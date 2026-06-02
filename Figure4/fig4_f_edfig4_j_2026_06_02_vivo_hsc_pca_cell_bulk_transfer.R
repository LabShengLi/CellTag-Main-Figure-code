#!/usr/bin/env Rscript
#
# vivo HSC PCA — INTEGRATED cell-level + sample-level (pseudo-bulk) analysis
#
# Merges logic from:
#   - 2026_05_31_vivo_hsc_pca_density_contour.R      (single-cell reference PCA + cell stats/plots)
#   - 2026_05_31_vivo_hsc_pca_density_contour_v2.R  (pseudo-bulk PERMANOVA / PERMDISP / replicate plot)
#
# Key design: bulk samples are NOT re-run through an independent prcomp PCA.
#   1. Build reference PCA on all HSC cells (Seurat NormalizeData -> HVG -> ScaleData -> RunPCA, 30 PCs).
#   2. Aggregate pseudo-bulk profiles per mouse_id (AverageExpression, RNA data slot).
#   3. Project pseudo-bulk into the cell PCA space via FindTransferAnchors +
#      IntegrateEmbeddings; save cell reference Seurat RDS + PC variance table.
#   4. Figures first (cell density contour + bulk replicate PCA), then statistical tests.
#   5. Cell- and bulk-level PERMANOVA / PERMDISP both use PC1–30 (same axes as cell reference).
#
# Data: CrossAge(exp2)_vivo.RDS; CloneID != 0; celltype_final == "HSC"
# Groups: orig.ident OO / OY / YO / YY; mouse_id = orig.ident-Rep (~12 pseudo-bulk profiles)

rm(list = ls())

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
    library(readr)
    library(Seurat)
    library(ggplot2)
    library(ggrepel)
    library(vegan)
})

# =============================================================================
# 0. Configuration
# =============================================================================
wdir <- "/project2/sli68423_1316/projects/U01_aim2/results/2026_06_02_pca_seperation"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)
setwd(wdir)

infn2 <- "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge(exp2)_vivo.RDS"

PERMANOVA_PERMUTATIONS <- 9999L
N_PCS_REFERENCE <- 30L          # cell reference PCA dimensions (fixed)
N_HVG <- 2000L
USER_SELECTED_N_PCS_BULK <- NA_integer_  # NA = PC1–N_PCS_REFERENCE (30); override e.g. 8 if needed

# =============================================================================
# 1. Shared helper functions
# =============================================================================
format_p_sci <- function(p, permutations = PERMANOVA_PERMUTATIONS) {
    p_floor <- 1 / (permutations + 1)
    vapply(
        p,
        function(x) {
            if (is.na(x)) {
                return(NA_character_)
            }
            if (x <= p_floor) {
                paste0("< ", formatC(p_floor, format = "e", digits = 2))
            } else {
                formatC(x, format = "e", digits = 2)
            }
        },
        character(1),
        USE.NAMES = FALSE
    )
}

add_permanova_p_sci_cols <- function(df, permutations = PERMANOVA_PERMUTATIONS) {
    out <- df
    if ("p_value" %in% names(out)) {
        out$p_value_sci <- format_p_sci(out$p_value, permutations = permutations)
    }
    if ("p_adj_BH" %in% names(out)) {
        out$p_adj_BH_sci <- format_p_sci(out$p_adj_BH, permutations = permutations)
    }
    out
}

get_col_or_na <- function(df, nm) {
    if (nm %in% colnames(df)) {
        return(df[[nm]])
    }
    rep(NA_real_, nrow(df))
}

extract_permdisp_group_dispersion_df <- function(permdisp_obj, group_order, test_label) {
    group_levels <- as.character(permdisp_obj$group)
    obs_distances <- as.numeric(permdisp_obj$distances)
    group_mean <- stats::setNames(
        as.numeric(permdisp_obj$group.distances),
        names(permdisp_obj$group.distances)
    )
    if (is.null(names(group_mean)) || any(names(group_mean) == "")) {
        names(group_mean) <- levels(permdisp_obj$group)
    }

    purrr::map_dfr(group_order, function(g) {
        idx <- group_levels == g
        dists <- obs_distances[idx]
        tibble::tibble(
            orig.ident = g,
            mean_distance_to_centroid = unname(group_mean[[g]]),
            median_distance_to_centroid = stats::median(dists, na.rm = TRUE),
            sd_distance_to_centroid = stats::sd(dists, na.rm = TRUE),
            n = sum(idx)
        )
    }) %>%
        dplyr::mutate(
            orig.ident = factor(orig.ident, levels = group_order),
            test = test_label
        )
}

run_permanova_global <- function(
    emb_mat,
    meta_df,
    permutations = PERMANOVA_PERMUTATIONS
) {
    adonis2(
        dist(emb_mat) ~ orig.ident,
        data = meta_df,
        permutations = permutations
    )
}

run_permanova_pairwise <- function(
    emb_mat,
    meta_df,
    group_order,
    permutations = PERMANOVA_PERMUTATIONS
) {
    pairs <- utils::combn(group_order, 2, simplify = FALSE)
    meta_group <- meta_df$orig.ident

    purrr::map_dfr(pairs, function(groups) {
        idx <- as.character(meta_group) %in% groups
        emb_sub <- emb_mat[idx, , drop = FALSE]
        meta_sub <- meta_df[idx, , drop = FALSE]
        meta_sub$orig.ident <- factor(meta_sub$orig.ident, levels = groups)
        fit <- adonis2(
            dist(emb_sub) ~ orig.ident,
            data = meta_sub,
            permutations = permutations
        )
        tibble::tibble(
            group1 = groups[1],
            group2 = groups[2],
            comparison = paste(groups, collapse = " vs "),
            permutations = permutations,
            p_value_floor = 1 / (permutations + 1),
            df = fit$Df[1],
            sumOfSqs = fit$SumOfSqs[1],
            R2 = fit$R2[1],
            F = fit$F[1],
            p_value = fit$`Pr(>F)`[1]
        )
    }) %>%
        dplyr::mutate(
            p_adj_BH = stats::p.adjust(p_value, method = "BH")
        ) %>%
        add_permanova_p_sci_cols() %>%
        dplyr::mutate(
            cross_recipient = comparison %in% c("OO vs OY", "YO vs YY")
        ) %>%
        dplyr::arrange(p_value)
}

build_permanova_asymmetry_tables <- function(permanova_pairwise_df) {
    permanova_asymmetry_df <- permanova_pairwise_df %>%
        dplyr::filter(comparison %in% c("OO vs OY", "YO vs YY")) %>%
        dplyr::select(comparison, R2, F, p_value, p_value_sci, p_adj_BH, p_adj_BH_sci)

    oo_oy_r2 <- permanova_asymmetry_df$R2[permanova_asymmetry_df$comparison == "OO vs OY"]
    yo_yy_r2 <- permanova_asymmetry_df$R2[permanova_asymmetry_df$comparison == "YO vs YY"]

    permanova_asymmetry_summary_df <- tibble::tibble(
        comparison_old_donor = "OO vs OY",
        comparison_young_donor = "YO vs YY",
        R2_OO_vs_OY = oo_oy_r2,
        R2_YO_vs_YY = yo_yy_r2,
        p_value_OO_vs_OY = permanova_asymmetry_df$p_value[permanova_asymmetry_df$comparison == "OO vs OY"],
        p_value_YO_vs_YY = permanova_asymmetry_df$p_value[permanova_asymmetry_df$comparison == "YO vs YY"],
        p_adj_BH_OO_vs_OY = permanova_asymmetry_df$p_adj_BH[permanova_asymmetry_df$comparison == "OO vs OY"],
        p_adj_BH_YO_vs_YY = permanova_asymmetry_df$p_adj_BH[permanova_asymmetry_df$comparison == "YO vs YY"],
        p_value_sci_OO_vs_OY = permanova_asymmetry_df$p_value_sci[permanova_asymmetry_df$comparison == "OO vs OY"],
        p_value_sci_YO_vs_YY = permanova_asymmetry_df$p_value_sci[permanova_asymmetry_df$comparison == "YO vs YY"],
        p_adj_BH_sci_OO_vs_OY = permanova_asymmetry_df$p_adj_BH_sci[permanova_asymmetry_df$comparison == "OO vs OY"],
        p_adj_BH_sci_YO_vs_YY = permanova_asymmetry_df$p_adj_BH_sci[permanova_asymmetry_df$comparison == "YO vs YY"],
        R2_difference_OO_OY_minus_YO_YY = oo_oy_r2 - yo_yy_r2
    )

    list(
        asymmetry = permanova_asymmetry_df,
        asymmetry_summary = permanova_asymmetry_summary_df
    )
}

format_permdisp_global_df <- function(
    permdisp_global_test,
    test_label,
    extra_cols = list()
) {
    permdisp_global_df <- tibble::as_tibble(
        as.data.frame(permdisp_global_test$tab),
        rownames = "term"
    )
    colnames(permdisp_global_df) <- gsub("\\s+", ".", colnames(permdisp_global_df))

    if (!"Pr(>F)" %in% colnames(permdisp_global_df)) {
        p_col <- colnames(permdisp_global_df)[grepl("^Pr", colnames(permdisp_global_df))]
        if (length(p_col) > 0) {
            colnames(permdisp_global_df)[match(p_col[1], colnames(permdisp_global_df))] <- "Pr(>F)"
        }
    }

    permdisp_global_df %>%
        dplyr::mutate(
            df = get_col_or_na(., "Df"),
            sumOfSqs = dplyr::coalesce(get_col_or_na(., "Sum.Sq"), get_col_or_na(., "Sum.Sqs")),
            meanSqs = dplyr::coalesce(get_col_or_na(., "Mean.Sq"), get_col_or_na(., "Mean.Sqs")),
            F = get_col_or_na(., "F"),
            p_value = get_col_or_na(., "Pr(>F)")
        ) %>%
        dplyr::select(term, df, sumOfSqs, meanSqs, F, p_value) %>%
        dplyr::mutate(
            test = test_label,
            permutations = PERMANOVA_PERMUTATIONS,
            p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
            p_value_sci = format_p_sci(p_value),
            !!!extra_cols
        )
}

compute_centroid_distance <- function(emb_mat, meta_group, group_order) {
    centroid_list <- lapply(group_order, function(g) {
        cells_use <- names(meta_group)[meta_group == g]
        colMeans(emb_mat[cells_use, , drop = FALSE])
    })
    names(centroid_list) <- group_order

    centroid_mat <- do.call(rbind, centroid_list)
    dist_mat <- as.matrix(stats::dist(centroid_mat, method = "euclidean"))
    dist_mat <- dist_mat[group_order, group_order, drop = FALSE]

    dist_long <- as.data.frame(as.table(dist_mat))
    colnames(dist_long) <- c("group1", "group2", "euclidean_distance")
    dist_long <- dist_long[
        match(dist_long$group1, group_order) < match(dist_long$group2, group_order),
    ]

    list(
        centroid_mat = centroid_mat,
        dist_mat = dist_mat,
        dist_long = dist_long
    )
}

write_distance_matrix_tsv <- function(dist_mat, out_path) {
    out_df <- as.data.frame(dist_mat)
    out_df <- tibble::rownames_to_column(out_df, var = "group")
    readr::write_tsv(out_df, out_path)
}

# =============================================================================
# 2. Load data and define metadata (shared by cell + bulk)
# =============================================================================
message("=== Section 2: Load HSC cells and define metadata ===")

ss <- readRDS(infn2)
ss <- subset(ss, subset = CloneID != 0)
ss_hsc <- subset(ss, subset = celltype_final == "HSC")

group_order <- c("OO", "OY", "YO", "YY")
group_order <- group_order[group_order %in% unique(ss_hsc$orig.ident)]

ss_hsc$orig.ident <- as.character(ss_hsc$orig.ident)
ss_hsc$Rep <- as.character(ss_hsc$Rep)
ss_hsc$mouse_id <- gsub(
    "_", "-",
    paste(ss_hsc$orig.ident, ss_hsc$Rep, sep = "-"),
    fixed = TRUE
)
ss_hsc$donor_age <- ifelse(substr(ss_hsc$orig.ident, 1, 1) == "O", "Old", "Young")
ss_hsc$host_age <- ifelse(substr(ss_hsc$orig.ident, 2, 2) == "O", "Old", "Young")

sample_meta <- ss_hsc[[]] %>%
    dplyr::select(mouse_id, orig.ident, Rep, donor_age, host_age) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        orig.ident = factor(orig.ident, levels = group_order),
        donor_age = factor(donor_age, levels = c("Young", "Old")),
        host_age = factor(host_age, levels = c("Young", "Old"))
    ) %>%
    dplyr::arrange(orig.ident, Rep)

# =============================================================================
# 3. Cell-level reference PCA (Seurat; defines shared PC axes)
# =============================================================================
message("=== Section 3: Cell-level reference PCA (", N_PCS_REFERENCE, " PCs) ===")

DefaultAssay(ss_hsc) <- "RNA"
ss_hsc <- NormalizeData(ss_hsc, verbose = FALSE)
ss_hsc <- FindVariableFeatures(
    ss_hsc,
    selection.method = "vst",
    nfeatures = N_HVG,
    verbose = FALSE
)
ss_hsc <- ScaleData(ss_hsc, features = VariableFeatures(ss_hsc), verbose = FALSE)
ss_hsc <- RunPCA(
    ss_hsc,
    features = VariableFeatures(ss_hsc),
    npcs = N_PCS_REFERENCE,
    verbose = FALSE
)

pca_emb <- Embeddings(ss_hsc, reduction = "pca")
pcs_use <- seq_len(N_PCS_REFERENCE)
pcs_use <- pcs_use[pcs_use <= ncol(pca_emb)]
pca_emb <- pca_emb[, pcs_use, drop = FALSE]

ref_pca_stdev <- ss_hsc[["pca"]]@stdev
ref_var_explained <- (ref_pca_stdev^2) / sum(ref_pca_stdev^2)
ref_pca_variance_df <- tibble::tibble(
    PC = seq_along(ref_pca_stdev),
    variance_explained_reference_cells = ref_var_explained,
    cumulative_variance_explained_reference_cells = cumsum(ref_var_explained)
)

message(
    "Reference PCA complete: ",
    nrow(pca_emb), " cells x ",
    ncol(pca_emb), " PCs."
)

saveRDS(
    ss_hsc,
    file.path(wdir, "HSC_cell_reference_pca_seurat_bulk_transfer.RDS")
)

ref_pca_variance_out_df <- ref_pca_variance_df %>%
    dplyr::mutate(
        stdev = ref_pca_stdev,
        variance_explained_pct = 100 * variance_explained_reference_cells,
        cumulative_variance_explained_pct = 100 * cumulative_variance_explained_reference_cells,
        n_cells = nrow(pca_emb),
        n_hvg = length(VariableFeatures(ss_hsc)),
        npcs = N_PCS_REFERENCE
    )

write_tsv(
    ref_pca_variance_out_df,
    file.path(wdir, "HSC_cell_reference_PCA_variance_explained.tsv")
)

cell_pca_coords_df <- tibble::rownames_to_column(
    as.data.frame(pca_emb),
    var = "cell_id"
) %>%
    dplyr::left_join(
        ss_hsc[[]] %>%
            tibble::rownames_to_column("cell_id") %>%
            dplyr::select(cell_id, orig.ident, Rep, mouse_id, donor_age, host_age),
        by = "cell_id"
    ) %>%
    dplyr::select(
        cell_id,
        orig.ident,
        Rep,
        mouse_id,
        donor_age,
        host_age,
        dplyr::starts_with("PC")
    )

write_tsv(
    cell_pca_coords_df,
    file.path(wdir, "HSC_cell_PCA_coordinates_all_PCs_bulk_transfer.tsv")
)

message(
    "Saved cell reference Seurat object: HSC_cell_reference_pca_seurat_bulk_transfer.RDS; ",
    "PC variance: HSC_cell_reference_PCA_variance_explained.tsv; ",
    "cell PC coordinates (PC1–", ncol(pca_emb), "): HSC_cell_PCA_coordinates_all_PCs_bulk_transfer.tsv."
)

# =============================================================================
# 4. Pseudo-bulk expression + projection into reference PCA (FindTransferAnchors)
# =============================================================================
message("=== Section 4: Pseudo-bulk aggregation and PCA transfer (no independent prcomp) ===")

pb_expr <- Seurat::AverageExpression(
    ss_hsc,
    assays = "RNA",
    slot = "data",
    group.by = "mouse_id",
    verbose = FALSE
)$RNA

sample_meta$mouse_id <- gsub("_", "-", sample_meta$mouse_id, fixed = TRUE)
sample_meta <- sample_meta %>%
    dplyr::filter(mouse_id %in% colnames(pb_expr))
pb_expr <- pb_expr[, sample_meta$mouse_id, drop = FALSE]

# Query Seurat object: one pseudo-bulk profile per mouse replicate
pb_seu <- CreateSeuratObject(counts = pb_expr, assay = "RNA")
pb_seu$orig.ident <- sample_meta$orig.ident[match(colnames(pb_seu), sample_meta$mouse_id)]
pb_seu$mouse_id <- colnames(pb_seu)
pb_seu <- NormalizeData(pb_seu, verbose = FALSE)

transfer_dims <- seq_len(min(N_PCS_REFERENCE, length(ref_pca_stdev)))
transfer_features <- VariableFeatures(ss_hsc)
transfer_features <- transfer_features[transfer_features %in% rownames(pb_seu)]

if (length(transfer_features) < 50) {
    stop("Too few shared variable features between reference cells and pseudo-bulk for transfer.")
}

# Seurat defaults k.score = 30, which fails when query has only ~12 pseudo-bulk samples.
# Both k.score and k.anchor must be strictly less than ncol(query) and ncol(reference).
n_ref_cells <- ncol(ss_hsc)
n_query_samples <- ncol(pb_seu)
k_anchor_transfer <- max(5L, min(5L, n_query_samples - 1L, n_ref_cells - 1L))
k_score_transfer <- max(5L, min(30L, n_query_samples - 1L, n_ref_cells - 1L))

message(
    "FindTransferAnchors: reference = ", n_ref_cells, " cells, query = ",
    n_query_samples, " pseudo-bulk samples, dims = ",
    min(transfer_dims), ":", max(transfer_dims),
    ", k.anchor = ", k_anchor_transfer, ", k.score = ", k_score_transfer, "."
)

transfer_anchors <- FindTransferAnchors(
    reference = ss_hsc,
    query = pb_seu,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = transfer_dims,
    features = transfer_features,
    k.anchor = k_anchor_transfer,
    k.score = k_score_transfer,
    verbose = FALSE
)

# Project query into reference PCA space (same step MapQuery runs for reduction.model = "pca").
# MapQuery is avoided: it also runs TransferData with default k.weight = 100, which fails
# when ncol(query) is small (~12 pseudo-bulk samples). k.weight = FALSE disables neighbor
# weighting and is appropriate for sample-level query objects.
bulk_pca_reduction <- "ref.pca"

message(
    "IntegrateEmbeddings: new.reduction.name = ", bulk_pca_reduction,
    ", dims = ", min(transfer_dims), ":", max(transfer_dims),
    ", k.weight = FALSE (", n_query_samples, " query samples)."
)

pb_seu <- IntegrateEmbeddings(
    anchorset = transfer_anchors,
    reference = ss_hsc,
    query = pb_seu,
    new.reduction.name = bulk_pca_reduction,
    dims.to.integrate = transfer_dims,
    k.weight = FALSE,
    verbose = FALSE
)

if (!bulk_pca_reduction %in% names(pb_seu@reductions)) {
    stop("IntegrateEmbeddings did not create reduction: ", bulk_pca_reduction)
}

bulk_pca_all <- Embeddings(pb_seu, reduction = bulk_pca_reduction)
rownames(bulk_pca_all) <- colnames(pb_seu)

if (nrow(bulk_pca_all) < 3L || ncol(bulk_pca_all) < 2L) {
    stop("Too few pseudo-bulk samples or transferred PCs for PERMANOVA.")
}

# Bulk stats use the same PC1–30 axes as the cell reference (not capped at n − 1 samples).
if (is.na(USER_SELECTED_N_PCS_BULK)) {
    pc_final_bulk <- min(N_PCS_REFERENCE, ncol(bulk_pca_all))
    pc_source_bulk <- "reference_matched_PC1_30"
} else {
    pc_final_bulk <- max(
        2L,
        min(as.integer(USER_SELECTED_N_PCS_BULK), N_PCS_REFERENCE, ncol(bulk_pca_all))
    )
    pc_source_bulk <- "user_selected"
}

sample_pcs <- bulk_pca_all[, seq_len(pc_final_bulk), drop = FALSE]
rownames(sample_pcs) <- rownames(bulk_pca_all)

message(
    "Transferred pseudo-bulk PCA: reduction = ", bulk_pca_reduction,
    "; bulk PERMANOVA/PERMDISP use PC1–", ncol(sample_pcs),
    " (same as cell reference; n_samples = ", nrow(sample_pcs),
    ", source = ", pc_source_bulk, ")."
)

# Variance across pseudo-bulk samples in transferred PC space (descriptive only)
bulk_var_explained <- apply(bulk_pca_all, 2, stats::var, na.rm = TRUE)
bulk_var_explained <- bulk_var_explained / sum(bulk_var_explained, na.rm = TRUE)
bulk_pca_variance_df <- tibble::tibble(
    PC = seq_along(bulk_var_explained),
    variance_explained_across_samples = bulk_var_explained,
    cumulative_variance_explained_across_samples = cumsum(bulk_var_explained),
    pca_method = "FindTransferAnchors_IntegrateEmbeddings_from_cell_reference",
    reference_n_pcs = N_PCS_REFERENCE,
    n_pcs_used_for_permanova = ncol(sample_pcs)
) %>%
    dplyr::left_join(
        ref_pca_variance_df %>% dplyr::rename_with(~paste0(.x, "_ref_annotation"), -PC),
        by = "PC"
    )

# =============================================================================
# 5. Figures — cell density contour + bulk replicate PCA (before slow permutation tests)
# =============================================================================
message("=== Section 5: Figures (cell + bulk PCA plots) ===")

sample_meta_bulk <- sample_meta %>%
    dplyr::mutate(orig.ident = factor(orig.ident, levels = group_order))

# --- Plot helpers ---
extract_group_colors <- function(p, groups) {
    built <- ggplot_build(p)
    colour_scale <- Filter(function(s) "colour" %in% s$aesthetics, built$plot$scales$scales)[[1]]
    setNames(colour_scale$map(groups), groups)
}

get_contour_axis_limits <- function(plot_df, bins = 5, expand = 0.05) {
    p_tmp <- ggplot(plot_df, aes(x = PC1, y = PC2, color = orig.ident)) +
        geom_density_2d(bins = bins, linewidth = 0.9, alpha = 0.9)

    built <- ggplot_build(p_tmp)
    contour_df <- purrr::detect(built$data, ~all(c("x", "y", "level") %in% names(.x)))

    if (is.null(contour_df) || nrow(contour_df) == 0) {
        return(list(
            xlim = range(plot_df$PC1, na.rm = TRUE),
            ylim = range(plot_df$PC2, na.rm = TRUE)
        ))
    }

    xrng <- range(contour_df$x, na.rm = TRUE)
    yrng <- range(contour_df$y, na.rm = TRUE)
    xpad <- diff(xrng) * expand
    ypad <- diff(yrng) * expand

    list(
        xlim = c(xrng[1] - xpad, xrng[2] + xpad),
        ylim = c(yrng[1] - ypad, yrng[2] + ypad)
    )
}

compute_density_peak_pc12 <- function(
    plot_df,
    group_order,
    n = 100,
    search_lims = NULL
) {
    plot_df$orig.ident <- factor(plot_df$orig.ident, levels = group_order)

    p_kde <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = orig.ident)) +
        ggplot2::stat_density_2d(geom = "tile", contour = FALSE, n = n)

    grid <- ggplot2::ggplot_build(p_kde)$data[[1]]
    if (!"density" %in% names(grid)) {
        stop("Could not extract KDE grid from stat_density_2d.")
    }

    group_ids <- sort(unique(grid$group))

    purrr::imap_dfr(group_order, function(g, i) {
        gi <- group_ids[i]
        dg <- grid[grid$group == gi, , drop = FALSE]

        if (!is.null(search_lims)) {
            dg <- dg[
                dg$x >= search_lims$xlim[1] &
                    dg$x <= search_lims$xlim[2] &
                    dg$y >= search_lims$ylim[1] &
                    dg$y <= search_lims$ylim[2],
                ,
                drop = FALSE
            ]
        }

        if (nrow(dg) == 0) {
            return(tibble::tibble(
                orig.ident = g,
                PC1 = NA_real_,
                PC2 = NA_real_,
                kde_density = NA_real_
            ))
        }

        idx <- which.max(dg$density)
        tibble::tibble(
            orig.ident = g,
            PC1 = dg$x[idx],
            PC2 = dg$y[idx],
            kde_density = dg$density[idx]
        )
    }) %>%
        dplyr::mutate(orig.ident = factor(orig.ident, levels = group_order))
}

compute_pc12_euclidean_dist_matrix <- function(marker_df, group_order, id_col = "orig.ident") {
    g <- as.character(group_order)
    ids <- as.character(marker_df[[id_col]])

    pc1 <- vapply(g, function(gi) {
        ii <- ids == gi
        if (!any(ii)) {
            return(NA_real_)
        }
        marker_df$PC1[ii][[1L]]
    }, numeric(1), USE.NAMES = FALSE)
    pc2 <- vapply(g, function(gi) {
        ii <- ids == gi
        if (!any(ii)) {
            return(NA_real_)
        }
        marker_df$PC2[ii][[1L]]
    }, numeric(1), USE.NAMES = FALSE)

    n <- length(g)
    dist_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(g, g))
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            dist_mat[i, j] <- sqrt((pc1[i] - pc1[j])^2 + (pc2[i] - pc2[j])^2)
        }
    }
    dist_mat
}

add_group_label_nudge <- function(marker_df, contour_lim) {
    xspan <- diff(contour_lim$xlim)
    yspan <- diff(contour_lim$ylim)

    marker_df %>%
        dplyr::mutate(
            nudge_x = dplyr::case_when(
                as.character(orig.ident) == "YO" ~ -0.12 * xspan,
                as.character(orig.ident) == "YY" ~ 0.12 * xspan,
                TRUE ~ 0
            ),
            nudge_y = dplyr::case_when(
                as.character(orig.ident) == "YO" ~ 0.14 * yspan,
                as.character(orig.ident) == "YY" ~ -0.14 * yspan,
                TRUE ~ 0
            )
        )
}

make_density_contour_marker_plot <- function(
    plot_df,
    marker_df,
    title,
    group_colors,
    donor_shape_values,
    legend_fill_values,
    contour_lim
) {
    marker_label_df <- add_group_label_nudge(marker_df, contour_lim)

    ggplot(plot_df, aes(x = PC1, y = PC2, color = orig.ident)) +
        geom_density_2d(linewidth = 0.9, alpha = 0.9, bins = 5, show.legend = FALSE) +
        scale_color_manual(
            values = group_colors,
            limits = names(group_colors),
            drop = FALSE
        ) +
        geom_point(
            data = marker_df,
            aes(x = PC1, y = PC2, color = orig.ident, fill = orig.ident, shape = donor),
            inherit.aes = FALSE,
            size = 5,
            stroke = 0.6
        ) +
        scale_shape_manual(
            values = donor_shape_values,
            limits = names(donor_shape_values),
            drop = FALSE,
            name = "Donor"
        ) +
        scale_fill_manual(
            values = group_colors,
            limits = names(group_colors),
            drop = FALSE
        ) +
        geom_label_repel(
            data = marker_label_df,
            aes(
                x = PC1,
                y = PC2,
                label = orig.ident,
                nudge_x = nudge_x,
                nudge_y = nudge_y
            ),
            inherit.aes = FALSE,
            colour = "black",
            fill = "white",
            size = 5,
            label.size = 0.35,
            label.padding = grid::unit(0.2, "lines"),
            show.legend = FALSE,
            box.padding = 0.4,
            point.padding = 1.1,
            min.segment.length = 0,
            force = 2,
            max.overlaps = Inf
        ) +
        coord_cartesian(xlim = contour_lim$xlim, ylim = contour_lim$ylim) +
        guides(
            color = guide_legend(
                title = "Group",
                order = 1,
                override.aes = list(
                    shape = 16,
                    fill = legend_fill_values,
                    colour = legend_fill_values,
                    linetype = "blank",
                    linewidth = 0,
                    stroke = 0.6,
                    size = 4
                )
            ),
            shape = guide_legend(
                title = "Donor",
                order = 2,
                override.aes = list(
                    colour = "black",
                    fill = "black",
                    size = 4,
                    stroke = 0.6
                )
            ),
            fill = "none"
        ) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14)
        ) +
        labs(
            title = title,
            x = "PC1",
            y = "PC2",
            color = "Group",
            shape = "Donor"
        )
}

pca_plot_df <- data.frame(
    PC1 = pca_emb[, 1],
    PC2 = pca_emb[, 2],
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    row.names = rownames(pca_emb)
)

donor_shape_values <- c(Old = 17, Young = 16)

p_ref_cell <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = orig.ident)) +
    geom_point(size = 1.5, alpha = 0.8) +
    theme_classic()

group_colors <- extract_group_colors(p_ref_cell, group_order)
contour_lim <- get_contour_axis_limits(pca_plot_df, bins = 5)
legend_fill_values <- unname(group_colors[group_order])

pca_density_peak_df <- compute_density_peak_pc12(
    pca_plot_df,
    group_order,
    n = 100,
    search_lims = contour_lim
) %>%
    dplyr::mutate(
        donor = factor(
            dplyr::if_else(substr(as.character(orig.ident), 1, 1) == "O", "Old", "Young"),
            levels = c("Old", "Young")
        )
    )

pca_density_peak_dist_mat <- compute_pc12_euclidean_dist_matrix(
    pca_density_peak_df,
    group_order = group_order
)

p_peak <- make_density_contour_marker_plot(
    plot_df = pca_plot_df,
    marker_df = pca_density_peak_df,
    title = "HSC 2D density contours and peaks (cell reference PCA)",
    group_colors = group_colors,
    donor_shape_values = donor_shape_values,
    legend_fill_values = legend_fill_values,
    contour_lim = contour_lim
)

print(p_peak)

make_replicate_label <- function(group, rep) {
    paste0(group, "_rep", tolower(sub("^Rep", "", rep)))
}

sample_plot_df <- sample_meta_bulk %>%
    dplyr::mutate(
        PC1 = bulk_pca_all[mouse_id, 1],
        PC2 = bulk_pca_all[mouse_id, 2],
        replicate_id = make_replicate_label(as.character(orig.ident), Rep),
        donor = factor(donor_age, levels = c("Old", "Young"))
    )

make_sample_pca_plot <- function(
    plot_df,
    group_colors,
    donor_shape_values,
    legend_fill_values
) {
    ggplot(plot_df, aes(x = PC1, y = PC2, color = orig.ident, shape = donor)) +
        geom_point(size = 5, stroke = 0.8) +
        scale_color_manual(
            values = group_colors,
            limits = names(group_colors),
            drop = FALSE
        ) +
        scale_shape_manual(
            values = donor_shape_values,
            limits = names(donor_shape_values),
            drop = FALSE,
            name = "Donor"
        ) +
        geom_text_repel(
            aes(label = replicate_id),
            colour = "black",
            size = 3.5,
            show.legend = FALSE,
            box.padding = 0.4,
            point.padding = 0.8,
            min.segment.length = Inf,
            segment.color = NA,
            force = 2,
            max.overlaps = Inf
        ) +
        guides(
            color = guide_legend(
                title = "Group",
                order = 1,
                override.aes = list(
                    shape = 16,
                    colour = legend_fill_values,
                    linetype = "blank",
                    linewidth = 0,
                    stroke = 0.6,
                    size = 4
                )
            ),
            shape = guide_legend(
                title = "Donor",
                order = 2,
                override.aes = list(
                    colour = "black",
                    fill = "black",
                    size = 4,
                    stroke = 0.6
                )
            )
        ) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14)
        ) +
        labs(
            title = "HSC pseudo-bulk replicate PCA (transferred from cell reference)",
            x = "PC1",
            y = "PC2",
            color = "Group",
            shape = "Donor"
        )
}

p_sample_pca <- make_sample_pca_plot(
    sample_plot_df,
    group_colors = group_colors,
    donor_shape_values = donor_shape_values,
    legend_fill_values = legend_fill_values
)

print(p_sample_pca)

sample_pca_coords_df <- tibble::rownames_to_column(
    as.data.frame(bulk_pca_all),
    var = "mouse_id"
) %>%
    dplyr::left_join(
        sample_meta %>%
            dplyr::mutate(
                replicate_id = make_replicate_label(
                    as.character(orig.ident),
                    Rep
                )
            ),
        by = "mouse_id"
    ) %>%
    dplyr::select(
        mouse_id,
        replicate_id,
        orig.ident,
        Rep,
        donor_age,
        host_age,
        dplyr::starts_with("PC")
    )

pb_expr_df <- pb_expr %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene")

ggsave(
    file.path(wdir, "HSC_density_contours_density_peak_PC1_PC2_PCA_all_groups.pdf"),
    p_peak,
    width = 5,
    height = 4
)
ggsave(
    file.path(wdir, "HSC_pseudobulk_replicate_PC1_PC2_PCA_bulk_transfer.pdf"),
    p_sample_pca,
    width = 6,
    height = 5
)

write_tsv(pca_density_peak_df, file.path(wdir, "HSC_density_peak_PC1_PC2_PCA_all_groups.tsv"))
write_distance_matrix_tsv(
    pca_density_peak_dist_mat,
    file.path(wdir, "HSC_density_peak_PC1_PC2_euclidean_distance_matrix.tsv")
)
write_tsv(
    sample_plot_df %>%
        dplyr::select(
            mouse_id, replicate_id, orig.ident, Rep, donor, host_age,
            PC1, PC2
        ),
    file.path(wdir, "HSC_pseudobulk_replicate_PC1_PC2_coordinates_bulk_transfer.tsv")
)
write_tsv(
    sample_pca_coords_df,
    file.path(wdir, "HSC_pseudobulk_replicate_PCA_coordinates_all_PCs_bulk_transfer.tsv")
)
write_tsv(pb_expr_df, file.path(wdir, "HSC_pseudobulk_expression_matrix_bulk_transfer.tsv"))
write_tsv(
    bulk_pca_variance_df,
    file.path(wdir, "HSC_pseudobulk_PCA_variance_explained_bulk_transfer.tsv")
)

message(
    "Section 5 complete: cell + bulk PCA PDFs and coordinate tables written. ",
    "Starting permutation-based statistical tests..."
)

# =============================================================================
# 6. Cell-level statistical tests (PERMANOVA, PERMDISP, centroids)
# =============================================================================
message("=== Section 6: Cell-level PERMANOVA / PERMDISP / centroids ===")

cell_meta_df <- data.frame(
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    stringsAsFactors = FALSE
)

permanova_global_cell <- run_permanova_global(
    emb_mat = pca_emb,
    meta_df = cell_meta_df
)

permanova_global_cell_df <- tibble::tibble(
    test = "PERMANOVA global (cell-level reference PCA)",
    n_pcs = ncol(pca_emb),
    n_cells = nrow(pca_emb),
    permutations = PERMANOVA_PERMUTATIONS,
    p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
    df = permanova_global_cell$Df[1],
    sumOfSqs = permanova_global_cell$SumOfSqs[1],
    R2 = permanova_global_cell$R2[1],
    F = permanova_global_cell$F[1],
    p_value = permanova_global_cell$`Pr(>F)`[1]
) %>%
    add_permanova_p_sci_cols()

permanova_global_cell_full_df <- tibble::as_tibble(permanova_global_cell, rownames = "term") %>%
    dplyr::mutate(`Pr(>F)_sci` = format_p_sci(`Pr(>F)`))

permanova_pairwise_cell_df <- run_permanova_pairwise(
    emb_mat = pca_emb,
    meta_df = cell_meta_df,
    group_order = group_order
)

cell_asymmetry <- build_permanova_asymmetry_tables(permanova_pairwise_cell_df)

cell_meta_interaction <- data.frame(
    mouse_id = ss_hsc$mouse_id,
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    donor_age = factor(ss_hsc$donor_age, levels = c("Young", "Old")),
    host_age = factor(ss_hsc$host_age, levels = c("Young", "Old"))
)
cell_dist <- dist(pca_emb)

permdisp_global_cell <- vegan::betadisper(cell_dist, cell_meta_interaction$orig.ident)
permdisp_global_cell_test <- vegan::permutest(
    permdisp_global_cell,
    permutations = PERMANOVA_PERMUTATIONS
)

permdisp_global_cell_df <- format_permdisp_global_df(
    permdisp_global_cell_test,
    test_label = "PERMDISP global (cell-level reference PCA, orig.ident)",
    extra_cols = list(
        n_pcs = ncol(pca_emb),
        n_cells = nrow(pca_emb)
    )
)

permdisp_group_cell_df <- extract_permdisp_group_dispersion_df(
    permdisp_global_cell,
    group_order = group_order,
    test_label = "PERMDISP group dispersion (cell-level reference PCA, orig.ident)"
) %>%
    dplyr::mutate(
        n_pcs = ncol(pca_emb),
        permutations = PERMANOVA_PERMUTATIONS
    )

interaction_note_cell <- paste(
    "Cell-level interaction without strata (donor/host invariant within mouse_id).",
    "Formal replicate-level asymmetry: see bulk-transfer interaction table."
)

permanova_interaction_cell <- adonis2(
    cell_dist ~ donor_age * host_age,
    data = cell_meta_interaction,
    permutations = PERMANOVA_PERMUTATIONS,
    by = "margin"
)

permanova_interaction_cell_df <- tibble::as_tibble(permanova_interaction_cell, rownames = "term") %>%
    dplyr::mutate(
        test = "PERMANOVA donor_age * host_age (cell-level, margin, no strata)",
        n_pcs = ncol(pca_emb),
        n_cells = nrow(pca_emb),
        n_mice = length(unique(cell_meta_interaction$mouse_id)),
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(`Pr(>F)`),
        note = interaction_note_cell
    )

pca_centroid_res <- compute_centroid_distance(
    emb_mat = pca_emb,
    meta_group = ss_hsc$orig.ident,
    group_order = group_order
)

pca_group_centroid_df <- tibble::rownames_to_column(
    as.data.frame(pca_centroid_res$centroid_mat),
    var = "orig.ident"
) %>%
    dplyr::mutate(orig.ident = factor(orig.ident, levels = group_order))

# =============================================================================
# 7. Bulk-level statistical tests (transferred PCA)
# =============================================================================
message("=== Section 7: Bulk-level PERMANOVA / PERMDISP (transferred PCA) ===")

permanova_global_bulk <- run_permanova_global(
    emb_mat = sample_pcs,
    meta_df = sample_meta_bulk
)

permanova_global_bulk_df <- tibble::tibble(
    test = "PERMANOVA global (sample-level pseudo-bulk, transferred PCA)",
    pca_method = "FindTransferAnchors_IntegrateEmbeddings_from_cell_reference",
    bulk_pca_reduction = bulk_pca_reduction,
    n_pcs = ncol(sample_pcs),
    n_pcs_available = ncol(bulk_pca_all),
    n_pcs_reference = N_PCS_REFERENCE,
    n_pcs_user_setting = USER_SELECTED_N_PCS_BULK,
    n_pcs_selection_source = pc_source_bulk,
    n_samples = nrow(sample_pcs),
    permutations = PERMANOVA_PERMUTATIONS,
    p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
    df = permanova_global_bulk$Df[1],
    sumOfSqs = permanova_global_bulk$SumOfSqs[1],
    R2 = permanova_global_bulk$R2[1],
    F = permanova_global_bulk$F[1],
    p_value = permanova_global_bulk$`Pr(>F)`[1]
) %>%
    add_permanova_p_sci_cols()

permanova_global_bulk_full_df <- tibble::as_tibble(permanova_global_bulk, rownames = "term") %>%
    dplyr::mutate(`Pr(>F)_sci` = format_p_sci(`Pr(>F)`))

permanova_pairwise_bulk_df <- run_permanova_pairwise(
    emb_mat = sample_pcs,
    meta_df = sample_meta_bulk,
    group_order = group_order
)

bulk_asymmetry <- build_permanova_asymmetry_tables(permanova_pairwise_bulk_df)

sample_dist_bulk <- dist(sample_pcs)

permdisp_global_bulk <- vegan::betadisper(sample_dist_bulk, sample_meta_bulk$orig.ident)
permdisp_global_bulk_test <- vegan::permutest(
    permdisp_global_bulk,
    permutations = PERMANOVA_PERMUTATIONS
)

permdisp_global_bulk_df <- format_permdisp_global_df(
    permdisp_global_bulk_test,
    test_label = "PERMDISP global (sample-level transferred PCA, orig.ident)",
    extra_cols = list(
        n_pcs = ncol(sample_pcs),
        n_samples = nrow(sample_pcs)
    )
)

permdisp_group_bulk_df <- extract_permdisp_group_dispersion_df(
    permdisp_global_bulk,
    group_order = group_order,
    test_label = "PERMDISP group dispersion (sample-level transferred PCA, orig.ident)"
) %>%
    dplyr::mutate(
        n_pcs = ncol(sample_pcs),
        n_samples = nrow(sample_pcs),
        permutations = PERMANOVA_PERMUTATIONS
    )

if (any(duplicated(sample_meta_bulk$mouse_id))) {
    permanova_interaction_bulk <- adonis2(
        sample_dist_bulk ~ donor_age * host_age,
        data = sample_meta_bulk,
        permutations = PERMANOVA_PERMUTATIONS,
        strata = sample_meta_bulk$mouse_id,
        by = "margin"
    )
    interaction_note_bulk <- "Interaction PERMANOVA used strata = mouse_id."
} else {
    permanova_interaction_bulk <- adonis2(
        sample_dist_bulk ~ donor_age * host_age,
        data = sample_meta_bulk,
        permutations = PERMANOVA_PERMUTATIONS,
        by = "margin"
    )
    interaction_note_bulk <- "Each mouse_id appears once at sample-level; strata = mouse_id not applied."
}

permanova_interaction_bulk_df <- tibble::as_tibble(permanova_interaction_bulk, rownames = "term") %>%
    dplyr::mutate(
        test = "PERMANOVA donor_age * host_age (sample-level transferred PCA)",
        pca_method = "FindTransferAnchors_IntegrateEmbeddings_from_cell_reference",
        n_pcs = ncol(sample_pcs),
        n_samples = nrow(sample_pcs),
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(`Pr(>F)`),
        note = interaction_note_bulk
    )

# =============================================================================
# 8. Write statistical result tables
# =============================================================================
message("=== Section 8: Write statistical tables ===")

write_tsv(pca_group_centroid_df, file.path(wdir, "HSC_group_centroid_PC1_30_PCA.tsv"))
write_distance_matrix_tsv(
    pca_centroid_res$dist_mat,
    file.path(wdir, "HSC_group_centroid_PC1_30_euclidean_distance_matrix.tsv")
)
write_tsv(
    pca_centroid_res$dist_long,
    file.path(wdir, "HSC_group_centroid_PC1_30_euclidean_distance_pairs.tsv")
)
write_tsv(permanova_global_cell_df, file.path(wdir, "HSC_PERMANOVA_global_PC_scores.tsv"))
write_tsv(permanova_global_cell_full_df, file.path(wdir, "HSC_PERMANOVA_global_PC_scores_full.tsv"))
write_tsv(permanova_pairwise_cell_df, file.path(wdir, "HSC_PERMANOVA_pairwise_PC_scores_BH.tsv"))
write_tsv(cell_asymmetry$asymmetry, file.path(wdir, "HSC_PERMANOVA_cross_recipient_pairwise.tsv"))
write_tsv(
    cell_asymmetry$asymmetry_summary,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_asymmetry.tsv")
)
write_tsv(permdisp_global_cell_df, file.path(wdir, "HSC_PERMDISP_global_PC_scores.tsv"))
write_tsv(permdisp_group_cell_df, file.path(wdir, "HSC_PERMDISP_group_dispersion_PC_scores.tsv"))
write_tsv(
    permanova_interaction_cell_df,
    file.path(wdir, "HSC_PERMANOVA_donor_host_interaction_PC_scores.tsv")
)

write_tsv(
    permanova_global_bulk_df,
    file.path(wdir, "HSC_PERMANOVA_global_sample_level_PC_scores_bulk_transfer.tsv")
)
write_tsv(
    permanova_global_bulk_full_df,
    file.path(wdir, "HSC_PERMANOVA_global_sample_level_PC_scores_full_bulk_transfer.tsv")
)
write_tsv(
    permanova_pairwise_bulk_df,
    file.path(wdir, "HSC_PERMANOVA_pairwise_sample_level_PC_scores_BH_bulk_transfer.tsv")
)
write_tsv(
    bulk_asymmetry$asymmetry,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_pairwise_sample_level_bulk_transfer.tsv")
)
write_tsv(
    bulk_asymmetry$asymmetry_summary,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_asymmetry_sample_level_bulk_transfer.tsv")
)
write_tsv(
    permdisp_global_bulk_df,
    file.path(wdir, "HSC_PERMDISP_global_sample_level_bulk_transfer.tsv")
)
write_tsv(
    permdisp_group_bulk_df,
    file.path(wdir, "HSC_PERMDISP_group_dispersion_sample_level_bulk_transfer.tsv")
)
write_tsv(
    permanova_interaction_bulk_df,
    file.path(wdir, "HSC_PERMANOVA_donor_host_interaction_sample_level_bulk_transfer.tsv")
)

save.image("2026_06_02_vivo_hsc_pca_cell_bulk_transfer.RData")

message(
    "Finished integrated analysis:\n",
    "  Section 3: HSC_cell_reference_pca_seurat_bulk_transfer.RDS, PC variance + cell PC1–30 TSV.\n",
    "  Section 5: cell + bulk PCA figures and coordinate tables.\n",
    "  Section 6–7: PERMANOVA/PERMDISP (cell + bulk, ", PERMANOVA_PERMUTATIONS, " permutations).\n",
    "  Bulk PCA: FindTransferAnchors/IntegrateEmbeddings (", ncol(sample_pcs), " PCs)."
)
