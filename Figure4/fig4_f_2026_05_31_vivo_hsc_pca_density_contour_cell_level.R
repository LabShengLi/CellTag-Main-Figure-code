#!/usr/bin/env Rscript
`#
# vivo HSC PCA — SINGLE-CELL level (pair with *_v2.R for pseudo-bulk / sample level)
#
# Data: CrossAge(exp2)_vivo.RDS; CloneID != 0; celltype_final == "HSC"
# Groups: orig.ident OO / OY / YO / YY (4 transplant combos); Rep -> mouse_id (metadata)
#
# PCA: Seurat NormalizeData -> HVG (2000) -> ScaleData -> RunPCA on all HSC cells (30 PCs)
# Unit of analysis for PERMANOVA / PERMDISP: one row per cell (cell-level PC distance matrix)
# Cross-recipient asymmetry: adonis2(d ~ donor_age * host_age, by = margin) without strata
#   (strata = mouse_id invalid here: donor/host constant within mouse; see v2 for replicate-level test)
#
# Outputs (wdir):
#   - PERMANOVA global + pairwise (BH) on cell PC scores
#   - PERMANOVA donor_age * host_age interaction (cell-level, margin; no strata)
#   - PERMDISP global + per-group dispersion on cell PC scores (betadisper + permutest)
#   - Group centroids (colMeans over cells, PC1-30) + Euclidean distance matrix / pairs
#   - PC1/PC2 density contours; KDE density peaks + peak distance matrix (PC1-PC2 only)
#
# See: fig4_f_2026_05_31_vivo_hsc_pca_density_contour_v2_bulk_level.R for mouse-level pseudo-bulk analysis

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

# -------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------
wdir <- "figure4"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Input data (filename only; read before setwd(wdir))
# -------------------------------------------------------------------------
infn2 <- "CrossAge(exp2)_vivo.RDS"

ss <- readRDS(infn2)
setwd(wdir)
ss <- subset(ss, subset = CloneID != 0)
ss_hsc <- subset(ss, subset = celltype_final == "HSC")

group_order <- c("OO", "OY", "YO", "YY")
group_order <- group_order[group_order %in% unique(ss_hsc$orig.ident)]

ss_hsc$orig.ident <- as.character(ss_hsc$orig.ident)
ss_hsc$Rep <- as.character(ss_hsc$Rep)
ss_hsc$mouse_id <- gsub("_", "-", paste(ss_hsc$orig.ident, ss_hsc$Rep, sep = "-"), fixed = TRUE)
ss_hsc$donor_age <- ifelse(substr(ss_hsc$orig.ident, 1, 1) == "O", "Old", "Young")
ss_hsc$host_age <- ifelse(substr(ss_hsc$orig.ident, 2, 2) == "O", "Old", "Young")

# PERMANOVA / PERMDISP: cell-level; min resolvable p = 1 / (permutations + 1)
PERMANOVA_PERMUTATIONS <- 9999

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

# -------------------------------------------------------------------------
# PCA
# -------------------------------------------------------------------------
DefaultAssay(ss_hsc) <- "RNA"
ss_hsc <- NormalizeData(ss_hsc, verbose = FALSE)
ss_hsc <- FindVariableFeatures(ss_hsc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ss_hsc <- ScaleData(ss_hsc, features = VariableFeatures(ss_hsc), verbose = FALSE)
ss_hsc <- RunPCA(ss_hsc, features = VariableFeatures(ss_hsc), npcs = 30, verbose = FALSE)

pca_emb <- Embeddings(ss_hsc, reduction = "pca")
pcs_use <- 1:30
pcs_use <- pcs_use[pcs_use <= ncol(pca_emb)]
pca_emb <- pca_emb[, pcs_use, drop = FALSE]

# -------------------------------------------------------------------------
# PERMANOVA on cell-level PC scores (vegan::adonis2)
# -------------------------------------------------------------------------
run_permanova_global <- function(
    emb_mat,
    meta_group,
    group_order,
    permutations = PERMANOVA_PERMUTATIONS
) {
    meta_df <- data.frame(
        orig.ident = factor(meta_group, levels = group_order)
    )
    adonis2(
        dist(emb_mat) ~ orig.ident,
        data = meta_df,
        permutations = permutations
    )
}

run_permanova_pairwise <- function(
    emb_mat,
    meta_group,
    group_order,
    permutations = PERMANOVA_PERMUTATIONS
) {
    pairs <- utils::combn(group_order, 2, simplify = FALSE)

    purrr::map_dfr(pairs, function(groups) {
        idx <- meta_group %in% groups
        emb_sub <- emb_mat[idx, , drop = FALSE]
        meta_sub <- data.frame(
            orig.ident = factor(meta_group[idx], levels = groups)
        )
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

permanova_global <- run_permanova_global(
    emb_mat = pca_emb,
    meta_group = ss_hsc$orig.ident,
    group_order = group_order
)

permanova_global_df <- tibble::tibble(
    test = "PERMANOVA global",
    n_pcs = ncol(pca_emb),
    n_cells = nrow(pca_emb),
    permutations = PERMANOVA_PERMUTATIONS,
    p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
    df = permanova_global$Df[1],
    sumOfSqs = permanova_global$SumOfSqs[1],
    R2 = permanova_global$R2[1],
    F = permanova_global$F[1],
    p_value = permanova_global$`Pr(>F)`[1]
) %>%
    add_permanova_p_sci_cols()

permanova_global_full_df <- tibble::as_tibble(permanova_global, rownames = "term") %>%
    dplyr::mutate(
        `Pr(>F)_sci` = format_p_sci(`Pr(>F)`)
    )

permanova_pairwise_df <- run_permanova_pairwise(
    emb_mat = pca_emb,
    meta_group = ss_hsc$orig.ident,
    group_order = group_order
)

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

message("PERMANOVA global (", ncol(pca_emb), " PCs, cell-level):")
print(permanova_global)
print(permanova_global_df)

message("PERMANOVA pairwise comparisons (BH-adjusted):")
print(permanova_pairwise_df)

message("Cross-recipient asymmetry — descriptive pairwise R2 (OO vs OY vs YO vs YY):")
print(permanova_asymmetry_df)
print(permanova_asymmetry_summary_df)

message(
    "PERMANOVA permutations = ", PERMANOVA_PERMUTATIONS,
    " (vegan default = 999); min resolvable p = ",
    format_p_sci(0)[[1]]
)

# -------------------------------------------------------------------------
# Cell-level distance matrix + metadata (PERMDISP + donor×host interaction)
# -------------------------------------------------------------------------
cell_meta_interaction <- data.frame(
    mouse_id = ss_hsc$mouse_id,
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    donor_age = factor(ss_hsc$donor_age, levels = c("Young", "Old")),
    host_age = factor(ss_hsc$host_age, levels = c("Young", "Old"))
)
cell_dist <- dist(pca_emb)

# -------------------------------------------------------------------------
# PERMDISP on cell-level PC scores (vegan::betadisper + permutest)
# -------------------------------------------------------------------------
permdisp_global <- vegan::betadisper(cell_dist, cell_meta_interaction$orig.ident)
permdisp_global_test <- vegan::permutest(
    permdisp_global,
    permutations = PERMANOVA_PERMUTATIONS
)

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

permdisp_global_df <- permdisp_global_df %>%
    dplyr::mutate(
        df = get_col_or_na(., "Df"),
        sumOfSqs = dplyr::coalesce(get_col_or_na(., "Sum.Sq"), get_col_or_na(., "Sum.Sqs")),
        meanSqs = dplyr::coalesce(get_col_or_na(., "Mean.Sq"), get_col_or_na(., "Mean.Sqs")),
        F = get_col_or_na(., "F"),
        p_value = get_col_or_na(., "Pr(>F)")
    ) %>%
    dplyr::select(term, df, sumOfSqs, meanSqs, F, p_value) %>%
    dplyr::mutate(
        test = "PERMDISP global (cell-level PC scores, orig.ident)",
        n_pcs = ncol(pca_emb),
        n_cells = nrow(pca_emb),
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(p_value)
    )

permdisp_group_df <- extract_permdisp_group_dispersion_df(
    permdisp_global,
    group_order = group_order,
    test_label = "PERMDISP group dispersion (cell-level PC scores, orig.ident)"
) %>%
    dplyr::mutate(
        n_pcs = ncol(pca_emb),
        permutations = PERMANOVA_PERMUTATIONS
    )

message("PERMDISP global (cell-level, ", ncol(pca_emb), " PCs):")
print(permdisp_global_df)

message("PERMDISP per-group dispersion (mean distance to group centroid):")
print(permdisp_group_df)

# -------------------------------------------------------------------------
# donor_age * host_age interaction PERMANOVA (cell-level, no strata)
# strata = mouse_id omitted: donor/host labels are constant within each mouse, so
# within-stratum permutations cannot test interaction (Pr(>F) -> 1). Labels are
# permuted across all cells; pair with *_v2.R for replicate-level formal inference.
# -------------------------------------------------------------------------
permanova_interaction <- adonis2(
    cell_dist ~ donor_age * host_age,
    data = cell_meta_interaction,
    permutations = PERMANOVA_PERMUTATIONS,
    by = "margin"
)
interaction_note <- paste(
    "Cell-level interaction without strata (donor/host invariant within mouse_id).",
    "Formal replicate-level asymmetry: see HSC_PERMANOVA_donor_host_interaction_sample_level_v2.tsv."
)

permanova_interaction_df <- tibble::as_tibble(permanova_interaction, rownames = "term") %>%
    dplyr::mutate(
        test = "PERMANOVA donor_age * host_age (cell-level, margin, no strata)",
        n_pcs = ncol(pca_emb),
        n_cells = nrow(pca_emb),
        n_mice = length(unique(cell_meta_interaction$mouse_id)),
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(`Pr(>F)`),
        note = interaction_note
    )

message("PERMANOVA donor_age * host_age interaction (cell-level, no strata):")
print(permanova_interaction_df)

# -------------------------------------------------------------------------
# Group centroids in PC1–30 and Euclidean distances
# -------------------------------------------------------------------------
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

message("Group centroids in PCA space (PC1–", ncol(pca_emb), ", colMeans per group):")
print(pca_group_centroid_df)

message("Group centroid Euclidean distance matrix (PC1–", ncol(pca_emb), "):")
print(round(pca_centroid_res$dist_mat, 4))

message("Group centroid pairwise Euclidean distances (PC1–", ncol(pca_emb), "):")
print(pca_centroid_res$dist_long)

pca_plot_df <- data.frame(
    PC1 = pca_emb[, 1],
    PC2 = pca_emb[, 2],
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    row.names = rownames(pca_emb)
)

# Donor shapes for density peaks: Old = triangle, Young = circle
donor_shape_values <- c(Old = 17, Young = 16)

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

# Reference plot to recover default ggplot group colors (same as 2026_05_28)
p_ref <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = orig.ident)) +
    geom_point(size = 1.5, alpha = 0.8) +
    theme_classic()

group_colors <- extract_group_colors(p_ref, group_order)
contour_lim <- get_contour_axis_limits(pca_plot_df, bins = 5)

# Density peak = argmax(density) on ggplot2::stat_density_2d grid (same KDE as geom_density_2d)
# Per-group bandwidth; panel-scale lims; optional search within plot coord limits
compute_density_peak_pc12 <- function(
    plot_df,
    group_order,
    n = 100,
    search_lims = NULL
) {
    plot_df <- plot_df
    plot_df$orig.ident <- factor(plot_df$orig.ident, levels = group_order)

    p_kde <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = orig.ident)) +
        ggplot2::stat_density_2d(geom = "tile", contour = FALSE, n = n)

    grid <- ggplot2::ggplot_build(p_kde)$data[[1]]
    if (!"density" %in% names(grid)) {
        stop("Could not extract KDE grid from stat_density_2d.")
    }

    group_ids <- sort(unique(grid$group))
    if (length(group_ids) != length(group_order)) {
        warning(
            "KDE grid group count (", length(group_ids), ") != group_order length (",
            length(group_order), "); matching by index."
        )
    }

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

# Euclidean distance matrix (PC1, PC2) between group markers
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

write_distance_matrix_tsv <- function(dist_mat, out_path) {
    out_df <- as.data.frame(dist_mat)
    out_df <- tibble::rownames_to_column(out_df, var = "group")
    readr::write_tsv(out_df, out_path)
}

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

message("Density peaks (mode) from ggplot2 KDE grid:")
print(pca_density_peak_df)

pca_density_peak_dist_mat <- compute_pc12_euclidean_dist_matrix(
    pca_density_peak_df,
    group_order = group_order
)

message("Density peak Euclidean distance matrix (PC1, PC2):")
print(round(pca_density_peak_dist_mat, 4))

# Legend keys: Group colors; Donor shapes shown in separate legend
legend_fill_values <- unname(group_colors[group_order])

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

plot_args <- list(
    plot_df = pca_plot_df,
    group_colors = group_colors,
    donor_shape_values = donor_shape_values,
    legend_fill_values = legend_fill_values,
    contour_lim = contour_lim
)

p_peak <- do.call(
    make_density_contour_marker_plot,
    c(
        plot_args,
        list(
            marker_df = pca_density_peak_df,
            title = "HSC 2D density contours and peaks"
        )
    )
)

print(p_peak)

ggsave(
    file.path(wdir, "HSC_density_contours_density_peak_PC1_PC2_PCA_all_groups.pdf"),
    p_peak,
    width = 5,
    height = 4
)

write_tsv(
    pca_density_peak_df,
    file.path(wdir, "HSC_density_peak_PC1_PC2_PCA_all_groups.tsv")
)

write_distance_matrix_tsv(
    pca_density_peak_dist_mat,
    file.path(wdir, "HSC_density_peak_PC1_PC2_euclidean_distance_matrix.tsv")
)

write_tsv(
    pca_group_centroid_df,
    file.path(wdir, "HSC_group_centroid_PC1_30_PCA.tsv")
)

write_distance_matrix_tsv(
    pca_centroid_res$dist_mat,
    file.path(wdir, "HSC_group_centroid_PC1_30_euclidean_distance_matrix.tsv")
)

write_tsv(
    pca_centroid_res$dist_long,
    file.path(wdir, "HSC_group_centroid_PC1_30_euclidean_distance_pairs.tsv")
)

write_tsv(
    permanova_global_df,
    file.path(wdir, "HSC_PERMANOVA_global_PC_scores.tsv")
)

write_tsv(
    permanova_global_full_df,
    file.path(wdir, "HSC_PERMANOVA_global_PC_scores_full.tsv")
)

write_tsv(
    permanova_pairwise_df,
    file.path(wdir, "HSC_PERMANOVA_pairwise_PC_scores_BH.tsv")
)

write_tsv(
    permanova_asymmetry_df,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_pairwise.tsv")
)

write_tsv(
    permanova_asymmetry_summary_df,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_asymmetry.tsv")
)

write_tsv(
    permdisp_global_df,
    file.path(wdir, "HSC_PERMDISP_global_PC_scores.tsv")
)

write_tsv(
    permdisp_group_df,
    file.path(wdir, "HSC_PERMDISP_group_dispersion_PC_scores.tsv")
)

write_tsv(
    permanova_interaction_df,
    file.path(wdir, "HSC_PERMANOVA_donor_host_interaction_PC_scores.tsv")
)

save.image("2026_05_31_vivo_hsc_pca_density_contour.RData")

message(
    "Finished: density contour peak plot, peak distance matrix, ",
    "group centroids (PC1–30) + Euclidean distances, ",
    "PERMANOVA (global + pairwise BH + donor×host interaction), ",
    "and PERMDISP (global + per-group dispersion)."
)
