#!/usr/bin/env Rscript
#
# vivo HSC PCA — SAMPLE-LEVEL pseudo-bulk (pair with fig4_f_2026_05_31_vivo_hsc_pca_density_contour_cell_level.R)
#
# Data: CrossAge(exp2)_vivo.RDS; CloneID != 0; celltype_final == "HSC"
# Groups: orig.ident OO / OY / YO / YY; sample = mouse_id (orig.ident-Rep), ~12 pseudo-bulk profiles
#
# Expression: Seurat::AverageExpression(RNA, data slot, group.by = mouse_id)
# PCA: prcomp on scaled HVG pseudo-bulk matrix (samples x genes); PCs by 80% cum. variance
#      (override with USER_SELECTED_N_PCS) or all available if threshold not reached
#
# Unit of analysis for PERMANOVA / PERMDISP: one row per mouse replicate (sample-level distances)
#   - Global + pairwise PERMANOVA (BH); PERMDISP global + per-group dispersion; donor_age * host_age interaction
#   - Cross-recipient asymmetry summary (OO vs OY, YO vs YY)
# Plot: PC1 vs PC2 of pseudo-bulk replicates (not single cells)
#
# Outputs (wdir): tables suffixed *_v2.tsv; PDF HSC_pseudobulk_replicate_PC1_PC2_PCA_v2.pdf

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

# PERMANOVA / PERMDISP: sample-level; min resolvable p = 1 / (permutations + 1)
PERMANOVA_PERMUTATIONS <- 9999
PCA_VAR_THRESHOLD <- 0.80   # used when USER_SELECTED_N_PCS is NA
USER_SELECTED_N_PCS <- NA_integer_  # e.g. 3, 5, 8 to fix PC count; NA = use 80% variance rule

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

# -------------------------------------------------------------------------
# Sample-level pseudo-bulk PCA + PERMANOVA / PERMDISP
# -------------------------------------------------------------------------
ss_hsc$orig.ident <- as.character(ss_hsc$orig.ident)
ss_hsc$Rep <- as.character(ss_hsc$Rep)
ss_hsc$mouse_id <- paste(ss_hsc$orig.ident, ss_hsc$Rep, sep = "-")
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

pb_expr <- Seurat::AverageExpression(
    ss_hsc,
    assays = "RNA",
    slot = "data",
    group.by = "mouse_id",
    verbose = FALSE
)$RNA

# Ensure sample IDs use the same separator convention as pseudo-bulk output.
sample_meta$mouse_id <- gsub("_", "-", sample_meta$mouse_id, fixed = TRUE)

sample_meta <- sample_meta %>%
    dplyr::filter(mouse_id %in% colnames(pb_expr))
pb_expr <- pb_expr[, sample_meta$mouse_id, drop = FALSE]

DefaultAssay(ss_hsc) <- "RNA"
if (length(Seurat::VariableFeatures(ss_hsc)) == 0) {
    ss_hsc <- Seurat::NormalizeData(ss_hsc, verbose = FALSE)
    ss_hsc <- Seurat::FindVariableFeatures(
        ss_hsc,
        selection.method = "vst",
        nfeatures = 2000,
        verbose = FALSE
    )
}

var_feats <- VariableFeatures(ss_hsc)
var_feats <- var_feats[var_feats %in% rownames(pb_expr)]
if (length(var_feats) < 50) {
    stop("Too few variable features available for pseudo-bulk PCA.")
}

pb_var <- apply(pb_expr[var_feats, , drop = FALSE], 1, stats::var, na.rm = TRUE)
var_feats_nonzero <- names(pb_var)[is.finite(pb_var) & pb_var > 0]
if (length(var_feats_nonzero) < 2) {
    stop("Too few non-zero-variance pseudo-bulk features for PCA.")
}
if (length(var_feats_nonzero) < length(var_feats)) {
    message(
        "Filtered zero-variance pseudo-bulk features: ",
        length(var_feats) - length(var_feats_nonzero),
        " removed; ",
        length(var_feats_nonzero),
        " kept."
    )
}

pb_mat <- t(pb_expr[var_feats_nonzero, , drop = FALSE])
sample_pca <- prcomp(pb_mat, center = TRUE, scale. = TRUE)

pc_max_available <- min(ncol(sample_pca$x), nrow(sample_pca$x) - 1)
if (pc_max_available < 2) {
    stop("Too few pseudo-bulk samples for PCA/PERMANOVA.")
}
sample_pcs_all <- sample_pca$x[, seq_len(pc_max_available), drop = FALSE]

var_explained <- (sample_pca$sdev^2) / sum(sample_pca$sdev^2)
pca_variance_df <- tibble::tibble(
    PC = seq_along(var_explained),
    variance_explained = var_explained,
    cumulative_variance_explained = cumsum(var_explained)
)

pc_recommended <- which(pca_variance_df$cumulative_variance_explained >= PCA_VAR_THRESHOLD)[1]
if (is.na(pc_recommended)) {
    pc_recommended <- pc_max_available
}

if (is.na(USER_SELECTED_N_PCS)) {
    pc_final <- pc_recommended
    pc_source <- "recommended"
} else {
    pc_final <- as.integer(USER_SELECTED_N_PCS)
    pc_final <- max(2L, min(pc_final, pc_max_available))
    pc_source <- "user_selected"
}

sample_pcs <- sample_pca$x[, seq_len(pc_final), drop = FALSE]
pc_final_var <- pca_variance_df$cumulative_variance_explained[pc_final]
pc_recommended_var <- pca_variance_df$cumulative_variance_explained[pc_recommended]

message("Pseudo-bulk PCA variance explained by each available PC:")
print(
    pca_variance_df %>%
        dplyr::mutate(
            variance_explained_pct = 100 * variance_explained,
            cumulative_variance_explained_pct = 100 * cumulative_variance_explained
        )
)
message(
    "Recommended PCs for ", sprintf("%.0f", 100 * PCA_VAR_THRESHOLD),
    "% cumulative variance: ", pc_recommended,
    " (cumulative = ", sprintf("%.2f%%", 100 * pc_recommended_var), ")."
)
message(
    "Using ", pc_final, " PCs for PERMANOVA (source = ", pc_source,
    ", cumulative variance = ", sprintf("%.2f%%", 100 * pc_final_var), ")."
)

run_permanova_global <- function(
    emb_mat,
    sample_df,
    permutations = PERMANOVA_PERMUTATIONS
) {
    adonis2(
        dist(emb_mat) ~ orig.ident,
        data = sample_df,
        permutations = permutations
    )
}

run_permanova_pairwise <- function(
    emb_mat,
    sample_df,
    group_order,
    permutations = PERMANOVA_PERMUTATIONS
) {
    pairs <- utils::combn(group_order, 2, simplify = FALSE)

    purrr::map_dfr(pairs, function(groups) {
        idx <- as.character(sample_df$orig.ident) %in% groups
        emb_sub <- emb_mat[idx, , drop = FALSE]
        meta_sub <- sample_df[idx, , drop = FALSE]
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

permanova_global <- run_permanova_global(
    emb_mat = sample_pcs,
    sample_df = sample_meta
)

permanova_global_df <- tibble::tibble(
    test = "PERMANOVA global (sample-level pseudo-bulk)",
    n_pcs = ncol(sample_pcs),
    n_pcs_available = pc_max_available,
    n_pcs_recommended_80pct = pc_recommended,
    n_pcs_user_setting = USER_SELECTED_N_PCS,
    n_pcs_selection_source = pc_source,
    n_samples = nrow(sample_pcs),
    permutations = PERMANOVA_PERMUTATIONS,
    p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
    cumulative_variance_recommended = pc_recommended_var,
    cumulative_variance_used = pc_final_var,
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
    emb_mat = sample_pcs,
    sample_df = sample_meta,
    group_order = group_order
)

sample_dist <- dist(sample_pcs)

permdisp_global <- vegan::betadisper(sample_dist, sample_meta$orig.ident)
permdisp_global_test <- vegan::permutest(permdisp_global, permutations = PERMANOVA_PERMUTATIONS)
permdisp_global_df <- tibble::as_tibble(
    as.data.frame(permdisp_global_test$tab),
    rownames = "term"
)

# Standardize column names from vegan outputs across versions
colnames(permdisp_global_df) <- gsub("\\s+", ".", colnames(permdisp_global_df))

if (!"Pr(>F)" %in% colnames(permdisp_global_df)) {
    p_col <- colnames(permdisp_global_df)[grepl("^Pr", colnames(permdisp_global_df))]
    if (length(p_col) > 0) {
        colnames(permdisp_global_df)[match(p_col[1], colnames(permdisp_global_df))] <- "Pr(>F)"
    }
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
        test = "PERMDISP global (orig.ident)",
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(p_value)
    )

permdisp_group_df <- extract_permdisp_group_dispersion_df(
    permdisp_global,
    group_order = group_order,
    test_label = "PERMDISP group dispersion (sample-level pseudo-bulk, orig.ident)"
) %>%
    dplyr::mutate(
        n_pcs = ncol(sample_pcs),
        n_samples = nrow(sample_pcs),
        permutations = PERMANOVA_PERMUTATIONS
    )

sample_meta_interaction <- sample_meta
if (any(duplicated(sample_meta_interaction$mouse_id))) {
    permanova_interaction <- adonis2(
        sample_dist ~ donor_age * host_age,
        data = sample_meta_interaction,
        permutations = PERMANOVA_PERMUTATIONS,
        strata = sample_meta_interaction$mouse_id,
        by = "margin"
    )
    interaction_note <- "Interaction PERMANOVA used strata = mouse_id."
} else {
    permanova_interaction <- adonis2(
        sample_dist ~ donor_age * host_age,
        data = sample_meta_interaction,
        permutations = PERMANOVA_PERMUTATIONS,
        by = "margin"
    )
    interaction_note <- "Each mouse_id appears once at sample-level; strata = mouse_id not applied."
}

permanova_interaction_df <- tibble::as_tibble(permanova_interaction, rownames = "term") %>%
    dplyr::mutate(
        permutations = PERMANOVA_PERMUTATIONS,
        p_value_floor = 1 / (PERMANOVA_PERMUTATIONS + 1),
        p_value_sci = format_p_sci(`Pr(>F)`),
        note = interaction_note
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

message("PERMANOVA global (pseudo-bulk sample-level, ", ncol(sample_pcs), " PCs):")
print(permanova_global)
print(permanova_global_df)

message("PERMANOVA pairwise comparisons (BH-adjusted):")
print(permanova_pairwise_df)

message("PERMDISP global:")
print(permdisp_global_df)

message("PERMDISP per-group dispersion (mean distance to group centroid):")
print(permdisp_group_df)

message("PERMANOVA donor_age * host_age interaction:")
print(permanova_interaction_df)

message("Cross-recipient asymmetry (OO vs OY vs YO vs YY):")
print(permanova_asymmetry_df)
print(permanova_asymmetry_summary_df)

message(
    "PERMANOVA permutations = ", PERMANOVA_PERMUTATIONS,
    " (vegan default = 999); min resolvable p = ",
    format_p_sci(0)[[1]]
)

# -------------------------------------------------------------------------
# Plot pseudo-bulk replicate points on PC1-PC2 (12 samples)
# -------------------------------------------------------------------------
make_replicate_label <- function(group, rep) {
    paste0(group, "_rep", tolower(sub("^Rep", "", rep)))
}

sample_plot_df <- sample_meta %>%
    dplyr::mutate(
        PC1 = sample_pca$x[, 1],
        PC2 = sample_pca$x[, 2],
        replicate_id = make_replicate_label(as.character(orig.ident), Rep),
        donor = factor(donor_age, levels = c("Old", "Young"))
    )

donor_shape_values <- c(Old = 17, Young = 16)

extract_group_colors <- function(p, groups) {
    built <- ggplot_build(p)
    colour_scale <- Filter(function(s) "colour" %in% s$aesthetics, built$plot$scales$scales)[[1]]
    setNames(colour_scale$map(groups), groups)
}

p_ref <- ggplot(sample_plot_df, aes(x = PC1, y = PC2, color = orig.ident)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_classic()

group_colors <- extract_group_colors(p_ref, group_order)
legend_fill_values <- unname(group_colors[group_order])

pc1_pct <- round(100 * var_explained[1], 1)
pc2_pct <- round(100 * var_explained[2], 1)

make_sample_pca_plot <- function(plot_df, group_colors, donor_shape_values, legend_fill_values) {
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
            title = "HSC pseudo-bulk replicate PCA",
            x = paste0("PC1 (", pc1_pct, "%)"),
            y = paste0("PC2 (", pc2_pct, "%)"),
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

ggsave(
    file.path(wdir, "HSC_pseudobulk_replicate_PC1_PC2_PCA_v2.pdf"),
    p_sample_pca,
    width = 6,
    height = 5
)

write_tsv(
    sample_plot_df %>%
        dplyr::select(
            mouse_id, replicate_id, orig.ident, Rep, donor, host_age,
            PC1, PC2
        ),
    file.path(wdir, "HSC_pseudobulk_replicate_PC1_PC2_coordinates_v2.tsv")
)

write_tsv(
    permanova_global_df,
    file.path(wdir, "HSC_PERMANOVA_global_sample_level_PC_scores_v2.tsv")
)

write_tsv(
    permanova_global_full_df,
    file.path(wdir, "HSC_PERMANOVA_global_sample_level_PC_scores_full_v2.tsv")
)

write_tsv(
    permanova_pairwise_df,
    file.path(wdir, "HSC_PERMANOVA_pairwise_sample_level_PC_scores_BH_v2.tsv")
)

write_tsv(
    permanova_asymmetry_df,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_pairwise_sample_level_v2.tsv")
)

write_tsv(
    permanova_asymmetry_summary_df,
    file.path(wdir, "HSC_PERMANOVA_cross_recipient_asymmetry_sample_level_v2.tsv")
)

write_tsv(
    permdisp_global_df,
    file.path(wdir, "HSC_PERMDISP_global_sample_level_v2.tsv")
)

write_tsv(
    permdisp_group_df,
    file.path(wdir, "HSC_PERMDISP_group_dispersion_sample_level_v2.tsv")
)

write_tsv(
    permanova_interaction_df,
    file.path(wdir, "HSC_PERMANOVA_donor_host_interaction_sample_level_v2.tsv")
)

write_tsv(
    pca_variance_df,
    file.path(wdir, "HSC_pseudobulk_PCA_variance_explained_v2.tsv")
)

save.image("2026_05_31_vivo_hsc_pca_density_contour_v2.RData")

message("Finished: pseudo-bulk statistics, PC1-PC2 replicate plot, PERMDISP per-group dispersion, and output tables.")
