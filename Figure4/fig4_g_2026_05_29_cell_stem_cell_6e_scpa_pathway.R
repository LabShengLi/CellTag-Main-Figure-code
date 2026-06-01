#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
    library(SCPA)
    library(dplyr)
    library(purrr)
    library(ComplexHeatmap)
    library(circlize)
    library(tibble)
    library(stringr)
    library(readr)
})

# -------------------------------------------------------------------------
# SCPA qvalue reference:
# https://github.com/jackbibby1/SCPA/discussions/79
#
# qval = sqrt(-log10(adjPval)); signed_qval = sign(FC) * qval
# Heatmap values are signed -log10Q (same scale as 2026_05_28 plots).
# -------------------------------------------------------------------------

wdir <- "figure4"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)

infn2 <- "CrossAge_vivo.RDS"

# -------------------------------------------------------------------------
# Pathway definitions and SCPA helper
# -------------------------------------------------------------------------
hallmark_pathways <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>%
    SCPA::format_pathways()

run_scpa_and_save <- function(obj, label, group1_population) {
    res <- SCPA::compare_seurat(
        obj,
        group1 = "orig.ident",
        group1_population = group1_population,
        group2 = "celltype_final",
        group2_population = "HSC",
        pathways = c(hallmark_pathways),
        downsample = 500,
        parallel = TRUE,
        cores = 10,
        min_genes = 9
    ) %>%
        mutate(
            signed_qval = case_when(
                FC > 0 ~ qval,
                FC < 0 ~ -qval,
                TRUE ~ 0
            )
        )

    saveRDS(res, file = file.path(wdir, paste0(label, "_pathway_result.rds")))
    write_tsv(res, file = file.path(wdir, paste0(label, "_pathway_result.tsv")))
    res
}

load_or_run_scpa <- function(obj, label, group1_population) {
    rds_path <- file.path(wdir, paste0(label, "_pathway_result.rds"))
    if (file.exists(rds_path)) {
        message("Loading cached SCPA result: ", rds_path)
        readRDS(rds_path)
    } else {
        message("Running SCPA: ", label)
        run_scpa_and_save(obj, label, group1_population)
    }
}

# -------------------------------------------------------------------------
# Run / load SCPA comparisons
# -------------------------------------------------------------------------
ss <- readRDS(infn2)
setwd(wdir)
ss <- subset(ss, subset = CloneID != 0)

OY_vs_OO_path <- load_or_run_scpa(ss, label = "OY_vs_OO", group1_population = c("OY", "OO"))
YO_vs_YY_path <- load_or_run_scpa(ss, label = "YO_vs_YY", group1_population = c("YO", "YY"))

# -------------------------------------------------------------------------
# Merged pathway result table
# -------------------------------------------------------------------------
selected_pathways <- c(
    "G2M_CHECKPOINT",
    "E2F_TARGETS",
    "MITOTIC_SPINDLE",
    "MTORC1_SIGNALING",
    "ESTROGEN_RESPONSE_LATE",
    "MYC_TARGETS_V1",
    "OXIDATIVE_PHOSPHORYLATION",
    "UV_RESPONSE_DN",
    "PI3K_AKT_MTOR_SIGNALING",
    "ESTROGEN_RESPONSE_EARLY"
)

scpa_merged_tbl <- list(
    OY_vs_OO = OY_vs_OO_path,
    YO_vs_YY = YO_vs_YY_path
) %>%
    purrr::imap(function(df, comp) {
        df %>%
            dplyr::select(Pathway, Pval, adjPval, qval, FC, signed_qval) %>%
            dplyr::rename(
                !!paste0(comp, "_Pval") := Pval,
                !!paste0(comp, "_adjPval") := adjPval,
                !!paste0(comp, "_qval") := qval,
                !!paste0(comp, "_FC") := FC,
                !!paste0(comp, "_signed_qval") := signed_qval
            )
    }) %>%
    purrr::reduce(dplyr::full_join, by = "Pathway") %>%
    dplyr::mutate(Pathway = stringr::str_remove(Pathway, "^HALLMARK_")) %>%
    dplyr::arrange(Pathway)

selected_pathways <- selected_pathways[selected_pathways %in% scpa_merged_tbl$Pathway]

write_tsv(scpa_merged_tbl, file = file.path(wdir, "2026_05_29_SCPA_merged_table_OY_OO_YO_YY.tsv"))
saveRDS(scpa_merged_tbl, file = file.path(wdir, "2026_05_29_SCPA_merged_table_OY_OO_YO_YY.rds"))

scpa_merged_tbl_selected <- scpa_merged_tbl %>%
    dplyr::filter(Pathway %in% selected_pathways)

write_tsv(
    scpa_merged_tbl_selected,
    file = file.path(wdir, "2026_05_29_SCPA_merged_table_OY_OO_YO_YY_selected.tsv")
)
saveRDS(
    scpa_merged_tbl_selected,
    file = file.path(wdir, "2026_05_29_SCPA_merged_table_OY_OO_YO_YY_selected.rds")
)

# -------------------------------------------------------------------------
# Signed -log10Q matrices (all pathways + heatmap-filtered + selected)
# -------------------------------------------------------------------------
pathways_mat_all <- list(
    OY_vs_OO = OY_vs_OO_path,
    YO_vs_YY = YO_vs_YY_path
) %>%
    purrr::map(\(x) dplyr::select(x, Pathway, signed_qval)) %>%
    purrr::imap(\(x, y) dplyr::rename(x, !!y := signed_qval)) %>%
    purrr::reduce(dplyr::full_join, by = "Pathway") %>%
    dplyr::mutate(Pathway = stringr::str_remove(Pathway, "^HALLMARK_")) %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()

pathways_mat <- pathways_mat_all[
    apply(pathways_mat_all, 1, function(x) any(abs(x) > 2, na.rm = TRUE)),
    ,
    drop = FALSE
]

pathways_mat_selected <- pathways_mat_all[selected_pathways, , drop = FALSE]

saveRDS(pathways_mat, file = file.path(wdir, "2026_05_29_pathways_OY_OO_and_YO_YY_mat.rds"))
saveRDS(pathways_mat_selected, file = file.path(wdir, "2026_05_29_pathways_OY_OO_and_YO_YY_selected_mat.rds"))

write.table(
    pathways_mat,
    file = file.path(wdir, "2026_05_29_pathways_OY_OO_and_YO_YY_mat.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)
write.table(
    pathways_mat_selected,
    file = file.path(wdir, "2026_05_29_pathways_OY_OO_and_YO_YY_selected_mat.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

# -------------------------------------------------------------------------
# Heatmaps: all pathways and selected pathways
# -------------------------------------------------------------------------
heatmap_col <- circlize::colorRamp2(
    c(-3, 0, 3),
    c("#2166AC", "white", "#B2182B")
)

make_pathway_heatmap <- function(mat, height_inch = 6) {
    Heatmap(
        mat,
        name = "signed -log10Qvalue\n(based on FC)",
        col = heatmap_col,
        show_row_names = TRUE,
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        border = TRUE,
        row_names_gp = grid::gpar(fontsize = 9),
        column_names_gp = grid::gpar(fontsize = 16),
        width = grid::unit(2.2, "inch"),
        height = grid::unit(height_inch, "inch")
    )
}

hm <- make_pathway_heatmap(pathways_mat, height_inch = 6)
draw(hm)

pdf(file.path(wdir, "2026_05_29_OY_OO_YO_YY_SCPA-pathways-enrichment.pdf"), height = 9, width = 7.5)
draw(hm)
invisible(dev.off())

heatmap_col_v2 <- circlize::colorRamp2(
    c(-5, 0, 5),
    c("#2166AC", "white", "#B2182B")
)

hm_selected <- Heatmap(
    pathways_mat_selected,
    name = "-log10Qvalue\n(signed based on FC)",
    col = heatmap_col_v2,
    show_row_names = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    border = TRUE,
    show_heatmap_legend = FALSE,
    row_names_gp = grid::gpar(fontsize = 16),
    column_names_gp = grid::gpar(fontsize = 16),
    width = grid::unit(1.2, "inch"),
    height = grid::unit(max(3, 0.35 * nrow(pathways_mat_selected)), "inch")
)

lgd_selected <- Legend(
    col_fun = heatmap_col_v2,
    title = "-log10Qvalue\n(signed based on FC)",
    at = c(-5, 0, 5),
    labels = c("-5", "0", "5"),
    direction = "horizontal",
    title_position = "topcenter",
    legend_width = grid::unit(4, "cm"),
    legend_height = grid::unit(0.4, "cm"),
    labels_gp = grid::gpar(fontsize = 14),
    title_gp = grid::gpar(fontsize = 14),
    border = FALSE
)

draw_hm_selected <- function() {
    draw(hm_selected, padding = grid::unit(c(14, 2, 2, 2), "mm"))
    draw(
        lgd_selected,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.98, "npc"),
        just = c("center", "top")
    )
}

draw_hm_selected()

pdf(file.path(wdir, "2026_05_29_OY_OO_YO_YY_SCPA-pathways-enrichment_selected.pdf"), height = 7.5, width = 7)
draw_hm_selected()
invisible(dev.off())

# -------------------------------------------------------------------------
# Paired Wilcoxon test: qval (not signed_qval), OY_vs_OO vs YO_vs_YY
# -------------------------------------------------------------------------
col_oy_oo <- "OY_vs_OO"
col_yo_yy <- "YO_vs_YY"

pathways_qval_mat_all <- list(
    OY_vs_OO = OY_vs_OO_path,
    YO_vs_YY = YO_vs_YY_path
) %>%
    purrr::map(\(x) dplyr::select(x, Pathway, qval)) %>%
    purrr::imap(\(x, y) dplyr::rename(x, !!y := qval)) %>%
    purrr::reduce(dplyr::full_join, by = "Pathway") %>%
    dplyr::mutate(Pathway = stringr::str_remove(Pathway, "^HALLMARK_")) %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()

pathways_qval_mat_selected <- pathways_qval_mat_all[selected_pathways, , drop = FALSE]

build_paired_qval_df <- function(mat, label) {
    tibble::tibble(
        pathway_set = label,
        Pathway = rownames(mat),
        OY_vs_OO_qval = mat[, col_oy_oo],
        YO_vs_YY_qval = mat[, col_yo_yy],
        diff_OY_OO_minus_YO_YY = mat[, col_oy_oo] - mat[, col_yo_yy]
    )
}

run_paired_wilcox_qval <- function(mat, pathway_set_label) {
    x <- mat[, col_oy_oo]
    y <- mat[, col_yo_yy]
    ok <- stats::complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]

    test_res <- stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)

    list(
        pathway_set = pathway_set_label,
        n_pathways = length(x),
        median_OY_vs_OO = stats::median(x, na.rm = TRUE),
        median_YO_vs_YY = stats::median(y, na.rm = TRUE),
        median_diff = stats::median(x - y, na.rm = TRUE),
        statistic = unname(test_res$statistic),
        p_value = test_res$p.value,
        alternative = test_res$alternative,
        method = test_res$method
    )
}

paired_qval_all <- build_paired_qval_df(pathways_qval_mat_all, "all_pathways")
paired_qval_selected <- build_paired_qval_df(pathways_qval_mat_selected, "selected_pathways")

wilcox_summary <- bind_rows(
    run_paired_wilcox_qval(pathways_qval_mat_all, "all_pathways"),
    run_paired_wilcox_qval(pathways_qval_mat_selected, "selected_pathways")
)

wilcox_summary_print <- wilcox_summary %>%
    mutate(
        p_value_fmt = format.pval(p_value, digits = 4, eps = 1e-16),
        significant_0.05 = p_value < 0.05
    )

print(wilcox_summary_print)

write_tsv(paired_qval_all, file.path(wdir, "2026_05_29_paired_qval_all_pathways.tsv"))
write_tsv(paired_qval_selected, file.path(wdir, "2026_05_29_paired_qval_selected_pathways.tsv"))
write_tsv(wilcox_summary_print, file.path(wdir, "2026_05_29_paired_wilcox_qval_summary.tsv"))

wilcox_txt <- file.path(wdir, "2026_05_29_paired_wilcox_qval.txt")
sink(wilcox_txt)
cat("Paired Wilcoxon test: qval (OY_vs_OO vs YO_vs_YY)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working directory:", wdir, "\n\n")
cat("Metric: SCPA qval = sqrt(-log10(adjPval)); unsigned (not signed_qval).\n")
cat("Test: paired Wilcoxon rank-sum test.\n")
cat("Comparison columns:", col_oy_oo, "vs", col_yo_yy, "\n\n")

cat("=== Summary table ===\n")
print(as.data.frame(wilcox_summary_print), row.names = FALSE)

for (df_label in c("all_pathways", "selected_pathways")) {
    df_use <- switch(
        df_label,
        all_pathways = paired_qval_all,
        selected_pathways = paired_qval_selected
    )
    cat("\n=== Per-pathway qval (", df_label, ") ===\n", sep = "")
    print(as.data.frame(df_use), row.names = FALSE)
}

sink()

message("Saved paired Wilcox results: ", wilcox_txt)
message("Finished pathway results, heatmaps, and paired Wilcoxon tests.")
