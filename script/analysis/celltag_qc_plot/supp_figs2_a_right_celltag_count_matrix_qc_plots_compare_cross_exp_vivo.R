#!/usr/bin/env Rscript
rm(list = ls())
library(dplyr)
library(ggplot2)

dsname <- "figs2_a_right_panel_cross_exp_vivo_4cases"
wdir <- "/project2/sli68423_1316/projects/U01_aim2/results/2026_01_24_celltag_qc"

cell_tag_matrix_fn_list <- c(
    "/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/ExtractBarcodes/Out/YY/hf1.d15.v1.YY.celltag.matrix.Rds",
    "/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/ExtractBarcodes/Out/YO/hf1.d15.v1.YO.celltag.matrix.Rds",
    "/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/ExtractBarcodes/Out/OY/hf1.d15.v1.OY.celltag.matrix.Rds",
    "/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/ExtractBarcodes/Out/OO/hf1.d15.v1.OO.celltag.matrix.Rds"
)

dsname_list <- c("YY", "YO", "OY", "OO")
umi_cutoff <- 2

dir.create(wdir, recursive = TRUE, showWarnings = FALSE)
setwd(wdir)

# -------------------------
# Helper: load + sanitize
# -------------------------
load_ctm <- function(fn) {
    ctm <- readRDS(fn)
    ctm <- as.data.frame(ctm)
    if ("Cell.BC" %in% colnames(ctm)) {
        rownames(ctm) <- ctm$Cell.BC
        ctm <- ctm[, setdiff(colnames(ctm), "Cell.BC"), drop = FALSE]
    }
    ctm
}

# -------------------------
# 1) Combine n_tags per cell for all datasets
# -------------------------
df_n_tags <- bind_rows(lapply(seq_along(cell_tag_matrix_fn_list), function(i) {
    ctm <- load_ctm(cell_tag_matrix_fn_list[i])
    n_tags <- rowSums(ctm >= umi_cutoff, na.rm = TRUE)
    data.frame(dataset = dsname_list[i], n_tags = n_tags)
}))

# (optional) for log10 boxplot, keep >=1 (or >=2 if you prefer)
df_n_tags_f <- df_n_tags %>% filter(n_tags >= 1)

df_n_tags_f$dataset <- factor(df_n_tags_f$dataset, levels = dsname_list)

bp1_3 <- ggplot(df_n_tags_f, aes(x = dataset, y = n_tags)) +
    geom_boxplot(width = 0.55,
                 outlier.shape = 1, outlier.size = 0.8, outlier.stroke = 0.5) +
    scale_y_log10(breaks = c(1, 2, 5, 10, 20)) +
    labs(x = NULL, y = "CellTags per cell") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.8))

bp1_3
ggsave(sprintf("%s_CellTag_per_cell_boxplot_4datasets.pdf", dsname),
       plot = bp1_3, width = 3, height = 3)

summary(df_n_tags_f[df_n_tags_f$dataset == 'YY', 2])
summary(df_n_tags_f[df_n_tags_f$dataset == 'YO', 2])

summary(df_n_tags_f[df_n_tags_f$dataset == 'OY', 2])
summary(df_n_tags_f[df_n_tags_f$dataset == 'OO', 2])

# > summary(df_n_tags_f[df_n_tags_f$dataset == 'YY', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   1.000   2.043   2.000  17.000
# > summary(df_n_tags_f[df_n_tags_f$dataset == 'YO', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   1.000   2.288   3.000  17.000
# >
# > summary(df_n_tags_f[df_n_tags_f$dataset == 'OY', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   2.000   2.298   3.000  15.000
# > summary(df_n_tags_f[df_n_tags_f$dataset == 'OO', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   1.000   1.721   2.000  12.000

# -------------------------
# 2) Combine tag_ncells (cells per tag) for all datasets
# -------------------------
df_tagfreq <- bind_rows(lapply(seq_along(cell_tag_matrix_fn_list), function(i) {
    ctm <- load_ctm(cell_tag_matrix_fn_list[i])
    tag_ncells <- colSums(ctm >= umi_cutoff, na.rm = TRUE)
    tag_ncells_pos <- tag_ncells[tag_ncells > 0]
    data.frame(dataset = dsname_list[i], tag_ncells = as.numeric(tag_ncells_pos))
}))

df_tagfreq$dataset <- factor(df_tagfreq$dataset, levels = dsname_list)

bp2_3 <- ggplot(df_tagfreq, aes(x = dataset, y = tag_ncells)) +
    geom_boxplot(width = 0.55,
                 outlier.shape = 1, outlier.size = 0.6, outlier.stroke = 0.3) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_y_continuous(breaks = seq(0, 20, by = 2)) +
    labs(x = NULL, y = "CellTag Frequency") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.8))

bp2_3
ggsave(sprintf("%s_CellTag_freqency_boxplot_4datasets.pdf", dsname),
       plot = bp2_3, width = 3, height = 3)


summary(df_tagfreq[df_tagfreq$dataset == 'YY', 2])
summary(df_tagfreq[df_tagfreq$dataset == 'YO', 2])

summary(df_tagfreq[df_tagfreq$dataset == 'OY', 2])
summary(df_tagfreq[df_tagfreq$dataset == 'OO', 2])


# > summary(df_tagfreq[df_tagfreq$dataset == 'YY', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1.00    1.00    4.00   14.14   15.25  211.00
# > summary(df_tagfreq[df_tagfreq$dataset == 'YO', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     1.0     1.0     3.0    10.3     9.0   181.0
# >
# > summary(df_tagfreq[df_tagfreq$dataset == 'OY', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   1.000   2.384   2.000  42.000
# > summary(df_tagfreq[df_tagfreq$dataset == 'OO', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1.00    1.00    2.00    5.31    6.00  112.00

##
# -------------------------
# 2) Combine plot histogram for all datasets
# -------------------------
df_hist_ntags <- bind_rows(lapply(seq_along(cell_tag_matrix_fn_list), function(i) {
    ctm <- load_ctm(cell_tag_matrix_fn_list[i])
    cell_ntags <- rowSums(ctm >= umi_cutoff, na.rm = TRUE)
    data.frame(dataset = dsname_list[i], cell_ntags = as.numeric(cell_ntags))
}))

topk <- 5
df2 <- df_hist_ntags %>% dplyr::filter(cell_ntags <= topk & cell_ntags > 0)

df2$dataset <- factor(df2$dataset, levels = dsname_list)


dodge <- position_dodge2(width = 0.9, preserve = "single")

hp1 <- ggplot(df2, aes(x = cell_ntags, fill = dataset)) +
    geom_bar(width = 0.9, position = dodge) +
    geom_text(
        stat = "count",
        aes(label = scales::comma(after_stat(count))),
        position = dodge,
        angle = 90,
        hjust = -0.1,
        vjust = 0.4,
        size = 3
    ) +
    scale_x_continuous(breaks = 0:max(df2$cell_ntags)) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.15))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        plot.margin = margin(5.5, 20, 5.5, 5.5),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "CellTag number", y = "Cell number", fill = NULL)

hp1
ggsave(sprintf("%s_CellTag_distribution_histogrm.pdf", dsname),
       plot = hp1, width = 4, height = 3.5)

cat("Done\n")
