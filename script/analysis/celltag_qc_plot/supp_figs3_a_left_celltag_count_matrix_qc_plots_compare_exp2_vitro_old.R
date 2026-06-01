#!/usr/bin/env Rscript
rm(list = ls())
library(dplyr)
library(ggplot2)

dsname <- "figs3_a_left_exp2_vitro_old"
wdir <- "celltag_qc_plot"

cell_tag_matrix_fn_list <- c(
    "4_Oa_hf1.d15.v1.celltag.matrix.Rds",
    "5_Oa_hf1.d15.v1.celltag.matrix.Rds",
    "6_Oa_hf1.d15.v1.celltag.matrix.Rds"
)

dsname_list <- c("Old1", "Old2", "Old3")
umi_cutoff <- 2

dir.create(wdir, recursive = TRUE, showWarnings = FALSE)

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

bp1_3 <- ggplot(df_n_tags_f, aes(x = dataset, y = n_tags)) +
    geom_boxplot(width = 0.55,
                 outlier.shape = 1, outlier.size = 0.8, outlier.stroke = 0.5) +
    scale_y_log10(breaks = c(1, 2, 5, 10, 20)) +
    labs(x = NULL, y = "CellTags per cell") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.8))

bp1_3
setwd(wdir)
ggsave(sprintf("%s_CellTag_per_cell_boxplot_3datasets.pdf", dsname),
       plot = bp1_3, width = 3, height = 3)

summary(df_n_tags_f[df_n_tags_f$dataset == 'Old1', 2])
summary(df_n_tags_f[df_n_tags_f$dataset == 'Old2', 2])
summary(df_n_tags_f[df_n_tags_f$dataset == 'Old3', 2])

# > summary(df_n_tags_f[df_n_tags_f$dataset == 'Old1', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1.00    1.00    3.00    4.44    6.00   38.00
# > summary(df_n_tags_f[df_n_tags_f$dataset == 'Old2', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   2.000   4.000   4.736   7.000  31.000
# > summary(df_n_tags_f[df_n_tags_f$dataset == 'Old3', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   3.000   5.000   5.864   8.000  34.000

# -------------------------
# 2) Combine tag_ncells (cells per tag) for all datasets
# -------------------------
df_tagfreq <- bind_rows(lapply(seq_along(cell_tag_matrix_fn_list), function(i) {
    ctm <- load_ctm(cell_tag_matrix_fn_list[i])
    tag_ncells <- colSums(ctm >= umi_cutoff, na.rm = TRUE)
    tag_ncells_pos <- tag_ncells[tag_ncells > 0]
    data.frame(dataset = dsname_list[i], tag_ncells = as.numeric(tag_ncells_pos))
}))

bp2_3 <- ggplot(df_tagfreq, aes(x = dataset, y = tag_ncells)) +
    geom_boxplot(width = 0.55,
                 outlier.shape = 1, outlier.size = 0.6, outlier.stroke = 0.3) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_y_continuous(breaks = seq(0, 20, by = 2)) +
    labs(x = NULL, y = "CellTag Frequency") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.8))

bp2_3
ggsave(sprintf("%s_CellTag_cells_per_tag_boxplot_3datasets.pdf", dsname),
       plot = bp2_3, width = 3, height = 3)


summary(df_tagfreq[df_tagfreq$dataset == 'Old1', 2])
summary(df_tagfreq[df_tagfreq$dataset == 'Old2', 2])
summary(df_tagfreq[df_tagfreq$dataset == 'Old3', 2])


# > summary(df_tagfreq[df_tagfreq$dataset == 'Old1', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   2.000   4.000   7.729  10.000 172.000
# > summary(df_tagfreq[df_tagfreq$dataset == 'Old2', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1.00    1.00    3.00    6.12    7.00  109.00
# > summary(df_tagfreq[df_tagfreq$dataset == 'Old3', 2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   2.000   4.000   7.176   9.000 136.000

##
# -------------------------
# 2) Combine plot histogram for all datasets
# -------------------------
df_hist_ntags <- bind_rows(lapply(seq_along(cell_tag_matrix_fn_list), function(i) {
    ctm <- load_ctm(cell_tag_matrix_fn_list[i])
    cell_ntags <- rowSums(ctm >= umi_cutoff, na.rm = TRUE)
    data.frame(dataset = dsname_list[i], cell_ntags = as.numeric(cell_ntags))
}))


topk <- 15
df2 <- df_hist_ntags %>% dplyr::filter(cell_ntags <= topk)

# Histogram (all), same as report from Sheng
hp1 <- ggplot(df2, aes(x = cell_ntags)) +
    geom_bar(width = 0.9) +
    geom_text(
        stat = "count",
        aes(label = ifelse(after_stat(count) >= 100,
                           scales::comma(after_stat(count)),
                           "")),
        angle = 90,
        hjust = -0.1,
        vjust = 0.4,
        size = 3
    ) +
    scale_x_continuous(breaks = 0:max(df2$cell_ntags)) +
    scale_y_continuous(
        labels = scales::comma,
        expand = expansion(mult = c(0, 0.15))   # Leave 15% headroom at the top
    ) +
    coord_cartesian(clip = "off") +           # Do not clip
    theme_classic() +
    theme(plot.margin = margin(5.5, 20, 5.5, 5.5)) +  # Right/top plot margin
    labs(x = "CellTag number", y = "Cell number")

hp1

ggsave(sprintf("%s_CellTag_distribution_histogrm_hp1.pdf", dsname),
       plot = hp1, width = 4, height = 3)


topk <- 5
df3 <- df_hist_ntags %>% dplyr::filter(cell_ntags <= topk & cell_ntags > 0)

dodge <- position_dodge2(width = 0.9, preserve = "single")


hp2 <- ggplot(df3, aes(x = cell_ntags, fill = dataset)) +
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
        legend.position = c(0.98, 0.98),   # Place inside top-right
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "CellTag number", y = "Cell number", fill = NULL)

hp2

ggsave(sprintf("%s_CellTag_distribution_histogrm_hp2.pdf", dsname),
       plot = hp2, width = 4, height = 3)

cat("Done\n")
