#!/usr/bin/env Rscript

# ============================================================
# Purpose:
# This script reads Jaccard similarity matrices (RDS format)
# from multiple datasets, extracts upper-triangular values,
# and visualizes their distributions using scatter plots.
#
# NOTE:
# - Output directory (wdir) is kept as the original results folder
# - Input files are read from a separate "Out" directory
# ============================================================

rm(list = ls())

library(Matrix)
library(ggplot2)

# ============================================================
# 1) Output working directory
#    All figures will be saved here
# ============================================================
wdir <- "celltag_qc_plot"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2) Input files (mapped filenames from copy_celltag_qc_inputs.sh)
# ============================================================
exp_name <- "figs2_d_crossexp4_vivo"
dsname_list <- c("OO", "OY", "YO", "YY")
jaccard_fn_list <- c(
    "OO_mef_Jaccard_mtx.RDS",
    "OY_mef_Jaccard_mtx.RDS",
    "YO_mef_Jaccard_mtx.RDS",
    "YY_mef_Jaccard_mtx.RDS"
)

# ============================================================
# 6) Function to read and process Jaccard matrix
#    - Reads RDS file
#    - Converts sparse Matrix to dense matrix if needed
#    - Extracts upper triangle (excluding diagonal)
#    - Removes NA/Inf values
# ============================================================
read_jaccard_vec <- function(fn) {
    jcm <- readRDS(fn)
    if (inherits(jcm, "Matrix")) jcm <- as.matrix(jcm)
    v <- jcm[upper.tri(jcm, diag = FALSE)]
    v <- v[is.finite(v)]
    as.numeric(v)
}

# ============================================================
# 7) Load all datasets into a list of numeric vectors
# ============================================================
vec_list <- lapply(jaccard_fn_list, read_jaccard_vec)
names(vec_list) <- dsname_list

# ============================================================
# 8) Convert to long-format data.frame for ggplot
#    - dsname: dataset label
#    - idx: index within dataset
#    - value: Jaccard similarity
# ============================================================
df <- do.call(rbind, lapply(names(vec_list), function(nm) {
    v <- vec_list[[nm]]
    data.frame(
        dsname = nm,
        idx = seq_along(v),
        value = v
    )
}))
df$dsname <- factor(df$dsname, levels = dsname_list)

# ============================================================
# 9) Plot
#    - Each dataset is shown in a separate facet panel
#    - X-axis is hidden for clarity
# ============================================================
p <- ggplot(df, aes(x = idx, y = value)) +
    geom_point(size = 0.2) +
    facet_wrap(~dsname, ncol = 2, scales = "free_x") +
    labs(x = NULL, y = "Jaccard Similarity") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

# ============================================================
# 10) Save outputs
# ============================================================
setwd(wdir)

# PNG
outfn_png <- sprintf("%s_Jaccard_4datasets_onepanel.png", exp_name)
ggsave(outfn_png, plot = p, width = 3, height = 2.5, dpi = 300)

# PDF
outfn_pdf <- sprintf("%s_Jaccard_4datasets_onepanel.pdf", exp_name)
ggsave(outfn_pdf, plot = p, width = 3, height = 2.5)

cat("Done\n")
