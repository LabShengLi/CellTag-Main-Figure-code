#!/usr/bin/env Rscript
rm(list = ls())
library(Matrix)
library(ggplot2)

wdir <- "celltag_qc_plot"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)

exp_name <- "figs3_d_exp2_vitro_old"
dsname_list <- c("Old1", "Old2", "Old3")
jaccard_fn_list <- c(
    "4_Oa_mef_Jaccard_mtx.RDS",
    "5_Oa_mef_Jaccard_mtx.RDS",
    "6_Oa_mef_Jaccard_mtx.RDS"
)

read_jaccard_vec <- function(fn) {
    jcm <- readRDS(fn)
    if (inherits(jcm, "Matrix")) jcm <- as.matrix(jcm)
    v <- jcm[upper.tri(jcm, diag = FALSE)]
    v <- v[is.finite(v)]
    as.numeric(v)
}

vec_list <- lapply(jaccard_fn_list, read_jaccard_vec)
names(vec_list) <- dsname_list

# auto colors
cols <- setNames(rainbow(length(vec_list)), names(vec_list))

# ---- build long data.frame: value + dsname + index (within each dataset)
df <- do.call(rbind, lapply(names(vec_list), function(nm) {
    v <- vec_list[[nm]]
    data.frame(
        dsname = nm,
        idx = seq_along(v),
        value = v
    )
}))
df$dsname <- factor(df$dsname, levels = dsname_list)

# ---- plot: one figure with 3 columns (facets)
p <- ggplot(df, aes(x = idx, y = value)) +
    geom_point(size = 0.2) +
    facet_wrap(~dsname, ncol = 3, scales = "free_x") +
    labs(x = NULL, y = "Jaccard Similarity") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

# p

setwd(wdir)

# PDF (3 inch)
outfn <- sprintf("%s_Jaccard_3datasets_onepanel.pdf", exp_name)
# pdf(outfn, width = 3, height = 3)
# print(p)
# dev.off()

ggsave(
    outfn,
    plot = p,
    width = 2.5,
    height = 2.5,
    dpi = 300
)

# PNG (3 inch @300 dpi)
dpi <- 300
outfn <- sprintf("%s_Jaccard_3datasets_onepanel.png", exp_name)
# png(outfn, width = 5 * dpi, height = 4 * dpi, res = dpi)
# print(p)
# dev.off()

ggsave(
    outfn,
    plot = p,
    width = 2.5,
    height = 2.5,
    dpi = 300
)

cat("Done\n")
