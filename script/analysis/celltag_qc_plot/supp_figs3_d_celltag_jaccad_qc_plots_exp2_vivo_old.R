#!/usr/bin/env Rscript
rm(list = ls())
library(Matrix)
library(ggplot2)

wdir <- "/project2/sli68423_1316/projects/U01_aim2/results/2026_01_24_celltag_qc"
input_dir <- file.path(wdir, "exp2_celltag_input")

dir.create(wdir, showWarnings = FALSE, recursive = TRUE)
setwd(wdir)

exp_name <- "figs3_d_exp2_vivo_old"

# Set these three names to match your current filename
dsname_list <- c("Old1B", "Old2B", "Old3B")

# Automatically read RDS files from the input directory
jaccard_fn_list <- sort(list.files(
    input_dir,
    pattern = "^O[123]B.*\\.RDS$",
    full.names = TRUE
))

if (length(jaccard_fn_list) != 3) {
    stop("Expected 3 RDS files, but found: ", length(jaccard_fn_list))
}

read_jaccard_vec <- function(fn) {
    jcm <- readRDS(fn)
    if (inherits(jcm, "Matrix")) jcm <- as.matrix(jcm)
    v <- jcm[upper.tri(jcm, diag = FALSE)]
    v <- v[is.finite(v)]
    as.numeric(v)
}

vec_list <- lapply(jaccard_fn_list, read_jaccard_vec)
names(vec_list) <- dsname_list

# long data.frame
df <- do.call(rbind, lapply(names(vec_list), function(nm) {
    v <- vec_list[[nm]]
    data.frame(
        dsname = nm,
        idx = seq_along(v),
        value = v
    )
}))
df$dsname <- factor(df$dsname, levels = dsname_list)

# plot
p <- ggplot(df, aes(x = idx, y = value)) +
    geom_point(size = 0.2) +
    facet_wrap(~dsname, ncol = 3, scales = "free_x") +
    labs(x = NULL, y = "Jaccard Similarity") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

# print(p)

# PNG
dpi <- 300
outfn <- sprintf("%s_Jaccard_3datasets_onepanel.png", exp_name)
# png(outfn, width = 3 * dpi, height = 3 * dpi, res = dpi)
# print(p)
# dev.off()

ggsave(
    outfn,
    plot = p,
    width = 2.5,
    height = 2.5,
    dpi = 300
)

# PDF
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

cat("Done\n")
