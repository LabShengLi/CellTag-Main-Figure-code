#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
    library(readr)
    library(Seurat)
    library(ggplot2)
    library(ggrepel)
})

# -------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------
wdir <- "/project2/sli68423_1316/projects/U01_aim2/results/2026_05_17_scpa"
dir.create(wdir, recursive = TRUE, showWarnings = FALSE)
setwd(wdir)

# -------------------------------------------------------------------------
# Input data
# -------------------------------------------------------------------------
infn2 <- "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"

ss <- readRDS(infn2)
ss <- subset(ss, subset = CloneID != 0)
ss_hsc <- subset(ss, subset = celltype_final == "HSC")

group_order <- c("OO", "OY", "YO", "YY")
group_order <- group_order[group_order %in% unique(ss_hsc$orig.ident)]

rep_order <- c("Rep1", "Rep2", "Rep3")
rep_order <- rep_order[rep_order %in% unique(ss_hsc$Rep)]

# -------------------------------------------------------------------------
# PCA (group centroid = per-PC median)
# -------------------------------------------------------------------------
compute_centroid_distance <- function(emb_mat, meta_group, group_order) {
    centroid_list <- lapply(group_order, function(g) {
        cells_use <- names(meta_group)[meta_group == g]
        apply(emb_mat[cells_use, , drop = FALSE], 2, stats::median)
    })
    names(centroid_list) <- group_order
    centroid_mat <- do.call(rbind, centroid_list)
    list(centroid_mat = centroid_mat)
}

DefaultAssay(ss_hsc) <- "RNA"
ss_hsc <- NormalizeData(ss_hsc, verbose = FALSE)
ss_hsc <- FindVariableFeatures(ss_hsc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ss_hsc <- ScaleData(ss_hsc, features = VariableFeatures(ss_hsc), verbose = FALSE)
ss_hsc <- RunPCA(ss_hsc, features = VariableFeatures(ss_hsc), npcs = 30, verbose = FALSE)

pca_emb <- Embeddings(ss_hsc, reduction = "pca")
pcs_use <- 1:30
pcs_use <- pcs_use[pcs_use <= ncol(pca_emb)]
pca_emb <- pca_emb[, pcs_use, drop = FALSE]

pca_res <- compute_centroid_distance(
    emb_mat = pca_emb,
    meta_group = ss_hsc$orig.ident,
    group_order = group_order
)

pca_plot_df <- data.frame(
    PC1 = pca_emb[, 1],
    PC2 = pca_emb[, 2],
    orig.ident = factor(ss_hsc$orig.ident, levels = group_order),
    Rep = factor(ss_hsc$Rep, levels = rep_order),
    row.names = rownames(pca_emb)
)

pca_centroid_df <- as.data.frame(pca_res$centroid_mat[, 1:2, drop = FALSE])
colnames(pca_centroid_df) <- c("PC1", "PC2")
pca_centroid_df$orig.ident <- factor(rownames(pca_centroid_df), levels = group_order)

# Match 2026_05_28 group colors; marker shapes: OO/OY solid triangle, YO/YY solid circle
hollow_shape_values <- c(Rep1 = 1, Rep2 = 0, Rep3 = 2)
hollow_shape_values <- hollow_shape_values[names(hollow_shape_values) %in% rep_order]

centroid_shape_values <- c(OO = 17, OY = 17, YO = 16, YY = 16)
centroid_shape_values <- centroid_shape_values[names(centroid_shape_values) %in% group_order]

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
p_ref <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = orig.ident, shape = Rep)) +
    geom_point(size = 1.5, stroke = 0.4, fill = NA, alpha = 0.8) +
    scale_shape_manual(values = hollow_shape_values) +
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
                dg$x >= search_lims$xlim[1] & dg$x <= search_lims$xlim[2] &
                    dg$y >= search_lims$ylim[1] & dg$y <= search_lims$ylim[2],
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
)

message("Density peaks (mode) from ggplot2 KDE grid:")
print(pca_density_peak_df)

pca_density_peak_dist_mat <- compute_pc12_euclidean_dist_matrix(
    pca_density_peak_df,
    group_order = group_order
)
pca_median_dist_mat <- compute_pc12_euclidean_dist_matrix(
    pca_centroid_df,
    group_order = group_order
)

message("Density peak Euclidean distance matrix (PC1, PC2):")
print(round(pca_density_peak_dist_mat, 4))
message("Median Euclidean distance matrix (PC1, PC2):")
print(round(pca_median_dist_mat, 4))

# Legend keys: OO/OY triangle (17), YO/YY circle (16); contours excluded from legend
legend_shape_values <- unname(centroid_shape_values[group_order])
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
    marker_shape_values,
    legend_shape_values,
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
            aes(x = PC1, y = PC2, color = orig.ident, fill = orig.ident, shape = orig.ident),
            inherit.aes = FALSE,
            size = 5,
            stroke = 0.6
        ) +
        scale_shape_manual(
            values = marker_shape_values,
            limits = names(marker_shape_values),
            drop = FALSE
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
                override.aes = list(
                    shape = legend_shape_values,
                    fill = legend_fill_values,
                    colour = legend_fill_values,
                    linetype = "blank",
                    linewidth = 0,
                    stroke = 0.6,
                    size = 4
                )
            ),
            fill = "none",
            shape = "none"
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
            color = "Group"
        )
}

plot_args <- list(
    plot_df = pca_plot_df,
    group_colors = group_colors,
    marker_shape_values = centroid_shape_values,
    legend_shape_values = legend_shape_values,
    legend_fill_values = legend_fill_values,
    contour_lim = contour_lim
)

p_median <- do.call(
    make_density_contour_marker_plot,
    c(
        plot_args,
        list(
            marker_df = pca_centroid_df,
            title = "HSC 2D density contours and median"
        )
    )
)

p_mode <- do.call(
    make_density_contour_marker_plot,
    c(
        plot_args,
        list(
            marker_df = pca_density_peak_df,
            title = "HSC 2D density contours and peaks"
        )
    )
)

print(p_median)
print(p_mode)

ggsave(
    file.path(wdir, "HSC_density_contours_median_PC1_PC2_PCA_all_groups.pdf"),
    p_median,
    width = 4,
    height = 3
)

ggsave(
    file.path(wdir, "HSC_density_contours_density_peak_PC1_PC2_PCA_all_groups.pdf"),
    p_mode,
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
write_distance_matrix_tsv(
    pca_median_dist_mat,
    file.path(wdir, "HSC_median_PC1_PC2_euclidean_distance_matrix.tsv")
)

save.image("2026_05_29_vivo_hsc_pca_density_contour.RData")

message("Finished: density contour plots with median centroids, density peaks (mode), and PC1/PC2 distance matrices.")
