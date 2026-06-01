###################################################################################################


### Compute Differentially variance Genes (DVGs) in Cross(Day0), Unmanipulated and exp2 data


###################################################################################################

# Maintainer: Chris Chen 
# Last updated: 03/03/2026

##########################

## Load packages ##

load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","scMayoMap",
    "readxl","pheatmap","Matrix","openxlsx","gt","glue","openxlsx"
  )
  
  suppressPackageStartupMessages(
    lapply(pkgs, require, character.only = TRUE)
  )
  message("Allpackages loaded.")
}
load_all_packages()

# set working directory
setwd('Main_figures/Figure3/Figures/')

################################

# Read in the 3 dataset 

seurat_in_vitro <- readRDS("data/CrossAge(exp2)_vitro.RDS")
seurat_unmanipulated <- readRDS("data/Unmanipulated_vitro.rds")
seurat_vitro_exp2 <- readRDS("data/Exp2(exp1)_vitro.RDS")

# inspect data

table(seurat_in_vitro$sampleName) # O_vitro: 6337 ; Y_vitro: 10971
seurat_unmanipulated  # 18046 genes and 48140 cells 
table(seurat_unmanipulated$AGE) # mid: 10877 ; old:11966 ; vold:12604 ; young: 12693
seurat_vitro_exp2 # 16611 genes and 55499 cells 

############################################

## Here we have two types of function for each of the DVGs analysis 

# 1. Compute differential variance gene function 
# 2. Plot boxplot and volcano plot function 

############################################

# Prepare input data 

############################################

# Create a function to remove unrelated/technical genes 
# ------------------------------------------------
# Remove unrelated / technical genes
# ------------------------------------------------

filter_unrelated_genes <- function(seurat_obj){
  all_genes <- rownames(seurat_obj) 
  genes_filtered <- all_genes[
    !grepl("^mt-", all_genes, ignore.case = TRUE) &
      !grepl("^Rps", all_genes) &
      !grepl("^Rpl", all_genes) &
      !grepl("^Hist", all_genes) &
      !grepl("^Gm[0-9]", all_genes) &
      !grepl("Rik$", all_genes)]
  cat("Original genes:", length(all_genes), "\n")
  cat("Filtered genes:", length(genes_filtered), "\n")
  seurat_obj <- subset(seurat_obj, features = genes_filtered)
  return(seurat_obj)}

#######################

## Cross Day0 DVGs 

#######################

# Subset the HSPC populations

# seurat_in_vitro <- seurat_celltag_in_vitro

table(seurat_in_vitro$celltype)
seurat_in_vitro_HSC <- subset(seurat_in_vitro,subset = celltype %in% c("LT-HSC"))
# Check the counts
table(seurat_in_vitro_HSC$celltype) # HSC: 9421
seurat_in_vitro_HSC # 33270 genes and 9421 cells 
DefaultAssay(seurat_in_vitro_HSC) <- "RNA"
seurat_in_vitro_HSC$donor_age <- ifelse(grepl("^O", seurat_in_vitro_HSC$sampleName), "O", "Y")
seurat_in_vitro_HSC$donor_age <- factor(seurat_in_vitro_HSC$donor_age, levels = c("Y", "O"))
table(seurat_in_vitro_HSC$celltype)
table(seurat_in_vitro_HSC$donor_age) # Young: 5040; Old: 4381

# Remove unrelated genes (technical)
all_genes <- rownames(seurat_in_vitro_HSC)
genes_filtered <- all_genes[!grepl("^mt-", all_genes, ignore.case = TRUE) &!grepl("^Rps", all_genes) & 
                            !grepl("^Rpl", all_genes) & !grepl("^Hist", all_genes) & !grepl("^Gm[0-9]", all_genes) & !grepl("Rik$", all_genes)]
length(all_genes)
length(genes_filtered) # 13366 genes 
seurat_in_vitro_HSC <- subset(seurat_in_vitro_HSC,features = genes_filtered)
seurat_in_vitro_HSC # 26732 genes and 9421 cells
seurat_in_vitro_HSC <- filter_unrelated_genes(seurat_in_vitro_HSC)
seurat_in_vitro_HSC

##### CrossAge Exp DVG #####

# Function: takes a seurat object and compute variance tests (Brown-Forsythe, Levene, Bartlett) and Coefficient of Variation (CV)
compute_variance_and_cv <- function(
    seu,
    assay = "RNA",
    slot = "data",
    group_col = "donor_age",
    n_hvg = 2000,
    selection.method = "vst",
    n_cores = 1,
    chunk_size = 250,
    tests = c("brown_forsythe", "levene", "bartlett"),
    min_pct = 0.1,
    min_mean = 0.1
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(purrr)
    library(Matrix)
    library(progressr)
    library(tidyr)
  })
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(group_col %in% colnames(seu@meta.data))
  groups_raw <- seu[[group_col]][, 1]
  groups_raw <- groups_raw[!is.na(groups_raw)]
  if (is.factor(groups_raw)) {
    groups <- levels(droplevels(groups_raw))
  } else {
    groups <- unique(as.character(groups_raw))
  }
  if (length(groups) != 2) {
    stop("group_col must contain exactly TWO groups.")
  }
  g1 <- groups[1]
  g2 <- groups[2]
  message("Comparing groups: ", g1, " vs ", g2)
  # ---------------- HVG selection ----------------
  hvgs <- VariableFeatures(seu)
  if (length(hvgs) == 0) {
    seu <- FindVariableFeatures(
      seu,
      assay = assay,
      selection.method = selection.method,
      nfeatures = n_hvg
    )
    hvgs <- VariableFeatures(seu)
  } else if (length(hvgs) > n_hvg) {
    hvgs <- hvgs[1:n_hvg]
  }
  expr <- tryCatch(
    GetAssayData(seu, assay = assay, layer = slot),
    error = function(e) GetAssayData(seu, assay = assay, slot = slot)
  )
  expr <- expr[intersect(rownames(expr), hvgs), , drop = FALSE]
  group_factor <- factor(seu[[group_col]][, 1], levels = groups)
  detect_pct <- Matrix::rowMeans(expr > 0)
  mean_expr  <- Matrix::rowMeans(expr)
  keep_genes <- names(which(detect_pct >= min_pct & mean_expr >= min_mean))
  expr <- expr[keep_genes, , drop = FALSE]
  message("Keeping ", length(keep_genes), " genes after filtering.")
  # ---------------- Variance tests ----------------
  do_tests <- function(x, group_factor) {
    res <- list(bf = NA, lev = NA, bart = NA)
    df <- data.frame(expr = x, grp = group_factor)
    if (length(unique(df$grp)) < 2) return(res)
    # Brown–Forsythe
    med <- tapply(df$expr, df$grp, median)
    dev <- abs(df$expr - med[df$grp])
    fit <- try(aov(dev ~ grp, data = df), silent = TRUE)
    if (!inherits(fit, "try-error"))
      res$bf <- summary(fit)[[1]]["grp", "Pr(>F)"]
    # Levene (mean-centered)
    mu <- tapply(df$expr, df$grp, mean)
    dev <- abs(df$expr - mu[df$grp])
    fit <- try(aov(dev ~ grp, data = df), silent = TRUE)
    if (!inherits(fit, "try-error"))
      res$lev <- summary(fit)[[1]]["grp", "Pr(>F)"]
    bt <- try(bartlett.test(expr ~ grp, data = df), silent = TRUE)
    if (!inherits(bt, "try-error"))
      res$bart <- bt$p.value
    res
  }
  message("Running variance tests ...")
  handlers(global = TRUE)
  chunks <- split(rownames(expr),
                  ceiling(seq_along(rownames(expr)) / chunk_size))
  results_list <- list()
  with_progress({
    p <- progressor(steps = length(chunks))
    for (i in seq_along(chunks)) {
      p(sprintf("Chunk %d / %d", i, length(chunks)))
      genes <- chunks[[i]]
      sub_expr <- expr[genes, , drop = FALSE]
      mat <- apply(sub_expr, 1, do_tests, group_factor = group_factor)
      bf_p   <- sapply(mat, `[[`, "bf")
      lev_p  <- sapply(mat, `[[`, "lev")
      bart_p <- sapply(mat, `[[`, "bart")
      var_tbl <- apply(sub_expr, 1, function(x) {
        df <- data.frame(expr = x, grp = group_factor)
        var_g1 <- var(df$expr[df$grp == g1], na.rm = TRUE)
        var_g2 <- var(df$expr[df$grp == g2], na.rm = TRUE)
        log2FC_variance <- ifelse(var_g1 > 0,
                                  log2(var_g2 / var_g1),
                                  NA_real_)
        c(var_g1, var_g2, log2FC_variance)
      })
      var_tbl <- as.data.frame(t(var_tbl))
      colnames(var_tbl) <- c(
        paste0("var_", g1),
        paste0("var_", g2),
        "log2FC_variance"
      )
      results_list[[i]] <- data.frame(
        gene = genes,
        p_brown_forsythe = bf_p,
        p_levene = lev_p,
        p_bartlett = bart_p,
        var_tbl
      )
    }
  })
  var_results <- bind_rows(results_list) %>%
    mutate(
      fdr_brown_forsythe = p.adjust(p_brown_forsythe, "BH"),
      fdr_levene = p.adjust(p_levene, "BH"),
      fdr_bartlett = p.adjust(p_bartlett, "BH")
    )
  # ---------------- Mean-adjusted variance ----------------
  message("Computing mean-adjusted variance and SD ...")
  mean_g1 <- apply(expr[, group_factor == g1], 1, mean, na.rm = TRUE)
  mean_g2 <- apply(expr[, group_factor == g2], 1, mean, na.rm = TRUE)
  mean_var_df <- data.frame(
    gene = rownames(expr),
    group = rep(c(g1, g2), each = nrow(expr)),
    mean = c(mean_g1, mean_g2),
    var = c(
      var_results[[paste0("var_", g1)]],
      var_results[[paste0("var_", g2)]]
    )
  )
  mean_var_df <- mean_var_df %>%
    filter(!is.na(mean) & mean > 0 & !is.na(var) & var > 0)
  fit_loess <- loess(log10(var) ~ log10(mean),
                     data = mean_var_df,
                     span = 0.75)
  mean_var_df$expected_logVar <- predict(fit_loess, newdata = mean_var_df)
  mean_var_df$residual_logVar <- log10(mean_var_df$var) -
    mean_var_df$expected_logVar
  mean_var_df$mean_adjusted_var <- 10 ^ mean_var_df$residual_logVar
  mean_var_df$mean_adjusted_sd  <- sqrt(mean_var_df$mean_adjusted_var)
  adj_tbl <- mean_var_df %>%
    dplyr::select(gene, group, mean_adjusted_var, mean_adjusted_sd) %>%
    pivot_wider(names_from = group,
                values_from = c(mean_adjusted_var,
                                mean_adjusted_sd))
  adj_tbl <- adj_tbl %>%
    mutate(
      log2FC_mean_adjusted_variance =
        log2((.data[[paste0("mean_adjusted_var_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_var_", g1)]] + 1e-8)),
      log2FC_mean_adjusted_SD =
        log2((.data[[paste0("mean_adjusted_sd_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_sd_", g1)]] + 1e-8))
    )
  # ---------------- CV ----------------
  expr_full <- GetAssayData(seu, assay = assay, slot = slot)
  pooled_mean <- Matrix::rowMeans(expr_full, na.rm = TRUE)
  cv_list <- purrr::map(groups, function(g) {
    cells <- colnames(seu)[seu[[group_col]][, 1] == g]
    mat <- expr_full[, cells, drop = FALSE]
    sd_vals <- apply(mat, 1, sd, na.rm = TRUE)
    sd_vals / (pooled_mean + 1e-8)
  })
  cv_tbl <- as.data.frame(cv_list)
  colnames(cv_tbl) <- paste0("CV_", groups)
  cv_tbl$gene <- rownames(expr_full)
  cv_tbl <- cv_tbl %>%
    mutate(
      log2_CV_ratio =
        log2((.data[[paste0("CV_", g2)]] + 1e-8) /
               (.data[[paste0("CV_", g1)]] + 1e-8))
    )
  merged <- var_results %>%
    inner_join(cv_tbl, by = "gene") %>%
    left_join(adj_tbl, by = "gene")
  message("Done. Generated table with variance tests, CV, and mean-adjusted variance metrics.")
  return(merged)
}

#### Run the function for Day0 YvsO DVGs #### 

Day0_hsc_dvgs_results <- compute_variance_and_cv(
  seu = seurat_in_vitro_HSC,
  group_col = "donor_age",
  assay = "RNA",
  slot = "data",
  n_hvg = 2000,
  min_pct = 0.1,
  min_mean = 0.1
) # 789 genes after filtering 

View(Day0_hsc_dvgs_results)

# save the Crossage Day0 DVGs table 
write.xlsx(Day0_hsc_dvgs_results,file = "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure3/Tables/Cross_vitro_HSC_DVG_results.xlsx",rowNames = FALSE)

#########

#### Function for generating the plots ####

plot_variance_cv_summaries <- function(
    df,
    fdr_col = "fdr_brown_forsythe",
    fc_col = "log2FC_mean_adjusted_variance",
    fc_col_var = "log2FC_variance",
    fdr_cutoff = 0.05,
    fc_cutoff = 0.25,
    color_palette = c("#2F7BAA", "#E58B1C"),
    prefix = "Day0",
    group_order = NULL   # optional parameter
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggsignif)
    library(ggrepel)
    library(gt)
    library(purrr)
  })
  # ---------------------------------------------------------
  # Detect group names from variance columns
  # ---------------------------------------------------------
  var_cols <- grep("^var_", colnames(df), value = TRUE)
  groups_detected <- gsub("^var_", "", var_cols)
  if (length(groups_detected) != 2)
    stop("Exactly two groups must be present.")
  if (is.null(group_order)) {
    g1 <- groups_detected[1]
    g2 <- groups_detected[2]
  } else {
    if (!all(group_order %in% groups_detected))
      stop("group_order must match detected groups: ",
           paste(groups_detected, collapse = ", "))
    g1 <- group_order[1]
    g2 <- group_order[2]
  }
  message("Plotting comparison: ", g2, " vs ", g1)
  get_wilcox_label <- function(p)
    if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
  # ==========================================================
  # Variance Boxplot
  # ==========================================================
  long_var <- df %>%
    dplyr::select(gene, all_of(var_cols)) %>%
    pivot_longer(-gene,
                 names_to = "Group",
                 values_to = "Variance") %>%
    mutate(
      Group = gsub("^var_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  p_wilcox <- wilcox.test(Variance ~ Group, data = long_var)$p.value
  p_var <- ggplot(long_var, aes(Group, Variance, color = Group)) +
    geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, linewidth = 1) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.6) +
    geom_signif(comparisons = list(c(g1, g2)),
                annotations = get_wilcox_label(p_wilcox),
                y_position = max(long_var$Variance, na.rm = TRUE) * 1.05) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    labs(title = sprintf("%s: Per-Gene Variance", prefix),
         subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                            p_wilcox, nrow(df)),
         y = "Variance", x = NULL)
  # ==========================================================
  # CV Boxplot
  # ==========================================================
  cv_cols <- grep("^CV_", colnames(df), value = TRUE)
  long_cv <- df %>%
    dplyr::select(gene, all_of(cv_cols)) %>%
    pivot_longer(-gene,
                 names_to = "Group",
                 values_to = "CV") %>%
    mutate(
      Group = gsub("^CV_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  p_cv <- wilcox.test(CV ~ Group, data = long_cv)$p.value
  p_cv_box <- ggplot(long_cv, aes(Group, CV, color = Group)) +
    geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, linewidth = 1) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.6) +
    geom_signif(comparisons = list(c(g1, g2)),
                annotations = get_wilcox_label(p_cv),
                y_position = max(long_cv$CV, na.rm = TRUE) * 1.05) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    labs(title = sprintf("%s: Coefficient of Variation", prefix),
         subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                            p_cv, nrow(df)),
         y = "Coefficient of Variation", x = NULL)
  # ==========================================================
  # Mean-adjusted Variance Boxplot
  # ==========================================================
  mvar_cols <- grep("^mean_adjusted_var_", colnames(df), value = TRUE)
  long_mvar <- df %>%
    dplyr::select(gene, all_of(mvar_cols)) %>%
    pivot_longer(-gene,
                 names_to = "Group",
                 values_to = "MeanAdjustedVar") %>%
    mutate(
      Group = gsub("^mean_adjusted_var_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  p_mv <- wilcox.test(MeanAdjustedVar ~ Group, data = long_mvar)$p.value
  p_mv_box <- ggplot(long_mvar, aes(Group, MeanAdjustedVar, color = Group)) +
    geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, linewidth = 1) +
    geom_jitter(aes(fill = Group),
                width = 0.15,
                alpha = 0.6,
                size = 2,
                shape = 21,
                stroke = 0.4,
                color = "black") +
    geom_signif(comparisons = list(c(g1, g2)),
                annotations = get_wilcox_label(p_mv),
                y_position = max(long_mvar$MeanAdjustedVar, na.rm = TRUE) * 1.05) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) + 
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    labs(title = sprintf("%s: Mean-Adjusted Variance", prefix),
         subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                            p_mv, nrow(df)),
         y = "Mean-Adj Var", x = NULL)
  # ==========================================================
  # Density Plot
  # ==========================================================
  p_density <- ggplot(long_mvar,
                      aes(x = MeanAdjustedVar, color = Group)) +
    geom_density(linewidth = 1.2, adjust = 1.1) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    labs(title = sprintf("%s: Distribution of Mean-Adjusted Variance", prefix),
         x = "Mean-Adj Var",
         y = "Density")
  
  # ==========================================================
  # Density Plot
  # ==========================================================
  
  df_sig <- df %>%
    filter(
      .data[[fdr_col]] < fdr_cutoff &
        abs(.data[[fc_col]]) > fc_cutoff
    )
  
  long_mvar_sig <- df_sig %>%
    dplyr::select(gene, all_of(mvar_cols)) %>%
    pivot_longer(
      -gene,
      names_to = "Group",
      values_to = "MeanAdjustedVar"
    ) %>%
    mutate(
      Group = gsub("^mean_adjusted_var_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  
  p_density_sig <- ggplot(long_mvar_sig,
                          aes(x = MeanAdjustedVar, color = Group)) +
    geom_density(linewidth = 1.2, adjust = 1.1) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    labs(
      #title = sprintf("%s: Distribution of Mean-Adjusted Variance (Sig)", prefix),
      #subtitle = sprintf("n = %d significant genes", nrow(df_sig)),
      x = "Mean-Adj Var",
      y = "Density"
    )
  
  # ==========================================================
  # Summary Counts 
  # ==========================================================
  fc_thresholds <- c(0, 0.1, 0.25, 0.5)
  summary_counts_mvar <- purrr::map_dfr(fc_thresholds, function(fc_thr) {
    df %>%
      filter(.data[[fdr_col]] < 0.5) %>%
      summarise(
        log2FC_cutoff = fc_thr,
        up_in_g2 = sum(.data[[fc_col]] >  fc_thr, na.rm = TRUE),
        up_in_g1 = sum(.data[[fc_col]] < -fc_thr, na.rm = TRUE)
      ) %>%
      mutate(metric = "Mean-Adjusted Variance")
  })
  summary_counts_var <- purrr::map_dfr(fc_thresholds, function(fc_thr) {
    df %>%
      filter(.data[[fdr_col]] < 0.5) %>%
      summarise(
        log2FC_cutoff = fc_thr,
        up_in_g2 = sum(.data[[fc_col_var]] >  fc_thr, na.rm = TRUE),
        up_in_g1 = sum(.data[[fc_col_var]] < -fc_thr, na.rm = TRUE)
      ) %>%
      mutate(metric = "Raw Variance")
  })
  summary_counts <- bind_rows(summary_counts_mvar, summary_counts_var)
  # ==========================================================
  # Volcano (Mean-adjusted)
  # ==========================================================
  
  df_vol_mvar <- df %>%
    mutate(
      neg_log10_fdr = -log10(pmax(.data[[fdr_col]], 1e-300)),
      significance = case_when(
        .data[[fdr_col]] < fdr_cutoff & .data[[fc_col]] >  fc_cutoff ~ paste0("↑", g2),
        .data[[fdr_col]] < fdr_cutoff & .data[[fc_col]] < -fc_cutoff ~ paste0("↑", g1),
        TRUE ~ "Not significant"
      )
    )
  
  # Top genes increased in g2
  top_up_g2 <- df_vol_mvar %>%
    filter(significance == paste0("↑", g2)) %>%
    arrange(.data[[fdr_col]], desc(.data[[fc_col]])) %>%   
    slice_head(n = 10)
  
  # Top genes increased in g1
  top_up_g1 <- df_vol_mvar %>%
    filter(significance == paste0("↑", g1)) %>%
    arrange(.data[[fdr_col]], .data[[fc_col]]) %>%         
    slice_head(n = 10)
  
  top_genes_mvar <- bind_rows(top_up_g2, top_up_g1) %>%
    mutate(gene = paste0("italic('", gene, "')"))
  
  volcano_colors <- c("#5271AE", "#D85B59", "grey80")
  names(volcano_colors) <- c(
    paste0("↑", g1),
    paste0("↑", g2),
    "Not significant"
  )
  # Significance threshold
  sig_line <- -log10(fdr_cutoff)
  
  # Add some padding room
  y_min <- sig_line - 10
  y_max <- max(df_vol_mvar$neg_log10_fdr, na.rm = TRUE)
  

  p_volcano_mvar <- ggplot(df_vol_mvar,
                           aes(x = .data[[fc_col]],
                               y = neg_log10_fdr)) +
    geom_point(aes(color = significance),
               alpha = 0.85,
               size = 3) +
    geom_text_repel(
      data = top_genes_mvar,
      aes(label = gene),
      parse = TRUE, 
      size = 5,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.4,
      force = 2,
      segment.size = 0.3
    ) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed") +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    xlim(-1.5, 1.5) +
    theme_classic(base_size = 20) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black"),
      legend.position = "none"
    ) +
    labs(
      x = paste0("log2 Mean-Adj Var (", g2, " / ", g1, ")"),
      y = expression(bold(-log[10](FDR)))
    )
  # -------------------------------------------------------------------------
  # Volcano Plot 2: Raw variance
  # -------------------------------------------------------------------------
  
  df_vol_var <- df %>%
    mutate(
      neg_log10_fdr = -log10(pmax(.data[[fdr_col]], 1e-300)),
      significance_var = case_when(
        .data[[fdr_col]] < fdr_cutoff & .data[[fc_col_var]] >  fc_cutoff ~ paste0("↑", g2),
        .data[[fdr_col]] < fdr_cutoff & .data[[fc_col_var]] < -fc_cutoff ~ paste0("↑", g1),
        TRUE ~ "Not significant"
      )
    )
  
  top_up_g2_var <- df_vol_var %>%
    filter(significance_var == paste0("↑", g2)) %>%
    arrange(.data[[fdr_col]], desc(.data[[fc_col]])) %>% 
    slice_head(n = 10)
  
  top_up_g1_var <- df_vol_var %>%
    filter(significance_var == paste0("↑", g1)) %>%
    arrange(.data[[fdr_col]], .data[[fc_col]]) %>%
    slice_head(n = 10)
  
  top_genes_var <- bind_rows(top_up_g2_var, top_up_g1_var) %>%
    mutate(gene = paste0("italic('", gene, "')"))
  
  volcano_colors <- c("#5271AE", "#D85B59", "grey80")
  names(volcano_colors) <- c(
    paste0("↑", g1),
    paste0("↑", g2),
    "Not significant"
  )
  y_min <- sig_line - 10
  y_max <- max(df_vol_mvar$neg_log10_fdr)
  
  p_volcano_var <- ggplot(df_vol_var,
                          aes(x = .data[[fc_col_var]],
                              y = neg_log10_fdr)) +
    geom_point(aes(color = significance_var),
               alpha = 0.85,
               size = 3) +
    geom_text_repel(
      data = top_genes_var,
      aes(label = gene),
      parse = TRUE, 
      size = 5,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.4,
      force = 2,
      segment.size = 0.3
    ) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed") +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    xlim(-3.5, 3.5) +
    theme_classic(base_size = 20) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      x = paste0("log2 Variance (", g2, "/", g1, ")"),
      y = expression(bold(-log[10](FDR)))
    )
  # ==========================================================
  # Return
  # ==========================================================
  list(
    variance_box = p_var,
    cv_box = p_cv_box,
    mean_adj_var_box = p_mv_box,
    density_mean_adj = p_density,
    density_mean_adj_sig = p_density_sig,
    volcano_mean_adj = p_volcano_mvar,
    volcano_variance = p_volcano_var,
    summary_counts = summary_counts
  )
}

####

day0_hsc_plots <- plot_variance_cv_summaries(
  df = Day0_hsc_dvgs_results,
  fc_cutoff = 0.1,
  group_order = c("Y","O"),
  prefix = "HSC (Day 0)"
)
day0_hsc_plots$summary_counts
day0_hsc_plots$mean_adj_var_box
day0_hsc_plots$density_mean_adj
day0_hsc_plots$density_mean_adj_sig
day0_hsc_plots$volcano_mean_adj
day0_hsc_plots$volcano_variance


###########################################################


##########################################

## Compute DVGs (Y vs O) in unmanipulated data

##########################################

# function 8.2: compute mean adjusted variance for unmanipulated data (multiple agegroups)

compute_variance_and_cv_multigroup_with_mean_adjustment <- function(
    seu,
    assay = "RNA",
    slot = "data",
    group_col = "AGE",
    n_hvg = 2000,
    selection.method = "vst",
    min_pct = 0.1,
    min_mean = 0.1,
    chunk_size = 250
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(progressr)
  })
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(group_col %in% colnames(seu@meta.data))
  seu[[group_col]][, 1] <- factor(
    seu[[group_col]][, 1],
    levels = c("young", "mid", "old", "vold"),
    ordered = TRUE
  )
  age <- seu[[group_col]][, 1]
  message("Age levels detected: ", paste(levels(age), collapse = ", "))
  # --- HVGs ---
  hvgs <- VariableFeatures(seu)
  if (length(hvgs) == 0) {
    seu <- FindVariableFeatures(seu, assay = assay, selection.method = selection.method, nfeatures = n_hvg)
    hvgs <- VariableFeatures(seu)
  } else if (length(hvgs) > n_hvg) hvgs <- hvgs[1:n_hvg]
  expr <- tryCatch(
    GetAssayData(seu, assay = assay, layer = slot),
    error = function(e) GetAssayData(seu, assay = assay, slot = slot)
  )
  expr <- expr[intersect(rownames(expr), hvgs), , drop = FALSE]
  # --- Filtering ---
  detect_pct <- Matrix::rowMeans(expr > 0)
  mean_expr <- Matrix::rowMeans(expr)
  keep_genes <- names(which(detect_pct >= min_pct & mean_expr >= min_mean))
  expr <- expr[keep_genes, , drop = FALSE]
  message("Keeping ", length(keep_genes), " genes after filtering.")
  # --- Variance tests (Brown–Forsythe, Levene, Bartlett) ---
  handlers(global = TRUE)
  message("Running variance tests ...")
  do_tests <- function(x, age) {
    res <- list(bf = NA, lev = NA, bart = NA)
    df <- data.frame(expr = x, age = age)
    if (length(unique(age)) < 2) return(res)
    # Brown–Forsythe
    med <- tapply(df$expr, df$age, median, na.rm = TRUE)
    dev <- abs(df$expr - med[df$age])
    fit <- try(aov(dev ~ age, data = df), silent = TRUE)
    if (!inherits(fit, "try-error"))
      res$bf <- summary(fit)[[1]]["age", "Pr(>F)"]
    # Levene
    mu <- tapply(df$expr, df$age, mean, na.rm = TRUE)
    dev <- abs(df$expr - mu[df$age])
    fit <- try(aov(dev ~ age, data = df), silent = TRUE)
    if (!inherits(fit, "try-error"))
      res$lev <- summary(fit)[[1]]["age", "Pr(>F)"]
    # Bartlett
    bt <- try(bartlett.test(expr ~ age, data = df), silent = TRUE)
    if (!inherits(bt, "try-error"))
      res$bart <- bt$p.value
    res
  }
  # Chunking for memory efficiency
  chunks <- split(rownames(expr), ceiling(seq_along(rownames(expr)) / chunk_size))
  results_list <- list()
  with_progress({
    p <- progressor(steps = length(chunks))
    for (i in seq_along(chunks)) {
      p(sprintf("Chunk %d / %d", i, length(chunks)))
      genes <- chunks[[i]]
      sub_expr <- expr[genes, , drop = FALSE]
      mat <- apply(sub_expr, 1, do_tests, age = age)
      bf_p   <- sapply(mat, `[[`, "bf")
      lev_p  <- sapply(mat, `[[`, "lev")
      bart_p <- sapply(mat, `[[`, "bart")
      results_list[[i]] <- data.frame(
        gene = genes,
        p_brown_forsythe = bf_p,
        p_levene = lev_p,
        p_bartlett = bart_p
      )
      gc()
    }
  })
  var_test_results <- bind_rows(results_list) %>%
    mutate(
      fdr_brown_forsythe = p.adjust(p_brown_forsythe, "BH"),
      fdr_levene = p.adjust(p_levene, "BH"),
      fdr_bartlett = p.adjust(p_bartlett, "BH")
    )
  # --- Group-wise mean, variance, and CV ---
  message("Computing group-wise mean, variance, and CV ...")
  mean_tbl <- apply(expr, 1, function(x) tapply(x, age, mean, na.rm = TRUE)) %>% t() %>% as.data.frame()
  var_tbl  <- apply(expr, 1, function(x) tapply(x, age, var,  na.rm = TRUE)) %>% t() %>% as.data.frame()
  cv_tbl   <- apply(expr, 1, function(x) tapply(x, age, function(v) sd(v, na.rm = TRUE) / (mean(v, na.rm = TRUE) + 1e-8))) %>%
    t() %>% as.data.frame()
  colnames(mean_tbl) <- paste0("mean_", colnames(mean_tbl))
  colnames(var_tbl)  <- paste0("var_",  colnames(var_tbl))
  colnames(cv_tbl)   <- paste0("CV_",   colnames(cv_tbl))
  gene_stats <- bind_cols(
    gene = rownames(expr),
    mean_tbl,
    var_tbl,
    cv_tbl
  )
  # --- Mean–variance adjustment ---
  long_stats <- gene_stats %>%
    pivot_longer(cols = starts_with("var_"), names_to = "var_group", values_to = "var") %>%
    mutate(mean = rep(as.numeric(t(mean_tbl)), each = 1)) %>%
    mutate(log_mean = log10(mean + 1e-8), log_var = log10(var + 1e-8)) %>%
    filter(is.finite(log_mean) & is.finite(log_var))
  fit_loess <- loess(log_var ~ log_mean, data = long_stats, span = 0.75)
  pred <- predict(fit_loess, newdata = long_stats$log_mean)
  long_stats$residual_logVar <- long_stats$log_var - pred
  long_stats$residual_SD <- sqrt(10 ^ long_stats$residual_logVar)
  resid_summary <- long_stats %>%
    mutate(group = sub("var_", "", var_group)) %>%
    group_by(gene, group) %>%
    summarise(
      mean_adjusted_var = mean(10 ^ residual_logVar),
      mean_adjusted_sd  = mean(residual_SD),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = group,
      values_from = c(mean_adjusted_var, mean_adjusted_sd),
      names_glue = "{.value}_{group}"
    )
  # --- Merge all ---
  merged <- gene_stats %>%
    left_join(resid_summary, by = "gene") %>%
    left_join(var_test_results, by = "gene") %>%
    mutate(
      # log2FCs
      log2_var_ratio_old_young  = log2((var_old  + 1e-8) / (var_young  + 1e-8)),
      log2_var_ratio_vold_young = log2((var_vold + 1e-8) / (var_young  + 1e-8)),
      log2_CV_ratio_old_young   = log2((CV_old   + 1e-8) / (CV_young   + 1e-8)),
      log2_CV_ratio_vold_young  = log2((CV_vold  + 1e-8) / (CV_young   + 1e-8)),
      log2FC_mean_adjusted_var_old_young  = log2((mean_adjusted_var_old  + 1e-8) / (mean_adjusted_var_young  + 1e-8)),
      log2FC_mean_adjusted_var_vold_young = log2((mean_adjusted_var_vold + 1e-8) / (mean_adjusted_var_young + 1e-8))
    )
  message("Done. Computed variance tests, CV, and mean-adjusted variance with FDR corrections.")
  return(merged)
}
### Plotting function ###

plot_variance_cv_summaries_multigroup_with_mean_adjustment <- function(
    df,
    fdr_col = "fdr_brown_forsythe",
    fdr_cutoff = 0.05,
    fc_cutoff = 0.25,
    color_palette = c("#D95A9A", "#4A9EB0", "#2F7BAA", "#E58B1C"),
    prefix = "Day0 HSC"
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggrepel)
    library(gt)
    library(ggsignif)
  })
  age_levels <- c("young", "mid", "old", "vold")
  n_genes <- nrow(df)
  wilcox_pairs <- combn(age_levels, 2, simplify = FALSE)
  # Helper to format significance stars
  signif_label <- function(p) {
    if (is.na(p)) return("ns")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else "ns"
  }
  # ==== Boxplot: Variance ====
  long_var <- df %>%
    select(gene, starts_with("var_")) %>%
    pivot_longer(-gene, names_to = "Group", values_to = "Variance") %>%
    mutate(Group = sub("var_", "", Group),
           Group = factor(Group, levels = age_levels, ordered = TRUE))
  # Compare young vs old
  if (all(c("var_young", "var_old") %in% colnames(df))) {
    p_var_test <- wilcox.test(df$var_young, df$var_old)$p.value
  } else {
    p_var_test <- NA
  }
  p_var <- ggplot(long_var, aes(x = Group, y = Variance, color = Group)) +
    geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.9) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    geom_signif(
      comparisons = wilcox_pairs,
      test = "wilcox.test",
      test.args = list(paired = FALSE),
      step_increase = 0.1,
      textsize = 3.5
    ) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("%s: Variance Across Ages", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes", p_var_test, n_genes),
      y = "Variance"
    )
  # ==== Boxplot: CV ====
  long_cv <- df %>%
    select(gene, starts_with("CV_")) %>%
    pivot_longer(-gene, names_to = "Group", values_to = "CV") %>%
    mutate(Group = sub("CV_", "", Group),
           Group = factor(Group, levels = age_levels, ordered = TRUE))
  if (all(c("CV_young", "CV_old") %in% colnames(df))) {
    p_cv_test <- wilcox.test(df$CV_young, df$CV_old)$p.value
  } else {
    p_cv_test <- NA
  }
  p_cv <- ggplot(long_cv, aes(x = Group, y = CV, color = Group)) +
    geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.9) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    geom_signif(
      comparisons = wilcox_pairs,
      test = "wilcox.test",
      test.args = list(paired = FALSE),
      step_increase = 0.1,
      textsize = 3.5
    ) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("%s: CV Across Ages", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes", p_cv_test, n_genes),
      y = "Coefficient of Variation"
    )
  # ==== Boxplot: Mean-adjusted variance (Y vs O) ====
  long_mvar <- df %>%
    select(gene, mean_adjusted_var_young, mean_adjusted_var_old) %>%
    pivot_longer(-gene, names_to = "Group", values_to = "MeanAdjustedVar") %>%
    mutate(Group = recode(Group,
                          mean_adjusted_var_young = "Young (Y)",
                          mean_adjusted_var_old = "Old (O)"))
  p_mvar_test <- wilcox.test(
    df$mean_adjusted_var_young,
    df$mean_adjusted_var_old
  )$p.value
  p_label_mvar <- signif_label(p_mvar_test)
  y_max_mv <- max(long_mvar$MeanAdjustedVar, na.rm = TRUE)
  p_mv_box <- ggplot(long_mvar, aes(x = Group, y = MeanAdjustedVar, color = Group)) +
    geom_boxplot(outlier.shape = NA, fill = "white", linewidth = 0.9) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    geom_signif(
      comparisons = list(c("Young (Y)", "Old (O)")),
      annotations = p_label_mvar,
      y_position = y_max_mv * 1.05,
      color = "black"
    ) +
    scale_color_manual(values = c("#67A5CC", "#D780AA")) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("%s: Mean-Adjusted Variance", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes", p_mvar_test, n_genes),
      y = "Mean-Adjusted Variance"
    )
  # ==== Density: Mean-adjusted variance ====
  p_density <- ggplot(long_mvar, aes(x = MeanAdjustedVar, color = Group)) +
    geom_density(linewidth = 1.2, adjust = 1.1) +
    scale_color_manual(values = c("#67A5CC", "#D780AA")) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = sprintf("%s: Mean-Adjusted Variance Distribution", prefix),
      x = "Mean-Adjusted Variance", y = "Density"
    )
  # ==== Volcano plots ====
  make_volcano <- function(df, fc_col, label) {
    df_plot <- df %>%
      mutate(
        neg_log10_fdr = -log10(pmax(.data[[fdr_col]], 1e-300)),
        significance = case_when(
          .data[[fdr_col]] < fdr_cutoff & .data[[fc_col]] >  fc_cutoff ~ "↑O (Higher in Old)",
          .data[[fdr_col]] < fdr_cutoff & .data[[fc_col]] < -fc_cutoff ~ "↑Y (Higher in Young)",
          TRUE ~ "Not significant"
        )
      )
    n_total <- nrow(df_plot)
    n_upO <- sum(df_plot$significance == "↑O (Higher in Old)", na.rm = TRUE)
    n_upY <- sum(df_plot$significance == "↑Y (Higher in Young)", na.rm = TRUE)
    n_sig <- n_upO + n_upY
    # ---- Label top10 up and down ----
    top_up <- df_plot %>%
      filter(significance == "↑O (Higher in Old)") %>%
      arrange(.data[[fdr_col]], desc(.data[[fc_col]])) %>%
      slice_head(n = 10)
    top_down <- df_plot %>%
      filter(significance == "↑Y (Higher in Young)") %>%
      arrange(.data[[fdr_col]], .data[[fc_col]]) %>%
      slice_head(n = 10)
    top_genes <- bind_rows(top_up, top_down) %>%
      mutate(gene_label = paste0("italic('", gene, "')"))
    ggplot(df_plot, aes(x = .data[[fc_col]], y = neg_log10_fdr)) +
      geom_point(aes(color = significance), alpha = 0.85, size = 3) +
      geom_text_repel(
        data = top_genes,
        aes(label = gene_label),
        parse = TRUE,
        color = "black", size = 5,
        box.padding = 0.4,
        point.padding = 0.4, force = 3,
        max.overlaps = Inf,
        segment.size = 0.3
      ) +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
      geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
      scale_color_manual(values = c(
        "↑O (Higher in Old)" = "#D85B59",
        "↑Y (Higher in Young)" = "#5271AE",
        "Not significant" = "grey80"
      )) +
      coord_cartesian(xlim = c(-0.5, 0.5)) +
      theme_classic(base_size = 20) +
      theme(
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text  = element_text(color = "black"),
        legend.position = "none"
      ) +
      labs(
        #title = sprintf("%s: Volcano Plot (%s)", prefix, label),
        #subtitle = sprintf(
        #  "n = %d | %d significant | ↑O: %d | ↑Y: %d | FDR < %.2f & |log₂FC| > %.2f",
        #  n_total, n_sig, n_upO, n_upY, fdr_cutoff, fc_cutoff
        #),
        x = expression(bold(log[2]("Mean-Adj Var (O/Y)"))),
        y = expression(bold(-log[10](FDR)))
      )
  }
  p_volcano_old  <- make_volcano(df, "log2FC_mean_adjusted_var_old_young",  "Old vs Young")
  p_volcano_vold <- make_volcano(df, "log2FC_mean_adjusted_var_vold_young", "Vold vs Young")
  # ==== Summary Table ====
  summary_tbl <- df %>%
    mutate(direction = case_when(
      .data[[fdr_col]] < fdr_cutoff & .data[["log2FC_mean_adjusted_var_old_young"]]  >  fc_cutoff ~ "↑O (Higher in Old)",
      .data[[fdr_col]] < fdr_cutoff & .data[["log2FC_mean_adjusted_var_old_young"]]  < -fc_cutoff ~ "↑Y (Higher in Young)",
      TRUE ~ "No significant change"
    )) %>%
    dplyr::count(direction, name = "n_genes") %>%
    mutate(Threshold = sprintf("FDR<%.2f & |log₂FC|>%.2f", fdr_cutoff, fc_cutoff))
  gt_summary <- summary_tbl %>%
    gt(rowname_col = "direction", groupname_col = "Threshold") %>%
    fmt_number(columns = n_genes, decimals = 0) %>%
    cols_label(n_genes = "Gene Count") %>%
    tab_header(
      title = md(sprintf("**%s: Summary of Mean-Adjusted DVGs**", prefix)),
      subtitle = md("Counts of genes with higher mean-adjusted variance in Old vs Young")
    )
  # ==== Return ====
  list(
    variance_box = p_var,
    cv_box = p_cv,
    mean_adj_box = p_mv_box,
    mean_adj_density = p_density,
    volcano_old = p_volcano_old,
    volcano_vold = p_volcano_vold,
    summary_table = gt_summary
  )
}

##### Unmanipulated HSC ######

seurat_unmanipulated
seurat_unmanipulated$AGE <- factor(seurat_unmanipulated$AGE,levels = c("young", "mid", "old", "vold"),ordered = TRUE)
table(seurat_unmanipulated$AGE)
table(seurat_unmanipulated$celltype)
seurat_unmanipulated_HSC <- subset(seurat_unmanipulated,subset = celltype == c("LT-HSC","ST-HSC"))
seurat_unmanipulated_HSC <- filter_unrelated_genes(seurat_unmanipulated_HSC) # filter out unrelated/ technical genes
seurat_unmanipulated_HSC # 14323 genes and 2797 cells 
table(seurat_unmanipulated_HSC$celltype)
table(seurat_unmanipulated_HSC$AGE)


#### Run the two functions ####

unmanipulated_results <- compute_variance_and_cv_multigroup_with_mean_adjustment(
  seu = seurat_unmanipulated_HSC,
  group_col = "AGE",
  assay = "RNA",
  slot = "data",
  min_pct = 0.1,
  min_mean = 0.1
) # 687 genes 

plots_unmanipulated <- plot_variance_cv_summaries_multigroup_with_mean_adjustment(
  df = unmanipulated_results,
  fdr_col = "fdr_brown_forsythe",
  fdr_cutoff = 0.05,
  fc_cutoff = 0.1,
  prefix = "Unmanipulated(HSC)"
)

plots_unmanipulated$variance_box
plots_unmanipulated$cv_box
plots_unmanipulated$mean_adj_box
plots_unmanipulated$mean_adj_density
plots_unmanipulated$volcano_old
plots_unmanipulated$volcano_vold
plots_unmanipulated$summary_table

##### Top Aging-up Unmanipylated DVGs #####

# ---- Select top 100 aging-up DVGs (higher mean-adjusted variance in old) ----
unmanipulated_HSC_results_ranked_top100 <- unmanipulated_results %>%
  filter(
    fdr_brown_forsythe < 0.05,
    log2FC_mean_adjusted_var_old_young > 0
  ) %>%
  arrange(fdr_brown_forsythe, desc(log2FC_mean_adjusted_var_old_young)) %>%
  slice_head(n = 100)

# ---- Prepare data for plotting (Young vs Old only) ----
mvar_long <- unmanipulated_HSC_results_ranked_top100 %>%
  select(mean_adjusted_var_young, mean_adjusted_var_old) %>%
  pivot_longer(
    cols = everything(),
    names_to = "AgeGroup",
    values_to = "mean_adjusted_variance"
  ) %>%
  mutate(
    AgeGroup = sub("mean_adjusted_var_", "", AgeGroup),
    AgeGroup = factor(AgeGroup, levels = c("young", "old"), ordered = TRUE)
  )

# ---- Define color palette ----

color_palette <- c(
  "young" = "cornflowerblue",
  "old"   = "orange2"
)


# ---- Wilcoxon test pairs ----
wilcox_pairs <- combn(levels(mvar_long$AgeGroup), 2, simplify = FALSE)

# ---- Boxplot ----
p_mvar_box <- ggplot(mvar_long, aes(x = AgeGroup, y = mean_adjusted_variance, color = AgeGroup)) +
  geom_boxplot(width = 0.5, fill = NA, outlier.shape = NA, linewidth = 0.9) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  geom_signif(
    comparisons = wilcox_pairs,
    test = "wilcox.test",
    test.args = list(paired = FALSE),
    step_increase = 0.1,
    textsize = 3.5
  ) +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text  = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Unmanipulated HSC: Mean-Adjusted Variance Across Ages",
    subtitle = sprintf(
      "Top 100 aging-up DVGs | n = %d genes",
      length(unique(unmanipulated_HSC_results_ranked_top100$gene))
    ),
    x = "Age Group",
    y = "Mean-Adjusted Variance"
  )
p_mvar_box

### Randomization by genes ###
# true aging-up genes
top100_genes <- unmanipulated_HSC_results_ranked_top100$gene
obs_df <- unmanipulated_results %>%
  filter(gene %in% top100_genes) %>%
  transmute(
    gene,
    young = mean_adjusted_var_young,
    old   = mean_adjusted_var_old,
    diff  = old - young
  )
obs_mean_diff <- mean(obs_df$diff, na.rm = TRUE)
obs_mean_diff # 0.1034
set.seed(42)
n_perm <- 100
all_genes <- setdiff(unmanipulated_results$gene, top100_genes)
perm_diffs <- numeric(n_perm)
for (i in 1:n_perm) {
  random_genes <- sample(all_genes, 100)
  perm_df <- unmanipulated_results %>%
    filter(gene %in% random_genes) %>%
    transmute(
      gene,
      young = mean_adjusted_var_young,
      old   = mean_adjusted_var_old,
      diff  = old - young
    )
  perm_diffs[i] <- mean(perm_df$diff, na.rm = TRUE)
}
emp_p <- (1 + sum(perm_diffs >= obs_mean_diff)) / (n_perm + 1)
emp_p
df_plot <- data.frame(
  perm_diffs = perm_diffs
)
p_null <- ggplot(df_plot, aes(x = perm_diffs)) +
  geom_histogram(bins = 40, color = "black", fill = "grey80") +
  geom_vline(xintercept = obs_mean_diff, color = "red", linewidth = 1.2) +
  theme_classic(base_size = 16) +
  labs(
    title = "Null Distribution of Old–Young Variance Difference\n(Random 100 Genes per Iteration)",
    subtitle = sprintf("Observed mean diff (red) = %.4f | Empirical p = %.4g", 
                       obs_mean_diff, emp_p),
    x = "Mean(Old − Young) variance difference",
    y = "Count"
  )
p_null

##### Compare randomized and observed DVG difference #####

# ---- Observed: per-gene values from top100
obs_plot <- unmanipulated_results %>% filter(gene %in% top100_genes) %>%
  select(
    gene,
    young = mean_adjusted_var_young,
    old   = mean_adjusted_var_old
  ) %>%
  pivot_longer(
    cols = c(young, old),
    names_to = "group",
    values_to = "mean_adj_var"
  ) %>%
  mutate(
    group = ifelse(group == "young", "Y", "O"),
    condition = "Observed"
  )
# ---- Randomized: use permutation MEAN across iterations
set.seed(42)
rand_plot <- map_dfr(1:n_perm, function(i) {
  random_genes <- sample(all_genes, 100)
  df <- unmanipulated_results %>%
    filter(gene %in% random_genes)
  tibble(
    iter = i,
    group = c("Y", "O"),
    mean_adj_var = c(
      mean(df$mean_adjusted_var_young, na.rm = TRUE),
      mean(df$mean_adjusted_var_old,  na.rm = TRUE)
    ),
    condition = "Randomized"
  )
})

### Combine observed + randomized
plot_df <- bind_rows(obs_plot, rand_plot) %>%
  mutate(
    group = recode(group, Y = "Y", O = "O"),
    group = factor(group, levels = c("Y", "O")),
    condition = factor(condition, levels = c("Observed", "Randomized"))
  )
### Final observed vs randomized boxplot

color_palette <- c("Y" = "cornflowerblue", "O" = "#E58B1C")
base_plot <- function(df) {
  ggplot(df, aes(x = group, y = mean_adj_var)) +
    geom_boxplot(
      aes(color = group),
      width = 0.55,
      fill = "white",
      outlier.shape = NA,
      linewidth = 1
    ) +
    geom_jitter(
      aes(fill = group),
      width = 0.15,
      alpha = 0.8,
      size = 3,
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_classic(base_size = 20) +
    theme(
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = "Mean-Adj Var"
    )
}
# observed 
p_cond1 <- base_plot(
  subset(plot_df, condition == unique(plot_df$condition)[1])
)
# randomized
p_cond2 <- base_plot(
  subset(plot_df, condition == unique(plot_df$condition)[2])
)
p_cond1
p_cond2


# ---- Observed Δ per gene ----
obs_delta <- unmanipulated_results %>%
  filter(gene %in% top100_genes) %>%
  transmute(
    gene,
    delta = mean_adjusted_var_old - mean_adjusted_var_young,
    condition = "Obs"
  )
set.seed(42)
perm_delta <- map_dfr(1:n_perm, function(i) {
  
  random_genes <- sample(all_genes, 100)
  
  df <- unmanipulated_results %>%
    filter(gene %in% random_genes)
  
  tibble(
    iter = i,
    delta = mean(df$mean_adjusted_var_old, na.rm = TRUE) -
      mean(df$mean_adjusted_var_young, na.rm = TRUE),
    condition = "Ran"
  )
})
delta_plot_df <- bind_rows(obs_delta, perm_delta) %>%
  mutate(condition = factor(condition, levels = c("Obs", "Ran")))
delta_plot_df 

p_delta <- ggplot(
  delta_plot_df,
  aes(x = condition, y = delta, color = condition)
) +
  geom_boxplot(
    width = 0.55,         
    fill = "white",
    outlier.shape = NA,
    linewidth = 1
  ) +
  geom_jitter(
    aes(fill = condition),
    width = 0.12,         
    alpha = 0.8,
    size = 3,
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_color_manual(
    values = c(
      "Obs" = "#51C2BB",
      "Ran" = "#F57075"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Obs"   = "#51C2BB",
      "Ran" = "#F57075"
    )
  ) +
  theme_classic(base_size = 20) +
  labs(
    x = NULL,
    y = "Δ Mean-Adj Var (O − Y)"
  ) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    plot.margin  = margin(6, 6, 6, 6),
    plot.title   = element_text(size = 22, hjust = 0.5),
    strip.text   = element_text(size = 16),
    legend.position = "none"  
  )
p_delta

## Patchwork to combined p_cond1 and p_delta
p_cond1
p_delta
p1 <- p_cond1 + facet_wrap(~"Observed") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.2),
    strip.text = element_text(face = "bold", size = 13, margin = margin(t = 6, b = 6, l = 40, r = 40)),
    axis.text.x = element_text(face = "bold")
  )

p2 <- p_delta + facet_wrap(~"Randomized") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.2),
    strip.text = element_text(face = "bold", size = 13, margin = margin(t = 6, b = 6, l = 40, r = 40)),
    axis.text.x = element_text(face = "bold")
  )
p_combined <- p1 + p2 + plot_layout(ncol = 2)
p_combined


#####################################

## Exp2 Y vs O DVGs 

#####################################

## Modify the function to compute DVGs for Exp2 

compute_variance_and_cv_YO_pairs_exp2 <- function(
    seu,
    assay = "RNA",
    slot = "data",
    group_col = "sampleName",
    n_hvg = 2000,
    selection.method = "vst",
    min_pct = 0.1,
    min_mean = 0.1
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(car)
  })
  # --- Align metadata ---
  if (!identical(colnames(seu), rownames(seu@meta.data))) {
    seu@meta.data <- seu@meta.data[colnames(seu), , drop = FALSE]
  }
  meta <- seu@meta.data
  # --- Aggregate replicates into Y vs O ---
  meta$AgeGroup <- ifelse(grepl("_Ya", meta[[group_col]]), "Y",
                          ifelse(grepl("_Oa", meta[[group_col]]), "O", NA))
  meta <- meta[!is.na(meta$AgeGroup), ]
  meta$AgeGroup <- factor(meta$AgeGroup, levels = c("Y", "O"))
  
  message("Aggregated AgeGroups: ", paste(levels(meta$AgeGroup), collapse = ", "))
  print(table(meta$AgeGroup))
  # --- HVG selection ---
  hvgs <- VariableFeatures(seu)
  if (length(hvgs) == 0) {
    seu <- FindVariableFeatures(seu, assay = assay,
                                selection.method = selection.method,
                                nfeatures = n_hvg)
    hvgs <- VariableFeatures(seu)
  } else if (length(hvgs) > n_hvg) hvgs <- hvgs[1:n_hvg]
  # --- Expression matrix ---
  expr <- tryCatch(
    GetAssayData(seu, assay = assay, layer = slot),
    error = function(e) GetAssayData(seu, assay = assay, slot = slot)
  )
  expr <- expr[intersect(rownames(expr), hvgs), , drop = FALSE]
  expr <- expr[, rownames(meta), drop = FALSE]
  # --- Filter low-expression genes ---
  detect_pct <- Matrix::rowMeans(expr > 0)
  mean_expr  <- Matrix::rowMeans(expr)
  keep_genes <- names(which(detect_pct >= min_pct & mean_expr >= min_mean))
  expr <- expr[keep_genes, , drop = FALSE]
  message("Keeping ", length(keep_genes), " genes after filtering.")
  # --- Mean & variance per group ---
  groups <- levels(meta$AgeGroup)
  mean_tbl <- map_dfc(groups, function(g) {
    cells <- rownames(meta)[meta$AgeGroup == g]
    apply(expr[, cells, drop = FALSE], 1, mean, na.rm = TRUE)
  })
  var_tbl <- map_dfc(groups, function(g) {
    cells <- rownames(meta)[meta$AgeGroup == g]
    apply(expr[, cells, drop = FALSE], 1, var, na.rm = TRUE)
  })
  colnames(mean_tbl) <- paste0("mean_", groups)
  colnames(var_tbl)  <- paste0("var_",  groups)
  gene_stats <- bind_cols(gene = rownames(expr), mean_tbl, var_tbl)
  # --- LOESS mean–variance correction (bias-stabilized) ---
  long_stats <- gene_stats %>%
    pivot_longer(cols = starts_with("mean_"), names_to = "mean_group", values_to = "mean") %>%
    mutate(group = sub("^mean_", "", mean_group)) %>%
    left_join(
      gene_stats %>%
        pivot_longer(cols = starts_with("var_"), names_to = "var_group", values_to = "var") %>%
        mutate(group = sub("^var_", "", var_group)),
      by = c("gene", "group")
    ) %>%
    mutate(
      log_mean = log10(pmax(mean, 1e-3)),   # floor to stabilize near-zero means
      log_var  = log10(pmax(var,  1e-4))    # floor to stabilize near-zero vars
    ) %>%
    filter(is.finite(log_mean), is.finite(log_var))
  # --- Fit LOESS on 5–95% quantiles only ---
  fit_loess <- loess(log_var ~ log_mean,
                     data = long_stats %>%
                       filter(between(log_mean,
                                      quantile(log_mean, 0.05, na.rm = TRUE),
                                      quantile(log_mean, 0.95, na.rm = TRUE))),
                     span = 0.75)
  # --- LOESS prediction ---
  predicted_logVar <- predict(fit_loess, newdata = data.frame(log_mean = long_stats$log_mean))
  # Replace NA predictions (outside fit range) with boundary values
  if (anyNA(predicted_logVar)) {
    range_fit <- range(fit_loess$x, na.rm = TRUE)
    low_val  <- predict(fit_loess, newdata = data.frame(log_mean = range_fit[1]))
    high_val <- predict(fit_loess, newdata = data.frame(log_mean = range_fit[2]))
    predicted_logVar[long_stats$log_mean < range_fit[1]] <- low_val
    predicted_logVar[long_stats$log_mean > range_fit[2]] <- high_val
  }
  # Compute residuals safely
  long_stats$residual_logVar <- long_stats$log_var - predicted_logVar
  # --- Summarize adjusted variance ---
  resid_summary <- long_stats %>%
    group_by(gene, group) %>%
    summarise(mean_adjusted_var = mean(10 ^ residual_logVar, na.rm = TRUE), .groups = "drop") %>%
    complete(gene, group = c("Y", "O")) %>%
    pivot_wider(names_from = group, values_from = mean_adjusted_var,
                names_glue = "mean_adjusted_var_{group}") %>%
    filter(!(is.na(mean_adjusted_var_Y) & is.na(mean_adjusted_var_O)))
  # --- Add log2FC (O/Y) ---
  resid_summary <- resid_summary %>%
    mutate(log2FC_mean_adjusted_variance_YO =
             log2((mean_adjusted_var_O + 1e-8) / (mean_adjusted_var_Y + 1e-8)))
  # --- Brown–Forsythe test (Y vs O) ---
  message("Running Brown–Forsythe variance test (Y vs O)...")
  res_df <- map_dfr(rownames(expr), function(g) {
    df_test <- data.frame(val = expr[g, ], grp = meta$AgeGroup)
    df_test <- df_test[!is.na(df_test$val), ]
    if (length(unique(df_test$grp)) < 2) return(tibble(gene = g, p_bf = NA))
    p_bf <- tryCatch(car::leveneTest(val ~ grp, data = df_test, center = median)[1, "Pr(>F)"], error = function(e) NA)
    tibble(gene = g, p_brown_forsythe_All_YO = p_bf)
  })
  res_df$fdr_brown_forsythe_All_YO <- p.adjust(res_df$p_brown_forsythe_All_YO, "BH")
  # --- Merge all results ---
  merged <- gene_stats %>%
    left_join(resid_summary, by = "gene") %>%
    left_join(res_df, by = "gene") %>%
    mutate(log2_var_ratio_OY = log2((var_O + 1e-8) / (var_Y + 1e-8)))
  # --- Direction check ---
  dir_check <- mean(sign(merged$mean_adjusted_var_O - merged$mean_adjusted_var_Y) ==
                      sign(merged$log2FC_mean_adjusted_variance_YO), na.rm = TRUE)
  message(sprintf("🔍 Direction consistency: %.1f%% genes match expected O/Y direction", dir_check * 100))
  message("Finished Exp2 Y–O variance computation (bias-stabilized, NA-safe).")
  return(merged)
}

# Plot function for Exp2 

plot_variance_cv_summaries_exp2 <- function(
    df,
    fdr_col = "fdr_brown_forsythe_All_YO",
    fdr_cutoff = 0.05,
    fc_cutoff = 0.25,
    color_palette = c("#67A5CC", "#D780AA"),
    prefix = "Exp2 (Y–O Comparison)"
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggrepel)
    library(ggsignif)
  })
  age_levels <- c("Y", "O")
  n_genes <- n_distinct(df$gene)
  # --- Boxplot: mean-adjusted variance ---
  long_adjvar <- df %>%
    dplyr::select(gene, starts_with("mean_adjusted_var_")) %>%
    pivot_longer(-gene, names_to = "Group", values_to = "MeanAdjVar") %>%
    mutate(Group = sub("mean_adjusted_var_", "", Group),
           Group = factor(Group, levels = age_levels))
  p_adjvar <- ggplot(long_adjvar, aes(Group, MeanAdjVar, color = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
    geom_signif(comparisons = list(c("Y", "O")), test = "wilcox.test", textsize = 3.5) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 15) +
    labs(
      title = sprintf("%s: Mean-Adjusted Variance (Y vs O)", prefix),
      subtitle = sprintf("n = %d genes | Wilcoxon test", n_genes),
      y = "Mean-Adjusted Variance"
    )
  # --- Volcano Plot ---
  df_volcano <- df %>%
    mutate(
      log2FC = log2FC_mean_adjusted_variance_YO,
      neg_log10_fdr = -log10(pmax(.data[[fdr_col]], 1e-300)),
      significance = case_when(
        .data[[fdr_col]] < fdr_cutoff & log2FC >  fc_cutoff ~ "↑O",
        .data[[fdr_col]] < fdr_cutoff & log2FC < -fc_cutoff ~ "↑Y",
        TRUE ~ "ns"))
  top_genes <- bind_rows(
    df_volcano %>%
      filter(significance == "↑O") %>%
      arrange(.data[[fdr_col]], desc(log2FC)) %>%
      slice_head(n = 10),
    df_volcano %>%
      filter(significance == "↑Y") %>%
      arrange(.data[[fdr_col]], log2FC) %>%
      slice_head(n = 10)) %>%
    mutate(gene_label = paste0("italic('", gene, "')"))
  p_volcano <- ggplot(df_volcano,
                      aes(x = log2FC,
                          y = neg_log10_fdr,
                          color = significance)) +
    geom_point(alpha = 0.85, size = 3) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene_label),
      parse = TRUE,
      color = "black",
      size = 5,
      box.padding = 0.5,
      point.padding = 0.5,
      force = 2
    ) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed") +
    coord_cartesian(xlim = c(-2.3, 2.3)) +
    scale_color_manual(
      values = c(
        "↑O" = "#D85B59",
        "↑Y" = "#5271AE",
        "ns" = "grey80"
      ),
      guide = "none"     # removes legend
    ) +
    theme_classic(base_size = 20) +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 22, face = "bold", color = "black"),
      axis.title.y = element_text(size = 22, face = "bold", color = "black"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    labs(
      x = expression(bold(log[2]("Mean-Adj Var (O / Y)"))),
      y = expression(bold(-log[10]("FDR"))))
  # --- Density Plot: Mean-Adjusted Variance (Y vs O) ---
  p_density <- ggplot(long_adjvar, aes(x = MeanAdjVar, fill = Group)) +
    geom_density(alpha = 0.4, adjust = 1.2, color = NA) +
    scale_fill_manual(values = color_palette, name = "Group", labels = c("Y", "O")) +
    theme_classic(base_size = 16) +
    labs(
      title = sprintf("%s – Mean-Adjusted Variance Distribution (Y vs O)", prefix),
      subtitle = sprintf("n = %d genes", n_genes),
      x = "Mean-Adjusted Variance",
      y = "Density"
    )
  # --- Return plots ---
  list(
    mean_adjusted_variance_box = p_adjvar,
    volcano = p_volcano,
    density = p_density
  )
}

### Compute DVG analysis for Exp2 ###

# Subset LT/ST-HSCs first
hsc_exp2 <- subset(seurat_vitro_exp2, celltype %in% c("LT-HSC", "ST-HSC"))
table(hsc_exp2$sampleName)
# filter out unrelated/ technical genes
hsc_exp2 <- filter_unrelated_genes(hsc_exp2) # 13617 genes 
# Compute aggregated variance–CV table (Y vs O)
exp2_HSC_results <- compute_variance_and_cv_YO_pairs_exp2(hsc_exp2) # 782 genes
View(exp2_HSC_results)

# Pick aging-up HVGs (positive FC, O > Y)
aging_up_genes <- exp2_HSC_results %>%
  filter(log2FC_mean_adjusted_variance_YO > 0, fdr_brown_forsythe_All_YO < 0.05)
nrow(aging_up_genes) # 256
# Select Top200 / Top100 aging-up genes
hsc_results_exp2_up_top200 <- aging_up_genes %>%
  arrange(fdr_brown_forsythe_All_YO) %>%
  slice_head(n = 200)
hsc_results_exp2_up_top100 <- aging_up_genes %>%
  arrange(fdr_brown_forsythe_All_YO) %>%
  slice_head(n = 100)

# Plot 
plots_exp2_all_cv <- plot_variance_cv_summaries_exp2(
  df = exp2_HSC_results,
  fc_cutoff = 0.1,
  prefix = "Exp2 LT/ST-HSC – All Genes (Aggregated Y vs O)"
)
plots_exp2_all_cv$density
plots_exp2_all_cv$volcano
# Top 100 DVGs plot 
plots_exp2_top100_cv <- plot_variance_cv_summaries_exp2(
  df = hsc_results_exp2_up_top200,
  prefix = "Exp2_HSC – Top 100 Aging-Up DVGs"
)
plots_exp2_top100_cv$mean_adjusted_variance_box
plots_exp2_top100_cv$volcano
plots_exp2_top100_cv <- plot_variance_cv_summaries_exp2(
  df = hsc_results_exp2_up_top100,
  prefix = "Exp2_HSC – Top 100 Aging-Up HVGs"
)
plots_exp2_top100_cv$mean_adjusted_variance_box
plots_exp2_top100_cv$volcano
plots_exp2_top100_cv$density

#### Save DVG tables from the 3 Exp #### 
# Cross Vitro DVGs 
Day0_hsc_dvgs_results$fdr_brown_forsythe
Day0_hsc_dvgs_results$log2FC_mean_adjusted_variance
# Unmanipulated DVGs
unmanipulated_results$fdr_brown_forsythe
unmanipulated_results$log2FC_mean_adjusted_var_old_young
# Exp2 DVGs 
exp2_HSC_results$fdr_brown_forsythe_All_YO
exp2_HSC_results$log2FC_mean_adjusted_variance_YO

############################################################
# 1 Standardize tables
############################################################

cross_tbl <- Day0_hsc_dvgs_results %>% dplyr::select(gene, fdr_brown_forsythe, log2FC_mean_adjusted_variance) %>%
  dplyr::rename(log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_variance)
unmanip_tbl <- unmanipulated_results %>% dplyr::select(gene, fdr_brown_forsythe,
                                                       log2FC_mean_adjusted_var_old_young) %>% dplyr::rename(log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_var_old_young)
exp2_tbl <- exp2_HSC_results %>% dplyr::select(gene,fdr_brown_forsythe_All_YO, log2FC_mean_adjusted_variance_YO) %>%
  dplyr::rename(fdr_brown_forsythe = fdr_brown_forsythe_All_YO, log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_variance_YO)

############################################################
# 2 Significant aging-up genes
############################################################

cross_sig <- cross_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)
unmanip_sig <- unmanip_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)
exp2_sig <- exp2_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)

############################################################
# 3 Gene sets
############################################################

cross_genes <- cross_sig$gene # 117
unmanip_genes <- unmanip_sig$gene # 114
exp2_genes <- exp2_sig$gene # 92

############################################################
# 4 Overlap sets
############################################################

cross_unmanip <- intersect(cross_genes, unmanip_genes)
cross_exp2 <- intersect(cross_genes, exp2_genes)
unmanip_exp2 <- intersect(unmanip_genes, exp2_genes)
consensus_genes <- Reduce(intersect,list(cross_genes,unmanip_genes,exp2_genes))

############################################################
# 5 Summary table with gene lists
############################################################

summary_table <- tibble(Comparison = c("Cross","Unmanipulated","Exp2","Cross ∩ Unmanipulated",
                                       "Cross ∩ Exp2","Unmanipulated ∩ Exp2","Consensus (All Three)"),
                        Gene_Count = c(length(cross_genes),length(unmanip_genes),length(exp2_genes),
                                       length(cross_unmanip),length(cross_exp2),length(unmanip_exp2),length(consensus_genes)),
                        Gene_List = c(paste(cross_genes, collapse = ", "),paste(unmanip_genes, collapse = ", "),
                                      paste(exp2_genes, collapse = ", "),paste(cross_unmanip, collapse = ", "),
                                      paste(cross_exp2, collapse = ", "),paste(unmanip_exp2, collapse = ", "),
                                      paste(consensus_genes, collapse = ", ")))

############################################################
# 6 Meta table (genes shared by ≥2 experiments)
############################################################

gene_list <- list(Cross = cross_genes,Unmanipulated = unmanip_genes,Exp2 = exp2_genes)
gene_counts <- table(unlist(gene_list))
shared_genes <- names(gene_counts[gene_counts >= 2])
meta_table <- bind_rows(cross_sig %>% filter(gene %in% shared_genes) %>%
                          mutate(Experiment = "Cross"), unmanip_sig %>% filter(gene %in% shared_genes) %>%
                          mutate(Experiment = "Unmanipulated"),
                        exp2_sig %>% filter(gene %in% shared_genes) %>% mutate(Experiment = "Exp2")) %>% arrange(gene)

############################################################
# 7 Write Excel
############################################################

wb <- createWorkbook()
addWorksheet(wb, "Cross_All")
writeData(wb, "Cross_All", cross_tbl)
addWorksheet(wb, "Unmanipulated_All")
writeData(wb, "Unmanipulated_All", unmanip_tbl)
addWorksheet(wb, "Exp2_All")
writeData(wb, "Exp2_All", exp2_tbl)
addWorksheet(wb, "Cross_Significant")
writeData(wb, "Cross_Significant", cross_sig)
addWorksheet(wb, "Unmanipulated_Significant")
writeData(wb, "Unmanipulated_Significant", unmanip_sig)
addWorksheet(wb, "Exp2_Significant")
writeData(wb, "Exp2_Significant", exp2_sig)
addWorksheet(wb, "Summary_Shared_Genes")
writeData(wb, "Summary_Shared_Genes", summary_table)
addWorksheet(wb, "Meta_Shared_Genes")
writeData(wb, "Meta_Shared_Genes", meta_table)
saveWorkbook(wb,"/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure3/Tables/Fig3_DVG_Master_Table.xlsx",overwrite = TRUE)

##################################################################


##################################################################

# Visualized per-gene experssion and mean-adj variance 

# 80 gene list 
genes_of_interest <- c(
  "Aldh1a1","Nupr1","Slamf1","Cdkn2c","Cd48","Slc14a1","Otos","Krt18","Apoe","Hmgb2",
  "Cavin3","Nkx2-3","Prkca","Ran","Tubb5","Mki67","Birc5","Tespa1","Top2a","Gng11",
  "Cdkn1a","Txnip","Ncl","Zmat3","Serpinb1a","Xbp1","Ccnd1","Tubb4b","Slc22a3","Pclaf",
  "Ybx1","Cks1b","Cdca8","Mid1","Cenpv","Lgals1","Lgals9","Hacd4","H1fx","Nrgn",
  "Smc4","Nfia","Ptms","Tuba1b","Rora","Smc2","Myl10","Cdk1","Fbxo5","Ccna2",
  "Pbk","Plek","Alcam","Cks2","Anln","Esco2","Cd53","Ube2s","Knl1","Hbb-bt",
  "Incenp","Frat2","Kif20b","Procr","Neat1","Dlgap5","Ccng1","Hspe1","Tnip3","Adgrg3",
  "Ckap2l","Slfn9","Pbx1","Sgo1","Racgap1","Brip1","Snhg3","Cycs","Retreg1","Clec4d")

# 33 gene list 
genes_of_interest <- c(
  "Nupr1","Slamf1","Cd48","Otos","Krt18","Apoe","Hmgb2","Cavin3",
  "Nkx2-3","Txnip","Ran","Tubb5","Mki67","Birc5","Tespa1","Top2a",
  "Gng11","Cdkn1a","Hacd4","Smc4","Serpinb1a","Xbp1","Ccnd1",
  "Tubb4b","Slc22a3","Ybx1","Cdca8","H1fx","Ptms","Smc2","Cks2",
  "Esco2","Cd53")

plot_density_gene <- function(
    seurat_obj,
    gene,
    group_col,
    young_labels,
    old_labels,
    title_prefix,
    dvgs_df = NULL,
    mean_adj_var_Y_col = "mean_adjusted_var_Y",
    mean_adj_var_O_col = "mean_adjusted_var_O"
){
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
  })
  # Skip gene if not present
  if(!(gene %in% rownames(seurat_obj))){
    message("Skipping ", gene, " (not found)")
    return(NULL)
  }
  
  expr_df <- FetchData(seurat_obj, vars = c(gene, group_col)) %>%
    mutate(GroupRaw = .data[[group_col]]) %>%
    mutate(
      Group = case_when(
        GroupRaw %in% young_labels ~ "Y",
        GroupRaw %in% old_labels   ~ "O",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Group)) %>%
    mutate(Expression = .data[[gene]]) %>%
    filter(Expression > 0) %>%
    mutate(
      Group = factor(Group, levels = c("Y", "O"), ordered = TRUE) 
    )
  color_palette <- c(
    "Y" = "#89AEEB",
    "O" = "#F5A36C"
  )
  
  p_density <- ggplot(expr_df, aes(x = Expression, color = Group, fill = Group)) +
    geom_density(alpha = 0.25, linewidth = 1.1) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 18) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    ) +
    labs(
      title = paste0(title_prefix, " — ", gene),
      x = "Normalized Expression",
      y = "Density"
    )
  
  # variance bar
  if(!is.null(dvgs_df) && gene %in% dvgs_df$gene){
    var_vals <- dvgs_df %>%
      dplyr::filter(gene == !!gene)
    
    mean_var_df <- tibble(
      Group = factor(c("Y", "O"),
                     levels = c("Y", "O")), 
      MeanAdjVar = c(
        var_vals[[mean_adj_var_Y_col]][1],
        var_vals[[mean_adj_var_O_col]][1]
      )
    )
    
    p_bar <- ggplot(mean_var_df,
                    aes(x = Group, y = MeanAdjVar, fill = Group)) +
      geom_bar(stat = "identity", width = 0.6, color = "black") +
      scale_fill_manual(values = color_palette) +
      theme_classic(base_size = 18) +
      theme(
        axis.title = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(
        title = "Mean-Adj Var",
        x = NULL,
        y = "Mean-Adj Var"
      ) +
      scale_y_continuous(expand = c(0,0), limits = c(0,NA))
    p_final <- p_density + p_bar + plot_layout(widths = c(2.2,1))
  } else {
    p_final <- p_density
  }
  return(p_final)
}

make_4panel_for_gene <- function(gene){
  p_unman <- plot_density_gene(
    seurat_obj   = seurat_unmanipulated_HSC,
    gene         = gene,
    group_col    = "AGE",
    young_labels = c("young","Y"),
    old_labels   = c("old","O"),
    title_prefix = "Unmanipulated HSCs",
    dvgs_df = unmanipulated_results,
    mean_adj_var_Y_col = 'mean_adjusted_var_young', 
    mean_adj_var_O_col = 'mean_adjusted_var_old'
  )
  p_vitro <- plot_density_gene(
    seurat_obj   = seurat_in_vitro_HSC,
    gene         = gene,
    group_col    = "sampleName",
    young_labels = c("Y_vitro"),
    old_labels   = c("O_vitro"),
    title_prefix = "Day0 In Vitro HSCs",
    dvgs_df = Day0_hsc_dvgs_results
  )
  p_exp1 <- plot_density_gene(
    seurat_obj   = hsc_exp1,
    gene         = gene,
    group_col    = "sampleName",
    young_labels = c("Y2","Y3","Y4"),
    old_labels   = c("O2","O3","O4"),
    title_prefix = "Exp1 HSCs",
    dvgs_df      = hsc_results_exp1
  )
  p_exp2 <- plot_density_gene(
    seurat_obj   = hsc_exp2,
    gene         = gene,
    group_col    = "sampleName",
    young_labels = c("1_Ya","2_Ya","3_Ya"),
    old_labels   = c("4_Oa","5_Oa","6_Oa"),
    title_prefix = "Exp2 HSCs",
    dvgs_df      = exp2_HSC_results
  )
  combined_panel <- (p_unman | p_vitro) /
    (p_exp1 | p_exp2)
  
  return(combined_panel)
}

panel_list <- lapply(genes_of_interest, make_4panel_for_gene)
names(panel_list) <- genes_of_interest

pdf("Figure3/Figures/DensityPanels_DVG_0.05FC_80_Genes.pdf", width = 15, height = 12)
for(g in genes_of_interest){
  if(!is.null(panel_list[[g]])){
    print(panel_list[[g]])
  }
}
dev.off()

#####################################################################################

# DEG for DVGs

# Cross
                    
table(seurat_in_vitro_HSC$sampleName)
deg_vitro <- run_clonewise_DEG_suite(
  seurat_obj = seurat_in_vitro_HSC,
  clone_set1_ids = "O_vitro",   # Group1 = Old
  clone_set2_ids = "Y_vitro",   # Group2 = Young
  clone_col = "sampleName",
  group1_label = "Old_vitro",
  group2_label = "Young_vitro",
  fc_cutoff  = 0.25,
  fdr_cutoff = 0.05,
  excel_name  = "Cross__Old_vs_Young_vitro_DEG.xlsx",
  volcano_name = "Cross__Old_vs_Young_vitro_DEG_Volcano.png",
  output_dir = "Cross_Vitro_DEG_Results"
)

View(deg_vitro$results$DEG2_TopHVGs_MAST)
deg_vitro$volcano_plot

# Exp2 

table(hsc_exp2$AgeGroup)
deg_exp2_vitro <- run_clonewise_DEG_suite(
  seurat_obj = hsc_exp2,
  clone_set1_ids = "Old",
  clone_set2_ids = "Young",
  clone_col = "AgeGroup",
  group1_label = "Old_vitro",
  group2_label = "Young_vitro",
  fc_cutoff  = 0.25,
  fdr_cutoff = 0.05,
  excel_name  = "Exp2_Old_vs_Young_vitro_DEG.xlsx",
  volcano_name = "Exp2_Old_vs_Young_vitro_DEG_Volcano.png",
  output_dir = "Exp2_Vitro_DEG_Results"
)

View(deg_exp2_vitro$results$DEG2_TopHVGs_MAST)
deg_exp2_vitro$volcano_plot

####################################################################################
# End of this script 
