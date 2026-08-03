###################################################################################################


### Compute Differentially variance Genes (DVGs) in Cross(Day0), wildetype, exp1 and exp2 data


###################################################################################################

# Maintainer: Chris Chen
# Last updated: 03/03/2026

##########################

setwd('...')

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

# Read in the 4 dataset

seurat_in_vitro <- readRDS(".../Data/CrossAge(exp2)_vitro.RDS")

seurat_unmanipulated <- readRDS('.../Data/Unmanipulated_vitro.rds')

seurat_vitro_exp2 <- readRDS('.../Data/Exp2(exp1)_vitro.RDS')

# inspect data

table(seurat_in_vitro$sampleName) # O_vitro: 6337 ; Y_vitro: 10971
seurat_unmanipulated  # 18046 genes and 48140 cells
table(seurat_unmanipulated$AGE) # mid: 10877 ; old:11966 ; vold:12604 ; young: 12693
seurat_vitro_exp1 # 18389 genes and 28386 cells
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
  
  idx_g1 <- which(group_factor == g1)
  idx_g2 <- which(group_factor == g2)
  
  # ---------------- Mean-variance trend (fit BEFORE testing) ----------------
  mean_g1 <- Matrix::rowMeans(expr[, idx_g1, drop = FALSE])
  mean_g2 <- Matrix::rowMeans(expr[, idx_g2, drop = FALSE])
  var_g1_raw <- apply(expr[, idx_g1, drop = FALSE], 1, var)
  var_g2_raw <- apply(expr[, idx_g2, drop = FALSE], 1, var)
  
  trend_df <- data.frame(
    mean = c(mean_g1, mean_g2),
    var  = c(var_g1_raw, var_g2_raw)
  )
  trend_df <- trend_df[
    !is.na(trend_df$mean) & trend_df$mean > 0 &
      !is.na(trend_df$var)  & trend_df$var  > 0,
  ]
  
  fit_loess <- loess(log10(var) ~ log10(mean), data = trend_df, span = 0.75,
                     control = loess.control(surface = "direct"))
  
  expected_var <- function(mean_vec) {
    out <- rep(NA_real_, length(mean_vec))
    ok <- !is.na(mean_vec) & mean_vec > 0
    out[ok] <- 10 ^ predict(fit_loess, newdata = data.frame(mean = mean_vec[ok]))
    out
  }
  
  scale_g1 <- sqrt(expected_var(mean_g1))
  scale_g2 <- sqrt(expected_var(mean_g2))
  testable <- !is.na(scale_g1) & !is.na(scale_g2) & scale_g1 > 0 & scale_g2 > 0
  
  adj_expr <- expr
  adj_expr[, idx_g1] <- expr[, idx_g1, drop = FALSE] / scale_g1
  adj_expr[, idx_g2] <- expr[, idx_g2, drop = FALSE] / scale_g2
  
  message(sum(!testable), " gene(s) excluded: mean-variance trend not evaluable.")
  
  # ---------------- Variance tests (run on trend-adjusted expression) ----------------
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
  test_genes <- rownames(expr)[testable]
  chunks <- split(test_genes, ceiling(seq_along(test_genes) / chunk_size))
  results_list <- list()
  with_progress({
    p <- progressor(steps = length(chunks))
    for (i in seq_along(chunks)) {
      p(sprintf("Chunk %d / %d", i, length(chunks)))
      genes <- chunks[[i]]
      sub_adj <- adj_expr[genes, , drop = FALSE]
      sub_raw <- expr[genes, , drop = FALSE]
      mat <- apply(sub_adj, 1, do_tests, group_factor = group_factor)
      bf_p   <- sapply(mat, `[[`, "bf")
      lev_p  <- sapply(mat, `[[`, "lev")
      bart_p <- sapply(mat, `[[`, "bart")
      var_tbl <- apply(sub_raw, 1, function(x) {
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
  # ---------------- Mean-adjusted variance (from the same adjusted values used for testing) ----------------
  message("Computing mean-adjusted variance and SD ...")
  mean_adjusted_var_g1 <- apply(adj_expr[, idx_g1, drop = FALSE], 1, var, na.rm = TRUE)
  mean_adjusted_var_g2 <- apply(adj_expr[, idx_g2, drop = FALSE], 1, var, na.rm = TRUE)
  adj_tbl <- data.frame(gene = rownames(expr))
  adj_tbl[[paste0("mean_adjusted_var_", g1)]] <- mean_adjusted_var_g1
  adj_tbl[[paste0("mean_adjusted_var_", g2)]] <- mean_adjusted_var_g2
  adj_tbl[[paste0("mean_adjusted_sd_", g1)]]  <- sqrt(mean_adjusted_var_g1)
  adj_tbl[[paste0("mean_adjusted_sd_", g2)]]  <- sqrt(mean_adjusted_var_g2)
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
write.xlsx(Day0_hsc_dvgs_results,file = ".../Cross_vitro_HSC_DVG_results.xlsx",rowNames = FALSE)

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
  # Density Plot (Significant Genes Only)
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
  # Summary Counts (log2FC already fixed upstream)
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
    arrange(.data[[fdr_col]], desc(.data[[fc_col]])) %>%   # significance first
    slice_head(n = 8)
  
  # Top genes increased in g1
  top_up_g1 <- df_vol_mvar %>%
    filter(significance == paste0("↑", g1)) %>%
    arrange(.data[[fdr_col]], .data[[fc_col]]) %>%         # significance first
    slice_head(n = 8)
  
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
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
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

# Figure 3 b inter-cell variance density plot
day0_hsc_plots$density_mean_adj_sig <- day0_hsc_plots$density_mean_adj_sig +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  scale_x_continuous(breaks = c(1.0, 1.5, 2.0)) +
  coord_cartesian(xlim = c(NA, 2.2)) +
  labs(title = "Inter-cell\nvariance") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22, lineheight = 1.1),
    plot.title.position = "plot"
  )

day0_hsc_plots$density_mean_adj_sig$layers[[1]]$aes_params$linewidth <- 1.6

day0_hsc_plots$density_mean_adj_sig

# Figure 3 c Differentially variable genes (DVG)
day0_hsc_plots$volcano_mean_adj <- day0_hsc_plots$volcano_mean_adj +
  labs(
    title = "Differentially variable\ngenes (DVG)",
    y = expression(-log[10](FDR)) 
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22, lineheight = 1.1),
    plot.title.position = "plot",
    axis.title.x = element_text(face = "plain", size = 22),
    axis.title.y = element_text(face = "plain", size = 22)
  )

day0_hsc_plots$volcano_mean_adj


###########################################################


##########################################

## Compute DVGs (Y vs O) in unmanipulated data

##########################################

# function: compute mean adjusted variance for unmanipulated data (multiple agegroups)

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
  age_levels <- levels(age)
  message("Age levels detected: ", paste(age_levels, collapse = ", "))
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
  
  idx_by_group <- lapply(age_levels, function(g) which(age == g))
  names(idx_by_group) <- age_levels
  
  # ---------------- Mean-variance trend (fit BEFORE testing) ----------------
  mean_by_group <- lapply(idx_by_group, function(idx) Matrix::rowMeans(expr[, idx, drop = FALSE]))
  var_by_group  <- lapply(idx_by_group, function(idx) apply(expr[, idx, drop = FALSE], 1, var))
  
  trend_df <- data.frame(
    mean = unlist(mean_by_group),
    var  = unlist(var_by_group)
  )
  trend_df <- trend_df[
    !is.na(trend_df$mean) & trend_df$mean > 0 &
      !is.na(trend_df$var)  & trend_df$var  > 0,
  ]
  
  fit_loess <- loess(log10(var) ~ log10(mean), data = trend_df, span = 0.75,
                     control = loess.control(surface = "direct"))
  
  expected_var <- function(mean_vec) {
    out <- rep(NA_real_, length(mean_vec))
    ok <- !is.na(mean_vec) & mean_vec > 0
    out[ok] <- 10 ^ predict(fit_loess, newdata = data.frame(mean = mean_vec[ok]))
    out
  }
  
  scale_by_group <- lapply(mean_by_group, function(m) sqrt(expected_var(m)))
  testable <- Reduce(`&`, lapply(scale_by_group, function(s) !is.na(s) & s > 0))
  
  adj_expr <- expr
  for (g in age_levels) {
    adj_expr[, idx_by_group[[g]]] <- expr[, idx_by_group[[g]], drop = FALSE] / scale_by_group[[g]]
  }
  
  message(sum(!testable), " gene(s) excluded: mean-variance trend not evaluable.")
  
  # --- Variance tests (Brown–Forsythe, Levene, Bartlett), run on trend-adjusted expression ---
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
  test_genes <- rownames(expr)[testable]
  chunks <- split(test_genes, ceiling(seq_along(test_genes) / chunk_size))
  results_list <- list()
  with_progress({
    p <- progressor(steps = length(chunks))
    for (i in seq_along(chunks)) {
      p(sprintf("Chunk %d / %d", i, length(chunks)))
      genes <- chunks[[i]]
      sub_adj <- adj_expr[genes, , drop = FALSE]
      mat <- apply(sub_adj, 1, do_tests, age = age)
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
  # --- Mean-adjusted variance (from the same adjusted expression used for testing) ---
  mean_adjusted_var <- lapply(age_levels, function(g) {
    apply(adj_expr[, idx_by_group[[g]], drop = FALSE], 1, var, na.rm = TRUE)
  })
  names(mean_adjusted_var) <- age_levels
  resid_summary <- data.frame(gene = rownames(expr))
  for (g in age_levels) {
    resid_summary[[paste0("mean_adjusted_var_", g)]] <- mean_adjusted_var[[g]]
    resid_summary[[paste0("mean_adjusted_sd_", g)]]  <- sqrt(mean_adjusted_var[[g]])
  }
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
        x = expression(log[2]("Mean-Adj Var (O/Y)")),
        y = expression(-log[10](FDR))
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

# figure3 h DVG unmanipulated HSC 
plots_unmanipulated$volcano_old <- plots_unmanipulated$volcano_old +
  labs(title = "DVG (unmanipulated HSC)") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22, lineheight = 1.1),
    plot.title.position = "plot"
  )

plots_unmanipulated$volcano_old


##########################

# Cross aging up DVG in EXP2 

gene_vec <- Day0_hsc_dvgs_results %>% dplyr::filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance > 0) %>% dplyr::pull(gene)
df_exp2 <- exp2_HSC_results %>% dplyr::filter(gene %in% gene_vec) %>% dplyr::select(gene, mean_adjusted_var_Y, mean_adjusted_var_O) %>%
  tidyr::pivot_longer(cols = c(mean_adjusted_var_Y, mean_adjusted_var_O), names_to = "group", values_to = "mean_adj_var") %>% dplyr::mutate(
    group = sub("mean_adjusted_var_", "", group))
df_exp2$group <- factor(df_exp2$group, levels = c("Y", "O"))
df_exp2_wide <- df_exp2 %>% tidyr::pivot_wider(names_from = group, values_from = mean_adj_var)
p_val <- wilcox.test(df_exp2_wide$Y, df_exp2_wide$O, paired = TRUE)$p.value # 3.6e-6
p_val
if (p_val < 0.001) {
  sig_label <- "***"
} else if (p_val < 0.01) {
  sig_label <- "**"
} else if (p_val < 0.05) {
  sig_label <- "*"
} else {
  sig_label <- "ns"
}
y_bracket_prop <- 0.9
y_star_prop    <- 0.95
y_range_obs <- c(0, 4)
y_bracket_obs <- y_range_obs[1] + diff(y_range_obs) * y_bracket_prop
y_star_obs    <- y_range_obs[1] + diff(y_range_obs) * y_star_prop
p_obs <- base_plot(df_exp2) +
  ggtitle("Old vs \nYoung") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, lineheight = 1.1),
    axis.title.y = element_text(face = "plain")
  ) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  annotate("segment", x = 1, xend = 2,
           y = 3.7, yend = 3.7,
           linewidth = 1.2) +
  annotate("text",
           x = 1.5,
           y = 3.85,
           label = sig_label,
           size = 7)
p_obs
## Randomized delta plot
obs_delta <- exp2_HSC_results %>%
  dplyr::filter(gene %in% gene_vec) %>%
  dplyr::transmute(
    gene,
    delta = mean_adjusted_var_O - mean_adjusted_var_Y,
    condition = "DVG")
set.seed(42)
gene_vec_exp2 <- intersect(gene_vec, exp2_HSC_results$gene)
n_genes <- length(gene_vec_exp2)
n_perm <- 200
perm_all <- purrr::map_dfr(1:n_perm, function(i) {
  rand_genes <- sample(all_genes, n_genes)
  exp2_HSC_results %>%
    dplyr::filter(gene %in% rand_genes) %>%
    dplyr::transmute(
      iter = i,
      gene,
      delta = mean_adjusted_var_O - mean_adjusted_var_Y
    )
})
perm_delta_avg <- perm_all %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    delta = mean(delta, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::slice_sample(n = n_genes) %>%
  dplyr::mutate(condition = "Ran")
delta_plot_df <- dplyr::bind_rows(obs_delta, perm_delta_avg)
delta_plot_df$condition <- factor(
  delta_plot_df$condition,
  levels = c("DVG", "Ran"))
p_rand <- ggplot(delta_plot_df,
                 aes(x = condition, y = delta, color = condition)) +
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
  scale_color_manual(values = c("DVG" = "#51C2BB", "Ran" = "#F57075")) +
  scale_fill_manual(values = c("DVG" = "#51C2BB", "Ran" = "#F57075")) +
  theme_classic(base_size = 20) +
  ggtitle("Observed vs \nRandomized") +
  labs(
    x = NULL,
    y = "Δ Mean-Adj Var (O − Y)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, lineheight = 1.1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "plain"),
    legend.position = "none"
  )
p_rand
levels(delta_plot_df$condition)
p_val_rand <- wilcox.test(delta ~ condition, data = delta_plot_df,alternative = "greater")$p.value
p_val_rand # 3.2e-8
if (p_val_rand < 0.001) {
  sig_label_rand <- "***"
} else if (p_val_rand < 0.01) {
  sig_label_rand <- "**"
} else if (p_val_rand < 0.05) {
  sig_label_rand <- "*"
} else {
  sig_label_rand <- "ns"
}
y_max_rand <- max(delta_plot_df$delta, na.rm = TRUE)
p_rand <- p_rand +
  # bracket
  annotate("segment", x = 1, xend = 2,
           y = y_max_rand * 1.15, yend = y_max_rand * 1.15,
           linewidth = 1.2) +
  annotate("text",
           x = 1.5,
           y = y_max_rand * 1.25,
           label = sig_label_rand,
           size = 7)
p_rand
### Combine 2 plots; ; figure 3 e
p_obs
p_rand
p_combined <- p_obs + p_rand + plot_layout(ncol = 2)
p_combined

##################################

# Cross DVG in Unmanipulated data
df_unmanip <- unmanipulated_results %>% dplyr::filter(gene %in% gene_vec) %>%
  dplyr::select(gene, mean_adjusted_var_young,mean_adjusted_var_old) %>%
  tidyr::pivot_longer(cols = c(mean_adjusted_var_young, mean_adjusted_var_old), names_to = "group", values_to = "mean_adj_var") %>%
  dplyr::mutate(group = sub("mean_adjusted_var_", "", group),
                group = dplyr::recode(group, young = "Young",old   = "Old")) %>%
  dplyr::mutate(group = factor(group, levels = c("Young", "Old")))
# paired test
df_unmanip_wide <- df_unmanip %>% tidyr::pivot_wider(names_from = group, values_from = mean_adj_var)
p_val <- wilcox.test(df_unmanip_wide$Young,df_unmanip_wide$Old,paired = TRUE)$p.value # 0.009
p_val
if (p_val < 0.001) {
  sig_label <- "***"
} else if (p_val < 0.01) {
  sig_label <- "**"
} else if (p_val < 0.05) {
  sig_label <- "*"
} else {
  sig_label <- "ns"
}
color_palette <- c("Young" = "cornflowerblue","Old"   = "#E58B1C")
# observed plot
p_obs_unmanip <- base_plot(df_unmanip) + ggtitle("Old vs \nYoung") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, lineheight = 1.1),
    axis.title.y = element_text(face = "plain")
  ) +
  scale_y_continuous(limits = c(0.5, 2.9), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Young" = "Y", "Old" = "O")) +
  annotate("segment", x = 1, xend = 2, y = 2.7, yend = 2.7,linewidth = 1.2) +
  annotate("text", x = 1.5, y = 2.8, label = sig_label,size = 7)
p_obs_unmanip
### Randomized
## Randomized delta plot (UNMANIPULATED)
obs_delta <- unmanipulated_results %>%
  dplyr::filter(gene %in% gene_vec) %>%
  dplyr::transmute(
    gene,
    delta = mean_adjusted_var_old - mean_adjusted_var_young,
    condition = "Obs"
  )
set.seed(42)
gene_vec_unmanip <- intersect(gene_vec, unmanipulated_results$gene)
n_genes <- length(gene_vec_unmanip)
n_perm <- 200
perm_all <- purrr::map_dfr(1:n_perm, function(i) {
  
  rand_genes <- sample(all_genes, n_genes)
  
  unmanipulated_results %>%
    dplyr::filter(gene %in% rand_genes) %>%
    dplyr::transmute(
      iter = i,
      gene,
      delta = mean_adjusted_var_old - mean_adjusted_var_young
    )
})
perm_delta_avg <- perm_all %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    delta = mean(delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::slice_sample(n = n_genes) %>%
  dplyr::mutate(condition = "Ran")
delta_plot_df <- dplyr::bind_rows(obs_delta, perm_delta_avg)
delta_plot_df$condition <- factor(
  delta_plot_df$condition,
  levels = c("Obs", "Ran")
)
nrow(obs_delta)
nrow(perm_delta_avg)
p_rand_unmanip <- ggplot(delta_plot_df,
                         aes(x = condition, y = delta, color = condition)) +
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
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Obs" = "DVG", "Ran" = "Ran")) +
  scale_color_manual(values = c("Obs" = "#51C2BB", "Ran" = "#F57075")) +
  scale_fill_manual(values = c("Obs" = "#51C2BB", "Ran" = "#F57075")) +
  theme_classic(base_size = 20) +
  ggtitle("Observed vs \nRandomized") +
  labs(
    x = NULL,
    y = "Δ Mean-Adj Var (O − Y)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, lineheight = 1.1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "plain"),
    legend.position = "none"
  )
p_val_rand <- wilcox.test(
  delta ~ condition,
  data = delta_plot_df,
  alternative = "greater"
)$p.value # 0.026
p_val_rand 
if (p_val_rand < 0.001) {
  sig_label_rand <- "***"
} else if (p_val_rand < 0.01) {
  sig_label_rand <- "**"
} else if (p_val_rand < 0.05) {
  sig_label_rand <- "*"
} else {
  sig_label_rand <- "ns"
}
y_max_rand <- max(delta_plot_df$delta, na.rm = TRUE)
p_rand_unmanip <- p_rand_unmanip +
  annotate("segment", x = 1, xend = 2,
           y = 0.88, yend = 0.88,
           linewidth = 1.2) +
  annotate("text",
           x = 1.5,
           y = 0.93,
           label = sig_label_rand,
           size = 7)
p_rand_unmanip
# combined the 2; figure 3 g 
p_combined_unmanip <- (p_obs_unmanip + p_rand_unmanip + plot_layout(ncol = 2))
p_combined_unmanip


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
  
  groups <- levels(meta$AgeGroup)
  idx_by_group <- lapply(groups, function(g) which(meta$AgeGroup == g))
  names(idx_by_group) <- groups
  
  # --- Mean & variance per group ---
  mean_by_group <- lapply(idx_by_group, function(idx) apply(expr[, idx, drop = FALSE], 1, mean, na.rm = TRUE))
  var_by_group  <- lapply(idx_by_group, function(idx) apply(expr[, idx, drop = FALSE], 1, var, na.rm = TRUE))
  mean_tbl <- as.data.frame(mean_by_group); colnames(mean_tbl) <- paste0("mean_", groups)
  var_tbl  <- as.data.frame(var_by_group);  colnames(var_tbl)  <- paste0("var_",  groups)
  gene_stats <- bind_cols(gene = rownames(expr), mean_tbl, var_tbl)
  
  # --- Mean-variance trend (fit BEFORE testing, bias-stabilized LOESS on 5-95% quantiles) ---
  trend_df <- data.frame(
    mean = unlist(mean_by_group),
    var  = unlist(var_by_group)
  ) %>%
    mutate(
      log_mean = log10(pmax(mean, 1e-3)),   # floor to stabilize near-zero means
      log_var  = log10(pmax(var,  1e-4))    # floor to stabilize near-zero vars
    ) %>%
    filter(is.finite(log_mean), is.finite(log_var))
  
  fit_loess <- loess(log_var ~ log_mean,
                     data = trend_df %>%
                       filter(between(log_mean,
                                      quantile(log_mean, 0.05, na.rm = TRUE),
                                      quantile(log_mean, 0.95, na.rm = TRUE))),
                     span = 0.75)
  
  expected_logvar <- function(mean_vec) {
    log_mean_vec <- log10(pmax(mean_vec, 1e-3))
    pred <- predict(fit_loess, newdata = data.frame(log_mean = log_mean_vec))
    if (anyNA(pred)) {
      range_fit <- range(fit_loess$x, na.rm = TRUE)
      low_val  <- predict(fit_loess, newdata = data.frame(log_mean = range_fit[1]))
      high_val <- predict(fit_loess, newdata = data.frame(log_mean = range_fit[2]))
      pred[log_mean_vec < range_fit[1]] <- low_val
      pred[log_mean_vec > range_fit[2]] <- high_val
    }
    pred
  }
  
  scale_by_group <- lapply(mean_by_group, function(m) sqrt(10 ^ expected_logvar(m)))
  
  adj_expr <- expr
  for (g in groups) {
    adj_expr[, idx_by_group[[g]]] <- expr[, idx_by_group[[g]], drop = FALSE] / scale_by_group[[g]]
  }
  
  # --- Mean-adjusted variance (from the same adjusted values used for testing) ---
  mean_adjusted_var <- lapply(groups, function(g) {
    apply(adj_expr[, idx_by_group[[g]], drop = FALSE], 1, var, na.rm = TRUE)
  })
  names(mean_adjusted_var) <- groups
  resid_summary <- data.frame(gene = rownames(expr))
  for (g in groups) {
    resid_summary[[paste0("mean_adjusted_var_", g)]] <- mean_adjusted_var[[g]]
  }
  # --- Add log2FC (O/Y) ---
  resid_summary <- resid_summary %>%
    mutate(log2FC_mean_adjusted_variance_YO =
             log2((mean_adjusted_var_O + 1e-8) / (mean_adjusted_var_Y + 1e-8)))
  # --- Brown–Forsythe test (Y vs O), run on trend-adjusted expression ---
  message("Running Brown–Forsythe variance test (Y vs O)...")
  res_df <- map_dfr(rownames(expr), function(g) {
    df_test <- data.frame(val = adj_expr[g, ], grp = meta$AgeGroup)
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
nrow(aging_up_genes) # 183
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
#  Standardize tables
############################################################

cross_tbl <- Day0_hsc_dvgs_results %>% dplyr::select(gene, fdr_brown_forsythe, log2FC_mean_adjusted_variance) %>%
  dplyr::rename(log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_variance)
unmanip_tbl <- unmanipulated_results %>% dplyr::select(gene, fdr_brown_forsythe,
                                                       log2FC_mean_adjusted_var_old_young) %>% dplyr::rename(log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_var_old_young)
exp2_tbl <- exp2_HSC_results %>% dplyr::select(gene,fdr_brown_forsythe_All_YO, log2FC_mean_adjusted_variance_YO) %>%
  dplyr::rename(fdr_brown_forsythe = fdr_brown_forsythe_All_YO, log2FC_mean_adjusted_variance_OY = log2FC_mean_adjusted_variance_YO)

############################################################
# Significant aging-up genes
############################################################

cross_sig <- cross_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)
unmanip_sig <- unmanip_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)
exp2_sig <- exp2_tbl %>% filter(fdr_brown_forsythe < 0.05, log2FC_mean_adjusted_variance_OY >= 0.1)

############################################################
# Gene sets
############################################################

cross_genes <- cross_sig$gene # 117
unmanip_genes <- unmanip_sig$gene # 114
exp2_genes <- exp2_sig$gene # 92

############################################################
# Overlap sets
############################################################

cross_unmanip <- intersect(cross_genes, unmanip_genes)
cross_exp2 <- intersect(cross_genes, exp2_genes)
unmanip_exp2 <- intersect(unmanip_genes, exp2_genes)
consensus_genes <- Reduce(intersect,list(cross_genes,unmanip_genes,exp2_genes))

############################################################
# Summary table with gene lists
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
# Meta table (genes shared by ≥2 experiments)
############################################################

gene_list <- list(Cross = cross_genes,Unmanipulated = unmanip_genes,Exp2 = exp2_genes)
gene_counts <- table(unlist(gene_list))
shared_genes <- names(gene_counts[gene_counts >= 2])
meta_table <- bind_rows(cross_sig %>% filter(gene %in% shared_genes) %>%
                          mutate(Experiment = "Cross"), unmanip_sig %>% filter(gene %in% shared_genes) %>%
                          mutate(Experiment = "Unmanipulated"),
                        exp2_sig %>% filter(gene %in% shared_genes) %>% mutate(Experiment = "Exp2")) %>% arrange(gene)

############################################################
# Consensus DVGs (significant in >=2 of 3 experiments) at two FC thresholds
############################################################

build_consensus_table <- function(fc_cutoff) {
  cross_sig_fc   <- cross_tbl   %>% filter(fdr_brown_forsythe        < 0.05, log2FC_mean_adjusted_variance_OY >= fc_cutoff)
  unmanip_sig_fc <- unmanip_tbl %>% filter(fdr_brown_forsythe        < 0.05, log2FC_mean_adjusted_variance_OY >= fc_cutoff)
  exp2_sig_fc    <- exp2_tbl    %>% filter(fdr_brown_forsythe        < 0.05, log2FC_mean_adjusted_variance_OY >= fc_cutoff)
  
  gene_list_fc <- list(Cross = cross_sig_fc$gene, Unmanipulated = unmanip_sig_fc$gene, Exp2 = exp2_sig_fc$gene)
  gene_counts_fc <- table(unlist(gene_list_fc))
  consensus_genes_fc <- names(gene_counts_fc[gene_counts_fc >= 2])
  
  cross_wide <- cross_tbl %>%
    filter(gene %in% consensus_genes_fc) %>%
    dplyr::rename(fdr_Cross = fdr_brown_forsythe, log2FC_Cross = log2FC_mean_adjusted_variance_OY)
  unmanip_wide <- unmanip_tbl %>%
    filter(gene %in% consensus_genes_fc) %>%
    dplyr::rename(fdr_Unmanipulated = fdr_brown_forsythe, log2FC_Unmanipulated = log2FC_mean_adjusted_variance_OY)
  exp2_wide <- exp2_tbl %>%
    filter(gene %in% consensus_genes_fc) %>%
    dplyr::rename(fdr_Exp2 = fdr_brown_forsythe, log2FC_Exp2 = log2FC_mean_adjusted_variance_OY)
  
  purrr::reduce(
    list(cross_wide, unmanip_wide, exp2_wide),
    dplyr::full_join, by = "gene"
  ) %>%
    mutate(
      n_experiments_significant =
        (!is.na(fdr_Cross) & fdr_Cross < 0.05 & log2FC_Cross >= fc_cutoff) +
        (!is.na(fdr_Unmanipulated) & fdr_Unmanipulated < 0.05 & log2FC_Unmanipulated >= fc_cutoff) +
        (!is.na(fdr_Exp2) & fdr_Exp2 < 0.05 & log2FC_Exp2 >= fc_cutoff)
    ) %>%
    arrange(desc(n_experiments_significant), gene)
}

consensus_2of3_fc0.05 <- build_consensus_table(0.05)
consensus_2of3_fc0.1  <- build_consensus_table(0.1)

############################################################
# Write Excel
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
addWorksheet(wb, "Consensus_2of3_FC0.05")
writeData(wb, "Consensus_2of3_FC0.05", consensus_2of3_fc0.05)
addWorksheet(wb, "Consensus_2of3_FC0.1")
writeData(wb, "Consensus_2of3_FC0.1", consensus_2of3_fc0.1)
saveWorkbook(wb,".../Fig3_DVG_Master_Table.xlsx",overwrite = TRUE)

############################################################


##################################################################

##################################################################

# Visualized per-gene experssion and mean-adj variance

# gene list: consensus DVGs (significant in >=2 of 3 experiments), FC >= 0.05
genes_of_interest_fc0.05 <- consensus_2of3_fc0.05$gene

# gene list: consensus DVGs (significant in >=2 of 3 experiments), FC >= 0.1
genes_of_interest_fc0.1 <- consensus_2of3_fc0.1$gene

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
    p_exp2
  
  return(combined_panel)
}

panel_list_fc0.05 <- lapply(genes_of_interest_fc0.05, make_4panel_for_gene)
names(panel_list_fc0.05) <- genes_of_interest_fc0.05

pdf(".../DensityPanels_DVG_FC0.05_Consensus2of3.pdf", width = 15, height = 12)
for(g in genes_of_interest_fc0.05){
  if(!is.null(panel_list_fc0.05[[g]])){
    print(panel_list_fc0.05[[g]])
  }
}
dev.off()

panel_list_fc0.1 <- lapply(genes_of_interest_fc0.1, make_4panel_for_gene)
names(panel_list_fc0.1) <- genes_of_interest_fc0.1

pdf(".../DensityPanels_DVG_FC0.1_Consensus2of3.pdf", width = 15, height = 12)
for(g in genes_of_interest_fc0.1){
  if(!is.null(panel_list_fc0.1[[g]])){
    print(panel_list_fc0.1[[g]])
  }
}
dev.off()

#####################################################################################

#####################################################################################


# DEG for DVGs

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
  output_dir = ".../Cross_Vitro_DEG_Results"
)

View(deg_vitro$results$DEG2_TopHVGs_MAST)
deg_vitro$volcano_plot

# DEG for Unmanipulated (Old vs Young)
deg_unmanip <- run_clonewise_DEG_suite(
  seurat_obj = seurat_unmanipulated_HSC,
  clone_set1_ids = "old",     # Group1 = Old
  clone_set2_ids = "young",   # Group2 = Young
  clone_col = "AGE",
  group1_label = "Old_unmanip",
  group2_label = "Young_unmanip",
  fc_cutoff  = 0.25,
  fdr_cutoff = 0.05,
  excel_name  = "Unmanipulated__Old_vs_Young_DEG.xlsx",
  volcano_name = "Unmanipulated__Old_vs_Young_DEG_Volcano.png",
  output_dir = ".../Unmanipulated_DEG_Results"
)

View(deg_unmanip$results$DEG2_TopHVGs_MAST)
deg_unmanip$volcano_plot

# DEG for Exp2 (Old vs Young)
deg_exp2 <- run_clonewise_DEG_suite(
  seurat_obj = hsc_exp2,
  clone_set1_ids = c("4_Oa","5_Oa","6_Oa"),   # Group1 = Old
  clone_set2_ids = c("1_Ya","2_Ya","3_Ya"),   # Group2 = Young
  clone_col = "sampleName",
  group1_label = "Old_exp2",
  group2_label = "Young_exp2",
  fc_cutoff  = 0.25,
  fdr_cutoff = 0.05,
  excel_name  = "Exp2__Old_vs_Young_DEG.xlsx",
  volcano_name = "Exp2__Old_vs_Young_DEG_Volcano.png",
  output_dir = ".../Exp2_DEG_Results"
)

View(deg_exp2$results$DEG2_TopHVGs_MAST)
deg_exp2$volcano_plot

#####################################################################################

#####################################################################################

# DEG vs DVG overlap 

extract_deg_genes <- function(deg_df, fdr_cutoff = 0.05, fc_cutoff = 0.25) {
  df <- deg_df
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  fdr_col <- intersect(c("p_val_adj", "padj", "FDR", "fdr"), colnames(df))[1]
  fc_col  <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), colnames(df))[1]
  stopifnot(!is.na(fdr_col), !is.na(fc_col))
  df %>%
    filter(.data[[fdr_col]] < fdr_cutoff, .data[[fc_col]] >= fc_cutoff) %>%
    pull(gene)
}

deg_genes_cross    <- extract_deg_genes(deg_vitro$results$DEG3_FilteredHVGs_MAST)
deg_genes_unmanip  <- extract_deg_genes(deg_unmanip$results$DEG3_FilteredHVGs_MAST)
deg_genes_exp2     <- extract_deg_genes(deg_exp2$results$DEG3_FilteredHVGs_MAST)

# Significant, aging-up DVGs per experiment 

build_overlap_counts <- function(deg_genes, dvg_genes, experiment_label) {
  shared   <- intersect(deg_genes, dvg_genes)
  deg_only <- setdiff(deg_genes, dvg_genes)
  dvg_only <- setdiff(dvg_genes, deg_genes)
  tibble(
    Experiment = experiment_label,
    Category = c("DEG", "Shared", "DVG"),
    Count = c(length(deg_only), length(shared), length(dvg_only))
  )
}

overlap_counts <- bind_rows(
  build_overlap_counts(deg_genes_cross,   dvg_genes_cross,   "Exp1"),
  build_overlap_counts(deg_genes_unmanip, dvg_genes_unmanip, "Unmanipulated"),
  build_overlap_counts(deg_genes_exp2,    dvg_genes_exp2,    "Exp2")
) %>%
  mutate(
    Experiment = factor(Experiment, levels = c("Exp1", "Exp2", "Unmanipulated")),
    # Level order below + default position_stack() puts Shared at the top of
    # the bar, DVG in the middle, DEG at the bottom.
    Category = factor(Category, levels = c("Shared", "DVG", "DEG"))
  )

p_deg_dvg_overlap <- ggplot(overlap_counts, aes(x = Experiment, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "DEG" = "#D95F02",
    "Shared" = "#1B9E77",
    "DVG" = "#7570B3"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text  = element_text(color = "black"),
    legend.title = element_blank()
  ) +
  labs(
    title = "Aging-Up DVG vs DEG Overlap",
    x = NULL,
    y = "Number of genes"
  )

p_deg_dvg_overlap


####################################################################################
# End of this script
