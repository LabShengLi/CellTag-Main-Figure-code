library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(grid)
library(purrr)
library(Seurat)
library(dplyr)
library(tidyr)

# ============================================================
# 0. Load 数据
# ============================================================

seurat.vivo.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)
samples_to_plot <- c("OY", "OO")
names(samples_to_plot) <- c("OY", "OO")

filename <- "FateHeatmap_OY_vs_OO_log2sorted_grouped_split.pdf"

# ============================================================
# 1. 全局预处理
# ============================================================

meta <- FetchData(seurat.vivo.cross,
                  c("celltype_final", "CloneID", "orig.ident"))

meta <- meta %>%
  filter(
    tolower(celltype_final) %in% tolower(c("HSC", "MPP", "GMP", "MkP", "EryP")),
    !is.na(CloneID),
    CloneID != 0
  ) %>%
  select(orig.ident, celltype_final, CloneID)
seurat.vivo.cross@meta.data <- meta

# ============================================================
# 2. Fate matrix 函数
# ============================================================
compute_fate_matrix <- function(seurat, cluster_col, larry_col, sample_col,
                                samples_to_filt, min_cells,
                                bcs_to_filt=NULL,
                                normalize_clusters=TRUE,
                                normalize_intraclone=TRUE,
                                normalize_pre_subset=FALSE) {
  df <- FetchData(seurat, c(cluster_col, larry_col, sample_col))
  df <- df[df[[sample_col]] %in% samples_to_filt, , drop=FALSE]
  if (nrow(df) == 0) return(matrix(0,0,0))
  
  clone_counts <- table(df[[larry_col]])
  big_clones <- names(clone_counts[clone_counts >= min_cells])
  df <- df[df[[larry_col]] %in% big_clones, , drop=FALSE]
  if (nrow(df) == 0) return(matrix(0,0,0))
  
  clonal_mat <- as.matrix(table(df[[larry_col]], df[[cluster_col]]))
  
  if (!is.null(bcs_to_filt) & normalize_pre_subset) {
    clonal_mat <- clonal_mat[rownames(clonal_mat) %in% bcs_to_filt,, drop=FALSE]
    clonal_mat <- clonal_mat[bcs_to_filt,, drop=FALSE]
  }
  
  if (normalize_clusters) {
    col_sums <- colSums(clonal_mat); col_sums[col_sums==0] <- 1
    clonal_mat <- sweep(clonal_mat, 2, col_sums, "/")
  }
  
  if (normalize_intraclone) {
    row_sums <- rowSums(clonal_mat); row_sums[row_sums==0] <- 1
    clonal_mat <- sweep(clonal_mat, 1, row_sums, "/")
  }
  
  if (!is.null(bcs_to_filt) & !normalize_pre_subset) {
    clonal_mat <- clonal_mat[rownames(clonal_mat) %in% bcs_to_filt,, drop=FALSE]
    clonal_mat <- clonal_mat[bcs_to_filt,, drop=FALSE]
  }
  return(clonal_mat)
}

get_fate_matrices <- function(seurat, samples_to_select, ...) {
  fate_mat_l <- lapply(unname(samples_to_select),
                       function(sample) compute_fate_matrix(seurat, samples_to_filt=sample, ...))
  nm <- names(samples_to_select)
  if (!is.null(nm)) names(fate_mat_l) <- nm else names(fate_mat_l) <- samples_to_select
  return(fate_mat_l)
}

# ---------------------- 2.1 fate heatmap 函数 (只在第一个 sample 显示标签) ----------------------
generate_fate_heatmaps <- function(
    fate_matrices,
    top_barplot    = TRUE,
    color_function = NULL,
    scale_quantile = 1,
    ...
) {
  if (top_barplot) {
    percent_bcs_in_fate <- map(fate_matrices, \(m) colSums(m > 0)/nrow(m)*100)
    
    top_heatmap_annot_l <- imap(
      percent_bcs_in_fate,
      \(x, nm) {
        if (nm == names(fate_matrices)[1]) {
          HeatmapAnnotation(
            `% bc in fate` = anno_barplot(
              x, beside = TRUE, gp = gpar(fill = "black"), height = unit(1, "cm")
            ),
            show_annotation_name = TRUE,
            annotation_name_side = "left"
          )
        } else {
          HeatmapAnnotation(
            `% bc in fate` = anno_barplot(
              x, beside = TRUE, gp = gpar(fill = "black"), height = unit(1, "cm")
            ),
            show_annotation_name = FALSE
          )
        }
      }
    )
  }
  
  if (is.null(color_function)) {
    sample_color_scale <- map(
      fate_matrices,
      \(x) circlize::colorRamp2(c(min(x), quantile(x, scale_quantile)),
                                hcl_palette = "Rocket", reverse = TRUE)
    )
  } else {
    sample_color_scale <- map(fate_matrices, color_function)
  }
  
  if (top_barplot) {
    fate_heatmap_l <- imap(
      fate_matrices,
      \(fate_mat, sample) Heatmap(
        fate_mat, 
        name             = "Intra-clone fraction", 
        col              = sample_color_scale[[sample]],
        column_title     = sample,
        top_annotation   = top_heatmap_annot_l[[sample]],
        ...
      )
    )
  } else {
    fate_heatmap_l <- imap(
      fate_matrices,
      \(fate_mat, sample) Heatmap(
        fate_mat, 
        name             = "Intra-clone fraction", 
        col              = sample_color_scale[[sample]],
        column_title     = sample,
        ...
      )
    )
  }
  
  h_l <- NULL
  for (hm in fate_heatmap_l) { h_l <- h_l + hm }
  return(h_l)
}
generate_cloneSize_heatmaps <- function(
    size_matrices,
    color_function = NULL,
    scale_quantile = 1,
    ...
) {
  
  #----- Define color scales if not specified
  if (is.null(color_function)) {
    
    sample_color_scale <- map(
      size_matrices,
      \(x) circlize::colorRamp2(c(min(x), quantile(x, scale_quantile)), hcl_palette = "Rocket", reverse = TRUE)
    )
    
  } else { sample_color_scale <- map(size_matrices, color_function) }
  
  
  #----- Prepare heatmaps
  size_heatmap_l <- imap(
    size_matrices,
    \(fate_mat, sample) Heatmap(fate_mat, 
                                name             = paste0(sample, " clone size"), 
                                col              = sample_color_scale[[sample]],
                                ...
    )
  )
  
  # Merge all heatmaps into 1 list and finish
  h_l <- NULL
  for (hm in size_heatmap_l) { h_l <- h_l + hm }
  
  return(h_l)
  
}
# ============================================================
# 3. Fate 矩阵
# ============================================================

l_fate_matrices <- get_fate_matrices(
  seurat = seurat.vivo.cross,
  cluster_col = "celltype_final",
  larry_col = "CloneID",   # 改这里
  sample_col = "orig.ident",
  samples_to_select = samples_to_plot,
  min_cells = 1
)

common_clones <- intersect(
  rownames(l_fate_matrices[[1]]),
  rownames(l_fate_matrices[[2]])
)

common_clones <- setdiff(common_clones, "0")

l_fate_matrices <- lapply(
  l_fate_matrices,
  function(x) x[common_clones,, drop=FALSE]
)

keep_types <- c("HSC", "MPP", "GMP", "MkP", "EryP")

l_fate_matrices <- lapply(l_fate_matrices, function(m) {
  cols_to_keep <- colnames(m) %in% keep_types
  m[, cols_to_keep, drop = FALSE]
})

# ============================================================
# 4. log2FC (HSC)
# ============================================================

log2fc_vec <- log2(
  (l_fate_matrices[["OO"]][,"HSC"]+0.01) /
    (l_fate_matrices[["OY"]][,"HSC"]+0.01)
)

groups_raw <- cut(log2fc_vec,
                  breaks = c(-Inf, -1, 1, Inf),
                  labels = c("lose HSCs", "unchanged", "gain HSCs"))

groups <- factor(groups_raw,
                 levels = c("gain HSCs", "unchanged", "lose HSCs"))

df_order <- data.frame(
  CloneID = names(log2fc_vec),
  log2fc = log2fc_vec,
  group = groups
)

df_order <- df_order[order(df_order$group, -df_order$log2fc), ]

log2fc_vec <- df_order$log2fc
names(log2fc_vec) <- df_order$CloneID
groups <- df_order$group

# 同步矩阵
l_fate_matrices <- lapply(
  l_fate_matrices,
  function(m) m[df_order$CloneID, , drop = FALSE]
)

# ============================================================
# 5. Clone size
# ============================================================

l_size_matrices <- get_fate_matrices(
  seurat=seurat.vivo.cross,
  cluster_col="celltype_final",
  larry_col="CloneID",
  sample_col="orig.ident",
  samples_to_select=samples_to_plot,
  bcs_to_filt=common_clones,
  min_cells=1,
  normalize_clusters=FALSE,
  normalize_intraclone=FALSE
) %>% lapply(rowSums)

l_size_matrices <- lapply(
  l_size_matrices,
  function(v) v[df_order$CloneID]
)

# ============================================================
# 6. Replicate 统计
# ============================================================

seurat.vivo.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)
md <- seurat.vivo.cross@meta.data %>%
  filter(CloneID %in% df_order$CloneID,
         orig.ident %in% c("OO","OY"))

df_counts <- md %>%
  group_by(CloneID, orig.ident, Rep) %>%
  summarise(n_cells = n(), .groups = "drop")

clone_rep_metrics <- df_counts %>%
  unite("SampleRep", orig.ident, Rep, sep = "_") %>%
  pivot_wider(names_from = SampleRep,
              values_from = n_cells,
              values_fill = 0) %>%
  slice(match(df_order$CloneID, CloneID))

choose_major_rep <- function(counts) {
  reps <- c("Rep1","Rep2","Rep3")
  reps[which.max(counts)]
}

clone_rep_metrics$OO_majorRep <- apply(
  clone_rep_metrics[, c("OO_Rep1","OO_Rep2","OO_Rep3")],
  1, choose_major_rep
)

clone_rep_metrics$OY_majorRep <- apply(
  clone_rep_metrics[, c("OY_Rep1","OY_Rep2","OY_Rep3")],
  1, choose_major_rep
)

# ============================================================
# 7. Heatmap
# ============================================================

col_log2 <- circlize::colorRamp2(c(-4,0,4), c("blue","white","red"))
rep_col <- c("Rep1"="#1f77b4",
             "Rep2"="#2ca02c",
             "Rep3"="#ff7f0e")

fate_heatmaps <- generate_fate_heatmaps(
  l_fate_matrices,
  top_barplot=TRUE,
  show_row_names=FALSE,
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  border=TRUE,
  row_split=groups,
  row_title_gp=gpar(fontsize=12, fontface="bold"),
  row_title_rot=0,
  row_title_side="left"
)

size_heatmaps <- generate_cloneSize_heatmaps(
  l_size_matrices,
  scale_quantile=0.95,
  show_row_names=FALSE,
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  border=TRUE,
  row_split=groups,
  row_title=NULL
)

hsc_log2fc <- Heatmap(
  matrix(log2fc_vec, ncol=1,
         dimnames=list(names(log2fc_vec), "Log2-HSC")),
  name="Log2-HSC",
  col=col_log2,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

oy_rep_mat <- matrix(clone_rep_metrics$OY_majorRep, ncol=1,
                     dimnames=list(df_order$CloneID_clean, "OY Rep"))

oo_rep_mat <- matrix(clone_rep_metrics$OO_majorRep, ncol=1,
                     dimnames=list(df_order$CloneID_clean, "OO Rep"))

oy_rep_heatmap <- Heatmap(
  oy_rep_mat,
  name="OY Rep",
  col=rep_col,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

oo_rep_heatmap <- Heatmap(
  oo_rep_mat,
  name="OO Rep",
  col=rep_col,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

# ============================================================
# 8. 输出
# ============================================================

ht_list <- fate_heatmaps +
  size_heatmaps +
  hsc_log2fc

pdf(file.path("/project2/sli68423_1316/users/Kailiang/Test_Rcode/3.17",
              filename),
    height=8, width=6)

draw(ht_list, ht_gap=unit(6,"mm"))
dev.off()
##############################
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(grid)
library(purrr)
library(Seurat)
library(dplyr)
library(tidyr)

# ---------------------- 0. Load 数据 ----------------------
seurat.vivo.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)

samples_to_plot <- c("YY", "YO")
names(samples_to_plot) <- c("YY", "YO")

filename <- "FateHeatmap_YY_vs_YO_log2sorted_grouped_split.pdf"

# ---------------------- 1. 全局预处理 ----------------------
meta <- FetchData(seurat.vivo.cross,
                  c("celltype_final", "CloneID", "orig.ident"))
meta <- meta %>%
  filter(
    tolower(celltype_final) %in% tolower(c("HSC", "MPP", "GMP", "MkP", "EryP")),
    !is.na(CloneID),
    CloneID != 0
  ) %>%
  select(orig.ident, celltype_final, CloneID)
seurat.vivo.cross@meta.data <- meta


# ---------------------- 2. Fate matrix ----------------------
l_fate_matrices <- get_fate_matrices(
  seurat=seurat.vivo.cross,
  cluster_col="celltype_final",
  larry_col="CloneID",
  sample_col="orig.ident",
  samples_to_select=samples_to_plot,
  min_cells=1
)

# 公共克隆
common_clones <- intersect(rownames(l_fate_matrices[[1]]), rownames(l_fate_matrices[[2]]))
common_clones <- setdiff(common_clones, "0")
#common_clones <- setdiff(common_clones, ambiguous_young)

l_fate_matrices <- lapply(l_fate_matrices, function(x) x[common_clones,, drop=FALSE])
keep_types <- c("HSC", "MPP", "GMP", "MkP", "EryP")

l_fate_matrices <- lapply(l_fate_matrices, function(m) {
  cols_to_keep <- colnames(m) %in% keep_types
  m[, cols_to_keep, drop = FALSE]
})

# ---------------------- log2FC ----------------------
log2fc_vec <- log2(
  (l_fate_matrices[["YO"]][,"HSC"]+0.01) /
    (l_fate_matrices[["YY"]][,"HSC"]+0.01)
)

groups_raw <- cut(log2fc_vec,
                  breaks = c(-Inf, -1, 1, Inf),
                  labels = c("lose HSCs", "unchanged", "gain HSCs"))

groups <- factor(groups_raw, levels = c("gain HSCs", "unchanged", "lose HSCs"))

df_order <- data.frame(
  CloneID = names(log2fc_vec),
  log2fc  = log2fc_vec,
  group   = groups
)

df_order <- df_order[order(df_order$group, -df_order$log2fc), ]

log2fc_vec <- df_order$log2fc
names(log2fc_vec) <- df_order$CloneID
groups <- df_order$group

# ---------------------- 同步矩阵 ----------------------
l_fate_matrices <- lapply(l_fate_matrices,
                          function(m) m[df_order$CloneID, , drop = FALSE])

# ---------------------- 重新加载数据 ----------------------
seurat.vivo.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)

md <- seurat.vivo.cross@meta.data %>%
  filter(CloneID %in% df_order$CloneID,
         orig.ident %in% c("YO","YY"))

# ---------------------- clone × replicate 统计 ----------------------
df_counts <- md %>%
  group_by(CloneID, orig.ident, Rep) %>%
  summarise(n_cells = n(), .groups = "drop")

clone_rep_metrics <- df_counts %>%
  unite("SampleRep", orig.ident, Rep, sep = "_") %>%
  pivot_wider(names_from = SampleRep,
              values_from = n_cells,
              values_fill = 0) %>%
  slice(match(df_order$CloneID, CloneID))

choose_major_rep <- function(counts) {
  reps <- c("Rep1","Rep2","Rep3")
  max_val <- max(counts)
  max_reps <- reps[counts == max_val]
  
  if (length(max_reps) == 1) {
    return(max_reps)
  } else if (all(c("Rep1","Rep2") %in% max_reps)) {
    return("Rep1")
  } else if (all(c("Rep2","Rep3") %in% max_reps)) {
    return("Rep2")
  } else if (all(c("Rep1","Rep3") %in% max_reps)) {
    return("Rep1")
  } else {
    return(max_reps[1])  
  }
}

clone_rep_metrics$YO_majorRep <- apply(
  clone_rep_metrics[, c("YO_Rep1","YO_Rep2","YO_Rep3")],
  1, choose_major_rep
)

clone_rep_metrics$YY_majorRep <- apply(
  clone_rep_metrics[, c("YY_Rep1","YY_Rep2","YY_Rep3")],
  1, choose_major_rep
)

# ---------------------- size matrix ----------------------
l_size_matrices <- get_fate_matrices(
  seurat=seurat.vivo.cross,
  cluster_col="celltype_final",
  larry_col="CloneID",
  sample_col="orig.ident",
  samples_to_select=samples_to_plot,
  bcs_to_filt=common_clones,
  min_cells=1,
  normalize_clusters=FALSE,
  normalize_intraclone=FALSE
) %>% lapply(rowSums)

l_size_matrices <- lapply(l_size_matrices,
                          function(v) v[df_order$CloneID])

# ---------------------- Heatmap ----------------------
col_log2 <- circlize::colorRamp2(c(-4,0,4), c("blue","white","red"))
rep_col <- c("Rep1" = "#1f77b4", "Rep2" = "#2ca02c", "Rep3" = "#ff7f0e")

fate_heatmaps <- generate_fate_heatmaps(
  l_fate_matrices,
  top_barplot=TRUE,
  show_row_names=FALSE,
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  border=TRUE,
  row_split=groups,
  row_title_gp=gpar(fontsize=12, fontface="bold"),
  row_title_rot=0,
  row_title_side="left"
)

size_heatmaps <- generate_cloneSize_heatmaps(
  l_size_matrices,
  scale_quantile=0.95,
  show_row_names=FALSE,
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  border=TRUE,
  row_split=groups,
  row_title=NULL
)

hsc_log2fc <- Heatmap(
  matrix(log2fc_vec, ncol=1,
         dimnames=list(names(log2fc_vec), "Log2-HSC")),
  name="Log2-HSC",
  col=col_log2,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

yy_rep_vec <- clone_rep_metrics$YY_majorRep[
  match(df_order$CloneID, clone_rep_metrics$CloneID)
]

yo_rep_vec <- clone_rep_metrics$YO_majorRep[
  match(df_order$CloneID, clone_rep_metrics$CloneID)
]

yy_rep_mat <- matrix(yy_rep_vec, ncol=1,
                     dimnames=list(df_order$CloneID, "YY Rep"))

yo_rep_mat <- matrix(yo_rep_vec, ncol=1,
                     dimnames=list(df_order$CloneID, "YO Rep"))

yy_rep_heatmap <- Heatmap(
  yy_rep_mat,
  name="YY Rep",
  col=rep_col,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

yo_rep_heatmap <- Heatmap(
  yo_rep_mat,
  name="YO Rep",
  col=rep_col,
  show_row_names=FALSE,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  border=TRUE,
  row_split=groups
)

# ---------------------- 输出 ----------------------
ht_list <- fate_heatmaps + size_heatmaps + hsc_log2fc

pdf(file.path("/project2/sli68423_1316/users/Kailiang/Test_Rcode/2.27",
              filename),
    height=8, width=6)

draw(ht_list, ht_gap=unit(6,"mm"))
dev.off()
