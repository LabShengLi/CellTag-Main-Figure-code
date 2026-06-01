seurat.vivo.exp2 <- readRDS(
  "data/Exp2(exp1)_vivo.RDS"
)
seurat.vitro.exp2 <- readRDS(
  "data/Exp2(exp1)_vitro.RDS"
)
library(Seurat)
library(purrr)
library(ComplexHeatmap)
table(seurat.vivo.exp2$celltype_final)
table(seurat.vitro.exp2$celltype_final)
lineage_map <- c(
  # ====================
  # HSC
  # ====================
  "LT-HSC" = "HSC",
  "ST-HSC" = "HSC",
  "HSC"    = "HSC",
  
  # ====================
  # MPP
  # ====================
  "MPP"  = "MPP",
  "MPP3" = "MPP",
  "MPP4" = "MPP",
  "LMPP-1" = "MPP",
  "LMPP-2" = "MPP",
  # ====================
  # CMP
  # ====================
  "CMP" = "CMP",
  
  # ====================
  # GMP branch
  # ====================
  "GMP"   = "GMP",
  "GMP1"  = "GMP",
  "GMP-1" = "GMP",
  "GMP-2" = "GMP",
  "GMP-3" = "GMP",
  
  # ====================
  # MDP (GMP downstream)
  # ====================
  "MDP" = "MDP",
  
  # ====================
  # MEP branch
  # ====================
  "MEP"   = "MEP",
  "MEP1"  = "MEP",
  "MEP-1" = "MEP",
  "MEP-2" = "MEP",
  "MEP-3" = "MEP",
  "MEP/EryP" = "MEP/EryP",
  
  # ====================
  # MEP downstream
  # ====================
  "MkP"   = "MkP",
  "MKP"   = "MkP",
  "MKP-1" = "MkP",
  "MkP1"  = "MkP",
  
  "EryP" = "EryP"
)
meta_vitro <- seurat.vitro.exp2@meta.data
meta_vitro$celltype_final <- trimws(as.character(meta_vitro$celltype_final))
meta_vitro$lineage <- NA_character_

meta_vitro$lineage[meta_vitro$celltype_final %in% names(lineage_map)] <-
  lineage_map[meta_vitro$celltype_final[meta_vitro$celltype_final %in% names(lineage_map)]]

seurat.vitro.exp2@meta.data <- meta_vitro
table(seurat.vitro.exp2$lineage, useNA = "ifany")
meta_vivo <- seurat.vivo.exp2@meta.data
meta_vivo$celltype_final <- trimws(as.character(meta_vivo$celltype_final))
meta_vivo$lineage <- NA_character_

meta_vivo$lineage[meta_vivo$celltype_final %in% names(lineage_map)] <-
  lineage_map[meta_vivo$celltype_final[meta_vivo$celltype_final %in% names(lineage_map)]]

seurat.vivo.exp2@meta.data <- meta_vivo
table(seurat.vivo.exp2$lineage, useNA = "ifany")
seurat.vitro.exp2.sub <- subset(seurat.vitro.exp2, subset = !is.na(lineage))
seurat.vivo.exp2.sub  <- subset(seurat.vivo.exp2,  subset = !is.na(lineage))

table(seurat.vitro.exp2.sub$lineage)
table(seurat.vivo.exp2.sub$lineage)
common_clones_exp2 <- intersect(
  unique(seurat.vitro.exp2.sub$CloneID_V2),
  unique(seurat.vivo.exp2.sub$CloneID_V2)
)

message("exp2 total clone number = ", length(common_clones_exp2))

common_clones_exp2_young <- common_clones_exp2[grepl("_Y", common_clones_exp2)]
common_clones_exp2_old   <- common_clones_exp2[grepl("_O", common_clones_exp2)]
length(common_clones_exp2_young)
length(common_clones_exp2_old)

head(common_clones_exp2_young)
head(common_clones_exp2_old)
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

# ---------------------- 2.1 fate heatmap  ----------------------
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

vitro_fate_exp2_young <- compute_fate_matrix(
  seurat = seurat.vitro.exp2.sub,
  cluster_col = "lineage",
  larry_col   = "CloneID_V2",
  sample_col  = "orig.ident",
  samples_to_filt = unique(seurat.vitro.exp2.sub$orig.ident),
  min_cells   = 0,
  bcs_to_filt = common_clones_exp2_young,
  normalize_clusters   = TRUE,
  normalize_intraclone = TRUE
)

vivo_fate_exp2_young <- compute_fate_matrix(
  seurat = seurat.vivo.exp2.sub,
  cluster_col = "lineage",
  larry_col   = "CloneID_V2",
  sample_col  = "orig.ident",
  samples_to_filt = unique(seurat.vivo.exp2.sub$orig.ident),
  min_cells   = 0,
  bcs_to_filt = common_clones_exp2_young,
  normalize_clusters   = TRUE,
  normalize_intraclone = TRUE
)
desired_order <- c("HSC", "MPP", "GMP", "MEP", "MkP")
for (fate in list(vitro_fate_exp2_young, vivo_fate_exp2_young)) {
  missing <- setdiff(desired_order, colnames(fate))
  for (m in missing) fate[[m]] <- 0
}

vitro_fate_exp2_young <- vitro_fate_exp2_young[, desired_order, drop = FALSE]
desired_order <- c("HSC", "MPP","CMP", "GMP", "EryP","MkP")
vivo_fate_exp2_young  <- vivo_fate_exp2_young[,  desired_order, drop = FALSE]
fate_matrices <- list(
  Vitro = as.matrix(vitro_fate_exp2_young),
  Vivo  = as.matrix(vivo_fate_exp2_young)
)
stopifnot(
  identical(
    rownames(fate_matrices$Vitro),
    rownames(fate_matrices$Vivo)
  )
)
mat_joint <- cbind(
  Vitro = fate_matrices$Vitro,
  Vivo  = fate_matrices$Vivo
)
dist_mat <- dist(mat_joint, method = "euclidean")


hc_joint <- hclust(dist_mat, method = "ward.D2")
clone_order <- rownames(mat_joint)[hc_joint$order]
Vitro_ord <- fate_matrices$Vitro[clone_order, ]
Vivo_ord  <- fate_matrices$Vivo[clone_order, ]
fate_matrices_ordered <- list(
  Vitro = Vitro_ord,
  Vivo  = Vivo_ord
)
h_exp2_young <- generate_fate_heatmaps(
  fate_matrices = fate_matrices_ordered,
  top_barplot   = TRUE,
  scale_quantile = 0.95,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
h_exp2_young <- generate_fate_heatmaps(
  fate_matrices = fate_matrices_ordered,
  top_barplot   = TRUE,
  scale_quantile = 0.95,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # 🔥 ========== font size = 6 ==========
  column_names_gp = grid::gpar(fontsize = 11),
  row_names_gp    = grid::gpar(fontsize = 11),
  
  # 🔥 ========== theme size = 18 ==========
  column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
  
  # 🔥 legend
  heatmap_legend_param = list(
    title_gp  = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 11)
  )
)
draw(
  h_exp2_young,
  gap = unit(8, "mm")   
)
