#!/usr/bin/env Rscript
rm(list = ls())
library(circlize)
library("Seurat")
library(ComplexHeatmap)
library(ggplot2)

wdir <- '/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/PreProcessing/Harmony_integration_Cutoff_1K/'
DataDir <- '/project2/sli68423_1316/projects/U01_aim2/Cross_Expirement/PreProcessing/Harmony_integration_Cutoff_1K/'

infn <- "cross.exp.combined.vivo.hto.celltag.rds"

############ load Seurat object
ss <- readRDS(paste0(DataDir, "/", infn))
ss

table(ss$sampleName)

#   OO   OY   YO   YY
# 8093 1514 3083 5151


dp1 <- DimPlot(ss, group.by = "seurat_clusters", label = T, label.size = 5)
print(dp1)

outfn = sprintf("%s_umap_seurat_clusters.pdf", "cross_vivo")
graphics.off()
pdf(file.path(wdir, outfn), width = 6, height = 6)
print(dp1)
dev.off()

####################
####################
# read marker file
marker_fn <- '/project2/sli68423_1316/users/yang/workspace/nanome/related_project/single_cell/chip_mouse_hsc_cell_types_v4.tsv'
marker_df <- read.table(marker_fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(marker_df)
dim(marker_df)


# heatmap of groups
# Compute average expression by cell type (CT)
marker_list <- marker_df$Marker

plot_heatmap_on_markers <- function(ss, clusterName, marker_list) {

	avg_exp <- AverageExpression(
		ss,
		features = marker_list,
		group.by = clusterName,
		assays = "RNA",
		slot = "scale.data"
	)$RNA

	common_genes <- intersect(marker_list, rownames(avg_exp))
	avg_exp <- avg_exp[common_genes,]

	# Scale rows (genes) to z-scores for better visualization
	avg_exp_scaled <- t(scale(t(avg_exp)))

	# Set color scale
	col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

	ht1 <- Heatmap(
		avg_exp_scaled,
		name = "Scaled Exp",
		col = col_fun,
		cluster_rows = F,
		cluster_columns = F,
		row_names_gp = gpar(fontsize = 12, fontface = "italic"),
		column_names_gp = gpar(fontsize = 12),
		heatmap_legend_param = list(title = "Z-score"),
		column_title = "Average Expression Group By Cell Type",
		row_title = "Known markers"
	)
	return(ht1)
}

ht1 <- plot_heatmap_on_markers(ss, 'seurat_clusters', marker_list)

draw(ht1)

outfn = sprintf("%s_heatmap_seurat_clusters_known_markers.pdf", "cross_vivo")
graphics.off()
pdf(file.path(wdir, outfn), width = 6, height = 10)
draw(ht1)
dev.off()


ss$celltype <- "Unknown"
ss$celltype[which(ss$seurat_clusters == "0")] <- "LT-HSC"
ss$celltype[which(ss$seurat_clusters == "1")] <- "ST-HSC"
ss$celltype[which(ss$seurat_clusters == "2")] <- "LT-HSC"
ss$celltype[which(ss$seurat_clusters == "3")] <- "ST-HSC"
ss$celltype[which(ss$seurat_clusters == "4")] <- "GMP"
ss$celltype[which(ss$seurat_clusters == "5")] <- "MPP"
ss$celltype[which(ss$seurat_clusters == "6")] <- "Granulocyte"
ss$celltype[which(ss$seurat_clusters == "7")] <- "ST-HSC"
ss$celltype[which(ss$seurat_clusters == "8")] <- "EryP"
ss$celltype[which(ss$seurat_clusters == "9")] <- "MEP"
ss$celltype[which(ss$seurat_clusters == "10")] <- "ST-HSC"
ss$celltype[which(ss$seurat_clusters == "11")] <- "DC"
ss$celltype[which(ss$seurat_clusters == "12")] <- "GMP"
ss$celltype[which(ss$seurat_clusters == "13")] <- "MAC"
ss$celltype[which(ss$seurat_clusters == "14")] <- "B-Cells"
ss$celltype[which(ss$seurat_clusters == "15")] <- "T-Cells"

colnames(ss[[]])
#  [1] "orig.ident"          "nCount_RNA"
#  [3] "nFeature_RNA"        "barcode"
#  [5] "Size_Factor"         "phases"
#  [7] "G1_score"            "G2M_score"
#  [9] "S_score"             "num_genes_expressed"
# [11] "Total_mRNAs"         "percentMito"
# [13] "percent.mt"          "percent.rb"
# [15] "Cells"               "orig.ident.y"
# [17] "nCount_RNA.y"        "nFeature_RNA.y"
# [19] "Cells.y"             "Library"
# [21] "percent.mt.y"        "percent.rb.y"
# [23] "RNA_snn_res.0.5"     "seurat_clusters"
# [25] "pANNcomputed"        "Doublet_LowConf"
# [27] "Doublet"             "DoubletFinder"
# [29] "SequencedCells"      "QCpass"
# [31] "PassPercent"         "RNA_snn_res.0.1"
# [33] "RNA_snn_res.0.2"     "numclust"
# [35] "predicted"           "sampleName"
# [37] "HTO_tag"             "HTO_rep"
# [39] "cloneTraceCombined"  "celltype"

lineage_order <- c(
	"LT-HSC", "ST-HSC", "MPP", "GMP",
	"MEP", "EryP", "DC", "MAC",
	"T-Cells", "B-Cells", "Granulocyte"
)

ss$celltype <- factor(ss$celltype, levels = lineage_order)

ss$Rep <- ss$HTO_rep
ss$CloneID <- ss$cloneTraceCombined

saveRDS(ss, file = paste0(wdir, "/cross.exp.combined.vivo.hto.celltag.annotation.rds"))

# save celltype annotation table
an_df <- ss@meta.data[, c('barcode', 'celltype')]

# Save as TSV (no row names)
write.table(
	an_df,
	file = "cross_exp_vivo_annotation_an_df.tsv",
	sep = "\t",
	quote = FALSE,
	row.names = T
)

# Save as RDS
saveRDS(an_df, file = "cross_exp_vivo_annotation_an_df.RDS")

graphics.off()
pdf(paste0(wdir, "UMAP_celltype_vitro.pdf"))
DimPlot(ss, group.by = "celltype", label = T, label.size = 5)
dev.off()


####################
####################
# read marker file
marker_fn <- '/project2/sli68423_1316/users/yang/workspace/nanome/related_project/single_cell/chip_mouse_hsc_cell_types_v4.tsv'
marker_df <- read.table(marker_fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(marker_df)
dim(marker_df)

marker_list <- marker_df$Marker
marker_list <- marker_list[1:39]

# dot plot for markers
dp1 <- DotPlot(object = ss, features = marker_list, group.by = 'celltype') +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"), plot.title = element_text(hjust = 0.5)) +
	ggtitle("scRNA")
dp1


outfn = sprintf("%s_dotplot_known_markers.pdf", "cross_vivo")
graphics.off()
pdf(file.path(outdir, outfn), width = 10, height = 5)
print(dp1)
dev.off()


# heatmap of groups
# Compute average expression by cell type (CT)
marker_list <- marker_df$Marker

avg_exp <- AverageExpression(
	ss,
	features = marker_list,
	group.by = "celltype",
	assays = "RNA",
	slot = "scale.data"
)$RNA

common_genes <- intersect(marker_list, rownames(avg_exp))
avg_exp <- avg_exp[common_genes,]

# Scale rows (genes) to z-scores for better visualization
avg_exp_scaled <- t(scale(t(avg_exp)))

# Set color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ht1 <- Heatmap(
	avg_exp_scaled,
	name = "Scaled Exp",
	col = col_fun,
	cluster_rows = F,
	cluster_columns = F,
	row_names_gp = gpar(fontsize = 12, fontface = "italic"),
	column_names_gp = gpar(fontsize = 12),
	heatmap_legend_param = list(title = "Z-score"),
	column_title = "Average Expression Group By Cell Type",
	row_title = "Known markers"
)
draw(ht1)


outfn = sprintf("%s_heatmap_celltype_known_markers.pdf", "cross_vivo")
graphics.off()
pdf(file.path(wdir, outfn), width = 6, height = 10)
draw(ht1)
dev.off()

cat("### Done \n")
