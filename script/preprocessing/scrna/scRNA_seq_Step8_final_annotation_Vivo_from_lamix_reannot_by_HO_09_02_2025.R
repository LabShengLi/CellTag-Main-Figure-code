#!/usr/bin/env Rscript
rm(list = ls())
library(circlize)
library("Seurat")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)

wdir <- 'Harmony_integration_Cutoff_1K/'
DataDir <- 'Harmony_integration_Cutoff_1K/'

infn <- "cross.exp.combined.vivo.hto.celltag.annotation_by_lamis_object.rds"

############ load Seurat object
ss <- readRDS(paste0(DataDir, "/", infn))
ss

colnames(ss@meta.data)

table(ss$sampleName)

#   OO   OY   YO   YY
# 8093 1514 3083 5151

table(ss$celltype)

#     LT-HSC      ST-HSC       LMPP1         MPP
#       4031        2529        1415        5530
#        GMP         MkP        EryP          DC
#       1772         466         510         202
#        MAC     T-Cells     B-cells Granulocyte
#         92          75          88         931
# Mast-Cells
#        200


# comments from HO:  Since Marco uses "MPP4" for his presentations, maybe we should use "MPP4" instead of "LMPP". In this case, the current "MPP" on the heatmap can be "MPP3"
ss$celltype <- as.character(ss$celltype)
ss$celltype[ss$celltype == "MPP"] <- "MPP3"
ss$celltype[ss$celltype == "LMPP1"] <- "MPP4"

table(ss$celltype)

lineage_order <- c(
	"LT-HSC", "ST-HSC", "MPP3", "MPP4",
	"GMP", "MkP", "EryP",
	"DC", "MAC", "T-Cells", "B-cells",
	"Granulocyte", "Mast-Cells"
)

# Make sure the metadata column is a factor with correct order
ss$celltype <- factor(
	ss$celltype,
	levels = lineage_order
)


dp3 <- DimPlot(ss, group.by = "celltype", label = T, label.size = 5)
print(dp3)


outfn = sprintf("%s_umap_seurat_clusters_corrected_by_HO_09_02_2025.pdf", "cross_vivo")
graphics.off()
pdf(file.path(wdir, outfn), width = 6, height = 6)
print(dp3)
dev.off()

####################
####################
# read marker file
marker_fn <- 'chip_mouse_hsc_cell_types_v4.tsv'
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


# include Cd34 genes
assay_to_use <- DefaultAssay(ss)
ss <- ScaleData(ss, assay = assay_to_use, features = rownames(ss[[assay_to_use]]), verbose = FALSE)
ht1 <- plot_heatmap_on_markers(ss, 'celltype', marker_list)
draw(ht1)

outfn = sprintf("%s_heatmap_seurat_clusters_known_markers_corrected_by_HO_09_02_2025.pdf", "cross_vivo")
graphics.off()
pdf(file.path(wdir, outfn), width = 6, height = 10)
draw(ht1)
dev.off()

ss$celltype_HO_09_02_2025 <- ss$celltype


ss$celltype_group <- as.character(ss$celltype)


# Create new grouping column based on celltype
ss$celltype_group <- dplyr::recode(ss$celltype_group,
								   "LT-HSC" = "HSC",
								   "ST-HSC" = "HSC",
								   "DC" = "Differentiated_cells",
								   "MAC" = "Differentiated_cells",
								   "T-Cells" = "Differentiated_cells",
								   "B-cells" = "Differentiated_cells",
								   "Granulocyte" = "Differentiated_cells",
								   "Mast-Cells" = "Differentiated_cells")


lineage_order <- c(
	"HSC", "MPP3", "MPP4",
	"GMP", "MkP", "EryP",
	"Differentiated_cells"
)

# Make sure the metadata column is a factor with correct order
ss$celltype_group <- factor(
	ss$celltype_group,
	levels = lineage_order
)

# Check the result
table(ss$celltype_group)

#                HSC               MPP3               MPP4
#               6560               5530               1415
#                GMP                MkP               EryP
#               1772                466                510
# Differentiated_cells
#               1588

saveRDS(ss, file = paste0(wdir, "/cross.exp.combined.vivo.hto.celltag.annotation_by_lamis_object.correct_by_HO_09_02_2025.rds"))

cat("### Done \n")


quit()


ss_130 <- subset(ss, subset = CloneID == 'Old_C_130')


ss_oy <- subset(ss, subset = orig.ident == 'OY')
ss_oo <- subset(ss, subset = orig.ident == 'OO')

ss_oy_clones<- ss_oy$CloneID
ss_oo_clones<- ss_oo$CloneID

old_common_clones <- intersect(ss_oy_clones, ss_oo_clones)


intersect(old_common_clones, 'Old_C_152')