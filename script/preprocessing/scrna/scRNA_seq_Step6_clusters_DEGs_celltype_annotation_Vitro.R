#!/usr/bin/env Rscript
rm(list = ls())

library(EnhancedVolcano)

lnWD <- 'PreProcessing/'

DataDir <- 'Harmony_integration_Cutoff_1K/'

RDSname <- "cross.exp.combined.vitro.hto.celltag.rds"

outdir <- paste0(lnWD, "biology"); dir.create(outdir, showWarnings = FALSE)

library(ggplot2)

source("common/funForLoading_PK.R")
source("common/data_preparation.R")
source("common/castle.R")
source("common/Useful_Functions.R")
source("common/funForSeurat.R")
source("common/plot_figures.R")

############ load Seurat object
hspc.combined <- readRDS(paste0(DataDir, "/", RDSname))
hspc.combined

library("Seurat")

graphics.off()
pdf(paste0(lnWD, "UMPAP_clusters_vitro.pdf"))
DimPlot(hspc.combined, group.by = "seurat_clusters", cols = ClusPalette, label = T, label.size = 5)
dev.off()

graphics.off()
pdf(paste0(lnWD, "UMPAP_ages_vitro.pdf"))
DimPlot(hspc.combined, group.by = "sampleName", cols = Light.Pallette)
dev.off()


################ Composition 
library(scales)
minPropCellExp = 0.001 #"character", "proportion minimal of cells that expressed the genes kept for the analaysis 0.001 by default",
logfc_threshold = 0.25 #"numeric", "logfc threshold for finding cluster markers (1 cluster vs all deg) 0.25 by default",
norm_method = "logNorm"

hspc.combined@meta.data$AGE <- hspc.combined@meta.data$sampleName
table(hspc.combined@meta.data$AGE, hspc.combined@meta.data$sampleName)


table(hspc.combined@meta.data$AGE)

library("plyr")
age <- getAGEPropPerClustBarplot(hspc.combined)


getSamplePropPerClustBarplot(hspc.combined)
getPhasePropPerClustBarplot(hspc.combined)
getPredictedPropPerClustBarplot(hspc.combined)
getAGEPropPerClustBarplot(hspc.combined)


############## Marker genes
Idents(hspc.combined) <- "seurat_clusters"

markers <- FindAllMarkers(hspc.combined, only.pos = T, logfc.threshold = logfc_threshold)
markers <- markers[which(markers$p_val_adj < 0.05),]
write.table(x = markers, paste(outdir, "/markers_vitro.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)
markers <- read.table(paste0(outdir, "/markers_vitro.tsv"), header = T)

library("dplyr")
top <- markers %>%
	group_by(cluster) %>%
	top_n(n = 10, wt = avg_log2FC); top <- top[!duplicated(top$gene),]; dim(top);
write.csv(top, file = paste0(lnWD, "top10_vitro.csv"))

Idents(hspc.combined) <- "seurat_clusters"
SCdata.temp.Heatmap <- subset(hspc.combined, downsample = 100)
table(SCdata.temp.Heatmap@meta.data[, "seurat_clusters"])
SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = top$gene, scale.max = 2)
dtype = "scale.data"
nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
nor.exp <- nor.exp[top$gene,]

SCdata.meta <- SCdata.temp.Heatmap@meta.data[, c("seurat_clusters"), drop = FALSE]
SCdata.meta$seurat_clusters <- factor(SCdata.meta$seurat_clusters, levels = sort(unique(SCdata.meta$seurat_clusters)))
SCdata.meta.order <- SCdata.meta[order(SCdata.meta$seurat_clusters), , drop = FALSE]
nor.exp <- nor.exp[, rownames(SCdata.meta.order)]

seurat_clusters.FULL = ClusPalette; names(seurat_clusters.FULL) <- sort(unique(SCdata.meta$seurat_clusters)); Cluster.FULL <- list(Cluster = seurat_clusters.FULL[!is.na(names(seurat_clusters.FULL))]); Cluster.FULL

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(viridis)
colors <- c(seq(-2, 2, by = 0.01))
my_palette <- c(colorRampPalette(colors = c("darkblue", "#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred"))(n = length(colors)))
ann_colors = Cluster.FULL

library("pheatmap")

graphics.off()
pdf(file = paste0(lnWD, "top10_vitro.pdf"), height = 14)

print(pheatmap(nor.exp[top$gene,], annotation_col = SCdata.meta, show_colnames = F, show_rownames = T, color = my_palette, breaks = colors,
			   annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap"), fontsize = 13))

dev.off()


############### cell type identification ##########

### this long list from Parveen
markers <- read.table(file = "markers_orderbyCT.txt", sep = "\t", header = F)
markers <- markers$V1

### this long list from chenx
markers <- read.table(file = "common/markers_orderbyCT_8S.txt", sep = "\t", header = F)[, 2]

marker_fn <- 'chip_mouse_hsc_cell_types_v4.tsv'
markers <- unique(read.table(file = marker_fn, sep = "\t", header = F)[, 3])
markers

### problem with gene number 59
# markers<-markers$V2[1:58]

DefaultAssay(hspc.combined) <- "RNA"
hspc.combined.data <- GetAssayData(hspc.combined, slot = "data", assay = "RNA")
markers <- intersect(markers, rownames(hspc.combined.data))
hspc.combined.data.markers <- hspc.combined.data[markers,]
hspc.combined.data.markers <- as.matrix(hspc.combined.data.markers)

res = "seurat_clusters"
cts <- sort(unique(hspc.combined@meta.data[, res])); cts
cts <- as.character(cts)
markers_meanexp <- NULL

for (i in cts) {
	#i="C01"
	print(paste0("Processing ", i))
	Idents(hspc.combined) <- res
	tmp_cells <- subset(hspc.combined, idents = c(i))
	print(table(tmp_cells@meta.data[, res]))
	print(dim(tmp_cells))
	tmp_cells <- colnames(tmp_cells)
	tmp_data <- hspc.combined.data.markers[, tmp_cells]
	tmp_mean <- apply(tmp_data, 1, mean)
	markers_meanexp <- cbind(markers_meanexp, tmp_mean)
}

colnames(markers_meanexp) <- cts

markers_meanexp_df <- as.data.frame(markers_meanexp)

# min-max normalization method
normalized_markers_meanexp <- apply(markers_meanexp_df, 1, function(x) { (x - min(x)) / max(x) })
library("grid")
mat = t(normalized_markers_meanexp)
my_palette <- c(colorRampPalette(colors = c("#6483ec", "#f3bf5a"))(n = length(colors)))

graphics.off()
pdf(file = paste0(lnWD, "celltypeMarkersHeatVitro.pdf"), width = 7, height = 22)
print(ComplexHeatmap::Heatmap(mat, name = "Expression", cluster_rows = F, cluster_columns = F, col = my_palette, row_names_gp = grid::gpar(fontsize = 20), heatmap_legend_param = list(legend_gp = gpar(fontsize = 50))))
dev.off()

save(normalized_markers_meanexp, file = paste0(lnWD, "heatmap_input.rda"))

hspc.combined$celltype <- "Delete"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "0")] <- "UN"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "1")] <- "MPP"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "2")] <- "ST-HSC"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "3")] <- "MEP/MKP"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "4")] <- "LT-HSC"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "5")] <- "MPP"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "6")] <- "MkP"
hspc.combined$celltype[which(hspc.combined$seurat_clusters == "7")] <- "Mast-cells"


saveRDS(hspc.combined, file = paste0(lnWD, "/annotatedCombinedVitro.rds"))

cat("### Done \n")