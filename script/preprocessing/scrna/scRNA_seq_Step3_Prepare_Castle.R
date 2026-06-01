#!/usr/bin/env Rscript
rm(list = ls())
library(Seurat)
library(Biobase)
library("lubripack")
lubripack(c("monocle", "scran", "Seurat"))
lubripack(c("DoubletFinder", "stringr", "useful", "ggplot2", "gridExtra", "grid", "gtable", "data.table", "cowplot", "plyr", "reshape2", "tidyr", "pheatmap", "RColorBrewer", "tidyverse", "colorspace", "splines", "AnnotationDbi", "dplyr", "gridBase", "ggrepel", "calibrate", "ggplotify"))
lubripack(c("RColorBrewer", "ggpubr", "harmony", "SingleCellExperiment", "SummarizedExperiment", "ComplexHeatmap"))
library("Seurat")
library("monocle")


# may not needed, create supervised cell annotation data, run only once to build the objects

opaque <- function(color, opacity, max = 255) {
	rgb(t(col2rgb(color)), alpha = opacity, max = max)
}

wdir <- 'projects/U01_aim2/Cross_Expirement/PreProcessing'

inputData <- "from_jax/Lamis/U01_Projects/U01/Chenx/CloneTracing_project/data/publicData/Rodriguez/rawCount/"

source("scripts/common/funForLoading_PK.R")
source("scripts/common/data_preparation.R")

opt = list('inputData' = inputData,
		   'outdir' = paste0(wdir, "/RodriguezPreparation"),
		   'minPropCellExp' = 0.001)

# if help was asked, print a friendly message
# and exit with a non-zero error code

if (is.null(opt$outdir)) {
	opt$outdir = "./"
}

dir.create(opt$outdir, recursive = T, showWarnings = F)


if (is.null(opt$minPropCellExp)) {
	opt$minPropCellExp <- 0.001
	print(opt$minPropCellExp)
}

#setwd(opt$outdir)

#get Seurat3 Rodriguez
rawCountsPaths <- list.files(opt$inputData, pattern = "raw_umifm_counts.csv", full.names = T)

#make sure we gate the correct file (no MPP4)
rawCountsPaths <- rawCountsPaths[which(grepl(pattern = 'LTHSC|MPP2|MPP3|MPP4|STHSC', rawCountsPaths))]


getSeuratRodriguez <- function(rawCountsPath) {

	rawTable <- read.csv(rawCountsPath)

	phenoData <- rawTable[, c(1:5)]
	rownames(phenoData) <- rawTable$cell_id
	data <- rawTable[, c(6:ncol(rawTable))]
	rownames(data) <- rawTable$cell_id
	seurat <- CreateSeuratObject(counts = t(data), meta.data = phenoData)

	#Use of Rodriguez QC to filter poor qualisty cells
	seurat <- subset(x = seurat, subset = pass_filter > 0.1)
	print(unique(seurat@meta.data$seq_run_id))
	return(seurat)

}


seuratObjects <- lapply(rawCountsPaths, getSeuratRodriguez)

# We do not take the MPP4, BUT I TAKE
seuratAll <- merge(x = seuratObjects[[1]],
				   y = c(seuratObjects[[2]], seuratObjects[[3]],
						 seuratObjects[[4]], seuratObjects[[5]]))

seuratAll@meta.data$library_id <- factor(seuratAll@meta.data$library_id, levels = c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4"))


## Switch to Monocle object to more easily discard unexpressed genes

pd <- new("AnnotatedDataFrame", data = seuratAll@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(seuratAll)))
rownames(fd) <- fd$gene_short_name


object <- JoinLayers(object = seuratAll, layers = NULL)
monocleAll <- newCellDataSet(GetAssayData(object = object, slot = "counts"),
							 phenoData = pd,
							 featureData = fd,
							 lowerDetectionLimit = 0.1,
							 expressionFamily = negbinomial.size())


pData(monocleAll)$Total_mRNAs <- Matrix::colSums(exprs(monocleAll))
monocleAll <- detectGenes(monocleAll, min_expr = 0.1)

# Filter out genes expressed in too few cells
fData(monocleAll)$use_for_castle <- fData(monocleAll)$num_cells_expressed > opt$minPropCellExp * ncol(monocleAll)
# 
gbm_to_seurat <- monocleAll[fData(monocleAll)$use_for_castle == T,]

seurat <- CreateSeuratObject(counts = exprs(gbm_to_seurat), meta.data = pData(gbm_to_seurat))

# Need normalized data for castle
seurat <- NormalizeData(object = seurat)

saveRDS(object = seurat, file = paste0(opt$outdir, "/seuratForCastle.rds"))
