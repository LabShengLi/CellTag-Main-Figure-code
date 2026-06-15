#!/usr/bin/env Rscript
rm(list = ls())
library("plyr")
library("ggplotify")
library("stringr")
library(gridExtra)
library(grid)
library("Biobase")
library("Seurat")
library("DoubletFinder")
library(doParallel)
library(foreach)

# create seurat object for each library

# register cores (use all available or set manually)
registerDoParallel(cores = parallel::detectCores())

wdir <- 'PreProcessing'

min.genes = 500
min.cells = 3
mtpercent = 20
rbpercent = 50
FeatureNum = 2000
RUNDoubletAgain = "YES"

setwd(wdir)

################## Convert Monocle processed data to Seurat Object
RUNDoubletAgain == "YES"
# RUNDoubletAgain=="NO"
cellRangerDir <- "Out"

sourceDir <- "common/"

source(paste0(sourceDir, "funForLoading_PK.R"))
source(paste0(sourceDir, "data_preparation.R"))
source(paste0(sourceDir, "Useful_Functions.R"))

sample_list <- c("OO", "O_vitro", "OY", "Y_vitro", "YY", "YO")
sample <- "OO"
foreach(sample = sample_list) %dopar% {
	# for (sample in sample_list) {
	print(paste0("Processing age: ", sample))

	setwd(wdir)

	outdir <- paste0(wdir, "/", sample); outdir
	dir.create(outdir, showWarnings = FALSE, recursive = T)
	cellCycleDir <- paste0(outdir, "/cellCycle")

	opt = list(
		'inputRDS' = paste0(cellCycleDir, "/", sample, "_CellCyclePhase_Processed.rds"),
		'outdir' = paste0(wdir, "/", sample, "/Seurat3"),   #,
		"num_dim" = "30",
		"correction" = "G2M_score+S_score+G1_score", #"Covariable to use as blocking factor (eg one or several columns of pData: betch, cell cycle phases... separated by +)",
		"minPropCellExp" = 0.001, #"numeric", "minimal proportion of cells that expressed the genes kept for the analaysis 0.001 by default",
		"norm_method" = "logNorm", #"character", "normalisation method, logNorm seurat (by default) or sctransform",
		"resolution" = 0.2, #"numeric", "resolution for Seurat clustering 0.9 by default",
		"identRemoved" = NULL, #"character", "Optionnal cluster to remove only work if input is a seurat object",
		"nonExpressedGenesRemoved" = FALSE, #"logical", "non expressed gene already removed default to FALSE",
		"gprofiler" = TRUE, # "logical", "if true doing gprofiler default to true",
		"logfc_threshold" = 0.25) #"numeric", "logfc threshold for finding cluster markers (1 cluster vs all deg) 0.25 by default",
	print(opt)

	dir.create(opt$outdir, recursive = TRUE); opt$outdir

	print(opt$gprofiler)
	print(opt$nonExpressedGenesRemoved)


	# get correction vector
	corrections <- strsplit(x = opt$correction, split = "\\+")[[1]]

	# create blank for the plots
	library("grid")
	blank <- grid.rect(gp = gpar(col = "white"))

	# input is a monocle object in the workflow
	# Code to use this script with seurat input (eg to re analyse)
	print("input monocle")
	gbm_cds <- readRDS(opt$inputRDS)

	if (opt$nonExpressedGenesRemoved == F) {
		print("remove non expressed genes (non expressed in at least X% of the cells X user option in monocle dp feature 5% in seurat tutorial 0,1%)")
		fData(gbm_cds)$use_for_seurat <- fData(gbm_cds)$num_cells_expressed > opt$minPropCellExp * ncol(gbm_cds)

		gbm_to_seurat <- gbm_cds[fData(gbm_cds)$use_for_seurat == T,]
	} else {
		gbm_to_seurat <- gbm_cds
	}

	if (is.element("Cluster", colnames(pData(gbm_to_seurat)))) {
		colnames(pData(gbm_to_seurat))[which(colnames(pData(gbm_to_seurat)) == "Cluster")] <- "Cluster_monocle"
	}

	# Only needed if ensemble id (it is the case in the workflow)

	dupGeneNames <- fData(gbm_to_seurat)[which((duplicated(featureData(gbm_to_seurat)$gene_short_name))), "gene_short_name"]

	if (length(dupGeneNames) == 0) {
		rownames(gbm_to_seurat) <- fData(gbm_to_seurat)$gene_short_name

	} else {
		write.csv(fData(gbm_to_seurat)[which(is.element(fData(gbm_to_seurat)$gene_short_name, dupGeneNames)),], paste(opt$outdir, "/dupGenesName.csv", sep = ""))
		print("Dup gene short names existing, making them unique...")
		rownames(gbm_to_seurat) <- make.unique(fData(gbm_to_seurat)$gene_short_name, sep = "--")

	}

	## seurat object. contains full data


	SequencedCells <- nrow(pData(gbm_to_seurat))
	seurat.FullData <- CreateSeuratObject(counts = exprs(gbm_to_seurat), meta.data = pData(gbm_to_seurat), min.cells = min.cells, min.features = min.genes)
	BasicQC <- SequencedCells - nrow(seurat.FullData@meta.data); BasicQC; BasicQC / SequencedCells

	seurat.FullData[["percent.mt"]] <- PercentageFeatureSet(seurat.FullData, pattern = "^mt-")
	seurat.FullData[["percent.rb"]] <- PercentageFeatureSet(seurat.FullData, pattern = "^Rp[sl]")

	#pdf(file = paste0(opt$outdir,"/Basic_QC_",sample,"_mincells_",min.cells,"_mingenes_",min.genes,"_mitoPerc_",mtpercent,"_riboPerc_",rbpercent,".pdf"), height = 6, width = 8)

	SeuratObject <- subset(seurat.FullData, subset = percent.mt < mtpercent & percent.rb < rbpercent)
	SeuratObject$Cells <- rownames(SeuratObject@meta.data)
	SeuratObject$orig.ident <- sample
	MitoRiboQC <- SequencedCells -
		BasicQC -
		nrow(SeuratObject@meta.data); MitoRiboQC; MitoRiboQC / SequencedCells
	QCpassCells <- nrow(SeuratObject@meta.data)

	cutoff.df <- t(data.frame(SequencedCells = SequencedCells, Min.Cells.Genes.Filter = BasicQC, Mt.Ribo.Filter = MitoRiboQC, QCpassCells = QCpassCells))
	print(cutoff.df)
	cutoff.df <- as.data.frame(cutoff.df)
	colnames(cutoff.df) <- "Cells"
	cutoff.df$Percentage <- round((cutoff.df$Cells / SequencedCells * 100), 1)
	cutoff.df[1, 2] <- ""

	# install.packages("ggpubr")
	library("ggpubr")
	# install.packages("ggplot2")
	ggtexttable(cutoff.df, theme = ttheme("mViolet"))
	print(VlnPlot(seurat.FullData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.1, ncol = 4, cols = "cyan4"))
	print(VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.1, ncol = 4, cols = "cyan4"))

	#dev.off()

	# install.packages("DescTools")
	library("DescTools")

	matrix_dir <- paste0(cellRangerDir, "/", sample); matrix_dir
	sample.Dir <- paste0(matrix_dir, "/outs/filtered_feature_bc_matrix")
	print(paste0("sample Dir is: ", sample.Dir))
	saveDIR <- opt$outdir; saveDIR

	head(SeuratObject@meta.data)
	saveRDS(SeuratObject, file = paste0(opt$outdir, "/SCdata_Initial.rds"))


	if (RUNDoubletAgain == "YES") {
		print(paste0("Already completed Doublet Detection, Skipping for Now"))

		Temp <- SeuratObject
		### Duplicate Detection
		print(paste0("Running Doublet Detection for ", sample))

		library("cowplot")
		library("tidyverse")
		DoubletSCdata <- Doublet_Detection_DF(sample.Dir, saveDIR, sample, FeatureUseCount = FeatureNum, PCAnum = 10, Species = "mmu", plotCCgene = TRUE,
											  mincells = min.cells, mingenes = min.genes, mtpercent = mtpercent)

		head(DoubletSCdata@meta.data)
		head(Temp@meta.data)
		Temp@meta.data <- merge(Temp@meta.data, DoubletSCdata@meta.data, by = 0, all.x = TRUE)  ## by rownames
		Temp@meta.data <- Temp@meta.data %>%
			dplyr::select(names(Temp@meta.data)[!names(Temp@meta.data) %like% "\\.y$"]) %>%
			rename_with(~str_remove(., "\\.x$")) %>%
			column_to_rownames(var = "Row.names")

		#head(DoubletSCdata@meta.data)
		#head(Temp@meta.data)

		#tail(DoubletSCdata@meta.data)
		#tail(Temp@meta.data)

		print(all(rownames(Temp@meta.data) == Temp$Cells))
		print(all(rownames(Temp@meta.data) == Temp$barcode))

		Idents(Temp) <- "DoubletFinder"
		SCdata <- subset(Temp, idents = c("Singlet")); SCdata

	} else {

		print(paste0("Reading the detected doublets:"))
		DDdir <- paste0(saveDIR, "/Doublet_Detection_", sample, "/DoubletFinder_", sample)
		setwd(DDdir)
		SingletsCells <- read.table(file = paste0(DDdir, "/Singlets_DoubletFinder_Calls_", sample, "_using_PCA_10_res_0.5.txt"), header = T); head(SingletsCells); nrow(SingletsCells)
		DoubletsCells <- read.table(file = paste0(DDdir, "/Doublets_DoubletFinder_Calls_", sample, "_using_PCA_10_res_0.5.txt"), header = T); head(DoubletsCells); nrow(DoubletsCells)

		SCdata <- subset(SeuratObject, cells = rownames(SingletsCells)); DoubletSCdata

		QCpassCells <- nrow(SCdata@meta.data)
		cutoff.df <- t(data.frame(SequencedCells = SequencedCells, Min.Cells.Genes.Filter = BasicQC, Mt.Ribo.Filter = MitoRiboQC, Doublets = nrow(DoubletsCells), QCpassCells = QCpassCells))
		print(cutoff.df)
		cutoff.df <- as.data.frame(cutoff.df)
		colnames(cutoff.df) <- "Cells"
		cutoff.df$Percentage <- round((cutoff.df$Cells / SequencedCells * 100), 1)
		cutoff.df[1, 2] <- ""

		setwd(opt$outdir)
		#pdf(file = paste0("Basic_QC_",sample,".pdf"), height = 6, width = 8)
		pdf(file = paste0("Basic_QC_", sample, "_mincells_", min.cells, "_mingenes_", min.genes, "_mitoPerc_", mtpercent, "_riboPerc_", rbpercent, "_Doublets_Removed_Sep.pdf"), height = 6, width = 8)

		PlotTable.withRowNames.New(df = cutoff.df, MainTitle = sample)
		print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.1, ncol = 4, cols = "cyan4"))
		dev.off()
	}


	print(SCdata)
	head(SCdata@meta.data)

	saveRDS(SCdata, file = paste0(opt$outdir, "/SCdata_Doublet_Removed.rds"))


	if (opt$norm_method == "sctransform") {
		paste0("Performing sctransform")
		SCdata <- SCTransform(object = SCdata, vars.to.regress = c("G2M_score", "S_score", "G1_score"))

	} else {

		paste0("Performing logNormalization")
		SCdata <- NormalizeData(object = SCdata)
		SCdata <- FindVariableFeatures(object = SCdata, selection.method = "vst", nfeatures = FeatureNum, verbose = T)
		SCdata <- ScaleData(object = SCdata, vars.to.regress = corrections)

	}

	SCdata <- RunPCA(object = SCdata)

	#png(paste0(opt$outdir,"/ElbowPlot.png"))
	ElbowPlot(object = SCdata, ndims = 30)
	#dev.off()


	print("Clustering...")

	SCdata <- FindNeighbors(object = SCdata, dims = c(1:opt$num_dim), k.param = 20)
	#SCdata <- FindClusters(object = SCdata,resolution = c(0.1,0.2,0.25,0.3,0.4,0.5,0.6))
	SCdata <- FindClusters(object = SCdata, resolution = c(0.1, 0.2))

	#print("Running TSNE...")
	#SCdata <- RunTSNE(SCdata,dims = c(1:opt$num_dim))

	print("Running UMAP...")
	SCdata <- RunUMAP(SCdata, dims = c(1:opt$num_dim))

	print("UMAP ok")
	print(colnames(SCdata@meta.data))

	colPrefix <- "RNA_snn_res."
	if (opt$norm_method == "sctransform") {
		colPrefix <- "SCT_snn_res."
	}

	umapListRes <- list()
	#for (r in c(0.6,0.8,1,1.2)) {
	for (r in c(0.1, 0.2)) {
		umapListRes[[as.character(r)]] <- DimPlot(SCdata,
												  reduction = "umap",
												  label = T,
												  group.by = paste(colPrefix, r, sep = "")) +
			NoLegend() +
			ggtitle(paste("res", r))

	}

	#png(paste(opt$outdir,"/umap_different_res.png",sep = ""),height = 800,width = 800)
	plot_grid(plotlist = umapListRes, nrow = 2)
	#dev.off()


	Idents(SCdata) <- paste(colPrefix, opt$resolution, sep = "")
	SCdata@meta.data$numclust <- SCdata@meta.data[, paste(colPrefix, opt$resolution, sep = "")]


	#Check for unwanted source of variation
	pPhases <- DimPlot(SCdata, group.by = "phases")
	if (!is.null(SCdata@meta.data$predicted)) {
		pPred <- DimPlot(SCdata, group.by = "predicted")
	} else {
		pPred <- blank
	}

	pUMI <- FeaturePlot(SCdata, "Total_mRNAs")
	pMito <- FeaturePlot(SCdata, "percentMito")

	#png(paste(opt$outdir,"/umap_factors.png",sep = ""),height = 800,width = 800)
	library(grid)
	grid.arrange(pPhases, pPred, pUMI, pMito)
	#dev.off()

	gene_list <- c("Pdzk1ip1", "Mllt3")

	gene_list <- gene_list[which(is.element(set = rownames(SCdata), el = gene_list))]

	dir.create(paste(opt$outdir, "/genes/umap/", sep = ""), recursive = T)
	for (g in gene_list) {
		png(paste(opt$outdir, "/genes/umap/", g, ".png", sep = ""))
		plot(FeaturePlot(SCdata, features = g))
		dev.off()
	}

	markers <- FindAllMarkers(SCdata, only.pos = T, logfc.threshold = opt$logfc_threshold)
	markers <- markers[which(markers$p_val_adj < 0.05),]

	write.table(x = markers, paste(opt$outdir, "/markers.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

	dir.create(paste(opt$outdir, "/markers/", sep = ""))

	ylab <- "LogNormalized UMI counts"
	if (opt$norm_method == "sctransform") {
		ylab <- "Expression level"
	}

	for (numClust in unique(markers$cluster)) {
		print(head(markers[which(markers$cluster == numClust),], n = 9))
		png(paste(opt$outdir, "/markers/Cluster_", numClust, "_topGenesVlnPlot.png", sep = ""), width = 1000, height = 1000)
		plot(VlnPlot(object = SCdata, features = head(markers[which(markers$cluster == numClust), "gene"], n = 9), pt.size = 0.5) +
				 labs(x = "Clusters", y = ylab, colour = "black") +
				 theme(axis.text = element_text(size = 20),
					   plot.title = element_text(size = 25)))
		dev.off()
	}

	saveRDS(SCdata, file = paste0(opt$outdir, "/SCdata.rds"))
	getwd()

	oo = readRDS(paste0(opt$outdir, "/SCdata.rds"))
	print(oo)
}



