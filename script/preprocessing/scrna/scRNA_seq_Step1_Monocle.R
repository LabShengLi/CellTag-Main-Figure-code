#!/usr/bin/env Rscript
rm(list = ls())
library(cli)
library(Matrix)
library(lubripack)
library("monocle")
library("scran")
library("ggpubr")
library(doParallel)
library(foreach)

# register cores (use all available or set manually)
registerDoParallel(cores = parallel::detectCores())

# may not needed

### Read in Raw data from 10X folders
wdir <- 'projects/U01_aim2/Cross_Expirement/PreProcessing'
cellRangerDir <- "projects/U01_aim2/Cross_Expirement/CellRanger/Out/"

sourceDir <- "scripts/common/"
source(paste0(sourceDir, "funForLoading_PK.R"))
source(paste0(sourceDir, "data_preparation.R"))

dir.create(wdir, showWarnings = FALSE, recursive = T)


sample_list = c("OO", "O_vitro", "OY", "Y_vitro", "YY", "YO")

#sample="YY"
foreach(sample = sample_list) %dopar% {
	# for (sample in c("OO", "O_vitro", "OY", "Y_vitro", "YY", "YO")) {

	print(paste0("Processing age: ", sample))

	setwd(wdir)

	outdir <- paste0(wdir, "/", sample); outdir
	dir.create(outdir, showWarnings = FALSE, recursive = T)

	opt = list('inputMatrixDir' = paste0(cellRangerDir, sample, "/outs/filtered_feature_bc_matrix"), #,
			   'outfile' = paste0(sample, "_10X_data_newCellDataSet.rds"),
			   'cellranger' = "3",
			   'subsample' = NULL,
			   "sampleName" = sample) #)
	print(opt)

	## -----------------------------------------------------------------------------
	## Processing data
	## -----------------------------------------------------------------------------

	if (opt$cellranger == "3") {
		crmat <- loadCellRangerMatrix_cellranger3(opt$inputMatrixDir, sample_name = opt$sampleName)
	} else {
		crmat <- loadCellRangerMatrix(opt$inputMatrixDir, sample_name = opt$sampleName)
	}

	fd <- crmat$fd
	pd <- crmat$pd

	# Column 'symbol' is the one (from cellRanger workflow) that corresponds to featureData's gene short names.

	colnames(fd)[which(colnames(fd) == "symbol")] <- "gene_short_name"

	cds <- newCellDataSet(crmat$exprs,
						  phenoData = new("AnnotatedDataFrame", data = pd),
						  featureData = new("AnnotatedDataFrame", data = fd),
						  expressionFamily = negbinomial.size())

	## Subsample
	if (is.null(opt$subsample) == F) {
		opt$subsample <- as.numeric(opt$subsample)
		print(opt$subsample)

		cellSubset <- sample(rownames(pData(cds)), size = opt$subsample * length(rownames(pData(cds))))
		cds <- cds[, cellSubset]
	}
	saveRDS(cds, file = paste0(outdir, "/", opt$outfile))


	# print some stuff
	# head of the expression matrix
	print(exprs(cds[c(1:5), c(1:5)]))

	# cells information
	print(head(pData(cds)))

	# gene information
	print(head(Biobase::fData(cds)))

	# First getting cell cycle phases with cyclone function from package scrann (see data_preparation.R))
	cellCycleDir <- paste0(outdir, "/cellCycle")

	########## ------>>>>>>>>>>>> TIME CONSUMING STEP ***************
	system(paste0("mkdir ", cellCycleDir))

	cds <- getCellCyclePhases(cds, outdir = cellCycleDir)

	# Filtering cells on UMI counts as explained in monocle doc
	# print num cells before filtering

	print("Matrix dim before cell filtering :")
	dim(pData(cds))

	########## ------>>>>>>>>>>>> TIME CONSUMING STEP ***************
	cds <- filterCells(cds, paste0(cellCycleDir, "/cellsFiltering"))
	# Filter cells with to few expressed genes to assign a cell cycle phases with cyclone
	# not the the case with cellranger2 workflow

	valid_cells <- row.names(subset(pData(cds),
									is.na(G2M_score) == FALSE &
										is.na(S_score) == FALSE &
										is.na(G1_score) == FALSE &
										is.na(phases) == FALSE))

	nonValidCellNum <- length(row.names(pData(cds))) - length(valid_cells)
	print(paste(nonValidCellNum, "cells were filtered out because they express to few genes to assign a cell cycle phase with cyclone to them. "))

	df <- t(data.frame(SequenceCells = length(valid_cells) + nonValidCellNum, ValidCells = length(valid_cells), DiscardedCells = nonValidCellNum))
	tab <- ggpubr::ggtexttable(df, theme = ttheme("mViolet"))

	#cells are row in pData(monocle) but column in monocle object
	cds <- cds[, valid_cells]

	print("Matrix dim after cell filtering :")
	print(dim(pData(cds)))

	saveRDS(object = cds, file = paste0(cellCycleDir, "/", sample, "_CellCyclePhase_Processed.rds"))
}
