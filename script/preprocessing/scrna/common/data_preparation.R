# R function to prepare single cell data for monocle & Seurat

getCellCyclePhases <- function(gbm_cds,outdir = "./") {
  
  dir.create(path = outdir,recursive = T,showWarnings = F)
  gene_count_matrix <- as.matrix(exprs(gbm_cds))
  set.seed(100)
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  assignments <- cyclone(gene_count_matrix, mm.pairs, gene.names=rownames(gene_count_matrix))

  ## add cell cycle phase to pData
  pData(gbm_cds)$phases <- assignments$phases
  pData(gbm_cds)$G1_score <- assignments$scores$G1
  pData(gbm_cds)$G2M_score <- assignments$scores$G2M
  pData(gbm_cds)$S_score <- assignments$scores$S
  pData(gbm_cds)[which(pData(gbm_cds)$phases=="G1"),"phases"] <- "G1_G0"
  pData(gbm_cds)[which(pData(gbm_cds)$phases=="G2M"),"phases"] <- "G2_M"
  
  
  png(paste(outdir,"/cell_cycle_phases_assignment_on_umi.png",sep = ""))
  plot(assignments$score$G1, assignments$score$G2M, 
       xlab="G1 score", ylab="G2/M score", pch=16)
  dev.off()
  
  print("Cell cycle phases added..")
  return(gbm_cds)
}

addPercentMitoch <- function(gbm_cds) {      #in fact this is fraction of mitochondrial transcripts
  mitoGenes <- grep(pattern = "mt-",x=fData(gbm_cds)$gene_short_name,value = T)
  mitoGenesEnsembl <- rownames(gbm_cds)[which(is.element(el = fData(gbm_cds)$gene_short_name,set = mitoGenes))]
  percentMito <- Matrix::colSums(as.matrix(exprs(gbm_cds))[mitoGenesEnsembl, ])/Matrix::colSums(as.matrix(exprs(gbm_cds)))
  pData(gbm_cds)$percentMito <- percentMito
  return(gbm_cds)
}


#Filtering low-quality cells

filterCells <- function(gbm_cds,outdir="./",num_cells_expressed=10,min_expr=0.1) {
  
  dir.create(path = outdir,recursive = T,showWarnings = F)
  
  #Set detection threshold
  gbm_cds <- detectGenes(gbm_cds, min_expr = min_expr)
  
  #Define gene expressed as gene detected in more than n cell (only use for plotting distribution)
  expressed_genes <- row.names(subset(fData(gbm_cds),            
                                    num_cells_expressed >= num_cells_expressed))

  #filtering cells as explained in monocle doc, cutting distribution tails
  #eg. remove cells with no RNA or too much RNA

  pData(gbm_cds)$Total_mRNAs <- Matrix::colSums(exprs(gbm_cds))
  
  png(paste(outdir,"/densityTotalmRNA_raw.png",sep =""))
  qplot(Total_mRNAs, data = pData(gbm_cds), geom ="density")
  dev.off()
  
  upper_bound <- 10^(mean(log10(pData(gbm_cds)$Total_mRNAs)) +
                     2*sd(log10(pData(gbm_cds)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(gbm_cds)$Total_mRNAs)) -
                     2*sd(log10(pData(gbm_cds)$Total_mRNAs)))

  png(paste(outdir,"/densityWithFirstFilter.png",sep =""))
  print(qplot(Total_mRNAs, data = pData(gbm_cds), geom ="density") + 
          geom_vline(xintercept = lower_bound) + 
          geom_vline(xintercept = upper_bound))
  dev.off()
  
  gbm_cds <- gbm_cds[,pData(gbm_cds)$Total_mRNAs > lower_bound &
                       pData(gbm_cds)$Total_mRNAs < upper_bound]

  # redefining expressed gene after cells filtering
  gbm_cds <- detectGenes(gbm_cds, min_expr = min_expr)
  expressed_genes <- row.names(subset(fData(gbm_cds),num_cells_expressed >= num_cells_expressed)) #genes expressed in at least 10 cells of the data set.
  
  # Log-transform each value in the expression matrix.
  L <- log(exprs(gbm_cds[expressed_genes,]))

  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  melted_dens_df <- reshape2::melt(Matrix::t(scale(Matrix::t(L))))

  # Plot the distribution of the standardized gene expression values.
#  png(paste(outdir,"/standardized_distribution.png",sep =""))
#  print(qplot(value, geom = "density", data = melted_dens_df) +
#        stat_function(fun = dnorm, size = 0.5, color = 'red') +
#        xlab("Standardized log(UMIcounts)") +
#        ylab("Density"))
#  dev.off()
  
  ## Visualize percentage of mitochondiral genes
  gbm_cds <- addPercentMitoch(gbm_cds)
  
  png(paste(outdir,"/percent_mito.png",sep =""))
  plot(hist(gbm_cds$percentMito*100,breaks = 100),main = "percentage of mitochondrial transcripts")
  dev.off()
  
  return(gbm_cds)

}


