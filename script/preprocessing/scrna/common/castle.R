castle <- function(source_seurat,target_seurat,outdir = "./") {
  set.seed(2018)
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  dir.create(outdir,showWarnings = F)

  # 1. Load datasets in scater format: loaded files expected to contain "Large SingleCellExperiment" object
  
  
  # 1.1 convert from seurat (be careful seurat 2) to sce object if seurat object as input
  source <- Convert(from = source_seurat, to = "sce")
  target <- Convert(from = target_seurat, to = "sce")


  
  
  ds1 = t(exprs(source)) 
  ds2 = t(exprs(target)) 
  sourceCellTypes = as.factor(colData(source)[,"library_id"])

  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(exprs(source), 1, function(x) { sum(x > 0) } )
  print("source OK")
  target_n_cells_counts = apply(exprs(target), 1, function(x) { sum(x > 0) } )
  common_genes = intersect( rownames(source)[source_n_cells_counts>10], 
                            rownames(target)[target_n_cells_counts>10]
  )
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  # 4. Highest mutual information in source
  topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),sourceCellTypes,method = "nmi") }), decreasing = T))
  
  # 5. Top n genes that appear in both mi and avg
  selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
  
  # 6. remove correlated features
  tmp = cor(as.matrix(ds[,selectedFeatures]), method = "pearson")
  tmp[!lower.tri(tmp)] = 0
  selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
  remove(tmp)
  
  # 7,8. Convert data from continous to binned dummy vars
  # break datasets to bins
  dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
  # use only bins with more than one value
  nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
  # convert to dummy vars
  ds = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
  remove(dsBins, nUniq)
  
  
  # 9. Classify
  train = runif(nrow(ds[isSource,]))<0.8
  #slightly different setup for multiclass and binary classification
  if (length(unique(sourceCellTypes)) > 2) {
    xg=xgboost(data=ds[isSource,][train, ] ,
               label=as.numeric(sourceCellTypes[train])-1,
               objective="multi:softmax", num_class=length(unique(sourceCellTypes)),
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  } else {
    xg=xgboost(data=ds[isSource,][train, ] ,
               label=as.numeric(sourceCellTypes[train])-1,
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  }
  
  
  
  # 10. Predict
  predictedClasses = predict(xg, ds[!isSource, ])
  
  target_seurat@meta.data$predicted <- predictedClasses
 
  
  conv <- data.frame(cellType = unique(sourceCellTypes),numCellTypes=unique(as.numeric(sourceCellTypes[train]))-1)
  
  
  for (num in conv$numCellTypes) {
    target_seurat@meta.data$predicted[target_seurat@meta.data$predicted==num] <- as.character(conv[which(conv$numCellType== num),"cellType"])
  }    
  
  
  
  target_seurat@meta.data$predicted <- factor(target_seurat@meta.data$predicted,levels = c("LTHSC","STHSC", "MPP2", "MPP3","MPP4"))
  png(paste(outdir,"/cellTypeLearned_tsne.png", sep = ""),height = 800,width = 800)
  TSNEPlot(target_seurat,group.by="predicted")
  dev.off()
  
  
  
  return(target_seurat)
}

castleSeurat3 <- function(source_seurat,target_seurat,outdir = "./") {
  set.seed(2018)
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  dir.create(outdir,showWarnings = F)
  
  # 1. Load datasets in scater format: loaded files expected to contain "Large SingleCellExperiment" object
  
  
  # 1.1 convert from seurat (be careful seurat 3) to sce object if seurat object as input
  source <- as.SingleCellExperiment(source_seurat, to = "sce",data = "logcounts")
  target <- as.SingleCellExperiment(target_seurat, to = "sce",data = "logcounts")
  
  
  
  
  ds1 = t(exprs(source)) 
  ds2 = t(exprs(target)) 
  sourceCellTypes = as.factor(colData(source)[,"library_id"])
  
  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(exprs(source), 1, function(x) { sum(x > 0) } )
  print("source OK")
  target_n_cells_counts = apply(exprs(target), 1, function(x) { sum(x > 0) } )
  common_genes = intersect( rownames(source)[source_n_cells_counts>10], 
                            rownames(target)[target_n_cells_counts>10]
  )
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  # 4. Highest mutual information in source
  topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),sourceCellTypes,method = "nmi") }), decreasing = T))
  
  # 5. Top n genes that appear in both mi and avg
  selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
  
  # 6. remove correlated features
  tmp = cor(as.matrix(ds[,selectedFeatures]), method = "pearson")
  tmp[!lower.tri(tmp)] = 0
  selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
  remove(tmp)
  
  # 7,8. Convert data from continous to binned dummy vars
  # break datasets to bins
  dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
  # use only bins with more than one value
  nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
  # convert to dummy vars
  ds = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
  remove(dsBins, nUniq)
  
  
  # 9. Classify
  train = runif(nrow(ds[isSource,]))<0.8
  #slightly different setup for multiclass and binary classification
  if (length(unique(sourceCellTypes)) > 2) {
    xg=xgboost(data=ds[isSource,][train, ] ,
               label=as.numeric(sourceCellTypes[train])-1,
               objective="multi:softmax", num_class=length(unique(sourceCellTypes)),
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  } else {
    xg=xgboost(data=ds[isSource,][train, ] ,
               label=as.numeric(sourceCellTypes[train])-1,
               eta=0.7 , nthread=5, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
  }
  
  
  
  # 10. Predict
  predictedClasses = predict(xg, ds[!isSource, ])
  
  target_seurat@meta.data$predicted <- predictedClasses
 
  
  conv <- data.frame(cellType = unique(sourceCellTypes),numCellTypes=unique(as.numeric(sourceCellTypes[train]))-1)
  
  
  for (num in conv$numCellTypes) {
    target_seurat@meta.data$predicted[target_seurat@meta.data$predicted==num] <- as.character(conv[which(conv$numCellType== num),"cellType"])
  }    
  
  
  
  target_seurat@meta.data$predicted <- factor(target_seurat@meta.data$predicted,levels = c("LTHSC","STHSC", "MPP2", "MPP3","MPP4"))
  png(paste(outdir,"/cellTypeLearned_tsne.png", sep = ""),height = 800,width = 800)
  DimPlot(target_seurat,group.by ="predicted")
  dev.off()
  
  
  
  return(target_seurat)
}


transfertToMonocle <- function(target_seurat,target_monocle,outdir = "./") {
  pData(target_monocle)$predicted <- target_seurat@meta.data$predicted
  png(paste(outdir,"/cellTypeLearned_ddrtree.png", sep = ""),height = 800,width = 800)
  print(plot_cell_trajectory(target_monocle,color_by = "predicted"))
  dev.off()
  return(target_monocle)
}

