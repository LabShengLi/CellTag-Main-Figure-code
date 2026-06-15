
library("lubripack")
lubripack(c("RColorBrewer", "ggsci"))
    mypal_npg = pal_npg("nrc", alpha = 0.7)(10)
    mypal_npg
    mypal_nejm = pal_nejm("default", alpha = 0.7)(8)
    mypal_nejm
    mypal_jco = pal_jco("default", alpha = 0.7)(10)
    mypal_jco
    mypal_ucscgb = pal_ucscgb("default", alpha = 0.7)(26)
    mypal_ucscgb
    mypal_startrek = pal_startrek("uniform", alpha = 0.7)(7)
    mypal_startrek
    mypal_d3 = pal_d3("category20")(20)
    mypal_d3
    ClusPallette <- c(mypal_ucscgb, mypal_npg, mypal_jco, mypal_nejm)
    ClusPallette
    PairedPalette12 <- brewer.pal(n = 12, name = "Paired")
    OtherPalette <- c(mypal_jco[1:2], mypal_jco[4:5], pal_d3("category20")(20)[2:19],
        pal_lancet("lanonc")(9))
    OtherPalette <- OtherPalette[!OtherPalette %in% "#FF9896FF"]
    #show_col(mypal)

ClusPalette <- ClusPallette

#ClusPallette <- c("#0000ee", "#ff0000", "#00ff00", "#E69F00", "#ffff00", "#bf3eff", "#00ffff", "#FF1493", "#ff8247", "#836fff", "#999999", 
#                  "#7FFFD4", "#b5d4ff", "#ff80bf", "#CAFF70",
#                  "#27408B", "#8b0000", "#008000", "#DEB887", "#8B008B", "#F0E442", "#56B4E9", "#CD5C5C", "#000000", "#D3D3D3"
#)
Dark.Pallette <- c("#27408B", "#8b0000", "#008000", "#D2691E", "#F0E442",
                   "#8B008B", "#7ac5cd", "#D2B48C","#66cdaa", "#8A2BE2",
                   "#DB7093", "#9ACD32", "#cdad00", "#eecbad", "#d8bfd8",
                   "#787878", "#000000")

Light.Pallette <- c("#0000ee", "#ff0000", "#00ff00", "#ff8247", "#ffff00",
                    "#bf3eff", "#98f5ff", "#E5D3B3", "#7fffd4","#836fff",
                    "#ff3e96", "#CAFF70", "#ffb90f", "#ffdab9", "#ffe1ff",
                    "#737373")


Doublet_Detection_DF <- function(matrix.DIR, saveDIR, Sample, Species="hsa", FeatureUseCount=2500, PCAnum=10, resClus = 0.5, mincells=3, mingenes=500, mtpercent=20, rbpercent=50, ClusPallette=ClusPallette, DoubletFinder = TRUE, plotCCgene = TRUE){
  
  #matrix.DIR=matrix.DIR
  #saveDIR=saveDIR
  #Sample=Sample
  #Species="hsa"
  #FeatureUseCount=2500
  #PCAnum=10
  #resClus = 0.5
  #mincells=3
  #mingenes=500
  #mtpercent=20
  #rbpercent=50
  #ClusPallette=ClusPallette
  #DoubletFinder = TRUE
  #plotCCgene = TRUE
  
  #library(DoubletDecon)
  library(DoubletFinder)
  library(ggpubr)

  ClusOrder <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60")
  
  setwd(saveDIR)
  DDdir <- paste(getwd(),paste0("Doublet_Detection_",Sample),sep="/"); print(DDdir)
  dir.create(file.path(getwd(),paste0("Doublet_Detection_",Sample)), showWarnings = FALSE)
 
  print("**** STEP1")
  
  print(paste0("ClusOrder is:"))
  print(paste0(ClusOrder))
 
  print("Creating Seurat Object")
  print(paste0("min.cells:",mincells, ", min.features:",mingenes))
  data <- Read10X(data.dir = matrix.DIR)
  BeforeFilter=ncol(data); BeforeFilter
  SequencedCells=BeforeFilter
  SeuratObject <- CreateSeuratObject(counts = data, project = Sample, min.cells = mincells, min.features = mingenes)
  SeuratObject@meta.data$Cells <- rownames(SeuratObject@meta.data)
  SeuratObject@meta.data$Library <- Sample
  head(SeuratObject@meta.data)
  
  if(Species=="hsa"){
    print("Counting MT % for Human")
    SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")
    SeuratObject[["percent.rb"]] <- PercentageFeatureSet(SeuratObject, pattern = "^RP[SL]")
    
  } else {
    print("Counting MT % for Mouse")
    SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^mt-")
    SeuratObject[["percent.rb"]] <- PercentageFeatureSet(SeuratObject, pattern = "^Rp[sl]")
  }
  
  print("**** STEP2")
  
  print("Filtering Data based on specified filters")
  Min.Cells.Genes.Filter <- SequencedCells-nrow(SeuratObject@meta.data)
  t1 <- VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.1, ncol = 4, cols = "cyan4") #+ ggtitle(paste0("rb%  (",nrow(SeuratObject@meta.data)," cells)\nfiltered ",Min.Cells.Genes.Filter, " cells"))
  SeuratObject <- subset(SeuratObject, subset = percent.mt < mtpercent & percent.rb < rbpercent)
  Mt.Ribo.Filter=SequencedCells-nrow(SeuratObject@meta.data)-Min.Cells.Genes.Filter
  t2 <- VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.1, ncol = 4, cols = "cyan4") #+ ggtitle(paste0("rb%\n  (",nrow(SeuratObject@meta.data)," cells)\nfiltered ",Mt.Ribo.Filter, " cells"))
  QCpassCells=nrow(SeuratObject@meta.data)
  p1 <- plot_grid(t1, t2, ncol = 1)
  
  print(paste0(Sample, " cells BEFORE filtering:",BeforeFilter))
  print(paste0("Cells Discarded during min.cell and min.genes QC filtering: ",Min.Cells.Genes.Filter, " cells, ",round((Min.Cells.Genes.Filter/BeforeFilter*100),1),"%"))
  print(paste0("Cells Discarded during mito and Ribo QC filtering: ",Mt.Ribo.Filter, " cells, ",round((Mt.Ribo.Filter/BeforeFilter*100),1),"%"))
  Discardedcells=Min.Cells.Genes.Filter+Mt.Ribo.Filter
  print(paste0("Total Cells Discarded during min.cell,  min.genes, mito and Ribo QC filtering: ",Discardedcells, " cells, ",round((Discardedcells/BeforeFilter*100),1),"%"))
  print(paste0("Cells After QC filtering: ",QCpassCells, " cells, ",round((QCpassCells/BeforeFilter*100),1),"%"))
  
  cutoff.df <- t(data.frame(SequencedCells = BeforeFilter, Min.Cells.Genes.Filter=Min.Cells.Genes.Filter, Mt.Ribo.Filter=Mt.Ribo.Filter,QCpassCells=QCpassCells)); print(cutoff.df)
  cutoff.df <- as.data.frame(cutoff.df)
  colnames(cutoff.df) <- "Cells"
  cutoff.df$Percentage <- round((cutoff.df$Cells / BeforeFilter*100),1)
  cutoff.df[1,2] <- ""
  #TableDF <- cutoff.df
  #FontsDF <- c(2.5,2.5,2.5)
  #titleDF <- paste0(Sample,": Stats")
  #func.PlotTable.withRowNames.General(TableDF, FontsDF, titleDF, 20)
  
  #https://rpkgs.datanovia.com/ggpubr/files/ggtexttable-theme.pdf
  p2 <- ggtexttable(cutoff.df,theme = ttheme("mViolet"))
  #ggtexttable(cutoff.df, theme = ttheme("minimal"))
  #ggtexttable(cutoff.df, theme = ttheme("classic"))
  
  print("**** STEP3")
  
  SeuratObject <- NormalizeData(SeuratObject) %>%
    FindVariableFeatures(nfeatures = FeatureNum, verbose = FALSE) %>%
    ScaleData() %>%
    RunPCA(verbose=FALSE) %>% 
    FindNeighbors(dims = 1:PCAnum, verbose=FALSE) %>%
    FindClusters(resolution = resClus) %>% 
    RunUMAP(reduction='pca', dims=1:PCAnum, verbose=FALSE)
  
  setwd(saveDIR)
  pdf(file=paste0("QC_Regular_",Sample,".pdf"),height = 10,width = 12)
  
  if(plotCCgene==TRUE){
    if(Species=="hsa"){
      p3 <- FeaturePlot(SeuratObject, features = "TOP2A", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1)
    } else {
      p3 <- FeaturePlot(SeuratObject, features = "Top2a", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1)
    }
    
    p4.1 <- plot_grid(p2,p3, ncol = 1)
    print(plot_grid(p1,p4.1, ncol = 2, rel_widths = c(1.5,1)))
  } else if(plotCCgene==FALSE){
    
    print(plot_grid(p1,p2, ncol = 2, rel_widths = c(1.5,1)))
  }
  
  dev.off()
  
  print("**** STEP4")
  
  print(paste0("PrePorocessing data before running Doublet Detection Algorithms for: ",Sample))
  RUNProcessData="YES"
  if(RUNProcessData == "YES"){
    
    
    setwd(DDdir)
    pdf(file=paste0("PreProcess_DoubletFinder_Doublet_Detection_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    print(paste0("PreProcess Seurat object for: ",Sample))
    print(ElbowPlot(SeuratObject, ndims = 50))
    print(DimHeatmap(SeuratObject, dims = 1:12, cells = 500, balanced = TRUE))
    dev.off()
    
    SeuratObject.temp = SeuratObject
    
  }
  
  
  print("**** STEP5")
  
  if (DoubletFinder == TRUE) {
    
    setwd(DDdir)
    DoubletFinderDir <- paste(getwd(),paste0("DoubletFinder_",Sample),sep="/"); print(DoubletFinderDir)
    dir.create(file.path(getwd(),paste0("DoubletFinder_",Sample)), showWarnings = FALSE)
    
    setwd(DoubletFinderDir)
    
    print(paste0("Started DoubletFinder process for ",Sample, " using PCA ",PCAnum))
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_SeuratObject <- paramSweep_v3(SeuratObject, PCs = 1:PCAnum, sct = FALSE)
    sweep.stats_SeuratObject <- summarizeSweep(sweep.res.list_SeuratObject, GT = FALSE)
    bcmvn_SeuratObject <- find.pK(sweep.stats_SeuratObject)
    
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(SeuratObject$seurat_clusters)           ## ex: annotations <- SeuratObject@meta.data$ClusteringResults
    nExp_poi <- round(0.075*length(colnames(x = SeuratObject))); nExp_poi  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)); nExp_poi.adj
    dev.off()
    print("Finished Doublet Finder steps")
    
    print("**** STEP6.1")
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:PCAnum, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    head(SeuratObject@meta.data)
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'pANN_'] <- "pANNcomputed"
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'DF.classifications_'] <- "Doublet_LowConf"
    head(SeuratObject@meta.data)
    
    SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:PCAnum, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANNcomputed", sct = FALSE)
    head(SeuratObject@meta.data)
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'DF.classifications_'] <- "Doublet"
    head(SeuratObject@meta.data)
    #table(SeuratObject@meta.data$Doublet_LowConf)
    table(SeuratObject@meta.data$Doublet)
    #table(SeuratObject@meta.data$Doublet_LowConf, SeuratObject@meta.data$Doublet_HighConf)
    
    SeuratObject@meta.data$DoubletFinder <- SeuratObject@meta.data$Doublet
    DoubletCells <- SeuratObject@meta.data[SeuratObject@meta.data$DoubletFinder != "Singlet",]; dim(DoubletCells)
    SingletCells <- SeuratObject@meta.data[SeuratObject@meta.data$DoubletFinder == "Singlet",]; dim(SingletCells)
    #SeuratObject@meta.data[SeuratObject@meta.data$Doublet_LowConf == "Doublet", "DoubletFinder"] <- "Doublet_LowConf"
    #table(SeuratObject@meta.data$DoubletFinder)
    #SeuratObject@meta.data[SeuratObject@meta.data$Doublet_HighConf == "Doublet", "DoubletFinder"] <- "Doublet_HighConf"
    print(table(SeuratObject@meta.data$DoubletFinder))
    
    SeuratObject@meta.data$seurat_clusters <- MakeClustersFrom1(SeuratObject@meta.data$seurat_clusters, ClusOrder)
    
    #head(SeuratObject@meta.data)
    #cutoff.df <- data.frame(Doublets = table(SeuratObject@meta.data$DoubletFinder)); print(cutoff.df)
    #TableDF <- cutoff.df
    #FontsDF <- c(2.5,2.5,2.5)
    #titleDF <- paste0(Sample,": Doublets Detected")
    #func.PlotTable.General(TableDF, FontsDF, titleDF, 20)
    
    print("**** STEP6.2")
    
    cutoff.df.doub <- data.frame(cutoff.df, Doublets=nrow(DoubletCells))
    cutoff.df.doub <- rbind(cutoff.df, Doublets=c(nrow(DoubletCells), round((nrow(DoubletCells)/BeforeFilter*100),1)))
    cutoff.df.doub["QCpassCells","Cells"] = nrow(SingletCells)
    cutoff.df.doub["QCpassCells","Percentage"] = round((nrow(SingletCells)/BeforeFilter*100),1)
    cutoff.df.doub <- cutoff.df.doub[c(1,4,2,3,5),]
    
    print("**** STEP6.2.1")
    d1 <- DimPlot(SeuratObject, group.by = "DoubletFinder", cols = c("purple", "cyan3"), label = F, label.size = 7, pt.size = 0.1)
    print("**** STEP6.2.2")
    d2 <- VlnPlot(SeuratObject, features = "nFeature_RNA", pt.size = 0.1, group.by = "seurat_clusters", cols = c("purple", "cyan3"), split.by = "DoubletFinder") + ggtitle(paste0("nFeature_RNA, Doublets = ",nrow(DoubletCells),", (",round((nrow(DoubletCells)/BeforeFilter*100),1),"%)"))
    print("**** STEP6.2.3")
    #Idents(SeuratObject) <- "seurat_clusters"
    print(table(SeuratObject$seurat_clusters))
    #d3 <- DimPlot(SeuratObject, cols = ClusPallette, label = T, label.size = 7, pt.size = 0.1)
    print("**** STEP6.2.4")
    d4 <- FeaturePlot(SeuratObject, features = "percent.mt", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1)
    print("**** STEP6.2.5")
    d5 <- FeaturePlot(SeuratObject, features = "percent.rb", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1)
    print("**** STEP6.2.6")
    d6 <- FeaturePlot(SeuratObject, features = "nCount_RNA", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1)
    print("**** STEP6.2.7")
    d7 <- FeaturePlot(SeuratObject, features = "nFeature_RNA") + scale_colour_viridis_c(option = "plasma", direction = -1)
    print("**** STEP6.2.8")
    d8 <- FeaturePlot(SeuratObject, features = "pANNcomputed", pt.size = 0.1) + scale_colour_viridis_c(option = "plasma", direction = -1) + ggtitle(paste0("DoubletFinder Score"))
    print("**** STEP6.2.9")
    d9 <- ggtexttable(cutoff.df.doub,theme = ttheme("mViolet"))
    print("**** STEP6.2.10")
    print("**** STEP6.3")
    
    setwd(DoubletFinderDir)
    pdf(file=paste0("Plots_DoubletFinder_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    m1 <- plot_grid(d9, d2, nrow=1, rel_widths = c(1,2))
    #m2 <- plot_grid(d1,d3,d4, nrow=1)
    m2 <- plot_grid(d1,NULL,d4, nrow=1)
    m3 <- plot_grid(d5, d6,d7,d8, nrow=1)
    print(plot_grid(m1, m2, m3, nrow = 3))
    dev.off()
    
    print("**** STEP6.4")
    
    print("Plotted Doublet Information")
    
    #p1 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = ClusPallette)
    #p2 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "DoubletFinder", pt.size=1.5, cols=c("red", "dodgerblue", "black"))
    #print(plot_grid(p1, p2, NULL))
    
    write.table(SeuratObject@meta.data,file=paste0("DoubletFinder_Calls_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    write.table(DoubletCells,file=paste0("Doublets_DoubletFinder_Calls_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    write.table(SingletCells,file=paste0("Singlets_DoubletFinder_Calls_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
    print("**** STEP6.5")
    
  }
  
  SeuratObject@meta.data$SequencedCells <- BeforeFilter
  SeuratObject@meta.data$QCpass <- nrow(SingletCells)
  SeuratObject@meta.data$PassPercent <- round(nrow(SingletCells)/BeforeFilter*100,1)
  head(SeuratObject@meta.data)
  
  print("**** STEP7")
  
  setwd(saveDIR)
  saveRDS(SeuratObject, file = paste0(Sample,"_After_Doublets_using_PCA_",PCAnum,"_res_",resClus,".rds"))
  
  return(SeuratObject)
  
  print("Done")
  print(Sys.time())
}




MakeClustersFrom1 <- function(temp, ClusOrder){
print("Before Changing")
print(table(temp))
temp <- as.character(as.numeric(temp))
temp <- factor(temp, levels = ClusOrder[ClusOrder %in% temp])
print("After Changing")
print(table(temp))

return(temp)

}
