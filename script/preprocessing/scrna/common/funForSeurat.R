FindAgingMarkers3 <- function(cluster,hspc.combined) {
  hspc.combined$cluster.AGE <- paste(Idents(object = hspc.combined), hspc.combined$AGE, sep = "_")
  Idents(object = hspc.combined) <- "cluster.AGE"
  agingMarkers <- FindMarkers(object = hspc.combined, 
                              ident.1 = paste(cluster,"_Old",sep = ""), 
                              ident.2 = paste(cluster,"_Young",sep = ""), 
                              verbose = FALSE)
  agingMarkers$Cluster <- paste(cluster,"_Old_up",sep ="")
  agingMarkers$Cluster[which(agingMarkers$avg_logFC < 0)] <- paste(cluster,"_Old_down",sep ="")
  agingMarkers$Gene <- rownames(agingMarkers)
  agingMarkers <- agingMarkers[order(agingMarkers$avg_logFC),]
  agingMarkers <- agingMarkers[which(agingMarkers$p_val_adj < 0.05),]
  return(agingMarkers)
}




getAGEPropPerClustBarplot <- function(hspc.combined,age_levels=c("young","mid","old","vold")) {
  clusterAge <- ddply(hspc.combined@meta.data,~numclust + AGE,nrow)
  
  propExpect <- table(hspc.combined@meta.data$AGE)/length(hspc.combined@meta.data$AGE)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$AGE)[1]]]
  
  #clusterAGE$numclust <- factor(v$clusterNature , levels = c(""))
  clusterAge$AGE <- factor(clusterAge$AGE , levels = age_levels) 
  
  
  AGE <- ggplot(data.frame(clusterAge), aes(fill = AGE,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$AGE)))))+
    scale_y_continuous(name = "Age (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() + geom_hline(yintercept = propYoungExp)+
    theme(legend.title=element_blank())
  return(AGE)
  
}


getSamplePropPerClustBarplot <- function(hspc.combined) {
  hspc.combined@meta.data$dataset <- hspc.combined@meta.data$sampleName
  clustersampleName <- ddply(hspc.combined@meta.data,~numclust + dataset,nrow)
  
  propExpect <- table(hspc.combined@meta.data$dataset)/length(hspc.combined@meta.data$dataset)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$dataset)[1]]]
  
  #clustersampleName$numclust <- factor(v$clusterNature , levels = c(""))
  #clustersampleName$sampleName <- factor(clustersampleName$predicted , levels = c("")) 
  
  
  sampleName <- ggplot(data.frame(clustersampleName), aes(fill = dataset,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$dataset)))))+
    scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank()) 
  return(sampleName)
  
}




getPredictedPropPerClustBarplot <- function(hspc.combined) {
  clusterpredicted <- ddply(hspc.combined@meta.data,~numclust + predicted,nrow)
  
  propExpect <- table(hspc.combined@meta.data$predicted)/length(hspc.combined@meta.data$predicted)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$predicted)[1]]]
  
  #clusterpredicted$numclust <- factor(v$clusterNature , levels = c(""))
  #clusterpredicted$predicted <- factor(clusterpredicted$predicted , levels = c("")) 
  
  
  predicted <- ggplot(data.frame(clusterpredicted), aes(fill = predicted,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= hue_pal()(length(unique(hspc.combined@meta.data$predicted))))+
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank())  
  return(predicted)
  
}


getPhasePropPerClustBarplot <- function(hspc.combined) {
  clusterphases <- ddply(hspc.combined@meta.data,~numclust + phases,nrow)
  
  propExpect <- table(hspc.combined@meta.data$phases)/length(hspc.combined@meta.data$phases)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$phases)[1]]]
  
  #clusterphases$numclust <- factor(v$clusterNature , levels = c(""))
  #clusterphases$phases <- factor(clusterphases$predicted , levels = c("")) 
  
  
  phases <- ggplot(data.frame(clusterphases), aes(fill = phases,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$phases)))))+
    scale_y_continuous(name = "Cell cycle phase (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() +
    theme(legend.title=element_blank())  
  return(phases)
  
}
  
