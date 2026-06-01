# =========================
# Load libraries
# =========================
library(Seurat)
library(ggplot2)

# =========================
# Load data
# =========================
seurat.vivo.exp2 <- readRDS(
  "data/Exp2(exp1)_vivo.RDS"
)

seurat.vitro.exp2 <- readRDS(
  "data/Exp2(exp1)_vitro.RDS"
)

seurat.vivo.cross <- readRDS(
  "data/CrossAge(exp2)_vivo.RDS"
)

seurat.vitro.cross <- readRDS(
  "data/CrossAge(exp2)_vitro.RDS"
)

# =========================
# Celltype reference
# =========================
celltype_ref <- c(
  "LT-HSC","ST-HSC","HSC",
  "MPP","MPP3","MPP4",
  "LMPP-1","LMPP-2",
  "CMP",
  "GMP","GMP-1","GMP-2","MDP",
  "MEP","MEP-1","MEP-2","EryP","MEP/EryP",
  "MkP","MKP","MKP-1"
)

# =========================
# Gene sets
# =========================
HSPC_genes <- c(
  "Mecom","Mllt3","Pdzk1ip1","Vwf","Sult1a1","Angpt1","Ifitm1","Meis1",
  "Esam","Hlf","Ly6a","Procr","Cdkn1c","Zfp367","Prkcdbp",
  "Itga2b","Mpl","Cd34","Flt3","Dntt","Ccl9",
  "Ctla2","Adgrl4","Rgs1","Hmga2",
  "Gcnt2","Wfdc17","Ifit1",
  "Elane","Cebpa","Csf1r","Irf8","Mpo","Prtn3","Calr","Ctsz","Fcgr3","Lgals1","Csrp3",
  "Car1","Car2","Mt1","Mt2","Trfc","Abcb7","Rhd","Casp3",
  "Epor","Gata1","Gata2","Klf1",
  "Mfsd2b","Vamp5","Gp5","Pf4","Gp1bb","Sdpr","Treml1"
)

genes_young_vivo_exp2 <- c(
  "Mecom","Mllt3","Meis1","Hlf","Ly6a","Procr",
  "Cd34","Flt3","Hmga2",
  "Elane","Cebpa","Mpo",
  "Gata1","Klf1","Epor",
  "Pf4","Gp5","Gp1bb","Treml1"
)

genes_young_vitro_exp2 <- c(
  "Mecom","Mllt3","Hlf","Ly6a","Procr","Cdkn1c",
  "Flt3","Hmga2","Cd34",
  "Elane","Cebpa","Mpo","Prtn3",
  "Car2","Klf1","Epor",
  "Pf4","Gp5","Gp1bb","Treml1"
)
genes_vivo_cross <- c(
  "Mecom","Mllt3","Hlf","Ly6a",
  "Cebpa","Mpo","Prtn3","Fcgr3",
  "Gata1","Klf1","Epor",
  "Pf4","Gp1bb"
)
# =========================
# Helper functions
# =========================
apply_celltype_ref <- function(obj, ref){
  keep_levels <- intersect(ref, unique(obj$celltype_final))
  obj <- subset(obj, subset = celltype_final %in% keep_levels)
  obj$celltype_final <- factor(obj$celltype_final, levels = keep_levels)
  Idents(obj) <- obj$celltype_final
  return(obj)
}

subset_by_ids <- function(obj, ids){
  subset(obj, subset = orig.ident %in% ids)
}

plot_hspc_heatmap <- function(obj, genes, title = NULL){
  
  genes_present <- intersect(genes, rownames(obj))
  
  p <- DoHeatmap(
    object   = obj,
    features = genes_present,
    group.by = "celltype_final",
    disp.min = -1.5,
    disp.max = 1.5,
    size     = 4
  ) +
    theme(
      legend.box = "vertical",
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

plot_hspc_dotplot <- function(obj, genes = HSPC_genes){
  DotPlot(
    obj,
    features = genes,
    group.by = "celltype_final",
    dot.scale = 6
  ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
      axis.text.y = element_text(size = 14),
      text = element_text(size = 14)
    ) +
    ylab("Cell Type") +
    xlab("") +
    scale_color_gradientn(
      colors = c(
        "blue4","lightskyblue","white",
        "lightpink","indianred1","brown3","brown4"
      )
    ) +
    guides(
      colour = guide_colourbar(title = "Average Expression"),
      size   = guide_legend(title = "Percent Expressed")
    )
}

save_plot <- function(name, plot, w = 10, h = 12){
  ggsave(name, plot, width = w, height = h, dpi = 300)
}

# =========================
# Apply celltype order
# =========================
seurat.vivo.exp2   <- apply_celltype_ref(seurat.vivo.exp2, celltype_ref)
seurat.vitro.exp2  <- apply_celltype_ref(seurat.vitro.exp2, celltype_ref)
seurat.vivo.cross  <- apply_celltype_ref(seurat.vivo.cross, celltype_ref)
seurat.vitro.cross <- apply_celltype_ref(seurat.vitro.cross, celltype_ref)

# =========================
# Subset Old / Young
# =========================
seurat.vitro.exp2.Old   <- subset_by_ids(seurat.vitro.exp2, c("4_Oa","5_Oa","6_Oa"))
seurat.vitro.exp2.Young <- subset_by_ids(seurat.vitro.exp2, c("1_Ya","2_Ya","3_Ya"))

seurat.vivo.exp2.Old    <- subset_by_ids(seurat.vivo.exp2, c("O1B","O2B","O3B"))
seurat.vivo.exp2.Young  <- subset_by_ids(seurat.vivo.exp2, c("Y1B","Y2B","Y3B"))

# =========================
# Heatmaps
# =========================
h_vivo_exp2        <- plot_hspc_heatmap(seurat.vivo.exp2, HSPC_genes, "vivo exp2")
h_vitro_exp2       <- plot_hspc_heatmap(seurat.vitro.exp2, HSPC_genes, "vitro exp2")

h_vivo_exp2_Old    <- plot_hspc_heatmap(seurat.vivo.exp2.Old, HSPC_genes, "vivo exp2 Old")
h_vitro_exp2_Old   <- plot_hspc_heatmap(seurat.vitro.exp2.Old, HSPC_genes, "vitro exp2 Old")

h_vivo_exp2_Young  <- plot_hspc_heatmap(seurat.vivo.exp2.Young, genes_young_vivo_exp2, "vivo exp2 Young")
h_vitro_exp2_Young <- plot_hspc_heatmap(seurat.vitro.exp2.Young, genes_young_vitro_exp2, "vitro exp2 Young")

h_vivo_cross  <- plot_hspc_heatmap(seurat.vivo.cross, genes_vivo_cross)
h_vitro_cross <- plot_hspc_heatmap(seurat.vitro.cross, HSPC_genes, "vitro cross")

# =========================
# DotPlots
# =========================
p_vivo_exp2        <- plot_hspc_dotplot(seurat.vivo.exp2)
p_vitro_exp2       <- plot_hspc_dotplot(seurat.vitro.exp2)

p_vivo_exp2_Old    <- plot_hspc_dotplot(seurat.vivo.exp2.Old)
p_vitro_exp2_Old   <- plot_hspc_dotplot(seurat.vitro.exp2.Old)

p_vivo_exp2_Young  <- plot_hspc_dotplot(seurat.vivo.exp2.Young, genes_young_vivo_exp2)
p_vitro_exp2_Young <- plot_hspc_dotplot(seurat.vitro.exp2.Young, genes_young_vitro_exp2)

p_vivo_cross  <- plot_hspc_dotplot(seurat.vivo.cross)
p_vitro_cross <- plot_hspc_dotplot(seurat.vitro.cross)

# =========================
# Save outputs
# =========================
save_plot("heatmap_vivo_exp2.pdf", h_vivo_exp2)
save_plot("heatmap_vitro_exp2.pdf", h_vitro_exp2)

save_plot("heatmap_vivo_cross.pdf", h_vivo_cross)
save_plot("heatmap_vitro_cross.pdf", h_vitro_cross)

save_plot("heatmap_vitro_exp2.png", h_vitro_exp2, 8, 6)
save_plot("heatmap_vivo_cross.png", h_vivo_cross, 8, 6)
