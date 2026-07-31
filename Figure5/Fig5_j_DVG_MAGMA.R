################################################################################
## Prepare DVG gene signatures for MAGMA and BeatAML somatic mutation analysis
###############################################################################
load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","scMayoMap",
    "readxl","pheatmap","Matrix","openxlsx","gt","glue","writexl","patchwork"
  )
  
  suppressPackageStartupMessages(
    lapply(pkgs, require, character.only = TRUE)
  )
  message("Allpackages loaded.")
}
load_all_packages()
library(org.Hs.eg.db)
library(AnnotationDbi)

# set working directory
setwd('.../Data')

############################################
# read in data 
Mouse_Human_DVGs <- 'file_input/Fig3_DVG_Master_Table.xlsx'
excel_sheets(Mouse_Human_DVGs)
file <- Mouse_Human_DVGs
mouse_fc005 <- read_excel(file, sheet = "Mouse_gt_0.05FC_2_exp")
mouse_fc01  <- read_excel(file, sheet = "Mouse_gt_0.1FC_2_exp")
human_fc005 <- read_excel(file, sheet = "Human_Overlap_FC>=0.05")
human_fc01  <- read_excel(file, sheet = "Human_Overlap_FC>=0.1")
ortholog_map <- read_excel("file_input/mouse_human_ortholog_map.xlsx")
colnames(ortholog_map)

ortholog_map <- ortholog_map %>%
  tidyr::separate_rows(human_gene, sep = "\\s*;\\s*") %>%
  mutate(human_gene = trimws(human_gene)) %>%
  filter(human_gene != "") %>%
  distinct()

mouse_fc005_human <- mouse_fc005 %>%
  left_join(ortholog_map, by = c("gene" = "mouse_gene")) %>%
  filter(!is.na(human_gene)) %>%
  pull(human_gene) %>%
  unique() # 69 genes 

mouse_fc01_human <- mouse_fc01 %>%
  left_join(ortholog_map, by = c("gene" = "mouse_gene")) %>%
  filter(!is.na(human_gene)) %>%
  pull(human_gene) %>%
  unique() # 31 genes 
human_fc005_genes <- unique(human_fc005$gene) # 37 genes 
human_fc01_genes  <- unique(human_fc01$gene) # 26 genes 

length(mouse_fc005_human) # 69 
length(mouse_fc01_human) # 31
length(human_fc005_genes) # 37 
length(human_fc01_genes) # 26

# Shared (aging-overlap) DVGs 
aging_overlap <- unique(c(
  "IGFBP7","ANGPT1","APP","HSPA5","PIM2","ID1","PRDX2","HIGD1A","SELENOP",
  "SELENOM","SMC4","HSPD1","STMN1","CD63","JUN","PDLIM1","MYC","HSPE1",
  "NCL","C1QBP","CCT6A","IER2","SMC2","NDUFAF8","RAN","NUCKS1","HNRNPAB",
  "ZFP36L2","PBX1","NAA38","TMSB10","CSF3R","HMGB2","YBX1","HSPH1","SET",
  "UBB","JUNB","PCLAF","LSP1","MCM5","UQCRQ","NME1","NOP10","PNP","PDIA6",
  "TXNIP","CAVIN2","FLT3","RGS1","NFIA","UNG","CKS2","EIF5A","BOLA3",
  "CKS1B","KIF22","PLEK","CLU","TUBB","MIF","CD69","ALDH1A1","CD9","CKAP2",
  "CPA3","HDGF","MEG3","MMRN1","SNRPG","TUBA1B","EGR1","MLLT3","NFKBIA",
  "RAB27B","NDUFS6"
))
aging_overlap
length(aging_overlap) # 76 genes

# Convert gene symbols to Entrez id 
convert_to_entrez <- function(genes){
  entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = genes,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  )
  entrez %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
}
ent_mouse_fc005 <- convert_to_entrez(mouse_fc005_human) # 84 
ent_mouse_fc01  <- convert_to_entrez(mouse_fc01_human) # 32 
ent_human_fc005 <- convert_to_entrez(human_fc005_genes) # 48 
ent_human_fc01  <- convert_to_entrez(human_fc01_genes) # 37 
ent_overlap     <- convert_to_entrez(aging_overlap) # 61  

# Filter to MAGMA gene background 
magma_genes <- read.table("file_input/all_genes.entrez",
                          comment.char = "#")
colnames(magma_genes) <- "ENTREZID"
filter_magma <- function(entrez_ids){
  entrez_ids[entrez_ids %in% magma_genes$ENTREZID]
}
ent_mouse_fc005 <- filter_magma(ent_mouse_fc005)
ent_mouse_fc01  <- filter_magma(ent_mouse_fc01)
ent_human_fc005 <- filter_magma(ent_human_fc005)
ent_human_fc01  <- filter_magma(ent_human_fc01)
ent_overlap     <- filter_magma(ent_overlap)
length(ent_mouse_fc005) # 67 
length(ent_mouse_fc01) # 29 
length(ent_human_fc005) # 34 
length(ent_human_fc01) # 25
length(ent_overlap) # 73 

# Write MAGMA input 
write_ids <- function(ids, file){
  write.table(ids,
              file,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}

write_ids(ent_mouse_fc005,"DVG_Mouse_FC0.05_2exp_human_Entrez.txt")
write_ids(ent_mouse_fc01,"DVG_Mouse_FC0.1_2exp_human_Entrez.txt")
write_ids(ent_human_fc005,"DVG_Human_FC0.05_human_Entrez.txt")
write_ids(ent_human_fc01,"DVG_Human_FC0.1_human_Entrez.txt")
write_ids(ent_overlap,"DVG_MouseHuman_AgingOverlap_human_Entrez.txt")

####################################################################################
magma_DVG_summary <- read.delim(
  ".../MAGMA/magma_DVG_summary_ALLGWAS_0130.tsv",
  stringsAsFactors = FALSE
)
# remove duplicated header row
magma_DVG_summary <- magma_DVG_summary %>%
  filter(X.e.GWAS != "GeneSet")
magma_DVG_summary <- magma_DVG_summary %>%
  filter(!NGENES %in% c("NGENES"))
magma_DVG_summary <- magma_DVG_summary %>%
  filter(!BETA %in% c("BETA"))
magma_DVG_summary <- magma_DVG_summary %>%
  filter(!SE %in% c("SE"))
magma_DVG_summary <- magma_DVG_summary %>%
  filter(!P %in% c("P"))
# rename first column
colnames(magma_DVG_summary)[1] <- "GWAS"
magma_DVG_summary
# convert numeric columns
magma_DVG_summary <- magma_DVG_summary %>%
  mutate(
    NGENES = as.numeric(NGENES),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P)
  )
magma_DVG_summary
write_xlsx(magma_DVG_summary,".../MAGMA/magma_DVG_summary.xlsx")

###############################
# Enrichment heatmap 
library(dplyr)
library(ggplot2)
magma <- magma_DVG_summary
magma
magma <- magma %>%
  dplyr::filter(!GeneSet %in% c("DVG_Mouse_FC0.1_2exp_human_Entrez", "DVG_Human_FC0.1_human_Entrez"))
# shorten gene set names
magma$GeneSet <- magma$GeneSet %>%
  gsub("_human_Entrez", "", .) %>%
  gsub("DVG_", "", .) %>%
  gsub("_2exp", "", .)
table(magma$GeneSet)
magma <- magma %>%
  dplyr::mutate(
    GeneSet = dplyr::recode(
      GeneSet,
      "MouseHuman_AgingOverlap" = "Shared DVG",
      "Mouse_FC0.05"              = "Mouse DVG",
      "Human_FC0.05"              = "Human DVG"
    )
  )
magma$GeneSet <- factor(
  magma$GeneSet,
  levels = c(
    "Human DVG",
    "Shared DVG",
    "Mouse DVG"
  )
)
magma <- magma %>%
  mutate(GWAS = dplyr::recode(GWAS,
                              "Overall_CH_hg19" = "CHIP ALL",
                              "DNMT3A_CH_hg19" = "CHIP DNMT3A",
                              "TET2_CH_hg19" = "CHIP TET2",
                              "AML_UKB_hg19" = "AML",
                              "MDS_UKB_hg19" = "MDS"))

# significance annotation
magma <- magma %>%
  mutate(sig = case_when(
    P < 0.01 ~ "***",
    P < 0.05  ~ "**",
    P < 0.1  ~ "*",
    TRUE ~ ""
  ))
View(magma)

# ordering
magma$GWAS <- factor(
  magma$GWAS,
  levels = c(
    "CHIP ALL",
    "CHIP DNMT3A",
    "CHIP TET2",
    "AML",
    "MDS"
  )
)

library(ggtext)

# GWAS sample sizes — confirm these match your current GWAS inputs
gwas_labels <- c(
  "CHIP ALL"    = "CHIP ALL<br><span style='font-size:16pt'>n = 10,203</span>",
  "CHIP DNMT3A" = "CHIP *DNMT3A*<br><span style='font-size:16pt'>n = 5,185</span>",
  "CHIP TET2"   = "CHIP *TET2*<br><span style='font-size:16pt'>n = 2,041</span>",
  "AML"         = "AML<br><span style='font-size:16pt'>n = 312</span>",
  "MDS"         = "MDS<br><span style='font-size:16pt'>n = 521</span>"
)

# heatmap
p <- ggplot(magma, aes(GeneSet, GWAS, fill = BETA)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = sig), size = 6, fontface = "bold", color = "black") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(labels = gwas_labels, expand = c(0,0)) +
  scale_fill_gradientn(
    colors = c("#95CEBE","#DAEAD7","#EDBF97","#D4897E"),
    values = scales::rescale(c(-0.3, 0, 0.3)),
    limits = c(-0.3, 0.3),
    breaks = c(-0.3, 0, 0.3),
    labels = c("-0.3", "0", "0.3"),
    name = "MAGMA\nBETA"
  ) +
  labs(title = "GWAS enrichment (MAGMA)") +
  theme_classic(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black"),
    axis.text.y = ggtext::element_markdown(size = 18, color = "black", lineheight = 1.1),
    axis.title = element_blank(),
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black"),
    panel.grid = element_blank()
  )
p

###########################################
# Check on key genes for all combos significant at p < 0.1
library(data.table)
# use gene set: human_fc005_genes; aging_overlap; mouse_fc005_human
aging_overlap
human_fc005_genes
mouse_fc005_human

setwd('.../MAGMA')
TET2_CH_hg19_DVG_Human_FC0.05 <- fread("TET2_CH_hg19_DVG_Human_FC0.05_human_Entrez.gsa.genes.out")
AML_UKB_hg19_DVG_Human_FC0.05 <- fread("AML_UKB_hg19_DVG_Human_FC0.05_human_Entrez.gsa.genes.out")
MDS_UKB_hg19_DVG_Mouse_FC0.05 <- fread("MDS_UKB_hg19_DVG_Mouse_FC0.05_2exp_human_Entrez.gsa.genes.out")


# convert the entrez id back to the gene names 
convert_symbol <- function(df){
  df$SYMBOL <- mapIds(
    org.Hs.eg.db,
    keys = as.character(df$GENE),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  return(df)
}
TET2_CH_hg19_DVG_Human_FC0.05             <- convert_symbol(TET2_CH_hg19_DVG_Human_FC0.05)
AML_UKB_hg19_DVG_Human_FC0.05             <- convert_symbol(AML_UKB_hg19_DVG_Human_FC0.05)
MDS_UKB_hg19_DVG_Mouse_FC0.05             <- convert_symbol(MDS_UKB_hg19_DVG_Mouse_FC0.05)


rank_genes <- function(df){
  df %>%
    filter(ZSTAT > 0) %>%
    mutate(absZ = abs(ZSTAT)) %>%
    arrange(desc(absZ))
}

get_drivers <- function(df, geneset){
  df %>%
    dplyr::filter(SYMBOL %in% geneset) %>%
    dplyr::filter(ZSTAT > 0) %>%
    dplyr::mutate(absZ = abs(ZSTAT)) %>%
    dplyr::arrange(desc(absZ)) %>%
    dplyr::select(SYMBOL, CHR, START, STOP, ZSTAT, NSNPS)
}

# CHIP TET2 × Human DVG (p = 0.041)
drivers_tet2_human <- get_drivers(TET2_CH_hg19_DVG_Human_FC0.05, human_fc005_genes)
View(head(drivers_tet2_human, 20))
# AML × Human DVG (p = 0.041)
drivers_aml_human <- get_drivers(AML_UKB_hg19_DVG_Human_FC0.05, human_fc005_genes)
View(head(drivers_aml_human, 20))
# MDS × Mouse DVG (p = 0.018)
drivers_mds_mouse <- get_drivers(MDS_UKB_hg19_DVG_Mouse_FC0.05, mouse_fc005_human)
View(head(drivers_mds_mouse, 20))
# MDS × Human DVG (p = 0.081)
drivers_mds_human <- get_drivers(MDS_UKB_hg19_DVG_Mouse_FC0.05, human_fc005_genes)
View(head(drivers_mds_human, 20))


#############################################
# Plot Z-score density plots 
plot_z_density <- function(df, geneset, title){
  plot_df <- df %>%
    mutate(group = ifelse(SYMBOL %in% geneset, "Gene set", "Other genes"))
  mean_set <- mean(plot_df$ZSTAT[plot_df$group == "Gene set"], na.rm = TRUE)
  mean_other <- mean(plot_df$ZSTAT[plot_df$group == "Other genes"], na.rm = TRUE)
  
  ggplot(plot_df, aes(x = ZSTAT, fill = group)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = mean_set, linetype = "dashed", color = "firebrick", size = 1) +
    geom_vline(xintercept = mean_other, linetype = "dashed", color = "navy", size = 1) +
    scale_fill_manual(values = c("Gene set" = "firebrick", "Other genes" = "grey70")) +
    labs(x = "MAGMA Z score", y = "Density", fill = "", title = title) +
    theme_classic(base_size = 20) + 
    theme(
      plot.title = element_text(size = 18)
    )
}

p1 <- plot_z_density(TET2_CH_hg19_DVG_Human_FC0.05, human_fc005_genes, "CHIP TET2 — Human DVG")
p2 <- plot_z_density(AML_UKB_hg19_DVG_Human_FC0.05, human_fc005_genes, "AML — Human DVG")
p3 <- plot_z_density(MDS_UKB_hg19_DVG_Mouse_FC0.05, mouse_fc005_human, "MDS — Mouse DVG")
p4 <- plot_z_density(MDS_UKB_hg19_DVG_Mouse_FC0.05, human_fc005_genes, "MDS — Human DVG")

combined_z_density <- (p1 | p2) / (p3 | p4) + patchwork::plot_layout(guides = "collect")
combined_z_density

##########################################################
# End of the script