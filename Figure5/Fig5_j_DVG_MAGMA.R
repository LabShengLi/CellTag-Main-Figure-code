################################################################################

## Prepare DVG gene signatures for MAGMA and BeatAML somatic mutation analysis

###############################################################################

load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","scMayoMap",
    "readxl","pheatmap","Matrix","openxlsx","gt","glue","writexl"
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
setwd('...')
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

mouse_fc005_human <- mouse_fc005 %>%
  left_join(ortholog_map, by = c("gene" = "mouse_gene")) %>%
  filter(!is.na(human_gene)) %>%
  pull(human_gene) %>%
  unique() # 83 genes 

mouse_fc01_human <- mouse_fc01 %>%
  left_join(ortholog_map, by = c("gene" = "mouse_gene")) %>%
  filter(!is.na(human_gene)) %>%
  pull(human_gene) %>%
  unique() # 32 genes 

human_fc005_genes <- unique(human_fc005$gene) # 49 genes 
human_fc01_genes  <- unique(human_fc01$gene) # 37 genes 

# Ovarlapped genes in human form 
listA <- strsplit(
  "TSC22D1; UQCR11; CSF3R; HSPA8; HMGB2; HSPH1; SET; UBB; JUNB; MED12L; PCLAF; LDHA; LSP1; MCM5; MINPP1; MT2A; NOP10; PNP; SAT1; TXNIP; ZEB2; CAVIN2; FLT3; RGS1; NFIA; UNG; CKS2; TMEM14B; TMEM14C; SDF2L1; NASP; BOLA3; CKS1B; CITED2; HLF; KIF22; TIMM10; CD69; CD9; CKAP2; HDGF; LAT2; MEG3; MMRN1",
  ";")[[1]]

listB <- strsplit(
  "HSP90AA1; PCNA; KLF6; PRDX2; SMC4; TUBA1B; JUN; PDLIM1; EGR1; HNRNPAB; ZFP36L2; JUNB; LYZ; RAB27B; NDUFS6; PLEK; TUBB; MIF; CD69; MMRN1",
  ";")[[1]]

aging_overlap <- unique(c(listA, listB))
aging_overlap <- trimws(unique(c(listA, listB)))
aging_overlap # 61 genes 

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

length(ent_mouse_fc005) # 78 
length(ent_mouse_fc01) # 31 
length(ent_human_fc005) # 43 
length(ent_human_fc01) # 33
length(ent_overlap) # 58 

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
  "magma_DVG_summary_ALLGWAS_0130.tsv",
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

write_xlsx(magma_DVG_summary,"magma_DVG_summary.xlsx")

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
  mutate(sig = ifelse(P < 0.05, "**", ""))

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

# heatmap
p <- ggplot(magma, aes(GeneSet, GWAS, fill = BETA)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = sig), size = 6, fontface = "bold") +
  scale_x_discrete(
    labels = c(
      "CHIP ALL" = expression(italic("CHIP ALL")),
      "CHIP DNMT3A" = expression(italic("CHIP DNMT3A")),
      "CHIP TET2" = expression(italic("CHIP TET2")),
      "AML" = expression(italic("AML")),
      "MDS" = expression(italic("MDS"))
    ),
    expand = c(0,0)
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('", x, "')")),
    expand = c(0,0)
  ) +
  scale_fill_gradientn(
    colors = c("#95CEBE","#DAEAD7","#EDBF97","#D4897E"),
    values = scales::rescale(c(-0.3, 0, 0.3)),
    limits = c(-0.3, 0.3),
    breaks = c(-0.3, 0, 0.3),
    labels = c("-0.3", "0", "0.3"),
    name = expression(bold("MAGMA BETA"))
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(size = 20, face = "bold", color = "black"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text  = element_text(color = "black"),
    panel.grid = element_blank()
  )
p

###########################################

# Check on key genes on: TET2_CH_hg19; MDS_UKB_hg19
library(data.table)
# use gene set: human_fc005_genes; aging_overlap; mouse_fc005_human

aging_overlap
human_fc005_genes
mouse_fc005_human

TET2_CH_hg19_DVG_MouseHuman_AgingOverlap <- fread("TET2_CH_hg19_DVG_MouseHuman_AgingOverlap_human_Entrez.gsa.genes.out")
head(TET2_CH_hg19_DVG_MouseHuman_AgingOverlap)

TET2_CH_hg19_DVG_Human_FC0.05 <- fread("TET2_CH_hg19_DVG_Human_FC0.05_human_Entrez.gsa.genes.out")
head(TET2_CH_hg19_DVG_Human_FC0.05)

MDS_UKB_hg19_DVG_Mouse_FC0.05 <- fread('MDS_UKB_hg19_DVG_Mouse_FC0.05_2exp_human_Entrez.gsa.genes.out')
head(MDS_UKB_hg19_DVG_Mouse_FC0.05)

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

TET2_CH_hg19_DVG_MouseHuman_AgingOverlap <- convert_symbol(TET2_CH_hg19_DVG_MouseHuman_AgingOverlap)
TET2_CH_hg19_DVG_MouseHuman_AgingOverlap
TET2_CH_hg19_DVG_Human_FC0.05 <- convert_symbol(TET2_CH_hg19_DVG_Human_FC0.05)
TET2_CH_hg19_DVG_Human_FC0.05
MDS_UKB_hg19_DVG_Mouse_FC0.05 <- convert_symbol(MDS_UKB_hg19_DVG_Mouse_FC0.05)
MDS_UKB_hg19_DVG_Mouse_FC0.05

rank_genes <- function(df){
  df %>%
    filter(ZSTAT > 0) %>%
    mutate(absZ = abs(ZSTAT)) %>%
    arrange(desc(absZ))}

rank_genes(TET2_CH_hg19_DVG_MouseHuman_AgingOverlap) %>%
  dplyr::select(SYMBOL, ZSTAT, absZ) %>%
  head(20)

get_drivers <- function(df, geneset){
  df %>%
    dplyr::filter(SYMBOL %in% geneset) %>%
    dplyr::filter(ZSTAT > 0) %>%
    dplyr::mutate(absZ = abs(ZSTAT)) %>%
    dplyr::arrange(desc(absZ)) %>%
    dplyr::select(SYMBOL, CHR, START, STOP, ZSTAT, NSNPS)
}

# Aging overlap in TET2 CH 
drivers_aging <- get_drivers(
  TET2_CH_hg19_DVG_MouseHuman_AgingOverlap,
  aging_overlap
)
View(head(drivers_aging,20))

# Human 0.05 in TET2 CH
drivers_human_fc005 <- get_drivers(
  TET2_CH_hg19_DVG_Human_FC0.05,
  human_fc005_genes
)

View(head(drivers_human_fc005, 20))

# Mouse 0.05 in MDS 
drivers_mouse_fc005_mds <- 
  get_drivers(
    MDS_UKB_hg19_DVG_Mouse_FC0.05,
    mouse_fc005_human
  )

rank_genes(MDS_UKB_hg19_DVG_Mouse_FC0.05) %>%
  dplyr::select(SYMBOL, ZSTAT, absZ) %>%
  head(20)
View(head(drivers_mouse_fc005_mds, 20))

#############################################

# Plot Z-score density plots 


plot_z_density <- function(df, geneset, title){
  plot_df <- df %>%
    mutate(group = ifelse(SYMBOL %in% geneset, "Gene set", "Other genes"))
  mean_set <- mean(plot_df$ZSTAT[plot_df$group == "Gene set"], na.rm = TRUE)
  mean_other <- mean(plot_df$ZSTAT[plot_df$group == "Other genes"], na.rm = TRUE)
  
  ggplot(plot_df, aes(x = ZSTAT, fill = group)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = mean_set,
               linetype = "dashed",
               color = "firebrick",
               size = 1) +
    geom_vline(xintercept = mean_other,
               linetype = "dashed",
               color = "navy",
               size = 1) +
    scale_fill_manual(
      values = c("Gene set" = "firebrick", "Other genes" = "grey70")
    ) +
    labs(
      x = "MAGMA Z score",
      y = "Density",
      fill = "",
      title = title
    ) +
    theme_classic(base_size = 20)
}

plot_z_density(
  TET2_CH_hg19_DVG_MouseHuman_AgingOverlap,
  aging_overlap,
  "TET2 CHIP — Aging overlap genes"
)

plot_z_density(
  TET2_CH_hg19_DVG_Human_FC0.05,
  human_fc005_genes,
  "TET2 CHIP — Human FC0.05 genes"
)
plot_z_density(
  MDS_UKB_hg19_DVG_Mouse_FC0.05,
  mouse_fc005_human,
  "MDS — Mouse FC0.05 genes"
)

##########################################################

# End of the script 


