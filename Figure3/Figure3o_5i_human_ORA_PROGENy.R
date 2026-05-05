###################################################################################################


### Compute Differentially variance Genes (DVGs) in Healthy Human Bone marrow data 


###################################################################################################

# Maintainer: Chris Chen 
# Last updated: 03/20/2026

##########################

setwd('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts')

## Load packages ##

load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","scMayoMap",
    "readxl","pheatmap","Matrix","openxlsx","gt","glue","openxlsx"
  )
  
  suppressPackageStartupMessages(
    lapply(pkgs, require, character.only = TRUE)
  )
  message("Allpackages loaded.")
}
load_all_packages()

# set working directory
setwd('Main_figures/Figure3/Figures/')

#########################

# GSE180298

Ainciburu_et_al <- readRDS('GSE180298_harmonized.rds')
Ainciburu_et_al

#### Calculate heterogeneity between young and old ####

Ainciburu_et_al <- combined
table(Ainciburu_et_al$RNA_snn_res.0.2) # cluster 0 is HSC 
Idents(Ainciburu_et_al) <- "RNA_snn_res.0.2"
colnames(Ainciburu_et_al@meta.data)
table(Ainciburu_et_al$AgeGroup) # Old: 43987 and Young: 36905
DimPlot(Ainciburu_et_al, group.by = 'AgeGroup')
# Subset HSC 
Ainciburu_HSC <- subset(Ainciburu_et_al, idents = "0")
Ainciburu_HSC # 23278 genes and 26642 cells 
table(Ainciburu_HSC$AgeGroup) # Old: 20075 cells and Young: 6567 cells 
# Add donor age 
Ainciburu_HSC$donor_age <- ifelse(Ainciburu_HSC$AgeGroup == "Young","Y","O")
Ainciburu_HSC$donor_age <- factor(Ainciburu_HSC$donor_age,levels = c("Y", "O"))
table(Ainciburu_HSC$donor_age)
### Remove ribosomal and mt genes ###
# Define technical gene patterns
exclude_pattern <- "^MT-|^RPL|^RPS|^HBB|^HBA|^MALAT1"
genes_to_exclude <- grep(pattern = exclude_pattern,x = rownames(Ainciburu_HSC),value = TRUE)
length(genes_to_exclude) # 120 genes 
# Subset object to remove them
Ainciburu_HSC <- subset(Ainciburu_HSC,features = setdiff(rownames(Ainciburu_HSC), genes_to_exclude))
Ainciburu_HSC # 23158 genes 
## Run the compute variance and cv function 

External_HSC_DVG_results <- compute_variance_and_cv(
  seu = Ainciburu_HSC,
  group_col = "donor_age",
  assay = "RNA",
  slot = "data",
  n_hvg = 2000,
  min_pct = 0.1,
  min_mean = 0.1
)
# 520 genes 
head(External_HSC_DVG_results)

GSE180298_HSC_DVG_res <- External_HSC_DVG_results
External_HSC_DVG_results_clean <- External_HSC_DVG_results %>%
  dplyr::filter(
    !grepl("^AC", gene),
    !grepl("^LINC", gene))

external_hsc_plots <- plot_variance_cv_summaries(
  df = External_HSC_DVG_results_clean,
  fdr_cutoff = 0.1,
  prefix = "GSE180298 HSC")

external_hsc_plots$summary_counts
external_hsc_plots$volcano_mean_adj
external_hsc_plots$volcano_variance
external_hsc_plots$mean_adj_var_box
external_hsc_plots$density_mean_adj

########################################

# GSE189161

#### Calculate heterogeneity between young and old ####

bm_obj <- readRDS('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltype_annotation/Other_human_single_cell_bone_marrow_data/GSE189161_data/GSE189161_25yr_to_77yr_BM_harmonized.rds')
Idents(bm_obj) <- 'RNA_snn_res.0.2'
table(Idents(bm_obj))
bm_HSC <- subset(
  bm_obj,
  subset = RNA_snn_res.0.2 == "0")
bm_HSC # 40040 genes and 6884 cells 
table(bm_HSC$AgeGroup) # Old: 3667 and Young: 3217
# create donor_age 
bm_HSC$donor_age <- ifelse(
  bm_HSC$AgeGroup == "Young",
  "Y",
  "O"
bm_HSC$donor_age <- factor(
  bm_HSC$donor_age,
  levels = c("Y", "O"))
table(bm_HSC$donor_age)

# remove technical genes 
exclude_pattern <- "^MT-|^MT\\.|^RPL|^RPS|^HBB|^HBA|^MALAT1|^RP11"
genes_to_exclude <- grep(
  pattern = exclude_pattern,
  x = rownames(bm_HSC),
  value = TRUE)
length(genes_to_exclude) # 7817 genes removed 

bm_HSC <- subset(
  bm_HSC,
  features = setdiff(rownames(bm_HSC), genes_to_exclude))

bm_HSC # 32082 genes and 6884 cells 

# run heterogeneity analysis 
bm_HSC_DVG_results <- compute_variance_and_cv(
  seu = bm_HSC,
  group_col = "donor_age",
  assay = "RNA",
  slot = "data",
  n_hvg = 2000,
  min_pct = 0.1,
  min_mean = 0.1
)

nrow(bm_HSC_DVG_results) # 265 genes after filtering 
View(bm_HSC_DVG_results)

GSE189161_HSC_DVG_res <- bm_HSC_DVG_results

bm_HSC_DVG_results <- bm_HSC_DVG_results %>%
  dplyr::filter(!grepl("^AC07", gene))

# plot 
bm_hsc_plots <- plot_variance_cv_summaries(
  df = bm_HSC_DVG_results,
  fc_cutoff = 0.1,
  prefix = "GSE189161 Adult BM HSC"
)

bm_hsc_plots$summary_counts
bm_hsc_plots$volcano_mean_adj
bm_hsc_plots$volcano_variance
bm_hsc_plots$mean_adj_var_box
bm_hsc_plots$density_mean_adj

#################################################


#################################################

# Plot density plot and heterogeneity bar plot 

genes_of_interest_human <- c(
  "HSP90AA1","PCNA","KLF6","PRDX2","SMC4","TUBA1B","JUN","PDLIM1",
  "EGR1","HNRNPAB","ZFP36L2","JUNB","LYZ","RAB27B","NDUFS6",
  "PLEK","TUBB","MIF","CD69","MMRN1")
genes_of_interest_human <- intersect(
  genes_of_interest_human,
  rownames(bm_HSC))
plot_density_gene_human <- function(
    seurat_obj,
    gene,
    group_col = "donor_age",
    title_prefix = "Human BM HSCs",
    dvgs_df = bm_HSC_DVG_results,
    mean_adj_var_Y_col = "mean_adjusted_var_Y",
    mean_adj_var_O_col = "mean_adjusted_var_O"
){
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  expr_df <- FetchData(seurat_obj, vars = c(gene, group_col)) %>%
    mutate(Group = ifelse(.data[[group_col]] == "Y","Young (Y)","Old (O)")) %>%
    mutate(Expression = .data[[gene]]) %>%
    filter(Expression > 0)
  color_palette <- c(
    "Young (Y)" = "#89AEEB",
    "Old (O)"   = "#F5A36C"
  )
  p_density <- ggplot(expr_df,
                      aes(x = Expression, color = Group, fill = Group)) +
    geom_density(alpha = 0.25, linewidth = 1.1) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_classic(base_size = 15) +
    labs(
      title = paste0(title_prefix," — ",gene),
      x = "Normalized Expression",
      y = "Density"
    )
  if(gene %in% dvgs_df$gene){
    var_vals <- dvgs_df %>% filter(gene == !!gene)
    mean_var_df <- tibble(
      Group = c("Young (Y)","Old (O)"),
      MeanAdjVar = c(
        var_vals[[mean_adj_var_Y_col]][1],
        var_vals[[mean_adj_var_O_col]][1]
      )
    )
    p_bar <- ggplot(mean_var_df,
                    aes(x = Group, y = MeanAdjVar, fill = Group)) +
      geom_bar(stat="identity",width=0.6,color="black") +
      scale_fill_manual(values=color_palette) +
      theme_classic(base_size=15) +
      labs(title="Mean-Adj Var",x=NULL,y="Mean-Adj Var")
    p_density + p_bar + patchwork::plot_layout(widths=c(2.2,1))
  } else {
    p_density
  }
}

panel_list_human <- lapply(
  genes_of_interest_human,
  function(g){
    plot_density_gene_human(
      seurat_obj = bm_HSC,
      gene = g
    )
  }
)
names(panel_list_human) <- genes_of_interest_human
pdf("/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure5/Human_DVG_overlap/Human_HSC_cross_species_overlap_genes.pdf",
  width = 12,
  height = 6)
for(g in genes_of_interest_human){
  if(!is.null(panel_list_human[[g]])){
    print(panel_list_human[[g]])
  }
}
dev.off()

##############################################


##############################################

# Human DVG Over-representation analysis (ORA)

file_name <- "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure3/Tables/Human_HSC_DVG__dataset_summary.xlsx"
excel_sheets(file_name)
dvg_FC_0.5_df <- read_excel(file_name, sheet = "Overlap_FC>=0.05")
colnames(dvg_FC_0.5_df)
genes <- dvg_FC_0.5_df %>%
  pull(gene) %>%
  unique() %>%
  na.omit()
length(genes) # 49 genes 
head(genes)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
gene_df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db)
# GO ontologies
ontologies <- c("BP", "CC", "MF")
entrez_genes <- gene_df$ENTREZID # 48 
go_list <- lapply(ontologies, function(ont) {
  enrichGO(
    gene          = entrez_genes,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_list) <- ontologies
# KEGG
kegg_res <- enrichKEGG(
  gene         = entrez_genes,
  organism     = "hsa",
  pvalueCutoff = 0.05)
# function to clean enrich results
process_enrich <- function(res, source_name) {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df %>%
    mutate(
      Source = source_name,
      EnrichmentScore = -log10(p.adjust)
    )
}
# GO combined
go_df <- bind_rows(
  lapply(names(go_list), function(name) {
    process_enrich(go_list[[name]], paste0("GO_", name))
  })
)
# KEGG
kegg_df <- process_enrich(kegg_res, "KEGG")
# combine all
all_df <- bind_rows(go_df, kegg_df)
top_n_terms <- 10
plot_df <- all_df %>%
  group_by(Source) %>%
  slice_max(order_by = EnrichmentScore, n = top_n_terms) %>%
  ungroup()
plot_df <- plot_df %>%
  mutate(
    Description = str_wrap(Description, width = 40),
    Description = fct_reorder(Description, EnrichmentScore))
selected_terms <- c(
  "tumor necrosis factor-mediated signaling pathway",
  "regulation of hemopoiesis",
  "negative regulation of hematopoietic progenitor cell differentiation",
  "homeostasis of number of cells",
  "leukocyte homeostasis",
  "lymphocyte homeostasis",
  "lymphocyte differentiation",
  "DNA-binding transcription activator activity, RNA polymerase II-specific")

selected_ids <- c(
  "GO:0033209",  # TNF signaling pathway
  "GO:1903706",  # regulation of hemopoiesis
  "GO:1901533",  # negative regulation of hematopoietic progenitor cell differentiation
  "GO:0048872",  # homeostasis of number of cells
  "GO:0001776",  # leukocyte homeostasis
  "GO:0002260",  # lymphocyte homeostasis
  "GO:0030098",  # lymphocyte differentiation
  "GO:0001228"   # DNA-binding transcription activator activity (RNA pol II)
)
plot_df_2 <- plot_df
plot_df_2 <- plot_df_2 %>%
  dplyr::filter(
    ID %in% selected_ids,
    qvalue < 0.05
  ) %>%
  dplyr::mutate(
    log10padj = -log10(p.adjust)
  )
plot_df_2

library(forcats)
library(stringr)

# Prepare data

plot_df_2 <- plot_df_2 %>%
  dplyr::filter(
    ID %in% selected_ids,
    qvalue < 0.05
  ) %>%
  dplyr::mutate(
    log10padj = -log10(p.adjust)
  )
# scaling for dual axis
plot_df_2 <- plot_df_2 %>%
  mutate(
    Description = case_when(
      str_detect(Description, "DNA-binding transcription activator") ~
        "TF activator activity (RNA Pol II)",
      str_detect(Description, "negative regulation of hematopoietic") ~
        "neg. reg. of HPC differentiation",
      str_detect(Description, "tumor necrosis factor") ~
        "TNF signaling pathway",
      TRUE ~ Description))

# labels + ordering
plot_df_2 <- plot_df_2 %>%
  mutate(
    Description = str_wrap(Description, 40),
    Description = fct_reorder(Description, zScore, .desc = TRUE))

range_z <- range(plot_df_2$zScore)
range_p <- range(plot_df_2$log10padj)

plot_df_2 <- plot_df_2 %>%
  mutate(
    log10padj_scaled =
      scales::rescale(log10padj, to = range_z))

p_human_ora_bar <- ggplot(plot_df_2, aes(y = Description)) +
  # zScore + RichFactor color
  geom_col(aes(x = zScore, fill = RichFactor), width = 0.7) +
  # log10(adj p)
  geom_point(
    data = plot_df_2,
    aes(x = log10padj_scaled, y = Description),
    inherit.aes = FALSE,
    shape = 21,
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.5
  ) +
  # RichFactor color scale
  scale_fill_gradientn(
    colors = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
    values = scales::rescale(c(
      min(plot_df_2$RichFactor),
      quantile(plot_df_2$RichFactor, 0.25),
      median(plot_df_2$RichFactor),
      quantile(plot_df_2$RichFactor, 0.75),
      max(plot_df_2$RichFactor)
    )),
    oob = scales::squish,
    name = "RichFactor"
  ) +
  # Dual axis
  scale_x_continuous(
    name = "z-score",
    expand = expansion(mult = c(0, 0.08)),
    sec.axis = sec_axis(
      ~ scales::rescale(., to = range_p, from = range_z),
      name = "-log10(P.Adj)"
    )
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 22) +
  theme(
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x.top = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    legend.title = element_text(size = 20, color = "black", face = "bold"),
    legend.text = element_text(size = 18, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

p_human_ora_bar

####################

# PROGENy pathway scores 

bm_obj <- readRDS('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltype_annotation/Other_human_single_cell_bone_marrow_data/GSE189161_data/GSE189161_25yr_to_77yr_BM_harmonized.rds')
Idents(bm_obj) <- 'RNA_snn_res.0.2'
table(Idents(bm_obj))
bm_HSC <- subset(bm_obj,subset = RNA_snn_res.0.2 == "0")
bm_HSC # 40040 genes and 6884 cells 
table(bm_HSC$AgeGroup) # Old: 3667 and Young: 3217
if (!requireNamespace("progeny", quietly = TRUE)) {
  BiocManager::install("progeny")}
library(progeny)
# extract expression matrix 
expr_mat <- GetAssayData(bm_HSC, slot = "data")  # log-normalized
expr_mat <- as.matrix(expr_mat)
# run progeny 
progeny_scores <- progeny(
  expr_mat,
  scale = TRUE,
  organism = "Human",
  top = 500,
  perm = 1
)
bm_HSC <- AddMetaData(bm_HSC, progeny_scores)
colnames(bm_HSC@meta.data)
# format the data
progeny_df <- bm_HSC@meta.data %>%
  dplyr::select(
    AgeGroup,
    Androgen, EGFR, Estrogen, Hypoxia,
    `JAK-STAT`, MAPK, NFkB, p53,
    PI3K, TGFb, TNFa, Trail, VEGF, WNT
  ) %>%
  pivot_longer(
    cols = -AgeGroup,
    names_to = "Pathway",
    values_to = "Activity"
  )
#  Mean + delta
summary_df <- progeny_df %>%
  group_by(Pathway, AgeGroup) %>%
  summarise(mean_activity = mean(Activity), .groups = "drop") %>%
  pivot_wider(
    names_from = AgeGroup,
    values_from = mean_activity
  ) %>%
  mutate(delta = Old - Young)

#  Wilcoxon test
stat_df <- progeny_df %>%
  group_by(Pathway) %>%
  summarise(
    p_value = wilcox.test(Activity ~ AgeGroup)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

#  summary table 
final_df <- summary_df %>%
  left_join(stat_df, by = "Pathway") %>%
  arrange(p_adj)

### Link human DVGs with pathway scores 

pathways_use <- c(
  "EGFR",
  "p53",
  "NFkB",
  "TNFa",
  "Hypoxia",
  "VEGF",
  "WNT"
)
dvg_005 <- c(
  "STMN1","MPC2","ZFP36L2","PLEK","CXCR4","NR4A2","GATA2","SMC4","AREG",
  "PLAC8","MMRN1","DUSP1","HNRNPAB","HIST1H4C","TUBB","HSPA1A","TNFAIP3",
  "CITED2","TMSB4X","ITM2A","HMGN5","TSC22D3","BIRC3","REXO2","KLF6","SRGN",
  "DNAJC9","PDLIM1","PDZD8","OAT","GAPDH","CD69","TUBA1A","NR4A1","DUSP6",
  "HSP90B1","HMGB1","TNFSF13B","HSP90AA1","SLC12A6","CPPED1","HLF","PCNA",
  "SMIM24","JUNB","PRDX2","ZFP36","FOSB","MIF"
)
dvg_01 <- c(
  "STMN1","MPC2","ZFP36L2","CXCR4","NR4A2","SMC4","PLAC8","MMRN1","DUSP1",
  "HNRNPAB","HSPA1A","ITM2A","HMGN5","TSC22D3","BIRC3","REXO2","KLF6","SRGN",
  "DNAJC9","PDLIM1","PDZD8","OAT","CD69","TUBA1A","DUSP6","HSP90B1","HMGB1",
  "TNFSF13B","HSP90AA1","SLC12A6","CPPED1","PCNA","SMIM24","JUNB","PRDX2",
  "ZFP36","FOSB"
)
dvg_025 <- c(
  "PLAC8","KLF6","SRGN","CD69","TUBA1A","HMGB1","FOSB"
)
dvg_005 <- intersect(dvg_005, rownames(bm_HSC)) # 49 DVGs
dvg_01  <- intersect(dvg_01,  rownames(bm_HSC)) # 37 DVGs
dvg_025 <- intersect(dvg_025, rownames(bm_HSC)) # 7 DVGs
# add module score 
bm_HSC <- AddModuleScore(bm_HSC, features = list(dvg_005), name = "DVG005")
bm_HSC <- AddModuleScore(bm_HSC, features = list(dvg_01),  name = "DVG01")
bm_HSC <- AddModuleScore(bm_HSC, features = list(dvg_025), name = "DVG025")
score_names <- c("DVG0051", "DVG011", "DVG0251")

cor_df <- data.frame()

meta <- bm_HSC@meta.data
for (score in score_names) {
  for (p in pathways_use) {
    x <- as.numeric(meta[[score]])
    y <- as.numeric(meta[[p]])
    # remove NA
    keep <- complete.cases(x, y)
    tmp <- cor.test(
      x[keep],
      y[keep],
      method = "spearman"
    )
    cor_df <- rbind(
      cor_df,
      data.frame(
        DVG_set = score,
        Pathway = p,
        rho = tmp$estimate,
        p = tmp$p.value
      )
    )
  }
}
cor_df$p_adj <- p.adjust(cor_df$p, method = "BH")

# addlabels
cor_df$DVG_set <- recode(cor_df$DVG_set,
                         "DVG0051" = "FC ≥ 0.05",
                         "DVG011"  = "FC ≥ 0.1",
                         "DVG0251" = "FC ≥ 0.25")
cor_df



###################################################


# End of this script 


