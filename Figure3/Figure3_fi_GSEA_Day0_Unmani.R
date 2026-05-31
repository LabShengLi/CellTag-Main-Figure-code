###################################################################################################


### Functional Analysis: GSEA and ORA ###


###################################################################################################

# Maintainer: Chris Chen 
# Last updated: 02/10/2026

##########################

## Load packages ##

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

# set working directory
setwd('...')

########################################

### Perform GSEA for CrossAge Day0 data ###

# Read in DVG file

file_name <- "Fig3_DVG_Master_Table.xlsx"
excel_sheets(file_name)
Day0_All <- read_excel(file_name, sheet = 1) # 789 genes
Unmanipulated_All <- read_excel(file_name, sheet = 2) # 687 genes
# Exp1_All <- read_excel(file_name, sheet = 3)
Exp2_All <- read_excel(file_name, sheet = 3) # 782 genes 

########################################

## Prepare DVG gene list ##
Exp2_All
dvg_list <- list(Day0 = Day0_All, Unmanipulated = Unmanipulated_All, Exp2 = Exp2_All)
get_sig_table <- function(df, logfc_cut = 0) {df %>% 
    dplyr::filter(log2FC_mean_adjusted_variance_OY > logfc_cut,fdr_brown_forsythe < 0.05) %>%
    dplyr::select(gene, log2FC_mean_adjusted_variance_OY)}

combined_table <- imap_dfr(dvg_list, ~ get_sig_table(.x, logfc_cut = 0) %>%
    mutate(experiment = .y))
head(combined_table)
union_df <- combined_table %>% group_by(gene) %>%
  summarise(mean_log2FC = mean(log2FC_mean_adjusted_variance_OY), n_exp = n(),
    .groups = "drop") %>% arrange(desc(mean_log2FC))

get_genes_2exp <- function(logfc_cut) {df <- imap_dfr(dvg_list,
    ~ get_sig_table(.x, logfc_cut) %>% mutate(experiment = .y))
  df %>% group_by(gene) %>% summarise(mean_log2FC = mean(log2FC_mean_adjusted_variance_OY),
      n_exp = n(), .groups = "drop") %>% filter(n_exp >= 2) %>% arrange(desc(mean_log2FC))}

exp2_005_df <- get_genes_2exp(0.05)   # ~80 genes
exp2_01_df  <- get_genes_2exp(0.1)    # ~33 genes
exp2_025_df <- get_genes_2exp(0.25)   # 3 genes

head(exp2_005_df)
head(exp2_01_df)
head(exp2_025_df)

count_table <- data.frame(Gene_List = c(
    "Union aging-up (log2FC > 0)",
    ">=2 Exp aging-up (log2FC > 0.05)",
    ">=2 Exp aging-up (log2FC > 0.1)",
    ">=2 Exp aging-up (log2FC > 0.25)"),
  Gene_Count = c(nrow(union_df), nrow(exp2_005_df),
    nrow(exp2_01_df), nrow(exp2_025_df)))
count_table

# save to an excel 
write_xlsx(list(Summary_Counts = count_table, Union_All_AgingUp = union_df,
    Replicated_2Exp_logFC005 = exp2_005_df,Replicated_2Exp_logFC01 = exp2_01_df,
    Replicated_2Exp_logFC025 = exp2_025_df),"Fig3_DVG_AgingUp_GeneLists.xlsx")

#### GSEA function ####

run_msigdb_gsea_all_collections <- function(
    hsc_results,
    species = "Mus musculus",
    background_genes = NULL,
    msig_collections = c(
      "H",
      "C2:CGP", "C2:BIOCARTA", "C2:KEGG_MEDICUS", "C2:PID", "C2:REACTOME",
      "C2:WIKIPATHWAYS", "C2:KEGG_LEGACY",
      "C3:MIRDB", "C3:MIR_LEGACY", "C3:GTRD", "C3:TFT_LEGACY",
      "C4:3CA", "C4:CGN", "C4:CM",
      "C5:BP", "C5:CC", "C5:MF", "C5:HPO",
      "C6",
      "C7:IMMUNESIGDB", "C7:VAX",
      "C8"
    ),
    rank_var = "log2FC_mean_adjusted_variance_YO",
    top_n = 6,
    title_prefix = "HSC"
) {
  suppressPackageStartupMessages({
    library(msigdbr)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Mm.eg.db)
    library(ggplot2)
    library(patchwork)
    library(ggplotify)
    library(stringr)
    library(dplyr)
  })
  # sanity checks
  stopifnot("gene" %in% colnames(hsc_results))
  stopifnot(rank_var %in% colnames(hsc_results))
  set.seed(42)
  # SYMBOL → ENTREZ
  gene_map <- bitr(
    unique(hsc_results$gene),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  hsc_results <- hsc_results %>%
    left_join(gene_map, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))
  # lookup table for later
  entrez2symbol <- gene_map %>%
    dplyr::select(ENTREZID, SYMBOL)
  # ranked gene list
  ranked_tbl <- hsc_results %>%
    arrange(desc(.data[[rank_var]])) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  geneList <- ranked_tbl[[rank_var]]
  names(geneList) <- ranked_tbl$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  n_genes_used <- length(geneList)
  message("Input genes: ", length(unique(hsc_results$gene)))
  message("Mapped ENTREZ: ", length(unique(entrez2symbol$ENTREZID)))
  message("Final ranked genes: ", n_genes_used)
  # convert core_enrichment to SYMBOLs
  convert_core_to_symbol <- function(core_string) {
    if (is.na(core_string) || core_string == "") return(NA_character_)
    entrez_ids <- unlist(strsplit(core_string, "/"))
    symbols <- entrez2symbol %>%
      filter(ENTREZID %in% entrez_ids) %>%
      pull(SYMBOL) %>%
      unique()
    paste(symbols, collapse = "; ")
  }
  # run one MSigDB collection
  run_one_collection <- function(msig_collection) {
    message("\n MSigDB collection: ", msig_collection)
    if (grepl(":", msig_collection)) {
      parts <- strsplit(msig_collection, ":")[[1]]
      m_df <- msigdbr(
        species = species,
        collection = parts[1],
        subcollection = parts[2]
      )
    } else {
      m_df <- msigdbr(
        species = species,
        collection = msig_collection
      )
    }
    if (nrow(m_df) == 0) {
      message("No gene sets found.")
      return(NULL)
    }
    term2gene <- m_df %>% dplyr::select(gs_name, ncbi_gene)
    term2name <- m_df %>% dplyr::select(gs_name, gs_description)
    gsea_res <- tryCatch({
      GSEA(
        geneList     = geneList,
        TERM2GENE    = term2gene,
        TERM2NAME    = term2name,
        pvalueCutoff = 1,
        verbose      = FALSE
      )
    }, error = function(e) NULL)
    if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
      return(list(
        gsea_result  = gsea_res,
        gsea_table   = NULL,
        panel_up     = NULL,
        panel_down   = NULL
      ))
    }
    # summary table with genes
    gsea_tbl <- gsea_res@result %>%
      mutate(
        collection = msig_collection,
        rank_variable = rank_var,
        condition = title_prefix,
        n_genes_used = n_genes_used,
        leading_edge_genes = vapply(
          core_enrichment,
          convert_core_to_symbol,
          FUN.VALUE = character(1)
        )
      )
    # FDR-filtered plotting
    sig_res <- gsea_res@result %>% filter(p.adjust < 0.05)
    top_up <- sig_res %>%
      filter(NES > 0) %>%
      arrange(desc(NES)) %>%
      slice_head(n = top_n)
    top_down <- sig_res %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      slice_head(n = top_n)
    make_panel <- function(df, label) {
      if (nrow(df) == 0) return(NULL)
      plots <- lapply(seq_len(nrow(df)), function(i) {
        tryCatch(
          ggplotify::as.ggplot(
            gseaplot2(
              gsea_res,
              geneSetID = df$ID[i],
              title = paste0(
                str_wrap(df$ID[i], 50),
                "\nFDR=", signif(df$p.adjust[i], 3),
                ", NES=", round(df$NES[i], 2)
              )
            )
          ),
          error = function(e) NULL
        )
      })
      plots <- plots[!sapply(plots, is.null)]
      if (length(plots) == 0) return(NULL)
      wrap_plots(plots)
    }
    list(
      gsea_result = gsea_res,
      gsea_table  = gsea_tbl,
      panel_up    = make_panel(top_up, "Up"),
      panel_down  = make_panel(top_down, "Down")
    )
  }
  # run all collections
  results <- lapply(msig_collections, run_one_collection)
  names(results) <- msig_collections
  # combine plots
  combined_up_panel <- wrap_plots(
    lapply(results, `[[`, "panel_up")[!sapply(lapply(results, `[[`, "panel_up"), is.null)]
  )
  combined_down_panel <- wrap_plots(
    lapply(results, `[[`, "panel_down")[!sapply(lapply(results, `[[`, "panel_down"), is.null)]
  )
  # combine summary tables
  gsea_summary_table <- bind_rows(
    lapply(results, `[[`, "gsea_table")
  ) %>%
    arrange(collection, desc(NES))
  
  return(list(
    per_collection      = results,
    combined_up_panel   = combined_up_panel,
    combined_down_panel = combined_down_panel,
    gsea_summary_table  = gsea_summary_table
  ))
}

gsea_Day0_all <- run_msigdb_gsea_all_collections(
  hsc_results  = Day0_All,
  rank_var     = "log2FC_mean_adjusted_variance_OY",
  title_prefix = "Day0_All_Genes"
)
Day0_GSEA_results <- gsea_Day0_all$gsea_summary_table

out_dir <- ".../GSEA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.xlsx(gsea_Day0_all$gsea_summary_table, file = file.path(out_dir, "GSEA_Day0_all_genes_Enriched_Pathways.xlsx"), overwrite = TRUE)

######################################################

### Run GSEA analysis using shared DVGs (Dayu0 and Unmanipulate) as a validation ###

shared_genes_Day0_Un <- intersect(unique(Day0_All$gene), unique(Unmanipulated_All$gene))
Day0_shared_filtered <- Day0_All %>% filter(gene %in% shared_genes_Day0_Un)
Day0_shared_filtered # 191 genes 

gsea_Day0_shared_validation <- run_msigdb_gsea_all_collections(
  hsc_results  = Day0_shared_filtered,
  rank_var     = "log2FC_mean_adjusted_variance_OY",
  title_prefix = "Day0_SharedGene_Validation"
)
Day0_Unmani_shared_genes_GSEA_results <- gsea_Day0_shared_validation$gsea_summary_table
View(Day0_Unmani_shared_genes_GSEA_results)

# Save the excel
out_dir <- ".../GSEA"
write.xlsx(gsea_Day0_shared_validation$gsea_summary_table, file = file.path(out_dir, "GSEA_Day0_Unmani_shared_genes_Enriched_Pathways.xlsx"), overwrite = TRUE)


### Meta GSEA using the average fold change as ranking variable ### 

shared_genes <- intersect(Day0_All$gene, Unmanipulated_All$gene)
length(shared_genes) # 202 genes found in both experiments 
length(shared_genes) / length(Day0_All$gene) # 24%
length(shared_genes) / length(Unmanipulated_All$gene) # 28%

day0_fc <- Day0_All %>% dplyr::select(gene, log2FC_mean_adjusted_variance_YO) %>% rename(log2FC_Day0 = log2FC_mean_adjusted_variance_YO)
unm_fc <- Unmanipulated_All %>% dplyr::select(gene, log2FC_mean_adjusted_variance_YO) %>% rename(log2FC_Unmanipulated = log2FC_mean_adjusted_variance_YO)
meta_gsea_input <- day0_fc %>% filter(gene %in% shared_genes) %>% inner_join(unm_fc %>% filter(gene %in% shared_genes), by = "gene")
nrow(meta_gsea_input)
# compute the combined fold change 
meta_gsea_input <- meta_gsea_input %>% mutate(log2FC_meta_mean = (log2FC_Day0 + log2FC_Unmanipulated) / 2)
cor(meta_gsea_input$log2FC_Day0,meta_gsea_input$log2FC_Unmanipulated,method = "spearman") # very weak correlation

gsea_meta_Day0_Unman_combined <- run_msigdb_gsea_all_collections(
  hsc_results  = meta_gsea_input,
  rank_var     = "log2FC_meta_mean",
  title_prefix = "Meta_Day0_Unmanipulated"
)
Day0_Unmani_shared_genes_combined_GSEA_results <- gsea_meta_Day0_Unman_combined$gsea_summary_table
Day0_Unmani_shared_genes_combined_GSEA_results

#######################################################

# Read in the two GSEA tables 


Day0_GSEA <- '.../GSEA_Day0_all_genes_Enriched_Pathways.xlsx'
Day0_GSEA_results <- read_excel(Day0_GSEA, sheet = 1)
Day0_GSEA_results # 8800 rows 

Day0_Unmani_shared_genes_GSEA <- '.../GSEA_Day0_Unmani_shared_genes_Enriched_Pathways.xlsx'
Day0_Unmani_shared_genes_GSEA_results <- read_excel(Day0_Unmani_shared_genes_GSEA, sheet = 1)
Day0_Unmani_shared_genes_GSEA_results # 1338 rows 


#######################################################

## GSEA plots ##

Day0_GSEA_results
colnames(Day0_GSEA_results)
# Related pathways: 
# HALLMARK_E2F_TARGETS:E2F cell-cycle program
# FISCHER_DREAM_TARGETS: DREAM cell-cycle program 

# CROONQUIST_IL6_DEPRIVATION_DN: IL6 signaling 
# CROONQUIST_NRAS_SIGNALING_DN: NRAS signaling 
# REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3: Rho GTPase signaling
# PID_CMYB_PATHWAY: c-MYB regulatory program
# KEGG_P53_SIGNALING_PATHWAY: p53 stress-response pathway
# MARTENS_TRETINOIN_RESPONSE_DN: ATRA response to APL program
# GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP: Quiescent leukemic stem cell program

# MORF_NPM1 # GSEA PLOT : NPM1 regulatory program

Day0_Unmani_shared_genes_GSEA_results

# HALLMARK_E2F_TARGETS : E2F cell-cycle program
# HALLMARK_G2M_CHECKPOINT: G2/M cell-cycle checkpoint
# REACTOME_INFECTIOUS_DISEASE: Reactome infectious disease 
# KIM_WT1_TARGETS_DN: WT1 regulatory program
# BYSTRYKH_HEMATOPOIESIS_STEM_CELL_QTL_TRANS: Hematopoietic stem cell trans-QTL
# HALLMARK_APOPTOSIS: Apoptotic signaling program
# REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM: Cytokine signaling pathway

# # HALLMARK_E2F_TARGETS : E2F cell-cycle program # GSEA plot 


# GSEA plots # 
library(enrichplot)
library(ggplotify)

# get the collection

# Day0_GSEA_results

plot_gsea_pathway <- function(
    pathway_id,
    gsea_obj,
    gsea_summary,
    color = "green",
    base_size = 18
){
  library(dplyr)
  library(enrichplot)
  library(aplot)
  library(ggplot2)
  stats <- gsea_summary %>%
    dplyr::filter(ID == pathway_id) %>%
    dplyr::slice(1)
  if(nrow(stats) == 0){
    stop(paste("Pathway not found:", pathway_id))
  }
  collection <- stats$collection
  gsea_res <- gsea_obj$per_collection[[collection]]$gsea_result
  # Generate GSEA plot
  p <- gseaplot2(
    x = gsea_res,
    geneSetID = pathway_id,
    title = paste0(
      pathway_id,
      "\nqval=", signif(stats$qvalue,3),
      " NES=", round(stats$NES,2)
    ),
    color = color,
    base_size = base_size,
    rel_heights = c(1.5,0.5,1),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom = "line"
  )
  rm_left <- theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
  p1 <- p[[1]] + rm_left + theme(plot.title = element_text(size = 22))
  p2 <- p[[2]] + rm_left
  p3 <- p[[3]] + rm_left
  plot_clean <- aplot::plot_list(
    p1, p2, p3,
    ncol = 1,
    heights = c(1.5,0.5,1)
  )
  return(list(
    pathway = pathway_id,
    stats = stats,
    plot = plot_clean
  ))
}

### Plot GSEA plot using the function ###

pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "FISCHER_DREAM_TARGETS",
  "CROONQUIST_IL6_DEPRIVATION_DN",
  "CROONQUIST_NRAS_SIGNALING_DN",
  "REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3",
  "PID_CMYB_PATHWAY",
  "KEGG_P53_SIGNALING_PATHWAY",
  "MARTENS_TRETINOIN_RESPONSE_DN",
  "GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP",
  "MORF_NPM1")
Day0_GSEA_plots <- lapply(
  pathways,
  plot_gsea_pathway,
  gsea_obj = gsea_Day0_all,
  gsea_summary = gsea_Day0_all$gsea_summary_table)
names(Day0_GSEA_plots) <- pathways
Day0_GSEA_plots$MORF_NPM1 # NPM1
Day0_GSEA_plots$HALLMARK_E2F_TARGETS # E2F targets
Day0_GSEA_plots$FISCHER_DREAM_TARGETS # FREAM targets
Day0_GSEA_plots$CROONQUIST_IL6_DEPRIVATION_DN
Day0_GSEA_plots$CROONQUIST_NRAS_SIGNALING_DN
Day0_GSEA_plots$REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3
Day0_GSEA_plots$PID_CMYB_PATHWAY
Day0_GSEA_plots$KEGG_P53_SIGNALING_PATHWAY
Day0_GSEA_plots$MARTENS_TRETINOIN_RESPONSE_DN
Day0_GSEA_plots$GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP

######################

### Dot plot ###

pathways_keep_day0 <- c(
  "HALLMARK_E2F_TARGETS",
  "FISCHER_DREAM_TARGETS",
  "CROONQUIST_IL6_DEPRIVATION_DN",
  "CROONQUIST_NRAS_SIGNALING_DN",
  "REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3",
  "PID_CMYB_PATHWAY",
  "KEGG_P53_SIGNALING_PATHWAY",
  "MARTENS_TRETINOIN_RESPONSE_DN",
  "GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP"
)
gsea_plot_df_day0 <- gsea_Day0_all$gsea_summary_table %>%
  dplyr::filter(ID %in% pathways_keep_day0) %>%
  dplyr::mutate(
    neg_log10_FDR = -log10(p.adjust))

pathway_labels <- c(
  "HALLMARK_E2F_TARGETS" = "E2F cell-cycle program",
  "FISCHER_DREAM_TARGETS" = "DREAM cell-cycle program",
  "CROONQUIST_IL6_DEPRIVATION_DN" = "IL6 signaling",
  "CROONQUIST_NRAS_SIGNALING_DN" = "NRAS signaling",
  "REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3" = "Rho GTPase signaling",
  "PID_CMYB_PATHWAY" = "c-MYB regulatory program",
  "KEGG_P53_SIGNALING_PATHWAY" = "p53 stress-response pathway",
  "MARTENS_TRETINOIN_RESPONSE_DN" = "ATRA response to APL program",
  "GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP" = "Quiescent leukemic stem cell program"
)

gsea_plot_df_day0 <- gsea_plot_df_day0 %>%
  dplyr::mutate(
    pathway_readable = pathway_labels[ID]
  ) %>%
  dplyr::arrange(NES)

Day0_GSEA_dotplot <- ggplot(
  gsea_plot_df_day0,
  aes(
    x = NES,
    y = reorder(pathway_readable, NES),
    size = setSize,
    fill = qvalue
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 1,
    alpha = 0.95
  ) +
  scale_fill_gradient(
    high = "#D73027",  
    low = "#4575B4",  
    name = "qvalue",
    trans = "reverse"
  ) +
  scale_size_continuous(
    range = c(6, 15),
    breaks = c(25, 50, 100),
    name = "Gene set size"
  ) +
  labs(
    x = "NES",
    y = NULL,
    title = "Day0 DVG GSEA"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text  = element_text(color = "black"),
    plot.margin = margin(5, 25, 5, 5)
  )

Day0_GSEA_dotplot

## Dot plot for Shared Gene GSEA (Day0 and Unmanipulated) ##

pathways_keep_shared <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "REACTOME_INFECTIOUS_DISEASE",
  "KIM_WT1_TARGETS_DN",
  "BYSTRYKH_HEMATOPOIESIS_STEM_CELL_QTL_TRANS",
  "HALLMARK_APOPTOSIS",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"
)

pathway_labels_shared <- c(
  "HALLMARK_E2F_TARGETS" = "E2F cell-cycle program",
  "HALLMARK_G2M_CHECKPOINT" = "G2/M cell-cycle checkpoint",
  "REACTOME_INFECTIOUS_DISEASE" = "Infectious disease response",
  "KIM_WT1_TARGETS_DN" = "WT1 regulatory program",
  "BYSTRYKH_HEMATOPOIESIS_STEM_CELL_QTL_TRANS" = "Hematopoietic stem cell trans-QTL",
  "HALLMARK_APOPTOSIS" = "Apoptotic signaling program",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" = "Cytokine signaling pathway"
)
gsea_plot_df_shared <- gsea_Day0_shared_validation$gsea_summary_table %>%
  dplyr::filter(ID %in% pathways_keep_shared) %>%
  dplyr::mutate(
    pathway_readable = pathway_labels_shared[ID]
  ) %>%
  dplyr::arrange(NES)


Day0_Unmani_shared_GSEA_dotplot <- ggplot(
  gsea_plot_df_shared,
  aes(
    x = NES,
    y = reorder(pathway_readable, NES),
    size = setSize,
    fill = qvalue
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 1,       
    alpha = 0.95
  ) +
  scale_fill_gradient(
    high = "#D73027",    
    low = "#4575B4",    
    name = "qvalue",
    trans = "reverse",
    # breaks = c(0.05, 0.1, 0.15),
    limits = range(gsea_plot_df_shared$qvalue)
  ) +
  scale_size_continuous(
    range = c(6, 15),
    breaks = c(10, 15, 20),
    name = "Gene set size"
  ) +
  labs(
    x = "NES",
    y = NULL,
    title = "Day0 Unmanipulated Shared DVGs GSEA"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text  = element_text(color = "black"),
    plot.margin = margin(5, 25, 5, 5)
  ) +
  guides(
    fill = guide_colorbar(order = 1), 
    size = guide_legend(order = 2)  
  )

Day0_Unmani_shared_GSEA_dotplot

##################################################


### End of this script 


















