#######################################################################################


# Compute Output Activity (OA) for CrossAge, Exp1 and Exp2 vitro data


#######################################################################################

# Maintainer: Chris Chen 
# Last updated: 02/03/2026

##########################

setwd('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts')

# save in 

# /project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/test_folder/Updated_main_figures_0223


## Load packages ##

load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","scMayoMap",
    "readxl","pheatmap","Matrix","openxlsx","gt","glue"
  )
  
  suppressPackageStartupMessages(
    lapply(pkgs, require, character.only = TRUE)
  )
  message("Allpackages loaded.")
}
load_all_packages()

# set working directory
setwd('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure1/OA_pipeline_out/')

################################

### Define low and high output markers ###

# From Rodriguez-Fraticelli et al., Nature 2020 and Singh et al., Cell Stem Cell, 2025

low_output_markers <- c(
  "Mpl","Ifitm3","Ifitm1","Tgm2","H2-K1","Socs2","Mycn","Nupr1","Hacd4",
  "Mllt3","Gda","B2m","Procr","Txnip","Clu","Sult1a1","S100a6","Rpl36a",
  "Tsc22d1","Gng11","Ccnd2","mt-Co3","Tmem176b","Lmo2","Ly6a","Mmrn1",
  "Trim47","St3gal1","Mecom","Pik3r1","Adgrg1","Esam","Ryk","Rps21","Hoxb8",
  "Ccnd1","Uba7","Rps28","Serpina3g","Rpl37a","Fkbp1a","Pdzk1ip1","Selp",
  "Eif4a2","Tmem176a","Bex4","Grb10","Jam3","Ptk2b","Gimap8","Csgalnact1",
  "H2-Q7","App","Mylk","Casp12","Rhof","Laptm4a","Tcf15","Gabarapl1","Ctsl",
  "Bex1","Rpl21","D630039A03Rik","Myof","Glul","Ppic","Tpt1","Rpl5","Cdkn1c",
  "mt-Nd5","Abcg3","Aplp2","Art4","Tbxas1","Cish","Vwf","Scarf1","mt-Atp6",
  "H2-Eb1","H2-D1","Cd74","Iigp1","Aldh1a1","Gstm1","mt-Nd4","mt-Cytb","mt-Co1"
)
high_output_markers <- c(
  "Plac8","H2afy","Cdk6","Cd34","Nkg7","Ptma","Mpo","Cd48","Stmn1","Fam117a",
  "Slc22a3","Adgrg3","Ppia","Car1","Ctsg","Flt3","Lgals1","Muc4","Gpx1",
  "Hmgb2","Ndufa4","Serpinb1a","Ccl9","Oaz1","H2afz","Crip1","Mcm7","Cpa3",
  "Vim","Ybx1","Sell","Sh3bgrl3","H3f3a","Dut","Atpif1","Ran","Hnrnpa2b1",
  "Hdgf","Mcm4","Elane","Rabgap1l","Cmtm7","Rpsa","Mcm6","Plek","Set","Atp5g3",
  "Myc","Taldo1","Tuba1b","Cks2","Slc25a5","Fam65a","Cebpa","Tmsb10","Cd52",
  "Klf1","Anp32b","Hn1","Parvg","Ffar2","Bex6","Emilin2","Itgal","Cox5a",
  "Hnrnpab","Ighm","BC035044","Lmnb1","Golm1","Bin1","Igfbp4","Tyrobp","E2f8",
  "Banf1","Cst7","Sh2d5","Aqp1","Ptbp3","Snrpf","Rfc2","Ramp1","Hmgn5","Nme1"
)

# All genes detetced in Rodriguez-Fraticelli et al. (background genes for hypergeometric test)

cellstem_path <- "/project2/sli68423_1316/users/chris/Celltag/HSC_heterogenecity/DEG_clonal_behaviors/Figure1_low_high_OA_Exp1_Exp2/Cell_Stem_Cell_and_nature_DEG_list/cell_stem_cell_deg_list.xlsx"

cellstem_df <- openxlsx::read.xlsx(cellstem_path, sheet = 2)

head(cellstem_df)
cellstem_genes <-  unique(cellstem_df$gene) # 4789 genes 

###############################################################################

# Read in CrossAge, Unmanipulated, Exp1 and Exp2 data 

#### Read in CrossAge data ####

seurat_in_vivo_day60 <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)

seurat_in_vitro_day0 <- readRDS("/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vitro.RDS")

#### Read in Exp1 and Exp2 data ####

# In vitro 

# seurat_vitro_exp1 <- readRDS('/project2/sli68423_1316/users/chris/Celltag/HSC_heterogenecity/DEG_clonal_behaviors/Figure1_low_high_OA_Exp1_Exp2/seurat_vitro_exp1.rds')

seurat_vitro_exp2 <- readRDS('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/Exp2_vitro_final.RDS')

# In vivo 

# seurat_vivo_exp1 <- readRDS('/project2/sli68423_1316/users/chris/Celltag/HSC_heterogenecity/DEG_clonal_behaviors/Figure1_low_high_OA_Exp1_Exp2/seurat_vivo_exp1.rds')

seurat_vivo_exp2 <- readRDS('//project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/Exp2_vivo_final.RDS')

######################################################

# Define the cell type color 

celltype_colors <- c(
  ## HSCs 
  "LT-HSC" = "#E76F51",
  "ST-HSC" = "#F4A261",
  ## Early progenitors 
  "MPP"  = "#5BC374",
  "MPP3" = "#4CAF50",
  "MPP4" = "#66CDAA",
  "CMP"  = "#C49A6C",
  "EryP" = "#F08080",
  "LMPP" = "#4DB6AC",
  ## Myeloid differentiation
  "GMP"   = "#1FA187",
  "GMP-3" = "#2BBBAD",
  "MDP"   = "#008B8B",
  "Granulocyte" = "#6C8EBF",
  ## Megakaryocyte / Platelet
  "MkP"   = "#5FBAC8",
  "MKP"   = "#4DA8B5",
  "MKP-1" = "#78C2D0",
  ## Erythroid / Megakaryocyte
  "MEP"     = "#8C973E",
  "MEP-1"   = "#7A8735",
  "MEP-2"   = "#9AAE4C",
  "MEP/EryP"= "#A4B95A",
  ## Mast / Basophil 
  "Ba/mast"   = "#C77DFF",
  "Mast-cell" = "#B565D9",
  ## Immune / Other 
  "DC"  = "#9E6E4A",
  "Mac" = "#7F5539",
  "Mast-Cells" = "#D98BB8",
  ## T cell 
  "T-cells" = "#6CA3D1",
  "T-cell"  = "#4F8AC9",
  "T-Cell"  = "#83B8E6",
  "T-Cells"  = "#83B8E6",
  ## Undefined
  "UN" = "grey65")


celltype_colors <- c(
  
  # Stem / early
  "HSC"  = "#E76F51",   # deep red
  "MPP"  = "#F4A261",   # teal-green
  "CMP"  = "#8C973E",   # warm brown
  
  # Myeloid / progenitors
  "GMP"  = "#007F5F",   # dark green
  "MEP"  = "#6A994E",   # olive green
  "MkP"  = "#277DA1",   # steel blue
  "EryP" = "#F08080",   # strong red
  
  # Immune / differentiated
  "DC"           = "#9C6644",  # brown
  "Mac"          = "#6D4C41",  # dark brown
  "Granulocyte"  = "#6C8EBF",  # strong blue
  
  # Lymphoid
  "B_cell" = "#B565D9",  # purple
  "T_cell" = "#4361EE",  # royal blue
  
  # Mast
  "Mast" = "#D98BB8"     # magenta
)


## Define the clone wise DEG function ##

run_clonewise_DEG_suite <- function(
    seurat_obj,
    clone_set1_ids,
    clone_set2_ids,
    clone_col   = "CloneID",
    assay_name  = "RNA",
    top_hvgs    = 2000,
    min_detect_prop = 0.10,
    min_mean_norm  = 0.10,
    min_mean_raw   = 2,
    excel_name  = "Clonewise_DEG.xlsx",
    volcano_name = "Volcano_MAST_TopHVGs.png",
    group1_label = "Low_OA",
    group2_label = "High_OA",
    fc_cutoff  = 0.25,
    fdr_cutoff = 0.05,
    output_dir = "Clonewise_DEG_Results"
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(glue)
    library(openxlsx)
    library(ggplot2)
    library(ggrepel)
  })
  DefaultAssay(seurat_obj) <- assay_name
  # ============================================================
  # ️ Create output directory for this DEG run
  # ============================================================
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_folder <- file.path(output_dir,
                          glue("{group1_label}_vs_{group2_label}_{timestamp}"))
  dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)
  excel_path   <- file.path(run_folder, excel_name)
  volcano_path <- file.path(run_folder, volcano_name)
  message(glue("📁 Output directory created: {run_folder}"))
  # ============================================================
  # 0️⃣ Assign group labels based on CloneID
  # ============================================================
  seurat_obj$group_label <- NA_character_
  seurat_obj$group_label[seurat_obj[[clone_col]][, 1] %in% clone_set1_ids] <- "Group1"
  seurat_obj$group_label[seurat_obj[[clone_col]][, 1] %in% clone_set2_ids] <- "Group2"
  seurat_sub <- subset(seurat_obj, subset = !is.na(group_label))
  Idents(seurat_sub) <- seurat_sub$group_label
  n1 <- sum(seurat_sub$group_label == "Group1")
  n2 <- sum(seurat_sub$group_label == "Group2")
  message(glue("🧬 Clone-wise DEG: {n1} {group1_label} cells vs {n2} {group2_label} cells"))
  # ============================================================
  # 1️⃣ Helper: adjust FDR & append avg expression
  # ============================================================
  adjust_and_append_avg <- function(deg_df, obj, features) {
    if (is.null(deg_df) || nrow(deg_df) == 0) return(tibble())
    deg_df <- deg_df %>%
      tibble::rownames_to_column("gene") %>%
      mutate(
        p_val_adj = p.adjust(p_val, method = "BH"),
        FDR_BY_manual = p.adjust(p_val, method = "BY")
      )
    ae_raw <- AverageExpression(obj, group.by = "group_label", features = features)
    ae_df <- if (is.list(ae_raw) && "RNA" %in% names(ae_raw)) as.data.frame(ae_raw$RNA) else as.data.frame(ae_raw)
    ae_df <- ae_df %>% tibble::rownames_to_column("gene")
    col_g1 <- grep("Group1", colnames(ae_df), value = TRUE)[1]
    col_g2 <- grep("Group2", colnames(ae_df), value = TRUE)[1]
    ae_df <- ae_df %>%
      rename_with(~ "avg_Group1", all_of(col_g1)) %>%
      rename_with(~ "avg_Group2", all_of(col_g2))
    left_join(deg_df, ae_df, by = "gene")
  }
  # ============================================================
  # 2️⃣ Helper: run one DEG method
  # ============================================================
  do_one_deg <- function(obj, features, method = "wilcox") {
    if (length(features) == 0) return(tibble())
    FindMarkers(
      obj,
      ident.1 = "Group1",
      ident.2 = "Group2",
      features = features,
      test.use = method,
      logfc.threshold = 0,
      min.pct = 0
    )
  }
  # ============================================================
  # 3️⃣ Gene universes
  # ============================================================
  all_genes <- rownames(seurat_sub)
  seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst",
                                     nfeatures = top_hvgs, verbose = FALSE)
  hvgs <- VariableFeatures(seurat_sub)
  data_mat <- GetAssayData(seurat_sub, slot = "data")[hvgs, , drop = FALSE]
  det_prop <- Matrix::rowMeans(data_mat > 0)
  mean_norm <- Matrix::rowMeans(data_mat)
  raw_mat <- tryCatch(GetAssayData(seurat_sub, slot = "counts")[hvgs, , drop = FALSE],
                      error = function(e) NULL)
  mean_raw <- if (!is.null(raw_mat)) Matrix::rowMeans(raw_mat) else rep(0, length(hvgs))
  names(mean_raw) <- hvgs
  filt_genes <- hvgs[
    det_prop[hvgs] >= min_detect_prop &
      mean_norm[hvgs] >= min_mean_norm &
      mean_raw[hvgs] >= min_mean_raw
  ]
  # ============================================================
  # 4️⃣ DEG results loop (3 gene sets × 2 methods)
  # ============================================================
  combos <- tribble(
    ~setting, ~features, ~label,
    "DEG1_AllGenes",  all_genes,       "All genes",
    "DEG2_TopHVGs",   hvgs,            glue("Top {length(hvgs)} HVGs"),
    "DEG3_FilteredHVGs", filt_genes,   glue("{length(filt_genes)} filtered HVGs")
  )
  methods <- c("wilcox", "MAST")
  deg_results <- list()
  summary_rows <- list()
  for (i in seq_len(nrow(combos))) {
    setting <- combos$setting[i]
    feats   <- combos$features[[i]]
    label   <- combos$label[i]
    for (m in methods) {
      message(glue("🔹 Running {m} on {setting} ({length(feats)} genes)"))
      res <- do_one_deg(seurat_sub, feats, m)
      res2 <- adjust_and_append_avg(res, seurat_sub, feats)
      key <- glue("{setting}_{toupper(m)}")
      deg_results[[key]] <- res2
      n_sig <- if (nrow(res2)) sum(res2$p_val_adj < fdr_cutoff) else 0
      n_up1 <- if (nrow(res2)) sum(res2$avg_log2FC >  fc_cutoff & res2$p_val_adj < fdr_cutoff) else 0
      n_up2 <- if (nrow(res2)) sum(res2$avg_log2FC < -fc_cutoff & res2$p_val_adj < fdr_cutoff) else 0
      summary_rows[[key]] <- tibble(
        Setting = setting,
        Method  = toupper(m),
        Description = label,
        Cells_Group1 = n1,
        Cells_Group2 = n2,
        Genes_Test = length(feats),
        Sig_FDR = n_sig,
        Up_in_Group1 = n_up1,
        Up_in_Group2 = n_up2
      )
    }
  }
  summary_tab <- bind_rows(summary_rows)
  # ============================================================
  # 5️⃣ Save Excel
  # ============================================================
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_tab)
  
  write_one <- function(name, df) {
    addWorksheet(wb, name)
    if (is.null(df) || nrow(df) == 0) {
      writeData(wb, name, data.frame(note = "No results."))
    } else {
      writeData(wb, name, df)
    }
  }
  purrr::iwalk(deg_results, ~ write_one(.y, .x))
  saveWorkbook(wb, excel_path, overwrite = TRUE)
  message(glue("💾 Excel saved: {excel_path}"))
  # ============================================================
  # 6️⃣ Volcano plot (MAST × Top HVGs)
  # ============================================================
  mast_key <- "DEG2_TopHVGs_MAST"
  volc_plot <- NULL
  mast_df <- deg_results[[mast_key]]
  if (!is.null(mast_df) && nrow(mast_df) > 0) {
    df <- mast_df %>%
      mutate(
        neg_log10_fdr = -log10(pmax(p_val_adj, 1e-300)),
        significance = case_when(
          p_val_adj < fdr_cutoff & avg_log2FC >  fc_cutoff ~ "↑Group1",
          p_val_adj < fdr_cutoff & avg_log2FC < -fc_cutoff ~ "↑Group2",
          TRUE ~ "Not significant"
        )
      )
    top_left <- df %>%
      filter(significance == "↑Group2") %>%
      arrange(avg_log2FC, p_val_adj) %>%  
      slice_head(n = 10)
    top_right <- df %>%
      filter(significance == "↑Group1") %>%
      arrange(desc(avg_log2FC), p_val_adj) %>%  
      slice_head(n = 10)
    top_genes <- bind_rows(top_left, top_right)
    top_genes$gene <- paste0("italic('", top_genes$gene, "')")
    # Color palette
    volcano_colors <- c("#5271AE", "#D85B59", "grey80")
    names(volcano_colors) <- c("↑Group2", "↑Group1", "Not significant")
    # Significance threshold
    sig_line <- -log10(fdr_cutoff)
    # Padding
    y_min <- sig_line - 10
    y_max <- max(df$neg_log10_fdr)
    
    volc_plot <- ggplot(df,
                        aes(x = avg_log2FC,
                            y = neg_log10_fdr)) +
      geom_point(
        aes(color = significance),
        alpha = 0.85,
        size = 3
      ) +
      geom_text_repel(
        data = top_genes,
        aes(label = gene),
        parse = TRUE,
        size = 5,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.4,
        force = 2,
        segment.size = 0.3
      ) +
      scale_color_manual(values = volcano_colors) +
      geom_vline(
        xintercept = c(-fc_cutoff, fc_cutoff),
        linetype = "dashed"
      ) +
      geom_hline(
        yintercept = sig_line,
        linetype = "dashed"
      ) +
      coord_cartesian(
        ylim = c(y_min, y_max),
        clip = "off"
      ) +
      xlim(-4, 4) +
      theme_classic(base_size = 20) +
      theme(
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size = 20, face = "bold", color = "black"),
        axis.text  = element_text(size = 18, color = "black"),
        legend.position = "none"
      ) +
      labs(
        title = glue("Clone-wise DEG Volcano (MAST, Top {top_hvgs} HVGs)"),
        x = expression(log[2]("FC (Group1 / Group2)")),
        y = expression(-log[10]("FDR"))
      )
    
    ggsave(volcano_path, volc_plot, width = 8.5, height = 7, dpi = 300)
    
    message(glue("🌋 Volcano saved: {volcano_path}"))
  }
  return(invisible(list(
    summary = summary_tab,
    results = deg_results,
    volcano_plot = volc_plot,
    folder = run_folder
  )))
}


## OA function (per-sample input)

run_low_high_OA_analysis_for_a_single_sample <- function(
    seurat_obj,
    sample_label = "Sample",
    output_dir = "OA_Analysis",
    clone_col = "CloneID",
    HSC_types = c("HSC"),
    nonHSC_types = NULL,
    deg_celltypes = c("HSC","MPP"),
    reduction = "umap",
    top_hvgs = 2000,
    min_detect_prop = 0.10,
    min_mean_norm  = 0.10,
    min_mean_raw   = 2,
    fc_cutoff      = 0.25,
    fdr_cutoff     = 0.05,
    low_output_markers,
    high_output_markers,
    celltype_colors,
    full_ref_deg_genes_list
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggrepel)
    library(glue)
    library(openxlsx)
    library(ggvenn)
    library(ggrastr)
    library(scales)
  })
  message("\n==============================")
  message(glue("🔷 Starting OA Analysis: {sample_label}"))
  message("==============================\n")
  # -----------------------------------------
  # 0. Create output folder
  # -----------------------------------------
  sample_dir <- file.path(output_dir, sample_label)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  message(glue("📁 Output folder created: {sample_dir}\n"))
  # -----------------------------------------
  # 1. DimPlot
  # -----------------------------------------
  message("🔹 Step 1: Creating DimPlot …")
  p_dim <- DimPlot(
    seurat_obj,
    reduction = reduction,
    cols      = celltype_colors,
    pt.size   = 0.4,
    label     = TRUE,
    label.size = 6
  ) +
    ggtitle(glue("{sample_label} — Celltype UMAP")) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text  = element_text(size = 16)
    )
  ggsave(
    filename = file.path(sample_dir, glue("{sample_label}_DimPlot.pdf")),
    plot = p_dim, width = 8.5, height = 7.5
  )
  message("Done: DimPlot saved.\n")
  # -----------------------------------------
  # 2. Compute OA per clone
  # -----------------------------------------
  message("🔹 Step 2: Computing Output Activity (OA)…")
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::rename(CloneID = !!clone_col)
  # define non-HSC types (user-specified or automatic)
  if (is.null(nonHSC_types)) {
    nonHSC_types <- setdiff(unique(meta$celltype), HSC_types)
  }
  message(glue("   HSC types: {paste(HSC_types, collapse=', ')}"))
  message(glue("   non-HSC types: {paste(nonHSC_types, collapse=', ')}"))
  # remove cloneID 0 or NA
  meta <- meta %>% filter(CloneID != "0", !is.na(CloneID))
  # clone counts
  clone_counts <- meta %>%
    mutate(group = ifelse(celltype %in% HSC_types, "HSC", "nonHSC")) %>%
    dplyr::count(CloneID, group) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = n,
      values_fill = 0
    )
  clone_freq <- clone_counts %>%
    mutate(
      total_HSC_all = sum(HSC),
      total_nonHSC_all = sum(nonHSC),
      HSC_freq = HSC / total_HSC_all,
      nonHSC_freq = nonHSC / total_nonHSC_all,
      OA = (nonHSC_freq + 1e-6) / (HSC_freq + 1e-6),
      log2OA = log2(OA)
    )
  # threshold low / high
  q_low  <- quantile(clone_freq$OA, 0.30)
  q_high <- quantile(clone_freq$OA, 0.70)
  low_clones  <- clone_freq %>% filter(OA <= q_low)  %>% pull(CloneID)
  high_clones <- clone_freq %>% filter(OA >= q_high) %>% pull(CloneID)
  message(glue("   - Low Output clones: {length(low_clones)}"))
  message(glue("   - High Output clones: {length(high_clones)}"))
  # Save OA tables
  write.xlsx(
    clone_freq,
    file = file.path(sample_dir, glue("{sample_label}_OA_by_clone.xlsx")),
    overwrite = TRUE
  )
  cell_OA_df <- meta %>%
    left_join(clone_freq %>% dplyr::select(CloneID, OA), by = "CloneID") %>%
    dplyr::select(cell_ID, CloneID, OA)
  write.xlsx(
    cell_OA_df,
    file = file.path(sample_dir, glue("{sample_label}_OA_by_cell.xlsx")),
    overwrite = TRUE
  )
  message("Done: OA tables saved.\n")
  # -----------------------------------------
  # 3. OA UMAP plot
  # -----------------------------------------
  message("🔹 Step 3: Creating OA UMAP …")
  emb <- Embeddings(seurat_obj, reduction)
  
  umap_df <- emb[, 1:2] %>%
    as.data.frame() %>%
    setNames(c("UMAP_1", "UMAP_2")) %>%
    tibble::rownames_to_column("cell_ID")
  meta_umap <- meta %>%
    left_join(umap_df, by = "cell_ID") %>%
    left_join(clone_freq %>% dplyr::select(CloneID, OA), by = "CloneID") %>%
    filter(!is.na(OA))
  meta_umap$OA <- squish(meta_umap$OA, c(0, 2))
  cell_centroids <- meta_umap %>%
    group_by(celltype) %>%
    summarize(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  p_oa <- ggplot(meta_umap, aes(UMAP_1, UMAP_2, color = OA)) +
    geom_point(size = 1.4, alpha = 0.55) +
    scale_color_gradientn(
      colors = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
      values = scales::rescale(c(0, 0.7, 1, 1.3, 2)),
      limits = c(0, 2),
      name = "OA"
    ) +
    geom_text(
      data = cell_centroids,
      aes(label = celltype),
      size = 5.5, fontface = "bold", color = "grey10"
    ) +
    coord_equal() +
    theme_void(base_size = 16) +
    ggtitle(glue("{sample_label} — OA (Ai) UMAP"))
  ggsave(
    filename = file.path(sample_dir, glue("{sample_label}_OA_UMAP.pdf")),
    plot = p_oa, width = 6.5, height = 5.5
  )
  message("Done: OA UMAP saved.\n")
  # -----------------------------------------
  # 4. DEG (subset for HSC / MPP)
  # -----------------------------------------
  message("🔹 Step 4: Running clonewise DEG (HSC-only)…")
  # Subset object for DEG cell types
  seurat_deg <- subset(seurat_obj, subset = celltype %in% deg_celltypes)
  # Run DEG module
  deg_out <- run_clonewise_DEG_suite(
    seurat_obj       = seurat_deg,
    clone_set1_ids   = low_clones,
    clone_set2_ids   = high_clones,
    clone_col        = clone_col,
    group1_label     = "Low_output",
    group2_label     = "High_output",
    top_hvgs         = top_hvgs,
    min_detect_prop  = min_detect_prop,
    min_mean_norm    = min_mean_norm,
    min_mean_raw     = min_mean_raw,
    fc_cutoff        = fc_cutoff,
    fdr_cutoff       = fdr_cutoff,
    excel_name       = glue("{sample_label}_DEG.xlsx"),
    volcano_name     = glue("{sample_label}_volcano.pdf"),
    output_dir       = sample_dir
  )
  message("Done: DEG analysis completed.\n")
  message("   Notes: Adding reference-marker volcano plot ...")
  # Ensure plots list exists
  if (!exists("plots")) plots <- list()
  # -----------------------------
  # Extract DEG2 Top HVGs (MAST)
  # -----------------------------
  mast_top <- deg_out$results$DEG2_TopHVGs_MAST
  if (!is.null(mast_top) && nrow(mast_top) > 0) {
    # Combined reference markers
    reference_genes <- unique(c(low_output_markers, high_output_markers))
    df_ref <- mast_top %>%
      mutate(
        neg_log10_fdr = -log10(pmax(p_val_adj, 1e-300)),
        significance = case_when(
          p_val_adj < fdr_cutoff & avg_log2FC >  fc_cutoff ~ "↑Low_OA",
          p_val_adj < fdr_cutoff & avg_log2FC < -fc_cutoff ~ "↑High_OA",
          TRUE ~ "Not significant"
        ),
        is_reference = gene %in% reference_genes
      )
    label_genes <- df_ref %>%
      filter(is_reference & significance != "Not significant") %>%
      mutate(gene = paste0("italic('", gene, "')"))
    p_ref_volcano <- ggplot(df_ref, aes(avg_log2FC, neg_log10_fdr)) +
      geom_point(aes(color = significance), alpha = 0.85, size = 3) +
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        parse = TRUE,
        size = 6,
        color = "black",
        box.padding = 0.6,
        point.padding = 0.6,
        force = 3,
        segment.size = 0.3
      ) +
      scale_color_manual(values = c(
        "↑Low_OA"  = "#D85B59",
        "↑High_OA" = "#5271AE",
        "Not significant" = "grey80"
      )) +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                 linetype = "dashed", linewidth = 0.6) +
      geom_hline(yintercept = -log10(fdr_cutoff),
                 linetype = "dashed", linewidth = 0.6) +
      coord_cartesian(xlim = c(-5, 5)) +
      theme_classic(base_size = 20) +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 22, face = "bold", color = "black"),
        axis.text.x  = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y  = element_text(size = 20, color = "black"),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black")
      ) +
      labs(
        x = expression(log[2]("Fold change")),
        y = expression(-log[10]("P.adj"))
      )
    # Save
    ref_vol_path <- file.path(
      sample_dir,
      glue("{sample_label}_reference_marker_volcano.pdf")
    )
    ggsave(ref_vol_path, p_ref_volcano,
           width = 8, height = 6, dpi = 300)
    message(glue("      ✔ Reference volcano saved: {ref_vol_path}"))
    plots$reference_volcano <- p_ref_volcano
  } else {
    message("  No MAST Top-HVGs DEG available — skipping reference volcano.")
  }
  # -----------------------------------------
  # 5. Venn diagram (Universe-filtered only)
  # -----------------------------------------
  message("🔹 Step 5: Creating universe-filtered Venn diagram …")
  mast_top <- deg_out$results$DEG2_TopHVGs_MAST
  # Significant DEGs
  sig_low_DEGs <- mast_top %>%
    filter(avg_log2FC > fc_cutoff, p_val_adj < fdr_cutoff) %>%
    pull(gene)
  # ======================================================
  # 📌 Background universe
  # ======================================================
  mast_genes <- mast_top$gene
  ref_genes  <- full_ref_deg_genes_list
  gene_universe <- intersect(mast_genes, ref_genes)
  message(glue("   • MAST genes: {length(mast_genes)}"))
  message(glue("   • External reference DEG list: {length(ref_genes)}"))
  message(glue("   • Background universe: {length(gene_universe)} genes"))
  # Restrict both sets to the same universe
  reference_set_universe <- intersect(low_output_markers, gene_universe)
  deg_set_universe       <- intersect(sig_low_DEGs, gene_universe)
  # Overlap
  overlap_genes <- intersect(reference_set_universe, deg_set_universe)
  # Sizes for hypergeometric test
  k <- length(overlap_genes)
  m <- length(reference_set_universe)
  q <- length(deg_set_universe)
  N <- length(gene_universe)
  p_hyper <- phyper(
    q = k - 1,
    m = m,
    n = N - m,
    k = q,
    lower.tail = FALSE
  )
  # -------------------------
  # Messages
  # -------------------------
  message(glue("   • Reference markers in universe: {m}"))
  message(glue("   • DEGs in universe: {q}"))
  message(glue("   • Overlap size: {k}"))
  message(glue("   • Hypergeometric P-value: {signif(p_hyper, 3)}"))
  if (p_hyper < 0.05) {
    message("   ✔ Overlap statistically enriched (p < 0.05)")
  } else {
    message("   ⚠ Overlap is NOT statistically enriched")
  }
  message("   • Overlapping genes:")
  if (k == 0) {
    message("     (None)")
  } else {
    message("     ", paste(overlap_genes, collapse = ", "))
  }
  # -------------------------
  # Filtered Venn diagram
  # -------------------------
  venn_list_filtered <- list(
    low_output_markers = reference_set_universe,
    low_output_DEG     = deg_set_universe
  )
  p_venn_filtered <- ggvenn(
    venn_list_filtered,
    fill_color = c("#4E9AC7", "#F79A63"),
    text_size = 12,
    stroke_size = 0.6,
    show_percentage = FALSE
  ) +
    ggtitle(glue("{sample_label}: Low Output — Universe-filtered")) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  ggsave(
    file.path(sample_dir, glue("{sample_label}_LowOutput_venn_filtered.pdf")),
    p_venn_filtered,
    width = 6.5,
    height = 5.5
  )
  message("Done: Filtered Venn diagram saved.\n")
  # -----------------------------------------
  # 6. Contour UMAP
  # -----------------------------------------
  message("🔹 Step 6: Creating contour UMAP …")
  meta_df <- seurat_obj@meta.data %>%
    dplyr::select(celltype) %>%
    cbind(
      Embeddings(seurat_obj, reduction)[,1:2] %>%
        as.data.frame() %>%
        setNames(c("UMAP_1","UMAP_2"))
    )
  p_contour <- ggplot(meta_df, aes(UMAP_1, UMAP_2)) +
    rasterise(
      geom_point(
        aes(fill = celltype),
        shape = 21, color = "black",
        size = 1.2, stroke = 0.15, alpha = 0.8
      ),
      dpi = 300
    ) +
    geom_density_2d(color = "black", linewidth = 0.2) +
    scale_fill_manual(values = celltype_colors) +
    coord_equal() +
    theme_void(base_size = 18) +
    ggtitle(glue("{sample_label} — Contour UMAP"))
  ggsave(
    filename = file.path(sample_dir, glue("{sample_label}_contour_umap.pdf")),
    plot = p_contour,
    width = 9, height = 7, device = cairo_pdf
  )
  message("Done: Contour UMAP saved.\n")
  # -----------------------------------------
  # 7. Save results object
  # -----------------------------------------  
  results <- list(
    OA_by_clone = clone_freq,
    OA_by_cell = cell_OA_df,
    low_clones = low_clones,
    high_clones = high_clones,
    DEG_results = deg_out,
    folder = sample_dir,
    # Add plots here
    plot_dim = p_dim,
    plot_OA = p_oa,
    plot_DEG_volcano = deg_out$volcano_plot,
    plot_ref_volcano = if (exists("p_ref_volcano")) p_ref_volcano else NULL,
    plot_venn_filtered = p_venn_filtered,
    plot_contour = p_contour
  )
  saveRDS(results, file = file.path(sample_dir, "results_object.rds"))
  message(glue("\n COMPLETE! All outputs saved to: {sample_dir}\n"))
  return(results)
}

############################

## Prepare Vitro input ##

############################

# Prepare seurat_vitro_exp1_Y

young_samples <- c("1_Y2", "2_Y3", "3_Y4")
seurat_vitro_exp1_Y <- subset(seurat_vitro_exp1,subset = sampleName %in% young_samples)
seurat_vitro_exp1_Y$CloneID <- seurat_vitro_exp1_Y$uniqueClonesTraced
table(seurat_vitro_exp1_Y$celltype)
# Ready: seurat_vitro_exp1_Y

# Prepare seurat_vitro_exp1_O

seurat_vitro_exp1
table(seurat_vitro_exp1$sampleName)
old_samples <- c('4_O2','5_O3','6_O4')
seurat_vitro_exp1_O <- subset(seurat_vitro_exp1,subset = sampleName %in% old_samples)
seurat_vitro_exp1_O # 15420 cells 
table(seurat_vitro_exp1_O$celltype)
table(Idents(seurat_vitro_exp1_O))
seurat_vitro_exp1_O$CloneID <- seurat_vitro_exp1_O$uniqueClonesTraced
# Ready: seurat_vitro_exp1_O

#######

# Prepare seurat_vitro_exp2_Y
seurat_vitro_exp2$celltype <- seurat_vitro_exp2$celltype_final
table(seurat_vitro_exp2$sampleName)
young_samples_exp2 <- c('1_Ya','2_Ya','3_Ya')
seurat_vitro_exp2_Y <- subset(seurat_vitro_exp2,subset = sampleName %in% young_samples_exp2)
seurat_vitro_exp2_Y # 25322 cells 
table(seurat_vitro_exp2_Y$celltype)
table(Idents(seurat_vitro_exp2_Y))
# Ready: seurat_vitro_exp2_Y

#######

# Prepare seurat_vitro_exp2_O
old_samples_exp2 <- c('4_Oa','5_Oa','6_Oa')
seurat_vitro_exp2_O <- subset(seurat_vitro_exp2,subset = sampleName %in% old_samples_exp2)
seurat_vitro_exp2_O # 30177 cells 
table(seurat_vitro_exp2_O$celltype)
# Ready: seurat_vitro_exp2_O

#######

# Prepare seurat_in_vitro_day0_Y
seurat_in_vitro_day0$celltype <- seurat_in_vitro_day0$celltype_final

table(seurat_in_vitro_day0$sampleName)
seurat_in_vitro_day0_Y <-subset(seurat_in_vitro_day0, 
                                subset = sampleName %in% "Y_vitro")
table(seurat_in_vitro_day0_Y$celltype)
Idents(seurat_in_vitro_day0_Y) <- seurat_in_vitro_day0_Y$celltype
table(Idents(seurat_in_vitro_day0_Y)) # Y_vitro

#######

# Prepare seurat_in_vitro_day0_O
seurat_in_vitro_day0_O <-subset(seurat_in_vitro_day0, 
                                subset = sampleName %in% "O_vitro")
seurat_in_vitro_day0_O
Idents(seurat_in_vitro_day0_O) <- seurat_in_vitro_day0_O$celltype
table(Idents(seurat_in_vitro_day0_O)) # O_vitro

##############################

## Prepare Vivo input ##

##############################

# Prepare seurat_vivo_exp1_Y
table(seurat_vivo_exp1$sampleName)
young_vivo_exp1_sample <- "1_YOUNG"
seurat_vivo_exp1_Y <- subset(seurat_vivo_exp1,subset = sampleName %in% young_vivo_exp1_sample)
seurat_vivo_exp1_Y # 3547 cells
table(seurat_vivo_exp1_Y$celltype)
# Ready: seurat_vivo_exp1_Y

########

# Prepare seurat_vivo_exp1_O
old_vivo_exp1_sample <- '2_Old'
seurat_vivo_exp1_O <- subset(seurat_vivo_exp1,subset = sampleName %in% old_vivo_exp1_sample)
seurat_vivo_exp1_O # 3703 cells 
table(seurat_vivo_exp1_O$celltype)
# Ready: seurat_vivo_exp1_O

#########

# Prepare seurat_vivo_exp2_Y
table(seurat_vivo_exp2$sampleName)
young_vivo_exp2_sample <- c('Y1B','Y2B','Y3B')
seurat_vivo_exp2_Y <- subset(seurat_vivo_exp2,subset = sampleName %in% young_vivo_exp2_sample)
seurat_vivo_exp2_Y # 10108 cells 
seurat_vivo_exp2_Y$CloneID
table(seurat_vivo_exp2_Y$celltype)
# Ready: seurat_vivo_exp2_Y

#########
# Prepare seurat_vivo_exp2_O
table(seurat_vivo_exp2$sampleName)
old_vivo_exp2_sample <- c('O1B','O2B','O3B')
seurat_vivo_exp2_O <- subset(seurat_vivo_exp2,subset = sampleName %in% old_vivo_exp2_sample)
seurat_vivo_exp2_O # 10598 cells 
table(seurat_vivo_exp2_O$celltype)
# Ready: seurat_vivo_exp2_O

#########

# Prepare seurat_in_vivo_day60_Y
seurat_in_vivo_day60_Y <-subset(seurat_in_vivo_day60,subset = sampleName %in% 'YY')
seurat_in_vivo_day60_Y # 3161 cells
Idents(seurat_in_vivo_day60_Y) <- seurat_in_vivo_day60_Y$celltype_final
seurat_in_vivo_day60_Y$celltype <- seurat_in_vivo_day60_Y$celltype_final
table(Idents(seurat_in_vivo_day60_Y))

#########

# Prepare seurat_in_vivo_day60_O
seurat_in_vivo_day60_O <-subset(seurat_in_vivo_day60,subset = sampleName %in% 'OO')
seurat_in_vivo_day60_O # 7643 cells
Idents(seurat_in_vivo_day60_O) <- seurat_in_vivo_day60_O$celltype_final
seurat_in_vivo_day60_O$celltype <- seurat_in_vivo_day60_O$celltype_final
table(Idents(seurat_in_vivo_day60_O))

###########################################

### Run EXP1 and EXP2 OA pipeline ###

setwd('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure1/OA_pipeline_out')

# Run both Exp1 and Exp2 
sample_list <- list(
  Exp1_vitro_Y = seurat_vitro_exp1_Y,
  Exp1_vitro_O = seurat_vitro_exp1_O,
  Exp2_vitro_Y = seurat_vitro_exp2_Y,
  Exp2_vitro_O = seurat_vitro_exp2_O,
  
  Exp1_vivo_Y  = seurat_vivo_exp1_Y,
  Exp1_vivo_O  = seurat_vivo_exp1_O,
  Exp2_vivo_Y  = seurat_vivo_exp2_Y,
  Exp2_vivo_O  = seurat_vivo_exp2_O
)

# Run Exp2 only
sample_list <- list(
  Exp2_vitro_Y = seurat_vitro_exp2_Y,
  Exp2_vitro_O = seurat_vitro_exp2_O,
  Exp2_vivo_Y  = seurat_vivo_exp2_Y,
  Exp2_vivo_O  = seurat_vivo_exp2_O
)

# Run Exp2 Y vivo only 
sample_list <- list(Exp2_vivo_Y = seurat_vivo_exp2_Y)

# Define global parameters

HSC_types       <- c("HSC")
deg_celltypes   <- c("HSC", "MPP")
low_markers     <- low_output_markers
high_markers    <- high_output_markers
full_ref_deg_genes_list <- cellstem_genes

# DEG parameters:
top_hvgs        <- 2000
min_detect_prop <- 0.10
min_mean_norm   <- 0.10
min_mean_raw    <- 2
fc_cutoff       <- 0.25
fdr_cutoff      <- 0.05
output_dir      <- "OA_Analysis_Exp2_vivo_Y"

results_list <- list()

for (sample_name in names(sample_list)) {
  
  message("\n==========================================")
  message(glue("Starting OA analysis for: {sample_name}"))
  message("==========================================\n")
  
  obj <- sample_list[[sample_name]]
  
  # ensure CloneID exists
  if (!"CloneID" %in% colnames(obj@meta.data)) {
    obj$CloneID <- obj$uniqueClonesTraced
  }
  
  # run the full OA pipeline
  res <- run_low_high_OA_analysis_for_a_single_sample(
    seurat_obj           = obj,
    sample_label         = sample_name,
    output_dir      = output_dir,
    
    HSC_types            = HSC_types,
    deg_celltypes        = deg_celltypes,
    
    low_output_markers   = low_markers,
    high_output_markers  = high_markers,
    
    clone_col            = "CloneID",
    
    # DEG params
    top_hvgs             = top_hvgs,
    min_detect_prop      = min_detect_prop,
    min_mean_norm        = min_mean_norm,
    min_mean_raw         = min_mean_raw,
    fc_cutoff            = fc_cutoff,
    fdr_cutoff           = fdr_cutoff,
    celltype_colors = celltype_colors,
    full_ref_deg_genes_list = cellstem_genes
  )
  
  # store results
  results_list[[sample_name]] <- res
  
  message(glue("\n Completed OA Analysis for {sample_name}\n"))
}

###############

### Run CrossAge (Day0 and Day60) OA pipeline ###

# define global paramters 

HSC_types     <- c("HSC")
deg_celltypes <- c("HSC", "MPP")
low_markers   <- low_output_markers
high_markers  <- high_output_markers
full_ref_deg_genes_list <- cellstem_genes
# DEG parameters
top_hvgs        <- 2000
min_detect_prop <- 0.10
min_mean_norm   <- 0.10
min_mean_raw    <- 2
fc_cutoff       <- 0.25
fdr_cutoff      <- 0.05
output_dir <- "OA_Analysis_CrossAge_Day0_Day60"

### Run OA function for CrossAge Vitro (Y&O) and Vivo (YY&OO) ###

# Y_vitro
names(seurat_in_vitro_day0_Y@reductions)

res_day0_vitro_Y <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vitro_day0_Y,
  sample_label = "day0_vitro_Y",
  output_dir = "OA_Analysis_Day0_Day60",
  clone_col = "CloneID",
  HSC_types = c("HSC"),
  deg_celltypes = c("HSC", "MPP"),
  reduction = 'umap',
  low_output_markers = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors,
  full_ref_deg_genes_list = cellstem_genes
)
res_day0_vitro_Y$plot_contour
res_day0_vitro_Y$plot_ref_volcano

# O_vitro

res_day0_vitro_O <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vitro_day0_O,
  sample_label = "day0_vitro_O",
  output_dir = "OA_Analysis_Day0_Day60",
  clone_col = "CloneID",
  HSC_types = c("HSC"),
  deg_celltypes = c("HSC", "MPP"),
  low_output_markers = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors,
  full_ref_deg_genes_list = cellstem_genes
)
res_day0_vitro_O$plot_contour
res_day0_vitro_O$plot_ref_volcano

# YY
res_day60_Y <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vivo_day60_Y,
  sample_label = "day60_Y_(YY)",
  output_dir = "OA_Analysis_Day0_Day60",
  clone_col = "CloneID",
  HSC_types = c("HSC"),
  deg_celltypes = c("HSC","MPP"),
  low_output_markers = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors, 
  full_ref_deg_genes_list = cellstem_genes
)
res_day60_Y$plot_contour

# OO
res_day60_vivo_O <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vivo_day60_O,
  sample_label = "day60_vivo_O_(OO)",
  output_dir = "OA_Analysis_Day0_Day60",
  clone_col = "CloneID",
  HSC_types = c("HSC"),
  deg_celltypes = c("HSC", "MPP"),
  low_output_markers = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors,
  full_ref_deg_genes_list = cellstem_genes
)
res_day60_vivo_O$plot_contour
res_day60_vivo_O$plot_venn_filtered
res_day60_vivo_O$plot_ref_volcano

###################

### RUN OA function for OY and YO ###

table(seurat_in_vivo_day60$sampleName)
seurat_in_vivo_day60$celltype <- seurat_in_vivo_day60$celltype_final
# -------------------------
# Day60 cross-age subsets
# -------------------------
seurat_in_vivo_day60_OY <- subset(
  seurat_in_vivo_day60,
  subset = sampleName == "OY"
)

seurat_in_vivo_day60_YO <- subset(
  seurat_in_vivo_day60,
  subset = sampleName == "YO"
)

Idents(seurat_in_vivo_day60_OY) <- seurat_in_vivo_day60_OY$celltype
Idents(seurat_in_vivo_day60_YO) <- seurat_in_vivo_day60_YO$celltype

table(Idents(seurat_in_vivo_day60_OY))
table(Idents(seurat_in_vivo_day60_YO))

# OY 
res_day60_OY <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vivo_day60_OY,
  sample_label = "Day60_OY",
  output_dir = "OA_Analysis_Day60_CrossAge",
  clone_col = "CloneID",
  HSC_types = c("HSC","MPP"),
  deg_celltypes = c("HSC","MPP"),
  low_output_markers  = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors,
  full_ref_deg_genes_list = cellstem_genes
)
# check results OY 
res_day60_OY$plot_DEG_volcano
res_day60_OY$plot_dim
res_day60_OY$plot_contour
res_day60_OY$plot_DEG_volcano


# YO
res_day60_YO <- run_low_high_OA_analysis_for_a_single_sample(
  seurat_obj = seurat_in_vivo_day60_YO,
  sample_label = "Day60_YO",
  output_dir = "OA_Analysis_Day60_CrossAge",
  clone_col = "CloneID",
  HSC_types = c("HSC","MPP"),
  deg_celltypes = c("HSC","MPP"),
  low_output_markers  = low_output_markers,
  high_output_markers = high_output_markers,
  celltype_colors = celltype_colors,
  full_ref_deg_genes_list = cellstem_genes
)
# check results YO
res_day60_YO$plot_dim
res_day60_YO$plot_contour
res_day60_YO$plot_ref_volcano

########################################################









