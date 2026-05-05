###################################################################################################


### # Calculate HSC heterogenecity using PCA distance in Cross(Day0) data


###################################################################################################

# Maintainer: Chris Chen 
# Last updated: 02/06/2026

##########################

setwd('/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts')

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
setwd('Main_figures/Figure3/Figures/')

#############################

# read in seurat object
# seurat_celltag_in_vitro <- readRDS('/project2/sli68423_1316/users/Kailiang/U1_celltag/data/seurat.vitro.cross.no_ambiguous_cloneID.rds')
seurat_celltag_in_vitro <- readRDS("/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/test_folder/Data_objects/crossage_vitro.3_2.rds")

colnames(seurat_celltag_in_vitro@meta.data)
meta <- seurat_celltag_in_vitro@meta.data  
table(meta$celltype) # check celltype 
nrow(meta) # 17308
colnames(meta)
meta_HSC_only <- meta[meta$celltype %in% c("LT-HSC"),]
nrow(meta_HSC_only) # 9421
colnames(meta_HSC_only)
# check number of unique clones 
meta_HSC_only
length(unique(meta_HSC_only$CloneID)) # 1935 
table(meta_HSC_only$sampleName) # O_vitro: 4381 ; Y_vitro: 5040
# filtered out CloneID_clean = 0
meta_HSC_only <- meta_HSC_only[meta_HSC_only$CloneID != '0',]
nrow(meta_HSC_only) # 3278
table(meta_HSC_only$sampleName) # O_vitro: 1218 ; Y_vitro: 2060

##############################

## Aging-up DVGs, we will use these genes as feature to compute clonal PCA distance ##

aging_up_DVGs <- c(
  "Actb","Lmnb1","Cd9","Anp32e","Cct6a","Ywhaz","Eif4g2","Nap1l1","Alyref","Hnrnpd",
  "Pa2g4","Bax","Ube2m","Tubb5","Ybx3","Lyar","Myb","Hnrnpu","Fam136a","Mrpl54",
  "Cdv3","Timm8a1","Lamtor4","Dph3","Ncl","Nucks1","Hnrnpab","Banf1","Atf4","Hdgf",
  "Hspd1","Mybbp1a","Cacybp","Lmo2","Hmgn2","Aldoa","Mrpl12","Gclm","Ptges3","Tagln2",
  "Dynll2","Svil","Gapdh","Rgs18","Lsm2","H2afy","Skp1a","Pgam1","Cenpv","Nhp2",
  "Tceal8","Dok2","Exoc4","Cdkn2d","Rap1b","Ldha","Tuba1b","Unc119","Idi1","Npm3",
  "Pmvk","Hsph1","Timm13","Smc2","Calr","Eno1","Sox4","Smc4","Tomm5","Plscr1",
  "Cks1b","Ubald2","Tuba1a","Asns","Hnrnpa3","Cdca3","Sap30","Tespa1","Sod1",
  "Rangap1","Myh9","Pbk","Msi2","Srm","Fut8","Pdcd4","0610010K14Rik","Tpm1","Sdf2l1",
  "Mef2c","Bub3","Hspa8","Hoxb2","Jund","Hirip3","Jpt1","Mvd","Rgcc","A430005L14Rik",
  "Fam3c","Hmgb3","Nt5dc2","Nolc1","Tacc3","Myef2","Clec4d","H1f0","Crip1","Hnrnpul2",
  "Gng11","Mctp1","Xbp1","Qdpr","Id1","Rrm2","Cenpb","Rangrf","Kpna2","Knstrn",
  "Dbf4","Set","Nuf2","Rad21","Tuba8","Ddt","Terf1","Il21r","Nrgn","Itga6","Mthfd2",
  "Ube2s","Bin1","Birc5","Ggct","Hmgcs1","B3gnt5","Phlda3","Fabp5","Rpl31","Mycn",
  "Myl10","Serpinb1a","Clec4e","Fbxo5","Egfl7","Atp5g1","Fchsd2","Ttc32","Ccng1",
  "Frat2","Plek","Pimreg","Eif2b3","Psrc1","Fli1","Lmo1","Bnip3","Gatm","Cdk1",
  "Uhrf1","Cited2","Kif22","Cdca8","Runx1","Cycs","Ydjc","Afp","Sgo1","Zbed3",
  "Plxdc2","Car2","Igf1r","Tcf15","Racgap1","Snhg3","Mc5r","Fzr1","Calm2","Rpl23a",
  "Cox7b","Kif15","Enah","Cish","Pitpnc1","Aurkb","Cip2a","Lgals9","Cdc20","Diaph3",
  "Hjurp","Zmynd19","Ppp1r14b","Rps17","Cbfa2t3","Kif2c","Cfl1","Bub1b","Nkg7",
  "Cep55","Rad51ap1"
)

# updated in 03/07/2026 aging-up gene list (gt_0.05FC_2_exp)

aging_up_DVGs <- c(
  "Aldh1a1","Nupr1","Slamf1","Cdkn2c","Cd48","Slc14a1","Otos","Krt18","Apoe","Hmgb2",
  "Cavin3","Nkx2-3","Prkca","Ran","Tubb5","Mki67","Birc5","Tespa1","Top2a","Gng11",
  "Cdkn1a","Txnip","Ncl","Zmat3","Serpinb1a","Xbp1","Ccnd1","Tubb4b","Slc22a3","Pclaf",
  "Ybx1","Cks1b","Cdca8","Mid1","Cenpv","Lgals1","Lgals9","Hacd4","H1fx","Nrgn",
  "Smc4","Nfia","Ptms","Tuba1b","Rora","Smc2","Myl10","Cdk1","Fbxo5","Ccna2",
  "Pbk","Plek","Alcam","Cks2","Anln","Esco2","Cd53","Ube2s","Knl1","Hbb-bt",
  "Incenp","Frat2","Kif20b","Procr","Neat1","Dlgap5","Ccng1","Hspe1","Tnip3","Adgrg3",
  "Ckap2l","Slfn9","Pbx1","Sgo1","Racgap1","Brip1","Snhg3","Cycs","Retreg1","Clec4d"
)

###############################

#### Compute PCA ####

seurat_in_vitro_HSC <- subset(seurat_celltag_in_vitro,subset = celltype %in% c("LT-HSC") & CloneID != "0")
seurat_in_vitro_HSC
# seurat_in_vitro_HSC$CloneID <- seurat_in_vitro_HSC$CloneID_clean # use CloneID as the downstream name

DefaultAssay(seurat_in_vitro_HSC) <- "RNA"
valid_genes <- intersect(aging_up_DVGs, rownames(seurat_in_vitro_HSC))
length(valid_genes) # 80 DVGs 

#seurat_in_vitro_HSC <- FindVariableFeatures(seurat_in_vitro_HSC,selection.method = "vst",nfeatures = 2000)
expr_data <- GetAssayData(seurat_in_vitro_HSC, layer = "data")
hvgs <- VariableFeatures(seurat_in_vitro_HSC)

gene_pct <- Matrix::rowMeans(expr_data[hvgs, ] > 0)
gene_avg <- Matrix::rowMeans(expr_data[hvgs, ])
genes_pass <- names(gene_pct)[gene_pct >= 0.1 & gene_avg >= 0.1]

seurat_in_vitro_HSC <- subset(seurat_in_vitro_HSC, features = genes_pass)
seurat_in_vitro_HSC # 1934 genes and 3278 cells 

seurat_in_vitro_HSC <- ScaleData(seurat_in_vitro_HSC,features = rownames(seurat_in_vitro_HSC),verbose = FALSE)
seurat_in_vitro_HSC <- RunPCA(seurat_in_vitro_HSC,features = rownames(seurat_in_vitro_HSC),npcs = 30,verbose = FALSE)

ElbowPlot(seurat_in_vitro_HSC, ndims = 30) +
  ggtitle("Elbow Plot for PCA (Filtered HVGs)")
pcs <- Embeddings(seurat_in_vitro_HSC, reduction = "pca")[, 1:30]
pcs

# Map PCs to cells 
meta <- seurat_in_vitro_HSC@meta.data %>% mutate(cell_id = rownames(seurat_in_vitro_HSC@meta.data))
pc_df <- as.data.frame(pcs) %>% mutate(cell_id = rownames(pcs)) %>%
  inner_join(meta[, c("cell_id", "CloneID", "sampleName")], by = "cell_id")
table(pc_df$sampleName) # O_vitro: 1968 ; Y_vitro: 3403
length(unique(pc_df$CloneID)) # 2513 clones 
head(pc_df[, 1:6])

# calculate centroid per clone
compute_clone_intra <- function(clone_id, df_pc) {
  # subset to the clone
  clone_df <- df_pc %>% filter(CloneID == clone_id)
  # skip clones with <2 cells (no variance)
  if (nrow(clone_df) < 2) {
    return(tibble(
      CloneID = clone_id,
      sampleName = clone_df$sampleName[1],
      n_cells = nrow(clone_df),
      intra_dist = NA_real_,
      centroid = list(NA)
    ))
  }
  # extract PC coordinates
  pc_mat <- as.matrix(clone_df[, grep("^PC", colnames(clone_df))])
  # compute centroid
  centroid <- colMeans(pc_mat, na.rm = TRUE)
  # compute Euclidean distances to centroid
  centered <- sweep(pc_mat, 2, centroid, "-")
  dists <- sqrt(rowSums(centered^2))
  # return metrics
  tibble(
    CloneID    = clone_id,
    sampleName = clone_df$sampleName[1],
    n_cells    = nrow(clone_df),
    intra_dist = mean(dists, na.rm = TRUE),
    centroid   = list(centroid)
  )
}
clone_list <- unique(pc_df$CloneID)
intra_results <- map_dfr(clone_list, ~compute_clone_intra(.x, pc_df))
CrossAge_Day0_intra_results_saved <- intra_results
CrossAge_Day0_intra_results_saved$centroid <- NULL # drop the centroid column 

#### Compute inter clone distance #####

# Function to compute inter clone distance #

compute_inter_clone_distance <- function(intra_results, df_pc,
                                              min_cells = 2,
                                              pc_prefix = "^PC") {
  
  message("STEP 1: Filtering valid clones")
  intra_results_filtered <- intra_results %>%
    dplyr::filter(!is.na(intra_dist), n_cells >= min_cells)
  clone_ids <- intra_results_filtered$CloneID
  message("Clones included: ", length(clone_ids))
  message("STEP 2: Extract PC matrix once")
  pc_cols <- grep(pc_prefix, colnames(df_pc), value = TRUE)
  df_pc <- df_pc %>%
    dplyr::filter(CloneID %in% clone_ids)
  message("PC dimensions: ", length(pc_cols))
  message("STEP 3: Precompute centroid list")
  centroid_list <- intra_results %>%
    dplyr::filter(CloneID %in% clone_ids) %>%
    dplyr::select(CloneID, centroid) %>%
    tibble::deframe()
  message("STEP 4: Split cell PC matrices by clone")
  clone_cells <- split(df_pc[, pc_cols], df_pc$CloneID)
  message("STEP 5: Generate clone pairs")
  pairwise_df <- expand.grid(
    CloneA = clone_ids,
    CloneB = clone_ids,
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(CloneA != CloneB)
  message("Total comparisons: ", nrow(pairwise_df))
  message("STEP 6: Compute distances")
  pb <- txtProgressBar(min = 0, max = nrow(pairwise_df), style = 3)
  inter_dist <- numeric(nrow(pairwise_df))
  for (i in seq_len(nrow(pairwise_df))) {
    clone_A <- pairwise_df$CloneA[i]
    clone_B <- pairwise_df$CloneB[i]
    pcs_A <- as.matrix(clone_cells[[clone_A]])
    centroid_B <- centroid_list[[clone_B]]
    if (is.null(centroid_B) || nrow(pcs_A) == 0) {
      inter_dist[i] <- NA_real_
    } else {
      centered <- sweep(pcs_A, 2, centroid_B, "-")
      dists <- sqrt(rowSums(centered^2))
      inter_dist[i] <- mean(dists)
    }
    if (i %% 1000 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  pairwise_df$inter_dist <- inter_dist
  message("\nSTEP 7: Compute mean inter-clone distance")
  avg_inter <- pairwise_df %>%
    dplyr::group_by(CloneA) %>%
    dplyr::summarise(
      mean_inter_dist = mean(inter_dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(CloneID = CloneA)
  message("STEP 8: Merge results")
  merged <- intra_results %>%
    dplyr::left_join(avg_inter, by = "CloneID")
  message("DONE")
  list(
    summary  = merged %>% dplyr::select(-centroid),
    pairwise = pairwise_df
  )
}

intra_results

res_inter <- compute_inter_clone_distance(intra_results, pc_df, min_cells = 2) # 815 clones 
head(res_inter)
clone_heterogeneity_summary_day0_in_vitro <- res_inter$summary
head(clone_heterogeneity_summary_day0_in_vitro)

# Save 
write.csv(clone_heterogeneity_summary_day0_in_vitro,file = "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Figure3/Tables/clonal_PCS_dist_summary_day0_in_vitro.csv",row.names = FALSE)  
  
##################################################

### Create a violin plot for inter - intra clonal PCA distance ###

clone_df <- clone_heterogeneity_summary_day0_in_vitro

clone_long <- clone_df %>% select(sampleName, intra_dist, mean_inter_dist) %>%
  pivot_longer(cols = c(intra_dist, mean_inter_dist), names_to = "DistanceType",values_to = "Distance") %>%
  mutate(DistanceType = recode(DistanceType,"intra_dist" = "Intra-clone","mean_inter_dist" = "Inter-clone"),
    sampleName = recode(sampleName,"O_vitro" = "O","Y_vitro" = "Y")) %>%
  drop_na(Distance) 
clone_long
clone_long <- clone_long %>%
  mutate(sampleName = factor(sampleName, levels = c("Y", "O")),         
    DistanceType = factor(DistanceType, levels = c("Intra-clone", "Inter-clone")))

ggplot(clone_long, aes(x = sampleName, y = Distance, fill = DistanceType)) +
  geom_violin(trim = TRUE, alpha = 0.7, scale = "width", position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Intra-clone" = "cornflowerblue","Inter-clone" = "orange3")) +
  labs(x = "", y = "PCA distance", fill = NULL,
    title = "Day 0 HSC Heterogeneity") +
  theme_classic(base_size = 18) + 
  theme(
    plot.title  = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18),
    #legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = "top",
    strip.text = element_text(size = 18, face = "bold")
  )
#### Compute within group test ####
# (Young: intra vs inter ; Old: intra vs inter)
within_group_tests <- clone_long %>% group_by(sampleName) %>%
  summarise(p_value = wilcox.test(Distance ~ DistanceType, exact = FALSE)$p.value,intra_mean = mean(Distance[DistanceType == "Intra-clone"], na.rm = TRUE),
    inter_mean = mean(Distance[DistanceType == "Inter-clone"], na.rm = TRUE), n_intra = sum(DistanceType == "Intra-clone"), n_inter = sum(DistanceType == "Inter-clone")) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"),y_pos = max(clone_long$Distance, na.rm = TRUE) * 1.05)
within_group_tests

#### Between group test ####
# (Intra: Y vs O ; Inter: Y vs O)
between_group_tests <- clone_long %>% group_by(DistanceType) %>%
  summarise(p_value = wilcox.test(Distance ~ sampleName, exact = FALSE)$p.value,young_mean = mean(Distance[sampleName == "Young"], na.rm = TRUE),
    old_mean   = mean(Distance[sampleName == "Old"],   na.rm = TRUE), n_young = sum(sampleName == "Young"), n_old = sum(sampleName == "Old")) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ), y_pos = max(clone_long$Distance, na.rm = TRUE) * 1.12)
between_group_tests

### Manually annotated the significance ###
manual_annot <- tibble(
  x_start = c(
    1 - 0.2,   # Young: Intra
    2 - 0.2,   # Old: Intra
    1 - 0.2,   # Between-group Intra
    1 + 0.2    # Between-group Inter
  ),
  x_end = c(
    1 + 0.2,   # Young: Inter
    2 + 0.2,   # Old: Inter
    2 - 0.2,   # Between-group Intra
    2 + 0.2    # Between-group Inter
  ),
  y_pos = c(
    max(clone_long$Distance) * 1.05,
    max(clone_long$Distance) * 1.05,
    max(clone_long$Distance) * 1.15,
    max(clone_long$Distance) * 1.25
  ),
  label = c("***", "***", "***", "***")
)
manual_annot

p <- ggplot(clone_long, aes(x = sampleName, y = Distance, fill = DistanceType)) +
  geom_violin(trim = TRUE, alpha = 0.7, scale = "width",
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(
    "Intra-clone" = "cornflowerblue",
    "Inter-clone" = "orange3"
  )) +
  labs(
    x = "",
    y = "PCA distance",
    fill = NULL,
    title = "Day 0 HSC Heterogeneity"
  ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title  = element_text(size = 20, hjust = 0.5),
    axis.title.y = element_text(size = 21, face = "bold"),
    axis.text.x  = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.position = "top"
  ) +
geom_segment(
  data = manual_annot,
  aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
  inherit.aes = FALSE,
  linewidth = 0.9
) +
geom_text(
  data = manual_annot,
  aes(
    x = (x_start + x_end) / 2,
    y = y_pos + max(clone_long$Distance) * 0.015,
    label = label
  ),
  inherit.aes = FALSE,
  size = 7, fontface = "bold"
)
p # Intra and Inter clonal PCA distance for Day0 HSC 
###############################################################

## End of this script 