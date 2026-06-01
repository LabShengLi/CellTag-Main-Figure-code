library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# =========================
# 1. Load data
# =========================
seurat.vivo.cross <- readRDS("data/CrossAge(exp2)_vivo.RDS")
seurat.vitro.cross <- readRDS("data/CrossAge(exp2)_vitro.RDS")

file_path <- "file_input/Fig3_DVG_Master_Table.xlsx"


# =========================
# 2. Define functions
# =========================
calc_clone_pool_freq <- function(meta,
                                 clone_col = "CloneID",
                                 hsc_types = c("LT-HSC", "ST-HSC","HSC","MPP3", "MPP4","MPP"), 
                                 celltype_col = "celltype_final") {
  
  meta_sub <- meta %>%
    dplyr::filter(
      !is.na(.data[[clone_col]]),
      .data[[clone_col]] != "0"
    ) %>%
    dplyr::mutate(
      is_HSC = .data[[celltype_col]] %in% hsc_types
    )
  
  total_hsc    <- sum(meta_sub$is_HSC)
  total_nonhsc <- sum(!meta_sub$is_HSC)
  
  meta_sub %>%
    dplyr::group_by(.data[[clone_col]]) %>%
    dplyr::summarise(
      n_HSC = sum(is_HSC),
      n_nonHSC = sum(!is_HSC),
      HSC_freq = ifelse(total_hsc == 0, NA_real_, n_HSC / total_hsc),
      nonHSC_freq = ifelse(total_nonhsc == 0, NA_real_, n_nonHSC / total_nonhsc),
      .groups = "drop"
    )
}

calc_clone_pool_freq_seurat <- function(
    seurat_obj,
    clone_col = "CloneID",
    celltype_col = "celltype_final",
    hsc_types = c("LT-HSC", "ST-HSC","HSC","MPP3", "MPP4","MPP"), 
    min_total_cells = 0
) {
  
  meta <- seurat_obj@meta.data
  
  freq <- calc_clone_pool_freq(
    meta = meta,
    clone_col = clone_col,
    celltype_col = celltype_col,
    hsc_types = hsc_types
  )
  
  freq %>%
    dplyr::mutate(
      total_cells = n_HSC + n_nonHSC
    ) %>%
    dplyr::filter(total_cells >= min_total_cells) %>%
    dplyr::select(-total_cells)
}

# =========================
# 3. Compute HSC frequency and define groups
# =========================
seurat.vivo.cross.oo <- subset(seurat.vivo.cross, orig.ident == "OO")

HSC_freq_OO <- calc_clone_pool_freq_seurat(seurat.vivo.cross.oo)

HSC_freq_OO2 <- HSC_freq_OO %>% filter(HSC_freq > 0.003)

q30 <- quantile(HSC_freq_OO2$HSC_freq, 0.3)
q70 <- quantile(HSC_freq_OO2$HSC_freq, 0.7)

top_ids <- HSC_freq_OO2 %>% filter(HSC_freq >= q70) %>% pull(CloneID)
bottom_ids <- HSC_freq_OO2 %>% filter(HSC_freq <= q30) %>% pull(CloneID)

# =========================
# 4. Load DVG gene list
# =========================
read_genes <- function(sheet_id){
  read_excel(file_path, sheet = sheet_id)$gene
}

genes_use <- unique(c(
  read_genes(10),
  read_genes(11),
  read_genes(12),
  read_genes(13)
)) %>%
  na.omit() %>%
  trimws() %>%
  .[. != ""]

# =========================
# 5. Expression matrix and clone-level mean
# =========================
seurat.vitro.cross.o <- subset(seurat.vitro.cross, orig.ident == "O_vitro")

HSC_obj <- subset(
  seurat.vitro.cross.o,
  subset = celltype_final %in% c("LT-HSC","ST-HSC","HSC","MPP") &
    !is.na(CloneID) & CloneID != 0
)

expr_mat <- GetAssayData(HSC_obj, layer = "data")

genes_use_filtered <- intersect(genes_use, rownames(expr_mat))
expr_mat <- expr_mat[genes_use_filtered, ]

expr_df <- as.data.frame(t(as.matrix(expr_mat)))
expr_df$CloneID <- HSC_obj$CloneID

clone_mean <- expr_df %>%
  group_by(CloneID) %>%
  summarise(across(all_of(genes_use_filtered), mean, na.rm = TRUE))

# =========================
# 6. Compute summary statistics
# =========================
get_summary <- function(top_ids, bottom_ids, clone_mean){
  
  top_mean <- clone_mean %>% filter(CloneID %in% top_ids)
  bottom_mean <- clone_mean %>% filter(CloneID %in% bottom_ids)
  
  Mean_Top <- colMeans(top_mean[,-1])
  Mean_Bottom <- colMeans(bottom_mean[,-1])
  Mean_All <- colMeans(clone_mean[,-1])
  
  df <- data.frame(
    gene = names(Mean_All),
    Mean_All = Mean_All,
    Mean_Top = Mean_Top,
    Mean_Bottom = Mean_Bottom
  )
  
  df %>%
    mutate(
      Top_minus_All = Mean_Top - Mean_All,
      Bottom_minus_All = Mean_Bottom - Mean_All,
      Top_minus_Bottom = Mean_Top - Mean_Bottom,
      Abs_Top = abs(Top_minus_All),
      Abs_Bottom = abs(Bottom_minus_All)
    )
}

summary_df <- get_summary(top_ids, bottom_ids, clone_mean)

# =========================
# 7. Visualization
# =========================
long_df <- summary_df %>%
  select(Abs_Top, Abs_Bottom) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Value") %>%
  mutate(Group = recode(Group,
                        Abs_Top = "Top 30%",
                        Abs_Bottom = "Bottom 30%"))

p_1 = ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = c("Top 30%" = "tomato",
                               "Bottom 30%" = "turquoise3")) +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold")
  )

# Wilcoxon test
wilcox_res <- wilcox.test(
  Value ~ Group,
  data = long_df
)

wilcox_res

long_df <- summary_df %>%
  dplyr::select(Abs_Top, Abs_Bottom) %>%
  tidyr::pivot_longer(everything(), names_to = "Group", values_to = "Value") %>%
  dplyr::mutate(Group = dplyr::recode(Group,
                                      "Abs_Top" = "Top",
                                      "Abs_Bottom" = "Bottom"),
                Group = factor(Group, levels = c("Bottom", "Top")))
p_2 = ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("Top" = "tomato",
                               "Bottom" = "turquoise3")) +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold")
  )

p_2
# =========================
# 8. Random sampling control
# =========================
set.seed(135)

pool_ids <- HSC_freq_OO2$CloneID

n_top <- length(top_ids)
n_bottom <- length(bottom_ids)

random_top_ids <- sample(pool_ids, n_top)

random_bottom_ids <- sample(
  setdiff(pool_ids, random_top_ids),
  n_bottom
)

# Compute random clone means
random_top_clone_mean <- clone_mean %>%
  filter(CloneID %in% random_top_ids)

random_bottom_clone_mean <- clone_mean %>%
  filter(CloneID %in% random_bottom_ids)

Mean_Top_r <- colMeans(random_top_clone_mean[,-1], na.rm = TRUE)
Mean_Bottom_r <- colMeans(random_bottom_clone_mean[,-1], na.rm = TRUE)
Mean_All <- colMeans(clone_mean[,-1], na.rm = TRUE)

# Align gene order
Mean_Top_r <- Mean_Top_r[names(Mean_All)]
Mean_Bottom_r <- Mean_Bottom_r[names(Mean_All)]

summary_df_random <- data.frame(
  gene = names(Mean_All),
  Mean_All = Mean_All,
  Mean_Top = Mean_Top_r,
  Mean_Bottom = Mean_Bottom_r
)

summary_df_random$Abs_Top <- abs(summary_df_random$Mean_Top - summary_df_random$Mean_All)
summary_df_random$Abs_Bottom <- abs(summary_df_random$Mean_Bottom - summary_df_random$Mean_All)

# Convert to long format
long_df_random <- summary_df_random %>%
  select(Abs_Top, Abs_Bottom) %>%
  pivot_longer(cols = everything(),
               names_to = "Group",
               values_to = "Value") %>%
  mutate(Group = recode(Group,
                        Abs_Top = "Random Top",
                        Abs_Bottom = "Random Bottom"))

# Plot random control
p_random <- ggplot(long_df_random, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = c("Random Top" = "tomato",
                               "Random Bottom" = "turquoise3")) +
  coord_cartesian(ylim = c(0, 0.4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold")
  )

p_random

# Paired Wilcoxon test
wilcox.test(summary_df_random$Abs_Top,
            summary_df_random$Abs_Bottom,
            paired = TRUE)

library(patchwork)
p_1 + p_random

# =========================
# 9. Permutation test
# =========================
set.seed(1)

pool_ids <- intersect(HSC_freq_OO2$CloneID, clone_mean$CloneID)
n_top <- length(top_ids)
n_bottom <- length(bottom_ids)

n_perm <- 1000
wilcox_stats <- numeric(n_perm)

Mean_All <- colMeans(clone_mean[,-1])

for (i in seq_len(n_perm)) {
  
  shuffled <- sample(pool_ids)
  rand_top <- shuffled[1:n_top]
  rand_bottom <- shuffled[(n_top+1):(n_top+n_bottom)]
  
  rand_top_mean <- colMeans(clone_mean[clone_mean$CloneID %in% rand_top, -1])
  rand_bottom_mean <- colMeans(clone_mean[clone_mean$CloneID %in% rand_bottom, -1])
  
  wilcox_stats[i] <- wilcox.test(
    abs(rand_top_mean - Mean_All),
    abs(rand_bottom_mean - Mean_All),
    paired = TRUE
  )$statistic
}

real_stat <- wilcox.test(summary_df$Abs_Top,
                         summary_df$Abs_Bottom,
                         paired = TRUE)$statistic

p_empirical <- mean(wilcox_stats >= real_stat)

# =========================
# 10. Permutation plot
# =========================
hist(wilcox_stats, breaks = 30, col = "grey",
     main = "Permutation Test",
     xlab = "Wilcoxon statistic")

abline(v = real_stat, col = "red", lwd = 3)

mean_perm <- mean(wilcox_stats)
sd_perm <- sd(wilcox_stats)
z_score <- (real_stat - mean_perm) / sd_perm
