# ===============================
# Libraries
# ===============================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

seurat.vivo.cross <- readRDS(
  "data/CrossAge(exp2)_vivo.RDS"
)
seurat.vitro.cross <- readRDS(
  "data/CrossAge(exp2)_vitro.RDS"
)

seurat.vivo.cross.oo <- subset(seurat.vivo.cross, orig.ident == "OO")
seurat.vivo.cross.oy <- subset(seurat.vivo.cross, orig.ident == "OY")
seurat.vivo.cross.yy <- subset(seurat.vivo.cross, orig.ident == "YY")
seurat.vivo.cross.yo <- subset(seurat.vivo.cross, orig.ident == "YO")
# ==========================================================
# OO vs OY
# ==========================================================

HSC_types <- c("LT-HSC", "ST-HSC", "HSC")
epsilon   <- 1e-4
cutoff    <- 5
clone_col <- "CloneID"

# -------- Subset --------
meta_oo <- subset(seurat.vivo.cross, orig.ident == "OO")@meta.data
meta_oy <- subset(seurat.vivo.cross, orig.ident == "OY")@meta.data

# -------- Shared clones --------
shared_clones_old <- intersect(
  unique(meta_oo[[clone_col]][meta_oo[[clone_col]] != "0"]),
  unique(meta_oy[[clone_col]][meta_oy[[clone_col]] != "0"])
)

meta_oo_sub <- meta_oo %>% filter(CloneID %in% shared_clones_old)
meta_oy_sub <- meta_oy %>% filter(CloneID %in% shared_clones_old)

# -------- Clone filter (total ≥ 5) --------
count_oo <- table(meta_oo_sub[[clone_col]])
count_oy <- table(meta_oy_sub[[clone_col]])

all_clones <- shared_clones_old

count_oo_vec <- setNames(as.numeric(count_oo[all_clones]), all_clones)
count_oy_vec <- setNames(as.numeric(count_oy[all_clones]), all_clones)

count_oo_vec[is.na(count_oo_vec)] <- 0
count_oy_vec[is.na(count_oy_vec)] <- 0

total_count_old <- count_oo_vec + count_oy_vec

kept_clones_old <- names(total_count_old)[total_count_old >= 5]
filtered_clones_old <- names(total_count_old)[total_count_old < 5]

# -------- Apply filter --------
meta_oo_kept <- meta_oo %>% filter(CloneID %in% kept_clones_old)
meta_oy_kept <- meta_oy %>% filter(CloneID %in% kept_clones_old)

# -------- Total HSC --------
total_hsc_oo <- sum(meta_oo_kept$celltype_final %in% HSC_types)
total_hsc_oy <- sum(meta_oy_kept$celltype_final %in% HSC_types)

# -------- HSC freq per clone --------
freq_oo <- meta_oo_kept %>%
  group_by(CloneID) %>%
  summarise(
    HSC_freq_OO = sum(celltype_final %in% HSC_types) / total_hsc_oo,
    .groups = "drop"
  )

freq_oy <- meta_oy_kept %>%
  group_by(CloneID) %>%
  summarise(
    HSC_freq_OY = sum(celltype_final %in% HSC_types) / total_hsc_oy,
    .groups = "drop"
  )

freq_old <- full_join(freq_oo, freq_oy, by = "CloneID") %>%
  mutate(
    HSC_freq_OO = replace_na(HSC_freq_OO, 0),
    HSC_freq_OY = replace_na(HSC_freq_OY, 0),
    log2FC_old = log2((HSC_freq_OO + epsilon) /
                        (HSC_freq_OY + epsilon))
  )

# ==========================================================
# YY vs YO
# ==========================================================

meta_yy <- subset(seurat.vivo.cross, orig.ident == "YY")@meta.data
meta_yo <- subset(seurat.vivo.cross, orig.ident == "YO")@meta.data

shared_clones_young <- intersect(
  unique(meta_yy$CloneID[meta_yy$CloneID != "0"]),
  unique(meta_yo$CloneID[meta_yo$CloneID != "0"])
)

meta_yy_sub <- meta_yy %>% filter(CloneID %in% shared_clones_young)
meta_yo_sub <- meta_yo %>% filter(CloneID %in% shared_clones_young)

# -------- Clone filter --------
count_yy <- table(meta_yy_sub$CloneID)
count_yo <- table(meta_yo_sub$CloneID)

all_clones_y <- shared_clones_young

count_yy_vec <- setNames(as.numeric(count_yy[all_clones_y]), all_clones_y)
count_yo_vec <- setNames(as.numeric(count_yo[all_clones_y]), all_clones_y)

count_yy_vec[is.na(count_yy_vec)] <- 0
count_yo_vec[is.na(count_yo_vec)] <- 0

total_count_young <- count_yy_vec + count_yo_vec

kept_clones_young <- names(total_count_young)[total_count_young >= 5]
filtered_clones_young <- names(total_count_young)[total_count_young < 5]

# -------- Apply filter --------
meta_yy_kept <- meta_yy %>% filter(CloneID %in% kept_clones_young)
meta_yo_kept <- meta_yo %>% filter(CloneID %in% kept_clones_young)

# -------- Total HSC --------
total_hsc_yy <- sum(meta_yy_kept$celltype_final %in% HSC_types)
total_hsc_yo <- sum(meta_yo_kept$celltype_final %in% HSC_types)

# -------- HSC freq --------
freq_yy <- meta_yy_kept %>%
  group_by(CloneID) %>%
  summarise(
    HSC_freq_YY = sum(celltype_final %in% HSC_types) / total_hsc_yy,
    .groups = "drop"
  )

freq_yo <- meta_yo_kept %>%
  group_by(CloneID) %>%
  summarise(
    HSC_freq_YO = sum(celltype_final %in% HSC_types) / total_hsc_yo,
    .groups = "drop"
  )

freq_young <- full_join(freq_yy, freq_yo, by = "CloneID") %>%
  mutate(
    HSC_freq_YY = replace_na(HSC_freq_YY, 0),
    HSC_freq_YO = replace_na(HSC_freq_YO, 0),
    log2FC_young = log2((HSC_freq_YO + epsilon) /
                          (HSC_freq_YY + epsilon))
  )
high_old  <- sum(freq_old$log2FC_old > cutoff)
low_old   <- sum(freq_old$log2FC_old < -cutoff)
unch_old  <- nrow(freq_old) - high_old - low_old

high_young <- sum(freq_young$log2FC_young > cutoff)
low_young  <- sum(freq_young$log2FC_young < -cutoff)
unch_young <- nrow(freq_young) - high_young - low_young

tab <- matrix(
  c(high_old, unch_old, low_old,
    high_young, unch_young, low_young),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Donor = c("Old", "Young"),
    Fate  = c("High", "Unchanged", "Low")
  )
)

tab
chisq.test(tab)


# ==========================================================
# Bubble plot — OO vs OY
# ==========================================================

df_plot_old <- freq_old %>%
  arrange(log2FC_old) %>%
  mutate(
    rank = row_number(),
    max_hsc_freq = pmax(HSC_freq_OO, HSC_freq_OY, na.rm = TRUE)
  )

up_pct_old   <- round(mean(df_plot_old$log2FC_old > cutoff) * 100, 1)
down_pct_old <- round(mean(df_plot_old$log2FC_old < -cutoff) * 100, 1)

# ==========================================================
# Bubble plot — YY vs YO
# ==========================================================

df_plot_young <- freq_young %>%
  arrange(log2FC_young) %>%
  mutate(
    rank = row_number(),
    max_hsc_freq = pmax(HSC_freq_YY, HSC_freq_YO, na.rm = TRUE)
  )

up_pct_young   <- round(mean(df_plot_young$log2FC_young > cutoff) * 100, 1)
down_pct_young <- round(mean(df_plot_young$log2FC_young < -cutoff) * 100, 1)



global_max <- max(
  df_plot_old$max_hsc_freq,
  df_plot_young$max_hsc_freq,
  na.rm = TRUE
)
size_scale <- scale_size_continuous(
  name   = "HSC freq",
  limits = c(0, global_max),
  breaks = c(0.00, 0.05, 0.1, round(global_max, 2)),
  range  = c(2, 12)
)

fill_scale <- scale_fill_gradient2(
  low      = "#2ecc71",
  mid      = "white",
  high     = "#f39c12",
  midpoint = 0,
  limits   = c(-3, 3),
  oob      = scales::squish,
  name     = expression(bold(Log[2]~"FC"))
)
p_old <- ggplot(df_plot_old, aes(
  x = rank,
  y = log2FC_old,
  size = max_hsc_freq,
  fill = log2FC_old
)) +
  geom_point(alpha = 0.75, colour = "black", pch = 24) +
  
  geom_hline(yintercept = c(cutoff, -cutoff),
             linetype = "dashed",
             linewidth = 0.6,
             colour = "grey40") +
  
  size_scale +
  fill_scale +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  ) +
  coord_cartesian(ylim = c(-15, 15)) +
  
  labs(
    title = "OO vs OY",
    x     = "Rank",
    y     = NULL) +
  
  annotate(
    "text",
    x = min(df_plot_old$rank) + 1,
    y = 15 * 0.95,
    label = paste0(up_pct_old, "% clones\nAging adaptation"),
    hjust = 0,
    vjust = 1,
    size  = 5,
    color = "#f39c12",
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = max(df_plot_old$rank) - 1,
    y = -15 * 0.95,
    label = paste0(down_pct_old, "% clones\nAging resistance"),
    hjust = 1,
    vjust = 0,
    size  = 5,
    color = "#2ecc71",
    fontface = "bold"
  ) +
  
  theme_bw(base_size = 20)  +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),

    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 16, face = "bold"),
    
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

p_old

p_young <- ggplot(df_plot_young, aes(
  x = rank,
  y = log2FC_young,
  size = max_hsc_freq,
  fill = log2FC_young
)) +
  geom_point(alpha = 0.75, colour = "black", pch = 24) +
  
  geom_hline(yintercept = c(cutoff, -cutoff),
             linetype = "dashed",
             linewidth = 0.6,
             colour = "grey40") +
  
  scale_size_continuous(
    name   = "HSC freq",
    limits = c(0, global_max),
    breaks = c(0.00, 0.05, 0.1, round(global_max, 2)),
    range  = c(2, 12)
  ) +
  
  scale_fill_gradient2(
    low      = "#2ecc71",
    mid      = "white",
    high     = "#f39c12",
    midpoint = 0,
    limits   = c(-3, 3),
    oob      = scales::squish,
    name     = expression(bold(Log[2]~"FC"))
  ) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  ) +
  coord_cartesian(ylim = c(-15, 15)) +
  
  labs(
    title = "YO vs YY",
    x     = "Rank",
    y     = bquote(bold(Log[2]~"HSC clonal fitness"))
  ) +
  
  annotate(
    "text",
    x = min(df_plot_young$rank) + 1,
    y = 15 * 0.95,
    label = paste0(up_pct_young, "% clones\nAging adaptation"),
    hjust = 0,
    vjust = 1,
    size  = 5,
    color = "#f39c12",
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = max(df_plot_young$rank) - 1,
    y = -15 * 0.95,
    label = paste0(down_pct_young, "% clones\nAging resistance"),
    hjust = 1,
    vjust = 0,
    size  = 5,
    color = "#2ecc71",
    fontface = "bold"
  ) +
  
  theme_bw(base_size = 20)  +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    
    # 🔥 轴
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 16, face = "bold"),
    # 🔥 legend
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold")
  )

p_young
#600 365 svg

p_old_clean <- p_old +
  theme(
    axis.title.y = element_blank()
  )
library(patchwork)

p_CFI_final <- (p_young | p_old_clean)

p_CFI_final + plot_annotation(
  theme = theme(
    axis.title.y = element_text(size = 14)
  )
)

