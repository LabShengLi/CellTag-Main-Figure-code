library(Seurat)
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(glue)
library(tidyr)

seurat.vivo.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vivo.RDS"
)
seurat.vitro.cross <- readRDS(
  "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vitro.RDS"
)
eps <- 1e-4
seurat.vivo.cross.oo <- subset(seurat.vivo.cross, orig.ident == "OO")
seurat.vivo.cross.oy <- subset(seurat.vivo.cross, orig.ident == "OY")
seurat.vivo.cross.yy <- subset(seurat.vivo.cross, orig.ident == "YY")
seurat.vivo.cross.yo <- subset(seurat.vivo.cross, orig.ident == "YO")

calc_clone_pool_freq <- function(meta,
                                 clone_col = "CloneID",
                                 hsc_types = c("HSC"),
                                 celltype_col = "celltype_final") {
  
  meta_sub <- meta %>%
    filter(
      !is.na(.data[[clone_col]]),
      .data[[clone_col]] != "0"
    ) %>%
    mutate(
      is_HSC = .data[[celltype_col]] %in% hsc_types
    )
  
  total_hsc    <- sum(meta_sub$is_HSC)
  total_nonhsc <- sum(!meta_sub$is_HSC)
  
  meta_sub %>%
    group_by(.data[[clone_col]]) %>%
    summarise(
      n_HSC = sum(is_HSC),
      n_nonHSC = sum(!is_HSC),
      HSC_freq = ifelse(total_hsc == 0, NA_real_, n_HSC / total_hsc),
      nonHSC_freq = ifelse(total_nonhsc == 0, NA_real_, n_nonHSC / total_nonhsc),
      .groups = "drop"
    )
}

calc_sister_clone_pool_freq <- function(
    seurat_A,
    seurat_B,
    label_A = "A",
    label_B = "B",
    clone_col = "CloneID",
    celltype_col = "celltype_final",
    hsc_types = c("HSC"),
    min_total_cells = 5
) {
  
  meta_A <- seurat_A@meta.data
  meta_B <- seurat_B@meta.data
  
  for (col in c(clone_col, celltype_col)) {
    if (!col %in% colnames(meta_A) || !col %in% colnames(meta_B)) {
      stop(glue("Column `{col}` must exist in both Seurat objects"))
    }
  }
  
  # ---------------------------
  # 1. 找 sister clones
  # ---------------------------
  sister_clones <- bind_rows(
    meta_A %>% select(.data[[clone_col]]) %>% mutate(donor = label_A),
    meta_B %>% select(.data[[clone_col]]) %>% mutate(donor = label_B)
  ) %>%
    filter(
      !is.na(.data[[clone_col]]),
      .data[[clone_col]] != "0"
    ) %>%
    distinct(.data[[clone_col]], donor) %>%
    group_by(.data[[clone_col]]) %>%
    summarise(n_donor = n_distinct(donor), .groups = "drop") %>%
    filter(n_donor == 2) %>%
    pull(.data[[clone_col]])
  
  # ---------------------------
  # 2. 先 subset sister clones
  # ---------------------------
  meta_A_sister <- meta_A %>%
    filter(.data[[clone_col]] %in% sister_clones)
  
  meta_B_sister <- meta_B %>%
    filter(.data[[clone_col]] %in% sister_clones)
  
  # ---------------------------
  # 3. 🔥 关键：先算 clone size 再 filter
  # ---------------------------
  clone_sizes <- bind_rows(
    meta_A_sister %>% count(.data[[clone_col]]),
    meta_B_sister %>% count(.data[[clone_col]])
  ) %>%
    group_by(.data[[clone_col]]) %>%
    summarise(total_cells = sum(n), .groups = "drop")
  
  kept_clones <- clone_sizes %>%
    filter(total_cells >= min_total_cells) %>%
    pull(.data[[clone_col]])
  
  # ---------------------------
  # 4. 🔥 用过滤后的 meta 再算 freq
  # ---------------------------
  meta_A_filtered <- meta_A_sister %>%
    filter(.data[[clone_col]] %in% kept_clones)
  
  meta_B_filtered <- meta_B_sister %>%
    filter(.data[[clone_col]] %in% kept_clones)
  
  freq_A <- calc_clone_pool_freq(meta_A_filtered,
                                 clone_col,
                                 hsc_types,
                                 celltype_col) %>%
    rename_with(~ paste0(.x, "_", label_A), -all_of(clone_col))
  
  freq_B <- calc_clone_pool_freq(meta_B_filtered,
                                 clone_col,
                                 hsc_types,
                                 celltype_col) %>%
    rename_with(~ paste0(.x, "_", label_B), -all_of(clone_col))
  
  # ---------------------------
  # 5. merge（不再需要再 filter）
  # ---------------------------
  full_join(freq_A, freq_B, by = clone_col)
}

freq_OO_OY <- calc_sister_clone_pool_freq(
  seurat.vivo.cross.oo,
  seurat.vivo.cross.oy,
  label_A = "OO",
  label_B = "OY"
)

freq_YO_YY <- calc_sister_clone_pool_freq(
  seurat.vivo.cross.yo,
  seurat.vivo.cross.yy,
  label_A = "YO",
  label_B = "YY"
)

global_max <- max(
  freq_OO_OY$max_hsc_freq,
  freq_YO_YY$max_hsc_freq,
  na.rm = TRUE
)

## ---- Old ----
freq_OO_OY <- freq_OO_OY %>%
  mutate(
    OA_OO = (nonHSC_freq_OO + eps) / (HSC_freq_OO + eps),
    OA_OY = (nonHSC_freq_OY + eps) / (HSC_freq_OY + eps),
    log2OA_OO = log2(OA_OO),
    log2OA_OY = log2(OA_OY),
    OA_ratio = (OA_OY + eps) / (OA_OO + eps),
    log2_OA_ratio = log2(OA_ratio),
    max_hsc_freq = pmax(HSC_freq_OO, HSC_freq_OY, na.rm = TRUE)
  )

## ---- Young ----
freq_YO_YY <- freq_YO_YY %>%
  mutate(
    OA_YO = (nonHSC_freq_YO + eps) / (HSC_freq_YO + eps),
    OA_YY = (nonHSC_freq_YY + eps) / (HSC_freq_YY + eps),
    log2OA_YO = log2(OA_YO),
    log2OA_YY = log2(OA_YY),
    OA_ratio = (OA_YY + eps) / (OA_YO + eps),
    log2_OA_ratio = log2(OA_ratio),
    max_hsc_freq = pmax(HSC_freq_YO, HSC_freq_YY, na.rm = TRUE)
  )

plot_clone_bias <- function(df, title_text, up_label, down_label) {
  
  df_plot <- df %>%
    mutate(
      rank = rank(log2_OA_ratio, ties.method = "first")
    )
  
  up_pct   <- round(mean(df_plot$log2_OA_ratio >=  2) * 100, 1)
  down_pct <- round(mean(df_plot$log2_OA_ratio <= -2) * 100, 1)
  
  ggplot(df_plot, aes(
    x    = rank,
    y    = log2_OA_ratio,
    size = max_hsc_freq,
    fill = log2_OA_ratio
  )) +
    geom_point(alpha = 0.6, colour = "black", pch = 21) +
    geom_hline(yintercept = c(2, -2), linetype = "dashed") +
    scale_size_continuous(
      name   = "HSC freq",
      limits = c(0, global_max),
      range  = c(2, 10),   # 推荐 2–10
      breaks = scales::pretty_breaks(n = 4)
    ) +
    colorspace::scale_fill_continuous_diverging(
      name   = expression(Log[2]~"FC"),
      limits = c(-3, 3),
      oob    = scales::squish
    ) +
    labs(
      title = title_text,
      x     = "Rank",
      y     = expression(Log[2]~"HSC self-renewal")
    ) +
    annotate(
      "text",
      x = min(df_plot$rank) + 3,
      y = max(df_plot$log2_OA_ratio, na.rm = TRUE) * 0.85,
      label = paste0(up_pct, "% clones\n", up_label),
      hjust = 0, vjust = 1, size = 5, color = "red"
    ) +
    annotate(
      "text",
      x = max(df_plot$rank) - 3,
      y = min(df_plot$log2_OA_ratio, na.rm = TRUE) * 0.85,
      label = paste0(down_pct, "% clones\n", down_label),
      hjust = 1, vjust = 0, size = 5, color = "blue"
    ) +
    theme_bw(base_size = 20)  +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}
plot_clone_bias <- function(df, title_text, up_label, down_label) {
  
  df_plot <- df %>%
    mutate(
      rank = rank(log2_OA_ratio, ties.method = "first")
    )
  
  up_pct   <- round(mean(df_plot$log2_OA_ratio >=  2) * 100, 1)
  down_pct <- round(mean(df_plot$log2_OA_ratio <= -2) * 100, 1)
  
  ggplot(df_plot, aes(
    x    = rank,
    y    = log2_OA_ratio,
    size = max_hsc_freq,
    fill = log2_OA_ratio
  )) +
    geom_point(alpha = 0.6, colour = "black", pch = 21) +
    geom_hline(yintercept = c(2, -2), linetype = "dashed") +
    scale_size_continuous(
      name   = "HSC freq",
      limits = c(0, global_max),
      range  = c(2, 10),   # 推荐 2–10
      breaks = scales::pretty_breaks(n = 4)
    ) +
    colorspace::scale_fill_continuous_diverging(
      name   = expression(bold(Log[2]~"FC")),
      limits = c(-3, 3),
      oob    = scales::squish
    ) +
    labs(
      title = title_text,
      x     = "Rank",
      y = expression(bold(Log[2]~"HSC self-renewal"))
    ) +
    annotate(
      "text",
      x = min(df_plot$rank) + 1,
      y = max(df_plot$log2_OA_ratio, na.rm = TRUE) * 0.95,
      label = paste0(up_pct, "% clones\n", up_label),
      hjust = 0, vjust = 1, size = 5, color = "red",
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = max(df_plot$rank) - 1,
      y = min(df_plot$log2_OA_ratio, na.rm = TRUE) * 0.95,
      label = paste0(down_pct, "% clones\n", down_label),
      hjust = 1, vjust = 0, size = 5, color = "blue",
      fontface = "bold"
    ) +
    theme_bw(base_size = 20)  +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      
      # 🔥 轴
      axis.title = element_text(size = 20, face = "bold"),
      axis.text  = element_text(size = 16, face = "bold"),
      
      # 🔥 legend
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14, face = "bold")
    )
}
p_old <- plot_clone_bias(
  freq_OO_OY,
  "OO vs OY",
  "HSC bias increased",
  "HSC bias decreased"
)

p_young <- plot_clone_bias(
  freq_YO_YY,
  "YO vs YY",
  "HSC bias increased",
  "HSC bias decreased"
)

print(p_old)
print(p_young)
# 520 370

library(patchwork)
p_old_clean <- p_old +
  theme(
    axis.title.y = element_blank()
  )
p_final <- (p_young | p_old_clean) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

p_final + plot_annotation(
  theme = theme(
    axis.title.y = element_text(size = 14)
  )
)


classify_clone_fate <- function(df) {
  df %>%
    mutate(
      Fate = case_when(
        log2_OA_ratio >=  2 ~ "High",
        log2_OA_ratio <= -2 ~ "Low",
        TRUE                ~ "Unchanged"
      )
    ) %>%
    count(Fate) %>%
    complete(
      Fate = c("High", "Unchanged", "Low"),
      fill = list(n = 0)
    ) %>%
    arrange(match(Fate, c("High", "Unchanged", "Low")))
}
old_counts   <- classify_clone_fate(freq_OO_OY)
young_counts <- classify_clone_fate(freq_YO_YY)
tab <- rbind(
  Old   = old_counts$n,
  Young = young_counts$n
)

colnames(tab) <- old_counts$Fate

tab
chisq.test(tab)
