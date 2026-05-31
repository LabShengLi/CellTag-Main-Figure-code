#######################################################################################


# Compute Relative Clone Size Distribution Plots and Replicate HSC clonal fraction plots 


#######################################################################################

# Maintainer: Chris Chen 
# Last updated: 02/04/2026

##########################

## Load packages ##

load_all_packages <- function() {
  pkgs <- c(
    "dplyr","tidyr","vegan","Seurat","ggplot2","tibble","stringr",
    "cowplot","purrr","ggrepel","harmony","patchwork","RColorBrewer",
    "scales","SingleR","celldex","EnhancedVolcano","gridExtra", 
    "readxl","Matrix","openxlsx","gt","glue"
  )
  
  suppressPackageStartupMessages(
    lapply(pkgs, require, character.only = TRUE)
  )
  message("Allpackages loaded.")
}
load_all_packages()

# set working directory
setwd('...')

################################

# read in CrossAge data (in vitro and in vivo)

seurat_in_vivo <- readRDS("data/CrossAge(exp2)_vivo.RDS")
seurat_in_vitro <- readRDS("data/CrossAge(exp2)_vitro.RDS")

# ------------------------------------------------------------
#  Define cell type groups (HSC/MPP and myeloid)
# ------------------------------------------------------------
simpleCTmap <- c(
  "LT-HSC" = "HSC/MPP",
  "ST-HSC" = "HSC/MPP",
  "MPP3"   = "HSC/MPP",
  "MPP4"   = "HSC/MPP"
)
myeloid_types <- c("GMP", "MkP", "EryP", "DC", "MAC", "Granulocyte")

# ------------------------------------------------------------
#  Extra metadata for OO and OY (sister clones only)
# ------------------------------------------------------------
meta_vivo <- seurat_in_vivo@meta.data %>%
  filter(sampleName %in% c("OO","OY")) %>%
  mutate(
    CloneID = as.character(CloneID),
    CellType = simpleCTmap[as.character(celltype)])

# ------------------------------------------------------------
#  Identify HSC/MPP clones
# ------------------------------------------------------------
hscmpp_clones <- meta_vivo %>%
  filter(CloneID != "0") %>%
  dplyr::count(sampleName, CloneID, name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from = sampleName,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  filter(OO > 1, OY > 1) %>%
  pull(CloneID)

# ------------------------------------------------------------
#  Compute relative clone frequencies for HSC/MPP 
# ------------------------------------------------------------
immature_OO <- meta_vivo %>% filter(sampleName == "OO") %>% mutate(CellType = simpleCTmap[as.character(celltype)])
immature_OY <- meta_vivo %>% filter(sampleName == "OY") %>% mutate(CellType = simpleCTmap[as.character(celltype)])

df_OO <- immature_OO %>% filter(CellType == "HSC/MPP", CloneID %in% hscmpp_clones) %>%
  dplyr::count(CloneID, name = "freq") %>% mutate(RelFreq = freq / sum(freq), Sample = "OO")

df_OY <- immature_OY %>% filter(CellType == "HSC/MPP", CloneID %in% hscmpp_clones) %>%
  dplyr::count(CloneID, name = "freq") %>% mutate(RelFreq = freq / sum(freq), Sample = "OY")

# ------------------------------------------------------------
#  Drop bottom 20% lowest-frequency clones
# ------------------------------------------------------------
drop_bottom_frac <- 0.20
filter_clones <- function(df, drop_bottom_frac = 0.2) {
  n_total <- nrow(df)
  keep_n  <- ceiling((1 - drop_bottom_frac) * n_total)
  df %>%
    arrange(desc(RelFreq)) %>%
    slice_head(n = keep_n) %>%
    mutate(CloneID = factor(CloneID, levels = CloneID[order(freq, decreasing = TRUE)]))
}
df_OO <- filter_clones(df_OO)
df_OY <- filter_clones(df_OY)

# ------------------------------------------------------------
#  Access difference in median, mean and max relative freq (OO vs OY)
# ------------------------------------------------------------

meta_hscmpp <- meta_vivo %>% mutate(CellType = simpleCTmap[as.character(celltype)],
  CloneID  = as.character(CloneID)) %>% 
  filter(CellType == "HSC/MPP", 
         CloneID %in% hscmpp_clones,
        sampleName %in% c("OO", "OY"),Rep != "UnMapped") %>%
  mutate(Sample = sampleName)

freq_hscmpp <- meta_hscmpp %>% dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>% mutate(RelFreq = freq / sum(freq)) %>% ungroup()

freq_hscmpp_ranked <- freq_hscmpp %>% group_by(Sample, Rep) %>%
  arrange(desc(RelFreq)) %>% mutate(Rank = row_number()) %>%
  filter(Rank <= 20) %>% ungroup()
freq_hscmpp_ranked

meta_myeloid <- meta_vivo %>% mutate(CloneID = as.character(CloneID)) %>%
  filter(
    celltype %in% myeloid_types,
    CloneID %in% hscmpp_clones,
    sampleName %in% c("OO", "OY"),
    Rep != "UnMapped"
  ) %>%
  mutate(Sample = sampleName)

freq_myeloid <- meta_myeloid %>% dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>% mutate(RelFreq = freq / sum(freq)) %>% ungroup()

freq_myeloid_ranked <- freq_myeloid %>%
  group_by(Sample, Rep) %>% arrange(desc(RelFreq)) %>%
  mutate(Rank = row_number()) %>% filter(Rank <= 20) %>% ungroup()
freq_myeloid_ranked

## Replicate level summary 

replicate_hsc_summary_OO_OY <- freq_hscmpp_ranked %>%
  group_by(Sample, Rep) %>%
  summarise(
    median_relfreq = median(RelFreq),
    mean_relfreq   = mean(RelFreq),
    max_relfreq    = max(RelFreq),
    .groups = "drop"
  )
replicate_hsc_summary_OO_OY

replicate_myeloid_summary_OO_OY <- freq_myeloid_ranked %>%
  group_by(Sample, Rep) %>%
  summarise(
    median_relfreq = median(RelFreq),
    mean_relfreq   = mean(RelFreq),
    max_relfreq    = max(RelFreq),
    .groups = "drop"
  )
replicate_myeloid_summary_OO_OY

# Wilcoxn rank sum 
wilcox.test(
  median_relfreq ~ Sample,
  data = replicate_hsc_summary_OO_OY,
  exact = FALSE
)
# max rel freq
wilcox.test(
  max_relfreq ~ Sample,
  data = replicate_hsc_summary_OO_OY,
  exact = FALSE
)

wilcox.test(
  median_relfreq ~ Sample,
  data = replicate_myeloid_summary_OO_OY,
  exact = FALSE
)
wilcox.test(
  mean_relfreq ~ Sample,
  data = replicate_myeloid_summary_OO_OY,
  exact = FALSE
)
#########################################

# ------------------------------------------------------------
# Plot per-sample distributions (OO & OY )
# ------------------------------------------------------------
# ------------------------------------------------------------
# Shared y-axis limits
# ------------------------------------------------------------
y_limit <- c(0, 0.15)
shared_clones <- intersect(df_OO$CloneID, df_OY$CloneID)
df_OO <- df_OO %>% filter(CloneID %in% shared_clones)
df_OY <- df_OY %>% filter(CloneID %in% shared_clones)

# ------------------------------------------------------------
# OO plot
# ------------------------------------------------------------
p_OO <- ggplot(df_OO, aes(x = CloneID, y = RelFreq)) +
  geom_col(fill = "#6A9FB5", width = 0.7) +
  scale_y_continuous(limits = y_limit, expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) + 
  theme_classic(base_size = 18) +   
  labs(
    x = "Clones",
    y = "Relative frequency",
    title = paste0("OO, n = ", nrow(df_OO))
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 16),
    plot.title   = element_text(hjust = 0.5, size = 22, face = "bold")
  )

# ------------------------------------------------------------
# OY plot
# ------------------------------------------------------------
p_OY <- ggplot(df_OY, aes(x = CloneID, y = RelFreq)) +
  geom_col(fill = "#E6AA68", width = 0.7) +
  scale_y_continuous(limits = y_limit, expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) + 
  theme_classic(base_size = 18) +
  labs(
    x = "Clones",
    y = "Relative frequency",
    title = paste0("OY, n = ", nrow(df_OY))
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 16),
    plot.title   = element_text(hjust = 0.5, size = 22, face = "bold")
  )
combined_plot_OO_OY <- grid.arrange(p_OO, p_OY, nrow = 1)

ggsave(
  filename = "/project2/sli68423_1316/users/chris/Celltag/HSC_heterogenecity/Publication_plots/Figure2_Clone_size_distribution/Supp/OO_OY_sister_clone_HSC_distribution.pdf",
  plot = combined_plot_OO_OY,
  width = 12,
  height = 5,
  units = "in"
)
###################

#### Replicate-aware ranked clone frequencies plot (HSC/MPP + Myeloid) ####

# HSC and MPP #
meta_hscmpp <- meta_vivo %>%
  mutate(CellType = simpleCTmap[as.character(celltype)]) %>%
  filter(CellType == "HSC/MPP",
         CloneID %in% hscmpp_clones,
         Rep != "UnMapped") %>%
  mutate(Sample = ifelse(sampleName == "OO", "OO", "OY")) 

freq_hscmpp <- meta_hscmpp %>%
  dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>%
  mutate(RelFreq = freq / sum(freq)) %>%
  ungroup()
top_n <- 20 # ranked top20 clones 
rep_levels <- c("OO_Rep1", "OO_Rep2", "OO_Rep3",
                "OY_Rep1", "OY_Rep2", "OY_Rep3")

freq_ranked <- freq_hscmpp %>%
  group_by(Sample, Rep) %>%
  arrange(desc(RelFreq)) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank <= top_n) %>%
  ungroup() %>%
  mutate(RepGroup = paste(Sample, Rep, sep = "_"),
         RepGroup = factor(RepGroup, levels = rep_levels),
         x_position = (Rank - 1) * 6 + as.numeric(RepGroup))
freq_ranked

# Myeloid #
meta_myeloid <- meta_vivo %>%
  filter(celltype %in% myeloid_types,
         CloneID %in% hscmpp_clones,
         Rep != "UnMapped") %>%
  mutate(Sample = ifelse(sampleName == "OO", "OO", "OY")) 

freq_myeloid <- meta_myeloid %>%
  dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>%
  mutate(RelFreq = freq / sum(freq)) %>%
  ungroup()

freq_myeloid_ranked <- freq_myeloid %>%
  group_by(Sample, Rep) %>%
  arrange(desc(RelFreq)) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank <= top_n) %>%
  ungroup() %>%
  mutate(RepGroup = paste(Sample, Rep, sep = "_"),
         RepGroup = factor(RepGroup, levels = rep_levels),
         x_position = (Rank - 1) * 6 + as.numeric(RepGroup))

## HSC_MPP ranked clonal fraction plot ##

p_HSCMPP_OO_OY <- ggplot(freq_ranked,
                   aes(x = x_position, y = RelFreq,
                       fill = Sample, group = interaction(Rank, Sample, Rep))) +
  geom_col(color = "black", width = 1.1, linewidth = 0.25) +
  scale_y_continuous(limits = c(0, 0.5),expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = c("OO" = "#cc246e", "OY" = "#661237"),
    name = "Sample"  
  ) +
  labs(x = "Clones by HSC and MPP output",
       y = "Fraction of\nall HSCs") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 12)),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    legend.position = "right",                
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.8, "cm")
  )
p_HSCMPP_OO_OY  

p_MYELOID_OO_OY <- ggplot(freq_myeloid_ranked,
                    aes(x = x_position, y = RelFreq,
                        fill = Sample, group = interaction(Rank, Sample, Rep))) +
  geom_col(color = "black", width = 1.1, linewidth = 0.25) +
  scale_y_continuous(limits = c(0, 0.5),expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = c("OO" = "#d9d98e", "OY" = "#808054"),
    name = "Sample"  
  ) +
  labs(x = "Clones by myeloid output",
       y = "Fraction of\nmyeloid cells") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 12)),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(0.8, "cm")
  )    
p_MYELOID_OO_OY 

combined_plot_OO_OY <- p_HSCMPP_OO_OY / p_MYELOID_OO_OY +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10))) &
  theme(legend.position = "right")

combined_plot_OO_OY

# save figure 2f
ggsave(
  filename = "OO_OY_HSC_MPP_output_ranked_plot.pdf",
  plot = combined_plot_OO_OY,
  width = 10,
  height = 6,
  units = "in"
)

##### Replicate level plots for YY and YO ######

# ------------------------------------------------------------
# 1. Define HSC/MPP group and subset for YY / YO
# ------------------------------------------------------------
meta_hscmpp <- seurat_in_vivo@meta.data %>%
  mutate(
    CellGroup = recode(celltype,
                       "LT-HSC" = "HSC/MPP",
                       "ST-HSC" = "HSC/MPP",
                       "MPP3"   = "HSC/MPP",
                       "MPP4"   = "HSC/MPP"),
    CloneID = as.character(CloneID)
  ) %>%
  filter(
    CellGroup == "HSC/MPP",
    sampleName %in% c("YY", "YO"),
    Rep != "UnMapped",
    CloneID != "0"
  )

sister_clones <- seurat_in_vivo@meta.data %>%
  filter(sampleName %in% c("YY", "YO"),
         Rep != "UnMapped",
         CloneID != "0") %>%
  mutate(CloneID = as.character(CloneID)) %>%
  dplyr::count(sampleName, CloneID, name = "n_cells") %>%
  pivot_wider(names_from = sampleName, values_from = n_cells, values_fill = 0) %>%
  filter(YY >= 1, YO >= 1) %>%
  pull(CloneID)

length(sister_clones) # 93 sister clones 

simpleCTmap <- c(
  "LT-HSC" = "HSC/MPP",
  "ST-HSC" = "HSC/MPP",
  "MPP3"   = "HSC/MPP",
  "MPP4"   = "HSC/MPP"
)
myeloid_types <- c("GMP", "MkP", "EryP", "DC", "MAC", "Granulocyte")

meta_vivo <- seurat_in_vivo@meta.data %>%
  filter(sampleName %in% c("YY", "YO")) %>%
  mutate(
    CloneID  = as.character(CloneID),
    CellType = simpleCTmap[as.character(celltype)]
  ) %>%
  filter(CloneID %in% sister_clones)

# Identify HSC/MPP clones present in both YY & YO

hscmpp_clones <- meta_vivo %>%
  filter(CellType == "HSC/MPP") %>%
  dplyr::count(sampleName, CloneID, name = "n_cells") %>%
  pivot_wider(
    names_from  = sampleName,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  filter(YY >= 1, YO >= 1) %>%
  pull(CloneID)


# Compute clone frequencies for YY and YO

immature_YY <- meta_vivo %>%
  filter(sampleName == "YY") %>%
  mutate(CellType = simpleCTmap[as.character(celltype)])

immature_YO <- meta_vivo %>%
  filter(sampleName == "YO") %>%
  mutate(CellType = simpleCTmap[as.character(celltype)])

df_YY <- immature_YY %>%
  filter(CellType == "HSC/MPP", CloneID %in% hscmpp_clones) %>%
  dplyr::count(CloneID, name = "freq") %>%
  mutate(RelFreq = freq / sum(freq), Sample = "YY")

df_YO <- immature_YO %>%
  filter(CellType == "HSC/MPP", CloneID %in% hscmpp_clones) %>%
  dplyr::count(CloneID, name = "freq") %>%
  mutate(RelFreq = freq / sum(freq), Sample = "YO")

filter_clones <- function(df) {
  df %>%
    arrange(desc(RelFreq)) %>%
    mutate(CloneID = factor(CloneID, levels = CloneID[order(RelFreq, decreasing = TRUE)]))
}

df_YY <- filter_clones(df_YY)
df_YO <- filter_clones(df_YO)

# ------------------------------------------------------------
# YY plot
# ------------------------------------------------------------
p_YY <- ggplot(df_YY, aes(x = CloneID, y = RelFreq)) +
  geom_col(fill = "#6A9FB5", width = 0.7) +
  scale_y_continuous(limits = c(0, 0.12), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 18) +
  labs(
    x = "Clones",
    y = "Relative frequency",
    title = paste0("YY, n = ", nrow(df_YY))
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 16),
    plot.title   = element_text(hjust = 0.5, size = 22, face = "bold")
  )

# ------------------------------------------------------------
# YO plot
# ------------------------------------------------------------
p_YO <- ggplot(df_YO, aes(x = CloneID, y = RelFreq)) +
  geom_col(fill = "#E6AA68", width = 0.7) +
  scale_y_continuous(limits = c(0, 0.12), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 18) +
  labs(
    x = "Clones",
    y = "Relative frequency",
    title = paste0("YO, n = ", nrow(df_YO))
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 16),
    plot.title   = element_text(hjust = 0.5, size = 22, face = "bold")
  )
combined_plot_YY_YO <- grid.arrange(p_YY, p_YO, nrow = 1)


meta_hscmpp <- meta_vivo %>%
  mutate(CellType = simpleCTmap[as.character(celltype)]) %>%
  filter(CellType == "HSC/MPP",
         CloneID %in% hscmpp_clones,
         Rep != "UnMapped") %>%
  mutate(Sample = sampleName)   # YY / YO

freq_hscmpp <- meta_hscmpp %>%
  dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>%
  mutate(RelFreq = freq / sum(freq)) %>%
  ungroup()

# Extract the top 20 clones 
top_n <- 20
rep_levels <- c("YY_Rep1","YY_Rep2","YY_Rep3",
                "YO_Rep1","YO_Rep2","YO_Rep3")

freq_ranked <- freq_hscmpp %>%
  group_by(Sample, Rep) %>%
  arrange(desc(RelFreq)) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank <= top_n) %>%
  ungroup() %>%
  mutate(
    RepGroup = paste(Sample, Rep, sep = "_"),
    RepGroup = factor(RepGroup, levels = rep_levels),
    x_position = (Rank - 1) * 6 + as.numeric(RepGroup)
  )

meta_myeloid <- meta_vivo %>%
  filter(celltype %in% myeloid_types,
         CloneID %in% hscmpp_clones,
         Rep != "UnMapped") %>%
  mutate(Sample = sampleName)

freq_myeloid <- meta_myeloid %>%
  dplyr::count(Sample, Rep, CloneID, name = "freq") %>%
  group_by(Sample, Rep) %>%
  mutate(RelFreq = freq / sum(freq)) %>%
  ungroup()

freq_myeloid_ranked <- freq_myeloid %>%
  group_by(Sample, Rep) %>%
  arrange(desc(RelFreq)) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank <= top_n) %>%
  ungroup() %>%
  mutate(
    RepGroup = paste(Sample, Rep, sep = "_"),
    RepGroup = factor(RepGroup, levels = rep_levels),
    x_position = (Rank - 1) * 6 + as.numeric(RepGroup)
  )

# Create clonal fraction plots #
p_HSCMPP_YY_YO <- ggplot(freq_ranked,
                   aes(x = x_position, y = RelFreq,
                       fill = Sample, group = interaction(Rank, Sample, Rep))) +
  geom_col(color = 'black', width = 1.1, linewidth = 0.25) +
  scale_fill_manual(values = c("YY"="#cc246e", "YO"="#661237"), breaks = c("YY", "YO") ) +
  scale_y_continuous(limits = c(0, 0.5),expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Clones by HSC and MPP output", 
       y = "Fraction of\nall HSCs") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 12)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12)
  )

p_MYELOID_YY_YO <- ggplot(freq_myeloid_ranked,
                    aes(x = x_position, y = RelFreq,
                        fill = Sample, group = interaction(Rank, Sample, Rep))) +
  geom_col(color = 'black', width = 1.1, linewidth = 0.25) +
  scale_fill_manual(values = c("YY"="#d9d98e", "YO"="#808054"), breaks = c("YY", "YO") ) +
  scale_y_continuous(limits = c(0, 0.5),expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Clones by myeloid output", 
       y = "Fraction of\nmyeloid cells") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 12)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12)
  )

combined_plot_YY_YO <- p_HSCMPP_YY_YO / p_MYELOID_YY_YO +
  plot_annotation(theme = theme(plot.margin = margin(10,10,10,10))) &
  theme(legend.position = "right")

combined_plot_YY_YO

### Compute statistical p values ###

replicate_hsc_summary <- freq_ranked %>%
  group_by(Sample, Rep) %>%
  summarise(
    median_relfreq = median(RelFreq),
    mean_relfreq   = mean(RelFreq),
    max_relfreq    = max(RelFreq),
    .groups = "drop"
  )
replicate_hsc_summary
wilcox.test(
  median_relfreq ~ Sample,
  data = replicate_hsc_summary,
  exact = FALSE
) # p = 0.08
wilcox.test(
  max_relfreq ~ Sample,
  data = replicate_hsc_summary,
  exact = FALSE
) # p = 0.19
####### Myeloid #######

replicate_myeloid_summary <- freq_myeloid_ranked %>%
  group_by(Sample, Rep) %>%
  summarise(
    median_relfreq = median(RelFreq),
    mean_relfreq   = mean(RelFreq),
    max_relfreq    = max(RelFreq),
    .groups = "drop"
  )

replicate_myeloid_summary

wilcox.test(
  median_relfreq ~ Sample,
  data = replicate_myeloid_summary,
  exact = FALSE
) # p = 0.66
wilcox.test(
  max_relfreq ~ Sample,
  data = replicate_myeloid_summary,
  exact = FALSE
) # p = 0.66

#############################################################

#### Figure 2e relative clonal size stacked bar chart ####

# count and create count labels for Y and O vitro clones 

label_counts_vitro <- seurat_in_vitro@meta.data %>%
  dplyr::filter(sampleName %in% c("Y_vitro","O_vitro"),Rep != "UnMapped", CloneID != "0"
  ) %>% dplyr::count(sampleName, Rep, CloneID, name = "n_cells") %>%
  dplyr::filter(n_cells >= 1) %>% dplyr::count(sampleName, Rep, name = "n_clones") %>%
  dplyr::mutate(Sample = sampleName, SampleLabel = paste0(sampleName, " (", n_clones, ")")
  ) %>% dplyr::select(Replicate = Rep, Sample, SampleLabel) %>%
  dplyr::arrange(Sample, Replicate)
label_counts_vitro

#### Y_Vitro and O_vitro ####

# helper function to create equal width bar 
add_dummy_facets <- function(df, facet_col, target_levels) {
  missing_levels <- setdiff(target_levels, unique(df[[facet_col]]))
  if (length(missing_levels) == 0) {
    return(df)
  }
  dummy <- df[rep(1, length(missing_levels)), , drop = FALSE]
  dummy[,] <- NA
  dummy[[facet_col]] <- missing_levels
  dplyr::bind_rows(df, dummy)
}

plot_clone_distribution_exact <- function(
    seurat_obj,
    sample_name,
    label_counts,
    rep_col = "Rep",
    clone_col = "CloneID",
    relfreq_threshold = 0.003,
    facet_levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6") 
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(ggsci)
    library(colorspace)
  })
  other_lab <- paste0("Other (<", relfreq_threshold * 100, "%)")
  meta_df <- seurat_obj@meta.data %>%
    filter(
      sampleName == sample_name,
      !!sym(rep_col) != "UnMapped",
      !!sym(clone_col) != "0"
    ) %>%
    mutate(
      RepLabel = !!sym(rep_col),
      CloneID  = as.character(!!sym(clone_col))
    )
  clone_freq <- meta_df %>%
    dplyr::count(RepLabel, CloneID, name = "n_cells") %>%
    dplyr::filter(n_cells >= 1) %>%
    group_by(RepLabel) %>%
    mutate(RelFreq = n_cells / sum(n_cells)) %>%
    ungroup()
  clone_freq <- clone_freq %>%
    mutate(CloneLabel = ifelse(RelFreq < relfreq_threshold, other_lab, CloneID)) %>%
    group_by(RepLabel, CloneLabel) %>%
    summarise(RelFreq = sum(RelFreq), .groups = "drop")
  clone_freq <- clone_freq %>%
    left_join(
      label_counts %>%
        filter(Sample == sample_name) %>%
        dplyr::select(Replicate, SampleLabel),
      by = c("RepLabel" = "Replicate")
    ) %>%
    mutate(x_bar = SampleLabel)
  clone_order <- clone_freq %>%
    group_by(CloneLabel) %>%
    summarise(total = sum(RelFreq)) %>%
    arrange(desc(total)) %>%
    pull(CloneLabel)
  clone_freq$CloneLabel <- factor(
    clone_freq$CloneLabel,
    levels = c(other_lab, setdiff(clone_order, other_lab))
  )
  set.seed(42)
  pal_raw <- c(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Accent"),
    pal_npg("nrc")(10),
    pal_d3("category20")(20)
  )
  hsv_vals <- coords(as(hex2RGB(pal_raw), "HSV"))
  pal_raw <- pal_raw[hsv_vals[, "S"] > 0.4]
  pal <- rep(pal_raw, length.out = nlevels(clone_freq$CloneLabel))
  names(pal) <- levels(clone_freq$CloneLabel)
  pal[other_lab] <- "grey80"
  clone_freq <- add_dummy_facets(
    clone_freq,
    facet_col = "RepLabel",
    target_levels = facet_levels
  )
  clone_freq$RepLabel <- factor(clone_freq$RepLabel, levels = facet_levels)
  ggplot(
    clone_freq,
    aes(x = x_bar, y = RelFreq, fill = CloneLabel)
  ) +
    geom_bar(
      stat = "identity",
      width = 0.9,
      position = position_stack(reverse = TRUE)
    ) +
    facet_wrap(~RepLabel, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Relative Clone Size Distribution: ", sample_name),
      x = NULL,
      y = "Relative clone size"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      axis.ticks.x = element_blank(),
      strip.text   = element_text(size = 20),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title   = element_text(size = 18, hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}
# Apply the function 
p_Y_vitro <- plot_clone_distribution_exact(
  seurat_obj = seurat_in_vitro,
  sample_name = "Y_vitro",
  label_counts = label_counts_vitro
)
p_Y_vitro
p_O_vitro <- plot_clone_distribution_exact(
  seurat_obj = seurat_in_vitro,
  sample_name = "O_vitro",
  label_counts = label_counts_vitro
)
p_O_vitro

#################

##### Repeat for in Vivo OY vs OO and YY vs YO ######

plot_relative_clone_size_distribution <- function(
    seurat_obj,
    sample_names = c("YY", "YO"),
    facet_samples = c("YY", "YO"),    
    facet_reps = c("Rep1", "Rep2", "Rep3"), 
    rep_col = "Rep",
    clone_col = "CloneID",
    relfreq_threshold = 0.01,
    label_counts,
    title = "Relative Clone Size Distribution"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(ggsci)
    library(colorspace)
  })
  set.seed(42)
  other_lab <- paste0("Other (<", relfreq_threshold * 100, "%)")
  meta_df <- seurat_obj@meta.data %>%
    filter(
      sampleName %in% sample_names,
      !!sym(rep_col) != "UnMapped",
      !!sym(clone_col) != "0"
    ) %>%
    mutate(
      Sample    = factor(sampleName, levels = sample_names),
      Replicate = !!sym(rep_col),
      CloneID   = as.character(!!sym(clone_col))
    )
  clone_freq <- meta_df %>%
    dplyr::count(Sample, Replicate, CloneID, name = "n_cells") %>%
    group_by(Sample, Replicate) %>%
    mutate(RelFreq = n_cells / sum(n_cells)) %>%
    ungroup()
  clone_freq <- clone_freq %>%
    mutate(CloneLabel = ifelse(RelFreq < relfreq_threshold, other_lab, CloneID)) %>%
    group_by(Sample, Replicate, CloneLabel) %>%
    summarise(RelFreq = sum(RelFreq), .groups = "drop")
  facet_levels <- unlist(
    lapply(facet_samples, function(s) paste(s, facet_reps, sep = "_"))
  )
  clone_freq <- clone_freq %>%
    left_join(label_counts, by = c("Replicate", "Sample")) %>%
    mutate(
      FacetOrder = factor(
        paste(Sample, Replicate, sep = "_"),
        levels = facet_levels
      )
    )
  clone_order <- clone_freq %>%
    group_by(CloneLabel) %>%
    summarise(total = sum(RelFreq)) %>%
    arrange(desc(total)) %>%
    pull(CloneLabel)
  clone_freq$CloneLabel <- factor(
    clone_freq$CloneLabel,
    levels = c(other_lab, setdiff(clone_order, other_lab))
  )
  pal_raw <- c(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Accent"),
    pal_npg("nrc")(10),
    pal_d3("category20")(20)
  )
  hsv_vals <- coords(as(hex2RGB(pal_raw), "HSV"))
  pal_raw <- pal_raw[hsv_vals[, "S"] > 0.4]
  pal <- rep(pal_raw, length.out = nlevels(clone_freq$CloneLabel))
  names(pal) <- levels(clone_freq$CloneLabel)
  pal[other_lab] <- "grey80"
  ggplot(
    clone_freq,
    aes(
      x = SampleLabel,
      y = RelFreq,
      fill = CloneLabel
    )
  ) +
    geom_bar(
      stat = "identity",
      width = 0.9,
      position = position_stack(reverse = TRUE)
    ) +
    facet_wrap(
      ~ FacetOrder,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(
        FacetOrder = function(x) sub(".*_", "", x)  
      )
    ) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Relative clone size"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x  = element_text(size = 16),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 18),
      axis.title.y = element_text(size = 20),
      strip.text   = element_text(size = 20),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

## OO vs OY ##
label_counts_OY_OO <- seurat_in_vivo@meta.data %>%
  dplyr::filter(sampleName %in% c("OY","OO"),Rep != "UnMapped",CloneID != "0"
  ) %>% dplyr::count(sampleName, Rep, CloneID, name = "n_cells") %>%
  dplyr::filter(n_cells >= 2) %>% dplyr::count(sampleName, Rep, name = "n_clones") %>%
  dplyr::mutate(Sample = sampleName, SampleLabel = paste0(sampleName, " (", n_clones, ")")
  ) %>% dplyr::select(Replicate = Rep, Sample, SampleLabel) %>%
  dplyr::arrange(Sample, Replicate)

label_counts_OY_OO

p_OY_OO <- plot_relative_clone_size_distribution(
  seurat_obj = seurat_in_vivo,
  sample_names  = c("OY", "OO"),
  facet_samples = c("OY", "OO"),   
  facet_reps    = c("Rep1", "Rep2", "Rep3"),
  rep_col = "Rep",
  clone_col = "CloneID",
  relfreq_threshold = 0.01,
  label_counts = label_counts_OY_OO,
  title = "Relative Clone Size Distribution (OY vs OO)"
)
p_OY_OO

## YY vs YO ##
label_counts_YY_YO <- seurat_in_vivo@meta.data %>%
  dplyr::filter(sampleName %in% c("YY","YO"), Rep != "UnMapped",CloneID != "0"
  ) %>% dplyr::count(sampleName, Rep, CloneID, name = "n_cells") %>%
  dplyr::filter(n_cells >= 2) %>% dplyr::count(sampleName, Rep, name = "n_clones") %>%
  dplyr::mutate(Sample = sampleName, SampleLabel = paste0(sampleName, " (", n_clones, ")")
  ) %>% dplyr::select(Replicate = Rep, Sample, SampleLabel) %>%
  dplyr::arrange(Sample, Replicate)

label_counts_YY_YO


p_YY_YO <- plot_relative_clone_size_distribution(
  seurat_obj = seurat_in_vivo,
  sample_names = c("YY","YO"),
  facet_samples = c("YY","YO"),   
  label_counts = label_counts_YY_YO,
  title = "Relative Clone Size Distribution (YY vs YO)"
)
p_YY_YO

### Now make them the same scale 

p_Y_vitro
p_O_vitro
p_OY_OO
p_YY_YO

combined_panel <-
  (p_Y_vitro | p_O_vitro) /
  (p_YY_YO   | p_OY_OO) +
  plot_layout(
    widths  = c(1, 1),
    heights = c(1, 1),
    guides  = "collect"
  ) &
  theme(
    legend.position = "none"
  )
combined_panel
getwd()
ggsave(
  filename = "Figure2_Relative_clone_size_distributions.pdf",
  plot     = combined_panel,
  device   = cairo_pdf,
  width    = 18,
  height   = 10,
  units    = "in",
  dpi      = 300
)
#############################################################

# Clonal distribution stats warpper function 

run_clone_distribution_engine <- function(
    comparison_list,
    output_file = NULL,
    rep_col = "Rep",
    clone_col = "CloneID"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(ggsci)
    library(colorspace)
    library(patchwork)
  })
  # --------------------------
  # Palette builder
  # --------------------------
  build_palette <- function(levels_vec, other_lab) {
    pal_raw <- c(
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Dark2"),
      brewer.pal(8, "Accent"),
      pal_npg("nrc")(10),
      pal_d3("category20")(20)
    )
    hsv_vals <- coords(as(hex2RGB(pal_raw), "HSV"))
    pal_raw <- pal_raw[hsv_vals[, "S"] > 0.4]
    pal <- rep(pal_raw, length.out = length(levels_vec))
    names(pal) <- levels_vec
    pal[other_lab] <- "grey80"
    pal
  }
  # --------------------------
  # Stats computation
  # --------------------------
  compute_stats <- function(meta_df) {
    meta_df %>%
      filter(!!sym(clone_col) != "0") %>%
      group_by(sampleName, !!sym(clone_col)) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(sampleName) %>%
      mutate(freq = n / sum(n)) %>%
      summarise(
        n_clones = n(),
        min_freq = min(freq),
        max_freq = max(freq),
        mean_freq = mean(freq),
        .groups = "drop"
      )
  }
  plot_list <- list()
  stats_list <- list()
  # ======================================================
  # LOOP THROUGH USER-DEFINED COMPARISONS
  # ======================================================
  for (comp in comparison_list) {
    seurat_obj  <- comp$seurat_obj
    sample_names <- comp$samples
    label_counts <- comp$label_counts
    threshold    <- comp$threshold
    title        <- comp$title
    comp_name    <- comp$name
    other_lab <- paste0("Other (<", threshold * 100, "%)")
    meta_df <- seurat_obj@meta.data %>%
      filter(sampleName %in% sample_names,
             !!sym(rep_col) != "UnMapped",
             !!sym(clone_col) != "0") %>%
      mutate(
        Sample = sampleName,
        Replicate = !!sym(rep_col),
        CloneID = as.character(!!sym(clone_col))
      )
    clone_freq <- meta_df %>%
      count(Sample, Replicate, CloneID, name = "n_cells") %>%
      group_by(Sample, Replicate) %>%
      mutate(RelFreq = n_cells / sum(n_cells)) %>%
      ungroup() %>%
      mutate(CloneLabel = ifelse(RelFreq < threshold, other_lab, CloneID)) %>%
      group_by(Sample, Replicate, CloneLabel) %>%
      summarise(RelFreq = sum(RelFreq), .groups = "drop") %>%
      left_join(label_counts, by = c("Replicate","Sample")) %>%
      mutate(FacetOrder = paste(Sample, Replicate, sep = "_"))
    clone_order <- clone_freq %>%
      group_by(CloneLabel) %>%
      summarise(total = sum(RelFreq)) %>%
      arrange(desc(total)) %>%
      pull(CloneLabel)
    clone_freq$CloneLabel <- factor(
      clone_freq$CloneLabel,
      levels = c(other_lab, setdiff(clone_order, other_lab))
    )
    pal <- build_palette(levels(clone_freq$CloneLabel), other_lab)
    p <- ggplot(clone_freq,
                aes(x = SampleLabel,
                    y = RelFreq,
                    fill = CloneLabel)) +
      geom_bar(stat = "identity",
               width = 0.9,
               position = position_stack(reverse = TRUE)) +
      facet_wrap(~FacetOrder, nrow = 1) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = title,
           y = "Relative clone size",
           x = NULL) +
      theme_classic() +
      theme(legend.position = "none")
    plot_list[[comp_name]] <- p
    stats_list[[comp_name]] <- compute_stats(seurat_obj@meta.data)
  }
  # --------------------------
  # Combine all panels
  # --------------------------
  combined_panel <- wrap_plots(plot_list)
  if (!is.null(output_file)) {
    ggsave(output_file,
           combined_panel,
           width = 18,
           height = 10,
           device = cairo_pdf)
  }
  return(list(
    plots = plot_list,
    combined_plot = combined_panel,
    descriptive_stats = stats_list
  ))
}

### Example usage ###

comparison_list <- list(
  
  list(
    name = "Vitro_Y_vs_O",
    seurat_obj = seurat_in_vitro,
    samples = c("Y_vitro","O_vitro"),
    label_counts = label_counts_vitro,
    threshold = 0.002,
    title = "Y vs O (In Vitro)"
  ),
  
  list(
    name = "OY_vs_OO",
    seurat_obj = seurat_in_vivo,
    samples = c("OY","OO"),
    label_counts = label_counts_OY_OO,
    threshold = 0.01,
    title = "OY vs OO"
  ),
  
  list(
    name = "YY_vs_YO",
    seurat_obj = seurat_in_vivo,
    samples = c("YY","YO"),
    label_counts = label_counts_YY_YO,
    threshold = 0.01,
    title = "YY vs YO"
  )
)

res <- run_clone_distribution_engine(
  comparison_list = comparison_list,
  output_file = "Figure2_clone_distribution.pdf"
)

res$combined_plot
res$descriptive_stats



############################################################

# End of this script 
