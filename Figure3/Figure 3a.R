# =========================
# Load libraries
# =========================
library(Seurat)
library(ggplot2)
library(scuttle)
library(ggpubr)

# =========================
# Step 1: Load data
# =========================
infn <- "/project2/sli68423_1316/users/Qiuyang/Qiuyang_Zhang/cell_tag/Celltag_main_scripts/Main_figures/Data_objects/CrossAge_vitro.RDS"
s <- readRDS(infn)

# =========================
# Step 2: Extract raw counts
# =========================
raw_counts <- GetAssayData(s, layer = "counts")

# =========================
# Step 3: Determine target library size (peak of density)
# =========================
lib_size <- Matrix::colSums(raw_counts)
d <- density(lib_size)
target_library_size <- round(d$x[which.max(d$y)])

cat("Target library size:", target_library_size, "\n")

# =========================
# Step 4: Downsample
# =========================
set.seed(123)

prop <- pmin(1, target_library_size / lib_size)
down_counts <- downsampleMatrix(raw_counts, prop = prop)

# =========================
# Step 5: Compute Shannon entropy
# =========================
entropy_raw <- apply(down_counts, 2, function(cell) {
  p <- cell / sum(cell)
  p <- p[p > 0]  
  -sum(p * log(p))
})

# =========================
# Step 6: Normalize entropy
# =========================
max_entropy <- log(colSums(down_counts > 0))
scaled_entropy <- entropy_raw / max_entropy

# Add to metadata
s$shannon_entropy <- scaled_entropy

# =========================
# Step 7: Define groups
# =========================
s$Group <- ifelse(grepl("Y_vitro", s$orig.ident), "Young",
                  ifelse(grepl("O_vitro", s$orig.ident), "Old", NA))

# =========================
# Step 8: Subset HSC cells
# =========================
s_LTHSC <- subset(s, subset = celltype_final == "HSC" & !is.na(Group))

# =========================
# Step 9: Prepare plotting data
# =========================
df_entropy <- s_LTHSC@meta.data[, c("Group", "shannon_entropy")]

df_entropy$Group <- factor(df_entropy$Group, levels = c("Young", "Old"))

cell_counts <- table(df_entropy$Group)

# =========================
# Step 10: Plot
# =========================
p <- ggplot(df_entropy, aes(x = Group, y = shannon_entropy, fill = Group)) +
  
  geom_boxplot(width = 0.4, color = "black", alpha = 0.9) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "black"
  ) +
  
  stat_compare_means(
    comparisons = list(c("Young", "Old")),
    method = "t.test",
    label = "p.signif",
    size = 6,
    tip.length = 0.01,
    label.y = max(df_entropy$shannon_entropy) * 0.98
  ) +
  
  scale_fill_manual(values = c(
    "Young" = "#88CCEE",
    "Old"   = "#DDCC77"
  )) +
  
  labs(x = "Age", y = "Normalized Shannon Entropy") +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  
  coord_cartesian(ylim = c(0.90, 0.96))

print(p)

# =========================
# Step 11: Save
# =========================
#ggsave("shannon_entropy_vitro_HSC.pdf", p, width = 4, height = 6, dpi = 600)