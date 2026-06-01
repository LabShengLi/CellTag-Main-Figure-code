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
infn <- "data/CrossAge(exp2)_vivo.RDS"
s <- readRDS(infn)

# =========================
# Step 2: Define Group (4 groups)
# =========================
s$Group <- ifelse(s$orig.ident %in% c("YY","YO","OY","OO"), s$orig.ident, NA)
s$Group <- factor(s$Group, levels = c("YY", "YO", "OY", "OO"))

# =========================
# Step 3: Downsample counts
# =========================
raw_counts <- GetAssayData(s, layer = "counts")
lib_size <- Matrix::colSums(raw_counts)

d <- density(lib_size)
target_library_size <- round(d$x[which.max(d$y)])
cat("Target library size:", target_library_size, "\n")

# Downsample
set.seed(123)
prop <- pmin(1, target_library_size / lib_size)
down_counts <- downsampleMatrix(raw_counts, prop = prop)

# =========================
# Step 4: Shannon entropy
# =========================
entropy_raw <- apply(down_counts, 2, function(cell) {
  p <- cell / sum(cell)
  p <- p[p > 0]
  -sum(p * log(p))
})

max_entropy <- log(colSums(down_counts > 0))
scaled_entropy <- entropy_raw / max_entropy

s$shannon_entropy <- scaled_entropy

# =========================
# Step 5: Subset HSC
# =========================
s_HSC <- subset(s, subset = celltype_final == "HSC" & !is.na(Group))

# =========================
# Step 6: Prepare plotting data
# =========================
df_entropy <- s_HSC@meta.data[, c("Group", "shannon_entropy")]
df_entropy$Group <- factor(df_entropy$Group, levels = c("YY","YO","OY","OO"))

cell_counts <- table(df_entropy$Group)

# =========================
# Step 7: Plot
# =========================
comparisons <- list(
  c("OO", "OY"),
  c("OO", "YO"),
  c("OO", "YY"),
  c("OY", "YO"),
  c("OY", "YY")
)

p <- ggplot(df_entropy, aes(x = Group, y = shannon_entropy, fill = Group)) +
  
  geom_boxplot(width = 0.5, color = "black", alpha = 0.9, outlier.shape = NA) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "black"
  ) +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    label = "p.signif",
    size = 8,
    tip.length = 0.01,
    label.y = 0.92,
    step.increase = 0.05
  ) +
  
  scale_fill_manual(values = c(
    "YY" = "#88CCEE",
    "YO" = "#44AA99",
    "OY" = "#DDCC77",
    "OO" = "#882255"
  )) +
  
  labs(
    x = "Donor-Age & Host Age",
    y = "Single Cell Shannon Entropy"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  
  coord_cartesian(ylim = c(0.86, 0.96))

print(p)

# =========================
# Step 8: Save
# =========================
#ggsave("shannon_entropy_invivo_HSC.pdf", p, width = 5, height = 5, dpi = 600)
