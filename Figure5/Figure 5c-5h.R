############################
# 1. Load packages
############################
library(readr)
library(dplyr)
library(ggplot2)

############################
# 2. Load data
############################

# ML results
ml_res <- read_csv("data/ml_results_summary.csv")

# Feature importance (vitro)
DVG_33_O_vitro_OO <- read_csv("data/DVG_33_Old_SF_OO_Top50_Importance.csv")
DVG_33_Y_vitro_YY <- read_csv("data/DVG_33_Young_SF_YY_Top50_Importance.csv")
DVG_80_O_vitro_OO <- read_csv("data/DVG_80_Old_SF_OO_Top50_Importance.csv")
DVG_union_O_vitro_OO <- read_csv("data/DVG_union_Old_SF_OO_Top50_Importance.csv")
DEG_O_vitro_OO    <- read_csv("data/LowOutput_Old_SF_OO_Top50_Importance.csv")

# Feature importance (vivo)
DVG_33_OO <- read_csv("data/DVG_33_OO_Top50_Importance.csv")
DVG_33_YY <- read_csv("data/DVG_33_YY_Top50_Importance.csv")
DVG_80_OO <- read_csv("data/DVG_80_OO_Top50_Importance.csv")
DVG_union_OO <- read_csv("data/DVG_union_OO_Top50_Importance.csv")
DEG_OO    <- read_csv("data/LowOutput_OO_Top50_Importance.csv")

# Additional datasets
DVG_33_Up_O_vitro_OO <- read_csv("data/DVG_33_Old_SF_OO_Top10_Up_in_High.csv")
DVG_80_OO <- read_csv("data/DVG_80_OO_Top50_Importance.csv")

############################
# 3. Function: F1 bar plot
############################
plot_f1_bar <- function(data, dataset_type, y_label,
                        gene_subset = c("LowOutput", "DVG_33", "DVG_80", "DVG_union"),
                        base_size = 14, show_legend = TRUE) {
  
  df <- data %>%
    filter(dataset == dataset_type,
           gene_list %in% gene_subset) %>%
    mutate(
      gene_list = recode(gene_list,
                         LowOutput = "DEG",
                         DVG_33 = "DVG33",
                         DVG_80 = "DVG80",
                         DVG_union = "DVGunion")
    )
  
  p <- ggplot(df, aes(
    x = condition,
    y = macro_F1,
    fill = gene_list
  )) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.75),
      width = 0.65
    ) +
    scale_fill_manual(values = c(
      DEG = "#1f77b4",
      DVG33 = "#e74c3c",
      DVG80 = "#2ecc71",
      DVGunion = "#9b59b6"
    )) +
    labs(
      x = "Condition",
      y = y_label,
      fill = "Gene List"
    ) +
    theme_classic(base_size = base_size)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

############################
# 4. Function: Top20 importance plot
############################
plot_top20_as_bar <- function(dat, title = NULL, font_size = 14) {
  
  colnames(dat)[1:2] <- c("Genes", "Overall")
  
  dat <- dat %>%
    arrange(desc(Overall)) %>%
    slice(1:20)
  
  ggplot(dat, aes(
    x = reorder(Genes, -Overall),
    y = Overall,
    fill = Genes
  )) +
    geom_col(width = 0.7) +
    scale_y_continuous(expand = c(0, 0))+
    labs(
      x = "Gene",
      y = "Feature Importance",
      title = title
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        size = font_size - 2,
        angle = 45,
        hjust = 1,
        vjust = 0.5,
        margin = margin(t = -12)
      ),
      axis.text.y = element_text(size = font_size + 2),
      axis.title  = element_text(size = font_size),
      plot.title  = element_text(size = font_size + 2, hjust = 0.5)
    )
}
plot_metric_bar <- function(data, dataset_type, y_label,
                            gene_subset = c("LowOutput", "DVG_33", "DVG_80", "DVG_union"),
                            base_size = 14, show_legend = TRUE) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  df <- data %>%
    filter(dataset == dataset_type,
           gene_list %in% gene_subset) %>%
    mutate(
      gene_list = recode(gene_list,
                         LowOutput = "DEG",
                         DVG_33 = "DVG33",
                         DVG_80 = "DVG80",
                         DVG_union = "DVGunion")
    ) %>%
  
    pivot_longer(
      cols = c(Accuracy),#macro_F1, 
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric,
                      macro_F1 = "F1 Score",
                      Accuracy = "Accuracy")
    )
  
  p <- ggplot(df, aes(
    x = condition,
    y = value,
    fill = metric   
  )) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.75),
      width = 0.65
    ) +
 
    scale_fill_manual(values = c(
      "F1 Score" = "#1f77b4",
      "Accuracy" = "#e74c3c"
    ))+
   scale_y_continuous(expand = c(0, 0))+
    labs(
      x = NULL,
      y = y_label,
      fill = NULL
    ) +
    theme_classic(base_size = base_size)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}
############################
# 5. Generate plots
############################

# F1 plots
p_vitro_main <- plot_f1_bar(
  ml_res,
  dataset_type = "vitro",
  y_label = "Macro F1 (Early Prediction)",
  gene_subset = c("DVG_33"),
  base_size = 18,
  show_legend = TRUE
)
p_vitro_full <- plot_f1_bar(
  ml_res,
  dataset_type = "vitro",
  y_label = "Macro F1 (Early Prediction)"
)

p_vivo_main <- plot_f1_bar(
  ml_res,
  dataset_type = "vivo",
  y_label = "Macro F1 (Late Prediction)",
  gene_subset = c("DVG_33"),
  base_size = 18,
  show_legend = TRUE
)
p_vivo_full <- plot_f1_bar(
  ml_res,
  dataset_type = "vivo",
  y_label = "Macro F1 (Late Prediction)"
)
# F1 + accuracy
p_vitro_early <- plot_metric_bar(
  ml_res,
  dataset_type = "vitro",
  y_label = "Early Prediction",
  gene_subset = c("DVG_33"),
  base_size = 18,
  show_legend = TRUE
)
p_vitro_late <- plot_metric_bar(
  ml_res,
  dataset_type = "vivo",
  y_label = "Late Prediction",
  gene_subset = c("DVG_33"),
  base_size = 18,
  show_legend = TRUE
)
# Feature importance plots
DVG_33_vitro_OO <- plot_top20_as_bar(
  DVG_33_O_vitro_OO[, 1:2],
  title = "Early Prediction(O >> OO)"
)
DVG_33_vivo_OO <- plot_top20_as_bar(
  DVG_33_OO[, 1:2],
  title = "Late Prediction(O >> OO)"
)
DVG_33_vitro_YY <- plot_top20_as_bar(
  DVG_33_Y_vitro_YY[, 1:2],
  title = "Early Prediction(Y >> YY)"
)
DVG_33_vivo_YY <- plot_top20_as_bar(
  DVG_33_YY[, 1:2],
  title = "Late Prediction(Y >> YY)"
)
DVG_80_vitro_OO <- plot_top20_as_bar(
  DVG_80_O_vitro_OO[, 1:2],
  title = "DVG Gene 80 O-vitro OO Self-renewal"
)
DVG_80_vivo_OO <- plot_top20_as_bar(
  DVG_80_OO[, 1:2],
  title = "DVG Gene 80 OO Self-renewal"
)
DVG_union_vitro_OO <- plot_top20_as_bar(
  DVG_union_O_vitro_OO[, 1:2],
  title = "DVG O-vitro OO Self-renewal"
)
DVG_union_vivo_OO <- plot_top20_as_bar(
  DVG_union_OO[, 1:2],
  title = "DVG OO Self-renewal"
)

DEG_vitro_OO <- plot_top20_as_bar(
  DEG_O_vitro_OO[, 1:2],
  title = "DEG O-vitro OO Self-renewal"
)
DEG_vivo_OO <- plot_top20_as_bar(
  DEG_OO[, 1:2],
  title = "DEG OO Self-renewal"
)
############################
# 6. Print plots
############################

