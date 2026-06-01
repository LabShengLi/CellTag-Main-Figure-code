## ==========================
## 0. LOAD LIBRARIES
## ==========================
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(edgeR)
  library(dplyr)
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(survival)
  library(survminer)
  library(readr)
  library(tibble)
})

## ==========================
## 1. PREPARE TCGA-LAML DATA
## ==========================
prepare_tcga_laml <- function() {
  
  query <- GDCquery(
    project = "TCGA-LAML",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query)
  aml_raw <- GDCprepare(query)
  
  counts <- assay(aml_raw, "unstranded")
  logCPM <- cpm(counts, log = TRUE)
  
  # ---- Gene annotation ----
  ensembl_ids <- gsub("\\..*", "", rownames(logCPM))
  
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  expr <- logCPM[!is.na(symbols), ]
  rownames(expr) <- toupper(symbols[!is.na(symbols)])
  
  # ---- Remove duplicated genes ----
  expr <- expr %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(avg = rowMeans(dplyr::across(where(is.numeric)))) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(avg, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-avg) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  # ---- Clinical ----
  clinical <- GDCquery_clinic("TCGA-LAML", "clinical") %>%
    mutate(
      sample = submitter_id,
      OS_time = ifelse(is.na(days_to_death),
                       days_to_last_follow_up,
                       days_to_death),
      OS_status = ifelse(vital_status == "Dead", 1, 0),
      age = as.numeric(age_at_index)
    )
  
  # ---- SAFE MATCHING----
  samples <- substr(colnames(expr), 1, 12)
  common <- intersect(samples, clinical$sample)
  
  expr <- expr[, samples %in% common]
  samples <- samples[samples %in% common]
  
  clinical <- clinical[match(samples, clinical$sample), ]
  
  stopifnot(all(samples == clinical$sample))  
  
  # ---- Age group ----
  clinical <- clinical %>% filter(!is.na(age))
  
  clinical$age_group <- ifelse(
    clinical$age > median(clinical$age, na.rm = TRUE),
    "Old", "Young"
  )
  
  list(expr = expr, clinical = clinical)
}

tcga <- prepare_tcga_laml()
expr <- tcga$expr
clinical <- tcga$clinical


## ==========================
## 2. SIGNATURE FUNCTION
## ==========================
compute_signature <- function(expr, genes_up, genes_down = NULL) {
  
  zscore <- function(x) {
    if (nrow(x) == 0) return(NULL)
    sd_vec <- apply(x, 1, sd)
    sd_vec[sd_vec == 0] <- 1
    (x - rowMeans(x)) / sd_vec
  }
  
  genes_up <- intersect(toupper(genes_up), rownames(expr))
  genes_up_input <- toupper(genes_up)
  genes_up <- intersect(genes_up_input, rownames(expr))
  
  cat("Genes used:", length(genes_up), "/", length(genes_up_input), "\n")
  z_up <- zscore(expr[genes_up, , drop = FALSE])
  score_up <- colMeans(z_up, na.rm = TRUE)
  
  if (is.null(genes_down)) return(scale(score_up)[,1])
  
  genes_down <- intersect(toupper(genes_down), rownames(expr))
  z_down <- zscore(expr[genes_down, , drop = FALSE])
  score_down <- colMeans(z_down, na.rm = TRUE)
  
  scale(score_up - score_down)[,1]
}


## ==========================
## 3. SURVIVAL FUNCTION
## ==========================
run_survival <- function(clinical, score, title, split_by = NULL) {
  
  df <- clinical
  df$score <- as.numeric(score)
  
  # subgroup
  if (!is.null(split_by)) {
    df <- df[df[[split_by$var]] == split_by$value, ]
  }
  
  # clean NA
  df <- df[complete.cases(df[, c("OS_time", "OS_status", "score")]), ]
  
  df$group <- ifelse(
    df$score > median(df$score, na.rm = TRUE),
    "High", "Low"
  )
  
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = df)
  
  plot <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = FALSE,
    ggtheme = theme_classic(),
    title = title
  )
  
  plot$plot <- plot$plot +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )
  
  cox <- summary(
    coxph(Surv(OS_time, OS_status) ~ score, data = df)
  )
  
  list(plot = plot, cox = cox)
}


## ==========================
## 4. TOP20 SIGNATURE
## ==========================
# ==========================
# MULTI-SIGNATURE (Top10/20/30)
# ==========================

top_tbl <- read_csv("data/DVG_union_OO_Top50_Importance.csv") %>%
  arrange(desc(Overall))
length(top_tbl$Gene)

top_sizes <- c(5, 10, 15, 20, 25, 30, 35)

results <- list()

for (n in top_sizes) {
  
  genes <- toupper(top_tbl$Gene[1:n])
  
  score <- compute_signature(expr, genes)
  
  res <- run_survival(
    clinical,
    score,
    paste0("Top", n, " Signature")
  )
  
  results[[paste0("Top", n)]] <- res
  
  print(res$plot) 
  print(res$cox)    
}
results$Top5$plot
results$Top10$plot
results$Top15$plot
results$Top20$plot
results$Top25$plot
results$Top30$plot
results$Top35$plot
results$Top5$cox
results$Top10$cox
results$Top15$cox
results$Top20$cox
results$Top25$cox
results$Top30$cox
results$Top35$cox
## ==========================
## 5. AGE STRATIFIED
## ==========================
res_old <- run_survival(
  clinical, score_top,
  "Top20 Signature - Old",
  split_by = list(var = "age_group", value = "Old")
)

res_young <- run_survival(
  clinical, score_top,
  "Top20 Signature - Young",
  split_by = list(var = "age_group", value = "Young")
)

res_old$plot
res_young$plot


## ==========================
## 6. SINGLE GENE
## ==========================
run_single_gene <- function(expr, clinical, gene) {
  
  gene <- toupper(gene)
  
  if (!gene %in% rownames(expr)) {
    stop(paste("Gene not found:", gene))
  }
  
  score <- as.numeric(expr[gene, ])
  
  run_survival(clinical, score, paste("AML OS ~", gene))
}
top20_genes <- toupper(top_tbl$Gene[1:20])
res_gene <- run_single_gene(expr, clinical, "HLF")
res_gene$plot
res_gene$cox


## ==========================
## 7. UP / DOWN SIGNATURE
## ==========================
out_dir <- "output/"

genes_up <- read.csv(
  file.path(out_dir, "LowOutput_OO_Top10_Up_in_High.csv")
)$Gene

genes_down <- read.csv(
  file.path(out_dir, "LowOutput_OO_Top10_Up_in_Low.csv")
)$Gene

score_up <- compute_signature(expr, genes_up)
score_down <- compute_signature(expr, genes_down)
score_combined <- compute_signature(expr, genes_up, genes_down)

res_up <- run_survival(clinical, score_up, "Up Signature")
res_down <- run_survival(clinical, score_down, "Down Signature")
res_combined <- run_survival(clinical, score_combined, "Up-Down Signature")

res_up$plot
res_down$plot
res_combined$plot

res_up$cox
res_down$cox

res_combined$cox
