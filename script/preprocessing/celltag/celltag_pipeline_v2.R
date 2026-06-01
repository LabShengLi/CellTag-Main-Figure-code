#!/usr/bin/env Rscript
# CellTagR clone calling (multi-BAM folder) + Barcode.Aggregate
# Save the object after each major step with meaningful filenames
library(argparser)
library(CellTagR)
library(Matrix)

# -----------------------------
# Parse args
# -----------------------------
p <- arg_parser("CellTagR clone calling (multi-BAM folder + Barcode.Aggregate + save each step)")

p <- add_argument(p, "--dsname", help = "Dataset/sample name for outputs", default = NULL)
p <- add_argument(p, "--wdir", help = "Working directory", default = NULL)
p <- add_argument(p, "--bam_dir", help = "Folder containing ONLY BAM files (passed to fastq.bam.directory)", default
    = NULL)
p <- add_argument(p, "--barcode_fns", help = "Comma-separated barcodes files in SAME ORDER as BAM order", default =
    NULL)
p <- add_argument(p, "--barcode_agg_fn", help = "Aggregated barcode output file (e.g., ./barcodes_all.tsv)", default
    = "barcodes_all.tsv")
p <- add_argument(p, "--celltag_version", help = "celltag version", default = "v1")
# A 90% percentile cut-off in terms of reads reported for each CellTag was used to select CellTags for inclusion on the whitelist of cell barcodes. ref: https://www-nature-com.libproxy2.usc.edu/articles/s41586-018-0744-4#Sec7
p <- add_argument(p, "--cutoff_Whitelist_fn", help = "Whitelist cutoff used in filename", type = "str", default =
    "whitelist/v1_whitelist_Cutoff_0.9.csv")
# The suggested cutoff that marks presence or absence is at least 2 counts per CellTag per Cell. For details regarding cutoff choice, please refer to the paper - https://www.nature.com/articles/s41586-018-0744-4. ref: https://github.com/morris-lab/CellTagR?tab=readme-ov-file#2-binarize-the-single-cell-celltag-umi-count-matrix
p <- add_argument(p, "--binarize_k", help = "k used in SingleCellDataBinarization()", type = "integer", default = 2)
# Cells expressing more than 20 CellTags (likely to correspond to cell multiplets), and less than 2 CellTags per cell were filtered out.  ref: https://www-nature-com.libproxy2.usc.edu/articles/s41586-018-0744-4#Sec7
p <- add_argument(p, "--metric_less", help = "MetricBasedFiltering threshold with comparison='less'", type =
    "numeric", default = 20)
p <- add_argument(p, "--metric_greater", help = "MetricBasedFiltering threshold with comparison='greater'", type =
    "numeric", default = 1)
p <- add_argument(p, "--correlation_cutoff", help = "CloneCalling correlation.cutoff", type = "numeric", default = 0.7)
p <- add_argument(p, "--skip_raw_read", help = "Skip read raw input", flag = TRUE)
p <- add_argument(p, "--count_rds", help = "Skip read raw input", type = "str", default = "CellTagMatrixCount.RDS")
p <- add_argument(p, "--collapsing", help = "Use collapsing step", flag = TRUE)
p <- add_argument(p, "--is_debug", help = "Is debug or not", flag = TRUE)
p <- add_argument(p, "--stop_at_count_matrix", help = "If stop at count matrix", flag = TRUE)
p <- add_argument(p, "--plot_heatmap", help = "If plot heatmap or not", flag = TRUE)
p <- add_argument(p, "--outdir", help = "Output directory", default = '.')
argv <- parse_args(p)

argv

# -----------------------------
dsname <- argv$dsname
wdir <- argv$wdir
bam_dir <- argv$bam_dir
barcode_fns <- argv$barcode_fns
barcode_agg_fn <- argv$barcode_agg_fn
celltag_version <- argv$celltag_version
cutoff_Whitelist_fn <- argv$cutoff_Whitelist_fn
binarize_k <- argv$binarize_k
metric_less <- argv$metric_less
metric_greater <- argv$metric_greater
correlation_cutoff <- argv$correlation_cutoff
outdir <- argv$outdir
skip_raw_read <- argv$skip_raw_read
collapsing <- argv$collapsing
count_rds <- argv$count_rds
is_debug <- argv$is_debug
plot_heatmap <- argv$plot_heatmap
stop_at_count_matrix <- argv$stop_at_count_matrix

if (is_debug) {
    dsname = "exp2_y1a_y1b"
    wdir <- "projects/U01_aim2/results/2026_01_27_celltag_collapsing_test"
    bam_dir <- "bam_in"
    barcode_fns <- "barcode/barcode_1.tsv.gz,barcode/barcode_2.tsv.gz"
    skip_raw_read <- T
}

# -----------------------------
# set wdir and outdir
# -----------------------------
dir.create(wdir)
setwd(wdir)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Utilities
# -----------------------------
bam_dir_norm <- bam_dir
if (!grepl("/$", bam_dir_norm)) bam_dir_norm <- paste0(bam_dir_norm, "/")

barcode_list <- trimws(unlist(strsplit(barcode_fns, ",")))

dt.mtx.whitelist.path <- cutoff_Whitelist_fn
print(dt.mtx.whitelist.path)

prefix <- paste0(outdir, '/')

save_obj <- function(obj, tag) {
    saveRDS(obj, paste0(prefix, tag, ".RDS"))
    invisible(obj)
}

if (!skip_raw_read) {
    # -----------------------------
    # 0) Aggregate barcodes
    # -----------------------------

    if (length(barcode_list) > 1) {
        final_barcode_agg_fn <- file.path(outdir, barcode_agg_fn)
        Barcode.Aggregate(barcode_list, final_barcode_agg_fn)

    } else {
        final_barcode_agg_fn <- barcode_list
    }
    print(final_barcode_agg_fn)

    saveRDS(
        list(barcode_inputs = barcode_list, barcode_agg_fn = final_barcode_agg_fn),
        paste0(prefix, "BarcodeAggregate.meta.RDS")
    )

    # -----------------------------
    # 1) CellTagR pipeline + save after each step
    # -----------------------------
    ct.obj <- CellTagObject(
        object.name = paste0(dsname, ".celltag.obj"),
        fastq.bam.directory = bam_dir_norm
    )
    save_obj(ct.obj, "CellTagObject")

    # Take time.......
    # Update: CellTagR now enables read-in of multiple BAM files at a time. When multiple BAM files need to be
    # processed, please use a folder that contains ONLY BAM files and put the fastq.bam.directory as the path of the
    # folder. For instance, two bam files need to be processed named as bam1.bam and bam2.bam. They will be put into a
    # folder named as beautiful_bams in the Desktop. Then, the input will be fastq.bam
    # .directory="~/Desktop/beautiful_bams/" as below.
    # name as bam1.bam, bam2.bam, etc.
    ct.obj <- CellTagExtraction(ct.obj, celltag.version = celltag_version)
    save_obj(ct.obj, "CellTagExtraction")

    ct.obj <- CellTagMatrixCount(
        celltag.obj = ct.obj,
        barcodes.file = final_barcode_agg_fn,
        replace.option = TRUE
    )
    save_obj(ct.obj, "CellTagMatrixCount")

    # save to tsv gz
    mat <- ct.obj@raw.count
    # may be large
    mat_dense <- as.matrix(mat)

    out <- sprintf("%s/celltag.raw.count.tsv.gz", outdir)  # Recommended to gzip directly
    write.table(
        mat_dense,
        file = gzfile(out),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
    )

    message(sprintf("raw.count dim: %s x %s", nrow(ct.obj@raw.count), ncol(ct.obj@raw.count)))
} else {
    ct.obj <- readRDS(count_rds)
}

if (stop_at_count_matrix) {
    quit()
}

if (collapsing) {
    # Generating the collapsing file
    collapsing_fn <- sprintf("%s/collapsing.txt", outdir)
    ct.obj <- CellTagDataForCollapsing(celltag.obj = ct.obj,
                                       output.file = collapsing_fn)
    save_obj(ct.obj, "CellTagDataForCollapsing")

    collapsed.rslt.dir <- sprintf("%s/star_collapse", outdir)
    dir.create(collapsed.rslt.dir)

    # starcode -s --print-clusters collapsing_Sample-1.txt > star_collapse/collapsing_result_Sample-1.txt
    # starcode -s --print-clusters collapsing_Sample-2.txt > star_collapse/collapsing_result_Sample-2.txt

    # Generating the collapsing file
    # Recount and generate collapsed matrix

    if (length(barcode_list) > 1) {
        for (k in 1:length(barcode_list)) {
            infn1 <- sprintf("%s/collapsing_Sample-%d.txt", outdir, k)
            outfn1 <- sprintf("%s/collapsing_result_Sample-%d.txt", collapsed.rslt.dir, k)
            print(infn1)
            print(outfn1)
            res <- system2(
                command = "starcode",
                args = c("-s", "--print-clusters", infn1),
                stdout = outfn1,
                stderr = paste0(outfn1, ".err")
            )
            # print(res)
        }
        ct.obj <- CellTagDataPostCollapsing(celltag.obj = ct.obj,
                                            collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
    } else {
        infn1 <- sprintf("%s/collapsing.txt", outdir)
        outfn1 <- sprintf("%s/collapsing_rslt.txt", collapsed.rslt.dir)
        print(infn1)
        print(outfn1)
        res <- system2(
            command = "starcode",
            args = c("-s", "--print-clusters", infn1),
            stdout = outfn1,
            stderr = paste0(outfn1, ".err")
        )
        ct.obj <- CellTagDataPostCollapsing(celltag.obj = ct.obj,
                                            collapsed.rslt.file = outfn1)
    }
    save_obj(ct.obj, "CellTagDataPostCollapsing")

    # save matrix
    # save to tsv gz
    mat <- ct.obj@collapsed.count
    # may be large
    mat_dense <- as.matrix(mat)

    out <- sprintf("%s/celltag.collapsed.count.tsv.gz", outdir)  # Recommended to gzip directly
    write.table(
        mat_dense,
        file = gzfile(out),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
    )
}

ct.obj <- SingleCellDataBinarization(ct.obj, binarize_k, replace.option = T)
save_obj(ct.obj, "SingleCellDataBinarization")
dim(ct.obj@binary.mtx)

out_png <- sprintf("%s/metric_plots_before_all_filtering_%s.png", outdir, dsname)
png(out_png, width = 14, height = 10, units = "in", res = 300)
MetricPlots(ct.obj)
dev.off()

ct.obj <- SingleCellDataWhitelist(ct.obj, dt.mtx.whitelist.path, replace.option = T)
save_obj(ct.obj, "SingleCellDataWhitelist")

# Workaround for some CellTagR versions
# ct.obj@metric.filtered.count <- as(matrix(NA, 0, 0), "dgCMatrix")

ct.obj <- MetricBasedFiltering(ct.obj, metric_less, comparison = "less", replace.option = T)
save_obj(ct.obj, "MetricBasedFiltering_less")

ct.obj <- MetricBasedFiltering(ct.obj, metric_greater, comparison = "greater", replace.option = T)
save_obj(ct.obj, "MetricBasedFiltering_greater")

out_png <- sprintf("%s/metric_plots_after_all_filtering_%s.png", outdir, dsname)
png(out_png, width = 14, height = 10, units = "in", res = 300)
MetricPlots(ct.obj)
dev.off()

ct.obj <- JaccardAnalysis(ct.obj, plot.corr = F, fast = TRUE)
save_obj(ct.obj, "JaccardAnalysis_fast")
dim(ct.obj@jaccard.mtx)

heatmap_jaccard <- function(ct.obj, out_png) {
    library(corrplot)
    Jac <- ct.obj@jaccard.mtx
    Jac <- as.matrix(Jac)
    diag(Jac) <- 1
    png(out_png, width = 10, height = 10, units = "in", res = 300)
    corrplot(
        Jac,
        is.corr = T,
        method = "color",
        order = "hclust",
        hclust.method = "ward.D2",
        col.lim = c(0, 1),
        tl.cex = 0.5
    )
    dev.off()
}

if (plot_heatmap) {
    heatmap_jaccard(ct.obj, sprintf("%s/heatmap_jaccard_%s.png", outdir, dsname))
}

ct.obj <- CloneCalling(ct.obj, correlation.cutoff = correlation_cutoff)
save_obj(ct.obj, "CloneCalling")

# -----------------------------
# Final tables
# -----------------------------
write.csv(ct.obj@clone.composition[[celltag_version]],
          paste0(prefix, "clone_composition.csv"),
          row.names = FALSE)
write.csv(ct.obj@clone.size.info[[celltag_version]],
          paste0(prefix, "clone_size_info.csv"),
          row.names = FALSE)

save.image(file = sprintf("%s/Run_workspace_celltag_v2_workflow_%s.RData", outdir, dsname))
message("Done.")
