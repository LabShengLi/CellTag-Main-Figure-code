# CellTag-Main-Figure-code

## Overview

This repository contains the code used to reproduce Figures 1–5 from the CellTag manuscript.

## Repository structure

Figure1/
Code used to generate Figure 1.

Figure2/
Code used to generate Figure 2.

Figure3/
Code used to generate Figure 3.

Figure4/
Code used to generate Figure 4.

Figure5/
Code used to generate Figure 5.

script/preprocessing/
Preprocessing scripts for scRNA-seq, CellTag, and hashtag analyses.

script/analysis/
Downstream analysis scripts.

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.

## Installation

No installation is required. Users only need access to an R environment with the required packages installed (see Environment below).

Typical setup time: less than 10 minutes.

## Script-to-figure mapping

| Script | Produces |
|---|---|
| `Figure1/Figure 1c and 1d.R` | Fig 1c, 1d |
| `Figure1/figure 1h.R` | Fig 1h |
| `Figure1/Figure1_befg_figure2_cd_OA_function.R` | Fig 1b, 1e, 1f, 1g; Fig 2c, 2d (shared function) |
| `Figure2/Figure2_ef_clonal_size_distributions.R` | Fig 2e, 2f |
| `Figure3/Figure 3a.R` | Fig 3a |
| `Figure3/Figure3_d_Day0_PCA_inter_intra_dist.R` | Fig 3d |
| `Figure3/Figure3_bcegh_DVGs_3_data_function.R` | Fig 3b, 3c, 3e, 3g, 3h |
| `Figure3/Figure3_fi_GSEA_Day0_Unmani.R` | Fig 3f, 3i |
| `Figure3/Figure 3k.R` | Fig 3k |
| `Figure3/Figure3_lmno_human_ORA_PROGENy.R` | Fig 3l, 3m, 3n, 3o |
| `Figure4/Panel A.R` | Fig 4a |
| `Figure4/Panel B.R` | Fig 4b |
| `Figure4/Panel C.R` | Fig 4c |
| `Figure4/Figure4_f_edfig5_c.R` | Fig 4f; Extended Data Fig 5c |
| `Figure4/Figure4_g.R` | Fig 4g |
| `Figure5/Figure 5a.R` | Fig 5a |
| `Figure5/Figure 5c-5h.R` | Fig 5c–5h |
| `Figure5/Figure 5i.R` | Fig 5i |
| `Figure5/Fig5_j_DVG_MAGMA.R` | Fig 5j |

Expected runtime per script: approximately 5–30 minutes depending on the computational environment and dataset size.

## Environment

Analyses were run under **R 4.4.1** on macOS (x86_64-apple-darwin20).

Key packages:

| Package | Version |
|---|---|
| ggplot2 | 4.0.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| data.table | 1.17.0 |
| clusterProfiler | 4.12.6 |
| enrichplot | 1.24.4 |
| msigdbr | 25.1.1 |
| biomaRt | 2.60.1 |
| org.Hs.eg.db | 3.19.1 |
| susieR | 0.14.2 |
| pheatmap | 1.0.12 |
| patchwork | 1.3.0 |

Full session details available via `sessionInfo()`.

## Data and code availability

Analysis code is archived on Zenodo: [10.5281/zenodo.20709794](https://doi.org/10.5281/zenodo.20709794), corresponding to release `v1.0-natsubmission` of this repository.

Processed data are deposited on Zenodo: [10.5281/zenodo.20480987](https://doi.org/10.5281/zenodo.20480987).
