#!/bin/bash
#SBATCH --job-name=cr_multi
#SBATCH -p main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16         # <— bump this to what your node can give
#SBATCH --mem=60G                 # <— and give it plenty of RAM
#SBATCH --time=2-00:00:00
#SBATCH --output=log/%x.%j.log

# # cellranger preprocessing of raw single cell with hashtag data

set -euo pipefail

module purge
module load cellranger/9.0.1

# --------------------- user inputs ---------------------
ID=${1:? "Usage: sbatch this.sh <RunID>"}             # e.g., CrossExp_YO_YY
CSV=multi_config_${ID}_ly.csv
REF=refdata-gex-mm10-2020-A
#OUTBASE=Out
OUTBASE=2026_02_26_cellranger_multi
# ------------------------------------------------------

CPUS=${SLURM_CPUS_PER_TASK:-16}
LOCALMEM_GB=60                      # keep a buffer below --mem

mkdir -p "${OUTBASE}"
cd "${OUTBASE}"

# clean any previous attempt
rm -rf "${ID}"

cellranger multi \
  --id="${ID}" \
  --csv="${CSV}" \
  --localcores="${CPUS}" \
  --localmem="${LOCALMEM_GB}"

echo "### Done"
