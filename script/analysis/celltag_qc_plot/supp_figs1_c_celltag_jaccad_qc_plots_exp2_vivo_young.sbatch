#!/bin/bash
#SBATCH --job-name=celltag_qc_jac
#SBATCH -p largemem
#SBATCH -n 4
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=/project2/sli68423_1316/users/yang/log/%x.%j.log
#SBATCH --mail-user=yliu8962@usc.edu
#SBATCH --mail-type=END

module purge
module load gcc/13.3.0 openblas/0.3.28 r/4.4.3

date
hostname
pwd
set -e

R --vanilla -f \
    supp_figs1_c_celltag_jaccad_qc_plots_exp2_vivo_young.R

echo "### Done"

