#!/usr/bin/env bash
#SBATCH --mem=100G
#SBATCH --output=cut_height.log
#SBATCH --error=cut_height.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=6
#SBATCH --ntasks=1
set -euo pipefail

module load gcc/7.4.0 openblas/0.3.6-openmp intel-tbb/2019.4 mvapich2/2.3.1 r boost/1.67.0-mpi

Rscript ../R/height_cut_S.R --beta_file ../data/betas/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.rds  --sigma_file ../data/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.csvLong --chain 1 --out ../output/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.rds
