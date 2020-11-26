#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=/work/ext-unil-ctgg/tdaniel/genomicarchitecture/logs/unscale_betas_${1}.log
#SBATCH --error=/work/ext-unil-ctgg/tdaniel/genomicarchitecture/logs/unscale_betas_${1}.err
module load intel intel-mkl r

Rscript ../R/unscale_betas.R --chain ../data/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_$1_unrelated_148covadjusted_w35520NA.betMap --start 500 --end 1200 --mafsd ../data/maf_sd.rds --out ../data/betas/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_$1_unrelated_148covadjusted_w35520NA.rds
