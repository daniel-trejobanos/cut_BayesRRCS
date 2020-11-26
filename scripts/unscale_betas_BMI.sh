#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=/work/ext-unil-ctgg/tdaniel/genomicarchitecture/logs/unscale_betas_BMI_${1}.log
#SBATCH --error=/work/ext-unil-ctgg/tdaniel/genomicarchitecture/logs/unscale_betas_BMI_${1}.err
module load intel intel-mkl r

BETMAPDIR=/work/ext-unil-ctgg/marion/annot0820/post_processing/ukb_BMI_78groups
BETMAPFILE=groups78_mix4_cpus4_tasks8_nodes10_mpisync10_BMI_$1_unrelated_148covadjusted_w35520NA.betMap

Rscript ../R/unscale_betas.R\
    --chain ${BETMAPDIR}/${BETMAPFILE}\
    --start 500 --end 1200\
    --mafsd ../data/maf_sd.rds\
    --out ../data/betas/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_BMI_$1_unrelated_148covadjusted_w35520NA.rds
