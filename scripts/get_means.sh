#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=get_means.log
#SBATCH --error=get_means.err
module load intel intel-mkl r
Rscript get_means.R
