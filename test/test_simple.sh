#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --time=00:20:00
#SBATCH --output=sim_height.log
#SBATCH --error=sim_height.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
module load gcc/7.4.0 openblas/0.3.6-openmp intel-tbb/2019.4 mvapich2/2.3.1 r boost/1.67.0-mpi
source /home/tdaniel/tensorflow/activate
Rscript ../R/height.R
