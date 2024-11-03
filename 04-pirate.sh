#!/bin/sh

#SBATCH --job-name="pirate"
#SBATCH --partition=iob_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e pirate.err-%N
#SBATCH -o pirate.out-%N

module load Miniconda3
module load R-bundle-Bioconductor
source activate /home/ad14556/conda-envs/envs/pirate

PIRATE -i /scratch/ad14556/pangenome-large/results/prokka/ \
       -o /scratch/ad14556/pangenome-large/results/pirate/ \
       -a -r -t 8
