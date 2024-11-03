#!/bin/sh

#SBATCH --job-name="bigscape"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e bigscape.err-%N
#SBATCH -o bigscape.out-%N

#Path Variables
module load Miniconda3

INPUT=/home/ad14556/pseudoalteromonas-pangenome/scripts/07-multismash.yaml

multismash "$INPUT"
