#!/bin/sh

#SBATCH --job-name="rrap"
#SBATCH --partition=iob_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e rrap.err-%N
#SBATCH -o rrap.out-%N

#Path Variables
module load Miniconda3
source activate /home/ad14556/conda-envs/envs/rrap

INPUT=/scratch/ad14556/nematode-microbiome-final
DATA=/home/ad14556/pseudoalteromonas-pangenome/data/genomes
OUTPUT=/scratch/ad14556/pangenome-large/results/rrap-bac

rrap -i "$INPUT"/path-metagenomes.txt -rg "$DATA" \
	-o "$OUTPUT" -n nematode_metagenomes -suffix "_R1.fastq.gz" --threads 25
