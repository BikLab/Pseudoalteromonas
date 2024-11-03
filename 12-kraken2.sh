#!/bin/sh

#SBATCH --job-name="kraken2"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e kraken2.err-%N
#SBATCH -o kraken2.out-%N

#Path Variables
module load Kraken2

INPUT=/scratch/ad14556/nematode-microbiome-final/raw-data-jgi
OUTPUT=/scratch/ad14556/pangenome-large/results/kraken
DB=

#kraken2-build --download-library bacteria --db "$OUTPUT"/bacteria --threads 24

#for file in "$INPUT"/*_R1.fastq.gz; do
#  base=$(basename $file _R1.fastq.gz)
#  kraken2 --db /db/kraken2/20240906/pluspf \
#  --paired --classified-out "$OUTPUT"/samples/"$base"_R#.fq "$INPUT"/"$base"_R1.fastq.gz "$INPUT"/"$base"_R2.fastq.gz --gzip-compressed --threads 24
#done

pigz --processes 24 /scratch/ad14556/pangenome-large/results/kraken/samples/*
