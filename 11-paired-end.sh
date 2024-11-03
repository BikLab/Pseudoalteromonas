#!/bin/sh

#SBATCH --job-name="PE"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e PE.err-%N
#SBATCH -o PE.out-%N

#Path Variables
module load BBMap

INPUT=/scratch/ad14556/nematode-microbiome-final/raw-data-jgi

for file in "$INPUT"/*; do
  base=$(basename $file .fastq.gz)
  reformat.sh in="$file" out1="$INPUT"/"$base"_R1.fastq.gz out2="$INPUT"/"$base"_R2.fastq.gz
done

