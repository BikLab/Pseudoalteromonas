#!/bin/sh

#SBATCH --job-name="metabolic"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e metabolic.err-%N
#SBATCH -o metabolic.out-%N

#Path Variables
module load HMMER
module load METABOLIC
module load GTDB-Tk

INPUT=/scratch/ad14556/pangenome-large/results/metabolic/input
OUTPUT=/scratch/ad14556/pangenome-large/results/metabolic/output

#cp /scratch/ad14556/pangenome-large/results/prokka/*.fna "$INPUT"

METABOLIC-G.pl -p meta -t 25 -in-gn "$INPUT" -o "$OUTPUT"
