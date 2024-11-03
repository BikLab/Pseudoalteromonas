#!/bin/sh

#SBATCH --job-name="antismash"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e antismash.err-%N
#SBATCH -o antismash.out-%N

#Path Variables
OUTPUT=/scratch/ad14556/pangenome-large/results/antismash

mkdir -p "$OUTPUT"

for file in /home/ad14556/pseudoalteromonas-pangenome/data/genomes/*; do \
 sample=$(basename "$file" .fna)
 singularity exec /apps/singularity-images/antismash_7.0.0--c222958.sif antismash \
 --output-basename "${sample}" --output-dir "$OUTPUT"/"${sample}" \
 --genefinding-tool prodigal "${file}"
done

# --fullhmmer --cassis --clusterhmmer --tigrfam --asf --cc-mibig --pfam2go --smcog-tree ${file} \
