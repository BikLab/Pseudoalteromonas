#!/bin/bash

#SBATCH --job-name="checkm"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e checkm.err-%N
#SBATCH -o checkm.out-%N

module load CheckM-Database
module load CheckM

INPUT=/home/ad14556/pangenome-large/data/genomes
OUTPUT=/home/ad14556/pangenome-large/results/checkm

#checkm lineage_wf -t 8 -x fna "$INPUT" "$OUTPUT"
checkm qa -o 2 --tab_table -f "$OUTPUT"/checkm-statistics.txt "$OUTPUT"/lineage.ms "$OUTPUT"
