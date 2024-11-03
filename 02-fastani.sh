#!/bin/sh

#SBATCH --job-name="fastANI"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e fastANI.err-%N
#SBATCH -o fastANI.out-%N

# load module
module load dRep

OUTPUT=/home/ad14556/pangenome-large/results/02-drep
INPUT=/home/ad14556/pangenome-large/data/genomes

#dRep dereplicate "$OUTPUT" -g "$INPUT"/* --S_algorithm fastANI --S_ani 0.999 --SkipMash  
dRep compare "$OUTPUT" -g "$INPUT"/* --S_algorithm fastANI --S_ani 0.95 --SkipMash

