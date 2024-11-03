#!/bin/sh

#SBATCH --job-name="grapple"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e grapple.err-%N
#SBATCH -o grapple.out-%N

# load module
module load Miniconda3
module load pandas

python /home/ad14556/pseudoalteromonas-pangenome/scripts/07-grapple/grapple/pw_similarity.py \
  -i /scratch/ad14556/pangenome-large/results/pirate/PIRATE_90_01.tsv -o /scratch/ad14556/pangenome-large/results/grapple-90-01/ -r "both" -s "jaccard"
