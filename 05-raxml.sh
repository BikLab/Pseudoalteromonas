#!/bin/sh

#SBATCH --job-name="raxml"
#SBATCH --partition=iob_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=10G
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e raxml.err-%N
#SBATCH -o raxml.out-%N

#Path Variables
#module load RAxML/8.2.12-foss-2022a-pthreads-avx
module load RAxML/8.2.12-gompi-2022a-hybrid-avx

INPUT=/scratch/ad14556/pangenome-large/results/pirate/core_alignment.fasta
OUTPUT=/scratch/ad14556/pangenome-large/results/raxml
MOLLUSCA=/scratch/ad14556/pangenome-large/results/pirate/mollusca_core_alignment.fasta
BIVALVIA=/scratch/ad14556/pangenome-large/results/pirate/mollusca_bivalvia_core_alignment.fasta
PORIFERA=/scratch/ad14556/pangenome-large/results/pirate/porifera_core_alignment.fasta
CNIDARIA=/scratch/ad14556/pangenome-large/results/pirate/cnidaria_core_alignment.fasta
SCLERACTINIA=/scratch/ad14556/pangenome-large/results/pirate/cnidaria_scleractinia_core_alignment.fasta
FAVIIDAE=/scratch/ad14556/pangenome-large/results/pirate/cnidaria_faviidae_core_alignment.fasta

mkdir -p "$OUTPUT"

#srun raxmlHPC -f a -T 20 -s "$INPUT" -n pseudoalteromonas-pangenome -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
#srun raxmlHPC -f a -T 25 -s "$MOLLUSCA" -n pseudoalteromonas-mollusca -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
#srun raxmlHPC -f a -T 25 -s "$PORIFERA" -n pseudoalteromonas-porifera -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
#srun raxmlHPC -f a -T 20 -s "$BIVALVIA" -n pseudoalteromonas-bivalvia -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log

#srun raxmlHPC -f a -T 20 -s "$CNIDARIA" -n pseudoalteromonas-cnidaria -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
srun raxmlHPC -f a -T 20 -s "$FAVIIDAE" -n pseudoalteromonas-faviidae -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
srun raxmlHPC -f a -T 20 -s "$SCLERACTINIA" -n pseudoalteromonas-scleractinia -m GTRGAMMA -p 1234 -x 500 -#100 -w "$OUTPUT" > job_${SLURM_JOB_ID}.log
