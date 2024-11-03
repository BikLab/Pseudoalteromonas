#!/bin/sh

#SBATCH --job-name="prokka"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e prokka.err-%N
#SBATCH -o prokka.out-%N


module load prokka
module load barrnap

# index hmm
hmmpress /scratch/ad14556/pangenome_psalt/database/hmm/Pfam-A.hmm

# run prokka for every pseudoalteromonas genome
for file in /home/ad14556/pangenome-large/data/genomes/*.fna
	do
		base=$(basename $file .fna)
		prokka  --outdir /scratch/ad14556/pangenome-large/results/prokka \
			--proteins /scratch/ad14556/pangenome_psalt/database/fasta/uniprot_sprot.fasta \
			--genus Pseudoalteromonas --kingdom Bacteria --gcode 11 --usegenus \
                        --prefix "$base" \
			--cpus 24 \
			--force --compliant --addgenes $file
	done

#                       --hmms /scratch/ad14556/pangenome_psalt/database/hmm/Pfam-A.hmm \
