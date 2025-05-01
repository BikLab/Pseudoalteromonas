#!/bin/sh

#SBATCH --job-name="genomes"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -e genomes.err-%N
#SBATCH -o genomes.out-%N

module load Miniconda3
source activate /home/ad14556/conda-envs/envs/ncbi_datasets/

cd /home/ad14556/pangenome-large

while read line
    do
        id=$(echo $line | cut -d ' ' -f1) # get accession id from metadata
        name=$(echo $line | cut -d ' ' -f2) # get sample id from metadata
        datasets download genome accession "$id" --include genome # download genome
        unzip -o ncbi_dataset.zip # unzip file
        mv ncbi_dataset/data/"$id"/GCF*.fna data/genomes/"$name".fna # move and rename genome fasta file
#        mv ncbi_dataset/data/"$id"/cds*.fna data/coding-sequences/"$name"-CDS.fna
#        mv ncbi_dataset/data/"$id"/*.gff data/general-feature-files/"$name".gff
#        mv ncbi_dataset/data/"$id"/*.faa data/proteins/"$name"-PROTEINS.faa
        rm ncbi_dataset.zip # remove zipped file
        rm -r ncbi_dataset/ # removed directory
    done < metadata/pseudoalteromonas-large-dataset.txt

