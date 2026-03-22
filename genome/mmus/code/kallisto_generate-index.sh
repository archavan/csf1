#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p ~/csf1/genome/mmus/data/ensembl/v111
cd ~/csf1/genome/mmus/data/ensembl/v111
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

mkdir -p ~/csf1/genome/mmus/data/ensembl/v111/kallisto-index
cd ~/csf1/genome/mmus/data/ensembl/v111/kallisto-index

kallisto index -i Mus_musculus.GRCm39.cdna.all.kallisto.idx \
~/csf1/genome/mmus/data/ensembl/v111/Mus_musculus.GRCm39.cdna.all.fa.gz
