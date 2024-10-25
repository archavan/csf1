#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

cd /home/arc78/csf1/genome/drer/GRCz11/results
wget https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz

kallisto index -i Danio_rerio.GRCz11.cdna.all.kallisto.idx Danio_rerio.GRCz11.cdna.all.fa.gz

