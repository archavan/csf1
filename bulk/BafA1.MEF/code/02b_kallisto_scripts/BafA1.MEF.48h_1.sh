#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_BafA1.MEF.48h_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p ~/csf1/bulk/BafA1.MEF/data/02_kallisto/BafA1.MEF.48h_1

cd ~/palmer_scratch/data/bulk-rna/csf1/BafA1.MEF/01_fastq

kallisto quant \
-i ~/csf1/genome/mmus/data/ensembl/v111/kallisto-index/Mus_musculus.GRCm39.cdna.all.kallisto.idx \
-t 6 -b 100 -o ~/csf1/bulk/BafA1.MEF/data/02_kallisto/BafA1.MEF.48h_1 --single -l 200 -s 20 \
RBTQJ6_11_48h_1.fastq.gz

