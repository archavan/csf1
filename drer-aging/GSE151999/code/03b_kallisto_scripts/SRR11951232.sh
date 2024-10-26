#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_SRR11951232
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load Python
module load kallisto

mkdir -p /home/arc78/csf1/drer-aging/GSE151999/data/02_kallisto/SRR11951232

cd ~/palmer_scratch/geo-data/GSE151999

kallisto quant \
-i ~/csf1/genome/drer/GRCz11/results/Danio_rerio.GRCz11.cdna.all.kallisto.idx \
-t 6 -b 100 -o /home/arc78/csf1/drer-aging/GSE151999/data/02_kallisto/SRR11951232 --single -l 200 -s 20 \
~/palmer_scratch/geo-data/GSE151999/SRR11951232/SRR11951232.fastq

