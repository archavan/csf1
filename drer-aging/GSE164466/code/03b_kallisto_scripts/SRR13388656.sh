#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_SRR13388656
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load Python
module load kallisto

mkdir -p /home/arc78/csf1/drer-aging/GSE164466/data/02_kallisto/SRR13388656

cd ~/palmer_scratch/geo-data/GSE164466

kallisto quant \
-i ~/csf1/genome/drer/GRCz11/results/Danio_rerio.GRCz11.cdna.all.kallisto.idx \
-t 6 -b 100 -o /home/arc78/csf1/drer-aging/GSE164466/data/02_kallisto/SRR13388656 \
~/palmer_scratch/geo-data/GSE164466/SRR13388656/SRR13388656_1.fastq ~/palmer_scratch/geo-data/GSE164466/SRR13388656/SRR13388656_2.fastq

