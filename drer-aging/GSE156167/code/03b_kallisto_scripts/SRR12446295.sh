#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_SRR12446295
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load Python
module load kallisto

mkdir -p /home/arc78/csf1/drer-aging/GSE156167/data/02_kallisto/SRR12446295

cd ~/palmer_scratch/geo-data/GSE156167

kallisto quant \
-i ~/csf1/genome/drer/GRCz11/results/Danio_rerio.GRCz11.cdna.all.kallisto.idx \
-t 6 -b 100 -o /home/arc78/csf1/drer-aging/GSE156167/data/02_kallisto/SRR12446295 \
~/palmer_scratch/geo-data/GSE156167/SRR12446295/SRR12446295_1.fastq ~/palmer_scratch/geo-data/GSE156167/SRR12446295/SRR12446295_2.fastq

