#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=download_fastq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu


module load SRA-Toolkit

cd ~/palmer_scratch/geo-data/PRJEB24858

for Sample in ERR2304198 ERR2304199 ERR2304200 ERR2304201 ERR2304202 ERR2304203 ERR2304204 ERR2304205 ERR2304206 ERR2304207 ERR2304208 ERR2304209
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

