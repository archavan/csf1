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

cd ~/palmer_scratch/geo-data/GSE126096

for Sample in SRR8530750 SRR8530751 SRR8530752 SRR8530753 SRR8530754 SRR8530755 SRR8530756 SRR8530757
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

