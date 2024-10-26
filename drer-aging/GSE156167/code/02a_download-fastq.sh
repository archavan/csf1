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

cd ~/palmer_scratch/geo-data/GSE156167

for Sample in SRR12446284 SRR12446285 SRR12446286 SRR12446287 SRR12446288 SRR12446289 SRR12446290 SRR12446291 SRR12446292 SRR12446293 SRR12446294 SRR12446295
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

