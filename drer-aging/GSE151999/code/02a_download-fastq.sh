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

cd ~/palmer_scratch/geo-data/GSE151999

for Sample in SRR11951228 SRR11951229 SRR11951230 SRR11951231 SRR11951232 SRR11951233 SRR11951234 SRR11951235 SRR11951236 SRR11951237 SRR11951238 SRR11951239 SRR11951240 SRR11951241 SRR11951242 SRR11951243 SRR11951244 SRR11951245 SRR11951246 SRR11951247 SRR11951248 SRR11951249 SRR11951250 SRR11951251
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

