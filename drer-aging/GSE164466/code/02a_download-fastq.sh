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

cd ~/palmer_scratch/geo-data/GSE164466

for Sample in SRR13388639 SRR13388640 SRR13388641 SRR13388642 SRR13388643 SRR13388644 SRR13388645 SRR13388646 SRR13388647 SRR13388648 SRR13388649 SRR13388650 SRR13388651 SRR13388652 SRR13388653 SRR13388654 SRR13388655 SRR13388656 SRR13388657 SRR13388658 SRR13388659 SRR13388660 SRR13388661 SRR13388662
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

