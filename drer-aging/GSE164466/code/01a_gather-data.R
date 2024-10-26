# Arun Chavan
# Started: 2024-10-24

# background ==================================================================

# setup =======================================================================
library(GEOfastq)
library(here)
library(tidyverse)
library(glue)

gse <- "GSE164466" 

basedir <- here("drer-aging", gse)

# get samples and sample metadata =============================================
srp_meta <- crawl_gse(gse) %>% 
  extract_gsms() %>% 
  crawl_gsms()

# generate script for downloading fastq files from SRA ========================

run_ids <- srp_meta$run

fs::dir_create(here("~/palmer_scratch/geo-data", gse))

sra_script <- glue(
  "
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

cd ~/palmer_scratch/geo-data/{gse}

for Sample in {paste0(run_ids, collapse = \" \")}
do
	prefetch $Sample
	fasterq-dump ./$Sample --split-3 -O ./$Sample
	rm ./$Sample/$Sample.sra
done

  "
)


# write files =================================================================
fs::dir_create(here(basedir, "data"))
write_tsv(srp_meta, here(basedir, "data/sample-metadata.tsv"))
writeLines(sra_script, here(basedir, "code/02a_download-fastq.sh"))

# end =========================================================================
