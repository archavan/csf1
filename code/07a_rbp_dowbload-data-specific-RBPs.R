# background ==================================================================

# Download full data for specific RBPs.

# setup =======================================================================
library(tidyverse)
library(here)
library(parallel)
library(parallelly)

resdir <- here("results/rbp/tpm")
fs::dir_create(resdir)

# data ========================================================================
ids <- read_tsv(here("results/rbp/expt-and-tx-quant-accessions.tsv"))

rbps <- c("PPIL4", "MTPAP")

# add replicate id
ids <- ids %>% 
  group_by(rbp, cell_line, expt_type, expt_accession) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup()


# read data ===================================================================

ids <- ids %>% 
  filter(rbp %in% rbps)

read_tpm <- function(tsv_url) {
  kal <- read_tsv(tsv_url) 
  kal <- separate(kal, 
    target_id, 
    into = c("ensembl_transcript_id", "ensembl_gene_id", NA, NA, "transcript_name", "gene_name", rep(NA, 3)), 
    sep = "\\|")
  kal <- kal[, c("ensembl_transcript_id", "ensembl_gene_id", "transcript_name", "gene_name", "tpm")]
  kal <- arrange(kal, ensembl_gene_id)
  return(kal)
}

dat <- ids %>% 
  mutate(tpm = mclapply(href, read_tpm, mc.cores = availableCores()))

# write data ==================================================================
for (i in seq_len(nrow(dat))) {
  write_tsv(dat$tpm[[i]], 
            here(resdir, paste0(dat$rbp[i], 
                                "_", 
                                dat$cell_line[i], 
                                "_", 
                                dat$expt_type[i], 
                                "_", 
                                dat$rep[i], 
                                ".tsv")))
}

# end =========================================================================
