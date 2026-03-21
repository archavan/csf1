# background ==================================================================

# We don't need the data for all transcripts, and we have ~1800 accessions with
# each file of about 25MB, which is about 45GB. Instead, we will only read in
# the data for the CSF1 transcripts.

# setup =======================================================================
library(tidyverse)
library(here)
library(parallel)
library(parallelly)

basedir <- here("rbp")
datadir <- here(basedir, "data")
resdir  <- here(basedir, "results")

# data ========================================================================
ids <- read_tsv(here(resdir, "expt-and-tx-quant-accessions.tsv"))

# add replicate id
ids <- ids %>% 
  group_by(rbp, cell_line, expt_type, expt_accession) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup()


# read data ===================================================================
csf1_tx_ids <- c("ENST00000527192.5", 
                 "ENST00000525659.5", 
                 "ENST00000357302.8", 
                 "ENST00000329608.10", 
                 "ENST00000488198.5", 
                 "ENST00000369802.7", 
                 "ENST00000420111.6", 
                 "ENST00000369801.1", 
                 "ENST00000526001.1")

read_csf1_tpm <- function(tsv_url) {
  kal <- read_tsv(tsv_url, n_max = 10000) # CSF1 transcripts appear before row 10000
  kal_csf1 <- filter(kal, str_detect(target_id, "\\|CSF1\\|"))
  kal_csf1 <- separate(kal_csf1, 
                       target_id, 
                       into = c("ensembl_transcript_id", "ensembl_gene_id", rep(NA, 7)), 
                       sep = "\\|")
  kal_csf1 <- kal_csf1[, c("ensembl_transcript_id", "ensembl_gene_id", "tpm")]
  stopifnot(all(csf1_tx_ids %in% kal_csf1$ensembl_transcript_id))
  return(kal_csf1)
}

dat <- ids %>% 
  mutate(tpm = mclapply(href, read_csf1_tpm, mc.cores = availableCores()))

dat <- dat %>% 
  unnest(cols = tpm)

# write data ==================================================================
write_tsv(dat, here(resdir, "tpm_csf1-transcripts.tsv"))

# end =========================================================================
