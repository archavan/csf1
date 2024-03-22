# background ==================================================================

# We have knockdown and control RNA-seq data for ~250 RBPs from
# https://www.nature.com/articles/s41586-020-2077-3. Here we try to create a
# table of accession ids associated with each RBP experiment in order to
# eventually get accession ids for all transcript quatification tsv files from
# the dataset.

# setup =======================================================================
library(tidyverse)
library(here)
library(jsonlite)
library(httr)
library(readxl)
library(parallel)
library(parallelly)

resdir <- here("results/rbp")
fs::dir_create(resdir)

# data ========================================================================
expt_ids <- read_xlsx("data/rbp/41586_2020_2077_MOESM4_ESM.xlsx", 
                     sheet = "KD-RNA-seq") %>% 
  janitor::clean_names() %>% 
  select(rbp, cell_line, rna_seq_exp, control_exp) %>% 
  rename(kd = rna_seq_exp, 
         ctrl = control_exp) %>% 
  pivot_longer(cols = c(kd, ctrl), 
               names_to = "expt_type", values_to = "expt_accession")

# get tx quant file accessions ================================================

get_tx_quant_file_accessions <- function(.exp_id) {
  
  baseurl <- "https://www.encodeproject.org/experiments/"
  exp_url <- paste0(baseurl, .exp_id, "/?format=json")
  
  response <- GET(exp_url)
  x <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  y <- x$files
  
  z <- y %>% 
    filter(output_type == "transcript quantifications",
           assembly == "GRCh38",
           genome_annotation == "V29") %>% 
    select(accession, href, analysis_step_version)
  
  z$analysis_title <- z$analysis_step_version$analysis_step$title
  z <- z[which(grepl("kallisto", z$analysis_title)), ]
  
  z$analysis_step_version <- NULL
  z$analysis_title <- NULL
  
  z$href <- paste0("https://www.encodeproject.org", z$href)
  names(z)[names(z) == "accession"] <- "tx_quant_accession"
  return(z)
}

expt_ids <- expt_ids %>% 
  mutate(tx_acc_df = mclapply(expt_accession, 
                              get_tx_quant_file_accessions, 
                              mc.cores = availableCores()))

expt_ids <- expt_ids %>% 
  unnest(cols = tx_acc_df)

# write files =================================================================

write_tsv(expt_ids, here(resdir, "expt-and-tx-quant-accessions.tsv"))

# end =========================================================================
