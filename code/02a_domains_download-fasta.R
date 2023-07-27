# Arun Chavan
# Started: 2023-07-26

# background ==================================================================

# Download multi-fasta file containing the sequences of all protein ids that
# appear in the gene tree.

# setup =======================================================================
library(httr)
library(xml2)
library(jsonlite)
library(tidyverse)
library(bioseq)

# data ========================================================================
tips <- read_tsv("results/csf1-tree-1_tips.tsv")

# rest api ====================================================================
get_protein_fasta <- function(pid) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", pid, "?")
  url <- paste(server, ext, sep = "")
  
  response <- GET(url, content_type("text/plain")) # ca also use text/x-fasta

  if (http_status(response)$category == "Success") {
    return(content(response, "text"))
  } else {
    cat("Failed to fetch sequences. Status code:", http_status(response)$message, "\n")
    stop()
  }
}

tips <- tips %>% 
  mutate(seq = map_chr(pid, ~ get_protein_fasta(.x)))

# write fasta =================================================================
fa <- set_names(aa(tips$seq), tips$newlab)
write_fasta(fa, "results/csf1-tree-1_tips.fasta")

# end =========================================================================
