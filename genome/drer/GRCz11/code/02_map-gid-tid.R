# transcript id to gene id mapping
# Arun Chavan

# background ==================================================================
# We use the ensembl cDNA fastq file's headers to extract transcript id and gene
# id mapping.

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)

genomedir <- here("genome/drer/GRCz11")

# data ========================================================================
fq <- read.FASTA("https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz")

# extract info
ids <- tibble(headers = names(fq)) %>% 
  mutate(
    ensembl_gene_id = str_replace(headers, ".* gene:(ENSDARG[^\\s\\.]+).*", "\\1"),
    external_gene_name = ifelse(str_detect(headers, "gene_symbol"),
                                str_replace(headers, ".* gene_symbol:([^\\s]+).*", "\\1"),
                                ""),
    ensembl_transcript_id_version = str_replace(headers, "(^ENSDART.+?) .*", "\\1")
  ) %>% 
  select(-headers)

write_csv(ids, here(genomedir, "results/gene-id_gene-name_transcript-id.csv"))

# end =========================================================================
