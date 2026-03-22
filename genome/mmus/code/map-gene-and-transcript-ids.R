# transcript id to gene id mapping
# Arun Chavan
# Started: 2024-04-30

# background ==================================================================
# We use the ensembl cDNA fastq file's headers to extract transcript id and gene
# id mapping.

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)

genomedir <- here("genome/mmus/data/ensembl/v111")
fs::dir_create(genomedir)

# data ========================================================================
fq <- read.FASTA(
  "https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
)

# extract info
ids <- tibble(headers = names(fq)) %>%
  mutate(
    ensembl_gene_id = str_replace(
      headers,
      ".* gene:(ENSMUSG[^\\s\\.]+).*",
      "\\1"
    ),
    external_gene_name = ifelse(
      str_detect(headers, "gene_symbol"),
      str_replace(headers, ".* gene_symbol:([^\\s]+).*", "\\1"),
      ""
    ),
    ensembl_transcript_id_version = str_replace(
      headers,
      "(^ENSMUST.+?) .*",
      "\\1"
    )
  ) %>%
  select(-headers)

write_csv(ids, here(genomedir, "mmus_gene-id_gene-name_transcript-id.csv"))

# end =========================================================================
