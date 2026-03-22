# description =================================================================
# 1. Compile kallisto output for all samples into one file
# 2. Aggregate count and TPM values by gene ids

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

basedir <- here("bulk/BafA1.MEF")
genomedir <- here("genome/mmus/data/ensembl/v111")

# metadata and genomic data ===================================================
metadata <- read_csv(here(basedir, "data/samples.csv"))
genes <- read_csv(here(genomedir, "mmus_gene-id_gene-name_transcript-id.csv"))

# compile kallisto output =====================================================
kallisto_read_abundance <- function(path, ids) {
  x <- list()
  for (i in ids) {
    x[[i]] <- readr::read_delim(
      glue("{path}/{i}/abundance.tsv"),
      delim = "\t",
      quote = "",
      col_names = TRUE
    )
  }
  return(x)
}

kal <- kallisto_read_abundance(
  path = here(basedir, "data/02_kallisto"),
  ids = as.character(metadata$sample_id)
)

# extract tpm and counts data (at transcript level) ===========================
tpm.tx <- imap(
  kal,
  ~ select(.x, target_id, tpm) %>%
    rename({{ .y }} := tpm, transcript_id = target_id)
) %>%
  reduce(left_join, by = "transcript_id")

counts.tx <- imap(
  kal,
  ~ select(.x, target_id, est_counts) %>%
    mutate(est_counts = round(est_counts)) %>%
    rename({{ .y }} := est_counts, transcript_id = target_id)
) %>%
  reduce(left_join, by = "transcript_id")

# summarize to gene level =====================================================
summarize_to_gene_level <- function(tx) {
  tx %>%
    inner_join(
      select(genes, ensembl_gene_id, ensembl_transcript_id_version),
      by = c("transcript_id" = "ensembl_transcript_id_version")
    ) %>%
    relocate(ensembl_gene_id) %>%
    pivot_longer(cols = where(is.double)) %>%
    group_by(ensembl_gene_id, name) %>%
    summarize(value = sum(value)) %>%
    ungroup() %>%
    pivot_wider(id_cols = "ensembl_gene_id")
}

tpm.gn <- summarize_to_gene_level(tpm.tx)
counts.gn <- summarize_to_gene_level(counts.tx)

# add gene names ==============================================================
add_gene_names <- function(x) {
  left_join(
    x,
    (genes %>%
      select(ensembl_gene_id, external_gene_name) %>%
      distinct()),
    by = "ensembl_gene_id"
  ) %>%
    relocate(external_gene_name, .after = ensembl_gene_id)
}

tpm.gn <- add_gene_names(tpm.gn)
counts.gn <- add_gene_names(counts.gn)

# calculate mean tpm ==========================================================
tpm.gn.mean <- tpm.gn %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = c("sample", "rep_id"),
    names_sep = "_",
    values_to = "tpm"
  ) %>%
  group_by(ensembl_gene_id, sample) %>%
  summarize(tpm = mean(tpm)) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "tpm") %>%
  add_gene_names()

# write files =================================================================
fs::dir_create(here(basedir, "data/03_compiled"))
write_csv(tpm.tx, here(basedir, "data/03_compiled/tpm_by-transcript-id.csv"))
write_csv(
  counts.tx,
  here(basedir, "data/03_compiled/counts_by-transcript-id.csv")
)
write_csv(tpm.gn, here(basedir, "data/03_compiled/tpm.csv"))
write_csv(counts.gn, here(basedir, "data/03_compiled/counts.csv"))
write_csv(tpm.gn.mean, here(basedir, "data/03_compiled/tpm_mean.csv"))

### end =======================================================================
