# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(janitor)

gse <- "GSE66362" # ***MAKE SURE THIS IS CORRECT***

basedir <- here("drer-aging", gse)
genomedir <- here("genome/drer/GRCz11")

# metadata and genomic data ===================================================
srp_meta <- read_tsv(here(basedir, "data/sample-metadata.tsv")) %>% 
  select(run, title, library_layout)

genes <- read_csv(here(genomedir, "results/gene-id_gene-name_transcript-id.csv"))

# fix sample names ============================================================
# *** THIS NEEDS TO DONE MANUALLY FOR EACH GSE ***
srp_meta <- srp_meta %>% 
  mutate(name_new = make_clean_names(title)) %>% 
  mutate(name_new = str_replace(name_new, "3_75n_m", "3.75nM"))

# compile kallisto output =====================================================
kallisto_read_abundance <- function(path, ids) {
  x <- list()
  for (i in ids) {
    x[[i]] <- readr::read_delim(glue("{path}/{i}/abundance.tsv"), 
                                delim = "\t", 
                                quote = "", 
                                col_names = TRUE)
  }
  return(x)
}

kal <- kallisto_read_abundance(path = here(basedir, "data/02_kallisto"), 
                               ids = as.character(srp_meta$run))

# extract tpm and counts data (at transcript level) ===========================
tpm.tx <- imap(kal, 
               ~ select(.x, target_id, tpm) %>% 
                 rename({{ .y }} := tpm,
                        transcript_id = target_id)) %>% 
  reduce(left_join, by = "transcript_id")

counts.tx <- imap(kal,
                  ~ select(.x, target_id, est_counts) %>%
                    mutate(est_counts = round(est_counts)) %>% 
                    rename({{ .y }} := est_counts,
                           transcript_id = target_id)) %>% 
  reduce(left_join, by = "transcript_id")

# summarize to gene level =====================================================
summarize_to_gene_level <- function(tx) {
  tx %>% 
    inner_join(select(genes,
                      ensembl_gene_id, 
                      ensembl_transcript_id_version), 
               by = c("transcript_id" = "ensembl_transcript_id_version")) %>% 
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
  left_join(x, 
            (genes %>% 
               select(ensembl_gene_id, external_gene_name) %>% 
               distinct()),
            by = "ensembl_gene_id") %>% 
    relocate(external_gene_name, .after = ensembl_gene_id)
}

tpm.gn <- add_gene_names(tpm.gn)
counts.gn <- add_gene_names(counts.gn)

# add gene ids and gene names to transcript-level files =======================
add_gid_gname <- function(x) {
  left_join(x,
            genes,
            by = c("transcript_id" = "ensembl_transcript_id_version")) %>%
    relocate(external_gene_name) %>% 
    relocate(ensembl_gene_id)
}

tpm.tx <- add_gid_gname(tpm.tx)
counts.tx <- add_gid_gname(counts.tx)

# change names ================================================================
names_dict <- srp_meta %>% 
  select(run, name_new) %>% 
  deframe()

replace_srr_with_name <- function(df) {
  df %>% 
    rename_with(function(x) names_dict[x], .cols = starts_with("SRR"))
}

tpm.gn <- replace_srr_with_name(tpm.gn)
counts.gn <- replace_srr_with_name(counts.gn)
tpm.tx <- replace_srr_with_name(tpm.tx)
counts.tx <- replace_srr_with_name(counts.tx)

# subset to csf1 genes ========================================================
subset_to_csf1 <- function(df) {
  df %>% 
    filter(external_gene_name %in% c("csf1a", "csf1b"))
}

tpm.gn.csf1 <- subset_to_csf1(tpm.gn)

# write files =================================================================
fs::dir_create(here(basedir, "data/03_compiled"))

write_csv(tpm.gn,      here(basedir, "data/03_compiled/tpm.csv"))
write_csv(tpm.gn.csf1, here(basedir, "data/03_compiled/tpm_csf1.csv"))

write_tsv(srp_meta, here(basedir, "data/03_compiled/sample-names-dictionary.tsv"))

# write_csv(tpm.tx,          here(basedir, "data/03_compiled/tpm_by-transcript-id.csv"))
# write_csv(counts.tx,       here(basedir, "data/03_compiled/counts_by-transcript-id.csv"))
# write_csv(counts.gn,       here(basedir, "data/03_compiled/counts.csv"))

### end =======================================================================
