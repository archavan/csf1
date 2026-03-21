# Arun Chavan
# Started: 2023-07-27

# background ==================================================================

# Here we download:
# 1. all CSF1 gene ids
# 2. all isoforms of those gene ids
# 3. domain information for all those isoforms

# setup =======================================================================
library(tidyverse)
library(here)
library(httr)
library(xml2)
library(jsonlite)

basedir <- here("gene-trees")
datadir <- here(basedir, "data")
resdir  <- here(basedir, "results")

# data ========================================================================
tr <- ape::read.tree(here(datadir, "ensembl", "CSF1_gene-tree.nh"))

# get gene ids from protein ids ===============================================
ids <- tibble(pid = tr$tip.label)
ids <- ids %>% 
  filter(str_detect(pid, "^ENS")) %>% 
  mutate(pid = str_replace_all(pid, "_\\w{4}_", "")) 

get_parent_id <- function(id) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/lookup/id/", id, "?")
  res <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(res)
  return(content(res)$Parent)
}

ids <- ids %>% 
  mutate(tid = map_chr(pid, get_parent_id, .progress = TRUE)) %>% 
  mutate(gid = map_chr(tid, get_parent_id, .progress = TRUE))

ids <- ids %>% 
  select(gid) %>% 
  distinct()

# download id info ============================================================
get_id_info <- function(id) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/lookup/id/", id, "?;expand=TRUE")
  res <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(res)
  return(content(res))
}

ids <- ids %>% 
  mutate(res = map(gid, get_id_info, .progress = TRUE)) 

fs::dir_create(here(resdir, "isoforms"))
write_rds(ids, here(resdir, "isoforms", "gid-info.rds"))

# pull id info from GET results ===============================================
pull_species <- function(x) x$species

pull_protein_variants <- function(x) {
  n_transcripts <- length(x$Transcript)
  res <- vector("list", n_transcripts)
  for (i in seq_len(n_transcripts)) {
    res[[i]]$tid <- x$Transcript[[i]]$id
    res[[i]]$display_name <- x$Transcript[[i]]$display_name
    res[[i]]$biotype <- x$Transcript[[i]]$biotype
    res[[i]]$pid <- x$Transcript[[i]]$Translation$id
    res[[i]]$length_aa <- x$Transcript[[i]]$Translation$length
    res[[i]] <- as.data.frame(res[[i]])
  }
  bind_rows(res)
} 

ids <- ids %>% 
  mutate(species = map_chr(res, pull_species)) %>% 
  mutate(vars = map(res, pull_protein_variants, .progress = TRUE))

ids <- ids %>% 
  select(-res) %>% 
  relocate(species) %>% 
  arrange(species) %>% 
  unnest(vars) 

ids <- filter(ids, !is.na(pid))

write_tsv(ids, here(resdir, "isoforms", "isoform-ids.tsv"))

# get domain info =============================================================
get_domains <- function(pid) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/overlap/translation/", pid, "?")
  response <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(response)
  res <- fromJSON(toJSON(content(response)))
  return(res)
}

ids <- ids %>% 
  filter(biotype == "protein_coding") %>% 
  mutate(domains = map(pid, get_domains, .progress = TRUE))

write_rds(ids, here(resdir, "isoforms", "domains.rds"))

# end =========================================================================

