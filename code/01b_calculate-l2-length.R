# Arun Chavan
# Started: 2023-07-28

# background ==================================================================

# Calculate L2 length for all isoforms.

# setup =======================================================================
library(tidyverse)

# data ========================================================================
ids <- read_rds("results/isoforms/domains.rds")

# process domain info =========================================================
# clean up columns
ids <- ids %>% 
  mutate(domains = map(domains, ~ select(.x, type, id, start, end)))

# filter to keep pids for which both CSF1 and TM domains were found
ids <- ids %>% 
  filter(map_lgl(domains, 
                 function(x) {
                   csf1_found <- any(x$id %in% c("PF05337"))
                   tm_found <- any(x$id %in% c("TMhelix"))
                   return(csf1_found & tm_found)
                 })) 

# keep only those two domains
ids <- ids %>% 
  mutate(domains = map(domains, ~ filter(.x, id %in% c("PF05337", "TMhelix"))))

# make sure only one CSF1 domain 
ids <- ids %>% 
  filter(map_int(domains, ~ sum(str_detect(.x$id, "PF05337"))) == 1)

# some pids have more than one TM domains. We want to keep only the distal TM,
# the one that comes after the CSF1 domain
check_if_tm_distal <- function(df) {
  csf1_end <- df$end[[which(df$id == "PF05337")]]
  stopifnot(length(csf1_end) == 1)
  df %>% 
    mutate(tm_distal = case_when(
      id == "TMhelix" & (start >  csf1_end) ~ TRUE,
      id == "TMhelix" & (start <= csf1_end) ~ FALSE,
      TRUE ~ NA
    ))
}

ids <- ids %>% 
  mutate(domains = map(domains, check_if_tm_distal)) %>% 
  mutate(domains = map(domains, ~ .x %>% 
                         filter(id == "PF05337" | tm_distal))) 

# After the filtering steps, ensure that all pids have a CSF and a TM domain
ids <- ids %>% 
  filter(map_lgl(domains, 
                 function(x) {
                   csf1_found <- any(str_detect(x$id, "PF05337"))
                   tm_found <- any(str_detect(x$id, "TMhelix"))
                   return(csf1_found & tm_found)
                 }))

# make sure only one CSF1 domain and only one TM domain
ids <- ids %>% 
  filter(map_int(domains, ~ sum(str_detect(.x$id, "PF05337"))) == 1) %>% 
  filter(map_int(domains, ~ sum(str_detect(.x$id, "TMhelix"))) == 1)

write_rds(ids, "results/isoforms/domains_filtered.rds")

# calculate the length of L2 linker domain ====================================
calc_l2_length <- function(df) {
  csf1_end <- df$end[[which(df$id == "PF05337")]]
  tm_start <- df$start[[which(df$id == "TMhelix")]]
  return(tm_start - csf1_end + 1)
}

ids <- ids %>% 
  mutate(l2_len = map_int(domains, calc_l2_length))

write_tsv(select(ids, -domains), "results/isoforms/l2-lengths.tsv")

# end =========================================================================
