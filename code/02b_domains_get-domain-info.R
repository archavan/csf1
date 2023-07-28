# Arun Chavan
# Started: 2023-07-27

# background ==================================================================

# Download domain information from Ensembl REST API

# setup =======================================================================
library(tidyverse)
library(httr)
library(xml2)
library(jsonlite)

# data ========================================================================
tips <- read_tsv("results/csf1-tree-1_tips.tsv")

# get domains =================================================================
get_domains <- function(pid) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/overlap/translation/", pid, "?")
  response <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  res <- fromJSON(toJSON(content(response)))
    
  if (http_status(response)$category == "Success") {
    return(res)
  } else {
    cat("Failed to fetch sequences. Status code:", http_status(response)$message, "\n")
    return(NULL)
  }
}

tips <- tips %>% 
  mutate(domains = map(pid, get_domains))

# filter to pids for which both CSF and TM domains were found
domains <- tips %>% 
  select(pid, domains) %>% 
  unnest_longer(domains) %>% 
  unpack(-pid) %>% 
  select(pid, type, id, start, end) %>% 
  unnest(-pid) %>% 
  nest(.by = pid) %>% 
  filter(map_lgl(data, 
                 function(x) {
                   csf1_found <- any(x$id %in% c("PF05337"))
                   tm_found <- any(x$id %in% c("TMhelix"))
                   return(csf1_found & tm_found)
                 })) %>% 
  mutate(data = map(data, 
                    ~ filter(.x, id %in% c("PF05337", "TMhelix"))))

# make sure only one CSF1 domain 
domains <- domains %>% 
  filter(map_int(data, ~ sum(str_detect(.x$id, "PF05337"))) == 1)

# some pids have more than one TM domains. We want to keep only the distal TM,
# the one that comes after the CSF1 domain
check_if_tm_distal <- function(df) {
  csf1_end <- df$end[which(df$id == "PF05337")]
  stopifnot(length(csf1_end) == 1)
  df %>% 
    mutate(tm_distal = case_when(
      id == "TMhelix" & start >  csf1_end ~ TRUE,
      id == "TMhelix" & start <= csf1_end ~ FALSE,
      TRUE ~ NA
    ))
}

domains <- domains %>% 
  mutate(data = map(data, check_if_tm_distal)) %>% 
  mutate(data = map(data, ~ .x %>% 
                      filter(id == "PF05337" | tm_distal))) 

# After the filtering steps, ensure that all pids have a CSF and a TM domain
domains <- domains %>% 
  filter(map_lgl(data, 
                 function(x) {
                   csf1_found <- any(str_detect(x$id, "PF05337"))
                   tm_found <- any(str_detect(x$id, "TMhelix"))
                   return(csf1_found & tm_found)
                 }))

domains <- domains %>% 
  filter(map_int(data, ~ sum(str_detect(.x$id, "PF05337"))) == 1) %>% 
  filter(map_int(data, ~ sum(str_detect(.x$id, "TMhelix"))) == 1)

# caculate the lenght of linker 2 domain ======================================
calc_l2_lenght <- function(df) {
  csf1_end <- df$end[which(df$id == "PF05337")]
  tm_start <- df$start[which(df$id == "TMhelix")]
  return(tm_start - csf1_end + 1)
}

l2_len <- domains %>% 
  mutate(l2_length = map_int(data, calc_l2_lenght)) %>% 
  select(pid, l2_length) %>% 
  left_join(select(tips, pid, newlab), by = "pid")
  

# write files =================================================================
fs::dir_create("results/linker2")
write_tsv(l2_len, "results/linker2/l2-lengths.tsv")
write_rds(tips, "results/linker2/domains_raw.rds")
write_rds(domains, "results/linker2/domains_filtered.rds")

# end =========================================================================
