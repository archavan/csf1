# 1. Extract kallisto mapping stats from json output
# 2. Plot mapping summary stats

# setup =======================================================================
library(tidyverse)
library(jsonlite)
library(glue)
library(ggrepel)
library(here)

basedir <- here("bulk/BafA1.MEF")

source(here(basedir, "code/utilities_ggplot.R"))

# functions ===================================================================
#' Read run_info.json from kallisto output
#'
#' @param dir Name of the directory in which run_info.json is found
#' @param path Path to the seqfile_id directory, without the trailing slash
#'
#' @return A dataframe with one row, with json fields as columns

read_kallisto_json <- function(dir, path = here(basedir, "data/02_kallisto")) {
  read_json(path = glue("{path}/{dir}/run_info.json")) %>%
    as.data.frame() %>%
    select(-call, -start_time)
}

# metadata ====================================================================
metadata <- read_csv(here(basedir, "data/samples.csv"))

# compile kallisto stats ======================================================
kalstats <- metadata %>%
  mutate(json = map(sample_id, ~ read_kallisto_json(.x))) %>%
  unnest_wider(col = json) %>%
  arrange(sample_id)

# write to file
fs::dir_create(here(basedir, "results/02_kallisto"))
write_csv(
  kalstats,
  here(basedir, "results/02_kallisto/kallisto_summary-stats.csv")
)

# plot kallisto statistics ====================================================
## numbers --------------------------------------------------------------------
n_reads <- kalstats %>%
  select(n_processed, n_pseudoaligned, n_unique, sample_id) %>%
  arrange(n_processed) %>%
  mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "statistic",
    names_prefix = "n_",
    values_to = "number_of_reads"
  ) %>%
  ggplot(aes(x = number_of_reads, y = sample_id, colour = statistic)) +
  geom_point(size = 2.5) +
  theme_half_open_grid()

ggsave(
  here(basedir, "results/02_kallisto/n_reads_aligned.pdf"),
  n_reads,
  width = 6.5,
  height = 3.2,
  units = "in"
)

## percentages -----------------------------------------------------------------
p_reads <- kalstats %>%
  select(p_pseudoaligned, p_unique, sample_id) %>%
  arrange(p_pseudoaligned) %>%
  mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "statistic",
    names_prefix = "p_",
    values_to = "percent_reads"
  ) %>%
  ggplot(aes(x = percent_reads, y = sample_id, color = statistic)) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(10, 100, 10)) +
  theme_half_open_grid()

ggsave(
  here(basedir, "results/02_kallisto/percent_reads_aligned.pdf"),
  p_reads,
  width = 6.5,
  height = 3.2,
  units = "in"
)

# end =========================================================================
