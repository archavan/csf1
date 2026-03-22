# Arun Chavan
# Started: 2024-05-01

# Setup =======================================================================
library(tidyverse)
library(here)
library(ggrepel)
library(pvclust)

basedir <- here("bulk/BafA1.MEF")
genomedir <- here("genome/mmus/data/ensembl/v111")

resdir <- here(basedir, "results/05_clustering_pca")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"), local = TRUE)

# Data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))
metadata <- read_csv(here(basedir, "data/samples.csv"))

# Clustering ==================================================================

pv <- tpm %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), sqrt)) %>%
  pvclust(
    method.hclust = "complete",
    method.dist = "corr",
    nboot = 1000,
    parallel = TRUE
  )

pdf(here(resdir, "pvclust.pdf"), 8, 6)
plot(pv, print.num = FALSE, print.pv = "bp")
dev.off()

# PCA =========================================================================

pca_out <- tpm %>%
  select(-external_gene_name) %>%
  column_to_rownames("ensembl_gene_id") %>%
  sqrt() %>%
  .[which(rowSums(.) > 0), ] %>%
  t() %>%
  prcomp(x = ., center = TRUE, scale. = TRUE)

pc_space <- pca_out$x %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  right_join(y = metadata, by = "sample_id") |>
  mutate(
    sample = fct(sample, paste0("BafA1.MEF.", c(0, 4, 8, 12, 24, 48), "h"))
  )

var_explained <- round(
  summary(pca_out)$`importance`["Proportion of Variance", ] * 100
)

pc_space |>
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = sample), color = "black", size = 4, shape = 21) +
  scale_fill_viridis_d(direction = -1) +
  labs(
    x = paste0("PC1 (", var_explained["PC1"], "%)"),
    y = paste0("PC2 (", var_explained["PC2"], "%)")
  ) +
  coord_equal() +
  theme_bw_2()

ggsave(here(resdir, "pca_1-2.pdf"), plot = last_plot(), width = 6, height = 4)

# end =========================================================================
