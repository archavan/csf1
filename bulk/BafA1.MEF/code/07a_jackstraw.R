# Arun Chavan
# Started: 2026-03-22

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(factoextra)
library(ggrepel)
library(pheatmap)
library(conflicted)

basedir <- here("bulk/BafA1.MEF")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_jackstraw.R"))

resdir <- here(basedir, "results/07a_jackstraw")
fs::dir_create(resdir)

conflicts_prefer(
  dplyr::filter,
  purrr::map,
  dplyr::select,
  dplyr::count,
  base::as.factor
)

###############################################################################
# data ========================================================================
###############################################################################

md <- read_csv(here(basedir, "data/samples.csv"))
tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))

genes <- tpm |>
  select(ensembl_gene_id, external_gene_name)

###############################################################################
# run PCA =====================================================================
###############################################################################

pca_dat <- tpm |>
  select(-external_gene_name) |>
  column_to_rownames("ensembl_gene_id") |>
  sqrt() |>
  (function(x) x[rowSums(x) > 3, ])()

js <- js_run_pca_with_jackstraw(pca_dat, 10, 1000, ndim = 9, scale. = TRUE)

## tidy output ================================================================
ldngs <- js$pca_out$rotation |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading")

pvals <- js$pval |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "pval")

js_tidy <- left_join(ldngs, pvals, by = c("PC", "ensembl_gene_id")) |>
  group_by(PC) |>
  mutate(padj = p.adjust(pval, method = "BH")) |>
  ungroup() |>
  mutate(PC = fct(PC, paste0("PC", sort(parse_number(unique(PC)))))) |>
  mutate(direction = ifelse(loading < 0, "negative", "positive")) |>
  arrange(PC, direction, -abs(loading)) |>
  left_join(genes, by = "ensembl_gene_id")


add_metadata_to_pc_coords <- function(.pca_out) {
  .pca_out$x |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    left_join(y = md, by = "sample_id")
}

pc_space <- add_metadata_to_pc_coords(js$pca_out) |>
  mutate(duration = fct(duration, paste0(c(0, 4, 8, 12, 24, 48), "h")))

## get sig genes ==============================================================
sig_genes <- js_tidy |>
  filter(padj < 0.05) |>
  group_by(PC, direction) |>
  nest() |>
  ungroup()

## plot variance explained ====================================================

var_explained <- get_eig(js$pca_out) |>
  rownames_to_column("PC") |>
  mutate(PC = gsub("Dim\\.", "PC", PC)) |>
  mutate(variance.percent = round(variance.percent, 1))

var_explained |>
  mutate(PC = str_remove(PC, "PC")) |>
  mutate(PC = as.integer(PC)) |>
  ggplot(aes(PC, variance.percent)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(name = "Variance Explained (%)") +
  theme_bw_2()
ggsave(here(resdir, "pca_var-explained.pdf"), width = 4, height = 3.5)

## plot PC score against segment ==============================================
pc_space |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_score"
  ) |>
  ggplot(aes(duration, PC_score, color = as.factor(rep))) +
  geom_line(aes(group = rep), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "Replicate") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(here(resdir, "pca_pc-score-vs-duration.pdf"), width = 6.5, height = 5.5)

## write significant genes ====================================================
fs::dir_create(here(resdir, "jackstraw-genes"))

for (i in seq_len(nrow(sig_genes))) {
  write_csv(
    sig_genes$data[[i]],
    here(
      resdir,
      "jackstraw-genes",
      glue("{sig_genes$PC[[i]]}_{sig_genes$direction[[i]]}.csv")
    )
  )
}

# end =========================================================================
