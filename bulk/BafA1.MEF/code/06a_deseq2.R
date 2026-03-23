# Arun Chavan
# Started: 2026-03-22

# Setup =======================================================================
library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(patchwork)

basedir <- here("bulk/BafA1.MEF")
resdir <- here(basedir, "results/06a_deseq")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

# Data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))
counts <- read_csv(here(basedir, "data/03_compiled/counts.csv"))
samples <- read_csv(here(basedir, "data/samples.csv"))

samples <- rename(samples, tx = duration)

genes <- select(tpm, ensembl_gene_id, external_gene_name)

# DESeq: pairwise =============================================================
# set up receptacle
pw <- tibble(
  x = c("0h", "0h", "0h", "0h", "0h"),
  y = c("4h", "8h", "12h", "24h", "48h"),
  coef = paste0("tx_", y, "_vs_", x)
)

p_cutoff <- 0.05

add_gene_names <- function(x) {
  left_join(x, genes, by = "ensembl_gene_id") %>%
    relocate(external_gene_name, .after = ensembl_gene_id)
}

# subset metadata and counts
pw <- pw %>%
  mutate(md = map2(x, y, ~ filter(samples, tx %in% c(.x, .y)))) %>%
  mutate(
    cts = map(
      md,
      ~ select(counts, all_of(c("ensembl_gene_id", .x$sample_id)))
    )
  )

# set up factor levels
for (i in seq_len(nrow(pw))) {
  pw$md[[i]] <- pw$md[[i]] %>%
    mutate(tx = as.character(tx)) %>%
    mutate(tx = fct(tx, c(pw$x[i], pw$y[i])))
}

pw <- pw %>%
  mutate(
    dds = map2(
      cts,
      md,
      ~ DESeqDataSetFromMatrix(
        countData = as.data.frame(.x),
        colData = as.data.frame(.y),
        design = ~tx,
        tidy = TRUE
      )
    )
  ) %>%
  mutate(dds = map(dds, ~ DESeq(.x, test = "LRT", reduced = ~1)))

pw <- pw %>%
  mutate(
    res_shrunken_raw = map2(
      dds,
      coef,
      ~ lfcShrink(.x, coef = .y, type = "apeglm")
    )
  )

pw <- pw %>%
  mutate(
    res_shrunken = map(
      res_shrunken_raw,
      ~ .x %>%
        as.data.frame() %>%
        filter(!is.na(padj)) %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        add_gene_names()
    )
  )

pw <- pw %>%
  mutate(
    res_raw = pmap(., function(x, y, dds, ...) {
      results(dds, contrast = c("tx", y, x), alpha = p_cutoff)
    })
  )

pw <- pw %>%
  mutate(
    res = map(
      res_raw,
      ~ as.data.frame(.x) %>%
        filter(!is.na(padj)) %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        add_gene_names()
    )
  )

pw <- pw %>%
  mutate(
    res_de = map(
      res,
      ~ .x %>%
        filter(padj < p_cutoff) %>%
        arrange(desc(abs(log2FoldChange)))
    )
  )

pw <- pw %>%
  mutate(
    res_up = map(res_de, ~ filter(.x, log2FoldChange > 0)),
    res_down = map(res_de, ~ filter(.x, log2FoldChange < 0))
  )

## p value distributions ====

pw <- pw |>
  mutate(
    pval_hist = map2(
      res,
      coef,
      ~ ggplot(.x, aes(pvalue)) +
        geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
        labs(title = .y) +
        theme_bw_2()
    )
  )


## extract normalized counts ====

pw <- pw %>%
  mutate(
    norm_cts = map(
      dds,
      ~ counts(.x, normalized = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        add_gene_names()
    )
  )

## volcano plots ====

plot_volcano <- function(.data, .title) {
  pdat <- .data
  ptitle <- .title

  pdat$DE <- pdat$padj < p_cutoff
  pdat$de_score <- abs(pdat$log2FoldChange) * -log10(pdat$padj)
  pdat <- arrange(pdat, desc(de_score))
  labelgenes <- pdat %>%
    filter(padj < p_cutoff) %>%
    slice_max(de_score, n = 500) %>%
    pull(external_gene_name)
  pdat$label <- ifelse(
    pdat$external_gene_name %in% labelgenes,
    pdat$external_gene_name,
    NA_character_
  )

  ggplot(pdat, aes(log2FoldChange, -log10(padj))) +
    geom_vline(xintercept = 0, linetype = 2, color = "red") +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    geom_point(aes(color = DE), size = 1) +
    geom_text_repel(aes(label = label)) +
    scale_color_manual(values = c("grey", "black")) +
    ggtitle(ptitle) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    )
}

pw <- mutate(pw, volcano = map2(res_shrunken, coef, plot_volcano))


## write files ====

fs::dir_create(here(resdir, "pairwise"))

for (i in seq_len(nrow(pw))) {
  pairdir <- here(resdir, "pairwise", paste0(pw$y[[i]], "_vs_", pw$x[[i]]))
  fs::dir_create(pairdir)

  ggsave(
    here(pairdir, "pvalue-distribution.pdf"),
    pw$pval_hist[[i]],
    width = 4,
    height = 3.5
  )

  ggsave(
    here(pairdir, "volcano-plot.pdf"),
    pw$volcano[[i]],
    width = 9,
    height = 7
  )

  write_csv(pw$norm_cts[[i]], here(pairdir, "normalized-counts.csv"))

  write_csv(pw$res_shrunken[[i]], here(pairdir, "deseq_out_shrunken.csv"))
  write_csv(pw$res[[i]], here(pairdir, "deseq_out.csv"))
  write_csv(pw$res_de[[i]], here(pairdir, "genes_de.csv"))
  write_csv(
    pw$res_up[[i]],
    here(pairdir, paste0("genes_", pw$y[[i]], "-up.csv"))
  )
  write_csv(
    pw$res_down[[i]],
    here(pairdir, paste0("genes_", pw$x[[i]], "-up.csv"))
  )

  rm(pairdir)
}

# end =========================================================================
