library(tidyverse)
library(here)
library(bioseq)
library(ape)
library(ggtree)
library(ggmsa)
library(Biostrings)
library(aplot)

basedir <- here("gene-trees")
datadir <- here(basedir, "data")
resdir  <- here(basedir, "results")

tr <- read.tree(here(resdir, "gene-tree", "csf1-tree-2_only-tips-with-l2-lengths.nwk"))
ids <- read_tsv(here(resdir, "isoforms", "isoform-ids.tsv"))
ids <- ids %>% 
  filter(pid %in% tr$tip.label) %>% 
  mutate(species = str_replace_all(species, "_", " ") %>% 
           str_to_sentence())
  
spp <- read_tsv(here(resdir, "l2-length", "species-filtered.tsv")) %>%
  select(common_name, scientific_name)

ids <- left_join(ids, spp, by = c("species" = "scientific_name")) %>% 
  relocate(common_name)

ids <- filter(ids, !is.na(common_name))
ids <- mutate(ids, newlab = paste0(pid, " - ", common_name))

ids <- select(ids, pid, newlab)
pid_newlab_mapping <- deframe(ids)

tr <- keep.tip(tr, names(pid_newlab_mapping))
tr$tip.label <- pid_newlab_mapping[tr$tip.label]

alfasta <- read_fasta(here(datadir, "ensembl", "CSF1_alignment.fa"), type = "AA")
alfasta <- alfasta[str_detect(names(alfasta), "^ENS")]
names(alfasta) <- str_replace(names(alfasta), "(^ENS.+?)_.*", "\\1")
alfasta <- alfasta[which(names(alfasta) %in% ids$pid)]
names(alfasta) <- pid_newlab_mapping[names(alfasta)]
alfasta <- alfasta[tr$tip.label]

fs::dir_create(here(resdir, "msa"))
write_fasta(alfasta, here(resdir, "msa", "msa_filtered.fasta"))


al <- here(resdir, "msa", "msa_filtered.fasta")
alf <- read.FASTA(al, type = "AA")


p_al <- ggmsa(al, border = NULL, seq_name = TRUE, ignore_gaps = TRUE, none_bg = TRUE, use_dot = TRUE) +
  scale_y_discrete(position = "right") +
  theme(axis.text.y.right = element_text(size = 2))

ggsave(here(resdir, "msa", "msa-plot.pdf"), p_al, width = 40, height = 20)

renv::install("bioc::msa")
library(msa)

msa::msaPrettyPrint(al)
 
p <- ggtree(tr) + geom_tiplab(size = 3)
msaplot(p, alf, offset = 3, width = 2)


library(seqvisr)
doms <- read_rds(here(resdir, "isoforms", "domains_filtered.rds"))
msavisr(al, myref = "ENSP00000327513 - Human", refontop = FALSE, basecolors = c("grey", "black", "white"))

library(msa)

# setup =======================================================================
library(httr)
library(xml2)
library(jsonlite)
library(tidyverse)
library(bioseq)

# rest api ====================================================================
get_protein_fasta <- function(pid) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", pid, "?")
  url <- paste(server, ext, sep = "")
  response <- GET(url, content_type("text/plain")) # ca also use text/x-fasta
  stop_for_status(response)
  return(content(response, "text"))
}

ids <- ids %>% 
  mutate(seq = map_chr(pid, get_protein_fasta))

fa <- set_names(aa(ids$seq), ids$newlab)
write_fasta(fa, here(resdir, "msa", "protein-aa-sequences.fasta"))

seqs <- readAAStringSet(here(resdir, "msa", "protein-aa-sequences.fasta"))
seqs <- seqs[tr$tip.label]

aln <- msa(seqs)

print(aln, show = "complete")

msaPrettyPrint(aln, 
               output = "pdf", 
               showNames = "right",
               showNumbering = "none",
               alFile = here(resdir, "msa", "msa.fasta"),
               showLogo = "bottom", 
               showConsensus = "top",
               askForOverwrite = FALSE, 
               verbose = FALSE, 
               consensusColors = "ColdHot", 
               showLegend = FALSE, 
               shadingMode = "similar",
               file = here(resdir, "msa", "msa.pdf"),
               furtherCode = c("\\defconsensus{.}{lower}{upper}",
                               "\\showruler{1}{top}"))

# write fasta =================================================================

