# Plot CSF1 gene tree
# Arun Chavan
# Started: 2023-07-25

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)
library(ggtree)

basedir <- here("gene-trees")
datadir <- here(basedir, "data")
resdir  <- here(basedir, "results")

# data ========================================================================
tr <- read.tree(here(datadir, "ensembl", "CSF1_gene-tree.nh"))
sp <- read_csv(here(datadir, "ensembl", "ensembl-species.csv")) %>%
  janitor::clean_names() %>%
  select(common_name, scientific_name)

# organize names ==============================================================
## tip labesl from the treefile ----
tips <- tibble(tiplab = tr$tip.label)

str_replace_if <- function(x, pattern, replacement, nomatch = NA_character_) {
  ifelse(str_detect(x, pattern),
         str_replace(x, pattern, replacement),
         nomatch)
} # there are some tips that don't follow the ENS...abbr pattern

tips <- tips %>% 
  mutate(pid  = str_replace_if(tiplab, "(ENS\\w+\\d+)_(\\w{4})_", "\\1"),
         abbr = str_replace_if(tiplab, "(ENS\\w+\\d+)_(\\w{4})_", "\\2")) 

# there is another abbreviation in the ensembl protein id, which we will extract
extract_ensembl_abbr <- function(x) {
  if (str_detect(x, "^ENS\\w{3}")) {
    substr(x, 4, 6)
  } else {
    NA_character_
  }
}

tips <- tips %>% 
  mutate(ens_abbr = map_chr(tiplab, extract_ensembl_abbr))

## species names from the species file ----
construct_abbr <- function(x) {
  # the species name abbr in gene tree tip labels are the first letter of genus,
  # 1-3 letters of species, 1-2 letters of subspecies, and 1 letter of
  # subsubspecies depending
  
  # e.g. Mus musculus = Mmus, Cervus hanglu yarkandensis = Chya
  
  x <- str_split_1(x, "\\s+")
  x[1] <- substr(x[1], 1, 1)
  x[2] <- case_when(length(x) == 2 ~ substr(x[2], 1, 3),
                    length(x) > 2 ~ substr(x[2], 1, 1))
  x[3] <- case_when(length(x) == 2 ~ "",
                    length(x) == 3 ~ substr(x[3], 1, 2),
                    length(x) > 3 ~ substr(x[3], 1, 1))
  x[4] <- case_when(length(x) == 2 ~ "",
                    length(x) == 3 ~ "",
                    length(x) >= 4 ~ substr(x[4], 1, 1))
  x <- head(x, 4)
  paste0(x, collapse = "")
}

sp <- sp %>% 
  mutate(abbr = map_chr(scientific_name, construct_abbr))

sp <- filter(sp, !duplicated(abbr)) # some species have multiple strains

## deal with duplicate abbreviations ---- 

# There are certain species for which the abbr within the tip label are the
# same. Since we will join tip labels and species labels dataframes based on
# abbr, we need to make sure we fix this problem.

ambiguous_abbr <- tips %>% 
  group_by(abbr) %>% 
  summarize(ens_abbr = list(ens_abbr)) %>% 
  mutate(ens_abbr = map(ens_abbr, unique)) %>% 
  filter(map_int(ens_abbr, length) > 1) %>% 
  unnest(cols = c(ens_abbr))

# Manually look these up and correct the entries
disambiguate_abbr <- function(x) {
  case_when(
    x == "CPR" ~ list(common_name = "Australian saltwater crocodile",
                      scientific_name = "Crocodylus porosus"),
    x == "CPO" ~ list(common_name = "Guinea Pig",
                      scientific_name = "Cavia porcellus"),
    x == "MMD" ~ list(common_name = "Pinecone soldierfish",
                      scientific_name = "Myripristis murdjan"),
    x == "MIC" ~ list(common_name = "Mouse Lemur",
                      scientific_name = "Microcebus murinus"),
    x == "PSI" ~ list(common_name = "Chinese softshell turtle",
                      scientific_name = "Pelodiscus sinensis"),
    x == "PSN" ~ list(common_name = "Vaquita",
                      scientific_name = "Phocoena sinus")
  )
}

## join tip labels data and species names data ----

# First separate the unambiguous and ambiguous tips, join, then rbind

ambiguous_tips <- filter(tips, ens_abbr %in% ambiguous_abbr$ens_abbr) %>% 
  mutate(temp = map(ens_abbr, disambiguate_abbr)) %>% 
  unnest_wider(col = "temp")

unambiguous_tips <- filter(tips, !(ens_abbr %in% ambiguous_abbr$ens_abbr)) %>% 
  left_join(sp, by = "abbr")


tips <- bind_rows(ambiguous_tips, unambiguous_tips)

## create display names ----
tips <- tips %>% 
  mutate(newlab = paste0(pid, " - ", common_name)) 

# drop tips ===================================================================
keep_species <- c("Red-bellied piranha", 
                  "Electric eel",
                  # "Channel catfish",
                  "Goldfish",
                  "Common carp",
                  "Zebrafish",
                  "Atlantic herring",
                  "Rainbow trout",
                  "Atlantic salmon",
                  "Northern pike",
                  "Tiger tail seahorse",
                  "Atlantic cod",
                  "Asian bonytongue",
                  # "Barramundi perch",
                  # "Climbing perch",
                  "Pike-perch",
                  "Stickleback",
                  # "European seabass",
                  "Orange clownfish",
                  # "Makobe Island cichlid",
                  # "Nile tilapia",
                  "Turquoise killifish",
                  # "Amazon molly",
                  # "Guppy",
                  # "Indian medaka",
                  # "Fugu",
                  "Spotted gar",
                  "Reedfish",
                  "Chicken",
                  "Zebra finch",
                  # "Australian saltwater crocodile",
                  "Painted turtle",
                  "Green anole",
                  "Golden Hamster",
                  # "Northern American deer mouse",
                  "Mouse",
                  "Rat",
                  "Kangaroo rat",
                  "Squirrel",
                  # "Eurasian red squirrel",
                  "Guinea Pig",
                  "Pika",
                  "Macaque",
                  "Human",
                  "Gibbon",
                  # "Chimpanzee",
                  "Bonobo",
                  "Gorilla",
                  # "Panamanian white-faced capuchin",
                  "Tarsier",
                  # "Mouse Lemur",
                  "Rabbit",
                  "Blue whale",
                  # "American bison",
                  "Goat",
                  # "Siberian musk deer",
                  "Cow",
                  "Pig",
                  # "Arabian camel",
                  "American mink",
                  # "American black bear",
                  "Dog",
                  "Cat",
                  "Horse",
                  "Microbat",
                  "Armadillo",
                  "Elephant",
                  "Opossum",
                  "Platypus")

keep_tips <- tips %>% 
  filter(common_name %in% keep_species) %>% 
  pull(tiplab) %>% 
  unique()

tr1 <- keep.tip(tr, keep_tips)

# replace tip labels
new_tip_labels <- filter(tips, tiplab %in% keep_tips) %>% 
  select(tiplab, newlab) %>% 
  deframe()

tr1$tip.label <- new_tip_labels[tr1$tip.label]

summary(tr1)

plot.phylo(ladderize(tr1), cex = 0.5)

# annoatate duplication nodes based on ensembl tree ===========================
dup_nodes <- c()
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr1, tr1$tip.label[grepl("Zebrafish", tr1$tip.label)]))
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr1, c("ENSSSAP00000150975 - Atlantic salmon",
                                   "ENSSSAP00000108598 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr1, c("ENSSSAP00000154249 - Atlantic salmon",
                                   "ENSSSAP00000149216 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr1, c("ENSCARP00000058590 - Goldfish",
                                   "ENSCARP00000062662 - Goldfish")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr1, c("ENSCARP00000091986 - Goldfish",
                                   "ENSCARP00000075789 - Goldfish")))

tr1.nodecols <- rep('NA', tr1$Nnode)
tr1.nodecols[dup_nodes - length(tr1$tip.label)] <- 'red' 

# plots =======================================================================
# without annotation
plot.phylo(ladderize(tr1), 
           cex = 0.7, 
           use.edge.length = TRUE, 
           edge.width = 1.5, 
           no.margin = TRUE, 
           align.tip.label = TRUE,
           font = 1, 
           label.offset = 0.02)
nodelabels(pch = 19, col = tr1.nodecols, cex = 1.25)

nodes_to_label <- c(`Mammals` = getMRCA(tr1, c("ENSOANP00000011057 - Platypus",
                                               "ENSLAFP00000026014 - Elephant")),
                    `Reptiles` = getMRCA(tr1, c("ENSGALP00010041133 - Chicken",
                                                "ENSACAP00000011757 - Green anole")),
                    `Teleosts` = getMRCA(tr1, c("ENSSFOP00015001355 - Asian bonytongue",
                                                "ENSGACP00000003727 - Stickleback")),
                    `Ray-finned fish` = getMRCA(tr1, c("ENSGACP00000003727 - Stickleback",
                                                       "ENSECRP00000005354 - Reedfish")))

p1 <- ggtree(tr1, ladderize = TRUE) +
  geom_tiplab(align = TRUE, geom = "text", size = 3, offset = 0.01, linesize = 0.25) +
  ggtree::xlim(c(0, 2.25)) +
  geom_nodepoint(colour = tr1.nodecols, size = 3, shape = 19) +
  # geom_hilight(node = nodes_to_label[["Mammals"]], fill = "#7fc97f", alpha = 0.2) +
  # geom_hilight(node = nodes_to_label[["Reptiles"]], fill = "#beaed4", alpha = 0.2) +
  # geom_hilight(node = nodes_to_label[["Ray-finned fish"]], fill = "#ffff99", alpha = 0.2) +
  geom_hilight(node = nodes_to_label[["Teleosts"]], fill = "#fdc086", alpha = 0.2) +
  geom_cladelabel(node = nodes_to_label[["Mammals"]], 
                  label = "Mammals", 
                  align = TRUE, color = "grey30", 
                  hjust = 0.5, offset.text = 0.02,
                  offset = 0.7, angle = 270, barsize = 1, fontsize = 3) +
  geom_cladelabel(node = nodes_to_label[["Reptiles"]], 
                  label = "Reptiles", 
                  align = TRUE, color = "grey30", 
                  hjust = 0.5, offset.text = 0.02,
                  offset = 0.7, angle = 270, barsize = 1, fontsize = 3) +
  geom_cladelabel(node = nodes_to_label[["Teleosts"]], 
                  label = "Teleosts", 
                  align = TRUE, color = "grey30", 
                  hjust = 0.5, offset.text = 0.02,
                  offset = 0.7, angle = 270, barsize = 1, fontsize = 3) +
  geom_cladelabel(node = nodes_to_label[["Ray-finned fish"]], 
                  label = "Ray-finned fish", 
                  align = TRUE, color = "grey30", 
                  hjust = 0.5, offset.text = 0.02,
                  offset = 0.8, angle = 270, barsize = 1, fontsize = 3) +
  labs(caption = "\nRed circles indicate gene duplication nodes.")

fs::dir_create(here(resdir, "gene-tree"))
ggsave(here(resdir, "gene-tree", "csf1-tree-1.pdf"), p1, width = 10, height = 11)

# write files =================================================================
write.tree(tr1, file = here(resdir, "gene-tree", "csf1-tree-1.nwk"))

tips %>%
  filter(newlab %in% tr1$tip.label) %>%
  write_tsv(here(resdir, "gene-tree", "csf1-tree-1_tips.tsv"))

# end =========================================================================
