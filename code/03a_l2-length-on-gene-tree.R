# Arun Chavan
# Started: 2023-07-29

# background ==================================================================

# The gene tree has protein ids on it. Here we
# 1. get the gene ids from those protein ids
# 2. get the lengths of linker2 domains for all isoforms of those genes
# 3. plot the gene tree along with the lengths of the L2 domains

# setup =======================================================================
library(tidyverse)
library(ape)
library(ggtree)
library(aplot)

# data ========================================================================
tr <- read.tree("results/gene-tree/csf1-tree-1.nwk")
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")

ids <- read_tsv("results/isoforms/l2-lengths.tsv")

spp <- read_csv("data/ensembl/ensembl-species.csv") %>% 
  janitor::clean_names() %>% 
  select(common_name, scientific_name) %>% 
  mutate(common_name = case_when(
    str_detect(scientific_name, "Mus musculus") ~ "Mouse",
    str_detect(scientific_name, "Sus scrofa") ~ "Pig",
    str_detect(scientific_name, "Rattus norvegicus") ~ "Rat",
    str_detect(scientific_name, "Capra hircus") ~ "Goat",
    str_detect(scientific_name, "Canis lupus familiaris") ~ "Dog",
    str_detect(scientific_name, "Gallus gallus") ~ "Chicken",
    TRUE ~ common_name
  )) %>% 
  distinct()

# match tip names =============================================================
tr$tip.label <- str_replace(tr$tip.label, "(^ENS.+?) .*", "\\1")

# For not all proteins that are tips in this tree we could calculate the L2
# length. We need to drop them.
tr <- keep.tip(tr, tr$tip.label[tr$tip.label %in% ids$pid])
write.tree(tr, "results/gene-tree/csf1-tree-2_only-tips-with-l2-lengths.nwk")

pid_gid_mapping <- ids %>% 
  select(pid, gid) %>% 
  distinct() %>% 
  deframe()

tr$tip.label <- pid_gid_mapping[tr$tip.label]

ids <- ids %>% 
  filter(gid %in% tr$tip.label)

# add species name to tip labels
ids <- ids %>% 
  mutate(species = str_replace_all(species, "_", " ")) %>% 
  mutate(species = str_to_sentence(species)) %>% 
  rename(scientific_name = species) %>% 
  left_join(spp, by = "scientific_name") %>% 
  relocate(common_name) %>% 
  mutate(newlab = paste0(gid, " - ", common_name))

gid_newlab_mapping <- ids %>% 
  select(gid, newlab) %>% 
  distinct() %>% 
  deframe()

tr$tip.label <- gid_newlab_mapping[tr$tip.label]

# add isoform info ===========================================================
ids <- ids %>% 
  group_by(gid) %>% 
  mutate(isoform = as.integer(as.factor(pid))) %>% 
  mutate(isoform = as.character(isoform)) %>% 
  ungroup()

# plot L2 lengths =============================================================
p_l2 <- ggplot(ids, aes(l2_len, newlab)) +
  geom_point(shape = 21, fill = "grey20", stroke = 0.5, size = 2.5, color = "white") +
  scale_y_discrete(position = "right") +
  scale_fill_discrete(na.translate = FALSE) +
  labs(x = "amino acids", 
       y = NULL,
       title = "Length of Linker2 domain") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.5),
        axis.text.y.right = element_text(size = 9),
        plot.title = element_text(size = 9))

# plot gene tree ==============================================================
# find duplication nodes
dup_nodes <- c()
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr, tr$tip.label[grepl("Zebrafish", tr$tip.label)]))
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr, c("ENSSSAG00000078129 - Atlantic salmon",
                                  "ENSSSAG00000077389 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr, c("ENSSSAG00000081356 - Atlantic salmon",
                                  "ENSSSAG00000081033 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr, c("ENSCARG00000025558 - Goldfish",
                                  "ENSCARG00000027791 - Goldfish")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr, c("ENSCARG00000035280 - Goldfish",
                                  "ENSCARG00000045570 - Goldfish")))

tr.nodecols <- rep('NA', tr$Nnode)
tr.nodecols[dup_nodes - length(tr$tip.label)] <- 'red' 

# nodes to label ==============================================================
nodes_to_label <- c(
  `Mammals` = getMRCA(tr, c("ENSOANG00000006943 - Platypus",
                            "ENSLAFG00000032354 - Elephant")),
  `Reptiles` = getMRCA(tr, c("ENSGALG00010027734 - Chicken",
                             "ENSACAG00000012025 - Green anole")),
  `Teleosts` = getMRCA(tr, c("ENSSFOG00015000970 - Asian bonytongue",
                             "ENSDARG00000068263 - Zebrafish")),
  `Ray-finned fish` = getMRCA(tr, c("ENSDARG00000068263 - Zebrafish",
                                    "ENSECRG00000003618 - Reedfish"))
)

# plot ========================================================================
p_tr <- ggtree(tr, ladderize = TRUE) +
  geom_tiplab(align = TRUE, geom = "text", size = 0, offset = 0, linesize = 0.25) +
  geom_nodepoint(colour = tr.nodecols, size = 3, shape = 19) +
  geom_hilight(node = nodes_to_label[["Teleosts"]], fill = "#fdc086", alpha = 0.2) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Teleosts"]]), 
                   label = "Teleosts"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -0.05, nudge_y = 0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Ray-finned fish"]]), 
                   label = "Ray-finned fish"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -0.05, nudge_y = 0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Mammals"]]), 
                   label = "Mammals"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -0.05, nudge_y = 0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Reptiles"]]), 
                   label = "Reptiles"), 
               color = "grey40",
               hjust = 1, vjust = 1, 
               nudge_x = -0.05, nudge_y = -0.5, 
               size = 3) 

# combine the plots ===========================================================
pp <- p_l2 %>% 
  insert_left(p_tr, width = 2)

fs::dir_create("results/l2-length")
ggsave("results/l2-length/l2-length-on-gene-tree.pdf", pp,
       width = 12, height = 10)

# end =========================================================================

