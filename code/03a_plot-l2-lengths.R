# Arun Chavan
# Started: 2023-07-27

# background ==================================================================

# Here we plot L2 lenghts against the phylogeny.

# setup =======================================================================
library(tidyverse)
library(ape)
library(ggtree)
library(aplot)

# data ========================================================================
tr1 <- read.tree("results/csf1-tree-1.nwk")
l2 <- read_tsv("results/linker2/l2-lengths.tsv")
tips <- read_tsv("results/csf1-tree-1_tips.tsv")

# filter data =================================================================

# For not all proteins that are tips in this tree we could calculate the L2
# length. We need to drop them.

tr1$tip.label <- str_replace_all(tr1$tip.label, "_", " ") 
tr2 <- keep.tip(tr1, l2$newlab)


# find duplication nodes ======================================================
dup_nodes <- c()
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr2, tr2$tip.label[grepl("Zebrafish", tr2$tip.label)]))
dup_nodes <- append(dup_nodes, 
                    getMRCA(tr2, c("ENSSSAP00000150975 - Atlantic salmon",
                                   "ENSSSAP00000108598 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr2, c("ENSSSAP00000154249 - Atlantic salmon",
                                   "ENSSSAP00000149216 - Atlantic salmon")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr2, c("ENSCARP00000058590 - Goldfish",
                                   "ENSCARP00000062662 - Goldfish")))
dup_nodes <- append(dup_nodes,
                    getMRCA(tr2, c("ENSCARP00000091986 - Goldfish",
                                   "ENSCARP00000075789 - Goldfish")))

tr2.nodecols <- rep('NA', tr2$Nnode)
tr2.nodecols[dup_nodes - length(tr2$tip.label)] <- 'red' 

# nodes to label ==============================================================
nodes_to_label <- c(`Mammals` = getMRCA(tr2, c("ENSOANP00000011057 - Platypus",
                                               "ENSLAFP00000026014 - Elephant")),
                    `Reptiles` = getMRCA(tr2, c("ENSGALP00010041133 - Chicken",
                                                "ENSACAP00000011757 - Green anole")),
                    `Teleosts` = getMRCA(tr2, c("ENSSFOP00015001355 - Asian bonytongue",
                                                "ENSDARP00000089407 - Zebrafish")),
                    `Ray-finned fish` = getMRCA(tr2, c("ENSDARP00000089407 - Zebrafish",
                                                       "ENSECRP00000005354 - Reedfish")))

# plot ========================================================================
p_tr <- ggtree(tr2, ladderize = TRUE) +
  geom_tiplab(align = TRUE, geom = "text", size = 0, offset = 0, linesize = 0.25) +
  geom_nodepoint(colour = tr2.nodecols, size = 3, shape = 19) +
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

p_l2 <- ggplot(l2, aes(l2_length, newlab)) +
  geom_point(shape = 19, size = 2) +
  scale_y_discrete(position = "right") +
  labs(x = "amino acids", 
       y = NULL,
       title = "Length of Linker2 domain in the longest isoform") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.5),
        axis.text.y.right = element_text(size = 9),
        plot.title = element_text(size = 9))

pp <- p_l2 %>% 
  insert_left(p_tr, width = 2)

ggsave("results/linker2/l2-length-on-gene-tree.pdf", pp, width = 12, height = 10)

# write files =================================================================
write.tree(tr2, "results/linker2/csf-tree-2.nwk")

# end =========================================================================
