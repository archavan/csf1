# Arun Chavan
# Started: 2023-07-28

# background ==================================================================

# Plot L2 lengths against the species tree. 

# setup =======================================================================
library(tidyverse)
library(ape)
library(ggtree)
library(aplot)

# data ========================================================================
ids <- read_tsv("results/isoforms/l2-lengths.tsv")
spp <- read_tsv("results/l2-length/species-filtered.tsv")
tr <- read.tree("results/l2-length/species-tree.nwk")
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")

# match species names across datasets =========================================
ids <- ids %>% 
  mutate(species = str_replace_all(species, "_", " ")) %>% 
  mutate(species = str_to_sentence(species))

spp <- select(spp, common_name, scientific_name)

ids <- left_join(spp, ids, by = c("scientific_name" = "species"))

# add paralogy info ===========================================================
ids <- ids %>% 
  group_by(scientific_name) %>% 
  mutate(paralog = as.integer(as.factor(gid))) %>% 
  mutate(paralog = as.character(paralog)) %>% 
  ungroup()

# plot length data ============================================================
p_l2 <- ggplot(ids, aes(l2_len, common_name, fill = paralog)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, size = 3) +
  scale_y_discrete(position = "right") +
  scale_fill_discrete(na.translate = FALSE) +
  labs(x = "amino acids", 
       y = NULL,
       title = "Length of Linker2 domain in all isoforms") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.5),
        axis.text.y.right = element_text(size = 9),
        plot.title = element_text(size = 9))

# plot tree ===================================================================
nodes_to_label <- c(
  `Mammals` = getMRCA(tr, c("Platypus", "Elephant")),
  `Reptiles` = getMRCA(tr, c("Painted turtle", "Green anole")),
  `Teleosts` = getMRCA(tr, c("Asian bonytongue", "Zebrafish")),
  `Ray-finned fish` = getMRCA(tr, c("Zebrafish", "Reedfish")),
  `Vertebrates` = getMRCA(tr, c("Human", "Hagfish")),
  `Lobe-finned fish` = getMRCA(tr, c("Human", "Coelacanth"))
)

p_tr <- ggtree(tr, ladderize = TRUE) +
  geom_tiplab(align = TRUE, geom = "text", size = 0, offset = 0, linesize = 0.25) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Teleosts"]]), 
                   label = "Teleosts"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -7, nudge_y = 0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Ray-finned fish"]]), 
                   label = "Ray-finned\nfish"), 
               color = "grey40",
               hjust = 1, vjust = 1, 
               nudge_x = -7, nudge_y = -0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Mammals"]]), 
                   label = "Mammals"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -7, nudge_y = 0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Reptiles"]]), 
                   label = "Reptiles"), 
               color = "grey40",
               hjust = 1, vjust = 1, 
               nudge_x = -7, nudge_y = -0.5, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Vertebrates"]]), 
                   label = "Vertebrates"), 
               color = "grey40", angle = 270,
               hjust = 0.5, vjust = 0, 
               nudge_x = 1, nudge_y = 0, 
               size = 3) +
  geom_nodelab(aes(subset = (node == nodes_to_label[["Lobe-finned fish"]]), 
                   label = "Lobe-finned fish"), 
               color = "grey40",
               hjust = 1, vjust = 0, 
               nudge_x = -7, nudge_y = 0.5, 
               size = 3) 

# combine plots ===============================================================
pp <- p_l2 %>% 
  insert_left(p_tr, width = 2)

ggsave("results/l2-length/l2-length-on-species-tree.pdf", pp,
       width = 13, height = 10)

# end =========================================================================
