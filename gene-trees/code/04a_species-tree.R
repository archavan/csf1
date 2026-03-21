# Arun Chavan
# Started: 2023-07-29

# background ==================================================================

# Generate a list of species to include in the species tree. Upload this files
# to timetree.org to get a newick file.

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)

basedir <- here("gene-trees")
datadir <- here(basedir, "data")
resdir  <- here(basedir, "results")

# data ========================================================================
spp <- read_csv(here(datadir, "ensembl", "ensembl-species.csv")) %>%
  janitor::clean_names()

# keep species ================================================================
keep_species <- c("Red-bellied piranha", 
                  "Electric eel",
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
                  "Orange clownfish",
                  "Turquoise killifish",
                  "Spotted gar",
                  "Reedfish",
                  "Chicken",
                  "Zebra finch",
                  "Painted turtle",
                  "Green anole",
                  "Golden Hamster",
                  "Mouse",
                  "Rat",
                  "Kangaroo rat",
                  "Squirrel",
                  "Guinea Pig",
                  "Pika",
                  "Macaque",
                  "Human",
                  "Gibbon",
                  "Bonobo",
                  "Tarsier",
                  "Rabbit",
                  "Blue whale",
                  "Goat",
                  "Cow",
                  "Pig",
                  "American mink",
                  "Dog",
                  "Cat",
                  "Horse",
                  "Microbat",
                  "Armadillo",
                  "Elephant",
                  "Platypus",
                  # the following don't have CSF1 but bracket the different clades
                  "Elephant shark",
                  "Coelacanth",
                  "Tropical clawed frog",
                  "Hagfish",
                  "Lamprey"
                  )

spp_filtered <- filter(spp, common_name %in% keep_species)

spp_filtered %>% 
  select(scientific_name) %>% 
  write_tsv(here(resdir, "l2-length", "species-for-tree.tsv"), col_names = FALSE)

# read tree downloaded from timetree.org ======================================
tr <- read.tree(here(datadir, "species-tree_timetree-org.nwk"))

# make sure tip labels match ==================================================
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")

spp_filtered <- spp_filtered %>% 
  filter(scientific_name %in% tr$tip.label)

spp_names <- spp_filtered %>% 
  select(scientific_name, common_name) %>% 
  deframe()

tr <- keep.tip(tr, names(spp_names))

# Use common names ============================================================
tr$tip.label <- spp_names[tr$tip.label]

# write files =================================================================
spp_filtered %>% 
  write_tsv(here(resdir, "l2-length", "species-filtered.tsv"))

write.tree(tr, here(resdir, "l2-length", "species-tree.nwk"))
# end =========================================================================
