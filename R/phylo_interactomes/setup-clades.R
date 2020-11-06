# Write a table denoting the species within each clade of the phylogenetic tree.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ape)

# read phylogenetic tree
tree = read.tree("data/resources/TimeTree/species.nwk")

# first, get all 'tips' (individual species)
tips = tree$tip.label

# get all possible 'internal' nodes (subtrees)
subtrees = subtrees(tree)
subtree_names = map_chr(subtrees, ~ .$node.label[1])
stopifnot(all(subtree_names == tree$node.label))
subtree_species = map(subtrees, 'tip.label') %>%
  setNames(subtree_names)

# add tips to the list of subtrees
tip_species = as.list(tips) %>% setNames(tips)
subtree_species %<>% c(tip_species)

# recode the species to match those in CF-MS-searches
subtree_species %<>% map(~ chartr('_', ' ', .) %>%
                           gsub("oleracea", "oleracea var. italica", .))
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")
setdiff(unlist(subtree_species), expts$Species)

# set up a table
grid = data.frame(node_id = names(subtree_species),
                  species = map_chr(subtree_species, paste0,
                                    collapse = ';')) %>%
  # rename empty string ('') to root
  mutate(node_id = ifelse(node_id == '', 'root', node_id))

# write
write.csv(grid, "data/analysis/phylo_interactomes/phylo_grid.csv",
          row.names = FALSE)
