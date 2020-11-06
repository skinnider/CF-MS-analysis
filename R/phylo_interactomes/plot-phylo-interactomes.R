setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ape)
library(ggtree)
library(phylobase)
library(treeio)
library(broom)
source("R/theme.R")

# read phylogenetic tree
tree = read.tree("data/resources/TimeTree/species.nwk")

# read clade grid
grid = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv")

# set up labelled internal nodes
internal_nodes = c('Euarchontoglires' = '14',
                   'Tetrapoda' = '13',
                   'Deuterostomia' = '11',
                   'Opiskothonta' = '22',
                   'BOP clade' = '6',
                   'Mesangiospermae' = '29',
                   'Viridiplantae' = '27',
                   'Plasmodium' = '40',
                   'Eukaryota' = '',
                   'Eukaryota' = 'root') %>%
  data.frame(label = ., internal_label = names(.))

# read results
PR = readRDS("data/analysis/phylo_interactomes/PR.rds")
n_ppi = readRDS("data/analysis/phylo_interactomes/n_ppis.rds") %>%
  filter(precision %in% c(0.5, 0.75, 0.8, 0.9, 0.95))
res = grid %>%
  left_join(n_ppi, by = 'node_id') %>%
  dplyr::rename(label = node_id) %>%
  mutate(label = replace(label, label == 'root', ''))

# add to tree
meta = res_sum %>%
  spread(precision, n_ppis) %>%
  dplyr::rename(n_50 = `0.5`) %>%
  mutate(label_clean = chartr('_', ' ', label))
tbl = as_tibble(tree) %>%
  left_join(meta, by = 'label') %>%
  left_join(internal_nodes, by = 'label') %>%
  # plot only interactions for tips or labelled internal nodes
  mutate(label_clean = ifelse(!is.na(as.integer(label_clean)), NA, label_clean)) %>%
  mutate(n_50 = ifelse(!is.na(label_clean) | !is.na(internal_label), n_50, NA))
# count the total number of interactions  
sum(tbl$n_50, na.rm = TRUE)

# convert back to tree for ggtree
tree2 = as.treedata(tbl)

# plot
pal = brewer.pal(4, 'Blues') %>% extract(-1)
p = ggtree(tree2, size = 0.3, layout = 'rectangular') +
  geom_rootedge(rootedge = 50, size = 0.3) +
  #              color = 'grey50') +
  geom_point(aes(size = n_50), color = darken(pal[1], 1.1), fill = pal[1],
             alpha = 0.7) +
  geom_tiplab(size = 2, offset = 100, aes(label = label_clean)) +
  geom_nodelab(aes(label = internal_label), size = 2) +
  scale_size_continuous(name = expression(Interactions~(10^3)), 
                        range = c(0.1, 3), labels = function(x) x / 1e3) +
  scale_x_continuous(limits = c(-50, 2800)) +
  boxed_theme() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p
ggsave("fig/analysis/phylo_interactomes/tree.pdf", p,
       width = 8.9, height = 7.5, units = 'cm')

# print number in each species or clade
res %>% 
  left_join(internal_nodes) %>% 
  mutate(label = ifelse(!is.na(internal_label), internal_label, label)) %>%
  filter(precision == 0.5) %$%
  setNames(n_ppis, label) %>%
  sort()

# plot correlation between # of experiments and # of interactions
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  dplyr::rename(species = Species) %>%
  dplyr::select(species, Accession, Replicate)
xy = tbl %>%
  dplyr::select(label, n_50, label_clean, internal_label) %>%
  mutate(label_clean = coalesce(label_clean, internal_label)) %>%
  drop_na(label_clean) %>%
  dplyr::select(-label_clean, -internal_label) %>%
  # merge in species
  dplyr::rename(node_id = label) %>%
  left_join(grid %>%
              mutate(node_id = replace(node_id, node_id == 'root', ''))) %>%
  # merge in experiments
  mutate(species = strsplit(species, ';')) %>%
  unnest(species) %>%
  left_join(expts, by = 'species') %>%
  # count experiments per species
  dplyr::count(node_id, n_50)
col = colours.cafe425[1]
p2 = xy %>%
  ggplot(aes(x = n, y = n_50)) +
  geom_point(size = 0.8, shape = 1) +
  geom_smooth(alpha = 0.15, size = 0.4, color = col) +
  scale_y_continuous('Interactions (thousands)', label = function(x) x / 1e3) +
  scale_x_continuous('# of CF-MS experiments') +
  boxed_theme() +
  theme(aspect.ratio = 1)
p2
ggsave("fig/analysis/phylo_interactomes/saturation.pdf", p2,
       width = 5.5, height = 5.5, units = "cm", useDingbats = FALSE)
