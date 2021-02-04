# Analyze the evolutionary conservation of interactions
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggtree)
library(phylobase)
library(treeio)
library(broom)
library(Matrix)
source("R/theme.R")

# list interactome files
network_files = list.files("~/Foster lab/CF-MS meta-analysis/Zenodo/Interactomes/50% precision",
                           pattern = "*.tsv", full.names = TRUE)
networks = map(network_files, read.delim) %>%
  setNames(gsub("\\..*$", "", basename(network_files)))

# calculate Jaccard index between all pairs
pairs = tidyr::crossing(network1 = names(networks), 
                        network2 = names(networks))
jaccard = map2_dfr(pairs$network1, pairs$network2, ~ {
  network1 = networks[[.x]]
  network2 = networks[[.y]]
  intersection = inner_join(network1, network2,
                            by = c('protein_A', 'protein_B')) %>% nrow() 
  union = full_join(network1, network2, by = c('protein_A', 'protein_B')) %>%
    nrow()
  sum = nrow(network1) + nrow(network2)
  data.frame(network1 = .x, network2 = .y, 
             intersection = intersection,
             union = union,
             jaccard = intersection / union,
             dice = 2 * intersection / sum,
             percentage1 = intersection / nrow(network1),
             percentage2 = intersection / nrow(network2))
})

# plot human, mouse, and Arabidopsis conservation up the tree
ancestors = c('Euarchontoglires', 
              'Tetrapoda',
              'Deuterostomia',
              'Opiskothonta',
              'Eukaryota')
filter(jaccard, network1 == 'Homo_sapiens' & network2 %in% ancestors) %>%
  arrange(jaccard)
filter(jaccard, network1 == 'Mus_musculus' & network2 %in% ancestors) %>%
  arrange(jaccard)
filter(jaccard, network1 == 'Arabidopsis_thaliana' & 
         network2 %in% c('Mesangiospermae', 'Viridiplantae', 'Eukaryota')) %>%
  arrange(jaccard)

# read phylogenetic tree
tree = read.tree("data/resources/TimeTree/species.nwk")

# read clade grid
grid = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv") %>%
  dplyr::rename(label = node_id)

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
  data.frame(label = ., clean_name = names(.))

# merge in % conserved interactions
## human
jaccard_hs = jaccard %>%
  filter(network1 == 'Homo_sapiens' & network2 %in% ancestors) %>%
  dplyr::rename(clean_name = network2)
res_hs = grid %>%
  left_join(internal_nodes, by = 'label') %>%
  left_join(jaccard_hs, by = 'clean_name') %>%
  drop_na(percentage1) %>%
  mutate(label = replace(label, label == 'root', ''))
## mouse
jaccard_mm = jaccard %>%
  filter(network1 == 'Mus_musculus' & network2 %in% ancestors) %>%
  dplyr::rename(clean_name = network2)
res_mm = grid %>%
  left_join(internal_nodes, by = 'label') %>%
  left_join(jaccard_mm, by = 'clean_name') %>%
  drop_na(percentage1) %>%
  mutate(label = replace(label, label == 'root', ''))
## arabidopsis
jaccard_at = jaccard %>%
  filter(network1 == 'Arabidopsis_thaliana' & 
           network2 %in% c('Mesangiospermae', 'Viridiplantae', 'Eukaryota')) %>%
  dplyr::rename(clean_name = network2)
res_at = grid %>%
  left_join(internal_nodes, by = 'label') %>%
  left_join(jaccard_at, by = 'clean_name') %>%
  drop_na(percentage1) %>%
  mutate(label = replace(label, label == 'root', ''))

# plot human
tbl_hs = as_tibble(tree) %>%
  left_join(res_hs, by = 'label') %>%
  mutate(text = ifelse(!is.na(percentage1), 
                       paste0(clean_name, ':\n', 
                              format(100 * percentage1, format = 'f', 
                                     digits = 3), '%'),
                       ifelse(label == 'Homo_sapiens', 'Homo sapiens', NA)),
         size = ifelse(is.na(text) & label != 'Homo_sapiens', NA, 1))
tree_hs = as.treedata(tbl_hs)
pal = pals::kovesi.diverging_isoluminant_cjo_70_c25(100)
range = range(tbl_hs$percentage1, na.rm = TRUE)
rounded = round(100 * range, digits = 1) %>% format(digits = 3)
p1 = ggtree(tree_hs, size = 0.3, layout = 'rectangular') +
  geom_rootedge(rootedge = 100, size = 0.3) +
  geom_point(aes(fill = percentage1, size = size), shape = 21, color = 'grey20') +
  geom_tiplab(size = 2, offset = 100, aes(label = text)) +
  geom_nodelab(aes(label = text), size = 2, lineheight = 1) +
  scale_size_continuous(range = c(0.1, 2.4), guide = FALSE) +
  scale_fill_gradientn(colours = pal, 
                       name = '% interactions conserved',
                       breaks = range,
                       limits = range,
                       labels = rounded) +
  scale_x_continuous(limits = c(-100, 2800)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.key.height = unit(0.25, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p1

# plot mouse
tbl_mm = as_tibble(tree) %>%
  left_join(res_mm, by = 'label') %>%
  mutate(text = ifelse(!is.na(percentage1), 
                       paste0(clean_name, ':\n', 
                              format(100 * percentage1, format = 'f', 
                                     digits = 3), '%'),
                       ifelse(label == 'Mus_musculus', 'Mus musculus', NA)),
         size = ifelse(is.na(text) & label != 'Mus_musculus', NA, 1))
tree_mm = as.treedata(tbl_mm)
range = range(tbl_mm$percentage1, na.rm = TRUE)
rounded = round(100 * range, digits = 1) %>% format(digits = 3)
p2 = ggtree(tree_mm, size = 0.3, layout = 'rectangular') +
  geom_rootedge(rootedge = 100, size = 0.3) +
  geom_point(aes(fill = percentage1, size = size), shape = 21, color = 'grey20') +
  geom_tiplab(size = 2, offset = 100, aes(label = text)) +
  geom_nodelab(aes(label = text), size = 2, lineheight = 1) +
  scale_size_continuous(range = c(0.1, 2.4), guide = FALSE) +
  scale_fill_gradientn(colours = pal, 
                       name = '% interactions conserved',
                       breaks = range,
                       limits = range,
                       labels = rounded) +
  scale_x_continuous(limits = c(-100, 2800)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.key.height = unit(0.25, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p2

# plot A. thaliana
tbl_at = as_tibble(tree) %>%
  left_join(res_at, by = 'label') %>%
  mutate(text = ifelse(!is.na(percentage1), 
                       paste0(clean_name, ':\n', 
                              format(100 * percentage1, format = 'f', 
                                     digits = 3), '%'),
                       ifelse(label == 'Arabidopsis_thaliana', 'Arabidopsis thaliana', NA)),
         size = ifelse(is.na(text) & label != 'Arabidopsis_thaliana', NA, 1))
tree_at = as.treedata(tbl_at)
range = range(tbl_at$percentage1, na.rm = TRUE)
rounded = round(100 * range, digits = 1) %>% format(digits = 3)
p3 = ggtree(tree_at, size = 0.3, layout = 'rectangular') +
  geom_rootedge(rootedge = 100, size = 0.3) +
  geom_point(aes(fill = percentage1, size = size), shape = 21, color = 'grey20') +
  geom_tiplab(size = 2, offset = 100, aes(label = text)) +
  geom_nodelab(aes(label = text), size = 2, lineheight = 1) +
  scale_size_continuous(range = c(0.1, 2.4), guide = FALSE) +
  scale_fill_gradientn(colours = pal, 
                       name = '% interactions conserved',
                       breaks = range,
                       limits = range,
                       labels = rounded) +
  scale_x_continuous(limits = c(-100, 2800)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.key.height = unit(0.25, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p3

# combine
p = p1 | p2 | p3
ggsave("fig/analysis/phylo_interactomes/parent-conservation.pdf", 
       width = 14, height = 5, units = "cm", useDingbats = FALSE)
