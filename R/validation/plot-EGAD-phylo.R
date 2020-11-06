setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
library(ggtree)
library(phylobase)
library(treeio)
source("R/theme.R")

# read EGAD data
egad = readRDS("data/analysis/validation/EGAD-phylo.rds") %>% 
  # recode nodes
  mutate(node_clean = fct_recode(node, 
                                 'Euarchontoglires' = '14',
                                 'Tetrapoda' = '13',
                                 'Deuterostomia' = '11',
                                 'Opiskothonta' = '22',
                                 'BOP clade' = '6',
                                 'Mesangiospermae' = '29',
                                 'Viridiplantae' = '27',
                                 'Plasmodium' = '40',
                                 'Eukaryota' = '',
                                 'Eukaryota' = 'root') %>%
           as.character()) %>%
  filter(is.na(as.integer(node_clean)))

# prepare to plot
dat = egad %>%
  # breadth filter
  filter(between(n_proteins, 10, 100)) %>%
  # tag y-value of each
  group_by(node) %>%
  arrange(auroc) %>%
  mutate(y = row_number(),
         y = y / max(y)) %>%
  ungroup() %>%
  # clean up network names
  mutate(node_clean = chartr('_', ' ', node_clean),
         network = node_clean)

# set up function
plot_spectrum = function(dat, range = NULL) {
  dat$network[grepl(" ", dat$network)] = 
    with(dat, strsplit(network, ' ') %>%
           map_chr(~ paste0(substr(.[1], 1, 1), '. ', .[2])))[
             grepl(" ", dat$network)]
  
  spectrum = dat %>%
    # tag y-value of each
    group_by(network) %>%
    arrange(auroc) %>%
    mutate(y = row_number(),
           y = y / max(y),
           y = rescale(y, c(0, 1))) %>%
    ungroup()
  spectrum_neg = spectrum %>%
    group_by(network) %>%
    summarise(y = mean(auroc < 0.5), median = median(auroc),
              mean = mean(auroc)) %>%
    arrange(y)
  
  max_dist = 3
  digits = 3
  lvls = spectrum_neg$network
  
  if (is.null(range))
    range = range(spectrum$auroc) %>% round(2)
  
  p1a = spectrum %>%
    mutate(network = factor(network, levels = lvls)) %>%
    mutate(auroc = winsorize(auroc, range)) %>%
    ggplot(aes(x = network, y = y)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_raster(aes(fill = auroc)) +
    geom_errorbar(data = spectrum_neg %>%
                    mutate(network = factor(network, levels = lvls)),
                  aes(ymin = y, ymax = y), width = 1, size = 0.3,
                  color = 'grey20') +
    scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),
                         name = 'AUC',
                         limits = range,
                         breaks = range,
                         labels = range) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous('% of GO terms', 
                       expand = c(0, 0), 
                       limits = c(-0.001, 1),
                       labels = function(x) paste0(x * 100, '%'),
                       breaks = seq(0, 1, 0.25)
    ) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                 title.pos = 'top', title.hjust = 0.5)) +
    coord_flip() +
    boxed_theme() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          plot.margin = margin(rep(0, 4)),
          legend.key.width = unit(0.25, 'lines'),
          legend.key.height = unit(0.25, 'lines'),
          legend.position = 'top')
  p1a
  
  # median as the second column
  p1b = spectrum_neg %>%
    mutate(network = factor(network, levels = lvls)) %>%
    arrange(network) %>%
    ggplot(aes(y = network, x = 0)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_text(aes(label = format(median, digits = 2)), 
              size = 2, hjust = 0, color = 'grey30') +
    scale_y_discrete(labels = function(x) gsub("^.*\\|", "", x),
                     expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) +
    boxed_theme() +
    theme(strip.text = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(rep(0, 4)),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'none')
  p1b
  
  p = p1a + p1b + plot_layout(widths = c(1, 0.2))
  return(p)
}

# spectrum
p1 = plot_spectrum(dat, range = c(0.45, 0.9))
p1
ggsave("fig/analysis/validation/EGAD-phylo.pdf", p1,
       width = 5.8, height = 10, units = "cm", useDingbats = FALSE)

# now, overlay median AUC onto phylogenetic tree
tree = read.tree("data/resources/TimeTree/species.nwk")
grid = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv")
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
median_auc = dat %>%
  group_by(node, node_clean) %>%
  summarise(auc = median(auroc), pct_neg = mean(auroc < 0.5)) %>%
  ungroup() %>%
  mutate(label = replace(node, node == 'root', '')) %>%
  dplyr::select(-node) %>%
  ## abbreviate species
  mutate(node_clean2 = strsplit(node_clean, ' ') %>%
           map_chr(~ paste0(substr(.[1], 1, 1), '. ', .[2])))
tbl = as_tibble(tree) %>%
  left_join(median_auc, by = 'label') %>%
  left_join(internal_nodes, by = 'label') %>%
  # plot only interactions for tips or labelled internal nodes
  mutate(auc = ifelse(is.na(node_clean), NA, auc),
         pct_neg = ifelse(is.na(node_clean), NA, pct_neg),
         size = ifelse(is.na(node_clean), NA, 1))
tree2 = as.treedata(tbl)
pal = pals::kovesi.isoluminant_cm_70_c39(100)
p2 = ggtree(tree2, size = 0.3, layout = 'rectangular') +
  geom_rootedge(rootedge = 100, size = 0.3) +
  geom_point(aes(fill = auc, size = size), shape = 21, color = 'grey20') +
  geom_tiplab(size = 2, offset = 100, aes(label = node_clean2)) +
  geom_nodelab(aes(label = node_clean), size = 2) +
  scale_size_continuous(range = c(0.1, 2.4), guide = FALSE) +
  scale_fill_gradientn(colours = pal, 'AUC (median)',
                       breaks = c(0.525, 0.616),
                       limits = c(0.525, 0.616),
                       labels = c(0.53, 0.62)) +
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
ggsave("fig/analysis/validation/EGAD-phylo-tree.pdf", p2,
       width = 6.2, height = 6.5, units = 'cm', useDingbats = FALSE)
