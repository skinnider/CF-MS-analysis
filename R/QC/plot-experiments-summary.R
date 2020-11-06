# Plot a summary figure of all experiments analyzed in this study.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

# list all summary files
files = list.files('~/projects/rrg-ljfoster-ab/CF-MS/maxquant',
                   pattern = 'summary.txt', full.names = T, recursive = T)

# count number of input files and fractions per experiment
summary = map_dfr(files, ~ {
  tab = fread(.)
  n_files = n_distinct(tab$`Raw file`)
  n_fractions = n_distinct(tab$Experiment)
  data.frame(file = ., raw_files = n_files, fractions = n_fractions)
})

# flag accession, replicate, and search
summary %<>%
  mutate(file = gsub("^.*maxquant/", "", file)) %>%
  separate(file, c("accession", "replicate", "search", "x"), '/') %>%
  dplyr::select(-x)

# merge in species
experiments = read.csv('~/git/CF-MS-searches/data/experiments.csv') %>%
  # drop channels
  mutate(Replicate = gsub("-heavy|-medium", "", Replicate))
colnames(experiments) %<>% tolower()
summary %<>% left_join(experiments, ., by = c('accession', 'replicate'))

# create data frame for plotting
df = summary %>%
  # recode one-off species
  mutate(species = fct_recode(species,
                              'Other\nplant' = 'Solanum lycopersicum',
                              'Other\nplant' = 'Selaginella moellendorffii',
                              'Other\nplant' = 'Glycine max',
                              'Other\nplant' = 'Brassica oleracea var. italica',
                              'Other\nplant' = 'Oryza sativa',
                              'Other\nplant' = 'Zea mays',
                              'Other\neukaryote' = 'Saccharomyces cerevisiae',
                              'Other\neukaryote' = 'Chaetomium thermophilum',
                              'Prokaryote' = 'Cyanothece ATCC 51142',
                              'Plasmodium\nspp.' = 'Plasmodium berghei',
                              'Plasmodium\nspp.' = 'Plasmodium falciparum',
                              'Plasmodium\nspp.' = 'Plasmodium knowlesi',
                              'H. sapiens' = 'Homo sapiens',
                              'A. thaliana' = 'Arabidopsis thaliana',
                              'T. brucei' = 'Trypanosoma brucei',
                              'M. musculus' = 'Mus musculus',
                              'C. elegans' = 'Caenorhabditis elegans',
                              'X. laevis' = 'Xenopus laevis',
                              'T. aestivum' = 'Triticum aestivum',
                              'Other\neukaryote' = 'Strongylocentrotus purpuratus',
                              'Other\neukaryote' = 'Dictyostelium discoideum',
                              'Other\nplant' = 'Chlamydomonas reinhardtii',
                              'Other\neukaryote' = 'Drosophila melanogaster',
                              'Other\neukaryote' = 'Nematostella vectensis'
                              )) %>%
  mutate(species = fct_recode(species, 
                              'Other\neukaryote' = 'Plasmodium\nspp.')) %>%
  group_by(species, accession) %>%
  arrange(fractions) %>%
  mutate(yval = factor(as.character(row_number()),
                       levels = as.character(seq_len(99))))
  
# plot
res = df %>%
  # arrange facets by species
  group_by(species) %>%
  mutate(n_expts = n(),
         n_xvals = n_distinct(accession)) %>%
  ungroup() %>%
  arrange(desc(n_xvals), desc(n_expts)) %>%
  mutate(species = factor(species, levels = unique(species))) %>%
  # let's just do it manually to get it right
  mutate(species = fct_relevel(species, 'H. sapiens', 'A. thaliana',
                               'M. musculus', 'X. laevis', 'T. aestivum', 
                               'C. elegans', 'T. brucei',
                               'Other\neukaryote', 'Other\nplant',
                               'Prokaryote')) %>%
  # arrange accessions within facets by experiments
  group_by(species, accession) %>%
  mutate(n_yvals = n()) %>%
  ungroup() %>%
  arrange(desc(n_yvals)) %>%
  mutate(accession = factor(accession, levels = unique(.$accession))) %>%
  # bin fractions
  mutate(bin = cut(fractions, breaks = c(0, 20, 30, 40, 50, 75, 100, 999)),
         bin = fct_recode(bin, 
                          '<20' = '(0,20]',
                          '<25' = '(0,25]',
                          '21-30' = '(20,30]',
                          '31-40' = '(30,40]',
                          '41-50' = '(40,50]',
                          '26-50' = '(25,50]',
                          '51-75' = '(50,75]',
                          '76-100' = '(75,100]',
                          '>100' = '(100,999]'))
labels = res %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(label = paste0('n = ', n))
pal = viridis::cividis(n = 20) %>% extract(-20)
p = res %>%
  mutate(fractions = winsorize(fractions, c(NA, 100))) %>%
  ggplot(aes(x = accession, y = yval)) + 
  facet_grid(~ species, scales = 'free_x', space = 'free',
             labeller = as_labeller(function(x) gsub("\n", " ", x))) +
  geom_tile(color = 'white', aes(fill = fractions)) + 
  geom_label(data = labels, aes(x = +Inf, y = +Inf, label = label),
             hjust = 1, vjust = 1, label.size = NA, size = 1.75,
             label.padding = unit(0.3, 'lines')) +
  scale_fill_gradientn('Fractions ', colors = pal, na.value = 'grey90',
                       breaks = c(20, 60, 100), labels = c(20, 60, '>100')) +
  scale_x_reordered() +
  scale_y_discrete('Experiments', breaks = as.character(seq(0, 18, 2))) +
  guides(fill = guide_colorbar(raster = F, nbin = 9, ticks = F)) +
  boxed_theme(size_sm = 5, size_lg = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.5, 'lines'),
        legend.key.height = unit(0.3, 'lines'),
        legend.key.width = unit(0.45, 'lines')
  )
p
ggsave("fig/QC/dataset-overview.pdf", p, width = 18, height = 5.3, units = "cm",
       useDingbats = F)
