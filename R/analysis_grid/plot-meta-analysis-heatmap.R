# Demonstrate the importance of meta-analysis by plotting the mean rank of 
# each coefficient across all 206 datasets.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# complexes
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # restore 'missing=NA'
  replace_na(list(missing = 'NA')) %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))

# GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))
# for each combination, filter by breadth, and get median AUROC across GO terms
dat2 %<>%
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# get mean per dataset
df1 = dat1 %>%
  group_by(analysis, accession, experiment, metric) %>%
  summarise(mean = mean(auroc)) %>%
  ungroup() %>%
  group_by(analysis, accession, experiment) %>%
  mutate(rank = rank(mean),
         rank_pct = (rank / max(rank)),
         rank_pct = rescale(rank_pct, c(0, 1))) %>%
  ungroup()
df2 = dat2 %>%
  group_by(analysis, accession, experiment, metric) %>%
  summarise(mean = mean(auroc)) %>%
  ungroup() %>%
  group_by(analysis, accession, experiment) %>%
  mutate(rank = rank(mean),
         rank_pct = (rank / max(rank)),
         rank_pct = rescale(rank_pct, c(0, 1))) %>%
  ungroup()
df = bind_rows(df1, df2)

# plot
p1 = df %>%
  mutate(analysis = fct_recode(analysis, 'Protein complexes' = 'complexes')) %>%
  mutate(metric = clean_metric(metric)) %>%
  unite(xval, accession, experiment) %>%
  ggplot(aes(x = xval, y = reorder(metric, mean, base::mean), 
             fill = rank_pct)) +
  facet_grid(~ analysis, space = 'free', scales = 'free') +
  # geom_tile(color = 'white', size = 0.1) +
  geom_tile(color = 'white', size = 0.08) +
  # scale_fill_paletteer_c("pals::ocean.deep", name = 'Rank (%) ',
  #                        labels = function(x) x * 100, 
  #                        breaks = seq(0, 1, 0.5)) +
  scale_fill_scico(palette = "davos", name = 'Rank (%) ',
                   labels = function(x) x * 100, 
                   breaks = seq(0, 1, 0.5)) +
  scale_x_discrete('Experiment', expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # guides(fill = guide_colorbar(nbin = 9, ticks = F, raster = F)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  boxed_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.3, 'lines'),
        legend.key.width = unit(0.25, 'lines'),
        legend.position = 'right',
        legend.justification = 'bottom')
p1
ggsave("fig/analysis/analysis_grid/AUCs-per-experiment-rank.pdf", p1,
       width = 18.5, height = 6.5, units = 'cm', useDingbats = F)
