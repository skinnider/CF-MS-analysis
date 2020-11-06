setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read data
complexes = readRDS("data/analysis/fractions_vs_replicates/complexes.rds") %>%
  filter(quant_mode == 'iBAQ')

# calculate mean over samples within each dataset
complex_means = complexes %>%
  # calculate mean over samples
  group_by_at(vars(-sample_idx, -auroc, -n_chroms, -n_pairs)) %>%
  summarise(auroc = mean(auroc), n = n()) %>%
  ungroup()

# calculate the mean over multiple datasets
complex_mmeans = complex_means %>% 
  filter(n == 10) %>%
  filter(n_fractions <= 50) %>%
  group_by_at(vars(-n_fractions, -n_replicates, -auroc, -n)) %>%
  mutate(auroc = rescale(auroc, c(0, 1))) %>%
  ungroup() %>%
  group_by(metric, n_fractions, n_replicates) %>%
  summarise(mean = mean(auroc), sd = sd(auroc), sem = sd / sqrt(n()),
            n1 = sum(n), n2 = n()) %>%
  ungroup() %>%
  type_convert()

# plot overview
p1 = complex_mmeans %>%
  filter(n_replicates <= 3,
         n_fractions <= 50,
         metric == 'pearson') %>%
  ggplot(aes(y = factor(n_fractions), x = factor(n_replicates), fill = mean)) + 
  facet_grid(~ "Complexes") +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC,\nscaled',
                         breaks = c(0.2, 0.8), labels = c('min', 'max')) +
  scale_y_discrete('Fractions', expand = c(0, 0)) +
  scale_x_discrete('Replicates', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  coord_fixed() + 
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        strip.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p1
ggsave("fig/analysis/fractions_vs_replicates/complex-AUC-average.pdf", p1,
       width = 3.5, height = 4.25, units = 'cm', useDingbats = FALSE)
## also plot with legend flipped
p1v2 = complex_mmeans %>%
  filter(n_replicates <= 3,
         n_fractions <= 50,
         metric == 'pearson') %>%
  ggplot(aes(y = factor(n_fractions), x = factor(n_replicates), fill = mean)) + 
  facet_grid(~ "Complexes") +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC, scaled',
                         breaks = c(0.2, 0.8), labels = c('min', 'max')) +
  scale_y_discrete('Fractions', expand = c(0, 0)) +
  scale_x_discrete('Replicates', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  coord_fixed() + 
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.position = 'top',
        legend.justification = 'right',
        legend.key.width = unit(0.3, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        strip.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p1v2
ggsave("fig/analysis/fractions_vs_replicates/complex-AUC-average-legend.pdf", p1v2,
       width = 3.5, height = 4.25, units = 'cm', useDingbats = FALSE)

# also try as a line chart
pos = position_dodge(width = 2)
colors = brewer.pal(n = 4, name = "Blues") %>% tail(-1)
p1b = complex_mmeans %>%
  filter(n_replicates <= 4, 
         n_fractions <= 50, 
         metric == 'pearson') %>%
  ggplot(aes(x = n_fractions, y = mean, color = factor(n_replicates))) +
  facet_grid(~ "Complexes") +
  geom_point(size = 0.8, position = pos) + 
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0, position = pos) +
  geom_line(position = pos) +
  scale_color_manual(name = '# of replicates', values = colors) +
  scale_x_continuous('# of fractions') +
  scale_y_continuous('AUC, scaled', limits = c(0, NA)) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.direction = 'horizontal',
        strip.text = element_text(size = 7))
p1b
ggsave("fig/analysis/fractions_vs_replicates/complex-AUC-average-line.pdf",
       p1b, width = 6, height = 6, units = "cm", useDingbats = FALSE)

# plot each individual dataset
cdf = readRDS("data/QC/protein-groups-CDF.rds")
max_fractions = cdf %>%
  filter(quant_mode == 'iBAQ') %>%
  group_by(accession, experiment) %>%
  summarise(max_fractions = max(n_fractions)) %>%
  ungroup() 

grey = sequential_hcl('Greens 2', n = 100)[100]
p2 = complex_means %>%
  filter(metric == 'pearson', n == 10) %>%
  type_convert() %>%
  group_by(accession, pattern) %>%
  mutate(scaled = rescale(auroc, to = c(0, 1))) %>%
  ungroup() %>%
  # clean up facets
  mutate(facet = fct_recode(paste(accession, pattern),
                            'PXD001220' = 'PXD001220 all',
                            'PXD002892\nBN-PAGE\nstimulated' = 'PXD002892 BN_heavy',
                            'PXD002892\nBN-PAGE\ncontrol' = 'PXD002892 BN_medium',
                            'PXD002892\nSEC\nstimulated' = 'PXD002892 SEC_heavy',
                            'PXD002892\nSEC\ncontrol' = 'PXD002892 SEC_medium')) %>%
  # add missing values
  tidyr::complete(nesting(facet, n_replicates), n_fractions,
                  fill = list(scaled = NA)) %>%
  # filter(n_fractions <= 55) %>%
  ggplot(aes(y = factor(n_fractions), x = factor(n_replicates),
             fill = scaled)) + 
  facet_grid(~ facet, scales = 'free', space = 'free') +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC, scaled',
                         breaks = c(0, 1), labels = c('min', 'max'),
                         na.value = grey) +
  scale_y_discrete('Fractions', expand = c(0, 0)) +
  scale_x_discrete('Replicates', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  # coord_fixed() + 
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.3, 'lines'),
        aspect.ratio = length(seq(10, 100, 5)) / 5,
        panel.spacing = unit(0.8, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p2
ggsave("fig/analysis/fractions_vs_replicates/complex-AUC-individual-datasets.pdf", p2,
       width = 15, height = 7, units = 'cm', useDingbats = FALSE)

# read data
go = readRDS("data/analysis/fractions_vs_replicates/GO.rds") %>%
  filter(quant_mode == 'iBAQ')

# calculate mean over GO terms per dataset (all GO terms)
go_means = go %>%
  # filter GO terms by breadth
  filter(between(n_annotated, 10, 100)) %>%
  # calculate mean AUC over all GO terms
  group_by_at(vars(-go_term, -auroc,
                   -n_proteins, -n_annotated, -n_chroms, -n_pairs)) %>%
  summarise(auroc = mean(auroc), n_terms = n()) %>%
  ungroup() %>%
  # calculate mean of means, over ten samples
  group_by_at(vars(-sample_idx, -auroc, -n_terms)) %>%
  summarise(auroc = mean(auroc), n = n()) %>%
  ungroup() %>%
  # filter to complete datasets
  filter(n == 10)

# calculate the mean over multiple datasets
go_mmeans = go_means %>%
  # complete datasets only
  filter(n == 10) %>%
  # filter by fractions, too, before scaling
  filter(n_fractions <= 50) %>%
  group_by_at(vars(-n_fractions, -n_replicates, -auroc, -n)) %>%
  mutate(auroc = rescale(auroc, c(0, 1))) %>%
  ungroup() %>%
  group_by(metric, n_fractions, n_replicates) %>%
  summarise(mean = mean(auroc), sd = sd(auroc), sem = sd / sqrt(n()),
            n1 = sum(n), n2 = n()) %>%
  ungroup() %>%
  type_convert()

# plot the mean over all datasets
p3 = go_mmeans %>%
  filter(n_fractions <= 50, metric == 'pearson') %>%
  ggplot(aes(y = factor(n_fractions), x = factor(n_replicates), fill = mean)) + 
  facet_grid(~ "GO") +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC,\nscaled',
                         breaks = c(0.2, 0.8), 
                         labels = c("min", 'max')) +
  scale_y_discrete('Fractions', expand = c(0, 0)) +
  scale_x_discrete('Replicates', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  coord_fixed() + 
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p3
ggsave("fig/analysis/fractions_vs_replicates/GO-AUC-average.pdf", p3,
       width = 4.5, height = 4.25, units = 'cm', useDingbats = FALSE)

# what about as a line chart?
pos = position_dodge(width = 4)
colors = brewer.pal(n = 6, name = "Blues") %>% tail(-1)
p3b = go_mmeans %>%
  filter(n_fractions <= 50, metric == 'pearson') %>%
  ggplot(aes(x = n_fractions, y = mean, color = factor(n_replicates))) +
  facet_grid(~ 'GO')+
  geom_point(size = 0.8, position = pos) + 
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0, position = pos) +
  geom_line(position = pos) +
  scale_color_manual(name = '# of replicates', values = colors) +
  scale_x_continuous('# of fractions') +
  scale_y_continuous('AUC, scaled') +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.direction = 'horizontal',
        strip.text = element_text(size = 7))
p3b
ggsave("fig/analysis/fractions_vs_replicates/GO-AUC-average-line.pdf",
       p3b, width = 6, height = 6, units = "cm", useDingbats = FALSE)

# plot individual datasets
p4 = go_means %>%
  filter(n == 10) %>%
  type_convert() %>%
  # scale within each 
  group_by(metric, accession, pattern) %>%
  mutate(scaled = rescale(auroc, to = c(0, 1))) %>%
  ungroup() %>%
  # rename facets
  mutate(facet = fct_recode(paste(accession, pattern),
                            'PXD001220' = 'PXD001220 all',
                            'PXD002319\nBead IEX' = 'PXD002319 Ce_Bead',
                            'PXD002324\nIEX' = 'PXD002324',
                            'PXD002325\nBead IEX' = 'XPD002325',
                            'PXD002892\nBN-PAGE\nstimulated' = 'PXD002892 BN_heavy',
                            'PXD002892\nBN-PAGE\ncontrol' = 'PXD002892 BN_medium',
                            'PXD002892\nSEC\nstimulated' = 'PXD002892 SEC_heavy',
                            'PXD002892\nSEC\ncontrol' = 'PXD002892 SEC_medium',
                            'PXD003754' = 'PXD003754 all',
                            'PXD011304' = 'PXD011304 2D_IEF')) %>%
  # add missing values
  tidyr::complete(nesting(facet, n_replicates), n_fractions, 
                  fill = list(scaled = NA)) %>%
  # plot 
  ggplot(aes(y = factor(n_fractions), x = factor(n_replicates), fill = scaled)) + 
  facet_grid(~ facet, scales = 'free', space = 'free') +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC, scaled',
                         breaks = c(0, 1), labels = c('min', 'max'), 
                         na.value = grey) +
  scale_y_discrete('Fractions', expand = c(0, 0)) +
  scale_x_discrete('Replicates', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.35, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        aspect.ratio = length(seq(10, 100, 5)) / 5,
        panel.spacing = unit(0.8, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
p4
ggsave("fig/analysis/fractions_vs_replicates/GO-AUC-individual-datasets.pdf", p4,
       width = 16, height = 7, units = 'cm', useDingbats = FALSE)
