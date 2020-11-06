setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
source("R/theme.R")

# read EGAD data
egad = readRDS("data/analysis/validation/EGAD.rds")

# manually recode networks
dat = egad %>%
  mutate(network = gsub("^.*\\|", "", network)) %>%
  # breadth filter
  filter(between(n_proteins, 10, 100)) %>%
  # tag y-value of each
  group_by(network) %>%
  arrange(auroc) %>%
  mutate(y = row_number(),
         y = y / max(y)) %>%
  ungroup() %>%
  # clean up network names
  mutate(network = fct_recode(network,
                              'HI-II-14' = 'Rolland2014',
                              'HuRI' = 'Luck2019',
                              'BioPlex' = 'Huttlin2015',
                              'BioPlex 2' = 'Huttlin2017',
                              'QUBIC' = 'Hein2015') %>% as.character()) %>%
  filter(!network %in% c('Havugimana2012', 'Wan2015'))

# set up for plotting
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

range = c(0.45, 0.9)
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
ggsave("fig/analysis/validation/EGAD.pdf", p, width = 4.75, height = 4.75, 
       units = "cm", useDingbats = FALSE)
