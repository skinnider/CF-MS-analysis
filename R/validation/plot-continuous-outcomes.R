setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
source("R/theme.R")

# read data
outcomes = readRDS("data/analysis/validation/continuous_outcomes.rds")

# extract random sample ('spectrum')
rnd = outcomes$'random sample' %>%
  # tag y-value of each
  group_by(network, matrix) %>%
  arrange(cor) %>%
  mutate(y = row_number()) %>%
  ungroup()

# set up spectrum
spectrum = rnd %>%
  mutate(network = gsub("^.*\\|", "", network)) %>%
  filter(!network %in% c("Havugimana2012", "Wan2015")) %>%
  type_convert() %>%
  # manually recode networks
  mutate(network = fct_recode(network,
                              'HI-II-14' = 'Rolland2014',
                              'HuRI' = 'Luck2019',
                              'BioPlex' = 'Huttlin2015',
                              'BioPlex 2' = 'Huttlin2017',
                              'QUBIC' = 'Hein2015') %>% as.character())
spectrum_neg = spectrum %>%
  group_by(network, matrix) %>%
  summarise(y = sum(cor < 0), median = median(cor), mean = mean(cor)) %>%
  arrange(y)
spectrum_neg

# set up function
plot_spectrum = function(spectrum, grep_for = NULL, range = NULL,
                         digits = digits) {
  if (!is.null(grep_for))
    spectrum %<>% filter(grepl(grep_for, matrix))
  
  spectrum_neg = spectrum %>%
    group_by(network) %>%
    summarise(y = sum(cor < 0), median = median(cor), mean = mean(cor)) %>%
    arrange(y)
  
  lvls = spectrum_neg$network
  
  p1a = spectrum %>%
    mutate(network = factor(network, levels = lvls)) %>%
    ggplot(aes(x = network, y = y)) +
    facet_grid(network ~ ., scales = 'free') +
    geom_raster(aes(fill = cor)) +
    geom_errorbar(data = spectrum_neg %>%
                    mutate(network = factor(network, levels = lvls)),
                  aes(ymin = y, ymax = y), width = 1, size = 0.3,
                  color = 'grey20') +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                           name = 'Correlation', breaks = c(-1, 1),
                           limits = c(-1, 1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous('% of interactions', 
                       expand = c(0, 0), 
                       labels = function(x) paste0(x / 10, '%'),
                       limits = c(0, 1000),
                       breaks = seq(0, 1000, 250)) +
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
    geom_text(aes(label = format(median, digits = digits)), 
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

## ProteomeHD
p1 = plot_spectrum(spectrum, grep_for = "Kustatscher", digits = 1)
p1
ggsave("fig/analysis/validation/continuous_outcomes/spectrum-Kustatscher2019.pdf", 
       p1, width = 4.75, height = 4.75, units = "cm", useDingbats = FALSE)

## SubCellBarCode
p2 = plot_spectrum(spectrum, grep_for = "Orre", digits = 2)
p2
ggsave("fig/analysis/validation/continuous_outcomes/spectrum-Orre2019.pdf", 
       p2, width = 4.75, height = 4.75, units = "cm", useDingbats = FALSE)

## hyperLOPIT
p3 = plot_spectrum(spectrum, grep_for = "Geladaki", digits = 2)
p3
ggsave("fig/analysis/validation/continuous_outcomes/spectrum-Geladaki2018.pdf", 
       p3, width = 4.75, height = 4.75, units = "cm", useDingbats = FALSE)

## Lapek et al. 2017
p4 = plot_spectrum(spectrum, grep_for = "Lapek", digits = 1)
p4
ggsave("fig/analysis/validation/continuous_outcomes/spectrum-Lapek2017.pdf", 
       p4, width = 4.75, height = 4.75, units = "cm", useDingbats = FALSE)
