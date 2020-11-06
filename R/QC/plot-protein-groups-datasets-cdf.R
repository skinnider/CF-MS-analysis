# Plot a CDF of the number of _experiments_ containing at least n 
# protein groups.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")

# read CDF
cdf = readRDS('data/QC/protein-groups-CDF.rds')

# reshape to # datasets/# of proteins
dataset_cdf = tidyr::crossing(n_proteins = seq(0, 8e3, 50),
                              n_fractions = c(seq_len(10), 15, 20, 25))
dataset_cdf$n_datasets = map2_int(
  dataset_cdf$n_proteins, dataset_cdf$n_fractions,
  ~ cdf %>%
    filter(n_fractions == .y) %>%
    filter(n_proteins > .x) %>%
    distinct(accession, experiment) %>%
    nrow()
)

# plot
pal = BuenColors::jdb_palette("wolfgang_extra") %>% colorRampPalette()
p1 = dataset_cdf %>%
  filter(n_fractions %in% c(1, 5, 10, 15, 20, 25)) %>%
  ggplot(aes(x = n_proteins, y = n_datasets, color = factor(n_fractions))) +
  # geom_hline(aes(yintercept = 209), size = 0.3, linetype = 'dotted') +
  geom_path() + 
  scale_color_manual("min. # of fractions", values = pal(7)[-1]) +
  scale_x_continuous("# of proteins") +
  scale_y_continuous("# of datasets") +
  coord_cartesian(xlim = c(0, 8e3)) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              nrow = 1)) +
  boxed_theme()
p1
ggsave("fig/QC/protein-groups-cdf-datasets.pdf", p1,
       width = 6, height = 5.5, units = "cm", useDingbats = F)
