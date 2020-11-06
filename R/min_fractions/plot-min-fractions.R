setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

## complexes
dat1 = readRDS("data/analysis/min_fractions/complexes.rds") %>%
  filter(quant_mode == 'iBAQ')

# join with CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds")
## calculate maximum % of fractions per experiment
max_fractions = cdf %>%
  group_by(accession, experiment) %>%
  summarise(max_fractions = max(n_fractions)) %>%
  ungroup()
## calculate maximum # of 
cdf1 = filter(cdf, n_fractions == 1) %>%
  dplyr::select(-n_fractions, -coverage) %>%
  dplyr::rename(n_1 = n_proteins)

# calculate mean over # fractions
mean1 = dat1 %>%
  # fix quant mode and merge in CDF
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  left_join(cdf, by = c('accession', 'experiment', 'quant_mode', 
                        'n_fractions')) %>%
  # also merge in the total number of fractions
  left_join(max_fractions, by = c('accession', 'experiment')) %>%
  # drop experiments without at least 20 fractions
  filter(max_fractions >= 20) %>%
  # merge in number of proteins found in one fraction
  left_join(cdf1, by = c("accession", "experiment", "quant_mode")) %>%
  # calculate mean AUC and # of proteins
  group_by(metric, quant_mode, n_fractions) %>%
  summarise(mean = mean(auroc, na.rm = T),
            sd = sd(auroc, na.rm = T),
            median = median(auroc, na.rm = T),
            q1 = quantile(auroc, na.rm = T, probs = 0.25),
            q3 = median(auroc, na.rm = T, probs = 0.75),
            pct_proteins = mean(n_proteins / n_1),
            n_proteins = mean(n_proteins),
            n = n()) %>%
  ungroup()

# plot mean and s.d.: AUC
pal = brewer.pal(8, 'Set2') %>% extract(c(3, 2, 4, 1, 7, 8))
lab1a = mean1 %>%
  filter(n_fractions == 4, metric == 'pearson') %>%
  mutate(label = paste0(n_fractions, ' fractions\nAUC = ', 
                        format(mean, format = 'f', digits = 2)))
p1a = mean1 %>%
  filter(metric == 'pearson') %>%
  ggplot(aes(x = n_fractions, y = mean, color = quant_mode)) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  geom_label_repel(data = lab1a, aes(label = label), size = 2, fill = NA,
                   color = 'black', label.padding = 0.4, hjust = 0,
                   nudge_x = 0.5, nudge_y = -0.03, label.size = NA, 
                   segment.size = 0.25) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions') +
  scale_y_continuous('AUC', # limits = c(0.6575, 0.6815),
                     breaks = seq(0.66, 0.68, 0.005)) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none')
p1a

# plot mean and s.d.: proteins quantified, as a % of max
lab1b = mean1 %>%
  filter(n_fractions == 4, metric == 'pearson') %>%
  mutate(label = paste0(n_fractions, ' fractions\n', 
                        format(100 * pct_proteins, format = 'f', digits = 3),
                        '% of proteins'))
p1b = mean1 %>%
  ggplot(aes(x = n_fractions, y = pct_proteins, color = quant_mode)) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  geom_label_repel(data = lab1b, aes(label = label), size = 2, 
                   color = 'black', label.padding = 0.4, hjust = 0, fill = NA,
                   nudge_x = 8, nudge_y = 0.1, label.size = NA, 
                   segment.size = 0.25) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions', breaks = c(1, seq(5, 20, 5))) +
  scale_y_continuous('% of proteins', labels = function(x) x * 100,
                     limits = c(0.4, 1)
                     ) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none')
p1b

# combine
p1 = (p1a + theme(axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   plot.margin = margin(rep(0, 4)))) +
  p1b + 
  plot_layout(ncol = 1, heights = c(1.15, 0.85))
p1
ggsave("fig/analysis/min_fractions/min_fractions-CORUM.pdf", p1,
       width = 5, height = 6, units = "cm", useDingbats = F)

# print stats at n=3-4-5 fractions
filter(mean1, between(n_fractions, 3, 5))
