setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# complexes
dat = readRDS("data/analysis/downsample_fractions/complexes.rds") %>%
  type_convert() %>%
  filter(quant_mode == 'iBAQ') %>%
  # summarise over ten iterations
  group_by(accession, experiment, analysis, metric, transform, missing, 
           n_fractions, quant_mode) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), mean) %>%
  ungroup()

# filter to experiments with at least 50 fractions
cdf = readRDS("data/QC/protein-groups-CDF.rds")
n_fractions = cdf %>%
  group_by(accession, experiment, quant_mode) %>%
  summarise(fraction_count = max(n_fractions)) %>%
  ungroup()
expts = filter(n_fractions, fraction_count >= 50)
dat %<>% inner_join(expts)

# merge in fractionation
fractionation = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  dplyr::select(Accession, Replicate, Fractionation) %>%
  set_colnames(c("accession", "experiment", "fractionation"))
dat %<>% left_join(fractionation)

# calculate delta-auroc
delta = dat %>% 
  filter(quant_mode == 'iBAQ') %>%
  group_by(accession, experiment, fractionation, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc - lag(auroc)) %>%
  ungroup()
mean_delta = delta %>% 
  group_by(fractionation, n_fractions, metric) %>%
  summarise(mean = mean(delta_auroc), 
            median = median(delta_auroc),
            sd = sd(delta_auroc),
            sem = sd / sqrt(n()),
            q1 = quantile(delta_auroc, probs = 0.25, na.rm = TRUE),
            q3 = quantile(delta_auroc, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(color = 'Change in AUC')
mean_delta

# also get mean
mean_auroc = dat %>% 
  group_by(fractionation, n_fractions, metric) %>%
  summarise(mean = mean(auroc), 
            median = median(auroc),
            sd = sd(auroc),
            sem = sd / sqrt(n()),
            q1 = quantile(auroc, probs = 0.25, na.rm = TRUE),
            q3 = quantile(auroc, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(color = 'AUC')

# plot
# pal = colours.cafe322[c(1, 3)]
pal = pals::stepped3()[c(1, 5, 10, 15, 19) + 1] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
p1a = mean_auroc %>%
  filter(fractionation != 'IEF')%>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  mutate(quant_mode = 'iBAQ',
         color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_ribbon(aes(x = n_fractions, ymin = mean - sem,
                  ymax = mean + sem, fill = fractionation),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(x = n_fractions, y = mean, color = fractionation)) +
  geom_point(aes(x = n_fractions, y = mean, color = fractionation),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC') +
  scale_color_manual('', values = pal, drop = FALSE) +
  scale_fill_manual('', values = pal, drop = FALSE) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p1a
ggsave("fig/analysis/downsample_fractions/separation-AUC-CORUM-pearson.pdf",
       p1a, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

p1b = mean_delta %>%
  filter(fractionation != 'IEF')%>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  drop_na() %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_ribbon(aes(x = n_fractions, ymin = mean - sem,
                  ymax = mean + sem, fill = fractionation),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(x = n_fractions, y = mean, color = fractionation)) +
  geom_point(aes(x = n_fractions, y = mean, color = fractionation),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC') +
  scale_color_manual('', values = pal) +
  scale_fill_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p1b
ggsave("fig/analysis/downsample_fractions/separation-deltaAUC-CORUM-pearson.pdf",
      p1b, width = 5.5, height = 4.75, units = "cm", useDingbats = F)
