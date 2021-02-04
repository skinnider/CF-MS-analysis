setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

dat = readRDS("data/analysis/sliding_window/complexes.rds") %>%
  type_convert() %>% 
  dplyr::select(-transform, -missing) %>%
  drop_na(auroc)

# filter to experiments with at least 50 fractions
cdf = readRDS("data/QC/protein-groups-CDF.rds")
n_fractions = cdf %>%
  group_by(accession, experiment, quant_mode) %>%
  summarise(fraction_count = max(n_fractions)) %>%
  ungroup()
expts = filter(n_fractions, fraction_count >= 50, quant_mode == 'iBAQ') %>%
  distinct(accession, experiment)
dat %<>% inner_join(expts) 

# calculate delta-auroc between adjacent sizes
delta = dat %>% 
  filter(quant_mode == 'iBAQ') %>%
  group_by(accession, experiment, metric) %>%
  arrange(window_size) %>%
  mutate(delta_auroc = auroc - lag(auroc)) %>%
  ungroup()
mean_delta = delta %>% 
  group_by(window_size, metric) %>%
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
  group_by(window_size, metric) %>%
  summarise(mean = mean(auroc), 
            median = median(auroc),
            sd = sd(auroc),
            sem = sd / sqrt(n()),
            q1 = quantile(auroc, probs = 0.25, na.rm = TRUE),
            q3 = quantile(auroc, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(color = 'AUC')

# plot
pal = colours.cafe322[c(1, 3)]
p1a = mean_auroc %>%
  filter(metric == 'pearson') %>%
  # filter(window_size <= 50) %>%
  mutate(quant_mode = 'iBAQ',
         color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_ribbon(aes(x = window_size, ymin = mean - sem,
                  ymax = mean + sem, fill = color),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(x = window_size, y = mean, color = color)) +
  geom_point(aes(x = window_size, y = mean, color = color),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', # limits = c(0.58, 0.67),
                     # breaks = seq(0.58, 0.68, 0.02)
                     ) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('AUC', 'Change in AUC')) +
  scale_fill_manual('', values = pal, drop = FALSE,
                    labels = c('AUC', 'Change in AUC')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p1a

p1b = mean_delta %>%
  filter(metric == 'pearson') %>%
  mutate(quant_mode = 'iBAQ') %>% 
  ggplot() + 
  geom_ribbon(aes(x = window_size, ymin = mean - sem,
                  ymax = mean + sem, fill = color),
              alpha = 0.3, show.legend = FALSE) +
  geom_line(aes(x = window_size, y = mean, color = quant_mode)) +
  geom_point(aes(x = window_size, y = mean, color = quant_mode),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', position = 'right') +
  scale_color_manual('', values = pal[-1]) +
  scale_fill_manual('', values = pal[-1]) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1)
p1b
# superimpose in cowplot
aligned_plots = align_plots(p1a, p1b, align = "hv", axis = "tblr")
p1 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p1
ggsave("fig/analysis/sliding_window/AUC-deltaAUC-CORUM-pearson.pdf",
       p1, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

# compare to downsampling randomly
dat_rnd = readRDS("data/analysis/downsample_fractions/complexes.rds") %>%
  type_convert() %>%
  filter(quant_mode == 'iBAQ') %>%
  # summarise over ten iterations
  group_by(accession, experiment, analysis, metric, transform, missing, 
           n_fractions, quant_mode) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), mean) %>%
  ungroup()
# filter to experiments with at least 50 fractions
dat_rnd %<>% inner_join(expts)
# calculate delta-auroc
delta_rnd = dat_rnd %>% 
  filter(quant_mode == 'iBAQ') %>%
  group_by(accession, experiment, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc - lag(auroc)) %>%
  ungroup()
mean_delta_rnd = delta_rnd %>% 
  group_by(n_fractions, metric) %>%
  summarise(mean = mean(delta_auroc), 
            median = median(delta_auroc),
            sd = sd(delta_auroc),
            sem = sd / sqrt(n()),
            q1 = quantile(delta_auroc, probs = 0.25, na.rm = TRUE),
            q3 = quantile(delta_auroc, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(color = 'Change in AUC')
# also get mean
mean_auroc_rnd = dat_rnd %>% 
  group_by(n_fractions, metric) %>%
  summarise(mean = mean(auroc), 
            median = median(auroc),
            sd = sd(auroc),
            sem = sd / sqrt(n()),
            q1 = quantile(auroc, probs = 0.25, na.rm = TRUE),
            q3 = quantile(auroc, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(color = 'AUC')

# combine
comp_mean = bind_rows(
  mutate(mean_auroc, downsampling = 'adjacent')%>% 
    dplyr::rename(n_fractions = window_size),
  mutate(mean_auroc_rnd, downsampling = 'random') 
)
pal = c('grey70', colours.cafe322[c(1)])
p2a = comp_mean %>%
  filter(n_fractions <= 50) %>% 
  filter(metric == 'pearson') %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_ribbon(aes(x = n_fractions, ymin = mean - sem,
                  ymax = mean + sem, fill = downsampling),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(x = n_fractions, y = mean, color = downsampling)) +
  geom_point(aes(x = n_fractions, y = mean, color = downsampling),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', limits = c(0.58, 0.67),
                     breaks = seq(0.58, 0.68, 0.02)
  ) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('Adjacent fractions', 'Random fractions')) +
  scale_fill_manual('', values = pal, drop = FALSE,
                    labels = c('Adjacent fractions', 'Random fractions')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p2a
ggsave("fig/analysis/sliding_window/AUC-vs-random-CORUM-pearson.pdf",
       p2a, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

# delta-AUC
comp_delta = bind_rows(
  mutate(mean_delta, downsampling = 'adjacent') %>% 
    dplyr::rename(n_fractions = window_size),
  mutate(mean_delta_rnd, downsampling = 'random') 
)
pal = colours.cafe322[c(2, 1)]
# pal = brewer.pal(8, 'PuBu')[c(6, 8)]
pal = c('grey70', colours.cafe322[c(1)])
p2b = comp_delta %>%
  drop_na() %>%
  filter(n_fractions <= 50) %>% 
  filter(metric == 'pearson') %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_ribbon(aes(x = n_fractions, ymin = mean - sem,
                  ymax = mean + sem, fill = downsampling),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(x = n_fractions, y = mean, color = downsampling)) +
  geom_point(aes(x = n_fractions, y = mean, color = downsampling),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', # limits = c(0.58, 0.67),
                     # breaks = seq(0.58, 0.68, 0.02)
  ) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('Adjacent fractions', 'Random fractions')) +
  scale_fill_manual('', values = pal, drop = FALSE,
                    labels = c('Adjacent fractions', 'Random fractions')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p2b
ggsave("fig/analysis/sliding_window/deltaAUC-vs-random-CORUM-pearson.pdf",
       p2b, width = 5.5, height = 4.75, units = "cm", useDingbats = F)
