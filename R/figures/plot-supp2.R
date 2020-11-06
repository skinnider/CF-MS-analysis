# Plot Supplementary Figure 2:
#' a. previously used metrics (elsewhere)
#' b. AUCs, metrics: single best preprocessing
#' c. AUCs, pipelines
#' d. min. fractions, MI
#' e. quant. mode, GO
#' f. AUCs, individual experiments
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

###############################################################################-
###### c. Min fractions ####
###############################################################################-

datC = readRDS("data/analysis/min_fractions/complexes.rds")

# join with CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds")
## calculate maximum % of fractions per experiment
max_fractions = cdf %>%
  group_by(accession, experiment) %>%
  summarise(max_fractions = max(n_fractions)) %>%
  ungroup()
## calculate maximum # of proteins per experiment
cdf1 = filter(cdf, n_fractions == 1) %>%
  dplyr::select(-n_fractions, -coverage) %>%
  dplyr::rename(n_1 = n_proteins)

# calculate mean over # fractions
meanC = datC %>%
  filter(metric == 'MI') %>%
  # fix quant mode and merge in CDF
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  left_join(cdf, by = c('accession', 'experiment', 'quant_mode',
                        'n_fractions')) %>%
  # also merge in the total number of fractions
  left_join(max_fractions, by = c('accession', 'experiment')) %>%
  # drop experiments without at least 20 fractions
  filter(max_fractions >= 20) %>%
  # merge in number of proteins found in one fraction
  left_join(cdf1, by = c("accession", "experiment", "search", "quant_mode")) %>%
  # calculate mean AUC and # of proteins
  group_by(quant_mode, n_fractions) %>%
  summarise(mean = mean(auroc, na.rm = T),
            sd = sd(auroc, na.rm = T),
            median = median(auroc, na.rm = T),
            q1 = quantile(auroc, na.rm = T, probs = 0.25),
            q3 = median(auroc, na.rm = T, probs = 0.75),
            pct_proteins = mean(n_proteins / n_1),
            n_proteins = mean(n_proteins),
            n = n()) %>%
  ungroup()

# plot AUC
pal = brewer.pal(8, 'Set2') %>% extract(c(3, 2, 4, 1, 7, 8))
labCa = meanC %>%
  filter(quant_mode == 'iBAQ', n_fractions == 4) %>%
  mutate(label = paste0(n_fractions, ' fractions\nAUC = ', 
                        format(mean, format = 'f', digits = 2)))
pCa = meanC %>%
  filter(quant_mode == 'iBAQ') %>%
  filter(!quant_mode %in% c('MS1 intensity', 'ratio')) %>%
  ggplot(aes(x = n_fractions, y = mean, color = '1')) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions') +
  scale_y_continuous('AUC', # limits = c(0.6575, 0.6815),
                     breaks = seq(0.658, 0.68, 0.002)) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none') + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(rep(0, 4)))
pCa

# plot # of proteins
labCb = meanC %>%
  filter(quant_mode == 'iBAQ', n_fractions == 4) %>%
  mutate(label = paste0(n_fractions, ' fractions\n', 
                        format(100 * pct_proteins, format = 'f', digits = 3),
                        '% of proteins'))
pCb = meanC %>%
  filter(quant_mode == 'iBAQ') %>%
  ggplot(aes(x = n_fractions, y = pct_proteins, color = '1')) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions', breaks = c(1, seq(5, 20, 5))) +
  scale_y_continuous('% of proteins', labels = function(x) x * 100,
                     limits = c(0.4, 1)
  ) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none')
pCb

# combine
pC = pCa + pCb + plot_layout(ncol = 1, heights = c(1.15, 0.85))
pC

###############################################################################-
###### d. LFQ ####
###############################################################################-

# read GO
go = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup()

# filter to GO terms that pass breadth cutoff in all quant. modes
filtered = go %>%
  # ignore MS1 intensity
  filter(quant_mode != 'MS1_intensity') %>%
  # filter by breadth
  filter(between(n_chromatograms, 10, 100)) %>%
  # keep the intersect of all three quant. modes
  dplyr::select(-n_proteins, -n_chromatograms, -auroc) %>%
  group_by_at(vars(-go_term, -quant_mode)) %>%
  filter(n_distinct(quant_mode) == 3) %>%
  ungroup()

# get median AUCs
datD = go %>%
  inner_join(filtered) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# read CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup() 

# first, plot CDF: mean # of protein quantitations from each method
max_quant = cdf %>%
  group_by(accession, experiment) %>%
  summarise(max_proteins = max(n_proteins),
            max_coverage = max(coverage)) %>%
  ungroup()
mean_cdf = cdf %>%
  group_by(quant_mode, n_fractions) %>%
  summarise(mean_proteins = mean(n_proteins),
            mean_coverage = mean(coverage)) %>%
  ungroup()
pal = brewer.pal(8, 'Set2') %>% extract(c(2, 1, 4))

## all datasets
pD1 = mean_cdf %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  filter(n_fractions <= 50) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(x = n_fractions, y = mean_proteins, color = quant_mode)) +
  geom_line(size = 0.5) +
  scale_y_continuous("# of proteins", breaks = seq(0, 3000, 500)) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme()
pD1

### GO
mediansD = datD %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  group_by(quant_mode) %>%
  summarise(median = median(auroc, na.rm = T),
            label = formatC(median, format = 'f', digits = 3))
pD2 = datD %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(y = reorder(quant_mode, auroc, stats::median),
             x = auroc, fill = quant_mode, color = quant_mode)) +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = mediansD, aes(x = 0.48, y = quant_mode, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC') +
  coord_cartesian(xlim = c(0.48, 0.55)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
pD2

# combine and save
pD = pD1 + pD2 + plot_layout(ncol = 1, heights = c(1, 0.4))
pD

###############################################################################-
###### combine all plots
###############################################################################-

p_bl = ggplot() + theme_void()
part2 = pCa + pCb + p_bl + pD1 + pD2 + 
  plot_layout(ncol = 1, heights = c(1.1, 0.75, 1.25, 1, 0.4))
ggsave("fig/final/supp2/part2.pdf", part2, height = 15, width = 5.6, 
       units = "cm", useDingbats = F)
