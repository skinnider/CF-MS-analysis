# Plot the impact of different quantitation strategies on datasets
# generated using label-free approaches.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read complexes
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  # restore 'missing=NA'
  replace_na(list(missing = 'NA')) %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile')) %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup()

# read GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))
# for each combination, filter by breadth, and get median AUROC across GO terms
dat2 %<>%
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup() %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
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
p1a = mean_cdf %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  filter(n_fractions <= 50) %>%
  ## remove MS1: same # quants. as iBAQ
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(x = n_fractions, y = mean_proteins, color = quant_mode)) +
  geom_line(size = 0.5) +
  # geom_point(size = 0.8) + 
  scale_y_continuous("# of proteins", breaks = seq(0, 3000, 500)) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme()
p1a

# repeat for human/mouse only
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species %in% c("Homo sapiens", "Mus musculus"))
mean_cdf2 = cdf %>%
  filter(paste(accession, experiment) %in% 
           with(expts, paste(Accession, Replicate))) %>%
  group_by(quant_mode, n_fractions) %>%
  summarise(mean_proteins = mean(n_proteins),
            mean_coverage = mean(coverage)) %>%
  ungroup()
pal = brewer.pal(8, 'Set2') %>% extract(c(2, 1, 4))
p1a2 = mean_cdf2 %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  filter(n_fractions <= 50) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(x = n_fractions, y = mean_proteins, color = quant_mode)) +
  geom_line(size = 0.5) +
  # geom_point(size = 0.8) + 
  scale_y_continuous("# of proteins", breaks = seq(0, 3000, 500)) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme()
p1a2

## also plot mean coverage
p1b = mean_cdf %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  filter(n_fractions <= 50) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(x = n_fractions, y = mean_coverage, color = quant_mode)) +
  geom_line(size = 0.5) +
  # geom_point(size = 0.8) + 
  scale_y_continuous("Proteome coverage (%)", labels = function(x) x * 100) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme()
p1b

## next, boxplots
### complexes: CORUM
medians2a = dat1 %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  group_by(quant_mode) %>%
  summarise(median = median(auroc, na.rm = T),
            label = formatC(median, format = 'f', digits = 3))
p2a = dat1 %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  ggplot(aes(y = reorder(quant_mode, auroc, stats::median),
             x = auroc, fill = quant_mode, color = quant_mode)) +
  # geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
  #            linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2a, aes(x = 0.47, y = quant_mode, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC') +
  coord_cartesian(xlim = c(0.475, 0.7)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2a

# combine and save
p2 = p1a2 + p2a + plot_layout(ncol = 1, heights = c(1, 0.4))
p2
ggsave("fig/analysis/analysis_grid/label-free-CORUM.pdf", p2,
       width = 6, height = 6.5, units = "cm", useDingbats = F)

### GO
medians2c = dat2 %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  group_by(quant_mode) %>%
  summarise(median = median(auroc, na.rm = T),
            label = formatC(median, format = 'f', digits = 3))
p2c = dat2 %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(y = reorder(quant_mode, auroc, stats::median),
             x = auroc, fill = quant_mode, color = quant_mode)) +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2c, aes(x = 0.478, y = quant_mode, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC') +
  coord_cartesian(xlim = c(0.48, 0.55)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2c

# combine and save
p3 = p1a + p2c + plot_layout(ncol = 1, heights = c(1, 0.4))
p3
ggsave("fig/analysis/analysis_grid/label-free-GO.pdf", p3,
       width = 6, height = 6.5, units = "cm", useDingbats = F)
