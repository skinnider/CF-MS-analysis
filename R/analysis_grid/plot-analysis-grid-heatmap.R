setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# 1. Complexes ####
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # restore 'missing=NA'
  replace_na(list(missing = 'NA')) %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))

# calculate median
medians1 = dat1 %>%
  group_by(metric, missing, transform) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup()

# plot 
grey = colorspace::sequential_hcl("Greens 2", n = 100)[100]
df1 = medians1 %>%
  tidyr::complete(metric, missing, transform, 
                  fill = list(median = NA)) %>%
  mutate(transform = factor(chartr(' ', '\n', transform), 
                            levels = c('none', 'log-transform',
                                       'quantile\nnormalization')),
         missing = factor(missing, levels = c('zero', 'noise', 'NA')),
         metric = clean_metric(metric),
         metric = reorder(metric, median, function(x) max(x, na.rm = TRUE)))
labels1 = df1 %>%
  drop_na() %>%
  mutate(label = formatC(median, digits = 2, format = 'f'))
p1 = df1 %>%
  ggplot(aes(x = missing, y = metric, fill = median)) +
  facet_grid(~ transform) +
  geom_tile(color = 'white') + 
  geom_text(data = labels1, aes(label = label), size = 1.75, color = 'white') +
  scale_x_discrete("Missing data", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_paletteer_c("pals::coolwarm", na.value = grey,
                         name = 'AUC, median', breaks = c(0.45, 0.65)) +
  ggtitle("Protein complexes") +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() + 
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 7),
        legend.key.height = unit(0.3, 'lines'),
        legend.key.width = unit(0.25, 'lines'),
        legend.position = 'right')
p1

# 2. GO ####
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

# calculate median
medians2 = dat2 %>%
  group_by(metric, missing, transform) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup()

# plot 
grey = colorspace::sequential_hcl("Greens 2", n = 100)[100]
df2 = medians2 %>%
  mutate(transform = factor(chartr(' ', '\n', transform), 
                            levels = c('none', 'log-transform',
                                       'quantile\nnormalization')),
         missing = factor(missing, levels = c('zero', 'noise', 'NA')),
         metric = clean_metric(metric)) 
labels2 = df2 %>%
  drop_na() %>%
  mutate(label = formatC(median, digits = 2, format = 'f'))
p2 = df2 %>%
  tidyr::complete(metric, missing, transform, 
                  fill = list(median = NA)) %>%
  ggplot(aes(x = missing, y = reorder(metric, median, 
                                      function(x) max(x, na.rm = T)),
             fill = median)) +
  facet_grid(~ transform) +
  geom_tile(color = 'white') + 
  geom_text(data = labels2, aes(label = label), size = 1.75, color = 'white') +
  scale_x_discrete("Missing data", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_paletteer_c("pals::coolwarm", na.value = grey,
                         name = 'AUC, median', breaks = c(0.5, 0.55)) +
  ggtitle("GO") +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() + 
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 7),
        legend.key.height = unit(0.3, 'lines'),
        legend.key.width = unit(0.25, 'lines'),
        legend.position = 'right')
p2

# combine and save
p = p1 | p2
ggsave("fig/analysis/analysis_grid/AUCs-heatmaps.pdf", p, 
       width = 20, height = 39, units = "cm", useDingbats = F)

# also save the medians for the replicate integration experiment
x1 = medians1 %>%
  arrange(desc(median)) %>%
  dplyr::select(metric, missing, transform, median)
x2 = medians2 %>%
  arrange(desc(median)) %>%
  dplyr::select(metric, missing, transform, median)
x3 = left_join(x1, x2, by = c('metric', 'missing', 'transform')) %>%
  mutate(median = median.x + median.y) %>%
  dplyr::select(-median.x, -median.y) %>%
  arrange(desc(median))
medians = list('complexes' = x1, 'GO' = x2, 'combined' = x3)
saveRDS(medians, "data/analysis/analysis_grid/best_features.rds")
