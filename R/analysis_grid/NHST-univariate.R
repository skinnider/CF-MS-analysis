setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
library(lawstat)
source("R/theme.R")

# read complexes
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile') %>% 
           as.character())

# print summary of the data
dplyr::count(dat1, metric, missing, transform)

# measure of association
## reorder metrics as in main text panel
metrics = with(dat1, reorder(metric, auroc, stats::median)) %>% levels()
pairs = tidyr::crossing(metric1 = metrics, metric2 = metrics) %>% 
  filter(metric1 != metric2)
df1 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(dat1, metric %in% c(current$metric1, current$metric2))
  x = filter(df, metric == current$metric1) %>% pull(auroc)
  y = filter(df, metric == current$metric2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(metric1 = factor(metric1, levels = metrics),
         metric2 = factor(metric2, levels = metrics)) 
sig1 = df1 %>% 
  filter(pval < 0.05)
# range = range(df1$delta)
range = c(-0.1, 0.1) ## cap range for better visualization
p1 = df1 %>% 
  mutate(delta = winsorize(delta, range)) %>% 
  ggplot(aes(x = metric1, y = metric2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = sig1, aes(label = '*'), size = 1.75, color = 'white',
            nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p1
ggsave("fig/analysis/analysis_grid/measure-of-association-BM.pdf", p1,
       width = 9.25, height = 9.25, units = 'cm', useDingbats = FALSE)

# repeat using the single best combination only
df2 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(dat1, metric %in% c(current$metric1, current$metric2))
  # extract the best combination for each metric
  best = df %>% 
    group_by(analysis, metric, missing, transform) %>% 
    summarise(median = median(auroc)) %>% 
    ungroup() %>% 
    group_by(metric) %>% 
    arrange(desc(median)) %>% 
    dplyr::slice(1) %>% 
    ungroup()
  df %<>% inner_join(best)
  x = filter(df, metric == current$metric1) %>% pull(auroc)
  y = filter(df, metric == current$metric2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(metric1 = factor(metric1, levels = metrics),
         metric2 = factor(metric2, levels = metrics)) 
sig2 = df2 %>% 
  filter(pval < 0.05)
range = range(df2$delta)
p2 = df2 %>% 
  ggplot(aes(x = metric1, y = metric2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  # geom_text(aes(label = label1), size = 1.75) +
  geom_text(data = sig2, aes(label = '*'), size = 1.75, color = 'white',
            nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p2
ggsave("fig/analysis/analysis_grid/measure-of-association-BM-bests.pdf", p2,
       width = 9.25, height = 9.25, units = 'cm', useDingbats = FALSE)

# repeat for GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile') %>% 
           as.character())
dat2 %<>%
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# measure of association
metrics = with(dat2, reorder(metric, auroc, stats::median)) %>% levels()
df3 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(dat2, metric %in% c(current$metric1, current$metric2))
  x = filter(df, metric == current$metric1) %>% pull(auroc)
  y = filter(df, metric == current$metric2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(metric1 = factor(metric1, levels = metrics),
         metric2 = factor(metric2, levels = metrics)) 
sig3 = df3 %>% 
  filter(pval < 0.05)
range = range(df3$delta)
p3 = df3 %>% 
  mutate(delta = winsorize(delta, range)) %>% 
  ggplot(aes(x = metric1, y = metric2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = sig3, aes(label = '*'), size = 1.75, color = 'white',
            nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p3
ggsave("fig/analysis/analysis_grid/measure-of-association-BM-GO.pdf", p3,
       width = 9.25, height = 9.25, units = 'cm', useDingbats = FALSE)

# repeat using the single best combination
df4 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(dat2, metric %in% c(current$metric1, current$metric2))
  best = df %>% 
    group_by(analysis, metric, missing, transform) %>% 
    summarise(median = median(auroc)) %>% 
    ungroup() %>% 
    group_by(metric) %>% 
    arrange(desc(median)) %>% 
    dplyr::slice(1) %>% 
    ungroup()
  df %<>% inner_join(best)
  x = filter(df, metric == current$metric1) %>% pull(auroc)
  y = filter(df, metric == current$metric2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(metric1 = factor(metric1, levels = metrics),
         metric2 = factor(metric2, levels = metrics)) 
sig4 = df4 %>% 
  filter(pval < 0.05)
range = range(df4$delta)
p4 = df4 %>% 
  ggplot(aes(x = metric1, y = metric2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  # geom_text(aes(label = label1), size = 1.75) +
  geom_text(data = sig4, aes(label = '*'), size = 1.75, color = 'white',
            nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p4
ggsave("fig/analysis/analysis_grid/measure-of-association-BM-bests-GO.pdf", p4,
       width = 9.25, height = 9.25, units = 'cm', useDingbats = FALSE)

# normalization
## remove co-occurrence first as preprocessing does not apply
x = filter(dat1, transform == 'none', !metric %in% co_occurrence)
y = filter(dat1, transform == 'log-transform', !metric %in% co_occurrence)
z = filter(dat1, transform == 'quantile normalization', 
           !metric %in% co_occurrence)

# check data completeness
dplyr::count(x, metric, missing) %$% table(n)
dplyr::count(y, metric, missing) %$% table(n)
dplyr::count(y, metric, missing) %>% filter(n < 67) ## a handful of treeClust with NA did not converge
dplyr::count(z, metric, missing) %$% table(n)

# run tests
library(lawstat)
pval1 = brunner.munzel.test(x$auroc, y$auroc)$p.value
delta1 = median(y$auroc) - median(x$auroc)
pval2 = brunner.munzel.test(x$auroc, z$auroc)$p.value
delta2 = median(z$auroc) - median(x$auroc)
pval3 = brunner.munzel.test(y$auroc, z$auroc)$p.value
delta3 = median(z$auroc) - median(y$auroc)

# create data frame and plot
df5 = data.frame(transformation1 = c('none', 'none', 'log-transform'),
                 transformation2 = c('log-transform', 'quantile normalization',
                                     'quantile normalization'),
                 delta = c(delta1, delta2, delta3),
                 pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df5 %<>% bind_rows(df5 %>% 
                     mutate(tmp = transformation2,
                            transformation2 = transformation1,
                            transformation1 = tmp,
                            delta = -delta))
range = range(df5$delta)
lvls = c('quantile normalization', 'log-transform', 'none') ## as in Fig. 2e
p5 = df5 %>% 
  ggplot(aes(x = factor(transformation1, levels = lvls),
             y = factor(transformation2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df5, sig), aes(label = '*'), size = 1.75, 
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p5
ggsave("fig/analysis/analysis_grid/normalization-BM.pdf", p5,
       width = 4, height = 4, units = 'cm', useDingbats = FALSE)

# repeat for GO
x = filter(dat2, transform == 'none', !metric %in% co_occurrence)
y = filter(dat2, transform == 'log-transform', !metric %in% co_occurrence)
z = filter(dat2, transform == 'quantile normalization', !metric %in% co_occurrence)
dplyr::count(x, metric, missing) %$% table(n)
dplyr::count(y, metric, missing) %$% table(n)
dplyr::count(y, metric, missing) %>% filter(n < 206) ## same handful of treeClust with NA
dplyr::count(z, metric, missing) %$% table(n)

# run tests
pval1 = brunner.munzel.test(x$auroc, y$auroc)$p.value
delta1 = median(y$auroc) - median(x$auroc)
pval2 = brunner.munzel.test(x$auroc, z$auroc)$p.value
delta2 = median(z$auroc) - median(x$auroc)
pval3 = brunner.munzel.test(y$auroc, z$auroc)$p.value
delta3 = median(z$auroc) - median(y$auroc)

# create data frame and plot
df6 = data.frame(transformation1 = c('none', 'none', 'log-transform'),
                 transformation2 = c('log-transform', 'quantile normalization',
                                     'quantile normalization'),
                 delta = c(delta1, delta2, delta3),
                 pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df6 %<>% bind_rows(df6 %>% 
                     mutate(tmp = transformation2,
                            transformation2 = transformation1,
                            transformation1 = tmp,
                            delta = -delta))
range = range(df6$delta)
p6 = df6 %>% 
  ggplot(aes(x = factor(transformation1, levels = lvls),
             y = factor(transformation2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df6, sig), aes(label = '*'), size = 1.75,
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 3)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p6
ggsave("fig/analysis/analysis_grid/normalization-BM-GO.pdf", p6,
       width = 4, height = 4, units = 'cm', useDingbats = FALSE)

# missing values
## again need to remove co-occurrence
x = filter(dat1, missing == 'zero', !metric %in% co_occurrence)
y = filter(dat1, missing == 'NA', !metric %in% co_occurrence)
z = filter(dat1, missing == 'noise', !metric %in% co_occurrence)

# check data completeness
dplyr::count(x, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %>% filter(n < 67) ## treeClust convergence
dplyr::count(z, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %>% filter(n < 67) ## treeClust convergence

# run tests
library(lawstat)
pval1 = brunner.munzel.test(x$auroc, y$auroc)$p.value
delta1 = median(y$auroc) - median(x$auroc)
pval2 = brunner.munzel.test(x$auroc, z$auroc)$p.value
delta2 = median(z$auroc) - median(x$auroc)
pval3 = brunner.munzel.test(y$auroc, z$auroc)$p.value
delta3 = median(z$auroc) - median(y$auroc)

# create data frame and plot
df7 = data.frame(missing1 = c('zero', 'zero', 'NA'),
                 missing2 = c('NA', 'noise', 'noise'),
                 delta = c(delta1, delta2, delta3),
                 pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df7 %<>% bind_rows(df7 %>% 
                     mutate(tmp = missing2,
                            missing2 = missing1,
                            missing1 = tmp,
                            delta = -delta))
range = range(df7$delta)
lvls = c('NA', 'noise', 'zero') ## as in Fig. 2e
p7 = df7 %>% 
  ggplot(aes(x = factor(missing1, levels = lvls),
             y = factor(missing2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df7, sig), aes(label = '*'), size = 1.75,
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p7
ggsave("fig/analysis/analysis_grid/missing-BM.pdf", p7,
       width = 2.9, height = 2.9, units = 'cm', useDingbats = FALSE)

## repeat for GO
x = filter(dat2, missing == 'zero', !metric %in% co_occurrence)
y = filter(dat2, missing == 'NA', !metric %in% co_occurrence)
z = filter(dat2, missing == 'noise', !metric %in% co_occurrence)
dplyr::count(x, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %>% filter(n < 206) ## treeClust NA
dplyr::count(z, metric, transform) %$% table(n)
dplyr::count(y, metric, transform) %>% filter(n < 206) ## treeClust NA

# run tests
library(lawstat)
pval1 = brunner.munzel.test(x$auroc, y$auroc)$p.value
delta1 = median(y$auroc) - median(x$auroc)
pval2 = brunner.munzel.test(x$auroc, z$auroc)$p.value
delta2 = median(z$auroc) - median(x$auroc)
pval3 = brunner.munzel.test(y$auroc, z$auroc)$p.value
delta3 = median(z$auroc) - median(y$auroc)

# create data frame and plot
df8 = data.frame(missing1 = c('zero', 'zero', 'NA'),
                 missing2 = c('NA', 'noise', 'noise'),
                 delta = c(delta1, delta2, delta3),
                 pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df8 %<>% bind_rows(df8 %>% 
                     mutate(tmp = missing2,
                            missing2 = missing1,
                            missing1 = tmp,
                            delta = -delta))
range = range(df8$delta)
lvls = c('NA', 'noise', 'zero') ## as in Fig. 2e
p8 = df8 %>% 
  ggplot(aes(x = factor(missing1, levels = lvls),
             y = factor(missing2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df8, sig), aes(label = '*'), size = 1.75,
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p8
ggsave("fig/analysis/analysis_grid/missing-BM-GO.pdf", p8,
       width = 2.9, height = 2.9, units = 'cm', useDingbats = FALSE)

## finally, compare label-free quantification
## here we can use the paired Brunner-Munzel test
# read complexes
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup() %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile') %>% 
           as.character())

# print summary of the data
dplyr::count(dat1, metric, missing, transform)

# run paired tests
x = filter(dat1, quant_mode == 'iBAQ')
y = filter(dat1, quant_mode == 'LFQ')
z = filter(dat1, quant_mode == 'spectral_counts')
# get delta-medians
delta1 = left_join(x, y, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
delta2 = left_join(x, z, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
delta3 = left_join(y, z, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
# paired BM test
pval1 = bind_rows(x, y) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
pval2 = bind_rows(x, z) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
pval3 = bind_rows(y, z) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
# double-check paired t-tests, as well
left_join(x, y, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)
left_join(x, z, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)
left_join(y, z, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)

# create data frame and plot
df9 = data.frame(quant_mode1 = c('iBAQ', 'iBAQ', 'MaxLFQ'),
                 quant_mode2 = c('MaxLFQ', 'spectral count', 'spectral count'),
                 delta = c(delta1, delta2, delta3),
                 pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df9 %<>% bind_rows(df9 %>% 
                     mutate(tmp = quant_mode2,
                            quant_mode2 = quant_mode1,
                            quant_mode1 = tmp,
                            delta = -delta))
range = range(df9$delta)
lvls = c('MaxLFQ', 'spectral count', 'iBAQ') ## as in Fig. 2j
p9 = df9 %>% 
  ggplot(aes(x = factor(quant_mode1, levels = lvls),
             y = factor(quant_mode2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df9, sig), aes(label = '*'), size = 1.75,
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 2)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p9
ggsave("fig/analysis/analysis_grid/LFQ-paired-BM.pdf", p9,
       width = 3.47, height = 3.47, units = 'cm', useDingbats = FALSE)

# read GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup() %>%
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile') %>% 
           as.character())

# filter to GO terms that pass breadth cutoff in all quant. modes
filtered = dat2 %>%
  # ignore MS1 intensity
  filter(quant_mode != 'MS1_intensity') %>%
  # filter by breadth
  filter(between(n_chromatograms, 10, 100)) %>%
  # keep the intersect of all three quant. modes
  dplyr::select(-n_proteins, -n_chromatograms, -auroc) %>%
  group_by_at(vars(-go_term, -quant_mode)) %>%
  filter(n_distinct(quant_mode) == 3) %>%
  ungroup()
# filter by breadth, and get median
means = dat2 %>%
  inner_join(filtered) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# run paired tests
x = filter(means, quant_mode == 'iBAQ')
y = filter(means, quant_mode == 'LFQ')
z = filter(means, quant_mode == 'spectral_counts')
# get delta-medians
delta1 = left_join(x, y, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
delta2 = left_join(x, z, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
delta3 = left_join(y, z, by = c('accession', 'experiment', 'analysis',
                                'metric', 'missing', 'transform')) %$%
  median(auroc.y - auroc.x, na.rm = TRUE)
# paired BM test
pval1 = bind_rows(x, y) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .,
                     nperm = 1, rounds = 1, plot.simci = FALSE) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
pval2 = bind_rows(x, z) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
pval3 = bind_rows(y, z) %>% 
  arrange(accession, experiment, metric, missing, transform) %>% 
  group_by(accession, experiment, metric, missing, transform) %>% 
  filter(n() == 2) %>% 
  ungroup() %>%
  npar.t.test.paired(auroc ~ quant_mode, data = .) %>%
  extract2('Analysis') %>% 
  extract('BM', 'p.value')
# double-check t-tests
left_join(x, y, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)
left_join(x, z, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)
left_join(y, z, by = c('accession', 'experiment', 'analysis',
                       'metric', 'missing', 'transform')) %$%
  t.test(auroc.y - auroc.x)

# create data frame and plot
df10 = data.frame(quant_mode1 = c('iBAQ', 'iBAQ', 'MaxLFQ'),
                  quant_mode2 = c('MaxLFQ', 'spectral count', 'spectral count'),
                  delta = c(delta1, delta2, delta3),
                  pval = c(pval1, pval2, pval3)) %>% 
  mutate(sig = pval < 0.05) 
## mirror
df10 %<>% bind_rows(df10 %>% 
                      mutate(tmp = quant_mode2,
                             quant_mode2 = quant_mode1,
                             quant_mode1 = tmp,
                             delta = -delta))
range = range(df10$delta)
lvls = c('MaxLFQ', 'spectral count', 'iBAQ') ## as in Fig. 2j
p10 = df10 %>% 
  ggplot(aes(x = factor(quant_mode1, levels = lvls),
             y = factor(quant_mode2, levels = lvls))) +
  geom_tile(aes(fill = delta), color = 'white') +
  geom_text(data = filter(df10, sig), aes(label = '*'), size = 1.75,
            color = 'white', nudge_y = -0.1) +
  scale_x_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_y_discrete(expand = c(0, 0), labels = clean_metric) +
  scale_fill_paletteer_c("pals::coolwarm",
                         name = expression(Delta~AUC~ "   "),
                         limits = range, breaks = range,
                         labels = round(range, digits = 3)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.2, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p10
ggsave("fig/analysis/analysis_grid/LFQ-paired-BM-GO.pdf", p10,
       width = 3.47, height = 3.47, units = 'cm', useDingbats = FALSE)

# write tables
library(openxlsx)
tables = list('3a. Measures of association' = df1 %>% 
                mutate_at(vars(metric1, metric2), clean_metric) %>% 
                set_colnames(c("Metric 1", "Metric 2", "Difference in medians", "P value")),
              '3b. Measures of association' = df2 %>% 
                mutate_at(vars(metric1, metric2), clean_metric) %>% 
                set_colnames(c("Metric 1", "Metric 2", "Difference in medians", "P value")),
              '3c. Missing values' = df7 %>% 
                dplyr::select(missing1, missing2, delta, pval) %>% 
                set_colnames(c("Missing values 1", "Missing values 2", "Difference in medians", "P value")),
              '3d. Missing values' = df8 %>% 
                dplyr::select(missing1, missing2, delta, pval) %>% 
                set_colnames(c("Missing values 1", "Missing values 2", "Difference in medians", "P value")),
              '3e. Normalization' = df5 %>% 
                dplyr::select(transformation1, transformation2, delta, pval) %>% 
                set_colnames(c("Normalization 1", "Normalization 2", "Difference in medians", "P value")),
              '3f. Normalization' = df6 %>% 
                dplyr::select(transformation1, transformation2, delta, pval) %>% 
                set_colnames(c("Normalization 1", "Normalization 2", "Difference in medians", "P value")),
              '3g. Measures of association' = df3 %>% 
                mutate_at(vars(metric1, metric2), clean_metric) %>% 
                set_colnames(c("Metric 1", "Metric 2", "Difference in medians", "P value")),
              '3h. Measures of association' = df4 %>% 
                mutate_at(vars(metric1, metric2), clean_metric) %>% 
                set_colnames(c("Metric 1", "Metric 2", "Difference in medians", "P value")),
              '3i. Label-free quantification' = df9 %>% 
                dplyr::select(quant_mode1, quant_mode2, delta, pval) %>% 
                set_colnames(c("Quant. mode 1", "Quant. mode 2", "Difference in medians", "P value")),
              '3j. Label-free quantification' = df10 %>% 
                dplyr::select(quant_mode1, quant_mode2, delta, pval) %>% 
                set_colnames(c("Quant. mode 1", "Quant. mode 2", "Difference in medians", "P value"))
              ) 
write.xlsx(tables, "data/analysis/analysis_grid/univariate-analysis.xlsx")
