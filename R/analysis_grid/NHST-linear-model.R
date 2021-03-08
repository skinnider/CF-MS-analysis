setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
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

# fit a linear model to identify significant preprocessing interactions
## use worst-performing metric as the baseline
metrics = with(dat1, reorder(metric, auroc, stats::median)) %>% levels()
## note that preprocessing is not applicable to co-occurrence
co_occurrence = c('hamming', 'binomial', 'dice', 'jaccard')
fit1 = dat1 %>% 
  mutate(missing = ifelse(metric %in% co_occurrence, 'co-occurrence', missing),
         transform = ifelse(metric %in% co_occurrence, 'co-occurrence', 
                            transform),
         missing = factor(missing, levels = c('NA', 'zero', 'noise',
                                              'co-occurrence')),
         metric = factor(metric, levels = metrics),
         transform = factor(transform, levels = c('quantile normalization',
                                                  'log-transform', 'none'))) %>% 
  distinct() %>% 
  do(tidy(lm(auroc ~ metric * missing + metric * transform, data = .))) %>% 
  as.data.frame() %>% 
  filter(term != '(Intercept)') %>% 
  # clean up terms
  mutate(type = ifelse(grepl(":", term), "Interaction terms", 
                       ifelse(grepl("metric", term), "Measure of association",
                              ifelse(grepl("missing", term), "Missing values",
                                     "Normalization"))),
         term = gsub("metric|missing|^transform", "", term),
         term = gsub(":transform", ":", term),
         term = gsub(":", " : ", term),
         term = ifelse(type == 'Measure of association', clean_metric(term),
                       ifelse(type == 'Interaction terms',
                              paste0(gsub(" :.*$", "", term) %>% clean_metric,
                                     " : ",
                                     gsub("^.*: ", "", term) %>% Hmisc::capitalize()),
                              Hmisc::capitalize(term))),
         sig = p.value < 0.05) %>% 
  # add the baselines back in 
  bind_rows(data.frame(term = c('treeClust (baseline)', 'NA (baseline)', 
                                'Quantile normalization (baseline)'),
                       estimate = c(0, 0, 0),
                       type = c('Measure of association', 'Missing values', 
                                'Normalization'))) %>% 
  mutate(type = fct_relevel(type, 
                            'Interaction terms',
                            'Normalization',
                            'Missing values',
                            'Measure of association'),
         term = gsub("metric|missing|transformation", "", term))

# plot
pal = c('FALSE' = 'grey70', 'TRUE' = 'black')
p1 = fit1 %>% 
  ggplot(aes(x = reorder_within(term, estimate, type), y = estimate, 
             color = sig)) +
  facet_grid(type ~ ., space = 'free', scales = 'free') +
  geom_point(size = 0.7) +
  geom_point(data = filter(fit1, is.na(sig)), shape = 1, color = 'black',
             size = NA) +
  geom_errorbar(aes(ymin = 0, ymax = estimate), width = 0) +
  scale_x_reordered() +
  scale_y_continuous(expression(beta)) +
  scale_color_manual('', values = pal, breaks = c('FALSE', 'TRUE'),
                     labels = c('FALSE' = 'p >= 0.05', 'TRUE' = 'p < 0.05')) +
  coord_flip() + 
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1
ggsave("fig/analysis/analysis_grid/linear-model.pdf", p1, 
       width = 9, height = 22, units = "cm", useDingbats = FALSE)

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

# fit linear model
metrics = with(dat2, reorder(metric, auroc, stats::median)) %>% levels()
fit2 = dat2 %>% 
  mutate(missing = ifelse(metric %in% co_occurrence, 'co-occurrence', missing),
         transform = ifelse(metric %in% co_occurrence, 'co-occurrence', 
                            transform),
         missing = factor(missing, levels = c('NA', 'zero', 'noise',
                                              'co-occurrence')),
         metric = factor(metric, levels = metrics),
         transform = factor(transform, levels = c('quantile normalization',
                                                  'log-transform', 'none'))) %>% 
  distinct() %>% 
  do(tidy(lm(auroc ~ metric * missing + metric * transform, data = .))) %>% 
  as.data.frame() %>% 
  filter(term != '(Intercept)') %>% 
  # clean up terms
  mutate(type = ifelse(grepl(":", term), "Interaction terms", 
                       ifelse(grepl("metric", term), "Measure of association",
                              ifelse(grepl("missing", term), "Missing values",
                                     "Normalization"))),
         term = gsub("metric|missing|^transform", "", term),
         term = gsub(":transform", ":", term),
         term = gsub(":", " : ", term),
         term = ifelse(type == 'Measure of association', clean_metric(term),
                       ifelse(type == 'Interaction terms',
                              paste0(gsub(" :.*$", "", term) %>% clean_metric,
                                     " : ",
                                     gsub("^.*: ", "", term) %>% Hmisc::capitalize()),
                              Hmisc::capitalize(term))),
         sig = p.value < 0.05) %>% 
  # add the baselines back in 
  bind_rows(data.frame(term = c('treeClust (baseline)', 'NA (baseline)', 
                                'Quantile normalization (baseline)'),
                       estimate = c(0, 0, 0),
                       type = c('Measure of association', 'Missing values', 
                                'Normalization'))) %>% 
  mutate(type = fct_relevel(type, 
                            'Interaction terms',
                            'Normalization',
                            'Missing values',
                            'Measure of association'),
         term = gsub("metric|missing|transformation", "", term))

# plot
p2 = fit2 %>% 
  ggplot(aes(x = reorder_within(term, estimate, type), y = estimate, 
             color = sig)) +
  facet_grid(type ~ ., space = 'free', scales = 'free') +
  geom_point(size = 0.7) +
  geom_point(data = filter(fit1, is.na(sig)), shape = 1, color = 'black',
             size = NA) +
  geom_errorbar(aes(ymin = 0, ymax = estimate), width = 0) +
  scale_x_reordered() +
  scale_y_continuous(expression(beta)) +
  scale_color_manual('', values = pal, breaks = c('FALSE', 'TRUE'),
                     labels = c('FALSE' = 'p >= 0.05', 'TRUE' = 'p < 0.05')) +
  coord_flip() + 
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2
ggsave("fig/analysis/analysis_grid/linear-model-GO.pdf", p2, 
       width = 9, height = 22, units = "cm", useDingbats = FALSE)

# save coefficients
library(openxlsx)
tables = list('4a. Complexes' = fit1,
              '4b. GO' = fit2) %>% 
  map(~ set_colnames(., c('Term', 'Estimate', 'Standard error', 'Statistic',
                          'P value', 'Type', 'sig')) %>% 
        dplyr::select(-sig) %>% 
        filter(!grepl("baseline", Term)))
write.xlsx(tables, "data/analysis/analysis_grid/linear-model.xlsx")
