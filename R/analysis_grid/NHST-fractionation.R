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

# merge in fractionation data and subset
expts = read.csv("~/git/PCPdb/data/experiments.csv") %>% 
  dplyr::rename(accession = Accession, experiment = Replicate,
                fractionation = Fractionation) %>% 
  dplyr::select(accession, experiment, fractionation)
input1 = dat1 %>% 
  left_join(expts, by = c('accession', 'experiment')) %>% 
  filter(fractionation != 'XL-SEC') %>% 
  distinct()

# fractionation
input1 %>% 
  group_by(fractionation) %>% 
  summarise(median = median(auroc))
fractionations = with(input1, reorder(fractionation, auroc, stats::median)) %>%
  levels()
pairs = tidyr::crossing(fractionation1 = fractionations,
                        fractionation2 = fractionations) %>% 
  filter(fractionation1 != fractionation2)
df1 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(input1, fractionation %in% c(current$fractionation1, 
                                           current$fractionation2))
  x = filter(df, fractionation == current$fractionation1) %>% pull(auroc)
  y = filter(df, fractionation == current$fractionation2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(fractionation1 = factor(fractionation1, levels = fractionations),
         fractionation2 = factor(fractionation2, levels = fractionations),
         padj = p.adjust(pval, 'BH'),
         label1 = ifelse(pval < 0.001, '***', 
                         ifelse(pval < 0.01, '**',
                                ifelse(pval < 0.05, '*', '')))) 
sig1 = df1 %>% 
  filter(pval < 0.05)
range = range(df1$delta)
p1 = df1 %>% 
  mutate(delta = winsorize(delta, range)) %>% 
  ggplot(aes(x = fractionation1, y = fractionation2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  # geom_text(aes(label = label1), size = 1.75) +
  geom_text(data = sig1, aes(label = '*'), size = 1.75, color = 'white', 
            nudge_y = -0.1) +
  # geom_tile(data = sig1, color = 'black', fill = NA) +
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
ggsave("fig/analysis/analysis_grid/fractionation-BM.pdf", p1,
       width = 3.25, height = 3.25, units = 'cm', useDingbats = FALSE)

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
input2 = dat2 %>% 
  left_join(expts, by = c('accession', 'experiment')) %>% 
  filter(fractionation != 'XL-SEC',
         fractionation != 'XL-N-PAGE',
         fractionation != 'Sucrose') %>% 
  distinct()

# fractionation
input2 %>% 
  group_by(fractionation) %>% 
  summarise(median = median(auroc))
fractionations = with(input2, reorder(fractionation, auroc, stats::median)) %>%
  levels()
df2 = pmap_dfr(pairs, function(...) {
  current = tibble(...)
  # print(current)
  df = filter(input2, fractionation %in% c(current$fractionation1, 
                                           current$fractionation2))
  x = filter(df, fractionation == current$fractionation1) %>% pull(auroc)
  y = filter(df, fractionation == current$fractionation2) %>% pull(auroc)
  delta = median(y) - median(x)
  pval = lawstat::brunner.munzel.test(x, y)$p.value
  cbind(current, delta = delta, pval = pval)
}) %>% 
  mutate(fractionation1 = factor(fractionation1, levels = fractionations),
         fractionation2 = factor(fractionation2, levels = fractionations),
         padj = p.adjust(pval, 'BH'),
         label1 = ifelse(pval < 0.001, '***', 
                         ifelse(pval < 0.01, '**',
                                ifelse(pval < 0.05, '*', '')))) 
sig2 = df2 %>% 
  filter(pval < 0.05)
range = range(df2$delta)
p2 = df2 %>% 
  mutate(delta = winsorize(delta, range)) %>% 
  ggplot(aes(x = fractionation1, y = fractionation2)) +
  geom_tile(aes(fill = delta), color = 'white') +
  # geom_text(aes(label = label1), size = 1.75) +
  geom_text(data = sig2, aes(label = '*'), size = 1.75, color = 'white', 
            nudge_y = -0.1) +
  # geom_tile(data = sig1, color = 'black', fill = NA) +
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
ggsave("fig/analysis/analysis_grid/fractionation-BM-GO.pdf", p2,
       width = 3.25, height = 3.25, units = 'cm', useDingbats = FALSE)

# multivariate analysis
# fit a linear model with interaction terms
fractionations = with(input1, reorder(fractionation, auroc, stats::median)) %>%
  levels()
## preprocessing is not applicable to co-occurrence
co_occurrence = c('hamming', 'binomial', 'dice', 'jaccard')
metrics = with(input1, reorder(metric, auroc, stats::median)) %>% levels()
fit1 = input1 %>% 
  mutate(fractionation = factor(fractionation, levels = fractionations),
         missing = ifelse(metric %in% co_occurrence, 'co-occurrence', missing),
         transform = ifelse(metric %in% co_occurrence, 'co-occurrence', 
                            transform),
         metric = factor(metric, levels = metrics),
         missing = factor(missing, levels = c('NA', 'zero', 'noise',
                                              'co-occurrence')),
         transform = factor(transform, levels = c('quantile normalization',
                                                  'log-transform', 'none'))) %>% 
  do(tidy(lm(auroc ~ metric * missing + metric * transform + fractionation, 
             data = .))) %>% 
  as.data.frame() %>% 
  filter(term != '(Intercept)') %>% 
  # clean up terms
  mutate(type = ifelse(grepl(":", term), "Interaction terms", 
                       ifelse(grepl("metric", term), "Measure of association",
                              ifelse(grepl("missing", term), "Missing values",
                                     ifelse(grepl("fractionation", term),
                                            "Fractionation",
                                            "Normalization")))),
         term = gsub("metric|missing|^transform|fractionation", "", term),
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
                                'Quantile normalization (baseline)',
                                'IEX (baseline)'),
                       estimate = c(0, 0, 0, 0),
                       type = c('Measure of association', 'Missing values', 
                                'Normalization', 'Fractionation')))
# plot the coefficients for fractionation only
pal = c('FALSE' = 'grey70', 'TRUE' = 'black')
p3 = fit1 %>% 
  filter(type == 'Fractionation') %>% 
  ggplot(aes(x = reorder(term, estimate), y = estimate, 
             color = sig)) +
  geom_point(size = 0.7) +
  geom_point(data = filter(fit1, is.na(sig), type == 'Fractionation'),
             shape = 1, color = 'black', size = NA) +
  geom_errorbar(aes(ymin = 0, ymax = estimate), width = 0) +
  scale_x_reordered() +
  scale_y_continuous(expression(beta), limits = c(0, 0.059)) +
  scale_color_manual('', values = pal, drop = FALSE, 
                     breaks = c('FALSE', 'TRUE'),
                     labels = c('FALSE' = 'p >= 0.05', 'TRUE' = 'p < 0.05')) +
  coord_flip() + 
  boxed_theme() +
  theme(axis.title.y = element_blank())
p3
ggsave("fig/analysis/analysis_grid/fractionation-lm.pdf", p3,
       width = 3.75, height = 3, ## height = 2.75, 
       units = "cm", useDingbats = FALSE)

# repeat for GO
fractionations = with(input2, reorder(fractionation, auroc, stats::median)) %>%
  levels()
metrics = with(input2, reorder(metric, auroc, stats::median)) %>% levels()
fit2 = input2 %>% 
  mutate(fractionation = factor(fractionation, levels = fractionations),
         missing = ifelse(metric %in% co_occurrence, 'co-occurrence', missing),
         transform = ifelse(metric %in% co_occurrence, 'co-occurrence', 
                            transform),
         metric = factor(metric, levels = metrics),
         missing = factor(missing, levels = c('NA', 'zero', 'noise',
                                              'co-occurrence')),
         transform = factor(transform, levels = c('quantile normalization',
                                                  'log-transform', 'none'))) %>% 
  do(tidy(lm(auroc ~ metric * missing + metric * transform + fractionation, 
             data = .))) %>% 
  as.data.frame() %>% 
  filter(term != '(Intercept)') %>% 
  # clean up terms
  mutate(type = ifelse(grepl(":", term), "Interaction terms", 
                       ifelse(grepl("metric", term), "Measure of association",
                              ifelse(grepl("missing", term), "Missing values",
                                     ifelse(grepl("fractionation", term),
                                            "Fractionation",
                                            "Normalization")))),
         term = gsub("metric|missing|^transform|fractionation", "", term),
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
                                'Quantile normalization (baseline)',
                                'IEF (baseline)'),
                       estimate = c(0, 0, 0, 0),
                       type = c('Measure of association', 'Missing values', 
                                'Normalization', 'Fractionation')))
# plot fractionation coefficients only
p4 = fit2 %>% 
  filter(type == 'Fractionation') %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, 
             color = sig)) +
  geom_point(size = 0.7) +
  geom_point(data = filter(fit2, is.na(sig), type == 'Fractionation'),
             shape = 1, color = 'black', size = NA) +
  geom_errorbar(aes(ymin = 0, ymax = estimate), width = 0) +
  scale_x_reordered() +
  scale_y_continuous(expression(beta), limits = c(0, 0.014),
                     breaks = seq(0, 0.012, 0.006)) +
  scale_color_manual('', values = pal, drop = FALSE, 
                     breaks = c('FALSE', 'TRUE'),
                     labels = c('FALSE' = 'p >= 0.05', 'TRUE' = 'p < 0.05')) +
  coord_flip() + 
  boxed_theme() +
  theme(axis.title.y = element_blank())
p4
ggsave("fig/analysis/analysis_grid/fractionation-lm-GO.pdf", p4,
       width = 3.75, height = 3, ## height = 2.75, 
       units = "cm", useDingbats = FALSE)

# write tables
library(openxlsx)
tables = list(
  df1 %>% 
    dplyr::select(fractionation1, fractionation2, delta, pval) %>% 
    set_colnames(c("Fractionation 1", "Fractionation 2", "Difference in medians", "P value")),
  df2 %>% 
    dplyr::select(fractionation1, fractionation2, delta, pval) %>% 
    set_colnames(c("Fractionation 1", "Fractionation 2", "Difference in medians", "P value")),
  fit1 %>% 
    dplyr::select(-sig) %>% 
    set_colnames(., c('Term', 'Estimate', 'Standard error', 'Statistic',
                      'P value', 'Type')) %>% 
    filter(!grepl("baseline", Term)),
  fit2 %>% 
    dplyr::select(-sig) %>% 
    set_colnames(., c('Term', 'Estimate', 'Standard error', 'Statistic',
                      'P value', 'Type')) %>% 
    filter(!grepl("baseline", Term))
) %>% 
  setNames(c('5a. Univariate analysis, complexes',
             '5b. Univariate analysis, GO',
             '5c. Multivariable analysis, complexes',
             '5d. Multivariable analysis, GO'))
write.xlsx(tables, "data/analysis/analysis_grid/fractionation.xlsx")

# how many coefficients in a four way model?
fit = lm(auroc ~ metric*transform*missing*fractionation, data=input2)
