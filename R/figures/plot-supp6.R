# Plot Supplementary Figure 6:
#' - split_by
#' - feature_select
#' - combine_features
#' - replace_missing_data
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read AUCs from the replicate_integration experiment
auc = readRDS("data/analysis/replicate_integration/AUC.rds") %>%
  # drop some (constant) columns
  dplyr::select(-min_fractions, -n_folds) %>%
  # pre-filter 1 dataset
  filter(n_datasets > 1)

################################################################################
###### a-b. split_by
################################################################################

# calculate delta AUC between split_by=proteins and pairs
delta_split = auc %>%
  mutate(split_by = fct_relevel(split_by, 'proteins', 'pairs')) %>%
  group_by_at(vars(-split_by, -auc, -n_1, -n_0, -n_proteins, -n_pairs)) %>%
  filter(n_distinct(split_by) == 2) %>%
  arrange(split_by) %>%
  summarise(delta = diff(auc), n = n()) %>%
  ungroup()

# filter that grid
df1 = filter(delta_split,
             classifier == 'NB',
             between(n_datasets, 2, 4),
             between(n_features, 2, 6),
             replace_missing_data == TRUE,
             combine_features == FALSE)

# plot
pal = brewer.pal(12, 'Set3')
cols = c('grey75', pal[4])

pB = df1 %>%
  mutate(n_datasets = fct_recode(as.character(n_datasets),
                                 '1 dataset' = '1',
                                 '2 datasets' = '2',
                                 '3 datasets' = '3',
                                 '4 datasets' = '4'),
         feature_select = fct_recode(feature_select,
                                     'Features:\nbest first' = 'best_first',
                                     'Features:\nrandom' = 'random'),
         complex_set = fct_recode(complex_set,
                                  'cross-validation' = 'held-in')) %>%
  ggplot(aes(x = factor(n_features), y = delta, fill = complex_set, 
             color = complex_set)) +
  facet_grid(feature_select ~ n_datasets) +
  geom_hline(aes(yintercept = 0), color = 'grey86') +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA,
               outlier.size = 0.3) + 
  scale_x_discrete('# of features') +
  scale_y_continuous(expression(AUC[pairs]~-~AUC[proteins])) +
  scale_color_manual('Evaluation', values = cols, 
                     labels = Hmisc::capitalize) +
  scale_fill_manual('Evaluation', values = cols,
                    labels = Hmisc::capitalize) +
  guides(color = guide_legend(title.pos = 'top', title.hjust = 0.5)) +
  coord_cartesian(ylim = c(-0.01, 0.05)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.45, 'lines'))
pB
ggsave("fig/final/supp6/b_split_by.pdf", pB,
       width = 11.5, height = 6.75, units = "cm", useDingbats = FALSE)

################################################################################
###### c-d. feature_select
################################################################################

# filter the grid
df2 = auc %>%
  filter(between(n_datasets, 2, 4),
         n_features <= 6,
         replace_missing_data == TRUE,
         combine_features == FALSE,
         split_by == 'proteins',
         complex_set == 'held-out')

# plot
pal = brewer.pal(10, 'Set3')
cols = c(pal[3], 'grey75', pal[6])
pC = df2 %>%
  mutate(n_datasets = fct_recode(as.character(n_datasets),
                                 '1 dataset' = '1',
                                 '2 datasets' = '2',
                                 '3 datasets' = '3',
                                 '4 datasets' = '4'),
         feature_select = fct_recode(feature_select,
                                     'best first' = 'best_first'),
         feature_select = fct_relevel(feature_select, 'best first',
                                      'random', 'PrInCE'),
         classifier = fct_recode(classifier,
                                 'Naive Bayes' = 'NB',
                                 'Random forest' = 'RF')) %>%
  mutate(xval = ifelse(feature_select == 'PrInCE', 5, n_features)) %>%
  ggplot(aes(x = factor(xval), y = auc, fill = feature_select, 
             color = feature_select)) +
  facet_grid(classifier ~ n_datasets) +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA,
               outlier.size = 0.3) + 
  scale_x_discrete('# of features') +
  scale_y_continuous('AUC') +
  scale_color_manual('Feature selection', values = cols,
                     labels = Hmisc::capitalize) +
  scale_fill_manual('Feature selection', values = cols,
                    labels = Hmisc::capitalize) +
  guides(color = guide_legend(title.pos = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(legend.position = 'top', # c(1, 0),
        # legend.justification = c(1, 0),
        legend.key.size = unit(0.45, 'lines'),
        aspect.ratio = 1)
pC
ggsave("fig/final/supp6/c_feature_select.pdf", pC,
       width = 11.5, height = 6.75, units = "cm", useDingbats = FALSE)

delta_feature_select = auc %>%
  filter(feature_select != 'PrInCE') %>%
  group_by_at(vars(-auc, -n_1, -n_0, -n_proteins, -n_pairs,
                   -feature_select)) %>%
  summarise(delta = auc[feature_select == 'best_first'] - 
              auc[feature_select == 'random'],
            n = n()) %>%
  ungroup()
df3 = delta_feature_select %>%
  filter(n_datasets <= 10,
         n_features <= 10,
         split_by == 'proteins',
         replace_missing_data,
         !combine_features) %>%
  mutate(classifier = fct_recode(classifier,
                                 'Naive Bayes' = 'NB',
                                 'Random forest' = 'RF'))
pD = df3 %>%
  # average over ten samples
  group_by_at(vars(-sample_idx, -delta)) %>%
  summarise(n = n(), auc = mean(delta)) %>%
  ungroup() %>%
  ggplot(aes(y = factor(n_features), x = factor(n_datasets))) +
  facet_grid(classifier ~ ., space = 'free', scales = 'free') +
  geom_tile(aes(fill = auc), color = 'white') +
  scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22", 
                         name = expression(AUC[best~first]~-~AUC[random]),
                         limits = c(-0.138, 0.138),
                         breaks = c(-0.12, 0.12)
  ) +
  scale_y_discrete('# of features', expand = c(0, 0)) +
  scale_x_discrete('# of datasets', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  # coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1,
        # legend.justification = 'middle',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
pD
ggsave("fig/final/supp6/d_delta_feature_select.pdf", pD,
       width = 6.4, height = 6.1, units = 'cm', useDingbats = FALSE)

################################################################################
###### e. combine_features
################################################################################

delta_combine = auc %>%
  group_by_at(vars(-auc, -n_1, -n_0, -n_proteins, -n_pairs,
                   -combine_features)) %>%
  summarise(delta = auc[combine_features] - auc[!combine_features],
            n = n()) %>%
  ungroup()
df4 = delta_combine %>%
  filter(feature_select != 'PrInCE',
         between(n_datasets, 2, 4),
         n_features <= 6,
         split_by == 'proteins',
         complex_set == 'held-out',
         replace_missing_data) %>%
  mutate(classifier = fct_recode(classifier,
                                 'Naive Bayes' = 'NB',
                                 'Random forest' = 'RF'))

pal = brewer.pal(4, 'Blues')[-1]
pE = df4 %>%
  mutate(feature_select = fct_recode(feature_select,
                                     'Features:\nbest first' = 'best_first',
                                     'Features:\nrandom' = 'random')) %>%
  ggplot(aes(x = factor(n_features), y = delta, fill = factor(n_datasets), 
             color = factor(n_datasets))) +
  facet_grid(classifier ~ feature_select) +
  geom_hline(aes(yintercept = 0), color = 'grey86') +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA, 
               outlier.size = 0.3) +
  scale_x_discrete('# of features') +
  scale_y_continuous(expression(AUC[combined~matrix]~-~AUC[separate~matrices])) +
  scale_color_manual('# of datasets', values = pal) +
  scale_fill_manual('# of datasets', values = pal) +
  guides(color = guide_legend(title.pos = 'top', title.hjust = 0.5)) +
  coord_cartesian(ylim = c(-0.2, 0.05)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.45, 'lines'))
pE
ggsave("fig/final/supp6/e_combine_features.pdf", pE,
       width = 9, height = 6.98, units = "cm", useDingbats = FALSE)

################################################################################
###### f. replace_missing_data
################################################################################

delta_missing = auc %>%
  # only NB can handle missing data
  filter(classifier == 'NB') %>%
  group_by_at(vars(-auc, -n_1, -n_0, -n_proteins, -n_pairs,
                   -replace_missing_data)) %>%
  summarise(delta = auc[replace_missing_data] - auc[!replace_missing_data],
            n = n()) %>%
  ungroup()

df5 = delta_missing %>%
  filter(feature_select != 'PrInCE',
         between(n_datasets, 2, 4),
         n_features <= 6,
         split_by == 'proteins',
         complex_set == 'held-out',
         !combine_features)

pal = brewer.pal(4, 'Blues')[-1]
pF = df5 %>%
  mutate(feature_select = fct_recode(feature_select,
                                     'Features:\nbest first' = 'best_first',
                                     'Features:\nrandom' = 'random')) %>%
  ggplot(aes(x = factor(n_features), y = delta, fill = factor(n_datasets), 
             color = factor(n_datasets))) +
  facet_grid(~ feature_select, scales = 'free') +
  geom_hline(aes(yintercept = 0), color = 'grey86') +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA, 
               outlier.size = 0.3) +
  scale_x_discrete('# of features') +
  scale_y_continuous(expression(AUC[replace~NAs]~-~AUC[keep~NAs])) +
  scale_color_manual('# of datasets', values = pal) +
  scale_fill_manual('# of datasets', values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.45, 'lines'))
pF
ggsave("fig/final/supp6/f_replace_missing_data.pdf", pF,
       width = 9, height = 5, units = "cm", useDingbats = FALSE)

################################################################################
###### g. n_features
################################################################################

df6 = auc %>%
  filter(feature_select != 'PrInCE',
         replace_missing_data,
         !combine_features,
         split_by == 'proteins',
         complex_set == 'held-out',
         n_datasets <= 10,
         n_features <= 10)
  
pG = df6 %>% 
  group_by_at(vars(-sample_idx, -auc, -n_1, -n_0, -n_proteins, -n_pairs)) %>%
  summarise(n = n(), auc = mean(auc)) %>%
  mutate(feature_select = fct_recode(feature_select,
                                     'Features:\nbest first' = 'best_first',
                                     'Features:\nrandom' = 'random'),
         classifier = fct_recode(classifier,
                                 'Random forest' = 'RF',
                                 'Naive Bayes' = 'NB')) %>%
  ggplot(aes(y = factor(n_features), x = factor(n_datasets), fill = auc)) +
  facet_grid(feature_select ~ classifier) +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC', breaks = c(0.6, 0.8)) +
  scale_y_discrete('# of features', expand = c(0, 0)) +
  scale_x_discrete('# of datasets', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
pG
ggsave("fig/final/supp6/g_n_features_n_datasets.pdf", pG,
       width = 6.48, height = 7, units = "cm", useDingbats = FALSE)

################################################################################
###### h. classifier (delta)
################################################################################

delta_clf = auc %>%
  group_by_at(vars(-auc, -n_1, -n_0, -n_proteins, -n_pairs, -classifier)) %>%
  summarise(delta = auc[classifier == 'RF'] - auc[classifier == 'NB'],
            n = n()) %>%
  ungroup()
df7 = delta_clf %>%
  filter(n_datasets <= 10,
         n_features <= 10,
         split_by == 'proteins',
         !combine_features,
         replace_missing_data,
         feature_select != 'PrInCE') %>%
  mutate(feature_select = fct_recode(feature_select,
                                     'Features:\nbest first' = 'best_first',
                                     'Features:\nrandom' = 'random'))

pH = df7 %>%
  # average over ten samples
  group_by_at(vars(-sample_idx, -delta)) %>%
  summarise(n = n(), auc = mean(delta)) %>%
  ungroup() %>%
  ggplot(aes(y = factor(n_features), x = factor(n_datasets))) +
  facet_grid(feature_select ~ .) +
  geom_tile(aes(fill = auc), color = 'white') +
  scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22", 
                         name = expression(AUC[RF]~-~AUC[NB]),
                         limits = c(-0.16, 0.16),
                         breaks = c(-0.16, 0.16)) +
  scale_y_discrete('# of features', expand = c(0, 0)) +
  scale_x_discrete('# of datasets', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
pH
ggsave("fig/final/supp6/h_delta_classifier.pdf", pH,
       width = 7, height = 5.31, units = "cm", useDingbats = FALSE)
