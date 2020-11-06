# Plot Figure 4.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read AUCs from the replicate_integration experiment
auc = readRDS("data/analysis/replicate_integration/AUC.rds") %>%
  # drop some (constant) columns
  dplyr::select(-min_fractions, -n_folds)

################################################################################
###### a. split_by
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
             n_datasets == 3,
             replace_missing_data == TRUE,
             feature_select == 'best_first',
             combine_features == FALSE)

# plot
pal = brewer.pal(12, 'Set3')
cols = c('grey75', pal[4])
pA = df1 %>%
  mutate(complex_set = fct_recode(complex_set,
                                  'cross-validation' = 'held-in')) %>%
  ggplot(aes(x = factor(n_features), y = delta, fill = complex_set, 
             color = complex_set)) +
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
  coord_cartesian(ylim = c(-0.005, 0.025)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.45, 'lines'))
pA

################################################################################
###### b. feature_select
################################################################################

# filter the grid
df2 = auc %>%
  filter(classifier == 'RF',
         n_datasets == 3,
         n_features <= 6,
         replace_missing_data == TRUE,
         combine_features == FALSE,
         split_by == 'proteins',
         complex_set == 'held-out')

# plot
pal = brewer.pal(10, 'Set3')
cols = c(pal[3], 'grey75', pal[6])
pB = df2 %>%
  mutate(feature_select = fct_recode(feature_select,
                                     'best first' = 'best_first'),
         feature_select = fct_relevel(feature_select, 'best first',
                                      'random', 'PrInCE'),
         classifier = fct_recode(classifier,
                                 'Naive Bayes' = 'NB',
                                 'Random forest' = 'RF')) %>%
  mutate(xval = ifelse(feature_select == 'PrInCE', 5, n_features)) %>%
  ggplot(aes(x = factor(xval), y = auc, fill = feature_select, 
             color = feature_select)) +
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
pB

################################################################################
###### Combine and plot row #1
################################################################################

row1 = pA + pB + plot_layout(nrow = 1)
ggsave("fig/final/fig4/row1.pdf", row1, width = 8.75, height = 10, 
       units = "cm", useDingbats = FALSE)

################################################################################
###### c. AUC, n_features x n_datasets (heatmap)
################################################################################

# filter the grid
df3 = filter(auc,
             feature_select == 'best_first',
             classifier == 'RF',
             between(n_datasets, 2, 10),
             n_features <= 10,
             replace_missing_data,
             !combine_features,
             split_by == 'proteins',
             complex_set == 'held-out') %>%
  # calculate mean over ten samples
  group_by_at(vars(-sample_idx, -auc, -n_1, -n_0, -n_proteins, -n_pairs)) %>%
  summarise(n = n(), auc = mean(auc)) %>%
  ungroup()

# plot
range = range(df3$auc)
pC = df3 %>%
  ggplot(aes(y = factor(n_features), x = factor(n_datasets), fill = auc)) +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC', 
                         limits = c(0.66, range[2]),
                         breaks = c(0.66, 0.82),
                         labels = format(range, format = 'f', digits = 2)) +
  scale_y_discrete('# of features', expand = c(0, 0)) +
  scale_x_discrete('# of datasets', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'top',
        legend.justification = 'center',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
pC

################################################################################
###### d. delta AUC, classifier, n_features x n_datasets (heatmap)
################################################################################

# calculate delta
delta_clf = auc %>%
  group_by_at(vars(-auc, -n_1, -n_0, -n_proteins, -n_pairs, -classifier)) %>%
  filter(n_distinct(classifier) == 2) %>%
  summarise(delta = auc[classifier == 'RF'] - auc[classifier == 'NB'],
            n = n()) %>%
  ungroup()

# filter the grid
df4 = delta_clf %>%
  filter(between(n_datasets, 2, 10),
         n_features <= 10,
         split_by == 'proteins',
         !combine_features,
         replace_missing_data,
         feature_select == 'best_first',
         complex_set == 'held-out') %>%
  # calculate mean over ten samples
  group_by_at(vars(-sample_idx, -delta)) %>%
  summarise(n = n(), auc = mean(delta)) %>%
  ungroup()

# plot
range = range(df4$auc)
limits = c(-max(abs(range)), max(abs(range)))
pD = df4 %>%
  ggplot(aes(y = factor(n_features), x = factor(n_datasets))) +
  geom_tile(aes(fill = auc), color = 'white') +
  scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22", 
                         name = expression(AUC[RF]~-~AUC[NB]),
                         limits = limits,
                         breaks = limits,
                         labels = round(limits, digits = 2)) +
  scale_y_discrete('# of features', expand = c(0, 0)) +
  scale_x_discrete('# of datasets', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                               title.position = 'top', title.hjust = 0.5)) +
  coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'top',
        legend.justification = 'center',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
pD

################################################################################
###### e. n_datasets
################################################################################

# subset the grid
df5 = auc %>%
  filter(complex_set == 'held-out',
         split_by == 'proteins',
         n_features == 1,
         n_datasets > 1,
         feature_select == 'best_first',
         replace_missing_data,
         !combine_features)

# plot
pal = brewer.pal(12, 'Set3') 
pal0 = pal[-c(3, 4, 6, 5, 9, 2)] # remove used or unusable colors
cols = pal0[c(1, 4)]
pE = df5 %>%
  ggplot(aes(x = factor(n_datasets), y = auc, fill = classifier, 
             color = classifier)) +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA, 
               outlier.size = 0.3) + 
  scale_y_continuous('AUC') +
  scale_x_discrete('# of datasets') +
  scale_color_manual('', values = cols,
                     labels = c('NB' = 'naive Bayes', 'RF' = 'random forest')) +
  scale_fill_manual('', values = cols,
                    labels = c('NB' = 'naive Bayes', 'RF' = 'random forest')) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  boxed_theme() +
  theme(legend.position = c(0.97, 0.03),
        legend.justification = c(1, 0),
        legend.key.size = unit(0.55, 'lines'),
        aspect.ratio = 1)
pE

################################################################################
###### f. downsample_complexes
################################################################################

# read the dataset
downs = readRDS("data/analysis/downsample_corum/AUC.rds") %>%
  # filter to the final grid
  filter(split_by == 'proteins', feature_select == 'best_first') %>%
  # filter to held-out complexes only
  filter(complex_set == 'held-out')

# convert outcomes to long
long = downs %>%
  mutate(corum_pairs = corum_TPs2 + corum_TNs2) %>%
  dplyr::select(split_by, classifier, n_features, n_datasets, downsample_pct,
                sample_idx, auc, corum_proteins2, corum_TPs2, corum_pairs) %>%
  gather('outcome', 'value', corum_proteins2:corum_pairs) %>%
  mutate(outcome = fct_recode(outcome,
                              '# of proteins' = 'corum_proteins2',
                              '# of true-positive PPIs' = 'corum_TPs2',
                              '# of total PPIs' = 'corum_pairs'))

# scale within sample_idx to % of maximum AUC
scaled = long %>%
  group_by_at(vars(-auc, -value, -downsample_pct)) %>%
  mutate(max = max(auc), 
         q99 = quantile(auc, probs = 0.99),
         q95 = quantile(auc, probs = 0.95),
         pct_of_max = auc / max,
         pct_of_q99 = auc / q99,
         pct_of_q95 = auc / q95,
         n = n()) %>%
  ungroup()

# plot
cols = pal0[c(1, 4)]
## highlight low-data regime for naive Bayes
cutoff = scaled %>%
  filter(pct_of_max < 0.98,
         classifier == 'NB',
         n_datasets == 3,
         outcome == '# of proteins') %>%
  arrange(desc(value)) %>%
  head(1)
print(cutoff)
pF = scaled %>%
  filter(n_datasets == 3, 
         outcome == '# of proteins') %>%
  mutate(classifier = fct_recode(classifier,
                                 'naive Bayes' = 'NB',
                                 'random forest' = 'RF')) %>%
  ggplot(aes(x = value, y = pct_of_max, color = classifier)) +
  geom_rect(data = cutoff, 
            aes(xmin = -Inf, xmax = value, ymin = -Inf, ymax = +Inf), 
            alpha = 0.3, fill = 'grey90', color = NA) +
  geom_vline(data = cutoff, aes(xintercept = value), linetype = 'dotted',
             color = 'grey20', size = 0.2) +
  geom_point(size = 0.2) + 
  scale_color_manual('', values = cols) +
  scale_x_continuous('# of complex proteins') +
  scale_y_continuous('% of maximum AUC', labels = function(x) x * 100) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              override.aes = list(size = 0.8))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.5, 'lines'),
        legend.position = c(0.97, 0.03), 
        legend.justification = c(1, 0)
        )
pF

################################################################################
###### Combine and plot row #2
################################################################################

row2 = pC + pD + plot_layout(nrow = 1)
ggsave("fig/final/fig4/row2.pdf", row2, width = 7.25, height = 10, 
       units = "cm", useDingbats = FALSE)

row3 = pE + pF + plot_layout(nrow = 1)
ggsave("fig/final/fig4/row3.pdf", row3, width = 8.9, height = 10, 
       units = "cm", useDingbats = FALSE)
