# Plot the impact of different quantitation strategies on datasets
# generated using isotopic labelling approaches.
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
  # filter to datasets with at least one ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(any(quant_mode == 'ratio')) %>%
  ungroup() %>%
  # flag label
  mutate(label = ifelse(accession == 'PXD006660', 'Dimethyl', 'SILAC'))

# read GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile')) %>%
  # filter to datasets with at least one ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(any(quant_mode == 'ratio')) %>%
  ungroup() %>%
  # flag label
  mutate(label = ifelse(accession == 'PXD006660', 'Dimethyl', 'SILAC'))

# read CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds")

### complexes
# calculate mean AUC per dataset
avg1 = dat1 %>%
  group_by(label, quant_mode, accession, experiment
           # group over metric, transformation, missing
  ) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup()
comp1 = avg1 %>%
  left_join(cdf, by = c('accession', 'experiment', 'quant_mode')) %>%
  filter(n_fractions == 1) %>%
  filter(quant_mode %in% c('iBAQ', 'ratio')) %>%
  distinct(accession, experiment, quant_mode, mean, median, n_proteins)
df1 = filter(comp1, quant_mode == 'iBAQ')
df2 = filter(comp1, quant_mode != 'iBAQ')
wide1 = left_join(df1, df2, by = c('accession', 'experiment'))
delta1 = wide1 %>%
  group_by(accession, experiment) %>%
  mutate(delta_mean = mean.y - mean.x,
         delta_median = median.y - median.x,
         delta_obs = n_proteins.y - n_proteins.x,
         delta_auc = delta_mean) %>%
  ungroup()
summary1 = delta1 %>%
  summarise(sd_mean = sd(delta_mean),
            sd_median = sd(delta_median),
            sd_obs = sd(delta_obs),
            delta_mean = mean(delta_mean),
            delta_median = mean(delta_median),
            delta_obs = mean(delta_obs),
            delta_auc = delta_mean,
            sd_auc = sd_mean) %>%
  mutate(label = paste0('# proteins: ', round(delta_obs), '\nAUC: +',
                        format(delta_median, format = 'f', digits = 2)))
red = Gpal[1]
p1 = delta1 %>%
  ggplot(aes(x = delta_obs, y = delta_auc)) +
  facet_grid(~ summary1$label) +
  geom_hline(aes(yintercept = 0), color = 'grey88', size = 0.65) +
  geom_vline(aes(xintercept = 0), color = 'grey88', size = 0.65) +
  geom_point(size = 0.8, shape = 1) + 
  geom_point(data = summary1, color = red, size = 0.8) +
  geom_errorbar(data = summary1, aes(ymin = delta_auc - sd_auc, 
                                     ymax = delta_auc + sd_auc), 
                width = 0, color = red, size = 0.4) +
  geom_errorbarh(data = summary1, aes(xmin = delta_obs - sd_obs, 
                                      xmax = delta_obs + sd_obs), 
                 width = 0, color = red, size = 0.4) +
  scale_x_continuous(expression(paste(Delta, '(Proteins)')),
                     limits = c(NA, 50)) +
  scale_y_continuous(expression(paste(Delta, '(AUC, SILAC'~-~'iBAQ)'))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        strip.text = element_text(hjust = 1))
p1
ggsave("fig/analysis/analysis_grid/SILAC-vs-iBAQ.pdf", p5a,
       width = 4.5, height = 5, units = 'cm', useDingbats = F)

# repeat for GO
### first: we need to filter to GO terms that pass the breadth filter in *both*
### iBAQ and ratio
filter1 = filter(dat2, quant_mode == 'ratio') %>%
  filter(between(n_chromatograms, 10, 100)) %>%
  dplyr::select(-quant_mode, -n_proteins, -n_chromatograms)
filter2 = filter(dat2, quant_mode == 'iBAQ') %>%
  filter(between(n_chromatograms, 10, 100)) %>%
  dplyr::select(-quant_mode, -n_proteins, -n_chromatograms)
intersect = inner_join(filter1, filter2,
                       by = setdiff(colnames(filter1), 'auroc')) %>%
  dplyr::select(-starts_with('auroc'))
go_means_filtered = dat2 %>%
  inner_join(intersect) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()
# now, calculate mean in each dataset
avg2 = go_means_filtered %>%
  group_by(label, quant_mode, accession, experiment
           # group over metric, transformation, missing
  ) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup()
comp2 = avg2 %>%
  left_join(cdf, by = c('accession', 'experiment', 'quant_mode')) %>%
  filter(n_fractions == 1) %>%
  filter(quant_mode %in% c('iBAQ', 'ratio')) %>%
  distinct(label, accession, experiment, quant_mode, mean, median, n_proteins)
df1 = filter(comp2, quant_mode == 'iBAQ')
df2 = filter(comp2, quant_mode != 'iBAQ')
wide2 = left_join(df1, df2, by = c('accession', 'experiment', 'label'))
delta2 = wide2 %>%
  group_by(label, accession, experiment) %>%
  mutate(delta_mean = mean.y - mean.x,
         delta_median = median.y - median.x,
         delta_obs = n_proteins.y - n_proteins.x,
         delta_auc = delta_mean) %>%
  ungroup()
summary2 = delta2 %>%
  group_by(label) %>%
  summarise(sd_mean = sd(delta_mean),
            sd_median = sd(delta_median),
            sd_obs = sd(delta_obs),
            delta_mean = mean(delta_mean),
            delta_median = mean(delta_median),
            delta_obs = mean(delta_obs),
            delta_auc = delta_mean,
            sd_auc = sd_mean) %>%
  mutate(text = paste0('# proteins: ', round(delta_obs), '\nAUC: ',
                        ifelse(delta_auc > 0, '+', ''),
                        format(delta_auc, format = 'f', digits = 2) %>% 
                          trimws()))
red = Gpal[1]
p2a = delta2 %>%
  ggplot(aes(x = delta_obs, y = delta_auc)) +
  facet_grid(~ label) +
  geom_hline(aes(yintercept = 0), color = 'grey88', size = 0.65) +
  geom_vline(aes(xintercept = 0), color = 'grey88', size = 0.65) +
  geom_label(data = summary2,
             aes(label = text, y = +Inf, x = +Inf),
             hjust = 1, vjust = 1,
             size = 2, hjust = 0.5, vjust = 0.5, color = red,
             label.size = NA, fill = NA, lineheight = 1,
             label.padding = unit(0.6, 'lines')) +
  geom_point(size = 0.8, shape = 1) + 
  geom_point(data = summary2, color = red, size = 0.8) +
  geom_errorbar(data = summary2, aes(ymin = delta_auc - sd_auc, 
                                     ymax = delta_auc + sd_auc), 
                width = 0, color = red, size = 0.4) +
  geom_errorbarh(data = summary2, aes(xmin = delta_obs - sd_obs, 
                                      xmax = delta_obs + sd_obs), 
                 width = 0, color = red, size = 0.4) +
  scale_x_continuous(expression(paste(Delta, '(# of proteins)')),
                     limits = c(NA, 50)) +
  scale_y_continuous(expression(paste(Delta, '(AUC, SILAC'~-~'iBAQ)'))) +
  boxed_theme() +
  theme(aspect.ratio = 1)
p2a
ggsave("fig/analysis/analysis_grid/SILAC-vs-iBAQ-GO-v1.pdf", p2a,
       width = 9, height = 4.75, units = 'cm', useDingbats = F)

## do SILAC only, same style as main-text
summary2b = filter(summary2, label == 'SILAC')
p2b = delta2 %>%
  filter(label == 'SILAC') %>%
  ggplot(aes(x = delta_obs, y = delta_auc)) +
  facet_grid(~ summary2b$text) +
  geom_hline(aes(yintercept = 0), color = 'grey88', size = 0.65) +
  geom_vline(aes(xintercept = 0), color = 'grey88', size = 0.65) +
  geom_point(size = 0.8, shape = 1) + 
  geom_point(data = summary2c, color = red, size = 0.8) +
  geom_errorbar(data = summary2c, aes(ymin = delta_auc - sd_auc, 
                                     ymax = delta_auc + sd_auc), 
                width = 0, color = red, size = 0.4) +
  geom_errorbarh(data = summary2c, aes(xmin = delta_obs - sd_obs, 
                                      xmax = delta_obs + sd_obs), 
                 width = 0, color = red, size = 0.4) +
  scale_x_continuous(expression(paste(Delta, '(# of proteins)')),
                     limits = c(NA, 50)) +
  scale_y_continuous(expression(paste(Delta, '(AUC, SILAC'~-~'iBAQ)'))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        strip.text = element_text(hjust = 1))
p2b
ggsave("fig/analysis/analysis_grid/SILAC-vs-iBAQ-GO-v2.pdf", p2b,
       width = 4.5, height = 5, units = 'cm', useDingbats = F)
