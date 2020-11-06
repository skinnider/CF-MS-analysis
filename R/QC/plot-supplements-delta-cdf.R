# Compare the difference in the number of proteins quantified between the 
# consistently searched data and the original paper supplements.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read protein quantitation CDF
cdf1 = readRDS("data/QC/protein-groups-CDF.rds") 

# now, generate an identical CDF for the supplements
cdf2_file = 'data/QC/supplements-CDF.rds'
if (file.exists(cdf2_file)) {
  cdf2 = readRDS(cdf2_file)
} else {
  # list all preprocessed supplement files
  files = list.files('~/git/CF-MS-searches/data/supplements', pattern = '*.rds', 
                     full.names = T)
  
  # extract accessions
  accessions = gsub("\\..*$", "", basename(files))
  
  # construct CDF for each file
  cdf2 = data.frame()
  for (file_idx in seq_along(files)) {
    file = files[file_idx]
    accession = accessions[file_idx]
    message('[', file_idx, '/', length(files), '] working on file: ', file)
    
    # extract CDF for this file
    cdf0 = data.frame()
    mats = readRDS(file)
    for (experiment in names(mats)) {
      mat = mats[[experiment]]      
      for (n_fractions in seq(0, ncol(mat))) {
        n_proteins = mat %>%
          replace(is.na(.), 0) %>%
          extract(rowSums(. > 0) >= n_fractions, , drop = F) %>%
          nrow()
        cdf0 %<>% bind_rows(data.frame(accession = accession, 
                                       experiment = experiment,
                                       n_fractions = n_fractions,
                                       n_proteins = n_proteins))
      }
    }
    
    # append to master CDF
    cdf2 %<>% bind_rows(cdf0)
  }
  
  # rearrange columns
  cdf2 %<>% 
    dplyr::select(accession, experiment, n_fractions, n_proteins)
  
  # save 
  saveRDS(cdf2, cdf2_file)
}

# take best number of proteins per fraction in the original data
best = cdf1 %>%
  group_by(accession, experiment, n_fractions) %>%
  # arrange(desc(n_proteins)) %>%
  # dplyr::slice(1) %>%
  summarise(n_proteins = max(n_proteins)) %>%
  ungroup()

# calculate difference relative to original
delta = best %>%
  left_join(cdf2, by = c('accession', 'experiment', 'n_fractions')) %>%
  mutate(delta = n_proteins.x - n_proteins.y) %>%
  # filter to a maximum of 50 fractions
  filter(between(n_fractions, 1, 50)) 

# remove a couple experiments where FDR >> 1% in the original pub
delta %<>%
  filter(accession != 'PXD002640')

# calculate the mean curve
mean = delta %>%
  group_by(n_fractions) %>%
  summarise(mean = mean(delta, na.rm = T),
            sd = sd(delta, na.rm = T),
            median = median(delta, na.rm = T),
            n = n()) %>%
  ungroup()

# plot each curve, plus the mean curve
## label n=10 fractions
lab = filter(mean, n_fractions == 10) %>%
  mutate(text = paste0('+', round(mean), 
                       ' protein groups\nin 10+ fractions'))
p1 = delta %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = delta)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.2) + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted', color = 'black') +
  geom_path(data = mean, aes(y = mean), color = kinney6[4], size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = mean), size = 2,
                   segment.size = 0.25, label.padding = 0.4, label.size = NA,
                   nudge_y = -2400, nudge_x = +12, hjust = 0) +
  scale_x_continuous('Fractions', expand = c(0, 0),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous(expression(Delta(protein~groups,~original~analysis))) +
  boxed_theme()
p1
ggsave("fig/QC/supplements-delta-cdf.pdf", p1,
       width = 6, height = 5, units = "cm", useDingbats = F)
