setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
}

# list files
input_dir = file.path(base_dir, "phylo_interactomes", "PR")
files = list.files(input_dir, recursive = T, pattern = 'rds',
                   full.names = TRUE)

# read data
dats = map(files, readRDS)

# extract PR curves
curves = map(dats, 'PR') %>%
  bind_rows()

# extract # of PPIs
rows = map(dats, 'PR') %>%
  map(~ extract(., 1, )) %>%
  bind_rows() %>%
  dplyr::select(-idx, -precision, -n_proteins) %>%
  type_convert()
n_ppis = map(seq_along(dats), ~ data.frame(precision = seq(0.5, 0.99, 0.01), 
                                           n_ppis = dats[[.x]]$n_ppis) %>% 
               cbind(rows[.x, ], .)) %>%
  bind_rows()

# save
output_dir = "data/analysis/phylo_interactomes"
if (!dir.exists(output_dir))
  dir.create(output_dir)
saveRDS(auc, file.path(output_dir, 'AUC.rds'))
saveRDS(n_ppis, file.path(output_dir, 'n_ppis.rds'))
