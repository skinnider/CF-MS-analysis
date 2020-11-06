setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
}

# set up IO
input_dir = file.path(base_dir, "downsample_corum")
output_dir = "~/git/CF-MS-analysis/data/analysis/downsample_corum"
if (!dir.exists(output_dir))
  dir.create(output_dir)

# list files
files = list.files(input_dir, recursive = T, pattern = 'rds')

# create data frame
df = data.frame(file = basename(files)) %>%
  separate(file, into = c('x', 'min_fractions', 'split_by', 'n_folds',
                          'classifier', 'replace_missing_data', 'n_features',
                          'feature_select', 'combine_features', 'n_datasets',
                          'downsample_pct', 'sample_idx'), 
           sep = '-') %>%
  dplyr::select(-x) %>%
  mutate_if(is.character, ~ gsub("^.*=|\\.rds$", "", .)) %>%
  type_convert()

# read all data
dats = map(files, ~ readRDS(file.path(input_dir, .x)))

# summarize AUC
auc = map(dats, 'AUC')
n_rows = map_int(auc, nrow)
## fix tables with four rows
aucs = map(auc, ~ {
  if (nrow(.) == 4) {
    head(., 2) %>%
      mutate(complex_set = c("held-in", "held-out")) %>%
      return()
  } else {
    return(.)
  }
})
df0 = df[rep(seq_len(nrow(df)), each = 2), ]
auc = bind_rows(aucs) %>%
  cbind(df0, .)

# save
output_dir = "data/analysis/downsample_corum"
if (!dir.exists(output_dir))
  dir.create(output_dir)
saveRDS(auc, file.path(output_dir, 'AUC.rds'))
