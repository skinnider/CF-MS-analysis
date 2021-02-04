# Identify synergistic or antagonistic interactions between pairs of features.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-feature-combinations.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# detect system
source("R/functions/detect_system.R")

# set up a list of metrics
metrics %<>% c('NA') ## one feature only

# get all combinations of metrics
feature_pairs = tidyr::crossing(feature1 = metrics, feature2 = metrics) %>%
  filter(feature1 < feature2 | feature2 == 'NA') %>%
  filter(feature1 != 'NA')

# set up the rest of the grid
opts = list(
  min_fractions = 4,
  ## ML setup
  split_by = 'proteins',
  n_folds = 5,
  classifier = c('NB', 'RF'),
  replace_missing_data = TRUE,
  ## datasets
  n_datasets = c(2:6),
  ## draw multiple samples
  sample_idx = seq_len(10)
)

# combine the two grids
grid = do.call(tidyr::crossing, opts) %>%
  tidyr::crossing(feature_pairs, .)

# either analyze combinations, or a larger range of datasets
grid %<>% 
  filter(feature2 == 'NA' | n_datasets %in% c(3, 6))

# write the raw array
grid_file = "sh/analysis/grids/feature_combinations_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "feature_combinations")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "feature_combinations"),
           output_filename = paste0('outcomes', 
                                    '-feature1=', feature1,
                                    '-feature2=', feature2,
                                    '-min_fractions=', min_fractions,
                                    '-split_by=', split_by,
                                    '-n_folds=', n_folds,
                                    '-classifier=', classifier,
                                    '-replace_missing_data=', replace_missing_data,
                                    '-n_datasets=', n_datasets,
                                    '-sample_idx=', sample_idx,
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file),
           idx = row_number()) %>%
    filter(!exists) %>%
    dplyr::select(-output_dir, -output_filename, -output_file, -exists,
                  -idx)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/feature_combinations.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/feature_combinations.sh",
                " ~/git/CF-MS-analysis/sh/analysis/feature_combinations.torque.sh")
submit_job(grid0, script, args$allocation, system)
