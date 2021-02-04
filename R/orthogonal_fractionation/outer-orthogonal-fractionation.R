# Perform network inference from mixtures of human SEC and IEX experiments in 
# varying ratios.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-orthogonal-fractionation.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# detect system
source("R/functions/detect_system.R")

# establish grid of analyses
opts = list(
  min_fractions = 4,
  ## ML setup
  split_by = 'proteins', 
  n_folds = 5,
  classifier = 'RF',
  replace_missing_data = TRUE,
  ## features
  n_features = 6,
  feature_select = 'best_first',
  combine_features = FALSE,
  ## datasets
  n_datasets = c(2, 3, 4, 5, 6, 7, 8),
  n_sec = seq(0, 8),
  ## draw multiple samples
  sample_idx = seq_len(20)
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  filter(n_sec <= n_datasets)

# write the raw array
grid_file = "sh/analysis/grids/orthogonal_fractionation_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "orthogonal_fractionation")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "orthogonal_fractionation"),
           output_filename = paste0('outcomes', 
                                    '-min_fractions=', min_fractions,
                                    '-split_by=', split_by,
                                    '-n_folds=', n_folds,
                                    '-classifier=', classifier,
                                    '-replace_missing_data=', replace_missing_data,
                                    '-n_features=', n_features,
                                    '-feature_select=', feature_select,
                                    '-combine_features=', combine_features,
                                    '-n_datasets=', n_datasets,
                                    '-n_sec=', n_sec,
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
write.table(grid0, "sh/analysis/grids/orthogonal_fractionation.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/orthogonal_fractionation.sh",
                " ~/git/CF-MS-analysis/sh/analysis/orthogonal_fractionation.torque.sh")
submit_job(grid0, script, args$allocation, system)
