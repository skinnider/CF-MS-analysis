# Define optimal strategies for integrating co-elution information across
# multiple replicates.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-replicate-integration.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# check that all of the dependent directories exist
if (!dir.exists("~/git/CF-MS-searches"))
  stop("repository `CF-MS-searches` does not exist")
if (!dir.exists("~/git/network-validation"))
  stop("repository `network-validation` does not exist")

# detect system
source("R/functions/detect_system.R")

# establish grid of analyses
opts = list(
  min_fractions = 4,
  ## ML setup
  split_by = c('proteins', 'pairs'),
  n_folds = 5,
  classifier = c('NB', 'RF'),
  replace_missing_data = c(TRUE, FALSE),
  ## features
  n_features = c(seq_len(6), seq(8, 20, 2)),
  feature_select = c('best_first', 'random', 'PrInCE'),
  combine_features = c(TRUE, FALSE), # calculate features on combined matrix?
  ## datasets
  n_datasets = c(seq_len(6), 8, 10, 15, 20, 30, 40),
  ## draw multiple samples?
  sample_idx = seq_len(10)
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  # ignore n_features with PrInCE
  filter(!(feature_select == 'PrInCE' & n_features > 1)) %>%
  # only NB can handle missing data
  filter(replace_missing_data | classifier == 'NB') %>%
  # ignore combine_features when n_datasets == 1
  filter(!(combine_features & n_datasets == 1)) %>%
  # only go above ten features in a single dataset
  filter(n_features <= 10 | n_datasets == 1) %>%
  # only go above ten datasets with a single feature
  filter(n_datasets <= 10 | n_features == 1) %>%
  ## note that PrInCE by default uses six features, not one
  filter(!(feature_select == 'PrInCE') & n_datasets > 10)

# write the raw array
grid_file = "sh/analysis/grids/replicate_integration_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "replicate_integration")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "replicate_integration"),
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
write.table(grid0, "sh/analysis/grids/replicate_integration.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/replicate_integration.sh",
                " ~/git/CF-MS-analysis/sh/analysis/replicate_integration.torque.sh")
submit_job(grid0, script, args$allocation, system)
