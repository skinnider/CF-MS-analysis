# Define the number of protein complexes needed for successful prediction.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-downsample-corum.R')
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
  split_by = 'proteins',
  n_folds = 5,
  classifier = c('NB', 'RF'),
  replace_missing_data = TRUE, 
  n_features = 6,
  feature_select = 'best_first',
  combine_features = FALSE,
  n_datasets = seq(2, 4),
  downsample_pct = c(0.001, seq(0.01, 0.05, 0.005), seq(0.1, 1, 0.05)),
  sample_idx = seq_len(10)
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  # minimum depth is different, depending on how we split the proteins
  filter(split_by == 'pairs' | downsample_pct >= 0.05)

# write the raw array
grid_file = "sh/analysis/grids/downsample_corum_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "downsample_corum")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "downsample_corum"),
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
                                    '-downsample_pct=', downsample_pct,
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
write.table(grid0, "sh/analysis/grids/downsample_corum.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/downsample_corum.sh",
                " ~/git/CF-MS-analysis/sh/analysis/downsample_corum.torque.sh")
submit_job(grid0, script, args$allocation, system)
