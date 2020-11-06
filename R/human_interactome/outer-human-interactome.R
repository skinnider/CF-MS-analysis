# Predict a 'consensus' human interactome using different numbers of replicates.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-human-interactome.R')
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
  classifier = c('NB', 'RF', 'LR', 'SVM'),
  n_features = c(1, 2, 4, 6),
  feature_select = 'best_first',
  n_datasets = c(seq(1, 10), seq(12, 20, 2), seq(25, 40, 5), 46),
  ## draw multiple samples?
  sample_idx = seq_len(46)
) 
grid = do.call(tidyr::crossing, opts) %>%
  # no point sampling with all 46 datasets
  filter(!(n_datasets == 46 & sample_idx > 1))

## fixed parameters:
#' min_fractions = 4
#' split_by = proteins
#' n_folds = 10
#' replace_missing_data = TRUE
#' combine_features = FALSE

# write the raw array
grid_file = "sh/analysis/grids/human_interactome_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "human_interactome"),
           output_filename = paste0('network', 
                                    '-classifier=', classifier,
                                    '-n_features=', n_features,
                                    '-feature_select=', feature_select,
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
write.table(grid0, "sh/analysis/grids/human_interactome.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/human_interactome.sh",
                " ~/git/CF-MS-analysis/sh/analysis/human_interactome.torque.sh")
submit_job(grid0, script, args$allocation, system)
