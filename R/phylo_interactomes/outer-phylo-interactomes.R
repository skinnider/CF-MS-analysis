setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-phylo-interactomes.R')
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
  classifier = 'RF',
  n_features = 1,
  feature_select = 'best_first',
  split_by = 'proteins',
  sample_idx = 1
) 
grid = do.call(tidyr::crossing, opts)

# rep this grid over each possible clade
clade_grid = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv") %>%
  dplyr::select(node_id)
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(clade_grid))) %>%
  mutate(node_id = rep(clade_grid$node_id, nrow(grid))) %>%
  left_join(clade_grid, by = 'node_id') %>%
  dplyr::select(node_id, everything())

# write the raw array
grid_file = "sh/analysis/grids/phylo_interactomes_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "phylo_interactomes"),
           output_filename = paste0(node_id, 
                                    '-classifier=', classifier,
                                    '-n_features=', n_features,
                                    '-feature_select=', feature_select,
                                    '-split_by=', split_by,
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
write.table(grid0, "sh/analysis/grids/phylo_interactomes.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/phylo_interactomes.sh",
                " ~/git/CF-MS-analysis/sh/analysis/phylo_interactomes.torque.sh")
submit_job(grid0, script, args$allocation, system)
