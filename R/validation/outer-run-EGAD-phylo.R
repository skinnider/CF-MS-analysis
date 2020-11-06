setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-EGAD-phylo.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "st-ljfoster-1")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# detect system
source("R/functions/detect_system.R")

# list network files
network_dir = file.path(base_dir, 'phylo_interactomes')
network_files = list.files(network_dir, pattern = '*rds$', full.names = TRUE)

# set up parameter grid
grid = data.frame(network_file = network_files)

# write the raw array
grid_file = "sh/validation/grids/run_EGAD.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# now, check for which parameters are already complete
grid0 = grid %>%
  mutate(output_dir = file.path(base_dir, "validation", "gene_set_outcomes",
                                "phylo"),
         output_filename = paste0('EGAD-', basename(network_file)),
         output_file = file.path(output_dir, output_filename),
         exists = file.exists(output_file),
         idx = row_number()) %>%
  filter(!exists) %>%
  dplyr::select(-output_dir, -output_filename, -output_file, -exists,
                -idx)

# write the grid that still needs to be run
write.table(grid0, "sh/validation/grids/run_EGAD_phylo.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = paste0("~/git/CF-MS-analysis/sh/validation/run_EGAD_phylo",
                fct_recode(system, '.sh' = 'cedar', '.torque.sh' = 'sockeye', 
                           '.elasti.sh' = 'elasti') %>% as.character())
submit_job(grid0, script, args$allocation, system)
