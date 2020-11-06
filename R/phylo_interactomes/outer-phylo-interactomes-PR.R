setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-phylo-interactomes-PR.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# detect system
source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, 'phylo_interactomes')
input_files = list.files(input_dir, pattern = 'rds', full.names = TRUE)

# establish grid of analyses
output_dir = file.path(input_dir, 'PR')
output_files = file.path(output_dir, basename(input_files))
grid = data.frame(input_file = input_files, output_file = output_files)

# now, check for which parameters are already complete
grid0 = filter(grid, !file.exists(output_file))
# analyze complete dataset first
grid0 %<>% arrange(!grepl("=46", input_file))

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/phylo_interactomes_PR.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/phylo_interactomes_PR.sh",
                " ~/git/CF-MS-analysis/sh/analysis/phylo_interactomes_PR.torque.sh")
submit_job(grid0, script, args$allocation, system)
