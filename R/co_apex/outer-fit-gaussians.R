setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-fit-gaussians.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/submit_job.R")

# what system are we on?
source("R/functions/detect_system.R")

# establish grid of analyses
chrom_dir = '~/git/CF-MS-searches/data/chromatograms'
chrom_files = list.files(chrom_dir, pattern = '*.rds', recursive = T) %>%
  # ignore the metadata files
  extract(!grepl("metadata", .))
split = strsplit(chrom_files, '/')
accessions = map_chr(split, 1)
experiments = map_chr(split, 2)
quant_modes = gsub("\\..*$", "", basename(chrom_files))
grid = data.frame(input_file = file.path(chrom_dir, chrom_files),
                  accession = accessions,
                  experiment = experiments,
                  quant_mode = quant_modes)

# write the raw array
grid_file = "sh/analysis/grids/fit_gaussians_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "fit_gaussians")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "fit_gaussians", accession, 
                                  experiment),
           output_filename = paste0('gaussians-', quant_mode, '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file),
           idx = row_number()) %>%
    filter(!exists) %>%
    dplyr::select(input_file, output_file)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/fit_gaussians.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/fit_gaussians.sh",
                " ~/git/CF-MS-analysis/sh/analysis/fit_gaussians.torque.sh")
submit_job(grid0, script, args$allocation, system)
