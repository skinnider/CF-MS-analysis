setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-write-matrices.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/detect_system.R")
source("R/functions.R")

# check that all of the dependent directories exist
if (!dir.exists("~/CF-MS-searches/CF-MS-searches"))
  stop("repository `CF-MS-searches` does not exist")

# read best features
best_features = readRDS("data/analysis/analysis_grid/best_features.rds") %>%
  extract2('combined') %>%
  # keep only one metric
  group_by(metric) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(desc(median))
# pre-generate matrices for up to three features
grid = head(best_features, 3)

# combine with input files
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") 
species = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv")$
  species[1] %>%
  strsplit(';') %>% 
  unlist()
expts %<>% filter(Species %in% species)
chrom_dirs = file.path("~/git/CF-MS-searches/data/chromatograms", 
                       expts$Accession, expts$Replicate)
chrom_files = file.path(chrom_dirs, 'iBAQ.rds')
split = chrom_files %>% gsub("^.*chromatograms/", "", .) %>% strsplit('/')
accessions = map_chr(split, 1)
experiments = map_chr(split, 2)
quant_modes = gsub("\\..*$", "", basename(chrom_files))
inputs = data.frame(file = chrom_files,
                    accession = accessions,
                    experiment = experiments,
                    quant_mode = quant_modes) %>%
  # merge these
  unite(input, accession, experiment, quant_mode, sep = '|')

# rep each analysis over each input
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input = rep(inputs$input, nrow(grid))) %>%
  left_join(inputs, by = 'input') %>%
  separate(input, into = c("accession", "experiment", "quant_mode"),
           sep = "\\|")

# clean up grid
grid %<>%
  dplyr::select(accession, experiment, quant_mode, metric, transform, missing)

# write the raw array
grid_file = "sh/analysis/grids/phylo_matrices_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "matrices_phylo")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "matrices_phylo", accession, 
                                  experiment),
           output_filename = paste0(quant_mode, 
                                    '-metric=', metric,
                                    '-transform=', transform,
                                    '-missing=', missing,
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file),
           idx = row_number()) %>%
    filter(!exists) %>%
    dplyr::select(-output_dir, -output_filename, -output_file, -exists,
                  -idx)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/phylo_matrices.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/phylo_matrices.sh",
                " ~/git/CF-MS-analysis/sh/analysis/phylo_matrices.torque.sh")
submit_job(grid0, script, args$allocation, system)
