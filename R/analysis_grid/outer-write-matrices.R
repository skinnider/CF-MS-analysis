# For a subset of metrics that take a very long time to compute, calculate
# correlation matrices and write them to a file.
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

# check that all of the dependent directories exist
if (!dir.exists("~/git/CF-MS-searches"))
  stop("repository `CF-MS-searches` does not exist")

# what system are we on?
source("R/functions/detect_system.R")

# establish grid of analyses
source("R/functions.R") ## contains metrics used to predict interactions
opts = list(
  analysis = c('complexes', 'GO'),
  metric = metrics,
  transform = c('none', 'quantile', 'log'),
  missing = c('NA', 'zero', 'noise')
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  # filter some combinations that cannot handle NAs
  filter(!(missing == 'NA' & 
             metric %in% c('phi_s', 'rho_p', 'GENIE3', 'distance_cor',
                           'cosine', 'wccor'))) %>%
  # filter some combinations for which zeros and NAs are identical
  filter(!(missing == 'NA' &
             metric %in% c('dice', 'hamming', 'jaccard', 'zi_kendall',
                           'binomial', 'bayes_cor'))) 

# combine with quantitation strategies (input files)
chrom_dir = file.path(base_dir, 'chromatograms')
chrom_files = list.files(chrom_dir, pattern = '*.rds', recursive = TRUE) %>%
  # ignore the metadata files
  extract(!grepl("metadata", .))
split = strsplit(chrom_files, '/')
accessions = map_chr(split, 1)
experiments = map_chr(split, 2)
quant_modes = gsub("\\..*$", "", basename(chrom_files))
inputs = data.frame(file = file.path(chrom_dir, chrom_files),
                    accession = accessions,
                    experiment = experiments,
                    quant_mode = quant_modes) %>%
  # process quant modes in 'inner' script
  distinct(accession, experiment, quant_mode) %>%
  # merge these
  unite(input, accession, experiment, quant_mode, sep = '|')

# rep each analysis over each input
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input = rep(inputs$input, nrow(grid))) %>%
  left_join(inputs, by = 'input') %>%
  separate(input, into = c("accession", "experiment",  "quant_mode"), 
           sep = "\\|")

# filter complex analysis where species is not human or mouse
experiments = read.csv("~/git/CF-MS-searches/data/experiments.csv") 
species = experiments %>% dplyr::select(Accession, Replicate, Species) %>%
  dplyr::rename(species = Species)
grid %<>%
  left_join(species, by = c('accession' = 'Accession', 
                            'experiment' = 'Replicate')) %>%
  filter(!(analysis == 'complexes' & 
             !species %in% c("Homo sapiens", "Mus musculus")))

# clean up grid
grid %<>%
  dplyr::select(accession, experiment, quant_mode, analysis, metric, transform, 
                missing)

# now, filter to the long metrics
grid %<>% 
  filter(metric %in% c("MI", "zi_kendall", "wccor", "GENIE3", 
                       "distance_cor")) %>%
  # complexes and GO are the same thing, in this outer
  filter(analysis == 'GO') %>%
  dplyr::select(-analysis)

# write the raw array
grid_file = "sh/analysis/grids/write_matrices_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "matrices")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "matrices", accession, 
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

# this is a big grid, do only part of it at a time if needed
if (nrow(grid0) >= 10000) {
  grid0 %<>% dplyr::slice(1:9900) ## allow for some other running jobs or sh
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/write_matrices.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/write_matrices.sh",
                " ~/git/CF-MS-analysis/sh/analysis/write_matrices.torque.sh")
submit_job(grid0, script, args$allocation, system)
