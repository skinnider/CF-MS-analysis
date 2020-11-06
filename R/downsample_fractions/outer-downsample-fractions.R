setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-downsample-fractions.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)

# check that all of the dependent directories exist
if (!dir.exists("~/git/CF-MS-searches"))
  stop("repository `CF-MS-searches` does not exist")
if (!dir.exists("~/git/network-validation"))
  stop("repository `network-validation` does not exist")

# what system are we on?
system = 'cedar'
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
  system = 'sockeye'
}

# establish grid of analyses
source("R/functions.R") ## contains metrics used to predict interactions
# first metric: mutual information
opts = list(
  analysis = c('complexes', 'GO'),
  metric = 'MI',
  transform = 'none',
  missing = 'zero',
  n_fractions = seq(5, 75, 5)
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F))
# add a second metric: Pearson correlation
grid2 = grid %>%
  mutate(metric = 'pearson')
grid %<>% bind_rows(grid2)

# combine with quantitation strategies (input files)
chrom_dir = '~/git/CF-MS-searches/data/chromatograms'
chrom_files = list.files(chrom_dir, pattern = '*.rds', recursive = T) %>%
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
  distinct(accession, experiment) %>%
  # merge these
  unite(input, accession, experiment, sep = '|')

# rep each analysis over each input
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input = rep(inputs$input, nrow(grid))) %>%
  left_join(inputs, by = 'input') %>%
  separate(input, into = c("accession", "experiment"), sep = "\\|")

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
  dplyr::select(accession, experiment, analysis, metric, transform, missing,
                n_fractions)

# write the raw array
grid_file = "sh/analysis/grids/downsample_fractions_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "downsample_fractions")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "downsample_fractions", accession, 
                                  experiment),
           output_filename = paste0(analysis, 
                                    '-metric=', metric,
                                    '-transform=', transform,
                                    '-missing=', missing,
                                    '-n_fractions=', n_fractions,
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file),
           idx = row_number()) %>%
    filter(!exists) %>%
    dplyr::select(-output_dir, -output_filename, -output_file, -exists,
                  -idx)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/downsample_fractions.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/downsample_fractions.sh",
                " ~/git/CF-MS-analysis/sh/analysis/downsample_fractions.torque.sh")
submit_job(grid0, script, args$allocation, system)
