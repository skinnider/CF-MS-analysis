setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-fractions-vs-replicates.R')
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
opts = list(
  analysis = c('complexes', 'GO'),
  metric = 'pearson',
  transform = 'none',
  missing = 'zero',
  min_fractions = 4,
  min_pairs = 0,
  n_fractions = seq(10, 100, 5),
  n_replicates = seq_len(5),
  sample_idx = seq_len(10)
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F))
## repeat with mutual information
grid2 = mutate(grid, metric = 'MI')
grid %<>% bind_rows(grid2)

# combine with input files
# manually set up a list of experiments with replicates
inputs = c('PXD001220' = 'all',
           'PXD002319' = 'Ce_Bead',
           'PXD002320' = 'Dd_HCW',
           'PXD002324' = 'Nv_IEX',
           'PXD002325' = 'Sp_Bead_IEX',
           'PXD002892' = 'BN_medium',
           'PXD002892' = 'BN_heavy',
           'PXD002892' = 'SEC_medium',
           'PXD002892' = 'SEC_heavy'
           )
# now, convert to data frame
inputs = data.frame(accession = names(inputs),
                    pattern = inputs)
# merge inputs
inputs %<>% unite(input, accession, pattern, sep = '|')

# rep each analysis over each input
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input = rep(inputs$input, nrow(grid))) %>%
  left_join(inputs, by = 'input') %>%
  separate(input, into = c("accession", "pattern"), sep = "\\|")

# ignore some combinations that can't be run
max_repls = c('PXD001220' = 3,
              'PXD002319' = 8,
              'PXD002320' = 3,
              'PXD002324' = 5,
              'PXD002325' = 3,
              'PXD002892' = 3) %>%
  data.frame(accession = names(.), max_replicates = .)
grid %<>% 
  left_join(max_repls, by = 'accession') %>%
  filter(n_replicates <= max_replicates) %>%
  dplyr::select(-max_replicates)

# filter complex analysis where species is not human or mouse
experiments = read.csv("~/git/CF-MS-searches/data/experiments.csv") 
species = experiments %>% dplyr::select(Accession, Species) %>%
  dplyr::rename(accession = Accession, species = Species) %>%
  distinct()
grid %<>%
  left_join(species, by = 'accession') %>%
  filter(!(analysis == 'complexes' & 
             !species %in% c("Homo sapiens", "Mus musculus")))

# clean up grid
grid %<>%
  dplyr::select(accession, pattern, analysis,
                metric, transform, missing, min_fractions, min_pairs,
                n_fractions, n_replicates, sample_idx)

# write the raw array
grid_file = "sh/analysis/grids/fractions_vs_replicates_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "fractions_vs_replicates")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    mutate(output_dir = file.path(base_dir, "fractions_vs_replicates",
                                  accession, pattern),
           output_filename = paste0(analysis, 
                                    '-metric=', metric,
                                    '-transform=', transform,
                                    '-missing=', missing,
                                    '-min_fractions=', min_fractions,
                                    '-min_pairs=', min_pairs,
                                    '-n_fractions=', n_fractions,
                                    '-n_replicates=', n_replicates,
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
write.table(grid0, "sh/analysis/grids/fractions_vs_replicates.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/analysis/fractions_vs_replicates.sh",
                " ~/git/CF-MS-analysis/sh/analysis/fractions_vs_replicates.torque.sh")
submit_job(grid0, script, args$allocation, system)
