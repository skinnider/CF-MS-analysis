setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-phylo-interactomes-PR.R')
parser$add_argument('--input_file', type = 'character', required = TRUE)
parser$add_argument('--output_file', type = 'character', required = TRUE)
args = parser$parse_args()

library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)
library(AUC)
library(PRROC)

# check output dir
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir)

# read the data 
ppi = readRDS(args$input_file)

# calculate precision
complexes = read.delim("data/analysis/phylo_interactomes/complexes_euk.txt") %>%
  as_annotation_list('eggNOG', 'complex')
complex_proteins = unique(unlist(complexes))
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(ppi)
precision = calculate_precision(labels)
ppi$label = labels
ppi$precision = precision

# calculate exact number of interactions at various precisions
thresholds = seq(0.5, 0.99, 0.01)
n_ppis = map_dbl(thresholds, ~ max(which(ppi$precision >= .x)))

# log every 100 pairs (to 200k) and every 1000 pairs thereafter to 5m
pair_idxs = c(seq(0, 200e3, 100), seq(200e3, 5e6, 1e3)) %>% unique()

# extract precision curves at these points
precisions = ppi$precision[pair_idxs]

# create data frame
PR = data.frame(filename = basename(args$input_file),
                idx = tail(pair_idxs, -1),
                precision = precisions) %>%
  # parse filename
  separate(filename, into = c('node_id', 'classifier', 'n_features',
                              'feature_select', 'split_by', 'sample_idx'),
           sep = '-') %>%
  mutate_at(vars(classifier, n_features, feature_select, split_by, sample_idx),
            ~ gsub("^.*=|\\.rds$", "", .))

# set up output
output = list(PR = PR, n_ppis = n_ppis)

# save the PR curve
saveRDS(output, args$output_file)
