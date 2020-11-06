setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-human-interactome-PR.R')
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

# label and calculate precision
complexes = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')
complex_proteins = unique(unlist(complexes))
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(ppi)
precision = calculate_precision(labels)
ppi$label = labels
ppi$precision = precision

# log PR curve every 100 pairs (to 200k) and every 1000 pairs thereafter to 5m
pair_idxs = c(seq(0, 200e3, 100), seq(200e3, 5e6, 1e3)) %>% unique()

# extract precision curves at these points
precisions = ppi$precision[pair_idxs]

# create data frame
PR = data.frame(filename = basename(args$input_file),
                idx = tail(pair_idxs, -1),
                precision = precisions) %>%
  # parse filename
  separate(filename, into = c('x', 'classifier', 'n_features', 'feature_select',
                              'n_datasets', 'sample_idx'), sep = '-') %>%
  dplyr::select(-x) %>%
  mutate_at(vars(classifier, n_features, feature_select, n_datasets,
                 sample_idx), ~ gsub("^.*=|\\.rds$", "", .))

# save the PR curve
output = list('PR' = PR)
saveRDS(output, args$output_file)
