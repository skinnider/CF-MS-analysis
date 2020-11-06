# Predict a 'consensus' human interactome using different numbers of replicates.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-human-interactome.R')
parser$add_argument('--output_dir', type = 'character', required = TRUE)
parser$add_argument('--classifier', type = 'character', 
                    choices = c("NB", "SVM", "RF", "RF2", "LR"),
                    required = TRUE)
parser$add_argument('--n_features', type = 'integer', required = TRUE)
parser$add_argument('--feature_select', type = 'character', required = TRUE,
                    choices = c('best_first', 'random', 'PrInCE'))
parser$add_argument('--n_datasets', type = 'integer', required = TRUE)
parser$add_argument('--sample_idx', type = 'integer', default = 0)
# fixed parameters
parser$add_argument('--split_by', type = 'character',
                    choices = c('pairs', 'proteins'), default = 'proteins')
parser$add_argument('--n_folds', type = 'integer', default = 10)
parser$add_argument('--min_fractions', type = 'integer', default = 4)
# parse arguments
args = parser$parse_args()

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(AUC)
library(ontologyIndex)
library(dismay)
source("R/functions/split_by_complex.R")
source("R/functions/get_features.R")
source("R/functions/predict_interactions.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
output_filename = paste0('network', 
                         '-classifier=', args$classifier,
                         '-n_features=', args$n_features,
                         '-feature_select=', args$feature_select,
                         '-n_datasets=', args$n_datasets,
                         '-sample_idx=', args$sample_idx,
                         '.rds')
output_file = file.path(args$output_dir, output_filename)

# read experiments and sample a certain set
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species == 'Homo sapiens')
if (args$n_datasets > 1) {
  set.seed(args$sample_idx)
  keep_rows = sample(seq_len(nrow(expts)), args$n_datasets)
  datasets = expts %>% extract(keep_rows, )
} else {
  datasets = expts %>% extract(args$sample_idx, )
}

idxs = numeric(0)
for (sample_idx in seq_len(46)) {
  set.seed(sample_idx)
  idxs[sample_idx] = sample(seq_len(nrow(expts)), 1)
}

# read CORUM
complexes = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')
complex_proteins = unique(unlist(complexes))
n_subunits = lengths(complexes)

# create pairwise interactions 
reference = complexes %>% 
  map_dfr(~ tidyr::crossing(protein_A = ., protein_B = .), .id = 'complex') %>%
  filter(protein_A < protein_B) %>%
  distinct() %>%
  # for interactions found in more than one complex, pick the smaller one
  mutate(n_subunits = n_subunits[complex] %>% unname()) %>%
  group_by(protein_A, protein_B) %>%
  arrange(n_subunits) %>%
  mutate(keep = row_number() == 1) %>%
  ungroup() %>%
  filter(keep)

# set up cross-validation folds in the reference set
percent_in_train = (args$n_folds - 1) / args$n_folds
set.seed(args$sample_idx)
if (args$split_by == 'proteins') {
  folds = cv_by_proteins(reference, n_folds = args$n_folds)
  splits = get_cv_labels(folds, mode = 'pairs')
} else if (args$split_by == 'pairs') {
  folds = cv_by_pairs(reference, n_folds = args$n_folds)
  splits = get_cv_labels(folds, mode = 'pairs')
} else {
  stop("not sure how to split by: ", args$split_by)
}

# now, read the input datasets
input_dirs = file.path("~/git/CF-MS-searches/data/chromatograms", 
                       datasets$Accession, datasets$Replicate)
input_files = file.path(input_dirs, 'iBAQ.rds')

# read chromatogram matrices and metadata, and map to gene names
mats = map(input_files, ~ {
  mat = readRDS(.)
  
  # handle infinite or NaN values
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  # filter by n_fractions
  keep = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) >= args$min_fractions
  mat %<>% extract(keep, )
  
  # also read the corresponding metadata file
  metadata_file = gsub("iBAQ", "metadata", .)
  meta = readRDS(metadata_file) %>%
    extract(keep, )
  
  # keep one row per gene
  gene_map = meta %>%
    dplyr::select(`Majority protein IDs`, `Gene names`) %>%
    set_colnames(c("protein_group", "gene")) %>%
    mutate(gene = strsplit(gene, ';')) %>%
    unnest(gene) %>%
    drop_na()
  genes = unique(gene_map$gene)
  gene_mat = matrix(NA, nrow = length(genes), ncol = ncol(mat),
                    dimnames = list(genes, colnames(mat)))
  n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
  out_map = data.frame()
  for (gene in genes) {
    protein_groups = gene_map$protein_group[gene_map$gene == gene]
    # pick the best protein for this replicate
    n_fractions0 = n_fractions[protein_groups]
    best = names(which(n_fractions0 == max(n_fractions0))) %>%
      dplyr::first()
    gene_mat[gene, ] = mat[best, ]
    # save the mapping
    out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
  }
  
  # set attribute on the gene matrix
  attr(gene_mat, 'gene_map') = out_map
  return(gene_mat)
})

# remove duplicate rows
mats = map(mats, ~ {
  mat = .
  gene_map = attr(mat, 'gene_map')
  order = rownames(mat) %>% order()
  mat %<>% extract(order, )
  gene_map %<>% extract(order, )
  order = rownames(mat) %in% complex_proteins %>% order(decreasing = TRUE)
  mat %<>% extract(order, )
  gene_map %<>% extract(order, )
  ## now drop duplicated rows
  drop = which(duplicated(mat))
  mat %<>% extract(-drop, )
  gene_map %<>% extract(-drop, )
  # reset attribute
  attr(mat, 'gene_map') = gene_map
  return(mat)
})

# calculate features on each matrix
feats = list()
for (mat_idx in seq_along(mats)) {
  message("calculating features for matrix ", mat_idx, " of ", length(mats), 
          " ...")
  mat = mats[[mat_idx]]
  
  # read matrices, if we can
  matrix_dir = "/scratch/st-ljfoster-1/CF-MS-analysis/matrices"
  if (!dir.exists(matrix_dir)) {
    # we must be on cedar
    matrix_dir = "~/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis/matrices"
  }
  chrom_dir = input_dirs[mat_idx] %>%
    gsub("^.*chromatograms\\/", "", .)
  matrix_dir %<>% file.path(chrom_dir)
  
  # calculate features
  if (args$feature_select == "best_first") {
    best_features = readRDS("data/analysis/analysis_grid/best_features.rds") %>%
      extract2('combined') %>%
      # keep only one pipeline per metric
      group_by(metric) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(desc(median))
    grid = head(best_features, args$n_features)
    features = get_features(mat,
                            feature_grid = grid,
                            n_features = args$n_features, 
                            matrix_dir = matrix_dir, 
                            quant_mode = 'iBAQ')
  } else if (args$feature_select == "random") {
    set.seed(args$sample_idx)
    features = get_features(mat,
                            feature_grid = NULL,
                            n_features = args$n_features, 
                            matrix_dir = matrix_dir,
                            quant_mode = 'iBAQ')
  } else if (args$feature_select == "PrInCE") {
    features = PrInCE::calculate_features(mat, gaussians = NULL, 
                                          co_apex = FALSE)
  }
  feats[[mat_idx]] = features
}

# combine features from each matrix
message("concatenating ", length(feats), " features ...")
features = concatenate_features(feats)
## clear feature matrices from memory
rm(feats)

# replace missing data
features %<>% replace_missing_data()

# predict interactions
message("predicting interactions ...")
network = predict_interactions(features, meta, splits,
                               classifier = args$classifier,
                               multi_map = FALSE)

# save the network itself
saveRDS(network, output_file)
