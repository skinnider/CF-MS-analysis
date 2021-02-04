setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-feature-combinations.R')
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--feature1', type = 'character', required = T)
parser$add_argument('--feature2', type = 'character', required = T)
parser$add_argument('--split_by', type = 'character',
                    choices = c('pairs', 'proteins'), required = T)
parser$add_argument('--n_folds', type = 'integer', required = T)
parser$add_argument('--classifier', type = 'character', 
                    choices = c("NB", "SVM", "RF", "LR"), required = T)
parser$add_argument('--replace_missing_data', type = 'character', required = T)
parser$add_argument('--n_datasets', type = 'integer', required = TRUE)
parser$add_argument('--min_fractions', type = 'integer', default = 4)
parser$add_argument('--sample_idx', type = 'integer', default = 0)
args = parser$parse_args()

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(AUC)
library(PRROC)
library(ontologyIndex)
library(dismay)
source("R/functions/split_by_complex.R")
source("R/functions/get_features.R")
source("R/functions/predict_interactions.R")

# convert args to logical
args$replace_missing_data %<>% as.logical()

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
output_filename = paste0('outcomes', 
                         '-feature1=', args$feature1,
                         '-feature2=', args$feature2,
                         '-min_fractions=', args$min_fractions,
                         '-split_by=', args$split_by,
                         '-n_folds=', args$n_folds,
                         '-classifier=', args$classifier,
                         '-replace_missing_data=', args$replace_missing_data,
                         '-n_datasets=', args$n_datasets,
                         '-sample_idx=', args$sample_idx,
                         '.rds')
output_file = file.path(args$output_dir, output_filename)

# read experiments and sample a certain set
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species == 'Homo sapiens')
set.seed(args$sample_idx)
keep_rows = sample(seq_len(nrow(expts)), args$n_datasets)
datasets = expts %>% extract(keep_rows, )

# read CORUM
complexes = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')
complex_proteins = unique(unlist(complexes))

# filter by maximum complex sizes
n_subunits = lengths(complexes)
n_ppis = map_dbl(n_subunits, ~ . * (. - 1) / 2) %>% as.integer()

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

# 70-30 split the 'held out' set, by proteins
set.seed(args$sample_idx)
split = split_by_proteins(reference, percent_in_train = 0.7)
held_out = split$test
# rename reference set
reference = split$train

# set up cross-validation folds in the reference set
## (n_folds)
percent_in_train = (args$n_folds - 1) / args$n_folds
## (split_by)
set.seed(args$sample_idx)
if (args$split_by == 'proteins') {
  # get CV folds
  folds = cv_by_proteins(reference, n_folds = args$n_folds)
  # get splits
  splits = get_cv_labels(folds, mode = 'pairs')
} else if (args$split_by == 'pairs') {
  # get CV folds
  folds = cv_by_pairs(reference, n_folds = args$n_folds)
  # get splits
  splits = get_cv_labels(folds, mode = 'pairs')
} else {
  stop("not sure how to split by: ", args$split_by)
}
# we now have our held-out set and cross-validation splits

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
    # matrixStats::rowSds(mat[proteins, ], na.rm = T)
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
  ## use the best preprocessing set up per feature
  best_features = readRDS("data/analysis/analysis_grid/best_features.rds") %>%
    extract2('combined')
  grid = best_features %>%
    filter(metric %in% c(args$feature1, args$feature2)) %>%
    group_by(metric) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    # spread 'transformation' back out
    mutate(transform = ifelse(transformation == 'log-transform',
                              'log', 'none'),
           normalize = ifelse(transformation == 'quantile\nnormalization',
                              'quantile', 'none'))
  features = get_features(mat,
                          feature_grid = grid,
                          n_features = 999, 
                          matrix_dir = matrix_dir, 
                          quant_mode = 'iBAQ')
  feats[[mat_idx]] = features
}

# combine features from each matrix
message("concatenating ", length(feats), " features ...")
features = concatenate_features(feats)
## clear features from memory
rm(feats)

# replace missing data
if (args$replace_missing_data) {
  features %<>% replace_missing_data()
} else {
  # we do need to handle infinite values
  # apply same fix as in SCT-MoA
  features = map_dfc(features, ~ {
    if (is.numeric(.)) {
      .[is.infinite(.)] = min(., na.rm = TRUE)
    }
    return(.)
  }) %>% as.data.frame()
  
  if (args$classifier != 'NB') {
    message("WARNING: not sure how to handle missing data with classifier: ", 
            args$classifier, "; removing features with missing data")
    # only NB can handle missing data;
    # for other classifiers, drop missing features
    features %<>% drop_na()
  }
}

# predict interactions
message("predicting interactions ...")
network = predict_interactions(features, meta, splits,
                               classifier = args$classifier,
                               multi_map = FALSE)

# now, calculate outcomes
## 1. CORUM AUC on the held-in set (cross-validation)
message("calculating complex AUC in the held-in set ...")
auc1 = network %>%
  drop_na(label) %$%
  auc(roc(predictions = score,
          labels = factor(label, levels = c('0', '1'))))

## 2. CORUM AUC on the held-out set
message("calculating complex AUC in the held-out set ...")
held_out_adj = held_out %>%
  distinct(protein_A, protein_B) %>%
  adjacency_matrix_from_data_frame()
held_out_labels = make_labels(held_out_adj, network)
auc2 = network %>%
  mutate(held_out_label = held_out_labels) %>%
  drop_na(held_out_label) %$%
  auc(roc(predictions = score,
          labels = factor(held_out_label, levels = c('0', '1'))))

# create data frame
n_proteins = with(network, n_distinct(c(protein_A, protein_B)))
results = data.frame(complex_set = c('held-in', 'held-out'),
                     auc = c(auc1, auc2),
                     n_1 = sum(network$label == 1, na.rm = TRUE),
                     n_0 = sum(network$label == 0, na.rm = TRUE),
                     n_proteins = n_proteins,
                     n_pairs = nrow(network))

# save all outputs
output = list('AUC' = results)
saveRDS(output, output_file)
