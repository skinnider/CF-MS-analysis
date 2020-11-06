# Define the number of protein complexes needed for successful network 
# inference.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-downsample-corum.R')
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--split_by', type = 'character',
                    choices = c('pairs', 'proteins'), required = T)
parser$add_argument('--n_folds', type = 'integer', required = T)
parser$add_argument('--classifier', type = 'character', 
                    choices = c("NB", "SVM", "RF", "LR"), required = T)
parser$add_argument('--replace_missing_data', type = 'character', required = T)
parser$add_argument('--n_features', type = 'integer', required = TRUE)
parser$add_argument('--feature_select', type = 'character', required = TRUE,
                    choices = c('best_first', 'random', 'PrInCE'))
parser$add_argument('--combine_features', type = 'character', required = T)
parser$add_argument('--n_datasets', type = 'integer', required = TRUE)
parser$add_argument('--min_fractions', type = 'integer', default = 4)
parser$add_argument('--downsample_pct', type = 'double', required = TRUE)
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
args$combine_features %<>% as.logical()

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
output_filename = paste0('outcomes', 
                         '-min_fractions=', args$min_fractions,
                         '-split_by=', args$split_by,
                         '-n_folds=', args$n_folds,
                         '-classifier=', args$classifier,
                         '-replace_missing_data=', args$replace_missing_data,
                         '-n_features=', args$n_features,
                         '-feature_select=', args$feature_select,
                         '-combine_features=', args$combine_features,
                         '-n_datasets=', args$n_datasets,
                         '-downsample_pct=', args$downsample_pct,
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
n_subunits = lengths(corum)

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
  # downsample the proteins
  set.seed(args$sample_idx)
  proteins = with(reference, unique(c(protein_A, protein_B)))
  sample_n = round(length(proteins) * args$downsample_pct)
  proteins = sample(proteins, sample_n)
  reference0 = reference %>% 
    filter(protein_A %in% proteins, protein_B %in% proteins)
  
  # get CV folds
  folds = cv_by_proteins(reference0, n_folds = args$n_folds)
  # get splits
  splits = get_cv_labels(folds, mode = 'pairs')
} else if (args$split_by == 'pairs') {
  # get CV folds
  folds = cv_by_pairs(reference, n_folds = args$n_folds)
  # get splits
  splits = get_cv_labels(folds, mode = 'pairs')
  
  # now, downsample them
  set.seed(args$sample_idx)
  splits %<>% map(~ {
    train = sample_frac(.$train, args$downsample_pct)
    test = sample_frac(.$test, args$downsample_pct)
    list(train = train, test = test)
  })
} else {
  stop("not sure how to split by: ", args$split_by)
}
# we now have our held-out set and cross-validation splits

# now, read the input datasets
input_dirs = file.path("~/git/CF-MS-searches/data/chromatograms", 
                       datasets$Accession, 
                       datasets$Replicate)
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

# optionally, combine across replicates
if (args$combine_features) {
  all_genes = map(mats, rownames) %>% Reduce(union, .)
  mats0 = map(mats, ~ {
    mat1 = .
    missing = setdiff(all_genes, rownames(mat1))
    mat2 = matrix(NA, nrow = length(missing), ncol = ncol(mat1),
                  dimnames = list(missing, colnames(mat1)))
    rbind(mat1, mat2) %>%
      extract(all_genes, )
  })
  mats = do.call(cbind, mats0)
  # pick one representative identifier for any 'duplicated' chromatogram
  # first, sort alphabetically
  order = rownames(mats) %>% order()
  mats %<>% extract(order, )
  # then, prefer CORUM proteins
  order = rownames(mats) %in% complex_proteins %>% order(decreasing = TRUE)
  mats %<>% extract(order, )
  ## now drop duplicated rows
  mats %<>% unique()
  # wrap in a list
  mats = list(mats)
} else {
  # we still need to remove duplicate rows
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
}

# calculate features on each matrix
feats = list()
for (mat_idx in seq_along(mats)) {
  message("calculating features for matrix ", mat_idx, " of ", length(mats), 
          " ...")
  mat = mats[[mat_idx]]
  
  # read matrices, if we can
  matrix_dir = NULL
  if (!args$combine_features) {
    matrix_dir = "/scratch/st-ljfoster-1/CF-MS-analysis/matrices"
    if (!dir.exists(matrix_dir)) {
      # we must be on cedar
      matrix_dir = "~/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis/matrices"
    }
    chrom_dir = input_dirs[mat_idx] %>%
      gsub("^.*chromatograms\\/", "", .)
    matrix_dir %<>% file.path(chrom_dir)
  }
  
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
    # skip Gaussians
    features = PrInCE::calculate_features(mat, gaussians = NULL, 
                                          co_apex = FALSE)
  }
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
    # for other classifiers, drop missing features if !replace_missing_data
    features %<>% drop_na()
  }
}

# predict interactions
message("predicting interactions ...")
network = predict_interactions(features, meta, splits,
                               classifier = args$classifier,
                               multi_map = FALSE)

# now, calculate outcomes
outcomes = data.frame()

# first, calculate the total number of proteins, TPs, and TNs in 
## i. downsampled CORUM
## i. overlapping with the dataset
corum = map(splits, bind_rows) %>% bind_rows()
corum_proteins1 = corum %$% unique(c(protein_A, protein_B)) %>% length()
corum_proteins2 = drop_na(network, label) %$% 
  unique(c(protein_A, protein_B)) %>% length()
corum_TPs1 = sum(corum$label == 1)
corum_TPs2 = sum(network$label == 1, na.rm = TRUE)
corum_TNs1 = sum(corum$label == 0)
corum_TNs2 = sum(network$label == 0, na.rm = TRUE)

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
auc_results = data.frame(complex_set = c('held-in', 'held-out'),
                         auc = c(auc1, auc2),
                         n_proteins = n_proteins,
                         n_pairs = nrow(network),
                         corum_proteins1 = corum_proteins1,
                         corum_proteins2 = corum_proteins2,
                         corum_TPs1 = corum_TPs1,
                         corum_TPs2 = corum_TPs2,
                         corum_TNs1 = corum_TNs1,
                         corum_TNs2 = corum_TNs2)

# save all outputs
output = list('AUC' = auc_results)
saveRDS(output, output_file)
