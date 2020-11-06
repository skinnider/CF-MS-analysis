setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-phylo-interactomes.R')
parser$add_argument('--output_dir', type = 'character', required = TRUE)
parser$add_argument('--node_id', type = 'character', required = TRUE)
parser$add_argument('--classifier', type = 'character', 
                    choices = c("NB", "SVM", "RF", "RF2", "LR"),
                    required = TRUE)
parser$add_argument('--n_features', type = 'integer', required = TRUE)
parser$add_argument('--feature_select', type = 'character', required = TRUE,
                    choices = c('best_first', 'random', 'PrInCE'))
parser$add_argument('--split_by', type = 'character', required = TRUE,
                    choices = c('pairs', 'proteins'), default = 'proteins')
parser$add_argument('--sample_idx', type = 'integer', default = 0)
# fixed parameters
parser$add_argument('--n_folds', type = 'integer', default = 10)
parser$add_argument('--min_fractions', type = 'integer', default = 4)
parser$add_argument('--n_threads', type = 'integer', default = 1)
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
output_filename = paste0(args$node_id, 
                         '-classifier=', args$classifier,
                         '-n_features=', args$n_features,
                         '-feature_select=', args$feature_select,
                         '-split_by=', args$split_by,
                         '-sample_idx=', args$sample_idx,
                         '.rds')
output_file = file.path(args$output_dir, output_filename)

# get clade species
clade_grid = read.csv("data/analysis/phylo_interactomes/phylo_grid.csv")
species = clade_grid %>%
  filter(node_id == args$node_id) %>% 
  pull(species) %>%
  strsplit(';') %>%
  unlist()

# read experiments
datasets = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species %in% species)

# read CORUM
complexes = read.delim("data/analysis/phylo_interactomes/complexes_euk.txt") %>%
  as_annotation_list('eggNOG', 'complex')
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

# read chromatogram matrices and metadata, and map to eggNOG IDs
mats = map(input_files, ~ {
  input_file = .x
  
  ## infer species
  split = gsub("^.*chromatograms\\/", "", input_file) %>% 
    strsplit('/') %>%
    unlist()
  accession = split[1]
  replicate = split[2]
  expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")
  species = filter(expts, Accession == accession, Replicate == replicate) %>%
    pull(Species)
  ## get the proteome string
  map = read.csv("data/GO/species-map.csv")
  proteome = map$proteome[map$species == species]
  fasta_files = list.files("~/git/CF-MS-searches/data/fasta/filtered")
  proteome_str = fasta_files %>% 
    extract(grepl(proteome, .)) %>%
    gsub("\\.fasta.*$", "", .)
  
  ## read map from UniProt IDs to euk identifiers
  eggnog_file = paste0("data/resources/eggNOG/", proteome_str,
                       ".emapper.annotations.gz")
  ortho = read.delim(eggnog_file,
                     comment.char = '#', header = FALSE,
                     col.names = c('query_name', 'seed_eggNOG_ortholog',
                                   'seed_ortholog_evalue', 'seed_ortholog_score',
                                   'predicted_gene_name', 'GO_terms', 
                                   'KEGG_KOs', 'BiGG_reactions',
                                   'Annotation_tax_scope', 'OGs', 
                                   'bestOG|evalue|score', 'COG cat', 
                                   'eggNOG annot')) %>%
    dplyr::select(query_name, OGs) %>%
    mutate(uniprot = map_chr(strsplit(query_name, '\\|'), 2)) %>%
    # get both KOGs and euk OGs
    mutate(euk = strsplit(OGs, ',') %>%
             map(~ extract(., grepl('euNOG', .) & !grepl("KOG|COG", .))),
           kog = strsplit(OGs, ',') %>%
             map(~ extract(., grepl('euNOG', .) & grepl("KOG", .)))) %>%
    # discard multi-mapping genes
    filter(lengths(euk) == 1 | lengths(kog) == 1) %>%
    # prefer KOG to euk OG
    mutate(eggNOG = map2(euk, kog, ~ ifelse(length(.y) > 0, .y, .x))) %>%
    # remove any remaining multi-mapping genes
    filter(lengths(eggNOG) == 1) %>%
    mutate(eggNOG = unlist(eggNOG)) %>%
    dplyr::select(-euk, -kog) %>%
    distinct(uniprot, eggNOG)
  
  # read matrix
  mat = readRDS(input_file)
  
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
  
  # map protein groups to eggNOG identifiers
  pg2nog = meta %>%
    dplyr::select(`Majority protein IDs`) %>%
    set_colnames("protein_group") %>%
    mutate(uniprot = strsplit(protein_group, ';')) %>%
    unnest(uniprot) %>% 
    left_join(ortho, by = 'uniprot') %>%
    drop_na(eggNOG) %>%
    # remove protein groups that map to >1 eggNOG ID
    group_by(protein_group) %>%
    filter(n_distinct(eggNOG) == 1) %>%
    ungroup() %>%
    distinct(protein_group, eggNOG)
  keep = which(meta$`Majority protein IDs` %in% pg2nog$protein_group)
  mat %<>% extract(keep, )
  meta %<>% extract(keep, )
  
  return(mat)
})

# calculate features on each matrix
feats = list()
for (mat_idx in seq_along(mats)) {
  message("calculating features for matrix ", mat_idx, " of ", length(mats), 
          " ...")
  mat = mats[[mat_idx]]
  
  # read matrices, if we can
  matrix_dir = "/scratch/st-ljfoster-1/CF-MS-analysis/matrices_phylo"
  if (!dir.exists(matrix_dir)) {
    # we must be on cedar
    matrix_dir = "~/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis/matrices_phylo"
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
                               multi_map = FALSE,
                               n_threads = args$n_threads)

# save the network itself
saveRDS(network, output_file)
