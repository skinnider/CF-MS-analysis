setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# read list of metrics from functions.R
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-write-matrices.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--metric', type = 'character', choices = metrics, 
                    required = T)
parser$add_argument('--transform', type = 'character', 
                    choices = c('none', 'log', 'quantile'), required = T)
parser$add_argument('--missing', type = 'character',
                    choices = c('NA', 'zero', 'noise'), required = T)
args = parser$parse_args()
print(args)

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
quant_mode = gsub("\\.rds$", "", basename(args$input_file))
output_filename = paste0(quant_mode,
                         '-metric=', args$metric,
                         '-transform=', args$transform,
                         '-missing=', args$missing,
                         '.rds')
output_file = file.path(args$output_dir, output_filename)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(AUC)
library(ontologyIndex)
library(dismay)
## other deps: preprocessCore; dismay

# read chromatogram matrix
mat = readRDS(args$input_file)
# also read metadata file
meta = readRDS(gsub("iBAQ", "metadata", args$input_file))

# handle infinite/NaN values
mat[!is.finite(mat)] = NA
mat[is.nan(mat)] = NA

# pre-filter any empty rows
empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) == 0
mat %<>% extract(!empty, )
meta %<>% extract(!empty, )

# step 0: map to eggNOG and _sum_ over euk identifiers
## infer species
split = gsub("^.*chromatograms\\/", "", args$input_file) %>% 
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

## map protein groups to eggNOG identifiers
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

# create eggNOG matrix
nogs = unique(pg2nog$eggNOG)
emat = matrix(NA, nrow = length(nogs), ncol = ncol(mat),
              dimnames = list(nogs, colnames(mat)))
for (nog in nogs) {
  protein_groups = pg2nog$protein_group[pg2nog$eggNOG == nog]
  emat[nog, ] = colSums(mat[protein_groups, , drop = FALSE], na.rm = TRUE)
}
# overwrite the original matrix with eggNOG matrix
mat = emat

# step 1: impute missing values (optional)
if (args$missing == 'zero') {
  mat[is.na(mat)] = 0
} else if (args$missing == 'noise') {
  mat[mat == 0] = NA
  mat %<>% clean_profiles(impute_NA = TRUE, smooth = FALSE)
} else if (args$missing == 'NA') {
  mat[mat == 0] = NA
}

# step 2: log-transform (optional)
if (args$transform == "log") {
  mat %<>% log()
  if (args$missing == 'zero') {
    mat[mat == -Inf] = 0
  }
}

# step 3: quantile normalize (optional)
if (args$transform == 'quantile') {
  dims = dimnames(mat)
  mat %<>% preprocessCore::normalize.quantiles()
  dimnames(mat) = dims
}

# step 4: score pairs
cor = get_coexpr(t(mat), args$metric)

# save matrix
saveRDS(cor, output_file)
