setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-EGAD-phylo.R')
parser$add_argument('--network_file', type = 'character', required = TRUE)
args = parser$parse_args()

library(tidyverse)
library(magrittr)
library(PrInCE)
library(igraph)
library(ontologyIndex)
library(flavin)
library(EGAD)
source("R/functions/detect_system.R")
source("R/functions/load_networks.R")
source("R/functions/gene_set_connectivity.R")

# read network
network = readRDS(args$network_file)
network_name = basename(args$network_file)

# read GO
ontology = get_ontology("~/git/network-validation/data/GO/go-basic.obo.gz")
goa = read_gaf("~/git/network-validation/data/GO/goa_human.gaf.gz",
               filter.NOT = T,
               filter.evidence = c("ND", "IPI", "IEA"), 
               ontology = ontology)
# remove roots 
rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
goa %<>% dplyr::filter(!GO_ID %in% rootNames)
# map GO to eggNOG
ortho = read.delim("data/resources/eggNOG/UP000005640-H.sapiens.emapper.annotations.gz",
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
goa %<>%
  dplyr::rename(uniprot = UNIPROT) %>%
  left_join(ortho, by = 'uniprot') %>%
  drop_na(eggNOG)

# create list of gene sets
gene_sets = as_annotation_list(goa, 'eggNOG', 'GO_ID') 
# save the original size of every gene set
gene_set_sizes = data.frame(gene_set = names(gene_sets),
                            n_genes = lengths(gene_sets))
# standardize column names
colnames(network)[c(1, 2)] = c('gene_A', 'gene_B')

# filter network to 50% precision
complexes = read.delim("data/analysis/phylo_interactomes/complexes_euk.txt") %>%
  as_annotation_list('eggNOG', 'complex')
complex_proteins = unique(unlist(complexes))
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(network)
precision = calculate_precision(labels)
network$label = labels
network$precision = precision
network %<>% threshold_precision(0.5)

# alphabetize interactors
sorted = t(apply(network[, 1:2], 1, sort))
network$gene_A = sorted[, 1]
network$gene_B = sorted[, 2]
network %<>% distinct()
# tibbles break EGAD
network %<>% as.data.frame()

# extract relevant gene sets
g = graph_from_data_frame(network, directed = FALSE)
nodes = names(V(g))
ann0 = map(gene_sets, ~ intersect(., nodes)) %>%
  # doesn't make sense to look for anything less than 3 nodes
  extract(lengths(.) >= 3)

## EGAD AUC
message("calculating EGAD AUCs ...")
# make EGAD network
genelist = make_genelist(as.data.frame(network))
gene_network = make_gene_network(network, genelist)
# make annotations
ann_long = map(ann0, ~ data.frame(gene = .)) %>%
  bind_rows(.id = 'gene_set') %>%
  dplyr::select(gene, gene_set)
annotations = make_annotations(ann_long, genelist, names(ann0))
# run GBA with no threshold
gba = run_GBA(gene_network, annotations, min = 3, max = 1e5)
# set up output data frame
aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(annotations)
n_terms = all_terms[names(aurocs)]
# get result
EGAD = data.frame(network = network_name,
                  term = names(aurocs), auroc = aurocs, 
                  n_proteins = n_terms, 
                  pct_proteins = n_terms / length(genelist))

# save output
output_dir = file.path(base_dir, 'validation/gene_set_outcomes/phylo')
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(EGAD, file.path(output_dir, paste0('EGAD-', network_name)))

# save gene set sizes, too
if (!file.exists(file.path(output_dir, 'gene_set_sizes.rds')))
  saveRDS(gene_set_sizes, file.path(output_dir, 'gene_set_sizes.rds'))
