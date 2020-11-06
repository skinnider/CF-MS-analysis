setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-EGAD.R')
parser$add_argument('--network_idx', type = 'integer', default = 1)
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

# read all networks
networks = load_networks()

# extract the network being analyzed in this run
network = networks[[args$network_idx]]
network_name = names(networks)[args$network_idx]

# read GO
ontology = get_ontology("~/git/network-validation/data/GO/go-basic.obo.gz")
goa = read_gaf("~/git/network-validation/data/GO/goa_human.gaf.gz",
               filter.NOT = T,
               filter.evidence = c("ND", "IPI", "IEA"), 
               ontology = ontology)
# remove roots 
rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
goa %<>% dplyr::filter(!GO_ID %in% rootNames)

# create list of gene sets
gene_sets = as_annotation_list(goa, 'DB_Object_Symbol', 'GO_ID') %>%
  setNames(paste0('GO|', names(.)))
# save the original size of every gene set
gene_set_sizes = data.frame(gene_set = names(gene_sets),
                            n_genes = lengths(gene_sets))
# standardize column names
colnames(network)[c(1, 2)] = c('gene_A', 'gene_B')

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
output_dir = file.path(base_dir, 'validation/gene_set_outcomes')
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(EGAD, file.path(output_dir, paste0('EGAD-network_idx=', 
                                           args$network_idx, '.rds')))

# save gene set sizes, too
if (!file.exists(file.path(output_dir, 'gene_set_sizes.rds')))
  saveRDS(gene_set_sizes, file.path(output_dir, 'gene_set_sizes.rds'))
