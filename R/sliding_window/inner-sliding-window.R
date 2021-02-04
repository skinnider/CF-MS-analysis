setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# read list of metrics from functions.R
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-sliding-window.R')
parser$add_argument('--input_dir', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--analysis', type = 'character',
                    choices = c('GO', 'complexes'), required = T)
parser$add_argument('--metric', type = 'character', choices = metrics, 
                    required = T)
parser$add_argument('--transform', type = 'character', 
                    choices = c('none', 'log', 'quantile'), required = T)
parser$add_argument('--missing', type = 'character',
                    choices = c('NA', 'zero', 'noise'), required = T)
parser$add_argument('--window', type = 'integer', required = T)
args = parser$parse_args()
print(args)

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
output_filename = paste0(args$analysis, 
                         '-metric=', args$metric,
                         '-transform=', args$transform,
                         '-missing=', args$missing,
                         '-window=', args$window,
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

# get species from experiments.csv
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")
split = strsplit(args$input_dir, '/') %>% unlist()
accession = split[length(split) - 2]
replicate = split[length(split) - 1]
species = expts %>%
  filter(Accession == accession, Replicate == replicate) %>%
  pull(Species)
if (length(species) != 1) 
  stop("couldn't get species for input directory: ", args$input_dir)

# list all files
input_files = list.files(args$input_dir, pattern = '*.rds', full.names = T) %>%
  extract(!grepl("metadata", .))
quant_modes = gsub("\\.rds$", "", basename(input_files))

# set up output data frame
results = data.frame()

# iterate over input files
for (input_file in input_files) {
  quant_mode = input_file %>%
    basename() %>%
    gsub("\\.rds$", "", .)
  message("analyzing input file #", which(input_files == input_file), " of ", 
          length(input_files), ": ", quant_mode)
  
  # read chromatogram matrix
  mat = readRDS(input_file)
  # handle infinite or NaN values
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  # abort if there aren't more fractions than the window size
  if (ncol(mat) < args$window)
    next
  
  # loop through windows
  windows = list()
  for (starting_pos in seq(1, ncol(mat) - args$window + 1)) {
    message("working on a window of ", args$window, 
            " fractions, starting at fraction ", starting_pos, " ...")
    
    # subset the matrix
    window_idxs = seq(starting_pos, starting_pos + args$window - 1)
    mat0 = mat[, window_idxs]
    
    # get correlation matrix for the subset matrix
    empty = rowSums(!is.na(mat0) & is.finite(mat0) & mat0 != 0) == 0
    mat0 %<>% extract(!empty, , drop = FALSE)
    # store number of paired observations in the original matrix
    n = crossprod(!is.na(t(mat0)) & t(mat0) != 0)
    
    # skip if empty 
    if (nrow(mat0) <= 1) 
      next
    
    # also read the corresponding metadata file
    metadata_file = file.path(args$input_dir, 'metadata.rds')
    meta = readRDS(metadata_file) %>%
      extract(!empty, )
    
    # step 1: impute missing values (optional)
    if (args$missing == 'zero') {
      mat0[is.na(mat0)] = 0
    } else if (args$missing == 'noise') {
      mat0[mat0 == 0] = NA
      mat0 %<>% clean_profiles(impute_NA = TRUE, smooth = FALSE)
    } else if (args$missing == 'NA') {
      mat0[mat0 == 0] = NA
    }
    
    # step 2: log-transform (optional)
    if (args$transform == "log") {
      mat0 %<>% log()
      if (args$missing == 'zero') {
        mat0[mat0 == -Inf] = 0
      }
    }
    
    # step 3: quantile normalize (optional)
    if (args$transform == 'quantile') {
      dims = dimnames(mat0)
      mat0 %<>% preprocessCore::normalize.quantiles()
      dimnames(mat0) = dims
    }
    
    # step 4: score pairs
    cor = try(get_coexpr(t(mat0), args$metric))
    
    # catch errors or exit
    if (class(cor) == 'try-error') {
      # allow only one particular exception - from treeClust
      if (grepl("Error in leaves|No tree produced anything", cor)) {
        next
      } else if (args$metric %in% c('rho_p', 'phi_s') & 
                 args$transform == 'log') {
        # ignore log-transformed ratios in proportionality
        next
      } else {
        stop(cor)
      }
    }
    
    # now, melt to data frame
    pairs = cor %>%
      reshape2::melt(varnames = c('protein_A', 'protein_B'), 
                     value.name = 'cor') %>%
      filter(as.integer(protein_A) < as.integer(protein_B)) %>%
      # flag window
      mutate(starting_pos = starting_pos)
    # append to list
    windows[[starting_pos]] = pairs
  }
  
  # coerce all matrices into a matrix of fixed dimension
  windows0 = map(windows, ~ {
    empty = matrix(NA, nrow = nrow(mat), ncol = nrow(mat),
                   dimnames = list(rownames(mat), rownames(mat)))
    idxing_mat = .x[, 1:2] %>% map_dfc(as.character) %>% as.matrix()
    empty[idxing_mat] = .x$cor
    return(empty)
  })
  
  # finally, pick out the highest value for each protein
  cor = do.call(pmax, c(windows0, na.rm = TRUE))
  
  # now, run one of two analyses
  if (args$analysis == 'complexes') {
    # throw an error if complexes is being run on a non-human/mouse species
    if (!species %in% c("Homo sapiens", "Mus musculus")) {
      stop("can't run complex analysis on species: ", species, " (file: ",
           input_file, ")")
    }
    
    # read complexes
    complex_species = fct_recode(species,
                                 'human' = 'Homo sapiens',
                                 'mouse' = 'Mus musculus')
    complex_file = paste0("~/git/network-validation/data/complex/CORUM", 
                          "/complexes_", complex_species, ".txt")
    adj = read.delim(complex_file) %>%
      as_annotation_list(., 'gene_name', 'complex') %>% 
      adjacency_matrix_from_list()
    
    # map proteins to genes
    map = with(meta, setNames(`Gene names`, `Majority protein IDs`))
    rownames(cor) = colnames(cor) = unname(map[rownames(cor)])
    # for this analysis, drop protein groups that map to >1 gene
    keep = !grepl(";", rownames(cor))
    cor %<>% extract(keep, keep)
    
    # convert correlations to pairwise data frame
    pairs = reshape2::melt(cor, varnames = c('gene1', 'gene2'),
                           value.name = 'cor') %>%
      drop_na() %>%
      filter(as.integer(gene1) < as.integer(gene2))
    
    # calculate AUC
    pairs0 = pairs %>%
      filter(gene1 %in% rownames(adj), gene2 %in% rownames(adj))
    x = pairs0$cor
    y = adj[as.matrix(pairs0[, 1:2])]
    ## catch cases where no complex pairs exist
    if (nrow(pairs0) == 0 | n_distinct(y) < 2) {
      # append empty row
      row = data.frame(quant_mode = quant_mode,
                       auroc = NA)
      results %<>% bind_rows(row)
    } else {
      auroc = auc(roc(predictions = x,
                      labels = factor(y, levels = c('0', '1'))))
      
      # append results
      row = data.frame(quant_mode = quant_mode,
                       auroc = auroc)
      results %<>% bind_rows(row)
    }
  } else if (args$analysis == 'GO') {
    # read GO
    go = get_ontology("~/git/network-validation/data/GO/goslim_generic.obo.gz")
    
    # read map from species to GOA file
    species_map = read.csv("data/GO/species-map.csv")
    goa_filename = species_map %>%
      filter(species == !!species) %>%
      pull(file)
    if (length(goa_filename) != 1) 
      stop("couldn't get GOA file for species: ", species,
           " (file: ", input_file, ")")
    
    # read GOA
    goa_file = paste0("data/GO/GOA/", goa_filename, ".gz")
    goa = read_tsv(goa_file, comment = '!', col_names = gaf_colnames) 
    
    # filter to GO slim terms
    goa_slim = goa %>%
      filter(`GO ID` %in% go$id) %>%
      # remove roots
      filter(!`GO ID` %in% c(BP = "GO:0008150", 
                             CC = "GO:0005575",
                             MF = "GO:0003674"))
    
    # convert to annotation list
    ann = as_annotation_list(goa_slim, 'DB Object ID', 'GO ID') %>%
      # remove any terms annotated to fewer than three proteins
      extract(lengths(.) >= 3)
    
    # keep only proteins that have at least one GO slim annotation
    keep = meta$`Majority protein IDs` %>%
      strsplit(';') %>%
      map_lgl(~ any(. %in% goa_slim$`DB Object ID`))
    cor %<>% extract(keep, keep)
    meta %<>% extract(keep, )
    
    # convert correlations to pairwise data frame
    pairs = reshape2::melt(cor, varnames = c('gene1', 'gene2'),
                           value.name = 'cor') %>%
      drop_na() %>%
      filter(as.integer(gene1) < as.integer(gene2))
    
    # map over GO slim terms
    for (term_idx in seq_along(ann)) {
      go_term = names(ann)[term_idx]
      message("[", term_idx, "/", length(ann), "] analyzing GO term: ",
              go_term, " ...")
      targets = ann[[term_idx]]
      
      # map protein groups to GO terms
      protein_groups = meta$`Majority protein IDs`
      annotations = protein_groups %>%
        strsplit(';') %>%
        map(~ . %in% targets)
      true_positives = protein_groups[map_lgl(annotations, ~ all(.) == TRUE)]
      # drop any protein groups with discordant mappings
      discordant = protein_groups[map_int(annotations, n_distinct) != 1]
      
      # create two vectors
      pairs0 = pairs %>%
        filter(!gene1 %in% discordant, !gene2 %in% discordant) %>%
        mutate(y = ifelse(gene1 %in% true_positives & gene2 %in% true_positives,
                          1, 0))
      x = pairs0$cor
      y = pairs0$y
      if (n_distinct(y) < 2) {
        message("  skipping GO term")
        next
      }
      auroc = auc(roc(predictions = x,
                      labels = factor(y, levels = c('0', '1'))))
      
      # append results
      row = data.frame(quant_mode = quant_mode,
                       go_term = go_term,
                       n_proteins = length(targets),
                       n_chromatograms = length(true_positives),
                       auroc = auroc)
      results %<>% bind_rows(row)
    }
  }
}

# save results
saveRDS(results, output_file)
