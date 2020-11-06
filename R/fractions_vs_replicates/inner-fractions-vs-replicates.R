setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# read list of metrics from functions.R
source("R/functions.R")

# parse arguments
parser = ArgumentParser(prog = 'inner-fractions-vs-replicates.R')
parser$add_argument('--input_dir', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--pattern', type = 'character', required = T)
parser$add_argument('--analysis', type = 'character',
                    choices = c('GO', 'complexes'), required = T)
parser$add_argument('--metric', type = 'character', choices = metrics, 
                    required = T)
parser$add_argument('--transform', type = 'character', 
                    choices = c('none', 'log'), required = T)
parser$add_argument('--missing', type = 'character',
                    choices = c('NA', 'zero', 'noise'), required = T)
parser$add_argument('--min_fractions', type = 'integer', required = T)
parser$add_argument('--min_pairs', type = 'integer', required = T,
                    default = 0)
parser$add_argument('--n_fractions', type = 'integer', required = T)
parser$add_argument('--n_replicates', type = 'integer', required = T)
parser$add_argument('--sample_idx', type = 'integer', default = 10)
args = parser$parse_args()
print(args)

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
output_filename = paste0(args$analysis, 
                         '-metric=', args$metric,
                         '-transform=', args$transform,
                         '-missing=', args$missing,
                         '-min_fractions=', args$min_fractions,
                         '-min_pairs=', args$min_pairs,
                         '-n_fractions=', args$n_fractions,
                         '-n_replicates=', args$n_replicates,
                         '-sample_idx=', args$sample_idx,
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
accession = basename(args$input_dir)
species = expts %>%
  filter(Accession == accession) %>%
  pull(Species) %>%
  unique()
if (length(species) != 1) 
  stop("couldn't get species for input directory: ", args$input_dir)

# list all _experiments_
experiments = list.dirs(args$input_dir, full.names = FALSE, recursive = FALSE)
## subset to those matching the pattern
if (args$pattern %in% c("BN_heavy", "BN_medium", "SEC_heavy", "SEC_medium")) {
  patt1 = gsub("_.*$", "", args$pattern)
  patt2 = gsub("^.*_", "", args$pattern)
  experiments %<>% extract(grepl(patt1, .) & grepl(patt2, .))
} else if (args$pattern == '2D_IEF') {
  experiments %<>% extract(grepl("2[dD]_IEF", .))
} else if (args$pattern != "all") {
  experiments %<>% extract(grepl(args$pattern, .))
}

# now, get files
input_dirs = file.path(args$input_dir, experiments)

# list all files
input_files = list.files(input_dirs[1], pattern = '*.rds', full.names = T) %>%
  extract(!grepl("metadata", .))
quant_modes = gsub("\\.rds$", "", basename(input_files))

# set up output data frame
results = data.frame()

# iterate over each quant. mode
for (quant_mode in quant_modes) {
  message("analyzing quant. mode #", which(quant_modes == quant_mode), " of ",
          length(quant_modes), ": ", quant_mode)
  
  # read inputs
  inputs = file.path(input_dirs, paste0(quant_mode, '.rds')) %>%
    map(readRDS) %>%
    # handle infinite or NaN values
    map(~ {
      mat = .
      mat[!is.finite(mat)] = NA
      mat[is.nan(mat)] = NA
      return(mat)
    })
  
  # read all metadata
  metas = file.path(input_dirs, 'metadata.rds') %>%
    map(readRDS) %>%
    map(`[`) ## remove attrs
  
  # split gene rows to match up across replicates
  exprs = map2(inputs, metas, ~ {
    protein_groups = dplyr::select(.y, `Majority protein IDs`) %>%
      rename(identifier = `Majority protein IDs`) %>%
      mutate(protein = strsplit(identifier, ';')) %>%
      unnest(protein)
    mat = matrix(NA, nrow = nrow(protein_groups), ncol = ncol(.x),
                 dimnames = list(protein_groups$protein, colnames(.x)))
    for (idx in seq_len(nrow(protein_groups))) {
      vec = .x[protein_groups$identifier[idx], ]
      mat[idx, ] = vec
    }
    # pre-filter any empty rows
    empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) == 0
    mat %<>% extract(!empty, )
    return(mat)
  })
  
  # seed the RNG, and sample replicates
  set.seed(args$sample_idx)
  exprs %<>% sample(args$n_replicates)
  
  # merge matrices
  proteins = map(exprs, rownames) %>% Reduce(intersect, .)
  expr = map(exprs, ~ extract(., proteins, )) %>% 
    do.call(cbind, .)
  
  # filter by n_fractions in this matrix
  n_fractions = rowSums(is.finite(expr) & expr > 0)
  expr %<>% extract(n_fractions >= args$min_fractions, )
  
  # pick one representative identifier for any 'duplicated' chromatogram
  # (arbitrarily, the least characters, to enrich for SwissProt)
  order = rownames(expr) %>% nchar() %>% order()
  expr %<>% extract(order, )
  ## now drop duplicated rows
  expr %<>% unique()
  
  # now that we have our combined matrix, we can sample fractions from it
  if (ncol(expr) >= args$n_fractions) {
    # set the seed again, and sample fractions
    set.seed(args$sample_idx)
    keep = sample(seq_len(ncol(expr)), args$n_fractions)
    expr %<>% extract(, keep)
    
    # store number of paired observations in the original matrix
    n = crossprod(is.finite(t(expr)) & t(expr) != 0)
    # also store number of fractions
    n_fractions = rowSums(is.finite(expr) & expr > 0)
    
    # catch infinite and NaN values
    expr[!is.finite(expr)] = NA
    expr[is.nan(expr)] = NA
    
    ## calculate matrix on the fly
    # step 1: impute missing values (optional)
    if (args$missing == 'zero') {
      expr[is.na(expr)] = 0
    } else if (args$missing == 'noise') {
      expr[expr == 0] = NA
      expr %<>% clean_profiles(impute_NA = TRUE, smooth = FALSE)
    } else if (args$missing == 'NA') {
      expr[expr == 0] = NA
    }
    
    # step 2: log-transform (optional)
    if (args$transform == "log") {
      expr %<>% log()
      if (args$missing == 'zero') {
        expr[expr == -Inf] = 0
      }
    }
    
    # step 3: quantile normalize (optional)
    if (args$transform == 'quantile') {
      dims = dimnames(expr)
      expr %<>% preprocessCore::normalize.quantiles()
      dimnames(expr) = dims
    }
    
    # step 4: score pairs
    cor = try(get_coexpr(t(expr), args$metric))
    
    # catch errors or exit
    if (class(cor) == 'try-error') {
      # allow only one particular exception - from treeClust
      if (grepl("Error in leaves", cor)) {
        next
      } else if (args$metric %in% c('rho_p', 'phi_s') & 
                 args$transform == 'log' &
                 quant_mode == 'ratio') {
        # ignore log-transformed ratios in proportionality
        next
      } else {
        stop(cor)
      }
    }
    
    # step 6: filter by minimum number of paired observations
    cor[n < args$min_pairs] = NA
    
    # finally, record the number of chromatograms and pairs
    n_chroms = sum(n_fractions > 0)
    n_pairs = sum(!is.na(cor))
    
    # run one of two analyses
    if (args$analysis == 'complexes') {
      # throw an error if complexes is being run on a non-human/mouse species
      if (!species %in% c("Homo sapiens", "Mus musculus")) {
        stop("can't run complex analysis on species: ", species, " (file: ",
             input_file, ")")
      }
      
      # read three sources of complexes
      complex_species = fct_recode(species,
                                   'human' = 'Homo sapiens',
                                   'mouse' = 'Mus musculus')
      complex_file = paste0("~/git/network-validation/data/complex/CORUM", 
                            "/complexes_", complex_species, ".txt")
      adj = read.delim(complex_file) %>%
        as_annotation_list(., 'gene_name', 'complex') %>%
        adjacency_matrix_from_list()
      
      # map proteins to genes
      map = metas %>%
        bind_rows() %>%
        dplyr::rename(gene = `Gene names`, protein = `Majority protein IDs`) %>%
        mutate(protein = strsplit(protein, ';')) %>%
        unnest(protein) %>%
        distinct(protein, gene) %$%
        setNames(gene, protein)
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
                         auroc = NA,
                         n_chroms = n_chroms,
                         n_pairs = n_pairs)
        results %<>% bind_rows(row)
      } else {
        auroc = auc(roc(predictions = x,
                        labels = factor(y, levels = c('0', '1'))))
        
        # append results
        row = data.frame(quant_mode = quant_mode,
                         auroc = auroc,
                         n_chroms = n_chroms,
                         n_pairs = n_pairs)
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
      meta = metas %>%
        bind_rows() %>%
        dplyr::rename(protein = `Majority protein IDs`) %>%
        mutate(protein = strsplit(protein, ';')) %>%
        unnest(protein) %>%
        distinct(protein)
      keep = intersect(rownames(expr), goa_slim$`DB Object ID`)
      cor %<>% extract(keep, keep)
      meta %<>% filter(protein %in% keep)
      
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
        true_positives = intersect(rownames(expr), targets)
        
        # create two vectors
        pairs0 = pairs %>%
          mutate(y = ifelse(gene1 %in% true_positives &
                              gene2 %in% true_positives,
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
                         n_annotated = length(true_positives),
                         auroc = auroc,
                         n_chroms = n_chroms,
                         n_pairs = n_pairs)
        results %<>% bind_rows(row)
      }
    }
  } else {
    # if there aren't enough fractions, ignore
    message("  warning: not enough fractions")
  }
}

# save results
saveRDS(results, output_file)
