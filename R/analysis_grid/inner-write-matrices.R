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
                    choices = c('none', 'quantile', 'log'), required = T)
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
# pre-filter any empty rows
empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) == 0
mat %<>% extract(!empty, )

# catch a single infinite value in one chromatogram
mat[is.infinite(mat)] = NA
# catch 'NaN' values in iBAQ/ratio
mat[is.nan(mat)] = NA

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
