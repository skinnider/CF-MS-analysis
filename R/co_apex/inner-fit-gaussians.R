setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-fit-gaussians.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_file', type = 'character', required = T)
args = parser$parse_args()
print(args)

# set up output filepath
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

library(tidyverse)
library(magrittr)
library(PrInCE)
source("R/functions.R")

# read chromatogram matrix
mat = readRDS(args$input_file)
# handle infinite or NaN values
mat[!is.finite(mat)] = NA
mat[is.nan(mat)] = NA

# pre-filter any empty rows
empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) == 0
mat %<>% extract(!empty, )

# fit Gaussians
gauss = PrInCE::build_gaussians(mat, min_points = 1, min_consecutive = 1)

# save Gaussians to file
saveRDS(gauss, args$output_file)
