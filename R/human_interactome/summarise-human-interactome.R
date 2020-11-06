setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
}

# list files
input_dir = file.path(base_dir, "human_interactome", "PR")
files = list.files(input_dir, recursive = T, pattern = 'rds', full.names = TRUE)

# read PR curves
dats = map(files, readRDS)
curves = map(dats, 'PR') %>%
  bind_rows()

# save
output_dir = "data/analysis/human_interactome"
if (!dir.exists(output_dir))
  dir.create(output_dir)
saveRDS(curves, file.path(output_dir, 'PR.rds'))
