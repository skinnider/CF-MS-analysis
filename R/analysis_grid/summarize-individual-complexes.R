setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
source("R/functions/detect_system.R")

# set up IO
input_dir = file.path(base_dir, "analysis_grid")
output_dir = "~/git/CF-MS-analysis/data/analysis/analysis_grid"
if (!dir.exists(output_dir))
  dir.create(output_dir)

# list files
files = list.files(input_dir, recursive = T, pattern = "individual_complexes")

# read each in turn and bind rows
dat = files %>%
  file.path(input_dir, .) %>%
  map(readRDS) %>%
  setNames(files) %>%
  bind_rows(.id = 'file')

# separate filename into columns
results = dat %>%
  separate(file, into = c("accession", "experiment", "filename"),
           sep = '/') %>%
  separate(filename, into = c("analysis", "metric", "transform", "missing"),
           sep = '-') %>%
  mutate_if(is.character, ~ gsub("^.*=|\\.rds$", "", .))

# keep only iBAQ
results %<>% filter(quant_mode == 'iBAQ') %>% dplyr::select(-quant_mode)

# write 
output_file = file.path(output_dir, 'individual_complexes.rds')
saveRDS(results, output_file)
