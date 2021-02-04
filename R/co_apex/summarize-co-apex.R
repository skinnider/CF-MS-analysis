setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
source("R/functions/detect_system.R")

# set up IO
input_dir = file.path(base_dir, "co_apex")
output_dir = "data/analysis/co_apex"
if (!dir.exists(output_dir))
  dir.create(output_dir)

for (analysis in c("complexes", "GO")) {
  # list files
  files = list.files(input_dir, recursive = TRUE, pattern = analysis)
  
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
  
  # write 
  output_file = file.path(output_dir, paste0(analysis, '.rds'))
  saveRDS(results, output_file)
}
