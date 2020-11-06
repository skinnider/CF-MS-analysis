setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
}

# set up IO
input_dir = file.path(base_dir, "fractions_vs_replicates")
output_dir = "~/git/CF-MS-analysis/data/analysis/fractions_vs_replicates"
if (!dir.exists(output_dir))
  dir.create(output_dir)

for (analysis in c("complexes", "GO")) {
  # list files
  files = list.files(input_dir, recursive = T, pattern = analysis)
  
  # read each in turn and bind rows
  dat = files %>%
    file.path(input_dir, .) %>%
    map(readRDS) %>%
    setNames(files) %>%
    bind_rows(.id = 'file')
  
  # separate filename into columns
  results = dat %>%
    separate(file, into = c("accession", "pattern", "filename"),
             sep = '/') %>%
    separate(filename, into = c("analysis", "metric", "transform", "missing",
                                "min_fractions", "min_pairs", 
                                "n_fractions", "n_replicates", "sample_idx"),
             sep = '-') %>%
    mutate_if(is.character, ~ gsub("^.*=|\\.rds$", "", .))
  
  # write 
  output_file = file.path(output_dir, paste0(analysis, '.rds'))
  saveRDS(results, output_file)
}
