# Collate all protein complex and GO AUCs from the min. pairs experiment
# into a single file.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
  system = 'sockeye'
}

# set up IO
input_dir = file.path(base_dir, "min_pairs")
output_dir = "~/git/CF-MS-analysis/data/analysis/min_pairs"
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
    separate(file, into = c("accession", "experiment", "filename"),
             sep = '/') %>%
    separate(filename, into = c("analysis", "metric", "transform",
                                "min_pairs", "missing", "min_fractions"),
             sep = '-') %>%
    mutate_if(is.character, ~ gsub("^.*=|\\.rds$", "", .)) %>%
    type_convert()
  
  # write 
  output_file = file.path(output_dir, paste0(analysis, '.rds'))
  saveRDS(results, output_file)
}
