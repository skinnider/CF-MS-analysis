setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# what system are we on?
system = 'cedar'
base_dir = "~skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
if (!dir.exists(base_dir)) {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
  system = 'sockeye'
}

# set up IO
input_dir = file.path(base_dir, "sliding_window")
output_dir = "~/git/CF-MS-analysis/data/analysis/sliding_window"
if (!dir.exists(output_dir))
  dir.create(output_dir)

# complexes only 
files = list.files(input_dir, recursive = T, pattern = "complexds")

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
  separate(filename, into = c("analysis", "metric", "transform", "missing", 
                              "window_size"),
           sep = '-') %>%
  mutate_if(is.character, ~ gsub("^.*=|\\.rds$", "", .))

# write 
output_file = file.path(output_dir, 'complexes.rds')
saveRDS(results, output_file)
