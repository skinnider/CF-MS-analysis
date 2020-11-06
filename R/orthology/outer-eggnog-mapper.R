setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-eggnog-mapper.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# check that all of the dependent directories exist
if (!dir.exists("~/git/CF-MS-searches"))
  stop("repository `CF-MS-searches` does not exist")

# detect system
source("R/functions/detect_system.R")

# list the input files
input_files = list.files("~/git/CF-MS-searches/data/fasta/filtered",
                         pattern = '*\\.fasta$', full.names = TRUE)
grid = data.frame(input_file = input_files) 

# check which files exist
grid0 = grid %>%
  mutate(output_dir = file.path(base_dir, 'eggNOG'),
         output_file = basename(input_files) %>% 
           gsub("\\.fasta$", "", .) %>%
           paste0(., '.emapper.annotations') %>%
           file.path(output_dir, .)) %>%
  filter(!file.exists(output_file)) %>%
  dplyr::select(input_file)

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/grids/eggnog_mapper.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = ifelse(system == 'cedar',
                " ~/git/CF-MS-analysis/sh/orthology/eggnog-mapper.sh",
                " ~/git/CF-MS-analysis/sh/orthology/eggnog-mapper.torque.sh")
submit_job(grid0, script, args$allocation, system)
