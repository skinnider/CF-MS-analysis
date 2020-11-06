# Generate a 'metrics file' that summarizes the experiment time associated with
# each Thermo raw file in the combined analysis.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-rawtools.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
parser$add_argument('--base_dir', type = 'character', 
                    default = "~/projects/rrg-ljfoster-ab/skinnim/CF-MS")
parser$add_argument('--output_dir', type = 'character', 
                    default = "data/QC/RawTools")
args = parser$parse_args()

# create the directory, if it doesn't exist
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir)

library(tidyverse)
library(magrittr)

# first, read list of all files
files = read.csv("~/git/CF-MS-searches/data/files.csv") %>%
  # filter to Thermo raw files
  filter(endsWith(tolower(File), '.raw')) %>%
  # expand full filepath
  mutate(filepath = file.path(args$base_dir, Accession, 
                              gsub("^.*\\\\", "", File)))

# make experiment grid; filter files that already exist
grid = files %>%
  mutate(input_exists = file.exists(filepath),
         output_file = file.path(args$output_dir, paste0(File,
                                                         '_Metrics.txt.gz')), 
         output_exists = file.exists(output_file)) %>%
  filter(input_exists, !output_exists) %>%
  dplyr::select(filepath)

# write the grid
grid_file = "sh/QC/grids/rawtools.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# finally, call the slurm script
n_jobs = nrow(grid)
## only 1,000 jobs per array
max_jobs = 1000
n_submissions = ifelse(n_jobs > max_jobs, ceiling(n_jobs / max_jobs), 1)
for (submission_idx in seq_len(n_submissions)) {
  job_start = (submission_idx - 1) * max_jobs + 1
  job_end = ifelse(submission_idx == n_submissions,
                   ifelse(n_jobs %% max_jobs == 0,
                          submission_idx * max_jobs,
                          job_start - 1 + n_jobs %% max_jobs),
                   submission_idx * max_jobs)
  system(paste0(
    "cd ~/project; ",
    "sbatch --account=", args$allocation, " --array=", job_start, "-", job_end,
    " ~/git/CF-MS-analysis/sh/QC/run_rawtools.sh"
  ))
}
