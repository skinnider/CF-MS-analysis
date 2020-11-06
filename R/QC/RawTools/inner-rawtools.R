setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-rawtools.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', 
                    default = "data/QC/RawTools")
parser$add_argument('--rawtools_dir', type = 'character', 
                    default = "~/RawTools-2.0.2")
args = parser$parse_args()

library(tidyverse)
library(magrittr)

# get instrument time from RawTools
cmd = paste("mono", file.path(args$rawtools_dir, "RawTools.exe"), 
            "-f", gsub("\\(", "\\\\(", 
                       gsub("\\)", "\\\\)",
                            gsub(" ", "\\\\ ", args$input_file))),
            "-o", args$output_dir, "-x")
out = system(cmd, intern = T)

# gzip the output file
output_file = paste0(file.path(args$output_dir, basename(args$input_file)),
                     '_Metrics.txt')
system(paste('gzip --force', shQuote(output_file)))