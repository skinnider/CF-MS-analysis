# Write species with 3+ CF-MS experiments to a line-delimited file for input to
# TimeTree.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read metadata
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")

# filter to species with at least three supporting experiments
species = expts %>%
  dplyr::count(Species) %>%
  filter(n >= 3) %>%
  arrange(desc(n))

# write these species as a line-delimited file
lines = species$Species
writeLines(lines, "data/resources/TimeTree/species.txt")

# upload to timetree.org
# export Newick file
# download