# Summarize key statistics about the consistently processed dataset.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

# list all summary files
files = list.files('~/projects/rrg-ljfoster-ab/skinnim/CF-MS/maxquant',
                   pattern = 'summary.txt', full.names = T, recursive = T)

# loop over files
summary = map_dfr(files, ~ {
  tab = fread(.)
  
  ## number of input raw files
  n_files = n_distinct(tab$`Raw file`)
  
  ## number of fractions
  n_fractions = n_distinct(tab$Experiment)
  
  ## number of MS/MS spectra
  msms = sum(tab$`MS/MS`)
  
  ## number of sequenced peptides
  peptides = sum(tab$`MS/MS Identified`)
  
  data.frame(file = ., raw_files = n_files, fractions = n_fractions,
             msms = msms, peptides = peptides)
})

# flag accession, replicate, and search
summary %<>%
  mutate(file = gsub("^.*maxquant/", "", file)) %>%
  separate(file, c("accession", "replicate", "search", "x"), '/') %>%
  dplyr::select(-x)

# print some statistics
message("number of experiments searched: ", 
        summary %>% unite(experiment, accession, replicate) %>% 
          pull(experiment) %>% n_distinct())
message("number of raw files searched: ", summary %>% pull(raw_files) %>% sum())
message("number of fractions searched: ", summary %>% pull(fractions) %>% sum())
message("number of MS/MS spectra analyzed: ", summary %>% pull(msms) %>% sum())
message("number of peptides sequenced: ", summary %>% pull(peptides) %>% sum())
