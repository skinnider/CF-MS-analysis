setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# detect system
args = list(); source("R/functions/detect_system.R")

# list Gaussian files
gauss_dir = file.path(base_dir, "fit_gaussians")
gauss_files = list.files(gauss_dir, recursive = TRUE, full.names = TRUE,
                         pattern = '*.rds')

# map over Gaussian files
all_gaussians = map_dfr(gauss_files, ~ {
  split = strsplit(.x, '/') %>% unlist()
  experiment = split[length(split) - 3]
  replicate = split[length(split) - 2]
  quant_mode = last(split) %>% gsub("^.*-|\\.rds", "", .)
  
  gauss = readRDS(.x)
  # n = map_int(gauss, 'n_gaussians')
  heights = map(gauss, ~ .$coefs$A)
  
  # create data frame
  df = data.frame(experiment = experiment, 
                  replicate = replicate, 
                  quant_mode = quant_mode,
                  protein = rep(names(heights), lengths(heights)),
                  height = unlist(heights)) %>%
    # normalize height as a % of max
    group_by(protein) %>%
    mutate(height_norm = height / max(height)) %>%
    ungroup()
  df
})

# save to git
saveRDS(all_gaussians, "data/analysis/co_apex/all_gaussians.rds")
