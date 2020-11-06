setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
source("R/functions/detect_system.R")
source("R/functions/load_networks.R")

# set up coexpression files
coexpr_names = c('protein/Kustatscher2019', 'protein/Lapek2017')
coexpr_files = file.path('~/git/network-validation/data/coexpression',
                         paste0(coexpr_names, '.txt.gz'))

# set up colocalization files
coloc_names = c('Geladaki2019', 'Orre2019')
coloc_files = file.path('~/git/network-validation/data/colocalization',
                        paste0(coloc_names, '.txt.gz'))

# read all matrices
mats1 = map(coexpr_files, ~ read.delim(.) %>% as.matrix()) %>%
  setNames(paste0('coexpression|', coexpr_names))
mats2 = map(coloc_files, ~ read.delim(.) %>% as.matrix()) %>%
  setNames(paste0('colocalization|', coloc_names))
mats = c(mats1, mats2)

# calculate coexpression
missing = map(mats, ~ !is.finite(.) | . == 0)
cors = map(mats, ~ cor(., use = 'p'))

# read all networks
networks = load_networks()

# set up analysis grid
grid = tidyr::crossing(network = names(networks),
                       matrix = names(mats))

# set up results containers
pct_neg = data.frame()
rnd = data.frame()
stats = data.frame()

# iterate through the grid
for (grid_idx in seq_len(nrow(grid))) {
  network_name = grid$network[grid_idx]
  matrix_name = grid$matrix[grid_idx]
  message("[", grid_idx, "/", nrow(grid), "] analyzing network: ",
          network_name, " / matrix: ", matrix_name, " ...")
  
  network = networks[[network_name]]
  cor = cors[[matrix_name]]
  
  # tag 
  colnames(network)[c(1, 2)] = c('protein_A', 'protein_B')
  net0 = network %>%
    filter(protein_A %in% rownames(cor), protein_B %in% rownames(cor))
  # coerce method to character
  if ("method" %in% colnames(net0))
    net0 %<>% mutate(method = as.character(method))
  net0$cor = cor[as.matrix(net0[, 1:2])]
  net0 %<>% drop_na(cor)
  
  ## 1. percent negatively coexpressed
  df1 = net0 %>%
    dplyr::summarise(pct = mean(cor < 0)) %>%
    mutate(network = network_name,
           matrix = matrix_name)
  pct_neg %<>% bind_rows(df1)
  
  ## 2. random sample of 1000 values ("spectrum")
  set.seed(0)
  sample_size = ifelse(nrow(net0) > 1e3, 1e3, nrow(net0))
  df2 = net0 %>%
    sample_n(sample_size) %>%
    mutate(network = network_name,
           matrix = matrix_name)
  rnd %<>% bind_rows(df2)
  
  ## 3. boxplot stats, mean, and s.d.
  df3 = net0 %>%
    do(tidy(summary(.$cor))) %>%
    mutate(sd = sd(net0$cor)) %>%
    mutate(network = network_name,
           matrix = matrix_name)
  stats %<>% bind_rows(df3)
}

# combine 
output = list('% negatively correlated' = pct_neg,
              'random sample' = rnd,
              'summary statistics' = stats)

# save
output_file = file.path(args$output_dir, 'continuous_outcomes.rds')
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = TRUE)
saveRDS(output, output_file)
