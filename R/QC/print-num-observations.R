# Print the total number of protein quantitations across all fractions.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read protein quantitation CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds")

# print total number of observations
cdf %>%
  group_by(quant_mode, search) %>%
  summarise(n_quants = sum(n_proteins)) %>%
  ungroup() %>%
  arrange(desc(n_quants))
