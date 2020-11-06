# Plot a CDF of protein groups vs. # of fractions in each experiment and 
# on average.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")

# read CDF file, or generate it if it doesn't exist
cdf_file = 'data/QC/protein-groups-CDF.rds'
if (file.exists(cdf_file)) {
  cdf = readRDS(cdf_file)
} else {
  # list all chromatogram files
  files = list.files('~/projects/rrg-ljfoster-ab/skinnim/CF-MS/chromatograms', 
                     pattern = '*.rds', full.names = T, recursive = T) %>%
    # ignore the metadata files
    extract(!grepl("metadata", .))
  
  # extract accession/experiment/search/quantitation
  split = strsplit(gsub("^.*CF-MS/", "", files), '/')
  accessions = map_chr(split, 3)
  experiments = map_chr(split, 4)
  searches = map_chr(split, 5)
  quant_modes = gsub("\\..*$", "", basename(files))
  
  # read all FASTA files
  fasta_files = list.files("~/git/CF-MS-searches/data/fasta/filtered", 
                           pattern = '*.gz', full.names = T)
  patts = with(expts, paste0(substr(Species, 1, 1), '.',
                             gsub("^.* ", "", Species))) %>%
    unique() %>%
    fct_recode('Cyanothece' = 'C.51142',
               'B.oleracea' = 'B.italica') %>%
    paste0(collapse = '|')
  fasta_files %<>% extract(grepl(patts, .))
  fastas = map(fasta_files, seqinr::read.fasta) %>%
    setNames(gsub("^.*-|\\.fasta.*$", "", basename(fasta_files)))
  proteome_sizes = map_int(fastas, length)
  
  # construct CDF for each file
  cdf = data.frame()
  for (file_idx in seq_along(files)) {
    file = files[file_idx]
    accession = accessions[file_idx]
    experiment = experiments[file_idx]
    search = searches[file_idx]
    quant_mode = quant_modes[file_idx]
    message('[', file_idx, '/', length(files), '] working on file: ', file)
    
    # extract CDF for this file
    mat = readRDS(file)
    cdf0 = data.frame()
    for (n_fractions in seq(0, ncol(mat))) {
      mat0 = mat %>%
        replace(is.na(.), 0) %>%
        extract(rowSums(. > 0) >= n_fractions, , drop = F)
      n_proteins = nrow(mat0)
      
      # also calculate coverage as a % of proteome size
      if (nrow(mat0) > 0) {
        protein_groups = mat0 %>%
          rownames() %>%
          strsplit(';') %>%
          unlist() %>%
          unique()
        species = expts %>%
          filter(Accession == accession, Replicate == experiment) %>%
          pull(Species) 
        key = paste0(substr(species, 1, 1), '.', gsub("^.* ", "", species)) %>%
          fct_recode('Cyanothece' = 'C.51142',
                     'B.oleracea' = 'B.italica') %>%
          as.character()
        proteome_size = proteome_sizes[key]
        coverage = length(protein_groups) / proteome_size
      } else {
        coverage = 0
      }
      
      # append
      cdf0 %<>% bind_rows(data.frame(n_fractions = n_fractions,
                                     n_proteins = n_proteins,
                                     coverage = coverage))
    }
    
    # add the rest of the metadata
    cdf0 %<>% mutate(accession = accession, experiment = experiment, 
                     search = search, quant_mode = quant_mode)
    
    # append to master CDF
    cdf %<>% bind_rows(cdf0)
  }
  
  # rearrange columns
  cdf %<>% 
    dplyr::select(accession, experiment, quant_mode, search,
                  n_fractions, n_proteins, coverage) %>%
    # clean up quantitation mode
    mutate(quant_mode = chartr('_', ' ', quant_mode))
  
  # save 
  saveRDS(cdf, cdf_file)
}

# take best number of proteins per fraction in the original data
best = cdf %>%
  group_by(accession, experiment, n_fractions) %>%
  summarise(n_proteins = max(n_proteins), coverage = max(coverage)) %>%
  ungroup()

# calculate the mean curve
mean = best %>%
  group_by(n_fractions) %>%
  summarise(mean = mean(n_proteins, na.rm = T),
            sd = sd(n_proteins, na.rm = T),
            median = median(n_proteins, na.rm = T),
            q1 = quantile(n_proteins, na.rm = T, probs = 0.25),
            q3 = quantile(n_proteins, na.rm = T, probs = 0.75),
            cov = mean(coverage, na.rm = T),
            n = n()) %>%
  ungroup()

# plot each curve, plus the mean curve
## label n=10 fractions
lab = filter(mean, n_fractions == 5) %>%
  mutate(text = paste0(format(round(mean), big.mark = ','), 
                       ' protein groups\nin 5+ fractions'))
col = BuenColors::jdb_palette("solar_extra") %>% extract(1)
p1 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = n_proteins)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.15) + 
  # geom_hline(aes(yintercept = 0), linetype = 'dotted', color = 'black') +
  geom_path(data = mean, aes(y = mean), color = col, size = 0.75) +
  # geom_path(data = mean, aes(y = q1), color = col, size = 0.35,
  #           linetype = 'dashed') +
  # geom_path(data = mean, aes(y = q3), color = col, size = 0.35,
  #           linetype = 'dashed') +
  geom_label_repel(data = lab, aes(label = text, y = mean), size = 2,
                   segment.size = 0.25, label.padding = 0.1, label.size = NA,
                   nudge_y = +4000, nudge_x = +50, hjust = 0, fill = NA) +
  geom_point(data = lab, aes(y = mean), size = 0.8, color = col) +
  scale_x_continuous('Fractions', expand = c(0, 0), limits = c(1, 50),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous(expression(Protein~groups~(10^3)), limits = c(0, 6400), 
                     expand = c(0, 0), labels = function(x) x / 1e3) +
  coord_cartesian(xlim = c(1, 50)) +
  boxed_theme()
p1
ggsave("fig/QC/protein-groups-cdf.pdf", p1,
       width = 4.6, height = 4.2, units = "cm", useDingbats = F)

# also plot as a % of the proteome size
lab = filter(mean, n_fractions == 5) %>%
  mutate(text = paste0(round(100 * cov), '% of proteome\nin 5+ fractions'))
# col = ggthemes::tableau_color_pal()(10)[2]
p2 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = coverage)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.2) + 
  # geom_hline(aes(yintercept = 0), linetype = 'dotted', color = 'black') +
  geom_path(data = mean, aes(y = cov), color = col, size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = cov), size = 2,
                   segment.size = 0.25, label.padding = 0.4, label.size = NA,
                   nudge_y = 0.8, nudge_x = +50, hjust = 0, fill = NA) +
  scale_x_continuous('Fractions', expand = c(0, 0),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous('% of proteome', labels = function(x) x * 100, 
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 50)) +
  boxed_theme()
p2
ggsave("fig/QC/protein-groups-cdf-percent.pdf", p2,
       width = 4.6, height = 4.2, units = "cm", useDingbats = F)
