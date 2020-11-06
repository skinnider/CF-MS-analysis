# Plot the second row of Figure 1:
#' fractionation techniques
#' protein groups CDF
#' CORUM coverage
#' abundance (Homo sapiens)
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

################################################################################
###### b. Fractionation methods
################################################################################

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  # recode fractionation slightly
  mutate(fract = fct_recode(Fractionation, 
                            'PAGE' = 'BN-PAGE',
                            'PAGE' = 'lpBN-PAGE',
                            'PAGE' = 'CN-PAGE',
                            'IEX' = 'SAX',
                            'IEX' = 'WWC'),
         fract = gsub("^XL-", "", fract)) %>%
  # reorder the columns
  mutate(quant = fct_relevel(Quantitation, 'Dimethyl', 'SILAC',
                             'MS1\nintensity', 'Spectral\ncounts',
                             'iBAQ', 'LFQ')) 

# make a separate labels data frame
labs = expts %>%
  arrange(desc(fract)) %>%
  mutate(idx = row_number()) %>%
  group_by(fract) %>%
  summarise(mean_y = mean(idx)) %>%
  ungroup()

# plot
pal = BuenColors::jdb_palette("solar_extra") %>% extract(-c(1, 9)) %>%
  colorRampPalette() %>% do.call(list(5))
df = expts %>%
  dplyr::count(fract) %>%
  mutate(fract = reorder(fract, -n))
p1 = df %>%
  ggplot(aes(x = fract, y = n, fill = fract, color = fract)) +
  geom_col(width = 0.65, alpha = 0.8) +
  geom_text(aes(y = n + 6, label = n), size = 2, color = 'black') +
  scale_x_discrete() +
  scale_y_continuous('Experiments', limits = c(0, 100), expand = c(0, 0),
                     breaks = seq(0, 90, 30)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p1

################################################################################
###### c. Protein groups CDF
################################################################################

cdf = readRDS("data/QC/protein-groups-CDF.rds")
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
## label n=5 fractions ("good quality" chromatogram)
lab = filter(mean, n_fractions == 5) %>%
  mutate(text = paste0(format(round(mean), big.mark = ','), 
                       ' protein groups\nin 5+ fractions'))
col = BuenColors::jdb_palette("solar_extra") %>% extract(1)
p2 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = n_proteins)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.15) + 
  geom_path(data = mean, aes(y = mean), color = col, size = 0.75) +
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
  boxed_theme(size_sm = 6, size_lg = 7) 
p2

################################################################################
###### d. CORUM coverage
################################################################################

cdf = readRDS("data/analysis/complexes/complex-cdf.rds")
pal = BuenColors::jdb_palette("solar_extra") %>% colorRampPalette()
p3 = cdf %>%
  filter(n_expts %in% c(1, 2, 3, 5, 7, 10, 20, 40, 60)) %>%
  ggplot(aes(x = min_fractions, y = mean, color = factor(n_expts))) +
  geom_path() +
  scale_x_continuous('Fractions', limits = c(1, NA),
                     breaks = c(1, seq(10, 50, 10))) +
  scale_y_continuous('Coverage (%)', labels = function(x) x * 100,
                     limits = c(0, NA)) +
  scale_color_manual('Experiments', values = pal(8)) +
  guides(color = guide_legend(ncol = 1, title.position = 'top',
                              title.hjust = 0,
                              override.aes = list(size = 0.5))) +
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.key.size = unit(0.5, 'lines'),
        legend.position = 'right',
        legend.title = element_text(size = 6))
p3

################################################################################
###### e. Protein abundance, PaxDb (human)
################################################################################

coverage = readRDS("data/analysis/abundance/coverage.rds")

# group by proteins
proteins = coverage %>%
  group_by(Species, gene) %>%
  summarise(detected = any(n_fractions > 0),
            fractions = sum(n_fractions),
            experiments = n_distinct(paste(Accession, Replicate)[
              n_fractions > 0])) %>%
  ungroup()

# merge in protein abundance data
paxdb_files = list.files("data/resources/PaxDB", 
                         pattern = '*WHOLE_ORGANISM*', full.names = T)
paxdbs = map(paxdb_files, read.delim, comment.char = '#',
             header = F, col.names = c('paxid', 'id', 'abundance')) %>%
  setNames(c("Mus musculus", "Homo sapiens"))
id_files = paste0("~/git/network-validation/data/identifier/", 
                  c("MOUSE_10090", "HUMAN_9606"), 
                  "_idmapping.dat.gz")
maps = map(id_files, read_tsv, col_names = c("uniprot", "db", "id")) %>%
  setNames(c("Mus musculus", "Homo sapiens"))
abundance = map2_dfr(paxdbs, maps, ~ {
  mutate(.x, id = gsub("^(10090|9606)\\.", "", id)) %>%
    left_join(.y, by = 'id') %>%
    drop_na(uniprot, abundance) %>%
    # keep only one protein ID per gene
    group_by(paxid) %>%
    arrange(desc(abundance)) %>%
    dplyr::slice(1) %>%
    ungroup() %>% 
    distinct(uniprot, abundance)
}) %>%
  dplyr::rename(gene = uniprot)
proteins %<>% left_join(abundance, by = 'gene')

# plot human only
prot1 = filter(proteins, Species == "Homo sapiens")
prot1 %>% filter(fractions > 0) %$%
  quantile(fractions, probs = seq(0, 1, 0.05))

fills = c('grey80', BuenColors::jdb_palette("solar_extra")[1]) %>% rev()
colors = c('grey50', darken(BuenColors::jdb_palette("solar_extra")[1])) %>% rev()
p4 = prot1 %>%
  filter(abundance > 0) %>%
  ggplot(aes(x = abundance, fill = detected, color = detected)) +
  geom_histogram(position = 'identity', alpha = 0.5, #color = 'grey40',
                 size = 0.3) + 
  scale_x_log10("Protein abundance (ppm)", labels = fancy_scientific) +
  scale_y_continuous(expression(Proteins~(10^2)), labels = function(x) x / 1e2,
                     expand = c(0, 0), limits = c(0, 1750)) +
  scale_fill_manual('', values = fills, breaks = c("TRUE", "FALSE"),
                    labels = c("Detected by CF-MS", "Never detected")) +
  scale_color_manual('', values = colors,
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Detected by CF-MS", "Never detected")) +
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.position = c(0.98, 1.04),
        legend.justification = c(1, 1),
        legend.key.size = unit(0.45, 'lines'),
        axis.line = element_blank())
p4

################################################################################
###### combine
################################################################################

p = p1 + p2 + p3 + p4 + plot_layout(ncol = 4, widths = c(0.7, 1, 0.85, 1))
p 
ggsave("fig/final/fig1/row2.pdf", p, width = 18.5, height = 4.93, units = "cm",
       useDingbats = F)
