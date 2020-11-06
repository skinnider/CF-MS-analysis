# Plot Supplementary Figure 1:
#' a. bar chart, quantitation methods
#' b. protein 
#' protein groups CDF
#' CORUM coverage
#' abundance (Homo sapiens)
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

################################################################################
###### a. Quantitation methods
################################################################################

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  # recode quantitation modes slightly
  mutate(quant = fct_recode(Quantitation, 
                            'MS1 intensity' = 'Intensity',
                            'Spectral counts' = 'Spectral counts')) %>%
  # reorder the columns
  mutate(quant = fct_relevel(quant, 'Dimethyl', 'SILAC', 'MS1 intensity', 
                             'Spectral counts', 'iBAQ', 'LFQ')) 

# make a separate labels data frame
labs = expts %>%
  arrange(desc(quant)) %>%
  mutate(idx = row_number()) %>%
  group_by(quant) %>%
  summarise(mean_y = mean(idx)) %>%
  ungroup()

# plot
pal = jdb_palette("wolfgang_basic") %>% colorRampPalette() %>% 
  do.call(list(7)) %>% rev()
df = expts %>%
  dplyr::count(quant) %>%
  mutate(quant = reorder(quant, -n))
p1 = df %>%
  ggplot(aes(x = quant, y = n, fill = quant, color = quant)) +
  geom_col(width = 0.65, alpha = 0.8) +
  geom_text(aes(y = n + 12, label = n), size = 2, color = 'black') +
  scale_x_discrete() +
  scale_y_continuous('Experiments', limits = c(0, 140), expand = c(0, 0),
                     breaks = seq(0, 150, 30)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  boxed_theme(size_sm = 6, size_lg = 7) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p1

################################################################################
###### b. Protein groups CDF, as a % of proteome
################################################################################

# read protein quantitation CDF
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

# plot as a % of the proteome size
lab = filter(mean, n_fractions == 5) %>%
  mutate(text = paste0(round(100 * cov), '% of proteome\nin 5+ fractions'))
col = ggthemes::tableau_color_pal()(10)[2]
col = jdb_palette("wolfgang_basic") %>% colorRampPalette() %>% 
  do.call(list(7)) %>% extract(6)
p2 = best %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = coverage)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.2) + 
  # geom_hline(aes(yintercept = 0), linetype = 'dotted', color = 'black') +
  geom_path(data = mean, aes(y = cov), color = col, size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = cov), size = 2,
                   segment.size = 0.25, label.padding = 0.1, label.size = NA,
                   nudge_y = 0.9, nudge_x = +50, hjust = 0, fill = NA) +
  geom_point(data = lab, aes(y = cov), size = 0.8, color = col) +
  scale_x_continuous('Fractions', expand = c(0, 0),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous('% of proteome', labels = function(x) x * 100, 
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 50)) +
  boxed_theme()
p2

################################################################################
###### b. Protein groups CDF #2: datasets vs. proteins
################################################################################

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")

# read CDF
cdf = readRDS('data/QC/protein-groups-CDF.rds')

# reshape to # datasets/# of proteins
dataset_cdf = tidyr::crossing(n_proteins = seq(0, 8e3, 50),
                              n_fractions = c(seq_len(10), 15, 20, 25))
dataset_cdf$n_datasets = map2_int(
  dataset_cdf$n_proteins, dataset_cdf$n_fractions,
  ~ cdf %>%
    filter(n_fractions == .y) %>%
    filter(n_proteins > .x) %>%
    distinct(accession, experiment) %>%
    nrow()
)

# plot
pal = jdb_palette("wolfgang_basic") %>% colorRampPalette() %>% 
  do.call(list(8)) %>% rev()
p3 = dataset_cdf %>%
  filter(n_fractions %in% c(1, 5, 10, 15, 20, 25)) %>%
  ggplot(aes(x = n_proteins, y = n_datasets, color = factor(n_fractions))) +
  geom_path() + 
  scale_color_manual("Min. fractions", values = pal) +
  scale_x_continuous("Proteins") +
  scale_y_continuous("Datasets") +
  coord_cartesian(xlim = c(0, 8e3)) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              ncol = 1)) +
  boxed_theme() +
  theme(legend.position = 'right')
p3

################################################################################
###### e. Protein abundance, PaxDb (mouse)
################################################################################

coverage = readRDS("data/analysis/abundance/coverage.rds")

# group by proteins
proteins = coverage %>%
  group_by(Species, gene) %>%
  summarise(detected = any(n_fractions > 0),
            fractions = sum(n_fractions),
            experiments = n_distinct(paste(Accession, Replicate)[n_fractions > 0])) %>%
  ungroup()

# merge in protein abundance data
paxdb_files = list.files("data/resources/PaxDB", 
                         pattern = '*WHOLE_ORGANISM*', full.names = T)
paxdbs = map(paxdb_files, read.delim, comment.char = '#',
             header = F, col.names = c('paxid', 'id', 'abundance')) %>%
  setNames(c("Mus musculus", "Homo sapiens"))
id_files = paste0("~/git/network-validation/data/identifier/", 
                  c("MOUSE_10090", "HUMAN_9606"), "_idmapping.dat.gz")
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

# plot mouse only
fills = c('grey80', BuenColors::jdb_palette("wolfgang_extra")[8]) %>% rev()
colors = c('grey50', darken(BuenColors::jdb_palette("wolfgang_extra")[8])) %>% rev()
p4 = proteins %>%
  filter(Species == 'Mus musculus') %>%
  filter(abundance > 0) %>%
  ggplot(aes(x = abundance, fill = detected, color = detected)) +
  geom_histogram(position = 'identity', alpha = 0.5, #color = 'grey40',
                 size = 0.3) + 
  scale_x_log10("Protein abundance (ppm)", labels = fancy_scientific) +
  scale_y_continuous(expression(Proteins~(10^2)), expand = c(0, 0), 
                     limits = c(0, 1750), labels = function(x) x / 1e2) +
  scale_fill_manual('', values = fills, breaks = c("TRUE", "FALSE"),
                    labels = c("Detected by CF-MS", "Never detected")) +
  scale_color_manual('', values = colors,
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Detected by CF-MS", "Never detected")) +
  boxed_theme() +
  theme(legend.position = c(0.98, 1.04),
        legend.justification = c(1, 1),
        legend.key.size = unit(0.45, 'lines'))
p4

################################################################################
###### combine row 1
################################################################################

row1 = p1 + p2 + p3 + p4 + plot_layout(ncol = 4, widths = c(0.7, 1, 0.93, 1))
row1
ggsave("fig/final/supp1/row1.pdf", row1, width = 18.5, height = 5.25, 
       units = "cm", useDingbats = F)

################################################################################
###### e. GO enrichment (CORUM, never detected)
################################################################################

enrs = readRDS("data/analysis/complexes/GO-never-detected.rds")

# plot
pal = get_pal("city haze")
enr = enrs %>%
  bind_rows() %>%
  filter(conditional == TRUE)

# plot a nice version 
enr %>%
  filter(## these were mistakenly flipped
         test == 'never detected') %>%
  group_by(test) %>%
  arrange(Pvalue) %>%
  pull(Term) %>%
  head(100)
terms_detected = c('intracellular organelle',
                   'RNA binding',
                   'cytoplasm',
                   'protein-containing complex',
                   'nuclear part',
                   'intracellular protein transport',
                   # 'focal adhesion',
                   'mitochondrion',
                   'cellular protein localization',
                   'cytoskeletal part')
enr %>%
  filter(## these were mistakenly flipped
         test != 'never detected') %>%
  group_by(test) %>%
  arrange(Pvalue) %>%
  pull(Term) %>%
  head(100)
terms_not = c('signaling receptor activity',
              'integral component of membrane',
              'G protein-coupled receptor activity',
              'ion channel activity',
              'neurotransmitter receptor activity',
              'cytokine activity',
              'plasma membrane part',
              'regulation of cell proliferation',
              'cell communication')

pal = get_pal("early meadow")[c(1, 3)]
p5a = enr %>%
  filter(test == 'never detected') %>%
  filter(Term %in% terms_detected) %>%
  ggplot(aes(x = reorder(Term, -Pvalue), y = -log10(Pvalue))) + 
  facet_grid(~ "Detected") +
  geom_col(width = 0.65, alpha = 0.7, size = 0.3, aes(color = '1', fill = '1')) +
  geom_hline(aes(yintercept = -log10(0.05)), color = 'red', size = 0.4) +
  scale_fill_manual(values = pal[2], guide = F) +
  scale_color_manual(values = pal[2], guide = F) +
  scale_x_reordered() +
  scale_y_continuous(expression(-log[10](P)), limits = c(0, 40), 
                     expand = c(0, 0)) +
  coord_flip() +
  clean_theme() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p5a  
p5b = enr %>%
  filter(test == 'detected') %>%
  filter(Term %in% terms_not) %>%
  ggplot(aes(x = reorder(Term, -Pvalue), y = -log10(Pvalue))) + 
  facet_grid(~ "Never detected") +
  geom_col(width = 0.65, alpha = 0.7, size = 0.3, aes(color = '1', fill = '1')) +
  geom_hline(aes(yintercept = -log10(0.05)), color = 'red', size = 0.4) +
  scale_fill_manual(values = pal[1], guide = F) +
  scale_color_manual(values = pal[1], guide = F) +
  scale_x_reordered() +
  scale_y_continuous(expression(-log[10](P)), limits = c(0, 40), 
                     expand = c(0, 0)) +
  coord_flip() +
  clean_theme() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p5b

# combine
p5 = p5b | p5a
p5
ggsave("fig/final/supp1/row2a.pdf", p5,
       width = 14, height = 4.75, units = "cm", useDingbats = F)

################################################################################
###### f. change over original CDF
################################################################################

# read protein quantitation CDF
cdf1 = readRDS("data/QC/protein-groups-CDF.rds")
# read supplement CDF
# now, generate an identical CDF for the supplements
cdf2 = readRDS('data/QC/supplements-CDF.rds')

# take best number of proteins per fraction in the original data
best = cdf1 %>%
  group_by(accession, experiment, n_fractions) %>%
  summarise(n_proteins = max(n_proteins)) %>%
  ungroup()

# calculate difference relative to original
delta = best %>%
  left_join(cdf2, by = c('accession', 'experiment', 'n_fractions')) %>%
  mutate(delta = n_proteins.x - n_proteins.y) %>%
  # filter to a maximum of 50 fractions
  filter(between(n_fractions, 1, 50)) 

# remove a couple experiments with FDR != 1%
delta %<>%
  filter(accession != 'PXD002640')

# calculate the mean curve
mean = delta %>%
  group_by(n_fractions) %>%
  summarise(mean = mean(delta, na.rm = T),
            sd = sd(delta, na.rm = T),
            median = median(delta, na.rm = T),
            n = n()) %>%
  ungroup()

# plot each curve, plus the mean curve
## label n=10 fractions
lab = filter(mean, n_fractions == 5) %>%
  mutate(text = paste0('+', round(mean), 
                       ' protein groups\nin 5+ fractions'))

# add percentage of the original to the label
delta2 = best %>%
  left_join(cdf2, by = c('accession', 'experiment', 'n_fractions')) %>%
  mutate(delta = n_proteins.x - n_proteins.y,
         delta_pct = delta / n_proteins.y) %>%
  # filter to a maximum of 50 fractions
  filter(between(n_fractions, 1, 50))
mean2 = delta2 %>%
  group_by(n_fractions) %>%
  filter(is.finite(delta_pct)) %>%
  summarise(mean = mean(delta_pct, na.rm = T),
            median = median(delta_pct, na.rm = T),
            sd = sd(delta_pct, na.rm = T),
            n = n()) %>%
  ungroup()
lab2 = filter(mean2, n_fractions == 5) %>%
  mutate(text = paste0(round(100 * mean), '% more protein groups\n',
                       'in 10+ fractions'))
lab$text %<>% paste0(' (+', round(100 * lab2$mean), '%)')

col = jdb_palette("wolfgang_basic") %>% colorRampPalette() %>% 
  do.call(list(7)) %>% extract(5)
p6 = delta %>%
  unite(group, accession, experiment, remove = FALSE) %>%
  ggplot(aes(x = n_fractions, y = delta)) +
  geom_path(aes(group = group), color = 'grey90', size = 0.2) + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted', color = 'black') +
  geom_path(data = mean, aes(y = mean), color = col, size = 0.75) +
  geom_label_repel(data = lab, aes(label = text, y = mean), size = 2,
                   segment.size = 0.25, label.padding = 0.2, label.size = NA,
                   nudge_y = -5000, nudge_x = +50, hjust = 0) +
  geom_point(data = lab, aes(y = mean), size = 0.8, color = col) +
  scale_x_continuous('Fractions', expand = c(0, 0),
                     breaks = c(1, seq(10, 50, 10)),
                     labels = c(1, seq(10, 40, 10), '>50')) +
  scale_y_continuous(expression(Delta(protein~groups))) +
  boxed_theme()
p6

################################################################################
###### combine row 2
################################################################################

# patchwork doesn't really work, so we have to manually save them
row2a = p5a + p5b
# row2a
ggsave("fig/final/supp1/row2-1.pdf", row2a, width = 13, height = 4.7,
       units = "cm", useDingbats = F)
ggsave("fig/final/supp1/row2-2.pdf", p6, width = 4.8, height = 4.2,
       units = "cm", useDingbats = F)

## layout reveals we also need to break the layout a bit for panel a
ggsave("fig/final/supp1/row1-1.pdf", p1, width = 3.8, height = 4.2,
       units = "cm", useDingbats = F)
