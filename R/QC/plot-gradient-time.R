setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(broom)
library(data.table)
source("R/theme.R")

# first, get a list of all human datasets
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species == "Homo sapiens")

# read human protein abundance from PaxDB
paxdb = read.delim("data/resources/PaxDB/9606-WHOLE_ORGANISM-integrated.txt.gz",
                   comment.char = '#', header = FALSE,
                   col.names = c('paxid', 'id', 'abundance'))

# map PaxDb identifiers to gene names
map = read_tsv(
  "~/git/network-validation/data/identifier/HUMAN_9606_idmapping.dat.gz", 
  col_names = c("uniprot", "db", "id"))
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
ensp = filter(map, db == 'Ensembl_PRO') %>% dplyr::select(-db)
ensp2gn = left_join(gn, ensp, by = 'uniprot') %>%
  dplyr::select(-uniprot) %>%
  set_colnames(c('gene_name', 'id')) %>%
  drop_na()
paxdb %<>% 
  mutate(id = gsub("^.*\\.", "", id)) %>%
  left_join(ensp2gn, by = 'id') %>%
  drop_na(gene_name, abundance) %>%
  # keep only one protein ID per gene
  group_by(gene_name) %>%
  mutate(n = n()) %>%
  arrange(desc(abundance)) %>%
  dplyr::slice(1) %>%
  ungroup() %>% 
  distinct(gene_name, abundance)

# compute bins
n_bins = 3
bins = paxdb %>%
  arrange(desc(abundance)) %>%
  mutate(bin = cut(row_number() / n(),
                   breaks = seq(0, n_bins) / n_bins),
         bin = as.integer(bin)) %>%
  split(.$bin) %>%
  map('gene_name') %>%
  set_names(c('highly expressed', 'moderately expressed', 'lowly expressed'))

# iterate through experiments
results = data.frame()
for (expt_idx in seq_len(nrow(expts))) {
  row = expts[expt_idx, ]
  expt_dir = with(row, file.path("~/git/CF-MS-searches/data/chromatograms", 
                                 Accession, Replicate))
  message("[", expt_idx, "/", nrow(expts), "] processing experiment: ", 
          expt_dir)
  
  # iterate through quant. modes
  chrom_files = list.files(expt_dir, pattern = '*rds', full.names = TRUE) %>%
    extract(!grepl("meta", .))
  for (chrom_file in chrom_files) {
    # read chromatogram matrix
    mat = readRDS(chrom_file)
    quant_mode = gsub("\\.rds", "", basename(chrom_file))
    
    # handle infinite or NaN values
    mat[!is.finite(mat)] = NA
    mat[is.nan(mat)] = NA
    
    # also read the corresponding metadata file
    metadata_file = gsub(quant_mode, "metadata", chrom_file)
    meta = readRDS(metadata_file)
    
    # keep one row per gene
    gene_map = meta %>%
      dplyr::select(`Majority protein IDs`, `Gene names`) %>%
      set_colnames(c("protein_group", "gene")) %>%
      mutate(gene = strsplit(gene, ';')) %>%
      unnest(gene) %>%
      drop_na()
    genes = unique(gene_map$gene)
    gene_mat = matrix(NA, nrow = length(genes), ncol = ncol(mat),
                      dimnames = list(genes, colnames(mat)))
    n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
    out_map = data.frame()
    for (gene in genes) {
      protein_groups = gene_map$protein_group[gene_map$gene == gene]
      # pick the best protein for this replicate
      n_fractions0 = n_fractions[protein_groups]
      best = names(which(n_fractions0 == max(n_fractions0))) %>%
        dplyr::first()
      # matrixStats::rowSds(mat[proteins, ], na.rm = T)
      gene_mat[gene, ] = mat[best, ]
      # save the mapping
      out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
    }
    
    # set attribute on the gene matrix
    attr(gene_mat, 'gene_map') = out_map
    
    # remove duplicate rows
    mat = gene_mat
    gene_map = attr(mat, 'gene_map')
    order = rownames(mat) %>% order()
    mat %<>% extract(order, )
    gene_map %<>% extract(order, )
    order = rownames(mat) %in% complex_proteins %>% order(decreasing = TRUE)
    mat %<>% extract(order, )
    gene_map %<>% extract(order, )
    ## now drop duplicated rows
    drop = which(duplicated(mat))
    mat %<>% extract(-drop, )
    gene_map %<>% extract(-drop, )
    # reset attribute
    attr(mat, 'gene_map') = gene_map
    
    # compute % of fractions each protein was detected in
    pct_fractions = rowMeans(is.finite(mat) & mat > 0)
    
    # now, compute mean coverage of highly, lowly, and moderately expressed
    # proteins
    mean_coverage = map(bins, ~ {
      map_dbl(.x, ~ { 
        if (.x %in% rownames(mat)) {
          pct_fractions[.x]
        } else {
          0
        }
      }) %>% setNames(.x)
    })
    
    # append to results
    df = map(mean_coverage, ~ data.frame(protein = names(.), coverage = .)) %>%
      bind_rows(.id = 'bin') %>%
      cbind(row, quant_mode = quant_mode, .)
    results %<>% bind_rows(df)
  }
}

# filter out proteins that were never detected in any experiment
keep = results %>%
  filter(quant_mode == 'iBAQ') %>%
  filter(coverage > 0) %$%
  unique(protein)

# for the remaining proteins, average over bins
avg = results %>%
  filter(quant_mode == 'iBAQ', protein %in% keep) %>%
  group_by_at(vars(-protein, -coverage)) %>%
  summarise(mean = mean(coverage)) %>%
  ungroup()

# plot overall coverage of highly - lowly - moderately expressed proteins
p1 = avg %>%
  mutate(bin = factor(bin, levels = rev(names(bins)))) %>%
  ggplot(aes(x = bin, y = mean, fill = '1', color = '1')) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = 'grey80') +
  scale_color_manual(values = 'grey72') +
  scale_x_discrete('Protein abundance', labels = function(x) gsub("ly.*$", "", x)) +
  scale_y_continuous('Coverage (%)', labels = function(x) x * 100,
                     limits = c(0, 0.6)) +
  boxed_theme() +
  theme(legend.position = 'none')
p1

# plot detection vs. gradient length
grad = read.csv("data/QC/instrument-time.csv") %>%
  filter(Species == 'Homo sapiens') %>%
  # compute the mode over all runs per experiment
  mutate(Time = round(Time)) %>%
  count(Accession, Replicate, Time) %>%
  group_by(Accession, Replicate) %>%
  filter(n == max(n)) %>%
  ungroup() 
## double Scott2017
grad %<>% rbind(
  grad %>% filter(Accession == 'PXD002892') %>% 
    mutate(Replicate = gsub("heavy", "medium", Replicate))
)
## confirm we aren't double-counting any
grad %>% group_by(Accession, Replicate) %>% filter(n() > 1)

# bin time
table(grad$Time) 
grad %<>%
  mutate(time_bin = fct_recode(as.character(Time),
                               ## 55-75 min
                               '0.5-1.5 h' = '55',
                               '0.5-1.5 h' = '55',
                               '0.5-1.5 h' = '60',
                               '0.5-1.5 h' = '65',
                               '0.5-1.5 h' = '75',
                               ## 100-146 min
                               '1.5-2.5 h' = '100',
                               '1.5-2.5 h' = '105',
                               '1.5-2.5 h' = '120',
                               '1.5-2.5 h' = '142',
                               '1.5-2.5 h' = '146',
                               ## 180 min +
                               '2.5+ h' = '180',
                               '2.5+ h' = '240'),
         time_bin = fct_recode(time_bin,
                               'Short (55-75 min, n = 11)' = '0.5-1.5 h',
                               'Intermediate (100-145 min, n = 18)' = '1.5-2.5 h',
                               'Long (180-240 min, n = 18)' = '2.5+ h') %>%
           fct_relevel('Short (55-75 min, n = 11)',
                       'Intermediate (100-145 min, n = 18)',
                       'Long (180-240 min, n = 18)'),
         # time_bin = fct_relevel(time_bin, '0.5-1.5 h', '1.5-2.5 h', '2.5+ h'),
         # also do a very coarse binning
         time_bin2 = fct_recode(as.character(Time),
                                '60' = '55', '60' = '65',
                                '100' = '105',
                                '140' = '142', '140' = '146')) 

# statistical tests
filter(avg, bin == 'lowly expressed') %>% 
  left_join(grad, by = c('Accession', 'Replicate')) %$%
  cor.test(.$Time, .$mean, method = 's')
tests = avg %>%
  left_join(grad, by = c('Accession', 'Replicate')) %>%
  group_by(bin) %>%
  do(tidy(cor.test(.$Time, .$mean, method = 's'))) %>%
  mutate(label = ifelse(p.value < 0.001, '***',
                        ifelse(p.value < 0.01, '**', 
                               ifelse(p.value < 0.05, '*', 'n.s.'))))

# plot
pal = brewer.pal(4, 'Blues')
p2 = avg %>%
  mutate(bin = factor(bin, levels = rev(names(bins)))) %>%
  left_join(grad, by = c('Accession', 'Replicate')) %>%
  ggplot(aes(x = bin, y = mean)) +
  geom_boxplot(aes(fill = time_bin, color = time_bin), 
               outlier.shape = NA, alpha = 0.6) +
  geom_text(data = tests, aes(label = label, y = 0.575), color = 'black',
            size = 2.25) +
  scale_x_discrete('Protein abundance', labels = function(x) gsub("ly.*$", "", x)) +
  scale_y_continuous('Coverage (%)', labels = function(x) x * 100,
                     limits = c(0, 0.6)) +
  scale_color_manual('MS time', values = pal) +
  scale_fill_manual('MS time', values = pal) +
  guides(fill = guide_legend(ncol = 1, nrow = 3, title.pos = 'top'),
         color = guide_legend(ncol = 1, nrow = 3, title.pos = 'top')) +
  boxed_theme()
p2

# combine plots
p = p1 + p2 + 
  plot_layout(nrow = 1, widths = c(0.8, 1.2)) & 
  theme(plot.background = element_blank(), 
        panel.background = element_blank())
p
ggsave("fig/QC/lowly-expressed.pdf", p, 
       width = 8, height = 6.5, units = "cm", useDingbats = FALSE)

# are protein complexes depleted in lowest bin?
corum = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt")
complex_proteins = with(corum, unique(gene_name))
pal = c(Gpal[4], jdb_palette("wolfgang_extra")[4])
p3 = paxdb %>%
  mutate(complex = gene_name %in% complex_proteins) %>%
  ggplot(aes(x = complex, y = abundance + 1, fill = complex, color = complex)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  scale_x_discrete(labels = c('Non-complex\nprotein', 'Complex\nprotein')) +
  scale_y_log10('Protein abundance (ppm)') +
  scale_fill_manual(values = pal, guide = FALSE) +
  scale_color_manual(values = darken(pal, 1.15), guide = FALSE) +
  annotation_logticks(sides = 'l', short = unit(0.04, 'cm'),
                      mid = unit(0.08, 'cm'), long = unit(0.12, 'cm'), 
                      size = 0.2) +
  # coord_cartesian(ylim = c(1, 1e5)) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p3
ggsave("fig/QC/complex-protein-abundance.pdf", p3, 
       width = 3, height = 4.5, units = "cm", useDingbats = FALSE)
