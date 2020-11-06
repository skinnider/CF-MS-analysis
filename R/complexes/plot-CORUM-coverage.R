# Calculate the number of CORUM complex proteins quantified in each 
# human and mouse CF-MS dataset.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)
source("R/theme.R")

# read CORUM 
corum_hs = read.delim(
  "~/git/network-validation/data/complex/CORUM/complexes_human.txt")
corum_mm = read.delim(
  "~/git/network-validation/data/complex/CORUM/complexes_mouse.txt")

# filter complexes to human genes
fa = read.fasta(
  "~/git/CF-MS-searches/data/fasta/filtered/UP000005640-H.sapiens.fasta.gz",
  as.string = T, seqtype = 'AA')
human_genes = map_chr(fa, getAnnot) %>%
  gsub("^.*GN=", "", .) %>%
  gsub(" .*$", "", .) %>%
  unique() %>%
  # remove proteins without gene names
  extract(!startsWith(., ">"))
corum_hs %<>% filter(gene_name %in% human_genes)

# create list of complexes
complexes = list(`Homo sapiens` = corum_hs, `Mus musculus` = corum_mm)

# read ortholog mapping
ortho = read.delim(
  "~/git/network-validation/data/ortholog/InParanoid-H.sapiens-M.musculus.txt.gz")

# read all human and mouse datasets
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species %in% c("Homo sapiens", "Mus musculus"))  

# calculate CORUM complex coverage in each dataset
coverage = data.frame()
for (idx in seq_len(nrow(expts))) {
  row = expts[idx, ]
  expt_dir = with(row, file.path("~/git/CF-MS-searches/data/chromatograms", 
                                 Accession, Replicate))
  message("[", idx, "/", nrow(expts), "] processing experiment: ", expt_dir)
  
  chrom_files = list.files(expt_dir, pattern = '*.rds', full.names = T) %>%
    extract(basename(.) != 'metadata.rds') %>%
    setNames(gsub("\\..*$", "", basename(.)))
  
  # read metadata (mapping to genes)
  meta_file = file.path(expt_dir, 'metadata.rds')
  if (!file.exists(meta_file)) {
    warning("file does not exist: ", meta_file)
    next
  }
  meta = readRDS(meta_file)
  
  # process each quantitation strategy independently
  for (quant_mode in names(chrom_files)) {
    chrom_file = chrom_files[quant_mode]
    chrom = readRDS(chrom_file)
    
    complex_genes = complexes[['Homo sapiens']]$gene_name %>%
      trimws() %>%
      unique()
    if (row$Species == "Homo sapiens") {
      # map protein complex genes to protein groups
      gene_meta = meta %>%
        mutate(gene = strsplit(`Gene names`, ';')) %>%
        unnest(gene) %>%
        # filter to protein complex genes
        filter(gene %in% complex_genes)
      protein_map = with(gene_meta, setNames(`Majority protein IDs`, gene))
    } else if (row$Species == "Mus musculus") {
      # in mouse, first map to human orthologs
      ortho_map = with(ortho, setNames(source_gene, target_gene))
      # subset complex genes to those with an ortholog in human
      complex_genes %<>% intersect(ortho_map)
      
      # map protein complex genes to protein groups
      gene_meta = meta %>%
        mutate(gene = strsplit(`Gene names`, ';')) %>%
        unnest(gene) %>%
        # map to ortholog
        mutate(ortho = ortho_map[gene]) %>%
        # filter to protein complex genes
        filter(ortho %in% complex_genes)
      protein_map = with(gene_meta, setNames(`Majority protein IDs`, ortho))
    }
    protein_groups = protein_map[complex_genes]
    
    # count number of fractions each complex gene was found in
    n_fractions = integer(length(protein_groups))
    complex_mat = chrom[na.omit(protein_groups), ]
    n_fractions[!is.na(protein_groups)] = rowSums(!is.na(complex_mat) & 
                                                    complex_mat > 0)
    
    # append 
    res = data.frame(quant_mode = quant_mode,
                     gene = complex_genes,
                     n_fractions = n_fractions) %>%
      cbind(row, .)
    coverage %<>% bind_rows(res)
  }
}

# save
saveRDS(coverage, "data/analysis/complexes/coverage.rds")

# visualize the shifting CDF:
# randomly pick one experiment, two experiments, three experiments... 
# (repeat 10x)
sampled = data.frame()
n_samples = 10
for (n_experiments in c(seq_len(4), seq(5, 60, 5), 64)) {
  message("drawing random samples of ", n_experiments, " experiments ...")
  for (sample_idx in seq_len(n_samples)) {
    message(".. random sample ", sample_idx, " of ", n_samples, " ...")
    set.seed(sample_idx)
    
    # draw a sample of n experiments
    query_expts = expts %>%
      sample_frac(1) %>%
      head(n_experiments)
    
    # calculate number of fractions each complex protein detected in,
    # summed across experiments
    df = coverage %>% 
      filter(quant_mode == 'iBAQ') %>%
      inner_join(query_expts) %>% 
      # sum all fractions
      group_by(gene) %>%
      summarise(sum_fractions = sum(n_fractions),
                expts_quantified = sum(n_fractions > 0)) %>%
      ungroup() %>%
      # flag sample index and n expts
      mutate(sample_idx = sample_idx, n_expts = n_experiments)
    
    # append to results
    sampled %<>% bind_rows(df)
  }
}

# now, calculate CDF
cdf = map_dfr(seq_len(50), ~ {
  min_fractions = .
  sampled %>%
    group_by(n_expts, sample_idx) %>%
    summarise(y = mean(sum_fractions >= min_fractions)) %>%
    ungroup() %>%
    # estimate s.d. over samples
    group_by(n_expts) %>%
    summarise(mean = mean(y), sd = sd(y)) %>%
    ungroup() %>%
    # flag minimum number of fractions
    mutate(min_fractions = min_fractions)
})

# save
saveRDS(cdf, "data/analysis/complexes/complex-cdf.rds")

# plot 
pal = BuenColors::jdb_palette("solar_extra") %>% colorRampPalette()
p = cdf %>%
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
                              override.aes = list(size = 0.7))) +
  boxed_theme() +
  theme(legend.key.size = unit(0.5, 'lines'),
        legend.position = 'right')
p
ggsave("fig/analysis/complexes/complex-cdf-CORUM.pdf", p, 
       width = 5.5, height = 4.2, units = 'cm', useDingbats = F)

# how many complex proteins were quantified in at least one experiment?
coverage %>%
  group_by(gene) %>%
  summarise(detected = sum(n_fractions) >= 1) %>%
  ungroup() %>%
  summarise(mean = mean(detected))
