# Run a GO enrichment analysis of complex proteins that were never quantified
# in any human or mouse CF-MS dataset.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)
library(GOstats)
source("R/theme.R")

# read complex protein coverage
coverage = readRDS("data/analysis/complexes/coverage.rds")

# group by proteins
proteins = coverage %>%
  group_by(gene) %>%
  summarise(detected = any(n_fractions > 0))

# how many proteins (and what %) were never detected?
proteins %>%
  summarise(n = sum(!detected), pct = mean(!detected))

# to do GOstats analysis, we need to map to Entrez IDs
map = read_tsv("~/git/network-validation/data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
id = filter(map, db == 'GeneID') %>% dplyr::select(-db)
gn2id = left_join(gn, id, by = 'uniprot') %>%
  dplyr::select(starts_with('id')) %>%
  set_colnames(c('gene', 'id'))
entrez = proteins %>%
  left_join(gn2id, by = 'gene') %>%
  drop_na(id) %>%
  distinct(id, detected)

# do GO analysis
enrs = list()
universe = entrez %>% pull(id)

# test overrepresentation among detected, never detected
for (test in c("detected", "never detected")) {
  if (test == "detected") {
    hits = entrez %>% filter(!detected) %>% pull(id)
  } else {
    hits = entrez %>% filter(detected) %>% pull(id)
  }
  
  # run GOstats enrichment analysis
  for (ontology in c("BP", "MF", "CC")) {
    message(".. analyzing enrichment of ", ontology, " terms ...")
    for (cond in c(T, F)) {
      message("... performing test with conditional parameter = ", cond)
      params = new("GOHyperGParams", 
                   geneIds = hits, universeGeneIds = universe,
                   annotation = "org.Hs.eg.db", ontology = ontology,
                   pvalueCutoff = 1, conditional = cond, 
                   testDirection = 'over')
      hgOver = hyperGTest(params)
      enr = summary(hgOver) %>%
        mutate(root = ontology, 
               conditional = cond,
               test = test)
      enrs[[length(enrs) + 1]] = enr
    }
  }
}

# save
saveRDS(enrs, "data/analysis/complexes/GO-never-detected.rds")

# plot
enr = enrs %>%
  bind_rows() %>%
  filter(conditional == TRUE)
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
              'integral component of plasma membrane',
              'G protein-coupled receptor activity',
              'ion channel activity',
              'neurotransmitter receptor activity',
              'cytokine activity',
              'plasma membrane part',
              'regulation of cell proliferation',
              'cell communication')
pal = get_pal("early meadow")[c(1, 3)]
p1a = enr %>%
  filter(test == 'never detected') %>%
  filter(Term %in% terms_detected) %>%
  ggplot(aes(x = reorder(Term, -Pvalue), y = -log10(Pvalue))) + 
  facet_grid(~ "Detected") +
  geom_col(width = 0.7, alpha = 0.7, size = 0.3, aes(color = '1', fill = '1')) +
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
p1a  
p1b = enr %>%
  filter(test == 'detected') %>%
  filter(Term %in% terms_not) %>%
  ggplot(aes(x = reorder(Term, -Pvalue), y = -log10(Pvalue))) + 
  facet_grid(~ "Never detected") +
  geom_col(width = 0.7, alpha = 0.7, size = 0.3, aes(color = '1', fill = '1')) +
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
p1b

# combine
p1 = p1b | p1a
p1
ggsave("fig/analysis/complexes/GO-never-detected.pdf", p1,
       width = 14, height = 4.75, units = "cm", useDingbats = F)
