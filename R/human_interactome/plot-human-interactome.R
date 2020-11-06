setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)

# read PR curves
curves = readRDS("data/analysis/human_interactome/PR.rds") %>%
  type_convert()

# plot saturation in the RF classifier
gradient = BuenColors::jdb_palette("brewer_celsius") %>% colorRampPalette() %>% 
  do.call(list(100))
p1 = curves %>%
  filter(precision >= 0.5, classifier == 'RF', feature_select == 'best_first',
         n_features == 1, sample_idx == 1) %>%
  ggplot(aes(x = idx, y = precision, color = n_datasets)) +
  geom_path(aes(group = n_datasets)) + 
  scale_color_gradientn(colours = gradient, name = '# of datasets',
                        breaks = c(1, 46)) +
  scale_y_continuous('Precision', breaks = seq(0.5, 1, 0.1)) +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3,
                     limits = c(0, 50e3)) +
  guides(color = guide_colorbar(ticks = FALSE, frame.colour = 'black',
                                title.position = 'top')) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        # legend.title.align = 1,
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'))
p1
ggsave("fig/analysis/human_interactome/PR-curve-saturation.pdf", p1, 
       width = 5, height = 4.25, units = 'cm', useDingbats = FALSE)

# compare classifiers with all 46 datasets
pal = brewer.pal(10, 'Set3')[-2]
p2 = curves %>%
  filter(precision >= 0.5, n_datasets == 46, classifier != 'RF2',
         feature_select != 'random', n_features == 1) %>%
  ggplot(aes(x = idx, y = precision, color = classifier)) +
  geom_path(aes(group = classifier)) + 
  scale_color_brewer(palette = 'Set2', name = 'Classifier') +
  scale_y_continuous('Precision', breaks = seq(0.5, 1, 0.1)) +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3,
                     limits = c(0, 50e3)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = c(0.94, 0.96),
        legend.justification = c(1, 1),
        legend.direction = 'vertical')
p2
ggsave("fig/analysis/human_interactome/PR-curve-classifier.pdf", p2, 
       width = 4.6, height = 4.6, units = 'cm', useDingbats = FALSE)

# plot RF against human high-throughput interactomes
final = filter(curves, n_datasets == 46, feature_select == 'best_first',
               classifier == 'RF', n_features == 1, sample_idx == 1)
hts_files = list.files("~/git/network-validation/data/high_throughput",
                       full.names = TRUE, pattern = '\\.gz$') %>%
  # ignore external data
  extract(!grepl("Havugimana|Wan", .))
hts_ppis = map(hts_files, read.delim) %>%
  setNames(gsub("\\..*$", "", basename(hts_files)))
# calculate precision and recall
complexes = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')
adj = adjacency_matrix_from_list(complexes)
precision = map_dbl(hts_ppis, ~ {
  make_labels(adj, .) %>% 
    mean(na.rm = TRUE)
})
recall = map_int(hts_ppis, nrow)
n_proteins = map_int(hts_ppis, ~ with(., n_distinct(c(gene_A, gene_B))))
hts = data.frame(interactome = names(hts_ppis),
                 precision = precision,
                 recall = recall,
                 n_proteins = n_proteins) %>%
  mutate(interactome = fct_recode(interactome,
                                  'BioPlex' = 'Huttlin2015',
                                  'BioPlex 2' = 'Huttlin2017',
                                  'QUBIC' = 'Hein2015',
                                  'HI-II-14' = 'Rolland2014',
                                  'HuRI' = 'Luck2019'))
pal = pals::stepped3()[seq(1, 19, 4)]
p3 = final %>%
  filter(idx <= 6e4) %>%
  ggplot(aes(x = idx, y = precision)) +
  geom_path(linetype = 'dotted', size = 0.4) +
  # geom_path(linetype = 'dotted', size = 0.3) +
  geom_point(data = hts, aes(x = recall, color = interactome), size = 0.9) +
  geom_label_repel(data = hts, 
                   aes(x = recall, color = interactome, label = interactome),
                   size = 2.25, label.padding = 0.25, label.size = NA,
                   fill = NA, min.segment.length = 0, segment.size = 0.4) +
  scale_y_continuous('Precision') +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = 'none')
p3
ggsave("fig/analysis/human_interactome/PR-curve-HTS.pdf", p3, 
       width = 5.5, height = 4.25, units = 'cm', useDingbats = FALSE)

# plot final PR curve against every individual-dataset PR curve
col = pals::stepped3()[2]
p4 = curves %>%
  filter(precision >= 0.5, classifier == 'RF', feature_select == 'best_first',
         n_features == 1, n_datasets %in% c(1, 46)) %>%
  mutate(color = fct_recode(as.character(n_datasets),
                            'single datasets' = '1',
                            'meta-analysis' = '46')) %>%
  ggplot(aes(x = idx, y = precision, color = color)) +
  geom_path(aes(group = sample_idx)) + 
  scale_y_continuous('Precision', breaks = seq(0.5, 1, 0.1)) +
  scale_x_continuous('Interactions (thousands)', labels = function(x) x / 1e3,
                     limits = c(0, 50e3)) +
  scale_color_manual('', values = c('grey80', col)) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.position = c(0.97, 1.02),
        legend.justification = c(1, 1),
        legend.direction = 'vertical',
        legend.title = element_blank(),
        # legend.title.align = 1,
        legend.key.size = unit(0.55, 'lines'))
p4
ggsave("fig/analysis/human_interactome/PR-curve-individual-datasets.pdf", p4,
       width = 4.6, height = 4.6, units = 'cm', useDingbats = FALSE)
