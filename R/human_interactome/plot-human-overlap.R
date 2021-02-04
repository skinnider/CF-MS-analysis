# Plot the overlap between the 'consensus' human interactome and:
# - Havugimana2012
# - hu.MAP
# - hu.MAP 2
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(VennDiagram)
library(CFTK)
source("R/theme.R")

# read the CF-MS meta-interactome 
source("R/functions/detect_system.R")
net = readRDS(file.path(base_dir, "human_interactome", 
                        "network-classifier=RF-n_features=1-feature_select=best_first-n_datasets=46-sample_idx=1.rds"))
# filter to 50% precision
net0 = head(net, 47575)

# read the other networks
## Havugimana2012
havugimana2012 = read.delim("~/git/network-validation/data/high_throughput/Havugimana2012.txt.gz")

## hu.MAP 1
humap1 = readLines("~/git/network-validation/data/hu.MAP/1.0/genename_clusters.txt") %>%
  strsplit("\t") %>%
  to_pairwise_df() %>%
  set_colnames(c("gene_A", "gene_B"))
  
## hu.MAP 2
humap2 = read.csv("~/git/network-validation/data/hu.MAP/2.0/humap2_complexes_20200809.txt") %>%
  pull(genenames) %>%
  strsplit(' ') %>%
  to_pairwise_df() %>%
  set_colnames(c("gene_A", "gene_B"))

# create list of networks
comparison_list = list('Havugimana et al., 2012' = havugimana2012,
                       'hu.MAP 1.0' = humap1,
                       'hu.MAP 2.0' = humap2)
for (comparison_network in names(comparison_list)) {
  target = comparison_list[[comparison_network]]
  
  # alphabetize target network interactors
  alphabetized = t(apply(target[, 1:2], 1, sort))
  target$gene_A = alphabetized[, 1]
  target$gene_B = alphabetized[, 2]
  target %<>% distinct(gene_A, gene_B)
  
  # get overlaps
  overlap1 = target %>%
    dplyr::rename(protein_A = gene_A, protein_B = gene_B) %>%
    left_join(net, by = c('protein_A', 'protein_B')) %>%
    mutate(flag = !is.na(precision))
  overlap2 = net0 %>%
    dplyr::rename(gene_A = protein_A, gene_B = protein_B) %>%
    left_join(target %>% mutate(flag = TRUE),
              by = c('gene_A', 'gene_B')) %>%
    replace_na(list(flag = FALSE))
  
  ## Venn diagram
  venn_file = paste0("fig/analysis/human_interactome/overlap/venn-",
                     fct_recode(comparison_network,
                                'Havugimana2012' = 'Havugimana et al., 2012'),
                     '.pdf')
  pdf(venn_file, width = 1.15, height = 1.15)
  draw.pairwise.venn(area1 = nrow(net),
                     area2 = nrow(target),
                     cross.area = sum(overlap1$flag),
                     category = c('Consensus CF-MS interactome',
                                  comparison_network),
                     scaled = TRUE,
                     lwd = c(0.5, 0.5),
                     col = c('grey20', 'grey20'),
                     fill = c('grey80', brewer.pal(8, 'Set2')[3]),
                     fontfamily = rep('sans', 3),
                     cat.fontfamily = rep('sans', 2),
                     cex = rep(0.514, 3),
                     cat.cex = rep(0.514, 2))
  dev.off()
  
  ## % of CF-MS interactions in the target network
  pal = brewer.pal(n = 6, name = 'Blues')[c(1, 5)]
  p2 = overlap2 %>%
    mutate(flag = factor(as.character(flag), levels = c('FALSE', 'TRUE'))) %>%
    ggplot(aes(x = '1', fill = flag)) +
    geom_bar(color = 'grey15', size = 0.15) +
    scale_y_continuous('CF-MS interactions (thousands)', expand = c(0, 0),
                       labels = function(x) x / 1e3) +
    scale_fill_manual(values = pal, name = '',
                      breaks = c('TRUE'),
                      labels = paste0('Overlap with ', comparison_network,
                                     ' n = ', format(sum(overlap2$flag), big.mark = ','), 
                                     ' (', round(100 * mean(overlap2$flag), digits = 1), '%)')) +
    coord_flip() +
    clean_theme() +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_line(size = 0.2),
          axis.line.y = element_blank(),
          axis.title.x = element_text(size = 6))
  p2
  
  ## % of target network interactions in CF-MS
  pal = brewer.pal(n = 6, name = 'Blues')[c(1, 5)]
  p3 = overlap1 %>%
    mutate(flag = factor(as.character(flag), levels = c('FALSE', 'TRUE'))) %>%
    ggplot(aes(x = '1', fill = flag)) +
    geom_bar(color = 'grey15', size = 0.15) +
    scale_y_continuous(paste(comparison_network, 'interactions (thousands)'),
                       expand = c(0, 0), labels = function(x) x / 1e3) +
    scale_fill_manual(values = pal, name = '',
                      breaks = c('TRUE'),
                      labels = paste0('Overlap with consensus CF-MS interactome',
                                      ' n = ', format(sum(overlap1$flag), big.mark = ','), 
                                      ' (', round(100 * mean(overlap1$flag), digits = 1), '%)')) +
    coord_flip() +
    clean_theme() +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_line(size = 0.2),
          axis.line.y = element_blank(),
          axis.title.x = element_text(size = 6))
  p3
  
  # combine 
  p = p2 / p3
  p
  bar_file = paste0("fig/analysis/human_interactome/overlap/bars-", 
                     fct_recode(comparison_network, 
                                'Havugimana2012' = 'Havugimana et al., 2012'),
                     '.pdf')
  ggsave(bar_file, p, width = 4, height = 4.9,
         units = "cm", useDingbats = FALSE)
}
