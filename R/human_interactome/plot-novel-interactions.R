setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read the CF-MS meta-interactome 
source("R/functions/detect_system.R")
net = readRDS(file.path(base_dir, "human_interactome", 
                        "network-classifier=RF-n_features=1-feature_select=best_first-n_datasets=46-sample_idx=1.rds"))
# filter to 50% precision
net0 = head(net, 47575)

# read all known interactions 
db_files = list.files("~/git/network-validation/data/database", 
                      pattern = '*-human\\.txt\\.gz$', full.names = TRUE)
dbs = map(db_files, read.delim)
db = bind_rows(dbs)
## alphabetize
sorted = t(apply(db[, 1:2], 1, sort))
db$gene_A = sorted[, 1]
db$gene_B = sorted[, 2]
db %<>% distinct(gene_A, gene_B, pmid)
# count PMIDs per interaction
known = db %>%
  dplyr::count(gene_A, gene_B) %>%
  dplyr::rename(protein_A = gene_A, protein_B = gene_B)

# flag number of PMIDs
pubs = net0 %>%
  left_join(known) %>%
  replace_na(list(n = 0))
 
# plot 
pal = brewer.pal(6, 'Blues')
p = pubs %>%
  mutate(fill = ifelse(n >= 5, 5, n) %>% 
           as.character() %>%
           fct_recode('0' = '0',
                      '5+' = '5')) %>%
  ggplot(aes(x = '1', fill = fill)) +
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0, yend = 1),
               color = 'grey50') +
  geom_bar(color = 'grey15', size = 0.15, position = 'fill') +
  scale_fill_manual('Publications', values = pal) +
  scale_y_reverse('% of interactions', # expand = c(0, 0),
                  breaks = seq(0, 1, 0.25), 
                  labels = function(x) 100 * (1 - x)) +
  coord_flip() +
  clean_theme() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        legend.position = 'top')
p
ggsave("fig/analysis/human_interactome/novel-interactions.pdf",
       p, width = 4.25, height = 2.5, units = 'cm', useDingbats = FALSE)
