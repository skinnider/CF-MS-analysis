# Plot the features that have been previously used, in published work, to 
# identify co-eluting proteins. 
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read review
dat = read.csv("data/analysis/analysis_grid/features.csv")

# convert to matrix
obs = dat %>%
  unite(AuthorYear, Author, Year, sep = "") %>%
  dplyr::select(-Accession, -X) %>%
  # convert some of the measures to 'other'
  mutate(Measure.of.association = fct_recode(Measure.of.association,
                                             'Other' = 'Hypergeometric',
                                             'Other' = 'Normalized cross-correlation')) %>%
  # drop 'noised' Pearson correlation - not fundamentally different
  filter(Measure.of.association != 'Noised pearson correlation')
rownames = unique(obs$AuthorYear)
colnames = unique(na.omit(obs$Measure.of.association))
mat = matrix(0, nrow = length(rownames), ncol = length(colnames),
             dimnames = list(rownames, colnames))
mat[as.matrix(obs)] = 1

# cluster
clust1 = hclust(dist(mat))
lvl1 = with(clust1, labels[order])
clust2 = hclust(dist(t(mat)))
lvl2 = with(clust2, labels[order])

# back to tidy
feat = reshape2::melt(mat, varnames = c('study', 'metric'),
                      value.name = 'fill', as.is = T) %>%
  mutate(study = factor(study, levels = lvl1),
         metric = factor(metric, levels = lvl2))

# plot
pal = brewer.pal(8, 'Dark2')
col = pal[8]
p1 = feat %>%
  mutate(facet = metric == 'External') %>%
  ggplot(aes(x = metric, y = study, fill = factor(fill))) + 
  facet_grid(~ facet, scales = 'free', space = 'free') + 
  geom_tile(color = 'white') +
  scale_x_discrete('Feature', expand = c(0, 0)) +
  scale_y_discrete('Study', expand = c(0, 0)) +
  scale_fill_manual(values = c('grey94', col), guide = F) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        strip.text = element_blank())
p1

p2 = feat %>%
  dplyr::count(metric, fill) %>%
  filter(fill != 0) %>%
  mutate(facet = metric == 'External') %>%
  ggplot(aes(x = metric, y = n)) +
  facet_grid(~ facet, scales = 'free', space = 'free') + 
  geom_col(aes(fill = '1'), color = 'grey94', width = 0.98, size = 0.15) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous('# of studies', expand = c(0, 0), limits = c(0, 16),
                     breaks = seq(0, 16, 4)) +
  scale_fill_manual(values = col, guide = F) +
  scale_color_manual(values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(rep(0, 4)),
        strip.text = element_blank(),
        panel.background = element_rect(color = NA, fill = 'grey94'))
p2

# combine
p = p2 + p1 + plot_layout(heights = c(0.2, 1), ncol = 1)
p
ggsave("fig/analysis/analysis_grid/previously-used-features.pdf", p,
       width = 7, height = 10.5, units = 'cm', useDingbats = F)
