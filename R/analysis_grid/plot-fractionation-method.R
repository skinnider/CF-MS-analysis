# Plot complex and GO AUCs as a function of the method used to separate
# co-eluting proteins.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read complexes
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # restore 'missing=NA'
  replace_na(list(missing = 'NA')) %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))

# read GO
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))
# for each combination, filter by breadth, and get median AUROC across GO terms
dat2 %<>%
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# read metadata
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")
fractionation = expts %>%
  dplyr::rename(accession = Accession, experiment = Replicate,
                fractionation = Fractionation) %>%
  dplyr::select(accession, experiment, fractionation)
dat1 %<>% left_join(fractionation, by = c('accession', 'experiment'))
dat2 %<>% left_join(fractionation, by = c('accession', 'experiment'))

# plot complex AUCs
dat1a = dat1 %>%
  # ignore cross-linked datasets
  filter(!grepl("^XL-", fractionation)) 
medians1a = dat1a %>%
  group_by(fractionation) %>%
  summarise(median = median(auroc, na.rm = TRUE),
            label = formatC(median, format = 'f', digits = 3))
arrange(medians1a, median)
pal = pals::stepped3()[c(1, 5, 10, 15, 19)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
p1 = dat1a %>%
  ggplot(aes(y = reorder(fractionation, auroc, median), x = auroc, 
             fill = fractionation, color = fractionation)) +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians1a, aes(x = 0.46, y = fractionation, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC', breaks = seq(0.5, 0.7, 0.1)) +
  coord_cartesian(xlim = c(0.46, 0.7)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1
ggsave("fig/analysis/analysis_grid/fractionation-AUC-CORUM.pdf", p1,
       width = 4.5, height = 3, units = "cm", useDingbats = FALSE)

# plot all metrics as a heatmap
summary1a = dat1a %>%
  group_by(fractionation, metric) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()
p2 = summary1a %>%
  mutate(metric = clean_metric(metric)) %>%
  ggplot(aes(y = reorder(fractionation, auroc, stats::median), 
             x = reorder(metric, -auroc, stats::median), 
             fill = auroc, color = auroc)) +
  geom_tile(color = 'white') +
  # geom_text(data = medians1a, aes(x = 0.475, y = fractionation, label = label),
  #           size = 2, color = 'grey20', hjust = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC',
                         breaks = c(0.45, 0.7), 
                         limits = c(0.45, 0.7)
                         ) +
  boxed_theme() +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'))
p2
ggsave("fig/analysis/analysis_grid/fractionation-AUC-CORUM-heatmap.pdf", p2,
       width = 9.5, height = 3.75, units = "cm", useDingbats = FALSE)

# plot GO AUCs
dat2a = dat2 %>%
  # drop cross-linked datasets
  filter(!grepl("^XL-", fractionation)) %>%
  # drop sucrose
  filter(!grepl("sucrose", fractionation, ignore.case = TRUE))
medians2a = dat2a %>%
  group_by(fractionation) %>%
  summarise(median = median(auroc, na.rm = TRUE),
            label = formatC(median, format = 'f', digits = 3))
arrange(medians2a, median)
pal = pals::stepped3()[c(1, 5, 10, 15)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF"))
p3 = dat2a %>%
  ggplot(aes(y = reorder(fractionation, auroc, median), x = auroc, 
             fill = fractionation, color = fractionation)) +
  geom_boxploth(width = 0.6, alpha = 0.6, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2a, aes(x = 0.48, y = fractionation, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC', breaks = seq(0.5, 0.6, 0.02)) +
  coord_cartesian(xlim = c(0.48, 0.555)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p3
ggsave("fig/analysis/analysis_grid/fractionation-AUC-GO.pdf", p3,
       width = 4.5, height = 3, units = "cm", useDingbats = FALSE)

# combine p1/p2
p4 = p1 + p3 + plot_layout(ncol = 1)#, heights = c(4, 5))
ggsave("fig/analysis/analysis_grid/fractionation-AUC.pdf", p4,
       width = 4.75, height = 5.3, units = "cm", useDingbats = FALSE)

# all metrics, as a heatmap
summary2a = dat2a %>%
  group_by(fractionation, metric) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()
p5 = summary2a %>%
  mutate(metric = clean_metric(metric)) %>%
  ggplot(aes(y = reorder(fractionation, auroc, stats::median), 
             x = reorder(metric, -auroc, stats::median), 
             fill = auroc, color = auroc)) +
  geom_tile(color = 'white') +
  # geom_text(data = medians1a, aes(x = 0.475, y = fractionation, label = label),
  #           size = 2, color = 'grey20', hjust = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC',
                         breaks = c(0.5, 0.55), 
                         # limits = c(0.45, 0.7)
  ) +
  boxed_theme() +
  guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.25, 'lines'),
        legend.key.height = unit(0.25, 'lines'))
p5
ggsave("fig/analysis/analysis_grid/fractionation-AUC-GO-heatmap.pdf", p5,
       width = 9.5, height = 3.75, units = "cm", useDingbats = FALSE)


# plot CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds") %>%
  filter(quant_mode == 'iBAQ')
mean_cdf = cdf %>%
  left_join(fractionation) %>%
  group_by(fractionation, n_fractions) %>%
  summarise_at(vars(n_proteins, coverage), mean) %>%
  ungroup()
pal = pals::stepped3()[c(1, 5, 10, 15)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF"))
p6 = mean_cdf %>%
  ## drop cross-linking datasets
  filter(!grepl("^XL-", fractionation)) %>%
  # drop sucrose
  filter(fractionation != 'Sucrose') %>%
  filter(n_fractions <= 50) %>%
  ggplot(aes(x = n_fractions, y = n_proteins, color = fractionation)) +
  geom_line(size = 0.5) +
  # geom_point(size = 0.8) + 
  scale_y_continuous("# of proteins", breaks = seq(0, 5000, 1000)) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.key.size = unit(0.5, 'lines'))
p6
ggsave("fig/analysis/analysis_grid/fractionation-CDF.pdf", p3, 
       width = 7, height = 4.5, units = "cm")

# combine with protein complex
p6b = p6 + 
  scale_y_continuous(expression(Proteins~(10^3)),
                     labels = function(x) x / 1e3,
                     breaks = seq(0, 5000, 1000)) +
  theme(legend.position = c(0.95, 1.1),
        legend.justification = c(1, 1))
p = p6b + p1 + plot_layout(heights = c(1.7, 1))
# p
ggsave("fig/analysis/analysis_grid/fractionation-CDF-combined.pdf", p, 
       width = 5, height = 6, units = "cm")
