# Plot the third column of Figure 2:
#' h. min fractions
#' i. min pairs 
#' j. label-free quantitation
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

###############################################################################-
###### h. Min fractions ####
###############################################################################-

datH = readRDS("data/analysis/min_fractions/complexes.rds")

# join with CDF
cdf = readRDS("data/QC/protein-groups-CDF.rds")
## calculate maximum % of fractions per experiment
max_fractions = cdf %>%
  group_by(accession, experiment) %>%
  summarise(max_fractions = max(n_fractions)) %>%
  ungroup()
## calculate maximum # of proteins per experiment
cdf1 = filter(cdf, n_fractions == 1) %>%
  dplyr::select(-n_fractions, -coverage) %>%
  dplyr::rename(n_1 = n_proteins)

# calculate mean over # fractions
meanH = datH %>%
  filter(metric == 'pearson') %>%
  # fix quant mode and merge in CDF
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  left_join(cdf, by = c('accession', 'experiment', 'quant_mode',
                        'n_fractions')) %>%
  # also merge in the total number of fractions
  left_join(max_fractions, by = c('accession', 'experiment')) %>%
  # drop experiments without at least 20 fractions
  filter(max_fractions >= 20) %>%
  # merge in number of proteins found in one fraction
  left_join(cdf1, by = c("accession", "experiment", "search", "quant_mode")) %>%
  # calculate mean AUC and # of proteins
  group_by(quant_mode, n_fractions) %>%
  summarise(mean = mean(auroc, na.rm = T),
            sd = sd(auroc, na.rm = T),
            median = median(auroc, na.rm = T),
            q1 = quantile(auroc, na.rm = T, probs = 0.25),
            q3 = median(auroc, na.rm = T, probs = 0.75),
            pct_proteins = mean(n_proteins / n_1),
            n_proteins = mean(n_proteins),
            n = n()) %>%
  ungroup()

# plot mean and s.d.: AUC (iBAQ only)
pal = brewer.pal(8, 'Set2') %>% extract(c(3, 2, 4, 1, 7, 8))
labHa = meanH %>%
  filter(quant_mode == 'iBAQ', n_fractions == 4) %>%
  mutate(label = paste0(n_fractions, ' fractions\nAUC = ', 
                        format(mean, format = 'f', digits = 2)))
pHa = meanH %>%
  filter(quant_mode == 'iBAQ') %>%
  ggplot(aes(x = n_fractions, y = mean, color = quant_mode)) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  # geom_label_repel(data = labHa, aes(label = label), size = 2, fill = NA,
  #                  color = 'black', label.padding = 0.4, hjust = 0,
  #                  nudge_x = 0.5, nudge_y = -0.03, label.size = NA, 
  #                  segment.size = 0.25) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions') +
  scale_y_continuous('AUC', # limits = c(0.6575, 0.6815),
                     breaks = seq(0.66, 0.68, 0.005)) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none') + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(rep(0, 4)))
pHa

# plot mean and s.d.: # of proteins (as a %, iBAQ only)
labHb = meanH %>%
  filter(quant_mode == 'iBAQ', n_fractions == 4) %>%
  mutate(label = paste0(n_fractions, ' fractions\n', 
                        format(100 * pct_proteins, format = 'f', digits = 3),
                        '% of proteins'))
pHb = meanH %>%
  filter(quant_mode == 'iBAQ') %>%
  ggplot(aes(x = n_fractions, y = pct_proteins, color = quant_mode)) +
  geom_rect(fill = 'grey93', color = NA, alpha = 0.1,
            aes(xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = +Inf)) +
  geom_line(size = 0.5) +
  # geom_label_repel(data = labHb, aes(label = label), size = 2, 
  #                  color = 'black', label.padding = 0.4, hjust = 0, fill = NA,
  #                  nudge_x = 8, nudge_y = 0.1, label.size = NA, 
  #                  segment.size = 0.25) +
  geom_point(size = 0.8) + 
  scale_x_continuous('Minimum # of fractions', breaks = c(1, seq(5, 20, 5))) +
  scale_y_continuous('% of proteins', labels = function(x) x * 100,
                     limits = c(0.4, 1)
  ) +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.position = 'none')
pHb

# combine
pH = pHa + pHb + plot_layout(ncol = 1, heights = c(1.15, 0.85))
pH

###############################################################################-
###### i. Min pairs ####
###############################################################################-

datIa = readRDS("data/analysis/min_pairs/complexes.rds")
# filter to those where min_fractions is not NA
datIa %<>% filter(!is.na(min_fractions))
# calculate mean AUC per min_fractions/min_pairs
meansIa = datIa %>%
  filter(metric == 'pearson') %>%
  filter(quant_mode == 'iBAQ') %>%
  group_by(min_pairs, min_fractions) %>%
  summarise(n = n(), auroc = mean(auroc, na.rm = T)) %>%
  ungroup()

# add in mean AUC from min_fractions experiment, where min_pairs == 0
datIb = readRDS("data/analysis/min_fractions/complexes.rds")
# calculate mean over # fractions
meansIb = datIb %>%
  filter(metric == 'pearson') %>%
  filter(quant_mode == 'iBAQ') %>%
  dplyr::rename(min_fractions = n_fractions) %>%
  mutate(min_pairs = 0) %>%
  group_by(min_pairs, min_fractions) %>%
  summarise(n = n(), auroc = mean(auroc, na.rm = T)) %>%
  ungroup()

# plot
pI = bind_rows(meansIa, meansIb) %>%
  filter(min_fractions <= 10) %>%
  filter(min_fractions >= min_pairs) %>%
  ggplot(aes(x = factor(min_fractions), y = factor(min_pairs), fill = auroc)) +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC ', 
                         breaks = seq(0.665, 0.675, 0.01)) +
  scale_x_discrete('Min. # of fractions', expand = c(0, 0)) +
  scale_y_discrete('Min. # of pairs', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.key.height = unit(0.25, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
pI

###############################################################################-
###### j. LFQ ####
###############################################################################-

# read complexes
datJ = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup()

# read CDF
cdfJ = readRDS("data/QC/protein-groups-CDF.rds") %>%
  # filter to datasets without ratiometric quantitation
  group_by(accession, experiment) %>%
  filter(!any(quant_mode == 'ratio')) %>%
  ungroup()

# first, plot CDF: mean # of protein quantitations from each method
# limit to human/mouse only
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  filter(Species %in% c("Homo sapiens", "Mus musculus"))
mean_cdfJ = cdfJ %>%
  filter(paste(accession, experiment) %in% 
           with(expts, paste(Accession, Replicate))) %>%
  group_by(quant_mode, n_fractions) %>%
  summarise(mean_proteins = mean(n_proteins),
            mean_coverage = mean(coverage)) %>%
  ungroup()
pal = brewer.pal(8, 'Set2') %>% extract(c(2, 1, 4))
pJ1 = mean_cdfJ %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  filter(n_fractions <= 50) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  ggplot(aes(x = n_fractions, y = mean_proteins, color = quant_mode)) +
  geom_line(size = 0.5) +
  # geom_point(size = 0.8) + 
  scale_y_continuous(expression(Proteins~(10^3)), breaks = seq(0, 3000, 500),
                     labels = function(x) x / 1e3) +
  scale_x_continuous('# of fractions') +
  scale_color_manual('', values = pal) +
  boxed_theme() +
  theme(legend.key.size = unit(0.6, "lines"))
pJ1

## next, boxplots
### complexes: CORUM
mediansJ = datJ %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  group_by(quant_mode) %>%
  summarise(median = median(auroc, na.rm = T),
            label = formatC(median, format = 'f', digits = 3))
pJ2 = datJ %>%
  mutate(quant_mode = gsub("_", " ", quant_mode)) %>%
  filter(quant_mode != 'MS1 intensity') %>%
  mutate(quant_mode = fct_recode(quant_mode, 'MaxLFQ' = 'LFQ')) %>%
  ggplot(aes(y = reorder(quant_mode, auroc, stats::median),
             x = auroc, fill = quant_mode, color = quant_mode)) +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = mediansJ, aes(x = 0.47, y = quant_mode, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC') +
  scale_y_discrete(labels = function(x) gsub(" ", "\n", x)) +
  coord_cartesian(xlim = c(0.475, 0.7)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
pJ2

# combine and save
pJ = pJ1 + pJ2 + plot_layout(ncol = 1, heights = c(1, 0.4))
pJ

###############################################################################-
###### combine all plots
###############################################################################-

p_bl = ggplot() + theme_void()
col3 = pHa + pHb + p_bl + pJ1 + pJ2 + 
  plot_layout(ncol = 1, heights = c(1.1, 0.75, 1.25, 1, 0.4))
ggsave("fig/final/fig2/col3.pdf", col3, height = 15, width = 5.6, 
       units = "cm", useDingbats = F)

# save heatmap separately (break grid)
ggsave("fig/final/fig2/col3-2.pdf", pI, height = 4.75, width = 4.75, 
       units = "cm", useDingbats = F)
