setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read all results from individual complexes
dat = readRDS("data/analysis/analysis_grid/individual_complexes.rds")

# tag fractionation
expts = read.csv("~/git/PCPdb/data/experiments.csv")
fract = expts %>%
  dplyr::select(Accession, Replicate, Fractionation) %>%
  set_colnames(c("accession", "experiment", "fractionation"))
dat %<>% left_join(fract)
# ignore cross-linked datasets
dat %<>% filter(!grepl("^XL-", fractionation)) 

# filter to complexes found by all methods, then compute mean
dat1 = dat %>%
  group_by(complex) %>%
  filter(n_distinct(fractionation) == 4) %>%
  ungroup() %>%
  group_by_at(vars(-complex, -n_proteins, -auroc)) %>% 
  summarise(median = median(auroc),
            mean = mean(auroc),
            n_complexes = n_distinct(complex)) %>% 
  ungroup()
medians1 = dat1 %>%
  group_by(fractionation) %>%
  summarise(median = median(mean, na.rm = TRUE),
            label = formatC(median, format = 'f', digits = 3))
arrange(medians1, median)
p1 = dat1 %>%
  ggplot(aes(y = reorder(fractionation, mean, stats::median), 
             x = mean, fill = fractionation, color = fractionation)) +
  # facet_grid(metric ~ ., scales = 'free_y') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians1, aes(x = 0.5, y = fractionation, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  scale_x_continuous('AUC') +
  scale_y_reordered() + 
  coord_cartesian(xlim = c(0.5, 0.73)) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1
ggsave("fig/analysis/analysis_grid/individual-complexes-AUC.pdf", p1,
       width = 4.4, height = 2.45, units = "cm", useDingbats = FALSE)

# identify complexes that are consistently resolved better by one method
complex_means = dat %>%
  group_by(accession, experiment, fractionation, complex) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>% 
  ungroup()
library(broom)
complex_tests = map(unique(dat$fractionation), ~ 
                      complex_means %>%
                      group_by(complex) %>%
                      filter(n_distinct(fractionation) >= 2) %>% 
                      mutate(class = (fractionation == .x)) %>% 
                      do(tidy(pairwise.t.test(x = .$mean, g = .$class,
                                              alternative = 'greater'))) %>% 
                      ungroup()) %>% 
  setNames(unique(dat$fractionation)) %>% 
  bind_rows(.id = 'fractionation') %>%
  mutate(padj = p.adjust(p.value, 'BH'))
labels = complex_tests %>% 
  filter(padj < 0.05) %>% 
  dplyr::count(fractionation)
p2 = complex_tests %>% 
  filter(padj < 0.05) %>% 
  ggplot(aes(x = fractionation, fill = fractionation, color = fractionation)) +
  geom_bar(alpha = 0.6, size = 0.3, width = 0.7) +
  geom_text(data = labels, aes(y = n + 10, label = n), size = 2,
            color = 'black') +
  scale_x_discrete(breaks = c('IEF', 'IEX', 'N-PAGE', 'SEC')) +
  scale_y_continuous('Complexes',
                     limits = c(0, 110), expand = c(0, 0),
                     breaks = seq(0, 110, 25)) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p2
ggsave("fig/analysis/analysis_grid/individual-complexes-n-best-by-method.pdf",
       p2, width = 2.75, height = 3.25, units = "cm", useDingbats = FALSE)

# how many complexes were only identified by one method?
only = dat %>%
  group_by(complex) %>%
  filter(n_distinct(fractionation) == 1) %>% 
  ungroup() %>%
  distinct(complex, fractionation)
labels = only %>%
  mutate(fractionation = factor(fractionation,
                                levels = c('IEF', 'IEX', 'N-PAGE', 'SEC'))) %>% 
  dplyr::count(fractionation, .drop = FALSE)
p3 = only %>%
  mutate(fractionation = factor(fractionation, levels = c('IEF', 'IEX',
                                                          'N-PAGE', 'SEC'))) %>% 
  ggplot(aes(x = fractionation, fill = fractionation, color = fractionation)) +
  # facet_grid(~ metric) +
  geom_bar(alpha = 0.6, size = 0.3, width = 0.7) +
  geom_text(data = labels, aes(y = n + 25, label = n), size = 2,
            color = 'black') +
  scale_x_discrete(breaks = c('IEF', 'IEX', 'N-PAGE', 'SEC'), drop = FALSE) +
  scale_y_continuous('Complexes',
                     limits = c(0, 320), expand = c(0, 0),
                     breaks = seq(0, 300, 100)) +
  scale_fill_manual('', values = pal, guide = F) +
  scale_color_manual('', values = pal, guide = F) +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p3
ggsave("fig/analysis/analysis_grid/individual-complexes-n-specific-by-method.pdf",
       p3, width = 2.75, height = 3.25, units = "cm", useDingbats = FALSE)

# boxplots of complexes detected best by each method
## IEF
filter(complex_tests, fractionation == 'IEF') %>% arrange(p.value)
ief_complex = 'H2AX complex I'
pal = pals::stepped3()[c(1, 5, 10, 15, 19)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
labels = complex_means %>%
  filter(complex == ief_complex) %>% 
  group_by(fractionation) %>% 
  summarise(mean = mean(mean),
            median = median(mean),
            label = round(median, digits = 3)) 
box = filter(labels, fractionation == 'IEF')
p4a = complex_means %>%
  filter(complex == ief_complex) %>% 
  mutate(facet = ief_complex) %>% 
  ggplot(aes(x = fractionation, y = mean, color = fractionation,
             fill = fractionation)) +
  facet_grid(~ facet) +
  geom_crossbar(data = box, aes(ymin = -Inf, ymax = +Inf),
                fill = 'grey90', color = NA, alpha = 0.33, width = 1.0) +
  geom_boxplot(alpha = 0.5, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.5) +
  geom_text(data = labels, aes(label = label, y = 0.76), color = 'black',
            vjust = 1, size = 1.75) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous('AUC') +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7))
p4a

## IEX
filter(complex_tests, fractionation == 'IEX') %>% arrange(p.value) %>% head(20)
iex_complex = "INO80 chromatin remodeling complex"
pal = pals::stepped3()[c(1, 5, 10, 15, 19)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
labels = complex_means %>%
  filter(complex == iex_complex) %>% 
  group_by(fractionation) %>% 
  summarise(mean = mean(mean),
            median = median(mean),
            label = round(median, digits = 3)) 
box = filter(labels, fractionation == 'IEX')
p4b = complex_means %>%
  mutate(fractionation = factor(fractionation,
                                levels = c("IEF", "IEX", "N-PAGE", "SEC"))) %>%
  filter(complex == iex_complex) %>% 
  mutate(facet = iex_complex) %>% 
  ggplot(aes(x = fractionation, y = mean, color = fractionation,
             fill = fractionation)) +
  facet_grid(~ facet) +
  geom_crossbar(data = box, aes(ymin = -Inf, ymax = +Inf),
                fill = 'grey90', color = NA, alpha = 0.33, width = 1) +
  geom_boxplot(alpha = 0.5, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.5) +
  geom_text(data = labels, aes(label = label, y = 0.74), color = 'black',
            vjust = 1, size = 1.75) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous('AUC') +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7))
p4b

## N-PAGE
filter(complex_tests, fractionation == 'N-PAGE') %>% arrange(p.value)
page_complex = 'Coat protein complex II (SAR1B, PREB, SEC23A , SEC24B, SEC23IP)'
pal = pals::stepped3()[c(1, 5, 10, 15, 19)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
labels = complex_means %>%
  filter(complex == page_complex) %>% 
  group_by(fractionation) %>% 
  summarise(mean = mean(mean),
            median = median(mean),
            label = round(median, digits = 3)) 
box = filter(labels, fractionation == 'N-PAGE')
p4c = complex_means %>%
  mutate(fractionation = factor(fractionation,
                                levels = c("IEF", "IEX", "N-PAGE", "SEC"))) %>% 
  filter(complex == page_complex) %>% 
  mutate(facet = page_complex %>% gsub(" \\(.*$", "", .)) %>% 
  ggplot(aes(x = fractionation, y = mean, color = fractionation,
             fill = fractionation)) +
  facet_grid(~ facet) +
  geom_crossbar(data = box, aes(ymin = -Inf, ymax = +Inf),
                fill = 'grey90', color = NA, alpha = 0.33, width = 1) +
  geom_boxplot(alpha = 0.5, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.5) +
  geom_text(data = labels, aes(label = label, y = 0.79), color = 'black',
            vjust = 1, size = 1.75) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous('AUC') +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7))
p4c

## SEC
filter(complex_tests, fractionation == 'SEC') %>% arrange(p.value)
# sec_complex = "CCR4-NOT complex"
# sec_complex = "PHAX-CBC complex (cap binding complex)"
# sec_complex = "Dynein complex"
# sec_complex = "CCR4-NOT-CNOT7-CNOT6 complex"
# sec_complex = "F1F0-ATP synthase, mitochondrial"
sec_complex = 'Ribosome, cytoplasmic'
pal = pals::stepped3()[c(1, 5, 10, 15, 19)] %>%
  setNames(c("SEC", "N-PAGE", "IEX", "IEF", "Sucrose"))
labels = complex_means %>%
  filter(complex == sec_complex) %>% 
  group_by(fractionation) %>% 
  summarise(mean = mean(mean),
            median = median(mean),
            label = round(median, digits = 3)) 
box = filter(labels, fractionation == 'SEC')
p4d = complex_means %>%
  filter(complex == sec_complex) %>% 
  mutate(facet = sec_complex) %>% 
  ggplot(aes(x = fractionation, y = mean, color = fractionation,
             fill = fractionation)) +
  facet_grid(~ facet) +
  geom_crossbar(data = box, aes(ymin = -Inf, ymax = +Inf),
                fill = 'grey90', color = NA, alpha = 0.33, width = 1) +
  geom_boxplot(alpha = 0.5, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.5) +
  geom_text(data = labels, aes(label = label, y = 0.82), color = 'black', # ifelse(fractionation %in% c('IEF', 'N-PAGE'), 0.815, 0.85)), color = 'black',
            vjust = 1, size = 1.75) +
  scale_y_continuous('AUC') +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7))
p4d

# combine all
p4 = p4a + p4b + p4c + p4d + plot_layout(ncol = 4)
p4
ggsave("fig/analysis/analysis_grid/individual-complexes-examples.pdf", p4,
       width = 14, height = 5.5, units = "cm", useDingbats = FALSE)
