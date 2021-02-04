setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read Gaussians
gauss = readRDS("data/analysis/co_apex/all_gaussians.rds") %>% 
  filter(quant_mode == 'iBAQ')

# plot number of Gaussians per protein group in each experiment/replicate
n_gauss = dplyr::count(gauss, experiment, replicate, quant_mode, protein)
pal = brewer.pal(5, 'Blues')
order = n_gauss %>%
  group_by(experiment, replicate) %>%
  summarise(pct_1 = mean(n == 1)) %>%
  ungroup() %>%
  arrange(desc(pct_1)) %$%
  paste(experiment, replicate)
p1 = n_gauss %>%
  unite(xval, experiment, replicate, sep = " ") %>%
  ggplot(aes(x = factor(xval, levels = order),
             fill = factor(n, levels = seq(5, 1, -1)))) +
  geom_segment(aes(x = -Inf, xend = -Inf, y = 0, yend = 1), color = 'grey50',
               size = 0.4) +
  geom_bar(color = 'grey15', size = 0.2, position = 'fill',
           width = 0.7) +
  scale_fill_manual('# of Gaussians', values = pal,
                    breaks = seq_len(5)) +
  scale_y_continuous('% of chromatograms', labels = function(x) x * 100) +
  guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  coord_flip() +
  boxed_theme() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.45, 'lines'))
p1
ggsave("fig/analysis/co_apex/n-gaussians.pdf", p1, 
       width = 6, height = 7.5, units = "cm", useDingbats = FALSE)
