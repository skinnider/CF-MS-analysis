setwd("~/git/CF-MS-analysis")
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read experiments
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
  # recode fractionation slightly
  mutate(fract = fct_recode(Fractionation, 
                            'N-PAGE' = 'XL-N-PAGE',
                            'SEC' = 'XL-SEC'))

# plot
pal = BuenColors::jdb_palette("solar_extra") %>% extract(-c(1, 9)) %>%
  colorRampPalette() %>% do.call(list(5))
df = expts %>%
  dplyr::count(fract) %>%
  mutate(fract = reorder(fract, -n))
p = df %>%
  ggplot(aes(x = fract, y = n, fill = fract, color = fract)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_text(aes(y = n + 6, label = n), size = 2, color = 'black') +
  scale_x_discrete() +
  scale_y_continuous('Experiments', limits = c(0, 100), expand = c(0, 0),
                     breaks = seq(0, 90, 30)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p
ggsave("fig/QC/fractionation-summary-bar.pdf", p, width = 3.5, height = 4.2,
       units = "cm", useDingbats = F)
