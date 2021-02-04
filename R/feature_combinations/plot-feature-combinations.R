setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(lme4)
library(lmerTest)
source("R/functions.R")
source("R/theme.R")

# read dataset
auc = readRDS("data/analysis/feature_combinations/AUC.rds") %>%
  # filter to held-out complexes only
  filter(complex_set == 'held-out') %>%
  # use camera ready feature names
  mutate(feature1 = clean_metric(feature1), 
         feature2 = clean_metric(feature2))

# extract individual features
indiv = filter(auc, is.na(feature2))

# extract feature pairs
pairs = filter(auc, !is.na(feature2))

# analyze RF and NB separately
classifiers = unique(auc$classifier)
all_interactions = list()
for (classifier in classifiers) {
  indiv0 = filter(indiv, classifier == !!classifier)
  pairs0 = filter(pairs, classifier == !!classifier)
  
  # order features
  means_indiv = indiv0 %>%
    group_by(classifier, n_datasets, feature1) %>%
    summarise(auc = mean(auc), n = n()) %>%
    ungroup() %>% 
    ## don't go above six datasets
    filter(n_datasets <= 6) 
  lvls = means_indiv %>% arrange(desc(auc)) %>% pull(feature1)
  
  # plot individual features
  range = range(means_indiv$auc)
  p1 = means_indiv %>%
    ggplot(aes(x = factor(n_datasets), y = reorder(feature1, auc),
               # y = factor(feature1, levels = rev(lvls)), 
               fill = auc)) + 
    facet_grid(~ classifier) +
    geom_tile(color = 'white') +
    scale_x_discrete('# of datasets', expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_paletteer_c("pals::coolwarm", name = 'AUC   ',
                           breaks = range, labels = format(range, digits = 2)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.height = unit(0.2, 'lines'),
          legend.key.width = unit(0.25, 'lines'))
  p1
  # save plot
  output_file1 = paste0("fig/analysis/feature_combinations/individual-features-",
                        classifier, ".pdf")
  ggsave(output_file1, p1, width = 8, height = 7.5, units = "cm",
         useDingbats = FALSE)
  
  # plot rank
  ranks_indiv = means_indiv %>%
    group_by(n_datasets) %>%
    arrange(auc) %>%
    mutate(rank = row_number(),
           rank_pct = rank / n(),
           rank_pct = rescale(rank_pct, c(0, 1))) %>%
    ungroup()
  lvls = means_indiv %$% reorder(feature1, auc) %>% levels()
  p1b = ranks_indiv %>%
    ggplot(aes(x = factor(n_datasets), y = factor(feature1, levels = lvls),
               # y = factor(feature1, levels = rev(lvls)), 
               fill = rank_pct)) + 
    facet_grid(~ classifier) +
    geom_tile(color = 'white') +
    scale_x_discrete('# of datasets', expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_viridis(option = "cividis", name = "Rank (%)  ", 
                       breaks = c(0, 1),
                       labels = function(x) 100 * x) +
    # scale_fill_paletteer_c("pals::coolwarm", name = 'AUC   ',
    #                        breaks = range, labels = format(range, digits = 2)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.height = unit(0.2, 'lines'),
          legend.key.width = unit(0.25, 'lines'))
  p1b
  # save plot
  output_file1b = paste0("fig/analysis/feature_combinations/individual-features-",
                        classifier, "-ranks.pdf")
  ggsave(output_file1b, p1b, width = 8, height = 7.5, units = "cm",
         useDingbats = FALSE)
  
  # analyze each n_datasets separately
  for (n_datasets in c(3, 6)) {
    # compute means for pairs
    means = pairs0 %>%
      group_by(classifier, n_datasets, feature1, feature2) %>%
      summarise(auc = mean(auc), n = n()) %>%
      ungroup()
    ## make symmetrical
    means %<>% rbind(
      means %>%
        dplyr::mutate(tmp = feature1, feature1 = feature2, feature2 = tmp) %>%
        dplyr::select(-tmp)
    )
    ## include individual features
    means %<>% rbind(
      means_indiv %>% mutate(feature2 = feature1) %>% 
        filter(n_datasets == !!n_datasets)
    )
    
    # filter
    means0 = filter(means, n_datasets == !!n_datasets) %>%
      mutate(classifier = fct_recode(classifier, 
                                     'Naive Bayes' = 'NB',
                                     'Random forest' = 'RF'),
             facet = paste0(classifier, ", ", n_datasets, " datasets"))
    
    # identify interactions
    # model = pairs0 %>%
    #   filter(n_datasets == !!n_datasets) %>%
    #   lm(auc ~ feature1 * feature2, data = .) %>%
    #   tidy()
    # interactions = filter(model, grepl(":", term)) %>%
    #   separate(term, into = c('feature1', 'feature2'), sep = ':') %>%
    #   mutate_at(vars(feature1, feature2), ~ gsub("feature[12]", "", .)) %>%
    #   mutate(padj = p.adjust(p.value, 'BH'))
    model = pairs0 %>%
      filter(n_datasets == !!n_datasets) %>% 
      mutate(sample_idx = factor(sample_idx)) %>% 
      lmer(auc ~ feature1 * feature2 + (1 | sample_idx), data = .,
           REML = TRUE, control = lmerControl(
             check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
    interactions = coef(summary(model)) %>%
      as.data.frame() %>%
      rownames_to_column('term') %>%
      filter(grepl(":", term)) %>%
      separate(term, into = c('feature1', 'feature2'), sep = ':') %>%
      mutate_at(vars(feature1, feature2), ~ gsub("feature[12]", "", .)) %>%
      mutate(pval = `Pr(>|t|)`,
             padj = p.adjust(pval, 'BH')) %>%
      arrange(padj)
    sig = filter(interactions, padj < 0.05)
    ## make symmetrical
    sig %<>% rbind(
      sig %>%
        dplyr::mutate(tmp = feature1, feature1 = feature2, feature2 = tmp) %>%
        dplyr::select(-tmp)
    )
    ## flag negative vs. positive interactions
    sig %<>%
      mutate(symbol = ifelse(Estimate < 0, '-', '+'))
    sig_lenient = filter(interactions, pval < 0.05) %>%
      arrange(Estimate)
    
    # append to a list
    # sigs[[length(sigs) + 1]] = mutate(sig, classifier = classifier, 
    #                                   n_datasets = n_datasets)
    all_interactions[[length(all_interactions) + 1]] = 
      mutate(interactions, classifier = classifier, n_datasets = n_datasets)
    
    # plot pairs
    range = range(means0$auc)
    lvls = means_indiv %>% 
      filter(n_datasets == !!n_datasets) %$%
      reorder(feature1, auc) %>%
      levels()
    p2 = means0 %>%
      mutate(feature1 = factor(feature1, levels = lvls),
             feature2 = factor(feature2, levels = lvls)) %>% 
      # filter(as.integer(feature2) > as.integer(feature1)) %>%
      ggplot(aes(x = feature1, y = feature2)) + 
      facet_grid(~ facet) +
      geom_tile(aes(fill = auc), color = 'white') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_paletteer_c("pals::coolwarm", name = 'AUC',
                             breaks = range, 
                             labels = formatC(range, format = 'f', digits = 2)) +
      guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
      coord_fixed() +
      boxed_theme() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'right',
            legend.justification = 'bottom',
            legend.key.height = unit(0.25, 'lines'),
            legend.key.width = unit(0.2, 'lines'))
    
    if (nrow(sig) > 0)
      p2 = p2 + 
      geom_text(data = sig, aes(label = symbol), color = 'white', size = 1.75, 
                # nudge_y = -0.15
                )
    p2
    # save plot
    output_file2 = paste0("fig/analysis/feature_combinations/pair-features-",
                          classifier, "-n_datasets=", n_datasets, ".pdf")
    ggsave(output_file2, p2, width = 9.25, height = 9, units = "cm",
           useDingbats = FALSE)
  }
}

# create a master plot of interactions
# master = bind_rows(sigs) %>%
#   dplyr::count(feature1, feature2, symbol)
master = bind_rows(all_interactions) %>%
  mutate(padj = p.adjust(pval, 'BH'),
         symbol = ifelse(Estimate < 0, '-', '+')) %>%
  filter(padj < 0.05) %>% 
  dplyr::count(feature1, feature2, symbol)
## make symmetrical
master %<>% rbind(
  master %>%
    dplyr::mutate(tmp = feature1, feature1 = feature2, feature2 = tmp) %>%
    dplyr::select(-tmp)
)

## are any discordant?
master %>% group_by(feature1, feature2) %>% filter(n_distinct(symbol) > 1)

# order features
lvls = indiv0 %>%
  group_by(n_datasets, feature1) %>%
  summarise(auc = mean(auc), n = n()) %>%
  ungroup() %>% 
  ## don't go above six datasets
  filter(n_datasets <= 6) %$%
  reorder(feature1, auc) %>%
  levels()

pal = pals::kovesi.diverging_cwm_80_100_c22(8) # [c(1, 3, 4, 6)]
p3 = master %>%
  mutate(feature1 = factor(feature1, levels = lvls),
         feature2 = factor(feature2, levels = lvls),
         n = ifelse(symbol == '-', n * -1, n)) %>% 
  ggplot(aes(x = feature1, y = feature2)) + 
  # facet_grid(~ facet) +
  geom_tile(aes(fill = n), color = 'white') +
  geom_text(aes(label = symbol), color = 'white', size = 1.75) +
  scale_x_discrete(expand = c(0, 0), drop = FALSE) +
  scale_y_discrete(expand = c(0, 0), drop = FALSE) +
  # scale_fill_manual(values = pal) +
  scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22",
                         breaks = seq(-4, 4, 1) %>% setdiff(0),
                         labels = abs,
                         name = "Synergistic\ninteractions") +
  guides(fill = guide_legend(title.pos = 'top', title.hjust = 0.5,
                             nrow = 1)) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(0.4, 'lines'),
        legend.key.width = unit(0.4, 'lines'))
p3
ggsave("fig/analysis/feature_combinations/interactions-master.pdf", p3,
       width = 9, height = 9, units = "cm", useDingbats = FALSE)
