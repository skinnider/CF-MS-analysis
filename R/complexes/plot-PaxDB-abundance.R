# Plot the abundance of human or Arabidopsis proteins detected once, or in 
# at least a certain number of fractions, compared to proteins that were
# never detected.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read or generate PaxDB file
abundance_file = "data/analysis/abundance/coverage.rds"
if (file.exists(abundance_file)) {
  # read complex protein coverage
  coverage = readRDS("data/analysis/abundance/coverage.rds")
} else {
  # read all human and mouse datasets
  expts = read.csv("~/git/CF-MS-searches/data/experiments.csv") %>%
    filter(Species %in% c("Homo sapiens", "Mus musculus"))
  
  # read PaxDB files
  paxdb_files = list.files("data/resources/PaxDB", 
                           pattern = '*WHOLE_ORGANISM*', full.names = T)
  paxdbs = map(paxdb_files, read.delim, comment.char = '#',
               header = F, col.names = c('paxid', 'id', 'abundance')) %>%
    setNames(c("Mus musculus", "Homo sapiens"))
  
  # map them to UniProt IDs
  id_files = paste0("~/git/network-validation/data/identifier/", 
                    c("MOUSE_10090", "HUMAN_9606"), "_idmapping.dat.gz")
  maps = map(id_files, read_tsv, col_names = c("uniprot", "db", "id")) %>%
    setNames(c("Mus musculus", "Homo sapiens"))
  
  coverage = data.frame()
  for (idx in seq_len(nrow(expts))) {
    row = expts[idx, ]
    expt_dir = with(row, file.path("~/git/CF-MS-searches/data/chromatograms", 
                                   Accession, Replicate))
    message("[", idx, "/", nrow(expts), "] processing experiment: ", expt_dir)
    
    chrom_files = list.files(expt_dir, pattern = '*.rds', full.names = T) %>%
      extract(basename(.) != 'metadata.rds') %>%
      setNames(gsub("\\..*$", "", basename(.)))
    
    # get PaxDB file for this species
    species = row$Species
    paxdb = paxdbs[[species]]
    
    # map to UniProt
    map = maps[[species]]
    paxdb %<>%
      mutate(id = gsub("^(10090|9606)\\.", "", id)) %>%
      left_join(map, by = 'id') %>%
      drop_na(uniprot, abundance) %>%
      # keep only one protein ID per gene
      group_by(paxid) %>%
      arrange(desc(abundance)) %>%
      dplyr::slice(1) %>%
      ungroup() %>% 
      distinct(uniprot, abundance)
    
    # read metadata (mapping to genes)
    meta_file = file.path(expt_dir, 'metadata.rds')
    if (!file.exists(meta_file)) {
      warning("file does not exist: ", meta_file)
      next
    }
    meta = readRDS(meta_file)
    
    # process each quantitation strategy independently
    for (quant_mode in names(chrom_files)) {
      chrom_file = chrom_files[quant_mode]
      chrom = readRDS(chrom_file)
      
      # map protein groups to PaxDB identifiers
      paxdb_proteins = meta %>%
        mutate(protein = strsplit(`Majority protein IDs`, ';')) %>%
        unnest(protein) %>%
        # filter to protein complex genes
        filter(protein %in% paxdb$uniprot)
      
      # count number of fractions each PaxDB protein was found in
      n_fractions = rep(NA, nrow(paxdb)) %>%
        setNames(paxdb$uniprot)
      complex_mat = chrom[paxdb_proteins$`Majority protein IDs`, ]
      n_fractions[paxdb_proteins$protein] = rowSums(!is.na(complex_mat) & 
                                                      complex_mat > 0)
      
      # append 
      res = data.frame(quant_mode = quant_mode,
                       gene = paxdb$uniprot,
                       n_fractions = n_fractions) %>%
        cbind(row, .)
      coverage %<>% bind_rows(res)
    }
  }
  
  # replace NAs with zeroes
  coverage %<>% replace_na(list(n_fractions = 0))
  
  # save
  saveRDS(coverage, "data/analysis/abundance/coverage.rds")
}

# group by proteins
proteins = coverage %>%
  group_by(Species, gene) %>%
  summarise(detected = any(n_fractions > 0),
            fractions = sum(n_fractions),
            experiments = n_distinct(paste(Accession, Replicate)[n_fractions > 0])) %>%
  ungroup()

# how many proteins (and what %) were never detected?
proteins %>%
  group_by(Species) %>%
  summarise(n = sum(!detected), pct = mean(!detected))

# merge in protein abundance data
paxdb_files = list.files("data/resources/PaxDB", 
                         pattern = '*WHOLE_ORGANISM*', full.names = T)
paxdbs = map(paxdb_files, read.delim, comment.char = '#',
             header = F, col.names = c('paxid', 'id', 'abundance')) %>%
  setNames(c("Mus musculus", "Homo sapiens"))
id_files = paste0("~/git/network-validation/data/identifier/", 
                  c("MOUSE_10090", "HUMAN_9606"), "_idmapping.dat.gz")
maps = map(id_files, read_tsv, col_names = c("uniprot", "db", "id")) %>%
  setNames(c("Mus musculus", "Homo sapiens"))
abundance = map2_dfr(paxdbs, maps, ~ {
  mutate(.x, id = gsub("^(3702|10090|9606)\\.", "", id)) %>%
    left_join(.y, by = 'id') %>%
    drop_na(uniprot, abundance) %>%
    # keep only one protein ID per gene
    group_by(paxid) %>%
    arrange(desc(abundance)) %>%
    dplyr::slice(1) %>%
    ungroup() %>% 
    distinct(uniprot, abundance)
}) %>%
  dplyr::rename(gene = uniprot)
proteins %<>% left_join(abundance, by = 'gene')

# plot detected vs. never detected
pal = Gpal[c(4, 1)]
pal = c('grey80', cividis(n = 10)[1])
fills = c('grey80', BuenColors::jdb_palette("solar_extra")[2]) %>% rev()
colors = c('grey50', BuenColors::jdb_palette("solar_extra")[1]) %>% rev()
p1 = prot1 %>%
  filter(abundance > 0) %>%
  ggplot(aes(x = abundance, fill = detected, color = detected)) +
  geom_histogram(position = 'identity', alpha = 0.5, size = 0.3) + 
  scale_x_log10("Protein abundance (ppm)", labels = fancy_scientific) +
  scale_y_continuous(expression(Proteins~(10^2)), labels = function(x) x / 1e2,
                     expand = c(0, 0), limits = c(0, 1750)) +
  scale_fill_manual('', values = fills, breaks = c("TRUE", "FALSE"),
                    labels = c("Detected by CF-MS", "Never detected")) +
  scale_color_manual('', values = colors,
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Detected by CF-MS", "Never detected")) +
  boxed_theme() +
  theme(legend.position = c(0.98, 1.04),
        legend.justification = c(1, 1),
        legend.key.size = unit(0.45, 'lines'),
        axis.line = element_blank())
p1
ggsave("fig/analysis/abundance/PaxDB-Homo-sapiens.pdf", p1,
       width = 4.75, height = 4.2, units = "cm", useDingbats = F)

## repeat for mouse
pal = c('grey80', plasma(n = 10)[1])
p2 = proteins %>%
  filter(Species == 'Mus musculus') %>%
  filter(abundance > 0) %>%
  ggplot(aes(x = abundance, fill = detected, color = detected)) +
  geom_histogram(position = 'identity', alpha = 0.5, size = 0.3) + 
  scale_x_log10("Protein abundance (ppm)", labels = fancy_scientific) +
  scale_y_continuous("Proteins", expand = c(0, 0), limits = c(0, 1750)) +
  scale_fill_manual('', values = rev(pal), breaks = c("TRUE", "FALSE"),
                    labels = c("Detected by CF-MS", "Never detected")) +
  scale_color_manual('', values = c(plasma(n = 10)[1], 'grey50'),
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Detected by CF-MS", "Never detected")) +
  boxed_theme() +
  theme(legend.position = c(0.98, 1.04),
        legend.justification = c(1, 1),
        legend.key.size = unit(0.45, 'lines'))
p2
ggsave("fig/analysis/abundance/PaxDB-Mus-musculus.pdf", p2,
       width = 6, height = 5, units = "cm", useDingbats = F)
