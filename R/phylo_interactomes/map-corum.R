# Map CORUM protein complexes to eggNOG (euk) accessions.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)

# read CORUM
complexes = read.delim("~/git/network-validation/data/complex/CORUM/coreComplexes.txt.gz") %>%
  dplyr::select(ComplexName, subunits.UniProt.IDs.) %>%
  set_colnames(c("complex", "uniprot")) %>%
  mutate(uniprot = strsplit(uniprot, ';')) %>%
  unnest(uniprot)

# read relevant eggNOG maps
eggnog_files = list.files("data/resources/eggNOG", full.names = TRUE,
                          pattern = "annotations\\.gz") %>%
  extract(grepl("musculus|norvegicus|sapiens", .))
orthos = map(eggnog_files, read.delim, comment.char = '#', header = FALSE,
             col.names = c('query_name', 'seed_eggNOG_ortholog',
                           'seed_ortholog_evalue', 'seed_ortholog_score',
                           'predicted_gene_name', 'GO_terms', 
                           'KEGG_KOs', 'BiGG_reactions',
                           'Annotation_tax_scope', 'OGs', 
                           'bestOG|evalue|score', 'COG cat',
                           'eggNOG annot')) 
ortho = orthos %>%
  bind_rows() %>%
  dplyr::select(query_name, OGs) %>%
  mutate(uniprot = map_chr(strsplit(query_name, '\\|'), 2)) %>%
  # get both KOGs and euk OGs
  mutate(euk = strsplit(OGs, ',') %>%
           map(~ extract(., grepl('euNOG', .) & !grepl("KOG|COG", .))),
         kog = strsplit(OGs, ',') %>%
           map(~ extract(., grepl('euNOG', .) & grepl("KOG", .)))) %>%
  # discard multi-mapping genes
  filter(lengths(euk) == 1 | lengths(kog) == 1) %>%
  # prefer KOG to euk OG
  mutate(eggNOG = map2(euk, kog, ~ ifelse(length(.y) > 0, .y, .x))) %>%
  # remove any remaining multi-mapping genes
  filter(lengths(eggNOG) == 1) %>%
  mutate(eggNOG = unlist(eggNOG)) %>%
  dplyr::select(-euk, -kog) %>%
  distinct(uniprot, eggNOG)

# map CORUM to eggNOG
mapped = complexes %>%
  left_join(ortho, by = 'uniprot') %>%
  drop_na() 
## confirm there are no multi-mapping proteins
mapped %>% 
  group_by(uniprot) %>%
  mutate(n_euk = n_distinct(eggNOG[!grepl("KOG", eggNOG)]),
         n_kog = n_distinct(eggNOG[grepl("KOG", eggNOG)])) %$%
  table(n_euk, n_kog)

# filter to complexes with at least two proteins
mapped %<>%
  distinct(complex, eggNOG) %>%
  group_by(complex) %>%
  filter(n_distinct(eggNOG) >= 2) %>%
  ungroup()

# save the mapped protein complexes
write.table(mapped, "data/analysis/phylo_interactomes/complexes_euk.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
