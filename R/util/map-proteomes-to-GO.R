setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read proteome2taxid map
map = read.delim("data/GO/proteome2taxid.gz", header = F,
                 col.names = c('species', 'taxid', 'file'))

# list proteomes
proteome_files = list.files("~/git/CF-MS-searches/data/fasta/filtered", 
                            pattern = "*.gz")
proteome_accns = gsub("-.*$", "", proteome_files)

# for each proteome, map to NCBI taxid using UniProt web server
taxid_file = 'data/util/uniprot2taxid.txt'
if (!file.exists(taxid_file)) {
  taxids = character(0)
  for (proteome_accn in proteome_accns) {
    url = paste0(
      "https://www.uniprot.org/proteomes/?sort=&desc=&compress=no&query=", 
      proteome_accn, 
      "%20redundant:no&fil=&limit=10&force=no&preview=true&format=tab",
      "&columns=id,organism-id")
    mapping = read_tsv(url)
    taxid = mapping$`Organism ID`[1]
    message('mapped ', proteome_accn, ' to ', nrow(mapping), ' taxid: ', taxid)
    taxids[proteome_accn] = taxid
  }
  # save this map
  taxid_map = data.frame(proteome = names(taxids), taxid = taxids)
  write.table(taxid_map, taxid_file, sep = '\t', quote = F, row.names = F)
} else {
  # read the map
  taxid_map = read.delim(taxid_file)
  taxids = with(taxid_map, setNames(taxid, proteome))
}

# now, filter GO map based on taxids
map0 = filter(map, taxid %in% taxids)
# download the appropriate GO files
for (file in map0$file) {
  message(file)
  remote = paste0("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/", file)
  local = paste0("data/GO/GOA/", file)
  if (!file.exists(paste0(local, '.gz'))) {
    download.file(url = remote, destfile = local, method = 'curl')
    # gzip it
    system(paste("gzip --force", local))
  }
}

# finally, write a 'master' map for the GO anlaysis:
## UniProt proteome
## NCBI taxid
## GO file name
## species clean name (experiments.csv)
expts = read.csv("~/git/CF-MS-searches/data/experiments.csv")
master_map = taxid_map %>%
  left_join(map, by = 'taxid')
master_map$species %<>%
  strsplit(" ") %>%
  map(head, 2) %>%
  map_chr(paste0, collapse = ' ') %>%
  # fix two that couldn't be auto-mapped
  fct_recode("Cyanothece ATCC 51142" = "Crocosphaera subtropica",
             "Brassica oleracea var. italica" = "Brassica oleracea") %>%
  as.character()
write.csv(master_map, "data/GO/species-map.csv", row.names = F)
