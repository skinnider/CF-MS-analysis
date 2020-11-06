load_networks = function() {
  library(flavin)
  library(PrInCE)
  
  # load CORUM
  complexes = read.delim("~/git/network-validation/data/complex/CORUM/complexes_human.txt") %>%
    as_annotation_list('gene_name', 'complex')
  adj = adjacency_matrix_from_list(complexes)
  
  # create master list
  networks = list()
  
  # 1. final CF-MS network #### 
  message("reading CF-MS consensus interactome ...")
  params = list(
    classifier = 'RF',
    n_features = 1,
    feature_select = 'best_first'
  )
  net1_file = file.path(base_dir, 'human_interactome',
                        paste0('network', 
                               '-classifier=', params$classifier,
                               '-n_features=', params$n_features,
                               '-feature_select=', params$feature_select,
                               '-n_datasets=', 46, '-sample_idx=', 1, '.rds'))
  net1 = readRDS(net1_file) %>%
    # threshold to 50% precision
    mutate(labels = make_labels(adj, .),
           precision = calculate_precision(labels)) %>%
    PrInCE::threshold_precision(0.5)
  networks[['CF-MS']] = net1
  
  # 2. human high-throughput interactomes (n=5) ####
  message("reading human high-throughput interactomes ...")
  hts_files = list.files("~/git/network-validation/data/high_throughput",
                         pattern = '*\\.gz$', full.names = TRUE)
  hts_nets = map(hts_files, read.delim) %>%
    setNames(gsub("\\..*$", "", basename(hts_files)) %>%
               paste0('high-throughput|', .))
  networks %<>% c(hts_nets)
  
  # return
  message("done")
  return(networks)
}
