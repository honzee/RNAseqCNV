#' get count diploid reference
get_count_ref <- function(base_matr){

  #create function for geometric mean
  gm_mean = function(a){prod(a)^(1/length(a))}

  #calculate per-gene geometric mean for normalization
  pseudo_ref <- apply(X = base_matr, MARGIN = 1, function(x) gm_mean(x))

  #calculate per-gene median for transformation (diploid centering)
  gene_med <- apply(X = base_matr, MARGIN = 1, function(x) median(x))

  return(list(pseudo_ref = pseudo_ref, gene_med = gene_med))
}
