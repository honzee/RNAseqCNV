#' Centromeric regions for genome build hg38.
#'
#' A table containg information on centromeric regions for each chromosomes from genome build hg38 extracted from Table Browser.
#' Group - Mapping and sequencing, Track - Centromeres, Table - centromeres. This produdes multiple centromeric regions overlapping/non-overlapping
#' on a single chromosome. Starting centromeric location was selected as the one with lowest number, end location as the one with the highest number.
#'
#' @format A data with 24 rows and 3 columns:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{cstart}{centromere starting genome location}
#'   \item{cend}{centromere end genome location}
#' }
"centromeres_hg38"
