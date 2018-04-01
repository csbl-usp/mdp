#' Expression data example
#'
#' Author expression data for GEO dataset GSE17156 of transcriptome blood samples from patients that were inoculated with the RSV virus that 
#' has been altered by collapsing for HGNC gene symbols.
#'
#' @format A data frame with 13838 rows and 40 variables:
#' @description 
#' \describe{
#'   \item{rownames}{HGNC gene names}
#'   \item{colnames}{sample expression data}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17156}
"example_data"

#' Phenotypic data example
#'
#' Subset of the annotation data for GEO dataset GSE17156, using only patients that have been inoculated with the RSV virus
#'
#' @format A data frame with 40 rows and 2 variables:
#' \describe{
#'   \item{Sample}{GSM identified}
#'   \item{Class}{Symtpomatic state}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17156}
"example_pheno"


#' Sample score results
#'
#' Resultant sample scores when the mdp is applied to example_data and example_pheno
#'
#' @format A data frame with 40 rows and 3 variables:
#' \describe{
#'   \item{Sample}{GSM identified}
#'   \item{Score}{Sample score}
#'   \item{Class}{Symtpomatic state}
#'   ...
#' }
#'
"sample_data"
