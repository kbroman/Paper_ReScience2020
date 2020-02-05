#' Number of insertion sites in each gene in M tb CDC1551
#'
#' Number of insertion sites in the initial 80\% of each gene in the
#' _Mycobacterium tuberculosis_ CDC1551 genome.
#'
#'
#' @name Mtb80
#'
#' @format A matrix with two columns.  Each row corresponds to a gene.  (The
#' row names are the MT numbers of the genes.)  The element in the first column
#' is the number of transposon insertion sites in the initial 80\% that appear
#' in the corresponding gene and in no other gene.  The element in the second
#' column is the number of transposon insertion sites in the initial 80\% of
#' both that gene and the following gene.  There are 4204 rows; the 46 genes
#' with no such site are not included.
#'
#' @seealso [negenes::negenes()], [negenes::sim.mutants()]
#'
#' @references Blades, N. J. and Broman, K. W. (2002) Estimating the number of
#' essential genes in a genome by random transposon mutagenesis.  Technical
#' Report MS02-20,Department of Biostatistics, Johns Hopkins University,
#' Baltimore, MD.
#' <https://www.biostat.wisc.edu/~kbroman/publications/ms0220.pdf>
#'
#' @source <http://www.tigr.org>
#'
#' @keywords datasets
#'
#' @examples
#' \dontrun{data(Mtb80)
#'
#' # simulate 44% of genes to be essential
#' essential <- rep(0,nrow(Mtb80))
#' essential[sample(1:nrow(Mtb80),ceiling(nrow(Mtb80)*0.44))] <- 1
#'
#' # simulate 759 mutants
#' counts <- sim.mutants(Mtb80[,1], essential, Mtb80[,2], 759)
#'
#' # run the Gibbs sampler
#' output <- negenes(Mtb80[,1], counts[,1], Mtb80[,2], counts[,2])}
#'
NULL
