#' Simulate data for a random transposon mutagenesis experiment
#'
#' Simulate data for a random transposon mutagenesis experiment.
#'
#'
#' @param n.sites A vector specifying the number of transposon insertion sites
#' in each gene.  All elements must by strictly positive.
#' @param essential A vector containing 1's (indicating that the corresponding
#' gene is essential) and 0's (indicating that the corresponding gene is not
#' essential). Must be the same length as `n.sites`.
#' @param n.sites2 A vector specfying the number of transposon insertion sites
#' shared by adjacent genes.  The *i*th element is the number of insertion
#' sites shared by genes *i* and *i*+1.  The last element is for
#' sites shared by genes *N* and 1.  If missing, these are assumed to be
#' all 0.
#' @param n.mutants Number of mutants to simulate.
#'
#' @return If `n.sites2` is missing or contains all 0's, a vector is
#' returned containing the number of mutants observed for each gene.
#'
#' If `n.sites2` is not missing and has some positive entries, a matrix
#' with two columns is returned.  The first column contains the number of
#' mutants observed for each gene alone; the second column contains the number
#' of mutants observed shared by adjacent genes.
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso [negenes::negenes()], [negenes::Mtb80()]
#'
#' @references Blades, N. J. and Broman, K. W. (2002) Estimating the number of
#' essential genes in a genome by random transposon mutagenesis.  Technical
#' Report MS02-20, Department of Biostatistics, Johns Hopkins University,
#' Baltimore, MD.
#' <https://www.biostat.wisc.edu/~kbroman/publications/ms0220.pdf>
#'
#' @keywords datagen
#'
#' @export
#'
#' @examples
#'
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
sim.mutants <-
function(n.sites, essential, n.sites2=NULL, n.mutants)
{
  n.genes <- length(n.sites)
  if(length(essential) != n.genes)
    stop("n.sites and essential must be the same length")

  if(is.null(n.sites2)) n.sites2 <- rep(0,n.genes)

  if(length(n.sites2) != n.genes)
    stop("n.sites and n.sites2 must be the same length")

  if(any(essential != 0 & essential != 1))
    stop("essential must contain only 0's and 1's.")

  if(n.mutants <= 0)
    stop("n.mutants must be positive")

  temp <- c(essential[-1],essential[1])
  p <- c(n.sites*(1-essential), n.sites2*(1-essential)*(1-temp))

  o <- table(factor(sample(1:(2*n.genes), n.mutants, replace=TRUE,
                           prob=p), levels=1:(2*n.genes)))
  names(o) <- NULL

  if(sum(n.sites2)==0) return(o[1:n.genes])
  else return(cbind(o[1:n.genes],o[-(1:n.genes)]))
}
