#' Estimate the number of essential genes in a genome
#'
#' Estimate, via a Gibbs sampler, the posterior distribution of the number of
#' essential genes in a genome with data from a random transposon mutagenesis
#' experiment. (See the technical report cited below.)
#'
#' @param n.sites A vector specifying the number of transposon insertion sites
#' in each gene (alone).  All elements must by strictly positive.
#' @param counts A vector specifying the number of mutants observed for each
#' gene (alone).  Must be the same length as `n.sites`, and all elements
#' must be non-negative integers.
#' @param n.sites2 A vector specfying the number of transposon insertion sites
#' shared by adjacent genes.  The *i*th element is the number of insertion
#' sites shared by genes *i* and *i*+1.  The last element is for
#' sites shared by genes *N* and 1. If NULL, assume all are 0.
#' @param counts2 A vector specfying the number of mutants shared by adjacent
#' gene (analogous to `n.sites2`). The *i*th element is the number of
#' mutants at sites shared by genes *i* and *i*+1.  The last element
#' is for sites shared by genes *N* and 1. If NULL, assume all are 0.
#' @param n.mcmc Number of Gibbs steps to perform.
#' @param skip An integer; only save every `skip` + 1st step.
#' @param burnin Number of initial Gibbs steps to run (output discarded).
#' @param startp Initial proportion of genes for which no mutant was observed
#' that will be assumed essential for the Gibbs sampler.  (Genes for which a
#' mutant was observed are assumed non-essential; other genes are assumed
#' essential independent with this probability.)
#' @param trace If TRUE, print iteration number occassionally.
#' @param calc.prob If TRUE, return the log posterior probability (up to an
#' additive constant) for each saved iteration.
#' @param return.output If TRUE, include detailed Gibbs results in the output.
#'
#' @return A list with components `n.essential` (containing the total
#' number of essential genes at each iteration of the Gibbs sampler)
#' `summary` (a vector containing the estimated mean, SD, 2.5 percentile
#' and 97.5 percentile of the posterior distribution of the number of essential
#' genes.
#'
#' The next component, `geneprob`, is a vector with one element for each
#' gene, containing the estimated posterior probability that each gene is
#' essential.  These are Rao-Blackwellized estimates.
#'
#' If the argument `calc.prob` was true, there will also be a component
#' `logprob` containing the log (base e) of the posterior probability (up
#' to an additive constant) at each Gibbs step.
#'
#' If the argument `return.output` was true, there will also be a matrix
#' with `n.mcmc` / (`skip` + 1) rows (corresponding to the Gibbs
#' steps) and a column for each gene The entries in the matrix are either 0
#' (essential gene) or 1 (non-essential gene) according to the state of that
#' gene at that step in the Gibbs sampler.
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @importFrom stats quantile sd
#' @export
#' @useDynLib negenes, .registration=TRUE
#'
#' @seealso [negenes::sim.mutants()], [negenes::Mtb80()]
#'
#' @references
#' - Blades, N. J. and Broman, K. W. (2002) Estimating the number of
#'   essential genes in a genome by random transposon mutagenesis.  Technical
#'   Report MS02-20, Department of Biostatistics, Johns Hopkins University,
#'   Baltimore, MD.
#'   <https://www.biostat.wisc.edu/~kbroman/publications/ms0220.pdf>
#'
#' - Lamichhane et al. (2003) A post-genomic method for predicting
#'   essential genes at subsaturation levels of mutagenesis:
#'   application to Mycobacterium, tuberculosis. Proc Natl Acad Sci USA
#'   100:7213-7218
#'   [doi:10.1073/pnas.1231432100](https://doi.org/10.1073/pnas.1231432100)
#'
#' @keywords models
#'
#' @examples
#' data(Mtb80)
#'
#' # simulate 44% of genes to be essential
#' essential <- rep(0,nrow(Mtb80))
#' essential[sample(1:nrow(Mtb80),ceiling(nrow(Mtb80)*0.44))] <- 1
#'
#' # simulate 759 mutants
#' counts <- sim.mutants(Mtb80[,1], essential, Mtb80[,2], 759)
#'
#' # run the Gibbs sampler without returning detailed output
#' \dontrun{output <- negenes(Mtb80[,1], counts[,1], Mtb80[,2], counts[,2])}
#' \dontshow{output <- negenes(Mtb80[,1], counts[,1], Mtb80[,2], counts[,2],
#'                             n.mcmc=100, skip=0, burnin=0)}
#' # run the Gibbs sampler, returning the detailed output
#' \dontrun{output2 <- negenes(Mtb80[,1], counts[,1], Mtb80[,2], counts[,2], return=TRUE)}
#' \dontshow{output2 <- negenes(Mtb80[,1], counts[,1], Mtb80[,2], counts[,2], return=TRUE,
#'                             n.mcmc=100, skip=0, burnin=0)}
#'
negenes <-
function(n.sites, counts, n.sites2=NULL, counts2=NULL,
         n.mcmc=5000, skip=49, burnin=500,
         startp=1, trace=TRUE,
         calc.prob=FALSE, return.output=FALSE)
{
  n.genes <- length(n.sites)

  # check for errors in the input
  if(length(counts) != n.genes)
    stop("n.sites and counts must be the same length")

  if(is.null(n.sites2)) n.sites2 <- rep(0,n.genes)
  if(is.null(counts2)) counts2 <- rep(0,n.genes)
  if(length(n.sites2) != n.genes)
    stop("n.sites2 and n.sites must be the same length")
  if(length(counts2) != n.genes)
    stop("counts2 and n.sites must be the same length")

  if(any(n.sites<0) || any(counts<0) || any(n.sites2<0) || any(counts2<0))
    stop("n.sites, counts, n.sites2, and counts2 must all be >= 0")

  # replace n.mcmc with number of iterations to be saved
  n.mcmc <- ceiling(n.mcmc/(skip+1))

  # update burnin, since using skip
  burnin <- ceiling(burnin/(skip+1))

  if(n.mcmc <= 0 || burnin < 0 || skip < 0)
    stop("n.mcmc, burnin or skip are incorrectly chosen")

  # genes known to be non-essential
  known <- as.numeric(((counts > 0) | (counts2 > 0) |
                       (c(counts2[n.genes],counts2[-n.genes])>0)))
  notknown <- ((1:n.genes)-1)[known==0]
  n.known <- sum(known)
  n.notknown <- n.genes-n.known

  if(startp < 0 || startp > 1)
    stop("startp must be between 0 and 1")

  temp <- as.numeric(return.output)

  n.mutants <- sum(counts) + sum(counts2)

  output <- .C("R_negenes",
               as.integer(n.genes),
               as.integer(n.mutants),
               as.integer(n.sites),
               # n.sites2 made so that first == last and last == first
               as.integer(c(n.sites2[length(n.sites2)],n.sites2,n.sites2[1])),
               as.integer(known),
               as.integer(n.mcmc),
               as.integer(burnin),
               as.integer(skip),
               output = as.integer(rep(0,(n.mcmc-1)*n.genes*temp+n.genes)),
               n.ess = as.integer(rep(0,n.mcmc)),
               geneprob = as.double(rep(0,n.genes)),
               as.integer(rep(0,n.genes+2)),
               as.integer(n.notknown),
               as.integer(notknown),
               as.integer(calc.prob),
               logprob = as.double(rep(0,n.mcmc)),
               as.integer(return.output),
               as.integer(trace),
               as.double(startp),
               PACKAGE="negenes")

  logprob <- output$logprob
  tot.ess <- output$n.ess
  geneprob <- output$geneprob
  if(return.output) {
    output <- matrix(output$output,ncol=n.genes)
    storage.mode(output) <- "integer"
  }
  summ <- c(mean=sum(geneprob),sd=sd(tot.ess),
            quantile(tot.ess,c(0.025,0.975)))

  if(return.output) {
    if(!calc.prob)
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,output=output))
    else
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,
                  logprob=logprob, output=output))
  }
  else {
    if(!calc.prob)
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob))
    else
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,logprob=logprob))
  }
}
