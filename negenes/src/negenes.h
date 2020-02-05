/**********************************************************************
 *
 * negenes.h
 *
 * copyright (c) 2002, Karl W Broman
 * last modified August, 2002
 * first written June, 2002
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at https://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/negenes package
 *
 * Contains: R_negenes, negenes, reorg_output, int_permute, random_int
 *
 **********************************************************************/

/**********************************************************************
 *
 * negenes
 *
 * Gibbs sampler to estimate the posterior distribution of the number of
 * essential genes in a genome with data from a random transposon
 * mutagenesis experiment.
 *
 * INPUT:
 *     n_genes = Number of different genes
 *     n_mutants = Number of observed viable mutants
 *     n_sites = Number of transposon insertion sites in each gene
 *               [vector of length n_genes]
 *     n_sites2 = Number of transposon insertion sites shared by adjacent genes.
 *               [vector of length n_genes]
 *     known = Vector of length n_genes; 1 = mutant observed; 0 = not
 *     n_mcmc  = Number of Gibbs iterations to perform (really, to save)
 *     burnin  = Number of initial Gibbs steps, with results discarded
 *               (really, burnin*(skip+1) initial steps are performed)
 *     skip    = Only return every skip+1st step
 *     output  = 1/0 matrix of dimension (n_genes x n_mcmc)
 *     n_ess   = vector of length n_mcmc
 *     geneprob = Rao-Blackwellized estimates of prob. essential for each
 *                gene [length = n_genes]
 *
 *     curstate  = Workspace of integers, of length n_genes+2
 *     n_w = n_genes - sum(known)
 *     w =         Workspace of length n_w, containing indices of
 *                 genes not known to be non-essential
 *     calcprob  = if 1, calculate the log posterior (up to a scalar multiple)
 *     logprob   = vector of length n.mcmc, to contain log posterior
 *     saveoutput = if 1, fill up all of output; if 0, don't.
 *
 *     trace     = if 1, print trace information
 *     startp    Initial proportion of genes for which no mutant was observed
 *               that will be assumed essential for the Gibbs sampler.
 *
 **********************************************************************/

/* Wrapper for R */
void R_negenes(int *n_genes, int *n_mutants, int *n_sites, int *n_sites2,
	       int *known, int *n_mcmc, int *burnin, int *skip, int *output,
	       int *n_ess, double *geneprob, int *curstate, int *n_w,
	       int *w, int *calcprob, double *logprob, int *saveoutput,
	       int *trace, double *startp);

/* The actual function */
void negenes(int n_genes, int n_mutants, int *n_sites, int *n_sites2,
	     int *known, int n_mcmc, int burnin, int skip, int **Output,
	     int *n_ess, double *geneprob, int *curstate, int n_w,
	     int *w, int calcprob, double *logprob, int saveoutput,
	     int trace, double startp);

/**********************************************************************
 *
 * gibbsProb: Calculate Gibbs step probability
 *
 * curstate      = current state at current position
 * curstate_next = current state at next position
 * curstate_prev = current state at prev position
 * n_sites       = no. sites for current position
 * n_sites2      = no. sites shared by current and next position
 * n_sites2_prev = no. sites shared by current and prev position
 * curnum        = current number of non-essential genes
 * cursum        = current number of viable sites
 * n_mutants     = number of mutants observed
 * n_genes       = number of genes
 **********************************************************************/
double gibbsProb(int curstate, int curstate_next, int curstate_prev,
		 int n_sites, int n_sites2, int n_sites2_prev,
		 int curnum, int cursum, int n_mutants, int n_genes);

/* reorganize output vector so it is a doubly indexed array rather
   than a single long vector; afterwards, indexed as
   Output[gene][iteration] */
void reorg_output(int n_mcmc, int n_notobs, int *output, int ***Output);

/* same as reorg_output, but for doubles
   than a single long vector; afterwards, indexed as
   Output[gene][iteration] */
/***** no longer needed
void reorg_output_dbl(int n_mcmc, int n_notobs, double *output, double ***Output);
*****/

/**********************************************************************
 *
 * int_permute
 *
 *   This function randomly permutes a vector of integers
 *
 * Input:
 *
 *   array = vector of ints; on output, it contains a random
 *           permutation of the input vector
 *
 *   len   = length of the vector
 *
 **********************************************************************/

void int_permute(int *array, int len);

/**********************************************************************
 *
 * random_int
 *
 * Generates a random int integer between "low" and "high", inclusive.
 *
 *  Input:
 *
 *    low
 *
 *    high
 *
 **********************************************************************/

int random_int(int low, int high);

/* end of negenes.h */



