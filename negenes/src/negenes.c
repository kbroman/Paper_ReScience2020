/**********************************************************************
 *
 * negenes.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "negenes.h"

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
 *               [vector of length n_genes+2]
 *               element 0 = element n_genes; element 1 = element n_genes+1
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
	       int *trace, double *startp)
{
  int **Output;

  reorg_output(*n_mcmc, *n_genes, output, &Output);

  GetRNGstate(); /* get random number seed */

  /* shift n_sites2 and curstate */
  n_sites2++;
  curstate++;

  negenes(*n_genes, *n_mutants, n_sites, n_sites2, known, *n_mcmc,
	  *burnin, *skip, Output, n_ess, geneprob, curstate, *n_w,
	  w, *calcprob, logprob, *saveoutput, *trace, *startp);

  PutRNGstate(); /* put random number seed */
}

/* The actual function */
void negenes(int n_genes, int n_mutants, int *n_sites, int *n_sites2,
	     int *known, int n_mcmc, int burnin, int skip, int **Output,
	     int *n_ess, double *geneprob, int *curstate, int n_w,
	     int *w, int calcprob, double *logprob, int saveoutput,
	     int trace, double startp)
{
  int i, j, s, oldstate;
  int cursum, curnum;
  double p;

  /* start point of Gibbs sampler */
  if(startp<1e-10) {
    for(i=0; i<n_genes; i++)
      curstate[i] = 1;  /* state = 1 -> non-essential */
  }
  else if(startp > 1.0 - 1e-10) {
    for(i=0; i<n_genes; i++) {
      if(known[i]==1) curstate[i] = 1;
      else curstate[i] = 0; /* state = 0 -> essential */
    }
  }
  else { /* random start */
    for(i=0; i<n_genes; i++) {
      if(known[i]==1 || unif_rand() > startp)
	curstate[i] = 1;
      else curstate[i] = 0;
    }
  }
  curstate[n_genes] = curstate[0];
  curstate[-1] = curstate[n_genes-1];

  /* get initial counts */
  cursum = curnum = 0;
  for(i=0; i<n_genes; i++) {
    curnum += curstate[i];
    cursum += (curstate[i]*n_sites[i] +
	       curstate[i]*curstate[i+1]*n_sites2[i]);
  }

  /* begin Gibbs iterations */
  for(i=-burnin; i<n_mcmc; i++) {
    for(s=0; s<=skip; s++) {

      /* consider genes in random order */
      int_permute(w, n_w);

      for(j=0; j<n_w; j++) {

	p = gibbsProb(curstate[w[j]], curstate[w[j]-1], curstate[w[j]+1],
		      n_sites[w[j]], n_sites2[w[j]], n_sites2[w[j]-1],
		      curnum, cursum, n_mutants, n_genes);

	oldstate = curstate[w[j]];

	/* simulate new state */
	if(unif_rand() < p)  curstate[w[j]] = 1;
	else curstate[w[j]] = 0;

	/* make curstate wrap around */
	if(w[j]==0) curstate[n_genes] = curstate[0];
	else if(w[j]==n_genes-1) curstate[-1] = curstate[n_genes-1];

	/* update curnum and cursum */
	curnum += (curstate[w[j]]-oldstate);

	cursum += (curstate[w[j]]-oldstate)*
	  (n_sites[w[j]] + curstate[w[j]+1]*n_sites2[w[j]] +
	   curstate[w[j]-1]*n_sites2[w[j]-1]);

      } /* end loop over genes */
    } /* end loop over skips */

    /* if after burnin period, save result */
    if(i >= 0) {
      n_ess[i] = 0;

      for(j=0; j<n_genes; j++) {
	/* Calculate Gibbs prob for Rao-Blackwellized ests */
	if(known[j]) p=1.0; /* mutant observed => non-essential */
	else p = gibbsProb(curstate[j],curstate[j+1],curstate[j-1],
			   n_sites[j], n_sites2[j], n_sites2[j-1],
			   curnum, cursum, n_mutants, n_genes);

	n_ess[i] += curstate[j];
	geneprob[j] += (1-p);
	if(saveoutput) Output[j][i] = curstate[j];
      }
      n_ess[i] = n_genes - n_ess[i]; /* current number essential */

      if(calcprob) {
	logprob[i] = -(double)n_mutants * log((double)cursum);
	logprob[i] += lgammafn((double)curnum+1) + lgammafn((double)(n_genes-curnum+1));
      }
    }

    if(trace && (i+1 == ((i+1)/20)*20))  /* print iteration */
      Rprintf("%4d\n", i+1);

  } /* end loop over gibbs iterations */

  /* update geneprob */
  for(j=0; j<n_genes; j++) geneprob[j] /= (double)n_mcmc;
}



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
		 int curnum, int cursum, int n_mutants, int n_genes)
{
  int Ai, Bi;
  double p;

  Ai = curnum - curstate;
  Bi = cursum - curstate*n_sites - curstate*curstate_next*n_sites2 -
    curstate_prev*curstate*n_sites2_prev;

  p = Bi + n_sites + n_sites2*curstate_next + n_sites2_prev*curstate_prev;

  p = exp(-(double)n_mutants * (log((double)p) - log((double)Bi)));

  /* p = Pr(gene is non-essential */
  return(p / (p + (double)(n_genes - Ai) / (double)(Ai + 1)));
}


/* reorganize output vector so it is a doubly indexed array rather
   than a single long vector; afterwards, indexed as
   Output[gene][iteration] */
void reorg_output(int n_mcmc, int n_notobs, int *output, int ***Output)
{
  int i;

  *Output = (int **)R_alloc(n_notobs, sizeof(int *));
  (*Output)[0] = output;
  for(i=1; i<n_notobs; i++)
    (*Output)[i] = (*Output)[i-1] + n_mcmc;

}

/* same as reorg_output, but for doubles
   than a single long vector; afterwards, indexed as
   Output[gene][iteration] */
/***************
** This is no longer needed **
void reorg_output_dbl(int n_mcmc, int n_notobs, double *output, double ***Output)
{
  int i;

  *Output = (double **)R_alloc(n_notobs, sizeof(double *));
  (*Output)[0] = output;
  for(i=1; i<n_notobs; i++)
    (*Output)[i] = (*Output)[i-1] + n_mcmc;

}
****************/

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

void int_permute(int *array, int len)
{
  int i, which;
  int tmp;

  for(i=0; i < len; i++) {
    which = random_int(i, len-1);
    tmp = array[which];
    array[which] = array[i];
    array[i] = tmp;
  }
}

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

int random_int(int low, int high)
{
  return((int)(unif_rand()*(double)(high - low + 1)) + low);
}


/* end of negenes.c */



