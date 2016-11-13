/*******************************************************************************
MODULE NAME
rdahmm_control

AUTHOR
Robert Granat

DESCRIPTION
See Robert Granat PhD thesis for details and notation.

COMPILE
make

NOTES
Not yet available

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id$";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log$
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* CP library */
#include "cp.h"

/* UT library */
#include "ut.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* local header files */
#include "rdahmm_arguments.h"
#include "rdahmm_gauss_output.h"
#include "rdahmm_generic.h"
#include "rdahmm_setup.h"
#include "rdahmm.h"

/* this program's header */
#include "rdahmm_control.h"

int test_gauss_output_hmm(double **cont_data, Gaussian_HMM hmm,
                          int *Q, double *L, int T, int viterbi, int verbose)
{
   int core_ok;

   core_ok = eval_gauss_output_hmm(cont_data, L, Q, hmm.pi, hmm.A, hmm.mu, 
                                   hmm.sigma, hmm.N, T, hmm.D, viterbi, 
                                   verbose);

   return(core_ok);
}

int train_gauss_output_hmm(double **cont_data, Gaussian_HMM hmm, 
                           RDAEM_Parameters params, int *Q, double *L,
                           int T, int viterbi, int verbose)
{
  int core_ok;
  int i, n;
  double best_L;
  int *best_Q;
  double *best_pi;
  double **best_A;
  double **best_mu;
  double ***best_sigma;

  for (n = 1; n <= params.ntries; n++)
  {
    /* seed the random number generator */

    if (params.seed == UT_FALSE)
    {
      set_rand_by_time_of_day();
    }
    else
    {
      set_rand(params.seed);
    }

    /* initialize model */
    init_gauss_output_hmm(params.init_type, cont_data, hmm.pi, hmm.A, hmm.mu, 
                          hmm.sigma, hmm.N, hmm.D, T);

    /* learn hmm model */
    core_ok = learn_gauss_output_hmm(cont_data, 
                                     params.regularize,
                                     params.reg_type,
                                     params.omega_Q1,
                                     params.omega_Q2,
                                     params.omega_Q3,
                                     params.omega_sigma,
                                     params.anneal,
                                     params.anneal_step,
                                     params.anneal_factor,
                                     params.beta_min,
                                     params.peps,
                                     params.thresh,
                                     params.covgraph,
                                     params.maxiters,
                                     L, Q, 
                                     hmm.pi, hmm.A, hmm.mu, hmm.sigma, 
                                     hmm.N, T, hmm.D,
                                     viterbi, verbose);

    if (params.ntries > 1)
    {
      if (n == 1) /* do initialization */
      {
        best_L = -DBL_MAX;
        best_Q = NR_ivector(1, T);
        best_pi = NR_dvector(1, hmm.N);
        best_A = NR_dmatrix(1, hmm.N, 1, hmm.N);
        best_mu = NR_dmatrix(1, hmm.N, 1, hmm.D);
        best_sigma = set_of_dmatrices(hmm.N, hmm.D, hmm.D);
      }

      if (*L > best_L)
      {
        best_L = *L;
        copy_ivec(Q, best_Q, T);
        copy_dvec(hmm.pi, best_pi, hmm.N);
        copy_dmat(hmm.A, best_A, hmm.N, hmm.N);
        copy_dmat(hmm.mu, best_mu, hmm.N, hmm.D);
        for (i = 1; i <= hmm.N; i++)
        {
          copy_dmat(hmm.sigma[i], best_sigma[i], hmm.D, hmm.D);
        }
      }

      if (n == params.ntries)
      {
        *L = best_L;
        copy_ivec(best_Q, Q, T);
        copy_dvec(best_pi, hmm.pi, hmm.N);
        copy_dmat(best_A, hmm.A, hmm.N, hmm.N);
        copy_dmat(best_mu, hmm.mu, hmm.N, hmm.D);
        for (i = 1; i <= hmm.N; i++)
        {
            copy_dmat(best_sigma[i], hmm.sigma[i], hmm.D, hmm.D);
        }
      }
    }
  }

  if (params.ntries > 1)
  {
    NR_free_ivector(best_Q, 1, T);
    NR_free_dvector(best_pi, 1, hmm.N);
    NR_free_dmatrix(best_A, 1, hmm.N, 1, hmm.N);
    NR_free_dmatrix(best_mu, 1, hmm.N, 1, hmm.D);
    free_set_of_dmatrices(best_sigma, hmm.N, hmm.D, hmm.D);
  }

  return(UT_OK);
}

