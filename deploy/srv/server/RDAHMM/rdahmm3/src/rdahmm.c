/*
Copyright 2006, by the California Institute of Technology. ALL RIGHTS RESERVED.
United States Government Sponsorship acknowledged. Any commercial use must be 
negotiated with the Office of Technology Transfer at the California Institute 
of Technology.

This software may be subject to U.S. export control laws. By accepting this 
software, the user agrees to comply with all applicable U.S. export laws and 
regulations. User has the responsibility to obtain export licenses, or other 
export authority as may be required before exporting such information to 
foreign countries or providing access to foreign persons.
*/

/*******************************************************************************
PROGRAM NAME
rdahmm

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
static char rcsid[] = "$Id: rdahmm.c,v 1.2 2003/08/27 18:10:37 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm.c,v $
 * Revision 1.2  2003/08/27 18:10:37  granat
 * removed hacks related to OSX symbol problems (solved by compiling libut
 * with libtool -c)
 *
 * Revision 1.1  2003/08/27 18:02:25  granat
 * Initial revision
 *
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
#include "rdahmm_setup.h"
#include "rdahmm_control.h"
#include "rdahmm_gauss_output.h"
#include "rdahmm_generic.h"

/* this program's header */
#include "rdahmm.h"

/*******************************************************************************
 MAIN
 Driver of program.
*******************************************************************************/
int main(int argc, char *argv[])

{
  /* arguments, used by libcp */
  arg     args[RDAHMM_NUM_ARGS+1];

  /* filenames */
  char    *data_file;
  char    *L_file;
  char    *Q_file;
  char    *pi_file;
  char    *A_file;
  char    *B_file;
  char    *min_val_file;
  char    *max_val_file;
  char    *range_file;
  char    *covars_weights_file;
  char    *covgraph_file;
  
  /* main program input */
  double  **cont_data;     /* for continuous input data */
  int     *discrete_data;  /* for discrete input data */
  int     T;               /* number of observations */
  int     D;               /* dimension of the observations */
  int     N;               /* number of states in the model (model size) */

  /* output distribution type */
  char    *output_type_names[RDAHMM_NUM_OUTPUT_TYPES+1];
  int     output_type;

  /* main program output */
  double  L;               /* log likelihood of the model */
  int     *Q;              /* state sequence */
  double  *pi;             /* vector of initial state probabilities */
  double  **A;             /* matrix of transition probabilities */

  /* main program output for Gaussian output distributions */
  double  **mu;
  double  ***sigma;

  /* secondary program input */
  double  thresh;          /* convergence threshold */
  double  peps;            /* perturbation epsilon */
  int     regularize;      /* flag to perform regularized learning */
  double  omega_Q1;        /* regularization weighting term */
  double  omega_Q2;        /* regularization weighting term */
  double  omega_Q3;        /* regularization weighting term */
  double  omega_sigma;     /* regularization weighting term (Gaussian output) */
  int     anneal;          /* flag to perform deterministic annealing */
  double  anneal_step;     /* temperature step for deterministic annealing */
  double  anneal_factor;   /* temperature factor for deterministic annealing */
  double  beta_min;        /* initial computational temperature for DA */
  int     eval_only;       /* flag for evaluating state assignments only */
  int     add_state;       /* flag to add a "blank" state to loaded model */
  int     weighted_covars; /* flag to weight covariance matrix in evaluation */
  int     use_covgraph;    /* covariance dimension connectivity */
  int     ntries;          /* number of times to retry the optimization */
  int     maxiters;        /* maximum number of EM iterations */
  int     seed;            /* seed for initializing random number generator */
  int     viterbi;         /* flag for performing viterbi state calculation */
  int     verbose;         /* flag for printing verbose output */

  /* storage for values over multiple retries */
  double  best_L;               /* log likelihood of the model */
  int     *best_Q;              /* state sequence */
  double  *best_pi;             /* vector of initial state probabilities */
  double  **best_A;             /* matrix of transition probabilities */

  /* storage over multiple retries for Gaussian output distributions */
  double  **best_mu;
  double  ***best_sigma;

  /* model initialization methods */
  char    *init_type_names[RDAHMM_NUM_INIT_TYPES+1];
  int     init_type;

  /* regularization methods */
  char    *reg_type_names[RDAHMM_NUM_REG_TYPES+1];
  int     reg_type;

  /* miscellaneous variables */
  int     continuous;
  double  *min_val;
  double  *max_val;
  double  *range;
  double  *covars_weights;
  int     *covgraph;
  FILE    *fp;
  int     t, i, j, d, d2, n;
  int     Nplus1;

  /* error checking */
  int     handle_args_ok;
  int     core_ok;

  /* structures */
  Gaussian_HMM hmm;
  RDAEM_Parameters params;

/*----------------------------------------------------------------------------*/
  /* Get command line arguments, set up data structures, and read input */

  rdahmm_setup(argc, argv, args, output_type_names, init_type_names, 
               reg_type_names, &data_file, &L_file, &Q_file, &pi_file, 
               &A_file, &B_file, &min_val_file, &max_val_file, &range_file,
               &T, &D, &N, &output_type, &init_type, &thresh, &peps, 
               &regularize, &reg_type, &omega_Q1, &omega_Q2, &omega_Q3, 
               &omega_sigma, &anneal, &anneal_step, &anneal_factor, &beta_min,
               &eval_only, &add_state, &viterbi, &weighted_covars, 
               &covars_weights_file, &use_covgraph, &covgraph_file, &covgraph, 
               &ntries, &maxiters, &seed, &verbose, &cont_data, &min_val, 
               &max_val, &range, &hmm, &params);

/*----------------------------------------------------------------------------*/
  /* Allocate memory for the underlying state sequence 
     (Note that this varies by data, and is not a model parameter) */

  Q = NR_ivector(1, T);

/*----------------------------------------------------------------------------*/
  /* if evaluating sequence based on a model */

  if (eval_only == UT_TRUE)
  {
    if (output_type == RDAHMM_GAUSS_OUTPUT)
    {
      test_gauss_output_hmm(cont_data, hmm, Q, &L, T, viterbi, verbose);
    }

    write_icol(Q_file, T, Q, "w");

    /* write the log likelihood file */
    fp = fopen(L_file, "w");
    fprintf(fp, "%g\n", L);
    fclose(fp);
  }

/*----------------------------------------------------------------------------*/
  /* otherwise learn the model and write model parameters */

  if (eval_only == UT_FALSE)
  {
    if (output_type == RDAHMM_GAUSS_OUTPUT)
    {
      train_gauss_output_hmm(cont_data, hmm, params, Q, &L, T,
                             viterbi, verbose);

      /* must unnormalize before writing model */
      range_unnormalize_dcols(hmm.mu, min_val, range, N, D);
      for (i = 1; i <= N; i++)
      {
        range_unnormalize_cov_dmatrix(hmm.sigma[i], D, range);
      }
    }

    /* write outputs */

    /* write hmm initial and transition probabilities */
    write_dcol(pi_file, N, hmm.pi, "w");
    write_dmatrix(A_file, N, N, hmm.A, "w");

    /* write state output distribution parameters */
    if (output_type == RDAHMM_GAUSS_OUTPUT)
    {
      write_gauss_output_params(B_file, hmm.mu, hmm.sigma, N, D);
    }

    /* write the state sequence file */
    write_icol(Q_file, T, Q, "w");

    /* write the log likelihood file */
    fp = fopen(L_file, "w");
    fprintf(fp, "%g\n", L);
    fclose(fp);

    /* write the range information for the data */
    write_drow(min_val_file, D, min_val, "w");
    write_drow(max_val_file, D, max_val, "w");
    write_drow(range_file, D, range, "w");
  }

/*----------------------------------------------------------------------------*/
  /* clean up memory */

  /* observations */
  if (continuous == UT_TRUE)
  {
    NR_free_dmatrix(cont_data, 1, T, 1, D);

    NR_free_dvector(min_val, 1, D);
    NR_free_dvector(max_val, 1, D);
    NR_free_dvector(range, 1, D);
  }

  /* optimization parameters */
  free_rdaem_params(params, hmm);

  /* model */
  if (output_type == RDAHMM_GAUSS_OUTPUT)
  {
    free_gauss_hmm(hmm);
  }

  /* state sequence */
  NR_free_ivector(Q, 1, T);

  /* filenames */
  free(data_file);
  free(L_file);
  free(Q_file);
  free(pi_file);
  free(A_file);
  free(B_file);

  return(UT_OK);
}
