/*******************************************************************************
MODULE NAME
rdahmm_setup

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
#include "rdahmm.h"

/* this program's header */
#include "rdahmm_setup.h"

int handle_rdahmm_input(int argc, char **argv, arg *args, 
                        char **output_type_names, char **init_type_names, 
                        char **reg_type_names, char **data_file, char **L_file,                         char **Q_file, char **pi_file, char **A_file, 
                        char **B_file, char **min_val_file, char **max_val_file,
                        char **range_file, int *T, int *D, int *N, 
                        int *output_type, int *init_type, double *thresh, 
                        double *peps, int *regularize, int *reg_type, 
                        double *omega_Q1, double *omega_Q2, double *omega_Q3, 
                        double *omega_sigma, int *anneal, double *anneal_step, 
                        double *anneal_factor, double *beta_min, int *eval_only,                        int *add_state, int *viterbi, int *weighted_covars, 
                        char **covars_weights_file, int *use_covgraph, 
                        char **covgraph_file, int *ntries, int *maxiters, 
                        int *seed, int *verbose,
                        double ***cont_data, double **min_val, double **max_val,
                        double **range)
{
  int handle_args_ok;
  int continuous;

  /* Get command line arguments */
  handle_args_ok = handle_args(argc, argv, args, output_type_names, 
                               init_type_names, reg_type_names, data_file, 
                               L_file, Q_file, pi_file, A_file, B_file,
                               min_val_file, max_val_file, range_file,
                               T, D, N, output_type, init_type, thresh, 
                               peps, regularize, reg_type, omega_Q1, 
                               omega_Q2, omega_Q3, omega_sigma, anneal,
                               anneal_step, anneal_factor, beta_min,
                               eval_only, add_state, viterbi,
                               weighted_covars, covars_weights_file,
                               use_covgraph, covgraph_file,
                               ntries, maxiters, seed, verbose);

  if (handle_args_ok != UT_OK)
  {
    exit(UT_ERROR);
  }

  /* read in the data matrix */
  if (*output_type == RDAHMM_GAUSS_OUTPUT)
  {
    read_dmatrix(*data_file, *T, *D, cont_data);
    continuous = UT_TRUE;
  }

  /* in the case of continuous data, normalize each dimension to a 0-1 scale */

  if (continuous == UT_TRUE)
  {
    if (*eval_only == UT_TRUE)
    {
      /* read range information from a file */
      read_drow(*min_val_file, *D, min_val);
      read_drow(*max_val_file, *D, max_val);
      read_drow(*range_file, *D, range);
    }
    else
    {
      /* create storage for the normalization */
      *min_val = NR_dvector(1, *D);
      *max_val = NR_dvector(1, *D);
      *range   = NR_dvector(1, *D);
      
      /* determine range based on data */
      range_of_dcols(*cont_data, *min_val, *max_val, *range, *T, *D);
    }

    range_normalize_dcols(*cont_data, *min_val, *range, *T, *D);
  }

  return(UT_OK);
}

int allocate_gauss_hmm(Gaussian_HMM *hmm, int N, int D)
{
  (*hmm).N = N;
  (*hmm).D = D;
  (*hmm).pi = NR_dvector(1, N);
  (*hmm).A = NR_dmatrix(1, N, 1, N);
  (*hmm).mu = NR_dmatrix(1, N, 1, D);
  (*hmm).sigma = set_of_dmatrices(N, D, D);

  return(UT_OK);
}
  
int setup_gauss_hmm(Gaussian_HMM *hmm, 
                    int add_state, int eval_only, int weighted_covars,
                    int use_covgraph,
                    double *min_val, double *range,
                    double **cont_data, int N, int D, int T,
                    char *pi_file, char *A_file, char *B_file,
                    char *covars_weights_file, char *covgraph_file,
                    int **covgraph)
{
  int i, j, t, d, d2;
  double *covars_weights;

  /* allocate memory for the HMM */
  if (add_state == UT_FALSE)
  {
    allocate_gauss_hmm(hmm, N, D);
  }
  else
  {
    allocate_gauss_hmm(hmm, N+1, D);
  }

  /* if evaluating sequence based on a model */
  if (eval_only == UT_TRUE)
  {
    /* load the model parameters */

    /* load initial and state-to-state transition probabilities */
    read_init_probs(pi_file, N, (*hmm).pi);
    read_transition_matrix(A_file, N, (*hmm).A);

    read_gauss_output_params(B_file, (*hmm).mu, (*hmm).sigma, N, D);
    if (weighted_covars == UT_TRUE)
    {
      read_dcol(covars_weights_file, D, &covars_weights);
    }

    /* normalize the output distributions */
    range_normalize_dcols((*hmm).mu, min_val, range, N, D);
    for (i = 1; i <= N; i++)
    {
      /* put this in a function later */
      for (d = 1; d <= D; d++)
      {
        for (d2 = 1; d2 <= D; d2++)
        {
          (*hmm).sigma[i][d][d2] /= (range[d] * range[d2]);
          if (weighted_covars == UT_TRUE)
          {
            (*hmm).sigma[i][d][d2] *= (covars_weights[d] * covars_weights[d2]);
          }
        }
      }
    }

    /* if adding a blank state, modify the model */
    if (add_state == UT_TRUE)
    {
      /* adjust existing values in initial and transition probabilities */
      for (i = 1; i <= N; i++)
      {
        (*hmm).pi[i] *= ((double) N) / ((double) (N+1));
        for (j = 1; j <= N; j++)
        {
          (*hmm).A[i][j] *= ((double) N) / ((double) (N+1));
        } 
      }

      /* add values for new state in initial and transition probilities */
      (*hmm).pi[N+1] = 1.0 / ((double) (N+1));
      for (i = 1; i <= N; i++)
      {
        (*hmm).A[i][N+1] = 1.0 / ((double) (N+1)); 
      }
      for (j = 1; j <= N; j++)
      {
        (*hmm).A[N+1][j] = 1.0 / ((double) (N+1)); 
      }
      (*hmm).A[N+1][N+1] = 1.0 / ((double) (N+1)); 
     
      /* add output distribution for new state */
      /* calculate mean and covariance of all the data */
      for (d = 1; d <= D; d++)
      {
        (*hmm).mu[N+1][d] = 0.0;
        for (t = 1; t <= T; t++)
        {
          (*hmm).mu[N+1][d] += cont_data[t][d];
        }
          (*hmm).mu[N+1][d] /= ((double) T);
      }
      for (d = 1; d <= D; d++)
      {
        for (d2 = 1; d2 <= D; d2++)
        {
          (*hmm).sigma[N+1][d][d2] = 0.0;
          for (t = 1; t <= T; t++)
          {
            (*hmm).sigma[N+1][d][d2] += 
              (cont_data[t][d] - (*hmm).mu[N+1][d]) * 
              (cont_data[t][d2] - (*hmm).mu[N+1][d2]);
          }
          (*hmm).sigma[N+1][d][d2] /= (range[d] * range[d2]);
        }
      }
    }
  }

  /* handle user-defined constraints on the output distribution */
  if (use_covgraph == UT_TRUE)
  {
    read_icol(covgraph_file, D, covgraph);
  }
  else
  {
    *covgraph = NR_ivector(1, D);
    (*covgraph)[1] = 0;
  }

  if (weighted_covars == UT_TRUE)
  {
    NR_free_dvector(covars_weights, 1, D);
  }

  return(UT_OK);
}
        
int setup_rdaem_params(RDAEM_Parameters *params, double thresh, double peps,
                       int regularize, int init_type, int reg_type,
                       double omega_Q1, double omega_Q2, double omega_Q3,
                       double omega_sigma, int anneal, double anneal_step,
                       double anneal_factor, double beta_min, int *covgraph,
                       int ntries, int maxiters, int seed)
{
  (*params).thresh = thresh;
  (*params).peps = peps;
  (*params).regularize = regularize;
  (*params).init_type = init_type;
  (*params).reg_type = reg_type;
  (*params).omega_Q1 = omega_Q1;
  (*params).omega_Q2 = omega_Q2;
  (*params).omega_Q3 = omega_Q3;
  (*params).omega_sigma = omega_sigma;
  (*params).anneal = anneal,
  (*params).anneal_step = anneal_step;
  (*params).anneal_factor = anneal_factor;
  (*params).beta_min = beta_min;
  (*params).covgraph = covgraph;
  (*params).ntries = ntries;
  (*params).maxiters = maxiters;
  (*params).seed = seed;
  
  return(UT_OK);
}

int rdahmm_setup(int argc, char **argv, arg *args, 
                 char **output_type_names, char **init_type_names, 
                 char **reg_type_names, char **data_file, char **L_file,
                 char **Q_file, char **pi_file, char **A_file, 
                 char **B_file, char **min_val_file, char **max_val_file,
                 char **range_file, int *T, int *D, int *N, 
                 int *output_type, int *init_type, double *thresh, 
                 double *peps, int *regularize, int *reg_type, 
                 double *omega_Q1, double *omega_Q2, double *omega_Q3, 
                 double *omega_sigma, int *anneal, double *anneal_step, 
                 double *anneal_factor, double *beta_min, int *eval_only,
                 int *add_state, int *viterbi, int *weighted_covars, 
                 char **covars_weights_file, int *use_covgraph, 
                 char **covgraph_file, int **covgraph, int *ntries, 
                 int *maxiters, int *seed, int *verbose,
                 double ***cont_data, double **min_val, double **max_val,
                 double **range,
                 Gaussian_HMM *hmm, 
                 RDAEM_Parameters *params)
{
  handle_rdahmm_input(argc, argv, args, output_type_names, init_type_names, 
                      reg_type_names, data_file, L_file,
                      Q_file, pi_file, A_file, B_file, min_val_file, 
                      max_val_file, range_file, T, D, N, 
                      output_type, init_type, thresh, 
                      peps, regularize, reg_type, 
                      omega_Q1, omega_Q2, omega_Q3, 
                      omega_sigma, anneal, anneal_step, 
                      anneal_factor, beta_min, eval_only,
                      add_state, viterbi, weighted_covars, 
                      covars_weights_file, use_covgraph, 
                      covgraph_file, ntries, maxiters, 
                      seed, verbose,
                      cont_data, min_val, max_val,
                      range);

  setup_gauss_hmm(hmm, *add_state, *eval_only, *weighted_covars, 
                  *use_covgraph, *min_val, *range, *cont_data, *N, *D, *T,
                  *pi_file, *A_file, *B_file, *covars_weights_file,
                  *covgraph_file, covgraph);

  setup_rdaem_params(params, *thresh, *peps, *regularize, *init_type,
                     *reg_type, *omega_Q1, *omega_Q2,
                     *omega_Q3, *omega_sigma, *anneal,
                     *anneal_step, *anneal_factor,
                     *beta_min, *covgraph, *ntries, *maxiters, *seed);

  return(UT_OK);
}

int free_gauss_hmm(Gaussian_HMM hmm)
{
  NR_free_dvector(hmm.pi, 1, hmm.N);
  NR_free_dmatrix(hmm.A, 1, hmm.N, 1, hmm.N);
  NR_free_dmatrix(hmm.mu, 1, hmm.N, 1, hmm.D);
  free_set_of_dmatrices(hmm.sigma, hmm.N, hmm.D, hmm.D);

  return(UT_OK);
}

int free_rdaem_params(RDAEM_Parameters params, Gaussian_HMM hmm)
{
  NR_free_ivector(params.covgraph, 1, hmm.D);
}
