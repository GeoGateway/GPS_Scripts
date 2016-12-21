/*******************************************************************************
MODULE NAME
rdahmm_generic

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
static char rcsid[] = "$Id: rdahmm_generic.c,v 1.1 2003/08/27 18:04:04 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_generic.c,v $
 * Revision 1.1  2003/08/27 18:04:04  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* CP library */
#include "cp.h"

/* UT library */
#include "ut.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* local header files */

/* this program's header */
#include "rdahmm_generic.h"

/*******************************************************************************
 READ_INIT_PROBS
 Read the initial state probabilities from a file.
 RG
*******************************************************************************/
int read_init_probs(char *pi_file, int N, double *pi)
{
  FILE   *fp;
  int    i;
  double doubleval;

  /* open file for reading */
  fp = fopen(pi_file, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector */
  for (i = 1; i <= N; i++)
  {
    fscanf(fp, "%lg ", &doubleval);
    pi[i] = doubleval;
  }

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
 READ_TRANSITION_MATRIX
 Read the state-to-state transition probabilities from a file.
 RG
*******************************************************************************/
int read_transition_matrix(char *A_file, int N, double **A)
{
  FILE   *fp;
  int    i, j;
  double doubleval;

  /* open file for reading */
  fp = fopen(A_file, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (i = 1; i <= N; i++)
  {
    for (j = 1; j <= N; j++)
    {
      fscanf(fp, "%lg ", &doubleval);
      A[i][j] = doubleval;
    }
  }

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
 CALC_LOG_LIKELIHOOD
*******************************************************************************/
double calc_log_likelihood(double *scale, int T)
{
  int    t;
  double L;

  L = 0.0;
  for (t = 1; t <= T; t++)
  {
    L -= log(scale[t] + DBL_MIN);
  }

  return(L);
}

/*******************************************************************************
 CALC_Q-FUNCTION
*******************************************************************************/
double calc_Q_function(double **tau_it, double ***tau_ijt, double *pi,
                       double **A, double **prob_data, int N, int T)
{
  int    i, j, t; 
  double Qfunc;
  double log_A_ij;

  Qfunc = 0.0;
  
  for (i = 1; i <= N; i++)
  {
    Qfunc += tau_it[i][1] * log(pi[i] + DBL_MIN); 
    for (t = 1; t <= T; t++)
    {
      Qfunc += tau_it[i][t] * log(prob_data[t][i] + DBL_MIN);
    }
  }

  for (i = 1; i <= N; i++)
  {
    for (j = 1; j <= N; j++)
    {
      log_A_ij = log(A[i][j] + DBL_MIN);
      for (t = 1; t < T; t++)
      {
        Qfunc += tau_ijt[i][j][t] * log_A_ij;
      }
    }
  }

  return(Qfunc);
}

/*******************************************************************************
 UPDATE_INIT_PROBS
*******************************************************************************/
int update_init_probs(double *pi, double **tau_it, double omega_Q1, 
                      int regularization, int N)
{
  int i;

  for (i = 1; i <= N; i++)
  {
    pi[i] = tau_it[i][1];
  }
  if (regularization == UT_TRUE)
  {
    for (i = 1; i <= N; i++)
    {
      pi[i] += omega_Q1;
    }
  }
  normalize_dvec(pi, N);

  return(UT_OK);
}

/*******************************************************************************
 UPDATE_TRANS_PROBS
*******************************************************************************/
int update_trans_probs(double **A, double **sum_t_tau_ijt, 
                       double *sum_t_tau_it, double omega_Q2, 
                       int regularization, int N)
{
  int i, j;

  if (regularization == UT_FALSE)
  {
    for (i = 1; i <= N; i++)
    {
      for (j = 1; j <= N; j++)
      {
        A[i][j] = sum_t_tau_ijt[i][j] / sum_t_tau_it[i];
      }
    }
  }
  else
  {
    for (i = 1; i <= N; i++)
    {
      for (j = 1; j <= N; j++)
      {
        A[i][j] = (sum_t_tau_ijt[i][j] + omega_Q2) /
                  (sum_t_tau_it[i] + ((double) N) * omega_Q2);
      }
    }
  }

  return(UT_OK);
}

/*******************************************************************************
 FORWARD_BACKWARD
*******************************************************************************/
int forward_backward(double *pi, double **A, double **prob_data, 
                     double **alpha, double **beta, double *scale, 
                     int N, int T)
{
  int    i, j, t;
  double sum;
  double beta_ti;
  double scale_t;
  double *alpha_t;
  double *alpha_t_minus;
  double *beta_t;
  double *beta_t_plus;
  double *prob_data_t;
  double *prob_data_t_plus;

  /* forward initialization */
  scale[1] = 0.0;

  for (i = 1; i <= N; i++)
  {
    alpha[1][i] = pi[i] * prob_data[1][i];
    scale[1] += alpha[1][i];
  }

  scale[1] = 1.0 / scale[1];

  for (i = 1; i <= N; i++)
  {
    alpha[1][i] *= scale[1];
  }

  /* calculate for t = 2,...,T */
  for (t = 2; t <= T; t++)
  {
    scale_t = 0.0;
    alpha_t = alpha[t];
    alpha_t_minus = alpha[t-1];
    prob_data_t = prob_data[t];

    for (i = 1; i <= N; i++)
    {
      sum = 0.0;
      for (j = 1; j <= N; j++)
      {
        sum += alpha_t_minus[j] * A[j][i];
      }
      alpha_t[i] = sum * prob_data_t[i];
 
      scale_t += alpha_t[i];
    }

    scale_t = 1.0 / scale_t;
    scale[t] = scale_t;

    for (i = 1; i <= N; i++)
    {
      alpha_t[i] *= scale_t;
    }
  }

  /* backwards initialization */
  for (i = 1; i <= N; i++)
  {
    beta[T][i] = 1.0;
  }

  /* calculate for t = T-1,...,1 */
  for (t = T-1; t >= 1; t--)
  {
    beta_t = beta[t];
    beta_t_plus = beta[t+1];
    prob_data_t_plus = prob_data[t+1];
    for (i = 1; i <= N; i++)
    {
      beta_ti = 0.0;
      for (j = 1; j <= N; j++)
      {
        beta_ti += A[i][j] * prob_data_t_plus[j] * beta_t_plus[j];
      }
      beta_ti *= scale[t+1];

      beta_t[i] = beta_ti;
    }
  }

  /* Less computationally efficient, but avoids (most) overflow/underflow issues
     generated by very large or very small scaling coefficients */
  for (t = 1; t <= T; t++)
  {
    beta_t = beta[t];
    for (i = 1; i <= N; i++)
    {
      beta_t[i] *= scale[t];
      if (beta_t[i] > DBL_MAX)
      {
        beta_t[i] = DBL_MAX;
      }
      else if (beta_t[i] < DBL_MIN)
      {
        beta_t[i] = DBL_MIN;
      }
    }
  }

  return(UT_OK);
}

/*******************************************************************************
 ESTIMATE_HIDDEN
*******************************************************************************/
int estimate_hidden(double **alpha, double **beta, double **A, 
                    double **prob_data, double **tau_it, double ***tau_ijt,
                    int N, int T)
{
  int    i, j, t;
  double sum, inv_sum;
  double tau_ijt_ijt;
  double alpha_ti;
  double *alpha_t;
  double *beta_t_plus;
  double *prob_data_t_plus;

  for (t = 1; t < T; t++)
  {
    alpha_t = alpha[t];
    beta_t_plus = beta[t+1];
    prob_data_t_plus = prob_data[t+1];

    /* calculate tau_ijt */
    sum = 0.0;

    for (i = 1; i <= N; i++)
    {
      alpha_ti = alpha_t[i];
      for (j = 1; j <= N; j++)
      {
        tau_ijt_ijt = alpha_ti * A[i][j] * prob_data_t_plus[j] *
                      beta_t_plus[j];
        sum += tau_ijt_ijt;
        tau_ijt[i][j][t] = tau_ijt_ijt;
      }
    }

    inv_sum = 1.0 / sum;

    for (i = 1; i <= N; i++)
    {
      for (j = 1; j <= N; j++)
      {
        tau_ijt[i][j][t] *= inv_sum;
      }
    }

    /* calculate tau_it */
    sum = 0.0;

    for (i = 1; i <= N; i++)
    {
      tau_it[i][t] = alpha[t][i] * beta[t][i];
      sum += tau_it[i][t];
    }

    for (i = 1; i <= N; i++)
    {
      tau_it[i][t] /= sum;
    }
  }

  /* finish tau_it for t = T */
  sum = 0.0;

  for (i = 1; i <= N; i++)
  {
    tau_it[i][T] = alpha[T][i] * beta[T][i];
    sum += tau_it[i][T];
  }

  for (i = 1; i <= N; i++)
  {
    tau_it[i][T] /= sum;
  }
 
  return(UT_OK);
}

/*******************************************************************************
 VITERBI_STATE_ASSIGNMENT
*******************************************************************************/
int viterbi_state_assignment(double **prob_data, double *pi, double **A, 
                            int *Q, int N, int T)
{
  int    i, j, t;
  double **delta;
  int    **psi;
  double *delta_t;
  double *delta_t_minus;
  int    *psi_t;
  double prob_data_tj;
  double max_delta;
  double max_delta_minus_a;
  double delta_minus_a;
  double delta_tj;
  double scale_t;

  /* compute the optimal state sequence using the Viterbi algorithm */
  delta = NR_dmatrix(1, T, 1, N);
  psi = NR_imatrix(1, T, 1, N);

  scale_t = 0.0;

  /* Viterbi initialization */
  for (i = 1; i <= N; i++)
  {
    delta[1][i] = pi[i] * prob_data[1][i];
    psi[1][i] = 0;
    scale_t += delta[1][i];
  }

  /* scale to avoid underflow/overflow issues */
  for (i = 1; i <= N; i++)
  {
    delta[1][i] /= scale_t;    
  }

  /* Viterbi forward calculation */
  for (t = 2; t <= T; t++)
  {
    delta_t = delta[t];
    delta_t_minus = delta[t-1];
    psi_t = psi[t];
    scale_t = 0.0;

    for (j = 1; j <= N; j++)
    {
      max_delta = 0.0;
      max_delta_minus_a = 0.0;
      prob_data_tj = prob_data[t][j];

      for (i = 1; i <= N; i++)
      {
        delta_minus_a = delta_t_minus[i] * A[i][j];
        if (delta_minus_a >= max_delta_minus_a)
        {
          psi_t[j] = i;
          max_delta_minus_a = delta_minus_a;
        }
        delta_tj = delta_minus_a * prob_data_tj;
        if (delta_tj >= max_delta)
        {
          delta_t[j] = delta_tj;
          max_delta = delta_tj;
        }
      }
      scale_t += delta_t[j];
    }

    /* scale to avoid underflow/overflow issues */
    for (j = 1; j <= N; j++)
    {
      delta_t[j] /= scale_t;
    }
  }
 
  /* determine final state */
  max_delta = 0.0;
  for (i = 1; i <= N; i++)
  {
    if (delta[T][i] >= max_delta)
    {
      Q[T] = i;
      max_delta = delta[T][i];
    }
  }

  /* backwards part of Viterbi */
  for (t = T-1; t >= 1; t--)
  {
    Q[t] = psi[t+1][Q[t+1]];
  }

  NR_free_dmatrix(delta, 1, T, 1, N);
  NR_free_imatrix(psi, 1, T, 1, N);

  return(UT_OK);
}
