/*******************************************************************************
MODULE NAME
rdahmm_gauss_output

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
static char rcsid[] = "$Id: rdahmm_gauss_output.c,v 1.6 2003/08/27 20:23:14 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_gauss_output.c,v $
 * Revision 1.6  2003/08/27 20:23:14  granat
 * reorganized operations in prob_data_gauss to improve efficiency
 *
 * Revision 1.5  2003/08/27 18:03:17  granat
 * various optimization methods tried, unrolling rejected as not cost effective
 *
 * Revision 1.4  2003/08/25 20:01:35  granat
 * eliminated (now commented) earlier attempt at unrolling
 *
 * Revision 1.3  2003/08/25 18:54:55  granat
 * some edits to change instruction ordering in inner loop
 *
 * Revision 1.2  2003/08/25 16:05:33  granat
 * unrolled inner loop of update_gauss_covars by a factor of two, made other
 * minor edits for computational efficiency as guided by Shark
 *
 * Revision 1.1  2003/08/24 21:51:57  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* BLAS library */
#include "cblas.h"

/* CP library */
#include "cp.h"

/* UT library */
#include "ut.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* local header files */
#include "rdahmm.h"
#include "rdahmm_generic.h"
#include "rdahmm_stat_utils.h"

/* this program's header */
#include "rdahmm_gauss_output.h"


/*******************************************************************************
 WRITE_GAUSS_OUTPUT_PARAMS
*******************************************************************************/
int write_gauss_output_params(char *B_file, double **mu, double ***sigma,
                              int N, int D)
{
  int i;

  for (i = 1; i <= N; i++)
  {
    if (i == 1)
    {
      write_drow(B_file, D, mu[i], "w");
    }
    else
    {
      write_drow(B_file, D, mu[i], "a");
    }
    write_dmatrix(B_file, D, D, sigma[i], "a");
  }

  return(UT_OK);
}

/*******************************************************************************
 INIT_GAUSS_OUTPUT_HMM
*******************************************************************************/
int init_gauss_output_hmm(int init_type, double *pi, double **A, double **mu,
                          double ***sigma, int N, int D)
{
  int i;

  if (init_type == RDAHMM_RANDOM_INIT)
  {
    /* initial and state transition probabilities */
    random_dvector(pi, N);
    normalize_dvec(pi, N);

    random_dmatrix(A, N, N);
    normalize_dmat_rows(A, N, N);

    /* output distribution */
    random_dmatrix(mu, N, D);
    for (i = 1; i <= N; i++)
    {
      generate_rand_covar_dmat(sigma[i], D);
    }
  }
    
  return(UT_OK);
}
 

/*******************************************************************************
 CALC_REGULARIZED_Q_FUNCTION_GAUSS
*******************************************************************************/
double calc_regularized_Q_function_gauss(double Qfunc, double *pi, double **A, 
                                         double **mu, double ***inv_chol_sigma,
                                         double *sum_T_tau_it, 
                                         double omega_Q1, double omega_Q2, 
                                         double omega_Q3, double omega_sigma, 
                                         int reg_type, int N, int D)
{
  int    i, j, d;
  double reg_Qfunc;

  reg_Qfunc = Qfunc;
  
  /* initial probabilities */
  if (omega_Q1 > 0.0)
  {
    for (i = 1; i <= N; i++)
    {
      reg_Qfunc += omega_Q1 * log(pi[i] + DBL_MIN);
    }
  }

  /* transition probabilities */
  if (omega_Q2 > 0.0)
  {
    for (i = 1; i <= N; i++)
    {
      for (j = 1; j <= N; j++)
      reg_Qfunc += omega_Q2 * log(A[i][j] + DBL_MIN);
    }
  }

  /* output distributions */
  if (omega_Q3 > 0.0)
  {
    if (reg_type == RDAHMM_EUCLID_REG)
    {
      for (i = 1; i <= N; i++)
      {
        for (j = 1; j <= N; j++)
        {
          reg_Qfunc += 0.5 * sum_T_tau_it[i] * omega_Q3 * 
                       euclid_dist_dvec(mu[i], mu[j], D);
        }
      }
    }
  }

  /* covariance matrices */
  if (omega_sigma > 0.0)
  {
    for (i = 1; i <= N; i++)
    {
      for (d = 1; d <= D; d++)
      {
        reg_Qfunc -= 0.5 * omega_sigma * NR_sqr(inv_chol_sigma[i][d][d]);
      }
    }
  }

  return(reg_Qfunc);
}

/*******************************************************************************
 PROB_DATA_GAUSS
*******************************************************************************/
int prob_data_gauss(double **data, double **mu, double ***inv_chol_sigma, 
                    double *sqrt_det_sigma, double **prob_data, int N,
                    int T, int D)
{
  int    i, t, d, d2;
  int    d2_minus;
  double two_pi;
  double c;
  double dot;
  double **diff;
  double *temp;
  double *inv_chol_sigma_id;
  double **inv_chol_sigma_i;
  double temp_d;
  double *data_t;
  double *mu_i;
  double *diff_i;

  /* allocate working memory */
  diff = NR_dmatrix(1, N, 1, D);
  temp = NR_dvector(1, D);

  /* pre-calculate constant */
  two_pi = 2 * M_PI;
  c = pow(two_pi, (double) D);
  c = sqrt(1.0 / c);

  /* calculate data probabilities */
  for (t = 1; t <= T; t++)
  {
    data_t = data[t];
    for (i = 1; i <= N; i++)
    {
      mu_i = mu[i];
      diff_i = diff[i];
      diff_dvec(data_t, mu_i, diff_i, D);
    }

    for (i = 1; i <= N; i++)
    {
      diff_i = diff[i];
      inv_chol_sigma_i = inv_chol_sigma[i];
      for (d = D; d >= 1; --d)
      {
        inv_chol_sigma_id = inv_chol_sigma_i[d];
        temp_d = 0.0;
        for (d2 = d; d2 >= 1; --d2)
        {
          temp_d += diff_i[d2] * inv_chol_sigma_id[d2];   
        }
        temp[d] = temp_d;
      }
/*
      copy_dvec(diff_i, temp, D); 
      cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, D,
                  &inv_chol_sigma_i[1][1], D, &temp[1], 1);
*/
                  
      /* inner product of temp vectors */
      dot = dot_dvec(temp, temp, D);

      /* Mahalanobis distance */
      prob_data[t][i] = c * exp(-0.5 * dot) / sqrt_det_sigma[i];

      /* This doesn't fix the problem but let's try it: */
      if (prob_data[t][i] < DBL_MIN)
      {
        prob_data[t][i] = DBL_MIN;
      }
    }
  }

  /* free working memory */
  NR_free_dmatrix(diff, 1, N, 1, D);
  NR_free_dvector(temp, 1, D);

  return(UT_OK);
}

/*******************************************************************************
 UPDATE_GAUSS_MEANS
*******************************************************************************/
int update_gauss_means(double **data, double **tau_it, double *sum_T_tau_it,
                       double **mu, int N, int T, int D)
{
  int i, t, d;
  double *mu_i;
  double *data_t;
  double tau_it_it;

  fast_zero_dmat(mu, N, D);

  for (i = 1; i <= N; i++)
  {
    mu_i = mu[i];
    for (t = 1; t <= T; t++)
    {
      tau_it_it = tau_it[i][t];
      data_t = data[t];
      for (d = 1; d <= D; d++)
      {
        mu_i[d] += tau_it_it * data_t[d];
      }
    }
    for (d = 1; d <= D; d++)
    {
      mu_i[d] /= sum_T_tau_it[i];
    }
  }

  return(UT_OK);
}

/*******************************************************************************
 UPDATE_GAUSS_COVARS
*******************************************************************************/
int update_gauss_covars(double **data, double **tau_it, double *sum_T_tau_it, 
                        double **mu, double ***sigma, int N, int T, int D)
{
  int    i, t, d, d2;
  int    d_plus;
  int    d2_plus;
  double tau_it_it;
  double *tau_it_i;
  double outer_product;
  double *diff;
/*
  double **diff;
*/
  double diff_d, diff_d2;
  double diff_id;
  double *diff_i;
  double *mu_i, **sigma_i, *sigma_id;
  double *data_t;
  double tau_diff;
  double tau_it_it_times_diff_id;
  double temp;

  diff = NR_dvector(1, D);
/*
  diff = NR_dmatrix(1, N, 1, D);
*/

  for (i = 1; i <= N; i++)
  {
    fast_zero_dmat(sigma[i], D, D);
  }

  for (t = 1; t <= T; t++)
  {
    data_t = data[t];
    for (i = 1; i <= N; i++)
    {
      mu_i = mu[i];
      sigma_i = sigma[i];
      tau_it_it = tau_it[i][t];

      diff_dvec(data_t, mu_i, diff, D);

      for (d = 1; d <= D; d++)
      {
        diff_d = diff[d];
        sigma_id = sigma_i[d];
        sigma_id[d] += tau_it_it * diff_d * diff_d;

        for (d2 = d-1; d2 >= 1; --d2)
        {
          sigma_id[d2] += tau_it_it * diff_d * diff[d2];
        }
      }
      /* unrolled secondmost inner loop */
/*
      if (D % 2 == 0)
      {
        for (d = 1; d <= D; d+=2)
        {
          d_plus = d + 1;

          diff_d = diff[d];
          sigma_id = sigma_i[d];
          sigma_id[d] += tau_it_it * diff_d * diff_d;

          update_gauss_covars_inner(sigma_id, diff, tau_it_it, diff_d, d);

          diff_d = diff[d_plus];
          sigma_id = sigma_i[d_plus];
          sigma_id[d_plus] += tau_it_it * diff_d * diff_d;

          update_gauss_covars_inner(sigma_id, diff, tau_it_it, diff_d, d_plus);
        }
      }
      else
      {
        diff_d = diff[1];
        sigma_id = sigma_i[1];
        sigma_id[1] += tau_it_it * diff_d * diff_d;

        for (d = 2; d <= D; d+=2)
        {
          d_plus = d + 1;

          diff_d = diff[d];
          sigma_id = sigma_i[d];
          sigma_id[d] += tau_it_it * diff_d * diff_d;

          update_gauss_covars_inner(sigma_id, diff, tau_it_it, diff_d, d);

          diff_d = diff[d_plus];
          sigma_id = sigma_i[d_plus];
          sigma_id[d_plus] += tau_it_it * diff_d * diff_d;

          update_gauss_covars_inner(sigma_id, diff, tau_it_it, diff_d, d_plus);
        }
      } 
*/
    }
  }

  for (i = 1; i <= N; i++)
  {
    for (d = 1; d <= D; d++)
    {
      for (d2 = 1; d2 <= d; d2++)
      {
        sigma[i][d][d2] /= sum_T_tau_it[i];
        sigma[i][d2][d] = sigma[i][d][d2];
      }
    }
  }

  NR_free_dvector(diff, 1, D);
/*
  NR_free_dmatrix(diff, 1, D, 1, D);
*/

  return(UT_OK);
}
  

/*******************************************************************************
 UPDATE_GAUSS_COVARS_INNER
*******************************************************************************/
int update_gauss_covars_inner(double *sigma_id, double *diff, double tau_it_it, 
                              double diff_d, int d)
{
  int    d2, d2_minus;

  for (d2 = d-1; d2 >= 1; --d2)
  {
    sigma_id[d2] += tau_it_it * diff_d * diff[d2];
  }

/*
  if (d % 2 == 0)
  {
    sigma_id[1] += tau_it_it * diff_d * diff[1];
    for (d2 = d-1; d2 >= 2; d2-=2)
    {
      d2_minus = d2 - 1;
      sigma_id[d2] += tau_it_it * diff_d * diff[d2];
      sigma_id[d2_minus] += tau_it_it * diff_d * diff[d2_minus];
    }
  }
  else if (d != 1)
  {
    sigma_id[1] += tau_it_it * diff_d * diff[1];
    sigma_id[2] += tau_it_it * diff_d * diff[2];

    for (d2 = d-1; d2 >= 3; d2-=2)
    {
      d2_minus = d2 - 1;
      sigma_id[d2] += tau_it_it * diff_d * diff[d2];
      sigma_id[d2_minus] += tau_it_it * diff_d * diff[d2_minus];
    }
  }
*/

  return(UT_OK);
}

/*******************************************************************************
 UPDATE_GAUSS_EUCLID
*******************************************************************************/
int update_gauss_euclid(double **data, double **tau_it, double *sum_T_tau_it, 
                        double **mu, double ***sigma, double omega_Q3,
                        double omega_sigma, int N, int T, int D, int verbose)
{
  int    i, j, d, d2;
  double norm_sigma;
  double **M;
  double ***inv_sigma;
  double **left_side, *right_side;
  double even_odd;
  int    *indx;

  /* allocate working space */
  M = NR_dmatrix(1, N, 1, D);
  inv_sigma = set_of_dmatrices(N, D, D);
  left_side = NR_dmatrix(1, N*D, 1, N*D);
  indx = NR_ivector(1, N*D);

  /* trick to save memory, use mu as working space */
  right_side = mu[1];

  /* check to make sure that omega_Q3 obayes concavity condition */
  for (i = 1; i <= N; i++)
  {
    norm_sigma = spec_norm_covar_dmat(sigma[i], D);
    if (omega_Q3 > 1 / (2 * N * norm_sigma))
    {
      if (verbose)
      {
        fflush(stdout);
        fprintf(stdout, "Corrected omega_Q3 for concavity condition\n");
      }
      omega_Q3 = 1 / (2 * N * norm_sigma);
    }
  }

  /* get initial approximations */
  update_gauss_means(data, tau_it, sum_T_tau_it, M, N, T, D);
  update_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, N, T, D);

  /* apply covariance regularization */
  if (omega_sigma > 0.0)
  {
    for (i = 1; i <= N; i++)
    {
      for (d = 1; d <= D; d++)
      {
        sigma[i][d][d] += omega_sigma;
      }
    }
  }

  /* calculate covariance inverses */
  for (i = 1; i <= N; i++)
  {
    copy_dmat(sigma[i], inv_sigma[i], D, D);
    invert_alloc_sym_dmat(sigma[i], D, inv_sigma[i]);    
  }

  /* fill in matrices that will form the linear system */
  fast_zero_dmat(left_side, N*D, N*D);
  fast_zero_dvec(right_side, N*D);

  for (i = 1; i <= N; i++)
  {
    for (d = 1; d <= D; d++)
    {
      for (d2 = 1; d2 <= D; d2++)
      {
        left_side[(i-1)*D+d][(i-1)*D+d2] += inv_sigma[i][d][d2];
      }
    }
    for (j = 1; j <= N; j++)
    {
      for (d = 1; d <= D; d++)
      {
        left_side[(i-1)*D+d][(j-1)*D+d] += omega_Q3;
      }
    }
    for (d = 1; d <= D; d++)
    {
      left_side[(i-1)*D+d][(i-1)*D+d] -= ((double) N) * omega_Q3;
    }
    for (d = 1; d <= D; d++)
    {
      for (d2 = 1; d2 <= D; d2++)
      {
        right_side[(i-1)*D+d] += inv_sigma[i][d][d2] * M[i][d2];
      }
    }
  }

  /* solve the system using LU decomposition */
  /* NOTE: makes no use of sparseness or symmetry */
  NR_dludcmp(left_side, N*D, indx, &even_odd);
  NR_dlubksb(left_side, N*D, indx, right_side);

  /* now that we have the new mu, resolve for new sigma */
  update_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, N, T, D);

  NR_free_dmatrix(M, 1, N, 1, D);
  free_set_of_dmatrices(inv_sigma, N, D, D);
  NR_free_dmatrix(left_side, 1, N*D, 1, N*D);
  NR_free_ivector(indx, 1, N*D);

  return(UT_OK);
}

/*******************************************************************************
 LEARN_GAUSS_OUTPUT_HMM
*******************************************************************************/
int learn_gauss_output_hmm(double **data, int regularize, int reg_type, 
                           double omega_Q1, double omega_Q2, double omega_Q3, 
                           double omega_sigma, int anneal, double anneal_step,
                           double peps, double thresh, double *L, int *Q, 
                           double *pi, double **A, double **mu, 
                           double ***sigma, int N, int T, int D,
                           int verbose)
{
  int    i, j, t, d;
  int    iters;
  int    converged;
  double Qfunc;
  double last_L, last_Qfunc;
  double omega_Q3_orig;
  double comp_temp;
  double norm_sigma;
  double ***inv_chol_sigma;
  double *sqrt_det_sigma;
  double **prob_data;
  double *pi_comp_temp, **A_comp_temp, **prob_data_comp_temp;
  double **alpha, **beta, *scale;
  double **tau_it, ***tau_ijt;
  double *sum_t_tau_it, *sum_T_tau_it, **sum_t_tau_ijt;
  double max_tau;

  /* allocate working memory */
  inv_chol_sigma = set_of_dmatrices(N, D, D);
  sqrt_det_sigma = NR_dvector(1, N);
  prob_data = NR_dmatrix(1, T, 1, N);
  pi_comp_temp = NR_dvector(1, N);
  A_comp_temp = NR_dmatrix(1, N, 1, N);
  prob_data_comp_temp = NR_dmatrix(1, T, 1, N);
  alpha = NR_dmatrix(1, T, 1, N);
  beta = NR_dmatrix(1, T, 1, N);
  scale = NR_dvector(1, T);
  tau_it = NR_dmatrix(1, N, 1, T);
  tau_ijt = set_of_dmatrices(N, N, T-1);
  sum_t_tau_it = NR_dvector(1, N);
  sum_T_tau_it = NR_dvector(1, N);
  sum_t_tau_ijt = NR_dmatrix(1, N, 1, N);

  /* Keep track of original omega_Q3 as we may be modifying omega_Q3 later */
  omega_Q3_orig = omega_Q3;

  /* Make sure that regularization weights are acceptable for the initial
     iteration */
  if ((regularize == UT_TRUE) && (reg_type == RDAHMM_EUCLID_REG))
  {
    for (i = 1; i <= N; i++)
    {
      norm_sigma = spec_norm_covar_dmat(sigma[i], D);
      if (omega_Q3 > 1 / (2 * N * norm_sigma)) 
      {
        if (verbose)
        {
          fflush(stdout);
          fprintf(stdout, "Corrected omega_Q3 for concavity condition\n");
        }
        omega_Q3 = 1 / (2 * N * norm_sigma);
      }
    }
  }

  /* set initial computational temperature */
  if (anneal == UT_TRUE)
  {
    comp_temp = 0.0;
  }
  else
  {
    comp_temp = 1.0;
  }

  /* main annealing loop */
  while (comp_temp <= 1.0)
  {
    if (verbose == UT_TRUE)
    {
      fflush(stdout);
      fprintf(stdout, "Computational temperature %g:\n", comp_temp);
    }

    /* initialization */
    iters = 1;
    converged = UT_FALSE;
    *L = -DBL_MAX;
    Qfunc = -DBL_MAX;

    /* perturb identical output distributions */
    if ((anneal == UT_TRUE) && (peps > 0.0))
    {
      for (i = 1; i <= N; i++)
      {
        for (j = 1; j < i; j++)
        {
          /* NOTE: this is a unprincipled threshold, but it seems to work 
             in practice */ 
          if (euclid_dist_dvec(mu[i], mu[j], D) < 1.0e-6 * peps * peps)
          {
            if (verbose == UT_TRUE)
            {
              fflush(stdout);
              fprintf(stdout, "Perturbing identical means %d and %d\n", i, j);
            }
            for (d = 1; d <= D; d++)
            {
              mu[i][d] += peps * (0.5 - (double) gen_rand());
            }
            break;
          }
        }
      }
    }

    while (converged == UT_FALSE)
    {
      last_L = *L;
      last_Qfunc = Qfunc;

      /* begin estimation step */

      /* calculate inverses of the covariance matrices */
      for (i = 1; i <= N; i++)
      {
        inv_chol_covar_dmat(sigma[i], inv_chol_sigma[i], &sqrt_det_sigma[i], D);
      }

      /* calculate P(O_t|B) for t = 1,...,T */
      prob_data_gauss(data, mu, inv_chol_sigma, sqrt_det_sigma, prob_data, 
                      N, T, D);

      /* calculate parameters to the power of the computational temperature */
      if (anneal == UT_TRUE)
      {
        for (i = 1; i <= N; i++)
        {
          pi_comp_temp[i] = pow(pi[i], comp_temp);
          for (j = 1; j <= N; j++)
          {
            A_comp_temp[i][j] = pow(A[i][j], comp_temp);
          }
          for (t = 1; t <= T; t++)
          {
            prob_data_comp_temp[t][i] = pow(prob_data[t][i], comp_temp);
          }
        }
      }
      else
      {
        copy_dvec(pi, pi_comp_temp, N);
        copy_dmat(A, A_comp_temp, N, N);
        copy_dmat(prob_data, prob_data_comp_temp, T, N);
      }

      /* calculate forward-backward parameters alpha and beta */
      forward_backward(pi_comp_temp, A_comp_temp, prob_data_comp_temp, 
                       alpha, beta, scale, N, T);

      /* calculate estimates of the hidden variable probabilities */
      estimate_hidden(alpha, beta, A_comp_temp, prob_data_comp_temp,
                      tau_it, tau_ijt, N, T);

      /* precalculate some quantities to save operations later */
      fast_zero_dvec(sum_t_tau_it, N);
      fast_zero_dmat(sum_t_tau_ijt, N, N);

      /* sums over t = 1,...,T-1 */
      for (t = 1; t < T; t++)
      {
        for (i = 1; i <= N; i++)
        {
          sum_t_tau_it[i] += tau_it[i][t];
          for (j = 1; j <= N; j++)
          {
            sum_t_tau_ijt[i][j] += tau_ijt[i][j][t];
          }
        }
      }

      /* we also want a quantity summed over t = 1,...,T */
      for (i = 1; i <= N; i++)
      {
        sum_T_tau_it[i] = sum_t_tau_it[i] + tau_it[i][T];
      }

      /* end estimation step */

      /* calculate log likelihood */
      *L = calc_log_likelihood(scale, T);
      if (verbose == UT_TRUE)
      {
        fflush(stdout);
        fprintf(stdout, "At iteration %d, the log likelihood is %g\n", 
                iters, *L);
      }

      /* evaluate the Q-function */
      Qfunc = calc_Q_function(tau_it, tau_ijt, pi, A, prob_data, N, T);

      if (regularize == UT_TRUE)
      {
        /* add in last bit of summation for t = T*/
        Qfunc = calc_regularized_Q_function_gauss(Qfunc, pi, A, mu, 
                                                  inv_chol_sigma, 
                                                  sum_T_tau_it, omega_Q1, 
                                                  omega_Q2, omega_Q3,
                                                  omega_sigma, reg_type, N, D);
      }

      if (verbose == UT_TRUE)
      {
        fflush(stdout);
        fprintf(stdout, "At iteration %d, the Q-function is %g\n", 
                iters, Qfunc);
      }

      /* check for convergence */
      if (iters > 1)
      {
        if (Qfunc < last_Qfunc - 1.0e-1)
        {
          if (verbose == UT_TRUE)
          {
            fflush(stdout);
            fprintf(stdout, "ERROR: decreasing Q-function\n");
          }
          break;
        }
        else 
        {
          if ((Qfunc - last_Qfunc) < thresh)
          {
            converged = UT_TRUE;
            break;
          }
          else
          {
            if (iters > 1000)
            {
              if (verbose == UT_TRUE)
              {
                fflush(stdout);
                fprintf(stdout, "ERROR: max iterations exceeded\n");
              }
              break;
            }
          }
        }
      }
      
      /* begin maximization step */ 

      /* initial state probabilities */
      update_init_probs(pi, tau_it, omega_Q1, regularize, N);

      /* transition probabilities */
      update_trans_probs(A, sum_t_tau_ijt, sum_t_tau_it, omega_Q2, 
                         regularize, N);

      /* output distributions */
      if (regularize == UT_FALSE) {
        update_gauss_means(data, tau_it, sum_T_tau_it, mu, N, T, D);
        update_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, N, T, D);
      }
      else
      {
        omega_Q3 = omega_Q3_orig;
        if (reg_type == RDAHMM_EUCLID_REG)
        {
          update_gauss_euclid(data, tau_it, sum_T_tau_it, mu, sigma, omega_Q3,
                              omega_sigma, N, T, D, verbose);
        }
      }

      /* apply covariance regularization */
      if ((omega_sigma > 0.0) || (regularize == UT_TRUE))
      {
        for (i = 1; i <= N; i++)
        {
          for (d = 1; d <= D; d++)
          {
            sigma[i][d][d] += omega_sigma;
          }
        }
      }


      /* end maximization step */

      iters++;
    } /* end inner hmm training loop for given temperature */

    /* account for small floating point variations */
    if ((comp_temp + anneal_step <= 1.0) || (comp_temp == 1.0))
    {
      comp_temp += anneal_step;
    }
    else
    {
      comp_temp = 1.0;
    }
  } /* end loop over temperatures */

  /* now compute individually most likely state assignemnts */
  for (t = 1; t <= T; t++)
  {
    max_tau = 0.0;
    for (i = 1; i <= N; i++)
    {
      if (tau_it[i][t] > max_tau)
      {
        max_tau = tau_it[i][t];
        Q[t] = i;
      }
    }
  }

  /* free working memory */
  free_set_of_dmatrices(inv_chol_sigma, N, D, D);
  NR_free_dvector(sqrt_det_sigma, 1, N);
  NR_free_dmatrix(prob_data, 1, T, 1, N);
  NR_free_dvector(pi_comp_temp, 1, N);
  NR_free_dmatrix(A_comp_temp, 1, N, 1, N);
  NR_free_dmatrix(prob_data_comp_temp, 1, T, 1, N);
  NR_free_dmatrix(alpha, 1, T, 1, N);
  NR_free_dmatrix(beta, 1, T, 1, N);
  NR_free_dvector(scale, 1, T);
  NR_free_dmatrix(tau_it, 1, N, 1, T);
  free_set_of_dmatrices(tau_ijt, N, N, T-1);
  NR_free_dvector(sum_t_tau_it, 1, N);
  NR_free_dvector(sum_T_tau_it, 1, N);
  NR_free_dmatrix(sum_t_tau_ijt, 1, N, 1, N);

  return(UT_OK);
}
