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
static char rcsid[] = "$Id: rdahmm_gauss_output.c,v 1.9 2005/07/20 21:07:19 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_gauss_output.c,v $
 * Revision 1.9  2005/07/20 21:07:19  granat
 * Added function for calculating state sequence given the model
 * This function is eval_gauss_output_hmm()
 *
 * Revision 1.8  2003/08/28 17:05:00  granat
 * more revisions for faster computation; little result
 *
 * Revision 1.7  2003/08/28 14:02:07  granat
 * minor changes to get some computational speedup
 *
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
 READ_GAUSS_OUTPUT_PARAMS
*******************************************************************************/
int read_gauss_output_params(char *B_file, double **mu, double ***sigma,
                             int N, int D)
{
  int i, d, d2;
  double **vals;

  read_dmatrix(B_file, N * (D + 1), D, &vals);

  for (i = 1; i <= N; i++)
  {
    for (d2 = 1; d2 <= D; d2++)
    {
      mu[i][d2] = vals[(i-1)*(D+1)+1][d2];
      for (d = 1; d <= D; d++)
      {
        sigma[i][d][d2] = vals[(i-1)*(D+1)+1+d][d2];
      }
    }
  }

  NR_free_dmatrix(vals, 1, N * (D + 1), 1, D);

  return(UT_OK);
}

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
int init_gauss_output_hmm(int init_type, double **cont_data, double *pi, 
                          double **A, double **mu, double ***sigma, int N, 
                          int D, int T)
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
  else if (init_type == RDAHMM_KMEANS_INIT)
  {
    kmeans_init_hmm(cont_data, pi, A, mu, sigma, N, D, T);    
  }
    
  return(UT_OK);
}
 

/*******************************************************************************
 KMEANS_INIT_HMM
*******************************************************************************/
int kmeans_init_hmm(double **cont_data, double *pi, double **A, double **mu,
                    double ***sigma, int N, int D, int T)
{
  int    i;
  int    k;
  int    n;
  int    num_tries;
  int    num_iters;
  int    num_empty_clust;
  double log_likelihood;
  double best_log_likelihood;
  double **init_means;
  double **kmeans_mu;
  double ***kmeans_sigma;
  double *kmeans_weights;
  double **probs;
  double *prob_data;
  double *sum_probs;
  double min_diag;
  int    *clust_size;
  int    *clust_label;

  /* Allocate working memory */
  init_means = NR_dmatrix(1, N, 1, D);
  kmeans_mu = NR_dmatrix(1, N, 1, D);
  kmeans_sigma = set_of_dmatrices(N, D, D);
  kmeans_weights = NR_dvector(1, N);
  probs = NR_dmatrix(1, N, 1, T);
  prob_data = NR_dvector(1, T);
  sum_probs = NR_dvector(1, N);
  clust_size = NR_ivector(1, N);
  clust_label = NR_ivector(1, T);

  /* For now these parameters are hardcoded */
  num_tries = 10;
  min_diag = 1.0e-6;

  best_log_likelihood = DBL_MIN;
  for (i = 1; i <= num_tries; i++)
  {
    /* choose K random points in the data as starting means for this run */
    random_means_from_data(init_means, N, D, cont_data, T);

    /* perform k-means */
    k_means(cont_data, T, D, kmeans_mu, init_means, clust_size,
            clust_label, N, &num_iters, &num_empty_clust, UT_FALSE);

    /* compute covariance matrices and mixture weights based on means and
       cluster memberships */
    k_means_mixture_estimate(cont_data, kmeans_mu, kmeans_sigma, 
                             kmeans_weights, clust_label, T, D, N);

    /* given the parameters, compute the likelihood of this model */
    log_likelihood = mixture_likelihood(N, kmeans_mu, kmeans_sigma, 
                                        kmeans_weights, cont_data, probs, 
                                        prob_data, sum_probs,
                                        T, D, "rows", min_diag);

    /* save the current best set of parameters */
    if (log_likelihood > best_log_likelihood)
    {
      best_log_likelihood = log_likelihood;
      for (k = 1; k <= N; k++)
      {
        copy_dvec(kmeans_mu[k], mu[k], D);
        copy_dmat(kmeans_sigma[k], sigma[k], D, D);
        copy_dvec(kmeans_weights, pi, N);
        for (n = 1; n <= N; n++)
        {
          copy_dvec(kmeans_weights, A[n], N);
        }
      }
    }
  }  /* end k-means loop over numtries */

  NR_free_dmatrix(init_means, 1, N, 1, D);
  NR_free_dmatrix(kmeans_mu, 1, N, 1, D);
  free_set_of_dmatrices(kmeans_sigma, N, D, D);
  NR_free_dvector(kmeans_weights, 1, N);
  NR_free_dmatrix(probs, 1, N, 1, T);
  NR_free_dvector(prob_data, 1, T);
  NR_free_dvector(sum_probs, 1, N);
  NR_free_ivector(clust_size, 1, N);
  NR_free_ivector(clust_label, 1, T);

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
  double data_td;
  double *mu_i;
  double *diff_i;
  double *inv_sqrt_det_sigma;
  double *prob_data_t;

  /* allocate working memory */
  inv_sqrt_det_sigma = NR_dvector(1, N);
  diff = NR_dmatrix(1, N, 1, D);
  temp = NR_dvector(1, D);

  /* pre-calculate constant */
  two_pi = 2 * M_PI;
  c = pow(two_pi, (double) D);
  c = sqrt(1.0 / c);

  /* pre-calculate determinant inverses */
  for (i = 1; i <= N; i++)
  {
    inv_sqrt_det_sigma[i] = 1.0 / sqrt_det_sigma[i];
  }

  /* calculate data probabilities */
  for (t = 1; t <= T; t++)
  {
    data_t = data[t];
    prob_data_t = prob_data[t];

    for (i = 1; i <= N; i++)
    {
      mu_i = mu[i];
      diff_i = diff[i];
      diff_dvec(data_t, mu_i, diff_i, D);
    }

    for (i = N; i >= 1; --i)
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

      /* inner product of temp vectors */
      dot = dot_dvec(temp, temp, D);

      /* Mahalanobis distance */
      prob_data_t[i] = c * exp(-0.5 * dot) * inv_sqrt_det_sigma[i];

      /* To prevent numerical problems make sure probabilities are > 0 */
      if (prob_data_t[i] < DBL_MIN)
      {
        prob_data_t[i] = DBL_MIN;
      }
    }
  }

  /* free working memory */
  NR_free_dvector(inv_sqrt_det_sigma, 1, N);
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
  int    d_plus, d_minus;
  int    d2_plus, d2_minus;
  double tau_it_it;
  double outer_product;
  double *diff;
  double diff_d, diff_d2;
  double *mu_i, **sigma_i, *sigma_id;
  double *data_t;
  double tau_diff;
  double **rank1;
  double *rank1_d;

  diff = NR_dvector(1, D);

  for (i = 1; i <= N; i++)
  {
    fast_zero_dmat(sigma[i], D, D);

    mu_i = mu[i];
    sigma_i = sigma[i];
    for (t = 1; t <= T; t++)
    {
      data_t = data[t];
      tau_it_it = tau_it[i][t];

      diff_dvec(data_t, mu_i, diff, D);

      sigma_i[1][1] += tau_it_it * NR_dsqr(diff[1]);
      for (d = 2; d <= D; d++)
      {
        diff_d = diff[d];
        sigma_id = sigma_i[d];
        tau_diff = tau_it_it * diff_d;
        for (d2 = 1; d2 < d; d2++)
        {
          sigma_id[d2] += tau_diff * diff[d2];
        }
        sigma_id[d] += tau_diff * diff_d;
      }
    }
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

  return(UT_OK);
}
  

/*******************************************************************************
 UPDATE_BLOCK_GAUSS_COVARS
*******************************************************************************/
int update_block_gauss_covars(double **data, double **tau_it, 
                              double *sum_T_tau_it, double **mu, 
                              double ***sigma, int *blocks, 
                              int N, int T, int D)
{
  int    i, t, d, d2;
  int    d_plus, d_minus;
  int    d2_plus, d2_minus;
  double tau_it_it;
  double outer_product;
  double *diff;
  double diff_d, diff_d2;
  double *mu_i, **sigma_i, *sigma_id;
  double *data_t;
  double tau_diff;
  double **rank1;
  double *rank1_d;

  diff = NR_dvector(1, D);

  for (i = 1; i <= N; i++)
  {
    fast_zero_dmat(sigma[i], D, D);

    mu_i = mu[i];
    sigma_i = sigma[i];
    for (t = 1; t <= T; t++)
    {
      data_t = data[t];
      tau_it_it = tau_it[i][t];

      diff_dvec(data_t, mu_i, diff, D);

      sigma_i[1][1] += tau_it_it * NR_dsqr(diff[1]);

      for (d = 2; d <= D; d++)
      {
        diff_d = diff[d];
        sigma_id = sigma_i[d];
        tau_diff = tau_it_it * diff_d;
        for (d2 = blocks[d]; d2 < d; d2++)
        {
          sigma_id[d2] += tau_diff * diff[d2];
        }
        sigma_id[d] += tau_diff * diff_d;
      }
    }
    for (d = 1; d <= D; d++)
    {
      for (d2 = blocks[d]; d2 <= d; d2++)
      {
        sigma[i][d][d2] /= sum_T_tau_it[i];
        sigma[i][d2][d] = sigma[i][d][d2];
      }
    }
  }

  NR_free_dvector(diff, 1, D);
 
  return(UT_OK);
}

/*******************************************************************************
 UPDATE_GRAPH_GAUSS_COVARS
*******************************************************************************/
int update_graph_gauss_covars(double **data, double **tau_it, 
                              double *sum_T_tau_it, double **mu, 
                              double ***sigma, int **graph, 
                              int N, int T, int D)
{
  int    i, t, d, d2;
  int    d_plus, d_minus;
  int    d2_plus, d2_minus;
  double tau_it_it;
  double outer_product;
  double *diff;
  double diff_d, diff_d2;
  double *mu_i, **sigma_i, *sigma_id;
  double *data_t;
  double tau_diff;
  double **rank1;
  double *rank1_d;

  diff = NR_dvector(1, D);

  for (i = 1; i <= N; i++)
  {
    fast_zero_dmat(sigma[i], D, D);

    mu_i = mu[i];
    sigma_i = sigma[i];
    for (t = 1; t <= T; t++)
    {
      data_t = data[t];
      tau_it_it = tau_it[i][t];

      diff_dvec(data_t, mu_i, diff, D);

      sigma_i[1][1] += tau_it_it * NR_dsqr(diff[1]);
      for (d = 2; d <= D; d++)
      {
        diff_d = diff[d];
        sigma_id = sigma_i[d];
        tau_diff = tau_it_it * diff_d;
        for (d2 = 1; d2 < d; d2++)
        {
          if (graph[d][d2] == 1)
          {
            sigma_id[d2] += tau_diff * diff[d2];
          }
        }
        sigma_id[d] += tau_diff * diff_d;
      }
    }
    for (d = 1; d <= D; d++)
    {
      for (d2 = 1; d2 <= d; d2++)
      {
        if (graph[d][d2] == 1)
        {
          sigma[i][d][d2] /= sum_T_tau_it[i];
        }
        sigma[i][d2][d] = sigma[i][d][d2];
      }
    }
  }

  NR_free_dvector(diff, 1, D);

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
                        double **mu, double ***sigma, double *omega_Q3,
                        double omega_sigma, int constrained, int block,
                        int *blocks, int **graph, int N, int T, int D, 
                        int verbose)
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

  /* check to make sure that omega_Q3 obeys concavity condition */
  for (i = 1; i <= N; i++)
  {
    norm_sigma = spec_norm_covar_dmat(sigma[i], D);
    if (*omega_Q3 > 1 / (2 * N * norm_sigma))
    {
      if (verbose)
      {
        fflush(stdout);
        fprintf(stdout, "Corrected omega_Q3 for concavity condition\n");
      }
      *omega_Q3 = 1 / (2 * N * norm_sigma);
    }
  }

  /* get initial approximations */
  update_gauss_means(data, tau_it, sum_T_tau_it, M, N, T, D);

  if (constrained == UT_FALSE)
  {
    update_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, N, T, D);
  }
  else
  {
    if (block == UT_TRUE)
    {
      update_block_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, blocks,
                                N, T, D);
    }
    else
    {
      update_graph_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, graph,
                                N, T, D);
    }
  }

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
        left_side[(i-1)*D+d][(j-1)*D+d] += *omega_Q3;
      }
    }
    for (d = 1; d <= D; d++)
    {
      left_side[(i-1)*D+d][(i-1)*D+d] -= ((double) N) * (*omega_Q3);
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
  if (constrained == UT_FALSE)
  {
    update_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, N, T, D);
  }
  else
  {
    if (block == UT_TRUE)
    {
      update_block_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, blocks,
                                N, T, D);
    }
    else
    {
      update_graph_gauss_covars(data, tau_it, sum_T_tau_it, M, sigma, graph,
                                N, T, D);
    }
  }

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
                           double anneal_factor, double beta_min, 
                           double peps, double thresh, int *covgraph,
                           int maxiters,
                           double *L, int *Q, double *pi, double **A, 
                           double **mu, double ***sigma, int N, int T, int D,
                           int viterbi, int verbose)
{
  int    i, j, t, d, d2;
  int    iters;
  int    converged;
  int    constrained;
  int    block;
  int    *blocks;
  int    **graph;
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
    comp_temp = beta_min;
  }
  else
  {
    comp_temp = 1.0;
  }

  /* see if the covariance matrix is being constrained by a connectivity
    graph */
  if (covgraph[1] > 0)
  {
    constrained = UT_TRUE;

    /* is it a block diagonal constraint? */
    block = UT_TRUE;
    for (d = 2; d <= D; d++)
    {
      if (covgraph[d] < covgraph[d-1])
      {
        block = UT_FALSE;
        break;
      }
    }

    if (block == UT_TRUE)
    {
      /* infer the first nonzero element of each row */
      blocks = NR_ivector(1, D);
      blocks[1] = 1;
      for (d = 2; d <= D; d++)
      {
        if (covgraph[d] == covgraph[d-1])
        {
          blocks[d] = blocks[d-1];
        }
        else
        {
          blocks[d] = d;
        }
      }
    }
    else
    {
      /* infer connection graph matrix */
      graph = NR_imatrix(1, D, 1, D);
      for (d = 1; d <= D; d++)
      {
        graph[d][d] = 1;
        for (d2 = 1; d2 < d; d2++)
        {
          if (covgraph[d] == covgraph[d2])
          {
            graph[d][d2] = 1;
            graph[d2][d] = 1;
          }
        }
      }
    }
  }
  else
  {
    constrained = UT_FALSE;
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
          if (euclid_dist_dvec(mu[i], mu[j], D) < 1.0e-3 * peps)
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
      /* test for NaNs entering the computation, looking at the likelihood
         is the easiest way */
      if (isnan(*L))
      {
        fprintf(stderr, "NaN discovered in computation, aborting.\n");
        return(UT_ERROR);
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
        if (Qfunc < (last_Qfunc - thresh))
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
            if (iters > maxiters)
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
        if (constrained == UT_FALSE)
        {
          update_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, N, T, D);
        }
        else
        {
          if (block == UT_TRUE)
          {
            update_block_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, 
                                      blocks, N, T, D);
          }
          else
          {
            update_graph_gauss_covars(data, tau_it, sum_T_tau_it, mu, sigma, 
                                      graph, N, T, D);
          }
        }
      }
      else
      {
        omega_Q3 = omega_Q3_orig;
        if (reg_type == RDAHMM_EUCLID_REG)
        {
          update_gauss_euclid(data, tau_it, sum_T_tau_it, mu, sigma, &omega_Q3,
                              omega_sigma, constrained, block, blocks, graph,
                              N, T, D, verbose);
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

    if (comp_temp == 1.0)
    {
      break;
    }

    if ((anneal_factor != 0) && (comp_temp * anneal_factor > comp_temp + anneal_step))
    {
      comp_temp *= anneal_factor;
    }
    else
    {
      comp_temp += anneal_step;
    }
   
    if (comp_temp >= 1.0)
    {
      comp_temp = 1.0;
    }
  } /* end loop over temperatures */

  if (viterbi == UT_FALSE)
  {
    /* compute individually most likely state assignemnts */
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
  }
  else
  {
    viterbi_state_assignment(prob_data, pi, A, Q, N, T);
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

  /* free optional working memory */
  if (constrained == UT_TRUE)
  {
    if (block == UT_TRUE)
    {
      NR_free_ivector(blocks, 1, D);
    }
    else
    {
      NR_free_imatrix(graph, 1, D, 1, D);
    }
  }

  return(UT_OK);
}

/*******************************************************************************
 EVAL_GAUSS_OUTPUT_HMM
*******************************************************************************/
int eval_gauss_output_hmm(double **data, double *L, int *Q, double *pi, 
                          double **A, double **mu, double ***sigma, int N, 
                          int T, int D, int viterbi, int verbose)
{
  int    i, j, t, d;
  double norm_sigma;
  double ***inv_chol_sigma;
  double *sqrt_det_sigma;
  double **prob_data;
  double **alpha, **beta, *scale;
  double **tau_it, ***tau_ijt;
  double max_tau;

  /* allocate working memory */
  inv_chol_sigma = set_of_dmatrices(N, D, D);
  sqrt_det_sigma = NR_dvector(1, N);
  prob_data = NR_dmatrix(1, T, 1, N);
  alpha = NR_dmatrix(1, T, 1, N);
  beta = NR_dmatrix(1, T, 1, N);
  scale = NR_dvector(1, T);
  tau_it = NR_dmatrix(1, N, 1, T);
  tau_ijt = set_of_dmatrices(N, N, T-1);

  /* calculate inverses of the covariance matrices */
  for (i = 1; i <= N; i++)
  {
    inv_chol_covar_dmat(sigma[i], inv_chol_sigma[i], &sqrt_det_sigma[i], D);
  }

  /* calculate P(O_t|B) for t = 1,...,T */
  prob_data_gauss(data, mu, inv_chol_sigma, sqrt_det_sigma, prob_data, 
                  N, T, D);

  /* calculate forward-backward parameters alpha and beta */
  forward_backward(pi, A, prob_data, alpha, beta, scale, N, T);

  if (viterbi == UT_FALSE)
  {
    /* calculate estimates of the hidden variable probabilities */
    estimate_hidden(alpha, beta, A, prob_data, tau_it, tau_ijt, N, T);

    /* compute individually most likely state assignemnts */
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
  }
  else
  {
    viterbi_state_assignment(prob_data, pi, A, Q, N, T);
  }

  /* compute the log likelihood of the model given the data */
  *L = calc_log_likelihood(scale, T);

  /* free working memory */
  free_set_of_dmatrices(inv_chol_sigma, N, D, D);
  NR_free_dvector(sqrt_det_sigma, 1, N);
  NR_free_dmatrix(prob_data, 1, T, 1, N);
  NR_free_dmatrix(alpha, 1, T, 1, N);
  NR_free_dmatrix(beta, 1, T, 1, N);
  NR_free_dvector(scale, 1, T);
  NR_free_dmatrix(tau_it, 1, N, 1, T);
  free_set_of_dmatrices(tau_ijt, N, N, T-1);

  return(UT_OK);
}
