/*******************************************************************************
MODULE NAME
da_timeseg

ONE-LINE SYNOPSIS
Functions related to time series segmentation, esp. using hidden Markov models.

SCOPE OF THIS MODULE
Any functions relating directly to time series segmentation should fall into 
this module.  Functions that are just as easily applicable to static clustering
should go in da_cluster.

SEE ALSO
da_cluster is highly related and follows the same general style of organization
of functionality, since HMM induction mirrors mixture model induction in a
natural way.

REFERENCE(S)
1. "A Tutorial on Hidden Markov Models and Selected Applications in Speech 
Recognition", L. R. Rabiner, 1989, Proc. of the IEEE, 77(2):257-286.  Notation
is based on this tutorial; some references to it are made in the code for
explanatory purposes.
2. "An Introduction to Hidden Markov Models", L. R. Rabiner and B. H. Juang,
IEEE ASSP Magazine, Jan. 1986.  A shorter version of the above.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_timeseg.c,v 1.10 2000/03/31 01:47:32 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_timeseg.c,v $
 * Revision 1.10  2000/03/31 01:47:32  granat
 * fixed typo bug
 *
 * Revision 1.9  1999/07/20 17:34:37  granat
 * changed to double precision
 *
 * Revision 1.8  1998/06/30 17:23:42  agray
 * changed includes.
 *
 * Revision 1.7  1998/05/08 00:40:51  granat
 * fixed function names to be consistent with name changes; added reference
 * to new da_util module, fixed old bugs
 *
 * Revision 1.6  1998/03/09 00:39:37  granat
 * fixed boolean bug
 *
 * Revision 1.5  1997/06/05 18:54:59  granat
 * edited to conform to changes in da_linalg
 *
 * Revision 1.4  1997/06/02 15:55:16  granat
 * changed to use new NR naming convention
 *
 * Revision 1.3  1997/05/13 23:49:33  agray
 * added typecasting to nulls.
 *
 * Revision 1.2  1997/01/29 21:52:50  agray
 * added init_hmm_model_type_names, cleaned up memory allocation
 * using ut_memory.
 *
 * Revision 1.1  1996/10/31 02:22:30  agray
 * Initial revision
 * 
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

/* UT library */
#include "ut_types.h"
#include "ut_string.h"
#include "ut_error.h"
#include "ut_output.h"
#include "ut_memory.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_util.h"
#include "da_io.h"
#include "da_random.h"
#include "da_linalg.h"
#include "da_cluster.h"
#include "da_probstat.h"

/* this module's header */
#include "da_timeseg.h"


/*******************************************************************************
HMM
Allocate an HMM structure and its parts, depending on the HMM model type
specified.
AG
*******************************************************************************/
HMM *hmm(int model_type, int num_states, int num_symbols, int num_comps,
         int num_dims)

{
  HMM* the_hmm;

  the_hmm = (HMM*) malloc_return_if_fail( sizeof(HMM), (HMM*)NULL );

  /* initial state probabilities, pi */
  the_hmm->prob_init = NR_dvector(1, num_states);

  /* state transition matrix, A */
  the_hmm->prob_trans = NR_dmatrix(1, num_states, 1, num_states); 

  /* output symbol probabilities, B */
  switch (model_type)
  {
    case HMM_DISCRETE_MODEL:
      the_hmm->prob_symbol = NR_dmatrix(1, num_states, 1, num_symbols);

      the_hmm->means = (double***)NULL;
      the_hmm->covars = (double****)NULL;
      the_hmm->weights = (double**)NULL;
      break;

    case HMM_GAUSSIAN_MODEL:
    case HMM_MIXTURE_MODEL:
      the_hmm->prob_symbol = (double**)NULL;

      /* allocate mean vectors for the mixture for each state */
      the_hmm->means = set_of_dmatrices(num_states, num_comps, num_dims);

      /* allocate covariance matrices for the mixture for each state */
      the_hmm->covars = set_of_sets_of_dmatrices(num_states, num_comps, 
                                                num_dims, num_dims);

      /* allocate mixture weights for the mixture for each state */
      the_hmm->weights = NR_dmatrix(1, num_states, 1, num_comps);
      break;

    default:
      err_printf();
      log_printf("Illegal model type specified: %d\n", model_type);
      break;
  }

  return (the_hmm);
}

/*******************************************************************************
FREE_HMM
De-allocate an HMM structure and its parts, depending on the HMM model type
specified.
AG
*******************************************************************************/
void free_hmm(HMM *the_hmm)

{
  int model_type, num_states, num_symbols, num_comps, num_dims;

  /* grab values from hmm for convenience */
  model_type =  the_hmm->model_type;
  num_states =  the_hmm->num_states;
  num_symbols = the_hmm->num_symbols;
  num_comps =   the_hmm->num_comps;
  num_dims =    the_hmm->num_dims;

  /* free the initial state probabilities */
  NR_free_dvector(the_hmm->prob_init, 1, num_states);

  /* free the state transition matrix */
  NR_free_dmatrix(the_hmm->prob_trans, 1, num_states, 1, num_states); 

  /* free the output probability values */
  switch (model_type)
  {
    case HMM_DISCRETE_MODEL:
      NR_free_dmatrix(the_hmm->prob_symbol, 1, num_states, 1, num_symbols);
      break;

    case HMM_GAUSSIAN_MODEL:
    case HMM_MIXTURE_MODEL:
      /* free mean vectors for the mixture for each state */
      free_set_of_dmatrices(the_hmm->means, num_states, num_comps, num_dims);

      /* free covariance matrices for the mixture for each state */
      free_set_of_sets_of_dmatrices(the_hmm->covars, num_states, num_comps, 
                                   num_dims, num_dims);

      /* free mixture weights for the mixture for each state */
      NR_free_dmatrix(the_hmm->weights, 1, num_states, 1, num_comps);
      break;

    default:
      err_printf();
      log_printf("Illegal model type specified: %d\n", model_type);
      break;
  }

  free(the_hmm);

  /* no error code returned, since this mirrors the NR free_...() functions */
}


/*******************************************************************************
LEARN_DISC_HMM_USING_EM
Given an empty hmm structure, learn the parameters of the hidden 
Markov model conditioned on a supplied observation sequence, using an EM
algorithm.  Returns the likelihood of the model found upon success, -1 (an
impossible value for the likelihood) otherwise.

Note that the discrete data must be a 1-dimensional array of ints.  To handle
multi-d discrete one day, a simple approach is to take the cartesian product 
of the alphabets of the various discrete dimensions and use that as a single
super discrete variable.

The second HMM structure passed to this function is used to hold intermediate
computations, namely, the hmm parameters corresponding to the last iteration
of the EM algorithm.  The first HMM passed in will end up containing the final
learned model.  The second HMM needs to have been allocated and have its values
identical to the first except that its array memory need not be initialized to
any meaningful values.
AG
*******************************************************************************/
double learn_disc_hmm_using_em(HMM *the_hmm, HMM *last_hmm, 
                              int *data, int num_data, 
                              double **alpha, double **beta, double *scale,
                              double **sum_t_xi, double **sum_t_gamma_by_symbol,
                              double *sum_t_gamma, double **zeta_t, 
                              double *gamma_t, double threshold)

{
  int      numiters;
  boolean  converged;
  double    log_likelihood, last_log_likelihood, diff;
  HMM      *tmp_hmm;

  /* perform EM iterations */
  
  /* set initial conditions for the convergence */
  numiters = 0;
  converged = UT_FALSE;
  log_likelihood = 10000;   /* should be the largest double possible */

  while (converged == UT_FALSE)
  {

    /* initialize for this iteration */
    last_log_likelihood = log_likelihood;

    /* E-step - "Expectation" */
    /* I. Compute the likelihood of the model given the data */
    log_likelihood = compute_disc_hmm_likelihood(the_hmm, data, num_data, 
                                                 alpha, beta, scale);

    /* see how much the likelihood improved */
    diff = (double) fabs( (double)( log_likelihood - last_log_likelihood ) );
                                   
    if (verbose_output())
    {
      log_printf("log-likelihood = %g\n",log_likelihood);
      log_printf("diff = %g\n",diff);
    }

    /* check for convergence */
      /*    if ((log_likelihood < 174.0) && (log_likelihood > 173.0)) */
    if (diff <= threshold) 
    {
      if (verbose_output())
        log_printf("Converged.\n\n");

      converged = UT_TRUE;
      continue;  /* this will jump back out to the while condition */
    }

    /* switch the two HMM structures so that the current one is now considered
       to be the last HMM */
    tmp_hmm = last_hmm; last_hmm = the_hmm; the_hmm = tmp_hmm;

    /* M-step - "Maximization" */
    /* II. Compute the model parameters given the probability of each datum
       with respect to the model */
    estimate_disc_hmm_params(last_hmm, the_hmm, data, num_data, alpha, beta, 
                             sum_t_xi, sum_t_gamma_by_symbol, sum_t_gamma,
                             zeta_t, gamma_t);

    /* keep track of num iterations */
    numiters++;
    if (verbose_output())
      log_printf("finished iteration# %d\n", numiters);

  } /* end of main em loop over iterations */

  if (verbose_output())
    log_printf("total iterations = %d\n\n",numiters);

  return(log_likelihood);
}


/*******************************************************************************
LEARN_CONT_HMM_USING_EM
Given an empty hmm structure, learn the parameters of the hidden 
Markov model conditioned on a supplied observation sequence, using an EM
algorithm.  Returns the likelihood of the model found upon success, -1 (an
impossible value for the likelihood) otherwise.

Note that the data is a 2-d array of doubles, so that multi-d (vector) data 
can be handled.

The second HMM structure passed to this function is used to hold intermediate
computations, namely, the hmm parameters corresponding to the last iteration
of the EM algorithm.  The first HMM passed in will end up containing the final
learned model.  The second HMM needs to have been allocated and have its values
identical to the first except that its array memory need not be initialized to
any meaningful values.
AG
*******************************************************************************/
double learn_cont_hmm_using_em(HMM *the_hmm, HMM *last_hmm,
                              double **data, int num_data,
                              double **alpha, double **beta, double *scale,
                              double **sum_t_xi, double **sum_t_gamma_by_comp,
                              double *sum_t_gamma, double **zeta_t,double *gamma_t,
                              double *mixture_probs, double **comp_wgtd_probs, 
                              double *diff_dvec, double threshold, double min_diag)

{
  int      numiters;
  boolean  converged;
  double    log_likelihood, last_log_likelihood, diff;
  HMM      *tmp_hmm;

  /* perform EM iterations */
  
  /* set initial conditions for the convergence */
  numiters = 0;
  converged = UT_FALSE;
  log_likelihood = 10000;   /* should be the largest double possible */

  while (converged == UT_FALSE)
  {

    /* initialize for this iteration */
    last_log_likelihood = log_likelihood;

    /* E-step - "Expectation" */
    /* I. Compute the likelihood of the model given the data */
    log_likelihood = compute_cont_hmm_likelihood(the_hmm, data, num_data, 
                                                 alpha, beta, scale, min_diag);

    /* see how much the likelihood improved */
    diff = (double) fabs( (double)( log_likelihood - last_log_likelihood ) );
                                   
    if (verbose_output())
    {
      log_printf("log-likelihood = %g\n",log_likelihood);
      log_printf("diff = %g\n",diff);
    }

    /* check for convergence */
    if (diff / fabs(log_likelihood) <= threshold)
    {
      if (verbose_output())
        log_printf("Converged.\n\n");

      converged = UT_TRUE;
      continue;  /* this will jump back out to the while condition */
    }

    /* switch the two HMM structures so that the current one is now considered
       to be the last HMM */
    tmp_hmm = last_hmm; last_hmm = the_hmm; the_hmm = tmp_hmm;

    /* M-step - "Maximization" */
    /* II. Compute the model parameters given the probability of each datum
       with respect to the model */
    estimate_cont_hmm_params(last_hmm, the_hmm, data, num_data, alpha, beta, 
                             sum_t_xi, sum_t_gamma_by_comp, sum_t_gamma,
                             zeta_t, gamma_t, mixture_probs, comp_wgtd_probs,
                             diff_dvec, min_diag);

    /* keep track of num iterations */
    numiters++;
    if (verbose_output())
      log_printf("finished iteration# %d\n", numiters);

  } /* end of main em loop over iterations */

  if (verbose_output())
    log_printf("total iterations = %d\n\n",numiters);

  return(log_likelihood);
}


/*******************************************************************************
COMPUTE_DISC_HMM_LIKELIHOOD
Given the parameters of a discrete hidden Markov model and an observation 
sequence, compute the log-likelihood of the model given the data, using the 
forward-backward algorithm.  Requires that allocated space for storing the 
matrices of intermediate probabilities alpha and beta be passed in, as well
as space for holding the intermediate scaling factors (one per datum).

Note that the discrete data must be a 1-dimensional array of ints.  To handle
multi-d discrete one day, a simple approach is to take the cartesian product 
of the alphabets of the various discrete dimensions and use that as a single
super discrete variable.

This procedure incorporates scaling as described in Rabiner, p. 282.  Scaling
affects the computation of alpha, beta, and the likelihood, and is confined
to this function only.
AG
*******************************************************************************/
double compute_disc_hmm_likelihood(HMM *the_hmm, int *data, int num_data, 
                                  double **alpha, double **beta, double *scale)

{
  /* indices */
  int      i, j, t;

  /* hmm parameters */
  int      num_states;
  double    **prob_trans, **prob_symbol, *prob_init;

  /* temporary references to places within matrices */
  double    *alpha_t, *alpha_next_t, *beta_t, *beta_next_t;
  int      curr_symbol, next_symbol;

  /* computed values */
  double    scale_t, scale_next_t;
  double    sum_i_alpha_times_a, beta_ti;
  double    log_likelihood;

  /* grab the hmm parameters stored in the hmm structure */
  num_states =  the_hmm->num_states;
  prob_trans =  the_hmm->prob_trans;
  prob_symbol = the_hmm->prob_symbol;
  prob_init =   the_hmm->prob_init;

  /* I. COMPUTE FORWARD VARIABLES (ALPHAS) */

  /* t = 1 case: */

  /* for each state, initialize alpha variables */
  t = 1;   /* set values for t = 1 */
  curr_symbol = data[t];  /* O_1 */
  alpha_t = alpha[t];
  scale_t = 0.0;

  for (i = 1; i <= num_states; i++)
  {
    /* compute the beginning alphas based on the initial probabilities in pi */

    /* alpha_t(i) = pi_i * b_i( O_1 ) */
    alpha_t[i] = prob_init[i] * prob_symbol[i][curr_symbol];

    /* accumulate the scale factor c_t = sum_i{ alpha_t(i) } */
    scale_t += alpha_t[i];
  }

  /* store the scale for this t, for use on the betas below */
  scale[t] = scale_t;

  /* scale the new set of alphas to be a number from 0 to 1 */
  scalar_div_dvec( alpha_t, num_states, scale_t);

  /* hack to deal with zeroed alphas */
  if (scale_t == 0.0)
    for (i = 1; i <= num_states; i++)
      alpha_t[i] = 1.0;

  /* t = 2 to T case:  (implementing using 1 to T-1) */

  /* for each observation, starting from the beginning */
  for (t = 1; t < num_data; t++)
  {
    /* set up references/values for this t */
    next_symbol = data[t+1];  /* O_t+1 */
    alpha_t = alpha[t];
    alpha_next_t = alpha[t+1];
    scale_next_t = 0.0;

    /* update the set of states for the next time step */
    for (j = 1; j <= num_states; j++)
    {
      /* compute next set of alphas based on current set */

      /* sum the inputs to this state */
      /* compute sum_i{ alpha_t(i) * a_ij } */
      sum_i_alpha_times_a = 0.0;
      for (i = 1; i <= num_states; i++)
        sum_i_alpha_times_a += alpha_t[i] * prob_trans[i][j];

      /* alpha_t+1(j) = [ sum_i{ alpha_t(i) * a_ij } ] * b_j( O_t+1 ) */
      alpha_next_t[j] = sum_i_alpha_times_a * prob_symbol[j][next_symbol];

      /* accumulate the scale factor c_t+1 = sum_i{ alpha_t+1(i) } */
      scale_next_t += alpha_next_t[j];
    }

    /* store the scale for this t, for use on the betas below */
    scale[t+1] = scale_next_t;

    /* scale the new set of alphas to be a number from 0 to 1 */
    scalar_div_dvec( alpha_next_t, num_states, scale_next_t );

    /* hack to deal with zeroed alphas */
    if (scale_next_t == 0.0)
      for (i = 1; i <= num_states; i++)
	alpha_next_t[i] = 1.0;
  }

  /* compute log-likelihood */
  /* p(O|model) = sum_i{ alpha_T(i) } without scaling */
  /* p(O|model) = - sum_t{ log(c_t) } with scaling */
  log_likelihood = -1.0 * sum_log_dvec(scale, num_data);

  /* II. COMPUTE BACKWARD VARIABLES (BETAS) */

  /* t = T case: */

  /* for each state, initialize beta variables */
  t = num_data;   /* set values for t = T */
  beta_t = beta[t];
  for (i = 1; i <= num_states; i++)
    beta_t[i] = 1;

  /* t = T-1 to 1 case: */

  /* for each observation, starting from the end */
  for (t = num_data - 1; t >= 1; t--)
  {
    /* set up references/values for this t */
    beta_t = beta[t];
    scale_t = scale[t];       /* use the same scale factor as used for the 
                                 same time step in the alpha's above */
    beta_next_t = beta[t+1];  /* "next" here means next in time, not next in 
                                 this loop; remember we are going backward */
    next_symbol = data[t+1];  /* O_t+1 */

    /* update the set of states for the next time step */
    for (i = 1; i <= num_states; i++)
    {
      /* beta_t(i) = sum_i{ a_ij * b_j( O_t+1 ) * beta_t+1(j) } */
      beta_ti = 0.0;
      for (j = 1; j <= num_states; j++)
      {
        beta_ti += prob_trans[i][j] * prob_symbol[j][next_symbol] *
                   beta_next_t[j];
      }

      /* then beta_t(i) /= c_t, the scaling factor computed above from alphas */
      beta_ti /= scale_t;

      /* hack to deal with zeroed betas */
      if (scale[t+1] == 0.0)
        beta_ti = 1.0;

      beta_t[i] = beta_ti;
    }
  }

  return( log_likelihood );
}


/*******************************************************************************
COMPUTE_CONT_HMM_LIKELIHOOD
Given the parameters of a continuous hidden Markov model and an observation 
sequence, compute the log-likelihood of the model given the data, using the 
forward-backward algorithm.  Requires that allocated space for storing the 
matrices of intermediate probabilities alpha and beta be passed in, as well
as space for holding the intermediate scaling factors (one per datum).

Note that the data is a 2-d array of doubles, so that multi-d (vector) data 
can be handled.

This procedure incorporates scaling as described in Rabiner, p. 282.  Scaling
affects the computation of alpha, beta, and the likelihood, and is confined
to this function only.
AG
*******************************************************************************/
double compute_cont_hmm_likelihood(HMM *the_hmm, double **data, int num_data, 
                                  double **alpha, double **beta, double *scale,
                                  double min_diag) 

{
  /* indices */
  int      i, j, t;

  /* hmm parameters */
  int      num_dims, num_states, num_comps;
  double    **prob_trans, *prob_init;

  /* temporary references to places within matrices */
  double    *alpha_t, *alpha_next_t, *beta_t, *beta_next_t;
  double    ***means, **means_i, ****covars, ***covars_i;
  double    **means_j, ***covars_j;
  double    **weights, *weights_i, *weights_j;
  double    *curr_value, *next_value;

  /* computed values */
  double    scale_t, scale_next_t;
  double    sum_i_alpha_times_a, beta_ti, prob_value_t, prob_value_next_t;
  double    log_likelihood;

  /* grab the hmm parameters stored in the hmm structure */
  num_dims =    the_hmm->num_dims;
  num_states =  the_hmm->num_states;
  num_comps =   the_hmm->num_comps;
  prob_trans =  the_hmm->prob_trans;
  prob_init =   the_hmm->prob_init;
  means =       the_hmm->means;
  covars =      the_hmm->covars;
  weights =     the_hmm->weights;

  /* I. COMPUTE FORWARD VARIABLES (ALPHAS) */

  /* t = 1 case: */

  /* for each state, initialize alpha variables */
  t = 1;   /* set values for t = 1 */
  curr_value = data[t];  /* O_1 */
  alpha_t = alpha[t];
  scale_t = 0.0;

  for (i = 1; i <= num_states; i++)
  {
    /* compute the beginning alphas based on the initial probabilities in pi */

    /* set up references/values for this i */
    means_i = means[i];
    covars_i = covars[i];
    weights_i = weights[i];

    /* compute b_i( O_1 ) */
    prob_value_t = prob_mixture(curr_value, means_i, covars_i, weights_i,
                                num_dims, num_comps, min_diag);

    /* alpha_t(i) = pi_i * b_i( O_1 ) */
    alpha_t[i] = prob_init[i] * prob_value_t;
/* printf("alpha_t[%d] = %g, t = %d; prob_init[%d] = %g; prob_value_t = %g\n", i, alpha_t[i], t, i, prob_init[i], prob_value_t); printf("DBL_MIN = %g\n", DBL_MIN); */
    /* accumulate the scale factor c_t = sum_i{ alpha_t(i) } */
    scale_t += alpha_t[i];
  }

  /* store the scale for this t, for use on the betas below */
  scale[t] = scale_t;

  /* scale the new set of alphas to be a number from 0 to 1 */
  scalar_div_dvec( alpha_t, num_states, scale_t);


  /* t = 2 to T case:  (implementing using 1 to T-1) */

  /* for each observation, starting from the beginning */
  for (t = 1; t < num_data; t++)
  {
    /* set up references/values for this t */
    next_value = data[t+1];  /* O_t+1 */
    alpha_t = alpha[t];
    alpha_next_t = alpha[t+1];
    scale_next_t = 0.0;

    /* update the set of states for the next time step */
    for (j = 1; j <= num_states; j++)
    {
      /* compute next set of alphas based on current set */

      /* set up references/values for this j */
      means_j = means[j];
      covars_j = covars[j];
      weights_j = weights[j];

      /* sum the inputs to this state */
      /* compute sum_i{ alpha_t(i) * a_ij } */
      sum_i_alpha_times_a = 0.0;
      for (i = 1; i <= num_states; i++)
        sum_i_alpha_times_a += alpha_t[i] * prob_trans[i][j];

      /* compute b_j( O_t+1 ) */
      prob_value_next_t = prob_mixture(next_value, means_j, covars_j, 
                                       weights_j, num_dims, num_comps,min_diag);
                                
      /* alpha_t+1(j) = [ sum_i{ alpha_t(i) * a_ij } ] * b_j( O_t+1 ) */
      alpha_next_t[j] = sum_i_alpha_times_a * prob_value_next_t;

      /* accumulate the scale factor c_t+1 = sum_i{ alpha_t+1(i) } */
      scale_next_t += alpha_next_t[j];
    }

    /* store the scale for this t, for use on the betas below */
    scale[t+1] = scale_next_t;

    /* scale the new set of alphas to be a number from 0 to 1 */
    scalar_div_dvec( alpha_next_t, num_states, scale_next_t );
  }
/* for(t=1;t<=num_data;t++) printf("scale[%d]=%g\n",t,scale[t]); */
  /* compute log-likelihood */
  /* p(O|model) = sum_i{ alpha_T(i) } without scaling */
  /* p(O|model) = - sum_t{ log(c_t) } with scaling */
  /*
  log_likelihood = -1.0 * sum_log_dvec(scale, num_data);
  */
  log_likelihood = sum_log_dvec(scale, num_data);


  /* II. COMPUTE BACKWARD VARIABLES (BETAS) */

  /* t = T case: */

  /* for each state, initialize beta variables */
  t = num_data;    /* set values for t = T */
  beta_t = beta[t];
  for (i = 1; i <= num_states; i++)
    beta_t[i] = 1;

  /* t = T-1 to 1 case: */

  /* for each observation, starting from the end */
  for (t = num_data - 1; t >= 1; t--)
  {
    /* set up references/values for this t */
    beta_t = beta[t];
    scale_t = scale[t];       /* use the same scale factor as used for the 
                                 same time step in the alpha's above */
    beta_next_t = beta[t+1];  /* "next" here means next in time, not next in 
                                 this loop; remember we are going backward */
    next_value = data[t+1];   /* O_t+1 */

    /* update the set of states for the next time step */
    for (i = 1; i <= num_states; i++)
    {
      /* beta_t(i) = sum_i{ a_ij * b_j( O_t+1 ) * beta_t+1(j) } */
      beta_ti = 0.0;
      for (j = 1; j <= num_states; j++)
      {
        /* set up references/values for this j */
        means_j = means[j];
        covars_j = covars[j];
        weights_j = weights[j];

        /* compute b_j( O_t+1 ) */
        prob_value_next_t = prob_mixture(next_value, means_j, covars_j, 
                                         weights_j, num_dims, num_comps, 
                                         min_diag);
                                  
        beta_ti += prob_trans[i][j] * prob_value_next_t * beta_next_t[j];
      }

      /* then beta_t(i) /= c_t, the scaling factor computed above from alphas */
      beta_ti /= scale_t;
      
      /* correct floating point overflow if it occurs */
      if (beta_ti > DBL_MAX) 
	beta_ti = DBL_MAX;

      beta_t[i] = beta_ti;
    }
  }

  return( log_likelihood );
}


/*******************************************************************************
ESTIMATE_DISC_HMM_PARAMS
Given the matrices of intermediate probabilities alpha and beta, and the 
previous values of the parameters of the discrete HMM, reestimate the HMM 
parameters.  Requires that allocated space for storing various intermediate 
probabilities be passed into this function.

Note that the discrete data must be a 1-dimensional array of ints.  To handle
multi-d discrete one day, a simple approach is to take the cartesian product 
of the alphabets of the various discrete dimensions and use that as a single
super discrete variable.

Note that this form of the reestimation equations is not yet designed to 
handle the case of multiple observation sequences.  p. 283 Rabiner.
AG
*******************************************************************************/
int estimate_disc_hmm_params(HMM *old_hmm, HMM *new_hmm, 
                             int *data, int num_data, 
                             double **alpha, double **beta, 
                             double **sum_t_xi, double **sum_t_gamma_by_symbol,
                             double *sum_t_gamma, double **zeta_t, double *gamma_t)

{
  /****** TODO ******
    - put in a check for near-zero state transitions, alphas, or betas
        (perhaps after some minimum number of iterations)
  ******************/

  /* indices */
  int      i, j, k, t;

  /* hmm parameters */
  int      num_states, num_symbols;
  double    **prob_trans, *prob_init, **prob_symbol;
  double    **new_prob_trans, *new_prob_init, **new_prob_symbol;

  /* temporary references to places within matrices */
  double    *prob_trans_i, *prob_symbol_j;
  double    *new_prob_trans_i, *new_prob_symbol_j;
  double    *alpha_t, *beta_t, *beta_next_t;
  double    *sum_t_xi_i, *sum_t_gamma_by_symbol_j;
  double    **xi_t;
  int      curr_symbol, next_symbol;

  /* intermediate/temporary computed values */
  double    alpha_ti, zeta_tij; 
  double    sum_t_gamma_i, sum_t_gamma_j;
  double    prob_data_to_t;


  /* grab the hmm parameters stored in the hmm structure */
  num_states =  old_hmm->num_states;
  num_symbols = old_hmm->num_symbols;

  prob_trans =  old_hmm->prob_trans;
  prob_symbol = old_hmm->prob_symbol;
  prob_init =   old_hmm->prob_init;

  new_prob_trans =  new_hmm->prob_trans;
  new_prob_symbol = new_hmm->prob_symbol;
  new_prob_init =   new_hmm->prob_init;

  /* initialize summands */
  set_dmat(sum_t_xi, num_states, num_states, 0.0);
  set_dmat(sum_t_gamma_by_symbol, num_states, num_symbols, 0.0);
  set_dvec(sum_t_gamma, num_states, 0.0);

  set_dmat(new_prob_trans, num_states, num_states, 0.0);
  set_dmat(new_prob_symbol, num_states, num_symbols, 0.0);
  set_dvec(new_prob_init, num_states, 0.0);

  /* computations for t=1 case, to start off the loop over t */

  t = 1;  /* set values for t=1 */
  alpha_t = alpha[t];
  beta_t  = beta[t];

  /* pi_i = gamma_t(i) where t=1, and
     gamma_t(i) = alpha_t(i) * beta_t(i) / p(O_to_t|model), where
     p(O_to_t|model) = sum_i{ alpha_t(i) * beta_t(i) } */
  for (i = 1; i <= num_states; i++)
    new_prob_init[i] = alpha_t[i] * beta_t[i];
  normalize_dvec(new_prob_init, num_states);


  /* computations for t = 1 to T-1; big loop over the data */

  for (t = 1; t < num_data; t++)
  {
    /* set up references/values for this t */
    alpha_t = alpha[t];
    beta_next_t = beta[t+1];
    curr_symbol = data[t];
    next_symbol = data[t+1];

    /* computations related to the state transition probabilities */

    /* a_ij = sum_t{ xi_t(i,j) } / sum_t{ gamma_t(i) }, where
       t goes from 1 to T-1 */

    /* xi_t(i,j) = zeta_t(i,j) / p(O_to_t|model), where
       zeta_t(i,j) = alpha_t(i) * a_ij * b_j( O_t+1 ) * beta_t+1(j), and
       p(O_to_t|model) = sum_i{ sum_j{ zeta_t(i,j) } }  */
    prob_data_to_t = 0.0;
    for (i = 1; i <= num_states; i++)
    {
      /* set up references/values for this i */
      prob_trans_i = prob_trans[i];
      alpha_ti = alpha_t[i];

      for (j = 1; j <= num_states; j++)
      {
        /* set up references/values for this j */
        prob_symbol_j = prob_symbol[j];

        /* compute zeta_t(i,j) */
        zeta_tij = alpha_ti * prob_trans_i[j] * prob_symbol_j[next_symbol]
                   * beta_next_t[j];
        zeta_t[i][j] = zeta_tij;

        /* accumulate p(O_to_t|model) */
        prob_data_to_t += zeta_tij;
      } 
    }   /* end of outer loop over states */
    
    /* compute each xi_t(i,j) */
    scalar_div_dmat(zeta_t, num_states, num_states, prob_data_to_t);
    xi_t = zeta_t;  /* after normalization by prob_data_to_t, we get xi_t */

    /* each sum_t_xi(i,j) += xi_t(i,j);
       note that sum_t_xi(i,j) is the numerator of the eqn. for a_ij */
    add_dmat(xi_t, sum_t_xi, sum_t_xi, num_states, num_states);

    /* gamma_t(i) = sum_j{ xi_t(i,j) } */
    sum_dmat_rows(xi_t, gamma_t, num_states, num_states);

    /* each sum_t_gamma(i) += gamma_t(i) */
    add_dvec(gamma_t, sum_t_gamma, sum_t_gamma, num_states);


    /* computations related to the output symbol probabilities */

    /* b_i(k) = sum_t{ gamma_t(i) s.t. O_t = v_k } / sum_t{ gamma_t(i) } */
    /* accumulate the numerator of this eqn. */
    for (i = 1; i <= num_states; i++)
      sum_t_gamma_by_symbol[i][curr_symbol] += gamma_t[i];

  }  /* end of loop over data */


  /* now divide by the denominator of eqn. for a_ij */
  for (i = 1; i <= num_states; i++)
  {
    /* set up references/values for this i */
    new_prob_trans_i = new_prob_trans[i];
    sum_t_xi_i = sum_t_xi[i];
    sum_t_gamma_i = sum_t_gamma[i];

    /* do the division */
    for (j = 1; j <= num_states; j++)
      new_prob_trans_i[j] = sum_t_xi_i[j] / sum_t_gamma_i;
  }


  /* more computations related to the output symbol probabilities */

  /* account for contributions from the last datum */

  /* compute gamma for the last time step */
  t = num_data;    /* set values for t=T */
  alpha_t = alpha[t];
  beta_t = beta[t];
  curr_symbol = data[t];  /* O_T */

  /* use the fact that gamma_t(i) = alpha_t(i) * beta_t(i) / p(O_to_t|model),
     where p(O_to_t|model) = sum_i{ alpha_t(i) * beta_t(i) } */
  for (j = 1; j <= num_states; j++)
    gamma_t[j] = alpha_t[j] * beta_t[j];
  normalize_dvec(gamma_t, num_states);

  /* add each gamma_T(i) to sum_t{ gamma_t(i) } */
  for (i = 1; i <= num_states; i++)
  {
    sum_t_gamma[i] += gamma_t[i];
    sum_t_gamma_by_symbol[i][curr_symbol] += gamma_t[i];
  }

  /* now divide output probabilities by the denominator, sum_t{ gamma_t(j) } */

  for (j = 1; j <= num_states; j++)
  {
    /* set up references/values for this j */
    new_prob_symbol_j = new_prob_symbol[j];
    sum_t_gamma_by_symbol_j = sum_t_gamma_by_symbol[j];
    sum_t_gamma_j = sum_t_gamma[j];

    for (k = 1; k <= num_symbols; k++)
      new_prob_symbol_j[k] = sum_t_gamma_by_symbol_j[k] / sum_t_gamma_j;
  }

  return (UT_OK);
}


/*******************************************************************************
ESTIMATE_CONT_HMM_PARAMS
Given the matrices of intermediate probabilities alpha and beta, and the 
previous values of the parameters of the continuous HMM, reestimate the HMM 
parameters.  Requires that allocated space for storing various intermediate 
probabilities be passed into this function.

Note that the data is a 2-d array of doubles, so that multi-d (vector) data 
can be handled.

Note that this form of the reestimation equations is not yet designed to 
handle the case of multiple observation sequences.  p. 283 Rabiner.
AG
*******************************************************************************/
int estimate_cont_hmm_params(HMM *old_hmm, HMM *new_hmm, 
                             double **data, int num_data, 
                             double **alpha, double **beta, 
                             double **sum_t_xi, double **sum_t_gamma_by_comp,
                             double *sum_t_gamma, double **zeta_t, double *gamma_t,
                             double *mixture_probs, double **comp_wgtd_probs, 
                             double *diff, double min_diag)

{
  /****** TODO ******
    - put in a check for near-zero state transitions, alphas, or betas
        (perhaps after some minimum number of iterations)
  ******************/

  /* indices */
  int      i, j, k, l, t, d, d2;

  /* hmm parameters */
  int      num_dims, num_states, num_comps;
  double    **prob_trans, *prob_init;
  double    **new_prob_trans, *new_prob_init;
  double    ***means, ****covars, **weights;
  double    ***new_means, ****new_covars, **new_weights;

  /* temporary references to places within matrices */
  double    *prob_trans_i, *new_prob_trans_i;
  double    *alpha_t, *beta_t, *beta_next_t;
  double    *sum_t_xi_i, *sum_t_gamma_by_comp_j;
  double    **xi_t;
  double    **means_j, *mean_jk, ***covars_j, *weights_j;
  double    **new_means_j, *new_mean_jk, ***new_covars_j, **new_covar_jk;
  double    *new_weights_j;
  double    *comp_wgtd_probs_j;
  double    *curr_value, *next_value;

  /* intermediate/temporary computed values */
  double    alpha_ti, zeta_tij;
  double    gamma_tj, sum_t_gamma_i, sum_t_gamma_j;
  double    sum_t_gamma_by_comp_jk, gamma_t_by_comp_jk;
  double    prob_data_to_t, prob_value_next_t, prob_comp_jk;
  double    mixture_probs_j;
  
  /* temp pointer to HMM */

  int  limit = 0;

  /* grab the hmm parameters stored in the hmm structures */
  num_dims =   old_hmm->num_dims;
  num_states = old_hmm->num_states;
  num_comps =  old_hmm->num_comps;

  prob_trans = old_hmm->prob_trans;
  prob_init =  old_hmm->prob_init;
  means =      old_hmm->means;
  covars =     old_hmm->covars;
  weights =    old_hmm->weights;

  new_prob_trans = new_hmm->prob_trans;
  new_prob_init =  new_hmm->prob_init;
  new_means =      new_hmm->means;
  new_covars =     new_hmm->covars;
  new_weights =    new_hmm->weights;

  /* check for parameters that alpha, beta, prob_trans close to zero.
     if some are close to zero; exit leaving all parameters the same.
     this should cause the loop to exit */
  for (i = 1; i <= num_states; i++)
    for (j = 1; j <= num_states; j++)
      if (prob_trans[i][j] < DBL_MIN) {
	copy_dmat(prob_trans, new_prob_trans, num_states, num_states);
	copy_dvec(prob_init, new_prob_init, num_states);
	for (k = 1; k <= num_states; k++) {
	  copy_dmat(means[k], new_means[k], num_comps, num_dims);
	  for (l = 1; l <= num_comps; l++)
	    copy_dmat(covars[k][l], new_covars[k][l], num_dims, num_dims);
	}
	copy_dmat(weights, new_weights, num_states, num_comps);

	return (UT_OK);
      }
  
  for (t = 1; t <= num_data; t++)
    for (i = 1; i <= num_states; i++) {
      if ((alpha[t][i] < DBL_MIN) || (beta[t][i] < DBL_MIN)) {
	copy_dmat(prob_trans, new_prob_trans, num_states, num_states);
	copy_dvec(prob_init, new_prob_init, num_states);
	for (k = 1; k <= num_states; k++) {
	  copy_dmat(means[k], new_means[k], num_comps, num_dims);
	  for (l = 1; l <= num_comps; l++)
	    copy_dmat(covars[k][l], new_covars[k][l], num_dims, num_dims);
	}
	copy_dmat(weights, new_weights, num_states, num_comps);

	return (UT_OK);
      }
    }

  /* initialize summands */
  set_dmat(sum_t_xi, num_states, num_states, 0.0);
  set_dmat(sum_t_gamma_by_comp, num_states, num_comps, 0.0);
  set_dvec(sum_t_gamma, num_states, 0.0);

  for (i = 1; i <= num_states; i++)
    set_dmat(new_means[i], num_comps, num_dims, 0.0);
  for (i = 1; i <= num_states; i++)
    for (k = 1; k <= num_comps; k++)
      set_dmat(new_covars[i][k], num_dims, num_dims, 0.0);
  set_dmat(new_weights, num_states, num_comps, 0.0);

  /* computations for t=1 case, to start off the loop over t */

  t = 1;  /* set values for t=1 */
  alpha_t = alpha[t];
  beta_t  = beta[t];
  curr_value = data[t];  /* O_1 */

  /* pi_i = gamma_t(i) where t=1, and  */
  /* gamma_t(i) = alpha_t(i) * beta_t(i) / p(O_to_t|model), where */
  /* p(O_to_t|model) = sum_i{ alpha_t(i) * beta_t(i) } */
  for (i = 1; i <= num_states; i++)
    new_prob_init[i] = alpha_t[i] * beta_t[i];
  normalize_dvec(new_prob_init, num_states);

  /* compute mixture probabilities for t=1 */
  for (j = 1; j <= num_states; j++)
  {
    /* set up references/values for this j */
    means_j = means[j];
    covars_j = covars[j];
    weights_j = weights[j];
    comp_wgtd_probs_j = comp_wgtd_probs[j];
        
    /* compute b_j( O_1 ) */
    mixture_probs[j] = prob_mixture_and_keep(curr_value, means_j, covars_j,
                                             weights_j, num_dims, num_comps,
                                             min_diag, comp_wgtd_probs_j);
  }

  /* computations for t = 1 to T-1; big loop over the data */
  
  for (t = 1; t < num_data; t++)
  {
    /* set up references/values for this t */
    alpha_t = alpha[t];
    beta_t = beta[t];
    beta_next_t = beta[t+1];
    curr_value = data[t];
    next_value = data[t+1];
    
    /* computations related to the state transition probabilities */

    /* a_ij = sum_t{ xi_t(i,j) } / sum_t{ gamma_t(i) }, where
       t goes from 1 to T-1 */

    /* xi_t(i,j) = zeta_t(i,j) / p(O_to_t|model), where
       zeta_t(i,j) = alpha_t(i) * a_ij * b_j( O_t+1 ) * beta_t+1(j), and
       p(O_to_t|model) = sum_i{ sum_j{ zeta_t(i,j) } }  */
    prob_data_to_t = 0.0;
    for (i = 1; i <= num_states; i++)
    {
      /* set up references/values for this i */
      prob_trans_i = prob_trans[i];
      alpha_ti = alpha_t[i];

      for (j = 1; j <= num_states; j++)
      {
        /* set up references/values for this j */
        means_j = means[j];
        covars_j = covars[j];
        weights_j = weights[j];
        comp_wgtd_probs_j = comp_wgtd_probs[j];
        
        /* compute b_j( O_t+1 ) */
        prob_value_next_t = prob_mixture_and_keep(next_value, means_j, 
                                                  covars_j, weights_j, 
                                                  num_dims, num_comps, 
                                                  min_diag, comp_wgtd_probs_j);
        mixture_probs[j] = prob_value_next_t;

        /* compute zeta_t(i,j) */
        zeta_tij = alpha_ti * prob_trans_i[j] * prob_value_next_t
                   * beta_next_t[j];
        zeta_t[i][j] = zeta_tij;

        /* accumulate p(O_to_t|model) */
        prob_data_to_t += zeta_tij;
      } 
    }   /* end of main loop over states */

    /* compute each xi_t(i,j) */
    scalar_div_dmat(zeta_t, num_states, num_states, prob_data_to_t);
    xi_t = zeta_t;  /* after normalization by prob_data_to_t, we get xi_t */

    /* each sum_t_xi(i,j) += xi_t(i,j);
       note that sum_t_xi(i,j) is the numerator of the eqn. for a_ij */
    add_dmat(xi_t, sum_t_xi, sum_t_xi, num_states, num_states);

    /* gamma_t(i) = sum_j{ xi_t(i,j) } */
    sum_dmat_rows(xi_t, gamma_t, num_states, num_states);

    /* each sum_t_gamma(i) += gamma_t(i) */
    add_dvec(gamma_t, sum_t_gamma, sum_t_gamma, num_states);


    /* computations related to the output symbol probabilities */

    /* update the mixture weights, means, and covariance matrices corresponding
       to each state's mixture;
       these all involve the value gamma_t(j,k) */

    /* compute gamma_t(j,k) = gamma_t(j) * prob_component(j,k), where
       prob_component(j,k) = c_jk * normal(datum_t, mean_jk, covar_jk) /
                      sum_k{ c_jk * normal(datum_t, mean_jk, covar_jk) } */
    for (j = 1; j <= num_states; j++)
    {
      /* set up references/values for this j */
      means_j = means[j];
      new_means_j = new_means[j];
      new_covars_j = new_covars[j];
      gamma_tj = gamma_t[j];
      sum_t_gamma_by_comp_j = sum_t_gamma_by_comp[j];
      comp_wgtd_probs_j = comp_wgtd_probs[j];
      mixture_probs_j = mixture_probs[j];

      for (k = 1; k <= num_comps; k++)
      {
        /* set up references/values for this k */
        mean_jk = means_j[k];
        new_mean_jk = new_means_j[k];
        new_covar_jk = new_covars_j[k];

        /* note that this uses the mixture probabilities for time t, which
           were computed previously; later in this big loop we'll compute
           the probabilities for time t+1, which will get used by the next
           iteration */
        prob_comp_jk = comp_wgtd_probs_j[k] / mixture_probs_j;
        gamma_t_by_comp_jk = gamma_tj * prob_comp_jk;
                             
        if ((debug_output()) && (t <= limit))
        {
          log_printf("state %d, comp %d: prob_comp_jk = %g\n",j,k,prob_comp_jk);
          log_printf("   comp wgtd probs_jk = %g\n", comp_wgtd_probs_j[k]);
          log_printf("   mixture_probs_j = %g\n", mixture_probs_j);
          log_printf("   gamma_t_by_comp_jk = %g\n", gamma_t_by_comp_jk);
        }

        sum_t_gamma_by_comp_j[k] += gamma_t_by_comp_jk;
      
        /* accumulate numerator of mu_jk, or sum_t{ gamma_t(j,k) * O_t } */
        for (d = 1; d <= num_dims; d++)
          new_mean_jk[d] += curr_value[d] * gamma_t_by_comp_jk;
	
        /* accumulate numerator of U_jk, or 
           sum_t{ gamma_t(j,k) * (O_t - mu_jk)(O_t - mu_jk)^T }  */
        subtract_dvec(curr_value, mean_jk, diff, num_dims);
        for (d = 1; d <= num_dims; d++)
          for (d2 = 1; d2 <= num_dims; d2++)
            new_covar_jk[d][d2] += diff[d] * diff[d2] * gamma_t_by_comp_jk;
      }  /* end of loop over mixture components */
    }  /* end of loop over states */

    if ((debug_output()) && (t <= limit))
    {
      log_printf("curr_value:\n");
      print_drow(ut_log_fp, num_dims, curr_value);
      log_printf("means:\n");
      for (i=1;i<=num_states;i++)
        for (j=1;j<=num_comps;j++)
          print_drow(ut_log_fp, num_dims, means[i][j]);
      log_printf("covars:\n");
      for (i=1;i<=num_states;i++)
        for (j=1;j<=num_comps;j++)
          print_dmatrix(ut_log_fp, num_dims, num_dims, covars[i][j]);
      log_printf("weights:\n");
      print_dmatrix(ut_log_fp, num_states, num_comps, weights);
      log_printf("\n");
      log_printf("mixture probs:\n");
      print_drow(ut_log_fp, num_states, mixture_probs);
      log_printf("comp wgtd probs:\n");
      print_dmatrix(ut_log_fp, num_states, num_comps, comp_wgtd_probs);
      log_printf("\n");
    }

  }  /* end of loop over data */


  /* now divide by the denominator of eqn. for a_ij */
  for (i = 1; i <= num_states; i++)
  {
    /* set up references/values for this i */
    new_prob_trans_i = new_prob_trans[i];
    sum_t_xi_i = sum_t_xi[i];
    sum_t_gamma_i = sum_t_gamma[i];

    /* do the division */
    for (j = 1; j <= num_states; j++)
      new_prob_trans_i[j] = sum_t_xi_i[j] / sum_t_gamma_i;
  }


  /* more computations related to the output symbol probabilities */

  /* account for contributions from the last datum */

  /* compute gamma for the last time step */
  t = num_data;    /* set values for t=T */
  alpha_t = alpha[t];
  beta_t = beta[t];
  curr_value = data[t];  /* O_T */

  /* compute gamma_t(j) = alpha_t(j) * beta_t(j) / p(O_to_t|model), where
     p(O_to_t|model) = sum_j{ alpha_t(j) * beta_t(j) } */
  for (j = 1; j <= num_states; j++)
    gamma_t[j] = alpha_t[j] * beta_t[j];
  normalize_dvec(gamma_t, num_states);

  if (debug_output())
  {
    log_printf("sum_t_gamma before adding gamma_T:\n");
    print_drow(ut_log_fp, num_states, sum_t_gamma);
  }

  /* add each gamma_T(i) to sum_t{ gamma_t(i) } */
  add_dvec(gamma_t, sum_t_gamma, sum_t_gamma, num_states);

  if (debug_output())
  {
    log_printf("gamma_T:\n");
    print_drow(ut_log_fp, num_states, gamma_t);
    log_printf("sum_t_gamma after adding gamma_T:\n");
    print_drow(ut_log_fp, num_states, sum_t_gamma);
  }

  /* this whole loop should be identical to a loop inside the big loop over 
     t, above */

  /* compute gamma_t(j,k) = gamma_t(j) * prob_component(j,k), where
     prob_component(j,k) = c_jk * normal(datum_t, mean_jk, covar_jk) /
                    sum_k{ c_jk * normal(datum_t, mean_jk, covar_jk) } 

     note that the mixture probabilities for t=T exist from the last iteration
        of the big loop above */
  for (j = 1; j <= num_states; j++)
  {
    /* set up references/values for this j */
    means_j = means[j];
    new_means_j = new_means[j];
    new_covars_j = new_covars[j];
    gamma_tj = gamma_t[j];
    sum_t_gamma_by_comp_j = sum_t_gamma_by_comp[j];
    comp_wgtd_probs_j = comp_wgtd_probs[j];
    mixture_probs_j = mixture_probs[j];

    for (k = 1; k <= num_comps; k++)
    {
      /* set up references/values for this k */
      mean_jk = means_j[k];
      new_mean_jk = new_means_j[k];
      new_covar_jk = new_covars_j[k];
      
      prob_comp_jk = comp_wgtd_probs_j[k] / mixture_probs_j;
      gamma_t_by_comp_jk = gamma_tj * prob_comp_jk;
                           
      if (debug_output())
      {
        log_printf("state %d, comp %d: prob_comp_jk = %g\n",j,k,prob_comp_jk);
        log_printf("   comp wgtd probs_jk = %g\n", comp_wgtd_probs_j[k]);
        log_printf("   mixture_probs_j = %g\n", mixture_probs_j);
        log_printf("   gamma_t_by_comp_jk = %g\n", gamma_t_by_comp_jk);
      }

      sum_t_gamma_by_comp_j[k] += gamma_t_by_comp_jk;
      
      /* accumulate numerator of mu_jk, or sum_t{ gamma_t(j,k) * O_t } */
      for (d = 1; d <= num_dims; d++)
        new_mean_jk[d] += curr_value[d] * gamma_t_by_comp_jk;
      
      /* accumulate numerator of U_jk, or 
         sum_t{ gamma_t(j,k) * (O_t - mu_jk)(O_t - mu_jk)^T }  */
      subtract_dvec(curr_value, mean_jk, diff, num_dims);
      for (d = 1; d <= num_dims; d++)
        for (d2 = 1; d2 <= num_dims; d2++)
          new_covar_jk[d][d2] += diff[d] * diff[d2] * gamma_t_by_comp_jk;
    }  /* end of loop over mixture components */
  }  /* end of loop over states */
  

  /* now divide mixture parameters by the denominator, sum_t{ gamma_t(j,k) } */

  for (j = 1; j <= num_states; j++)
  {
    /* set up references/values for this j */
    new_means_j = new_means[j];
    new_weights_j = new_weights[j];
    new_covars_j = new_covars[j];
    sum_t_gamma_j = sum_t_gamma[j];
    sum_t_gamma_by_comp_j = sum_t_gamma_by_comp[j];

    if (debug_output())
    {
      log_printf("sum_t_gamma_by_comp_j, j = %d:\n",j);
      print_drow(ut_log_fp, num_comps, sum_t_gamma_by_comp_j);
      log_printf("sum_t_gamma_j =  %g\n",sum_t_gamma_j);
    }

    for (k = 1; k <= num_comps; k++)
    {
      /* set up references/values for this k */      
      new_mean_jk = new_means_j[k];
      new_covar_jk = new_covars_j[k];
      sum_t_gamma_by_comp_jk = sum_t_gamma_by_comp_j[k];

      /* compute and normalize weights */
      new_weights_j[k] = sum_t_gamma_by_comp_j[k] / sum_t_gamma_j;

      /* normalize means */
      scalar_div_dvec( new_mean_jk, num_dims, sum_t_gamma_by_comp_jk );

      /* normalize covars */
      scalar_div_dmat( new_covar_jk, num_dims, num_dims, sum_t_gamma_by_comp_jk);
    }

    if (debug_output())
    {
      log_printf("new weights:\n");
      print_drow(ut_log_fp, num_comps, new_weights_j);
    }

  }

  return (UT_OK);
}


/*******************************************************************************
SIMULATE_DATA_USING_DISC_HMM
Use an instantiated discrete HMM to generate simulated data.

Sets the random number generator before generating data.
AG
*******************************************************************************/
int simulate_data_using_disc_hmm(HMM* the_hmm, int *data, int *states,
                                 int num_data, double *cum_prob_states, 
                                 double *cum_prob_symbols, int seed)

{
  /* indices */
  int t, curr_state;

  /* hmm parameters */
  int num_states, num_symbols;
  double *prob_init, **prob_trans, **prob_symbol;

  /* get hmm parameters for convenience */
  num_states  = the_hmm->num_states;
  num_symbols = the_hmm->num_symbols;
  prob_init   = the_hmm->prob_init;
  prob_trans  = the_hmm->prob_trans;
  prob_symbol = the_hmm->prob_symbol;

  /* set the random number generator */
  set_rand(seed);

  /* choose the first state based on the specified initial probabilities,
     then select a symbol based on the probabilities for that state */
  t = 1;
  curr_state = generate_disc_rv_value(prob_init, num_states, cum_prob_states);
  states[t]  = curr_state;
  data[t]    = generate_disc_rv_value(prob_symbol[curr_state], num_symbols, 
                                      cum_prob_symbols);

  /* for each datum to be created, transition to a randomly selected state;
     then select a symbol based on the probabilities for that state */
  for (t = 2; t <= num_data; t++)
  {
    curr_state = generate_disc_rv_value(prob_trans[curr_state], num_states, 
                                        cum_prob_states);
    states[t]  = curr_state;
    data[t]    = generate_disc_rv_value(prob_symbol[curr_state], num_symbols, 
                                        cum_prob_symbols);
  }
  
  return (UT_OK);
}


/*******************************************************************************
SIMULATE_DATA_USING_CONT_HMM
Use an instantiated continuous HMM to generate simulated data.

Sets the random number generator before generating data.
AG
*******************************************************************************/
int simulate_data_using_cont_hmm(HMM* the_hmm, double **data, int *states,
                                 int num_data, double *cum_prob_states,
                                 double *cum_prob_comps, double **chol_factor,
                                 double *diag, int seed)

{
  /* indices */
  int t, curr_state;

  /* hmm parameters */
  int   num_dims, num_comps, num_states;
  double *prob_init, **prob_trans;
  double ***means, ****covars, **weights;

  /* get hmm parameters for convenience */
  num_dims    = the_hmm->num_dims;
  num_states  = the_hmm->num_states;
  num_comps   = the_hmm->num_comps;
  prob_init   = the_hmm->prob_init;
  prob_trans  = the_hmm->prob_trans;
  means       = the_hmm->means;
  covars      = the_hmm->covars;
  weights     = the_hmm->weights;

  /* set the random number generator */
  set_rand(seed);

  /* choose the first state based on the specified initial probabilities,
     then select a symbol based on the probabilities for that state */
  t = 1;
  curr_state = generate_disc_rv_value(prob_init, num_states, cum_prob_states);
  states[t]  = curr_state;
  generate_mixture_rv_value(data[t], num_dims, num_comps, means[curr_state], 
                            covars[curr_state], weights[curr_state], 
                            cum_prob_comps, chol_factor, diag);
                            

  /* for each datum to be created, transition to a randomly selected state;
     then select a symbol based on the probabilities for that state */
  for (t = 2; t <= num_data; t++)
  {
    curr_state = generate_disc_rv_value(prob_trans[curr_state], num_states, 
                                        cum_prob_states);
    states[t] = curr_state;
    generate_mixture_rv_value(data[t], num_dims, num_comps, means[curr_state], 
                              covars[curr_state], weights[curr_state], 
                              cum_prob_comps, chol_factor, diag);
  }
  
  return (UT_OK);
}


/*******************************************************************************
ASSIGN_STATES_USING_IML
Assign a state label to each datum, using the individually most likely (IML)
state, i.e. the state with the highest probability at the time of the datum
given the observation sequence.
AG
*******************************************************************************/
int assign_states_using_iml(int *states, double **alpha, double **beta, 
                            int num_states, int num_data)

{
  int i, t;
  int arg_max_prob;
  int num_equal;
  double *alpha_t, *beta_t, gamma_ti, max_prob;

  for (t = 1; t <= num_data; t++)
  {
    /* set up references/values for this t */
    alpha_t = alpha[t];
    beta_t  = beta[t];
    gamma_ti = 0.0;
    max_prob = 0.0;
    arg_max_prob = -1;
    num_equal = 1;

    /* sum_i_gamma_t = 0.0; */

    for (i = 1; i <= num_states; i++)
    {
      gamma_ti = alpha_t[i] * beta_t[i];
      /* sum_i_gamma_t += gamma_ti; */
      if (gamma_ti > max_prob)
      {
        max_prob = gamma_ti;
        arg_max_prob = i;
        num_equal = 1;
      }
      else if (gamma_ti == max_prob)
	num_equal++;
    }

    if (arg_max_prob == -1)
    {
      err_printf();
      log_printf("No state was more probable than zero.\n");
    }

    /* note that computing gamma_ti properly requires normalizing by
       sum_i_gamma_t, but that the max of the numerator is all we need here */
    states[t] = arg_max_prob;

    if (num_equal == num_states)
      states[t] = -1;
  }

  return (UT_OK);
}


/*******************************************************************************
 INIT_HMM_MODEL_TYPE_NAMES
 Initialize the name strings corresponding to the available model type options.
*******************************************************************************/
int init_hmm_model_type_names ( char **model_type_names )
{
  model_type_names[HMM_DISCRETE_MODEL] = (char*) 
                                         strdup (HMM_DISCRETE_MODEL_NAME);
  model_type_names[HMM_GAUSSIAN_MODEL] = (char*) 
                                         strdup (HMM_GAUSSIAN_MODEL_NAME);
  model_type_names[HMM_MIXTURE_MODEL] = (char*) 
                                        strdup (HMM_MIXTURE_MODEL_NAME);

  return(UT_OK);
}

