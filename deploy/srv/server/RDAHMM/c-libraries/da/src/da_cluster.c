/*******************************************************************************
MODULE NAME
da_cluster

ONE-LINE SYNOPSIS
Functions related to clustering.

SCOPE OF THIS MODULE
Any functions relating directly to static clustering should fall into this 
module.  Functions relating to clustering in time series, or time series
segmentation, should instead go into the da_timeseg module.

SEE ALSO
Many functions which are analogous to those in this module are in da_timeseg.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/sc, AG.
2. /proj/cooltools/kmeans, AG.


NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_cluster.c,v 1.15 1998/06/26 02:04:31 agray Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_cluster.c,v $
 * Revision 1.15  1998/06/26 02:04:31  agray
 * updated some function names; updated includes
 *
 * Revision 1.14  1998/05/07 23:51:48  granat
 * made changes to include new da_util module
 *
 * Revision 1.13  1998/03/09 00:38:06  granat
 * fixed boolean bug
 *
 * Revision 1.12  1997/06/05 18:55:23  granat
 * edited to conform to changes in da_linalg
 *
 * Revision 1.11  1997/06/02 15:34:48  granat
 * changed to use new NR naming convention
 *
 * Revision 1.10  1997/04/05 19:07:59  granat
 * made adjustments to account for changes in da_linalg
 *
 * Revision 1.9  1997/01/29 21:42:41  agray
 * new formatting, cleaning up debugging output using ut_output.
 *
 * Revision 1.8  1996/10/31 02:03:03  agray
 * renamed from "da_clust" to "da_cluster";
 * added functions from HMM project;
 * changed .h and .c formats throughout library;
 * some reorganizing between modules.
 * 
 * Revision 1.7  1996/10/30 20:29:09  agray
 * changed naming from "mixmod" to "mixture", other naming changes.
 * 
 * Revision 1.6  1996/09/19 22:45:22  agray
 * changed name of prob_gauss() to gauss_eval(); similar for vector_prob_gauss().
 * 
 * Revision 1.5  1996/08/28 19:41:15  agray
 * changed random_means() and random_means_from_data() to wrappers around
 * functions in da_rand.
 * 
 * Revision 1.4  1996/07/17 20:41:48  agray
 * cosmetic.
 * 
 * Revision 1.3  1996/07/11 16:38:23  agray
 * major edits to kmeans(); changed randomize_means() to random_means();
 * added random_means_from_data(), mixmod_likelihood(), estimate_mixmod_
 * params().
 * 
 * Revision 1.2  1996/05/14 01:03:02  agray
 * minor compile problems.
 * 
 * Revision 1.1  1996/05/07 20:47:06  agray
 * Initial revision
 * 
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* UT library */
#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_util.h"
#include "da_io.h"
#include "da_linalg.h"
#include "da_probstat.h"
#include "da_random.h"

/* this module's header */
#include "da_cluster.h"


/*******************************************************************************
K_MEANS
Perform the k-means clustering algorithm on a dataset of continuous values,
returning the matrix of k mean vectors and a class labelling for each of the
data.  The input matrices and vectors are all assumed to have been allocated
before calling this function.  The argument initmeans should contain the
mean vectors which will start off the k-means iterations.  The data is assumed
to be stored in row vector form (i.e. each row is a datum).

NOTE: In practice, it is a good idea to normalize the columns of the data to
the same scale before performing k-means. (See the function range_normalize_
cols()).
AG
*******************************************************************************/
int k_means(double **data, int numrows, int numcols, double **means, 
            double **initmeans, int *clustsize, int *class, int K,
            int *numiters, int *num_empty_clusters, boolean empty_cluster_ok)

{
    double          **oldmeans, **tempmeans;
    int            i, j, k;
    int            bestclass;
    double          dist, mindist;
    boolean        stable;
    int            empty_clusters, iters;

    /* initialize before iterating */
    stable = UT_FALSE;
    empty_clusters = 0;
    iters = 0;
    
    /* create storage for temporary information */
    oldmeans = initmeans;

    /* perform k-means iterations */
    while (stable == UT_FALSE)
    {
      stable = UT_TRUE;

      /* zero out new means */
      for (k=1; k<=K; k++)
        for (j=1; j<=numcols; j++)
          means[k][j] = 0;

      /* zero out new cluster sizes */
      for (k=1; k<=K; k++)
        clustsize[k] = 0;

      /* compute distance of each datum to each mean */
      for (i=1; i<=numrows; i++)
      {
       
        mindist = euclid_dist_dvec(data[i],oldmeans[1],numcols);
        bestclass = 1;
          /* 
            printf("item %d to mean %d: %f\n",
                   i,1,euclid_dist_dvec(data[i],oldmeans[1],numcols));
             */
        for (k=2; k<=K; k++)
        {
          /*
            printf("item %d to mean %d: %f\n",
                   i,k,euclid_dist_dvec(data[i],oldmeans[k],numcols));
           */

          if ((dist = euclid_dist_dvec(data[i],oldmeans[k],numcols)) < mindist)
          {
            mindist = dist;
            bestclass = k;
          }
        }

        /* assign class and see if different from last iteration */
        if (class[i] != bestclass)
        {
          class[i] = bestclass;
          stable = UT_FALSE;
        }

        /* contribute to the new mean based on class assignment */
        for (j=1; j<=numcols; j++)
          means[class[i]][j] += data[i][j];
        clustsize[class[i]] += 1;

      }

      /* keep track of num iterations */
      iters++;
      if (verbose_output())
        log_printf("iteration %d\n", iters);

      /* check for eliminated clusters */
      for (k=1; k<=K; k++)
        if (clustsize[k] == 0)
          empty_clusters++;

      /* decide what to do if a cluster was eliminated, based on flag */
      if (empty_clusters > 0)
      {
        if (empty_cluster_ok == UT_FALSE)
        {
          /* set some values to be returned */
          *num_empty_clusters = empty_clusters;
          *numiters = iters;

          return (UT_ERROR);
        }
        else
        {
          K -= empty_clusters;  /* caller must also do this */
          /* also get rid of the corresponding means in the means matrix */
          printf("This option does not work yet.\n");

          /* set some values to be returned */
          *num_empty_clusters = empty_clusters;
          *numiters = iters;

          return (UT_ERROR);  /* for now... when this works, just continue
                                 with smaller K */
        }
      }

      /* finish computing new means */
      for (k=1; k<=K; k++)
        for (j=1; j<=numcols; j++)
          means[k][j] /= (double) clustsize[k];

      /* make the new means into the oldmeans for next iteration */
      tempmeans = oldmeans; oldmeans = means; means = tempmeans;

      /* 
      print_dmatrix(ut_log_fp, K, numcols, means);
       */

    }
    
    /* set some values to be returned */
    *num_empty_clusters = empty_clusters;
    *numiters = iters;

    return (UT_OK);
}


/*******************************************************************************
K_MEANS_MIXTURE_ESTIMATE
Compute covariance matrices and mixture weights based on the means and cluster
memberships resulting from a k-means run.
AG
*******************************************************************************/
int k_means_mixture_estimate(double **data, double **means, double ***covars, 
                             double *weights, int *clust_label, 
                             int num_data, int num_dims, int K)

{
  int i, j, j2, k;
  double **covar_k, *mean_k, *data_i;

  /* set up parameters for summing */
  set_dvec(weights, K, 0.0);
  for (k = 1; k <= K; k++)
    set_dmat(covars[k], num_dims, num_dims, 0.0);

  for (i = 1; i <= num_data; i++)
  {
    k = clust_label[i];
        
    /* get the parameters corresponding to the class of this datum */
    covar_k = covars[k];
    mean_k  = means[k];

    /* estimate mixture weights given members of k-means clusters */
    /* add the contribution of this datum to its class */
    weights[k]++;
        
    /* estimate initial covariance matrices given members of k-means 
       clusters */
    /* add the square of the deviation vectors to the right position in the
       covariance */
    for (j = 1; j <= num_dims; j++)
    {
      data_i = data[i];
      for (j2 = j; j2 <= num_dims; j2++)
        covar_k[j][j2] += (data_i[j] - mean_k[j]) * (data_i[j2] - mean_k[j2]);
    }
  }

  /* normalize weights by number of data */
  scalar_div_dvec(weights, K, num_data);

  /* normalize the covariance matrices by the number of data */
  for (k = 1; k <= K; k++)
    scalar_div_dmat(covars[k], num_dims, num_dims, num_data);

  /* fill in the lower triangle of each covariance matrix based on upper 
     triangle */
  for (k = 1; k <= K; k++)
  {
    covar_k = covars[k];
    for (j=1; j<=num_dims; j++)
      for (j2 = 1; j2 < j; j2++)
        covar_k[j][j2] = covar_k[j2][j];
  }

  return (UT_OK);
}


/*******************************************************************************
BEST_K_MEANS_MIXTURE_ESTIMATE
Performs the k-means algorithm a specified number of times on a dataset, each
time also estimating covariance matrices and mixture weights corresponding to
the resulting mean vectors and k-means cluster memberships, in order to learn
an estimate of a mixture model based on each k-means run.  The likelihood of
the data given the model is computed for each mixture estimate, and the para-
meters of the best mixture estimate chosen by that metric are returned.  The
return value of the function is the corresponding log-likelihood.

Seeds the random number generator before beginning.  Assumes dataset has been
range-normalized.
AG
*******************************************************************************/
double best_k_means_mixture_estimate(double **data, int num_data,
                                    int num_dims, double **means, 
                                    double ***covars, double *weights, 
                                    double **best_means, double ***best_covars,
                                    double *best_weights, double **init_means,
                                    double **probs, double *prob_data, 
                                    double *sum_probs, double num_tries, 
                                    double min_diag, double *min_val, 
                                    double *range, int K, int *clust_size,
                                    int *clust_label, int empty_clust_ok)

{
  int i, k, num_iters, num_empty_clust;
  double log_likelihood, best_log_likelihood;

  if (verbose_output())
    log_printf("%d tries of k-means mixture estimates:\n\n", num_tries);

  best_log_likelihood = 0.0;
  for (i = 1; i <= num_tries; i++)
  {
    if (verbose_output())
      log_printf("--> try # %d of k-means:\n", i);

    /* choose K random points in the data as starting means for this run */
    random_means_from_data(init_means, K, num_dims, data, num_data);

    if (verbose_output())
    {
      log_printf("initial means given to k-means\n");
      print_unnorm_dmatrix(ut_log_fp, K, num_dims, init_means, range, min_val);
    }

    /* perform k-means */
    k_means(data, num_data, num_dims, means, init_means, clust_size, 
            clust_label, K, &num_iters, &num_empty_clust, empty_clust_ok);
    
    if (verbose_output())
    {
      log_printf("means determined by k-means\n");
      print_unnorm_dmatrix(ut_log_fp, K, num_dims, means, range, min_val);
    }

    /* compute covariance matrices and mixture weights based on means and 
       cluster memberships */
    k_means_mixture_estimate(data, means, covars, weights, clust_label,
                             num_data, num_dims, K);

    if (debug_output())
    {
      log_printf("covariance matrices after computing\n");
      for (k = 1; k <= K; k++)
        print_unnorm_cov_dmatrix(ut_log_fp, num_dims, covars[k], range);
      log_printf("class weights after computing\n");
      print_drow(ut_log_fp, K, weights);
    }

    /* given the parameters, compute the likelihood of this model */
    log_likelihood = mixture_likelihood(K, means, covars, weights, data, 
                                        probs, prob_data, sum_probs, num_data, 
                                        num_dims, "rows", min_diag);

    /* save the current best set of parameters */
    if (log_likelihood > best_log_likelihood)
    {
      best_log_likelihood = log_likelihood;
      for (k = 1; k <= K; k++)
      {
        copy_dvec(means[k], best_means[k], num_dims);
        copy_dmat(covars[k], best_covars[k], num_dims, num_dims);
      }
    }
        
    /* report the parameters found */
    if (normal_output())
    {
      log_printf("log-likelihood = %g\n", log_likelihood);
      log_printf("\nk-means mixture estimate parameters\n");
      log_printf(  "-----------------------------------\n\n");

      print_unnorm_gauss_parms_set(ut_log_fp, num_dims, K, means, covars, range,
                                   min_val);
    }
  }  /* end k-means loop over numtries */

  return (best_log_likelihood);
}


/*******************************************************************************
MIXTURE_LIKELIHOOD
Given mixture model parameters, a dataset, and allocated space for storing
all the intermediate probabilities, compute the log-likelihood of the model 
given the data.

K is the number of components in the mixture.
means is the  matrix of mean vectors, where each row is a mean.
covars is the array of covariance matrices.
weights is the vector of class weights.
data is the dataset.
probs is the matrix to contain all intermediate probabilities.
prob_data is the vector to contain the posterior probability of each datum
given the model.
sum_probs is the vector to contain the posterior probability of each class
given the data.
numrows is the number of data.
numcols is the number of attributes.
rows_or_cols is the string specifying whether data are stored as rows or as
columns.
min_diag is a small number for perturbing ill-conditioned covariance matrices.
AG
*******************************************************************************/
double mixture_likelihood(int K, double **means, double ***covars, double *weights,
                         double **data, double **probs, double *prob_data, 
                         double *sum_probs, int numrows, int numcols, 
                         char *rows_or_cols, double min_diag)

{
  int      i,k;
  double    log_likelihood;
  double   *mean_k, **covar_k, *probs_k;

  /* initialize */
  set_dvec(prob_data, numrows, 0.0);
  set_dvec(sum_probs, K, 0.0);
  log_likelihood = 0.0;

  for (k=1; k<=K; k++)
  {
    /* get the parameters corresponding to this class; also get the place
       where the probabilities will be stored */
    mean_k  = means[k];
    covar_k = covars[k];
    probs_k = probs[k];
    
    /* note that the following 3 operations all cycle over the data; they
       should be consolidated for efficiency */

    /* A. first compute the LIKELIHOOD of the data given the parameters for
       this class; store it in the array probs[k] */
    vector_prob_gauss(data, mean_k, covar_k, numcols, numrows, probs_k,
                      min_diag, rows_or_cols);
    
    if (debug_output())
    {
      log_printf("First three density values from vector_prob_gauss\n");
      print_drow(ut_log_fp, 3, probs_k);
      log_printf("%f %f %f\n", probs_k[1],probs_k[2],probs_k[3]);
      log_printf("%g %g %g\n", probs_k[1],probs_k[2],probs_k[3]);
    }

    /* B. multiply by the class prior, or weight, to get the JOINT
       probability of the data and the parameters; store it in probs[k] */
    scalar_mult_dvec(probs_k, numrows, weights[k]);

    /* contribute to the probability of each datum across all models */
    for (i=1; i<=numrows; i++)
      prob_data[i] += probs_k[i];
    /*
    add_dvec(numrows, probs_k, prob_data, prob_data);
    */
  }

  if (debug_output())
  {
    log_printf("First three density values from prob_data\n");
    print_drow(ut_log_fp, 3, prob_data);
  }

  /* C. normalize by the probability of the data across all models to get 
     the POSTERIOR probability of each datum; store it in probs[k] */
  for (k=1; k<=K; k++)
  {
    probs_k = probs[k];
    for (i=1; i<=numrows; i++)
      if (prob_data[i] != 0)
        probs_k[i] /= prob_data[i];
  }
  /*
    div_dvec_elt(numrows, probs_k, prob_data, probs_k);
    */

  if (debug_output())
    for (k=1; k<=K; k++) 
    {
      log_printf("first three posterior probabilities for class %d\n", k);
      log_printf("    in true class A\n");
      print_drow(ut_log_fp, 3, probs[k]);
      log_printf("first three posterior probabilities for class %d\n", k);
      log_printf("    in true class B\n");
      for (i=251; i<=253; i++)
        log_printf("%g ",probs[k][i]);
      log_printf("\n");
    }

  /* compute the SUM OF THE POSTERIOR probabilities for each class, which
     is used in computing the parameters later */
  for (k=1; k<=K; k++)
    sum_probs[k] = sum_dvec(probs[k], numrows);

  if (debug_output()) {
    log_printf("sum of posteriors for all classes\n");
    print_drow(ut_log_fp, K, sum_probs);
  }

  /* compute the log-likelihood of all the data given the mixture model */
  for (i=1; i<=numrows; i++)
    log_likelihood += (double) log((double) (prob_data[i] + DBL_MIN));

  if (debug_output()) {
    log_printf("log-likelihood = %g\n", log_likelihood);
  }

  return(log_likelihood);
}


/*******************************************************************************
ESTIMATE_MIXTURE_PARAMS
Given the matrix containing the posterior probability of each datum belonging
to each class in the mixture model, a dataset, and allocated space for storing
all the parameters to be computed, estimate the mixture model parameters from
the class probabilities.
Assumes the data are stored in columns, not rows.
K is the number of components in the mixture.
means is the  matrix of mean vectors, where each row is a mean.
covars is the array of covariance matrices.
weights is the vector of class weights.
data is the dataset.
probs is the matrix to contain all intermediate probabilities.
sum_probs is the vector to contain the posterior probability of each class
given the data.
numrows is the number of data.
numcols is the number of attributes.
min_weight is the minimum class weight; it determines when a class is consi-
dered to have collapsed.
num_empty_clusters is the number of classes that have collapsed.
empty_cluster_ok is the flag indicating whether to continue in the empty
cluster case.
AG
*******************************************************************************/
int estimate_mixture_params(int K, double **means, double ***covars, 
                            double *weights, double **data, double **probs, 
                            double *sum_probs, int numrows, int numcols, 
                            double min_weight, int *num_empty_clusters, 
                            boolean empty_cluster_ok)

{
  int      i,j,j2,k;
  double    log_likelihood;
  double    *mean_k, **covar_k, *probs_k, *data_j, *data_j2;
  int      empty_clusters;

  /* A. compute the class WEIGHTS, or priors */
  for (k=1; k<=K; k++)
    weights[k] = sum_probs[k] / numrows;
  
  if (debug_output()) {
    log_printf("weights for all classes\n");
    print_drow(ut_log_fp, K, weights);
  }

  /* check for eliminated clusters */
  empty_clusters = 0;
  for (k=1; k<=K; k++)
    if (weights[k] <= min_weight)
      empty_clusters++;

  /* decide what to do if a cluster was eliminated, based on flag */
  if (empty_clusters > 0)
  {
    if (empty_cluster_ok == UT_FALSE)
    {
      /* set some values to be returned */
      *num_empty_clusters = empty_clusters;

      return (UT_ERROR);
    }
    else
    {
      K -= empty_clusters;  /* caller must also do this */
      /* also get rid of the corresponding means in the means matrix 
         and covar. matrices and weights */
      log_printf("This option does not work yet.\n");

      /* set some values to be returned */
      *num_empty_clusters = empty_clusters;

      return (UT_ERROR);  /* for now... when this works, just continue
                             with smaller K */
    }
  }
      
  /* B. compute the class MEAN vectors */
  /* can i just use mat_mult() here ???  I think so. */
  for (k=1; k<=K; k++)
  {
    mean_k  = means[k];
    probs_k = probs[k];
      
    /* first multiply the posteriors matrix by the data matrix */
    for (j=1; j<=numcols; j++)
    {
      data_j = data[j];

      mean_k[j] = 0.0;
      for (i=1; i<=numrows; i++)
        mean_k[j] += probs_k[i] * data_j[i];
    }

    /* then normalize the means by the sum of the posteriors */
    scalar_div_dvec(mean_k, numcols, sum_probs[k]);
  }

  /*
  if (debug_output()) {
    log_printf("means for all classes\n");
    print_unnorm_dmatrix(ut_log_fp, K, numcols, means, range, minval);
  }
  */

  /* C. compute the class COVARIANCE matrices */
  for (k=1; k<=K; k++)
  {
    covar_k = covars[k];
    mean_k  = means[k];
    probs_k = probs[k];

    /* multiply the deviation vectors by the posteriors;
       also normalize by the sum of the posteriors */
    for (j=1; j<=numcols; j++)
    {
      data_j = data[j];
      for (j2=j; j2<=numcols; j2++)
      {
        data_j2 = data[j2];
        
        covar_k[j][j2] = 0.0;
        for (i=1; i<=numrows; i++)
          covar_k[j][j2] += (data_j[i] - mean_k[j]) * (data_j2[i] - mean_k[j2])
                            * probs_k[i];

        covar_k[j][j2] /= sum_probs[k];
      }
    }

    /* fill in the lower triangle of the matrix based on upper triangle */
    for (j=1; j<=numcols; j++)
      for (j2=1; j2<j; j2++)
        covar_k[j][j2] = covar_k[j2][j];

    /*
    if (debug_output()) {
      log_printf("cov. matrix for class %d\n", k);
      print_unnorm_cov_dmatrix(ut_log_fp, numcols, numcols, covar_k, range);
    }
    */

  }

  /* set some values to be returned */
  *num_empty_clusters = empty_clusters;

  return (UT_OK);
}

/*******************************************************************************
RANDOM_MEANS
Randomly set the mean vectors, drawing from 0-1 uniform distribution for each
attribute value.  The matrix means contains the resulting row vectors.
AG
*******************************************************************************/
int random_means(double **means, int num_means, int num_cols)

{

  return ( random_dmatrix(means, num_means, num_cols) );
}


/*******************************************************************************
RANDOM_MEANS_FROM_DATA
Randomly set the mean vectors by drawing randomly from the row vectors of a
dataset.  The matrix means contains the resulting row vectors.
AG
*******************************************************************************/
int random_means_from_data(double **means, int num_means, int num_cols, 
                               double **data, int num_data)

{
  boolean *selected;

  /* allocate temporary storage */
  selected = (boolean*) NR_cvector(1, num_data);

  random_select_drows(data, num_data, num_cols, means, num_means, selected);

  NR_free_cvector(selected, 1, num_data);
  return (UT_OK);
}

