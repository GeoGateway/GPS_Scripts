/*******************************************************************************
MODULE NAME
da_cluster_pp

ONE-LINE SYNOPSIS
Parallel functions related to clustering.

SCOPE OF THIS MODULE
Analogous to da_cluster in scope, except that functions in this module are
written to run on a parallel machine.  

SEE ALSO
Most or all of the functions in this module are mirrors of their serial versions
in da_cluster.

REFERENCE(S)
-

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_cluster_pp.c,v 1.15 2000/02/24 18:50:54 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_cluster_pp.c,v $
 * Revision 1.15  2000/02/24 18:50:54  granat
 * changed communications to use MPI_Allreduce
 * this increases communication speed
 *
 * Revision 1.14  1998/06/30 17:23:23  agray
 * changed includes.
 *
 * Revision 1.13  1998/04/21 17:07:35  roden
 * Took out #include <pvm3.h> because I think it's not needed here, and it
 * doesn't conform to da_platform.h control over presence/absence of pvm.
 *
 * Revision 1.12  1998/03/09 00:38:36  granat
 * fixed boolean bug
 *
 * Revision 1.11  1997/06/05 18:55:46  granat
 * edited to conform to changes in da_linalg
 *
 * Revision 1.10  1997/06/02 15:35:17  granat
 * changed to use new NR naming convention
 *
 * Revision 1.9  1997/04/05 19:16:58  granat
 * made adjustments to account for changes in da_linalg
 *
 * Revision 1.8  1997/01/29 21:44:05  agray
 * new formatting, cleaning up debugging output using ut_output.
 *
 * Revision 1.7  1996/10/31 02:06:43  agray
 * renamed from "da_clust_pp" to "da_cluster_pp";
 * changed .h and .c formats throughout library;
 * 
 * Revision 1.6  1996/10/31 00:23:18  agray
 * changed naming from "mixmod" to "mixture"; other naming changes.
 * 
 * Revision 1.5  1996/09/23 23:57:23  agray
 * updated utRandomPick() to use char's for selected array.
 * 
 * Revision 1.4  1996/09/23 22:40:46  agray
 * minor changes.
 * 
 * Revision 1.3  1996/07/19 17:53:03  agray
 * added random_means_pp(), started random_means_from_data_pp()
 * 
 * Revision 1.2  1996/07/17 20:42:25  agray
 * cosmetic.
 * 
 * Revision 1.1  1996/07/11 16:42:35  agray
 * Initial revision
 * 
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* MPI library */
#include "mpi.h"

/* UT library */
#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_cluster.h"
#include "da_comm_pp.h"
#include "da_linalg.h"
#include "da_random.h"
#include "da_util.h"
#include "da_io.h"
#include "da_probstat.h"

/* this module's header */
#include "da_cluster_pp.h"


/*******************************************************************************
K_MEANS_PP
Parallel version of k_means().
Perform the k-means clustering algorithm on a dataset of continuous values,
returning the matrix of k mean vectors and a class labelling for each of the
data.  The input matrices and vectors are all assumed to have been allocated
before calling this function.  The argument initmeans should contain the
mean vectors which will start off the k-means iterations.  The data is assumed
to be stored in row vector form (i.e. each row is a datum).
When the function is finished, each of the processors has the final means.
Each processor will end up with its own chunk of the labels file.

NOTE: In practice, it is a good idea to normalize the columns of the data to
the same scale before performing k-means (see the function normalize_data()).
AG
*******************************************************************************/
int k_means_pp(double **data, int numrows, int numcols, double **means, 
               double **initmeans, int *clustsize, int *class, int K,
               int *numiters, int *num_empty_clusters, boolean empty_cluster_ok,
               int pe, int numPE)

{
    double          **oldmeans, **tempmeans, **othermeans;
    int            *otherclustsize;
    int            i, j, k, p, cc;
    int            bestclass;
    double          dist, mindist;
    int            stable, otherstable;
    int            empty_clusters, iters;

    if (verbose_output() && (pe == 0))
      log_printf("%d: beginning k_means_pp\n",pe);

    /* make storage related to message passing */
    otherclustsize = NR_ivector(1,K);
    othermeans = NR_dmatrix(1,K,1,numcols);

    /* initialize before iterating */
    stable = UT_FALSE;
    empty_clusters = 0;
    iters = 0;
    fast_set_ivec(class, numrows, 1);

    /* create storage for temporary information */
    oldmeans = initmeans;

    if (debug_output() && (pe == 0))
      log_printf("%d: beginning k_means_pp iterations\n",pe);

    /* perform k-means iterations */
    while (stable == UT_FALSE)
    {
      if (debug_output() && (pe == 0))
        log_printf("%d: iteration number %d\n",pe,iters);

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
            log_printf("item %d to mean %d: %f\n",
                   i,1,euclid_dist(data[i],oldmeans[1],numcols));
             */
        for (k=2; k<=K; k++)
        {
          /*
            log_printf("item %d to mean %d: %f\n",
                   i,k,euclid_dist(data[i],oldmeans[k],numcols));
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
      if (verbose_output() && (pe == 0))
        log_printf("iter# %d\n", iters);

      if (debug_output() && (pe == 0))
        log_printf("%d: before parallel comm. to collate clustsize\n", pe);

      /* Do parallel communication to collate cluster sizes */
      MPI_Allreduce(&clustsize[1], &otherclustsize[1], K, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      copy_ivec(otherclustsize, clustsize, K);

/*
      DA_barrier();
      if (pe != 0) 
      {
        cc = DA_send_msg(&clustsize[1],K,0,pe,DA_INT);

        if (debug_output() && (pe == 0))
          log_printf("%d: after sending local clustsize\n",pe);

        cc = DA_broadcast_msg(&clustsize[1],K,DA_ALL_PES,0,DA_INT);

        if (debug_output() && (pe == 0))
          log_printf("%d: after receiving global clustsize\n",pe);
      }
      else 
      {
        for(p=1;p<numPE;p++)
        {
          cc = DA_recv_msg(&otherclustsize[1],K,p,p,DA_INT);

          if (debug_output() && (pe == 0))
            log_printf("%d: after receiving local clustsize\n",pe);

          for(k=1;k<=K;k++)
            clustsize[k] += otherclustsize[k];
        }
        cc = DA_broadcast_msg(&clustsize[1],K,DA_ALL_PES,0,DA_INT);

        if (debug_output() && (pe == 0))
          log_printf("%d: after sending global clustsize\n",pe);
      }    
*/

      if (debug_output() && (pe == 0))
        log_printf("%d: after parallel comm. to collate clustsize\n",pe);

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

          /* free temporary storage */
          NR_free_ivector(otherclustsize,1,K);
          NR_free_dmatrix(othermeans,1,K,1,numcols);

          return (UT_ERROR);
        }
        else
        {
          /* K -= empty_clusters; */ /* caller must also do this */
          /* also get rid of the corresponding means in the means matrix */
	  /*
          log_printf("This option does not work yet.\n");
	  */

          /* set some values to be returned */
          *num_empty_clusters = empty_clusters;
          *numiters = iters;

          /* free temporary storage */
          NR_free_ivector(otherclustsize,1,K);
          NR_free_dmatrix(othermeans,1,K,1,numcols);

          return (UT_ERROR);  /* for now... when this works, just continue
                                 with smaller K */
        }
      }

      if (debug_output() && (pe == 0))
        log_printf("%d: after checking for empty clusters in k_means_pp\n",pe);

      /* normalize the local means by the global cluster sizes */
      for (k=1; k<=K; k++)
        for (j=1; j<=numcols; j++) {
          means[k][j] /= (double) clustsize[k];
	  /*
	  printf("%d: means[%d][%d] = %f\n",pe,k,j,means[k][j]);
	  fflush(stdout);
	  */
	}

      if (debug_output() && (pe == 0))
        log_printf("%d: after normalizing by clustsizes in k_means_pp\n",pe);

      if (debug_output() && (pe == 0))
        log_printf("%d: before doing parallel comm. to collate means\n",pe);

      /* Do parallel communication to collate means */
      MPI_Allreduce(&means[1][1], &othermeans[1][1], K * numcols, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      copy_dmat(othermeans, means, K, numcols);

      /* Do parallel communication to calculate global stability */
      MPI_Allreduce(&stable, &otherstable, 1, MPI_INT, MPI_LAND, 
                    MPI_COMM_WORLD);
      stable = otherstable;

/*
      DA_barrier();
      if (pe != 0)
      {
        cc = DA_send_msg(&means[1][1],numcols*K,0,pe,DA_DOUBLE);
        cc = DA_broadcast_msg(&means[1][1],numcols*K,DA_ALL_PES,0,DA_DOUBLE);
        cc = DA_send_msg(&stable,1,0,pe,DA_INT);
        cc = DA_broadcast_msg(&stable,1,DA_ALL_PES,0,DA_INT);
      }
      else
      {
        for(p=1;p<numPE;p++)
        {
          cc = DA_recv_msg(&othermeans[1][1],numcols*K,p,p,DA_DOUBLE);
  
          for(k=1;k<=K;k++)
            for(j=1;j<=numcols;j++)
              means[k][j] += othermeans[k][j];
        }

        cc = DA_broadcast_msg(&means[1][1],numcols*K,DA_ALL_PES,0,DA_DOUBLE);

        for (p=1;p<numPE;p++)
        {
          cc = DA_recv_msg(&otherstable,1,p,p,DA_INT);

          if ((stable != UT_TRUE) || (otherstable != UT_TRUE))
            stable = UT_FALSE;
        }

        cc = DA_broadcast_msg(&stable,1,DA_ALL_PES,0,DA_INT);
      }
*/

      if (debug_output() && (pe == 0))
        log_printf("%d: after doing parallel comm. to collate means\n",pe);

      /* make the new means into the oldmeans for next iteration */
      tempmeans = oldmeans; oldmeans = means; means = tempmeans;


      if (debug_output() && (pe == 0))
        print_dmatrix(ut_log_fp, K, numcols, means);
    }

    if (debug_output() && (pe == 0))
      log_printf("%d: after convergence loop in k_means_pp()\n",pe);

    if (debug_output() && (pe == 0))
      log_printf("%d: empty_clusters = %d in k_means_pp()\n",pe,empty_clusters);

    /* set some values to be returned */
    *num_empty_clusters = empty_clusters;
    *numiters = iters;

    if (verbose_output() && (pe == 0))
      log_printf("%d: num_empty_clusters = %d in k_means_pp()\n",pe,
                 *num_empty_clusters);

    /* free temporary storage */
    NR_free_ivector(otherclustsize,1,K);
    NR_free_dmatrix(othermeans,1,K,1,numcols);

    return (UT_OK);
}


/*******************************************************************************
MIXTURE_LIKELIHOOD_PP
Parallel version of mixture_likelihood().
Given mixture model parameters, a dataset, and allocated space for storing
all the intermediate probabilities, compute the log-likelihood of the model 
given the data.
 
When the function is finished, each processor has the global log-likelihood 
value.
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
pe is the current processing element's id.
numPE is the total number of processing elements.
AG
*******************************************************************************/
double mixture_likelihood_pp(int K, double **means, double ***covars, 
                           double *weights, double **data, double **probs, 
                           double *prob_data, double *sum_probs, int numrows, 
                           int numcols, char *rows_or_cols, double min_diag,
                           int pe, int numPE)

{
  int      i,k,p,cc;
  double    log_likelihood, other_log_likelihood;
  double   *mean_k, **covar_k, *probs_k;
  double   *other_sum_probs;

  double  a,b,c,d,e,f;

  /* make temporary storage related to message passing */
  other_sum_probs = NR_dvector(1,K);

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
    
    if (debug_output() && (pe == 0))
    {
      log_printf("First three density values from vector_gauss_eval\n");
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
    /* Faster than add_dvec(numrows, probs_k, prob_data, prob_data); */
  }

  if (debug_output() && (pe == 0)) 
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
  /* Faster than div_dvec_elt(numrows, probs_k, prob_data, probs_k); */

  if (debug_output() && (pe == 0))
  {
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
  }

  /* compute the SUM OF THE POSTERIOR probabilities for each class, which
     is used in computing the parameters later */
  for (k=1; k<=K; k++)
    sum_probs[k] = sum_dvec(probs[k], numrows);

  if (debug_output() && (pe == 0)) 
  {
    log_printf("sum of posteriors for all classes\n");
    print_drow(ut_log_fp, K, sum_probs);
    log_printf("first 10 cols of probs\n");
    print_dmatrix(ut_log_fp,K,10,probs);
    a=0; b=0;
    for (i=1;i<=numrows;i++)
      a += probs[1][i];
    for (i=1;i<=numrows;i++)
      b += probs[2][i];
    log_printf("sum of 1st row = %g, sum of 2nd row = %g\n",a,b);
  }

  /* compute the log-likelihood of all the data given the mixture model */
  for (i=1; i<=numrows; i++)
    log_likelihood += (double) log((double) (prob_data[i] + DBL_MIN));

  if (debug_output() && (pe == 0))
    log_printf("log-likelihood = %g\n", log_likelihood);

  /* Do parallel communication to collate sums of posteriors and the global
     log-likelihood */

  MPI_Allreduce(&sum_probs[1], &other_sum_probs[1], K, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  copy_dvec(other_sum_probs, sum_probs, K);

  MPI_Allreduce(&log_likelihood, &other_log_likelihood, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  log_likelihood = other_log_likelihood;

/*
  DA_barrier();
  if (pe != 0)
  {
    cc = DA_send_msg(&sum_probs[1],K,0,pe,DA_DOUBLE);
    cc = DA_broadcast_msg(&sum_probs[1],K,DA_ALL_PES,0,DA_DOUBLE);
    
    cc = DA_send_msg(&log_likelihood,1,0,pe,DA_DOUBLE);
    cc = DA_broadcast_msg(&log_likelihood,1,DA_ALL_PES,0,DA_DOUBLE);
  }
  else 
  {
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&other_sum_probs[1],K,p,p,DA_DOUBLE);
      add_dvec(sum_probs, other_sum_probs, sum_probs, K);
    }
    cc = DA_broadcast_msg(&sum_probs[1],K,DA_ALL_PES,0,DA_DOUBLE);
    
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&other_log_likelihood,1,p,p,DA_DOUBLE);
      log_likelihood += other_log_likelihood;
    }
    cc = DA_broadcast_msg(&log_likelihood,1,DA_ALL_PES,0,DA_DOUBLE);
  }     
*/

  if (debug_output() && (pe == 0))
  {
    log_printf("sum_probs in mixture_likelihood:\n");
    print_drow(ut_log_fp, K, sum_probs);
  }
  
  /* free temporary storage */
  NR_free_dvector(other_sum_probs,1,K);

  return(log_likelihood);
}


/*******************************************************************************
MIXTURE_LIKELIHOOD_ANNEAL_PP
Parallel version of mixture_likelihood(), with deterministic annealing.
Given mixture model parameters, a dataset, and allocated space for storing
all the intermediate probabilities, compute the log-likelihood of the model 
given the data.
 
When the function is finished, each processor has the global log-likelihood 
value.
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
pe is the current processing element's id.
numPE is the total number of processing elements.
RG
*******************************************************************************/
double mixture_likelihood_anneal_pp(int K, double **means, double ***covars, 
                                    double *weights, double **data, 
				    double **probs, double *prob_data, 
				    double *sum_probs, int numrows, 
				    int numcols, char *rows_or_cols, 
				    double min_diag, double temperature,
				    int pe, int numPE)

{
  int       i,k,p,cc;
  double    log_likelihood, other_log_likelihood;
  double   *mean_k, **covar_k, *probs_k;
  double   *other_sum_probs;

  double  a,b,c,d,e,f;

  /* make temporary storage related to message passing */
  other_sum_probs = NR_dvector(1,K);

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
    
    if (debug_output() && (pe == 0))
    {
      log_printf("First three density values from vector_gauss_eval\n");
      print_drow(ut_log_fp, 3, probs_k);
      log_printf("%f %f %f\n", probs_k[1],probs_k[2],probs_k[3]);
      log_printf("%g %g %g\n", probs_k[1],probs_k[2],probs_k[3]);
    }

    /* B. multiply by the class prior, or weight, to get the JOINT
       probability of the data and the parameters; store it in probs[k] */
    scalar_mult_dvec(probs_k, numrows, weights[k]);

    /* raise each joint probability to the power of the temperature
     * parameter to do the deterministic annealing version
     */
    for (i = 1; i <= numrows; i++)
      probs_k[i] = pow(probs_k[i], temperature);

    /* contribute to the probability of each datum across all models */
    for (i=1; i<=numrows; i++)
      prob_data[i] += probs_k[i];
    /* Faster than add_dvec(numrows, probs_k, prob_data, prob_data); */
  }

  if (debug_output() && (pe == 0)) 
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
  /* Faster than div_dvec_elt(numrows, probs_k, prob_data, probs_k); */

  if (debug_output() && (pe == 0))
  {
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
  }

  /* compute the SUM OF THE POSTERIOR probabilities for each class, which
     is used in computing the parameters later */
  for (k=1; k<=K; k++)
    sum_probs[k] = sum_dvec(probs[k], numrows);

  if (debug_output() && (pe == 0)) 
  {
    log_printf("sum of posteriors for all classes\n");
    print_drow(ut_log_fp, K, sum_probs);
    log_printf("first 10 cols of probs\n");
    print_dmatrix(ut_log_fp,K,10,probs);
    a=0; b=0;
    for (i=1;i<=numrows;i++)
      a += probs[1][i];
    for (i=1;i<=numrows;i++)
      b += probs[2][i];
    log_printf("sum of 1st row = %g, sum of 2nd row = %g\n",a,b);
  }

  /* compute the log-likelihood of all the data given the mixture model */
  for (i=1; i<=numrows; i++)
    log_likelihood += (double) log((double) (prob_data[i] + DBL_MIN));

  if (debug_output() && (pe == 0))
    log_printf("log-likelihood = %g\n", log_likelihood);

  /* Do parallel communication to collate sums of posteriors and the global
     log-likelihood */

  MPI_Allreduce(&sum_probs[1], &other_sum_probs[1], K, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  copy_dvec(other_sum_probs, sum_probs, K);

  MPI_Allreduce(&log_likelihood, &other_log_likelihood, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  log_likelihood = other_log_likelihood;

  if (debug_output() && (pe == 0))
  {
    log_printf("sum_probs in mixture_likelihood:\n");
    print_drow(ut_log_fp, K, sum_probs);
  }
  
  /* free temporary storage */
  NR_free_dvector(other_sum_probs,1,K);

  return(log_likelihood);
}


/*******************************************************************************
ESTIMATE_MIXTURE_PARAMS_PP
Parallel version of estimate_mixture_params().
Given the matrix containing the posterior probability of each datum belonging
to each class in the mixture model, a dataset, and allocated space for storing
all the parameters to be computed, estimate the mixture model parameters from
the class probabilities.
Assumes the data are stored in columns, not rows.
When the function is finished, each processor has the global values for the
parameters.
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
pe is the current processing element's id.
numPE is the total number of processing elements.
AG
*******************************************************************************/
int estimate_mixture_params_pp(int K, double **means, double ***covars, 
                              double *weights, double **data, double **probs, 
                              double *sum_probs, int total_numrows, int numrows, 
                              int numcols, double min_weight, 
                              int *num_empty_clusters, boolean empty_cluster_ok,
                              int pe, int numPE)

{
  int      i,j,j2,k,p,cc;
  double    log_likelihood;
  double    *mean_k, **covar_k, *probs_k, *data_j, *data_j2;
  double    **othermeans, **other_covar_k;
  int      empty_clusters;

  /* make storage related to message passing */
  othermeans = NR_dmatrix(1,K,1,numcols);
  other_covar_k = NR_dmatrix(1,numcols,1,numcols);

  /* A. compute the class WEIGHTS, or priors */
  for (k=1; k<=K; k++)
    weights[k] = sum_probs[k] / total_numrows;
  
  if (debug_output() && (pe == 0))
  {
    log_printf("weights for all classes\n");
    print_drow(ut_log_fp, K, weights);
  }

  if (debug_output() && (pe == 0))
  {
    log_printf("%d: after computing weights\n",pe);
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

      /* clean up allocated memory before exiting */
      NR_free_dmatrix(othermeans,1,K,1,numcols);
      NR_free_dmatrix(other_covar_k,1,numcols,1,numcols);

      return (UT_ERROR);
    }
    else
    {
      /* K -= empty_clusters; */ /* caller must also do this */
      /* also get rid of the corresponding means in the means matrix 
         and covar. matrices and weights */
      /*
      log_printf("This option does not work yet.\n");
      */

      /* set some values to be returned */
      *num_empty_clusters = empty_clusters;

      /* clean up allocated memory before exiting */
      NR_free_dmatrix(othermeans,1,K,1,numcols);
      NR_free_dmatrix(other_covar_k,1,numcols,1,numcols);

      return (UT_ERROR);  /* for now... when this works, just continue
                             with smaller K */
    }
  }
      
  if (debug_output() && (pe == 0))
    log_printf("%d: after checking for eliminated clusters\n",pe);

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

  if (debug_output() && (pe == 0))
    log_printf("%d: after computing means\n",pe);

  if (debug_output() && (pe == 0))
    log_printf("%d: before parallel comm. to collate means\n",pe);

  /* Do parallel communication to collate means */
  MPI_Allreduce(&means[1][1], &othermeans[1][1], K * numcols, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  copy_dmat(othermeans, means, K, numcols);

/*
  DA_barrier();
  if (pe != 0)
  {
    cc = DA_send_msg(&means[1][1],numcols*K,0,pe,DA_DOUBLE);
    
    if (debug_output() && (pe == 0))
      log_printf("%d: after sending local means\n",pe);

    cc = DA_broadcast_msg(&means[1][1],numcols*K,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output() && (pe == 0))
      log_printf("%d: after receiving global means\n",pe);

  }
  else 
  {
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&othermeans[1][1],numcols*K,p,p,DA_DOUBLE);
      add_dmat(means, othermeans, means, K, numcols);
    }

    if (debug_output() && (pe == 0))
      log_printf("%d: after receiving local means\n",pe);

    cc = DA_broadcast_msg(&means[1][1],numcols*K,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output() && (pe == 0))
      log_printf("%d: after broadcasting global means\n",pe);
  }
*/

  if (debug_output() && (pe == 0))
    log_printf("%d: after parallel comm. to collate means\n",pe);

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

  }

  if (debug_output() && (pe == 0))
    log_printf("%d: after computing cov. matrices\n",pe);

  if (debug_output() && (pe == 0))
    log_printf("%d: before parallel comm. to collate cov. matrices\n",pe);

  /* Do parallel communication to collate covariance matrices */
  for (k = 1; k <= K; k++) {
    covar_k = covars[k];
    MPI_Allreduce(&covar_k[1][1], &other_covar_k[1][1], numcols * numcols,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    copy_dmat(other_covar_k, covar_k, numcols, numcols);
  }

/*
  DA_barrier();
  if (pe != 0)
  {
    for (k=1; k<=K; k++)
    {
      covar_k = covars[k];
      cc = DA_send_msg(&covar_k[1][1],numcols*numcols,0,k,DA_DOUBLE);
    }

    for (k=1; k<=K; k++)
    {
      covar_k = covars[k];
      cc = DA_broadcast_msg(&covar_k[1][1],numcols*numcols,DA_ALL_PES,0,
		            DA_DOUBLE);
      DA_barrier();
    }
  }
  else 
  {
    for(p=1;p<numPE;p++)
    {
      for (k=1; k<=K; k++)
      {
        covar_k = covars[k];
        cc = DA_recv_msg(&other_covar_k[1][1],numcols*numcols,p,k,DA_DOUBLE);
        add_dmat(covar_k, other_covar_k, covar_k, numcols, numcols);
      }
    }

    if (debug_output() && (pe == 0))
      log_printf("%d: after receiving local cov. matrices\n",pe);

    for (k=1; k<=K; k++)
    {
      covar_k = covars[k];
      cc = DA_broadcast_msg(&covar_k[1][1],numcols*numcols,DA_ALL_PES,0,
                         DA_DOUBLE);
      DA_barrier();
    }

    if (debug_output() && (pe == 0))
      log_printf("%d: after broadcasting global cov. matrices\n",pe);

  }     
*/

  if (debug_output() && (pe == 0))
    log_printf("%d: after parallel comm. to collate cov. matrices\n",pe);

  /* free temporary storage space */
  NR_free_dmatrix(othermeans,1,K,1,numcols);
  NR_free_dmatrix(other_covar_k,1,numcols,1,numcols);

  /* set some values to be returned */
  *num_empty_clusters = empty_clusters;

  return (UT_OK);
}


/*******************************************************************************
RANDOM_MEANS_PP
Parallel version of random_means().
AG
*******************************************************************************/
int random_means_pp(means, num_means, num_cols, pe)

    double **means;
    int   num_means, num_cols;
    int   pe;
{
    int cc;

    if (debug_output() && (pe == 0))
      log_printf("%d: before message passing for choosing random means\n",pe);

    DA_barrier();
    if (pe != 0) /* This is not the controlling process */
    {	
      cc = DA_broadcast_msg(&means[1][1],num_means*num_cols,DA_ALL_PES,0,DA_DOUBLE);
    }
    else
    {
      /* choose initial values for means arbitrarily */
      random_means(means, num_means, num_cols);
      
      cc = DA_broadcast_msg(&means[1][1],num_means*num_cols,DA_ALL_PES,0,DA_DOUBLE);
    }

    if (debug_output() && (pe == 0))
    {
      log_printf("%d: means: \n",pe);
      print_dmatrix(ut_log_fp, num_means, num_cols, means);
    }

    return (UT_OK);
}


/*******************************************************************************
RANDOM_MEANS_FROM_DATA_PP
Parallel version of random_means_from_data().
AG
*******************************************************************************/
int random_means_from_data_pp(double **means, int num_means, int num_cols, 
		              double **data, int num_data, int pe, int numPE)
{
  int cc;
  int i, p;
  int num_local_means;
  int total_num_data;
  int startrow;
  int startmean;
  boolean *selected;

  selected = (boolean *) NR_cvector(1, num_data);

  split_data_pp(&startmean, &num_local_means, num_means, pe, numPE);
  
  if (num_local_means > 0)
    random_select_drows(data, num_data, num_cols, means, num_local_means, 
  		        selected);
  
  if (pe != 0) {
    if (num_local_means > 0)
      cc = DA_send_msg(&means[1][1], num_local_means * num_cols, 0, pe, 
  		       DA_DOUBLE);
    
    cc = DA_broadcast_msg(&means[1][1], num_means * num_cols, DA_ALL_PES, 0,
		          DA_DOUBLE);
  }
  else {
    for (p = 1; p < numPE; p++) {
      split_data_pp(&startmean, &num_local_means, num_means, p, numPE);
      
      if (num_local_means > 0)
        cc = DA_recv_msg(&means[startmean][1], num_local_means * num_cols, p, p,
  		         DA_DOUBLE);
    }
    
    cc = DA_broadcast_msg(&means[1][1], num_means * num_cols, DA_ALL_PES, 0,
		          DA_DOUBLE);
  }
  
  NR_free_cvector(selected, 1, num_data);

  return(UT_OK);
}
