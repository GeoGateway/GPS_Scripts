/*******************************************************************************
MODULE HEADER:
da_cluster.h
*******************************************************************************/

#ifndef _DA_CLUSTER_H_
#define _DA_CLUSTER_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_cluster_h_rcsid[] = "$Id: da_cluster.h,v 1.10 1998/03/09 00:40:28 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_cluster.h,v $
 * Revision 1.10  1998/03/09 00:40:28  granat
 * fixed boolean bug
 *
 * Revision 1.9  1997/01/29 21:29:20  agray
 * new format.
 *
 * Revision 1.8  1996/10/31 02:05:10  agray
 * renamed from "da_clust" to "da_cluster";
 * added functions from HMM project;
 * changed .h and .c formats throughout library;
 * some reorganizing between modules.
 *
 * Revision 1.7  1996/10/30 20:30:06  agray
 * changed naming from "mixmod" to "mixture", other naming changes.
 *
 * Revision 1.6  1996/09/27 17:57:57  agray
 * fixed log message.
 *
 * Revision 1.5  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.4  1996/08/28 19:41:52  agray
 * changed arguments to random_means_from_data().
 *
 * Revision 1.3  1996/07/11 16:41:47  agray
 * changed randomize_means() to random_means();
 * added random_means_from_data(), mixmod_likelihood(), estimate_mixmod_
 * params().
 * .
 *
 * Revision 1.2  1996/05/14 01:02:45  agray
 * minor compile problems.
 *
 * Revision 1.1  1996/05/07 20:47:19  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* k-means clustering */

int k_means(double **data, int numrows, int numcols, double **means, 
            double **initmeans, int *clustsize, int *class, int K,
            int *numiters, int *num_empty_clusters, boolean empty_cluster_ok);
int k_means_mixture_estimate(double **data, double **means, double ***covars, 
                             double *weights, int *clust_label, 
                             int num_data, int num_dims, int K);
double best_k_means_mixture_estimate(double **data, int num_data,
                                    int num_dims, double **means, 
                                    double ***covars, double *weights, 
                                    double **best_means, double ***best_covars,
                                    double *best_weights, double **init_means,
                                    double **probs, double *prob_data, 
                                    double *sum_probs, double num_tries, 
                                    double min_diag, double *min_val, 
                                    double *range, int K, int *clust_size,
                                    int *clust_label, int empty_clust_ok);

/* EM algorithm for mixtures */

double mixture_likelihood(int K, double **means, double ***covars, double *weights,
                         double **data, double **probs, double *prob_data, 
                         double *sum_probs, int numrows, int numcols, 
                         char *rows_or_cols, double min_diag);
int estimate_mixture_params(int K, double **means, double ***covars, 
                            double *weights, double **data, double **probs, 
                            double *sum_probs, int numrows, int numcols, 
                            double min_weight, int *num_empty_clusters, 
                            boolean empty_cluster_ok);

/* parameter initialization */

int random_means(double **means, int num_means, int num_cols);
int random_means_from_data(double **means, int num_means, int num_cols, 
                           double **data, int num_data);

#endif /* _DA_CLUSTER_H_ */

