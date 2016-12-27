/*******************************************************************************
MODULE HEADER:
da_cluster_pp.h
*******************************************************************************/

#ifndef _DA_CLUSTER_PP_H_
#define _DA_CLUSTER_PP_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_cluster_pp_h_rcsid[] = "$Id: da_cluster_pp.h,v 1.6 1998/03/09 00:40:54 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_cluster_pp.h,v $
 * Revision 1.6  1998/03/09 00:40:54  granat
 * fixed boolean bug
 *
 * Revision 1.5  1997/01/29 21:29:37  agray
 * new format.
 *
 * Revision 1.4  1996/10/31 02:06:43  agray
 * renamed from "da_clust_pp" to "da_cluster_pp";
 * changed .h and .c formats throughout library;
 *
 * Revision 1.3  1996/10/31 00:23:46  agray
 * changed naming from "mixmod" to "mixture"; other naming changes.
 *
 * Revision 1.2  1996/07/19 17:53:23  agray
 * added random_means_pp(), started random_means_from_data_pp()
 *
 * Revision 1.1  1996/07/11 16:44:01  agray
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

int k_means_pp(double **data, int numrows, int numcols, double **means, 
              double **initmeans, int *clustsize, int *class, int K,
              int *numiters, int *num_empty_clusters, boolean empty_cluster_ok,
              int pe, int numPE);

/* EM algorithm for mixtures */

double mixture_likelihood_pp(int K, double **means, double ***covars, 
                           double *weights, double **data, double **probs, 
                           double *prob_data, double *sum_probs, int numrows, 
                           int numcols, char *rows_or_cols, double min_diag,
                           int pe, int numPE);
int estimate_mixture_params_pp(int K, double **means, double ***covars, 
                              double *weights, double **data, double **probs, 
                              double *sum_probs, int total_numrows, int numrows, 
                              int numcols, double min_weight, 
                              int *num_empty_clusters, boolean empty_cluster_ok,
                              int pe, int numPE);

/* variations of EM algorithm for mixture of gaussians */
double mixture_likelihood_anneal_pp(int K, double **means, double ***covars, 
		                    double *weights, double **data, 
				    double **probs, double *prob_data, 
				    double *sum_probs, int numrows, 
				    int numcols, char *rows_or_cols, 
				    double min_diag, double temperature, 
				    int pe, int numPE);

/* parameter initialization */

int random_means_pp(double **means, int num_means, int num_cols, int pe);


#endif /* _DA_CLUSTER_PP_H_ */
