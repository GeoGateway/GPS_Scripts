/*******************************************************************************
MODULE HEADER:
rdahmm_gauss_output.h
*******************************************************************************/
#ifndef _RDAHMM_GAUSS_OUTPUT_H_
#define _RDAHMM_GAUSS_OUTPUT_H_ 1

#ifndef lint
static char rdahmm_gauss_output_hdr_rcsid[] = "$Id: rdahmm_gauss_output.h,v 1.2 2003/08/27 18:05:47 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_gauss_output.h,v $
 * Revision 1.2  2003/08/27 18:05:47  granat
 * modified function prototypes
 *
 * Revision 1.1  2003/08/25 20:10:31  granat
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
int read_gauss_output_params(char *B_file, double **mu, double ***sigma,
                             int N, int D);

int write_gauss_output_params(char *B_file, double **mu, double ***sigma,
                              int N, int D);

int init_gauss_output_hmm(int init_type, double **cont_data, double *pi,
                          double **A, double **mu, double ***sigma, int N,
                          int D, int T); 

int kmeans_init_hmm(double **cont_data, double *pi, double **A, double **mu,
                    double ***sigma, int N, int D, int T);

double calc_regularized_Q_function_gauss(double Qfunc, double *pi, double **A,
                                         double **mu, double ***inv_chol_sigma,
                                         double *sum_T_tau_it,
                                         double omega_Q1, double omega_Q2,
                                         double omega_Q3, double omega_sigma,
                                         int reg_type, int N, int D);

int prob_data_gauss(double **data, double **mu, double ***inv_chol_sigma,
                    double *sqrt_det_sigma, double **prob_data, int N,
                    int T, int D);

int update_gauss_means(double **data, double **tau_it, double *sum_T_tau_it,
                       double **mu, int N, int T, int D);

int update_gauss_covars(double **data, double **tau_it, double *sum_T_tau_it,   
                        double **mu, double ***sigma, int N, int T, int D);

int update_block_gauss_covars(double **data, double **tau_it,
                              double *sum_T_tau_it, double **mu,
                              double ***sigma, int *blocks,
                              int N, int T, int D);

int update_graph_gauss_covars(double **data, double **tau_it,
                              double *sum_T_tau_it, double **mu,
                              double ***sigma, int **graph,
                              int N, int T, int D);

int update_gauss_covars_inner(double *sigma_id, double *diff, double tau_it_it,
                              double diff_d, int d);

int update_gauss_euclid(double **data, double **tau_it, double *sum_T_tau_it,
                        double **mu, double ***sigma, double *omega_Q3,
                        double omega_sigma, int constrained, int block,
                        int *blocks, int **graph, int N, int T, int D, 
                        int verbose);

int learn_gauss_output_hmm(double **data, int regularize, int reg_type,
                           double omega_Q1, double omega_Q2, double omega_Q3,
                           double omega_sigma, int anneal, double anneal_step,
                           double anneal_factor, double beta_min,
                           double peps, double thresh, int *covgraph,
                           int maxiters,
                           double *L, int *Q, double *pi, double **A, 
                           double **mu, double ***sigma, int N, int T, int D,
                           int viterbi, int verbose);

int eval_gauss_output_hmm(double **data, double *L, int *Q, double *pi,
                          double **A, double **mu, double ***sigma, int N,
                          int T, int D, int viterbi, int verbose);

#endif
