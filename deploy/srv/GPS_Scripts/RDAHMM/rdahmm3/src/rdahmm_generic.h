/*******************************************************************************
MODULE HEADER:
rdahmm_generic.h
*******************************************************************************/
#ifndef _RDAHMM_GENERIC_H_
#define _RDAHMM_GENERIC_H_ 1

#ifndef lint
static char rdahmm_generic_hdr_rcsid[] = "$Id: rdahmm_generic.h,v 1.1 2003/08/27 18:06:15 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_generic.h,v $
 * Revision 1.1  2003/08/27 18:06:15  granat
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
int read_init_probs(char *pi_file, int N, double *pi);

int read_transition_matrix(char *A_file, int N, double **A);

double calc_log_likelihood(double *scale, int T);

double calc_Q_function(double **tau_it, double ***tau_ijt, double *pi,
                       double **A, double **prob_data, int N, int T);

int update_init_probs(double *pi, double **tau_it, double omega_Q1,
                      int regularization, int N);

int update_trans_probs(double **A, double **sum_t_tau_ijt,
                       double *sum_t_tau_it, double omega_Q2,
                       int regularization, int N);

int forward_backward(double *pi, double **A, double **prob_data,
                     double **alpha, double **beta, double *scale,
                     int N, int T);

int estimate_hidden(double **alpha, double **beta, double **A,
                    double **prob_data, double **tau_it, double ***tau_ijt,
                    int N, int T);

int viterbi_state_assignment(double **prob_data, double *pi, double **A, 
                            int *Q, int N, int T);

#endif
