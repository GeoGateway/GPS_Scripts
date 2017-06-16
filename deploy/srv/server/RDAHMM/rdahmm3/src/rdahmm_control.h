/*******************************************************************************
MODULE HEADER:
rdahmm_control.h
*******************************************************************************/
#ifndef _RDAHMM_CONTROL_H_
#define _RDAHMM_CONTROL_H_ 1

#ifndef lint
static char rdahmm_control_hdr_rcsid[] = "$Id$";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log$
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
int test_gauss_output_hmm(double **cont_data, Gaussian_HMM hmm,
                          int *Q, double *L, int T, int viterbi, int verbose);

int train_gauss_output_hmm(double **cont_data, Gaussian_HMM hmm,
                           RDAEM_Parameters params, int *Q, double *L,
                           int T, int viterbi, int verbose);

#endif
