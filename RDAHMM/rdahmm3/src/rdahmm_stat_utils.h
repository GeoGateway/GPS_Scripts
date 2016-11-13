/*******************************************************************************
MODULE HEADER:
rdahmm_stat_utils.h
*******************************************************************************/
#ifndef _RDAHMM_STAT_UTILS_H_
#define _RDAHMM_STAT_UTILS_H_ 1

#ifndef lint
static char rdahmm_stat_utils_hdr_rcsid[] = "$Id: rdahmm_stat_utils.h,v 1.1 2003/08/27 18:06:31 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_stat_utils.h,v $
 * Revision 1.1  2003/08/27 18:06:31  granat
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
int generate_rand_covar_dmat(double **sigma, int D);

double spec_norm_covar_dmat(double **sigma, int D);

int inv_chol_covar_dmat(double **sigma, double **inv_chol_sigma,
                        double *sqrt_det_sigma, int D);

int inv_covar_dmat(double **sigma, double **inv_sigma, int D);

int cholesky_dmat(double **mat, int size);

int invert_lower_dmat(double **l, int n);

int diff_dvec(double *v1, double *v2, double *diff, int D);

double dot_dvec(double *v1, double *v2, int D);

#endif
