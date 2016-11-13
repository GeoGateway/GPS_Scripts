/*******************************************************************************
MODULE HEADER:
da_nrhacks.h
*******************************************************************************/

#ifndef _DA_NRHACKS_H_
#define _DA_NRHACKS_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_nrhacks_h_rcsid[] = "$Id: da_nrhacks.h,v 1.3 1997/08/24 00:53:45 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_nrhacks.h,v $
 * Revision 1.3  1997/08/24 00:53:45  granat
 * changed function names to DA_ convention
 *
 * Revision 1.2  1997/06/20 22:15:25  granat
 * minor changes
 *
 * Revision 1.1  1997/05/13 20:58:31  granat
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/
typedef float (*AMOEBA_PF)( void *arg, int nargs, float *vals, int n );

/*==============================================================================
Defines and Macros
==============================================================================*/

/* Amoeba defines and macros...  Should be fixed later */

#define NMAX    2000

#define ALPHA   1.0
#define BETA    0.5
#define GAMMA   2.0

#define COST_THRESH 0.0

#define GET_PSUM        for( j = 1; j <= ndim; ++j ) {                  \
                            for( sum = 0, i = 1; i <= mpts; ++i ) {     \
                                sum += p[i][j];                         \
                                psum[j] = sum;                          \
                            }                                           \
                        }

/* End amoeba defines and macros */

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/
void DA_nrerror(char error_text[]);

int DA_svdcmp(float **a, int m, int n, float w[], float **v, int maxiter);
int DA_dsvdcmp(double **a, int m, int n, double w[], double **v, int maxiter);
float DA_amotry( float **p, float *y, float *psum, int ndim, AMOEBA_PF func,
               void *arg, int nargs, int ihi, int *nfunc, float fac );

int DA_amoeba( float **p, float *y, int ndim, float ftol, AMOEBA_PF func, 
               void *arg, int nargs, int *nfunc );

#endif /* _DA_NRHACKS_H_ */
