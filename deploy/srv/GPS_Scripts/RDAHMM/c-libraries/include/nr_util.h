/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.h.  Do not confuse this file with the same-named
   file nrutil.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char nr_utils_hdr_rcsid[] = "$Id: nr_util.h,v 1.7 1997/07/29 15:46:12 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: nr_util.h,v $
 * Revision 1.7  1997/07/29 15:46:12  granat
 * revised to use new "NR_" paradigm
 *
 * Revision 1.6  1996/11/19 22:31:01  agray
 * changed fmax(), fmin() to max(), min(), since there was a function
 * called fmin() in NR lib.
 *
 * Revision 1.5  1996/11/19 21:23:34  agray
 * New formatting, added lower-case versions of NR macros, moved
 * two #define's from nrutil.c.
 *
 * Revision 1.4  1996/09/20 21:09:50  granat
 * changed prototype of NR_free_cmatrix so that it uses unsigned chars instead
 * of chars
 *
 * Revision 1.3  1996/09/18 17:15:21  granat
 * Changed cmatrix prototype so that it makes matrices of unsigned chars;
 * Added prototype for csubmatrix function
 *
 * Revision 1.2  1996/07/25 18:15:11  granat
 * added prototype for nrerror2
 *
 * Revision 1.1  1996/07/25 17:28:17  agray
 * Initial revision
 *
 * I made some additions in the initial check-in - cvector(), cmatrix(), 
 * NR_free_cvector(), and NR_free_cmatrix().
 * AG
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

#define NR_END 1
#define NR_FREE_ARG char*

/* original NR macros */

static float sqrarg;
#define NR_SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define NR_DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define NR_DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define NR_DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define NR_FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define NR_FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define NR_LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define NR_LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define NR_IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define NR_IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* lower-case alternatives */

#define NR_sqr(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define NR_dsqr(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

#define NR_max(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define NR_min(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
#define NR_dmax(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))
#define NR_dmin(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))
#define NR_lmax(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))
#define NR_lmin(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))
#define NR_imax(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define NR_imin(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define NR_sign(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* error handling */
void NR_error(char error_text[]);

/* matrix and vector structures */

float *NR_vector(long nl, long nh);
int *NR_ivector(long nl, long nh);
unsigned char *NR_cvector(long nl, long nh);
unsigned long *NR_lvector(long nl, long nh);
double *NR_dvector(long nl, long nh);

float **NR_matrix(long nrl, long nrh, long ncl, long nch);
double **NR_dmatrix(long nrl, long nrh, long ncl, long nch);
int **NR_imatrix(long nrl, long nrh, long ncl, long nch);
unsigned char **NR_cmatrix(long nrl, long nrh, long ncl, long nch);

float **NR_submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
unsigned char **NR_csubmatrix(unsigned char **a, long oldrl, long oldrh, 
                          long oldcl, long oldch, long newrl, long newcl);
float **NR_convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);

float ***NR_f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void NR_free_vector(float *v, long nl, long nh);
void NR_free_ivector(int *v, long nl, long nh);
void NR_free_cvector(unsigned char *v, long nl, long nh);
void NR_free_lvector(unsigned long *v, long nl, long nh);
void NR_free_dvector(double *v, long nl, long nh);

void NR_free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void NR_free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void NR_free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void NR_free_cmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch);

void NR_free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void NR_free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);

void NR_free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);

#endif /* _NR_UTILS_H_ */
