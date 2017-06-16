/*******************************************************************************
MODULE NAME
rdahmm_stat_utils

AUTHOR
Robert Granat

DESCRIPTION
See Robert Granat PhD thesis for details and notation.

COMPILE
make

NOTES
Not yet available

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: rdahmm_stat_utils.c,v 1.2 2004/10/04 22:01:46 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: rdahmm_stat_utils.c,v $
 * Revision 1.2  2004/10/04 22:01:46  granat
 * *** empty log message ***
 *
 * Revision 1.1  2003/08/27 18:04:33  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* LAPACK library */
/*
#include <clapack.h>
*/

/* CP library */
#include "cp.h"

/* UT library */
#include "ut.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* local header files */

/* this program's header */
#include "rdahmm_stat_utils.h"

/*******************************************************************************
GENERATE_RAND_COVAR_DMAT
Generate a random sigmaiance matrix in double precision.
RG
*******************************************************************************/
int generate_rand_covar_dmat(double **sigma, int D)
{
  int     i, j, k;
  int     posdef;
  int     sing;
  double  sum, tau;
  double  *c, *d, *w;
  double  **a, **qt;

  c = NR_dvector(1, D);
  d = NR_dvector(1, D);
  w = NR_dvector(1, D);
  qt = NR_dmatrix(1, D, 1, D);

  a = sigma;

  do 
  {
    posdef = UT_FALSE;

    random_dmatrix(a, D, D);

    NR_dqrdcmp(a, D, c, d, &sing);

    if (sing)
      continue;

    fast_zero_dmat(qt, D, D);
    for (i = 1; i <= D; i++)
      qt[i][i] = 1.0;

    for (k = 1; k <= D; k++)
      for (j = 1; j < D; j++) {
        for (sum = 0.0, i = j; i <= D; i++)
          sum += a[i][j] * qt[i][k];
        tau = sum / c[j];
        for (i = j; i <= D; i++)
          qt[i][k] -= tau * a[i][j];
      }

    random_dvector(w, D);

    fast_zero_dmat(sigma, D, D);
    for (i = 1; i <= D; i++)
      for (j = 1; j <= D; j++)
        for (k = 1; k <= D; k++)
          sigma[i][j] += qt[k][i] * qt[k][j] * w[k];

    posdef_sym_alloc_dmat(sigma, D, &posdef);
  }
  while (posdef == UT_FALSE);

  NR_free_dvector(c, 1, D);
  NR_free_dvector(d, 1, D);
  NR_free_dvector(w, 1, D);
  NR_free_dmatrix(qt, 1, D, 1, D);

  return(UT_OK);
}

/*******************************************************************************
SPEC_NORM_COVAR_DMAT
Compute the spectral norm of a covariance matrix of doubles.
RG
*******************************************************************************/
double spec_norm_covar_dmat(double **sigma, int D)
{
  double spec_norm;
  double **a;
  double *w;
  char   jobz = 'N';
  char   uplo = 'U';
  char   range = 'I';
  double vl, vu;      /* not referenced */
  double abstol;      /* absolute error tolerance for eigenvalues */
  int    m;           /* number of eigenvalues returned, should be 1 */ 
  double *z;          /* not referenced */
  double *work;       /* doubles workspace vector */
  int    *iwork;      /* integer workspace vector */
  int    lwork;       /* size of workspace vector "work" */
  double workinfo[1]; /* holds value of optimal size of "work" */
  int    info;        /* return value of DSYEV */
  int    *ifail;      /* not referenced */

  a = NR_dmatrix(1, D, 1, D);
  w = NR_dvector(1, D);
  iwork = NR_ivector(1, 5*D);

  abstol = 0;

  copy_dmat(sigma, a, D, D);

  /* run once initially to determine optimal workspace size */
  lwork = -1;
/*
  dsyev_(&jobz, &uplo, &D, &a[1][1], &D, &w[1], workinfo, &lwork, &info);
*/
  dsyevx_(&jobz, &range, &uplo, &D, &a[1][1], &D, &vl, &vu, &D, &D, &abstol, 
          &m, &w[1], z, &D, workinfo, &lwork, &iwork[1], ifail, &info);

  lwork = workinfo[0];
  work = NR_dvector(1, lwork);
  
/*
  dsyev_(&jobz, &uplo, &D, &a[1][1], &D, &w[1], &work[1], &lwork, &info);
*/
  dsyevx_(&jobz, &range, &uplo, &D, &a[1][1], &D, &vl, &vu, &D, &D, &abstol, 
          &m, &w[1], z, &D, &work[1], &lwork, &iwork[1], ifail, &info);

  if (info != 0)
  {
    fprintf(stderr, "Error calculating spectral norm of the covariance \
            matrix.  DSYEVX returned %d\n", info);
  }

  spec_norm = w[1];

  NR_free_dmatrix(a, 1, D, 1, D);
  NR_free_dvector(w, 1, D);
  NR_free_dvector(work, 1, lwork);
  NR_free_ivector(iwork, 1, 5*D);

  return(spec_norm);
}

/*******************************************************************************
INV_CHOL_COVAR_DMAT
Invert a covariance matrix of doubles.
RG
*******************************************************************************/
int inv_chol_covar_dmat(double **sigma, double **inv_chol_sigma, 
                        double *sqrt_det_sigma, int D)
{
  int    d, d2;
  int    cholesky_ok;
  int    invert_ok;
  double sum;

  copy_dmat(sigma, inv_chol_sigma, D, D);

  /* calculate the Cholesky decomposition of the covariance */
  cholesky_ok = cholesky_dmat(inv_chol_sigma, D);

  if (cholesky_ok == UT_ERROR)
  {
    fflush(stderr);
    fprintf(stderr, "Error performing Cholesky decomposition\n");
    return(UT_ERROR);
  }
  
  /* calculate square root of the determinant */
  *sqrt_det_sigma = 1.0;
  for (d = 1; d <= D; d++)
  {
    *sqrt_det_sigma *= inv_chol_sigma[d][d];
  }

  /* check for numerical problems */
  if ((*sqrt_det_sigma >= DBL_MAX) || (*sqrt_det_sigma <= DBL_MIN))
  {
    /* recalculate determinant if necessary */
    *sqrt_det_sigma = 0.0;
    for (d = 1; d <= D; d++)
    {
      *sqrt_det_sigma += log(inv_chol_sigma[d][d]);
    }
    *sqrt_det_sigma = exp(*sqrt_det_sigma);
  }

  /* invert the Cholesky decomposition */
  invert_ok = invert_lower_dmat(inv_chol_sigma, D);

  if (invert_ok == UT_ERROR)
  {
    fflush(stderr);
    fprintf(stderr, "Error inverting Cholesky decomposition\n");
    return(UT_ERROR);
  }

  /* zero out the upper triangular part of the matrix */
  /* needed for later matrix operations */
/*
  for (d = D; d >= 1; --d)
  {
    for (d2 = D; d2 > d; --d2)
    {
      inv_chol_sigma[d][d2] = 0.0;
    }
  }
*/

  return(UT_OK);
}

/*******************************************************************************
INV_COVAR_DMAT
Invert a covariance matrix of doubles.
RG
*******************************************************************************/
int inv_covar_dmat(double **sigma, double **inv_sigma, int D)
{
  int    d, d2, d3;
  int    cholesky_ok;
  int    invert_ok;
  double diag;
  double **chol_sigma;

  copy_dmat(sigma, inv_sigma, D, D);

  chol_sigma = inv_sigma;

  /* calculate the Cholesky decomposition of the covariance */
  cholesky_ok = cholesky_dmat(chol_sigma, D);

  /* invert the Cholesky decomposition */
  invert_ok = invert_lower_dmat(chol_sigma, D);

  /* calculate full matrix */

  /* fill in strictly upper triangular portion */
  for (d = 1; d <= D; d++)
  {
    for (d2 = d+1; d2 <= D; d2++)
    {
      inv_sigma[d][d2] = 0.0;
      for (d3 = 1; d3 <= NR_min(d,d2); d3++)
      {
        inv_sigma[d][d2] += chol_sigma[d][d3] * chol_sigma[d2][d3];
      }
    }
  }

  /* calculate the values along the diagonal */
  for (d = 1; d <= D; d++)
  {
    diag = 0.0;
    for (d2 = 1; d2 <= d; d2++)
    {
      diag += chol_sigma[d][d2];
    }
    inv_sigma[d][d] = diag;
  }

  /* mirror the matrix to complete */
  for (d = 1; d <= D; d++)
  {
    for (d2 = 1; d2 < d; d2++)
    {
      inv_sigma[d][d2] = inv_sigma[d2][d];
    }
  }

  return(UT_OK);
}


/*******************************************************************************
CHOLESKY_DMAT
Perform an in-place Cholesky decomposition of a matrix of doubles.
RG
*******************************************************************************/
int cholesky_dmat(double **mat, int size)
{
  double fac;   
  int k;
  int j;
  int i;

  for (k = 1; k <= size; k++) 
  {
    if (mat[k][k] <= 0)
    {
      return(UT_ERROR);
    }
    mat[k][k] = sqrt(mat[k][k]);
    fac = 1/mat[k][k];
    for (j = k+1; j <= size; j++)
    {
      mat[j][k] *= fac;
    }
    for (j = k+1; j <= size; j++)
    {
      for (i = j; i <= size; i++)
      {
        mat[i][j] -= mat[i][k] * mat[j][k];
      }
    }
  }

  return(UT_OK);
}


/*******************************************************************************
INV_LOWER_DMAT
Invert a lower triangular matrix of doubles.
RG
*******************************************************************************/
int invert_lower_dmat(double **l, int n)
{
  int i,j,k;
  double sum;

  for (j = 1; j <= n; j++) 
  {
    if (l[j][j] == 0.0)
    {
      return(UT_ERROR);
    }
    l[j][j] = 1 / l[j][j];
    for (i = j+1; i <= n; i++) 
    {
      sum = 0.0;
      for (k = j; k < i; k++) 
      {
        sum -= l[i][k] * l[k][j];
      }
      l[i][j] = sum/l[i][i];
    }
  }

  return(UT_OK);
}

/*******************************************************************************
 DIFF_DVEC
*******************************************************************************/
int diff_dvec(double *v1, double *v2, double *diff, int D)
{
  int    d;
  
  for (d = D; d >= 1; --d)
  {
    diff[d] = v1[d] - v2[d];
  } 

  return(UT_OK);
}

/*******************************************************************************
 DOT_DVEC
*******************************************************************************/
double dot_dvec(double *v1, double *v2, int D)
{
  int    d;
  double dot;
  
  dot = 0.0;
  for (d = D; d >= 1; --d)
  {
    dot += v1[d] * v2[d];
  } 

  return(dot);
}
