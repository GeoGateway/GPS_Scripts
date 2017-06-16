/*******************************************************************************
MODULE NAME
da_util

ONE-LINE SYNOPSIS
General functions related to matrix and vector structures.

SCOPE OF THIS MODULE
Any functions relating directly to other concepts which have representative 
modules in this library should go in the appropriate module.  Functions that 
apply more generally are intended to go here.  For instance, reading or
writing a matrix from/to a file would more appropriately be placed in 
da_io, because it deals with I/O.

SEE ALSO
Because the definition of this module is quite broad, there is some potential
overlap with several other modules in this library.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/pca, AG.

NOTES
-

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_util.c,v 1.5 1999/07/20 17:35:21 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_util.c,v $
 * Revision 1.5  1999/07/20 17:35:21  granat
 * added double precision versions of some functions
 *
 * Revision 1.4  1998/07/02 01:21:04  granat
 * fixed fast_set functions so that they return ints
 * eliminated redundant variable in scaling functions
 *
 * Revision 1.3  1998/06/30 17:22:28  agray
 * fixed bugs related to 1-base indexing vs. 0-base indexing.
 *
 * Revision 1.2  1998/06/29 22:06:47  granat
 * aborted addition of functions related to sparse matrices
 *
 * Revision 1.1  1998/05/07 23:53:40  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* UT library */
#include "ut_error.h"
#include "ut_math.h"
#include "ut_memory.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr_util.h"

/* DA library */
#include "da_linalg.h"

/* this module's header */
#include "da_util.h"


/*******************************************************************************
RANGE_OF_COLS
For each of k columns in a matrix, find the minimum and maximum values, and
the range (difference between the two).  Pass in a size k vector for each of 
these latter sets of values, and they will be filled in by this function.  
AG
*******************************************************************************/
int range_of_cols(float **mat, float *minval, float *maxval, 
                  float *range, int numrows, int numcols)

{
  int i,j;

  /* for each attribute, find the min and max */
  for (j=1; j<=numcols; j++)
    minval[j] = maxval[j] = mat[1][j];

  for (i=2; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
    {
      if (mat[i][j] < minval[j])
        minval[j] = mat[i][j];
      if (mat[i][j] > maxval[j])
        maxval[j] = mat[i][j];
    }

  /* for each attribute, compute the range */
  for (j=1; j<=numcols; j++)
    range[j] = maxval[j] - minval[j];
  
  return (UT_OK);
}

/*******************************************************************************
RANGE_OF_DCOLS
For each of k columns in a matrix, find the minimum and maximum values, and
the range (difference between the two).  Pass in a size k vector for each of 
these latter sets of values, and they will be filled in by this function.  
AG
*******************************************************************************/
int range_of_dcols(double **mat, double *minval, double *maxval, 
                  double *range, int numrows, int numcols)

{
  int i,j;

  /* for each attribute, find the min and max */
  for (j=1; j<=numcols; j++)
    minval[j] = maxval[j] = mat[1][j];

  for (i=2; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
    {
      if (mat[i][j] < minval[j])
        minval[j] = mat[i][j];
      if (mat[i][j] > maxval[j])
        maxval[j] = mat[i][j];
    }

  /* for each attribute, compute the range */
  for (j=1; j<=numcols; j++)
    range[j] = maxval[j] - minval[j];
  
  return (UT_OK);
}

/*******************************************************************************
RANGE_NORMALIZE_COLS
Takes as input vectors containing the minimum values, and the range of values
(difference between the max and min), for each of the columns of a specified 
matrix.  Modifies the original matrix such that each value in each column is 
normalized by the range specified.  Thus, all entries in the matrix will end up
scaled to [0,1].
AG
*******************************************************************************/
int range_normalize_cols(float **mat, float *minval, float *range, 
                         int numrows, int numcols)

{
  int i,j;

  /* translate each attribute value by its minval and scale by its range */
  for (i=1; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
      mat[i][j] = (mat[i][j] - minval[j]) / range[j];

  return (UT_OK);
}

/*******************************************************************************
RANGE_NORMALIZE_DCOLS
Takes as input vectors containing the minimum values, and the range of values
(difference between the max and min), for each of the columns of a specified 
matrix.  Modifies the original matrix such that each value in each column is 
normalized by the range specified.  Thus, all entries in the matrix will end up
scaled to [0,1].
AG
*******************************************************************************/
int range_normalize_dcols(double **mat, double *minval, double *range, 
                         int numrows, int numcols)

{
  int i,j;

  /* translate each attribute value by its minval and scale by its range */
  for (i=1; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
      mat[i][j] = (mat[i][j] - minval[j]) / range[j];

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_COLS
Restore a matrix of vectors reflecting data which have been normalized by the 
range of values in each column.
AG
*******************************************************************************/
int range_unnormalize_cols(float **mat, float *minval, float *range, 
                         int numrows, int numcols)

{
  int i,j;

  /* scale each attribute value by its range and translate by its minval */
  for (i=1; i<=numrows; i++)
  {
    for (j=1;j<=numcols;j++)
      mat[i][j] = mat[i][j] * range[j] + minval[j];
  }

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_DCOLS
Restore a matrix of vectors reflecting data which have been normalized by the 
range of values in each column.
AG
*******************************************************************************/
int range_unnormalize_dcols(double **mat, double *minval, double *range, 
                         int numrows, int numcols)

{
  int i,j;

  /* scale each attribute value by its range and translate by its minval */
  for (i=1; i<=numrows; i++)
  {
    for (j=1;j<=numcols;j++)
      mat[i][j] = mat[i][j] * range[j] + minval[j];
  }

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_COV_MATRIX
Restore a covariance matrix derived from data which have been normalized by 
the range of values in each column.
AG
*******************************************************************************/
int range_unnormalize_cov_matrix(float **mat, int nc, float *range)


{
  int i,j;

  for (i=1; i<=nc; i++)
    for (j=1;j<=nc;j++)
      mat[i][j] = mat[i][j] * range[i] * range[j];

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_COV_DMATRIX
Restore a covariance matrix derived from data which have been normalized by 
the range of values in each column.
AG
*******************************************************************************/
int range_unnormalize_cov_dmatrix(double **mat, int nc, double *range)


{
  int i,j;

  for (i=1; i<=nc; i++)
    for (j=1;j<=nc;j++)
      mat[i][j] = mat[i][j] * range[i] * range[j];

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_VECTOR
Restore a vector reflecting data which have been normalized by the range of 
values in each attribute (column, often).
AG
*******************************************************************************/
int range_unnormalize_vector(float *vec, int nc, float *range, float *minval)

{
  int i;

  for (i=1; i<=nc; i++)
    vec[i] = vec[i] * range[i] + minval[i];

  return (UT_OK);
}

/*******************************************************************************
RANGE_UNNORMALIZE_DVECTOR
Restore a vector reflecting data which have been normalized by the range of 
values in each attribute (column, often).
AG
*******************************************************************************/
int range_unnormalize_dvector(double *vec, int nc, double *range, double *minval)

{
  int i;

  for (i=1; i<=nc; i++)
    vec[i] = vec[i] * range[i] + minval[i];

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_MATRIX
Print matrix of vectors which have been normalized by the range of values in
each column, in unnormalized form.
AG
******************************************************************************/
int print_unnorm_matrix(stream, nr, nc, mat, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  float  **mat;   /* matrix to be printed */
  float  *range;  /* range of values for each attribute */
  float  *minval; /* minimum value in each attribute */
{
  int i,j;

  for (i=1; i<=nr; i++)
  {
    for (j=1;j<=nc;j++)
      fprintf(stream, "%g ", mat[i][j] * range[j] + minval[j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_DMATRIX
Print matrix of vectors which have been normalized by the range of values in
each column, in unnormalized form.
AG
******************************************************************************/
int print_unnorm_dmatrix(stream, nr, nc, mat, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  double  **mat;   /* matrix to be printed */
  double  *range;  /* range of values for each attribute */
  double  *minval; /* minimum value in each attribute */
{
  int i,j;

  for (i=1; i<=nr; i++)
  {
    for (j=1;j<=nc;j++)
      fprintf(stream, "%g ", mat[i][j] * range[j] + minval[j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_COV_MATRIX
Print covariance matrix which has been normalized by the range of values in
each row/column, in unnormalized form.
AG
*******************************************************************************/
int print_unnorm_cov_matrix(stream, nc, mat, range)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* number of rows, columns of input matrix - should be same */
  float  **mat;   /* matrix to be printed */
  float  *range;  /* range of values for each attribute */
{
  int i,j;

  for (i=1; i<=nc; i++) {
    for (j=1;j<=nc;j++) {
      fprintf(stream, "%g ", mat[i][j] * range[i] * range[j]);
    }
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_COV_DMATRIX
Print covariance matrix which has been normalized by the range of values in
each row/column, in unnormalized form.
AG
*******************************************************************************/
int print_unnorm_cov_dmatrix(stream, nc, mat, range)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* number of rows, columns of input matrix - should be same */
  double  **mat;   /* matrix to be printed */
  double  *range;  /* range of values for each attribute */
{
  int i,j;

  for (i=1; i<=nc; i++) {
    for (j=1;j<=nc;j++) {
      fprintf(stream, "%g ", mat[i][j] * range[i] * range[j]);
    }
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_ROW
Print vector which has been normalized by the range of values in each 
column, in unnormalized form.  Prints in row form.
AG
*******************************************************************************/
int print_unnorm_row(stream, nc, v, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* length of the row vector to be printed */
  float  *v;      /* vector to be printed */
  float  *range;  /* range of values for each attribute */
  float  *minval; /* minimum value in each attribute */
{
  int     i;

  for (i = 1; i <= nc; i++) {
      fprintf(stream, "%g ", v[i] * range[i] + minval[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_DROW
Print vector which has been normalized by the range of values in each 
column, in unnormalized form.  Prints in row form.
AG
*******************************************************************************/
int print_unnorm_drow(stream, nc, v, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* length of the row vector to be printed */
  double  *v;      /* vector to be printed */
  double  *range;  /* range of values for each attribute */
  double  *minval; /* minimum value in each attribute */
{
  int     i;

  for (i = 1; i <= nc; i++) {
      fprintf(stream, "%g ", v[i] * range[i] + minval[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_COL
Print vector which has been normalized by the range of values in each 
row, in unnormalized form.  Prints in column form.
AG
*******************************************************************************/
int print_unnorm_col(stream, nc, v, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* length of the column vector to be printed */
  float  *v;      /* vector to be printed */
  float  *range;  /* range of values for each attribute */
  float  *minval; /* minimum value in each attribute */
{
  int     i;

  for (i = 1; i <= nc; i++)
    fprintf(stream, "%g\n", v[i] * range[i] + minval[i]);

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_DCOL
Print vector which has been normalized by the range of values in each 
row, in unnormalized form.  Prints in column form.
AG
*******************************************************************************/
int print_unnorm_dcol(stream, nc, v, range, minval)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nc;      /* length of the column vector to be printed */
  double  *v;      /* vector to be printed */
  double  *range;  /* range of values for each attribute */
  double  *minval; /* minimum value in each attribute */
{
  int     i;

  for (i = 1; i <= nc; i++)
    fprintf(stream, "%g\n", v[i] * range[i] + minval[i]);

  return (UT_OK);
}

/*******************************************************************************
SUM_VEC
Sum the values in a vector.
RG
*******************************************************************************/
float sum_vec(float *v, int n)
{
  float  *p;
  float  *p_end;
  float   sum = 0.0;;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    sum += *p;

  return (sum);
}


/*******************************************************************************
SUM_DVEC
Sum the values in a vector of doubles.
RG
*******************************************************************************/
double sum_dvec(double *v, int n)
{
  double  *p;
  double  *p_end;
  double   sum = 0.0;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    sum += *p;

  return (sum);
}


/*******************************************************************************
SUM_IVEC
Sum the values in a vector of integers, and return the result as a integer.
RG
*******************************************************************************/
int sum_ivec(int *v, int n)
{
  int    *p;
  int    *p_end;
  int    sum = 0;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    sum += (*p);

  return (sum);
}


/*******************************************************************************
SUM_CVEC
Sum the values in a vector of integers, and return the result as a integer.
RG
*******************************************************************************/
int sum_cvec(unsigned char *v, int n)
{
  unsigned char  *p;
  unsigned char  *p_end;
  int            sum = 0;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    sum += (int) (*p);

  return (sum);
}


/*******************************************************************************
SUM_MAT
Sum the values in a matrix.
RG
*******************************************************************************/
float sum_mat(float **m, int nr, int nc)
{
  float *p;
  float *p_end;
  float sum = 0.0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* the memory is one continuous strip, so it can be stepped through */
  /* in a single loop to reduce branching */

  for (p = &m[1][1]; p <= p_end; p++)
    sum += *p;

  return (sum);
}


/*******************************************************************************
SUM_DMAT
Sum the values in a matrix of doubles, and return the result as a double.
RG
*******************************************************************************/
double sum_dmat(double **m, int nr, int nc)
{
  double *p;
  double *p_end;
  double sum = 0.0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* the memory is one continuous strip, so it can be stepped through */
  /* in a single loop to reduce branching */

  for (p = &m[1][1]; p <= p_end; p++)
    sum += *p;

  return (sum);
}


/*******************************************************************************
SUM_IMAT
Sum the values in a matrix of integers, and return the result as a integer.
RG
*******************************************************************************/
int sum_imat(int **m, int nr, int nc)
{
  int   *p;
  int   *p_end;
  int   sum = 0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* the memory is one continuous strip, so it can be stepped through */
  /* in a single loop to reduce branching */

  for (p = &m[1][1]; p <= p_end; p++)
    sum += (float) (*p);

  return (sum);
}


/*******************************************************************************
SUM_CMAT
Sum the values in a matrix of unsigned chars, and return the result as a 
integer.
RG
*******************************************************************************/
int sum_cmat(unsigned char **m, int nr, int nc)
{
  unsigned char  *p;
  unsigned char  *p_end;
  int            sum = 0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* the memory is one continuous strip, so it can be stepped through */
  /* in a single loop to reduce branching */

  for (p = &m[1][1]; p <= p_end; p++)
    sum += (float) (*p);

  return (sum);
}


/*******************************************************************************
SUM_MAT_ROWS
Compute the vector which is the sum of the rows of a given matrix.
AG
*******************************************************************************/
int sum_mat_rows(float **m, float *v, int nr, int nc)
{
  int   i;
  float *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    v[i] = sum_vec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
SUM_DMAT_ROWS
Compute the vector which is the sum of the rows of a given matrix of doubles.
AG, RG
*******************************************************************************/
int sum_dmat_rows(double **m, double *v, int nr, int nc)
{
  int    i;
  double *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    v[i] = sum_dvec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
SUM_IMAT_ROWS
Compute the vector which is the sum of the rows of a given matrix of integers.
AG, RG
*******************************************************************************/
int sum_imat_rows(int **m, int *v, int nr, int nc)
{
  int    i;
  int    *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    v[i] = sum_ivec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
SUM_CMAT_ROWS
Compute the vector of integers which is the sum of the rows of a given 
matrix of unsigned chars.
AG, RG
*******************************************************************************/
int sum_cmat_rows(unsigned char **m, int *v, int nr, int nc)
{
  int            i;
  unsigned char  *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    v[i] = sum_cvec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
NORMALIZE_VEC
Given a vector, return the sum of the members, and divide each element by that
sum.
RG
*******************************************************************************/
float normalize_vec(float *v, int n)
{
  float sum;

  sum = sum_vec(v, n);
  scalar_div_vec(v, n, sum);

  return (sum);
}


/*******************************************************************************
NORMALIZE_DVEC
Given a vector of doubles, return the sum of the members, and divide each 
element by that sum.
RG
*******************************************************************************/
double normalize_dvec(double *v, int n)
{
  double sum;

  sum = sum_dvec(v, n);
  scalar_div_dvec(v, n, sum);

  return (sum);
}


/*******************************************************************************
NORMALIZE_MAT
Given a matrix, return the sum of the members, and divide each element by that
sum.
RG
*******************************************************************************/
float normalize_mat(float **m, int nr, int nc)
{
  float sum;

  sum = sum_mat(m, nr, nc);
  scalar_div_mat(m, nr, nc, sum);

  return (sum);
}


/*******************************************************************************
NORMALIZE_DMAT
Given a matrix of doubles, return the sum of the members, and divide each 
element by that sum.
RG
*******************************************************************************/
double normalize_dmat(double **m, int nr, int nc)
{
  double sum;

  sum = sum_dmat(m, nr, nc);
  scalar_div_dmat(m, nr, nc, sum);

  return (sum);
}


/*******************************************************************************
NORMALIZE_MAT_ROWS
For each row of a given matrix, divide each element by the sum of the elements
in the row.
 
Note:  Could later have norm_sum_mat_rows_and_keep(), which keeps a vector of
the row sums, returning it (analogous to prob_mixture_and_keep()).
AG
*******************************************************************************/
int normalize_mat_rows(float **m, int nr, int nc)
{
  int   i;
  float *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    normalize_vec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
NORMALIZE_DMAT_ROWS
For each row of a given matrix of doubles, divide each element by the sum of 
the elements in the row.
 
Note:  Could later have norm_sum_mat_rows_and_keep(), which keeps a vector of
the row sums, returning it (analogous to prob_mixture_and_keep()).
AG
*******************************************************************************/
int normalize_dmat_rows(double **m, int nr, int nc)
{
  int   i;
  double *curr_row;

  for (i = 1; i <= nr; i++)
  {
    curr_row = m[i];
    normalize_dvec(curr_row, nc);
  }

  return (UT_OK);
}


/*******************************************************************************
SUM_LOG_VEC
Given a vector, return the sum of the log() of each element.
RG
*******************************************************************************/
float sum_log_vec(float *v, int n)
{
  float *p;
  float *p_end;
  float sum = 0.0;

  /* assign pointer to last element */
  p_end = &v[n];

  for (p = &v[1]; p <= p_end; p++)
    sum += log((double) (*p));

  return (sum);
}


/*******************************************************************************
SUM_LOG_DVEC
Given a vector of doubles, return the sum of the log() of each element.
RG
*******************************************************************************/
double sum_log_dvec(double *v, int n)
{
  double *p;
  double *p_end;
  double sum = 0.0;

  /* assign pointer to last element */
  p_end = &v[n];

  for (p = &v[1]; p <= p_end; p++)
    sum += log(*p);

  return (sum);
}


/*******************************************************************************
SUM_LOG_MAT
Given a matrix, return the sum of the log() of each element.
RG
*******************************************************************************/
float sum_log_mat(float **m, int nr, int nc)
{
  float *p;
  float *p_end;
  float sum = 0.0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  for (p = &m[1][1]; p <= p_end; p++)
    sum += log((double) (*p));

  return (sum);
}


/*******************************************************************************
SUM_LOG_DMAT
Given a matrix of doubles, return the sum of the log() of each element.
RG
*******************************************************************************/
double sum_log_dmat(double **m, int nr, int nc)
{
  double *p;
  double *p_end;
  double sum = 0.0;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  for (p = &m[1][1]; p <= p_end; p++)
    sum += log(*p);

  return (sum);
}


/*******************************************************************************
MAX_VEC
Return the maximum value in a vector.
RG
*******************************************************************************/
float max_vec(float *v, int n)
{
  float  *p;
  float  *p_end;
  float   maxval;

  /* use first element as initial maximum */
  maxval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;

  return (maxval);
}


/*******************************************************************************
MAX_DVEC
Return the maximum value in a vector of doubles.
RG
*******************************************************************************/
double max_dvec(double *v, int n)
{
  double  *p;
  double  *p_end;
  double   maxval;

  /* use first element as initial maximum */
  maxval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;

  return (maxval);
}


/*******************************************************************************
MAX_IVEC
Return the maximum value in a vector of integers.
RG
*******************************************************************************/
int max_ivec(int *v, int n)
{
  int  *p;
  int  *p_end;
  int   maxval;

  /* use first element as initial maximum */
  maxval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;

  return (maxval);
}


/*******************************************************************************
MAX_CVEC
Return the maximum value in a vector of unsigned chars.
RG
*******************************************************************************/
unsigned char max_cvec(unsigned char *v, int n)
{
  unsigned char  *p;
  unsigned char  *p_end;
  unsigned char   maxval;

  /* use first element as initial maximum */
  maxval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;

  return (maxval);
}


/*******************************************************************************
MIN_VEC
Return the minimum value in a vector.
RG
*******************************************************************************/
float min_vec(float *v, int n)
{
  float  *p;
  float  *p_end;
  float   minval;

  /* use first element as initial maximum */
  minval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;

  return (minval);
}


/*******************************************************************************
MIN_DVEC
Return the minimum value in a vector of doubles.
RG
*******************************************************************************/
double min_dvec(double *v, int n)
{
  double  *p;
  double  *p_end;
  double   minval;

  /* use first element as initial maximum */
  minval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;

  return (minval);
}


/*******************************************************************************
MIN_IVEC
Return the minimum value in a vector of integers.
RG
*******************************************************************************/
int min_ivec(int *v, int n)
{
  int  *p;
  int  *p_end;
  int   minval;

  /* use first element as initial maximum */
  minval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;

  return (minval);
}


/*******************************************************************************
MIN_CVEC
Return the minimum value in a vector of unsigned chars.
RG
*******************************************************************************/
unsigned char min_cvec(unsigned char *v, int n)
{
  unsigned char  *p;
  unsigned char  *p_end;
  unsigned char   minval;

  /* use first element as initial maximum */
  minval = v[1];

  /* assign pointer to last element */
  p_end = &v[n];

  /* search for maximum element */
  for (p = &v[2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;

  return (minval);
}


/*******************************************************************************
MAX_MAT
Return the maximum value of a matrix.
RG
*******************************************************************************/
float max_mat(float **m, int nr, int nc)
{
  float *p;
  float *p_end;
  float maxval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  maxval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;
 
  return (maxval);
}


/*******************************************************************************
MAX_DMAT
Return the maximum value of a matrix of doubles.
RG
*******************************************************************************/
double max_dmat(double **m, int nr, int nc)
{
  double *p;
  double *p_end;
  double maxval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  maxval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;
 
  return (maxval);
}


/*******************************************************************************
MAX_IMAT
Return the maximum value of a matrix of integers.
RG
*******************************************************************************/
int max_imat(int **m, int nr, int nc)
{
  int *p;
  int *p_end;
  int maxval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  maxval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;
 
  return (maxval);
}


/*******************************************************************************
MAX_CMAT
Return the maximum value of a matrix of unsigned chars.
RG
*******************************************************************************/
unsigned char max_cmat(unsigned char **m, int nr, int nc)
{
  unsigned char *p;
  unsigned char *p_end;
  unsigned char maxval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  maxval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p > maxval)
      maxval = *p;
 
  return (maxval);
}


/*******************************************************************************
MIN_MAT
Return the minimum value of a matrix.
RG
*******************************************************************************/
float min_mat(float **m, int nr, int nc)
{
  float *p;
  float *p_end;
  float minval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  minval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;
 
  return (minval);
}


/*******************************************************************************
MIN_DMAT
Return the minimum value of a matrix of doubles.
RG
*******************************************************************************/
double min_dmat(double **m, int nr, int nc)
{
  double *p;
  double *p_end;
  double minval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  minval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;
 
  return (minval);
}


/*******************************************************************************
MIN_IMAT
Return the minimum value of a matrix of integers.
RG
*******************************************************************************/
int min_imat(int **m, int nr, int nc)
{
  int *p;
  int *p_end;
  int minval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  minval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;
 
  return (minval);
}


/*******************************************************************************
MIN_CMAT
Return the minimum value of a matrix of unsigned chars.
RG
*******************************************************************************/
int min_cmat(unsigned char **m, int nr, int nc)
{
  unsigned char *p;
  unsigned char *p_end;
  unsigned char minval;
 
  /* assign pointer to last element */
  p_end = &m[nr][nc];
 
  minval = m[1][1];
 
  for (p = &m[1][2]; p <= p_end; p++)
    if (*p < minval)
      minval = *p;
 
  return (minval);
}


/*******************************************************************************
ARG_MAX_VEC
Return the index of the maximum value in a vector.
RG
*******************************************************************************/
int arg_max_vec(float *v, int n)
{
  int     i;
  float  *p;
  float   maxval;
  int     arg;

  maxval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p > maxval) {
      maxval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MAX_DVEC
Return the index of the maximum value in a vector of doubles.
RG
*******************************************************************************/
int arg_max_dvec(double *v, int n)
{
  int     i;
  double *p;
  double  maxval;
  int     arg;

  maxval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p > maxval) {
      maxval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MAX_IVEC
Return the index of the maximum value in a vector of integers.
RG
*******************************************************************************/
int arg_max_ivec(int *v, int n)
{
  int     i;
  int    *p;
  int     maxval;
  int     arg;

  maxval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p > maxval) {
      maxval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MAX_CVEC
Return the index of the maximum value in a vector of unsigned chars.
RG
*******************************************************************************/
int arg_max_cvec(unsigned char *v, int n)
{
  int               i;
  unsigned char    *p;
  unsigned char     maxval;
  int               arg;

  maxval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p > maxval) {
      maxval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MIN_VEC
Return the index of the minimum value in a vector.
RG
*******************************************************************************/
int arg_min_vec(float *v, int n)
{
  int     i;
  float  *p;
  float   minval;
  int     arg;

  minval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p < minval) {
      minval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MIN_DVEC
Return the index of the minimum value in a vector of doubles.
RG
*******************************************************************************/
int arg_min_dvec(double *v, int n)
{
  int     i;
  double *p;
  double  minval;
  int     arg;

  minval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p < minval) {
      minval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MIN_IVEC
Return the index of the minimum value in a vector of integers.
RG
*******************************************************************************/
int arg_min_ivec(int *v, int n)
{
  int     i;
  int    *p;
  int     minval;
  int     arg;

  minval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p < minval) {
      minval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
ARG_MIN_CVEC
Return the index of the minimum value in a vector of unsigned chars.
RG
*******************************************************************************/
int arg_min_cvec(unsigned char *v, int n)
{
  int               i;
  unsigned char    *p;
  unsigned char     minval;
  int               arg;

  minval = v[1];
  arg = 1;

  for (i = 2, p = &v[2]; i <= n; i++, p++)
    if (*p < minval) {
      minval = *p;
      arg = i;
    }

  return (arg);
}


/*******************************************************************************
MINMAX_OF_COLS
For each of k columns in a matrix, find the minimum and maximum values, and
the range (difference between the two).  Pass in a size k vector for each of 
these latter sets of values, and they will be filled in by this function.  
AG
*******************************************************************************/
int minmax_of_cols(float **mat, int nr, int nc, float *minval, float *maxval, 
                   float *range)
{
  int i, j;

  /* for each attribute, find the min and max */
  for (j = 1; j <= nc; j++)
    minval[j] = maxval[j] = mat[1][j];

  for (i = 2; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      if (mat[i][j] < minval[j])
        minval[j] = mat[i][j];
      if (mat[i][j] > maxval[j])
        maxval[j] = mat[i][j];
    }

  /* for each attribute, compute the range */
  for (j = 1; j <= nc; j++)
    range[j] = maxval[j] - minval[j];
  
  return (UT_OK);
}


/*******************************************************************************
RANGE_VEC
Return the range of values (difference between max and min values) in a vector.
This value may be negative.
RG
*******************************************************************************/
float range_vec(float *v, int n)
{
  float  *p;
  float  *p_end;
  float   minval, maxval;

  /* use first element as initial minimum and maximum */
  minval = v[1];
  maxval = v[1];
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  return (maxval - minval);
}


/*******************************************************************************
RANGE_DVEC
Return the range of values (difference between max and min values) in a vector
of doubles.  This value may be negative.
RG
*******************************************************************************/
double range_dvec(double *v, int n)
{
  double  *p;
  double  *p_end;
  double   minval, maxval;

  /* use first element as initial minimum and maximum */
  minval = v[1];
  maxval = v[1];
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  return (maxval - minval);
}


/*******************************************************************************
RANGE_MAT
Return the range of values (difference between max and min values) in a matrix.
This value may be negative.
RG
*******************************************************************************/
float range_mat(float **m, int nr, int nc)
{
  float  *p;
  float  *p_end;
  float   minval, maxval;

  /* use first element as initial minimum and maximum */
  minval = m[1][1];
  maxval = m[1][1];
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  return (maxval - minval);
}


/*******************************************************************************
RANGE_DMAT
Return the range of values (difference between max and min values) in a matrix
of doubles.  This value may be negative.
RG
*******************************************************************************/
double range_dmat(double **m, int nr, int nc)
{
  double  *p;
  double  *p_end;
  double   minval, maxval;

  /* use first element as initial minimum and maximum */
  minval = m[1][1];
  maxval = m[1][1];
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  return (maxval - minval);
}


/*******************************************************************************
SCALE_VEC
Shift and scale the values of a vector so that they lie between 0 and 1.
RG
*******************************************************************************/
int scale_vec(float *v, int n)
{
  float  *p;
  float  *p_end;
  float   minval;
  float   maxval;
  float   scaling;
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* use first element as initial minimum and maximum */
  minval = v[1];
  maxval = v[1];
  
  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  /* calculate scaling factor */
  scaling = 1.0 / (maxval - minval);

  /* now rescale elements onto new range of values */
  for (p = &v[1]; p <= p_end; p++) {
    *p -= minval;
    *p *= scaling;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_DVEC
Shift and scale the values of a vector of doubles so that they lie between 0 
and 1.
RG
*******************************************************************************/
int scale_dvec(double *v, int n)
{
  double  *p;
  double  *p_end;
  double   minval;
  double   maxval;
  double   scaling;
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* use first element as initial minimum and maximum */
  minval = v[1];
  maxval = v[1];
  
  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  /* calculate scaling factor */
  scaling = 1.0 / (maxval - minval);

  /* now rescale elements onto new range of values */
  for (p = &v[1]; p <= p_end; p++) {
    *p -= minval;
    *p *= scaling;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_MAT
Shift and scale the values of a matrix so that they lie between 0 and 1.
RG
*******************************************************************************/
int scale_mat(float **m, int nr, int nc)
{
  float  *p;
  float  *p_end;
  float   minval;
  float   maxval;
  float   scaling;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* use first element as initial minimum and maximum */
  minval = m[1][1];
  maxval = m[1][1];
  
  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  /* calculate scaling factor */
  scaling = 1.0 / (maxval - minval);

  /* now rescale elements onto new range of values */
  for (p = &m[1][2]; p <= p_end; p++) {
    *p -= minval;
    *p *= scaling;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_DMAT
Shift and scale the values of a matrix of doubles so that they lie between 0 
and 1.
RG
*******************************************************************************/
int scale_dmat(double **m, int nr, int nc)
{
  double  *p;
  double  *p_end;
  double   minval;
  double   maxval;
  double   scaling;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* use first element as initial minimum and maximum */
  minval = m[1][1];
  maxval = m[1][1];
  
  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < minval)
      minval = *p;
    if (*p > maxval)
      maxval = *p;
  }

  /* calculate scaling factor */
  scaling = 1.0 / (maxval - minval);

  /* now rescale elements onto new range of values */
  for (p = &m[1][2]; p <= p_end; p++) {
    *p -= minval;
    *p *= scaling;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_RANGE_VEC
Shift and scale the values of a vector so that they lie in a specified
range of values.

minval is the new minimum value.
maxval is the new maximum value.
RG
*******************************************************************************/
int scale_range_vec(float *v, int n, float minval, float maxval)
{
  float  *p;
  float  *p_end;
  float   oldminval;
  float   oldmaxval;
  float   scaling;
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* use first element as initial minimum and maximum of original vector */
  oldminval = v[1];
  oldmaxval = v[1];
  
  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < oldminval)
      oldminval = *p;
    if (*p > oldmaxval)
      oldmaxval = *p;
  }

  /* calculate scaling factor */
  scaling = (maxval - minval) / (oldmaxval - oldminval);

  /* now rescale elements onto new range of values */
  for (p = &v[1]; p <= p_end; p++) {
    *p -= oldminval;
    *p *= scaling;
    *p += minval;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_RANGE_DVEC
Shift and scale the values of a vector of doubles so that they lie in a 
specified range of values.

minval is the new minimum value.
maxval is the new maximum value.
RG
*******************************************************************************/
int scale_range_dvec(double *v, int n, double minval, double maxval)
{
  double  *p;
  double  *p_end;
  double   oldminval;
  double   oldmaxval;
  double   scaling;
  
  /* assign pointer to last element */
  p_end = &v[n];

  /* use first element as initial minimum and maximum of original vector */
  oldminval = v[1];
  oldmaxval = v[1];
  
  /* search for minimum and maximum elements */
  for (p = &v[2]; p <= p_end; p++) {
    if (*p < oldminval)
      oldminval = *p;
    if (*p > oldmaxval)
      oldmaxval = *p;
  }

  /* calculate scaling factor */
  scaling = (maxval - minval) / (oldmaxval - oldminval);

  /* now rescale elements onto new range of values */
  for (p = &v[1]; p <= p_end; p++) {
    *p -= oldminval;
    *p *= scaling;
    *p += minval;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_RANGE_MAT
Shift and scale the values of a matrix so that they lie in a specified
range of values.

minval is the new minimum value.
maxval is the new maximum value.
RG
*******************************************************************************/
int scale_range_mat(float **m, int nr, int nc, float minval, float maxval)
{
  float  *p;
  float  *p_end;
  float   oldminval;
  float   oldmaxval;
  float   scaling;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* use first element as initial minimum and maximum of original vector */
  oldminval = m[1][1];
  oldmaxval = m[1][1];
  
  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < oldminval)
      oldminval = *p;
    if (*p > oldmaxval)
      oldmaxval = *p;
  }

  /* calculate scaling factor */
  scaling = (maxval - minval) / (oldmaxval - oldminval);

  /* now rescale elements onto new range of values */
  for (p = &m[1][2]; p <= p_end; p++) {
    *p -= oldminval;
    *p *= scaling;
    *p += minval;
  }

  return (UT_OK);
}


/*******************************************************************************
SCALE_RANGE_DMAT
Shift and scale the values of a matrix of doubles so that they lie in a 
specified range of values.

minval is the new minimum value.
maxval is the new maximum value.
RG
*******************************************************************************/
int scale_range_dmat(double **m, int nr, int nc, double minval, double maxval)
{
  double  *p;
  double  *p_end;
  double   oldminval;
  double   oldmaxval;
  double   scaling;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* use first element as initial minimum and maximum of original vector */
  oldminval = m[1][1];
  oldmaxval = m[1][1];
  
  /* search for minimum and maximum elements */
  for (p = &m[1][2]; p <= p_end; p++) {
    if (*p < oldminval)
      oldminval = *p;
    if (*p > oldmaxval)
      oldmaxval = *p;
  }

  /* calculate scaling factor */
  scaling = (maxval - minval) / (oldmaxval - oldminval);

  /* now rescale elements onto new range of values */
  for (p = &m[1][2]; p <= p_end; p++) {
    *p -= oldminval;
    *p *= scaling;
    *p += minval;
  }

  return (UT_OK);
}


/*******************************************************************************
FLIP_VECTOR
Change a vector so that it has its elements in reverse order.
RG
*******************************************************************************/
int flip_vector(float *v, int n)
{
  float *p, *p_flip, *p_half_n;
  float  temp;
  int    half_n;

  half_n = (int) (n / 2);
  p_half_n = &v[half_n];

  p = &v[1];
  p_flip = &v[n];

  for (p = &v[1], p_flip = &v[n]; p <= p_half_n; p++, p_flip--) {
    temp = *p;
    *p = *p_flip;
    *p_flip = temp;
  }

  return (UT_OK);
}


/*******************************************************************************
FLIP_DVECTOR
Change a vector of doubles so that it has its elements in reverse order.
RG
*******************************************************************************/
int flip_dvector(double *v, int n)
{
  double *p, *p_flip, *p_half_n;
  double  temp;
  int     half_n;

  half_n = (int) (n / 2);
  p_half_n = &v[half_n];

  p = &v[1];
  p_flip = &v[n];

  for (p = &v[1], p_flip = &v[n]; p <= p_half_n; p++, p_flip--) {
    temp = *p;
    *p = *p_flip;
    *p_flip = temp;
  }

  return (UT_OK);
}


/*******************************************************************************
FLIP_IVECTOR
Change a vector of integers so that it has its elements in reverse order.
RG
*******************************************************************************/
int flip_ivector(int *v, int n)
{
  int *p, *p_flip, *p_half_n;
  int  temp;
  int  half_n;

  half_n = (int) (n / 2);
  p_half_n = &v[half_n];

  p = &v[1];
  p_flip = &v[n];

  for (p = &v[1], p_flip = &v[n]; p <= p_half_n; p++, p_flip--) {
    temp = *p;
    *p = *p_flip;
    *p_flip = temp;
  }

  return (UT_OK);
}


/*******************************************************************************
FLIP_CVECTOR
Change a vector of unsigned chars so that it has its elements in reverse order.
RG
*******************************************************************************/
int flip_cvector(unsigned char *v, int n)
{
  unsigned char *p, *p_flip, *p_half_n;
  unsigned char  temp;
  int            half_n;

  half_n = (int) (n / 2);
  p_half_n = &v[half_n];

  p = &v[1];
  p_flip = &v[n];

  for (p = &v[1], p_flip = &v[n]; p <= p_half_n; p++, p_flip--) {
    temp = *p;
    *p = *p_flip;
    *p_flip = temp;
  }

  return (UT_OK);
}


/*******************************************************************************
FLIP_LEFT_RIGHT_MATRIX
Change a matrix so that its columns are in reverse order.
RG
*******************************************************************************/
int flip_left_right_matrix(float **a, int nr, int nc)
{
  int    i, j;
  float *p_left, *p_right;
  float  temp;
  int    half_nc;
 
  half_nc = (int) (nc / 2);

  p_left = &a[1][1];
 
  for (i = 1; i <= nr; i++) {

    /* start index pointers at the beginning and end of a row */

    p_left = a[i];
    p_right = &a[i][nc];

    /* switch the values and move each pointer towards the center */

    for (j = 1; j <= half_nc; j++) {
      temp = *p_left;
      *p_left = *p_right;
      *p_right = temp;

      p_left++; /* move the pointer to the right */
      p_right--; /* move the pointer to the left */
    }
  }
 
  return (UT_OK);
}
 

/*******************************************************************************
FLIP_LEFT_RIGHT_DMATRIX
Change a matrix of doubles so that its columns are in reverse order.
RG
*******************************************************************************/
int flip_left_right_dmatrix(double **a, int nr, int nc)
{
  int     i, j;
  double *p_left, *p_right;
  double  temp;
  int     half_nc;
 
  half_nc = (int) (nc / 2);

  p_left = &a[1][1];
 
  for (i = 1; i <= nr; i++) {

    /* start index pointers at the beginning and end of a row */

    p_left = a[i];
    p_right = &a[i][nc];

    /* switch the values and move each pointer towards the center */

    for (j = 1; j <= half_nc; j++) {
      temp = *p_left;
      *p_left = *p_right;
      *p_right = temp;

      p_left++; /* move the pointer to the right */
      p_right--; /* move the pointer to the left */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_LEFT_RIGHT_IMATRIX
Change a matrix of integers so that its columns are in reverse order.
RG
*******************************************************************************/
int flip_left_right_imatrix(int **a, int nr, int nc)
{
  int  i, j;
  int *p_left, *p_right;
  int  temp;
  int  half_nc;
 
  half_nc = (int) (nc / 2);

  p_left = &a[1][1];
 
  for (i = 1; i <= nr; i++) {

    /* start index pointers at the beginning and end of a row */

    p_left = a[i];
    p_right = &a[i][nc];

    /* switch the values and move each pointer towards the center */

    for (j = 1; j <= half_nc; j++) {
      temp = *p_left;
      *p_left = *p_right;
      *p_right = temp;

      p_left++; /* move the pointer to the right */
      p_right--; /* move the pointer to the left */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_LEFT_RIGHT_CMATRIX
Change a matrix of unsigned chars so that its columns are in reverse order.
RG
*******************************************************************************/
int flip_left_right_cmatrix(unsigned char **a, int nr, int nc)
{
  int            i, j;
  unsigned char *p_left, *p_right;
  unsigned char  temp;
  int            half_nc;
 
  half_nc = (int) (nc / 2);

  p_left = &a[1][1];
 
  for (i = 1; i <= nr; i++) {

    /* start index pointers at the beginning and end of a row */

    p_left = a[i];
    p_right = &a[i][nc];

    /* switch the values and move each pointer towards the center */

    for (j = 1; j <= half_nc; j++) {
      temp = *p_left;
      *p_left = *p_right;
      *p_right = temp;

      p_left++; /* move the pointer to the right */
      p_right--; /* move the pointer to the left */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_TOP_BOTTOM_MATRIX
Change a matrix so that its rows are in reverse order.
RG
*******************************************************************************/
int flip_top_bottom_matrix(float **a, int nr, int nc)
{
  int    i, j;
  float *p_top, *p_bottom;
  float  temp;
  int    half_nr;
 
  half_nr = (int) (nr / 2);

  p_top = &a[1][1];
 
  for (i = 1; i <= half_nr; i++) {

    /* one pointer naturally starts at the bottommost unflipped top row */
    /* start the second pointer at the topmost unflipped bottom row */

    p_bottom = a[nr - i + 1];

    /* switch the values and move each pointer down the row */

    for (j = 1; j <= nc; j++) {
      temp = *p_top;
      *p_top = *p_bottom;
      *p_bottom = temp;

      p_top++; /* move the pointer right */
      p_bottom++; /* move the pointer right */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_TOP_BOTTOM_DMATRIX
Change a matrix of doubles so that its rows are in reverse order.
RG
*******************************************************************************/
int flip_top_bottom_dmatrix(double **a, int nr, int nc)
{
  int     i, j;
  double *p_top, *p_bottom;
  double  temp;
  int     half_nr;
 
  half_nr = (int) (nr / 2);

  p_top = &a[1][1];
 
  for (i = 1; i <= half_nr; i++) {

    /* one pointer naturally starts at the bottommost unflipped top row */
    /* start the second pointer at the topmost unflipped bottom row */

    p_bottom = a[nr - i + 1];

    /* switch the values and move each pointer down the row */

    for (j = 1; j <= nc; j++) {
      temp = *p_top;
      *p_top = *p_bottom;
      *p_bottom = temp;

      p_top++; /* move the pointer right */
      p_bottom++; /* move the pointer right */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_TOP_BOTTOM_IMATRIX
Change a matrix of integers so that its rows are in reverse order.
RG
*******************************************************************************/
int flip_top_bottom_imatrix(int **a, int nr, int nc)
{
  int  i, j;
  int *p_top, *p_bottom;
  int  temp;
  int  half_nr;
 
  half_nr = (int) (nr / 2);

  p_top = &a[1][1];
 
  for (i = 1; i <= half_nr; i++) {

    /* one pointer naturally starts at the bottommost unflipped top row */
    /* start the second pointer at the topmost unflipped bottom row */

    p_bottom = a[nr - i + 1];

    /* switch the values and move each pointer down the row */

    for (j = 1; j <= nc; j++) {
      temp = *p_top;
      *p_top = *p_bottom;
      *p_bottom = temp;

      p_top++; /* move the pointer right */
      p_bottom++; /* move the pointer right */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
FLIP_TOP_BOTTOM_CMATRIX
Change a matrix of unsigned chars so that its rows are in reverse order.
RG
*******************************************************************************/
int flip_top_bottom_cmatrix(unsigned char **a, int nr, int nc)
{
  int            i, j;
  unsigned char *p_top, *p_bottom;
  unsigned char  temp;
  int            half_nr;
 
  half_nr = (int) (nr / 2);

  p_top = &a[1][1];
 
  for (i = 1; i <= half_nr; i++) {

    /* one pointer naturally starts at the bottommost unflipped top row */
    /* start the second pointer at the topmost unflipped bottom row */

    p_bottom = a[nr - i + 1];

    /* switch the values and move each pointer down the row */

    for (j = 1; j <= nc; j++) {
      temp = *p_top;
      *p_top = *p_bottom;
      *p_bottom = temp;

      p_top++; /* move the pointer right */
      p_bottom++; /* move the pointer right */
    }
  }
 
  return (UT_OK);
}


/*******************************************************************************
GRAB_ROW_MAT
Create the vector obtained from the specified row of a matrix.
RG
*******************************************************************************/
int grab_row_mat(float **mat, int index, int nc, float *vec)
{
  int    mem_size;
  float *p_vec;
  float *p_mat;

  mem_size = nc * sizeof(float);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_vec, p_mat, mem_size);

  return (UT_OK);
}


/*******************************************************************************
GRAB_ROW_DMAT
Create the vector obtained from the specified row of a matrix of doubles.
RG
*******************************************************************************/
int grab_row_dmat(double **mat, int index, int nc, double *vec)
{
  int     mem_size;
  double *p_vec;
  double *p_mat;

  mem_size = nc * sizeof(double);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_vec, p_mat, mem_size);

  return (UT_OK);
}


/*******************************************************************************
GRAB_ROW_IMAT
Create the vector obtained from the specified row of a matrix of integers.
RG
*******************************************************************************/
int grab_row_imat(int **mat, int index, int nc, int *vec)
{
  int  mem_size;
  int *p_vec;
  int *p_mat;

  mem_size = nc * sizeof(int);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_vec, p_mat, mem_size);

  return (UT_OK);
}


/*******************************************************************************
GRAB_ROW_CMAT
Create the vector obtained from the specified row of a matrix of unsigned chars.
RG
*******************************************************************************/
int grab_row_cmat(unsigned char **mat, int index, int nc, unsigned char *vec)
{
  int            mem_size;
  unsigned char *p_vec;
  unsigned char *p_mat;

  mem_size = nc * sizeof(unsigned char);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_vec, p_mat, mem_size);

  return (UT_OK);
}


/*******************************************************************************
GRAB_COL_MAT
Create the vector obtained from the specified column of a matrix.
AG
*******************************************************************************/
int grab_col_mat(float **mat, int index, int nr, float *vec)
{
  int   j;

  for (j=1; j <= nr; j++)
    vec[j] = mat[j][index];

  return (UT_OK);
}


/*******************************************************************************
GRAB_COL_DMAT
Create the vector obtained from the specified column of a matrix of doubles.
AG
*******************************************************************************/
int grab_col_dmat(double **mat, int index, int nr, double *vec)
{
  int   j;

  for (j=1; j <= nr; j++)
    vec[j] = mat[j][index];

  return (UT_OK);
}


/*******************************************************************************
GRAB_COL_IMAT
Create the vector obtained from the specified column of a matrix of integers.
AG
*******************************************************************************/
int grab_col_imat(int **mat, int index, int nr, int *vec)
{
  int   j;

  for (j=1; j <= nr; j++)
    vec[j] = mat[j][index];

  return (UT_OK);
}


/*******************************************************************************
GRAB_COL_CMAT
Create the vector obtained from the specified column of a matrix of unsigned
chars.
AG
*******************************************************************************/
int grab_col_cmat(unsigned char **mat, int index, int nr, unsigned char *vec)
{
  int   j;

  for (j=1; j <= nr; j++)
    vec[j] = mat[j][index];

  return (UT_OK);
}


/*******************************************************************************
COPY_MAT_SECTION
Copy a rectangular section of one matrix into a part of another 
matrix by means of direct memory copy.  (a_tlr, a_tlc) is the top
left corner of the section in the source matrix.  (b_tlr, b_tlc) is
the top left corner of the section in the destination matrix.  (nr, nc)
are the dimensions of the matrix section to be copied.
RG
*******************************************************************************/
int copy_mat_section(float **a, float **b, int a_tlr, 
                      int a_tlc, int b_tlr, int b_tlc, int nr, int nc)
{
  int   i;
  int   mem_size;
  float *p_a, *p_b;

  /* calculate the number of bytes in a row of the section */

  mem_size = nc * sizeof(float);

  /* copy a section of A into B */

  for (i = 1; i <= nr; i++) {

    /* set pointers to the beginning of each row */

    p_a = &a[a_tlr + i - 1][a_tlc];
    p_b = &b[b_tlr + i - 1][b_tlc];

    /* copy a row of A into B */

    memcpy(p_b, p_a, mem_size);
  }

  return (UT_OK);
}


/*******************************************************************************
COPY_IMAT_SECTION
Copy a rectangular section of one integer matrix into a part of another 
matrix by means of direct memory copy.  (a_tlr, a_tlc) is the top
left corner of the section in the source matrix.  (b_tlr, b_tlc) is
the top left corner of the section in the destination matrix.
matrix by means of direct memory copy.  (nr, nc) are the dimensions of
the matrix section to be copied.
RG
*******************************************************************************/
int copy_imat_section(int **a, int **b, int a_tlr, int a_tlc, int b_tlr, 
                       int b_tlc, int nr, int nc)
{
  int   i;
  int   mem_size;
  int   *p_a, *p_b;

  mem_size = nc * sizeof(int);

  /* copy a section of A into B */

  for (i = 1; i <= nr; i++) {
    p_a = &a[a_tlr + i - 1][a_tlc];
    p_b = &b[b_tlr + i - 1][b_tlc];

    memcpy(p_b, p_a, mem_size);
  }

  return (UT_OK);
}


/*******************************************************************************
COPY_CMAT_SECTION
Copy a rectangular section of one unsigned char matrix into a part of another 
matrix by means of direct memory copy.  (a_tlr, a_tlc) is the top
left corner of the section in the source matrix.  (b_tlr, b_tlc) is
the top left corner of the section in the destination matrix.  (nr, nc) are
the dimensions of the matrix section to be copied.
RG
*******************************************************************************/
int copy_cmat_section(unsigned char **a, unsigned char **b, int a_tlr, 
                       int a_tlc, int b_tlr, int b_tlc, int nr, int nc)
{
  int           i;
  int           mem_size;
  unsigned char *p_a, *p_b;

  mem_size = nc * sizeof(unsigned char);

  /* copy a section of A into B */

  for (i = 1; i <= nr; i++) {
    p_a = &a[a_tlr + i - 1][a_tlc];
    p_b = &b[b_tlr + i - 1][b_tlc];

    memcpy(p_b, p_a, mem_size);
  }

  return (UT_OK);
}


/*******************************************************************************
COPY_DMAT_SECTION
Copy a rectangular section of one double matrix into a part of another 
matrix by means of direct memory copy.  (a_tlr, a_tlc) is the top
left corner of the section in the source matrix.  (b_tlr, b_tlc) is
the top left corner of the section in the destination matrix.  (nr, nc)
are the dimensions of the matrix section to be copied.
RG
*******************************************************************************/
int copy_dmat_section(double **a, double **b, int a_tlr, 
                       int a_tlc, int b_tlr, int b_tlc, int nr, int nc)
{
  int    i;
  int    mem_size;
  double *p_a, *p_b;

  mem_size = nc * sizeof(double);

  /* copy a section of A into B */

  for (i = 1; i <= nr; i++) {
    p_a = &a[a_tlr + i - 1][a_tlc];
    p_b = &b[b_tlr + i - 1][b_tlc];

    memcpy(p_b, p_a, mem_size);
  }

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_ROW_MAT
Overwrite the specified row of the matrix with a given vector of values.
The inverse operation of grab_row().
RG
*******************************************************************************/
int overwrite_row_mat(float **mat, int index, int nc, float *vec)
{
  int   mem_size;
  float *p_vec;
  float *p_mat;

  mem_size = nc * sizeof(float);

  p_vec = &vec[1];
  p_mat = mat[index];

  memcpy(p_mat, p_vec, mem_size);

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_ROW_DMAT
Overwrite the specified row of a matrix of doubles with a given vector of 
values.  The inverse operation of grab_row().
RG
*******************************************************************************/
int overwrite_row_dmat(double **mat, int index, int nc, double *vec)
{
  int     mem_size;
  double *p_vec;
  double *p_mat;

  mem_size = nc * sizeof(double);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_mat, p_vec, mem_size);

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_ROW_IMAT
Overwrite the specified row of a matrix of integers with a given vector of 
values.  The inverse operation of grab_row().
RG
*******************************************************************************/
int overwrite_row_imat(int **mat, int index, int nc, int *vec)
{
  int  mem_size;
  int *p_vec;
  int *p_mat;

  mem_size = nc * sizeof(int);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_mat, p_vec, mem_size);

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_ROW_CMAT
Overwrite the specified row of a matrix of unsigned chars with a given vector 
of values.  The inverse operation of grab_row().
RG
*******************************************************************************/
int overwrite_row_cmat(unsigned char **mat, int index, int nc, 
		       unsigned char *vec)
{
  int  mem_size;
  unsigned char *p_vec;
  unsigned char *p_mat;

  mem_size = nc * sizeof(unsigned char);

  p_vec = &vec[1];
  p_mat = &mat[index][1];

  memcpy(p_mat, p_vec, mem_size);

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_COL_MAT
Overwrite the specified column of the matrix with a given vector of values.
The inverse operation of grab_col().
AG
*******************************************************************************/
int overwrite_col_mat(float **mat, int index, int nr, float *vec)
{
  int j;

  for (j=1; j <= nr; j++)
    mat[j][index] = vec[j];

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_COL_DMAT
Overwrite the specified column of a matrix of doubles with a given vector of 
values.  The inverse operation of grab_col().
AG,RG
*******************************************************************************/
int overwrite_col_dmat(double **mat, int index, int nr, double *vec)
{
  int j;

  for (j=1; j <= nr; j++)
    mat[j][index] = vec[j];

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_COL_IMAT
Overwrite the specified column of a matrix of integers with a given vector of 
values.  The inverse operation of grab_col().
AG,RG
*******************************************************************************/
int overwrite_col_imat(int **mat, int index, int nr, int *vec)
{
  int j;

  for (j=1; j <= nr; j++)
    mat[j][index] = vec[j];

  return (UT_OK);
}


/*******************************************************************************
OVERWRITE_COL_CMAT
Overwrite the specified column of a matrix of unsigned chars with a given 
vector of values.  The inverse operation of grab_col().
AG,RG
*******************************************************************************/
int overwrite_col_cmat(unsigned char **mat, int index, int nr, 
                       unsigned char *vec)
{
  int j;

  for (j=1; j <= nr; j++)
    mat[j][index] = vec[j];

  return (UT_OK);
}



/*******************************************************************************
SET_MAT
Sets each of the values of a matrix to a specified number.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_mat(float **m, int nr, int nc, float constant)
{
  int     i,j;

  /* change each element of the vector */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_DMAT
Sets each of the values of a matrix of doubles to a specified number.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_dmat(double **m, int nr, int nc, double constant)
{
  int     i,j;

  /* change each element of the vector */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_IMAT
Sets each of the values of a matrix of integers to a specified number.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_imat(int **m, int nr, int nc, int constant)
{
  int     i,j;

  /* change each element of the vector */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_CMAT
Sets each of the values of a matrix of unsigned chars to a specified number.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_cmat(unsigned char **m, int nr, int nc, unsigned char constant)
{
  int     i,j;

  /* change each element of the vector */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_MAT
Sets each of the values of a matrix of floats to a specified number by
setting memory rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
*******************************************************************************/
int fast_set_mat(float **m, int nr, int nc, float constant)
{
  memset(&m[1][1], constant, (nr * nc) * sizeof(float));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_DMAT
Sets each of the values of a matrix of doubles to a specified number by
setting memory rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
*******************************************************************************/
int fast_set_dmat(double **m, int nr, int nc, double constant)
{
  memset(&m[1][1], constant, (nr * nc) * sizeof(double));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_IMAT
Sets each of the values of a matrix of integers to a specified number by
setting memory rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
*******************************************************************************/
int fast_set_imat(int **m, int nr, int nc, int constant)
{
  memset(&m[1][1], constant, (nr * nc) * sizeof(int));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_CMAT
Sets each of the values of a matrix of unsigned chars to a specified number by
setting memory rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
*******************************************************************************/
int fast_set_cmat(unsigned char **m, int nr, int nc, unsigned char constant)
{
  memset(&m[1][1], constant, (nr * nc) * sizeof(unsigned char));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_MAT
Sets each of the values of a matrix of floats to zero by setting memory rather 
than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
AG
*******************************************************************************/
int fast_zero_mat(float **m, int nr, int nc)
{
  memset(&m[1][1], 0, (nr * nc) * sizeof(float));

  return (UT_OK);
}

/*******************************************************************************
FAST_ZERO_DMAT
Sets each of the values of a matrix of doubles to zero by setting memory 
rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
AG
*******************************************************************************/
int fast_zero_dmat(double **m, int nr, int nc)
{
  memset(&m[1][1], 0, (nr * nc) * sizeof(double));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_IMAT
Sets each of the values of a matrix of integers to zero by setting memory 
rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
AG
*******************************************************************************/
int fast_zero_imat(int **m, int nr, int nc)
{
  memset(&m[1][1], 0, (nr * nc) * sizeof(int));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_CMAT
Sets each of the values of a matrix of unsigned chars to zero by setting memory 
rather than a normal loop.

nr is the number of rows of the matrix.
nc is the number of cols of the matrix.
m is the input matrix.
AG
*******************************************************************************/
int fast_zero_cmat(unsigned char **m, int nr, int nc)
{
  memset(&m[1][1], 0, (nr * nc) * sizeof(unsigned char));

  return (UT_OK);
}


/*******************************************************************************
SET_VEC
Sets each of the values of a vector to a specified number.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_vec(float *v, int n,  float constant)
{
  int     i;

  /* change each element of the vector */
  for (i = 1; i <= n; i++) {
    v[i] = constant;
  }

  return (UT_OK);
}


/*******************************************************************************
SET_DVEC
Sets each of the values of a vector of doubles to a specified number.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_dvec(double *v, int n, double constant)
{
  int     i;

  /* change each element of the vector */
  for (i = 1; i <= n; i++) {
    v[i] = constant;
  }

  return (UT_OK);
}


/*******************************************************************************
SET_IVEC
Sets each of the values of a vector of integers to a specified number.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_ivec(int *v, int n, int constant)
{
  int     i;

  /* change each element of the vector */
  for (i = 1; i <= n; i++) {
    v[i] = constant;
  }

  return (UT_OK);
}


/*******************************************************************************
SET_CVEC
Sets each of the values of a vector of unsigned chars to a specified number.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int set_cvec(unsigned char *v, int n, unsigned char constant)
{
  int     i;

  /* change each element of the vector */
  for (i = 1; i <= n; i++) {
    v[i] = constant;
  }

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_VEC
Sets each of the values of a vector of floats to a specified number by
setting the memory values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int fast_set_vec(float *v, int n, float constant)
{
  memset(&v[1], constant, n * sizeof(float));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_DVEC
Sets each of the values of a vector of doubles to a specified number by
setting the memory values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int fast_set_dvec(double *v, int n, double constant)
{
  memset(&v[1], constant, n * sizeof(double));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_IVEC
Sets each of the values of a vector of integers to a specified number by
setting the memory values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int fast_set_ivec(int *v, int n, int constant)
{
  memset(&v[1], constant, n * sizeof(int));

  return (UT_OK);
}


/*******************************************************************************
FAST_SET_CVEC
Sets each of the values of a vector of unsigned chars to a specified number by
setting the memory values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
constant is the amount to set all the values to.
AG
*******************************************************************************/
int fast_set_cvec(unsigned char *v, int n, unsigned char constant)
{
  memset(&v[1], constant, n * sizeof(unsigned char));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_VEC
Sets each of the values of a vector of floats to zero by setting the memory 
values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
AG
*******************************************************************************/
int fast_zero_vec(float *v, int n)
{
  memset(&v[1], 0, n * sizeof(float));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_DVEC
Sets each of the values of a vector of doubles to zero by setting the memory 
values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
AG
*******************************************************************************/
int fast_zero_dvec(double *v, int n)
{
  memset(&v[1], 0, n * sizeof(double));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_IVEC
Sets each of the values of a vector of integers to zero by setting the memory 
values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
AG
*******************************************************************************/
int fast_zero_ivec(int *v, int n)
{
  memset(&v[1], 0, n * sizeof(int));

  return (UT_OK);
}


/*******************************************************************************
FAST_ZERO_CVEC
Sets each of the values of a vector of unsigned chars to zero by setting the 
memory values directly rather than by using a loop.

nc is the length of the vector (num. cols).
v is the input vector.
AG
*******************************************************************************/
int fast_zero_cvec(unsigned char *v, int n)
{
  memset(&v[1], 0, n * sizeof(unsigned char));

  return (UT_OK);
}


/*******************************************************************************
SET_DIAG_MAT
Set the diagonal of a square matrix to a given constant, and set the other
elements to zero.
RG
*******************************************************************************/
int set_diag_mat(float **m, int n, float constant)
{
  float *p, *p_end;
  int    n_plus_1;

  n_plus_1 = n + 1;

  fast_zero_mat(m, n, n);

  p_end = &m[n][n];

  for (p = &m[1][1]; p <= p_end; p += n_plus_1)
    *p = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_DIAG_DMAT
Set the diagonal of a square matrix of doubles to a given constant, and set 
the other elements to zero.
RG
*******************************************************************************/
int set_diag_dmat(double **m, int n, double constant)
{
  double *p, *p_end;
  int     n_plus_1;

  n_plus_1 = n + 1;

  fast_zero_dmat(m, n, n);

  p_end = &m[n][n];

  for (p = &m[1][1]; p <= p_end; p += n_plus_1)
    *p = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_DIAG_IMAT
Set the diagonal of a square matrix of integers to a given constant, and set 
the other elements to zero.
RG
*******************************************************************************/
int set_diag_imat(int **m, int n, int constant)
{
  int *p, *p_end;
  int  n_plus_1;

  n_plus_1 = n + 1;

  fast_zero_imat(m, n, n);

  p_end = &m[n][n];

  for (p = &m[1][1]; p <= p_end; p += n_plus_1)
    *p = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_DIAG_CMAT
Set the diagonal of a square matrix of unsigned chars to a given constant, and 
set the other elements to zero.
RG
*******************************************************************************/
int set_diag_cmat(unsigned char **m, int n, unsigned char constant)
{
  unsigned char *p, *p_end;
  int            n_plus_1;

  n_plus_1 = n + 1;

  fast_zero_cmat(m, n, n);

  p_end = &m[n][n];

  for (p = &m[1][1]; p <= p_end; p += n_plus_1)
    *p = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_BAND_DIAG_MAT
Set the diagonal band of a square matrix to a given constant, and set the other
elements to zero.  The diagonal band is defined as the diagonal plus a number
of strips of elements adjacent to it, on either side of it; this number is 
specified as an argument.
AG
*******************************************************************************/
int set_band_diag_mat(float **m, int n, int width, float constant)
{
  int i, j;

  fast_zero_mat(m, n, n);

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      if ((i >= j-width) && (i <= j+width))
        m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_BAND_DIAG_DMAT
Set the diagonal band of a square matrix of doubles to a given constant, and 
set the other elements to zero.  The diagonal band is defined as the diagonal 
plus a number of strips of elements adjacent to it, on either side of it; this 
number is specified as an argument.
AG, RG
*******************************************************************************/
int set_band_diag_dmat(double **m, int n, int width, double constant)
{
  int i, j;

  fast_zero_dmat(m, n, n);

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      if ((i >= j-width) && (i <= j+width))
        m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_BAND_DIAG_IMAT
Set the diagonal band of a square matrix of integers to a given constant, and 
set the other elements to zero.  The diagonal band is defined as the diagonal 
plus a number of strips of elements adjacent to it, on either side of it; this 
number is specified as an argument.
AG, RG
*******************************************************************************/
int set_band_diag_imat(int **m, int n, int width, int constant)
{
  int i, j;

  fast_zero_imat(m, n, n);

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      if ((i >= j-width) && (i <= j+width))
        m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
SET_BAND_DIAG_CMAT
Set the diagonal band of a square matrix of unsigned chars to a given constant, 
and set the other elements to zero.  The diagonal band is defined as the 
diagonal plus a number of strips of elements adjacent to it, on either side of 
it; this number is specified as an argument.
AG, RG
*******************************************************************************/
int set_band_diag_cmat(unsigned char **m, int n, int width, 
		       unsigned char constant)
{
  int i, j;

  fast_zero_cmat(m, n, n);

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      if ((i >= j-width) && (i <= j+width))
        m[i][j] = constant;

  return (UT_OK);
}


/*******************************************************************************
COPY_MAT
Copy a matrix "a" into a matrix "b".
RG, AG
*******************************************************************************/
int copy_mat(float **a, float **b, int nr, int nc)
{
  int   mem_size;
  float *p_a, *p_b;

  p_a = &a[1][1];
  p_b = &b[1][1];

  mem_size = nr * nc * sizeof(float);

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_IMAT
Copy a matrix of ints "a" into a matrix "b".
RG,AG
*******************************************************************************/
int copy_imat(int **a, int **b, int nr, int nc)
{
  int   mem_size;
  int   *p_a, *p_b;

  p_a = &a[1][1];
  p_b = &b[1][1];

  mem_size = nr * nc * sizeof(int);

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_DMAT
Copy a matrix of doubles "a" into a matrix "b".
RG,AG
*******************************************************************************/
int copy_dmat(double **a, double **b, int nr, int nc)
{
  int     mem_size;
  double  *p_a, *p_b;

  p_a = &a[1][1];
  p_b = &b[1][1];

  mem_size = nr * nc * sizeof(double);

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_DMAT_SLOW
Copy a matrix of doubles "a" into a matrix "b".
RG
*******************************************************************************/
int copy_dmat_slow(double **a, double **b, int nr, int nc)
{
  int     i, j;

  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
      b[i][j] = a[i][j];

  return (UT_OK);
}

/*******************************************************************************
COPY_CMAT
Copy a matrix of unsigned chars "a" into a matrix "b".
RG,AG
*******************************************************************************/
int copy_cmat(unsigned char **a, unsigned char **b, int nr, int nc)
{
  int            mem_size;
  unsigned char  *p_a, *p_b;

  p_a = &a[1][1];
  p_b = &b[1][1];

  mem_size = nr * nc * sizeof(unsigned char);

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_VEC
Copy a vector "a" into a vector "b".
RG,AG
*******************************************************************************/
int copy_vec(float *a, float *b, int nc)
{
  int    mem_size;
  float  *p_a, *p_b;

  mem_size = nc * sizeof(float);

  p_a = &a[1];
  p_b = &b[1];

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_IVEC
Copy a vector of ints "a" into a ints "b".
RG,AG
*******************************************************************************/
int copy_ivec(int *a, int *b, int nc)
{
  int    mem_size;
  int    *p_a, *p_b;

  mem_size = nc * sizeof(int);

  p_a = &a[1];
  p_b = &b[1];

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_DVEC
Copy a vector of doubles "a" into a vector "b".
RG,AG
*******************************************************************************/
int copy_dvec(double *a, double *b, int nc)
{
  int    mem_size;
  double *p_a, *p_b;

  mem_size = nc * sizeof(double);

  p_a = &a[1];
  p_b = &b[1];

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_CVEC
Copy a vector of unsigned chars "a" into a vector "b".
RG,AG
*******************************************************************************/
int copy_cvec(unsigned char *a, unsigned char *b, int nc)
{
  int           mem_size;
  unsigned char *p_a, *p_b;

  mem_size = nc * sizeof(unsigned char);

  p_a = &a[1];
  p_b = &b[1];

  memcpy(p_b, p_a, mem_size);

  return (UT_OK);
}

/*******************************************************************************
COPY_SUB_VEC_TO_VEC
Copy a vector into a larger vector, starting at a given position in the larger
vector.
RG
*******************************************************************************/
int copy_sub_vec_to_vec(float *sub_vec, float *vec, int sub_vec_size, 
                        int vec_size, int vec_start_pos)
{
  float *p;
  float *p_sub;
  float *p_end;

  /* assign pointer to last element of vec */
  p_end = &vec[vec_size];

  for (p = &vec[vec_start_pos], p_sub = &sub_vec[1]; p <= p_end; p++, p_sub++)
    *p = *p_sub;

  return (UT_OK);
}


/*******************************************************************************
COPY_SUB_DVEC_TO_DVEC
Copy a vector of doubles into a larger vector, starting at a given position in 
the larger vector.
RG
*******************************************************************************/
int copy_sub_dvec_to_dvec(double *sub_vec, double *vec, int sub_vec_size, 
                          int vec_size, int vec_start_pos)
{
  double *p;
  double *p_sub;
  double *p_end;

  /* assign pointer to last element of vec */
  p_end = &vec[vec_size];

  for (p = &vec[vec_start_pos], p_sub = &sub_vec[1]; p <= p_end; p++, p_sub++)
    *p = *p_sub;

  return (UT_OK);
}


/*******************************************************************************
COPY_SUB_IVEC_TO_IVEC
Copy a vector of integers into a larger vector, starting at a given position in 
the larger vector.
RG
*******************************************************************************/
int copy_sub_ivec_to_ivec(int *sub_vec, int *vec, int sub_vec_size, 
                          int vec_size, int vec_start_pos)
{
  int *p;
  int *p_sub;
  int *p_end;

  /* assign pointer to last element of vec */
  p_end = &vec[vec_size];

  for (p = &vec[vec_start_pos], p_sub = &sub_vec[1]; p <= p_end; p++, p_sub++)
    *p = *p_sub;

  return (UT_OK);
}


/*******************************************************************************
COPY_SUB_CVEC_TO_CVEC
Copy a vector of unsigned chars into a larger vector, starting at a given 
position in the larger vector.
RG
*******************************************************************************/
int copy_sub_cvec_to_cvec(unsigned char *sub_vec, unsigned char *vec, 
                          int sub_vec_size, int vec_size, int vec_start_pos)
{
  unsigned char *p;
  unsigned char *p_sub;
  unsigned char *p_end;

  /* assign pointer to last element of vec */
  p_end = &vec[vec_size];

  for (p = &vec[vec_start_pos], p_sub = &sub_vec[1]; p <= p_end; p++, p_sub++)
    *p = *p_sub;

  return (UT_OK);
}


/*******************************************************************************
SET_OF_MATRICES
Allocate a structure which holds many matrices.
AG
*******************************************************************************/
float ***set_of_matrices(int num_mats, int num_rows, int num_cols)
{
  int i;
  float ***mat_set;

  /* accomodate the NR indexing convention */
  mat_set = (float***) malloc_return_if_fail((num_mats+1) * sizeof(float**),
                                              (float***)NULL);

  mat_set[0] = (float**)NULL;
  for (i = 1; i <= num_mats; i++)
    mat_set[i] = NR_matrix(1, num_rows, 1, num_cols);

  return (mat_set);
}

/*******************************************************************************
SET_OF_IMATRICES
Allocate a structure which holds many matrices.
AG
*******************************************************************************/
int ***set_of_imatrices(int num_mats, int num_rows, int num_cols)
{
  int i;
  int ***mat_set;

  /* accomodate the NR indexing convention */
  mat_set = (int***) malloc_return_if_fail((num_mats+1) * sizeof(int**),
                                              (int***)NULL);

  mat_set[0] = (int**)NULL;
  for (i = 1; i <= num_mats; i++)
    mat_set[i] = NR_imatrix(1, num_rows, 1, num_cols);

  return (mat_set);
}

/*******************************************************************************
SET_OF_DMATRICES
Allocate a structure which holds many matrices.
AG
*******************************************************************************/
double ***set_of_dmatrices(int num_mats, int num_rows, int num_cols)
{
  int i;
  double ***mat_set;

  /* accomodate the NR indexing convention */
  mat_set = (double***) malloc_return_if_fail((num_mats+1) * sizeof(double**),
                                              (double***)NULL);

  mat_set[0] = (double**)NULL;
  for (i = 1; i <= num_mats; i++)
    mat_set[i] = NR_dmatrix(1, num_rows, 1, num_cols);

  return (mat_set);
}

/*******************************************************************************
FREE_SET_OF_MATRICES
De-allocate a structure which holds many matrices.
AG
*******************************************************************************/
void free_set_of_matrices(float ***mat_set, int num_mats, int num_rows, 
                          int num_cols)
{
  int i;

  for (i = 1; i <= num_mats; i++)
    NR_free_matrix(mat_set[i], 1, num_rows, 1, num_cols);

  free(mat_set);

  /* no error code returned, since this mirrors the NR NR_free_...() 
     functions */
}

/*******************************************************************************
FREE_SET_OF_IMATRICES
De-allocate a structure which holds many matrices.
AG
*******************************************************************************/
void free_set_of_imatrices(int ***mat_set, int num_mats, int num_rows, 
                           int num_cols)
{
  int i;

  for (i = 1; i <= num_mats; i++)
    NR_free_imatrix(mat_set[i], 1, num_rows, 1, num_cols);

  free(mat_set);

  /* no error code returned, since this mirrors the NR NR_free_...() 
     functions */
}


/*******************************************************************************
FREE_SET_OF_DMATRICES
De-allocate a structure which holds many matrices.
AG
*******************************************************************************/
void free_set_of_dmatrices(double ***mat_set, int num_mats, int num_rows, 
                          int num_cols)
{
  int i;

  for (i = 1; i <= num_mats; i++)
    NR_free_dmatrix(mat_set[i], 1, num_rows, 1, num_cols);

  free(mat_set);

  /* no error code returned, since this mirrors the NR NR_free_...() 
     functions */
}

/*******************************************************************************
SET_OF_SETS_OF_MATRICES
Allocate a structure which holds many sets of matrices.
AG
*******************************************************************************/
float ****set_of_sets_of_matrices(int num_sets, int num_mats, int num_rows, 
                                  int num_cols)
{
  int i;
  float ****mat_set_set;

  /* accomodate the NR indexing convention */
  mat_set_set = (float****) malloc_return_if_fail((num_sets+1) * 
                                                   sizeof(float***), 
                                                   (float****)NULL);

  mat_set_set[0] = (float***)NULL;
  for (i = 1; i <= num_sets; i++)
    mat_set_set[i] = set_of_matrices(num_mats, num_rows, num_cols);

  return (mat_set_set);
}

/*******************************************************************************
SET_OF_SETS_OF_DMATRICES
Allocate a structure which holds many sets of matrices.
AG
*******************************************************************************/
double ****set_of_sets_of_dmatrices(int num_sets, int num_mats, int num_rows, 
                                   int num_cols)
{
  int i;
  double ****mat_set_set;

  /* accomodate the NR indexing convention */
  mat_set_set = (double****) malloc_return_if_fail((num_sets+1) * 
                                                   sizeof(double***), 
                                                   (double****)NULL);

  mat_set_set[0] = (double***)NULL;
  for (i = 1; i <= num_sets; i++)
    mat_set_set[i] = set_of_dmatrices(num_mats, num_rows, num_cols);

  return (mat_set_set);
}

/*******************************************************************************
FREE_SET_OF_SETS_OF_MATRICES
De-allocate a structure which holds many sets of matrices.
AG
*******************************************************************************/
void free_set_of_sets_of_matrices(float ****mat_set_set, int num_sets, 
                                  int num_mats, int num_rows, int num_cols)
{
  int i;

  for (i = 1; i <= num_sets; i++)
    free_set_of_matrices(mat_set_set[i], num_mats, num_rows, num_cols);

  free(mat_set_set);

  /* no error code returned, since this mirrors the NR NR_free_...() functions */
}

/*******************************************************************************
FREE_SET_OF_SETS_OF_DMATRICES
De-allocate a structure which holds many sets of matrices.
AG
*******************************************************************************/
void free_set_of_sets_of_dmatrices(double ****mat_set_set, int num_sets, 
                                  int num_mats, int num_rows, int num_cols)
{
  int i;

  for (i = 1; i <= num_sets; i++)
    free_set_of_dmatrices(mat_set_set[i], num_mats, num_rows, num_cols);

  free(mat_set_set);

  /* no error code returned, since this mirrors the NR NR_free_...() functions */
}
