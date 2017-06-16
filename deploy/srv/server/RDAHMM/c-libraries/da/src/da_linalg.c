/*******************************************************************************
MODULE NAME
da_linalg

ONE-LINE SYNOPSIS
General functions related to linear algebra.

SCOPE OF THIS MODULE
Any functions relating directly to other linear algebra concepts which have 
representative modules in this library should go in the appropriate module.  
Functions that apply more generally are intended to go here.  For instance,
generating a matrix of random numbers is more directly related to random
numbers than to matrix operations, so it belongs in da_random.

SEE ALSO
Because the definition of this module is quite broad, there is some potential
overlap with several other modules in this library.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/pca, AG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_linalg.c,v 1.1.1.1 2001/05/21 23:36:14 diane Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_linalg.c,v $
 * Revision 1.1.1.1  2001/05/21 23:36:14  diane
 * Import of various mls & woldlab source packages
 *
 * Revision 1.45  2000/11/22 01:02:43  granat
 * added some extra matrix inverse routines for symmetric matrices
 *
 * Revision 1.44  2000/11/06 19:23:14  granat
 * added some routines for inversion of symmetric positive definite matrices
 * that use the cholesky decomposition
 *
 * Revision 1.43  2000/04/04 23:40:06  granat
 * added Jacobi preconditioned version of inverse function
 *
 * Revision 1.42  1998/07/02 01:15:29  granat
 * added fast matrix multiply
 *
 * Revision 1.41  1998/06/29 22:11:24  granat
 * added LSQR method of Paige and Saunders
 *
 * Revision 1.40  1998/05/07 23:56:43  granat
 * moved many functions to da_util.c
 *
 * Revision 1.39  1998/05/01 17:37:03  granat
 * added versions of many functions to handle other data types
 * added many new functions
 * changed some function names to match convention
 * some minor reformatting changes
 * edited many functions to make them faster
 *
 * Revision 1.38  1997/10/21 14:32:23  granat
 * fixed some bugs
 * added versions for different data types of some functions
 * changed some argument orders to agree with new standard
 *
 * Revision 1.37  1997/09/10 14:50:22  granat
 * changed many functions to their faster versions, added more comments
 *
 * Revision 1.36  1997/09/04 20:21:59  granat
 * added variants of sum_mat
 * added fast_copy_mat_section and variants
 *
 * Revision 1.35  1997/08/11 18:31:19  granat
 * added set_imat
 *
 * Revision 1.34  1997/06/20 22:09:41  granat
 * fixed transpose memory handling, some penny optimization, some cosmetic changes
 *
 * Revision 1.33  1997/06/05 18:53:57  granat
 * added flip_vector, changed prototypes to match with conventions, some cosmetic changes
 *
 * Revision 1.32  1997/06/02 15:37:02  granat
 * changed to use new NR naming convention
 *
 * Revision 1.31  1997/05/06 22:23:05  agray
 * added some things from dp cooltool.
 *
 * Revision 1.30  1997/05/06 19:04:15  granat
 * added function copy_mat_section
 *
 * Revision 1.29  1997/04/05 19:09:34  granat
 * added sum_mat() and norm_sum_mat()
 * adjusted functions so that they all follow input parameter conventions
 *
 * Revision 1.28  1997/04/04 23:09:12  granat
 * changed many functions so that they use pointer arithmatic
 *
 * Revision 1.27  1997/04/04 19:00:30  granat
 * added comments to transpose_in_situ_ functions and flip_ functions
 *
 * Revision 1.26  1997/03/27 18:09:48  granat
 * added transpose_in_situ_alloc_matrix(), transpose_in_situ_sqr_matrix(),
 * flip_left_right_matrix(), flip_top_bottom_matrix()
 *
 * Revision 1.25  1997/03/15 17:50:49  granat
 * Added flip_left_right_matrix() and flip_top_bottom_matrix()
 *
 * Revision 1.24  1997/03/14 19:26:06  agray
 * fixed minor bug in set_of_sets_of_matrices().
 *
 * Revision 1.23  1997/02/18 22:25:15  granat
 * fixed subtle error in transpose
 *
 * Revision 1.22  1997/02/14 15:00:37  granat
 * fixed mistake in transpose_matrix
 *
 * Revision 1.21  1997/02/14 00:25:15  granat
 * fixed transpose function in pca
 *
 * Revision 1.20  1997/02/14 00:06:59  granat
 * fixed typo
 *
 * Revision 1.19  1997/02/14 00:01:25  granat
 * Changed transpose funtion to transpose_matrix
 * speeded up transpose algorithm
 * added transpose_imatrix to transpose matrices of integers
 *
 * Revision 1.18  1997/01/29 21:45:45  agray
 * new formatting, cleaning up debugging output using ut_output,
 * cleaning up memory allocation with ut_memory.
 *
 * Revision 1.17  1996/10/31 02:16:35  agray
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.16  1996/09/24 17:54:42  agray
 * changed svdcmp2() to da_svdcmp2().
 *
 * Revision 1.15  1996/07/30 23:17:30  agray
 * changed heuristic for number of iterations in pca().
 *
 * Revision 1.14  1996/07/17 20:43:37  agray
 * cosmetic.
 *
 * Revision 1.13  1996/07/13 01:07:59  agray
 * moved out print/write_row/col/irow/icol() to da_data module
 *
 * Revision 1.12  1996/07/11 18:03:45  agray
 * moved out read_gauss_parms(), write_gauss_parms() to da_prob module;
 * added add_mat(), subtract_mat(), invert_mat_copy(), det_copy(),
 * restrict_illcond_matrix(), scalar_mult/div/add/subtract_mat/vec(), set_mat(),
 * set_vec(), mult/div_vec_elt(), sum_vec(), max/min_vec(), arg_max/min_vec(),
 * copy_vec().
 *
 * Revision 1.11  1996/04/09 02:49:36  agray
 * added print_irow(), print_icol(), write_irow(), write_icol().
 * ag
 *
 * Revision 1.10  1996/03/01  00:22:16  agray
 * moved read_data() and write_data() to da_data module.
 * changed %f to %g in all format strings.
 * ag
 *
 * Revision 1.9  1996/02/29  02:33:32  agray
 * moved write_bin_matrix() and read_bin_matrix() to da_data module
 * ag
 *
 * Revision 1.8  1996/02/29 00:54:12  agray
 * changed pca() to call svdcmp2() instead of svdcmp(); added svdcmp2()
 * to NR lib. to make the number of svd iterations a controllable
 * parameter.
 * ag
 *
 * Revision 1.7  1996/02/28 03:57:52  agray
 * put in checks for all memory allocations
 * AG
 *
 * Revision 1.6  1996/02/21 05:16:23  agray
 * NR_free_matrix() in pca().
 * ag
 *
 * Revision 1.5  1996/02/21  04:08:46  agray
 * added nr.h inclusion
 * ag
 *
 * Revision 1.4  1996/02/21  04:02:21  agray
 * moved some general notes on coding conventions into the file conventions.txt
 * in the /doc directory.
 * ag
 *
 * Revision 1.3  1996/02/21  03:51:36  agray
 * added pca().
 * ag
 *
 * Revision 1.2  1996/02/21  00:37:20  agray
 * added write_bin_matrix() and read_bin_matrix()
 * ag
 *
 * Revision 1.1  1996/02/06  03:27:47  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

/* UT library */
#include "ut_error.h"
#include "ut_math.h"
#include "ut_memory.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_util.h"
#include "da_nrhacks.h"

/* this module's header */
#include "da_linalg.h"


/*******************************************************************************
NORM_VEC
Returns the norm of a vector.
RG
*******************************************************************************/
float norm_vec(float *v, int n)
{
  float *p;
  float *p_end;
  float result = 0.0;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    result += NR_sqr(*p);

  return (sqrt((double) result));
}


/*******************************************************************************
NORM_DVEC
Returns the norm of a vector of doubles.
RG
*******************************************************************************/
double norm_dvec(double *v, int n)
{
  double *p;
  double *p_end;
  double result = 0.0;

  /* assign pointer to last element */
  p_end = &v[n];

  /* accumulate square of the vector */
  for (p = &v[1]; p <= p_end; p++)
    result += NR_sqr(*p);

  return (sqrt(result));
}


/*******************************************************************************
FROBENIUS_NORM_MAT
Returns the frobenius norm of a matrix.
RG
*******************************************************************************/
float frobenius_norm_mat(float **m, int nr, int nc)
{
  float *p;
  float *p_end;
  float result = 0.0;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];
  
  /* calculate the frobenius norm */
  for (p = &m[1][1]; p <= p_end; p++)
    result += NR_sqr(*p);
  
  return (sqrt(result));
}


/*******************************************************************************
FROBENIUS_NORM_DMAT
Returns the frobenius norm of a matrix.
RG
*******************************************************************************/
double frobenius_norm_dmat(double **m, int nr, int nc)
{
  double *p;
  double *p_end;
  double result = 0.0;
  
  /* assign pointer to last element */
  p_end = &m[nr][nc];
  
  /* calculate the frobenius norm */
  for (p = &m[1][1]; p <= p_end; p++)
    result += NR_sqr(*p);
  
  return (sqrt(result));
}


/*******************************************************************************
TRANSPOSE_MATRIX
Create the transpose of a matrix of floats in a separate matrix.
RG
*******************************************************************************/
int transpose_matrix(float **a, int nr, int nc, float **a_trans)
{
  float  *p;
  float  *p_trans, *p_end_trans;
  float  *p_col, *p_end_col;

  /* assign pointer to last element of transpose*/
  p_end_trans = &a_trans[nc][nr];
  
  /* assign pointer to last element of first transpose row */
  p_end_col = &a_trans[1][nr];

  /* start pointer at first element */
  p = &a[1][1];

  /* increment pointer to the original matrix element by element, but
     increment pointer to the transpose along each column */
  for (p_col = &a_trans[1][1]; p_col <= p_end_col; p_col++)
    for (p_trans = p_col; p_trans <= p_end_trans; p_trans += nr, p++)
      *p_trans = *p;

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_DMATRIX
Create the transpose of a matrix of doubles in a separate matrix.
RG
*******************************************************************************/
int transpose_dmatrix(double **a, int nr, int nc, double **a_trans)
{
  double  *p;
  double  *p_trans, *p_end_trans;
  double  *p_col, *p_end_col;

  /* assign pointer to last element of transpose*/
  p_end_trans = &a_trans[nc][nr];
  
  /* assign pointer to last element of first transpose row */
  p_end_col = &a_trans[1][nr];

  /* start pointer at first element */
  p = &a[1][1];

  /* increment pointer to the original matrix element by element, but
     increment pointer to the transpose along each column */
  for (p_col = &a_trans[1][1]; p_col <= p_end_col; p_col++)
    for (p_trans = p_col; p_trans <= p_end_trans; p_trans += nr, p++)
      *p_trans = *p;

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IMATRIX
Create the transpose of a matrix of doubles in a separate matrix.
RG
*******************************************************************************/
int transpose_imatrix(int **a, int nr, int nc, int **a_trans)
{
  int     *p;
  int     *p_trans, *p_end_trans;
  int     *p_col, *p_end_col;

  /* assign pointer to last element of transpose*/
  p_end_trans = &a_trans[nc][nr];
  
  /* assign pointer to last element of first transpose row */
  p_end_col = &a_trans[1][nr];

  /* start pointer at first element */
  p = &a[1][1];

  /* increment pointer to the original matrix element by element, but
     increment pointer to the transpose along each column */
  for (p_col = &a_trans[1][1]; p_col <= p_end_col; p_col++)
    for (p_trans = p_col; p_trans <= p_end_trans; p_trans += nr, p++)
      *p_trans = *p;

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_CMATRIX
Create the transpose of a matrix of unsigned chars in a separate matrix.
RG
*******************************************************************************/
int transpose_cmatrix(unsigned char **a, int nr, int nc, unsigned char **a_trans)
{
  unsigned char  *p;
  unsigned char  *p_trans, *p_end_trans;
  unsigned char  *p_col, *p_end_col;

  /* assign pointer to last element of transpose*/
  p_end_trans = &a_trans[nc][nr];
  
  /* assign pointer to last element of first transpose row */
  p_end_col = &a_trans[1][nr];

  /* start pointer at first element */
  p = &a[1][1];

  /* increment pointer to the original matrix element by element, but
     increment pointer to the transpose along each column */
  for (p_col = &a_trans[1][1]; p_col <= p_end_col; p_col++)
    for (p_trans = p_col; p_trans <= p_end_trans; p_trans += nr, p++)
      *p_trans = *p;

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_ALLOC_MATRIX
Create the transpose of a matrix of floats in the same memory space as the
original matrix.  Takes as one of its arguments a vector for temporary
storage space, which must be of length at least equal to the number of rows
in the input matrix.  Allocation and deallocation may be required to maintain
the appropriate matrix data structure.
RG
*******************************************************************************/
int transpose_in_situ_alloc_matrix(float ***a, int nr, int nc, float *temp_vect,
                                   char *mem_choice)
{
  int     i, j, k;   
  float  *p;          /* pointer into the input matrix */
  float  *p_temp;     /* pointer into temporary values */
  float **m;          /* pointer to new vector of pointers to rows */
  int     skip;
  int     nr_less1;
  int     nc_less1;

  /****
   * Algorithm outline:
   * For each row in the transpose, skip through the matrix, picking out
   * the appropriate element of each column, and copy them into temporary
   * storage.  Then compress the remain parts of the untransposed matrix
   * to fill in the gaps, leaving a space in memory.  In the space created 
   * copy the stored elements, thus creating a new row.
   ****/

  p = &(*a)[1][1];

  nr_less1 = nr-1;
  nc_less1 = nc-1;

  skip = nc;

  /* extract columns and compose into rows */

  for (i = 1; i <= nc_less1; i++) {
    p_temp = &temp_vect[1]; /* p_temp points into the temporary storage */

    /* move p along so the it points to every "skip"th element */
    /* (ie, the values in a column) and copy the values into the */
    /* temporary storage */

    for (j = 1; j <= nr; j++) {
      *p_temp = *p;
      p += skip;
      p_temp++;
    }

    /* move the pointer back into the array */

    p -= skip;

    /* account for effective smaller column size now that an element 
       is removed */

    skip--;

    /* pack the remaining untransposed data into the last part of memory */
    /* this leaves space in memory in front of the packed data */

    for (j = 1; j <= nr_less1; j++)
      for (k = 1; k <= skip; k++) {
        *p = *(p-j);
        p--;
      }

    /* set pointer back to the beginning of the new space */

    p -= nr_less1;
    p_temp = &temp_vect[1];

    /* copy the stored elements into the new space */

    for (j = 1; j <= nr; j++) {
      *p = *p_temp;
      p++;
      p_temp++;
    }
  }

  /* make necessary changes to new pointer scheme */

  if (!strcmp(mem_choice, "realloc") && (nr != nc)) {

    /* allocate more memory and reassign pointers */
 
    m=(float **) malloc((size_t) ((nc+NR_END)*sizeof(float*)));
    if (!m)
      NR_error("allocation failure 1 in NR_matrix()");
    m += NR_END;
    m--;

    m[1] = (*a)[1];
    for(i = 2; i <= nc; i++)
      m[i] = m[i-1]+nr;

    free(*a);

    *a = m;
  }
  else
    if (nr > nc) { 

      /* reassign pointers and free excess memory */

      for (i = 2; i <= nc; i++)
        (*a)[i] = (*a)[i-1]+nr;
      free(&(*a)[i]);
    }
    else if (nc > nr) { 

      /* allocate more memory and reassign pointers */

      m=(float **) malloc((size_t) ((nc+NR_END)*sizeof(float*)));
      if (!m) 
        NR_error("allocation failure 1 in NR_matrix()");
      m += NR_END;
      m--;

      m[1] = (*a)[1];
      for(i = 2; i <= nc; i++)
        m[i] = m[i-1]+nr;

      free(*a);

      *a = m;
    }

  /* nothing needs to be done if nr = nc */

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_ALLOC_DMATRIX
Create the transpose of a matrix of doubles in the same memory space as the
original matrix.  Takes as one of its arguments a vector for temporary
storage space, which must be of length at least equal to the number of rows
in the input matrix.  Allocation and deallocation may be required to maintain
the appropriate matrix data structure.
RG
*******************************************************************************/
int transpose_in_situ_alloc_dmatrix(double ***a, int nr, int nc, 
		                    double *temp_vect, char *mem_choice)
{
  int      i, j, k;   
  double  *p;          /* pointer into the input matrix */
  double  *p_temp;     /* pointer into temporary values */
  double **m;          /* pointer to new vector of pointers to rows */
  int      skip;
  int      nr_less1;
  int      nc_less1;

  /****
   * Algorithm outline:
   * For each row in the transpose, skip through the matrix, picking out
   * the appropriate element of each column, and copy them into temporary
   * storage.  Then compress the remain parts of the untransposed matrix
   * to fill in the gaps, leaving a space in memory.  In the space created 
   * copy the stored elements, thus creating a new row.
   ****/

  p = &(*a)[1][1];

  nr_less1 = nr-1;
  nc_less1 = nc-1;

  skip = nc;

  /* extract columns and compose into rows */

  for (i = 1; i <= nc_less1; i++) {
    p_temp = &temp_vect[1]; /* p_temp points into the temporary storage */

    /* move p along so the it points to every "skip"th element */
    /* (ie, the values in a column) and copy the values into the */
    /* temporary storage */

    for (j = 1; j <= nr; j++) {
      *p_temp = *p;
      p += skip;
      p_temp++;
    }

    /* move the pointer back into the array */

    p -= skip;

    /* account for effective smaller column size now that an element 
       is removed */

    skip--;

    /* pack the remaining untransposed data into the last part of memory */
    /* this leaves space in memory in front of the packed data */

    for (j = 1; j <= nr_less1; j++)
      for (k = 1; k <= skip; k++) {
        *p = *(p-j);
        p--;
      }

    /* set pointer back to the beginning of the new space */

    p -= nr_less1;
    p_temp = &temp_vect[1];

    /* copy the stored elements into the new space */

    for (j = 1; j <= nr; j++) {
      *p = *p_temp;
      p++;
      p_temp++;
    }
  }

  /* make necessary changes to new pointer scheme */

  if (!strcmp(mem_choice, "realloc") && (nr != nc)) {

    /* allocate more memory and reassign pointers */
 
    m=(double **) malloc((size_t) ((nc+NR_END)*sizeof(double*)));
    if (!m)
      NR_error("allocation failure 1 in NR_matrix()");
    m += NR_END;
    m--;

    m[1] = (*a)[1];
    for(i = 2; i <= nc; i++)
      m[i] = m[i-1]+nr;

    free(*a);

    *a = m;
  }
  else
    if (nr > nc) { 

      /* reassign pointers and free excess memory */

      for (i = 2; i <= nc; i++)
        (*a)[i] = (*a)[i-1]+nr;
      free(&(*a)[i]);
    }
    else if (nc > nr) { 

      /* allocate more memory and reassign pointers */

      m=(double **) malloc((size_t) ((nc+NR_END)*sizeof(double*)));
      if (!m) 
        NR_error("allocation failure 1 in NR_matrix()");
      m += NR_END;
      m--;

      m[1] = (*a)[1];
      for(i = 2; i <= nc; i++)
        m[i] = m[i-1]+nr;

      free(*a);

      *a = m;
    }

  /* nothing needs to be done if nr = nc */

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_ALLOC_IMATRIX
Create the transpose of a matrix of integers in the same memory space as the
original matrix.  Takes as one of its arguments a vector for temporary
storage space, which must be of length at least equal to the number of rows
in the input matrix.  Allocation and deallocation may be required to maintain
the appropriate matrix data structure.
RG
*******************************************************************************/
int transpose_in_situ_alloc_imatrix(int ***a, int nr, int nc, int *temp_vect, 
                                    char *mem_choice)
{
  int   i, j, k;   
  int  *p;          /* pointer into the input matrix */
  int  *p_temp;     /* pointer into temporary values */
  int **m;          /* pointer to new vector of pointers to rows */
  int   skip;
  int   nr_less1;
  int   nc_less1;

  /****
   * Algorithm outline:
   * For each row in the transpose, skip through the matrix, picking out
   * the appropriate element of each column, and copy them into temporary
   * storage.  Then compress the remain parts of the untransposed matrix
   * to fill in the gaps, leaving a space in memory.  In the space created 
   * copy the stored elements, thus creating a new row.
   ****/

  p = &(*a)[1][1];

  nr_less1 = nr-1;
  nc_less1 = nc-1;

  skip = nc;

  /* extract columns and compose into rows */

  for (i = 1; i <= nc_less1; i++) {
    p_temp = &temp_vect[1]; /* p_temp points into the temporary storage */

    /* move p along so the it points to every "skip"th element */
    /* (ie, the values in a column) and copy the values into the */
    /* temporary storage */

    for (j = 1; j <= nr; j++) {
      *p_temp = *p;
      p += skip;
      p_temp++;
    }

    /* move the pointer back into the array */

    p -= skip;

    /* account for effective smaller column size now that an element 
       is removed */

    skip--;

    /* pack the remaining untransposed data into the last part of memory */
    /* this leaves space in memory in front of the packed data */

    for (j = 1; j <= nr_less1; j++)
      for (k = 1; k <= skip; k++) {
        *p = *(p-j);
        p--;
      }

    /* set pointer back to the beginning of the new space */

    p -= nr_less1;
    p_temp = &temp_vect[1];

    /* copy the stored elements into the new space */

    for (j = 1; j <= nr; j++) {
      *p = *p_temp;
      p++;
      p_temp++;
    }
  }

  /* make necessary changes to new pointer scheme */

  if (!strcmp(mem_choice, "realloc") && (nr != nc)) {

    /* allocate more memory and reassign pointers */
 
    m=(int **) malloc((size_t) ((nc+NR_END)*sizeof(int*)));
    if (!m)
      NR_error("allocation failure 1 in NR_matrix()");
    m += NR_END;
    m--;

    m[1] = (*a)[1];
    for(i = 2; i <= nc; i++)
      m[i] = m[i-1]+nr;

    free(*a);

    *a = m;
  }
  else
    if (nr > nc) { 

      /* reassign pointers and free excess memory */

      for (i = 2; i <= nc; i++)
        (*a)[i] = (*a)[i-1]+nr;
      free(&(*a)[i]);
    }
    else if (nc > nr) { 

      /* allocate more memory and reassign pointers */

      m=(int **) malloc((size_t) ((nc+NR_END)*sizeof(int*)));
      if (!m) 
        NR_error("allocation failure 1 in NR_matrix()");
      m += NR_END;
      m--;

      m[1] = (*a)[1];
      for(i = 2; i <= nc; i++)
        m[i] = m[i-1]+nr;

      free(*a);

      *a = m;
    }

  /* nothing needs to be done if nr = nc */

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_ALLOC_CMATRIX
Create the transpose of a matrix of unsigned chars in the same memory space as 
the original matrix.  Takes as one of its arguments a vector for temporary
storage space, which must be of length at least equal to the number of rows
in the input matrix.  Allocation and deallocation may be required to maintain
the appropriate matrix data structure.
RG
*******************************************************************************/
int transpose_in_situ_alloc_cmatrix(unsigned char ***a, int nr, int nc, 
		                    unsigned char *temp_vect, char *mem_choice)
{
  int             i, j, k;   
  unsigned char  *p;          /* pointer into the input matrix */
  unsigned char  *p_temp;     /* pointer into temporary values */
  unsigned char **m;          /* pointer to new vector of pointers to rows */
  int             skip;
  int             nr_less1;
  int             nc_less1;

  /****
   * Algorithm outline:
   * For each row in the transpose, skip through the matrix, picking out
   * the appropriate element of each column, and copy them into temporary
   * storage.  Then compress the remain parts of the untransposed matrix
   * to fill in the gaps, leaving a space in memory.  In the space created 
   * copy the stored elements, thus creating a new row.
   ****/

  p = &(*a)[1][1];

  nr_less1 = nr-1;
  nc_less1 = nc-1;

  skip = nc;

  /* extract columns and compose into rows */

  for (i = 1; i <= nc_less1; i++) {
    p_temp = &temp_vect[1]; /* p_temp points into the temporary storage */

    /* move p along so the it points to every "skip"th element */
    /* (ie, the values in a column) and copy the values into the */
    /* temporary storage */

    for (j = 1; j <= nr; j++) {
      *p_temp = *p;
      p += skip;
      p_temp++;
    }

    /* move the pointer back into the array */

    p -= skip;

    /* account for effective smaller column size now that an element 
       is removed */

    skip--;

    /* pack the remaining untransposed data into the last part of memory */
    /* this leaves space in memory in front of the packed data */

    for (j = 1; j <= nr_less1; j++)
      for (k = 1; k <= skip; k++) {
        *p = *(p-j);
        p--;
      }

    /* set pointer back to the beginning of the new space */

    p -= nr_less1;
    p_temp = &temp_vect[1];

    /* copy the stored elements into the new space */

    for (j = 1; j <= nr; j++) {
      *p = *p_temp;
      p++;
      p_temp++;
    }
  }

  /* make necessary changes to new pointer scheme */

  if (!strcmp(mem_choice, "realloc") && (nr != nc)) {

    /* allocate more memory and reassign pointers */
 
    m=(unsigned char **) malloc((size_t) ((nc+NR_END)*sizeof(unsigned char*)));
    if (!m)
      NR_error("allocation failure 1 in NR_matrix()");
    m += NR_END;
    m--;

    m[1] = (*a)[1];
    for(i = 2; i <= nc; i++)
      m[i] = m[i-1]+nr;

    free(*a);

    *a = m;
  }
  else
    if (nr > nc) { 

      /* reassign pointers and free excess memory */

      for (i = 2; i <= nc; i++)
        (*a)[i] = (*a)[i-1]+nr;
      free(&(*a)[i]);
    }
    else if (nc > nr) { 

      /* allocate more memory and reassign pointers */

      m=(unsigned char **) malloc((size_t) ((nc+NR_END)*sizeof(unsigned char*)));
      if (!m) 
        NR_error("allocation failure 1 in NR_matrix()");
      m += NR_END;
      m--;

      m[1] = (*a)[1];
      for(i = 2; i <= nc; i++)
        m[i] = m[i-1]+nr;

      free(*a);

      *a = m;
    }

  /* nothing needs to be done if nr = nc */

  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_SQR_MATRIX
Create the transpose of a matrix of floats in the same memory space as the
original matrix.
RG
*******************************************************************************/
int transpose_in_situ_sqr_matrix(float **a, int n)
{
  int    i, j;
  float  *p_by_row, *p_by_col, *col_start, *row_start;
  float  temp;
 
  row_start = &a[1][1];
  col_start = &a[1][1];
 
  for (i = 1; i <= n; i++) {

    /* assign index pointers to the row and column being transposed */

    p_by_row = row_start;
    p_by_col = col_start;

    /* go through the row and column exchanging elements up to the diagonal */

    for (j = 1; j < i; j++) {

      /* switch the elements */

      temp = *p_by_row;
      *p_by_row = *p_by_col;
      *p_by_col = temp;

      p_by_row++; /* move to next element in the row */
      p_by_col += n; /* move to the next element in the column */
    }

    /* move pointers to the next row and column */

    row_start += n; /* move to the next row */
    col_start++; /* move to the next column */
  }
 
  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_SQR_DMATRIX
Create the transpose of a matrix of doubles in the same memory space as the
original matrix.
RG
*******************************************************************************/
int transpose_in_situ_sqr_dmatrix(double **a, int n)
{
  int      i, j;
  double  *p_by_row, *p_by_col, *col_start, *row_start;
  double   temp;
 
  row_start = &a[1][1];
  col_start = &a[1][1];
 
  for (i = 1; i <= n; i++) {

    /* assign index pointers to the row and column being transposed */

    p_by_row = row_start;
    p_by_col = col_start;

    /* go through the row and column exchanging elements up to the diagonal */

    for (j = 1; j < i; j++) {

      /* switch the elements */

      temp = *p_by_row;
      *p_by_row = *p_by_col;
      *p_by_col = temp;

      p_by_row++; /* move to next element in the row */
      p_by_col += n; /* move to the next element in the column */
    }

    /* move pointers to the next row and column */

    row_start += n; /* move to the next row */
    col_start++; /* move to the next column */
  }
 
  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_SQR_IMATRIX
Create the transpose of a matrix of integers in the same memory space as the
original matrix.
RG
*******************************************************************************/
int transpose_in_situ_sqr_imatrix(int **a, int n)
{
  int   i, j;
  int  *p_by_row, *p_by_col, *col_start, *row_start;
  int   temp;
 
  row_start = &a[1][1];
  col_start = &a[1][1];
 
  for (i = 1; i <= n; i++) {

    /* assign index pointers to the row and column being transposed */

    p_by_row = row_start;
    p_by_col = col_start;

    /* go through the row and column exchanging elements up to the diagonal */

    for (j = 1; j < i; j++) {

      /* switch the elements */

      temp = *p_by_row;
      *p_by_row = *p_by_col;
      *p_by_col = temp;

      p_by_row++; /* move to next element in the row */
      p_by_col += n; /* move to the next element in the column */
    }

    /* move pointers to the next row and column */

    row_start += n; /* move to the next row */
    col_start++; /* move to the next column */
  }
 
  return (UT_OK);
}


/*******************************************************************************
TRANSPOSE_IN_SITU_SQR_CMATRIX
Create the transpose of a matrix of unsigned chars in the same memory space as 
the original matrix.
RG
*******************************************************************************/
int transpose_in_situ_sqr_cmatrix(unsigned char **a, int n)
{
  int             i, j;
  unsigned char  *p_by_row, *p_by_col, *col_start, *row_start;
  unsigned char   temp;
 
  row_start = &a[1][1];
  col_start = &a[1][1];
 
  for (i = 1; i <= n; i++) {

    /* assign index pointers to the row and column being transposed */

    p_by_row = row_start;
    p_by_col = col_start;

    /* go through the row and column exchanging elements up to the diagonal */

    for (j = 1; j < i; j++) {

      /* switch the elements */

      temp = *p_by_row;
      *p_by_row = *p_by_col;
      *p_by_col = temp;

      p_by_row++; /* move to next element in the row */
      p_by_col += n; /* move to the next element in the column */
    }

    /* move pointers to the next row and column */

    row_start += n; /* move to the next row */
    col_start++; /* move to the next column */
  }
 
  return (UT_OK);
}


/*******************************************************************************
INVERT_ALLOC_MAT
Invert a matrix.  
Destroys the original.  
Allocates some memory for use in the calculation.

Ref: Numerical Recipes 2e, p. 48.
AG
*******************************************************************************/
int invert_alloc_mat(float **mat, int n, float **inv)
{
  int i, j, *index;
  float d, *col;

  /* first check for singleton matrix case */
  if (n == 1) {
    inv[1][1] = 1.0/mat[1][1];
    return (UT_OK);
  }

  /* allocate temporary storage */
  col = NR_vector(1,n);
  if (col == (float*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  index = NR_ivector(1,n);
  if (index == (int*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  
  /* invert, using LU decomposition followed by backsubstitution */
  NR_ludcmp(mat,n,index,&d);
  for(j=1; j<=n; j++) {
    for (i=1; i<=n; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    NR_lubksb(mat,n,index,col);
    for (i=1; i<=n; i++)
      inv[i][j] = col[i];
  }
  
  /* delete temporary storage */
  NR_free_vector(col, 1, n);
  NR_free_ivector(index, 1, n);
  
  return (UT_OK);
}


/*******************************************************************************
INVERT_ALLOC_DMAT
Invert a matrix.  
Destroys the original.  
Allocates some memory for use in the calculation.

Ref: Numerical Recipes 2e, p. 48.
AG
*******************************************************************************/
int invert_alloc_dmat(double **mat, int n, double **inv)
{
  int i, j, *index;
  double d, *col;

  /* first check for singleton matrix case */
  if (n == 1) {
    inv[1][1] = 1.0/mat[1][1];
    return (UT_OK);
  }

  /* allocate temporary storage */
  col = NR_dvector(1,n);
  if (col == (double*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  index = NR_ivector(1,n);
  if (index == (int*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  
  /* invert, using LU decomposition followed by backsubstitution */
  NR_dludcmp(mat,n,index,&d);
  for(j=1; j<=n; j++) {
    for (i=1; i<=n; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    NR_dlubksb(mat,n,index,col);
    for (i=1; i<=n; i++)
      inv[i][j] = col[i];
  }
  
  /* delete temporary storage */
  NR_free_dvector(col, 1, n);
  NR_free_ivector(index, 1, n);
  
  return (UT_OK);
}


/*******************************************************************************
INVERT_ALLOC_SYM_DMAT
Invert a symmetric matrix.  
Destroys the original.  
Allocates some memory for use in the calculation.

Ref: Numerical Recipes 2e, p. 97.
RG
*******************************************************************************/
int invert_alloc_sym_dmat(double **mat, int n, double **inv)
{
  int i, j;
  double d, *col, *p, *b;

  /* first check for singleton matrix case */
  if (n == 1) {
    inv[1][1] = 1.0/mat[1][1];
    return (UT_OK);
  }

  /* allocate temporary storage */
  col = NR_dvector(1,n);
  if (col == (double*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  p = NR_dvector(1,n);
  if (p == (double*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  b = NR_dvector(1,n);
  if (b == (double*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  
  /* invert, using cholesky decomposition followed by backsubstitution */
  NR_dcholdc(mat,n,p);
  for(j=1; j<=n; j++) {
    for (i=1; i<=n; i++)
      b[i] = 0.0;
    b[j] = 1.0;
    NR_dcholsl(mat,n,p,b,col);
    for (i=1; i<=n; i++)
      inv[i][j] = col[i];
  }
  
  /* delete temporary storage */
  NR_free_dvector(col, 1, n);
  NR_free_dvector(p, 1, n);
  NR_free_dvector(b, 1, n);
  
  return (UT_OK);
}


/*******************************************************************************
INVERT_ALLOC_SYM_MAT
Invert a symmetric matrix.
Destroys the original.
Allocates some memory for use in the calculation.

Ref: Numerical Recipes 2e, p. 97.
RG
*******************************************************************************/
int invert_alloc_sym_mat(float **mat, int n, float **inv)
{
  int i, j;
  float d, *col, *p, *b;

  /* first check for singleton matrix case */
  if (n == 1) {
    inv[1][1] = 1.0/mat[1][1];
    return (UT_OK);
  }

  /* allocate temporary storage */
  col = NR_vector(1,n);
  if (col == (float*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  p = NR_vector(1,n);
  if (p == (float*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  b = NR_vector(1,n);
  if (b == (float*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* invert, using cholesky decomposition followed by backsubstitution */
  NR_choldc(mat,n,p);
  for(j=1; j<=n; j++) {
    for (i=1; i<=n; i++)
      b[i] = 0.0;
    b[j] = 1.0;
    NR_cholsl(mat,n,p,b,col);
    for (i=1; i<=n; i++)
      inv[i][j] = col[i];
  }
 
  /* delete temporary storage */
  NR_free_vector(col, 1, n);
  NR_free_vector(p, 1, n);
  NR_free_vector(b, 1, n);
 
  return (UT_OK);
}

/*******************************************************************************
INVERT_COPY_ALLOC_MAT
Invert matrix copy.  
Like invert_mat(), but does not destroy the original; makes a copy first.
Allocates memory to hold a copy of the matrix.
AG
*******************************************************************************/
int invert_copy_alloc_mat(float **mat, int n, float **inv)
{
  float **mat_copy;

  /* make copy of matrix to be inverted */
  mat_copy = NR_matrix(1,n,1,n);
  copy_mat(mat, mat_copy, n, n);

  /* invert */
  invert_alloc_mat(mat_copy,n,inv);

  /* free the garbled matrix */
  NR_free_matrix(mat_copy,1,n,1,n);

  return (UT_OK);
}


/*******************************************************************************
INVERT_COPY_ALLOC_DMAT
Invert matrix copy.  
Like invert_dmat(), but does not destroy the original; makes a copy first.
Allocates memory to hold a copy of the matrix.
RG
*******************************************************************************/
int invert_copy_alloc_dmat(double **mat, int n, double **inv)
{
  double **mat_copy;

  /* make copy of matrix to be inverted */
  mat_copy = NR_dmatrix(1,n,1,n);
  copy_dmat(mat, mat_copy, n, n);

  /* invert */
  invert_alloc_dmat(mat_copy,n,inv);

  /* free the garbled matrix */
  NR_free_dmatrix(mat_copy,1,n,1,n);

  return (UT_OK);
}


/*******************************************************************************
INVERT_COPY_ALLOC_SYM_DMAT
Invert symmetric matrix copy.  
Like invert_sym_dmat(), but does not destroy the original; makes a copy first.
Allocates memory to hold a copy of the matrix.
RG
*******************************************************************************/
int invert_copy_alloc_sym_dmat(double **mat, int n, double **inv)
{
  double **mat_copy;

  /* make copy of matrix to be inverted */
  mat_copy = NR_dmatrix(1,n,1,n);
  copy_dmat(mat, mat_copy, n, n);

  /* invert */
  invert_alloc_sym_dmat(mat_copy,n,inv);

  /* free the garbled matrix */
  NR_free_dmatrix(mat_copy,1,n,1,n);

  return (UT_OK);
}


/*******************************************************************************
JACOBI_INVERT_COPY_ALLOC_DMAT
Invert matrix copy, using Jacobi preconditioning.  
Like invert_mat(), but does not destroy the original; makes a copy first.
Allocates memory to hold a copy of the matrix and store the preconditioner
values.
RG
*******************************************************************************/
int jacobi_invert_copy_alloc_dmat(double **mat, int n, double **inv)
{
  int    i, j;
  double *m;
  double **mat_copy;

  /* allocate storage for the preconditioner matrix */
  m = NR_dvector(1,n);

  /* make copy of matrix to be inverted */
  mat_copy = NR_dmatrix(1,n,1,n);
  copy_dmat(mat, mat_copy, n, n);

  /* precondition */
  for (i = 1; i <= n; i++) {
    m[i] = mat_copy[i][i];
    for (j = 1; j <= n; j++)
      mat_copy[i][j] /= m[i];
  }

  /* invert preconditioned matrix */
  invert_alloc_dmat(mat_copy,n,inv);

  /* recover the inverse */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      inv[i][j] /= m[j];

  /* free the garbled matrix */
  NR_free_dmatrix(mat_copy,1,n,1,n);

  /* free the preconditioner */
  NR_free_dvector(m, 1, n);

  return (UT_OK);
}


/*******************************************************************************
JACOBI_INVERT_COPY_ALLOC_SYM_DMAT
Invert matrix copy, using Jacobi preconditioning.  
Like invert_mat(), but does not destroy the original; makes a copy first.
Allocates memory to hold a copy of the matrix and store the preconditioner
values.
RG
*******************************************************************************/
int jacobi_invert_copy_alloc_sym_dmat(double **mat, int n, double **inv)
{
  int    i, j;
  double *m;
  double **mat_copy;

  /* allocate storage for the preconditioner matrix */
  m = NR_dvector(1,n);

  /* make copy of matrix to be inverted */
  mat_copy = NR_dmatrix(1,n,1,n);
  copy_dmat(mat, mat_copy, n, n);

  /* precondition */
  for (i = 1; i <= n; i++) {
    m[i] = mat_copy[i][i];
    for (j = 1; j <= n; j++)
      mat_copy[i][j] /= m[i];
  }

  /* invert preconditioned matrix */
  invert_alloc_sym_dmat(mat_copy,n,inv);

  /* recover the inverse */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      inv[i][j] /= m[j];

  /* free the garbled matrix */
  NR_free_dmatrix(mat_copy,1,n,1,n);

  /* free the preconditioner */
  NR_free_dvector(m, 1, n);

  return (UT_OK);
}


/*******************************************************************************
DET_ALLOC_MAT
Determinant of a matrix.
Destroys the original.  
Allocates memory for use in the computation.
Ref: Numerical Recipes 2e, p. 49.
AG
*******************************************************************************/
float det_alloc_mat(float **mat, int n)
{
  int   j;
  int  *index;
  float d;

  /* allocate temporary storage */
  index = NR_ivector(1,n);
  if (index == (int*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  
  /* use a side-effect of LU decomposition, which returns d as +1 or -1 */
  NR_ludcmp(mat,n,index,&d);

  /* the determinant of an LU decomposed matrix is the product of the diagonal
     elements */
  for (j=1; j<=n; j++)
    d *= mat[j][j];
  
  /* free temporary storage */
  NR_free_ivector(index, 1, n);

  return (d);
}


/*******************************************************************************
DET_ALLOC_DMAT
Determinant of a matrix.
Destroys the original.  
Allocates memory for use in the computation.
Ref: Numerical Recipes 2e, p. 49.
AG
*******************************************************************************/
double det_alloc_dmat(double **mat, int n)
{
  int   j;
  int  *index;
  double d;

  /* allocate temporary storage */
  index = NR_ivector(1,n);
  if (index == (int*)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
  
  /* use a side-effect of LU decomposition, which returns d as +1 or -1 */
  NR_dludcmp(mat,n,index,&d);

  /* the determinant of an LU decomposed matrix is the product of the diagonal
     elements */
  for (j=1; j<=n; j++)
    d *= mat[j][j];
  
  /* free temporary storage */
  NR_free_ivector(index, 1, n);

  return (d);
}


/*******************************************************************************
DET_COPY_ALLOC_MAT
Determinant of matrix copy.
Like det(), but does not destroy the original; makes a copy first.
Allocates memory to hold the copy.
AG
*******************************************************************************/
float det_copy_alloc_mat(float **mat, int n)
{
  float **mat_copy, d;
  
  /* make copy of matrix */
  mat_copy = NR_matrix(1,n,1,n);
  copy_mat(mat, mat_copy, n, n);

  /* find determinant */
  d = det_alloc_mat(mat_copy,n);

  /* free the garbled matrix */
  NR_free_matrix(mat_copy,1,n,1,n);

  return (d);
}


/*******************************************************************************
DET_COPY_ALLOC_DMAT
Determinant of matrix copy.
Like det(), but does not destroy the original; makes a copy first.
Allocates memory to hold the copy.
AG
*******************************************************************************/
double det_copy_alloc_dmat(double **mat, int n)
{
  double **mat_copy, d;
  
  /* make copy of matrix */
  mat_copy = NR_dmatrix(1,n,1,n);
  copy_dmat(mat, mat_copy, n, n);

  /* find determinant */
  d = det_alloc_dmat(mat_copy,n);

  /* free the garbled matrix */
  NR_free_dmatrix(mat_copy,1,n,1,n);

  return (d);
}


/*******************************************************************************
RESTRICT_ILLCOND_MATRIX
A hack procedure to restrict an ill-conditioned matrix such that it stays 
out of the realm of ill-conditioned matrices (hopefully), by ensuring that
the diagonal elements are above some small specified value.

mat is the input matrix.
dim is the size of the square matrix.
min_diag is the minimum value for a diagonal element of the matrix.
AG
*******************************************************************************/
int restrict_illcond_matrix(float **mat, int dim, float min_diag)
{
  int     i;

  /* NOTE: this should ideally compute the condition number of the matrix,
     and only do this if the condition number is less than some value */

  /* check each element of the diagonal */
  for (i = 1; i <= dim; i++) {
    if (mat[i][i] <= min_diag)
      mat[i][i] = min_diag;
    else /* add a small number to the diagonal */
      mat[i][i] += 1.0e-16;
  }

  return (UT_OK);
}


/*******************************************************************************
RESTRICT_ILLCOND_DMATRIX
A hack procedure to restrict an ill-conditioned matrix such that it stays 
out of the realm of ill-conditioned matrices (hopefully), by ensuring that
the diagonal elements are above some small specified value.

mat is the input matrix.
dim is the size of the square matrix.
min_diag is the minimum value for a diagonal element of the matrix.
AG
*******************************************************************************/
int restrict_illcond_dmatrix(double **mat, int dim, double min_diag)
{
  int     i;

  /* NOTE: this should ideally compute the condition number of the matrix,
     and only do this if the condition number is less than some value */

  /* check each element of the diagonal */
  for (i = 1; i <= dim; i++) {
    if (mat[i][i] <= min_diag)
      mat[i][i] = min_diag;
    else /* add a small number to the diagonal */
      mat[i][i] += 1.0e-16;
  }

  return (UT_OK);
}


/*******************************************************************************
EUCLID_DIST_VEC
Compute the Euclidean distance from one vector to another
RG
*******************************************************************************/
float euclid_dist_vec(float *v1, float *v2, int n)
{
  float *p1;
  float *p2;
  float *p1_end;
  float dist = 0.0;

  /* assign pointer to last element */
  p1_end = &v1[n];

  for (p1 = &v1[1], p2 = &v2[1]; p1 <= p1_end; p1++, p2++)
    dist += NR_sqr(*p1 - *p2);

  return (dist);
}


/*******************************************************************************
EUCLID_DIST_DVEC
Compute the Euclidean distance from one vector of doubles to another
RG
*******************************************************************************/
double euclid_dist_dvec(double *v1, double *v2, int n)
{
  double *p1;
  double *p2;
  double *p1_end;
  double dist = 0.0;

  /* assign pointer to last element */
  p1_end = &v1[n];

  for (p1 = &v1[1], p2 = &v2[1]; p1 <= p1_end; p1++, p2++)
    dist += NR_sqr(*p1 - *p2);

  return (dist);
}


/*******************************************************************************
RIGHT_MULT_MATRIX
Multiply a matrix by a vector on the right.
RG
*******************************************************************************/
int right_mult_matrix(float **m, int nr, int nc, float *v, float *result)
{
  float *p_m;
  float *p_v;
  float *p_result;
  float *p_end_v;
  float *p_end_result;
  
  /* initialize the result vector */
  fast_zero_vec(result, nr);
  
  /* assign pointers to last elements of u, result */
  p_end_v = &v[nc];
  p_end_result = &result[nr];
  
  /* assign pointer to first element of matrix */
  p_m = &m[1][1];
  
  /* perform the calculation */
  for (p_result = &result[1]; p_result <= p_end_result; p_result++)
    for (p_v = &v[1]; p_v <= p_end_v; p_v++) {
      *p_result += (*p_m) * (*p_v);
      p_m++;
    }

  return (UT_OK);
}


/*******************************************************************************
RIGHT_MULT_DMATRIX
Multiply a matrix of doubles by a vector of doubles on the right.
RG
*******************************************************************************/
int right_mult_dmatrix(double **m, int nr, int nc, double *v, double *result)
{
  double *p_m;
  double *p_v;
  double *p_result;
  double *p_end_v;
  double *p_end_result;
  
  /* initialize the result vector */
  fast_zero_dvec(result, nr);
  
  /* assign pointers to last elements of u, result */
  p_end_v = &v[nc];
  p_end_result = &result[nr];
  
  /* assign pointer to first element of matrix */
  p_m = &m[1][1];
  
  /* perform the calculation */
  for (p_result = &result[1]; p_result <= p_end_result; p_result++)
    for (p_v = &v[1]; p_v <= p_end_v; p_v++) {
      *p_result += (*p_m) * (*p_v);
      p_m++;
    }

  return (UT_OK);
}


/*******************************************************************************
LEFT_MULT_MATRIX
Multiply a matrix by a vector on the left.
RG
*******************************************************************************/
int left_mult_matrix(float **m, int nr, int nc, float *v, float *result)
{
  float *p_m;
  float *p_v;
  float *p_result;
  float *p_end_v;
  float *p_end_result;
  
  /* initialize the result vector */
  fast_zero_vec(result, nc);
  
  /* assign pointers to last elements of u, result */
  p_end_v = &v[nr];
  p_end_result = &result[nc];
  
  /* assign pointer to first element of matrix */
  p_m = &m[1][1];
  
  /* perform the calculation */
  for (p_v = &v[1]; p_v <= p_end_v; p_v++)
    for (p_result = &result[1]; p_result <= p_end_result; p_result++) {
      *p_result += (*p_m) * (*p_v);
      p_m++;
    }

  return (UT_OK);
}


/*******************************************************************************
LEFT_MULT_DMATRIX
Multiply a matrix of doubles by a vector of doubles on the left.
RG
*******************************************************************************/
int left_mult_dmatrix(double **m, int nr, int nc, double *v, double *result)
{
  double *p_m;
  double *p_v;
  double *p_result;
  double *p_end_v;
  double *p_end_result;
  
  /* initialize the result vector */
  fast_zero_dvec(result, nc);
  
  /* assign pointers to last elements of u, result */
  p_end_v = &v[nr];
  p_end_result = &result[nc];
  
  /* assign pointer to first element of matrix */
  p_m = &m[1][1];
  
  /* perform the calculation */
  for (p_v = &v[1]; p_v <= p_end_v; p_v++)
    for (p_result = &result[1]; p_result <= p_end_result; p_result++) {
      *p_result += (*p_m) * (*p_v);
      p_m++;
    }

  return (UT_OK);
}


/*******************************************************************************
MAT_MULT
Multiply two matrices.
AG
*******************************************************************************/
int mat_mult(float **a, int nra, int nca, float **b, int nrb, int ncb, 
             float **c)
{
  int    i, j, k;
  int    nrc, ncc;
  float *p_a;
  float *p_b;
  float *p_c;

  /* check for compatible dimensions */
  if (nca != nrb) {
    err_printf();
    log_printf("Matrix dimensions are incompatible\n");
    return (UT_ERROR);
  }

  nrc = nra;
  ncc = ncb;
 
  fast_zero_mat(c, nrc, ncc);

  p_c = &c[1][1];

  /* accumulate results */
  for (i = 1; i <= nrc; i++) {
    for (j = 1; j <= ncc; j++) {
      p_a = &a[i][1];
      p_b = &b[1][j];
      for (k = 1; k <= nca; k++) {
        *p_c += (*p_a) * (*p_b);
        p_a++;
        p_b += ncb;
      }
      p_c++;
    }
  }

  return (UT_OK);
}


/*******************************************************************************
DMAT_MULT
Multiply two matrices of doubles.
AG
*******************************************************************************/
int dmat_mult(double **a, int nra, int nca, double **b, int nrb, int ncb, 
              double **c)
{
  int    i, j, k;
  int    nrc, ncc;
  double *p_a;
  double *p_b;
  double *p_c;

  /* check for compatible dimensions */
  if (nca != nrb) {
    err_printf();
    log_printf("Matrix dimensions are incompatible\n");
    return (UT_ERROR);
  }

  nrc = nra;
  ncc = ncb;
 
  fast_zero_dmat(c, nrc, ncc);

  p_c = &c[1][1];

  /* accumulate results */
  for (i = 1; i <= nrc; i++) {
    for (j = 1; j <= ncc; j++) {
      p_a = &a[i][1];
      p_b = &b[1][j];
      for (k = 1; k <= nca; k++) {
        *p_c += (*p_a) * (*p_b);
        p_a++;
        p_b += ncb;
      }
      p_c++;
    }
  }

  return (UT_OK);
}


/*******************************************************************************
FAST_MAT_MULT
A more efficient routine for performing matrix multiplies.  Based off of code
originally written by Edward Keyes at M.I.T.  The general idea is to minimize
memory access by using the register variables as much as possible.  The
routine calculates a 4x4 block of answers each pass.  Note that the boundary
cases are not treated differently; answers are calculated with dummy values,
but just never stored (this leaves room for improvement).
RG
*******************************************************************************/
int fast_mat_mult(float **A, int nra, int nca, float **B, int nrb, int ncb, 
                  float **C)
{
  register float c00, c01, c02, c03;
  register float c10, c11, c12, c13;
  register float c20, c21, c22, c23;
  register float c30, c31, c32, c33;
  register float ax0, ax1, ax2, ax3;
  register float bx0, bx1, bx2, bx3;

  register float *baseC, *baseA, *baseB;
  register int   numr, numc, numi;
  register int   i, j, k;
  register int   border1, border2;
  register int   s1, s2, s3;
  register int   mod1, mod2, mod3;

  /* initialize the output array */
  fast_zero_mat(C, nra, ncb);

  /* round to next highest multiple of 4 */
  numr = (nra + 3) / 4;   
  numc = (ncb + 3) / 4;
  mod1 = nra % 4;
  mod2 = nca % 2;
  mod3 = ncb % 4;
  if (mod1 == 0) 
    mod1 = 4;
  if (mod3 == 0) 
    mod3 = 4;
  numi = nca / 2;
  s1 = nra; 
  s2 = nca; 
  s3 = ncb;
  i = 0;
  while (i < numr) { /*  loop over rows of C  */
    if (i == numr - 1)
      border1 = 1;
    else 
      border1 = 0;
    border2 = 0;
    j = 0;
    baseC = &C[1][1] + i * s3 * 4;
    while (j < numc) { /*  loop over columns of C  */
      if (j == numc - 1)
        border2 = 1;
      k = 0;
      baseA = &A[1][1] + i  *  4  *  s2;
      baseB = &B[1][1] + j  *  4;
       
      /*  increment early for the pipeline  */
      j++; 

      /*  copy registers */
      c00 = baseC[0];
      c01 = baseC[1];
      c02 = baseC[2];
      c03 = baseC[3];
      baseC += s3;
      c10 = baseC[0];
      c11 = baseC[1];
      c12 = baseC[2];
      c13 = baseC[3];
      baseC += s3;
      c20 = baseC[0];
      c21 = baseC[1];
      c22 = baseC[2];
      c23 = baseC[3];
      baseC += s3;
      c30 = baseC[0];
      c31 = baseC[1];
      c32 = baseC[2];
      c33 = baseC[3];
      while (k < numi) { /*  do the calculations  */
        k++;
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2 + s2];
        ax3 = baseA[s2 + s2 + s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2+s2];
        ax3 = baseA[s2+s2+s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
      }
      if (mod2 == 1) {
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2 + s2];
        ax3 = baseA[s2 + s2+ s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
      }
      if ((!(border1 + border2)) || (mod1 + mod3 == 8)) {
        baseC[0] = c30;
        baseC[1] = c31;
        baseC[2] = c32;
        baseC[3] = c33;
        baseC -= s3;
        baseC[0] = c20;
        baseC[1] = c21;
        baseC[2] = c22;
        baseC[3] = c23;
        baseC -= s3;
        baseC[0] = c10;
        baseC[1] = c11;
        baseC[2] = c12;
        baseC[3] = c13;
        baseC -= s3;
        baseC[0] = c00;
        baseC[1] = c01;
        baseC[2] = c02;
        baseC[3] = c03;
      }
      else if (border1 - border2 == 1) {
        if (mod1 > 3) {
          baseC[0] = c30;
          baseC[1] = c31;
          baseC[2] = c32;
          baseC[3] = c33;
        }
        baseC -= s3;
        if (mod1 > 2) {
          baseC[0] = c20;
          baseC[1] = c21;
          baseC[2] = c22;
          baseC[3] = c23;
        }
        baseC -= s3;
        if (mod1 > 1) {
          baseC[0] = c10;
          baseC[1] = c11;
          baseC[2] = c12;
          baseC[3] = c13;
        }
        baseC -= s3;
        baseC[0] = c00;
        baseC[1] = c01;
        baseC[2] = c02;
        baseC[3] = c03;
      }
      else if (border2 - border1 == 1) {
        baseC[0] = c30;
        if (mod3 > 1) 
	  baseC[1] = c31;
        if (mod3 > 2) 
	  baseC[2] = c32;
        if (mod3 > 3) 
	  baseC[3] = c33;
        baseC -= s3;
        baseC[0] = c20;
        if (mod3 > 1) 
	  baseC[1] = c21;
        if (mod3 > 2) 
	  baseC[2] = c22;
        if (mod3 > 3) 
	  baseC[3] = c23;
        baseC -= s3;
        baseC[0] = c10;
        if (mod3 > 1) 
	  baseC[1] = c11;
        if (mod3 > 2) 
	  baseC[2] = c12;
        if (mod3 > 3) 
	  baseC[3] = c13;
        baseC -= s3;
        baseC[0] = c00;
        if (mod3 > 1) 
	  baseC[1] = c01;
        if (mod3 > 2) 
          baseC[2] = c02;
        if (mod3 > 3) 
          baseC[3] = c03;
      }
      else if (border1 + border2 == 2) {
        if (mod1 > 3) {
          baseC[0] = c30;
          if (mod3 > 1) 
	    baseC[1] = c31;
          if (mod3 > 2) 
	    baseC[2] = c32;
          if (mod3 > 3) 
	    baseC[3] = c33;
        }
        baseC -= s3;
        if (mod1 > 2) {
          baseC[0] = c20;
          if (mod3 > 1) 
	    baseC[1] = c21;
          if (mod3 > 2) 
	    baseC[2] = c22;
          if (mod3 > 3) 
	    baseC[3] = c23;
        }
        baseC -= s3;
        if (mod1 > 1) {
          baseC[0] = c10;
          if (mod3 > 1) 
	    baseC[1] = c11;
          if (mod3 > 2) 
	    baseC[2] = c12;
          if (mod3 > 3) 
	    baseC[3] = c13;
        }
        baseC -= s3;
        baseC[0] = c00;
        if (mod3 > 1) 
	  baseC[1] = c01;
        if (mod3 > 2) 
	  baseC[2] = c02;
        if (mod3 > 3) 
	  baseC[3] = c03;
      }
      baseC += 4; /*  next group of results  */
    }
    i++;
  }
  
  return (UT_OK);
}


/*******************************************************************************
FAST_DMAT_MULT
A more efficient routine for performing matrix multiplies.  Based off of code
originally written by Edward Keyes at M.I.T.  The general idea is to minimize
memory access by using the register variables as much as possible.  The
routine calculates a 4x4 block of answers each pass.  Note that the boundary
cases are not treated differently; answers are calculated with dummy values,
but just never stored (this leaves room for improvement).
RG
*******************************************************************************/
int fast_dmat_mult(double **A, int nra, int nca, double **B, int nrb, int ncb, 
                   double **C)
{
  register double c00, c01, c02, c03;
  register double c10, c11, c12, c13;
  register double c20, c21, c22, c23;
  register double c30, c31, c32, c33;
  register double ax0, ax1, ax2, ax3;
  register double bx0, bx1, bx2, bx3;

  register double *baseC, *baseA, *baseB;
  register int   numr, numc, numi;
  register int   i, j, k;
  register int   border1, border2;
  register int   s1, s2, s3;
  register int   mod1, mod2, mod3;

  /* initialize the output array */
  fast_zero_dmat(C, nra, ncb);

  /* round to next highest multiple of 4 */
  numr = (nra + 3) / 4;   
  numc = (ncb + 3) / 4;
  mod1 = nra % 4;
  mod2 = nca % 2;
  mod3 = ncb % 4;
  if (mod1 == 0) 
    mod1 = 4;
  if (mod3 == 0) 
    mod3 = 4;
  numi = nca / 2;
  s1 = nra; 
  s2 = nca; 
  s3 = ncb;
  i = 0;
  while (i < numr) { /*  loop over rows of C  */
    if (i == numr - 1)
      border1 = 1;
    else 
      border1 = 0;
    border2 = 0;
    j = 0;
    baseC = &C[1][1] + i * s3 * 4;
    while (j < numc) { /*  loop over columns of C  */
      if (j == numc - 1)
        border2 = 1;
      k = 0;
      baseA = &A[1][1] + i  *  4  *  s2;
      baseB = &B[1][1] + j  *  4;
       
      /*  increment early for the pipeline  */
      j++; 

      /*  copy registers */
      c00 = baseC[0];
      c01 = baseC[1];
      c02 = baseC[2];
      c03 = baseC[3];
      baseC += s3;
      c10 = baseC[0];
      c11 = baseC[1];
      c12 = baseC[2];
      c13 = baseC[3];
      baseC += s3;
      c20 = baseC[0];
      c21 = baseC[1];
      c22 = baseC[2];
      c23 = baseC[3];
      baseC += s3;
      c30 = baseC[0];
      c31 = baseC[1];
      c32 = baseC[2];
      c33 = baseC[3];
      while (k < numi) { /*  do the calculations  */
        k++;
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2 + s2];
        ax3 = baseA[s2 + s2 + s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2+s2];
        ax3 = baseA[s2+s2+s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
      }
      if (mod2 == 1) {
        ax0 = baseA[0];
        ax1 = baseA[s2];
        ax2 = baseA[s2 + s2];
        ax3 = baseA[s2 + s2+ s2];
        bx0 = baseB[0];
        bx1 = baseB[1];
        bx2 = baseB[2];
        bx3 = baseB[3];
        baseA += 1;
        baseB += s3;
        c00 += ax0 * bx0;
        c01 += ax0 * bx1;
        c02 += ax0 * bx2;
        c03 += ax0 * bx3;
        c10 += ax1 * bx0;
        c11 += ax1 * bx1;
        c12 += ax1 * bx2;
        c13 += ax1 * bx3;
        c20 += ax2 * bx0;
        c21 += ax2 * bx1;
        c22 += ax2 * bx2;
        c23 += ax2 * bx3;
        c30 += ax3 * bx0;
        c31 += ax3 * bx1;
        c32 += ax3 * bx2;
        c33 += ax3 * bx3;
      }
      if ((!(border1 + border2)) || (mod1 + mod3 == 8)) {
        baseC[0] = c30;
        baseC[1] = c31;
        baseC[2] = c32;
        baseC[3] = c33;
        baseC -= s3;
        baseC[0] = c20;
        baseC[1] = c21;
        baseC[2] = c22;
        baseC[3] = c23;
        baseC -= s3;
        baseC[0] = c10;
        baseC[1] = c11;
        baseC[2] = c12;
        baseC[3] = c13;
        baseC -= s3;
        baseC[0] = c00;
        baseC[1] = c01;
        baseC[2] = c02;
        baseC[3] = c03;
      }
      else if (border1 - border2 == 1) {
        if (mod1 > 3) {
          baseC[0] = c30;
          baseC[1] = c31;
          baseC[2] = c32;
          baseC[3] = c33;
        }
        baseC -= s3;
        if (mod1 > 2) {
          baseC[0] = c20;
          baseC[1] = c21;
          baseC[2] = c22;
          baseC[3] = c23;
        }
        baseC -= s3;
        if (mod1 > 1) {
          baseC[0] = c10;
          baseC[1] = c11;
          baseC[2] = c12;
          baseC[3] = c13;
        }
        baseC -= s3;
        baseC[0] = c00;
        baseC[1] = c01;
        baseC[2] = c02;
        baseC[3] = c03;
      }
      else if (border2 - border1 == 1) {
        baseC[0] = c30;
        if (mod3 > 1) 
	  baseC[1] = c31;
        if (mod3 > 2) 
	  baseC[2] = c32;
        if (mod3 > 3) 
	  baseC[3] = c33;
        baseC -= s3;
        baseC[0] = c20;
        if (mod3 > 1) 
	  baseC[1] = c21;
        if (mod3 > 2) 
	  baseC[2] = c22;
        if (mod3 > 3) 
	  baseC[3] = c23;
        baseC -= s3;
        baseC[0] = c10;
        if (mod3 > 1) 
	  baseC[1] = c11;
        if (mod3 > 2) 
	  baseC[2] = c12;
        if (mod3 > 3) 
	  baseC[3] = c13;
        baseC -= s3;
        baseC[0] = c00;
        if (mod3 > 1) 
	  baseC[1] = c01;
        if (mod3 > 2) 
          baseC[2] = c02;
        if (mod3 > 3) 
          baseC[3] = c03;
      }
      else if (border1 + border2 == 2) {
        if (mod1 > 3) {
          baseC[0] = c30;
          if (mod3 > 1) 
	    baseC[1] = c31;
          if (mod3 > 2) 
	    baseC[2] = c32;
          if (mod3 > 3) 
	    baseC[3] = c33;
        }
        baseC -= s3;
        if (mod1 > 2) {
          baseC[0] = c20;
          if (mod3 > 1) 
	    baseC[1] = c21;
          if (mod3 > 2) 
	    baseC[2] = c22;
          if (mod3 > 3) 
	    baseC[3] = c23;
        }
        baseC -= s3;
        if (mod1 > 1) {
          baseC[0] = c10;
          if (mod3 > 1) 
	    baseC[1] = c11;
          if (mod3 > 2) 
	    baseC[2] = c12;
          if (mod3 > 3) 
	    baseC[3] = c13;
        }
        baseC -= s3;
        baseC[0] = c00;
        if (mod3 > 1) 
	  baseC[1] = c01;
        if (mod3 > 2) 
	  baseC[2] = c02;
        if (mod3 > 3) 
	  baseC[3] = c03;
      }
      baseC += 4; /*  next group of results  */
    }
    i++;
  }
  
  return (UT_OK);
}


/*******************************************************************************
ADD_MAT
Add one matrix element-wise to another.
RG
*******************************************************************************/
int add_mat(float **m1, float **m2, float **m3, int nr, int nc)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;
 
  p1 = &m1[1][1];
  p2 = &m2[1][1];

  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = *p1 + *p2;
    p1++;
    p2++;
  }

  return (UT_OK);
}

/*******************************************************************************
ADD_DMAT
Add one matrix of doubles element-wise to another.
RG
*******************************************************************************/
int add_dmat(double **m1, double **m2, double **m3, int nr, int nc)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;
 
  p1 = &m1[1][1];
  p2 = &m2[1][1];

  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = *p1 + *p2;
    p1++;
    p2++;
  }

  return (UT_OK);
}

/*******************************************************************************
SUBTRACT_MAT
Subtract one matrix element-wise from another.  The second matrix is
subtracted from the first matrix.
RG
*******************************************************************************/
int subtract_mat(float **m1, float **m2, float **m3, int nr, int nc)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;
 
  p1 = &m1[1][1];
  p2 = &m2[1][1];

  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = *p1 - *p2;
    p1++;
    p2++;
  }

  return (UT_OK);
}


/*******************************************************************************
SUBTRACT_DMAT
Subtract one matrix of doubles element-wise from another.  The second matrix is
subtracted from the first matrix.
RG
*******************************************************************************/
int subtract_dmat(double **m1, double **m2, double **m3, int nr, int nc)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;
 
  p1 = &m1[1][1];
  p2 = &m2[1][1];

  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = *p1 - *p2;
    p1++;
    p2++;
  }

  return (UT_OK);
}

/*******************************************************************************
ADD_VEC
Add one vector element-wise to another.
RG
*******************************************************************************/
int add_vec(float *v1, float *v2, float *v3, int n)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;
 
  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = *p1 + *p2;

  return (UT_OK);
}

/*******************************************************************************
ADD_DVEC
Add one vector of doubles element-wise to another.
RG
*******************************************************************************/
int add_dvec(double *v1, double *v2, double *v3, int n)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;
 
  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = *p1 + *p2;

  return (UT_OK);
}

/*******************************************************************************
SUBTRACT_VEC
Subtract one vector element-wise from another.  The second vector is 
subtracted from the first vector.
RG
*******************************************************************************/
int subtract_vec(float *v1, float *v2, float *v3, int n)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;
 
  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = *p1 - *p2;  

  return (UT_OK);
}

/*******************************************************************************
SUBTRACT_DVEC
Subtract one vector of doubles element-wise from another.  The second
vector is subtracted from the first vector.
RG
*******************************************************************************/
int subtract_dvec(double *v1, double *v2, double *v3, int n)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;
 
  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = *p1 - *p2;

  return (UT_OK);
}

/*******************************************************************************
DOT_PRODUCT_VEC
Compute the dot product of two vectors.
RG
*******************************************************************************/
float dot_product_vec(float *x, float *y, int n)
{
  float *p_x;
  float *p_y;
  float *p_end_x;
  float  result = 0.0;

  /* assign pointer to last element of x */
  p_end_x = &x[n];

  /* compute dot product */
  for (p_x = &x[1], p_y = &y[1]; p_x <= p_end_x; p_x++, p_y++)
    result += (*p_x) * (*p_y);

  return (result);
}

/*******************************************************************************
DOT_PRODUCT_DVEC
Compute the dot product of two vectors of doubles.
RG
*******************************************************************************/
double dot_product_dvec(double *x, double *y, int n)
{
  double *p_x;
  double *p_y;
  double *p_end_x;
  double  result = 0.0;

  /* assign pointer to last element of x */
  p_end_x = &x[n];

  /* compute dot product */
  for (p_x = &x[1], p_y = &y[1]; p_x <= p_end_x; p_x++, p_y++)
    result += (*p_x) * (*p_y);

  return (result);
}

/*******************************************************************************
OUTER_PRODUCT_VEC
Compute the outer product of two vectors.
RG
*******************************************************************************/
int outer_product_vec(float *x, float *y, int n, float **prod)
{
  float *p_x;
  float *p_y;
  float *p_end_x;
  float *p_end_y;
  float *p_prod;

  /* assign pointers to last elements of x and y */
  p_end_x = &x[n];
  p_end_y = &y[n];

  /* assign pointer to first element of outer product matrix */
  p_prod = &prod[1][1];

  /* calculate the outer product */
  for (p_x = &x[1]; p_x <= p_end_x; p_x++)
    for (p_y = &y[1]; p_y <= p_end_y; p_y++) {
      *p_prod = (*p_x) * (*p_y);
      p_prod++;
    }

  return (UT_OK);
}  

/*******************************************************************************
OUTER_PRODUCT_DVEC
Compute the outer product of two vectors of doubles.
RG
*******************************************************************************/
int outer_product_dvec(double *x, double *y, int n, double **prod)
{
  double *p_x;
  double *p_y;
  double *p_end_x;
  double *p_end_y;
  double *p_prod;

  /* assign pointers to last elements of x and y */
  p_end_x = &x[n];
  p_end_y = &y[n];

  /* assign pointer to first element of outer product matrix */
  p_prod = &prod[1][1];

  /* calculate the outer product */
  for (p_x = &x[1]; p_x <= p_end_x; p_x++)
    for (p_y = &y[1]; p_y <= p_end_y; p_y++) {
      *p_prod = (*p_x) * (*p_y);
      p_prod++;
    }

  return (UT_OK);
}

/*******************************************************************************
SCALAR_MULT_VEC
Scales a vector by a specified multiplicand.

n is the length of the vector.
v is the input vector.
constant is the multiplicand to scale the vector by.
RG,AG
*******************************************************************************/
int scalar_mult_vec(float *v, int n,  float constant)
{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p *= constant;

  return (UT_OK);
}

/*******************************************************************************
SCALAR_MULT_DVEC
Scales a vector of doubles by a specified multiplicand.

n is the length of the vector.
v is the input vector.
constant is the multiplicand to scale the vector by.
RG
*******************************************************************************/
int scalar_mult_dvec(double *v, int n, double constant)
{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p *= constant;

  return (UT_OK);
}

/*******************************************************************************
SCALAR_DIV_VEC
Scales a vector by a specified dividend.

n is the length of the vector.
v is the input vector.
constant is the dividend to scale the vector by.
RG
*******************************************************************************/
int scalar_div_vec(float *v, int n, float constant)
{
  float  *p;
  float  *p_end;
  float  inv_constant;

  /* assign pointer to last element */
  p_end = &v[n];

  /* calculate the inverse of the constant to save on divides */
  if (constant == 0.0)
    inv_constant = FLT_MAX;
  else
    inv_constant = 1.0 / constant;

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p *= inv_constant;

  return (UT_OK);
}

/*******************************************************************************
SCALAR_DIV_DVEC
Scales a vector of doubles by a specified dividend.

n is the length of the vector.
v is the input vector.
constant is the dividend to scale the vector by.
RG
*******************************************************************************/
int scalar_div_dvec(double *v, int n, double constant)
{
  double  *p;
  double  *p_end;
  double  inv_constant;

  /* assign pointer to last element */
  p_end = &v[n];

  /* calculate the inverse of the constant to save on divides */
  if (constant == 0.0)
    inv_constant = DBL_MAX;
  else
    inv_constant = 1.0 / constant;

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p *= inv_constant;

  return (UT_OK);
}

/*******************************************************************************
SCALAR_ADD_VEC
Translates a vector by a specified addend.

nc is the length of the vector.
v is the input vector.
constant is the addend to translate the vector by.
RG
*******************************************************************************/
int scalar_add_vec(float *v, int n, float constant)
{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p += constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_ADD_DVEC
Translates a vector of doubles by a specified addend.

nc is the length of the vector.
v is the input vector.
constant is the addend to translate the vector by.
RG
*******************************************************************************/
int scalar_add_dvec(double *v, int n, double constant)
{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p += constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_SUBTRACT_VEC
Translates a vector by a specified subtrahend.

nc is the length of the vector.
v is the input vector.
constant is the subtrahend to translate the vector by.
RG
*******************************************************************************/
int scalar_subtract_vec(float *v, int n, float constant)
{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p -= constant;

  return (UT_OK);
}

/*******************************************************************************
SCALAR_SUBTRACT_DVEC
Translates a vector of doubles by a specified subtrahend.

nc is the length of the vector.
v is the input vector.
constant is the subtrahend to translate the vector by.
RG
*******************************************************************************/
int scalar_subtract_dvec(double *v, int n, double constant)
{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &v[n];

  /* change each element of the vector */
  for (p = &v[1]; p <= p_end; p++)
    *p -= constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_MULT_MAT
Scales a matrix by a specified multiplicand.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the multiplicand to scale the matrix by.
RG
*******************************************************************************/
int scalar_mult_mat(float **m, int nr, int nc, float constant)
{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p *= constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_MULT_DMAT
Scales a matrix of doubles by a specified multiplicand.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the multiplicand to scale the matrix by.
RG
*******************************************************************************/
int scalar_mult_dmat(double **m, int nr, int nc, double constant)
{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p *= constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_DIV_MAT
Scales a matrix by a specified dividend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the dividend to scale the matrix by.
RG
*******************************************************************************/
int scalar_div_mat(float **m, int nr, int nc, float constant)
{
  float  *p;
  float  *p_end;
  float  inv_constant;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* calculate the inverse of the constant to save on divides */
  if (constant == 0.0)
    inv_constant = FLT_MAX;
  else
    inv_constant = 1.0 / constant;

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p *= inv_constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_DIV_DMAT
Scales a matrix of doubles by a specified dividend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the dividend to scale the matrix by.
RG
*******************************************************************************/
int scalar_div_dmat(double **m, int nr, int nc, double constant)
{
  double  *p;
  double  *p_end;
  double  inv_constant;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* calculate the inverse of the constant to save on divides */
  if (constant == 0.0)
    inv_constant = DBL_MAX;
  else
    inv_constant = 1.0 / constant;

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p *= inv_constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_ADD_MAT
Translates a matrix by a specified addend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the addend to translate the matrix by.
RG
*******************************************************************************/
int scalar_add_mat(float **m, int nr, int nc, float constant)

{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p += constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_ADD_DMAT
Translates a matrix of doubles by a specified addend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the addend to translate the matrix by.
RG
*******************************************************************************/
int scalar_add_dmat(double **m, int nr, int nc, double constant)

{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p += constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_SUBTRACT_MAT
Translates a matrix by a specified subtrahend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the subtrahend to translate the matrix by.
RG
*******************************************************************************/
int scalar_subtract_mat(float **m, int nr, int nc, float constant)
{
  float  *p;
  float  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p -= constant;

  return (UT_OK);
}


/*******************************************************************************
SCALAR_SUBTRACT_DMAT
Translates a matrix of doubles by a specified subtrahend.

nr and nc are the dimensions of the matrix (num. rows and num. cols).
m is the input matrix.
constant is the subtrahend to translate the matrix by.
RG
*******************************************************************************/
int scalar_subtract_dmat(double **m, int nr, int nc, double constant)
{
  double  *p;
  double  *p_end;

  /* assign pointer to last element */
  p_end = &m[nr][nc];

  /* change each element of the matrix */
  for (p = &m[1][1]; p <= p_end; p++)
    *p -= constant;

  return (UT_OK);
}


/*******************************************************************************
MULT_VEC_ELT
Multiply one vector element-wise by another.
RG
*******************************************************************************/
int mult_vec_elt(float *v1, float *v2, float *v3, int n)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;

  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = (*p1) * (*p2);

  return (UT_OK);
}


/*******************************************************************************
MULT_DVEC_ELT
Multiply one vector of doubles element-wise by another.
RG
*******************************************************************************/
int mult_dvec_elt(double *v1, double *v2, double *v3, int n)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;

  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++);
    *p3 = (*p1) * (*p2);

  return (UT_OK);
}


/*******************************************************************************
DIV_VEC_ELT
Divide one vector element-wise by another.
RG
*******************************************************************************/
int div_vec_elt(float *v1, float *v2, float *v3, int n)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;

  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = (*p1) / (*p2);

  return (UT_OK);
}


/*******************************************************************************
DIV_DVEC_ELT
Divide one vector element-wise by another.
RG
*******************************************************************************/
int div_dvec_elt(double *v1, double *v2, double *v3, int n)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;

  /* assign pointer to last element of output vector */
  p3_end = &v3[n];

  for (p1 = &v1[1], p2 = &v2[1], p3 = &v3[1]; p3 <= p3_end; p1++, p2++, p3++)
    *p3 = (*p1) / (*p2);

  return (UT_OK);
}


/*******************************************************************************
MULT_MAT_ELT
Multiply one matrix element-wise by another.
RG
*******************************************************************************/
int mult_mat_elt(float **m1, float **m2, float **m3, int nr, int nc)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;

  p1 = &m1[1][1];
  p2 = &m2[1][1];
  
  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = (*p1) * (*p2);
    p1++;
    p2++;
  }

  return (UT_OK);
}


/*******************************************************************************
MULT_DMAT_ELT
Multiply one matrix of doubles element-wise by another.
RG
*******************************************************************************/
int mult_dmat_elt(double **m1, double **m2, double **m3, int nr, int nc)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;

  p1 = &m1[1][1];
  p2 = &m2[1][1];
  
  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = (*p1) * (*p2);
    p1++;
    p2++;
  }

  return (UT_OK);
}


/*******************************************************************************
DIV_MAT_ELT
Divide one matrix element-wise by another.
RG
*******************************************************************************/
int div_mat_elt(float **m1, float **m2, float **m3, int nr, int nc)
{
  float  *p1;
  float  *p2;
  float  *p3;
  float  *p3_end;

  p1 = &m1[1][1];
  p2 = &m2[1][1];
  
  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = (*p1) / (*p2);
    p1++;
    p2++;
  }

  return (UT_OK);
}


/*******************************************************************************
DIV_DMAT_ELT
Divide one matrix of doubles element-wise by another.
RG
*******************************************************************************/
int div_dmat_elt(double **m1, double **m2, double **m3, int nr, int nc)
{
  double  *p1;
  double  *p2;
  double  *p3;
  double  *p3_end;

  p1 = &m1[1][1];
  p2 = &m2[1][1];
  
  /* assign pointer to last element of output vector */
  p3_end = &m3[nr][nc];

  for (p3 = &m3[1][1]; p3 <= p3_end; p3++) {
    *p3 = (*p1) / (*p2);
    p1++;
    p2++;
  }

  return (UT_OK);
}


/*******************************************************************************
PCA_MAT
Perform principal components analysis (PCA) on a matrix.  The output is three 
matrices, one containing the principal components for the data, one containing
the associated eigenvalues, and one containing the other orthogonal matrix in 
the decomposition.  These correspond to U, W, and V^T respectively in the 
decomposition of the input matrix A = U * W * V^T.

Done using singular value decomposition (SVD).

Notes: A is overwritten by U.  nr and nc are the dimensions of A (and thus U).
Use the vector w instead of the matrix W to hold the eigenvalues, of size nc.
VT has dimensions nc x nc.
AG, RG
*******************************************************************************/
int pca_mat(float **A, int nr, int nc, float *w, float **VT)
{
  int maxiter;

  /* heuristic to set the number of svd iterations */
  maxiter = NR_imax(30, (int) (nr * nc / 10000));

  /* compute the svd; A will be overwritten by U */
  DA_svdcmp(A, nr, nc, w, VT, maxiter);

  /* this returns V rather than the transpose V^T */
  transpose_in_situ_sqr_matrix(VT, nc);

  return (UT_OK);
}


/*******************************************************************************
PCA_ALLOC_MAT
Perform principal components analysis (PCA) on a matrix.  The output is three 
matrices, one containing the principal components for the data, one containing
the associated eigenvalues, and one containing the other orthogonal matrix in 
the decomposition.  These correspond to U, W, and V^T respectively in the 
decomposition of the input matrix A = U * W * V^T.

Done using singular value decomposition (SVD).

Notes: A is overwritten by U.  nr and nc are the dimensions of A (and thus U).
Use the vector w instead of the matrix W to hold the eigenvalues, of size nc.
VT has dimensions nc x nc.
AG
*******************************************************************************/
int pca_alloc_mat(float **A, int nr, int nc, float *w, float **VT)
{
  float **V;
  int maxiter;

  /* allocate the transpose of VT, V */
  V = NR_matrix(1, nc, 1, nc);
  if (V == (float**)NULL) {
    err_printf();
    log_printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* heuristic to set the number of svd iterations */
  maxiter = NR_imax(30, (int) (nr * nc / 10000));

  /* compute the svd; A will be overwritten by U */
  DA_svdcmp(A, nr, nc, w, V, maxiter);

  /* this returns V rather than the transpose V^T */
  transpose_matrix(V, nc, nc, VT);

  NR_free_matrix(V, 1, nc, 1, nc);
  return (UT_OK);
}


/*******************************************************************************
LSQR_MAT
An iterative method for solving least-squares/minimum norm problems of the
form Ax = b.  It requires three workspace vectors u, v and w of length nr, nc,
and nc respectively, where nr and nc are the dimensions of A.  The stopping
conditions are set by three input parameters, atol, btol, and conlim.  
Recommended values for these parameters depend on the machine value EPSILON:

atol = EPSILON
btol = EPSILON
conlim = (1.0 / (10 * sqrt(EPSILON).

Using floating point numbers may result in poorer convergence properties;
use the double-precision version of this function for better convergence 
properties.
 
Reference: CC Paige, MA Saunders.  "LSQR: An Algorithm for Sparse Linear
Equations and Sparse Least Squares."  ACM Trans Math Soft, Vol 8, No 1,
March 1982, pp43-71.
 
Most of the notation used in this function is derived from this paper.
RG
*******************************************************************************/
int lsqr_mat(float **A, int nr, int nc, float *b, float *x, float *u,
             float *v, float *w, float atol, float btol, float conlim)
{
  int   i, j;
  float alpha;
  float beta;
  float norm_A;
  float norm_b;
  float norm_B;
  float norm_D;
  float norm_Ar;
  float norm_r;
  float norm_x;
  float norm_xx;
  float norm_d;
  float cond_A;
  float phi, phi_bar;
  float rho, rho_bar;
  float gamma, gamma_bar;
  float z, z_bar;
  float delta;
  float c;
  float s;
  float theta;
  float sn;
  float cs;
  float rhs;
  int   iters;
  int   stop1, stop2, stop3;
 
  /* initialize */

  copy_vec(b, u, nr);
 
  beta = norm_vec(u, nr);
  scalar_div_vec(u, nr, beta);

  left_mult_matrix(A, nr, nc, u, v);
  alpha = norm_vec(v, nc);
  scalar_div_vec(v, nc, alpha);

  copy_vec(v, w, nc);;
  
  fast_zero_vec(x, nc);
  
  phi_bar = beta;
  rho_bar = alpha;
 
  /* quantities needed for stopping conditions */
 
  norm_A = frobenius_norm_mat(A, nr, nc);
  norm_b = norm_vec(b, nr);

  cond_A = 0.0;
  norm_B = 0.0;
  norm_D = 0.0;
  norm_x = 0.0;
  norm_Ar = 0.0;

  /* quantities used for estimating the norm of x */

  cs = -1.0;
  sn = 0.0;
  z = 0.0;
  norm_xx = 0.0;

  /* now perform the iterations */
 
  iters = 0;
  stop1 = 0;
  stop2 = 0;
  stop3 = 0;
 
  do {
    iters++;
 
    /* bidiagonalization */
 
    for (i = 1; i <= nr; i++) {
      u[i] *= -alpha;
      for (j = 1; j <= nc; j++)
        u[i] += A[i][j] * v[j];
    }
 
    beta = norm_vec(u, nr);
    scalar_div_vec(u, nr, beta);
 
    /* before we overwrite alpha(k) with alpha(k+1) */
 
    norm_B = norm_B + NR_sqr(alpha) + NR_sqr(beta);
 
    /* continue with bidiagonalization */
 
    for (j = 1; j <= nc; j++) {
      v[j] *= -beta;
      for (i = 1; i <= nr; i++)
        v[j] += A[i][j] * u[i];
    }
 
    alpha = norm_vec(v, nc);
    scalar_div_vec(v, nc, alpha);
 
    /* orthogonal transformation */
 
    rho = sqrt(NR_sqr(rho_bar) + NR_sqr(beta));
    c = rho_bar / rho;
    s = beta / rho;
    theta = s * alpha;
    rho_bar = -c * alpha;
    phi = c * phi_bar;
    phi_bar *= s;
 
    /* update x, w */
 
    for (i = 1; i <= nc; i++) {
      x[i] += phi * w[i] / rho;
      w[i] = v[i] - theta * w[i] / rho;
    }
 
    /* stopping criteria */
 
    norm_r = phi_bar;
 
    norm_Ar = phi_bar * alpha * fabs(c);
 
    delta = sn * rho;
    gamma_bar = -cs * rho;
    rhs = phi - delta * z;
    z_bar = rhs / gamma_bar;
    norm_x = sqrt(norm_xx + NR_sqr(z_bar));
    gamma = sqrt(NR_sqr(gamma_bar) + NR_sqr(theta));
    cs = gamma_bar / gamma;
    sn = theta / gamma;
    z = rhs / gamma;
    norm_xx += NR_sqr(z);
 
    norm_d = norm_vec(w, nc) / rho;
    norm_D = norm_D + norm_d;
    cond_A = norm_d * norm_D;
 
    if (norm_r <= btol * norm_b + atol * norm_A * norm_x)
      stop1 = 1;
 
    if (norm_Ar / (norm_A * norm_r) <= atol)
      stop2 = 1;
 
    if (cond_A >= conlim)
      stop3 = 1;
 
  } while (!stop1 && !stop2 && !stop3);

  return (UT_OK);
}


/*******************************************************************************
LSQR_DMAT
An iterative method for solving least-squares/minimum norm problems of the
form Ax = b.  It requires three workspace vectors u, v and w of length nr, nc,
and nc respectively, where nr and nc are the dimensions of A.  The stopping
conditions are set by three input parameters, atol, btol, and conlim.  
Recommended values for these parameters depend on the machine value
DBL_EPSILON:

atol = DBL_EPSILON
btol = DBL_EPSILON
conlim = (1.0 / (10 * sqrt(DBL_EPSILON).

For better convergence properties, LSQR takes its inputs in double precision.
 
Reference: CC Paige, MA Saunders.  "LSQR: An Algorithm for Sparse Linear
Equations and Sparse Least Squares."  ACM Trans Math Soft, Vol 8, No 1,
March 1982, pp43-71.
 
Most of the notation used in this function is derived from this paper.
RG
*******************************************************************************/
int lsqr_dmat(double **A, int nr, int nc, double *b, double *x, double *u,
              double *v, double *w, double atol, double btol, double conlim)
{
  int    i, j;
  double alpha;
  double beta;
  double norm_A;
  double norm_b;
  double norm_B;
  double norm_D;
  double norm_Ar;
  double norm_r;
  double norm_x;
  double norm_xx;
  double norm_d;
  double cond_A;
  double phi, phi_bar;
  double rho, rho_bar;
  double gamma, gamma_bar;
  double z, z_bar;
  double delta;
  double c;
  double s;
  double theta;
  double sn;
  double cs;
  double rhs;
  int    iters;
  int    stop1, stop2, stop3;
 
  /* initialize */

  copy_dvec(b, u, nr);
 
  beta = norm_dvec(u, nr);
  scalar_div_dvec(u, nr, beta);

  left_mult_dmatrix(A, nr, nc, u, v);
  alpha = norm_dvec(v, nc);
  scalar_div_dvec(v, nc, alpha);

  copy_dvec(v, w, nc);;
  
  fast_zero_dvec(x, nc);
  
  phi_bar = beta;
  rho_bar = alpha;
 
  /* quantities needed for stopping conditions */
 
  norm_A = frobenius_norm_dmat(A, nr, nc);
  norm_b = norm_dvec(b, nr);

  cond_A = 0.0;
  norm_B = 0.0;
  norm_D = 0.0;
  norm_x = 0.0;
  norm_Ar = 0.0;

  /* quantities used for estimating the norm of x */

  cs = -1.0;
  sn = 0.0;
  z = 0.0;
  norm_xx = 0.0;

  /* now perform the iterations */
 
  iters = 0;
  stop1 = 0;
  stop2 = 0;
  stop3 = 0;
 
  do {
    iters++;
 
    /* bidiagonalization */
 
    for (i = 1; i <= nr; i++) {
      u[i] *= -alpha;
      for (j = 1; j <= nc; j++)
        u[i] += A[i][j] * v[j];
    }
 
    beta = norm_dvec(u, nr);
    scalar_div_dvec(u, nr, beta);
 
    /* before we overwrite alpha(k) with alpha(k+1) */
 
    norm_B = norm_B + NR_sqr(alpha) + NR_sqr(beta);
 
    /* continue with bidiagonalization */
 
    for (j = 1; j <= nc; j++) {
      v[j] *= -beta;
      for (i = 1; i <= nr; i++)
        v[j] += A[i][j] * u[i];
    }
 
    alpha = norm_dvec(v, nc);
    scalar_div_dvec(v, nc, alpha);
 
    /* orthogonal transformation */
 
    rho = sqrt(NR_sqr(rho_bar) + NR_sqr(beta));
    c = rho_bar / rho;
    s = beta / rho;
    theta = s * alpha;
    rho_bar = -c * alpha;
    phi = c * phi_bar;
    phi_bar *= s;
 
    /* update x, w */
 
    for (i = 1; i <= nc; i++) {
      x[i] += phi * w[i] / rho;
      w[i] = v[i] - theta * w[i] / rho;
    }
 
    /* stopping criteria */
 
    norm_r = phi_bar;
 
    norm_Ar = phi_bar * alpha * fabs(c);
 
    delta = sn * rho;
    gamma_bar = -cs * rho;
    rhs = phi - delta * z;
    z_bar = rhs / gamma_bar;
    norm_x = sqrt(norm_xx + NR_sqr(z_bar));
    gamma = sqrt(NR_sqr(gamma_bar) + NR_sqr(theta));
    cs = gamma_bar / gamma;
    sn = theta / gamma;
    z = rhs / gamma;
    norm_xx += NR_sqr(z);
 
    norm_d = norm_dvec(w, nc) / rho;
    norm_D = norm_D + norm_d;
    cond_A = norm_d * norm_D;
 
    if (norm_r <= btol * norm_b + atol * norm_A * norm_x)
      stop1 = 1;
 
    if (norm_Ar / (norm_A * norm_r) <= atol)
      stop2 = 1;
 
    if (cond_A >= conlim)
      stop3 = 1;
 
  } while (!stop1 && !stop2 && !stop3);

  return (UT_OK);
}


/*******************************************************************************
POSDEF_SYM_ALLOC_DMAT
Checks whether a symmetric matrix of doubles is positive definite or not, by 
using a Cholesky decomposition.
RG
*******************************************************************************/
int posdef_sym_alloc_dmat(double **mat, int dim, int *posdef)
{
  double *p;
  double **mat_copy;
  int     i, j, k;
  double  sum;

  /* allocate some temporary storage */
  p = NR_dvector(1, dim);
  mat_copy = NR_dmatrix(1, dim, 1, dim);

  /* copy the matrix so that it isn't destroyed */
  copy_dmat(mat, mat_copy, dim, dim);

  /* check whether it is postive definite by performing the decomposition */
  *posdef = UT_TRUE;
  for (i = 1; i <= dim; i++) { 
    for (j = i; j <= dim; j++) { 
      for (sum = mat_copy[i][j], k = i - 1; k >= 1; k--) 
	sum -= mat_copy[i][k] * mat_copy[j][k]; 
      if (i == j) { 
	if (sum <= 0.0) {
	  *posdef = UT_FALSE;
	  break;
	}
	p[i] = sqrt(sum); 
      } else 
       mat_copy[j][i] = sum / p[i];
    }
    if (*posdef == UT_FALSE)
      break;
  }

  /* clean up allocated memory */
  NR_free_dvector(p, 1, dim);
  NR_free_dmatrix(mat_copy, 1, dim, 1, dim);

  return(UT_OK);
}

/*******************************************************************************
POSDEF_SYM_ALLOC_MAT
Checks whether a symmetric matrix of floats is positive definite or not, by
using a Cholesky decomposition.
RG
*******************************************************************************/
int posdef_sym_alloc_mat(float **mat, int dim, int *posdef)
{
  float *p;
  float **mat_copy;
  int     i, j, k;
  float  sum;

  /* allocate some temporary storage */
  p = NR_vector(1, dim);
  mat_copy = NR_matrix(1, dim, 1, dim);       

  /* copy the matrix so that it isn't destroyed */
  copy_mat(mat, mat_copy, dim, dim);

  /* check whether it is postive definite by performing the decomposition */
  *posdef = UT_TRUE;
  for (i = 1; i <= dim; i++) {
    for (j = i; j <= dim; j++) {      
      for (sum = mat_copy[i][j], k = i - 1; k >= 1; k--)
        sum -= mat_copy[i][k] * mat_copy[j][k];
      if (i == j) {
        if (sum <= 0.0) {
          *posdef = UT_FALSE;
          break;
        }
        p[i] = sqrt(sum);
      } else
       mat_copy[j][i] = sum / p[i];
    }
    if (*posdef == UT_FALSE)
      break;
  }

  /* clean up allocated memory */
  NR_free_vector(p, 1, dim);
  NR_free_matrix(mat_copy, 1, dim, 1, dim);

  return(UT_OK);
}
