/*******************************************************************************
MODULE NAME
da_memory

ONE-LINE SYNOPSIS
Utility functions related to allocating and de-allocating matrices and vectors.

SCOPE OF THIS MODULE
As stated.

SEE ALSO
-

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/hmm, AG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_memory.c,v 1.1 1997/07/30 23:31:06 agray Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_memory.c,v $
 * Revision 1.1  1997/07/30 23:31:06  agray
 * Initial revision
 *
 *
 * */

/* C library */
#include <stdio.h>
#include <stdlib.h>

/* UT library */
#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"
#include "ut_memory.h"

/* NR library */
#include "nr.h"

/* DA library */

/* this module's header */
#include "da_memory.h"


/*******************************************************************************
MATRIX_AND_TRACK, CMATRIX_AND_TRACK, IMATRIX_AND_TRACK, DMATRIX_AND_TRACK
Allocate a matrix, checking for invalid arguments and keeping track of memory
use.  Analogous to malloc_and_track().
AG
*******************************************************************************/
float** matrix_and_track(int nr, int nc)

{
  int   size;
  float **mat;

  /* compute size of matrix to be allocated */
  size  = (nr + 1) * sizeof(float*);
  size += (nr * nc + 1) * sizeof(float);

  /* check size argument for validity first */
  if ((nr <= 0) || (nc <= 0))
  {
    err_printf();
    log_printf("Tried to allocate matrix of size %d x %d.\n",nr,nc);
    return ( (float**) NULL );
  }

  /* check that allocation succeeded before tracking memory usage */
  if ( (mat = (float**) NR_matrix(1,nr,1,nc)) != (float**)NULL)
  {
    ut_mem_used += size;
    ut_mem_allocated += size;
  }

  return (mat);
}

unsigned char** cmatrix_and_track(int nr, int nc)

{
  int   size;
  unsigned char **mat;

  /* compute size of matrix to be allocated */
  size  = (nr + 1) * sizeof(unsigned char*);
  size += (nr * nc + 1) * sizeof(unsigned char);

  /* check size argument for validity first */
  if ((nr <= 0) || (nc <= 0))
  {
    err_printf();
    log_printf("Tried to allocate matrix of size %d x %d.\n",nr,nc);
    return ( (unsigned char**) NULL );
  }

  /* check that allocation succeeded before tracking memory usage */
  if ( (mat = (unsigned char**) NR_cmatrix(1,nr,1,nc)) != (unsigned char**)NULL)
  {
    ut_mem_used += size;
    ut_mem_allocated += size;
  }

  return (mat);
}

int** imatrix_and_track(int nr, int nc)

{
  int   size;
  int **mat;

  /* compute size of matrix to be allocated */
  size  = (nr + 1) * sizeof(int*);
  size += (nr * nc + 1) * sizeof(int);

  /* check size argument for validity first */
  if ((nr <= 0) || (nc <= 0))
  {
    err_printf();
    log_printf("Tried to allocate matrix of size %d x %d.\n",nr,nc);
    return ( (int**) NULL );
  }

  /* check that allocation succeeded before tracking memory usage */
  if ( (mat = (int**) NR_imatrix(1,nr,1,nc)) != (int**)NULL)
  {
    ut_mem_used += size;
    ut_mem_allocated += size;
  }

  return (mat);
}

double** dmatrix_and_track(int nr, int nc)

{
  int   size;
  double **mat;

  /* compute size of matrix to be allocated */
  size  = (nr + 1) * sizeof(double*);
  size += (nr * nc + 1) * sizeof(double);

  /* check size argument for validity first */
  if ((nr <= 0) || (nc <= 0))
  {
    err_printf();
    log_printf("Tried to allocate matrix of size %d x %d.\n",nr,nc);
    return ( (double**) NULL );
  }

  /* check that allocation succeeded before tracking memory usage */
  if ( (mat = (double**) NR_dmatrix(1,nr,1,nc)) != (double**)NULL)
  {
    ut_mem_used += size;
    ut_mem_allocated += size;
  }

  return (mat);
}


/*******************************************************************************
FREE_MATRIX_AND_TRACK, FREE_CMATRIX_AND_TRACK, FREE_IMATRIX_AND_TRACK, 
  FREE_DMATRIX_AND_TRACK
De-allocate a matrix, checking for invalid arguments and keeping track of memory
use.  Analogous to free_and_track().
AG
*******************************************************************************/
void free_matrix_and_track(float **mat, int nr, int nc)

{
  int size;

  /* compute size of matrix to be de-allocated */
  size  = (nr + 1) * sizeof(float*);
  size += (nr * nc + 1) * sizeof(float);

  /* check pointer argument for validity first */
  if (mat == (float**) NULL)
  {
    err_printf();
    log_printf("Tried to free a NULL matrix.\n");
    return;
  }
  else
  {
    NR_free_matrix(mat, 1, nr, 1, nc);
    ut_mem_used -= size;
  }
}

void free_cmatrix_and_track(unsigned char **mat, int nr, int nc)

{
  int size;

  /* compute size of matrix to be de-allocated */
  size  = (nr + 1) * sizeof(unsigned char*);
  size += (nr * nc + 1) * sizeof(unsigned char);

  /* check pointer argument for validity first */
  if (mat == (unsigned char**) NULL)
  {
    err_printf();
    log_printf("Tried to free a NULL matrix.\n");
    return;
  }
  else
  {
    NR_free_cmatrix(mat, 1, nr, 1, nc);
    ut_mem_used -= size;
  }
}

void free_imatrix_and_track(int **mat, int nr, int nc)

{
  int size;

  /* compute size of matrix to be de-allocated */
  size  = (nr + 1) * sizeof(int*);
  size += (nr * nc + 1) * sizeof(int);

  /* check pointer argument for validity first */
  if (mat == (int**) NULL)
  {
    err_printf();
    log_printf("Tried to free a NULL matrix.\n");
    return;
  }
  else
  {
    NR_free_imatrix(mat, 1, nr, 1, nc);
    ut_mem_used -= size;
  }
}

void free_dmatrix_and_track(double **mat, int nr, int nc)

{
  int size;

  /* compute size of matrix to be de-allocated */
  size  = (nr + 1) * sizeof(double*);
  size += (nr * nc + 1) * sizeof(double);

  /* check pointer argument for validity first */
  if (mat == (double**) NULL)
  {
    err_printf();
    log_printf("Tried to free a NULL matrix.\n");
    return;
  }
  else
  {
    NR_free_dmatrix(mat, 1, nr, 1, nc);
    ut_mem_used -= size;
  }
}

