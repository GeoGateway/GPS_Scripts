/*******************************************************************************
MODULE NAME
da_io

ONE-LINE SYNOPSIS
General functions related to data input and output.

SCOPE OF THIS MODULE
All functions in the library pertaining to reading, writing, and printing data
structures should be in this module.

SEE ALSO
-

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/sc, AG.
2. /proj/cooltools/qf, RG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_io.c,v 1.31 2000/03/31 01:46:24 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_io.c,v $
 * Revision 1.31  2000/03/31 01:46:24  granat
 * fixed some compile-time bugs
 *
 * Revision 1.30  2000/03/31 01:38:47  granat
 * added double precision sets of matrices functions
 *
 * Revision 1.29  1998/07/24 16:14:01  granat
 * fixed bugs in matrix input functions relating to types
 * eliminated some unused index variables
 *
 * Revision 1.28  1998/05/07 23:51:07  granat
 * made changes to include new da_util module
 *
 * Revision 1.27  1998/05/01 20:50:20  granat
 * added functions to read/write unsigned char matrices into/from matrices
 * of doubles and floats
 *
 * Revision 1.26  1997.20/21 14:52:24  granat
 * added some functions to read/write vectors of different data types
 *
 * Revision 1.25  1997/09/04 20:21:39  granat
 * fixed bug in skip_header
 *
 * Revision 1.24  1997/08/07 17:36:52  granat
 * fixed bugs from automated scripting of "nr_" convention
 *
 * Revision 1.23  1997/07/31 21:46:14  agray
 * added #include "da_memory.h"
 *
 * Revision 1.22  1997/06/20 22:12:47  granat
 * fixed bugs introduced by auto-editing scripts
 *
 * Revision 1.21  1997/06/02 15:35:42  granat
 * changed to use new NR naming convention
 *
 * Revision 1.20  1997/04/02 03:05:50  agray
 * added assert.h and string.h inclusions.  to fix compile problems introduced by
 * assert().
 *
 * Revision 1.19  1997/03/27 23:03:14  granat
 * added functions skip_header() and read_vicar()
 *
 * Revision 1.18  1997/01/29 21:54:01  agray
 * new formatting, cleaning up debugging output using ut_output,
 * added write_matrix_transpose().
 *
 * Revision 1.17  1996.20/31 02:13:49  agray
 * renamed from "da_data" to "da_io";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.16  1996/09/25 00:50:14  granat
 * Changed print_cmatrix and read_matrix to use integer conversion becausened character conversion is unavailable.
 *
 * Revision 1.15  1996/09/18 23:40:11  agray
 * cosmetic changes to comments.
 *
 * Revision 1.14  1996/09/18 16:33:41  granat
 * Fixed "char" matrix functions so that they use the "unsigned char" type
 *
 * Revision 1.13  1996/09/06 22:16:37  agray
 * fixed read_irow().
 *
 * Revision 1.12  1996/08/28 20:31:19  agray
 * fixed some bugs, removed duplicated copy of write_matrix().
 *
 * Revision 1.11  1996/08/28 19:52:30  agray
 * changed read/write_data() to read/write_lc_matrix(); added read/write/print_
 * indexed_matrix(); added write_matrix().
 *
 * Revision 1.20  1996/07/17 20:42:36  agray
 * cosmetic.
 *
 * Revision 1.9  1996/07/15 18:16:52  agray
 * added  write/read_bin_row/col/irow/icol();
 * renamed and tweaked read_bin_subset2matrix().
 *
 * Revision 1.8  1996/07/15 17:25:02  agray
 * rearranged order of functions; added read_row/col/irow/icol().
 *
 * Revision 1.7  1996/07/13 01:06:49  agray
 * moved in print/write_row/col/irow/icol() from da_linalg module
 *
 * Revision 1.6  1996/07/11 17:59:00  agray
 * added read_matrix_transpose(), write_comment().  some cosmetic changes.
 *
 * Revision 1.5  1996/05.20 17:16:15  granat
 * added read_subset2matrix and read_binary2matrix functions
 * rg
 *
 * Revision 1.4  1996/04/29 20:25:14  granat
 * cvector function.
 * rg
 *
 * Revision 1.3  1996/03/01 00:17:28  agray
 * added many i/o functions: read_matrix() (versions for floats, chars, ints, doubles),
 * write_matrix() (versions for f,c,i,d), print_matrix() (f,c,i,d), read_bin_matrix
 * (f,c,i,d), write_bin_matrix() (f,c,i,d), read_data() and write_data() (both moved
 * from da_linalg.c)
 * ag
 *
 * Revision 1.2  1996/02/29  02:32:08  agray
 * moved write_bin_matrix() and read_bin_matrix() to this module, added type
 * argument to them
 * ag
 *
 * Revision 1.1  1996/02/28  04:55:02  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* UT library */
#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"
#include "ut_file_io.h"
#include "ut_memory.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_util.h"
#include "da_memory.h"

/* this module's header */
#include "da_io.h"


/*******************************************************************************
READ_MATRIX
Read in ascii file containing a matrix.
Allocates the matrix of floats that will hold the data.
AG
*******************************************************************************/
int read_matrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (float**) matrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++)
      fscanf(fp, "%g ", (float*) &(*vals)[i][j]); 

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_CMATRIX
Read in ascii file containing a matrix.
Allocates the matrix of unsigned chars that will hold the data.
AG
*******************************************************************************/
int read_cmatrix(infile, nr, nc, vals)

  char           *infile;   /* name of file to read */
  int            nr;        /* number of rows in matrix */
  int            nc;        /* number of columns in matrix */
  unsigned char  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  int    intval;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (unsigned char**) cmatrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%d ", &intval);
      (*vals)[i][j] = (unsigned char) intval;
    }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_IMATRIX
Read in ascii file containing a matrix.
Allocates the matrix of ints that will hold the data.
AG
*******************************************************************************/
int read_imatrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  int    ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (int**) imatrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++)
      fscanf(fp, "%d ", (int*) &(*vals)[i][j]);

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_DMATRIX
Read in ascii file containing a matrix.
Allocates the matrix of doubles that will hold the data.
AG
*******************************************************************************/
int read_dmatrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  double ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  double doubleval;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (double**) dmatrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%lg ", &doubleval); 
      (*vals)[i][j] = doubleval;
    }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_CMATRIX_AS_MATRIX
Read in ascii file containing a matrix of unsigned chars into a matrix of
floats.  Allocates the matrix of floats that will hold the data.
AG
*******************************************************************************/
int read_cmatrix_as_matrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (float**) matrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++)
      fscanf(fp, "%g ", (float*) &(*vals)[i][j]); 

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_CMATRIX_AS_DMATRIX
Read in ascii file containing a matrix of unsigned chars into a matrix of
doubles.  Allocates the matrix of doubles that will hold the data.
AG
*******************************************************************************/
int read_cmatrix_as_dmatrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  double  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (double**) dmatrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++)
      fscanf(fp, "%lg ", (double*) &(*vals)[i][j]); 

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
PRINT_MATRIX
Print the contents of a matrix of floats.
AG
*******************************************************************************/
int print_matrix(stream, nr, nc, mat)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  float  **mat;   /* matrix to be printed */
{
  int i,j;

  for (i = 1; i <= nr; i++) 
  {    
	for (j = 1;j <= nc;j++)
      fprintf(stream, "%g ", mat[i][j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_CMATRIX
Print the contents of a matrix of unsigned chars.
AG
*******************************************************************************/
int print_cmatrix(stream, nr, nc, mat)

  FILE*          stream;  /* stream to print to, e.g. stdout, file pointer */
  int            nr, nc;  /* number of rows, columns of input matrix */
  unsigned char  **mat;   /* matrix to be printed */
{
  int i,j;

  for (i = 1; i <= nr; i++)
  {
    for (j = 1;j <= nc;j++)
      fprintf(stream, "%d ", mat[i][j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_IMATRIX
Print the contents of a matrix of ints.
AG
*******************************************************************************/
int print_imatrix(stream, nr, nc, mat)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  int    **mat;   /* matrix to be printed */
{
  int i,j;

  for (i = 1; i <= nr; i++)
  {
    for (j = 1;j <= nc;j++)
      fprintf(stream, "%d ", mat[i][j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_DMATRIX
Print the contents of a matrix of doubles.
AG
*******************************************************************************/
int print_dmatrix(stream, nr, nc, mat)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  double **mat;   /* matrix to be printed */
{
  int i,j;

  for (i = 1; i <= nr; i++)
  {
    for (j = 1;j <= nc;j++)
      fprintf(stream, "%.20g ", mat[i][j]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}


/*******************************************************************************
WRITE_MATRIX
Write a matrix of floats to a file.
AG
*******************************************************************************/
int write_matrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  float  **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);

  /* print matrix contents to file */
  print_matrix(fp, nr, nc, mat);

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
WRITE_CMATRIX
Write a matrix of unsigned chars to a file.
AG
*******************************************************************************/
int write_cmatrix(outfile, nr, nc, mat, mode)

  char           *outfile;   /* name of file to write to */
  int            nr, nc;     /* number of rows, columns in matrix */
  unsigned char  **mat;      /* matrix to write */
  char           *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);

  /* print matrix contents to file */
  print_cmatrix(fp, nr, nc, mat);

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
WRITE_IMATRIX
Write a matrix of ints to a file.
AG
*******************************************************************************/
int write_imatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  int    **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);

  /* print matrix contents to file */
  print_imatrix(fp, nr, nc, mat);

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
WRITE_DMATRIX
Write a matrix of doubles to a file.
AG
*******************************************************************************/
int write_dmatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  double **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);

  /* print matrix contents to file */
  print_dmatrix(fp, nr, nc, mat);

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
PRINT_MATRIX_AS_CMATRIX
Prints a matrix of floats as unsigned characters from 0 to 255.  Values
over 255 are thresholded to 255.  Negative values are thresholded to 0.
RG
*******************************************************************************/
int print_matrix_as_cmatrix( FILE *stream, int nr, int nc, float **m )
{
  int    i, j;
 
  for (i = 1; i <= nr; i++) {
    for (j = 1; j <= nc; j++) {
      if (m[i][j] > 255.0)
        fprintf( stream, "255 ");
      else if (m[i][j] < 0.0)
        fprintf( stream, "0 ");
      else
        fprintf( stream, "%d ", (int) m[i][j] );
    }
    fprintf( stream, "\n" );
  }
 
  return( UT_OK );
}


/*******************************************************************************
PRINT_DMATRIX_AS_CMATRIX
Prints a matrix of doubles as unsigned characters from 0 to 255.  Values
over 255 are thresholded to 255.  Negative values are thresholded to 0.
RG
*******************************************************************************/
int print_dmatrix_as_cmatrix( FILE *stream, int nr, int nc, double **m )
{
  int    i, j;
 
  for (i = 1; i <= nr; i++) {
    for (j = 1; j <= nc; j++) {
      if (m[i][j] > 255.0)
        fprintf( stream, "255 ");
      else if (m[i][j] < 0.0)
        fprintf( stream, "0 ");
      else
        fprintf( stream, "%d ", (int) m[i][j] );
    }
    fprintf( stream, "\n" );
  }
 
  return( UT_OK );
}


/*******************************************************************************
WRITE_MATRIX_AS_CMATRIX
Writes a matrix of floats to a file as unsigned characters from 0 to 255.
Values over 255 are thresholded to 255.  Negative values are thresholded to 0.
RG
*******************************************************************************/
int write_matrix_as_cmatrix( char *outfile, int nr, int nc, float **m,
                             char *mode )
{
  FILE *fp;
 
  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);
 
  /* print matrix contents to file */
  print_matrix_as_cmatrix(fp, nr, nc, m);
 
  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}
 
 
/*******************************************************************************
WRITE_DMATRIX_AS_CMATRIX
Writes a matrix of doubles to a file as unsigned characters from 0 to 255.
Values over 255 are thresholded to 255.  Negative values are thresholded to 0.
RG
*******************************************************************************/
int write_dmatrix_as_cmatrix( char *outfile, int nr, int nc, double **m,
                              char *mode )
{
  FILE *fp;
 
  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);
 
  /* print matrix contents to file */
  print_dmatrix_as_cmatrix(fp, nr, nc, m);
 
  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}
 
 
/*******************************************************************************
READ_BIN_MATRIX
Read in binary file of floats containing a matrix.
Allocates the matrix of floats that will hold the data.
AG
*******************************************************************************/
int read_bin_matrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (float**) matrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++) {
    fread( (float*) &((*vals)[i][1]), sizeof(float), nc, fp );
  }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_CMATRIX
Read in binary file of unsigned chars containing a matrix.
Allocates the matrix of unsigned chars that will hold the data.
AG
*******************************************************************************/
int read_bin_cmatrix(infile, nr, nc, vals)

  char          *infile;   /* name of file to read */
  int           nr;        /* number of rows in matrix */
  int           nc;        /* number of columns in matrix */
  unsigned char ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_cmatrix(1, nr, 1, nc);
  if (*vals == (unsigned char**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++) {
    fread( (unsigned char *) &((*vals)[i][1]), sizeof(unsigned char), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_IMATRIX
Read in binary file of ints containing a matrix.
Allocates the matrix of ints that will hold the data.
AG
*******************************************************************************/
int read_bin_imatrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  int    ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_imatrix(1, nr, 1, nc);
  if (*vals == (int**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++) {
    fread( (int*) &((*vals)[i][1]), sizeof(int), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_SMATRIX
Read in binary file of shorts containing a matrix.
Allocates the matrix of ints that will hold the data.
RG
*******************************************************************************/
int read_bin_smatrix(infile, nr, nc, vals)
 
  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  int    ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;
 
  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }
 
  /* allocate data matrix */
  *vals = NR_imatrix(1, nr, 1, nc);
  if (*vals == (int**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
 
  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++) {
    fread( (short*) &((*vals)[i][1]), sizeof(short), nc, fp );
  }
 
  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_DMATRIX
Read in binary file of doubles containing a matrix.
Allocates the matrix of doubles that will hold the data.
AG
*******************************************************************************/
int read_bin_dmatrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  double ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_dmatrix(1, nr, 1, nc);
  if (*vals == (double**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++) {
    fread( (double*) &((*vals)[i][1]), sizeof(double), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_CMATRIX_AS_MATRIX
Read in binary file of unsigned chars containing a matrix, and read it into
a matrix of floats.  Allocates the matrix of floats that will hold the data.
RG
*******************************************************************************/
int read_bin_cmatrix_as_matrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  unsigned char charval;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (float**) matrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++) {
      fread(&charval, sizeof(unsigned char), 1, fp);
      (*vals)[i][j] = (float) charval;
    }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_CMATRIX_AS_DMATRIX
Read in binary file of unsigned chars containing a matrix, and read it into
a matrix of doubles.  Allocates the matrix of doubles that will hold the data.
RG
*******************************************************************************/
int read_bin_cmatrix_as_dmatrix(infile, nr, nc, vals)

  char     *infile;   /* name of file to read */
  int      nr;        /* number of rows in matrix */
  int      nc;        /* number of columns in matrix */
  double   ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  unsigned char charval;

  /* open file for reading */
  fp = fopen_return_if_fail(infile, "r", UT_ERROR);

  /* allocate data matrix */
  *vals = (double**) dmatrix_return_if_fail(nr, nc, UT_ERROR);

  /* read from file into the matrix, in binary format */
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++) {
      fread(&charval, sizeof(unsigned char), 1, fp);
      (*vals)[i][j] = (double) charval;
    }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_MATRIX
Write a matrix to a file, in binary form (floats).
AG
*******************************************************************************/
int write_bin_matrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  float  **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) {
    fwrite( (float*) &mat[i][1], sizeof(float), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_CMATRIX
Write a matrix to a file, in binary form (unsigned chars).
AG
*******************************************************************************/
int write_bin_cmatrix(outfile, nr, nc, mat, mode)

  char          *outfile;   /* name of file to write to */
  int           nr, nc;     /* number of rows, columns in matrix */
  unsigned char **mat;      /* matrix to write */
  char          *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) {
    fwrite( (unsigned char *) &mat[i][1], sizeof(unsigned char), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_IMATRIX
Write a matrix to a file, in binary form (ints).
AG
*******************************************************************************/
int write_bin_imatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  int    **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) {
    fwrite( (int*) &mat[i][1], sizeof(int), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_DMATRIX
Write a matrix to a file, in binary form (doubles).
AG
*******************************************************************************/
int write_bin_dmatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  double **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) {
    fwrite( (double*) &mat[i][1], sizeof(double), nc, fp );
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_MATRIX_AS_BIN_CMATRIX
Write a matrix to a file in binarie form (unsigned chars).
RG
*******************************************************************************/
int write_matrix_as_bin_cmatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  float  **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i, j;
  unsigned char charval;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++) {
      if (mat[i][j] > 255.0)
	charval = 255;
      else if (mat[i][j] < 0.0)
	charval = 0;
      else
	charval = (unsigned char) mat[i][j];
      fwrite( &charval, sizeof(unsigned char), 1, fp );
    }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_DMATRIX_AS_BIN_CMATRIX
Write a matrix to a file in binarie form (unsigned chars).
RG
*******************************************************************************/
int write_dmatrix_as_bin_cmatrix(outfile, nr, nc, mat, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  double **mat;      /* matrix to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;
  int     i, j;
  unsigned char charval;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file, in binary */
  for (i = 1; i <= nr; i++) 
    for (j = 1; j <= nc; j++) {
      if (mat[i][j] > 255.0)
	charval = 255;
      else if (mat[i][j] < 0.0)
	charval = 0;
      else
	charval = (unsigned char) mat[i][j];
      fwrite( &charval, sizeof(unsigned char), 1, fp );
    }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_ROW
Read in ascii file containing a vector in row form.
Allocates the vector of floats that will hold the data.
AG
*******************************************************************************/
int read_row(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  float  **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_vector(1, dim);
  if (*vec == (float*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector */
  for (i = 1; i <= dim; i++) {
    fscanf(fp, "%g ", (float*) &(*vec)[i]);
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_COL
Read in ascii file containing a vector in column form.
Allocates the vector of floats that will hold the data.
Implemented exactly the same as read_row().
AG
*******************************************************************************/
int read_col(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  float  **vec;     /* pointer to vector to be filled in */
{
  return( read_row(infile, dim, vec) );
}


/*******************************************************************************
READ_IROW
Read in ascii file containing a vector of integers in row form.
Allocates the vector of integers that will hold the data.
AG
*******************************************************************************/
int read_irow(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  int    **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_ivector(1, dim);
  if (*vec == (int*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector */
  for (i = 1; i <= dim; i++) {
    fscanf(fp, "%d ", (int*) &(*vec)[i]);
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_ICOL
Read in ascii file containing a vector of integers in column form.
Allocates the vector of integers that will hold the data.
Implemented exactly the same as read_irow().
AG
*******************************************************************************/
int read_icol(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  int    **vec;     /* pointer to vector to be filled in */
{
  return( read_irow(infile, dim, vec) );
}


/*******************************************************************************
READ_DROW
Read in ascii file containing a vector of doubles in row form.
Allocates the vector of integers that will hold the data.
AG
*******************************************************************************/
int read_drow(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  double **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;
  int    i;
  double doubleval;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_dvector(1, dim);
  if (*vec == (double*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector */
  for (i = 1; i <= dim; i++) {
    fscanf(fp, "%lg ", &doubleval);
    (*vec)[i] = doubleval;
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_DCOL
Read in ascii file containing a vector of doubles in column form.
Allocates the vector of integers that will hold the data.
Implemented exactly the same as read_irow().
AG
*******************************************************************************/
int read_dcol(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  double **vec;     /* pointer to vector to be filled in */
{
  return( read_drow(infile, dim, vec) );
}


/*******************************************************************************
READ_CROW
Read in ascii file containing a vector of unsigned chars in row form.
Allocates the vector of integers that will hold the data.
AG
*******************************************************************************/
int read_crow(infile, dim, vec)

  char              *infile;   /* name of file to read */
  int               dim;       /* dimension of vector */
  unsigned char     **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;
  int    i;
  int    intval;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_cvector(1, dim);
  if (*vec == (unsigned char*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector */
  for (i = 1; i <= dim; i++) {
    fscanf(fp, "%d ", &intval);
    (*vec)[i] = intval;
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_CCOL
Read in ascii file containing a vector of unsigned chars in column form.
Allocates the vector of integers that will hold the data.
Implemented exactly the same as read_irow().
AG
*******************************************************************************/
int read_ccol(infile, dim, vec)

  char           *infile;   /* name of file to read */
  int            dim;       /* dimension of vector */
  unsigned char  **vec;     /* pointer to vector to be filled in */
{
  return( read_crow(infile, dim, vec) );
}


/*******************************************************************************
PRINT_ROW
Print the contents of a vector, in row form.
AG
*******************************************************************************/
int print_row(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the row vector to be printed */
  float  *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%g ", v[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}


/*******************************************************************************
PRINT_COL
Print the contents of a vector, in column form.
AG
*******************************************************************************/
int print_col(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the column vector to be printed */
  float  *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%g\n", v[i]);
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_IROW
Print the contents of an integer vector, in row form.
AG
*******************************************************************************/
int print_irow(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the row vector to be printed */
  int    *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%d ", v[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}


/*******************************************************************************
PRINT_ICOL
Print the contents of an integer vector, in column form.
AG
*******************************************************************************/
int print_icol(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the column vector to be printed */
  int    *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%d\n", v[i]);
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_DROW
Print the contents of an double vector, in row form.
AG
*******************************************************************************/
int print_drow(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the row vector to be printed */
  double *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%.20g ", v[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}


/*******************************************************************************
PRINT_DCOL
Print the contents of an double vector, in column form.
AG
*******************************************************************************/
int print_dcol(stream, dim, v)

  FILE*  stream; /* stream to print to, e.g. stdout, file pointer */
  int    dim;    /* length of the column vector to be printed */
  double *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%.20g\n", v[i]);
  }

  return (UT_OK);
}


/*******************************************************************************
PRINT_CROW
Print the contents of an unsigned char vector, in row form.
AG
*******************************************************************************/
int print_crow(stream, dim, v)

  FILE*            stream; /* stream to print to, e.g. stdout, file pointer */
  int              dim;    /* length of the row vector to be printed */
  unsigned char    *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%d ", (int) v[i]);
  }
  fprintf(stream, "\n");

  return (UT_OK);
}


/*******************************************************************************
PRINT_CCOL
Print the contents of an unsigned char vector, in column form.
AG
*******************************************************************************/
int print_ccol(stream, dim, v)

  FILE*            stream; /* stream to print to, e.g. stdout, file pointer */
  int              dim;    /* length of the column vector to be printed */
  unsigned char    *v;     /* vector to be printed */
{
  int     i;

  for (i = 1; i <= dim; i++) {
      fprintf(stream, "%d\n", (int) v[i]);
  }

  return (UT_OK);
}


/*******************************************************************************
WRITE_ROW
Write a vector to a file, in row form.
AG
*******************************************************************************/
int write_row(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row */
  float  *vec;      /* row vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_row(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_COL
Write a vector to a file, in column form.
AG
*******************************************************************************/
int write_col(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the column */
  float  *vec;      /* column vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_col(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_IROW
Write an integer vector to a file, in row form.
AG
*******************************************************************************/
int write_irow(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row */
  int    *vec;      /* row vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_irow(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_ICOL
Write an integer vector to a file, in column form.
AG
*******************************************************************************/
int write_icol(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the column */
  int    *vec;      /* column vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_icol(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_DROW
Write a vector of doubles to a file, in row form.
AG
*******************************************************************************/
int write_drow(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row */
  double *vec;      /* row vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_drow(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_DCOL
Write an integer vector to a file, in column form.
AG
*******************************************************************************/
int write_dcol(outfile, dim, vec, mode)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the column */
  double *vec;      /* column vector to be written */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_dcol(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_CROW
Write a vector of unsigned chars to a file, in row form.
AG
*******************************************************************************/
int write_crow(outfile, dim, vec, mode)

  char            *outfile;  /* name of file to write to */
  int             dim;       /* length of the row */
  unsigned char   *vec;      /* row vector to be written */
  char            *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_crow(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_CCOL
Write an unsigned char vector to a file, in column form.
AG
*******************************************************************************/
int write_ccol(outfile, dim, vec, mode)

  char             *outfile;  /* name of file to write to */
  int              dim;       /* length of the column */
  unsigned char    *vec;      /* column vector to be written */
  char             *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_ccol(fp, dim, vec);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_ROW
Read in binary file of floats containing a vector in row form.
Allocates the vector of floats that will hold the data.
AG
*******************************************************************************/
int read_bin_row(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  float  **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_vector(1, dim);
  if (*vec == (float*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector, in binary format */
  fread( (float*) &((*vec)[1]), sizeof(float), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_COL
Read in binary file of floats containing a vector in column form.
Allocates the vector of floats that will hold the data.
Implemented exactly the same as read_bin_row() since there is no
distinction in binary format.
AG
*******************************************************************************/
int read_bin_col(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  float  **vec;     /* pointer to vector to be filled in */
{
  return( read_bin_row(infile, dim, vec) );
}


/*******************************************************************************
READ_BIN_IROW
Read in binary file of integers containing a vector in row form.
Allocates the vector of integers that will hold the data.
AG
*******************************************************************************/
int read_bin_irow(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  int    **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_ivector(1, dim);
  if (*vec == (int*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector, in binary format */
  fread( (int*) &((*vec)[1]), sizeof(int), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_ICOL
Read in binary file of integers containing a vector in column form.
Allocates the vector of integers that will hold the data.
Implemented exactly the same as read_bin_irow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int read_bin_icol(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  int    **vec;     /* pointer to vector to be filled in */
{
  return( read_bin_irow(infile, dim, vec) );
}


/*******************************************************************************
READ_BIN_DROW
Read in binary file of doubles containing a vector in row form.
Allocates the vector of doubles that will hold the data.
AG
*******************************************************************************/
int read_bin_drow(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  double **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_dvector(1, dim);
  if (*vec == (double*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector, in binary format */
  fread( (double*) &((*vec)[1]), sizeof(double), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_DCOL
Read in binary file of doubles containing a vector in column form.
Allocates the vector of doubles that will hold the data.
Implemented exactly the same as read_bin_drow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int read_bin_dcol(infile, dim, vec)

  char   *infile;   /* name of file to read */
  int    dim;       /* dimension of vector */
  double **vec;     /* pointer to vector to be filled in */
{
  return( read_bin_drow(infile, dim, vec) );
}


/*******************************************************************************
READ_BIN_CROW
Read in binary file of unsigned chars containing a vector in row form.
Allocates the vector of unsigned chars that will hold the data.
AG
*******************************************************************************/
int read_bin_crow(infile, dim, vec)

  char           *infile;   /* name of file to read */
  int            dim;       /* dimension of vector */
  unsigned char  **vec;     /* pointer to vector to be filled in */
{
  FILE   *fp;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data vector */
  *vec = NR_cvector(1, dim);
  if (*vec == (unsigned char*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the vector, in binary format */
  fread( (unsigned char*) &((*vec)[1]), sizeof(unsigned char), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_BIN_CCOL
Read in binary file of unsigned chars containing a vector in column form.
Allocates the vector of unsigned chars that will hold the data.
Implemented exactly the same as read_bin_crow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int read_bin_ccol(infile, dim, vec)

  char          *infile;   /* name of file to read */
  int           dim;       /* dimension of vector */
  unsigned char **vec;     /* pointer to vector to be filled in */
{
  return( read_bin_crow(infile, dim, vec) );
}


/*******************************************************************************
WRITE_BIN_ROW
Write a row vector to a file, in binary form (floats).
AG
*******************************************************************************/
int write_bin_row(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  float  *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file, in binary */
  fwrite( (float*) &(vec[1]), sizeof(float), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_COL
Write a column vector to a file, in binary form (floats).
Implemented exactly the same as read_bin_row() since there is no
distinction in binary format.
AG
*******************************************************************************/
int write_bin_col(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  float  *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  return( write_bin_row(outfile, dim, vec, mode) );
}


/*******************************************************************************
WRITE_BIN_IROW
Write a row vector to a file, in binary form (integers).
AG
*******************************************************************************/
int write_bin_irow(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  int    *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file, in binary */
  fwrite( (int*) &(vec[1]), sizeof(int), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_ICOL
Write a column vector to a file, in binary form (integers).
Implemented exactly the same as read_bin_irow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int write_bin_icol(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  int    *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  return( write_bin_irow(outfile, dim, vec, mode) );
}


/*******************************************************************************
WRITE_BIN_DROW
Write a row vector to a file, in binary form (doubles).
AG
*******************************************************************************/
int write_bin_drow(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  double *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file, in binary */
  fwrite( (double*) &(vec[1]), sizeof(double), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_DCOL
Write a column vector to a file, in binary form (doubles).
Implemented exactly the same as read_bin_drow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int write_bin_dcol(outfile, dim, vec, mode)

  char   *outfile;   /* name of file to write to */
  int    dim;        /* dimension of vector */
  double *vec;       /* vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  return( write_bin_drow(outfile, dim, vec, mode) );
}


/*******************************************************************************
WRITE_BIN_CROW
Write a row vector to a file, in binary form (unsigned char).
AG
*******************************************************************************/
int write_bin_crow(outfile, dim, vec, mode)

  char          *outfile;   /* name of file to write to */
  int           dim;        /* dimension of vector */
  unsigned char *vec;       /* vector to write */
  char          *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file, in binary */
  fwrite( (unsigned char*) &(vec[1]), sizeof(unsigned char), dim, fp );

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_BIN_CCOL
Write a column vector to a file, in binary form (unsigned chars).
Implemented exactly the same as read_bin_crow() since there is no
distinction in binary format.
AG
*******************************************************************************/
int write_bin_ccol(outfile, dim, vec, mode)

  char          *outfile;   /* name of file to write to */
  int           dim;        /* dimension of vector */
  unsigned char *vec;       /* vector to write */
  char          *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  return( write_bin_crow(outfile, dim, vec, mode) );
}


/*******************************************************************************
READ_INDEXED_MATRIX
Read in ascii file containing a matrix having the first column as an index
column.
Allocates a matrix of floats that will hold the data and a vector of floats
which will hold the indices.
AG
*******************************************************************************/
int read_indexed_matrix(infile, nr, nc, vals, index)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix, excluding index column */
  float  ***vals;   /* ptr to matrix to be filled in with data */
  float  **index;   /* ptr to vector to be filled in with index column */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_matrix(1, nr, 1, nc);
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* allocate index vector */
  *index = NR_vector(1, nr);
  if (*index == (float*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    fscanf(fp, "%g ", (float*) &(*index)[i]);
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%g ", (float*) &(*vals)[i][j]);
    }
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
PRINT_INDEXED_MATRIX
Print the contents of an indexed matrix, i.e. a matrix of floats and a 
corresponding index vector of floats.
AG
*******************************************************************************/
int print_indexed_matrix(stream, nr, nc, mat, index)

  FILE*  stream;  /* stream to print to, e.g. stdout, file pointer */
  int    nr, nc;  /* number of rows, columns of input matrix */
  float  **mat;   /* matrix to be printed */
  float  *index;  /* matrix to be printed */
{
  int i,j;

  for (i = 1; i <= nr; i++) {
    fprintf(stream, "%g ", index[i]);
    for (j = 1;j <= nc;j++) {
      fprintf(stream, "%g ", mat[i][j]);
    }
    fprintf(stream, "\n");
  }

  return (UT_OK);
}


/*******************************************************************************
WRITE_INDEXED_MATRIX
Write the contents of an indexed matrix, i.e. a matrix of floats and a 
corresponding index vector of floats, to a file.
AG
*******************************************************************************/
int write_indexed_matrix(outfile, nr, nc, mat, index, mode)

  char   *outfile;   /* name of file to write to */
  int    nr, nc;     /* number of rows, columns in matrix */
  float  **mat;      /* matrix to write */
  float  *index;     /* index vector to write */
  char   *mode;      /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file */
  print_indexed_matrix(fp, nr, nc, mat, index);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_INDEXED_COL
Read in ascii file containing two columns, the first containing arbitrary 
index values, and the second constituting a vector of float values.

Allocates a vector of floats that will hold the data and a vector of floats
which will hold the indices.
AG
*******************************************************************************/
int read_indexed_col(char *infile, int dim, float **vals, float **index)

{
  FILE   *fp;
  int    i;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_vector(1, dim);
  if (*vals == (float*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* allocate index vector */
  *index = NR_vector(1, dim);
  if (*index == (float*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (i = 1; i <= dim; i++)
  {
    fscanf(fp, "%g ", (float*) &(*index)[i]);
    fscanf(fp, "%g ", (float*) &(*vals)[i]);
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
PRINT_INDEXED_COL
Print the contents of an indexed vector, i.e. a vector of floats and a 
corresponding index vector of floats, to a stream in column format.
AG
*******************************************************************************/
int print_indexed_col(FILE *stream, int dim, float *vec, float *index)

{
  int i;

  for (i = 1; i <= dim; i++) 
  {
    fprintf(stream, "%g ", index[i]);
    fprintf(stream, "%g\n", vec[i]);
  }

  return (UT_OK);
}


/*******************************************************************************
WRITE_INDEXED_COL
Write the contents of an indexed vector, i.e. a vector of floats and a 
corresponding index vector of floats, to a file in column format.
AG
*******************************************************************************/
int write_indexed_col(char *outfile, int dim, float *vec, float *index, 
                      char *mode)

{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(outfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print vector contents to file */
  print_indexed_col(fp, dim, vec, index);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_LC_MATRIX
Read in ascii file containing the rows of data.
Assumes that the file is in the standard dataset format, with line count at 
top, separated by spaces.
Allocates the matrix of floats that will hold the data.
AG
*******************************************************************************/
int read_lc_matrix(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    *nr;       /* ptr to number of rows in data, to be filled in */
  int    nc;        /* number of columns in data */
  float  ***vals;   /* ptr to matrix to be filled in with data */
{
  FILE    *fp;
  int     i, j;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* fills in numrows */
  fscanf(fp, "%d", nr);

  /* allocate data matrix */
  *vals = NR_matrix(1, *nr, 1, nc);
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (i = 1; i <= *nr; i++) {
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%g ", &(*vals)[i][j]);
    }
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_LC_MATRIX
Write a matrix of floats to a file in the standard data file format, with 
line count at top of file.
AG
*******************************************************************************/
int write_lc_matrix(outfile, nr, nc, mat)

  char   *outfile;  /* name of file to write to */
  int    nr, nc;    /* number of rows, columns of data */
  float  **mat;     /* matrix of data to be written */
{
  FILE    *fp;

  /* open file for writing */
  fp = fopen(outfile, "w");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print number of lines at top */
  fprintf(fp, "%d\n", nr);
  fclose(fp);

  /* open file for appending */
  fp = fopen(outfile, "a");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print matrix contents to file */
  print_matrix(fp, nr, nc, mat);
  fclose(fp);

  return (UT_OK);
}


/*******************************************************************************
READ_SET_OF_MATRICES
Read in ascii file containing a set of matrices, each separated by a carriage
return.
Allocates the vector of matrices of floats that will hold the data.
AG
*******************************************************************************/
int read_set_of_matrices(char *infile, int num_mats, int num_rows, 
                         int num_cols, float ****vals)
                                  
{
  FILE   *fp;
  int    i, j, m;
  float  **curr_mat;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = set_of_matrices(num_mats, num_rows, num_cols);
  if (*vals == (float***)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (m = 1; m <= num_mats; m++)
  {
    curr_mat = (*vals)[m];
    for (i = 1; i <= num_rows; i++)
      for (j = 1; j <= num_cols; j++)
        fscanf(fp, "%g ", (float*) &(curr_mat[i][j]));
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_SET_OF_DMATRICES
Read in ascii file containing a set of matrices, each separated by a carriage
return.
Allocates the vector of matrices of floats that will hold the data.
RG
*******************************************************************************/
int read_set_of_dmatrices(char *infile, int num_mats, int num_rows, 
                          int num_cols, double ****vals)
                                  
{
  FILE   *fp;
  int    i, j, m;
  double **curr_mat;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = set_of_dmatrices(num_mats, num_rows, num_cols);
  if (*vals == (double***)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (m = 1; m <= num_mats; m++)
  {
    curr_mat = (*vals)[m];
    for (i = 1; i <= num_rows; i++)
      for (j = 1; j <= num_cols; j++)
        fscanf(fp, "%lg ", (double*) &(curr_mat[i][j]));
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_SET_OF_SETS_OF_MATRICES
Read in ascii file containing a set of sets of matrices, where each of the 
sets of matrices is separated by two carriage returns, and each of the 
matrices in a set is separated by a single carriage return.

Allocates the vector of vectors of matrices of floats that will hold the data.
AG
*******************************************************************************/
int read_set_of_sets_of_matrices(char *infile, int num_sets, int num_mats, 
                                 int num_rows, int num_cols, float *****vals)
                                  
{
  FILE   *fp;
  int    i, j, m, s;
  float  ***curr_set, **curr_mat;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = set_of_sets_of_matrices(num_sets, num_mats, num_rows, num_cols);
  if (*vals == (float****)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (s = 1; s <= num_sets; s++)
  {
    curr_set = (*vals)[s];
    for (m = 1; m <= num_mats; m++)
    {
      curr_mat = curr_set[m];
      for (i = 1; i <= num_rows; i++)
        for (j = 1; j <= num_cols; j++)
          fscanf(fp, "%g", (float*) &(curr_mat[i][j]));
    }
  }

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
READ_SET_OF_SETS_OF_DMATRICES
Read in ascii file containing a set of sets of matrices, where each of the 
sets of matrices is separated by two carriage returns, and each of the 
matrices in a set is separated by a single carriage return.

Allocates the vector of vectors of matrices of floats that will hold the data.
RG
*******************************************************************************/
int read_set_of_sets_of_dmatrices(char *infile, int num_sets, int num_mats, 
                                 int num_rows, int num_cols, double *****vals)
                                  
{
  FILE   *fp;
  int    i, j, m, s;
  double  ***curr_set, **curr_mat;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = set_of_sets_of_dmatrices(num_sets, num_mats, num_rows, num_cols);
  if (*vals == (double****)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (s = 1; s <= num_sets; s++)
  {
    curr_set = (*vals)[s];
    for (m = 1; m <= num_mats; m++)
    {
      curr_mat = curr_set[m];
      for (i = 1; i <= num_rows; i++)
        for (j = 1; j <= num_cols; j++)
          fscanf(fp, "%lg", (double*) &(curr_mat[i][j]));
    }
  }

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
PRINT_SET_OF_MATRICES
Print the contents of a set-of-matrices structure, laying out the matrices
one after the other, separated by one carriage return.
AG
*******************************************************************************/
int print_set_of_matrices(FILE *stream, float ***mat_set, int num_mats, 
                          int num_rows, int num_cols)

{
  int i;

  for (i = 1; i <= num_mats; i++)
  {
    print_matrix(stream, num_rows, num_cols, mat_set[i]);
    if (i < num_mats)
      fprintf(stream, "\n");
  }
 
  return (UT_OK);
}


/*******************************************************************************
PRINT_SET_OF_DMATRICES
Print the contents of a set-of-matrices structure, laying out the matrices
one after the other, separated by one carriage return.
RG
*******************************************************************************/
int print_set_of_dmatrices(FILE *stream, double ***mat_set, int num_mats, 
                           int num_rows, int num_cols)

{
  int i;

  for (i = 1; i <= num_mats; i++)
  {
    print_dmatrix(stream, num_rows, num_cols, mat_set[i]);
    if (i < num_mats)
      fprintf(stream, "\n");
  }
 
  return (UT_OK);
}


/*******************************************************************************
PRINT_SET_OF_SETS_OF_MATRICES
Print the contents of a set-of-sets-of-matrices structure, laying out the 
matrices in the same set one after the other, separated by one carriage return,
and laying out the sets one after the other, separated by two carriage returns.
AG
*******************************************************************************/
int print_set_of_sets_of_matrices(FILE *stream, float ****mat_set_set, 
                                  int num_sets, int num_mats, int num_rows, 
                                  int num_cols)

{
  int i;

  for (i = 1; i <= num_sets; i++)
  {
    print_set_of_matrices(stream, mat_set_set[i], num_mats, num_rows, num_cols);
    fprintf(stream, "\n\n");
  }
 
  return (UT_OK);
}

/*******************************************************************************
PRINT_SET_OF_SETS_OF_DMATRICES
Print the contents of a set-of-sets-of-matrices structure, laying out the 
matrices in the same set one after the other, separated by one carriage return,
and laying out the sets one after the other, separated by two carriage returns.
RG
*******************************************************************************/
int print_set_of_sets_of_dmatrices(FILE *stream, double ****mat_set_set, 
                                   int num_sets, int num_mats, int num_rows, 
                                   int num_cols)

{
  int i;

  for (i = 1; i <= num_sets; i++)
  {
    print_set_of_dmatrices(stream, mat_set_set[i], num_mats, num_rows, num_cols);
    fprintf(stream, "\n\n");
  }
 
  return (UT_OK);
}


/*******************************************************************************
WRITE_SET_OF_MATRICES
Write the contents of a set-of-matrices structure to a file, laying out the 
matrices one after the other, separated by one carriage return.
AG
*******************************************************************************/
int write_set_of_matrices(char *out_file, float ***mat_set, int num_mats, 
                          int num_rows, int num_cols, char *mode)

{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(out_file, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print structure contents to file */
  print_set_of_matrices(fp, mat_set, num_mats, num_rows, num_cols);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_SET_OF_DMATRICES
Write the contents of a set-of-matrices structure to a file, laying out the 
matrices one after the other, separated by one carriage return.
RG
*******************************************************************************/
int write_set_of_dmatrices(char *out_file, double ***mat_set, int num_mats, 
                           int num_rows, int num_cols, char *mode)

{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(out_file, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print structure contents to file */
  print_set_of_dmatrices(fp, mat_set, num_mats, num_rows, num_cols);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_SET_OF_SETS_OF_MATRICES
Write the contents of a set-of-sets-of matrices structure to a file, laying 
out the  matrices one after the other, separated by one carriage return,
and laying out the sets one after the other, separated by two carriage returns.
AG
*******************************************************************************/
int write_set_of_sets_of_matrices(char *out_file, float ****mat_set_set, 
                                  int num_sets, int num_mats, int num_rows,
                                  int num_cols, char *mode)

{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(out_file, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print structure contents to file */
  print_set_of_sets_of_matrices(fp, mat_set_set, num_sets, num_mats, num_rows, 
                                num_cols);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
WRITE_SET_OF_SETS_OF_DMATRICES
Write the contents of a set-of-sets-of matrices structure to a file, laying 
out the  matrices one after the other, separated by one carriage return,
and laying out the sets one after the other, separated by two carriage returns.
RG
*******************************************************************************/
int write_set_of_sets_of_dmatrices(char *out_file, double ****mat_set_set, 
                                   int num_sets, int num_mats, int num_rows,
                                   int num_cols, char *mode)

{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(out_file, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* print structure contents to file */
  print_set_of_sets_of_dmatrices(fp, mat_set_set, num_sets, num_mats, num_rows, 
                                num_cols);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_GAUSS_PARMS
Read in a standard file containing mean vector and covariance matrix.
Allocates the structures that will hold the mean vector and covariance matrix.
AG
*******************************************************************************/
int read_gauss_parms(parmsfile, nc, mean, cov)

  char   *parmsfile;      /* name of file to read */
  int    nc;              /* number of dimensions of Gaussian */
  double  **mean, ***cov;  /* ptrs to mean vector and covariance matrix to be 
                             filled in */
{
  FILE    *fp;
  int     j, k;

  /* open file for reading */
  fp = fopen(parmsfile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate mean vector */
  *mean = NR_dvector(1, nc);
  if (*mean == (double*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* allocate covariance matrix */
  *cov = NR_dmatrix(1, nc, 1, nc);
  if (*cov == (double**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read mean vector from file */
  for (j = 1; j <= nc; j++) {
    fscanf(fp, "%lg ", &(*mean)[j]);
  }

  /* read covariance matrix from file */
  for (j = 1; j <= nc; j++) {
    for (k = 1; k <= nc; k++) {
      fscanf(fp, "%lg ", &(*cov)[j][k]);
    }
  }

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
PRINT_GAUSS_PARMS
Print the contents of a mean vector and of a covariance matrix.
AG
*******************************************************************************/
int print_gauss_parms(stream, nc, mean, cov)

  FILE   *stream;       /* stream to print to, e.g. stdout, file pointer */
  int    nc;            /* number of dimensions in Gaussian */
  double  *mean, **cov;  /* mean vector and covariance matrix to write */
{
  /* print mean vector */
  print_drow(stream, nc, mean);

  /* print covariance matrix */
  print_dmatrix(stream, nc, nc, cov);

  return (UT_OK);
}

/*******************************************************************************
PRINT_GAUSS_PARMS_SET
Print the contents of K mean vectors and covariance matrices.
AG
*******************************************************************************/
int print_gauss_parms_set(stream, nc, K, means, covars)

  FILE   *stream;       /* stream to print to, e.g. stdout, file pointer */
  int    nc;            /* number of dimensions in Gaussian */
  int    K;             /* number of Gaussians in the set */
  double  **means;       /* matrix containing K mean vectors to write */
  double  ***covars;     /* matrix containing K covariance matrices to write */
{
  int k;

  for (k = 1; k<=K; k++)
  {
    /* print identification number */
    fprintf(stream, "#%d\n\n",k);

    /* print mean vector to file */
    print_drow(stream, nc, means[k]);
    fprintf(stream, "\n");

    /* print covariance matrix to file */
    print_dmatrix(stream, nc, nc, covars[k]);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
PRINT_UNNORM_GAUSS_PARMS_SET
Print the contents of a set of K mean vectors and covariance matrices.
Takes data which has normalized by the range of values in each attribute,
and prints it in unnormalized form.
AG
*******************************************************************************/
int print_unnorm_gauss_parms_set(stream, nc, K, means, covars, range,
                                 minval)

  FILE   *stream;       /* stream to print to, e.g. stdout, file pointer */
  int    nc;            /* number of dimensions in Gaussian */
  int    K;             /* number of Gaussians in the set */
  double  **means;       /* matrix containing K mean vectors to write */
  double  ***covars;     /* matrix containing K covariance matrices to write */
  double  *range;        /* range of values for each attribute */
  double  *minval;       /* minimum value in each attribute */
{
  int k;

  for (k = 1; k<=K; k++)
  {
    /* print identification number */
    fprintf(stream, "#%d\n\n",k);

    /* print mean vector to file */
    print_unnorm_drow(stream, nc, means[k], range, minval);
    fprintf(stream, "\n");

    /* print covariance matrix to file */
    print_unnorm_cov_dmatrix(stream, nc, covars[k], range);
    fprintf(stream, "\n");
  }

  return (UT_OK);
}

/*******************************************************************************
WRITE_GAUSS_PARMS
Write out a standard file containing mean vector and covariance matrix.
AG
*******************************************************************************/
int write_gauss_parms(parmsfile, nc, mean, cov, mode)

  char   *parmsfile;    /* name of file to write to */
  int    nc;            /* number of dimensions in Gaussian */
  double  *mean, **cov;  /* mean vector and covariance matrix to write */
  char   *mode;         /* file i/o mode to use: e.g. "w" or "a" */
{
  FILE    *fp;

  /* open file in the specified mode */
  fp = fopen(parmsfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  print_gauss_parms(fp, nc, mean, cov);

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
WRITE_GAUSS_PARMS_SET
Write the contents of K mean vectors and covariance matrices to a standard
file.
AG
*******************************************************************************/
int write_gauss_parms_set(parmsfile, nc, K, means, covars, mode)

  char   *parmsfile;    /* name of file to write to */
  int    nc;            /* number of dimensions in Gaussian */
  int    K;             /* number of Gaussians in the set */
  double  **means;       /* matrix containing K mean vectors to write */
  double  ***covars;     /* matrix containing K covariance matrices to write */
  char   *mode;         /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE *fp;

  /* open file in the specified mode */
  fp = fopen(parmsfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  print_gauss_parms_set(fp, nc, K, means, covars);

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
WRITE_UNNORM_GAUSS_PARMS_SET
Write out a standard file containing a set of K mean vectors and covariance 
matrices.  Takes data which has normalized by the range of values in each 
attribute, and prints it in unnormalized form.
AG
*******************************************************************************/
int write_unnorm_gauss_parms_set(parmsfile, nc, K, means, covars, range,
                                 minval, mode)

  char   *parmsfile;    /* name of file to write to */
  int    nc;            /* number of dimensions in Gaussian */
  int    K;             /* number of Gaussians in the set */
  double  **means;       /* matrix containing K mean vectors to write */
  double  ***covars;     /* matrix containing K covariance matrices to write */
  double  *range;        /* range of values for each attribute */
  double  *minval;       /* minimum value in each attribute */
  char   *mode;         /* file i/o mode to use: e.g. "w", "r", or "a" */
{
  FILE *fp;

  /* open file in the specified mode */
  fp = fopen(parmsfile, mode);
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  print_unnorm_gauss_parms_set(fp, nc, K, means, covars, range, minval);

  fclose(fp);
  return (UT_OK);
}


/*******************************************************************************
READ_SUBSET2MATRIX
Read in part of an ascii data file into a matrix.
Allocates the matrix of floats that will hold the data.
Useful for data-parallel programs, where each processor reads in part of a
giant data file.  The starting row must be computed beforehand.
AG
*******************************************************************************/
int read_subset2matrix(infile, startrow, nr, nc, vals)
 
  char   *infile;   /* name of file to read */
  int    startrow;  /* first row of the file to start reading */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  int    skipsize;  /* number of elements to skip before starting read */
 
  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }
 
  /* allocate data matrix */
  *vals = NR_matrix(1, nr, 1, nc);
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
 
  /* set the file pointer */
  skipsize = (startrow-1)*nc;
  for (i = 1; i <= skipsize; i++)
    fscanf(fp, "%*g ");
 
  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    for (j = 1; j <= nc; j++) {
        fscanf(fp, "%g ", (float*) &(*vals)[i][j]);
    }
  }
 
  fclose(fp);
  return (UT_OK);
}
 
 
/*******************************************************************************
READ_SUBSET2DMATRIX
Read in part of an ascii data file into a matrix.
Allocates the matrix of floats that will hold the data.
Useful for data-parallel programs, where each processor reads in part of a
giant data file.  The starting row must be computed beforehand.
AG
*******************************************************************************/
int read_subset2dmatrix(infile, startrow, nr, nc, vals)
 
  char   *infile;   /* name of file to read */
  int    startrow;  /* first row of the file to start reading */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  double  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;
  double doubleval;
  int    skipsize;  /* number of elements to skip before starting read */
 
  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }
 
  /* allocate data matrix */
  *vals = NR_dmatrix(1, nr, 1, nc);
  if (*vals == (double**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
 
  /* set the file pointer */
  skipsize = (startrow-1)*nc;
  for (i = 1; i <= skipsize; i++)
    fscanf(fp, "%*g ");
 
  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%lg ", &doubleval);
      (*vals)[i][j] = doubleval;
    }
  }
 
  fclose(fp);
  return (UT_OK);
}
 
 
/*******************************************************************************
READ_BIN_SUBSET2MATRIX
Read in part of a binary data file (floats) into a matrix.
Allocates the matrix of floats that will hold the data.
Useful for data-parallel programs, where each processor reads in part of a
giant data file.  The starting row must be computed beforehand.
AG
*******************************************************************************/
int read_bin_subset2matrix(infile, startrow, nr, nc, vals)
 
  char   *infile;   /* name of file to read */
  int    startrow;  /* first row of the file to start reading */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;
 
  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }
 
  /* allocate data matrix */
  *vals = NR_matrix(1, nr, 1, nc);
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
 
  /* set the file pointer */
  fseek(fp,(startrow-1)*nc*sizeof(float),SEEK_SET);
 
  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    fread((float*) &(*vals)[i][1],sizeof(float),nc,fp);
  }
 
  fclose(fp);
  return (UT_OK);
}
 

/*******************************************************************************
READ_BIN_SUBSET2DMATRIX
Read in part of a binary data file (floats) into a matrix.
Allocates the matrix of floats that will hold the data.
Useful for data-parallel programs, where each processor reads in part of a
giant data file.  The starting row must be computed beforehand.
AG
*******************************************************************************/
int read_bin_subset2dmatrix(infile, startrow, nr, nc, vals)
 
  char   *infile;   /* name of file to read */
  int    startrow;  /* first row of the file to start reading */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  double  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i;
 
  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }
 
  /* allocate data matrix */
  *vals = NR_dmatrix(1, nr, 1, nc);
  if (*vals == (double**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }
 
  /* set the file pointer */
  fseek(fp,(startrow-1)*nc*sizeof(double),SEEK_SET);
 
  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    fread((double*) &(*vals)[i][1],sizeof(double),nc,fp);
  }
 
  fclose(fp);
  return (UT_OK);
}
 

/*******************************************************************************
READ_MATRIX_TRANSPOSE
Read in ascii file containing a matrix, into a form where the first (fast) 
index accesses the columns of the original data.  One way to think of it is
reading in the matrix 'sideways'.
Allocates the matrix of floats that will hold the data.
AG
*******************************************************************************/
int read_matrix_transpose(infile, nr, nc, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of rows in matrix */
  int    nc;        /* number of columns in matrix */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i, j;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix */
  *vals = NR_matrix(1, nc, 1, nr);  /* instead of NR_matrix(1, nr, 1, nc) */
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into the matrix */
  for (i = 1; i <= nr; i++) {
    for (j = 1; j <= nc; j++) {
      fscanf(fp, "%g ", (float*) &(*vals)[j][i]);
    }
  }

  fclose(fp);
  return (UT_OK);
}

/*******************************************************************************
WRITE_MATRIX_TRANSPOSE
Write out ascii file containing a matrix, into a transposed form, i.e the rows
become the columns and vice versa.  One way to think of it is writing out the 
matrix 'sideways'.
AG
*******************************************************************************/
int write_matrix_transpose(char *outfile, int nr, int nc, float **mat,
                           char *mode)

{
  FILE    *fp;
  int     i,j;

  /* open file in the specified mode */
  fp = fopen_return_if_fail(outfile, mode, UT_ERROR);

  /* print matrix contents to file */
  for (j = 1;j <= nc;j++)
  {    
    for (i = 1; i <= nr; i++) 
      fprintf(fp, "%g ", mat[i][j]);
    fprintf(fp, "\n");
  }

  fclose_return_if_fail(fp, UT_ERROR);
  return (UT_OK);
}


/*******************************************************************************
READ_LANDSAT
Read in file containing a LandSat image, which is in binary form (chars).
Note that each pixel is one byte (char), and it consists of nr x nc images, in
scan-line order, nb of them, one after the other.
It will be stored in a matrix such that each image becomes a single column
having all its nr x nc pixels strung out.
Allocates the matrix that will hold the data.
AG
*******************************************************************************/
int read_landsat(infile, nr, nc, nb, vals)

  char   *infile;   /* name of file to read */
  int    nr;        /* number of pixel rows in each image */
  int    nc;        /* number of pixel columns in each image */
  int    nb;        /* number of bands, or images */
  float  ***vals;   /* pointer to matrix to be filled in */
{
  FILE   *fp;
  int    i,j;
  unsigned char   *tempcol;

  /* open file for reading */
  fp = fopen(infile, "r");
  if (fp == (FILE*) NULL) {
    printf("Error in trying to open a file.\n");
    return (UT_ERROR);
  }

  /* allocate data matrix -
     each column corresponds to one of the images, one for each band;
     it will consist of the rows all strung together into one vertical string */
  *vals = NR_matrix(1, nr*nc, 1, nb);
  if (*vals == (float**)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* allocate vector that will hold each image as it is read in */
  tempcol = NR_cvector(1, nr*nc);
  if (tempcol == (unsigned char*)NULL) {
    printf("Error in memory allocation.\n");
    return (UT_ERROR);
  }

  /* read from file into temporary column, in binary format;
     then copy to matrix of floats */
  for (i = 1; i <= nb; i++) {
    /* more efficient to read in whole block... */
    fread( (unsigned char*) tempcol, sizeof(char), nr*nc, fp );  

    for (j = 1; j <= nr*nc; j++) {
      (*vals)[j][i] = (float) *(tempcol+j);  /* convert to float */
    }
  }

  /* clean up */
  NR_free_cvector(tempcol, 1, nr*nc);
  fclose(fp);

  return (UT_OK);
}

/*******************************************************************************
SKIP_HEADER
Move the file pointer forward in a file to skip the header
*******************************************************************************/
int skip_header(FILE *fp, char *file_type)
{
    char buf[64]; /* The number 64 is arbitrary */
    int  len;
    int  skip;
 
    fseek( fp, 0, SEEK_SET );  /* Make sure that the file pointer is at
                                   the beginning of the file */
 
    if (!strcmp( file_type, "none" )) /* The file has no header */
      skip = 0;
    else if (!strcmp( file_type, "vicar" )) { /* VICAR format */
 
      len = strlen( VIC_HEADER_LABEL );
      assert( sizeof( buf ) > len + 2 );
 
      fread( buf, 1, sizeof( buf ), fp );
 
      if( strncmp( buf, VIC_HEADER_LABEL, len ) != 0 ) /* not a VICAR file */
        return( 0 );
      else {
        buf[sizeof( buf ) - 1] = 0; /* Terminate the buffer string */
        skip = atoi( buf + len ); /* Parse the buffer string */
      }
    }
 
    fseek( fp, skip, SEEK_SET );
 
    return( UT_OK );
}

/*******************************************************************************
READ_VICAR
Read in an image in vicar format to a matrix of unsigned characters.
*******************************************************************************/
int read_vicar( char *infile, int nr, int nc, unsigned char ***vals )
{
    FILE   *fp;
    int    i;
   
    /* open file for reading */
    fp = fopen(infile, "r");
    if (fp == (FILE*) NULL) {
      printf("Error in trying to open a file.\n");
      return (UT_ERROR);
    }
   
    /* allocate data matrix */
    *vals = NR_cmatrix(1, nr, 1, nc);
    if (*vals == (unsigned char**)NULL) {
      printf("Error in memory allocation.\n");
      return (UT_ERROR);
    }
   
    /* skip the vicar file header */
    skip_header( fp, "vicar" );
   
    /* read in the data */
    for (i = 1; i <= nr; i++) {
      fread( (unsigned char *) &((*vals)[i][1]), sizeof(unsigned char), nc, fp );
    }
 
    fclose(fp);
 
    return( UT_OK );
}

