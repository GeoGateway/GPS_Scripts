/*******************************************************************************
MODULE NAME
da_io_pp

ONE-LINE SYNOPSIS
Parallel general functions related to data input and output.

SCOPE OF THIS MODULE
Analogous to da_io in scope, except that functions in this module are
written to run on a parallel machine.  

SEE ALSO
Many of the functions in this module will be mirrors of their serial versions
in da_io.  

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/sc_pp, AG.
2. /proj/cooltools/kmeans_pp, RG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_io_pp.c,v 1.13 1999/07/22 16:52:37 granat Exp $";
#endif
/* 
 * $Log: da_io_pp.c,v $
 * Revision 1.13  1999/07/22 16:52:37  granat
 * added double and unsigned char versions of write_col_cascade()
 *
 * Revision 1.12  1999/07/14 22:12:31  granat
 * added function to count total lines
 *
 * Revision 1.11  1998/06/29 22:12:05  granat
 * fixed bug in paralle control flow
 *
 * Revision 1.10  1998/05/15 18:16:11  granat
 * fixed bug in write_col_channel_cascade_pp and variants
 *
 * Revision 1.9  1998/04/21 17:37:20  roden
 * Took out #include <pvm3.h> because I think it's not needed here, and it
 * doesn't conform to da_platform.h control over presence/absence of pvm.
 *
 * Revision 1.8  1997/09/10 14:52:37  granat
 * revised many routines to reflect new experience
 *
 * Revision 1.7  1997/01/29 21:45:09  agray
 * added new includes.
 *
 * Revision 1.6  1996/10/31 02:15:43  agray
 * renamed from "da_data_pp" to "da_io_pp";
 * changed .h and .c formats throughout library;
 *
 * Revision 1.5  1996/09/23 22:45:42  agray
 * added pvm3.h inclusion.
 *
 * Revision 1.4  1996/09/23 22:42:14  agray
 * minor changes.
 *
 * Revision 1.3  1996/07/19 17:55:19  agray
 * added write_col_cascade_pp(), write_icol_cascade_pp(),
 * write_col_channel_cascade_pp(), write_icol_channel_cascade_pp()
 *
 * Revision 1.2  1996/07/16 00:52:47  agray
 * correction to comment.
 *
 * Revision 1.1  1996/07/16 00:43:28  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"
#include "ut_types.h"

/* DA library */
#include "da_io.h"
#include "da_comm_pp.h"

/* this module's header */
#include "da_io_pp.h"


/*******************************************************************************
SPLIT_DATA_PP
Computes which portion of a dataset to read in, depending on which processor 
you are, for the purpose of splitting a large dataset over several processors.
my_numrows is the number of data (rows) this processor gets; startrow is the
index of the row at which this processor should start reading.
AG
*******************************************************************************/
int split_data_pp(int *startrow, int *my_numrows, int total_numrows, int pe, 
                  int numPE)

{
  int extrarows;

  extrarows = total_numrows%numPE;
  if (extrarows > pe)
  {
    *my_numrows = (int)(total_numrows/numPE)+1;
    *startrow = pe*(*my_numrows) + 1;
  }
  else
  {
    *my_numrows = (int)(total_numrows/numPE);
    *startrow = pe*(*my_numrows) + extrarows + 1;
  }

  return (UT_OK);
}


/*******************************************************************************
TOTAL_DATA_PP
Get the total number of rows stored across all processors.
RG
*******************************************************************************/
int total_data_pp(int numrows, int *total_numrows, int pe, int numPE)
{
  int p;
  int cc;
  int other_numrows;
  
  if (pe != 0) {
    /* Send out the local rows */
    cc = DA_send_msg(&numrows, 1, 0, pe, DA_INT);
    
    /* Now block until receive broadcast message from the controlling
       process giving the total number of rows */
    cc = DA_broadcast_msg(total_numrows, 1, DA_ALL_PES, 0, DA_INT);
  }
  else {
    *total_numrows = numrows;

    /* Receive the other numbers of rows */
    for (p = 1; p < numPE; p++) {
      cc = DA_recv_msg(&other_numrows, 1, p, p, DA_INT);
      *total_numrows += other_numrows;
    }

    /* Now broadcast the total number of rows */
    cc = DA_broadcast_msg(total_numrows, 1, DA_ALL_PES, 0, DA_INT);
  }
  
  return (UT_OK);
}


/*******************************************************************************
READ_MAT_CHANNEL_CASCADE
Parallel version of read_bin_matrix(), which reads part of a large matrix so
that the whole of the matrix is distributed over several processors.
RG
*******************************************************************************/
int read_mat_channel_cascade( char *datafile, int numrows, int numcols, 
                              float ***data, int startrow, int channels,
                              int pe, int numPE )
{
  int i;
  
  for (i = 0; i < numPE; i += channels) {
    if ((pe >= i) && (pe < (i + channels)))
      read_bin_subset2matrix(datafile, startrow, numrows, numcols, data);
    DA_barrier();
  }

  return( UT_OK );
}


/*******************************************************************************
WRITE_COL_CASCADE_PP
Parallel version of write_col(), which assumes that a large vector is stored
on several processors, in the order dictated by their processor numbers;
the processors each write their portion of the vector in turn to the same
destination file.
RG
*******************************************************************************/
int write_col_cascade_pp(outfile, dim, v, mode, pe, numPE)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row vector to be printed */
  float  *v;        /* vector to be printed */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;        /* processor number */
{
  int     cc;
  int     dummy;

  /* check to see that the previous processors have already finished */

  DA_barrier();

  /* if this isn't the first processes, check to see if the previous
     process is finished writing */

  if (pe > 0)
    cc = DA_recv_msg( &dummy, 1, pe-1, numPE, DA_INT);

  if (pe == 0)
    write_col(outfile, dim, v, mode);
  else
    write_col(outfile, dim, v, "a");

  /* if this isn't the last process, send a message to the next process
     telling it to start writing */

  if (pe < (numPE - 1)) /* Not the last process */
    cc = DA_send_msg( &dummy, 1, pe+1, numPE, DA_INT);

  return (UT_OK);
}


/*******************************************************************************
WRITE_ICOL_CASCADE_PP
Similar to write_col_cascade_pp().
RG
*******************************************************************************/
int write_icol_cascade_pp(outfile, dim, v, mode, pe, numPE)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row vector to be printed */
  int    *v;        /* vector to be printed */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;        /* processor number */
{
  int     cc;
  int     dummy;

  /* check to see that the previous processors have already finished */

  DA_barrier();

  /* if this isn't the first processes, check to see if the previous
     process is finished writing */

  if (pe > 0)
    cc = DA_recv_msg( &dummy, 1, pe-1, numPE, DA_INT);

  if (pe == 0)
    write_icol(outfile, dim, v, mode);
  else {
    write_icol(outfile, dim, v, "a");
  }

  /* if this isn't the last process, send a message to the next process
     telling it to start writing */

  if (pe < (numPE - 1)) /* Not the last process */
    cc = DA_send_msg( &dummy, 1, pe+1, numPE, DA_INT);

  return (UT_OK);
}


/*******************************************************************************
WRITE_DCOL_CASCADE_PP
Similar to write_col_cascade_pp().
RG
*******************************************************************************/
int write_dcol_cascade_pp(outfile, dim, v, mode, pe, numPE)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row vector to be printed */
  double *v;        /* vector to be printed */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;        /* processor number */
{
  int     cc;
  int     dummy;

  /* check to see that the previous processors have already finished */

  DA_barrier();

  /* if this isn't the first processes, check to see if the previous
     process is finished writing */

  if (pe > 0)
    cc = DA_recv_msg( &dummy, 1, pe-1, numPE, DA_INT);

  if (pe == 0)
    write_dcol(outfile, dim, v, mode);
  else {
    write_dcol(outfile, dim, v, "a");
  }

  /* if this isn't the last process, send a message to the next process
     telling it to start writing */

  if (pe < (numPE - 1)) /* Not the last process */
    cc = DA_send_msg( &dummy, 1, pe+1, numPE, DA_INT);

  return (UT_OK);
}


/*******************************************************************************
WRITE_CCOL_CASCADE_PP
Similar to write_col_cascade_pp().
RG
*******************************************************************************/
int write_ccol_cascade_pp(outfile, dim, v, mode, pe, numPE)

  char   *outfile;  /* name of file to write to */
  int    dim;       /* length of the row vector to be printed */
  unsigned char *v; /* vector to be printed */
  char   *mode;     /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;        /* processor number */
{
  int     cc;
  int     dummy;

  /* check to see that the previous processors have already finished */

  DA_barrier();

  /* if this isn't the first processes, check to see if the previous
     process is finished writing */

  if (pe > 0)
    cc = DA_recv_msg( &dummy, 1, pe-1, numPE, DA_INT);

  if (pe == 0)
    write_ccol(outfile, dim, v, mode);
  else {
    write_ccol(outfile, dim, v, "a");
  }

  /* if this isn't the last process, send a message to the next process
     telling it to start writing */

  if (pe < (numPE - 1)) /* Not the last process */
    cc = DA_send_msg( &dummy, 1, pe+1, numPE, DA_INT);

  return (UT_OK);
}


/*******************************************************************************
WRITE_COL_CHANNEL_CASCADE_PP
Parallel version of write_col(), which assumes that a large vector is stored
on several processors, in the order dictated by their processor numbers.
Similar to write_col_cascade_pp(), except that this function will only use
up to the number of i/o channels specified, to avoid thrashing the i/o. 
It is generally faster than write_col_cascade_pp().
AG
******************************************************************************/
int write_col_channel_cascade_pp(outfile, dim, v, num_channels, mode, pe, numPE)

  char   *outfile;      /* name of file to write to */
  int    dim;           /* length of the row vector to be printed */
  float  *v;            /* vector to be printed */
  int    num_channels;  /* number of i/o channels to occupy */
  char   *mode;         /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;            /* processor number */
  int    numPE;         /* total number of processors */
{
  int     i;

  DA_barrier();
  for(i = 0; i < numPE; i += num_channels)
  {
    /* only print in groups of num_channels processors */
    if ((pe >= i) && (pe < (i + num_channels)))
    {
      /* print our portion of the column */
      if (pe == 0)
        write_col(outfile, dim, v, mode);
      else
        write_col(outfile, dim, v, "a");

      DA_barrier();
    }
  }

  return (UT_OK);
}


/*******************************************************************************
WRITE_ICOL_CHANNEL_CASCADE_PP
Similar to write_icol_channel_cascade_pp().
AG
*******************************************************************************/
int write_icol_channel_cascade_pp(outfile, dim, v, num_channels, mode, pe,numPE)

  char   *outfile;      /* name of file to write to */
  int    dim;           /* length of the row vector to be printed */
  int    *v;            /* vector to be printed */
  int    num_channels;  /* number of i/o channels to occupy */
  char   *mode;         /* file i/o mode to use: e.g. "w", "r", or "a" */
  int    pe;            /* processor number */
  int    numPE;         /* total number of processors */
{
  int     i;

  DA_barrier();
  for(i=0; i<numPE; i+=num_channels)
  {
    /* only print in groups of num_channels processors */
    if ((pe >= i) && (pe < (i + num_channels)))
    {
      /* print our portion of the column */
      if (pe==0)
        write_icol(outfile, dim, v, mode);
      else
        write_icol(outfile, dim, v, "a");

      DA_barrier();
    }
  }

  return (UT_OK);
}
