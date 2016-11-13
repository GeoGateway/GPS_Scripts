/*******************************************************************************
MODULE NAME
da_signal_pp

ONE-LINE SYNOPSIS
Parallel general functions related to signal processing.

SCOPE OF THIS MODULE
Analogous to da_signal in scope, except that functions in this module are
written to run on a parallel machine.  

SEE ALSO
Most or all of the functions in this module are mirrors of their serial versions
in da_signal.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/qf, RG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_signal_pp.c,v 1.7 1998/04/21 17:38:05 roden Exp granat $";
#endif
/* 
 * $Log: da_signal_pp.c,v $
 * Revision 1.7  1998/04/21 17:38:05  roden
 * Took out #include <pvm3.h> because I think it's not needed here, and it
 * doesn't conform to da_platform.h control over presence/absence of pvm.
 *
 * Revision 1.6  1997/06/20 22:14:36  granat
 * changed to match changes in NR library
 *
 * Revision 1.5  1997/01/29 21:52:15  agray
 * new formatting, cleaning up debugging output using ut_output.
 *
 * Revision 1.4  1996/10/31 02:20:44  agray
 * renamed from "da_dist_pp" to "da_signal_pp";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * 
 * Revision 1.3  1996/09/23 22:45:54  agray
 * added pvm3.h inclusion.
 * 
 * Revision 1.2  1996/09/13 01:23:02  agray
 * change name and comments for normalize_data_pp(). added some headers.
 * 
 * Revision 1.1  1996/07/19 17:59:42  agray
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

/* NR library */
#include "nr.h"

/* DA library */
#include "da_comm_pp.h"

/* this module's header */
#include "da_signal_pp.h"


/*******************************************************************************
RANGE_NORMALIZE_COLS_PP
For each of k columns in a matrix, find the minimum and maximum values, and
the range (difference between the two).  Pass in a size k vector for each of 
these latter sets of values, and they will be filled in by this function.  
Modifies the original matrix such that each value in each column is normalized
ized by the range of that column.  Thus, all entries in the matrix will end up
scaled from 0-1.  Does this in parallel.
AG
*******************************************************************************/
int range_normalize_cols_pp(float **mat, float *minval, float *maxval, 
                            float *range, int numrows, int numcols, int pe, 
                            int numPE)

{
  int   i,j, p, cc;
  float *otherminval, *othermaxval;

  /* create storage for the passed min/max values */
  otherminval = NR_vector(1,numcols);
  othermaxval = NR_vector(1,numcols);

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

  if (debug_output() && (pe == 0))
    log_printf("%d: after getting local ranges, minval = %g, maxval = %g\n",
               pe,minval[1],maxval[1]);
  /*
    print_row(ut_log_fp, numcols, minval);
    print_row(ut_log_fp, numcols, maxval);
    */

  if (debug_output() && (pe == 0))
    log_printf("%d: before message passing of ranges\n",pe);

  /* Pass minimum and maximum values in from each process to the controlling
     process */
  DA_barrier();
  if (pe != 0) /* This is not the controlling process */
  {	
    /* Do the minimum value */
    cc = DA_send_msg(&minval[1],numcols,0,pe,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after sending local minvals\n",pe);

    /* Now block until receive broadcast message from the controlling
       process giving global minimum values */
    cc = DA_broadcast_msg(&minval[1],numcols,DA_ALL_PES,0,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after receiving global minvals, minval=%g, maxval=%g\n",
                 pe,minval[1],maxval[1]);
    /*
      log_printf("%d: minval\n",pe);
      print_row(ut_log_fp, numcols, minval);
      */

    /* Do the maximum value */
    cc = DA_send_msg(&maxval[1],numcols,0,pe,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after sending local maxvals\n",pe);
    
    /* Now block until receive broadcast message from the controlling
       process giving global maximum values */
    cc = DA_broadcast_msg(&maxval[1],numcols,DA_ALL_PES,0,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after receiving global maxvals, minval=%g, maxval=%g\n",
                 pe,minval[1],maxval[1]);
    /*
      log_printf("%d: maxval\n",pe);
      print_row(ut_log_fp, numcols, maxval);
      */

  }
  else /* This is the controlling process */
  {
    /* Now receive the minimum values from each PE */
    /*      time1 = time_now(); */
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&otherminval[1],numcols,p,p,DA_FLOAT);

      if (debug_output())
        log_printf("%d: after receiving local minvals\n",pe);

      for(j=1; j<=numcols; j++)
        if (otherminval[j] < minval[j])
          minval[j] = otherminval[j];
    } 
    /* Now broadcast the global minimum values */
    cc = DA_broadcast_msg(&minval[1],numcols,DA_ALL_PES,0,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after broadcasting global minvals\n",pe);

    /* Now receive the maximum values from each PE */
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&othermaxval[1],numcols,p,p,DA_FLOAT);

      if (debug_output())
        log_printf("%d: after receiving local maxvals\n",pe);
 
      for(j=1; j<=numcols; j++)
        if (othermaxval[j] > maxval[j])
          maxval[j] = othermaxval[j];
    }
    /* Now broadcast the global minimum values */
    cc = DA_broadcast_msg(&maxval[1],numcols,DA_ALL_PES,0,DA_FLOAT);

    if (debug_output())
      log_printf("%d: after broadcasting global maxvals\n",pe);
    /* time2 = time_now(); */
    /* log_printf("%d: time to pass min/max values = %lf\n",pe,time2-time1); */
  }

  if (debug_output() && (pe == 0))
    log_printf("%d: after message passing of ranges\n",pe);

  if (debug_output())
    log_printf("%d: after message passing of ranges, minval=%g, maxval=%g\n",
               pe,minval[1],maxval[1]);
  /*
    log_printf("%d: minval, maxval\n",pe);
    print_row(ut_log_fp, numcols, minval);
    print_row(ut_log_fp, numcols, maxval);
    */

/*----------------------------------------------------------------------------*/
  /* finish normalizing */

  /* scale each attribute value by its range and translate by its minval */
  for (j=1; j<=numcols; j++)
    range[j] = maxval[j] - minval[j];

  for (i=1; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
      mat[i][j] = (mat[i][j] - minval[j]) / range[j];
  
  /* free storage related to message passing */
  NR_free_vector(othermaxval, 1, numcols);
  NR_free_vector(otherminval, 1, numcols);

  if (debug_output() && (pe == 0))
    log_printf("%d: after finish normalizing\n",pe);

  return (UT_OK);
}


/*******************************************************************************
RANGE_NORMALIZE_DCOLS_PP
For each of k columns in a matrix, find the minimum and maximum values, and
the range (difference between the two).  Pass in a size k vector for each of 
these latter sets of values, and they will be filled in by this function.  
Modifies the original matrix such that each value in each column is normalized
ized by the range of that column.  Thus, all entries in the matrix will end up
scaled from 0-1.  Does this in parallel.
AG
*******************************************************************************/
int range_normalize_dcols_pp(double **mat, double *minval, double *maxval, 
                            double *range, int numrows, int numcols, int pe, 
                            int numPE)

{
  int   i,j, p, cc;
  double *otherminval, *othermaxval;

  /* create storage for the passed min/max values */
  otherminval = NR_dvector(1,numcols);
  othermaxval = NR_dvector(1,numcols);

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

  if (debug_output() && (pe == 0))
    log_printf("%d: after getting local ranges, minval = %g, maxval = %g\n",
               pe,minval[1],maxval[1]);
  /*
    print_row(ut_log_fp, numcols, minval);
    print_row(ut_log_fp, numcols, maxval);
    */

  if (debug_output() && (pe == 0))
    log_printf("%d: before message passing of ranges\n",pe);

  /* Pass minimum and maximum values in from each process to the controlling
     process */
  DA_barrier();
  if (pe != 0) /* This is not the controlling process */
  {	
    /* Do the minimum value */
    cc = DA_send_msg(&minval[1],numcols,0,pe,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after sending local minvals\n",pe);

    /* Now block until receive broadcast message from the controlling
       process giving global minimum values */
    cc = DA_broadcast_msg(&minval[1],numcols,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after receiving global minvals, minval=%g, maxval=%g\n",
                 pe,minval[1],maxval[1]);
    /*
      log_printf("%d: minval\n",pe);
      print_row(ut_log_fp, numcols, minval);
      */

    /* Do the maximum value */
    cc = DA_send_msg(&maxval[1],numcols,0,pe,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after sending local maxvals\n",pe);
    
    /* Now block until receive broadcast message from the controlling
       process giving global maximum values */
    cc = DA_broadcast_msg(&maxval[1],numcols,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after receiving global maxvals, minval=%g, maxval=%g\n",
                 pe,minval[1],maxval[1]);
    /*
      log_printf("%d: maxval\n",pe);
      print_row(ut_log_fp, numcols, maxval);
      */

  }
  else /* This is the controlling process */
  {
    /* Now receive the minimum values from each PE */
    /*      time1 = time_now(); */
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&otherminval[1],numcols,p,p,DA_DOUBLE);

      if (debug_output())
        log_printf("%d: after receiving local minvals\n",pe);

      for(j=1; j<=numcols; j++)
        if (otherminval[j] < minval[j])
          minval[j] = otherminval[j];
    } 
    /* Now broadcast the global minimum values */
    cc = DA_broadcast_msg(&minval[1],numcols,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after broadcasting global minvals\n",pe);

    /* Now receive the maximum values from each PE */
    for(p=1;p<numPE;p++)
    {
      cc = DA_recv_msg(&othermaxval[1],numcols,p,p,DA_DOUBLE);

      if (debug_output())
        log_printf("%d: after receiving local maxvals\n",pe);
 
      for(j=1; j<=numcols; j++)
        if (othermaxval[j] > maxval[j])
          maxval[j] = othermaxval[j];
    }
    /* Now broadcast the global minimum values */
    cc = DA_broadcast_msg(&maxval[1],numcols,DA_ALL_PES,0,DA_DOUBLE);

    if (debug_output())
      log_printf("%d: after broadcasting global maxvals\n",pe);
    /* time2 = time_now(); */
    /* log_printf("%d: time to pass min/max values = %lf\n",pe,time2-time1); */
  }

  if (debug_output() && (pe == 0))
    log_printf("%d: after message passing of ranges\n",pe);

  if (debug_output())
    log_printf("%d: after message passing of ranges, minval=%g, maxval=%g\n",
               pe,minval[1],maxval[1]);
  /*
    log_printf("%d: minval, maxval\n",pe);
    print_row(ut_log_fp, numcols, minval);
    print_row(ut_log_fp, numcols, maxval);
    */

/*----------------------------------------------------------------------------*/
  /* finish normalizing */

  /* scale each attribute value by its range and translate by its minval */
  for (j=1; j<=numcols; j++)
    range[j] = maxval[j] - minval[j];

  for (i=1; i<=numrows; i++)
    for (j=1; j<=numcols; j++)
      mat[i][j] = (mat[i][j] - minval[j]) / range[j];
  
  /* free storage related to message passing */
  NR_free_dvector(othermaxval, 1, numcols);
  NR_free_dvector(otherminval, 1, numcols);

  if (debug_output() && (pe == 0))
    log_printf("%d: after finish normalizing\n",pe);

  return (UT_OK);
}
