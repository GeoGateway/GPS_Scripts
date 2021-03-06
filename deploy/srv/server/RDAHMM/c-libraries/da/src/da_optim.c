/*******************************************************************************
MODULE NAME
da_optim

ONE-LINE SYNOPSIS
General functions related to function optimization.

SCOPE OF THIS MODULE
Any functions that perform function optimization should go here, unless
they are more directly related to some other module (ie, clustering goes
in da_cluster).

SEE ALSO
For the most part there, should little overlap with other modules, although
da_cluster is an exception.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/qf, RG.

NOTES
-

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_optim.c,v 1.4 1999/07/20 16:34:59 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_optim.c,v $
 * Revision 1.4  1999/07/20 16:34:59  granat
 * added braces to remove ambiguous else
 *
 * Revision 1.3  1997/08/24 00:56:32  granat
 * conformed to changes in da_nrhacks
 *
 * Revision 1.2  1997/06/02 15:53:27  granat
 * changed to use new NR naming convention
 *
 * Revision 1.1  1997/05/15 16:46:28  granat
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
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_random.h"
#include "da_nrhacks.h"

/* this module's header */
#include "da_optim.h"

int simplex( float **p, float *y, int ndim, float ftol, AMOEBA_PF func,
             void *arg, int nargs, int *nfunc, float *best_p, float *scale )
{
  int   i, j;
  int   best;
  int   success;
  float best_sum;
 
  best = 1;
  best_sum = 2 * y[best];
 
  set_rand_by_time_of_day();
 
  while ((y[best]/best_sum < 0.99) && (y[best] >= 0)) {
 
    /* store values from last iteration */
 
    best_sum = y[best];
    for (i = 1; i <= ndim; i++)
      best_p[i] = p[best][i];
 
    /* do the simplex optimization */

    success = DA_amoeba( p, y, ndim, ftol, func, arg, nargs, nfunc );
 
    /* Handle major error */
 
    if (success == -1) {
      printf( "Warning: Function value went under threshold -- breaking out of \
               optimization.\n" );
      success = 0;
      break;
    }
 
    best = 1;
    for (i = 1; i <= (ndim+1); i++) {
      if (y[best] > y[i])
        best = i;
    }
 
    /* initialize the seed values for restart, adding random noise */
 
    for (i = 1; i <= (ndim+1); i++) {
      for (j = 1; j <= ndim; j++)
        if (i < best) {
          if (i == j)
            p[i][j] = p[best][j] + scale[j] * (gen_rand() + 0.5);
          else
            p[i][j] = p[best][j];
	} 
	else if (i > best) {
          if ((i - 1) == j)
            p[i][j] = p[best][j] + scale[j] * (gen_rand() + 0.5);
          else
            p[i][j] = p[best][j];
	}
 
      y[i] = func( arg, nargs, p[i], ndim );
    }
  }
 
  if (success)
    return( UT_OK );
  else
    return( UT_ERROR );
}
