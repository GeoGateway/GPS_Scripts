/*******************************************************************************
MODULE NAME
ut_math

ONE-LINE SYNOPSIS
Utility functions related to mathematical computations.

SCOPE OF THIS MODULE
The Data Analysis (DA) library is meant to cover just about all mathematical 
code.  Thus, this module is intended only to contain the simplest mathematical
concepts, for programs for which it is not sensible to invoke the DA library.

SEE ALSO
The da_random module is one which might be considered to overlap significantly
with the concept of this module.  However, since the portable random number 
generation routines we have selected are in the Numerical Recipes library, the
functions built using them lie in the DA library, which is meant to provide
higher-level concepts one level above the NR library.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/hmm, AG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: ut_math.c,v 1.3 1997/05/06 22:32:03 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_math.c,v $
 * Revision 1.3  1997/05/06 22:32:03  agray
 * added header includes.
 *
 * Revision 1.2  1997/05/06 22:26:18  agray
 * added safe_log().
 *
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <string.h>
#include <math.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"
 
/* this module's header */
#include "ut_math.h"
 

/*******************************************************************************
SAFE_LOG
A wrapper around log() which handles arguments outside the domain of the log
function more gracefully.  Note that this functions takes and returns floats.
AG
*******************************************************************************/
float safe_log (float x)
{
  /* report negative arguments */
  if (x < 0.0)
  {
    err_printf();
    log_printf("Bad argument to log() function: %f\n",x);
    return(UT_ERROR);
  }

  /* handle everything between zero and the smallest representable float */
  if ( (x >= 0.0) && (x < UT_EPSILON_FLOAT) )
    return ( (float) log(UT_EPSILON_FLOAT) );

  /* finally, the safe case */
  else
    return ( (float) log(x) );
}

