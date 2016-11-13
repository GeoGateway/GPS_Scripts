/*******************************************************************************
MODULE NAME
ut_error

ONE-LINE SYNOPSIS
Utility functions related to checking for and handling error conditions.

SCOPE OF THIS MODULE
As stated.

SEE ALSO
The reporting of error conditions falls in the realm of ut_output.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/hmm, AG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: ut_error.c,v 1.1 1997/01/29 23:45:59 agray Exp agray $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_error.c,v $
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <string.h>

/* UT library */
#include "ut_output.h"

/* this module's header */
#include "ut_error.h"

/* global variables */
int   ut_status;
char *ut_status_ptr;


