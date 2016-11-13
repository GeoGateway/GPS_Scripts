/*******************************************************************************
MODULE NAME
ut_output

ONE-LINE SYNOPSIS
Utility functions related to logging the output of a program.

SCOPE OF THIS MODULE
All functions regarding printing in various formats or to various places should
go here.

SEE ALSO
Since error reporting falls under the purview of this module, ut_error might
also be consulted, as it concerns other aspects of error checking and handling.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/hmm, AG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: ut_output.c,v 1.1 1997/01/29 23:45:59 agray Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_output.c,v $
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <stdarg.h>

/* UT library */
#include "ut_error.h"

/* this module's header */
#include "ut_output.h"

/* global variables */
/* gcc doesn't seem to like this under redhat 7 */
/*
FILE* ut_log_fp = stdout;
*/
FILE* ut_log_fp;
int   ut_log_level;
char  ut_err_msg[UT_MAX_ERROR_MSG_SIZE];


/*******************************************************************************
INIT_LOG_LEVEL_NAMES
Initialize the name strings corresponding to the available log level options.
AG
*******************************************************************************/
int init_log_level_names ( char **log_level_names )

{
  log_level_names[UT_SILENT_LOG] = (char*) strdup (UT_SILENT_LOG_NAME);
  log_level_names[UT_NORMAL_LOG] = (char*) strdup (UT_NORMAL_LOG_NAME);
  log_level_names[UT_VERBOSE_LOG] = (char*) strdup (UT_VERBOSE_LOG_NAME);
  log_level_names[UT_DEBUG_LOG] = (char*) strdup (UT_DEBUG_LOG_NAME);

  return(UT_OK);
}


/*******************************************************************************
LOG_PRINTF
A wrapper around fprintf() which automatically uses the stream ut_log_fp.
AG
*******************************************************************************/
void log_printf(char *format, ...)

{
  va_list all_args;

  /* start at beginning of argument list */
  va_start(all_args, format);

  /* print message */
  vfprintf(ut_log_fp, format, all_args);

  /* finish looking at argument list */
  va_end(all_args);
}
