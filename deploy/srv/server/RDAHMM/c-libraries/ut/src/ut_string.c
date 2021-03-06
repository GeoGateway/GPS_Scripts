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
static char rcsid[] = "$Id: ut_string.c,v 1.1 1997/01/29 23:45:59 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_string.c,v $
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <stdarg.h>

/* UT library */
#include "ut_error.h"
#include "ut_platform.h"
#include "ut_memory.h"
#include "ut_output.h"

/* this module's header */
#include "ut_string.h"


/*******************************************************************************
REMOVE_NEWLINES
Replace all instances of the newline character in the given string with the 
space character.
AG, based on JR
*******************************************************************************/
void remove_newlines (char *string)

{
  int   num_chars;
  int   i;
     
  num_chars = strlen(string);

  for (i=0; i < num_chars; i++)
    if (string[i] == '\n')
      string[i] = ' ';
}


/*******************************************************************************
STRDUP
This function determines the length of the given string. Then it creates an
allocated string buffer of that size, followed by the copying of the given
string into the new buffer.  Then it returns the newly created string.
AG, based on JT
*******************************************************************************/
/* issue 2 */
#if (UT_STRDUP_EXISTS == UT_DEFINE_STRDUP)
char *strdup (char *s)

{
  char	*p=NULL;

  /* +1 for '\0' */
  if ((p = (char *) malloc_and_track((unsigned)strlen(s) + 1)) != (char *)NULL)
    strcpy(p,s);

  return((char *)p);
}
#endif


/*******************************************************************************
STRCASECMP
This function compares the two given input strings and returns 0 if they are
equal and non-zero if not. This compare function is like strsmp, but not case
sensitive.  (Note: this uses just the strcmp for now).
AG, based on JT
*******************************************************************************/
/* issue 3 */
#if (UT_STRCASECMP_EXISTS == UT_DEFINE_STRCASECMP)
int strcasecmp (char *s1, char *s2)
	
{
  return ( (int)strcmp(s1,s2) );
}
#endif


