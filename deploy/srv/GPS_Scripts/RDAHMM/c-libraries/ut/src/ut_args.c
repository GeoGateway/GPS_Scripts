/*******************************************************************************
MODULE NAME
ut_args

ONE-LINE SYNOPSIS
Utility functions related to processing command-line program arguments.

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
static char rcsid[] = "$Id: ut_args.c,v 1.3 1997/03/06 20:07:11 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_args.c,v $
 * Revision 1.3  1997/03/06 20:07:11  agray
 * corrected malloc() bug in make_extended_name().
 *
 * Revision 1.2  1997/02/12 03:16:01  agray
 * fixed bug in match_option_name().
 *
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <string.h>

/* UT library */
#include "ut_output.h"
#include "ut_error.h"
#include "ut_string.h"
#include "ut_memory.h"

/* this module's header */
#include "ut_args.h"


/*******************************************************************************
GET_BASE_NAME
Return the subset of a given filename corresponding to the part of the filename
before the first period.  Returns the length of the basename.
AG
*******************************************************************************/
int get_base_name(char *file_name, char *base_name)

{
  int i, len;

  len = strlen(file_name);

  for (i = 0; i < len; i++)
    if (file_name[i] != '.')
      base_name[i] = file_name[i];
    else
      break;

  base_name[i] = (char)NULL;

  return( i );
}


/*******************************************************************************
MAKE_EXTENDED_NAME
Return a newly allocated string composed of a basename string followed by an
extension string.
AG
*******************************************************************************/
char *make_extended_name(char *base_name, char *ext_name)

{
  char *file_name;
  int  size;

  size = (strlen(base_name) + strlen(ext_name) + 1) * sizeof(char);
                                
  file_name = (char*) malloc_return_if_fail( size, (char*)NULL );

  strcpy(file_name, base_name);
  strcat(file_name, ext_name);

  return (file_name);
}


/*******************************************************************************
BUILD_OPTIONS_STRING
Build a string containing all the option names for a particular parameter,
given the array of option names.
AG
*******************************************************************************/
int build_options_string(char *string, char **option_names, int num_options,
                         char *option_descrip)

{
  int i;

  strcpy (string, option_descrip);
  strcat (string, " {");
  for (i = 1; i <= num_options; i++)
  {
    if (i>1)
      strcat (string, " | ");
    strcat (string, option_names[i]);
  }
  strcat (string, "}");

  return (UT_OK);
}


/*******************************************************************************
MATCH_OPTION_NAME
Given a string and the array of option names, return the index of name the
string matches, if any.
AG
*******************************************************************************/
int match_option_name(char *string, char **option_names, int num_options,
                      char *option_descrip)

{
  int i, index;

  index = -1;
  for (i = 1; i <= num_options; i++)
    if (streq(string, option_names[i]))
    {
      index = i;
      break; /* from for loop */
    }

  /* if an incorrect option was specified */
  if (index == -1)
  {
    err_printf();
    log_printf("Illegal %s specified: %s\n", option_descrip, 
               option_names[index]);

    /* show all legal options */
    log_printf ("Available %ss: ", option_descrip);
    for (i = 1; i <= num_options; i++)
      log_printf ("%s%s", ((i>1)?", ":""), option_names[i]);
    log_printf ("\n");
      
    return(UT_ERROR);
  }

  return( index );
}

