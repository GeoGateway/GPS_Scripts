/*******************************************************************************
MODULE NAME
ut_memory

ONE-LINE SYNOPSIS
Utility functions related to allocating and de-allocating memory.

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
static char rcsid[] = "$Id: ut_memory.c,v 1.3 1997/07/29 03:26:08 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_memory.c,v $
 * Revision 1.3  1997/07/29 03:26:08  agray
 * moved out all matrix/vector stuff to da_memory.
 *
 * Revision 1.2  1997/06/02 15:32:56  granat
 * changed to use new NR name convention
 *
 * Revision 1.1  1997/01/29 23:45:59  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdio.h>
#include <stdlib.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"

/* this module's header */
#include "ut_memory.h"

/* global variables */
int ut_mem_used = 0;
int ut_mem_allocated = 0;


/*******************************************************************************
MALLOC_AND_TRACK
Allocate a chunk of memory.  Checks for invalid size argument, and keeps track
of the amount of memory in use and the total amount allocated (this total
includes memory which might have subsequently been freed).  

Uses calloc() in debug mode, which also clears the allocated memory of garbage 
values.  Otherwise does the quicker malloc().
AG
*******************************************************************************/
void* malloc_and_track(int size)

{
  /* check size argument for validity first */
  if (size <= 0)
  {
    err_printf();
    log_printf("Tried to allocate %s bytes.\n",size);
    return ( (void*) NULL );
  }
  else
  {
    ut_mem_used += size;
    ut_mem_allocated += size;
    if (debug_output())
      return ( (void*) calloc(1,size) );
    else
      return ( (void*) malloc(size) );
  }
}


/*******************************************************************************
FREE_AND_TRACK
De-allocate a chunk of memory.  Checks for invalid pointer argument, and keeps 
track of the amount of memory in use and the total amount allocated (this total
includes memory which might have subsequently been freed).  

Since the free() function of the C library has no method of returning errors,
this function currently simply ignores the call if a NULL pointer was passed.
AG
*******************************************************************************/
void free_and_track(void *mem, int size)

{
  /* check pointer argument for validity first */
  if (mem == (void*) NULL)
  {
    err_printf();
    log_printf("Tried to free a NULL pointer.\n");
    return;
  }
  else
  {
    ut_mem_used -= size;
    free(mem);
  }
}

