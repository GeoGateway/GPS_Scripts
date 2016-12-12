/*******************************************************************************
MODULE HEADER:
ut_memory.h
*******************************************************************************/

#ifndef _UT_MEMORY_H_
#define _UT_MEMORY_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_memory_h_rcsid[] = "$Id: ut_memory.h,v 1.2 1997/07/29 03:26:43 agray Exp agray $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_memory.h,v $
 * Revision 1.2  1997/07/29 03:26:43  agray
 * moved out all matrix/vector stuff to da_memory.
 *
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

/*******************************************************************************
MALLOC_REPORT_IF_FAIL
Attempt to allocate a chunk of memory.  If the attempt fails, report an error.

Note that the caller can pretend that this macro returns the memory pointer
needed, i.e. one can say "array = (int*) malloc_report_if_fail(s)".
AG
*******************************************************************************/
#define malloc_report_if_fail(s) \
  ut_status_ptr = (char*) malloc_and_track(s); \
  sprintf(ut_err_msg, "Couldn't allocate memory chunk of size %d\n", s); \
  report_if_null(ut_status_ptr, ut_err_msg);

/*******************************************************************************
MALLOC_RETURN_IF_FAIL
Attempt to allocate a chunk of memory.  If the attempt fails, report an error
and return from the calling function.  (Only a macro can do this.)

Note that the caller can pretend that this macro returns the memory pointer
needed, i.e. one can say "array = (int*) malloc_return_if_fail(s)".
AG
*******************************************************************************/
#define malloc_return_if_fail(s,r) \
  ut_status_ptr = (char*) malloc_and_track(s); \
  sprintf(ut_err_msg, "Couldn't allocate memory chunk of size %d\n", s); \
  return_if_null(ut_status_ptr, ut_err_msg, r);

/*******************************************************************************
MALLOC_EXIT_IF_FAIL
Attempt to allocate a chunk of memory.  If the attempt fails, report an error
and exit from the calling program.  Closes the log file first.

Note that the caller can pretend that this macro returns the memory pointer
needed, i.e. one can say "array = (int*) malloc_exit_if_fail(s)".
AG
*******************************************************************************/
#define malloc_exit_if_fail(s) \
  ut_status_ptr = (char*) malloc_and_track(s); \
  sprintf(ut_err_msg, "Couldn't allocate memory chunk of size %d\n", s); \
  exit_if_null(ut_status_ptr, ut_err_msg);

/*==============================================================================
Variables
==============================================================================*/

extern int ut_mem_used;
extern int ut_mem_allocated;

/*==============================================================================
Function Declarations
==============================================================================*/

/* general memory allocation and de-allocation */

void* malloc_and_track(int size);
void free_and_track(void *mem, int size);

#endif /* _UT_MEMORY_H_ */
