/*******************************************************************************
MODULE HEADER:
ut_error.h
*******************************************************************************/

#ifndef _UT_ERROR_H_
#define _UT_ERROR_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_error_h_rcsid[] = "$Id: ut_error.h,v 1.1 1997/01/29 23:44:54 agray Exp agray $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_error.h,v $
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

#define UT_OK     0
#define UT_ERROR  1

/*******************************************************************************
REPORT_IF_BAD
Evaluates a given status integer.  If the status is not UT_OK, print a given
error message to the log file.
AG
*******************************************************************************/
#define report_if_bad(s,m)  if (s!=UT_OK) { err_printf(); log_printf(m); }

/*******************************************************************************
RETURN_IF_BAD
Evaluates a given status integer.  If the status is not UT_OK, print a given
error message to the log file and return from the calling function.  (Only a
macro can do this.)
AG
*******************************************************************************/
#define return_if_bad(s,m,r)  if (s!=UT_OK) { err_printf(); log_printf(m); \
                                              return(r); }

/*******************************************************************************
EXIT_IF_BAD
Evaluates a given status integer.  If the status is not UT_OK, print a given
error message to the log file and exit from the calling program.  Before 
exiting, indicates in the log file that we are exiting and closes the log file.
AG
*******************************************************************************/
#define exit_if_bad(s,m)    if (s!=UT_OK) { err_printf(); log_printf(m); \
                                            log_printf("Exiting.\n"); \
                                            fclose(ut_log_fp); \
                                            exit(UT_ERROR); }

/*******************************************************************************
REPORT_IF_NULL, RETURN_IF_NULL, EXIT_IF_NULL
Same as report_if_bad(), except that it evaluates a status pointer (i.e. a 
pointer cast to char*) rather than an integer.
AG
*******************************************************************************/
#define report_if_null(s,m) if ((void*)s==(void*)NULL) \
                                { err_printf(); log_printf(m); }

#define return_if_null(s,m,r) if ((void*)s==(void*)NULL) \
                                  { err_printf(); log_printf(m); \
                                    return(r); }

#define exit_if_null(s,m)  if ((void*)s==(void*)NULL) \
                               { err_printf(); log_printf(m); \
                               log_printf("Exiting.\n"); \
                               fclose(ut_log_fp); \
                               exit(UT_ERROR); }

/*******************************************************************************
REPORT_IF_NULL, RETURN_IF_NULL, EXIT_IF_NULL
Same as report_if_bad(), except that it evaluates a status pointer pointer
(i.e. a pointer cast to char**) rather than an integer.
AG
*******************************************************************************/
#define report_if_null(s,m) if ((void*)s==(void*)NULL) \
                                    { err_printf(); log_printf(m); }

#define return_if_null(s,m,r) if ((void*)s==(void*)NULL) \
                                      { err_printf(); log_printf(m); \
                                        return(r); }

#define exit_if_null(s,m)  if ((void*)s==(void*)NULL) \
                                   { err_printf(); log_printf(m); \
                                     log_printf("Exiting.\n"); \
                                     fclose(ut_log_fp); \
                                     exit(UT_ERROR); }

/*==============================================================================
Variables
==============================================================================*/

extern int   ut_status;
extern char *ut_status_ptr;

/*==============================================================================
Function Declarations
==============================================================================*/

#endif /* _UT_ERROR_H_ */

