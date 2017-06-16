/*******************************************************************************
MODULE HEADER:
ut_output.h
*******************************************************************************/

#ifndef _UT_OUTPUT_H_
#define _UT_OUTPUT_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_output_h_rcsid[] = "$Id: ut_output.h,v 1.1 1997/01/29 23:44:54 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_output.h,v $
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

#define UT_SILENT_LOG  1
#define UT_NORMAL_LOG  2
#define UT_VERBOSE_LOG 3
#define UT_DEBUG_LOG   4

#define UT_NUM_LOG_LEVELS  4

#define UT_SILENT_LOG_NAME  "silent"
#define UT_NORMAL_LOG_NAME  "normal"
#define UT_VERBOSE_LOG_NAME "verbose"
#define UT_DEBUG_LOG_NAME   "debug"

#define UT_DEFAULT_LOG_FP  ((FILE*)stderr)

#define UT_MAX_ERROR_MSG_SIZE  256

/*******************************************************************************
SILENT_OUTPUT, NORMAL_OUTPUT, VERBOSE_OUTPUT, DEBUG_OUTPUT
True if the log level is at least 1, 2, 3, or 4, respectively.
Useful as a condition for an 'if' statement, say.
AG
*******************************************************************************/
#define silent_output()  (ut_log_level >= UT_SILENT_LOG)
#define normal_output()  (ut_log_level >= UT_NORMAL_LOG)
#define verbose_output() (ut_log_level >= UT_VERBOSE_LOG)
#define debug_output()   (ut_log_level >= UT_DEBUG_LOG)

/*******************************************************************************
ERR_PRINTF
Prints a standard error-reporting string, meant to be printed to the log file
right before calling log_printf() to print an error message.
AG
*******************************************************************************/
#define err_printf() fprintf(ut_log_fp, \
                             "--> Error, from line %d of file %s:\n    ", \
                             __LINE__, __FILE__);

/*******************************************************************************
PRINT_LOG_HEADER
Prints a standard header for a log file which includes the date and time.
AG
*******************************************************************************/
#define print_log_header() log_printf("Runtime log -- log level = %d\n", \
                                      ut_log_level); \
                           log_printf("Version compiled on %s at %s\n\n", \
                                      __DATE__, __TIME__);

/*******************************************************************************
PRINT_LOG_FOOTER
Prints a standard footer for a log file which includes the date and time.
AG
*******************************************************************************/
#define print_log_footer() log_printf("Done.\n");


/*==============================================================================
Variables
==============================================================================*/

extern FILE* ut_log_fp;
extern int   ut_log_level;
extern char  ut_err_msg[];

/*==============================================================================
Function Declarations
==============================================================================*/

int init_log_level_names(char **log_level_names);

void log_printf(char *format, ...);

#endif /* _UT_OUTPUT_H_ */




