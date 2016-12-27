/*******************************************************************************
MODULE HEADER:
ut_args.h
*******************************************************************************/

#ifndef _UT_ARGS_H_
#define _UT_ARGS_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_args_h_rcsid[] = "$Id: ut_args.h,v 1.1 1997/01/29 23:44:54 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_args.h,v $
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

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* filenames */
int get_base_name(char *file_name, char *base_name);
char *make_extended_name(char *base_name, char *ext_name);

/* option names */
int build_options_string(char *string, char **option_names, int num_options,
                         char *option_descrip);
int match_option_name(char *string, char **option_names, int num_options,
                      char *option_descrip);

#endif /* _UT_ARGS_H_ */




