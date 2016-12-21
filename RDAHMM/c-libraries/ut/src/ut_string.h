/*******************************************************************************
MODULE HEADER:
ut_string.h
*******************************************************************************/

#ifndef _UT_STRING_H_
#define _UT_STRING_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_string_h_rcsid[] = "$Id: ut_string.h,v 1.2 1999/06/03 02:23:20 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_string.h,v $
 * Revision 1.2  1999/06/03 02:23:20  agray
 * added strcaseeq.
 *
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

#include "ut_platform.h"  /* to ensure this is accounted for */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

#define UT_MAX_STRING  256

/*******************************************************************************
STREQ
True if the two string arugments are equal to each other.
Useful as a condition for an 'if' statement, say.
Case-sensitive and -insensitive versions.
AG
*******************************************************************************/
#define streq(a, b) (strcmp((a), (b)) == 0)
#define strcaseeq(a, b) (strcasecmp((a), (b)) == 0)

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

void remove_newlines (char *string);

/* machine dependent declarations */

#if (UT_STRDUP_EXISTS == UT_DEFINE_STRDUP)
char *strdup (char *s);
#endif

#if (UT_STRCASECMP_EXISTS == UT_DEFINE_STRCASECMP)
int strcasecmp (char *s1, char *s2);
#endif


#endif /* _UT_STRING_H_ */
