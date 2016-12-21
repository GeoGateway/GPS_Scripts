/*******************************************************************************
MODULE HEADER:
ut_types.h
*******************************************************************************/

#ifndef _UT_TYPES_H_
#define _UT_TYPES_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_types_h_rcsid[] = "$Id: ut_types.h,v 1.3 1997/11/20 01:06:12 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_types.h,v $
 * Revision 1.3  1997/11/20 01:06:12  agray
 * removing bool, since it conflicts with something in C++/STL
 *
 * Revision 1.2  1997/11/19 23:36:26  agray
 * put #ifndef wrappers around typedefs.
 *
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures, Types
==============================================================================*/

#ifndef real
typedef float real;
#endif

#ifndef boolean
typedef unsigned char boolean;
#endif

/*==============================================================================
Constants, Macros
==============================================================================*/

#ifndef NULL
#define NULL 0
#endif

#ifndef UT_FALSE
#define UT_FALSE 0
#define UT_TRUE  1
#endif

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

#endif /* _UT_TYPES_H_ */

