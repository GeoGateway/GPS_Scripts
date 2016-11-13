/*******************************************************************************
MODULE HEADER:
ut_math.h
*******************************************************************************/

#ifndef _UT_MATH_H_
#define _UT_MATH_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_math_h_rcsid[] = "$Id: ut_math.h,v 1.3 1997/05/06 22:26:35 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_math.h,v $
 * Revision 1.3  1997/05/06 22:26:35  agray
 * added safe_log().
 *
 * Revision 1.2  1997/01/29 23:47:46  agray
 * improved UT_PI slightly.
 *
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

#include <float.h>
/* for FLT_MAX, FLT_MIN, FLT_EPSILON */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

#define UT_PI  3.14159265358979
#define UT_E   2.718281828

#define UT_MAX_FLOAT      FLT_MAX      /* from float.h */
#define UT_MIN_FLOAT      FLT_MIN      /* from float.h */
#define UT_EPSILON_FLOAT  FLT_EPSILON  /* from float.h */


/*******************************************************************************
SGN
An expression that evaluates to 1 if the argument is positive and -1 otherwise.
AG
*******************************************************************************/
#define sgn(a) ((a) >= 0.0 ? 1 : -1)

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* elementary unary operations */

float safe_log (float x);


#endif /* _UT_MATH_H_ */
