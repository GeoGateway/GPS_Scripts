/*******************************************************************************
MODULE HEADER:
da_optim.h
*******************************************************************************/

#ifndef _DA_OPTIM_H_
#define _DA_OPTIM_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_optim_h_rcsid[] = "$Id: da_optim.h,v 1.1 1997/05/15 16:46:53 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_optim.h,v $
 * Revision 1.1  1997/05/15 16:46:53  granat
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
int simplex( float **p, float *y, int ndim, float ftol, AMOEBA_PF func,
             void *arg, int nargs, int *nfunc, float *best_p, float *scale );

#endif /* _DA_OPTIM_H_ */
