/*******************************************************************************
MODULE HEADER:
da_signal_pp.h
*******************************************************************************/

#ifndef _DA_SIGNAL_PP_H_
#define _DA_SIGNAL_PP_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_signal_pp_hdr_rcsid[] = "$Id: da_signal_pp.h,v 1.4 1997/01/29 21:39:34 agray Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_signal_pp.h,v $
 * Revision 1.4  1997/01/29 21:39:34  agray
 * new format.
 *
 * Revision 1.3  1996/10/31 02:20:44  agray
 * renamed from "da_dist_pp" to "da_signal_pp";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 *
 * Revision 1.2  1996/09/13 01:23:49  agray
 * changed name for normalize_data_pp().
 *
 * Revision 1.1  1996/07/19 18:00:02  agray
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

/* normalization of values by their range */

int range_normalize_cols_pp(float **mat, float *minval, float *maxval, 
                            float *range, int numrows, int numcols, int pe, 
                            int numPE);

int range_normalize_dcols_pp(double **mat, double *minval, double *maxval, 
                            double *range, int numrows, int numcols, int pe, 
                            int numPE);

#endif /* _DA_SIGNAL_PP_H_ */
