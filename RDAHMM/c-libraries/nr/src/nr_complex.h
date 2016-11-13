/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file complex.h.  Do not confuse this file with the same-named
   file complex.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_


#ifndef lint
static char nr_complex_hdr_rcsid[] = "$Id: nr_complex.h,v 1.2 1997/07/29 15:48:10 granat Exp $";
#endif
/* 
 * $Log: nr_complex.h,v $
 * Revision 1.2  1997/07/29 15:48:10  granat
 * edited to use new "NR_" paradigm
 *
 * Revision 1.1  1996/07/25 17:26:54  agray
 * Initial revision
 *
 * */

#ifndef _NR_FCOMPLEX_DECLARE_T_
typedef struct NR_FCOMPLEX {float r,i;} nr_fcomplex;
#define _NR_FCOMPLEX_DECLARE_T_
#endif /* _NR_FCOMPLEX_DECLARE_T_ */

nr_fcomplex NR_Cadd(nr_fcomplex a, nr_fcomplex b);
nr_fcomplex NR_Csub(nr_fcomplex a, nr_fcomplex b);
nr_fcomplex NR_Cmul(nr_fcomplex a, nr_fcomplex b);
nr_fcomplex NR_Complex(float re, float im);
nr_fcomplex NR_Conjg(nr_fcomplex z);
nr_fcomplex NR_Cdiv(nr_fcomplex a, nr_fcomplex b);
float NR_Cabs(nr_fcomplex z);
nr_fcomplex NR_Csqrt(nr_fcomplex z);
nr_fcomplex NR_RCmul(float x, nr_fcomplex a);

#endif /* _NR_COMPLEX_H_ */
