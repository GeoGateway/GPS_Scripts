/*******************************************************************************
MODULE HEADER:
da_signal.h
*******************************************************************************/

#ifndef _DA_SIGNAL_H_
#define _DA_SIGNAL_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_signal_h_rcsid[] = "$Id: da_signal.h,v 1.13 1998/07/02 01:14:26 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_signal.h,v $
 * Revision 1.13  1998/07/02 01:14:26  granat
 * added median filtering prototypes
 *
 * Revision 1.12  1998/06/29 22:10:52  granat
 * added functions to calculate energies
 *  revised functions to perform convolution via the FFT
 *
 * Revision 1.11  1997/06/20 22:11:42  granat
 * filled out conv_/corr functions, cosmetic changes
 *
 * Revision 1.10  1997/06/05 18:57:20  granat
 * added prototype for slow_full_conv_vector
 *
 * Revision 1.9  1997/01/30 01:53:19  agray
 * added semicolon.
 *
 * Revision 1.8  1997/01/29 21:58:38  agray
 * added range_of_cols(), changed range_normalize_cols().
 *
 * Revision 1.7  1997/01/29 21:39:11  agray
 * new format.
 *
 * Revision 1.6  1996/10/31 02:19:54  agray
 * renamed from "da_dist" to "da_signal";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.5  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.4  1996/09/13 01:15:22  agray
 * change name and comments for normalize_data().
 *
 * Revision 1.3  1996/09/13 00:57:42  agray
 * reduced arguments of print_unnorm_cov_matrix().
 *
 * Revision 1.2  1996/07/11 18:01:13  agray
 * added print_unnorm_matrix(), print_unnorm_cov_matrix(), print_unnorm_row(),
 * print_unnorm_col().
 *
 * Revision 1.1  1996/05/06 23:22:16  agray
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

/* signal energy and normalization */

float energy_vec(float *v, int n);
double energy_dvec(double *v, int n);

float normalize_energy_vec(float *v, int n);
double normalize_energy_dvec(double *v, int n);

/* FFTs */

int real_fft_2D(float **data, float *speq, int nn2, int nn3, int isign);

/* convolution and correlation */

int fft_conv_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                    int nc_b, float **c, int nr_c, int nc_c, float **data1,
		    float **data2, float *speq1, float *speq2);
int fft_conv_alloc_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                          int nc_b, float **c, int nr_c, int nc_c);

int slow_full_conv_vector( float *a, int n_a, float *b, int n_b, float *c,
                           int n_c );
int slow_crop_conv_vector( float *a, int n_a, float *b, int n_b, float *c,
                            int n_c );

int slow_full_conv_matrix( float **a, int nr_a, int nc_a, float **b, int nr_b,
                           int nc_b, float **c, int nr_c, int nc_c );
int slow_full_norm_conv_matrix( float **a, int nr_a, int nc_a, float **b,
                                int nr_b, int nc_b, float **c, int nr_c,
                                int nc_c );
int slow_crop_conv_matrix( float **a, int nr_a, int nc_a, float **b, int nr_b,
                           int nc_b, float **c, int nr_c, int nc_c );
int slow_full_corr_vector( float *a, int n_a, float *b, int n_b, float *c,
                           int n_c );
int slow_full_corr_matrix( float **a, int nr_a, int nc_a, float **b, int nr_b,
                           int nc_b, float **c, int nr_c, int nc_c );
int slow_crop_corr_matrix( float **a, int nr_a, int nc_a, float **b, int nr_b,
                           int nc_b, float **c, int nr_c, int nc_c );

int median3_filt_matrix(float **A, int nr, int nc, float **B);
int median3_filt_dmatrix(double **A, int nr, int nc, double **B);

#endif /* _DA_SIGNAL_H_ */
