/*******************************************************************************
MODULE HEADER:
da_linalg.h
*******************************************************************************/

#ifndef _DA_LINALG_H_
#define _DA_LINALG_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_linalg_h_rcsid[] = "$Id: da_linalg.h,v 1.29 1998/07/02 01:15:48 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_linalg.h,v $
 * Revision 1.29  1998/07/02 01:15:48  granat
 * added prototypes for fast matrix multiply
 *
 * Revision 1.28  1998/06/29 22:11:37  granat
 * added prototype for LSQR method of Paige and Saunders
 *
 * Revision 1.27  1998/05/07 23:56:59  granat
 * removed prototypes for functions moved to da_util
 *
 * Revision 1.26  1998/05/01 17:42:52  granat
 * edited prototypes to match major changes in da_linalg.c
 * added comments and grouped functions
 *
 * Revision 1.25  1997/10/21 14:33:16  granat
 * added prototypes to match changes in da_linalg.c
 *
 * Revision 1.24  1997/09/10 14:50:50  granat
 * added fast_set_mat macros
 *
 * Revision 1.23  1997/09/04 20:22:54  granat
 * added protypes for variants of sum_mat
 * added prototypes for fast_copy_mat_section and variants
 *
 * Revision 1.22  1997/08/11 18:32:12  granat
 * added prototype for set_imat
 *
 * Revision 1.21  1997/06/20 22:10:24  granat
 * changed prototype for transpose function
 *
 * Revision 1.20  1997/06/05 18:54:37  granat
 * added flip_vector, changed prototypes to match with conventions, some cosmetic changes
 *
 * Revision 1.19  1997/05/06 23:05:49  granat
 * fixed bugs in memory setting functions, added fast_zero_mat and other
 * similar macros
 *
 * Revision 1.18  1997/05/06 22:23:51  agray
 * added some things from dp cooltool.
 *
 * Revision 1.17  1997/05/06 21:19:29  granat
 * added prototype for copy_mat_section()
 *
 * Revision 1.16  1997/04/05 19:10:43  granat
 * added prototypes for sum_mat() and norm_sum_mat()
 * adjust prototypes of functions so that they all follow input parameter conventions
 *
 * Revision 1.15  1997/03/27 18:10:54  granat
 * added prototypes for transpose_in_situ_alloc_matrix(),
 * transpose_in_situ_sqr_matrix(), flip_left_right_matrix(),
 * flip_top_bottom_matrix()
 *
 * Revision 1.14  1997/03/15 17:51:32  granat
 * Added prototypes for flip_left_right_matrix() and flip_top_bottom_matrix()
 *
 * Revision 1.13  1997/02/14 00:03:37  granat
 * changed prototype of transpose to transpose_matrix
 * added prototype for transpose_imatrix
 *
 * Revision 1.12  1997/01/29 21:31:04  agray
 * new format.
 * also added copy_sub_vec_to_vec().
 *
 * Revision 1.11  1996/10/31 02:16:35  agray
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.10  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.9  1996/07/15 18:51:02  agray
 * updated k&r prototypes.
 *
 * Revision 1.8  1996/07/13 01:25:31  agray
 * moved out print/write_row/col/irow/icol() to da_data module
 *
 * Revision 1.7  1996/07/11 18:12:16  agray
 * moved out read_gauss_parms(), write_gauss_parms() to da_prob module;
 * added add_mat(), subtract_mat(), invert_mat_copy(), det_copy(),
 * restrict_illcond_matrix(), scalar_mult/div/add/subtract_mat/vec(), set_mat(),
 * set_vec(), mult/div_vec_elt(), sum_vec(), max/min_vec(), arg_max/min_vec(),
 * copy_vec().
 *
 * Revision 1.6  1996/04/09 02:50:11  agray
 * added print_irow(), print_icol(), write_irow(), write_icol().
 * ag
 *
 * Revision 1.5  1996/04/09  02:44:04  agray
 * should have gotten checked in with da_linalg.c's changes,
 * long ago.
 * ag
 *
 * Revision 1.4  1996/02/29  02:33:55  agray
 * moved write_bin_matrix() and read_bin_matrix() to da_data module
 * ag
 *
 * Revision 1.3  1996/02/21 03:51:43  agray
 * added pca().
 * ag
 *
 * Revision 1.2  1996/02/21  00:37:41  agray
 * added write_bin_matrix() and read_bin_matrix()
 * ag
 *
 * Revision 1.1  1996/02/06  03:28:23  agray
 * Initial revision
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

/* norms of vectors and matrices */

float norm_vec(float *v, int n);
double norm_dvec(double *v, int n);

float frobenius_norm_mat(float **m, int nr, int nc);
double frobenius_norm_dmat(double **m, int nr, int nc);

/* transpose of matrices */

int transpose_matrix(float **a, int nr, int nc, float **a_trans);
int transpose_dmatrix(double **a, int nr, int nc, double **a_trans);
int transpose_imatrix(int **a, int nr, int nc, int **a_trans);
int transpose_cmatrix(unsigned char **a, int nr, int nc, unsigned char **a_trans);

int transpose_in_situ_alloc_matrix(float ***a, int nr, int nc, float *temp_vect,
                                   char *mem_choice);
int transpose_in_situ_alloc_dmatrix(double ***a, int nr, int nc, 
		                    double *temp_vect, char *mem_choice);
int transpose_in_situ_alloc_imatrix(int ***a, int nr, int nc, int *temp_vect, 
                                    char *mem_choice);
int transpose_in_situ_alloc_cmatrix(unsigned char ***a, int nr, int nc, 
		                    unsigned char *temp_vect, char *mem_choice);

int transpose_in_situ_sqr_matrix(float **a, int n);
int transpose_in_situ_sqr_dmatrix(double **a, int n);
int transpose_in_situ_sqr_imatrix(int **a, int n);
int transpose_in_situ_sqr_cmatrix(unsigned char **a, int n);

/* inversion and determinants of matrices */

int invert_alloc_mat(float **mat, int n, float **inv);
int invert_alloc_dmat(double **mat, int n, double **inv);
int invert_alloc_sym_mat(float **mat, int n, float **inv);
int invert_alloc_sym_dmat(double **mat, int n, double **inv);
int invert_copy_alloc_mat(float **mat, int n, float **inv);
int invert_copy_alloc_sym_dmat(double **mat, int n, double **inv);
int jacobi_invert_copy_alloc_dmat(double **mat, int n, double **inv);
int jacobi_invert_copy_alloc_sym_dmat(double **mat, int n, double **inv);
float det_alloc_mat(float **mat, int n);
double det_alloc_dmat(double **mat, int n);
float det_copy_alloc_mat(float **mat, int nc);
double det_copy_alloc_dmat(double **mat, int nc);

/* conditioning of matrices */

int restrict_illcond_matrix(float **mat, int dim, float min_diag);
int restrict_illcond_dmatrix(double **mat, int dim, double min_diag);

/* distance metrics */

float euclid_dist_vec(float *v1, float *v2, int n);
double euclid_dist_dvec(double *v1, double *v2, int n);

/* basic vector and matrix operations */

int right_mult_matrix(float **m, int nr, int nc, float *v, float *result);
int right_mult_dmatrix(double **m, int nr, int nc, double *v, double *result);

int left_mult_matrix(float **m, int nr, int nc, float *v, float *result);
int left_mult_dmatrix(double **m, int nr, int nc, double *v, double *result);

int mat_mult(float **a, int nra, int nca, float **b, int nrb, int ncb, 
             float **c);
int dmat_mult(double **a, int nra, int nca, double **b, int nrb, int ncb, 
              double **c);
int fast_mat_mult(float **A, int nra, int nca, float **B, int nrb, int ncb, 
                  float **C);
int fast_dmat_mult(double **A, int nra, int nca, double **B, int nrb, int ncb, 
                   double **C);

int add_mat(float **m1, float **m2, float **m3, int nr, int nc);
int add_dmat(double **m1, double **m2, double **m3, int nr, int nc);

int subtract_mat(float **m1, float **m2, float **m3, int nr, int nc);
int subtract_dmat(double **m1, double **m2, double **m3, int nr, int nc);

int add_vec(float *v1, float *v2, float *v3, int n);
int add_dvec(double *v1, double *v2, double *v3, int n);

int subtract_vec(float *v1, float *v2, float *v3, int n);
int subtract_dvec(double *v1, double *v2, double *v3, int n);

float dot_product_vec(float *x, float *y, int n);
double dot_product_dvec(double *x, double *y, int n);

int outer_product_vec(float *x, float *y, int n, float **prod);
int outer_product_dvec(double *x, double *y, int n, double **prod);

/* scalar operations on vectors and matrices */

int scalar_mult_vec(float *v, int n, float constant);
int scalar_mult_dvec(double *v, int n, double constant);

int scalar_div_vec(float *v, int n, float constant);
int scalar_div_dvec(double *v, int n, double constant);

int scalar_add_vec(float *v, int n, float constant);
int scalar_add_dvec(double *v, int n, double constant);

int scalar_subtract_vec(float *v, int n, float constant);
int scalar_subtract_dvec(double *v, int n, double constant);

int scalar_mult_mat(float **m, int nr, int nc, float constant);
int scalar_mult_dmat(double **m, int nr, int nc, double constant);

int scalar_div_mat(float **m, int nr, int nc, float constant);
int scalar_div_dmat(double **m, int nr, int nc, double constant);

int scalar_add_mat(float **m, int nr, int nc, float constant);
int scalar_add_dmat(double **m, int nr, int nc, double constant);

int scalar_subtract_mat(float **m, int nr, int nc, float constant);
int scalar_subtract_dmat(double **m, int nr, int nc, double constant);

/* element-wise operations on vectors and matrices */

int mult_vec_elt(float *v1, float *v2, float *v3, int n);
int mult_dvec_elt(double *v1, double *v2, double *v3, int n);

int div_vec_elt(float *v1, float *v2, float *v3, int n);
int div_dvec_elt(double *v1, double *v2, double *v3, int n);

int mult_mat_elt(float **m1, float **m2, float **m3, int nr, int nc);
int mult_dmat_elt(double **m1, double **m2, double **m3, int nr, int nc);

int div_mat_elt(float **m1, float **m2, float **m3, int nr, int nc);
int div_dmat_elt(double **m1, double **m2, double **m3, int nr, int nc);

/* other operations on matrices and vectors */

int pca_mat(float **A, int nr, int nc, float *w, float **VT);
int pca_alloc_mat(float **A, int nr, int nc, float *w, float **VT);

/* solving systems of linear equations */

int lsqr_mat(float **A, int nr, int nc, float *b, float *x, float *u,
             float *v, float *w, float atol, float btol, float conlim);

int lsqr_dmat(double **A, int nr, int nc, double *b, double *x, double *u,
              double *v, double *w, double atol, double btol, double conlim);

/* testing matrices */

int posdef_sym_alloc_mat(float **mat, int dim, int *posdef);
int posdef_sym_alloc_dmat(double **mat, int dim, int *posdef);

#endif /* _DA_LINALG_H_ */
