/*******************************************************************************
MODULE HEADER:
da_probstat.h
*******************************************************************************/

#ifndef _DA_PROBSTAT_H_
#define _DA_PROBSTAT_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_probstat_h_rcsid[] = "$Id: da_probstat.h,v 1.12 1997/09/10 14:55:38 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_probstat.h,v $
 * Revision 1.12  1997/09/10 14:55:38  granat
 * added prototypes for routines to extend functionality to matrices of ints,
 * chars
 *
 * Revision 1.11  1997/05/13 23:57:37  granat
 * added prototypes for mean_mat, stdev_mat, var_vec, var_mat
 * changed  name of mean to mean_vec, name of stdev to stdev_vec
 * changed order of input variables on some functions
 *
 * Revision 1.10  1997/05/06 22:49:38  granat
 * fixed interpolation prototypes
 *
 * Revision 1.9  1997/05/06 22:22:28  agray
 * added some things from dp cooltool.
 *
 * Revision 1.8  1997/02/24 22:59:57  granat
 * added prototypes for linint, vector_linint, vector_ratint, vector_polint
 *
 * Revision 1.7  1997/01/29 21:33:16  agray
 * new format.
 * also changed chop_data_equally() to partition_vector().
 *
 * Revision 1.6  1996/10/31 02:18:27  agray
 * renamed from "da_prob" to "da_probstat";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.5  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.4  1996/09/20 22:13:43  agray
 * changed names back.
 *
 * Revision 1.3  1996/09/19 22:45:08  agray
 * changed name of prob_gauss() to gauss_eval(); similar for vector_prob_gauss().
 *
 * Revision 1.2  1996/07/11 18:27:02  agray
 * added read/print/write_gauss_parms(), print/write_gauss_parms_set(),
 * print/write_unnorm_gauss_parms_set(), prob_gauss(), vector_prob_gauss().
 *
 * Revision 1.1  1996/02/21 05:18:29  agray
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

/* basic statistics */

float mean_vec(float *v, int nc);
float mean_mat(float **m, int nr, int nc);
float mean_imat(int **m, int nr, int nc);
float mean_cmat(unsigned char **m, int nr, int nc);
float stdev_vec(float *v, int nc, float mean);
float stdev_mat(float **m, int nr, int nc, float mean);
float var_vec(float *v, int nc, float mean);
float var_mat(float **m, int nr, int nc, float mean);
int est_covar_mat(float **covar, float **data, float *mean, int num_data,
                  int num_dims);
int stats_of_cols(float **mat, float *means, float *stdevs, int numrows, 
                  int numcols);

/* pdf's of distributions */

double prob_gauss(double *datum, double *mean, double **covar, int dim, 
                 double min_diag);
double* vector_prob_gauss(double **data, double *mean, double **covar, int dim, 
                         int numdata, double *probs, double min_diag,
                         char *rows_or_cols);

double prob_mixture(double *datum, double **means, double ***covars, double *weights,
                   int dim, int num_comps, double min_diag);
double prob_mixture_and_keep(double *datum, double **means, double ***covars, 
                            double *weights, int dim, int num_comps, 
                            double min_diag, double *wgtd_probs);

/* random variables */

int generate_disc_rv_value(double *prob_symbol, int num_symbols, 
                           double *cum_prob_symbol);
int generate_mixture_rv_value(double *value, int num_dims, int num_comps, 
                              double **means, double ***covars, double *weights,
                              double *cum_prob_comps, double **chol_factor,
                              double *diag);
int generate_gaussian_rv_value(double *value, int num_dims, double *mean, 
                               double **covar, double **chol_factor, double *diag);

/* computations using distributions */

int merf(double *mean, double **cov, double **data, int nr, int nc, double *probs);

/* interpolation */

int linint( float *xa, float *ya, int n, float x, float *y, float slop );
int vector_linint( float *xa, float *ya, int n, float *x, float *y, int m,
                   float slop );
int vector_ratint(float *xa, float *ya, int n, float *x, float *y, int m);
int vector_polint(float *xa, float *ya, int n, float *x, float *y, int m);

/* polynomials */

int fit_line_to_data(float *data, float *index, int num_data, float *slope, 
                     float *intercept);
int compute_line(float *data, float *index, int num_data, float slope, 
                 float intercept);

/* histograms */

int partition_vector(float *data, float *start_values, int num_data,
                   int num_intervals);
int count_in_intervals(float *data, int *count, float *start_value,
                       int num_data, int num_buckets);


#endif /* _DA_PROBSTAT_H_ */
