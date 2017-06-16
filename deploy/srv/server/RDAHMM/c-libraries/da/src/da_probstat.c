/*******************************************************************************
MODULE NAME
da_probstat

ONE-LINE SYNOPSIS
General functions related to probability and statistics.

SCOPE OF THIS MODULE
Any functions relating specifically to other probabilistic or statistical 
concepts which have representative modules in this library should go in the 
appropriate module.  Functions that apply more generally are intended to go 
here.  For instance, HMM induction, while entirely a probabilistic concept,
is handled in da_timeseg.

SEE ALSO
Because the definition of this module is quite broad, there is some potential
overlap with several other modules in this library.

REFERENCE(S)
-

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_probstat.c,v 1.27 2001/01/03 00:04:05 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_probstat.c,v $
 * Revision 1.27  2001/01/03 00:04:05  granat
 * changed vector_prob_gauss to use cholesky decomposition
 *
 * Revision 1.26  2000/11/22 00:58:33  granat
 * changed prob_gauss so that it uses cholesky decomposition
 * also eliminated extraneous decomposition for approximate doubling in
 * speed
 *
 * Revision 1.25  2000/11/06 19:22:37  granat
 * changed gaussian probability routines so that they called symmetric matrix
 * inverse functions
 *
 * Revision 1.24  2000/04/04 23:39:12  granat
 * changed prob_gauss so that it returns DBL_EPSILON instead of zero
 *
 * Revision 1.23  1998/06/30 17:24:12  agray
 * changed includes.
 *
 * Revision 1.22  1998/05/07 23:45:06  granat
 * Made edits to conform with function name changes in da_linalg and da_util
 *
 * Revision 1.21  1998/04/21 21:48:51  roden
 * Changed calls to conform to new naming conventions as follows:
 *   mat_vect_mult() -> right_mult_matrix()
 *   dot_product()   -> dot_product_vec()
 * Previous version, 1.20, was Granat's mod/fix to var_mat(); I stole lock
 * to make above changes.
 *
 * Revision 1.20  1998/04/21 21:08:11  roden
 * *** empty log message ***
 *
 * Revision 1.19  1998/03/09 00:44:25  granat
 * fixed 2nd matrix stdev bug
 *
 * Revision 1.18  1998/03/09 00:37:14  granat
 * fixed matrix stdev
 * fixed boolean bug
 * ,.
 *
 * Revision 1.17  1997/09/10 14:54:25  granat
 * added routines to extend functionality to matrices of ints, uchars
 *
 * Revision 1.16  1997/06/20 22:13:36  granat
 * changed to reflect bug fix in da_random
 *
 * Revision 1.15  1997/06/05 18:56:17  granat
 * edited to conform to changes in da_linalg
 *
 * Revision 1.14  1997/06/02 15:54:35  granat
 * changed to use new NR naming convention
 *
 * Revision 1.13  1997/05/13 23:56:30  granat
 * added mean_mat, stdev_mat, var_vec, var_mat
 * changed name of mean to mean_vec, name of stdev to stdev_vec
 * changed order of input variables on some functions
 *
 * Revision 1.12  1997/05/06 22:22:43  agray
 * added some things from dp cooltool.
 *
 * Revision 1.11  1997/04/05 19:01:51  granat
 * changed a couple of things to adjust for changes to da_linalg
 *
 * Revision 1.10  1997/02/24 22:59:13  granat
 * added linint, vector_linint, vector_ratint, vector_polint functions
 * for function interpolation
 *
 * Revision 1.9  1997/01/29 21:47:16  agray
 * new formatting, cleaning up debugging output using ut_output,
 * cleaning up random number generation using da_random.  also
 * changed chop_data_equally() to partition_vector().
 *
 * Revision 1.8  1996/10/31 02:18:27  agray
 * renamed from "da_prob" to "da_probstat";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.7  1996/09/20 22:19:03  agray
 * corrected comments.
 *
 * Revision 1.6  1996/09/20 22:13:31  agray
 * changed names back.
 *
 * Revision 1.5  1996/09/19 22:44:40  agray
 * changed name of prob_gauss() to gauss_eval(); similar for 
 * vector_prob_gauss().
 *
 * Revision 1.4  1996/09/16 21:04:40  agray
 * updated call to print_unnorm_cov_matrix().
 *
 * Revision 1.3  1996/07/17 20:45:14  agray
 * cosmetic.
 *
 * Revision 1.2  1996/07/11 18:25:25  agray
 * added read/print/write_gauss_parms(), print/write_gauss_parms_set(),
 * print/write_unnorm_gauss_parms_set(), prob_gauss(), vector_prob_gauss().
 *
 * Revision 1.1  1996/02/21 05:18:05  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* UT library */
#include "ut_error.h"
#include "ut_math.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"
#include "float.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_util.h"
#include "da_io.h"
#include "da_linalg.h"
#include "da_random.h"

/* this module's header */
#include "da_probstat.h"


/*******************************************************************************
MEAN_VEC
Returns the mean of a vector.
AG
*******************************************************************************/
float mean_vec(float *v, int nc)
{
  float result;

  result = sum_vec(v, nc);

  return ( result / (float) nc );
}

/*******************************************************************************
MEAN_MAT
Returns the mean of a matrix.
RG
*******************************************************************************/
float mean_mat( float **m, int nr, int nc )
{
  float result;

  result = sum_mat( m, nr, nc );

  return( result / ((float) nr * nc) );
}

/*******************************************************************************
MEAN_IMAT
Returns the mean of an integer matrix.
RG
*******************************************************************************/
float mean_imat( int **m, int nr, int nc )
{
  float result;

  result = sum_imat( m, nr, nc );

  return( result / ((float) nr * nc) );
}

/*******************************************************************************
MEAN_CMAT
Returns the mean of an unsigned char matrix.
RG
*******************************************************************************/
float mean_cmat( unsigned char **m, int nr, int nc )
{
  float result;

  result = sum_cmat( m, nr, nc );

  return( result / ((float) nr * nc) );
}

/*******************************************************************************
STDEV_VEC
Returns the standard deviation of a vector.
RG
*******************************************************************************/
float stdev_vec(float *v, int nc, float mean)

{
  float result;

  result = var_vec( v, nc, mean );

  result = (float) sqrt( (double) result / (double) (nc - 1) );

  return ( result );
}

/*******************************************************************************
STDEV_MAT
Returns the standard deviation of a matrix.
RG
*******************************************************************************/
float stdev_mat( float **m, int nr, int nc, float mean )
{
  float result;

  result = var_mat( m, nr, nc, mean );

  result = (float) sqrt( (double) result / ((double) (nr * nc - 1)) );

  return( result );
}

/*******************************************************************************
VAR_VEC
Returns the variance of a vector.
RG
*******************************************************************************/
float var_vec( float *v, int nc, float mean )
{
  int   i;
  float result = 0.0;
  float *p;

  p = &v[1];

  for (i = 1; i <= nc; i++) {
    result += NR_sqr(*p - mean);
    p++;
  }

  return( result );
}

/*******************************************************************************
VAR_MAT
Returns the variance of a matrix.
RG
*******************************************************************************/
float var_mat( float **m, int nr, int nc, float mean )
{
  int   i, j;
  float result = 0.0;
  float *p;

  p = &m[1][1];

  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++) {
      result += NR_sqr(*p - mean);
      p++;
    }

  return( result );
}


/*******************************************************************************
EST_COVAR_MAT
Estimate covariance matrix given a data set and the mean of the data set.
AG
*******************************************************************************/
int est_covar_mat(float **covar, float **data, float *mean, int num_data,
                  int num_dims)

{
  int i, j, j2;
  float *data_i;

  for (i = 1; i <= num_data; i++)
  {
    /* add the square of the deviation vectors to the right position in the
       covariance matrix */
    for (j = 1; j <= num_dims; j++)
    {
      data_i = data[i];
      for (j2 = j; j2 <= num_dims; j2++)
        covar[j][j2] += (data_i[j] - mean[j]) * (data_i[j2] - mean[j2]);
    }
  }

  /* normalize the covariance matrices by the number of data */
  scalar_div_mat(covar, num_dims, num_dims, num_data);

  /* fill in the lower triangle of each covariance matrix based on upper 
     triangle */
  for (j = 1; j <= num_dims; j++)
    for (j2 = 1; j2 < j; j2++)
      covar[j][j2] = covar[j2][j];

  return (UT_OK);
}


/*******************************************************************************
STATS_OF_COLS
For each of k columns in a matrix, find the mean and standard deviation.  Pass
in a size k vector for each of these sets of values, and they will be filled in
by this function.  

Note:  This performs two passes over the data.
AG
*******************************************************************************/
int stats_of_cols(float **mat, float *means, float *stdevs, int numrows, 
                  int numcols)

{
  int   i,j;
  float *mat_row_i;

  /* initialize */
  fast_zero_vec(means, numcols);
  fast_zero_vec(stdevs, numcols);

  /* for each attribute, find the sum */
  for (i=1; i<=numrows; i++)
  {
    mat_row_i = mat[i];

    for (j=1; j<=numcols; j++)
      means[j] += mat_row_i[j];
  }
    
  /* for each attribute, normalize the sum */
  for (j=1; j<=numcols; j++)
    means[j] = means[j] / numrows;

  /* for each attribute, find the deviation from the mean */
  for (i=1; i<=numrows; i++)
  {
    mat_row_i = mat[i];

    for (j=1; j<=numcols; j++)
      stdevs[j] += NR_sqr(mat_row_i[j] - means[j]);
  }
    
  /* for each attribute, normalize the deviations */
  for (j=1; j<=numcols; j++)
    stdevs[j] = (float) sqrt( (double) stdevs[j] / (double) (numrows - 1) );
  
  return(UT_OK);
}


/*******************************************************************************
PROB_GAUSS
Compute the conditional probability of the datum assuming it was drawn from a 
multivariate normal distribution with the specified mean and covariance matrix,
using the normal probability density function.
RG
*******************************************************************************/
double prob_gauss(double *datum, double *mean, double **covar, int dim, 
                  double min_diag)

{
    int    i, j, k;
    double mahalanobis, prob;
    double **covar_inv, **covar_copy;
    double *diff, *tempvec;
    double denom;
    double sum;

    /* make temporary storage */
    diff =       NR_dvector(1,dim);
    tempvec =    NR_dvector(1,dim);
    covar_copy = NR_dmatrix(1,dim,1,dim);

    /* try to ensure that the covariance matrix is not ill-conditioned,
       before we try to invert it */
    restrict_illcond_dmatrix(covar, dim, min_diag);

    /* copy the covariance matrix */
    copy_dmat(covar, covar_copy, dim, dim);

    /* compute the cholesky decomposition L*L^T of the covariance matrix */
    NR_dcholdc(covar_copy, dim, tempvec);

    /* compute the inverse cholesky decomposition L^-1 
     * reference: Numerical Recipies, p. 98
     */
    for (i = 1; i <= dim; i++) {
      covar_copy[i][i] = 1.0 / tempvec[i];
      for (j = i+1; j <= dim; j++) {
	sum = 0.0;
	for (k = i; k < j; k++)
	  sum -= covar_copy[j][k] * covar_copy[k][i];
	covar_copy[j][i] = sum / tempvec[j];
      }
    }
    covar_inv = covar_copy;

    /* compute the sqare root of the determinant of the covariance matrix */
    denom = 1.0;
    for (i = 1; i <= dim; i++)
      denom *= tempvec[i];

    /* compute the constant factor */
    denom = pow((double)(2.0*UT_PI), (double)(dim/2.0)) * denom;

    /* compute Mahalanobis distance */
    subtract_dvec(datum, mean, diff, dim);
    for (i = 1; i <= dim; i++) {
      tempvec[i] = 0.0;
      for (j = 1; j <= i; j++)
	tempvec[i] += diff[j] * covar_inv[i][j];
    }
    mahalanobis = dot_product_dvec(tempvec, tempvec, dim);

    /* compute probability */
    prob = (double) exp(-0.5 * mahalanobis) / (double) denom;
    
    /* protect against returning unexpected zeros */
    if (prob == 0.0)
      prob = DBL_MIN;
    
    /* protect against returning infinity */
    if (prob > DBL_MAX)
      prob = DBL_MAX;

    /* clean up */
    NR_free_dmatrix(covar_inv,1,dim,1,dim);
    NR_free_dvector(tempvec,1,dim);
    NR_free_dvector(diff,1,dim);

    return (prob);
}


/*******************************************************************************
VECTOR_PROB_GAUSS
Compute the conditional probability of the datum assuming it was drawn from a 
multivariate normal distribution with the specified mean and covariance matrix,
using the normal probability density function.  Returns the vector of proba-
bilities, one for each datum.
RG
*******************************************************************************/
double* vector_prob_gauss(double **data, double *mean, double **covar, int dim, 
                         int numdata, double *probs, double min_diag,
                         char *rows_or_cols)

{
    int    i, j, k;
    double mahalanobis;
    double **covar_inv, **covar_copy;
    double *datum, *diff, *tempvec;
    double denom;
    double sum;
    boolean row_flag;

    /* make temporary storage */
    datum =      NR_dvector(1,dim);
    diff =       NR_dvector(1,dim);
    tempvec =    NR_dvector(1,dim);
    covar_copy = NR_dmatrix(1,dim,1,dim);

    /* set flag to avoid string comparisons */
    if (streq(rows_or_cols, "rows"))
      row_flag = UT_TRUE;
    else
      row_flag = UT_FALSE;

    /* try to ensure that the covariance matrix is not ill-conditioned,
       before we try to invert it */
    restrict_illcond_dmatrix(covar, dim, min_diag);

    /* copy the covariance matrix */
    copy_dmat(covar, covar_copy, dim, dim);

    /* compute the cholesky decomposition L*L^T of the covariance matrix */
    NR_dcholdc(covar_copy, dim, tempvec);

    /* compute the inverse cholesky decomposition L^-1
     * reference: Numerical Recipies, p. 98
     */
    for (i = 1; i <= dim; i++) {
      covar_copy[i][i] = 1.0 / tempvec[i];
      for (j = i+1; j <= dim; j++) {
        sum = 0.0;
        for (k = i; k < j; k++)
          sum -= covar_copy[j][k] * covar_copy[k][i];
        covar_copy[j][i] = sum / tempvec[j];
      }
    }
    covar_inv = covar_copy;

    /* compute the sqare root of the determinant of the covariance matrix */
    denom = 1.0;
    for (i = 1; i <= dim; i++)
      denom *= tempvec[i];

    /* compute the constant factor */
    denom = pow((double)(2.0*UT_PI), (double)(dim/2.0)) * denom;

    /* compute probability for each datum in the dataset */
    for (k = 1; k <= numdata; k++)
    {
      /* get the k'th datum from the dataset */
      if (row_flag == UT_TRUE)
        grab_row_dmat(data, k, dim, datum);
      else
        grab_col_dmat(data, k, dim, datum);

      /* compute Mahalanobis distance */
      subtract_dvec(datum, mean, diff, dim);
      for (i = 1; i <= dim; i++) {
        tempvec[i] = 0.0;
        for (j = 1; j <= i; j++)
          tempvec[i] += diff[j] * covar_inv[i][j];
      }
      mahalanobis = dot_product_dvec(tempvec, tempvec, dim);

      /* compute probability */
      probs[k] = (double) exp(-0.5 * mahalanobis) / (double) denom;

      /* protect against returning unexpected zeros */
      if (probs[k] == 0.0)
	probs[k] = DBL_MIN;

      /* protect against returning infinity */
      if (probs[k] > DBL_MAX)
        probs[k] = DBL_MAX;
    }

    /* clean up */
    NR_free_dvector(datum,1,dim);
    NR_free_dmatrix(covar_inv,1,dim,1,dim);
    NR_free_dvector(tempvec,1,dim);
    NR_free_dvector(diff,1,dim);

    return (probs);
}

/*******************************************************************************
GENERATE_DISC_RV_VALUE
Given the parameters of a discrete distribution, randomly generate a value
from that distribution.  Takes a vector containing the probability of each
discrete symbol, such that they sum to 1, and returns the index of the value
selected.

Given an empty vector to contain the cumulative symbol probabilities, this
function will fill in the cumulative probabilities and use these to sample
from the set of values.  Note that the pseudo-random number generator is 
assumed to have already been seeded.
AG
*******************************************************************************/
int generate_disc_rv_value(double *prob_symbol, int num_symbols, 
                           double *cum_prob_symbol)

{
  int i, index;
  double r;

  /* compute cumulative symbol probabilities */
  cum_prob_symbol[1] = prob_symbol[1];
  for (i = 2; i <= num_symbols; i++)
    cum_prob_symbol[i] = cum_prob_symbol[i - 1] + prob_symbol[i];

  /* select a value from the distribution */
  r = (double) gen_rand();

  for (i = 1; i <= num_symbols; i++)
    if (r <= cum_prob_symbol[i])
    {
      index = i;
      break;
    }
  if (r > cum_prob_symbol[num_symbols])
  {
    log_printf("Invariant not met in generate_disc_rv_value.\n");
    return(UT_ERROR);
  }

  return( index );
}


/*******************************************************************************
GENERATE_MIXTURE_RV_VALUE
Given the parameters of a continuous mixture distribution, randomly generate 
a value from that distribution.  Takes mixture weights, means and covariance
matrices defining a mixture.

Given an empty vector to contain the cumulative component probabilities, this
function will fill in the cumulative probabilities and use these to sample
from the set of values.  Note that the pseudo-random number generator is 
assumed to have already been seeded.
AG
*******************************************************************************/
int generate_mixture_rv_value(double *value, int num_dims, int num_comps, 
                              double **means, double ***covars, double *weights,
                              double *cum_prob_comps, double **chol_factor,
                              double *diag)

{
  int index;

  /* first choose which component will generate the value */
  index = generate_disc_rv_value(weights, num_comps, cum_prob_comps);

  /* now generate a value from the selected component */
  generate_gaussian_rv_value(value, num_dims, means[index], covars[index],
                             chol_factor, diag);

  return (UT_OK);
}


/*******************************************************************************
GENERATE_GAUSSIAN_RV_VALUE
Given the parameters of a Gaussian distribution, randomly generate a value
from that distribution.  Takes a mean vector and covariance matrix defining a
multivariate Gaussian distribution, and space for the Cholesky factor that must
be computed internally, a matrix having the same dimensions as the covariance
matrix, as well as space for the diagonal of the Cholesky factor.

Note that the pseudo-random number generator is assumed to have already been 
seeded, e.g. by using set_rand().
AG
*******************************************************************************/
int generate_gaussian_rv_value(double *value, int num_dims, double *mean, 
                               double **covar, double **chol_factor, double *diag)

{
  int i,j;

  /* compute Cholesky decomposition of the covariance matrix */
  copy_dmat(covar, chol_factor, num_dims, num_dims);
  NR_dcholdc(chol_factor, num_dims, diag);

  /* due to the way the answer is stored by choldc(), we must transform the
     result to get the true lower diagonal Cholesky factor we want */
  for (i = 1; i <= num_dims; i++)
    for (j = i; j <= num_dims; j++)
    {
      if (i == j)
        chol_factor[i][j] = diag[i];
      else  
        chol_factor[i][j] = 0.0;
    }

  /* generate a vector of independent standard normals */
  for (i = 1; i <= num_dims; i++)
    diag[i] = NR_dgasdev(&da_curr_rand);

  /* multiply by Cholesky factor to obtain a vector with the proper 
     covariance */
  right_mult_dmatrix(chol_factor, num_dims, num_dims, diag, value);

  /* add the mean to obtain a vector with the proper mean */
  add_dvec(value, mean, value, num_dims);

  return (UT_OK);
}


/*******************************************************************************
MERF
Compute the 'multidimensional erf() (error function)' for set of vectors in 
a dataset.  The output is a probability for each datum, stored in the vector
probs, which has nr elements.
AG
*******************************************************************************/
int merf(mean, cov, data, nr, nc, probs)

  double  *mean;     /* mean vector of gaussian distribution */
  double  **cov;     /* covariance matrix of gaussian distribution */
  double  **data;    /* set of data vectors (rows) to be evaluated */
  int    nr;        /* number of rows in data */
  int    nc;        /* number of dimensions of data */
  double  *probs;    /* holds the resulting probabilities; filled in by this 
                       function */
{
  double  **inv_cov, *datum, *t1, *t2;
  double  p_r, r2;
  int    i;

  /* Compute inverse of covariance matrix */
  inv_cov = NR_dmatrix(1, nc, 1, nc);
  invert_copy_alloc_sym_dmat(cov, nc, inv_cov);

  /* Compute probabilities using gamma function */

  /* allocate intermediate vectors */
  datum = NR_dvector(1, nc);
  t1 = NR_dvector(1, nc);
  t2 = NR_dvector(1, nc);
  
  for (i = 1; i <= nr; i++) {
    
    /* grab datum from data matrix */
    grab_row_dmat(data, i, nc, datum);
    
    /* transform each vector to standard multivariate normal */
    subtract_dvec(datum, mean, t1, nc);
    right_mult_dmatrix(inv_cov, nc, nc, t1, t2);
    
    /* perform radius substitution */
    r2 = (double) dot_product_dvec(t1, t2, nc);

    /* compute lower-case gamma function */
    p_r = NR_dgammp( (double) nc/2, (double) r2 / 2.0 );
    
    probs[i] = p_r;
  }

  /* free */
  NR_free_dvector(datum, 1, nc);
  NR_free_dvector(t1, 1, nc);
  NR_free_dvector(t2, 1, nc);
  NR_free_dmatrix(inv_cov, 1, nc, 1, nc);

  return (UT_OK);
}

/*******************************************************************************
PROB_MIXTURE
Compute the conditional probability of the datum assuming it was drawn from a 
mixture of multivariate normal distributions having the specified mean vectors
and covariance matrices.  If the weights sum to 1, this is a stochastic value,
i.e. this is a probability density function.
AG
*******************************************************************************/
double prob_mixture(double *datum, double **means, double ***covars, double *weights,
                   int dim, int num_comps, double min_diag)

{
  int i;
  double prob;
    
  prob = 0.0;
  for (i = 1; i <= num_comps; i++)
    prob += weights[i] * prob_gauss(datum, means[i], covars[i], dim, min_diag);
  
  return (prob);
}


/*******************************************************************************
PROB_MIXTURE_AND_KEEP
Compute the conditional probability of the datum assuming it was drawn from a 
mixture of multivariate normal distributions having the specified mean vectors
and covariance matrices.  If the weights sum to 1, this is a stochastic value,
i.e. this is a probability density function.

Same as prob_mixture(), except that the individual weighted component proba-
bility values are returned in the vector comp_wgtd_probs (this is the 'keep' 
option).
AG
*******************************************************************************/
double prob_mixture_and_keep(double *datum, double **means, double ***covars, 
                            double *weights, int dim, int num_comps, 
                            double min_diag, double *comp_wgtd_probs)

{
  int i;
  double prob;

  prob = 0.0;
  for (i = 1; i <= num_comps; i++)
  {
    comp_wgtd_probs[i] = weights[i] * 
                         prob_gauss(datum, means[i], covars[i], dim, min_diag);
    prob += comp_wgtd_probs[i];
  }

  return (prob);
}


/******************************************************************************o
LININT
Given a vector of x locations xa and y locations ya representing a function, 
perform a linear interpolation to calculate the value y of that function for 
the independent variable x.  The slop measures how close the value x must be
to the endpoints of the interval covered by xa in order to be successfully
evaluated (necessary for floating point "equality").
*******************************************************************************/
int linint( float *xa, float *ya, int n, float x, float *y, float slop )
{
  int   i;
  float interval;
 
  for (i = 1; xa[i] < x; i++)
    if (i == n+1) /* interp point is above data range */
      return( UT_ERROR );
 
  if (i == 1) { /* interp point is below data range */
    if ((x - xa[1]) < slop)
      i = 2;
    else
      return( UT_ERROR );
  }
 
  interval = xa[i] - xa[i-1];
  *y = ya[i-1] * ( xa[i] - x ) + ya[i] * ( x - xa[i-1] );
  *y /= interval;
 
  return( UT_OK );
}

/*******************************************************************************
VECTOR_LININT
A wrapper around the function linint so that it is called for a series of
values contained in the vector x; the interpolated function values are
returned in y.  The slop measures how close the values x must be to the 
endpoints of the interval covered by xa in order to be successfully
evaluated (necessary for floating point "equality").  Note that this code
is somewhat more inefficient than a special function designed to do this
thing without including the linint function.
*******************************************************************************/
int vector_linint( float *xa, float *ya, int n, float *x, float *y, int m,
                   float slop )
{
  int i;
  int error_code;
 
  for (i = 1; i <= m; i++) {
    error_code = linint( xa, ya, n, x[i], &y[i], slop );
    if (error_code == UT_ERROR)
      return( UT_ERROR );
  }
 
  return( UT_OK );
}


/*******************************************************************************
VECTOR_RATINT
A wrapper around the function ratint so that it is called for a series of
values contained in the vector x; the interpolated function values are
returned in y.  
*******************************************************************************/
int vector_ratint(float *xa, float *ya, int n, float *x, float *y, int m)
{
  int i;
  float error;
 
  for (i = 1; i <= m; i++)
    NR_ratint( xa, ya, n, x[i], &y[i], &error );
 
  return( UT_OK );
}


/*******************************************************************************
VECTOR_POLINT
A wrapper around the function polint so that it is called for a series of
values contained in the vector x; the interpolated function values are
returned in y.  
*******************************************************************************/
int vector_polint(float *xa, float *ya, int n, float *x, float *y, int m)
{
  int i;
  float error;
 
  for (i = 1; i <= m; i++)
    NR_polint( xa, ya, n, x[i], &y[i], &error );
 
  return( UT_OK );
}


/*******************************************************************************
FIT_LINE_TO_DATA
Given an index vector x and a corresponding data vector f(x), fit a line to it
and return the line's slope and intercept.
AG
******************************************************************************/
int fit_line_to_data(float *data, float *index, int num_data, float *slope, 
                     float *intercept)

{
  int mwt;
  float sig, siga, sigb, chi2, q;

  mwt = 0;  /* this says to ignore option of weighting the different data
               differently; fit() will thus ignore the contents of sig */

  /* perform linear regression */
  NR_fit (index, data, num_data, &sig, mwt, intercept, slope, &siga, &sigb, &chi2,
       &q);

  /* note that here we ignore the additional returned values siga, sigb, chi2, 
     and q */

  return (UT_OK);
}


/*******************************************************************************
COMPUTE_LINE
Given an index vector x, compute the function l(x), which is specified by
a slope and intercept.
AG
*******************************************************************************/
int compute_line(float *data, float *index, int num_data, float slope, 
                 float intercept)
                 
{
  int i;

  for (i = 1; i <= num_data; i++)
    data[i] = ( slope * index[i] ) + intercept;

  return(UT_OK);
}

/*******************************************************************************
COUNT_IN_INTERVALS
Compute a histogram given a vector of data.  Takes two other vectors, one 
containing the start points for each bucket, and one to contain the number of 
data in each bucket, which will be filled in by this function.
AG
*******************************************************************************/
int count_in_intervals(float *data, int *count, float *start_value,
                       int num_data, int num_buckets)

{
  int i, b;
  float curr_value;

  set_ivec(count, num_buckets, 0);

  for (i = 1; i <= num_data; i++)
  {
    curr_value = data[i];

    /* each bucket is defined by its starting value */
    for (b = 1; b < num_buckets; b++)
      if ( (curr_value >= start_value[b]) && (curr_value < start_value[b+1]) )
      {
        count[b]++;
        break;  /* from the inner for loop */
      }
    /* the last bucket includes everything up to infinity */
    if (curr_value >= start_value[num_buckets])
      count[num_buckets]++;
  }
  
  return(UT_OK);
}

/*******************************************************************************
PARTITION_VECTOR
Given a vector of data, determine the start-points of a given number of
equally-spaced subintervals.
AG
*******************************************************************************/
int partition_vector(float *data, float *start_values, int num_data, 
                     int num_intervals)

{
  int i;
  float range_data, min_data, max_data, interval_len;

  /* get interval length using range of data */
  min_data = min_vec(data, num_data);
  max_data = max_vec(data, num_data);
  range_data = max_data - min_data;

  interval_len = range_data / (float) num_intervals;

  /* compute start points */
  start_values[1] = min_data;
  for (i = 2; i <= num_intervals; i++)
    start_values[i] = start_values[i-1] + interval_len;

  return (UT_OK);
}

