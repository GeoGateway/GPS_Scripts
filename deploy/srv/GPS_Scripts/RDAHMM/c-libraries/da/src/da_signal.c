/*******************************************************************************
MODULE NAME
da_signal

ONE-LINE SYNOPSIS
General functions related to signal processing.

SCOPE OF THIS MODULE
Any functions relating directly to other signal processing concepts which have 
representative modules in this library should go in the appropriate module.  
Functions that apply more generally are intended to go here.

SEE ALSO
Because the definition of this module is quite broad, there is some potential
overlap with several other modules in this library.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/qf, RG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_signal.c,v 1.17 1998/07/02 01:14:39 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_signal.c,v $
 * Revision 1.17  1998/07/02 01:14:39  granat
 * added median filtering
 *
 * Revision 1.16  1998/06/29 22:10:13  granat
 * added functions to calculate energies
 * revised functions to perform convolution via the FFT
 *
 * Revision 1.15  1997/10/21 14:53:11  granat
 * fixed corr/conv bug
 * made improvements to speed up corr/conv functions
 *
 * Revision 1.14  1997/06/20 22:11:22  granat
 * filled out conv/corr functions, cosmetic changes
 *
 * Revision 1.13  1997/06/05 18:56:51  granat
 * added slow_full_conv_vector
 *
 * Revision 1.12  1997/06/02 15:53:55  granat
 * changed to use new NR naming convention
 *
 * Revision 1.11  1997/01/30 02:19:50  agray
 * fixed bugs in range_normalize/unnormalize_cols().
 *
 * Revision 1.10  1997/01/29 22:01:53  agray
 * reordered functions.
 *
 * Revision 1.9  1997/01/29 21:49:28  agray
 * added range_of_cols(), changed range_normalize_cols().
 *
 * Revision 1.8  1997/01/29 21:11:15  agray
 * lots of source code formatting changes, some hacks to
 * range_normalize_cols(), changed ut lib. interface.
 *
 * Revision 1.7  1996/10/31 02:19:54  agray
 * renamed from "da_dist" to "da_signal";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 * 
 * Revision 1.6  1996/09/13 01:14:56  agray
 * change name and comments for normalize_data().
 * 
 * Revision 1.5  1996/09/13 00:57:33  agray
 * cosmetic; reduced arguments of print_unnorm_cov_matrix().
 * 
 * Revision 1.4  1996/08/28 20:17:52  agray
 * cosmetic.
 * 
 * Revision 1.3  1996/07/19 17:57:30  agray
 * cosmetic, minor efficiency hack
 * 
 * Revision 1.2  1996/07/11 18:00:18  agray
 * added print_unnorm_matrix(), print_unnorm_cov_matrix(), print_unnorm_row(),
 * print_unnorm_col().
 * 
 * Revision 1.1  1996/05/06 23:21:57  agray
 * Initial revision
 * 
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da.h"

/* this module's header */
#include "da_signal.h"


/*******************************************************************************
ENERGY_VEC
Return the total energy in a vector of floats.
RG
*******************************************************************************/
float energy_vec(float *v, int n)
{
  float *p;
  float *p_end;
  float energy = 0.0;
 
  /* assign pointer to last element */
  p_end = &v[n];
 
  /* calculate energy */
  for (p = &v[1]; p <= p_end; p++)
    energy += NR_sqr(*p);
 
  return (energy);
}


/*******************************************************************************
ENERGY_DVEC
Return the total energy in a vector of doubles.
RG
*******************************************************************************/
double energy_dvec(double *v, int n)
{
  double *p;
  double *p_end;
  double energy = 0.0;
 
  /* assign pointer to last element */
  p_end = &v[n];
 
  /* calculate energy */
  for (p = &v[1]; p <= p_end; p++)
    energy += NR_sqr(*p);
 
  return (energy);
}
 

/*******************************************************************************
NORMALIZE_ENERGY_VEC
Normalize a vector (signal) so that it has energy 1, and return the energy
in that vector before normalization.
RG
*******************************************************************************/
float normalize_energy_vec(float *v, int n)
{
  float energy;
  
  energy = energy_vec(v, n);
  scalar_div_vec(v, n, sqrt((double) energy));
  
  return (energy);
}


/*******************************************************************************
NORMALIZE_ENERGY_DVEC
Normalize a vector (signal) of doubles so that it has energy 1, and return the 
energy in that vector before normalization.
RG
*******************************************************************************/
double normalize_energy_dvec(double *v, int n)
{
  double energy;
  
  energy = energy_dvec(v, n);
  scalar_div_dvec(v, n, sqrt(energy));
  
  return (energy);
}


/*******************************************************************************
REAL_FFT_2D
Calculate the fast fourier transform of a matrix, or the inverse fast fourier
transform of a matrix in transform space.  Hacked from NR function nrlft3.
RG
*******************************************************************************/
int real_fft_2D(float **data, float *speq, int nn2, int nn3, int isign)
{
unsigned long i2,i3,j1,j2,j3,nn[4],ii3;
double theta,wi,wpi,wpr,wr,wtemp;
float c1,c2,h1r,h1i,h2r,h2i;

if (1+&data[nn2][nn3]-&data[1][1] != nn2*nn3)
  NR_error("rlft3: problem with dimensions or contiguity of data array\n");
  c1=0.5;
  c2 = -0.5*isign;
  theta=isign*(6.28318530717959/nn3);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  nn[1]=nn2;
  nn[2]=nn3 >> 1;
  if (isign == 1) {
    NR_fourn(&data[1][1]-1,nn,2,isign);
    for (i2=1,j2=0;i2<=nn2;i2++) {
      speq[++j2]=data[i2][1];
      speq[++j2]=data[i2][2];
    }
  }
  j1=1;
  wr=1.0;
  wi=0.0;
  for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
    for (i2=1;i2<=nn2;i2++) {
      if (i3 == 1) {
        j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
        h1r=c1*(data[i2][1]+speq[j2]);
        h1i=c1*(data[i2][2]-speq[j2+1]);
        h2i=c2*(data[i2][1]-speq[j2]);
        h2r= -c2*(data[i2][2]+speq[j2+1]);
        data[i2][1]=h1r+h2r;
        data[i2][2]=h1i+h2i;
        speq[j2]=h1r-h2r;
        speq[j2+1]=h2i-h1i;
      } else {
        j2=(i2 != 1 ? nn2-i2+2 : 1);
        j3=nn3+3-(i3<<1);
        h1r=c1*(data[i2][ii3]+data[j2][j3]);
        h1i=c1*(data[i2][ii3+1]-data[j2][j3+1]);
        h2i=c2*(data[i2][ii3]-data[j2][j3]);
        h2r= -c2*(data[i2][ii3+1]+data[j2][j3+1]);
        data[i2][ii3]=h1r+wr*h2r-wi*h2i;
        data[i2][ii3+1]=h1i+wr*h2i+wi*h2r;
        data[j2][j3]=h1r-wr*h2r+wi*h2i;
        data[j2][j3+1]= -h1i+wr*h2i+wi*h2r;
      }
    }
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == -1)
    NR_fourn(&data[1][1]-1,nn,2,isign);

  return (UT_OK);
}

/*******************************************************************************
FFT_CONV_MATRIX
Convolve two matrices using the fft.  The two input matrices may be of 
arbitrary size, but the function requires two workspace matrices and two 
workspace vectors be supplied.

If the input matrices are of dimension (M x N) and (U x V), the workspace
matrices must have minimum dimensions of the nearest powers of two greater
than or equal to (M + U - 1) x (N + V - 1).  The workspace vectors must
be of length equal to twice the number of rows of the workspace matrices.

The output matrix should be no larger than (M + U - 1) x (N + V - 1).
RG
*******************************************************************************/
int fft_conv_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                    int nc_b, float **c, int nr_c, int nc_c, float **data1,
		    float **data2, float *speq1, float *speq2)
{
  int    i;
  int    nr_radix2, nc_radix2;
  int    row_start, col_start;
  int    max_index;
  float  real, imag;
  float  scale;
  float *p1_real, *p1_imag;
  float *p2_real, *p2_imag;

  /* recalculate rows and columns of the temporary matrices */
  nr_radix2 = ((nr_a + nr_b - 1) >> 1) << 1;
  nc_radix2 = ((nc_a + nc_b - 1) >> 1) << 1;
  
  /* initialize the temporary matrices */
  fast_zero_mat(data1, nr_radix2, nc_radix2);
  fast_zero_mat(data2, nr_radix2, nc_radix2);

  /* transfer values of matrix a to temporary array with zero-padding */
  row_start = ((nr_radix2 - nr_a) >> 1) + 1;
  col_start = ((nc_radix2 - nc_a) >> 1) + 1;
  copy_mat_section(a, data1, 1, 1, row_start, col_start, nr_a, nc_a);

  /* transfer values of matrix b to temporary array with zero-padding */
  row_start = ((nr_radix2 - nr_b) >> 1) + 1;
  col_start = ((nc_radix2 - nc_b) >> 1) + 1;
  copy_mat_section(b, data2, 1, 1, row_start, col_start, nr_b, nc_b);
  
  /* perform the fft of the two data matrices */
  real_fft_2D(data1, speq1, nr_radix2, nc_radix2, 1);
  real_fft_2D(data2, speq2, nr_radix2, nc_radix2, 1);
  
  /* multiply in the frequency domain */
  max_index = nr_radix2 * nc_radix2 >> 1;
  scale = 2.0 / (nr_radix2 * nc_radix2);
  
  p1_real = &data1[1][1];
  p1_imag = &data1[1][2];

  p2_real = &data2[1][1];
  p2_imag = &data2[1][2];

  for (i = 1; i <= max_index; i++) {
    real = (*p1_real) * (*p2_real) - (*p1_imag) * (*p2_imag);
    imag = (*p1_real) * (*p2_imag) + (*p1_imag) * (*p2_real);
    *p1_real = scale * real;
    *p1_imag = scale * imag;
    p1_real += 2;
    p1_imag += 2;
    p2_real += 2;
    p2_imag += 2;
  }
  
  p1_real = &speq1[1];
  p1_imag = &speq1[2];

  p2_real = &speq2[1];
  p2_imag = &speq2[2];
  
  for (i = 1; i <= nr_radix2; i++) {
    real = (*p1_real) * (*p2_real) - (*p1_imag) * (*p2_imag);
    imag = (*p1_real) * (*p2_imag) + (*p1_imag) * (*p2_real);
    *p1_real = scale * real;
    *p1_imag = scale * imag;
    p1_real += 2;
    p1_imag += 2;
    p2_real += 2;
    p2_imag += 2;
  }
  
  /* perform the reverse transform */
  real_fft_2D(data1, speq1, nr_radix2, nc_radix2, 1);
  
  /* copy the data from the temporary matrix to the output matrix */
  row_start = ((nr_radix2 - nr_c) >> 1) + 1;
  col_start = ((nc_radix2 - nc_c) >> 1) + 1;
  copy_mat_section(data1, c, row_start, col_start, 1, 1, nr_c, nc_c);

  return (UT_OK);
}

/*******************************************************************************
FFT_CONV_ALLOC_MATRIX
Convolve two matrices using the fft.  This routine simply allocates workspace
and calls fft_conv_matrix().

See fft_conv_matrix().
*******************************************************************************/
int fft_conv_alloc_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                          int nc_b, float **c, int nr_c, int nc_c)
{
  int     nr_radix2, nc_radix2;
  float **data1, **data2;
  float  *speq1, *speq2;
  
  /* calculate the rows and columns of the temporary matrices */
  nr_radix2 = ((nr_a + nr_b - 1) >> 1) << 1;
  nc_radix2 = ((nc_a + nc_b - 1) >> 1) << 1;
  
  /* allocate the temporary matrices and vectors */
  data1 = NR_matrix(1, nr_radix2, 1, nc_radix2);
  data2 = NR_matrix(1, nr_radix2, 1, nc_radix2);
  
  speq1 = NR_vector(1, 2 * nr_radix2);
  speq2 = NR_vector(1, 2 * nr_radix2);
  
  /* perform the convolution */
  fft_conv_matrix(a, nr_a, nc_a, b, nr_b, nc_b, c, nr_c, nc_c, data1, data2,
                  speq1, speq2);
  
  /* free the workspace */
  NR_free_matrix(data1, 1, nr_radix2, 1, nc_radix2);
  NR_free_matrix(data2, 1, nr_radix2, 1, nc_radix2);
  
  NR_free_vector(speq1, 1, 2 * nr_radix2);
  NR_free_vector(speq2, 1, 2 * nr_radix2);
  
  return (UT_OK);
}

/*******************************************************************************
SLOW_FULL_CONV_VECTOR
Perform a 1D convolution of two vectors of length N and M, and return the 
result in a vector of size (N+M-1).  There are no approximations made to speed 
up the computation, but elements of the first vector are checked to make sure 
they are non-zero before operations are performed.  Therefore, the vector with 
the greater number of non-zero elements should be input as the first vector.
*******************************************************************************/
int slow_full_conv_vector(float *a, int n_a, float *b, int n_b, float *c, 
                          int n_c)
{
  int    i, j;
  float  *p_a, *p_b, *p_c;

  fast_zero_vec(c, n_c);

  p_a = &a[1];
  for (i = 1; i <= n_a; i++) {
    if (*p_a != 0.0) {
      p_b = &b[1];
      p_c = &c[i];
      for (j = 1; j <= n_b; j++) {
        *p_c += (*p_a) * (*p_b);
        p_b++;
        p_c++;
      }
    }
    p_a++;
  }

  return (UT_OK);
}


/*******************************************************************************
SLOW_CROP_CONV_VECTOR
Perform a 1D convolution of two vectors with dimensions N and U, and return
the result in a vector of length T, where T <= (N + U - 1).  See
slow_full_conv_vector.
RG
*******************************************************************************/
int slow_crop_conv_vector( float *a, int n_a, float *b, int n_b, float *c,
                           int n_c )
{
  int i, j, k;
  int n_full;
  int n_diff;
  int start;
 
  n_full = n_a + n_b - 1;
 
  n_diff = n_full - n_c;
 
  start = (n_diff >> 1) + 1;
 
  fast_zero_vec( c, n_c );
 
  for (i = 1; i <= n_a; i++)
    if (a[i] != 0.0)
      for (j = 1; j <= n_b; j++) {
        k = i + j - start;
        if ((k >= 1) && (k <= n_c))
          c[k] += a[i] * b[j];
      }
 
  return( UT_OK );
}

/*******************************************************************************
SLOW_FULL_CONV_MATRIX
Perform a 2D convolution of two matrices with dimensions NxM and UxV, and
return the result in a matrix of size (N+U-1)x(M+V-1).  There are no
approximations made to speed up the computation, but elements of the first
matrix are checked to make sure they are non-zero before operations are
performed.  Therefore, the matrix with the greater number of non-zero elements
should be input as the first matrix.
*******************************************************************************/
int slow_full_conv_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b, 
                          int nc_b, float **c, int nr_c, int nc_c) 
{ 
  int    i, j, k, l;
  float  *p_a, *p_b, *p_c;

  fast_zero_mat(c, nr_c, nc_c);

  p_a = &a[1][1];
  for (i = 1; i <= nr_a; i++)
    for (j = 1; j <= nc_a; j++) {
      if (*p_a != 0.0) {
        p_b = &b[1][1];
        for (k = 1; k <= nr_b; k++) {
          p_c = &c[i+k-1][j];
          for (l = 1; l <= nc_b; l++) {
            *p_c += (*p_a) * (*p_b);
            p_b++;
            p_c++;
          }
        }
      }
      p_a++;
    }

  return (UT_OK);
}

/*******************************************************************************
SLOW_FULL_NORM_CONV_MATRIX
Perform a 2D convolution of two matrices with dimensions NxM and UxV, and
return the normalized result in a matrix of size (N+U-1)x(M+V-1).  See 
slow_full_conv_matrix().
*******************************************************************************/
int slow_full_norm_conv_matrix(float **a, int nr_a, int nc_a, float **b, 
                               int nr_b, int nc_b, float **c, int nr_c, 
                               int nc_c)
{
  int    i, j, k, l;
  float  *p_a, *p_b, *p_c;
  float  mean_a, mean_b;
  float  var_a, var_b;
  float  norm;
 
  /* calculate statistics of a and b */

  mean_a = mean_mat(a, nr_a, nc_a);
  mean_b = mean_mat(b, nr_b, nc_b);

  var_a = var_mat(a, nr_a, nc_a, mean_a);
  var_b = var_mat(b, nr_b, nc_b, mean_b);

  norm = (float) (sqrt((double) var_a) * sqrt((double) var_b));

  /* subtract off the mean to save operations in the inner loop */

  scalar_subtract_mat(a, nr_a, nc_a, mean_a);
  scalar_subtract_mat(b, nr_b, nc_b, mean_b);

  /* zero out the new matrix */

  fast_zero_mat(c, nr_c, nc_c);

  p_a = &a[1][1];
  for (i = 1; i <= nr_a; i++)
    for (j = 1; j <= nc_a; j++) {
      if (*p_a != 0.0) {
        p_b = &b[1][1];
        for (k = 1; k <= nr_b; k++) {
          p_c = &c[i+k-1][j];
          for (l = 1; l <= nc_b; l++) {
            *p_c += (*p_a) * (*p_b);
            p_b++;
            p_c++;
          }
        }
      }
      p_a++;
    }

  /* normalize the result */

  scalar_div_mat(c, nr_c, nc_c, norm);

  /* restore the input matrices */

  scalar_add_mat(a, mean_a, nr_a, nc_a);
  scalar_add_mat(b, mean_b, nr_b, nc_b);

  return (UT_OK);
}


/*******************************************************************************
SLOW_CROP_CONV_MATRIX
Perform a 2D convolution of two matrices with dimensions NxM and UxV, and
return the center of the result in a matrix of size T x S, where 
T < (N + U - 1) and S < (M + V - 1).  See slow_full_conv_matrix().
*******************************************************************************/
int slow_crop_conv_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b, 
                          int nc_b, float **c, int nr_c, int nc_c)
{
  int    i, j, k, l;
  int    nr_full, nc_full;
  int    nr_diff, nc_diff;
  int    r_start, c_start;
  int    r_end, c_end;
  int    row, col;
  int    full_row, full_col;
  float  *p_a, *p_b, *p_c;
 
  nr_full = nr_a + nr_b - 1;
  nc_full = nc_a + nc_b - 1;

  nr_diff = nr_full - nr_c;
  nc_diff = nc_full - nc_c;

  r_start = (int) floor(((float) nr_diff) / 2.0) + 1;
  c_start = (int) floor(((float) nc_diff) / 2.0) + 1;

  r_end = nr_full - (int) ceil(((float) nr_diff) / 2.0);
  c_end = nr_full - (int) ceil(((float) nc_diff) / 2.0);
  
  fast_zero_mat(c, nr_c, nc_c);
 
  for (i = 1; i <= nr_a; i++)
    if (i <= r_end) {
      p_a = &a[i][1];
      for (j = 1; j <= nc_a; j++)
        if (j <= c_end) {
          if (*p_a != 0.0) {
            if (j > c_start)
              col = j - c_start;
            else
              col = 1;
            for (k = 1, row = i - r_start + 1, full_row = i; \
                 full_row <= r_end; k++, row++, full_row++)
              if (row >= 1) {
                p_b = &b[k][1];
                p_c = &c[row][col];
                for (l = 1, full_col = j; full_col <= c_end; l++, full_col++) {
                  if (full_col >= c_start) {
                    *p_c += (*p_a) * (*p_b);
                    p_c++;
                  }
                  p_b++;
                }
              }
          }
          p_a++;
        }
    }

  return (UT_OK);
}

/*******************************************************************************
SLOW_FULL_CORR_VECTOR
Perform a 1D correlation of two vectors of length N and M, and return the 
result in a vector of size (N+M-1).  No tricks are used to make the computation
faster, but there are no approximations made.
*******************************************************************************/
int slow_full_corr_vector(float *a, int n_a, float *b, int n_b, float *c, 
                          int n_c)
{
  flip_vector(a, n_a);

  slow_full_conv_vector(a, n_a, b, n_b, c, n_c);

  flip_vector(a, n_a);
  
  return (UT_OK);
}

/*******************************************************************************
SLOW_FULL_CORR_MATRIX
Perform a 2D correlation of two matrices with dimensions NxM and UxV, and
return the result in a matrix of size (N+U-1)x(M+V-1).  No tricks are used 
to make the computation faster, but there are no approximations made.
*******************************************************************************/
int slow_full_corr_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                          int nc_b, float **c, int nr_c, int nc_c)
{
  flip_top_bottom_matrix(a, nr_a, nc_a);
  flip_left_right_matrix(a, nr_a, nc_a);

  slow_full_conv_matrix(a, nr_a, nc_a, b, nr_b, nc_b, c, nr_c, nc_c);

  flip_left_right_matrix(a, nr_a, nc_a);
  flip_top_bottom_matrix(a, nr_a, nc_a);

  return (UT_OK);
}

/*******************************************************************************
SLOW_CROP_CORR_MATRIX
Perform a 2D correlation of two matrices with dimensions NxM and UxV, and
return the center of the result in a matrix of size T x S, where 
T < (N + U - 1) and S < (M + V - 1).  No tricks are used to make the 
computation faster, but there are no approximations made.
*******************************************************************************/
int slow_crop_corr_matrix(float **a, int nr_a, int nc_a, float **b, int nr_b,
                          int nc_b, float **c, int nr_c, int nc_c)
{
  flip_top_bottom_matrix(a, nr_a, nc_a);
  flip_left_right_matrix(a, nr_a, nc_a);

  slow_crop_conv_matrix(a, nr_a, nc_a, b, nr_b, nc_b, c, nr_c, nc_c);

  flip_left_right_matrix(a, nr_a, nc_a);
  flip_top_bottom_matrix(a, nr_a, nc_a);

  return (UT_OK);
}

/*******************************************************************************
SLOW_FULL_NORM_CORR_MATRIX
Perform a 2D correlation of two matrices with dimensions NxM and UxV, and
return the normalized result in a matrix of size (N+U-1)x(M+V-1).  No tricks
are used to make the computation faster, but there are no approximations made.
*******************************************************************************/
int slow_full_norm_corr_matrix(float **a, int nr_a, int nc_a, float **b,
                               int nr_b, int nc_b, float **c, int nr_c, 
                               int nc_c)
{
  flip_top_bottom_matrix(a, nr_a, nc_a);
  flip_left_right_matrix(a, nr_a, nc_a);
 
  slow_full_norm_conv_matrix(a, nr_a, nc_a, b, nr_b, nc_b, c, nr_c, nc_c);

  flip_left_right_matrix(a, nr_a, nc_a);
  flip_top_bottom_matrix(a, nr_a, nc_a);

  return (UT_OK);
}


/*******************************************************************************
MEDIAN3_FILT_MATRIX
Perform median filtering on a matrix using a 3x3 filter.  Ignores the boundary
elements.  Some gimmicks are here to maximize the computation done per memory
access.
RG
*******************************************************************************/
int median3_filt_matrix(float **A, int nr, int nc, float **B)
{
  register int    i, j;
  register float *baseA0, *baseA1, *baseA2;
  register float  ax00, ax01, ax02;
  register float  ax10, ax11, ax12;
  register float  ax20, ax21, ax22;
  register float  b0, b1, b2, b3, b4;

  if ((nr > 2) && (nc > 2)) {

    /* copy the borders of the matrix */

    memcpy(&B[1][1], &A[1][1], nc * sizeof(float));
    memcpy(&B[nr][1], &A[nr][1], nc * sizeof(float));
    for (i = 2; i < nr; i++) {
      B[i][1] = A[i][1];
      B[i][nc] = A[i][nc];
    }

    /* perform median filtering on the internal region */

    baseA0 = &A[1][1];
    baseA1 = baseA0 + nc;
    baseA2 = baseA1 + nc;
    
    for (i = 2; i < nr; i++) {
      
      /* assign values for the 3x3 window at the start of the row */

      ax00 = baseA0[0];
      ax01 = baseA0[1];
      ax02 = baseA0[2];
      ax10 = baseA1[0];
      ax11 = baseA1[1];
      ax12 = baseA1[2];
      ax20 = baseA2[0];
      ax21 = baseA2[1];
      ax22 = baseA2[2];

      for (j = 2; j < nc; j++) {

	/* sort the elements */

        b0 = ax00;

        if (ax01 >= b0)
          b1 = ax01;
        else {
          b1 = b0;
          b0 = ax00;
        }
    
	if (ax02 >= b1)
	  b2 = ax02;
	else if (ax02 >= b0) {
	  b2 = b1;
	  b1 = ax02;
	}
	else {
	  b2 = b1;
	  b1 = b0;
	  b0 = ax02;
	}

	if (ax10 >= b2)
	  b3 = ax10;
	else if (ax10 >= b1) {
          b3 = b2;
	  b2 = ax10;
	}
        else if (ax10 > b0) {
	  b3 = b2;
	  b2 = b1;
	  b1 = ax10;
        }
        else {
          b3 = b2;
          b2 = b1;
          b1 = b0;
          b0 = ax10;
        }
    
	if (ax11 >= b3) 
	  b4 = ax11;
	else if (ax11 >= b2) {
	  b4 = b3;
	  b3 = ax11;
	}
	else if (ax11 >= b1) {
	  b4 = b3;
	  b3 = b2;
	  b2 = ax11;
    	}
        else if (ax11 >= b0) {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
	  b1 = ax11;
        }
        else {
          b4 = b3;
          b3 = b2;
          b2 = b1;
          b1 = b0;
          b0 = ax11;
        }
    
	if (ax12 >= b4);
	else if (ax12 >= b3)
	  b4 = ax12;
	else if (ax12 >= b2) {
	  b4 = b3;
	  b3 = ax12;
	}
	else if (ax12 >= b1) {
	  b4 = b3;
	  b3 = b2;
	  b2 = ax12;
	}
        else if (ax12 >= b0) {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
	  b1 = ax12;
        }
        else {
          b4 = b3;
          b3 = b2;
          b2 = b1;
          b1 = b0;
        }
    
	if (ax20 >= b4);
	else if (ax20 >= b3)
	  b4 = ax20;
	else if (ax20 >= b2) {
	  b4 = b3;
	  b3 = ax20;
	}
	else if (ax20 >= b1) {
          b4 = b3;
	  b3 = b2;
	  b2 = ax20;
	}
        else {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
        }
    
	if (ax21 >= b4);
	else if (ax21 >= b3)
	  b4 = ax21;
	else if (ax21 >= b2) {
	  b4 = b3;
	  b3 = ax21;
	}
	else if (ax21 >= b1) {
	  b4 = b3;
	  b3 = b2;
        }
    
	/* assign the median element */

	if (ax22 >= b4)
	  B[i][j] = b4;
	else if (ax22 >= b3)
	  B[i][j] = ax22;
	else
	  B[i][j] = b3;
	
	/* assign new values for the shifted 3x3 window */

	if (j < nc - 1) {
	  baseA0++;
	  baseA1++;
	  baseA2++;

	  ax00 = ax01;
	  ax01 = ax02;
	  ax02 = baseA0[2];
	  ax10 = ax11;
	  ax11 = ax12;
	  ax12 = baseA1[2];
	  ax20 = ax21;
	  ax21 = ax22;
	  ax22 = baseA2[2];
	}
      }

      baseA0 += 3;
      baseA1 += 3;
      baseA2 += 3;
    }
  }
  else {
    printf("Matrix too small for median filtering\n");
    return (UT_ERROR);
  }
  
  return (UT_OK);
}


/*******************************************************************************
MEDIAN3_FILT_DMATRIX
Perform median filtering on a matrix of doubles using a 3x3 filter.  Ignores 
the boundary elements.  Some gimmicks are here to maximize the computation 
done per memory access.
RG
*******************************************************************************/
int median3_filt_dmatrix(double **A, int nr, int nc, double **B)
{
  register int     i, j;
  register double *baseA0, *baseA1, *baseA2;
  register double  ax00, ax01, ax02;
  register double  ax10, ax11, ax12;
  register double  ax20, ax21, ax22;
  register double  b0, b1, b2, b3, b4;

  if ((nr > 2) && (nc > 2)) {

    /* copy the borders of the matrix */

    memcpy(&B[1][1], &A[1][1], nc * sizeof(double));
    memcpy(&B[nr][1], &A[nr][1], nc * sizeof(double));
    for (i = 2; i < nr; i++) {
      B[i][1] = A[i][1];
      B[i][nc] = A[i][nc];
    }

    /* perform median filtering on the internal region */

    baseA0 = &A[1][1];
    baseA1 = baseA0 + nc;
    baseA2 = baseA1 + nc;
    
    for (i = 2; i < nr; i++) {
      
      /* assign values for the 3x3 window at the start of the row */

      ax00 = baseA0[0];
      ax01 = baseA0[1];
      ax02 = baseA0[2];
      ax10 = baseA1[0];
      ax11 = baseA1[1];
      ax12 = baseA1[2];
      ax20 = baseA2[0];
      ax21 = baseA2[1];
      ax22 = baseA2[2];

      for (j = 2; j < nc; j++) {

	/* sort the elements */

        b0 = ax00;

        if (ax01 >= b0)
          b1 = ax01;
        else {
          b1 = b0;
          b0 = ax00;
        }
    
	if (ax02 >= b1)
	  b2 = ax02;
	else if (ax02 >= b0) {
	  b2 = b1;
	  b1 = ax02;
	}
	else {
	  b2 = b1;
	  b1 = b0;
	  b0 = ax02;
	}

	if (ax10 >= b2)
	  b3 = ax10;
	else if (ax10 >= b1) {
          b3 = b2;
	  b2 = ax10;
	}
        else if (ax10 > b0) {
	  b3 = b2;
	  b2 = b1;
	  b1 = ax10;
        }
        else {
          b3 = b2;
          b2 = b1;
          b1 = b0;
          b0 = ax10;
        }
    
	if (ax11 >= b3) 
	  b4 = ax11;
	else if (ax11 >= b2) {
	  b4 = b3;
	  b3 = ax11;
	}
	else if (ax11 >= b1) {
	  b4 = b3;
	  b3 = b2;
	  b2 = ax11;
    	}
        else if (ax11 >= b0) {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
	  b1 = ax11;
        }
        else {
          b4 = b3;
          b3 = b2;
          b2 = b1;
          b1 = b0;
          b0 = ax11;
        }
    
	if (ax12 >= b4);
	else if (ax12 >= b3)
	  b4 = ax12;
	else if (ax12 >= b2) {
	  b4 = b3;
	  b3 = ax12;
	}
	else if (ax12 >= b1) {
	  b4 = b3;
	  b3 = b2;
	  b2 = ax12;
	}
        else if (ax12 >= b0) {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
	  b1 = ax12;
        }
        else {
          b4 = b3;
          b3 = b2;
          b2 = b1;
          b1 = b0;
        }
    
	if (ax20 >= b4);
	else if (ax20 >= b3)
	  b4 = ax20;
	else if (ax20 >= b2) {
	  b4 = b3;
	  b3 = ax20;
	}
	else if (ax20 >= b1) {
          b4 = b3;
	  b3 = b2;
	  b2 = ax20;
	}
        else {
	  b4 = b3;
	  b3 = b2;
	  b2 = b1;
        }
    
	if (ax21 >= b4);
	else if (ax21 >= b3)
	  b4 = ax21;
	else if (ax21 >= b2) {
	  b4 = b3;
	  b3 = ax21;
	}
	else if (ax21 >= b1) {
	  b4 = b3;
	  b3 = b2;
        }
    
	/* assign the median element */

	if (ax22 >= b4)
	  B[i][j] = b4;
	else if (ax22 >= b3)
	  B[i][j] = ax22;
	else
	  B[i][j] = b3;
	
	/* assign new values for the shifted 3x3 window */

	if (j < nc - 1) {
	  baseA0++;
	  baseA1++;
	  baseA2++;

	  ax00 = ax01;
	  ax01 = ax02;
	  ax02 = baseA0[2];
	  ax10 = ax11;
	  ax11 = ax12;
	  ax12 = baseA1[2];
	  ax20 = ax21;
	  ax21 = ax22;
	  ax22 = baseA2[2];
	}
      }

      baseA0 += 3;
      baseA1 += 3;
      baseA2 += 3;
    }
  }
  else {
    printf("Matrix too small for median filtering\n");
    return (UT_ERROR);
  }
  
  return (UT_OK);
}
