/*******************************************************************************
MODULE HEADER:
da_io.h
*******************************************************************************/

#ifndef _DA_IO_H_
#define _DA_IO_H_
/* Protects from multiple inclusion. */

#define VIC_HEADER_LABEL "LBLSIZE="

#ifndef lint
static char da_io_h_rcsid[] = "$Id: da_io.h,v 1.20 2000/03/31 01:38:17 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_io.h,v $
 * Revision 1.20  2000/03/31 01:38:17  granat
 * added prototypes for double precision sets of matrices functions
 *
 * Revision 1.19  1998/05/01 20:51:04  granat
 * added prototypes for unctions to read/write unsigned char matrices
 * into/from matrices of doubles and floats
 *
 * Revision 1.18  1997/10/21 14:52:47  granat
 * added prototypes to match changes in da_io.c
 *
 * Revision 1.17  1997/06/20 22:16:21  granat
 * fixed bugs introduced by auto-editing script
 *
 * Revision 1.16  1997/03/27 23:02:36  granat
 * added prototypes for skip_header() and read_vicar()
 *
 * Revision 1.15  1997/01/29 21:30:24  agray
 * new format.  also added write_matrix_transpose().
 *
 * Revision 1.14  1996/10/31 02:13:49  agray
 * renamed from "da_data" to "da_io";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.13  1996/09/27 17:57:45  agray
 * fixed log message.
 *
 * Revision 1.12  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.11  1996/09/18 16:34:37  granat
 * Fixed "char" matrix functions so that they use the "unsigned char" type
 *
 * Revision 1.10  1996/08/28 20:15:02  agray
 * changed read/write_data() to read/write_lc_matrix(); added read/write/print_
 * indexed_matrix().
 *
 * Revision 1.9  1996/07/15 18:56:06  agray
 * updated k&r prototypes.
 *
 * Revision 1.8  1996/07/15 18:17:57  agray
 * added  write/read_bin_row/col/irow/icol();
 * renamed read_bin_subset2matrix().
 *
 * Revision 1.7  1996/07/15 17:25:26  agray
 * rearranged order of functions; added read_row/col/irow/icol().
 *
 * Revision 1.6  1996/07/13 01:07:45  agray
 * moved in print/write_row/col/irow/icol() from da_linalg module
 *
 * Revision 1.5  1996/07/11 17:59:41  agray
 * added read_matrix_transpose(), write_comment().  some cosmetic changes.
 * /
 *
 * Revision 1.4  1996/05/10 21:06:39  granat
 * added prototypes for read_subset2matrix and read_binary2matrix
 * rg
 *
 * Revision 1.3  1996/03/01 00:20:02  agray
 * matched changes to da_data.c
 * ag
 *
 * Revision 1.2  1996/02/29  02:32:50  agray
 * moved write_bin_matrix() and read_bin_matrix() to this module, added type
 * argument to them
 * ag
 *
 * Revision 1.1  1996/02/28  04:55:24  agray
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

/* ascii matrices */

int read_matrix(char *infile, int nr, int nc, float ***vals);
int read_cmatrix(char *infile, int nr, int nc, unsigned char ***vals);
int read_imatrix(char *infile, int nr, int nc, int ***vals);
int read_dmatrix(char *infile, int nr, int nc, double ***vals);

int read_cmatrix_as_matrix(char *infile, int nr, int nc, float ***vals);
int read_cmatrix_as_dmatrix(char *infile, int nr, int nc, double ***vals);

int print_matrix(FILE* stream, int nr, int nc, float **mat);
int print_cmatrix(FILE* stream, int nr, int nc, unsigned char **mat);
int print_imatrix(FILE* stream, int nr, int nc, int **mat);
int print_dmatrix(FILE* stream, int nr, int nc, double **mat);

int write_matrix(char *outfile, int nr, int nc, float **mat, char *mode);
int write_cmatrix(char *outfile, int nr, int nc, unsigned char **mat, 
                  char *mode);
int write_imatrix(char *outfile, int nr, int nc, int **mat, char *mode);
int write_dmatrix(char *outfile, int nr, int nc, double **mat, char *mode);

int print_matrix_as_cmatrix( FILE *stream, int nr, int nc, float **m );
int print_dmatrix_as_cmatrix( FILE *stream, int nr, int nc, double **m );

int write_matrix_as_cmatrix( char *outfile, int nr, int nc, float **m,
                             char *mode );
int write_dmatrix_as_cmatrix( char *outfile, int nr, int nc, double **m,
                              char *mode );

/* binary-format matrices */

int read_bin_matrix(char *infile, int nr, int nc, float ***vals);
int read_bin_cmatrix(char *infile, int nr, int nc, unsigned char ***vals);
int read_bin_imatrix(char *infile, int nr, int nc, int ***vals);
int read_bin_smatrix(char *infile, int nr, int nc, int ***vals);
int read_bin_dmatrix(char *infile, int nr, int nc, double ***vals);

int read_bin_cmatrix_as_matrix(char *infile, int nr, int nc, float ***vals);
int read_bin_cmatrix_as_dmatrix(char *infile, int nr, int nc, double ***vals);

int write_bin_matrix(char *outfile, int nr, int nc, float **mat, char *mode);
int write_bin_cmatrix(char *outfile, int nr, int nc, unsigned char **mat, 
                      char *mode);
int write_bin_imatrix(char *outfile, int nr, int nc, int **mat, char *mode);
int write_bin_dmatrix(char *outfile, int nr, int nc, double **mat, char *mode);

int write_matrix_as_bin_cmatrix(char *outfile, int nr, int nc, float **mat, 
                                char *mode);
int write_dmatrix_as_bin_cmatrix(char *outfile, int nr, int nc, double **mat, 
                                 char *mode);

/* ascii vectors, row and column formats */

int read_row(char* infile, int dim, float **vec);
int read_col(char* infile, int dim, float **vec);
int read_irow(char* infile, int dim, int **vec);
int read_icol(char* infile, int dim, int **vec);
int read_drow(char* infile, int dim, double **vec);
int read_dcol(char* infile, int dim, double **vec);
int read_crow(char* infile, int dim, unsigned char **vec);
int read_ccol(char* infile, int dim, unsigned char **vec);

int print_row(FILE* stream, int dim, float *v);
int print_col(FILE* stream, int dim, float *v);
int print_irow(FILE* stream, int dim, int *v);
int print_icol(FILE* stream, int dim, int *v);
int print_drow(FILE* stream, int dim, double *v);
int print_dcol(FILE* stream, int dim, double *v);
int print_crow(FILE* stream, int dim, unsigned char *v);
int print_ccol(FILE* stream, int dim, unsigned char *v);

int write_row(char *outfile, int dim, float *vec, char *mode);
int write_col(char *outfile, int dim, float *vec, char *mode);
int write_irow(char *outfile, int dim, int *vec, char *mode);
int write_icol(char *outfile, int dim, int *vec, char *mode);
int write_drow(char *outfile, int dim, double *vec, char *mode);
int write_dcol(char *outfile, int dim, double *vec, char *mode);
int write_crow(char *outfile, int dim, unsigned char *vec, char *mode);
int write_ccol(char *outfile, int dim, unsigned char *vec, char *mode);

/* binary-format vectors */

int read_bin_row(char* infile, int dim, float **vec);
int read_bin_col(char* infile, int dim, float **vec);
int read_bin_irow(char* infile, int dim, int **vec);
int read_bin_icol(char* infile, int dim, int **vec);
int read_bin_drow(char* infile, int dim, double **vec);
int read_bin_dcol(char* infile, int dim, double **vec);
int read_bin_crow(char* infile, int dim, unsigned char **vec);
int read_bin_ccol(char* infile, int dim, unsigned char **vec);

int write_bin_row(char *outfile, int dim, float *vec, char *mode);
int write_bin_col(char *outfile, int dim, float *vec, char *mode);
int write_bin_irow(char *outfile, int dim, int *vec, char *mode);
int write_bin_icol(char *outfile, int dim, int *vec, char *mode);
int write_bin_drow(char *outfile, int dim, double *vec, char *mode);
int write_bin_dcol(char *outfile, int dim, double *vec, char *mode);
int write_bin_crow(char *outfile, int dim, unsigned char *vec, char *mode);
int write_bin_ccol(char *outfile, int dim, unsigned char *vec, char *mode);

/* ascii matrices with an index column */

int read_indexed_matrix(char *infile, int nr, int nc, float ***vals, 
                        float **index);
int print_indexed_matrix(FILE* stream, int nr, int nc, float **mat,
                         float *index);
int write_indexed_matrix(char *outfile, int nr, int nc, float **mat,
                         float *index, char *mode);

/* ascii vectors with an index column */

int read_indexed_col(char *infile, int dim, float **vals, float **index);
int print_indexed_col(FILE *stream, int dim, float *vec, float *index);
int write_indexed_col(char *outfile, int dim, float *vec, float *index, 
                      char *mode);

/* ascii matrices with the line count as the first line */

int read_lc_matrix(char *infile, int *nr, int nc, float ***vals);
int write_lc_matrix(char *outfile, int nr, int nc, float **mat);

/* ascii composite structures */

int read_set_of_matrices(char *infile, int num_mats, int num_rows, 
                         int num_cols, float ****vals);
int read_set_of_dmatrices(char *infile, int num_mats, int num_rows, 
                          int num_cols, double ****vals);
int read_set_of_sets_of_matrices(char *infile, int num_sets, int num_mats, 
                                 int num_rows, int num_cols, float *****vals);
int read_set_of_sets_of_dmatrices(char *infile, int num_sets, int num_mats, 
                                  int num_rows, int num_cols, double *****vals);

int print_set_of_matrices(FILE *stream, float ***mat_set, int num_mats, 
                          int num_rows, int num_cols);
int print_set_of_dmatrices(FILE *stream, double ***mat_set, int num_mats, 
                           int num_rows, int num_cols);
int print_set_of_sets_of_matrices(FILE *stream, float ****mat_set_set, 
                                  int num_sets, int num_mats, int num_rows, 
                                  int num_cols);
int print_set_of_sets_of_dmatrices(FILE *stream, double ****mat_set_set, 
                                  int num_sets, int num_mats, int num_rows, 
                                  int num_cols);

int write_set_of_matrices(char *out_file, float ***mat_set, int num_mats, 
                          int num_rows, int num_cols, char *mode);
int write_set_of_dmatrices(char *out_file, double ***mat_set, int num_mats, 
                           int num_rows, int num_cols, char *mode);
int write_set_of_sets_of_matrices(char *out_file, float ****mat_set_set, 
                                  int num_sets, int num_mats, int num_rows,
                                  int num_cols, char *mode);
int write_set_of_sets_of_dmatrices(char *out_file, double ****mat_set_set, 
                                   int num_sets, int num_mats, int num_rows,
                                  int num_cols, char *mode);

/* ascii files containing Gaussian parameters */

int read_gauss_parms(char *parmsfile, int nc, double **mean, double ***cov);

int print_gauss_parms(FILE* stream, int nc, double *mean, double **cov);
int print_gauss_parms_set(FILE *stream, int nc, int K, double **means, 
                          double ***covars);
int print_unnorm_gauss_parms_set(FILE *stream, int nc, int K, double **means, 
                                 double ***covars, double *range, double *minval);

int write_gauss_parms(char *parmsfile, int nc, double *mean, double **cov,
                      char *mode);
int write_gauss_parms_set(char *parmsfile, int nc, int K, double **means, 
                          double ***covars, char *mode);
int write_unnorm_gauss_parms_set(char *parmsfile, int nc, int K, double **means, 
                                 double ***covars, double *range, double *minval, 
                                 char *mode);

/* manipulating parts of matrices/vectors */

int read_subset2matrix(char *infile, int startrow, int nr, int nc, 
                       float ***vals);
int read_subset2dmatrix(char *infile, int startrow, int nr, int nc, 
                       double ***vals);
int read_bin_subset2matrix(char *infile, int startrow, int nr, int nc, 
                           float ***vals);
int read_bin_subset2dmatrix(char *infile, int startrow, int nr, int nc, 
                           double ***vals);

/* applying format transformations */

int read_matrix_transpose(char *infile, int nr, int nc, float ***vals);
int write_matrix_transpose(char *outfile, int nr, int nc, float **vals,
                           char *mode);

/* specialized data formats */

int read_landsat(char *infile, int nr, int nc, int nb, float ***vals);
int skip_header(FILE *fp, char *file_type);
int read_vicar( char *infile, int nr, int nc, unsigned char ***vals );

#endif /* _DA_IO_H_ */

