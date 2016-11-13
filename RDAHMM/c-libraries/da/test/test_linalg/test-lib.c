/*******************************************************************************

  Title:     test-lib
  Author:    Alexander Gray
  Function:  Test functions in the Data Analysis Library.
  Reference: -
  How used:  To test correctness of the library after modifications.

  Compile:   make
  Example:   test-lib
  Notes:     

  Revisions: 2/1/96 AG created
             2/5/96 AG modified to match da lib. changes, esp. in memory alloc-
                       ation method.

*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: test-lib.c,v 1.6 1996/09/24 22:52:51 granat Exp $";
#endif
/* 
 * $Log: test-lib.c,v $
 * Revision 1.6  1996/09/24 22:52:51  granat
 * Changed the character matrix variable to be a matrix of unsigned chars;
 * Comment out some things that broke the compile.
 * The current version seg faults on "invert a matrix"
 *
 * Revision 1.5  1996/03/01 00:23:25  agray
 * added tests for new i/o functions from da_data.c for different data types.
 * ag
 *
 * Revision 1.4  1996/02/29  02:34:49  agray
 * changed call to read_bin_matrix() and write_bin_matrix();
 * note about read_landsat()
 * ag
 *
 * Revision 1.3  1996/02/21 03:39:12  agray
 * note about tests for pca, merf.
 * ag
 *
 * Revision 1.2  1996/02/21  00:38:34  agray
 * added write_bin_matrix() and read_bin_matrix()
 * ag
 *
 * Revision 1.1  1996/02/06  03:31:50  agray
 * Initial revision
 * */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "nr.h"
#include "nrutil.h"
#include "da.h"
#include "test-lib.h"

/*******************************************************************************
 MAIN
 Driver of program.
*******************************************************************************/
main(argc, argv)
    int   argc;
    char  *argv[];
{
    int     dim, numdata;
    float   a;
    float   *mean, *row, *col, *sum, *diff;
    float   **data, **cov, **inv_cov;
    float   **matrix1, **matrix2, **matrix3, **matrix4, **matrix5;
    float   **mat_trans;
    unsigned char    **cmatrix1;
    int     **imatrix1;
    double  **dmatrix1;
    
    dim = 3;
/*----------------------------------------------------------------------------*/
    /* float matrix i/o functions */
      printf("START ================================================\n");
      printf("--> Reading ascii matrix into floats.\n");
    read_matrix("test.matrix", dim, dim, &matrix1);
      printf("--> Printing matrix of floats.\n");
    print_matrix(stdout, dim, dim, matrix1);
      printf("--> Writing matrix of floats.\n");
    write_matrix("out.float.matrix", dim, dim, matrix1, "w");
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* char matrix i/o functions */
      printf("--> Reading ascii matrix into chars.\n");
    read_cmatrix("test.matrix", dim, dim, &cmatrix1);
      printf("--> Printing matrix of chars.\n");
    print_cmatrix(stdout, dim, dim, cmatrix1);
      printf("--> Writing matrix of chars.\n");
    write_cmatrix("out.char.matrix", dim, dim, cmatrix1, "w");
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* int matrix i/o functions */
      printf("--> Reading ascii matrix into ints.\n");
    read_imatrix("test.matrix", dim, dim, &imatrix1);
      printf("--> Printing matrix of ints.\n");
    print_imatrix(stdout, dim, dim, imatrix1);
      printf("--> Writing matrix of ints.\n");
    write_imatrix("out.int.matrix", dim, dim, imatrix1, "w");
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* double matrix i/o functions */
      printf("--> Reading ascii matrix into doubles.\n");
    read_dmatrix("test.matrix", dim, dim, &dmatrix1);
      printf("--> Printing matrix of doubles.\n");
    print_dmatrix(stdout, dim, dim, dmatrix1);
      printf("--> Writing matrix of doubles.\n");
    write_dmatrix("out.double.matrix", dim, dim, dmatrix1, "w");
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* binary float matrix i/o functions */
      printf("--> Writing float matrix in binary.\n");
    write_bin_matrix("out.bin.float.matrix", dim, dim, matrix1, "w");
      printf("--> Reading binary float matrix.\n");
    free_matrix(matrix1, 1, dim, 1, dim); /* since it'll be allocated again
                                             by read_bin_matrix() */
    read_bin_matrix("out.bin.float.matrix", dim, dim, &matrix1);
      printf("--> Printing float matrix.\n");
    print_matrix(stdout, dim, dim, matrix1);
      printf("======================================================\n");
      
/*----------------------------------------------------------------------------*/
    /* binary char matrix i/o functions */
      printf("--> Writing char matrix in binary.\n");
    write_bin_cmatrix("out.bin.char.matrix", dim, dim, cmatrix1, "w");
      printf("--> Reading binary char matrix.\n");
    free_cmatrix(cmatrix1, 1, dim, 1, dim); /* since it'll be allocated again
                                               by read_bin_cmatrix() */
    read_bin_cmatrix("out.bin.char.matrix", dim, dim, &cmatrix1);
      printf("--> Printing char matrix.\n");
    print_cmatrix(stdout, dim, dim, cmatrix1);
      printf("======================================================\n");
      
/*----------------------------------------------------------------------------*/
    /* binary int matrix i/o functions */
      printf("--> Writing int matrix in binary.\n");
    write_bin_imatrix("out.bin.int.matrix", dim, dim, imatrix1, "w");
      printf("--> Reading binary int matrix.\n");
    free_imatrix(imatrix1, 1, dim, 1, dim); /* since it'll be allocated again
                                             by read_bin_imatrix() */
    read_bin_imatrix("out.bin.int.matrix", dim, dim, &imatrix1);
      printf("--> Printing int matrix.\n");
    print_imatrix(stdout, dim, dim, imatrix1);
      printf("======================================================\n");
      
/*----------------------------------------------------------------------------*/
    /* binary double matrix i/o functions */
      printf("--> Writing double matrix in binary.\n");
    write_bin_dmatrix("out.bin.double.matrix", dim, dim, dmatrix1, "w");
      printf("--> Reading binary double matrix.\n");
    free_dmatrix(dmatrix1, 1, dim, 1, dim); /* since it'll be allocated again
                                             by read_bin_dmatrix() */
    read_bin_dmatrix("out.bin.double.matrix", dim, dim, &dmatrix1);
      printf("--> Printing double matrix.\n");
    print_dmatrix(stdout, dim, dim, dmatrix1);
      printf("======================================================\n");
      
/*----------------------------------------------------------------------------*/
    /* gaussian parameters i/o functions */
    /* create structures using nrutil.h functions */
/*
      printf("--> Reading Gaussian parameters.\n");
    read_gauss_parms("test.params", dim, &mean, &cov);
      printf("--> Writing Gaussian parameters.\n");
    write_gauss_parms("out.params", dim, mean, cov);
      printf("======================================================\n");
*/

/*----------------------------------------------------------------------------*/
    /* data i/o functions */
/*
      printf("--> Reading data.\n");
    read_data("test.data", &numdata, dim, &data);
      printf("--> Writing data.\n");
    write_data("out.data", numdata, dim, data);
      printf("======================================================\n");
*/

/*----------------------------------------------------------------------------*/
    /* rows and columns */
    row = vector(1,dim); 
    col = vector(1,dim); 
      printf("--> Grabbing row.\n");
    grab_row(matrix1, 1, dim, row);
      printf("--> Printing row.\n");
    print_row(stdout, dim, row);
      printf("--> Grabbing column.\n");
    grab_col(matrix1, 3, dim, col);
      printf("--> Printing column.\n");
    print_col(stdout, dim, col);
      printf("--> Writing row.\n");
    write_row("out.row", dim, row, "w");
      printf("--> Writing column.\n");
    write_col("out.col", dim, col, "w");
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* some vector operations */
      printf("--> Adding vectors.\n");
    sum = vector(1, dim);
    add_vec(dim, row, col, sum);
      print_row(stdout, dim, sum);
      printf("--> Subtracting vectors.\n");
    diff = vector(1, dim);
    subtract_vec(dim, row, col, diff);
      print_row(stdout, dim, diff);
      printf("--> Dot product of vectors.\n");
    a = dot_product(row, dim, col);
      printf("%f\n", a);
      printf("--> Outer product of vectors.\n");
    matrix2 = matrix(1,dim,1,dim);
    outer_product(dim, col, row, matrix2);
      print_matrix(stdout, dim, dim, matrix2);
      printf("--> Norm of a vector.\n");
    a = norm(dim, row);
      printf("%f\n", a);
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* some matrix operations */
      printf("--> Transpose of a matrix.\n");
    mat_trans = matrix(1,dim,1,dim);
    transpose(matrix1, dim, dim, mat_trans);
      print_matrix(stdout, dim, dim, mat_trans);
      printf("--> Matrix multiply.\n");
    matrix3 = matrix(1,dim,1,dim);
    mat_mult(matrix1, dim, dim, matrix2, dim, dim, matrix3);
      print_matrix(stdout, dim, dim, matrix3);
    /* note that invert_mat() destroys the original matrix, cov */
    inv_cov = matrix(1, dim, 1, dim);
      printf("--> Invert a matrix.\n");
    invert_mat(cov, dim, inv_cov);  
      print_matrix(stdout, dim, dim, inv_cov);
    matrix4 = matrix(1,dim,1,dim);      
    copy_mat(matrix1,matrix4,dim,dim);  /* copy this before it's destroyed */
      printf("--> Determinant of a matrix.\n");
    a = det(matrix4,dim);
      printf("%f\n", a);
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
    /* example of calling a Numerical Recipes routine */
    /* Gauss-Jordan elimination - p. 30, NR 2e
       note that gaussj() destroys the original matrices */
      printf("--> Copying a matrix.\n");
    copy_mat(matrix1,matrix3,dim,dim);  /* this should give identity matrix as 
                                           solution matrix result */
      printf("--> Call a Numerical Recipes routine.\n");
    gaussj(matrix1, dim, matrix3, dim);
      printf("Inverse of matrix:\n");
      print_matrix(stdout, dim, dim, matrix1);
      printf("Solution vectors:\n");
      print_matrix(stdout, dim, dim, matrix3);
      printf("======================================================\n");

/*----------------------------------------------------------------------------*/
      /* pca() - see CoolTools */
      /* merf() - see CoolTools */
      /* read_landsat() - see pca in CoolTools */

/*----------------------------------------------------------------------------*/
    /* clean up memory */
      printf("--> Freeing vectors.\n");
    free_vector(mean, 1, dim);
    free_vector(row, 1, dim);
    free_vector(col, 1, dim);
    free_vector(sum, 1, dim);
    free_vector(diff, 1, dim);
      printf("--> Freeing matrices.\n");
    free_matrix(cov, 1, dim, 1, dim);
    free_matrix(inv_cov, 1, dim, 1, dim);
    free_matrix(data, 1, numdata, 1, dim);
    free_matrix(matrix1, 1, dim, 1, dim);
    free_matrix(matrix2, 1, dim, 1, dim);
    free_matrix(matrix3, 1, dim, 1, dim);
    free_matrix(matrix4, 1, dim, 1, dim);
    free_matrix(mat_trans, 1, dim, 1, dim);
    free_cmatrix(cmatrix1, 1, dim, 1, dim);
    free_imatrix(imatrix1, 1, dim, 1, dim);
    free_dmatrix(dmatrix1, 1, dim, 1, dim);
      printf("DONE =================================================\n");

    return(UT_OK);
}
