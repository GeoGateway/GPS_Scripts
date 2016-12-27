/*******************************************************************************
 
  Title:     test-new
  Author:    Robert Granat
  Function:  Test functions in the Data Analysis Library.
  Reference: -
  How used:  To test correctness of the library after modifications.
 
  Compile:   make
  Example:   test-lib
  Notes:
 
  Revisions: 3/26/97 RG created
 
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ut.h"
#include "nr.h"
#include "nrutil.h"
#include "da.h"

/*******************************************************************************
 MAIN
 Driver of program.
*******************************************************************************/
main(argc, argv)
    int   argc;
    char  *argv[];
{
    int     i, j;
    int     nr, nc, nn;
    float   **a, **b;
    float   *temp1, *temp2;
 
    nr = 5;
    nc = 3;
    nn = 5;
 
    a = matrix( 1, nr, 1, nc);
    b = matrix( 1, nn, 1, nn);
    temp1 = vector( 1, nr );
    temp2 = vector( 1, nc );
/*----------------------------------------------------------------------------*/
    printf("START ================================================\n");
 
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
        a[i][j] = i*nc + j;

    for (i = 1; i <= nn; i++)
      for (j = 1; j <= nn; j++)
        b[i][j] = i*nn + j;
 
    printf("--> Printing matrix.\n");
 
    print_matrix( stdout, nr, nc, a );
 
    printf("--> In situ transpose of matrix.\n");
 
    transpose_in_situ_alloc_matrix( &a, nr, nc, temp1 );
 
    print_matrix( stdout, nc, nr, a );
 
    printf("--> In situ transpose of matrix (again).\n");
 
    transpose_in_situ_alloc_matrix( &a, nc, nr, temp2 );
 
    print_matrix( stdout, nr, nc, a );

    printf("--> Flipping matrix (left/right).\n");

    flip_left_right_matrix( a, nr, nc );

    print_matrix( stdout, nr, nc, a );

    printf("--> Flipping matrix (top/bottom).\n");

    flip_top_bottom_matrix( a, nr, nc );

    print_matrix( stdout, nr, nc, a );

    printf("--> Printing square matrix.\n");

    print_matrix( stdout, nn, nn, b );

    printf("--> In situ transpose of square matrix.\n");

    transpose_in_situ_sqr_matrix( b, nn );

    print_matrix( stdout, nn, nn, b );

}

