/*******************************************************************************
 
  Title:     test-lib
  Author:    Robert Granat
  Function:  Test functions in the Data Analysis Library.
  Reference: -
  How used:  To test correctness of the library after modifications.
 
  Compile:   make
  Example:   test-lib
  Notes:
 
  Revisions: 4/02/97 RG created
 
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
    float   **a, **b, **c;
    float   *temp1, *temp2, *temp3;
 
    nr = 5;
    nc = 3;
    nn = 8;
 
    a = matrix( 1, nr, 1, nc);
    b = matrix( 1, nn, 1, nn);
    c = matrix( 1, nn, 1, nn);
    temp1 = vector( 1, nr );
    temp2 = vector( 1, nc );
    temp3 = vector( 1, 2 * nn );
/*----------------------------------------------------------------------------*/
    printf("START ================================================\n");
 
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
        a[i][j] = i*nc + j;

    for (i = 1; i <= nn; i++)
      for (j = 1; j <= nn; j++) {
        b[i][j] = i*nn + j;
        if (j % 2)
          c[i][j] = 0.5;
        else
          c[i][j] = -0.5;
      } 
 
    printf("--> Printing matrix.\n");
 
    print_matrix( stdout, nn, nn, c );
 
    printf("--> FFT resultant matrix.\n");

    real_fft_2D( c, temp3, nn, nn, 1);

    print_matrix( stdout, nn, nn, c );

    printf("--> FFT resulting spectrum vector.\n");

    print_col( stdout, 2 * nn, temp3 );

    printf("--> Matrix after reverse transform.\n");

    real_fft_2D( c, temp3, nn, nn, -1);

    print_matrix( stdout, nn, nn, c );
}

