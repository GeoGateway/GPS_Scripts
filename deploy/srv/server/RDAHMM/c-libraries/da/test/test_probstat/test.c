#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrutil.h"
#include "nr.h"

#include "util.h"
#include "da.h"

main(int argc, char *argv[])

{
  float **covar, *mean, *datum, px;

  mean =  vector(1,1);
  covar = matrix(1,1,1,1);
  datum = vector(1,1);

  mean[1] = 0;
  covar[1][1] = 1;

  datum[1] = -3;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = -2;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = -1;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = -0;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = 1;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = 2;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
  datum[1] = 3;
  px = gauss_eval(datum, mean, covar, 1, .1); 
  printf("px = %g\n", px);
}

