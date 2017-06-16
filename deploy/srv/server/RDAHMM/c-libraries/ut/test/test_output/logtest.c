#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "/proj/ml/mltt/include/cp.h"
#include "log.h"

FILE* ut_log_fp;
int   ut_log_level;

main(int argc, char *argv[])

{
  int val1 = 5;
  int val2 = 6;

  ut_log_fp = stdout;
  ut_log_level = 3;

  print_log_header();
  err_printf(); 
  log_printf("new value = %d\n",val1);
}

