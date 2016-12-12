#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"
#include "ut_file_io.h"

/*
gcc -o test1 -I/tools/code/c/include -L/tools/code/c/lib test.c -lut
*/

main()
{

  FILE* fp;
  int   s;
  char* sptr;

  s = UT_ERROR;
  sptr = NULL;
  ut_log_fp = stdout;

  report_if_bad(s,"s is bad\n");;
  report_if_null(sptr,"sptr is NULL\n");;

  fp = fopen_report_if_fail("bullshit", "r");
}
