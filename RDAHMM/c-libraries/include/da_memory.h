/*******************************************************************************
MODULE HEADER:
da_memory.h
*******************************************************************************/

#ifndef _DA_MEMORY_H_
#define _DA_MEMORY_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_memory_h_rcsid[] = "$Id: da_memory.h,v 1.1 1997/07/30 23:31:19 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_memory.h,v $
 * Revision 1.1  1997/07/30 23:31:19  agray
 * Initial revision
 *
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

/*******************************************************************************
MATRIX_REPORT_IF_FAIL, MATRIX_RETURN_IF_FAIL, MATRIX_EXIT_IF_FAIL
Attempt to allocate a matrix.  If the attempt fails, perform some error hand-
ling.  Analogous to the malloc_report_if_fail(), etc. macros.
*******************************************************************************/
#define matrix_report_if_fail(nr,nc) \
  ut_status_ptr = (void*) matrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d matrix\n", nr,nc); \
  report_if_null(ut_status_ptr, ut_err_msg);

#define matrix_return_if_fail(nr,nc,r) \
  ut_status_ptr = (void*) matrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d matrix\n", nr,nc); \
  return_if_null(ut_status_ptr, ut_err_msg, r);

#define matrix_exit_if_fail(nr,nc) \
  ut_status_ptr = (void*) matrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d matrix\n", nr,nc); \
  exit_if_null(ut_status_ptr, ut_err_msg);


/*******************************************************************************
DMATRIX_REPORT_IF_FAIL, DMATRIX_RETURN_IF_FAIL, DMATRIX_EXIT_IF_FAIL
Analogous to the matrix_report_if_fail(), etc. macros.
*******************************************************************************/
#define dmatrix_report_if_fail(nr,nc) \
  ut_status_ptr = (void*) dmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d dmatrix\n", nr,nc); \
  report_if_null(ut_status_ptr, ut_err_msg);

#define dmatrix_return_if_fail(nr,nc,r) \
  ut_status_ptr = (void*) dmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d dmatrix\n", nr,nc); \
  return_if_null(ut_status_ptr, ut_err_msg, r);

#define dmatrix_exit_if_fail(nr,nc) \
  ut_status_ptr = (void*) dmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d dmatrix\n", nr,nc); \
  exit_if_null(ut_status_ptr, ut_err_msg);

/*******************************************************************************
CMATRIX_REPORT_IF_FAIL, CMATRIX_RETURN_IF_FAIL, CMATRIX_EXIT_IF_FAIL
Analogous to the matrix_report_if_fail(), etc. macros.
*******************************************************************************/
#define cmatrix_report_if_fail(nr,nc) \
  ut_status_ptr = (void*) cmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d cmatrix\n", nr,nc); \
  report_if_null(ut_status_ptr, ut_err_msg);

#define cmatrix_return_if_fail(nr,nc,r) \
  ut_status_ptr = (void*) cmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d cmatrix\n", nr,nc); \
  return_if_null(ut_status_ptr, ut_err_msg, r);

#define cmatrix_exit_if_fail(nr,nc) \
  ut_status_ptr = (void*) cmatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d cmatrix\n", nr,nc); \
  exit_if_null(ut_status_ptr, ut_err_msg);

/*******************************************************************************
IMATRIX_REPORT_IF_FAIL, IMATRIX_RETURN_IF_FAIL, IMATRIX_EXIT_IF_FAIL
Analogous to the matrix_report_if_fail(), etc. macros.
*******************************************************************************/
#define imatrix_report_if_fail(nr,nc) \
  ut_status_ptr = (void*) imatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d imatrix\n", nr,nc); \
  report_if_null(ut_status_ptr, ut_err_msg);

#define imatrix_return_if_fail(nr,nc,r) \
  ut_status_ptr = (void*) imatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d imatrix\n", nr,nc); \
  return_if_null(ut_status_ptr, ut_err_msg, r);

#define imatrix_exit_if_fail(nr,nc) \
  ut_status_ptr = (void*) imatrix_and_track(nr,nc); \
  sprintf(ut_err_msg, "Couldn't allocate %d x %d imatrix\n", nr,nc); \
  exit_if_null(ut_status_ptr, ut_err_msg);

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* matrix allocation and de-allocation */

float** matrix_and_track(int nr, int nc);
int** imatrix_and_track(int nr, int nc);
unsigned char** cmatrix_and_track(int nr, int nc);
double** dmatrix_and_track(int nr, int nc);

void free_matrix_and_track(float **mat, int nr, int nc);
void free_cmatrix_and_track(unsigned char **mat, int nr, int nc);
void free_imatrix_and_track(int **mat, int nr, int nc);
void free_dmatrix_and_track(double **mat, int nr, int nc);

#endif /* _DA_MEMORY_H_ */
