/*******************************************************************************
MODULE HEADER:
da_random.h
*******************************************************************************/

#ifndef _DA_RANDOM_H_
#define _DA_RANDOM_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_random_h_rcsid[] = "$Id: da_random.h,v 1.7 1999/07/20 15:45:48 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_random.h,v $
 * Revision 1.7  1999/07/20 15:45:48  granat
 * added prototypes for random_dvector() and random_dmatrix().
 *
 * Revision 1.6  1998/03/09 00:40:40  granat
 * fixed boolean bug
 *
 * Revision 1.5  1997/06/18 19:08:08  agray
 * changed da_curr_rand from long* to long
 *
 * Revision 1.4  1997/01/29 21:34:31  agray
 * shifted all of ut_rand module into da_random.
 * this added da_curr_rand variable, set_rand_by_clock(), set_rand_by_time_of_day(),
 * and random_partition_vector().  changed random_select_rows() to
 * random_mark_vector().
 *
 * Revision 1.3  1996/10/31 02:19:10  agray
 * renamed from "da_rand" to "da_random";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 *
 * Revision 1.2  1996/09/27 17:55:46  agray
 * removed k&r prototypes so that this library can be linked with c++ code.
 *
 * Revision 1.1  1996/08/23 18:00:26  agray
 * Initial revision
 *
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

extern long  da_curr_rand;

/*==============================================================================
Function Declarations
==============================================================================*/

/* random number generation and seeding */

float set_rand(int seed);
float gen_rand();

float set_rand_by_clock();
float set_rand_by_time_of_day();

/* random structures */

int random_vector(float *vec, int dim);
int random_dvector(double *vec, int dim);
int random_matrix(float **mat, int num_rows, int num_cols);
int random_dmatrix(double **mat, int num_rows, int num_cols);

/* random selection */

int random_mark_vector (int pick_size, int num_elements, boolean *selected);
int random_select_rows(float **source_mat, int num_source_rows, int num_cols,
                       float **dest_mat, int num_dest_rows, boolean *selected);
int random_select_drows(double **source_mat, int num_source_rows, int num_cols,
                       double **dest_mat, int num_dest_rows, boolean *selected);
int random_partition_vector(int num_groups, int num_elements, int *labels,
                            boolean *marked, int *elements_left);


#endif /* _DA_RANDOM_H_ */
