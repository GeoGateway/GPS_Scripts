/*******************************************************************************
MODULE NAME
da_random

ONE-LINE SYNOPSIS
General functions related to randomness.

SCOPE OF THIS MODULE
Any functions relating directly to other randomness concepts which have 
representative modules in this library should go in the appropriate module.  
Functions that apply more generally are intended to go here.  For example,
simulation of data from HMM's belongs in da_hmm, since that is less general.

SEE ALSO
Generation of random variables from statistical distributions lies in 
da_probstat.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/sc, AG.

NOTES
-

AG
*******************************************************************************/

#ifndef lint
static char rcsid[] = "$Id: da_random.c,v 1.11 1999/07/20 15:46:17 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_random.c,v $
 * Revision 1.11  1999/07/20 15:46:17  granat
 * added random_dvector() and random_dmatrix()
 *
 * Revision 1.10  1999/07/12 22:27:33  granat
 * changed indexing to DA conventions
 *
 * Revision 1.9  1998/03/09 00:38:50  granat
 * fixed boolean bug
 *
 * Revision 1.8  1997/06/18 19:07:18  agray
 * changed da_curr_rand from long* to long, to avoid possible problem in set_rand()
 * where contents of pointer to first_rand get overwritten
 *
 * Revision 1.7  1997/06/02 15:54:50  granat
 * changed to use new NR naming convention
 *
 * Revision 1.6  1997/01/29 21:48:18  agray
 * new formatting, cleaning up debugging output using ut_output,
 * added functions from ut_rand - see da_random.h rcs comments.
 *
 * Revision 1.5  1996/10/31 02:19:10  agray
 * renamed from "da_rand" to "da_random";
 * changed .h and .c formats throughout library;
 * some reorganizing between modules;
 * added some functions from HMM project.
 * 
 * Revision 1.4  1996/09/06 22:18:08  agray
 * fixed random_matrix().
 * 
 * Revision 1.3  1996/08/28 23:18:28  agray
 * changed random number generation functions to the Sun form, from the Cray
 * form... later this should be made transparent.
 * 
 * Revision 1.2  1996/08/28 20:18:38  agray
 * removed some headers.
 * 
 * Revision 1.1  1996/08/23 18:00:13  agray
 * Initial revision
 * 
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* UT library */
#include "ut_types.h"
#include "ut_string.h"
#include "ut_error.h"
#include "ut_output.h"

/* NR library */
#include "nr.h"

/* DA library */
#include "da_platform.h"

/* machine-specific */
#if (DA_INCL_FILE_LOC == DA_INCL_FILE_NOT_IN_SYS)
  #include <types.h>
  #include <sys/time.h>
#else
  #include <sys/types.h>
  #include <sys/timeb.h>
  #include <sys/time.h>
#endif

/* this module's header */
#include "da_random.h"

long  da_curr_rand;

/*******************************************************************************
SET_RAND
Call this function to initialize the random number generator, to allow
successive calls to gen_rand().

We chose ran1() rather than, say, ran2(), because gasdev() uses ran1().
AG
*******************************************************************************/
float set_rand(int seed)

{
  long first_rand;

  /* make sure the first value given to ran1() is negative */
  first_rand = (long) ( -1 * abs(seed) );
  da_curr_rand = first_rand;

  return ( NR_ran1(&da_curr_rand) );
}


/*******************************************************************************
GEN_RAND
Assuming that the function set_rand() has been called to initialize the random
number generator, generate the next number in a pseudo-random sequence.

We chose ran1() rather than, say, ran2(), because gasdev() uses ran1().
AG
*******************************************************************************/
float gen_rand()

{
  return ( NR_ran1(&da_curr_rand) );
}


/*******************************************************************************
SET_RAND_BY_CLOCK
Sets the random seed by determining a seed number which is a function of the 
current time.  This function is useful for determining a seed when you don't 
want to pick one manually.
JR, mod. by AG
*******************************************************************************/
float set_rand_by_clock()

{
  time_t *tloc = NULL;
  time_t t;

  t = time(tloc);

  return (set_rand( (int)t ));
}


/*******************************************************************************
SET_RAND_BY_TIME_OF_DAY
Sets the random seed by determining a seed number which is a function of the 
current time.  This function is useful for determining a seed when you don't 
want to pick one manually.  Uses more possible seeds than set_rand_by_clock().
AG
*******************************************************************************/
float set_rand_by_time_of_day()

{
  struct timeval tp;

  gettimeofday(&tp,0);

  return (set_rand(tp.tv_usec));
}


/*******************************************************************************
RANDOM_VECTOR
Randomly set the elements of a vector, drawing from a 0-1 uniform distribution.
AG
*******************************************************************************/
int random_vector(float *vec, int dim)

{
  int i;

  for (i=1; i<=dim; i++)
    vec[i] = gen_rand();

  return (UT_OK);
}


/*******************************************************************************
RANDOM_DVECTOR
Randomly set the elements of a vector, drawing from a 0-1 uniform distribution.
AG
*******************************************************************************/
int random_dvector(double *vec, int dim)

{
  int i;

  for (i=1; i<=dim; i++)
    vec[i] = (double) gen_rand();

  return (UT_OK);
}


/*******************************************************************************
RANDOM_MATRIX
Randomly set the elements of a matrix, drawing from a 0-1 uniform distribution.
AG
*******************************************************************************/
int random_matrix(float **mat, int num_rows, int num_cols)

{
  int i,j;

  for (i=1; i<=num_rows; i++)
    for (j=1; j<=num_cols; j++)
      mat[i][j] = gen_rand();

  return (UT_OK);
}


/*******************************************************************************
RANDOM_DMATRIX
Randomly set the elements of a matrix, drawing from a 0-1 uniform distribution.
AG
*******************************************************************************/
int random_dmatrix(double **mat, int num_rows, int num_cols)

{
  int i,j;

  for (i=1; i<=num_rows; i++)
    for (j=1; j<=num_cols; j++)
      mat[i][j] = (double) gen_rand();

  return (UT_OK);
}


/*******************************************************************************
RANDOM_MARK_VECTOR
Mark pick_size randomly chosen elements of the pre-allocated array selected (of
size num_elements) as UT_TRUE (picked), and the remaining elements as UT_FALSE.
JR, mod. by AG
*******************************************************************************/
int random_mark_vector (int pick_size, int num_elements, boolean *selected)

{
  int  i, index;

  if (pick_size < num_elements/2)
  {
    for (i=1; i<=num_elements; i++)
      selected[i] = UT_FALSE;
    for (i=1; i<=pick_size; i++)
    {
      while (UT_TRUE)
      {
        index = (int) (gen_rand() * (float) num_elements) + 1;
        if (selected[index] == UT_FALSE)
        {
          selected[index] = UT_TRUE;
          break;
        }
      }
    }
  }
  else
  {
    for (i=1; i<=num_elements; i++)
      selected[i] = UT_TRUE;
    for (i=1; i<=num_elements-pick_size; i++)
    {
      while (UT_TRUE)
      {
        index = (int) (gen_rand() * (float) num_elements) + 1;
        if (selected[index] == UT_TRUE)
        {
          selected[index] = UT_FALSE;
          break;
        }
      }
    }
  }

  return (UT_OK);
}


/*******************************************************************************
RANDOM_SELECT_ROWS
Randomly select from the rows vectors of a source matrix, filling in a 
destination matrix with the selected rows.  The temporary vector selected must
be supplied; it is a vector of booleans of length num_source_rows.
AG
*******************************************************************************/
int random_select_rows(float **source_mat, int num_source_rows, int num_cols, 
                       float **dest_mat, int num_dest_rows, boolean *selected)

{
  int   i,j,k;

  /* randomly choose amongst the indices of the data */
  random_mark_vector(num_dest_rows, num_source_rows, selected);

  /* copy the chosen data to the mat matrix */
  k = 1;
  for (i=1; i<=num_source_rows; i++)
  {
    if (selected[i] == UT_TRUE)
    {
      for (j=1; j<=num_cols; j++)
        dest_mat[k][j] = source_mat[i][j];
      k++;
    }
  }
    
  return (UT_OK);
}


/*******************************************************************************
RANDOM_SELECT_DROWS
Randomly select from the rows vectors of a source matrix, filling in a 
destination matrix with the selected rows.  The temporary vector selected must
be supplied; it is a vector of booleans of length num_source_rows.
AG
*******************************************************************************/
int random_select_drows(double **source_mat, int num_source_rows, int num_cols, 
                       double **dest_mat, int num_dest_rows, boolean *selected)

{
  int   i,j,k;

  /* randomly choose amongst the indices of the data */
  random_mark_vector(num_dest_rows, num_source_rows, selected);

  /* copy the chosen data to the mat matrix */
  k = 1;
  for (i=1; i<=num_source_rows; i++)
  {
    if (selected[i] == UT_TRUE)
    {
      for (j=1; j<=num_cols; j++)
        dest_mat[k][j] = source_mat[i][j];
      k++;
    }
  }
    
  return (UT_OK);
}


/*******************************************************************************
RANDOM_PARTITION_VECTOR
Create a randomly-decided partitioning of an vector into k groups, by assigning 
each of the elements of the vector a number between 1 and k.  If N is the size 
of the vector, and N is not evenly divisible by k, there will be a remainder 
group, consisting of vector elements labelled k+1.  The input vector is assumed 
to have been pre-allocated.  Temporary storage marked, a vector of booleans of
length num_elements, and elements_left, a vector of ints of length num_elements,
must be passed in.
AG
*******************************************************************************/
int random_partition_vector(int num_groups, int num_elements, int *labels,
                            boolean *marked, int *elements_left)

{
  int  i, j, k, index, index2;
  int  curr_num_elements, group_size, remainder_size;

  /* initialize */
  for (i=1; i<=num_elements; i++)
  {
    marked[i] = UT_FALSE;
    elements_left[i] = i;
  }

  /* compute various sizes */
  curr_num_elements = num_elements;
  group_size = (int) (num_elements - (num_elements % num_groups)) / num_groups;
  remainder_size = (int) num_elements % num_groups;

  /* obtain the labellings for each group */
  for (k=1; k<=num_groups-1; k++)
  {
    /* only look at the elements that haven't been put into a group yet */
    for (j=1; j<=group_size; j++)
    {
      /* keep trying until we get a new random number */
      while(1)
      {
        index = (rand() % curr_num_elements);
        if (marked[elements_left[index]] == UT_FALSE)
        {
          marked[elements_left[index]] = UT_TRUE;
          break;
        }
      }
      labels[elements_left[index]] = k;
      /*      printf("item = %d\n",elements_left[index]); */
    }

    /* keep track of number of elements remaining to choose from */
    curr_num_elements -= group_size;

    /* update elements_left */
    index2 = 1;
    for (i=1; i<=num_elements; i++)
    {
      if (marked[i] == UT_FALSE)
      {
        elements_left[index2] = i;
        index2++;
      }
    }

    /* check that we have the number of elements we expected */
    if (index2 != curr_num_elements)
    {
      log_printf("Invariant not met in random_partition_vector().\n");
      return (UT_ERROR);
    }
  }

  /* now do the leftover elements, which make up the last group plus the
     'remainder' chunk; first do the last group: */
  for (j=1; j<=group_size; j++)
  {
    marked[elements_left[j]] = UT_TRUE;
    labels[elements_left[j]] = num_groups-1;
    /* printf("item = %d\n",elements_left[j]); */
  }

  /* now do the remainder chunk */
  for (j=group_size+1; j<=(group_size + remainder_size); j++)
  {
    marked[elements_left[j]] = UT_TRUE;
    labels[elements_left[j]] = num_groups;
    printf("item = %d\n",elements_left[j]);
  }

  return (UT_OK);
}

