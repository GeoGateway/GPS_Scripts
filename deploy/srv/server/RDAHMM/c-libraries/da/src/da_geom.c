/*******************************************************************************
MODULE NAME
da_geom

ONE-LINE SYNOPSIS
General functions related to geometry.

SCOPE OF THIS MODULE
Any functions that perform geometrical calculations should go here.

SEE ALSO
For the most part there, should little overlap with other modules.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
-

NOTES
-

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_geom.c,v 1.3 1998/03/09 00:39:12 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_geom.c,v $
 * Revision 1.3  1998/03/09 00:39:12  granat
 * fixed boolean bug
 *
 * Revision 1.2  1997/06/20 21:10:43  granat
 * cosmetic changes
 *
 * Revision 1.1  1997/06/20 20:54:31  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr.h"

/* DA library */

/* this module's header */
#include "da_geom.h"

/*******************************************************************************
DA_TABLE_SIN
Return the sine of the given angle by looking up an approximation in a pre-
calculated table.  If the table has not yet been calculated, generate it
(ie, on the first function call).
*******************************************************************************/
double DA_table_sin( double angle )
{
  static boolean  da_sin_init = UT_FALSE;
  static double   da_sin_table[DA_TRIG_TABLE_SIZE];

  if (da_sin_init) /* The table has already been created */
    return( da_sin_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  else {
    int i;
    da_sin_init = UT_TRUE;
    for (i = 0; i < DA_TRIG_TABLE_SIZE; i++)
      da_sin_table[i] = sin( i * DA_TRIG_SAMPLE_ANG );
    return( da_sin_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  }
}

/*******************************************************************************
DA_TABLE_COS
Return the cosine of the given angle by looking up an approximation in a pre-
calculated table.  If the table has not yet been calculated, generate it
(ie, on the first function call).
*******************************************************************************/
double DA_table_cos( double angle )
{
  static boolean  da_cos_init = UT_FALSE;
  static double   da_cos_table[DA_TRIG_TABLE_SIZE];
 
  if (da_cos_init) /* The table has already been created */
    return( da_cos_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  else {
    int i;
    da_cos_init = UT_TRUE;
    for (i = 0; i < DA_TRIG_TABLE_SIZE; i++)
      da_cos_table[i] = sin( i * DA_TRIG_SAMPLE_ANG );
    return( da_cos_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  }
}

/*******************************************************************************
DA_TABLE_TAN
Return the tangent of the given angle by looking up an approximation in a pre-
calculated table.  If the table has not yet been calculated, generate it
(ie, on the first function call).
*******************************************************************************/
double DA_table_tan( double angle )
{
  static boolean  da_tan_init = UT_FALSE;
  static double   da_tan_table[DA_TRIG_TABLE_SIZE];
 
  if (da_tan_init) /* The table has already been created */
    return( da_tan_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  else {
    int i;
    da_tan_init = UT_TRUE;
    for (i = 0; i < DA_TRIG_TABLE_SIZE; i++)
      da_tan_table[i] = sin( i * DA_TRIG_SAMPLE_ANG );
    return( da_tan_table[(int) (angle * DA_INV_TRIG_SAMPLE_ANG)] );
  }
}
