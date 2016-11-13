/*******************************************************************************
MODULE HEADER:
ut_time.h
*******************************************************************************/

#ifndef _UT_TIME_H_
#define _UT_TIME_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_time_h_rcsid[] = "$Id: ut_time.h,v 1.1 1997/01/29 23:44:54 agray Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_time.h,v $
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

#define UT_LEAP_DAY -29

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

int month_date_to_day (int month, int date);
int days_in_year (int year);
int days_in_year_range (int startyear, int stopyear);

#endif /* _UT_TIME_H_ */
