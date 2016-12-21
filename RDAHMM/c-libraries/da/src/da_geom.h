/*******************************************************************************
MODULE HEADER:
da_geom.h
*******************************************************************************/

#ifndef _DA_GEOM_H_
#define _DA_GEOM_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_geom_h_rcsid[] = "$Id: da_geom.h,v 1.1 1997/06/20 20:55:09 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_geom.h,v $
 * Revision 1.1  1997/06/20 20:55:09  granat
 * Initial revision
 *
 * Revision 1.1  1997/05/15 16:46:53  granat
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/
#define DA_TRIG_TABLE_SIZE 256
#define DA_TRIG_SAMPLE_ANG (2 * M_PI / DA_TRIG_TABLE_SIZE )
#define DA_INV_TRIG_SAMPLE_ANG (1 / DA_TRIG_SAMPLE_ANG)

/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/
double DA_table_sin( double angle );
double DA_table_cos( double angle );
double DA_table_tan( double angle );

#endif /* _DA_GEOM_H_ */
