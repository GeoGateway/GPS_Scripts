/*******************************************************************************
MODULE HEADER:
da_io_pp.h
*******************************************************************************/

#ifndef _DA_IO_PP_H_
#define _DA_IO_PP_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_io_pp_hdr_rcsid[] = "$Id: da_io_pp.h,v 1.8 1999/07/22 16:52:57 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_io_pp.h,v $
 * Revision 1.8  1999/07/22 16:52:57  granat
 * added prototypes for double and unsigned char versions of write_col_cascade()
 *
 * Revision 1.7  1999/07/14 22:12:07  granat
 * fixed #ifdef statement
 *
 * Revision 1.6  1997/09/10 14:56:40  granat
 * changed prototypes to reflect major revisions in the module
 *
 * Revision 1.5  1997/01/29 21:30:47  agray
 * new format.
 *
 * Revision 1.4  1996/10/31 02:15:43  agray
 * renamed from "da_data_pp" to "da_io_pp";
 * changed .h and .c formats throughout library;
 *
 * Revision 1.3  1996/09/23 22:43:19  agray
 * minor changes.
 *
 * Revision 1.2  1996/07/19 17:56:21  agray
 * added write_col_cascade_pp(), write_icol_cascade_pp(),
 * write_col_channel_cascade_pp(), write_icol_channel_cascade_pp()
 *
 * Revision 1.1  1996/07/16 00:43:46  agray
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

/*==============================================================================
Function Declarations
==============================================================================*/

/* input */

int split_data_pp( int *startrow, int *my_numrows, int total_numrows, int pe, 
                   int numPE);

int total_data_pp(int numrows, int *total_numrows, int pe, int numPE);

int read_mat_channel_cascade( char *datafile, int numrows, int numcols, 
                              float ***data, int startrow, int channels,
                              int pe, int numPE );

/* output */

int write_col_cascade_pp( char *outfile, int dim, float *v, char *mode, int pe,
                          int numPE);
int write_icol_cascade_pp( char *outfile, int dim, int *v, char *mode, int pe,
                           int numPE);
int write_dcol_cascade_pp( char *outfile, int dim, double *v, char *mode, 
		           int pe, int numPE);
int write_ccol_cascade_pp( char *outfile, int dim, unsigned char *v, char *mode,
                           int pe, int numPE);
int write_col_channel_cascade_pp( char *outfile, int dim, float *v, 
                                  int num_channels, char *mode, int pe, 
                                  int numPE);
int write_icol_channel_cascade_pp( char *outfile, int dim, int *v, 
                                   int num_channels, char *mode, int pe, 
                                   int numPE);

#endif /* _DA_IO_PP_H_ */
