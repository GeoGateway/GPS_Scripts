/*******************************************************************************
MODULE HEADER:
da_comm_pp.h
*******************************************************************************/

#ifndef _DA_COMM_PP_H_
#define _DA_COMM_PP_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_comm_pp_h_rcsid[] = "$Id: da_comm_pp.h,v 1.8 1997/06/25 23:54:36 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_comm_pp.h,v $
 * Revision 1.8  1997/06/25 23:54:36  granat
 * added prototypes for init, finalize, and barrier
 *
 * Revision 1.7  1997/06/24 23:30:21  granat
 * added mpi versions of functions, added recv_bcast_msg()
 *
 * Revision 1.6  1997/01/29 21:29:54  agray
 * new format.
 *
 * Revision 1.5  1996/10/31 02:10:15  agray
 * renamed from "da_msg_pp" to "da_comm_pp";
 * changed .h and .c formats throughout library.
 *
 * Revision 1.4  1996/10/31 01:23:19  agray
 * consolidated functions for different platforms.
 *
 * Revision 1.3  1996/09/12 23:37:00  granat
 * Added constant PVM_GROUP
 * Added prototypes for convert_pe_to_tid and convert_tid_to_pe
 *
 * Revision 1.2  1996/07/16 00:09:47  agray
 * added get_pe_info_cray().
 *
 * Revision 1.1  1996/07/11 18:13:31  agray
 * Initial revision
 *
 * */

/*==============================================================================
Data Structures
==============================================================================*/

/*==============================================================================
Constants, Macros
==============================================================================*/

#define DA_INT    'i'
#define DA_FLOAT  'f'
#define DA_DOUBLE 'd'
#define DA_UCHAR  'c'

#define DA_NO_MSGTAG  0
#define DA_ALL_PES    (char*)NULL

#define PVM_GROUP "meta_group"


/*==============================================================================
Variables
==============================================================================*/

/*==============================================================================
Function Declarations
==============================================================================*/

/* initialize and finalizing parallel programs */

int DA_init_pp( int *argc_p, char ***argv_p );
int DA_finalize_pp();

/* barrier */

int DA_barrier();  

/* identifying the processor */

int DA_get_pe_info(int *pe, int *numPE);

int DA_convert_pe_to_tid(int pe);
int DA_convert_tid_to_pe(int tid);

/* sending and receiving messages */

int DA_send_msg( void *data, int length, int target_pe, int msgtag, 
                 char datatype);
int DA_recv_msg( void *data, int length, int sending_pe, int msgtag, 
                 char datatype);
int DA_broadcast_msg( void *data, int length, char *target_group, int msgtag, 
                      char datatype);
int DA_recv_bcast_msg( void *data, int length, int sending_pe, int msgtag, 
                       char datatype);

#endif /* _DA_COMM_PP_H_ */
