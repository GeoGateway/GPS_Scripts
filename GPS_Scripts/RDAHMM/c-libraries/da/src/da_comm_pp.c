/*******************************************************************************
MODULE NAME
da_comm_pp

ONE-LINE SYNOPSIS
Functions related to communication on parallel machines.

SCOPE OF THIS MODULE
This module is concerned only with a substrate of functions for parallel compu-
ting, not functions which actually perform parallel computations.  

SEE ALSO
Functions which perform parallel computations are contained in other modules, 
denoted by the "_pp" suffix.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
1. /proj/cooltools/sc_pp, AG.
2. /proj/cooltools/kmeans_pp, RG.

NOTES
-

AG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_comm_pp.c,v 1.13 1997/09/10 14:51:36 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_comm_pp.c,v $
 * Revision 1.13  1997/09/10 14:51:36  granat
 * fixed some error handling
 *
 * Revision 1.12  1997/06/27 16:16:27  granat
 * added dummy functions for non-parallel case
 *
 * Revision 1.11  1997/06/26 15:19:37  granat
 * cosmetic changes
 *
 * Revision 1.10  1997/06/25 23:53:55  granat
 * fixed error handling, added init, finalize, and barrier
 * fixed platform handling
 *
 * Revision 1.9  1997/06/24 23:29:49  granat
 * added mpi versions of functions, added recv_bcast_msg()
 *
 * Revision 1.8  1997/06/19 14:18:47  granat
 * Fixed bug in Cray version of broadcast_msg()
 *
 * Revision 1.7  1997/01/29 21:44:36  agray
 * new formatting, cleaning up debugging output using ut_output.
 *
 * Revision 1.6  1996/10/31 02:10:15  agray
 * renamed from "da_msg_pp" to "da_comm_pp";
 * changed .h and .c formats throughout library.
 *
 * Revision 1.5  1996/10/31 01:23:51  agray
 * consolidated functions for different platforms.
 *
 * Revision 1.4  1996/09/12 23:37:53  granat
 * Added functions convert_pe_to_tid and convert_tid_to_pe
 * Modified send_msg and receive_msg to work on the suns
 *
 * Revision 1.3  1996/07/17 20:44:52  agray
 * cosmetic.
 *
 * Revision 1.2  1996/07/16 00:09:33  agray
 * added get_pe_info_cray().
 *
 * Revision 1.1  1996/07/11 18:13:17  agray
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* DA library */
#include "da_platform.h"

/* MPI library */
#if (DA_MPI == DA_MPI_EXISTS)
#include "mpi.h"
#endif

/* PVM library */
#if (DA_PVM == DA_PVM_EXISTS)
#include "pvm3.h"
#endif

/* UT library */
#include "ut_types.h"
#include "ut_error.h"
#include "ut_output.h"

/* this module's header */
#include "da_comm_pp.h"


/*******************************************************************************
DA_INIT_PP
Initialize a parallel program.  This should be the first thing done in main.
The inputs should be the addresses of argc and argv that are the arguments to 
main.  This function should not be called more than once.
RG
*******************************************************************************/
#if (DA_MPI == DA_MPI_EXISTS)

int DA_init_pp( int *argc_p, char ***argv_p )
{
  int cc;

  cc = MPI_Init( argc_p, argv_p );

  if (cc == MPI_SUCCESS)
    return( UT_OK );
  else {
    err_printf();
    log_printf( "MPI_Init failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
}

#elif (DA_PVM == DA_PVM_EXISTS)
#if (DA_PVM_SPAWN == DA_AUTO_SPAWN)

int DA_init_pp( int *argc_p, char ***argv_p )
{
  /* In this case there is no need to do any real initialization */

  return( UT_OK );
}

#elif (DA_PVM_SPAWN == DA_MANUAL_SPAWN)

int DA_init_pp( int *argc_p, char ***argv_p )
{
  int  *tids;
  int  id;
  int  nprocs;
  int  cc;
  int  i, n;
  char *cp;

 /****
  * First parse the number of processes from the command line arguments. 
  * This is done here to avoid forcing the user to explicitly parse 
  * parallel arguments.
  ****/

  nprocs = 1;

  for (i = 1; i <= (*argc_p); i++)
    if ( !strcmp( (*argv_p)[i], "-npes" ) ) {
      cp = 0;
      n = strtol( (*argv_p)[i + 1], &cp, 10 );
      if ( cp == (*argv_p)[i + 1] ) {
        err_printf();
        log_printf( "Expecting an integer arg to -npes, not \"%s\"\n", 
                    (*argv_p)[i + 1] );
        return( UT_ERROR );
      }
      else if( cp && *cp ) {
        err_printf();
        log_printf( "Ignoring extra args to -npes\n" );
        return( UT_ERROR );
      }
      else if( n <= 0 ) {
        err_printf();
        log_printf( "-npes must be > 0, not %d\n", n );
        return( UT_ERROR );
      }
      else
        nprocs = n;

      memmove( *argv_p + i, *argv_p + i + 2, 
               (*argc_p - i - 1) * sizeof( char * ) ); 
      *argc_p -= 2;
    }

  /* now identify this process and join the group */

  id = pvm_joingroup( PVM_GROUP );
  if (id < PvmOk) {
    printf( "pvm_joingroup failed with error #%d.\n", id );
    return( UT_ERROR );
  }

  /* spawn the processes if this is the first process */

  if ( (id == 0) && (nprocs > 1) ) {
    tids = malloc( (nprocs - 1) * sizeof( int ) );
    cc = pvm_spawn( (*argv_p)[0], &(*argv_p)[1], PvmTaskDefault, 0, nprocs - 1,
                    tids );
    free( tids );
    if (cc < PvmOk) {
      err_printf();
      log_printf( "Error in pvm_spawn (%d)\n", cc );
      return( UT_ERROR );
    }
    else if (cc < (nprocs - 1)) {
      err_printf();
      log_printf( "Insufficient processes (%d) spawned\n", cc );
      return( UT_ERROR );
    }
    else if (cc > (nprocs - 1)) {
      err_printf();
      log_printf( "Too many processes (%d) spawned\n", cc );
      return( UT_ERROR );
    }
  }

  return( UT_OK );
} 

#endif
#else /* No parallel library available */

int DA_init_pp( int *argc_p, char ***argv_p )
{
  return( UT_OK );
}

#endif

/*******************************************************************************
DA_FINALIZE_PP
Clean up after a parallel program.  This should be the last function called
before main returns.
RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_finalize_pp()
{
  int cc;

  cc = MPI_Finalize();

  if (cc == MPI_SUCCESS)
    return( UT_OK );
  else {
    err_printf();
    log_printf( "MPI_finalize failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
}

#elif (DA_PVM == DA_PVM_EXISTS)

int DA_finalize_pp()
{
  int cc;

  cc = pvm_exit();

  if (cc < PvmOk) {
    err_printf();
    log_printf( "pvm_exit failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
  else
    return( UT_OK );
}

#else /* No parallel library available */

int DA_finalize_pp()
{
  return( UT_OK );
}

#endif

/*******************************************************************************
DA_BARRIER
Synchronize parallel processes.
RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_barrier()
{
  int cc;

  cc = MPI_Barrier( MPI_COMM_WORLD );

  if (cc == MPI_SUCCESS)
    return( UT_OK );
  else {
    err_printf();
    log_printf( "MPI_Barrier failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
}

#elif (DA_PVM == DA_PVM_EXISTS)

int DA_barrier()
{
  int cc;
  int nprocs;

  nprocs = pvm_gsize( PVM_GROUP );

  if (nprocs < PvmOk) {
    err_printf();
    log_printf( "pvm_gsize failed with error #%d.\n", nprocs );
    return( UT_ERROR );
  }

  cc = pvm_barrier( PVM_GROUP, nprocs );

  if (cc < PvmOk) {
    err_printf();
    log_printf( "pvm_barrier failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
  else
    return( UT_OK );
} 

#else /* No parallel library available */

int DA_barrier()
{
  return( UT_OK );
}

#endif
    
/*******************************************************************************
DA_GET_PE_INFO
Get the virtual process number and the number of processes
AG, RG
*******************************************************************************/
#if (DA_MPI == DA_MPI_EXISTS)

int DA_get_pe_info( int *pe, int *numPE )
{
  int cc;

  cc = MPI_Comm_rank( MPI_COMM_WORLD, pe); 

  if (cc != MPI_SUCCESS) {
    err_printf();
    log_printf( "MPI_Comm_rank failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
 
  cc = MPI_Comm_size( MPI_COMM_WORLD, numPE); 

  if (cc == MPI_SUCCESS)
    return( UT_OK );
  else {
    err_printf();
    log_printf( "MPI_Comm_size failed with error #%d.\n", cc );
    return( UT_ERROR );
  }
}

#elif (DA_PVM == DA_PVM_EXISTS)
#if (DA_PVM_TID == DA_PE_ALIAS)

int DA_get_pe_info( int *pe, int *numPE )
{
  int pvm_tid;

  pvm_tid = pvm_mytid();     /*Find out which process number*/

  if (pvm_tid < PvmOk) {
    err_printf();
    log_printf( "pvm_mytid failed with error #%d.\n", pmv_tid );
    return( UT_ERROR );
  }

  *pe = pvm_get_PE( pvm_tid ); /*Get the virtual PE number*/

  if (*pe < PvmOk) {
    err_printf();
    log_printf( "pvm_get_PE faled with error #%d.\n", *pe );
    return( UT_ERROR );
  }

  *numPE = pvm_gsize( 0 );     /*Do the head count*/

  if (*numPE < PvmOk) {
    err_printf();
    log_printf( "pvm_gsize failed with error #%d.\n", *numPE );
    return( UT_ERROR );
  } 
  else
    return(UT_OK);
}

#elif (DA_PVM_TID == DA_NO_PE_ALIAS)

int DA_get_pe_info( int *pe, int *numPE )
{
  int pvm_tid;

  pvm_tid = pvm_mytid();

  if (pvm_tid < PvmOk) {
    err_printf();
    log_printf( "pvm_mytid failed with error #%d.\n", pvm_tid );
    return( UT_ERROR );
  }

  *pe = pvm_getinst( PVM_GROUP, pvm_tid );
 
  if (*pe < PvmOk) {
    err_printf();
    log_printf( "pvm_getinst failed with error #%d.\n", *pe );
    return( UT_ERROR );
  }

  *numPE = pvm_gsize( PVM_GROUP );     /*Do the head count*/

  if (*numPE < PvmOk) {
    err_printf();
    log_printf( "pvm_gsize failed with error #%d.\n", *numPE );
    return( UT_ERROR ); 
  }
  else
    return(UT_OK);
}

#endif
#else /* No parallel library available */

int DA_get_pe_info( int *pe, int *numPE )
{
  *pe = 0;
  *numPE = 1;

  return( UT_OK );
}

#endif

/*******************************************************************************
DA_CONVERT_PE_TO_TID
Convert the processing element number to the arbitrary pvm tid number of that
process.
AG, RG
*******************************************************************************/

#if ((DA_MPI == DA_NO_MPI) && (DA_PVM == DA_PVM_EXISTS))

int DA_convert_pe_to_tid( int pe )
{
  int tid = pvm_gettid( PVM_GROUP, pe );

  if (tid < PvmOk) {
    err_printf();
    log_printf( "pvm_gettid failed with error #%d.\n", tid );
    return( UT_ERROR );
  }

  return( tid );
}

#endif

/*******************************************************************************
DA_CONVERT_TID_TO_PE
Convert the arbitrary tid number of a process to the processing element number.
AG, RG
*******************************************************************************/

#if ((DA_MPI == DA_NO_MPI) && (DA_PVM == DA_PVM_EXISTS))

int DA_convert_tid_to_pe ( int tid )
{
  int pe = pvm_getinst( PVM_GROUP, tid );

  if (pe < PvmOk) {
    err_printf();
    log_printf( "pvm_getinst failed with error #%d.\n", pe );
    return( UT_ERROR );
  }

  return( pe );
}

#endif

/*******************************************************************************
DA_SEND_MSG
Send a message to another processor.
AG, RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_send_msg( void *data, int length, int target_pe, int msgtag, 
                 char datatype)
 
{
  int cc;
 
  switch(datatype)
  {
    case DA_INT:
      cc =
       MPI_Send(data, length, MPI_INT, target_pe, msgtag, MPI_COMM_WORLD);
      break;
    case DA_FLOAT:
      cc =
        MPI_Send(data, length, MPI_FLOAT, target_pe, msgtag, MPI_COMM_WORLD);
      break;
    case DA_DOUBLE:
      cc =
        MPI_Send(data, length, MPI_DOUBLE, target_pe, msgtag, MPI_COMM_WORLD);
      break;
    case DA_UCHAR:
      cc = MPI_Send(data, length, MPI_UNSIGNED_CHAR, target_pe, msgtag, 
	       	    MPI_COMM_WORLD);
      break;
    default:
      err_printf();
      log_printf("bad data type passed to send_msg.\n");
      return (UT_ERROR);
  }
  if (cc != MPI_SUCCESS) {
    err_printf();
    log_printf("mpi_send failed with error #%d.\n",cc);
    return( UT_ERROR );
  }
  else
    return( UT_OK );
}

#elif (DA_PVM == DA_PVM_EXISTS)
#if (DA_PVM_TID == DA_PE_ALIAS)

int DA_send_msg( void *data, int length, int target_pe, int msgtag, 
                 char datatype)

{
  int cc;
  int stride = 1;

  cc = pvm_initsend(PvmDataRaw);
  if (cc < PvmOk)
  {
    err_printf();
    log_printf("pvm_initsend failed with error #%d.\n",cc);
    return cc;
  }
  else
  {
    switch(datatype)
    {
      case DA_INT:
        cc = pvm_pkint((int*)data,length,stride);
        break;
      case DA_FLOAT:
        cc = pvm_pkfloat((float*)data,length,stride);
        break;
      case DA_DOUBLE:
        cc = pvm_pkdouble((double*)data,length,stride);
        break;
      default:
        err_printf();
        log_printf("bad data type passed to send_msg.\n");
        return (UT_ERROR); 
        /* should make pvm error codes correspond with UT */
    }
    if (cc < PvmOk)
    {
      err_printf();
      log_printf("pvm packing function failed with error #%d.\n",cc);
      return cc;
    }
    else
    {
      cc = pvm_send(target_pe,msgtag);
      if (cc < PvmOk)
      {
        err_printf();
        log_printf("pvm_send failed with error #%d.\n",cc);
        return cc;
      }
    }
  }
}

#elif (DA_PVM_TID == DA_NO_PE_ALIAS)

int DA_send_msg( void *data, int length, int target_pe, int msgtag, 
                 char datatype)
 
{
  int cc;
  int stride = 1;
 
  cc = pvm_initsend(PvmDataRaw);
  if (cc < PvmOk)
  {
    err_printf();
    log_printf("pvm_initsend failed with error #%d.\n",cc);
    return cc;
  }
  else
  {
    switch(datatype)
    {
      case DA_INT:
        cc = pvm_pkint((int*)data,length,stride);
        break;
      case DA_FLOAT:
        cc = pvm_pkfloat((float*)data,length,stride);
        break;
      case DA_DOUBLE:
        cc = pvm_pkdouble((double*)data,length,stride);
        break;
      default:
        log_printf("bad data type passed to send_msg.\n");
        return (UT_ERROR);
        /* should make pvm error codes correspond with UT */
    }
    if (cc < PvmOk)
    {
      err_printf();
      log_printf("pvm packing function failed with error #%d.\n",cc);
      return cc;
    }
    else
    {
      cc = pvm_send(convert_pe_to_tid(target_pe),msgtag);
      if (cc < PvmOk)
      {
        err_printf();
        log_printf("pvm_send failed with error #%d.\n",cc);
        return cc;
      }
    }
  }
}

#endif
#else /* No parallel library available */

int DA_send_msg( void *data, int length, int target_pe, int msgtag,
                 char datatype)
{
  return( UT_OK );
}

#endif

/*******************************************************************************
DA_RECV_MSG
Receive a message from another processor.
AG, RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_recv_msg( void *data, int length, int sending_pe, int msgtag,
                 char datatype)
 
{
  int cc;
  MPI_Status status;
 
  switch(datatype)
  {
    case DA_INT:
      cc = MPI_Recv( data, length, MPI_INT, sending_pe, msgtag,
                     MPI_COMM_WORLD, &status);
      break;
    case DA_FLOAT:
      cc = MPI_Recv( data, length, MPI_FLOAT, sending_pe, msgtag,
                     MPI_COMM_WORLD, &status);
      break;
    case DA_DOUBLE:
      cc = MPI_Recv( data, length, MPI_DOUBLE, sending_pe, msgtag,
                     MPI_COMM_WORLD, &status);
      break;
    case DA_UCHAR:
      cc = MPI_Recv( data, length, MPI_UNSIGNED_CHAR, sending_pe, msgtag,
                     MPI_COMM_WORLD, &status);
      break;
    default:
      err_printf();
      log_printf("bad data type passed to recv_msg.\n");
      return (UT_ERROR); /* should make pvm error codes correspond with UT */
  }
  if (cc != MPI_SUCCESS) {
    err_printf();
    log_printf("mpi_recv failed with error #%d.\n",cc);
    return UT_ERROR;
  }
  else
    return UT_OK;
}

#elif (DA_PVM == DA_PVM_EXISTS)
#if (DA_PVM_TID == DA_PE_ALIAS)

int DA_recv_msg( void *data, int length, int sending_pe, int msgtag, 
                 char datatype)

{
  int cc, buf, stride = 1;

  buf = pvm_recv(sending_pe,msgtag);
  switch(datatype)
  {
    case DA_INT:
      cc = pvm_upkint((int*)data,length,stride);
      break;
    case DA_FLOAT:
      cc = pvm_upkfloat((float*)data,length,stride);
      break;
    case DA_DOUBLE:
      cc = pvm_upkdouble((double*)data,length,stride);
      break;
    default:
      err_printf();
      log_printf("bad data type passed to recv_msg.\n");
      return (UT_ERROR); /* should make pvm error codes correspond with UT */
  }
  if (cc < PvmOk)
  {
    err_printf();
    log_printf("pvm unpacking function failed with error #%d.\n",cc);
    return cc;
  }
  else
  {
    cc = pvm_freebuf(buf);
    if (cc < PvmOk)
    {
      err_printf();
      log_printf("pvm_freebuf failed with error #%d.\n",cc);
      return cc;
    }
  }
}

#elif (DA_PVM_TID == DA_NO_PE_ALIAS)

int DA_recv_msg( void *data, int length, int sending_pe, int msgtag,
                 char datatype)
 
{
  int cc, buf, stride = 1;
 
  buf = pvm_recv(convert_pe_to_tid(sending_pe),msgtag);
  switch(datatype)
  {
    case DA_INT:
      cc = pvm_upkint((int*)data,length,stride);
      break;
    case DA_FLOAT:
      cc = pvm_upkfloat((float*)data,length,stride);
      break;
    case DA_DOUBLE:
      cc = pvm_upkdouble((double*)data,length,stride);
      break;
    default:
      err_printf();
      log_printf("bad data type passed to recv_msg.\n");
      return (UT_ERROR); /* should make pvm error codes correspond with UT */
  }
  if (cc < PvmOk)
  {
    err_printf();
    log_printf("pvm unpacking function failed with error #%d.\n",cc);
    return cc;
  }
  else
  {
    cc = pvm_freebuf(buf);
    if (cc < PvmOk)
    {
      err_printf();
      log_printf("pvm_freebuf failed with error #%d.\n",cc);
      return cc;
    }
  }
}

#endif
#else /* No parallel library available */

int DA_recv_msg( void *data, int length, int sending_pe, int msgtag,
                 char datatype)
{
  return( UT_OK );
}

#endif

/*******************************************************************************
DA_BROADCAST_MSG
Broadcast a message to a group of processors.
AG, RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_broadcast_msg( void *data, int length, char *target_group, int root,
                      char datatype)
{
  int cc;

  switch(datatype)
  {
    case DA_INT:
      cc = MPI_Bcast( data, length, MPI_INT, root, MPI_COMM_WORLD );
      break;
    case DA_FLOAT:
      cc = MPI_Bcast( data, length, MPI_FLOAT, root, MPI_COMM_WORLD );
      break;
    case DA_DOUBLE:
      cc = MPI_Bcast( data, length, MPI_DOUBLE, root, MPI_COMM_WORLD );
      break;
    case DA_UCHAR:
      cc = MPI_Bcast( data, length, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD );
      break;
    default:
      log_printf("bad data type passed to broadcast_msg.\n");
      return (UT_ERROR);
  }
  if (cc != MPI_SUCCESS) {
    err_printf();
    log_printf("mpi_bcast failed with error #%d.\n",cc);
    return UT_ERROR;
  }
  else
    return UT_OK;
}

#elif (DA_PVM == DA_PVM_EXISTS)

int DA_broadcast_msg( void *data, int length, char *target_group, int root, 
                      char datatype)

{
  int cc, stride = 1;
  int pe, numPE;
  int msgtag;
  struct timeval tp;
  
  gettimeofday(&tp,0);
  msgtag = (int) tp.tv_usec;

  DA_get_pid_info(&pe, &numPE);

  if (pe == root) {
    cc = pvm_initsend(PvmDataRaw);

    if (cc < PvmOk) {
      err_printf();
      log_printf( "pvm_initsend failed with error #%d.\n", cc );
      return( UT_ERROR );
    }

    switch(datatype)
    {
      case DA_INT:
        cc = pvm_pkint((int*)data,length,stride);
        break;
      case DA_FLOAT:
        cc = pvm_pkfloat((float*)data,length,stride);
        break;
      case DA_DOUBLE:
        cc = pvm_pkdouble((double*)data,length,stride);
        break;
      default:
        err_printf();
        log_printf("bad data type passed to broadcast_msg.\n");
        return (UT_ERROR); /* should make pvm error codes correspond with UT */
    }
    if (cc < PvmOk)
    {
      err_printf();
      log_printf("pvm packing function failed with error #%d.\n",cc);
      return( UT_ERROR );
    }
    else
    {
      cc = pvm_bcast(target_group,msgtag);
      if (cc < PvmOk)
      {
        err_printf();
        log_printf("pvm_bcast failed with error #%d.\n",cc);
        return( UT_ERROR );
      }
    }
  }
  else
    DA_recv_msg(data, length, root, msgtag, datatype);

  return( UT_OK );
}

#else /* No parallel library available */

int DA_broadcast_msg( void *data, int length, char *target_group, int msgtag,
                      char datatype)
{
  return( UT_OK );
}

#endif

/*******************************************************************************
DA_RECV_BCAST_MSG
Receive a message that was broadcast to a group of processors.
RG
*******************************************************************************/

#if (DA_MPI == DA_MPI_EXISTS)

int DA_recv_bcast_msg( void *data, int length, int sending_pe, int msgtag,
                       char datatype )
{
  int cc;
 
  switch(datatype)
  {
    case DA_INT:
      cc = MPI_Bcast( data, length, MPI_INT, sending_pe, MPI_COMM_WORLD );
      break;
    case DA_FLOAT:
      cc = MPI_Bcast( data, length, MPI_FLOAT, sending_pe, MPI_COMM_WORLD );
      break;
    case DA_DOUBLE:
      cc = MPI_Bcast( data, length, MPI_DOUBLE, sending_pe, MPI_COMM_WORLD );
      break;
    case DA_UCHAR:
      cc = MPI_Bcast( data, length, MPI_UNSIGNED_CHAR, sending_pe, MPI_COMM_WORLD );
      break;
    default:
      err_printf();
      log_printf("bad data type passed to broadcast_msg.\n");
      return( UT_ERROR ); /* should make pvm error codes correspond with UT */
  }
  if (cc != MPI_SUCCESS) {    
    err_printf();
    log_printf("mpi_bcast failed with error #%d.\n",cc);
    return( UT_ERROR );
  }
  else
    return( UT_OK );
}

#elif (DA_PVM == DA_PVM_EXISTS)

int DA_recv_bcast_msg( void *data, int length, int sending_pe, int msgtag,
                       char datatype )
{
  return( DA_recv_msg( data, length, sending_pe, msgtag, datatype ) );
}

#else /* No parallel library available */

int DA_recv_bcast_msg( void *data, int length, int sending_pe, int msgtag,
                       char datatype )
{
  return( UT_OK );
}

#endif
