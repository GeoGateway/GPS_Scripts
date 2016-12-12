/*******************************************************************************
MODULE HEADER:
da_platform.h
*******************************************************************************/

#ifndef _DA_PLATFORM_H_
#define _DA_PLATFORM_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_platform_h_rcsid[] = "$Id";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_platform.h,v $
 * Revision 1.4  1998/04/21 17:10:38  roden
 * Added definition of DA_LINUX supported platform and defined flag settings
 * for Linux to indicate no PVM present.  This is as per hyglac and cism
 * clusters.  If DA_LINUX is too general in the future we could define settings
 * on a per-system basis, e.g. DA_HYGLAC, DA_CISM, etc.
 *
 * Revision 1.3  1997/06/25 23:53:18  granat
 * added back some defs needed in da_random, some cosmetic changes
 *
 * Revision 1.2  1997/06/24 23:33:25  granat
 * added some parallel communications defs
 *
 * Revision 1.1  1997/01/29 21:31:49  agray
 * Initial revision
 *
 * */

/* defining table for the different supported systems */
#define DA_UNIX									0
#define DA_DOS									1
#define DA_MACOS								2
#define DA_CRAY								    3
#define DA_LINUX								    4

/* CHANGED BY USER FOR THE SPECIFIC PLATFORM DESIRED */
/* platform switch */  
#define DA_CURR_PLATFORM						DA_UNIX
/*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*/


/* values for decision of each issues */
#define DA_INCL_FILE_IN_SYS                     0
#define DA_INCL_FILE_NOT_IN_SYS                 1
 
#define DA_STRDUP_ALREADY_DEFINED               0
#define DA_DEFINE_STRDUP                        1
 
#define DA_STRCASECMP_ALREADY_DEFINED           0
#define DA_DEFINE_STRCASECMP                    1

#define DA_NO_PVM				0
#define DA_PVM_EXISTS				1

#define DA_NO_MPI				0
#define DA_MPI_EXISTS				1

#define DA_AUTO_SPAWN				0
#define DA_MANUAL_SPAWN				1

#define DA_NO_PE_ALIAS				0
#define DA_PE_ALIAS				1


/* default platform ( UNIX ) */
        /* issue 1 */
#define DA_INCL_FILE_LOC			DA_INCL_FILE_IN_SYS
        /* issue 2 */
#define DA_STRDUP_EXISTS			DA_STRDUP_ALREADY_DEFINED
	/* issue 3 */
#define DA_STRCASECMP_EXISTS			DA_STRCASECMP_ALREADY_DEFINED
        /* issue 4 */
#define DA_MPI					DA_MPI_EXISTS
        /* issue 5 */
#define DA_PVM					DA_NO_PVM
        /* issue 6 */
#define DA_PVM_SPAWN				DA_MANUAL_SPAWN
        /* issue 7 */
#define DA_PVM_TID				DA_NO_PE_ALIAS


/* platform flag setting */
#if DA_CURR_PLATFORM == DA_DOS
#elif DA_CURR_PLATFORM == DA_MACOS
  /* issue 1 */
  #undef DA_INCL_FILE_LOC
  #define DA_INCL_FILE_LOC				DA_INCL_FILE_NOT_IN_SYS
  /* issue 2 */
  #undef DA_STRDUP_EXISTS
  #define DA_STRDUP_EXISTS				DA_DEFINE_STRDUP
  /* issue 3 */
  #undef DA_STRCASECMP_EXISTS
  #define DA_STRCASECMP_EXISTS				DA_DEFINE_STRCASECMP
  /* issue 4 */
  #undef DA_MPI
  #define DA_MPI                                        DA_NO_MPI
  /* issue 5 */
  #undef DA_PVM
  #define DA_PVM                                        DA_NO_PVM

#elif DA_CURR_PLATFORM == DA_CRAY
  /* issue 4 */
  #undef DA_MPI
  #define DA_MPI					DA_MPI_EXISTS
  /* issue 5 */
  #undef DA_PVM
  #define DA_PVM					DA_PVM_EXISTS
  /* issue 6 */
  #undef DA_PVM_SPAWN
  #define DA_PVM_SPAWN					DA_AUTO_SPAWN
  /* issue 7 */
  #undef DA_PVM_TID
  #define DA_PVM_TID					DA_PE_ALIAS

#elif DA_CURR_PLATFORM == DA_LINUX
  /* issue 4 */
  #undef DA_MPI
  #define DA_MPI					DA_MPI_EXISTS
  /* issue 5 */
  #undef DA_PVM
  #define DA_PVM					DA_NO_PVM


#endif  /* #if DA_CURR_PLATFORM == whatever */

#endif  /* #ifndef _DA_PLATFORM_H_ */

