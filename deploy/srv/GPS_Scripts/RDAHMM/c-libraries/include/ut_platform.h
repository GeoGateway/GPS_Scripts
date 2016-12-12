/*******************************************************************************
MODULE HEADER:
ut_platform.h
*******************************************************************************/

#ifndef _UT_PLATFORM_H_
#define _UT_PLATFORM_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char ut_platform_h_rcsid[] = "$Id";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: ut_platform.h,v $
 * Revision 1.1  1997/01/29 23:44:54  agray
 * Initial revision
 *
 * */

/* defining table for the different supported systems */
#define UT_UNIX									0
#define UT_DOS									1
#define UT_MACOS								2

/* CHANGED BY USER FOR THE SPECIFIC PLATFORM DESIRED */
/* platform switch */  
#define UT_CURR_PLATFORM						UT_UNIX
/*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*/


/* values for decision of each issues */
#define UT_INCL_FILE_IN_SYS						0
#define UT_INCL_FILE_NOT_IN_SYS					1

#define UT_STRDUP_ALREADY_DEFINED				0
#define UT_DEFINE_STRDUP						1

#define UT_STRCASECMP_ALREADY_DEFINED			0
#define UT_DEFINE_STRCASECMP					1


/* default platform ( UNIX ) */
	/* issue 1 */
#define UT_INCL_FILE_LOC					UT_INCL_FILE_IN_SYS
	/* issue 2 */
#define UT_STRDUP_EXISTS					UT_STRDUP_ALREADY_DEFINED
	/* issue 3 */
#define UT_STRCASECMP_EXISTS				UT_STRCASECMP_ALREADY_DEFINED


/* platform flag setting */
#if UT_CURR_PLATFORM == UT_DOS
#elif UT_CURR_PLATFORM == UT_MACOS
	/* issue 1 */
	#undef UT_INCL_FILE_LOC
	#define UT_INCL_FILE_LOC					UT_INCL_FILE_NOT_IN_SYS
	/* issue 2 */
	#undef UT_STRDUP_EXISTS
	#define UT_STRDUP_EXISTS					UT_DEFINE_STRDUP
	/* issue 3 */
	#undef UT_STRCASECMP_EXISTS
	#define UT_STRCASECMP_EXISTS				UT_DEFINE_STRCASECMP
#endif


#endif  /* _UT_PLATFORM_H_ */

