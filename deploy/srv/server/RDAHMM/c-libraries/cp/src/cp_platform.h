/* cp_platform.h */

#ifndef CP_PLATFORM_HDR
#define CP_PLATFORM_HDR

#ifndef lint
static char cp_platform_hdr_rcsid[] = "$Id";
#endif
/* $Log: cp_platform.h,v $
 * Revision 1.2  1995/08/23 18:43:31  jctran
 * Cleaned up formating.
 * JCT
 *
 * Revision 1.1  1995/08/03  23:38:08  jctran
 * Initial revision
 * */

/* defining table for the different supported systems */
#define CP_UNIX									0
#define CP_DOS                            1
#define CP_MACOS                          2

/* CHANGED BY USER FOR THE SPECIFIC PLATFORM DESIARED */
/* platform switch */
#define CP_CURR_PLATFORM                  CP_UNIX

/*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*CHANGEME*/


/* values for decision of each issues */
#define CP_CMD_LINE_ARGS_ZERO					0
#define CP_CMD_LINE_ARGS_ONE					1

#define CP_STRDUP_IN_CP_STRING_H				0
#define CP_STRDUP_IN_STRING_H					1


/* default platform ( CP_UNIX ) */
	/* issue 1 */
#define CP_NUM_OF_CMD_LINE_ARGS				CP_CMD_LINE_ARGS_ONE
	/* issue 2 */
#define CP_STRDUP_DEFINED						CP_STRDUP_IN_STRING_H


/* platform flag setting */
#if CP_CURR_PLATFORM == CP_DOS
#elif CP_CURR_PLATFORM == CP_MACOS
	/* issue 1 */
	#undef CP_NUM_OF_CMD_LINE_ARGS
	#define CP_NUM_OF_CMD_LINE_ARGS			CP_CMD_LINE_ARGS_ZERO
	/* issue 2 */
	#undef CP_STRDUP_DEFINED
	#define CP_STRDUP_DEFINED					CP_STRDUP_IN_CP_STRING_H
#endif


#endif

