/* cp_version.h */

#ifndef CP_VERSION_HDR
#define CP_VERSION_HDR

static char cpVerNum[] = "$Id: cp_version.h,v 1.5 1997/01/30 02:46:08 agray Exp $";
static char cpVerLog[] = " log file ...";
/* 
 * $Log: cp_version.h,v $
 * Revision 1.5  1997/01/30 02:46:08  agray
 * made strings into static char's to fix a compile problem.
 *
 * Revision 1.4  1996/09/27 17:39:36  agray
 * commented out k&r prototypes so this library can be linked with c++ code.
 *
 * Revision 1.3  1995/09/06 00:05:11  jctran
 * Update for compiling with cc.
 * JCT
 *
 * Revision 1.2  1995/08/03  23:38:08  jctran
 * Added and modified the header comments, locations, and names for
 * all the .h and .c files.  Also, libcp will now support multi-platforms.
 * Changes can be made in cp_platform.h for the desired platforms.
 * Also, prototypes were added where needed.
 * JCT
 *
 * Revision 1.1  1995/08/03  23:26:10  jctran
 * Initial revision
 *";
 */
#define CP_MAX_COL                        79


/* prototypes */
#ifdef __STDC__
int cpPrintVersionLog();
#endif

/*
int cpPrintVersionLog();
*/


#endif
