/* cp_string.h */

#ifndef CP_STRING_HDR
#define CP_STRING_HDR

#ifndef lint
static char cp_string_hdr_rcsid[] = "$Id: cp_string.h,v 1.3 1996/05/01 22:34:43 agray Exp $";
#endif
/* $Log: cp_string.h,v $
 * Revision 1.3  1996/05/01 22:34:43  agray
 * removed "#ifndef strdup"; same with strcasecmp
 *
 * Revision 1.2  1996/01/30 01:02:05  roden
 * no change found as of today.
 * JR, 29 Jan 1996
 *
 * Revision 1.1  1995/08/03 23:38:08  jctran
 * Initial revision
 * */


/* constants */

/* typedefs */

/* macros */

/* function prototypes */
#ifdef __STDC__

#if (CP_STRDUP_DEFINED == CP_STRDUP_IN_CP_STRING_H)
char *strdup(char *s);
#endif

#endif

#if (CP_STRDUP_DEFINED == CP_STRDUP_IN_CP_STRING_H)
char *strdup();
#endif


#endif
