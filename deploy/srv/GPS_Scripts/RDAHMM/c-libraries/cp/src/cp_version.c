/* cp_version.c */

#ifndef lint
static char rcsid[] = "$Id: cp_version.c,v 1.2 1995/08/03 23:38:08 jctran Exp $";
#endif
/* $Log: cp_version.c,v $
 * Revision 1.2  1995/08/03 23:38:08  jctran
 * *** See cp_version.h to see the version update logs ***
 * JCT
 *
 * Revision 1.1  1995/08/03  23:26:10  jctran
 * Initial revision
 * */

/**************************************************************************
 * This is a module that generate the version log for this library
 * such that it complies with our standards for use under rcs.  The general
 * purpose of this is to keep an updatable log that allow the users to know
 * the specifics for the checked out version.
 **************************************************************************/

#include <stdio.h>
#include "cp_version.h"



int cpPrintVersionLog()
   {
   int i=0;
   int j=0;

   printf("Version %s\n", cpVerNum);
   printf("-----------------------\n");
   printf("%s", cpVerLog);
   printf("\n");
   return(0);
   } /* end cpPrintVersionLog */



