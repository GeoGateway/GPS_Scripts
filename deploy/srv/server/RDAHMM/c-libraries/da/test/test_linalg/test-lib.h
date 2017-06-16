/* test-lib.h */

#ifndef DA_TEST_LIB_HDR
#define DA_TEST_LIB_HDR

#ifndef lint
static char da_test_lib_hdr_rcsid[] = "$Id: test-lib.h,v 1.1 1996/02/06 03:32:22 agray Exp $";
#endif
/* $Log: test-lib.h,v $
 * Revision 1.1  1996/02/06  03:32:22  agray
 * Initial revision
 * */

/* needed for gaussj() function */
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#endif
