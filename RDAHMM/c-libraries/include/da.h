/*******************************************************************************
BLANKET HEADER:
da.h
*******************************************************************************/

#ifndef _DA_H_
#define _DA_H_
/* Protects from multiple inclusion. */

#ifndef lint
static char da_h_rcsid[] = "$Id: da.h,v 1.9 1998/06/29 22:12:37 granat Exp $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da.h,v $
 * Revision 1.9  1998/06/29 22:12:37  granat
 * added da_util
 *
 * Revision 1.8  1997/07/29 03:27:46  agray
 * added da_memory
 *
 * Revision 1.7  1997/06/20 20:46:05  granat
 * added header files for da_nrhacks, da_optim, da_geom
 *
 * Revision 1.6  1997/01/29 21:28:48  agray
 * new format.
 *
 * Revision 1.5  1996/10/31 02:17:16  agray
 * updated due to library overhaul involving renaming and reorganization.
 *
 * Revision 1.4  1996/09/23 22:50:13  agray
 * added new module.
 *
 * Revision 1.3  1996/07/11 16:35:51  agray
 * alphabetized include list
 *
 * Revision 1.2  1996/05/29 03:01:58  agray
 * added some headers
 *
 * Revision 1.1  1996/02/29 02:29:06  agray
 * Initial revision
 *
 * */

/* should include this first */
#include "da_platform.h"

#include "da_util.h"
#include "da_nrhacks.h"
#include "da_cluster.h"
#include "da_io.h"
#include "da_linalg.h"
#include "da_memory.h"
#include "da_probstat.h"
#include "da_random.h"
#include "da_signal.h"
#include "da_timeseg.h"
#include "da_optim.h"
#include "da_geom.h"

#endif /* _DA_H_ */
