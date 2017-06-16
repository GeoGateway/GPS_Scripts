# Makefile for Data Analysis library
#
# $Id: Makefile,v 1.23 1998/07/02 01:06:42 granat Exp granat $
#
# $Log: Makefile,v $
# Revision 1.23  1998/07/02 01:06:42  granat
# checking in temporary hack by agray
#
# Revision 1.22  1998/04/20 22:42:24  roden
# - Simplified $(LIBRARY) rule so that Gnu make doesn't call ar unnecessarily.
# - Removed %.o rule all together since implicit rule was being used mostly,
#   and explicit rule caused two compilations under certain circumstances,
#   not to mention $(OUTPUT_OPTION) was screwing things up under Linux.
#
# Revision 1.21  1998/04/20 22:41:37  roden
# *** empty log message ***
#
# Revision 1.20  1997/07/31 21:59:06  agray
# added da_memory, alphabetized list of sources and headers, added makedepend
#
# Revision 1.19  1997/06/26 15:09:33  granat
# took out special mpi includes
#
# Revision 1.18  1997/06/26 00:06:35  granat
# added mpi includes, change later when mpich gets set up properly
#
# Revision 1.17  1997/06/20 20:39:27  granat
# added da_nrhacks, da_optim, da_geom
#
# Revision 1.16  1997/05/13 23:50:26  agray
# tons of changes.
#
# Revision 1.15  1996/10/31 02:47:40  agray
# updated after overhaul of library, which included renaming of modules.
#
# Revision 1.14  1996/09/23 23:58:03  agray
# added PVM properly, added -Wall to standard debug compile options.
#
# Revision 1.13  1996/09/23 23:08:29  granat
# Added inclusion of pvm include directory
#
# Revision 1.12  1996/09/23 23:02:23  agray
# un-commented parallel modules.
#
# Revision 1.11  1996/08/28 20:37:16  agray
# udpated list of modules - .c's and .h's.
#
# Revision 1.10  1996/07/15 18:18:54  agray
# commented out parallel modules for now, until PVM works on Suns.
#
# Revision 1.9  1996/07/11 16:46:51  agray
# added parallel modules da_clust_pp and da_msg_pp.
#
# Revision 1.8  1996/05/23 23:23:03  granat
# cosmetic
# rg
#
# Revision 1.6  1996/05/06 23:18:49  agray
# updated for solaris port; added standard directory variables
#
# Revision 1.5  1996/04/30 21:34:44  agray
# added da_dist.c,h.
# ag
#
# Revision 1.4  1996/02/28 04:54:49  agray
# added da_data module
# ag
#
# Revision 1.3  1996/02/21 05:17:43  agray
# added da_prob.c,.h.
# ag
#
# Revision 1.2  1996/02/21  00:37:57  agray
# changed da.h to da_linalg.h
# ag
#
# Revision 1.1  1996/02/06  03:28:36  agray
# Initial revision
#

# Library, source code, and objects.

LIBRARY		= libda.a

DESTDIR		= ../../lib

INCDIR		= ../../include

SOURCES.c	= \
			da_cluster.c \
			da_geom.c \
			da_io.c \
			da_linalg.c \
			da_memory.c \
			da_nrhacks.c \
			da_optim.c \
			da_probstat.c \
			da_random.c \
			da_signal.c \
			da_timeseg.c \
			da_util.c
#			da_cluster_pp.c \
#			da_comm_pp.c \
#			da_io_pp.c \
#			da_signal_pp.c

SOURCES.h = \
			da_cluster.h \
			da_geom.h \
			da_io.h \
			da_linalg.h \
			da_memory.h \
			da_nrhacks.h \
			da_optim.h \
			da_platform.h \
			da_probstat.h \
			da_random.h \
			da_signal.h \
			da_timeseg.h \
			da_util.h
#			da_cluster_pp.h \
#			da_comm_pp.h \
#			da_io_pp.h \
#			da_signal_pp.h

SOURCES 	= $(SOURCES.h) $(SOURCES.c)

OBJECTS 	= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC			= gcc
#CC			= purify $(DA_CC)
CFLAGS 			= -m64
#CFLAGS		= -pg
#CFLAGS          = -g
CPPFLAGS	= -I../../include

# Special make commands.

.PRECIOUS:	$(LIBRARY)

# Standard targets.

#all:		depend $(SOURCES) $(OBJECTS) $(LIBRARY) install
all:		$(SOURCES) $(OBJECTS) $(LIBRARY) install
objects:	$(SOURCES) $(OBJECTS)
sources:	$(SOURCES)
test:		

install: 
			cp $(LIBRARY) $(DESTDIR)
			cp $(SOURCES.h) $(INCDIR)

clean:	
			/bin/rm -f core $(OBJECTS) $(LIBRARY)

depend:
			makedepend -I$(C_CODE_HOME)/include *.c

# Targets' dependencies.

$(LIBRARY):	$(OBJECTS)
#			ar rv $@ $?
			libtool -c -static -o $@ $?
#			libtool -static -o $@ $?

# Using implicit .c.o compile rule

# Specific objects' dependencies.
# DO NOT DELETE THIS LINE -- make depend depends on it.
