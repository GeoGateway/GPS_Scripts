# Makefile for command parse library

# $Id: Makefile,v 1.7 1998/04/20 23:43:39 roden Exp granat $
# $Log: Makefile,v $
# Revision 1.7  1998/04/20 23:43:39  roden
# Changed Makefile to conform to other /proj/code/c library Makefiles,
# including use of new CP_HOME, CP_CC and related variables from a new
# cp/setup.csh, and new cp directory structure like the other libraries.
#
# Revision 1.6  1998/04/20 22:51:52  roden
# - Simplified $(LIBRARY) rule so that Gnu make doesn't call ar unnecessarily.
# - Removed %.o rule all together since implicit rule was being used mostly,
#   and explicit rule caused two compilations under certain circumstances,
#   not to mention $(OUTPUT_OPTION) was screwing things up under Linux.
#
# Revision 1.5  1996/05/01 01:41:00  agray
# Removed ranlib line, which is unnecessary and possibly confusing.
#
# Revision 1.4  1996/05/01 01:35:44  agray
# Changed install target to match one of Joe's current Solaris makefile examples.
# Part of the general updating of libcp which is happening as part of this port
# from SunOS to Solaris.  (We also changed some .c filenames to match our
# general conventions more closely.)
#
# Revision 1.3  1995/09/05 18:27:06  jctran
# Commented out compiler flag to purify, using ML_CC instead.
# JCT
#
# Revision 1.2  1995/08/03  23:38:08  jctran
# Added files for compilation.  New files are:
# cp_string.c, cp_string.h, cp_platform.h, cp_version.h, cp_version.c
# JCT
#
# Revision 1.1  1995/07/27  22:07:29  agray
# Initial revision
#

# Library, source code, and objects.

LIBRARY		= libcp.a

DESTDIR		= ../../lib

INCDIR		= ../../include

SOURCES.c	= \
			cp_string.c \
			cp_args.c \
			cp_version.c

SOURCES.h	= \
			cp_platform.h \
			cp_string.h \
			cp_args.h \
			cp_version.h

SOURCES	= $(SOURCES.h) $(SOURCES.c)

OBJECTS	= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC			= gcc
#CC			= purify $(CP_CC)
CFLAGS			= -m64
#CFLAGS		= -pg $(CP_CFLAGS)
CPPFLAGS		= -I../../include

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

# Using implicit .c.o compile rule

# Specific objects' dependencies.
# DO NOT DELETE THIS LINE -- make depend depends on it.
