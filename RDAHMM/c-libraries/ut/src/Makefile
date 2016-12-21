# Makefile for UT library

# $Id: Makefile,v 1.5 1999/10/06 02:40:12 cataract Exp granat $
# $Log: Makefile,v $
# Revision 1.5  1999/10/06 02:40:12  cataract
# Added position-independent flags
#
# Revision 1.4  1998/04/20 22:40:27  roden
# - Simplified $(LIBRARY) rule so that Gnu make doesn't call ar unnecessarily.
# - Removed %.o rule all together since implicit rule was being used mostly,
#   and explicit rule caused two compilations under certain circumstances,
#   not to mention $(OUTPUT_OPTION) was screwing things up under Linux.
#
# Revision 1.3  1997/07/31 22:08:21  agray
# added makedepend facility.
#
# Revision 1.2  1997/07/31 22:04:33  granat
# added flags for profiler/purify
#
# Revision 1.1  1997/05/13 23:33:13  granat
# Initial revision
#
#

# Library, source code, and objects.

LIBRARY		= libut.a

DESTDIR		= ../../lib

INCDIR		= ../../include

SOURCES.c	= \
			ut_args.c \
			ut_error.c \
			ut_file_io.c \
			ut_math.c \
			ut_memory.c \
			ut_output.c \
			ut_platform.c \
			ut_string.c \
			ut_time.c \
			ut_types.c

SOURCES.h	= \
			ut_args.h \
			ut_error.h \
			ut_file_io.h \
			ut_math.h \
			ut_memory.h \
			ut_output.h \
			ut_platform.h \
			ut_string.h \
			ut_time.h \
			ut_types.h \
			ut.h

SOURCES 	= $(SOURCES.h) $(SOURCES.c)

OBJECTS 	= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC		= gcc
#CC		= purify $(UT_CC)
#CFLAGS		= -m64
CFLAGS		= -m64
#CFLAGS		= -pg $(UT_CFLAGS)
#CFLAGS		= $(UT_FLAGS) -fpic
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
			cp libut.dylib $(DESTDIR)/../dylib
			cp $(SOURCES.h) $(INCDIR)

clean:
			/bin/rm -f core $(OBJECTS) $(LIBRARY)

depend:
			makedepend -I$(C_CODE_HOME)/include *.c

# Targets' dependencies.

$(LIBRARY):	$(OBJECTS)
#			ar rv $@ $?
			libtool -c -static -o $@ $?
			libtool -v -arch_only x86_64 -dynamic -o libut.dylib $? -lc

# Using implicit .c.o compile rule

# Specific objects' dependencies.
# DO NOT DELETE THIS LINE -- make depend depends on it.
