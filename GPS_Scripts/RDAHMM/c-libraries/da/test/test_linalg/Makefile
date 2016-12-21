# Makefile for test-lib
#
# (Based on original merf makefile)
# $Id: Makefile,v 1.2 1996/02/19 21:34:42 agray Exp $
# $Log: Makefile,v $
# Revision 1.2  1996/02/19 21:34:42  agray
# clarified comment
# AG
#
# Revision 1.1  1996/02/06 03:31:22  agray
# Initial revision
#


BINDIR		= $(DA_HOME)/test

UTILLIBDIR	= $(ML_HOME)
NRLIBDIR    = $(NR_HOME)
DALIBDIR    = $(DA_HOME)

PROGRAM		= $(BINDIR)/test-lib

SOURCES.c	= \
		test-lib.c

SOURCES.h	= \
		test-lib.h

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

CC			= $(DA_CC)
#CC			= purify -windows=yes $(DA_CC)
CFLAGS		= $(DA_CFLAGS)
CPPFLAGS	= \
		-I$(UTILLIBDIR)/include \
		-I$(NRLIBDIR)/include \
		-I$(DALIBDIR)/include

LDFLAGS		= \
		-L$(UTILLIBDIR)/lib \
		-L$(NRLIBDIR)/lib \
		-L$(DALIBDIR)/lib

LDLIBS		= -lda -lrecipes_c -lutil -lm

all:		$(SOURCES) $(PROGRAM)
objects:	$(SOURCES.c) $(OBJECTS)
sources:	$(SOURCES)

$(PROGRAM):	$(SOURCES.c) $(SOURCES.h) $(OBJECTS) Makefile
		$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:		$(RM) $(OBJECTS) *~ core

new:		clean all


# specific objects dependencies (in case libraries and/or headers changing)

$(OBJECTS):	$(UTILLIBDIR)/lib/libutil.a \
			$(NRLIBDIR)/lib/librecipes_c.a \
			$(DALIBDIR)/lib/libda.a 
