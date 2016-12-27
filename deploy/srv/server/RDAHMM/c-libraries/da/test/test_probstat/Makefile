# Makefile for tmp
#
# (Based on hmm makefile)
# $Id$
# $Log$
#

BINDIR		= $(MY_BIN)

CPLIBDIR	= $(ML_HOME)
UTILLIBDIR	= $(ML_HOME)
NRLIBDIR    = $(NR_HOME)
DALIBDIR    = $(DA_HOME)

PROGRAM		= $(BINDIR)/tmp

SOURCES.c	= \
		tmp.c

SOURCES.h	=

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

CC			= $(CT_CC)
#CC			= purify -windows=yes $(CT_CC)
CC			= $(CT_CC) -pg
CFLAGS		= $(CT_CFLAGS)
CPPFLAGS	= \
		-I$(CPLIBDIR)/include \
		-I$(UTILLIBDIR)/include \
		-I$(NRLIBDIR)/include \
		-I$(DALIBDIR)/include

LDFLAGS		= \
		-L$(CPLIBDIR)/lib \
		-L$(UTILLIBDIR)/lib \
		-L$(NRLIBDIR)/lib \
		-L$(DALIBDIR)/lib

LDLIBS		= -lda -lrecipes_c -lutil -lcp -lm

all:		$(SOURCES) $(PROGRAM)
objects:	$(SOURCES.c) $(OBJECTS)
sources:	$(SOURCES)

$(PROGRAM):	$(SOURCES.c) $(SOURCES.h) $(OBJECTS) Makefile
		$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:		
			$(RM) $(OBJECTS) *~ core

new:		clean all


# specific objects dependencies (in case libraries and/or headers changing)

$(OBJECTS):	$(CPLIBDIR)/lib/libcp.a \
			$(UTILLIBDIR)/lib/libutil.a \
			$(NRLIBDIR)/lib/librecipes_c.a \
			$(DALIBDIR)/lib/libda.a 
