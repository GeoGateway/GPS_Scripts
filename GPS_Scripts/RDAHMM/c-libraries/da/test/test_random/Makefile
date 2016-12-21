# Makefile for kmeans
#
# (Based on pca in CoolTools)
# $Id$
# $Log$
#

BINDIR		= .

UTILLIBDIR	= $(ML_HOME)

PROGRAM		= $(BINDIR)/test_rand

SOURCES.c	= \
		test_rand.c

SOURCES.h	= \
		test_rand.h

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

CC			= gcc
#CC			= purify -windows=yes $(DA_CC)
CFLAGS		= $(CT_CFLAGS)
CPPFLAGS	= \
		-I$(UTILLIBDIR)/include

LDFLAGS		= \
		-L$(UTILLIBDIR)/lib 

LDLIBS		= -lutil -lm

all:		$(SOURCES) $(PROGRAM)
objects:	$(SOURCES.c) $(OBJECTS)
sources:	$(SOURCES)

$(PROGRAM):	$(SOURCES.c) $(SOURCES.h) $(OBJECTS) Makefile
		$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

info:
	echo sources = $(SOURCES)
	echo utildir = $(UTILLIBDIR)
	echo ldflags = $(LDFLAGS)
	echo ldlibs = $(LDLIBS)

clean:		
		$(RM) $(OBJECTS) *~ core

new:		clean all
