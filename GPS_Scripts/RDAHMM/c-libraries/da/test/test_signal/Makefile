# Makefile for test-new
#
# (Based on em makefile)
# To test new functions in da_linalg
# RG

PROGRAM		= test-signal

DESTDIR		= $(C_CODE_HOME)/bin

SOURCES.c	= \
			test-lib.c

SOURCES.h	=

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC		= $(CT_CC)
CFLAGS		= $(CT_CFLAGS)
CPPFLAGS	= -I$(C_CODE_HOME)/include
LDFLAGS		= -L$(C_CODE_HOME)/lib
LDLIBS		= -lda -lnr -lut -lcp -lm

# Special make commands.

.PRECIOUS: $(PROGRAM)

# Standard targets.

all:		$(SOURCES) $(PROGRAM) install
objects:	$(SOURCES.c) $(OBJECTS)
sources:	$(SOURCES)
test:		

install:
			cp $(PROGRAM) $(DESTDIR)

clean:	
			/bin/rm -f core $(OBJECTS) $(PROGRAM)

# Targets' dependencies.

$(PROGRAM):	$(OBJECTS) Makefile
			$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

$(OBJECTS):	\
			$(C_CODE_HOME)/lib/libcp.a \
			$(C_CODE_HOME)/lib/libut.a \
			$(C_CODE_HOME)/lib/libnr.a \
			$(C_CODE_HOME)/lib/libda.a 

# Specific objects' dependencies.  Use "makedepend -I/tools/code/c/include *.c".
# DO NOT DELETE THIS LINE -- make depend depends on it.
