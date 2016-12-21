# Makefile for hmm
#
# $Id$
# $Log$
#

C_CODE_HOME = /home/yuma/RDAHMM/c-libraries

PROGRAM		= rdahmm

DESTDIR		= ../bin

SOURCES.c	= \
			rdahmm_arguments.c \
			rdahmm_setup.c \
			rdahmm_control.c \
			rdahmm_gauss_output.c \
			rdahmm_generic.c \
			rdahmm_stat_utils.c \
			rdahmm.c

SOURCES.h	= \
			rdahmm_arguments.h \
			rdahmm_setup.h \
			rdahmm_control.h \
			rdahmm_gauss_output.h \
			rdahmm_generic.h \
			rdahmm_stat_utils.h \
			rdahmm.h

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC		= gcc
#CC		= /sw/bin/gcc-4
#CFLAGS		= -g -O3
CFLAGS		= -m64 -g
#CFLAGS		= -g -Wall
CPPFLAGS	= -I$(C_CODE_HOME)/include
LDFLAGS		= -L$(C_CODE_HOME)/lib
#LDLIBS		= -lda -lnr -lut -lcp -lm
#Note: -lcc_dynamic is to resolve gcc version linkage issue for g77 under OSX
#remove for other platforms
#LDLIBS		= -lda -lnr -lut -lcp -llapack -lf77blas -latlas -lm -lcc_dynamic
#LDLIBS		= -lda -lnr -lut -lcp -llapack_LINUX -lblas_LINUX -lm
LDLIBS		= -lda -lnr -lut -lcp -llapack -lblas -lm
#LDLIBS		= -lda -lnr -lut -lcp -lm -Wl,-framework -Wl,Accelerate 
#LDLIBS		= -lda_fast -lnr_fast -lut -lcp -lcblas -latlas -lm
#LDLIBS		= -lda_fast -lnr_fast -lut -lcp -lm

# Special make commands.

.PRECIOUS: $(PROGRAM)

# Standard targets.

all:		$(SOURCES) $(PROGRAM) install
objects:	$(SOURCES.c) $(OBJECTS)
sources:	$(SOURCES)
test:		

install:
			cp $(PROGRAM) $(DESTDIR)
mex:
			/Applications/MATLAB_R2010b.app/bin/mex -v -f ./mexopts.sh -DMACVERSION fit_gauss_hmm.c $(SOURCES.c) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -lnr_ut -lc 
			/Applications/MATLAB_R2010b.app/bin/mex -v -f ./mexopts.sh -DMACVERSION label_using_gauss_hmm.c $(SOURCES.c) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -lnr_ut -lc 

clean:	
			/bin/rm -f core $(OBJECTS) $(PROGRAM)

# Targets' dependencies.

$(PROGRAM):	$(OBJECTS) Makefile
#			$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)
			gcc -m64 -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

$(OBJECTS):	\
			$(C_CODE_HOME)/lib/libcp.a \
			$(C_CODE_HOME)/lib/libut.a \
			$(C_CODE_HOME)/lib/libnr.a \
			$(C_CODE_HOME)/lib/libda.a 

# Specific objects' dependencies.  Use "makedepend -I/tools/code/c/include *.c".
# DO NOT DELETE THIS LINE -- make depend depends on it.
