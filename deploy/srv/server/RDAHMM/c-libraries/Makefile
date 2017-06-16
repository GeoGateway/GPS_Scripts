#
# Makefile for Common Libraries
# 

# targets:
#   (default) -- make library and subscriber
#   all       -- same as default
#   install   -- same as default 
#              (NB: does not copy into /usr/local/bin or the like!)
#   clean     -- remove files not needed for execution
#   pristine  -- return to state as distributed

# convert an "all" target into "install"
# these will move the binary from the src/ dir to ../../bin/
# The rationale is that a make all at this level
# should mean installing the subscriber fully.
install: all

pristine: clean

# instructions for all and install, clean and pristine
all clean:
	cd nr/src; $(MAKE) $@
	cd ut/src; $(MAKE) $@
	cd cp/src; $(MAKE) $@
	cd da/src; $(MAKE) $@
	@ echo Common make done
