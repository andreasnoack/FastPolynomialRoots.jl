##################################################
# Makefile
##################################################
# Set environment in file environment
include environment

##################################################
TESTS := $(patsubst $(TODIR)/%.f90,$(TODIR)/%.o,$(wildcard $(TODIR)/*.f90))

all:
	make -C $(TESTDIR)

clean:
	make clean -C $(TODIR) &&\
	make clean -C $(SRCDIR) &&\
	make clean -C $(TESTDIR)	
