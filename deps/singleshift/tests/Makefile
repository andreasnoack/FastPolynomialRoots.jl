include ../environment

TOOBJS := $(patsubst $(TODIR)/%.f90,$(TODIR)/%.o,$(wildcard $(TODIR)/*.f90))
TOSRCS := $(wildcard $(TODIR)/*.f90)

OBJS := $(patsubst $(SRCDIR)/%.f90,$(SRCDIR)/%.o,$(wildcard $(SRCDIR)/*.f90))
SRCS := $(wildcard $(SRCDIR)/*.f90)

MPSRCS := $(wildcard $(MPDIR)/*.f)
MPOBJS := $(MPDIR)/mpfun90.o $(MPDIR)/mpmod90.o $(MPDIR)/mpmodm90.o $(MPDIR)/mpmodx90.o

all: zamvw2test race raceunit racejt racesp

zamvw2test: zamvw2test.f90 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o zamvw2test zamvw2test.f90 $(TOOBJS) $(OBJS) $(LIBS)

race: rootrace.f95 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o race rootrace.f95 $(TOOBJS) $(OBJS) $(LIBS)

raceunit: rootrace_unit.f95 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o raceunit rootrace_unit.f95 $(TOOBJS) $(OBJS) $(LIBS)

racejt: rootrace_jt.f95 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o racejt rootrace_jt.f95 $(TOOBJS) $(OBJS) $(LIBS)

racesp: rootrace_special_polynomials.f95 rootstocoeffs.f90 backward.f90 $(TOSRCS) $(SRCS) $(MPSRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o racesp rootrace_special_polynomials.f95 rootstocoeffs.f90 backward.f90  $(TOOBJS) $(OBJS) $(LIBS) $(MPOBJS) -I$(MPDIR)


$(SRCS):

$(TOSRCS):	

clean:
	-rm zamvw2test race raceunit racejt racesp
