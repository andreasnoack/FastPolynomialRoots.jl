include ../environment

TOOBJS := $(patsubst $(TODIR)/%.f90,$(TODIR)/%.o,$(wildcard $(TODIR)/*.f90))
TOSRCS := $(wildcard $(TODIR)/*.f90)

OBJS := $(patsubst $(SRCDIR)/%.f90,$(SRCDIR)/%.o,$(wildcard $(SRCDIR)/*.f90))
SRCS := $(wildcard $(SRCDIR)/*.f90)

MPSRCS := $(wildcard $(MPDIR)/*.f)
MPOBJS := $(MPDIR)/mpfun90.o $(MPDIR)/mpmod90.o $(MPDIR)/mpmodm90.o $(MPDIR)/mpmodx90.o

all: test race raceunit racejt racesp racebs

test: testDAMVW.f95 $(TOSRCS) $(SRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o test testDAMVW.f95 $(TOOBJS) $(OBJS) $(LIBS)

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

racesp: rootrace_special_polynomials.f95 rootstocoeffs.f90 backward.f95 $(TOSRCS) $(SRCS) $(MPSRCS)
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o racesp rootrace_special_polynomials.f95 rootstocoeffs.f90 backward.f95  $(TOOBJS) $(OBJS) $(LIBS) $(MPOBJS) -I$(MPDIR)

racebs: rootrace_backward_stability.f95 $(TOSRCS) $(SRCS) 
	make -C $(TODIR)
	make -C $(SRCDIR)
	$(FC) $(FFLAGS) -o racebs rootrace_backward_stability.f95 $(TOOBJS) $(OBJS) $(LIBS)


$(SRCS):

$(TOSRCS):

$(CSRCS):	

clean:
	-rm test race raceunit racejt racesp racebs damvw_errors.txt degrees.txt fig_damvw_times.txt fig_dgeev_times.txt fig_lapack_times.txt lapack_times.txt sp_table1.txt damvw_times.txt fig_damvw_errors.txt fig_dgeev_errors.txt fig_lapack_errors.txt lapack_errors.txt sp_degrees.txt sp_table2.txt

