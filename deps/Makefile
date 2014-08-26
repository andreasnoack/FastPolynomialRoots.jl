FC = gfortran
ifeq ($(OS), Windows_NT)
	SLIB = dll
else
	UNAME := $(shell uname)
	ifeq ($(UNAME), Darwin)		
		SLIB = dylib
	else
		SLIB = so
	endif
endif
LBLAS = -lblas
FFLAGS = -shared -fPIC -O3

OBJD = 	doubleshift/src/DAMVW.f90 \
		doubleshift/src/DCB.f90 \
		doubleshift/src/DCDB.f90 \
		doubleshift/src/DCFD.f90 \
		doubleshift/src/DCFT.f90 \
		doubleshift/src/DFCC.f90 \
		doubleshift/src/DFGR.f90 \
		doubleshift/src/DGR.f90 \
		doubleshift/src/DMQF.f90 \
		doubleshift/src/DNORMALPOLY.f90 \
		doubleshift/src/DRANDPOLYJT.f90 \
		doubleshift/src/RESCHECK.f90 \
		doubleshift/src/balance.f90 \
		doubleshift/src/init_random_seed.f90 \
		doubleshift/src/turnovers/DGTO2.f90 
OBJS =	singleshift/src/balance.f90 \
		singleshift/src/buildbulge.f90 \
		singleshift/src/chasebulge.f90 \
		singleshift/src/cnormalpoly.f90 \
		singleshift/src/crgivens.f90 \
		singleshift/src/deflation.f90 \
		singleshift/src/diagblock.f90 \
		singleshift/src/factor.f90 \
		singleshift/src/fuse.f90 \
		singleshift/src/init_random_seed.f90 \
		singleshift/src/modified_quadratic.f90 \
		singleshift/src/normalpoly.f90 \
		singleshift/src/rescheck.f90 \
		singleshift/src/throughdiag.f90 \
		singleshift/src/zamvw.f90 \
		singleshift/src/zamvw2.f90 \
		singleshift/src/turnovers/dto4.f90

all:
	$(FC) $(FFLAGS) $(OBJD) -o libamvwdouble.$(SLIB) $(LBLAS)
	$(FC) $(FFLAGS) $(OBJS) -o libamvwsingle.$(SLIB) $(LBLAS)

clean:
	rm *.$(SLIB)