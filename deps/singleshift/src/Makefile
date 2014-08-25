include ../environment

OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90))

all: $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(LIBS)

clean:
	-rm $(OBJS)	
