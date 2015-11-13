.SUFFIXES :
.SUFFIXES : .F90 .f90 .f .o

TARGETS = micore
RET_SRCS = micore_core.f90 micore_main.f90
RET_OBJS = $(RET_SRCS:.f90=.o)

FC = gfortran
FCFLAGS = -O -Wall

all: $(TARGETS)

.f90.o :
	$(FC) -c $< $(FCFLAGS)

micore : $(RET_OBJS)
	$(FC) -o micore $^ $(FCFLAGS)


micore_main.o : micore_core.o
micore_core.o : 

clean:
	rm -f *.o *.mod $(TARGETS) *~

