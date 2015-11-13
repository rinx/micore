.SUFFIXES :
.SUFFIXES : .F90 .f90 .f .o

TARGETS = retrieval
RET_SRCS = lib_retrieval.f90 retrieval.f90
RET_OBJS = $(RET_SRCS:.f90=.o)

FC = gfortran
FCFLAGS = -O -Wall

all: $(TARGETS)

.f90.o :
	$(FC) -c $< $(FCFLAGS)

retrieval : $(RET_OBJS)
	$(FC) -o retrieval $^ $(FCFLAGS)


retrieval.o : lib_retrieval.o
lib_retrieval.o : 


clean:
	rm -f *.o *.mod $(TARGETS) *~

