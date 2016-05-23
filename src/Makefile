.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

PROGRAM = main

OBJS = $(PROGRAM).o constants.o dmft.o \
	   ed_solver.o matsubara_grid.o mpi.o \
	   impurity_solver.o dmft_params.o

default: $(PROGRAM)
##################
### FDF Module ###
##################
FDF=fdf/libfdf.a

LIBS+=$(FDF)
FFLAGS+=-I./fdf

####################
### Dependencies ###
####################
$(OBJS): $(FDF)
$(PROGRAM).o: dmft.o mpi.o
dmft.o: dmft_params.o  matsubara_grid.o impurity_solver.o
impurity_solver.o: ed_solver.o
matsubara_grid.o: constants.o

$(FDF): 
	cd fdf ; $(MAKE);

$(PROGRAM): $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

%.o:%.mod
