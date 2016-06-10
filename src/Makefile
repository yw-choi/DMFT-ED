.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

PROGRAM = main

OBJS = $(PROGRAM).o constants.o dmft.o \
	   ed_solver.o matsubara_grid.o mpi.o \
	   impurity_solver.o dmft_params.o io_units.o numeric_utils.o \
	   dmft_lattice.o ed_projection.o ed_params.o ed_basis.o ed_full_diag.o \
	   ed_hamiltonian.o ed_operator.o

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
$(PROGRAM).o: mpi.o io_units.o dmft.o
dmft.o: dmft_params.o matsubara_grid.o impurity_solver.o dmft_lattice.o
impurity_solver.o: ed_solver.o 
matsubara_grid.o: constants.o
dmft_lattice.o: numeric_utils.o
ed_solver.o: ed_params.o ed_projection.o ed_full_diag.o ed_hamiltonian.o
ed_full_diag.o: ed_hamiltonian.o
ed_hamiltonian.o: ed_operator.o ed_basis.o 
ed_params.o: numeric_utils.o

$(FDF): 
	cd fdf ; $(MAKE);

$(PROGRAM): $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

%.o:%.mod
