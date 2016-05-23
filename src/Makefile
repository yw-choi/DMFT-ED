.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
include ../arch.make

PROGRAM = dmft_ed

OBJS = $(PROGRAM).o constants.o check_dmft_converged.o dmft.o dmft_finalize.o dmft_init.o \
	   dmft_post_processing.o ed_solver.o matsubara_grid.o mpi.o new_local_green_ftn.o \
	   update_weiss_field.o impurity_solver.o

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
dmft_init.o: mpi.o impurity_solver.o dmft.o
dmft_ed.o: dmft_init.o dmft_finalize.o matsubara_grid.o

matsubara_grid.o: constants.o
dmft.o: matsubara_grid.o
ed_solver.o: dmft.o matsubara_grid.o
impurity_solver.o: ed_solver.o


$(FDF): 
	cd fdf ; $(MAKE);

$(PROGRAM): $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

#%.o:%.mod
