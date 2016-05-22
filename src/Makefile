.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf

include ../arch.make

PROGRAM = dmft_ed

OBJS = $(PROGRAM).o


default: $(PROGRAM)

##################
### FDF Module ###
##################
FDF=./fdf/libfdf.a
$(FDF):
	(cd fdf ; $(MAKE))
LIBS+=$(FDF)
FFLAGS+=-I./fdf

$(PROGRAM): $(FDF) $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)

#%.o:%.mod
