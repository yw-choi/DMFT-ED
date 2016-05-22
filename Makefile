.PHONY: all clean

default: all

all:
	(cd src; $(MAKE); cd ..; cp src/dmft_ed.x bin)	
clean:
	(cd src; $(MAKE) clean)
