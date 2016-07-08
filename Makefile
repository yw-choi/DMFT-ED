.PHONY: all clean

default: all

all:
	(cd src; $(MAKE); cd ..; \
		if [ ! -f bin ]; then mkdir bin; fi;\
		cp src/main.x bin/)	
clean:
	(cd src; $(MAKE) clean)
