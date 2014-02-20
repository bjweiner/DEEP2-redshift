
include $(HOME)/programs/Make.defs


OBJS = getskysigma.o readskysigma.o readobjno.o

getskysigma: $(OBJS)
	$(FC) -O -o getskysigma $(OBJS) $(FITSIOWRAP) 
