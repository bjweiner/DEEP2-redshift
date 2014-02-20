
include $(HOME)/programs/Make.defs

SRC = plotpstamp.f

OBJS = readfield.o \
       medsmooth.o doerr.o 

plotpstamp: plotpstamp.f $(OBJS)
	$(FC) -O -o plotpstamp plotpstamp.f $(OBJS) $(PGPLOT) $(FITSIO) $(NUMREC)

