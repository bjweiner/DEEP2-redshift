
include $(HOME)/programs/Make.defs

SRC = deepconvert.f

OBJS = readfield.o get1ddeepspec.o linrebinerr.o \
       doerr.o

deepconvert: deepconvert.f $(OBJS)
	$(FC) -O -o deepconvert deepconvert.f $(OBJS) $(FITSIO) 