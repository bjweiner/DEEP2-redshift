
include $(HOME)/programs/Make.defs

SRC = plotdeepspec.f plotdeepspec.m \
       read1dspec.f readfield.f medsmooth.f

OBJS = plotdeepspec.o \
       read1dspec.o readfield.o medsmooth.o

plotdeepspec: $(OBJS)
	$(FC) -O -o plotdeepspec $(OBJS) $(FITSIOWRAP) $(PGPLOT) $(NUMREC)


