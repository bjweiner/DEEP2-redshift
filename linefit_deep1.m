
include $(HOME)/programs/Make.defs

SRC = linefit_deep1.f linefit.m \
       read1dspec_deep1.f lininterp.f findends.f

OBJS = linefit_deep1.o \
       read1dspec_deep1.o lininterp.o findends.o

FITGAUS = $(HOME)/programs/statprogs/fitgaus2.o
FITDOUB = $(HOME)/programs/statprogs/fitdoublet.o

linefit_deep1: $(OBJS)
	$(FC) -O -o linefit_deep1 $(OBJS) $(FITGAUS) $(PGPLOT) $(FITSIOWRAP) \
 $(BIWGT) $(MRQFIT) $(FITDOUB)
