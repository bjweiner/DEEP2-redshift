
include $(HOME)/programs/Make.defs

SRC = linefit.f linefit.m \
       read1dspec.f readfield.f calcsnrmed.f

OBJS = linefit.o \
       read1dspec.o readfield.o calcsnrmed.o

FITGAUS = $(HOME)/programs/statprogs/fitgaus2.o
FITDOUB = $(HOME)/programs/statprogs/fitdoublet.o
FITDOUB2 = $(HOME)/programs/statprogs/fitdoublet2a.o

linefit: $(OBJS)
	$(FC) -O -o linefit $(OBJS) $(FITGAUS) $(PGPLOT) $(FITSIOWRAP) \
  $(BIWGT) $(MRQFIT) $(FITDOUB) $(FITDOUB2) $(NUMREC)
