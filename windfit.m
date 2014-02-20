
include $(HOME)/programs/Make.defs

SRC = windfit.f windfit.m \
       readfield.f 
#       read1dspec.f readfield.f linrebinerr.f

OBJS = windfit.o \
       readfield.o
#       read1dspec.o readfield.o linrebinerr.o

#FITGAUS = $(HOME)/programs/statprogs/fitgaus2.o
FITGAUS = $(HOME)/programs/statprogs/fitgaus2special.o
FITDOUB = $(HOME)/programs/statprogs/fitdoublet.o

windfit: $(OBJS)
	 $(FC) -O -o windfit $(OBJS) $(FITSIOWRAP) $(FITGAUS) $(FITDOUB) $(MRQFIT) $(PGPLOT) $(NUMREC)
