
include $(HOME)/programs/Make.defs

SRC = zchoose3.f

OBJS = readfield.o get1ddeepspec.o showspecerr.o findends.o linrebinerr.o \
       medsmooth.o doerr.o showspec.o

zchoose3: zchoose3.f $(OBJS)
	$(FC) -O -o zchoose3 zchoose3.f $(OBJS) $(PGPLOT) $(FITSIO) $(NUMREC)

