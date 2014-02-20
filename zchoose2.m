
include $(HOME)/programs/Make.defs

SRC = zchoose2.f

OBJS = readfield.o get1ddeepspec.o showspec.o findends.o linrebinerr.o \
       medsmooth.o

zchoose2: zchoose2.f $(OBJS)
	$(FC) -O -o zchoose2 zchoose2.f $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

