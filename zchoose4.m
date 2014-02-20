
include $(HOME)/programs/Make.defs

SRC = zchoose4.f

OBJS = readfield.o get1ddeepspec.o showspecerr.o findends.o linrebinerr.o \
       medsmooth.o doerr.o weightsmooth.o boxcar.o showspec.o

zchoose4: zchoose4.f $(OBJS)
	$(FC) -O -o zchoose4 zchoose4.f $(OBJS) $(PGPLOT) $(FITSIO) $(NUMREC)

