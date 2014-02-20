

include $(HOME)/programs/Make.defs

OBJS = get1ddeepspec.o readfield.o linrebinerr.o doerr.o showspecerr.o \
       findends.o

writedeepspec: $(OBJS)
	$(FC) -O -o writedeepspec writedeepspec.f $(OBJS) $(FITSIO) $(PGPLOT)

