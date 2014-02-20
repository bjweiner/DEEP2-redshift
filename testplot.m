
include $(HOME)/programs/Make.defs

SRC = testplot.f testplot.m \
       readspec.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f contsubconst.f medsmooth.f \
       findends.f cleanends.f

OBJS = testplot.o \
       readspec.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o contsubconst.o medsmooth.o \
       findends.o cleanends.o

testplot: $(OBJS)
	$(FC) -O -o testplot $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm testplot *.o

tarball:
	tar cvf testplot.tar $(SRC)
