
include $(HOME)/programs/Make.defs

SRC = testplotsky.f testplotsky.m \
       readspec.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f contsubconst.f medsmooth.f \
       findends.f cleanends.f scalenoise.f

OBJS = testplotsky.o \
       readspec.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o contsubconst.o medsmooth.o \
       findends.o cleanends.o scalenoise.o

testplotsky: $(OBJS)
	$(FC) -O -o testplotsky $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm testplotsky *.o

tarball:
	tar cvf testplotsky.tar $(SRC)
