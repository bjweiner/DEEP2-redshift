
include $(HOME)/programs/Make.defs

SRC = pcadeeptwo.f pcadeeptwo.m \
       readspecdeeptwo.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f \
       scalenoise.f get1ddeepspec.f \
       doerr.f linrebin.f readfield.f

OBJS = pcadeeptwo.o \
       readspecdeeptwo.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o \
       scalenoise.o get1ddeepspec.o \
       doerr.o linrebin.o readfield.o

pcadeeptwo: $(OBJS)
	$(FC) -O -o pcadeeptwo $(OBJS) $(PGPLOT) $(FITSIO) $(NUMREC)

clean:
	rm pcadeeptwo $(OBJS)

tarball:
	tar cvf pcadeeptwo.tar $(SRC)
