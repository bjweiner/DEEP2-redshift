
include $(HOME)/programs/Make.defs

SRC = pcatest.f pcatest.m \
       readspec.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f

OBJS = pcatest.o \
       readspec.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o

pcatest: $(OBJS)
	$(FC) -O -o pcatest $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm pcatest $(OBJS)

tarball:
	tar cvf pcatest.tar $(SRC)
