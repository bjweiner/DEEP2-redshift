
include $(HOME)/programs/Make.defs

SRC = pcatestsky.f pcatestsky.m \
       readspec.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f \
       scalenoise.f

OBJS = pcatestsky.o \
       readspec.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o \
       scalenoise.o

pcatestsky: $(OBJS)
	$(FC) -O -o pcatestsky $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm pcatestsky $(OBJS)

tarball:
	tar cvf pcatestsky.tar $(SRC)
