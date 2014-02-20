
include $(HOME)/programs/Make.defs

SRC = fitz.f fitz.m \
       readspec.f showspec.f logrebin.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f \
       scalenoise.f readevects.f blankoutsky.f svdsubs.f findindex.f \
       findlocmin.f findmultmin.f

OBJS = fitz.o \
       readspec.o showspec.o logrebin.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o \
       scalenoise.o readevects.o blankoutsky.o svdsubs.o findindex.o \
       findlocmin.o findmultmin.o

fitz: $(OBJS)
	$(FC) -O -o fitz $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm fitz $(OBJS)

tarball:
	tar cvf fitz.tar $(SRC)
