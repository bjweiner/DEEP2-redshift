
include $(HOME)/programs/Make.defs

SRC = fitz3.f fitz3.m \
       readspec.f showspec.f logrebinerr.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f \
       scalenoise.f readevects.f blankoutsky.f svdsubs.f findindex.f \
       findlocmin.f findmultmin.f get1ddeepspec.f readfield.f linrebinerr.f \
       contsubpoly.f contsubleg.f localsigclip.f findbadends.f \
       centroidmin.f

OBJS = fitz3.o \
       readspec.o showspec.o logrebinerr.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o \
       scalenoise.o readevects.o blankoutsky.o svdsubs.o findindex.o \
       findlocmin.o findmultmin.o get1ddeepspec.o readfield.o linrebinerr.o \
       contsubpoly.o contsubleg.o localsigclip.o findbadends.o \
       centroidmin.o

fitz3: $(OBJS)
	$(FC) -O -o fitz3 $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm fitz3 $(OBJS)

tarball:
	tar cvf fitz3.tar $(SRC)
