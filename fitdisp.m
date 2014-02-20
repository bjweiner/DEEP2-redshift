
include $(HOME)/programs/Make.defs

SRC = fitdisp.f fitdisp.m \
       readspec.f showspec.f logrebinerr.f deshiftz2.f contsubmed.f \
       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
       contsubconst.f medsmooth.f findends.f cleanends.f blankout.f \
       scalenoise.f readevects.f blankoutsky2.f blankoutlines.f svdsubs.f findindex.f \
       findlocmin.f findmultmin.f get1ddeepspec.f readfield.f linrebinerr.f \
       contsubpoly.f contsubleg.f localsigclip.f findbadends.f \
       selectabsregions.f contdivpctile.f get1dcoaddspec.f

OBJS = fitdisp.o \
       readspec.o showspec.o logrebinerr.o deshiftz2.o contsubmed.o \
       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
       contsubconst.o medsmooth.o findends.o cleanends.o blankout.o \
       scalenoise.o readevects.o blankoutsky2.o blankoutlines.o svdsubs.o findindex.o \
       findlocmin.o findmultmin.o get1ddeepspec.o readfield.o linrebinerr.o \
       contsubpoly.o contsubleg.o localsigclip.o findbadends.o \
       selectabsregions.o contdivpctile.o get1dcoaddspec.o

fitdisp: $(OBJS)
	$(FC) -O -o fitdisp $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

#clean:
#	rm fitdisp $(OBJS)

#tarball:
#	tar cvf fitdisp.tar $(SRC)
