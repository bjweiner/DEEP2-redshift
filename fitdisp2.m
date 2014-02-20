
include $(HOME)/programs/Make.defs

SRC = fitdisp2.f fitdisp2.m \
       readspec.f showspec.f logrebinerr.f deshiftz2.f contsubmed.f \
       contsubconst.f medsmooth.f \
       readevects.f blankoutsky2.f blankoutlines.f svdsubs.f \
       get1ddeepspec.f readfield.f contsubpoly.f contsubleg.f \
       selectabsregions.f contdivpctile.f findends.f cleanends.f \
       linrebinerr.f localsigclip.f get1dcoaddspec.f 
# findbadends.f \
# blankout.f \
#       readz.f boxcar.f compmean.f submean.f normalize.f zerobad.f \
#       scalenoise.f \
# findindex.f \
#       findlocmin.f findmultmin.f

OBJS = fitdisp2.o \
       readspec.o showspec.o logrebinerr.o deshiftz2.o contsubmed.o \
       contsubconst.o medsmooth.o \
       readevects.o blankoutsky2.o blankoutlines.o svdsubs.o \
       get1ddeepspec.o readfield.o contsubpoly.o contsubleg.o \
       selectabsregions.o contdivpctile.o findends.o cleanends.o \
       linrebinerr.o localsigclip.o get1dcoaddspec.o
# findbadends.o \
# blankout.o \
#       readz.o boxcar.o compmean.o submean.o normalize.o zerobad.o \
#       scalenoise.o
# findindex.o \
#       findlocmin.o findmultmin.o

fitdisp2: $(OBJS)
	$(FC) -O -o fitdisp2 $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

#clean:
#	rm fitdisp2 $(OBJS)

#tarball:
#	tar cvf fitdisp2.tar $(SRC)
