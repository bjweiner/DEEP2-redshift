
include $(HOME)/programs/Make.defs

OBJS = readspec.o showspec.o logrebin.o deshiftz.o contsubmed.o

test1: test1.f $(OBJS)
	$(FC) -O -o test1 test1.f $(OBJS) $(PGPLOT) $(FITSIOWRAP) $(NUMREC)

clean:
	rm test1 *.o
