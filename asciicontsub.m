

include $(HOME)/programs/Make.defs

OBJS = contsubmed.o contsubpoly.o contsubleg.o svdsubs.o localsigclip.o

asciicontsub: $(OBJS)
	$(FC) -O -o asciicontsub asciicontsub.f $(OBJS) $(NUMREC) $(PGPLOT)

