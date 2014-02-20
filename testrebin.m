
include $(HOME)/programs/Make.defs

SRC = testrebin.f

OBJS = showspecerr.o findends.o linrebinerr.o logrebinerr.o\
       linrebin.o logrebin.o
#       medsmooth.o

testrebin: testrebin.f $(OBJS)
	$(FC) -O -o testrebin testrebin.f $(OBJS) $(PGPLOT) 

