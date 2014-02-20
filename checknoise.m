
include $(HOME)/programs/Make.defs

#OBJS = checknoise.o scalenoise.o showspec.o
OBJS = checknoise.o scalenoise.o findends.o

checknoise: $(OBJS)
	$(F77) -o checknoise $(OBJS) $(PGPLOT) $(FITSIOWRAP)
