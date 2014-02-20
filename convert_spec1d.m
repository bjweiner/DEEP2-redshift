
include $(HOME)/programs/Make.defs

SRC = convert_spec1d.f convert_spec1d.m \
       read1dspec.f readfield.f linrebinerr.f

OBJS = convert_spec1d.o \
       read1dspec.o readfield.o linrebinerr.o

convert_spec1d: $(OBJS)
	$(FC) -O -o convert_spec1d $(OBJS) $(FITSIOWRAP) 

