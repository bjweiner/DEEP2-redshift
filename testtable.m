
include $(HOME)/programs/Make.defs

testtable: readfield.o
	$(FC) -o testtable testtable.f readfield.o $(FITSIO)
