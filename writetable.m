
include $(HOME)/programs/Make.defs

writetable: readfield.o
	$(FC) -O -o writetable writetable.f readfield.o $(FITSIO)
