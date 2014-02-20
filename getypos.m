
include $(HOME)/programs/Make.defs

getypos: getypos.f readfield.o doerr.o
	$(FC) -O -o getypos getypos.f readfield.o doerr.o $(FITSIO)

