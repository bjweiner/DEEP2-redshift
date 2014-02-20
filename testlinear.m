
include $(HOME)/programs/Make.defs

testlinear: testlinear.f linrebin.o
	$(FC) -O -o testlinear testlinear.f linrebin.o

