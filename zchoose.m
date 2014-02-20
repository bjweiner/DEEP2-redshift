
include $(HOME)/programs/Make.defs

SRC = zchoose.f

zchoose: zchoose.f
	$(FC) -O -o zchoose zchoose.f $(PGPLOT) $(FITSIOWRAP)

