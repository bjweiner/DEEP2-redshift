
include $(HOME)/programs/Make.defs

compsigma: compsigma.f
	$(FC) -O -o compsigma compsigma.f
