
include $(HOME)/programs/Make.defs

SRC = plotpstamp_allslit.f

plotpstamp_allslit: plotpstamp_allslit.f
	$(FC) -O -o plotpstamp_allslit plotpstamp_allslit.f $(PGPLOT) $(FITSIOWRAP)

