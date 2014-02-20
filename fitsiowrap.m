
FC = f77
# FC = gfortran

fitsiowrap: fitsiowrap.f
  $(FC) -c -o fitsiowrap.o fitsiowrap.f

fitsiowrap.o: fitsiowrap.f
  $(FC) -c -o fitsiowrap.o fitsiowrap.f
