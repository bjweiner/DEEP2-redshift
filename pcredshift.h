c  common values to include in all appropriate files

c  ignore values less than the bad flag
c      bad = -10000.
      bad = -600.
c      bad = -200.

c  when we actually set a pixel to be bad, make sure
c  it is less than the flag.  Hopefully this also means
c  if a bad pixel and good get combined somehow, the 
c  new value will still be below "bad"
      badset = -1.0e6


