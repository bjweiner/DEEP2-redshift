
c  find the first/last points with data that is not zero or bad

c      subroutine findends(n,wave,spec,imin,imax,wmin,wmax)
      subroutine findends(n,spec,imin,imax)

c      real wave(n),spec(n)
      real spec(n)

      include 'pcredshift.h'

      small = 1.e-4

c buffer sizes - ignore the first ibuff1 and last ibuff2 pixels
      ibuff1= 0
      ibuff2= 0

c find first point with good data, i.e. data is not <bad or near zero
c      wmin = wave(1)
      imin = 1 + ibuff1
      i = 1
 100  continue
      if (spec(i) .gt. bad .and. 
     $     (spec(i) .gt. small .or. spec(i) .lt. -small)) then
c         wmin = wave(i)
         imin = i
      else if (i .eq. n) then
         go to 110
      else
         i = i+1
         go to 100
      end if
 110  continue
c find last point with good data by counting down
c      wmax = wave(n)
      imax = n - ibuff2
      i=n
 150  continue
      if (spec(i) .gt. bad .and.
     $     (spec(i) .gt. small .or. spec(i) .lt. -small)) then
c         wmax = wave(i)
         imax = i
      else if (i .eq. 1) then
         go to 160
      else
         i = i-1
         go to 150
      end if
 160  continue

      return
      end
