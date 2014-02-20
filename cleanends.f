
c  flag probably-bad values near the ends of the spectrum

      subroutine cleanends(n,spec,imin,imax)

      real spec(n)

      include 'pcredshift.h'

c  use a very tight low threshold and a loose high threshold
c  I don't think the high threshold will ever actually catch anything
c  but it shouldn't hurt.
c      cutlow = bad / 10.
c      cutlow = -30.
c      cuthigh = 1000.
      cutlow = -300.
      cuthigh = 10000.

c  distance from end to clean
      ibuff=16

      do i=imin,imin+ibuff
         if (spec(i) .lt. cutlow .or. spec(i) .gt. cuthigh)
     $        spec(i) = badset
      end do
      do i=imax-ibuff,imax
         if (spec(i) .lt. cutlow .or. spec(i) .gt. cuthigh)
     $        spec(i) = badset
      end do

      return
      end
