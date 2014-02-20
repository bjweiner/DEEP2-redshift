
c  de-redshift a spectrum by z, in log space.  this means
c  shifting wl scale by -log(1+z).

c do the shift in place, or not
c      subroutine deshiftz(spec,nsp,dw,z,specout)
      subroutine deshiftz(spec,nsp,dwlog,z)

      real spec(nsp)
      real specout(nsp)

      include 'pcredshift.h'

      dzlog = log10(1.+z)
      dindx = dzlog / dwlog
      ishift = int(dindx)
      fshift = dindx - ishift
c now ishift and fshift are the positive integer and fractional
c shifts (if z is positive).  we want to shift to lower w

      do i=1,nsp
         j1 = i + ishift
         if (spec(j1) .gt. bad .and. spec(j1+1) .gt. bad) then
            specout(i) = (1.-fshift)*spec(j1) + fshift*spec(j1+1)
         else if (spec(j1) .gt. bad) then
            specout(i) = spec(j1)
         else if (spec(j1+1) .gt. bad) then
            specout(i) = spec(j1+1)
         else
            specout(i) = bad
         end if
      end do
c to deshift in place
      do i=1,nsp
         spec(i) = specout(i)
      end do

      return
      end


      
