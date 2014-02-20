
c  de-redshift a spectrum by z, in log space.  this means
c  shifting wl scale by -log(1+z).

c do the shift into a new array with possibly different wl zeropoint
c but same dw/dpix
      subroutine deshiftz2(spec,nsp,w0,dwlog,z,specout,nspout,w0out)
c      subroutine deshiftz2(spec,nsp,dwlog,z)

      real spec(nsp)
      real specout(nspout)

      include 'pcredshift.h'

      dzlog = log10(1.+z)
c this is the shift due to deredshifting. positive moves high index to low
      dindx = dzlog / dwlog
c this is the shift due to zeropoint.  if w0out<w0, index for a given
c wl increases, so dindx2 is negative
      dindx2 = (w0out-w0) / dwlog
c total shift
      dindx = dindx + dindx2
      ishift = int(dindx)
      fshift = dindx - ishift

c now ishift and fshift are the positive integer and fractional
c shifts (if z is positive).  we want to shift to lower w
c  i indexes the new array
c  j1 indexes the old array
c      do i=1,nsp
      do i=1,nspout
         j1 = i + ishift
         if (j1 .ge. 1 .and. j1+1 .le. nsp) then
            if (spec(j1) .gt. bad .and. spec(j1+1) .gt. bad) then
               specout(i) = (1.-fshift)*spec(j1) + fshift*spec(j1+1)
            else if (spec(j1) .gt. bad) then
               specout(i) = spec(j1)
            else if (spec(j1+1) .gt. bad) then
               specout(i) = spec(j1+1)
            else
               specout(i) = badset
            end if
         else
            specout(i) = badset
         end if
      end do
c to deshift in place
c      do i=1,nsp
c         spec(i) = specout(i)
c      end do

      return
      end


      
