
c do 0th order continuum subtraction by just subtracting the mean
c of the good data points.  (Could also try the median?)

c do the subtraction in place, or not
c      subroutine contsubmed(spec,nsp,dwlog,wrad,specout)
      subroutine contsubconst(spec,nsp)

      real spec(nsp)

      include 'pcredshift.h'

      np=0
      sum=0.0
      do i=1,nsp
         if(spec(i) .gt. bad) then
            np=np+1
            sum = sum + spec(i)
         end if
      end do
      if (np .ge. 1) then
         rmean = sum/np
         do i=1,nsp
            if (spec(i) .gt. bad) then
               spec(i) = spec(i) - rmean
            end if
         end do
      end if

      return
      end

