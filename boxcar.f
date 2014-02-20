
c  boxcar smooth a spectrum with a window of radius = irad pixels
c  the original spectrum is preserved.

      subroutine boxcar(n,spin,irad,sout)

      real spin(n),sout(n)
      integer irad

      include 'pcredshift.h'

      i1 = 1+irad
      i2 = n-irad
      do i=1,i1-1
         sout(i) = spin(i)
      end do
      do i=i2+1,n
         sout(i) = spin(i)
      end do
      do i=i1,i2
         sum=0.0
         np =0
         do j=i-irad,i+irad
            if (spin(j) .gt. bad) then
               np=np+1
               sum=sum+spin(j)
            end if
         end do
         if (np .ge. 1) then
            sout(i) = sum/np
         else
            sout(i) = badset
         end if
      end do

      return
      end
