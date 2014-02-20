
c  bad-flag all of the spectrum outside some region of interest

      subroutine blankout(specarray,nspec,nspmax,nwave,wave,wmin,wmax)

      real specarray(nspmax,nwave)
      real wave(nwave)
      
      include 'pcredshift.h'

c  find the indexes corresponding to the region of interest
c  wave array is assumed monotonically increasing
      i=1
 100  continue
      if (i .gt. nwave .or. wave(i) .gt. wmin) then
         imin=i
         go to 110
      end if
      i=i+1
      go to 100
 110  continue

      i=nwave
 150  continue
      if (i .lt. 1 .or. wave(i) .lt. wmax) then
         imax=i
         go to 160
      end if
      i=i-1
      go to 150
 160  continue

      do j=1,nspec
         do i=1,imin-1
            specarray(j,i) = badset
         end do
         do i=imax+1,nwave
            specarray(j,i) = badset
         end do
      end do

      return
      end

