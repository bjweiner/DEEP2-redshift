c  Zero out all the "bad" flagged values in an array of spectra.
c  Something to do prior to SVD'ing.

      subroutine zerobad(specarray,nspec,nspmax,nwave,nwmax)

      real specarray(nspmax,nwmax)

      include 'pcredshift.h'

c      do i=1,nspec
      do i=1,nspmax
c         do j=1,nwave
         do j=1,nwmax
            if (specarray(i,j) .gt. bad) then
               continue
            else
               specarray(i,j) = 0.0
            end if
         end do
      end do

      return
      end
