c  subtract the mean spectrum (or any spectrum)

      subroutine submean(specarray,nspec,nspmax,nwave,specmean)

      real specarray(nspmax,nwave)
      real specmean(nwave)

      include 'pcredshift.h'

      do i=1,nspec
         do j=1,nwave
            if (specarray(i,j) .gt. bad) then
               specarray(i,j) = specarray(i,j) - specmean(j)
            end if
         end do
      end do

      return
      end
