c  normalize an array of spectra, in place.
c  normalize to the modulus (sum of squares) -
c  this isn't very robust.

      subroutine normalize(specarray,nspec,nspmax,nwave)

      real specarray(nspmax,nwave)

      include 'pcredshift.h'

      do i=1,nspec
         np = 0
         sumsq = 0.0
         do j=1,nwave
            if (specarray(i,j) .gt. bad) then
               np = np+1
               sumsq = sumsq + specarray(i,j)*specarray(i,j)
            end if
         end do
         smod = sqrt(sumsq)
         if (np .ge. 1) then
            do j=1,nwave
               if(specarray(i,j) .gt. bad) then
                  specarray(i,j) = specarray(i,j) / smod
               end if
            end do
         end if
      end do

      return
      end
