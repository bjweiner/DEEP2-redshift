c  compute the mean spectrum
c  nspec is the # of spectra, nspmax is physical dimension of array

      subroutine compmean(specarray,nspec,nspmax,nwave,specmean,specrms)

      real specarray(nspmax,nwave)
      real specmean(nwave),specrms(nwave)

      include 'pcredshift.h'

      do j=1,nwave
         sum = 0.0
         sumsq = 0.0
         np = 0
         do i=1,nspec
            if (specarray(i,j) .gt. bad) then
               sum = sum + specarray(i,j)
               sumsq = sumsq + specarray(i,j)*specarray(i,j)
               np = np + 1
            end if
         end do
         if (np .eq. 0) then
            specmean(j) = badset
            specrms(j) = badset
         else if (np .eq. 1) then
            specmean(j) = sum/np
            specrms(j) = badset
         else
            specmean(j) = sum/np
            specrms(j) = sqrt(sumsq/np - specmean(j)*specmean(j))
         end if
      end do

      return
      end
