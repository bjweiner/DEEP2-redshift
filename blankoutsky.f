
c  bad-flag specified spectral regions in the array

      subroutine blankoutsky(specarray,nspec,nspmax,nwave,wave)

      real specarray(nspmax,nwave)
      real wave(nwave)
      
      parameter (NTOBLANK=2)
      real blankreg(2*NTOBLANK)

      data blankreg /5574.,5580.,7590.,7650./
c 5577: 5574-5580, B-band: 6855-6870 ?, A-band: 7590-7650 ?

      include 'pcredshift.h'

      nblank=NTOBLANK

      do ii=1,nblank
      wmin = blankreg(2*ii-1)
      wmax = blankreg(2*ii)

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
         do i=imin,imax
            specarray(j,i) = badset
         end do
      end do

c end the do ii loop
      end do

      return
      end

