
c  bad-flag specified spectral regions in the array

c blankoutsky2 takes out the B-band as well as 5577 and the A-band

c could also blank locations of possible strong em lines, mostly
c 3727 and 5007, if we are using this for abs line fits.  Or use
c selectabsregions to exclude those.

      subroutine blankoutsky2(specarray,nspec,nspmax,nwave,wave)

      real specarray(nspmax,nwave)
      real wave(nwave)
      
      parameter (NTOBLANK=3)
      real blankreg(2*NTOBLANK)

c      data blankreg /5574.,5580.,7590.,7650./
c      data blankreg /5574.,5580.,6855.,6870.,7590.,7650./
      data blankreg /5574.,5580.,6840.,6890.,7590.,7650./
c 5577: 5574-5580, B-band: 6855-6870 ?, A-band: 7590-7650 ?
c or maybe 6840-6890 for B-band?

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

