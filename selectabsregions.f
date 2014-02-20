
c  bad-flag specified spectral regions in the array

c selectabsregions is derived from blankoutsky2 but it
c eliminates data in all but a small number of regions
c in restframe w.l. (where there are useful absorption lines).
c it can be called with the observed wavelengths,
c so you also pass it the z to transform the restframe regions 
c to observed

c can use this to excludde regions around 3728 and 5007, etc

c if iflog is 1, then the input array is in log wavelength in A

      subroutine selectabsregions(specarray,nspec,nspmax,nwave,wave,z,
     $     iflog)

      real specarray(nspmax,nwave)
      real tmparray(nspmax,nwave)
      real wave(nwave)

c this is derived from blankoutsky so the regions are named
c "blankreg" even though we are selecting them.
      parameter (NTOBLANK=7)
      real blankreg(2*NTOBLANK)
      real blankreg2(2*NTOBLANK)

c      data blankreg /5574.,5580.,7590.,7650./
c      data blankreg /5574.,5580.,6855.,6870.,7590.,7650./
c 5577: 5574-5580, B-band: 6855-6870 ?, A-band: 7590-7650 ?
      data blankreg /3735.,4000.,4070.,4130.,4270.,4370.,4800.,5002.,
     $     5012.,5380.,5840.,5960.,6520.,6600./

      include 'pcredshift.h'

      nblank=NTOBLANK
      do ib = 1,2*nblank
         blankreg2(ib) = blankreg(ib)* (1.0 + z)
         if (iflog .eq. 1) blankreg2(ib) = log10(blankreg2(ib))            
      end do

c set the tmparray to be all bad, then we write the good data into it
      do j=1,nspec
         do i=1,nwave
            tmparray(j,i) = badset
         end do
      end do

      do ii=1,nblank
      wmin = blankreg2(2*ii-1)
      wmax = blankreg2(2*ii)

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

c only copy data in the region into the tmparray
      do j=1,nspec
         do i=imin,imax
            tmparray(j,i) = specarray(j,i)
         end do
      end do

c end the do ii loop
      end do

c copy the tmparray back into the specarray
      do j=1,nspec
         do i=1,nwave
            specarray(j,i) = tmparray(j,i)
         end do
      end do
      return
      end

