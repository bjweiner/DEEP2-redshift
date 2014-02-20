
c  bad-flag specified spectral regions in the array

c blankoutsky2 takes out the B-band as well as 5577 and the A-band

c could also blank locations of possible strong em lines, mostly
c 3727 and 5007, if we are using this for abs line fits.  Or use
c selectabsregions to exclude those.

c blankoutlines is directly from blankoutsky2, but it also
c needs to have the z as input so it can blank in restframe

      subroutine blankoutlines(specarray,nspec,nspmax,nwave,wave,z)

      real specarray(nspmax,nwave)
      real wave(nwave)
      
      parameter (NTOBLANK=11)
      real blankreg(2*NTOBLANK)

c      data blankreg /5574.,5580.,7590.,7650./
c      data blankreg /5574.,5580.,6855.,6870.,7590.,7650./

c These might need to be in vacuum instead of air
c 3726+3729,4101,4340,4861,4959,5007,6548,6563,6583,6717,6731
cc medium size regions
c      data blankreg /3717.,3738.,4095.,4107.,4333.,4347.,4851.,4871.,
c     $     4954.,4964.,4997.,5017.,6545.,6551.,6550.,6576.,
c     $     6578.,6588.,6712.,6722.,6727.,6735./
c large size regions
      data blankreg /3707.,3747.,4091.,4111.,4325.,4355.,4841.,4881.,
     $     4949.,4969.,4977.,5037.,6540.,6556.,6523.,6603.,
     $     6563.,6603.,6707.,6727.,6721.,6741./

      include 'pcredshift.h'

      nblank=NTOBLANK

c I wrote the regions in air, so need to shift to vacuum?
c vacuum wavelengths are longer.
c Only do this if the spectrum was not converted to air on reading in.
c      do ii=1,nblank
c         blankreg(ii) = blankreg(ii) * 1.00029
c      end do

      do ii=1,nblank
      wmin = blankreg(2*ii-1) * (1.0+z)
      wmax = blankreg(2*ii) * (1.0+z)

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

