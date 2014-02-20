
c given a sky spectrum, convert it into a std.dev. spectrum
c nexp is the number of combined exposures, exptime the
c exp.time of each individual exp, dpix the number of 
c pixels averaged in the extracted spectrum
      subroutine scalenoise(skyspec,nw,nexp,exptime,dpix,espec)

      real skyspec(nw),espec(nw)

      include 'pcredshift.h'

c  Exposure time to which everything has been scaled
      reftime = 15.*60.
c  CCD parameters, gain in e-/DN, read noise in e-
      gain = 2.0
      rn = 7.0
      rnsq = rn*rn
c  Original pixels in the wavelength direction in one 
c  combined pixel.  For the DEEP 1 data, 600/5000 or 600/7500
c  gratings, there is 1.28 A/pix and in the 98-99 combined data
c  it is still 1.28 A/pix.
      wpixscale = 1.28 / 1.28

c Number of original pixels in one combined pixel
      rnorigpix = dpix * nexp * wpixscale
c 1 final DN = this many original total DN, takes out the scaling,
c the averaging of Nexp images, and averaging of rnorigpix pixels
      dnorig = exptime / reftime * rnorigpix
c 1 final DN = this many original electrons
      elecorig = dnorig * gain

c For nexp=2, exptime=1800, dpix=5, wpixscale=1, then
c  rnorigpix = 10
c  dnorig = 20
c  elecorig = 40
c  photons = 40*DN
c  shotnoisesq = 40*DN
c  readnoisesq = 490
c  totnoise = sqrt(40*DN + 490)
c  espec(i) = sqrt(40*DN + 490)/40 = sqrt(DN/40 + 0.306)
c which is much smaller than the real noise?

      do i=1,nw
         if (skyspec(i) .gt. bad) then
            photons = elecorig * skyspec(i)
            shotnoisesq = max(photons,1.e-4)
            readnoisesq = rnorigpix * rnsq
c totnoise is in original photons total
            totnoise = sqrt(shotnoisesq + readnoisesq)
c divide it down to final DN
            espec(i) = totnoise / elecorig
         else
            espec(i) = badset
         end if
      end do

      return
      end

