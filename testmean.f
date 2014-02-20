
c test various routines

      program testplot

      include 'specarray.h'
      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname(NSPECMAX,15)
      common /specname/ specname
      character specname(NSPECMAX)*60

      real spec1(NWAVEMAX),spec2(NWAVEMAX),spec3(NWAVEMAX)
      real wave1(NWAVEMAX),wave2(NWAVEMAX),wave3(NWAVEMAX)
      character sname*60

      include 'pcredshift.h'

      call pgbeg(0,'?',2,2)
      call pgscf(2)
      call pgsch(1.5)

      call readspec(nspec,nwave)

      nw = nwave
      nwlog = nwave

c loop through some spectra
      do k=1,nspec

      do i=1,nwave
         spec1(i) = specarr(k,i)
      end do

      sname = specname(k)
c      w0=6001.28
c      dw=1.28
c      w0log = log10(6000.)
      w0 = specw0(k)
      dw = specdw(k)
      write(*,'(a,", w0= ",f9.3,", dw= ",f6.3)') sname,w0,dw
      w0log = log10(w0)
      dwlog = 1.e-4
      do i=1,nw
         wave1(i) = w0+i*dw
      end do
      do i=1,nwlog
         wave2(i) = w0log+i*dwlog
      end do

c do a plot
      call showspec(nw,wave1,spec1)
      call pglabel("wavelength (A)","counts",sname)

      call logrebin(spec1,nw,w0,dw,nwlog,w0log,dwlog,spec2)

      call showspec(nwlog,wave2,spec2)
      call pglabel('log wavelength','counts',sname)

c      call deshiftz(spec2,nwlog,dwlog,0.43)
      contwrad = 100.*dwlog
      call contsubmed(spec2,nwlog,dwlog,contwrad)

      call showspec(nwlog,wave2,spec2)
      call pglabel('log wavelength','continuum-subtracted counts',sname)

      call pgpage()
c end the do k loop
      end do

      end

