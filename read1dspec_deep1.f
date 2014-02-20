
c This is a modification of read1dspec.f but is intended to
c work on DEEP 1 1-d data, and fudges the errors somehow.

c In this one the fname passed in is actually the directory
c and we look for <fname>/Be.fits and <fname>/Re.fits

c This is derived from get1ddeepspec.f but it returns
c the arrays of wavelength, flux, rms, rather than doing a rebinning.
c Also, make the boxcar/optimal switch an input parameter.

c  Open a Berkeley-format DEEP 1-d spectrum BINTABLE fits file,
c  retrieve the spectra and inverse variance, 
c  (don't put onto linear wavelength scale, 
c  convert variance to std dev.
c  return number of data points, wave, spec and error, err flag

c input params: fname, ifopt (0=boxcar, 1=optimal)
      subroutine read1dspec_deep1(fname,charspec,
     $     ifopt,n,wavetot,spectot,espectot,ifer)

      include 'specarray.h'

      parameter(maxhdu=12,maxfield=128)
c are we reading 1/variance or 1/rms from the spec1d file?
      parameter(IFVAR=1)

c      real specb(NWAVEMAX),specr(NWAVEMAX)
c      real waveb(NWAVEMAX),waver(NWAVEMAX)
c      real w0b,dwb,refpixb,w0r,dwr,refpixr
      real w0,dw,refpix
      real spec(NWAVEMAX),espec(NWAVEMAX)
      character charspec*1
      character fname*(*)
      character fname1*200,comment*80
      character skyname*200
      real w0sky,dwsky,refpixsky

      real wave(NWAVEMAX),var(NWAVEMAX)
      real wavetot(NWAVEMAX),spectot(NWAVEMAX),espectot(NWAVEMAX)
      integer laxes(7)
      logical flag

      include 'pcredshift.h'

      indx = index(fname,' ') -1
      if (charspec .eq. 'B') then
         fname1 = fname(1:indx) // '/Be.fits'
      else
         fname1 = fname(1:indx) // '/Re.fits'
      end if

      ierr=0
      ifer=0

c buffer of pixels to trim at each end, since there tends not to 
c be a sharp transition from good data to bad flag
      nbuffer = 10

c setting this here is a terrible kludge
c      skyname='/net/kubala/c/bjw/deepspec/skyspec/bothesky.fits'
      skyname='/net/asturias/a/kubala/' // 
     $     'c/bjw/deepspec/skyspec/bothesky.fits'

c open the files
      call ftgiou(im1,ierr)
      call doerr(ierr)
      call ftopen(im1,fname1,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         ifok = 0
         write(*,'("Error opening ",a)') fname1
         call ftfiou(im1,ierr)
         call doerr(ierr)
         ifer = 1
         go to 666
      else
         call getcalib(im1,w0,dw,iferr)
         if (iferr .eq. 0) then
            ifok = 1
            write(*,'(f9.3,2x,f6.3,2x,a60)') w0,dw,fname1
         else
            ifok = 0
            write(*,'("No wcal for ",a60)') fname1
c  this is not actually right for pre-97 data
            if (charspec .eq. 'B') then
               w0 = 3600.0
               dw = 1.28
            else
               w0 = 6001.28
               dw = 1.28
            end if
         end if
      end if

c      if (ifbok .ne. 0) then
         call ftgknj(im1,'NAXIS',1,7,laxes,idim,ierr)
         call doerr(ierr)
         nx = laxes(1)
c         write(*,*) "nx = ",nx
         call ftgpve(im1,0,1,nx,0.0,spec,flag,ierr)
         call doerr(ierr)
         do i=1,nx
            wave(i) = w0 + i*dw
         end do
c         call findends(nx,wave,spec,imin,imax,wmin,wmax)
         call findends(nx,spec,imin,imax)
c trim the buffer a bit more beyond these ends
         imin = min(imin+nbuffer,nx)
         wmin = wave(imin)
         imax = max(imax-nbuffer,1)
         wmax = wave(imax)
         write(*,*) "ends ",imin,imax,wmin,wmax

         call ftclos(im1,ierr)
         call doerr(ierr)
         call ftfiou(im1,ierr)
         call doerr(ierr)
c      end if

c note I've assumed dw is >0            

c copy data into output array
      np = 0
c      if (ifok .eq. 1) then
         ioffset = imin - (np+1)
         np = np + imax-imin+1
         do i=imin,imax
            j=i-ioffset
            wavetot(j) = wave(i)
            spectot(j) = spec(i)
         end do
c      end if

c      write(*,*) "Got ", np," points"

c need to fake up the error spectrum somehow,
c hopefully by reading a file, and maybe scaling
c with exposure time and # of rows extracted?
      do i=1,np
         espectot(i) = 1.0
      end do

      call ftgiou(im1,ierr)
      call doerr(ierr)
      call ftopen(im1,skyname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         write (*,'("Error opening sky spec: ",a)') skyname
         do i=1,np
            espectot(i) = 1.0
         end do
         go to 666
      end if
        
      call getcalib(im1,w0sky,dwsky,iferr)
      if (iferr .ne. 0) then
         w0sky = 3600.0
         dwsky = 1.28
      end if
      call ftgknj(im1,'NAXIS',1,7,laxes,idim,ierr)
      call doerr(ierr)
      nxsky = laxes(1)
c      write(*,*) "nxsky = ",nxsky
      call ftgpve(im1,0,1,nxsky,0.0,spec,flag,ierr)
      call doerr(ierr)
      do i=1,nxsky
         wave(i) = w0sky + i*dwsky
      end do
      call ftclos(im1,ierr)
      call doerr(ierr)
      call ftfiou(im1,ierr)
      call doerr(ierr)

c      write(*,*) "Sky image ",nxsky,wave(1),wave(nxsky)

c interpolate sky onto wavelength scale of spectrum
      call lininterp(nxsky,wave,spec,np,wavetot,espectot)
c don't want to use linrebin because it chooses its own w0 and dw
c      call linrebin(nxsky,wave,spec,nout,w0out,dwout,spectot)

c in principle we should calculate the noise before interpolating.
c probably not that big a deal esp. since read noise is not dominant.
c also, I found that when rebinning, you want to rebin the variance
c linearly (same way you rebin flux) otherwise things get ugly.

c calculate error from the sky values.
c sky file is in DN/1500sec, for 3600? sec total exp, 7 rows averaged
c we could try to scale this for each individual spectrum's 
c exp time and #rows averaged, but don't bother for now
c set a floor for error of 1.0 
c convert to photons/3600sec, compute error, then convert back
c to DN/1500 sec.
      gain = 2.0
      rnsq = 6.5**2
      elecsky = 3600./1500.*gain
      rowdiam = 7.0
      do i=1,np
         espectot(i) = max(1.0,(
     $        sqrt(espectot(i)*elecsky*rowdiam + rowdiam*rnsq)
     $        / rowdiam / elecsky) )
      end do

      n=np
c bad-flag any points with suspicious values
      do i=1,np
         if (spectot(i) .lt. -300. .or. spectot(i) .gt. 1.e4 .or.
     $        espectot(i) .lt. 0.1) then
            spectot(i) = badset
         end if
      end do

 666  continue

      return
      end

cccccccccccccccccccccccccccccc
c return wavelength of pixel 0 and dw
      subroutine getcalib(im,w0,dw,iferr)

      integer im,iferr
      real w0,dw,refpix
      character comment*80

      ierr=0

      call ftgkye(im,'CRPIX1',refpix,comment,ierr)
      call doerr(ierr)
      call ftgkye(im,'CRVAL1',w0,comment,ierr)
      call doerr(ierr)
      call ftgkye(im,'CDELT1',dw,comment,ierr)
      if (ierr .ne. 0) then
         ierr=0
         call ftgkye(im,'CD1_1',dw,comment,ierr)
      end if
      call doerr(ierr)
      if (w0 .gt. 0.0 .and. abs(dw) .gt. 1.0e-3) then
         w0 = w0 - refpix*dw
         iferr = 0
      else
         iferr = 1
      end if

      return
      end

