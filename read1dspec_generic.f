
c This is a modification of read1dspec_deep1.f that is intended to
c read a generic 1-d spectrum from a FITS image.  I set it up to read the
c flux from row 1 and the rms error from row 2; modify as needed for your
c data format.  It assumes that you have a linear wavelength
c calibration in the FITS header using CRVAL1, CRPIX1, and
c CDELT1 or CD1_1 - see the getcalib() subroutine at the end.

c This is derived from get1ddeepspec.f but it returns
c the arrays of wavelength, flux, rms, rather than doing a rebinning.
c Also, make the boxcar/optimal switch an input parameter.
c  (for generic input spectrum, the boxcar/optimal switch does nothing)

c  (orig) Open a Berkeley-format DEEP 1-d spectrum BINTABLE fits file,
c  retrieve the spectra and inverse variance, 
c  (don't put onto linear wavelength scale, 
c  convert variance to std dev.
c  return number of data points, wave, spec and error, err flag

c input params: fname, ifopt (0=boxcar, 1=optimal)
      subroutine read1dspec_generic(fname,
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

c charspec was used to indicate whether reading a DEEP1 B or R spectrum.

c      indx = index(fname,' ') -1
c      if (charspec .eq. 'B') then
c         fname1 = fname(1:indx) // '/Be.fits'
c      else
c         fname1 = fname(1:indx) // '/Re.fits'
c      end if

      ierr=0
      ifer=0

c buffer of pixels to trim at each end, if there is not 
c a sharp transition from good data to bad flag
      nbuffer = 5

c open the files
      call ftgiou(im1,ierr)
      call doerr(ierr)
      call ftopen(im1,fname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         ifok = 0
         write(*,'("Error opening ",a)') fname
         call ftfiou(im1,ierr)
         call doerr(ierr)
         ifer = 1
         go to 666
      else
         call getcalib(im1,w0,dw,iferr)
         if (iferr .eq. 0) then
            ifok = 1
            write(*,'(f9.3,2x,f6.3,2x,a60)') w0,dw,fname
         else
            ifok = 0
            write(*,'("No wavecal for ",a60)') fname
            write(*,'("enter wave(pixel 0) and dwave/dpix: ",$)')
            read(*,*) w0, dw
         end if
      end if

c      if (ifbok .ne. 0) then
         call ftgknj(im1,'NAXIS',1,7,laxes,idim,ierr)
         call doerr(ierr)
         nx = laxes(1)
         ny = laxes(2)
c         write(*,*) "nx = ",nx
c Here we read flux from the file from pixel 1 to nx
         call ftgpve(im1,0,1,nx,0.0,spec,flag,ierr)
         call doerr(ierr)
c And here we try to read fluxerr from the file from pixel nx+1 to 2*nx
c if we can't find it, fake it to 1.0
         if (ny .ge. 2) then
            call ftgpve(im1,0,nx+1,2*nx,0.0,espec,flag,ierr)
         else
            write(*,*) "Couldnt read error from 2nd row of spectrum,",
     $       " need to change read1dspec_generic.f, setting errors to 1"
            do i=1,nx
               espec(i) = 1.0
            end do
         end if
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
            espectot(j) = espec(i)
         end do
c      end if

c      write(*,*) "Got ", np," points"

c return n
      n=np
c bad-flag any points with suspicious values
      do i=1,np
         if (spectot(i) .lt. -300. .or. spectot(i) .gt. 1.e4 .or.
     $        espectot(i) .lt. 0.001) then
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

