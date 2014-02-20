c This is derived from get1ddeepspec.f but it returns
c the arrays of wavelength, flux, rms, rather than doing a rebinning.
c Also, make the boxcar/optimal switch an input parameter.

c  Open a Berkeley-format DEEP 1-d spectrum BINTABLE fits file,
c  retrieve the spectra and inverse variance, 
c  (don't put onto linear wavelength scale, 
c  convert variance to std dev.
c  return number of data points, wave, spec and error, err flag

c input params: fname, ifopt (0=boxcar, 1=optimal)
      subroutine read1dspec(fname,ifopt,n,wavetot,spectot,espectot,ifer)

      include 'specarray.h'

      parameter(maxhdu=12,maxfield=128)
c are we reading 1/variance or 1/rms from the spec1d file?
      parameter(IFVAR=1)
c option to read the sky spectrum rather than the sky-subtracted
      parameter(IFREADSKY=0)
c get HDUs with spectra by name rather than number?
      parameter(IFHDUNAME=1)
      parameter(IFDEBUG=1)

      real spec(NWAVEMAX),espec(NWAVEMAX)
      character fname*(*)
      character hduname(2)*12
c      character hduname1*12,hduname2*12

      real wave(NWAVEMAX),var(NWAVEMAX)
      real wavetot(NWAVEMAX),spectot(NWAVEMAX),espectot(NWAVEMAX)
      integer indhdu(maxhdu)

      include 'pcredshift.h'

      ierr=0
      ifer=0

c open the file
      call ftgiou(im,ierr)
      call doerr(ierr)
      call ftopen(im,fname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         ifer = 1
         write(*,'("Error opening ",a)') fname
         go to 666
      end if

c HDU's to get spectra from.  Blue and red chips are separate
c In typical file, HDU 1 is common, HDUs 2 and 3 are boxcar
c  4 and 5 are optimal.
c  Note this is a problem if the red or blue side data does not 
c  exist and the HDUs are omitted
      nhdu = 2
c ifopt switches boxcar/optimal
      if (ifopt .ne. 0) then
         indhdu(1) = 4
         indhdu(2) = 5
      else
         indhdu(1) = 2
         indhdu(2) = 3
      end if

c Try to get the HDUs by name rather than number.  
      if (IFHDUNAME .eq. 1) then
         if (ifopt .ne. 0) then
            hduname(1) = 'Horne-B'
            hduname(2) = 'Horne-R'
         else
            hduname(1) = 'Bxspf-B'
            hduname(2) = 'Bxspf-R'
         end if
         if (IFDEBUG .eq. 1) then
            write (*,'(a20,a20)') hduname(1),hduname(2)
         end if
c move to specified HDU.  the 2 specifies it's a binary table
         call ftmnhd(im,2,hduname(1),0,ierr)
         if (ierr .ne. 0) then
            indhdu(1) = 99
            call doerr(ierr)
         else
c retrieve number of current HDU             
            call ftghdn(im,indhdu(1))
         end if
         call ftmnhd(im,2,hduname(2),0,ierr)
         if (ierr .ne. 0) then
            indhdu(2) = 99
            call doerr(ierr)
         else
            call ftghdn(im,indhdu(2))
         end if
         if (IFDEBUG .eq. 1) then
            write(*,'(a12,2x,i3)') hduname(1),indhdu(1)
            write(*,'(a12,2x,i3)') hduname(2),indhdu(2)
         end if
      end if

c npread is total number of points read
      npread = 0
      do ihdu=1,nhdu
         if (IFREADSKY .eq. 0) then
            call readfield(im,indhdu(ihdu),'SPEC',itype1,isize1,
     $           naxis,laxes,spec)
         else
            call readfield(im,indhdu(ihdu),'SKYSPEC',itype1,isize1,
     $           naxis,laxes,spec)
         end if
         call readfield(im,indhdu(ihdu),'IVAR',itype2,isize2,
     $        naxis,laxes,espec)
         call readfield(im,indhdu(ihdu),'LAMBDA',itype3,isize3,
     $        naxis,laxes,wave)
c  We could also consider reading the mask and knocking out
c  bad pixels.

         if (isize1 .ne. isize3 .or. isize2 .ne. isize3) then
            write(*,*) "vector lengths dont match ",isize1,isize2,isize3
         end if
         if (itype1 .ne. 42 .or. itype2 .ne. 42 .or. 
     $        itype3 .ne. 42) then
            write(*,*) "Non-float type, problems? ",itype1,itype2,itype3
         end if
         nstart = npread+1
         nend = nstart+isize1
         npread = npread+isize1
c         write(*,*) indhdu(ihdu),nstart,nend,npread
c unify the spectra from all HDUs read
         do i=1,isize1
            spectot(i+nstart-1) = spec(i)
            espectot(i+nstart-1) = espec(i)
            wavetot(i+nstart-1) = wave(i)
         end do
      end do
      call ftclos(im,ierr)
      call doerr(ierr)
      call ftfiou(im,ierr)
      call doerr(ierr)

      n = npread
c      write(*,'("read1dspec read ",i6," pixels")') n

c I assume that the blue and red spectra both have wavelength increasing
c and don't overlap.  Could insert a sort on wavelength here.

c This is new and hopefully correct
      do i=1,n
         if (espectot(i) .gt. 1.0e-10) then
            if (IFVAR .ne. 0) then
c we read 1/variance
               espectot(i) = 1./sqrt(espectot(i))
            else
c we read 1/rms?
               espectot(i) = 1./espectot(i)
            end if
         else
            spectot(i) = badset
            espectot(i) = badset
         end if
      end do

 666  return
      end
