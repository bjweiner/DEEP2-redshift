c this version also retrieves fwhm and ypos, if it works
c it should replace the old version

c  Open a Berkeley-format DEEP 1-d spectrum BINTABLE fits file,
c  retrieve the spectra and inverse variance, put onto 
c  linear wavelength scale, convert variance to std dev.
c  return number of data points, spec and error, wave calib, err flag

      subroutine get1ddeepspec(fname,n,spec,espec,w0,dw,fwhm,ypos,ifer)

      include 'specarray.h'

      parameter(maxhdu=12,maxfield=128)
c are we reading 1/variance or 1/rms from the spec1d file?
      parameter(IFVAR=1)

      real spec(NWAVEMAX),espec(NWAVEMAX)
      character fname*(*)

      real wave(NWAVEMAX),var(NWAVEMAX)
      real wavetot(NWAVEMAX),spectot(NWAVEMAX),espectot(NWAVEMAX)
      integer indhdu(maxhdu)

      include 'pcredshift.h'

      ierr=0
      ifer=0
c HDU's to get spectra from.  Blue and red chips are separate
c In typical file, HDU 1 is common, HDUs 2 and 3 are boxcar?
c  4 and 5 are optimal?
      nhdu = 2
      indhdu(1) = 2
      indhdu(2) = 3

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

c npread is total number of points read
      npread = 0
      do ihdu=1,nhdu
c only retrieve fwhm and ypos from first hdu read
         if (ihdu .eq. 1) then 
            call readfield(im,indhdu(ihdu),'FWHM',itype3,isize3,
     $        naxis,laxes,fwhm)
            call readfield(im,indhdu(ihdu),'YPOS',itype3,isize3,
     $        naxis,laxes,ypos)
         end if
         call readfield(im,indhdu(ihdu),'SPEC',itype1,isize1,
     $        naxis,laxes,spec)
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
c      write(*,'("get1ddeepspec read ",i6," pixels")') n

c I assume that the blue and red spectra both have wavelength increasing
c and don't overlap.

c 10/29/02 - This was old and rebinning the variance is wrong!
cc  Going to deal with rebinning by assuming I can just 
cc  rebin the variance.
c      do i=1,n
c         if (espectot(i) .gt. bad) then
c            if (IFVAR .ne. 0) then
c               espectot(i) = 1. / espectot(i)
c            else
cc we read 1/rms from file
c               espectot(i) = 1. / espectot(i)**2
c            end if               
c         else
c            espectot(i) = badset
c         end if
c      end do
c
c      call linrebin(n,wavetot,spectot,n,w0,dw,spec)
c      call linrebin(n,wavetot,espectot,n,w0e,dwe,espec)

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

      call linrebinerr(n,wavetot,spectot,espectot,n,w0,dw,spec,espec)

c      write(*,*) "lin rebinned w0, dw = ",w0,dw

c could check to make sure that w0=w0e and dw=dwe, but that
c seems like overkill
cc 10/29/02 - no longer need to convert the variance to std dev      
c      do i=1,n
c         if(espec(i) .gt. bad) then
c            if (espec(i) .gt. 0.0) then
c               espec(i) = sqrt(espec(i))
c            else
cc  this shouldn't happen
c               espec(i) = badset
c               spec(i) = badset
c            end if
c         else
c            espec(i) = badset
c            spec(i) = badset
c         end if
c      end do

 666  return
      end
