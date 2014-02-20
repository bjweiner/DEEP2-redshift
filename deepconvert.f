
c convert a 1-d spectrum fits binary table
c into a fits image with spectrum and error
c linearized in wavelength

      program deepconvert

      parameter (nmax=20000)

      real flux1(nmax),ferr1(nmax)
      real w0,dw
      integer laxes(2)
      character sname*160,oname*160

      ierr=0

 100  continue
      write(*,'("spec1d table file [quit]: ",$)')
      read(*,'(a)') sname
      if (sname(1:3) .eq. '   ') go to 666

      call get1ddeepspec(sname,n1,flux1,ferr1,w0,dw,ifer)
      if (ifer .ne. 0) then
         write(*,*) "Error reading 1-d spectrum"
         go to 100
      end if
      
 200  continue
      write(*,'("name for new fits image: ",$)')
      read(*,'(a)') oname

      call ftgiou(im1,ierr)
      call doerr(ierr)
      iblock=1
      call ftinit(im1,oname,iblock,ierr)
      if (ierr .ne. 0) then
         write(*,*) "Error opening new file ",oname
         call doerr(ierr)
         go to 200
      end if
c bitpix=-32 for 32 bit floating point
      ibitpix=-32
      naxis=2
      laxes(1)=n1
      laxes(2)=2
      write(*,'("npix, w0, dw = ",i5,1x,f7.2,1x,f6.3)') n1,w0,dw
c      call ftphpr(im1,.true.,ibitpix,naxis,laxes,0,1,.true.,ierr)
      call ftphps(im1,ibitpix,naxis,laxes,ierr)
      call doerr(ierr)

c write spectrum to row 1 and error to row 2.
c This is ridiculously primitive
      nfpixel=1
      call ftppre(im1,1,nfpixel,n1,flux1,ierr)
      call doerr(ierr)
      nfpixel=1+n1
      call ftppre(im1,1,nfpixel,n1,ferr1,ierr)
      call doerr(ierr)

c write header values for wavelength calibration
      call ftpkyj(im1,'CRPIX1',0,'reference pixel',ierr)
      call doerr(ierr)
      call ftpkye(im1,'CRVAL1',w0,-8,'wavelength',ierr)
      call doerr(ierr)
      call ftpkye(im1,'CDELT1',dw,-7,'dw',ierr)
      call doerr(ierr)

c We should really copy the header out of the binary table
c and stick it in the fits image

c close output image and go around for another
      call ftclos(im1,ierr)
      call ftfiou(im1,ierr)
      call doerr(ierr)
      
      go to 100

 666  continue

      end

