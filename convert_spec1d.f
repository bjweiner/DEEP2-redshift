
c convert a DEEP 2 spec1d spectrum into a linearized-wavelength
c fits image that iraf can read

c now writes a 2xN spectrum, row 1 is flux, row 2 is rms

      program convert_spec1d

      include 'specarray.h'

      real wave(NWAVEMAX),spec(NWAVEMAX),espec(NWAVEMAX)
      real fout(NWAVEMAX),fouterr(NWAVEMAX)
      integer isize(7)
      character fname*160,sname*200,cline*160,objid*20
      character fname2*160,oname*200

      ifoptextract=1

 130  continue
      write(*,'("list of 1-d spectra: ",$)')
      read(*,'(a)') fname
      open(4,file=fname,status='old',err=130)
 140  continue
      write(*,'("list of output filenames: ",$)')
      read(*,'(a)') fname2
      open(2,file=fname2,status='old',err=140)

      nobj=1
 300  continue
      read(4,'(a)',err=666,end=666) sname
      read(2,'(a)',err=666,end=666) oname
      call read1dspec(sname,ifoptextract,np,wave,spec,espec,ifer)
      if (ifer .ne. 0) then
         write(*,'(a,a)') "Error getting spectrum for ",sname
         go to 300
      end if
      write(*,'(a,a)') "Read    ",sname(1:72)
      write(*,'(a,a)') "Writing ",oname(1:72)

      nout = np
      call linrebinerr(np,wave,spec,espec,nout,w0,dw,fout,fouterr)

c linrebinerr sets any pixels without data to -1.e6
c This is mildly annoying for plotting with autoscaling like in
c IRAF splot or implot, so as a hack let's reset them to -5000.

      do i=1,nout
         if (fout(i) .lt. -9.99e5) fout(i) = -5000.
         if (fouterr(i) .lt. -9.99e5) fouterr(i) = -5000.
      end do
      
c create new image, write data into it,
c make headers with w0, dw, DISPAXIS?
      isize(1) = nout
      isize(2) = 2
      do i=3,7
         isize(i) = 1
      end do
c      call imcrea(oname,isize,1,6,ierr)
      call imcrea(oname,isize,2,6,ierr)
      call doerr(ierr)
      call imopen(oname,3,im1,ierr)
      call doerr(ierr)
      call imakwr(im1,"CRPIX1",1.0," ",ierr)
      call doerr(ierr)
      call imakwr(im1,"CRVAL1",w0+dw," ",ierr)
      call doerr(ierr)
      call imakwr(im1,"CDELT1",dw," ",ierr)
      call doerr(ierr)
      call imakwi(im1,"DISPAXIS",1," ",ierr)
      call doerr(ierr)
c      call impl1r(im1,fout,ierr)
      call impl2r(im1,fout,1,ierr)
      call doerr(ierr)
      call impl2r(im1,fouterr,2,ierr)
      call doerr(ierr)
      call imclos(im1,ierr)
      call doerr(ierr)
      
      go to 300

 666  continue
      close(4)
      end


