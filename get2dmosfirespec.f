
c open a 2d spectrum file from the format Jon Trump put MOSFIRE
c data in, slightly different from DEIMOS
c this might work for any single hdu fits binary table with
c fields named FLUX, LAMBDA, IVAR

c for DEEP2 2d spec:
c open the two slit2d files for a given slit and return
c one 8192 or two 4096-long arrays of wave, flux, error, 
c depending on whether get2ddeepspec1 or get2ddeepspec is called

c ROTCURVE calls get2ddeepspec1 so leave that name unchanged
c and just use it as a wrapper to call get2dmosfirespec

c Here the second file name will be a dummy argument

      subroutine get2ddeepspec1(sfile1,sfile2,nx,ny,wave,
     $     flux,ferr)

      parameter (NMAXSLIT=256)
      parameter (NMAXTWO=8192)
      parameter (NMAXDISP=4096)

      parameter (IFDEBUG=1)

      real wave(NMAXTWO,NMAXSLIT),flux(NMAXTWO,NMAXSLIT)
      real ferr(NMAXTWO,NMAXSLIT)
c      real wave1(NMAXDISP,NMAXSLIT),wave2(NMAXDISP,NMAXSLIT)
c      real flux1(NMAXDISP,NMAXSLIT),flux2(NMAXDISP,NMAXSLIT)
c      real ferr1(NMAXDISP,NMAXSLIT),ferr2(NMAXDISP,NMAXSLIT)
      integer nx,ny
      character sfile1*120,sfile2*120

      call get2dmosfirespec(sfile1,nx,ny,wave,
     $     flux,ferr)
c      call get2ddeepspec(sfile1,sfile2,nx,ny,wave1,wave2,
c     $     flux1,flux2,ferr1,ferr2)
c check for errors 
      if (nx .le. 0 .or. nx .ge. 100000 
     +     .or. ny .le. 0 .or. ny .ge. 100000) then
         write(*,*) "get2ddeepspec1: bad image sizes ",nx,ny
         return
      end if

c      if (IFDEBUG .ge. 1) write(*,*) "get2ddeepspec1: combining images"

      return
      end

cccccccccccccccccccc
c unlike get1ddeepspec, doesn't linearize, or fill gap

      subroutine get2dmosfirespec(sfile1,nx,ny,wave1,
     $     flux1,ferr1)

c maximum number of pixels along slit
      parameter(NMAXSLIT=256)
      parameter(NMAXDISP=8192)

      parameter (IFDEBUG=1)
     
      real wave1(NMAXDISP,NMAXSLIT)
      real flux1(NMAXDISP,NMAXSLIT)
      real ferr1(NMAXDISP,NMAXSLIT)
      integer nx,ny
      character sfile1*120,sfile2*120

      real specdata(NMAXSLIT*NMAXDISP)
      real basewave1(NMAXDISP),basewave2(NMAXDISP)
      integer laxes1(5),laxes2(5)
      integer lwaxis1(5),lwaxis2(5)

      maxcol = NMAXDISP
      maxrow = NMAXSLIT
      ierr = 0

      if (IFDEBUG .ge. 2) write(*,*) "get2dmosfirespec: starting"
      nimag=0
      ifokblue=0
      call ftgiou(im1,ierr)
      if (IFDEBUG .ge. 2) write(*,*) "get2dmosfirespec image1 ",im1,ierr
      call doerr(ierr)
      call ftopen(im1,sfile1,0,block,ierr)
      if (ierr .ne. 0) then
         write(*,'(a,a)') "Couldnt open image ",sfile1
      else
         nimag=nimag+1
         ifokblue=1
         if (IFDEBUG .ge. 2) 
     $        write(*,*) "get2dmosfirespec: opened image 1",im1
      end if
      call doerr(ierr)
      ifokred=0
c      call ftgiou(im2,ierr)
c      if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec image2 ",im2,ierr
c      call doerr(ierr)
c      call ftopen(im2,sfile2,0,block,ierr)
c      if (ierr .ne. 0) then
c         write(*,'(a,a)') "Couldnt open image ",sfile2
c      else
c         nimag=nimag+1
c         ifokred=1
c         if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec: opened red",im2
c      end if
c      call doerr(ierr)

      if (nimag .eq. 0) then
c bail
         write(*,*) "Found no images!"
         nx = 0
         ny = 0
         return
c         go to 400
      else if (nimag .eq. 1) then
         write(*,*) "get2dmosfirespec got 1 image"
c should do something        
      else if (nimag .eq. 2 .and. IFDEBUG .ge.1) then
         write(*,*) "Got 2 images"
      end if

c read the 2-d data from the sfiles, for mosfire, hdu 2 which is the
c first extension. hdu 1 is the primary. see readfield.f
      ihdu = 2

      if (ifokblue .eq. 1) then
       call readfield(im1,ihdu,'FLUX',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,flux1)
       call readfield(im1,ihdu,'IVAR',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,ferr1)
c LAMDBDA not DLAMBDA here
       call readfield(im1,ihdu,'LAMBDA',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,wave1)
c       call readfield(im1,ihdu,'LAMBDA0',itype1,iw0size1,
c     $     nw0axis1,lwaxis1,basewave1)
       nx1=laxes1(1)
       ny1=laxes1(2)
       if (IFDEBUG .ge. 1) 
     +      write(*,*) "get2ddeepspec image 1 size ",nx1,ny1
       call ftclos(im1,ierr)
       call doerr(ierr)
       call ftfiou(im1,ierr)
       call doerr(ierr)
      else
         do j=1,NMAXSLIT
            do i=1,NMAXDISP
               wave1(i,j) = 0.0
               flux1(i,j) = 0.0
               ferr1(i,j) = 1.0e-8
            end do
         end do
      end if

c      if (ifokred .eq. 1) then
c       call readfield(im2,ihdu,'FLUX',itype2,isize2,
c     $     naxis2,laxes2,specdata)
c       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,flux2)
c       call readfield(im2,ihdu,'IVAR',itype2,isize2,
c     $     naxis2,laxes2,specdata)
c       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,ferr2)
c       call readfield(im2,ihdu,'DLAMBDA',itype2,isize2,
c     $     naxis2,laxes2,specdata)
c       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,wave2)
c       call readfield(im2,ihdu,'LAMBDA0',itype2,iw0size2,
c     $     nw0axis2,lwaxis2,basewave2)
c       nx2=laxes2(1)
c       ny2=laxes2(2)
c       call ftclos(im2,ierr)
c       call doerr(ierr)
c       call ftfiou(im2,ierr)
c       call doerr(ierr)
c       if (IFDEBUG .ge. 1) 
c     +      write(*,*) "get2ddeepspec image 2 size ",nx2,ny2
c      else
c         do j=1,NMAXSLIT
c            do i=1,NMAXDISP
c               wave2(i,j) = 0.0
c               flux2(i,j) = 0.0
c               ferr2(i,j) = 1.0e-8
c            end do
c         end do
c      end if

c      write(*,*) "flux arrays ",isize1,isize2,naxis1,laxes1(1),laxes1(2)

c      nx = min(nx1,nx2)
c      ny = min(ny1,ny2)
      if (nx1 .le. 0 .or. nx1 .ge. 100000 
     +     .or. ny1 .le. 0 .or. ny1 .ge. 100000) then
         write(*,*) "bad image sizes on first image ",nx1,ny1
         ifbad1 = 1
      else
         ifbad1 = 0
      end if
c      if (nx2 .le. 0 .or. nx2 .ge. 100000 
c     +     .or. ny2 .le. 0 .or. ny2 .ge. 100000) then
c         write(*,*) "bad image sizes on second image ",nx2,ny2
c         ifbad2 = 1
c      else
c         ifbad2 = 0
c      end if

      if (ifbad1 .eq. 0 .and. ifbad2 .eq. 0) then
         nx = max(nx1,nx2)
         ny = max(ny1,ny2)
      else if (ifbad1 .eq. 0) then
         nx = nx1
         ny = ny1
c      else if (ifbad2 .eq. 0) then
c         nx = nx2
c         ny = ny2
      else
         write (*,*) "Image has bad sizes, punting!"
         nx = 0
         ny = 0
         return
      end if

      if (IFDEBUG .ge. 2) 
     $     write(*,*) "get2dmosfirespec: converting errors"
c      open(8,file='testget2d.out',status='unknown')
      do j=1,ny1
         do i=1,nx1
c            if (j .eq. 1) write(8,*) wave1(i,j),basewave1(i)
c Don't need basewave1 for MOSFIRE since it is LAMBDA not DLAMBDA
c            wave1(i,j) = wave1(i,j) + basewave1(i)
c convert IVAR to error
            if (ferr1(i,j) .gt. 1.0e-12 .and.
     $           ferr1(i,j) .lt. 1.0e10) then
               ferr1(i,j) = 1.0 / sqrt(ferr1(i,j))
            else
               ferr1(i,j) = 1.0e6
            end if
c check FLUX for bad values
            if (flux1(i,j) .lt. -100. .or. 
     $           flux1(i,j) .gt. 1.0e5) then
               flux1(i,j) = 0.0
            end if
         end do
      end do
c      close(8)
c      do j=1,ny2
c         do i=1,nx2
c            wave2(i,j) = wave2(i,j) + basewave2(i)
c            if (ferr2(i,j) .gt. 0.0) then
c               ferr2(i,j) = 1.0 / sqrt(ferr2(i,j))
c            else
c               ferr2(i,j) = 1.0e6
c            end if
c         end do
c      end do

      return
      end

cccccccccccccccccccccccccccccc

c from zchoose2..4

c  reformat a 1-d data array into proper 2-d
      subroutine reformat(n1,data1d,naxis,laxes,nxphys,nyphys,data2d)

      real data1d(n1),data2d(nxphys,nyphys)
      integer laxes(naxis)

      if (naxis .le. 0) then
         write(*,*) "Error in reformat, naxis <= 0"
         return
      else if (naxis .eq. 1) then
         nx = laxes(1)
         ny = 1
      else
         nx = laxes(1)
         ny = laxes(2)
      end if
c bounds checking
      if (nx .gt. nxphys) then
         write(*,'("Warning! 2d image x is ",i6," max ",i6)') nx,nxphys
         nx=nxphys
      end if
      if (ny .gt. nyphys) then
         write(*,'("Warning! 2d image y is ",i6," max ",i6)') ny,nyphys
         ny=nyphys
      end if

      do j=1,ny
         ioffset = (j-1)*nx
         do i=1,nx
            data2d(i,j) = data1d(ioffset+i)
         end do
      end do
      return
      end

