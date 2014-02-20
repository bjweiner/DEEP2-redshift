
c open the two slit2d files for a given slit and return
c one 8192 or two 4096-long arrays of wave, flux, error, 
c depending on whether get2ddeepspec1 or get2ddeepspec is called

      subroutine get2ddeepspec1(sfile1,sfile2,nx,ny,wave,
     $     flux,ferr)

      parameter (NMAXSLIT=512)
      parameter (NMAXTWO=8192)
      parameter (NMAXDISP=4096)

      parameter (IFDEBUG=1)

      real wave(NMAXTWO,NMAXSLIT),flux(NMAXTWO,NMAXSLIT)
      real ferr(NMAXTWO,NMAXSLIT)
      real wave1(NMAXDISP,NMAXSLIT),wave2(NMAXDISP,NMAXSLIT)
      real flux1(NMAXDISP,NMAXSLIT),flux2(NMAXDISP,NMAXSLIT)
      real ferr1(NMAXDISP,NMAXSLIT),ferr2(NMAXDISP,NMAXSLIT)
      integer nx,ny
      character sfile1*120,sfile2*120
     
      call get2ddeepspec(sfile1,sfile2,nx,ny,wave1,wave2,
     $     flux1,flux2,ferr1,ferr2)
c check for errors 
      if (nx .le. 0 .or. nx .ge. 100000 
     +     .or. ny .le. 0 .or. ny .ge. 100000) then
         write(*,*) "get2ddeepspec1: bad image sizes ",nx,ny
         return
      end if

      if (IFDEBUG .ge. 1) write(*,*) "get2ddeepspec1: combining images"

      do j=1,ny
         do i=1,NMAXDISP
            wave(i,j) = wave1(i,j)
            flux(i,j) = flux1(i,j)
            ferr(i,j) = ferr1(i,j)
            wave(i+NMAXDISP,j) = wave2(i,j)
            flux(i+NMAXDISP,j) = flux2(i,j)
            ferr(i+NMAXDISP,j) = ferr2(i,j)
         end do
      end do
c this is sort of dumb since we are hardcoding 4096 everywhere
      nx = nx+NMAXDISP

      return
      end

cccccccccccccccccccc
c unlike get1ddeepspec, doesn't linearize, or fill gap

      subroutine get2ddeepspec(sfile1,sfile2,nx,ny,wave1,wave2,
     $     flux1,flux2,ferr1,ferr2)

c maximum number of pixels along slit
      parameter(NMAXSLIT=512)
      parameter(NMAXDISP=4096)

      parameter (IFDEBUG=1)
     
      real wave1(NMAXDISP,NMAXSLIT),wave2(NMAXDISP,NMAXSLIT)
      real flux1(NMAXDISP,NMAXSLIT),flux2(NMAXDISP,NMAXSLIT)
      real ferr1(NMAXDISP,NMAXSLIT),ferr2(NMAXDISP,NMAXSLIT)
      integer nx,ny
      character sfile1*120,sfile2*120

      real specdata(NMAXSLIT*NMAXDISP)
      real basewave1(NMAXDISP),basewave2(NMAXDISP)
      integer laxes1(5),laxes2(5)
      integer lwaxis1(5),lwaxis2(5)

      maxcol = NMAXDISP
      maxrow = NMAXSLIT
      ierr = 0

      if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec: starting"
      nimag=0
      ifokblue=0
      call ftgiou(im1,ierr)
      if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec image1 ",im1,ierr
      call doerr(ierr)
      call ftopen(im1,sfile1,0,block,ierr)
      if (ierr .ne. 0) then
         write(*,'(a,a)') "Couldnt open image ",sfile1
      else
         nimag=nimag+1
         ifokblue=1
         if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec: opened blue",im1
      end if
      call doerr(ierr)
      ifokred=0
      call ftgiou(im2,ierr)
      if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec image2 ",im2,ierr
      call doerr(ierr)
      call ftopen(im2,sfile2,0,block,ierr)
      if (ierr .ne. 0) then
         write(*,'(a,a)') "Couldnt open image ",sfile2
      else
         nimag=nimag+1
         ifokred=1
         if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec: opened red",im2
      end if
      call doerr(ierr)

      if (nimag .eq. 0) then
c bail
         write(*,*) "Found no images!"
         nx = 0
         ny = 0
         return
c         go to 400
      else if (nimag .eq. 1) then
         write(*,*) "Only got 1 image!"
c should do something        
      else if (nimag .eq. 2 .and. IFDEBUG .ge.1) then
         write(*,*) "Got 2 images"
      end if

c read the 2-d data from the sfiles, hdu 2
      ihdu = 2

      if (ifokblue .eq. 1) then
       call readfield(im1,ihdu,'FLUX',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,flux1)
       call readfield(im1,ihdu,'IVAR',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,ferr1)
       call readfield(im1,ihdu,'DLAMBDA',itype1,isize1,
     $     naxis1,laxes1,specdata)
       call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,wave1)
       call readfield(im1,ihdu,'LAMBDA0',itype1,iw0size1,
     $     nw0axis1,lwaxis1,basewave1)
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

      if (ifokred .eq. 1) then
       call readfield(im2,ihdu,'FLUX',itype2,isize2,
     $     naxis2,laxes2,specdata)
       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,flux2)
       call readfield(im2,ihdu,'IVAR',itype2,isize2,
     $     naxis2,laxes2,specdata)
       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,ferr2)
       call readfield(im2,ihdu,'DLAMBDA',itype2,isize2,
     $     naxis2,laxes2,specdata)
       call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,wave2)
       call readfield(im2,ihdu,'LAMBDA0',itype2,iw0size2,
     $     nw0axis2,lwaxis2,basewave2)
       nx2=laxes2(1)
       ny2=laxes2(2)
       call ftclos(im2,ierr)
       call doerr(ierr)
       call ftfiou(im2,ierr)
       call doerr(ierr)
       if (IFDEBUG .ge. 1) 
     +      write(*,*) "get2ddeepspec image 2 size ",nx2,ny2
      else
         do j=1,NMAXSLIT
            do i=1,NMAXDISP
               wave2(i,j) = 0.0
               flux2(i,j) = 0.0
               ferr2(i,j) = 1.0e-8
            end do
         end do
      end if

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
      if (nx2 .le. 0 .or. nx2 .ge. 100000 
     +     .or. ny2 .le. 0 .or. ny2 .ge. 100000) then
         write(*,*) "bad image sizes on second image ",nx2,ny2
         ifbad2 = 1
      else
         ifbad2 = 0
      end if

      if (ifbad1 .eq. 0 .and. ifbad2 .eq. 0) then
         nx = max(nx1,nx2)
         ny = max(ny1,ny2)
      else if (ifbad1 .eq. 0) then
         nx = nx1
         ny = ny1
      else if (ifbad2 .eq. 0) then
         nx = nx2
         ny = ny2
      else
         write (*,*) "Both images have bad sizes, punting!"
         nx = 0
         ny = 0
         return
      end if

      if (IFDEBUG .ge. 2) write(*,*) "get2ddeepspec: converting errors"
c      open(8,file='testget2d.out',status='unknown')
      do j=1,ny1
         do i=1,nx1
c            if (j .eq. 1) write(8,*) wave1(i,j),basewave1(i)
            wave1(i,j) = wave1(i,j) + basewave1(i)
            if (ferr1(i,j) .gt. 0.0) then
               ferr1(i,j) = 1.0 / sqrt(ferr1(i,j))
            else
               ferr1(i,j) = 1.0e6
            end if
         end do
      end do
c      close(8)
      do j=1,ny2
         do i=1,nx2
            wave2(i,j) = wave2(i,j) + basewave2(i)
            if (ferr2(i,j) .gt. 0.0) then
               ferr2(i,j) = 1.0 / sqrt(ferr2(i,j))
            else
               ferr2(i,j) = 1.0e6
            end if
         end do
      end do

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

