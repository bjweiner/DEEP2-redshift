
c  read a list of images, y coord and possible redshifts
c  and plot postage stamps around locations of features
c  allow the user to pick the correct z by typing a number
c  (or hitting the cursor?)

      program zchoose

      parameter (maxcol=10000,maxrow=128,maxlines=15,maxz=8)
      parameter(IFSTAT=0)

      real spec2d(maxcol,maxrow)
c      real spec1d(maxcol)
      real tmpdat(maxcol)

      real wave(maxcol),flux(maxcol)
      real wlem(maxlines),wlabs(maxlines),wluv(maxlines)
      real wlplot(maxlines)
      real zcand(maxz)
      real trans(6)
      real cheight
      integer isize(7)

      integer function pgopen
      integer function pgcurs

      character fname*80,sname*80,cline*80,answ*1
      character toplabel*80,xlabel*10
      
      data wlem/1215.67, 1550.0, 2800.0, 3728.0, 3868.74, 4340.5,
     $  4861.3, 4959.9, 5006.8, 6562.8, 6583.4, 6716.5,
     $     0.0, 0.0, 0.0/
      data wlabs/1550.0, 2800.0, 3798.6, 3835.6, 3889.1, 3933.7,
     $  3968.5, 4101.75, 4304.4, 4340.5, 0.0, 0.0,
     $     0.0, 0.0, 0.0/
      data wluv/1215.67, 1304.37, 1334.53, 1526.71, 1548.2, 1550.77,
     $          1608.45, 1670.79, 2344.21, 2374.46, 2382.77, 2586.65,
     $     2600.17, 2796.35, 2803.53/

c max number of candidate z's to read
      nzmax = 5
      if (nzmax .gt. maxz) nzmax=maxz
      
c number of lines to actually plot
      nlplot = 7
      wlplot(1) = 3933.7
      wlplot(2) = 3968.5
      wlplot(3) = 3728.0
      wlplot(4) = 4861.3
      wlplot(5) = 4959.9
      wlplot(6) = 5006.8
      wlplot(7) = 6562.8

c wavelength and y radius to plot
c      wrad = 10.0
      wrad = 15.0
      yrad = 25.0
      iyrad = int(yrad)
c      nyp = int(2*yrad)

c coord transf for image plotting
      trans(1) = 0.
      trans(2) = 1.
      trans(3) = 0.
      trans(4) = 0.
      trans(5) = 0.
      trans(6) = 1.

      call pgbeg(0,'?',nlplot,nzmax)
c      call pgbeg(0,'?',1,1)
      call pgscf(1)
      call pgsch(1.5)
      call pgask(.false.)

c      call pgqcir(ilo,ihi)
c      write(*,*) ilo,ihi
      call pgscir(16,254)

 100  continue

 110  write(*,'("File with files: ",$)') 
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') stop
      open(2,file=fname,status='old',err=110)
 120  write(*,'("File with y-coords: ",$)') 
      read(*,'(a)') fname
c      if (fname(1:3) .eq. '   ') stop
      open(3,file=fname,status='old',err=120)
 130  write(*,'("File with candidate redshifts: ",$)') 
      read(*,'(a)') fname
c      if (fname(1:3) .eq. '   ') stop
      open(4,file=fname,status='old',err=130)

      open(10,file='zchoose.out',status='unknown',err=140)
      go to 145
 140  continue
      write(*,'("couldnt open zchoose.out; name for output file: ",$)')
      read (*,'(a)') fname
      open(10,file=fname,status='new',err=140)
 145  continue

c main loop through objects
 200  continue
      read(2,'(a)',end=666) sname
      read(3,*,end=666) y
      iy = nint(y)
      read(4,'(a)',end=666) cline

      nz=nzmax
c  readz resets nz
c      call readz(cline,nz,zcand)
c  just use the version that requires nz in the file for now
      call readz2(cline,nz,zcand)

c open image
      call imopen(sname,1,img,ier)
      if (ier .ne. 0) then
         write(*,*) 'Couldnt open image ',sname
         go to 200
      end if
      call imgsiz(img,isize,idim,itype,ier)
      nx = isize(1)
      ny = isize(2)

c  try to get wavelength info from header
      dw=0.0
      call imgkwr(img,'CRPIX1',refpix,ier)
      if (ier .ne. 0) refpix=0.0
      call imgkwr(img,'CRVAL1',w0,ier)
      call imgkwr(img,'CDELT1',dw,ier)
      if (ier .ne. 0 .or. dw .lt. 0.001) then
         call imgkwr(img,'CD1_1',dw,ier)
      end if
      w0 = w0 - dw*refpix
      write(*,'(a40,2x,f7.2,2x,f5.3)') sname,w0,dw
      xrad = wrad /dw
      ixrad = int(xrad)

      iy1 = max(1,iy - iyrad)
      iy2 = min(ny,iy + iyrad)
      nyp = iy2-iy1+1
      iycen = iyrad+1
      do j=1,nyp
         jimage=j+iy1-1
c         call imgl2r(img,spec2d(j,1),jimage,ier)
         call imgl2r(img,tmpdat,jimage,ier)
         do i=1,nx
            spec2d(i,j) = tmpdat(i)
         end do
      end do
      call imclos(img,ier)

c      call findends(nx,spec2d(iycen,1),imin,imax)
c      wmin = w0+imin*dw
c      wmax = w0+imax*dw

      fmin=-3.
      fmax=10.

 300  continue
c make a grid of postage stamps

c      call pgsubp(nlplot,nzmax)
c      call pgpage()
      do jpan=1,nz
c         call pgpanl(1,jpan)
c         write(toplabel,'(i2,"  z= ",f7.4)') jpan,zcand(jpan)
c         call pgmtxt('T',1.,0.,0.,toplabel)
cc         ifirst=1
         do ipan=1,nlplot
c            call pgpanl(ipan,jpan)
            wcen = wlplot(ipan) * (1. + zcand(jpan))
            ixcen = int((wcen-w0)/dw)
            rix1=real(ixcen-ixrad)
            rix2=real(ixcen+ixrad)
            riy1=real(iycen-iyrad)
            riy2=real(iycen+iyrad)
c square pixels with JUST (5th arg) = 1
            call pgenv(rix1,rix2,riy1,riy2,1,-1)
            call pgqch(cheight)
            call pgsch(4.0)
c label number and z over first panel
            if (ipan .eq. 1) then
               write(toplabel,'(i2,", z= ",f7.4)') jpan,zcand(jpan)
               call pgmtxt('T',1.4,0.5,0.5,toplabel)
            end if            
c label wl under each panel
            write(xlabel,'(f7.1)') wlplot(ipan)
            call pgmtxt('B',1.4,0.5,0.5,xlabel)
c label y and image  over panel (2,1) and (4,1)
            if (ipan .eq. 2 .and .jpan .eq. 1) then
               write(toplabel,'("y = ",i5)') iy
               call pgmtxt('T',1.4,0.5,0.5,toplabel)
            end if               
            if (ipan .eq. nlplot .and .jpan .eq. 1) then
               call pgmtxt('T',1.4,1.0,1.0,sname)
            end if
            call pgsch(cheight)
c bounds checking
c  if totally off the plot, skip pggray
            if (rix1 .lt. real(nx) .and. rix2 .gt. 1.0 .and. 
     $           riy1 .lt. real(nyp) .and. riy2 .gt. 1.0) then

c            call pgpanl(ipan,jpan)
c            call pgenv(rix1,rix2,riy1,riy2,1,-1)

c            if (ifirst .eq. 1) then
c               write(toplabel,'(i2,"  z= ",f7.4)') jpan,zcand(jpan)
c               call pgmtxt('T',1.,0.,0.,toplabel)
c               ifirst=0
c            end if
c trim plot region
               ix1=max(rix1,1.)
               ix2=min(rix2,real(nx))
               iy1=max(riy1,1.)
               iy2=min(riy2,real(nyp))
c could determine flux limits for each plot
               if (IFSTAT .ne. 0) then
                  npix=0
                  sum=0.
                  sumsq=0.0
                  smin=1.e10
                  smax=-1.e10
                  do iy=iy1,iy2
                     do ix=ix1,ix2
                        npix=npix+1
                        d=spec2d(ix,iy)
                        sum = sum+d
                        sumsq=sumsq+d*d
                        smin=min(smin,d)
                        smax=max(smax,d)
                     end do
                  end do
                  write(*,*) smin,smax,sum/n,sqrt(sumsq-sum*sum/n)
               end if
               call pggray(spec2d,maxcol,maxrow,ix1,ix2,iy1,iy2,
     $              fmax,fmin,trans)
            end if
         end do
      end do
      
c  this is a total kludge, forced because pgsubp and pgpanl just 
c  do not work right.
      do j=nz+1,nzmax
         do i=1,nlplot
            call pgpage()
         end do
      end do

 400  continue
      write(*,'("# of right z (0-none,f-change scale): ",$)') 
      read(*,'(a)') cline
      if (cline(1:1) .eq. 'f') then
         write(*,'("fmin, fmax for new plot: ",$)')
         read(*,*) fmin, fmax
         go to 300
      else if (cline(1:1) .eq. 'z') then
         write(*,'("enter z: ",$)')
         read(*,*) zcand(nz)
         go to 300
      else if (cline(1:1) .eq. 't') then
         write(*,'(a)') "Tagging with z=-10"
         zcorr=-10.0
      else
         read(cline,*,err=400,end=400) iright
         if (iright .le. 0 .or. iright .gt. nz) then
            zcorr = -99.9
         else
            zcorr = zcand(iright)
         end if
      end if
      write(10,'(f8.4,3x,f7.1,3x,a58)') zcorr,y,sname

c      call pgsubp(1,1)
c      call pgpage()
c      call pgsubp(nlplot,nzmax)
c      call pgpage()

      go to 200

 666  continue
      close(2)
      close(3)
      close(4)
      close(10)

      end


c  somehow read an arbitrary number of z's out of a line
      subroutine readz(cline,nz,zarr)

      character cline*80
      real zarr(nz)

      do i=1,nz
         zarr(i)=0.0
      end do
      read(cline,*,err=700) (zarr(i),i=1,nz)
 700  continue
      i=1
 720  continue
      if (zarr(i) .ne. 0.) then
         go to 720
      else
         nz=i-1
      end if

      return
      end

c  read the number of z's and then read the z's
      subroutine readz2(cline,nz,zarr)
      character cline*80
      real zarr(nz)

      read(cline,*) nz
      read(cline,*,err=800,end=800) idummy,(zarr(i),i=1,nz)
 800  continue

      return
      end
