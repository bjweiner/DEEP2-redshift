
c take code from zchoose4 to make nice plots of postage stamps
c from spectra

      program plotpstamp

      parameter (maxcol=10000,maxrow=400)
      parameter (maxpstampx=1000, maxpstampy=500)
      parameter (maxaxis=5)
c finagle plot order (works only for columns varying faster than
c rows)
      parameter (IFCOLTRICK=0)
c skip plotting in bottom row?
      parameter (ISKIPBOTTOM=1)
c put filenames above imaging pstamp?
      parameter (IPRINTNAME=0)
c change stretch for each spectrum?
      parameter (ISPECSTRETCH=1)
c change stretch for each image?
      parameter (IIMSTRETCH=0)
c fatter, shorter spectra?
      parameter (IFATSPEC=1)

      real specdata(maxcol*maxrow)
      real spec2d(maxcol,maxrow)
      real data1d(maxcol)
      real dwave1(maxcol,maxrow)
      real wave1(maxcol)
      integer laxes1(maxaxis),laxes2(maxaxis)
      logical flag

      real pstamp(maxpstampx,maxpstampy)
c      real pstasmooth(maxpstampx,maxpstampy)
      real trans(6)

      character fname1*160,fname2*160,xlabel*30,toplabel*50

c      include 'pcredshift.h'

c change this if need to invert background/out of bounds color
      badflag = -10000.

c coord transf for image plotting
      trans(1) = 0.
      trans(2) = 1.
      trans(3) = 0.
      trans(4) = 0.
      trans(5) = 0.
      trans(6) = 1.

c image scales, arcsec/pixel
      deimospix = 0.119
      cfhtpix = 0.206

      nplot=0
      ierr=0

      idwin=pgopen('?')
      write(*,'("pgplot nx by ny: ",$)')
      read(*,*) nxplot,nyplot
      call pgsubp(nxplot,nyplot)
      nxplotold = nxplot
      nxplot = abs(nxplot)
      call pgpage()
      call pgscf(2)
      call pgsch(2.7)
      call pgscir(16,254)

      ixpl = 1
      iypl = 1
      call pgpanl(ixpl,iypl)

      write(*,'("y radius (Deimos pixels): ",$)')
      read (*,*) yrad
      write(*,'("km/sec radius: " ,$)')
      read(*,*) velrad
      write(*,'("Deimos flux min, max: ",$)')
      read(*,*) fmin,fmax
      write(*,'("CFHT flux min, max: ",$)')
      read(*,*) fmin_im,fmax_im

 100  continue
      write(*,'("2d spectrum file [quit]: ",$)')
      read(*,'(a)') fname1
      if (fname1(1:3) .eq. '   ') go to 666
      
      write(*,'("ypos (pixels): ",$)')
      read(*,*) ypos
      iypos = int(ypos)
      iyrad = int(yrad)
c      write(*,'("wavelength radius (A): ",$)')
c      read(*,*) wrad
      write(*,'("rest wavelength and z: ",$)')
      read(*,*) wrest,z
      wcen = wrest*(1.+z)

      if (ISPECSTRETCH .eq. 1) then
         write(*,'("spectrum fmin, fmax: ",$)')
         read(*,*) fmin,fmax
      end if

      call ftgiou(im1,ierr)
      call doerr(ierr)
      call ftopen(im1,fname1,0,block,ierr)
      if (ierr .ne. 0) then
         write(*,'(a,a)') "Couldnt open image ",fname1
         go to 100
      end if
      call doerr(ierr)

c retrieve 2-d data from hdu 2
      ihdu = 2
      
      call readfield(im1,ihdu,'FLUX',itype1,isize1,
     $     naxis1,laxes1,specdata)
      call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,spec2d)
      nx1=laxes1(1)
      ny1=laxes1(2)
      call readfield(im1,ihdu,'DLAMBDA',itype1,isize1,
     $     naxis1,laxes1,specdata)
      call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,dwave1)
      call readfield(im1,ihdu,'LAMBDA0',itype1,iw0size1,
     $     nw0axis1,lwaxis1,wave1)

      call ftclos(im1,ierr)
      call ftfiou(im1,ierr)
      call doerr(ierr)
      do i=1,nx1
         wave1(i) = wave1(i) + dwave1(i,iypos)
      end do
      dw1 = (wave1(nx1) - wave1(1)) / (nx1-1.)
      w01 = wave1(1) - dw1
c      write(*,'(a80,2x,f7.2,2x,f5.3)') fname1,w01,dw1

      call findwave(nx1,wave1,wcen,ixcen1)
c      xrad = wrad /dw1
      xrad = velrad / 3.0e5 * wcen / dw1
      ixrad = int(xrad)
      write(*,*) ixcen1,ixrad,iypos,iyrad

      ioffset = ixcen1-ixrad-1
      joffset = iypos-iyrad-1
      do j=iypos-iyrad,iypos+iyrad
         jpsta = j-joffset
         do i=ixcen1-ixrad,ixcen1+ixrad
            ipsta = i-ioffset
            if (i .ge. 1 .and. i .le. nx1 .and. 
     $           j .ge. 1 .and. j .le. ny1) then
               pstamp(ipsta,jpsta) = spec2d(i,j)
            else
               pstamp(ipsta,jpsta) = badflag
            end if
         end do
      end do

      nplotx = 2*ixrad+1
      nploty = 2*iyrad+1

c      call pgqvp(0,x1,x2,y1,y2)
c      write(*,*) x1,x2,y1,y2
c      call pgsvp(0.1,0.5,0.18,0.82)
      if (IFATSPEC .eq. 0) then
         call pgsvp(0.1,0.62,0.3,0.7)
      else
         call pgsvp(0.1,0.53,0.2,0.8)
      end if
      nplot=nplot+1

c      call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,1,0)
c square pixels
c      call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,1,-2)
c non-square pixels
c      call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,0,-2)
c      call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,0,-2)
      call pgswin(0.5,nplotx+0.5,0.5,nploty+0.5)
c arcsec per y pixel
      arcppix = deimospix
      yarcsecrad = yrad*arcppix
c      trans(1) = wcen-wrad
c      trans(2) = dw1
c      trans(3) = 0.0
c      trans(4) = -yarcsecrad
c      trans(5) = 0.0
c      trans(6) = arcppix
c      call pgenv(wcen-wrad,wcen+wrad,-yarcsecrad,yarcsecrad,0,0)
c 0 for linear
      call pgsitf(2)
      call pggray(pstamp,maxpstampx,maxpstampy,
c     $     ixpmin,ixpmax,iypmin,iypmax,
     $     1,nplotx,1,nploty,
     $     fmin,fmax,trans)
c      call pgswin(wcen-wrad,wcen+wrad,-yarcsecrad,yarcsecrad)
      call pgswin(-velrad,velrad,-yarcsecrad,yarcsecrad)
      write(xlabel,'("km/sec")')
      write(toplabel,'(i4," at z = ",f6.4)') nint(wrest),z
c      call pglabel(xlabel, 'arcsec', toplabel)
c      call pglabel(' ','arcsec',' ')
c if plot is on left hand column, write y label
      if (ixpl .eq. 1)
     $     call pgmtxt('L',2.0,0.5,0.5,"arcsec")
c if plot is on bottom, write x label
      if ((ISKIPBOTTOM .eq. 0 .and. iypl .eq. nyplot) .or.
     $     (ISKIPBOTTOM .eq. 1 .and. iypl .eq. nyplot-1)) 
     $     call pgmtxt('B',2.6,0.5,0.5,'km/sec')
      call pgmtxt('T',1.1,0.0,0.0,toplabel)

c ticks out
c      call pgbox('BCINST',0.0,0,'BCINST',0.0,0)
      call pgbox('BCNST',0.0,0,'BCNST',0.0,0)

c temp
c      go to 100

c now read in an image
 200  continue
      write(*,'("image fits file: ",$)')
      read(*,'(a)') fname1
            
      call ftgiou(im1,ierr)
      call doerr(ierr)
      call ftopen(im1,fname1,0,block,ierr)
      if (ierr .ne. 0) then
         write(*,'(a,a)') "Couldnt open image ",fname1
         go to 200
      end if
      call doerr(ierr)

      call ftgknj(im1,'NAXIS',1,7,laxes1,idim,ierr)
      call doerr(ierr)
      nx = laxes1(1)
      ny = laxes1(2)
c      write(*,'("xmin, xmax, ymin, ymax: ",$)')
c      read (*,*) xmin,xmax,ymin,ymax

c or compute it automatically, centering the image and making
c its size match the deimos slit
      yrad2 = yrad * deimospix/cfhtpix
      yarcsecrad2 = yrad2 * cfhtpix
      xrad2 = yrad2
      xmin = (nx+1)/2. -xrad2
      xmax = (nx+1)/2. +xrad2
      ymin = (ny+1)/2. -yrad2
      ymax = (ny+1)/2. +yrad2

      ixmin = int(xmin)
      ixmax = int(xmax)+1
      iymin = int(ymin)
      iymax = int(ymax)+1
      nplotx = ixmax-ixmin+1
      nploty = iymax-iymin+1

      if (IIMSTRETCH .eq. 1) then
         write(*,'("image fmin, fmax: ",$)')
         read(*,*) fmin_im,fmax_im
      end if

c      write(*,*) nx,ny,nplotx,nploty,yrad2,yarcsecrad2,
c     $     ixmin,ixmax,iymin,iymax

      do j=iymin,iymax
         jpsta = j-iymin+1
         if (j .lt. 1 .or. j. gt. ny) then
            do i=ixmin,ixmax
               ipsta = i-ixmin+1
               pstamp(ipsta,jpsta) = badflag
            end do
         else
            ipixel = (j-1)*nx + 1
            call ftgpve(im1,0,ipixel,nx,0.0,data1d,flag,ierr)
            call doerr(ierr)
            do i=ixmin,ixmax
               ipsta = i-ixmin+1
               if (i .lt. 1 .or. i .gt. nx) then
                  pstamp(ipsta,jpsta) = badflag
               else
                  pstamp(ipsta,jpsta) = data1d(i)
               end if
            end do
         end if
      end do

c      call pgsvp(0.6,0.9,0.18,0.82)
      if (IFATSPEC .eq. 0) then
         call pgsvp(0.68,0.92,0.3,0.7)
      else
         call pgsvp(0.58,0.92,0.2,0.8)
      end if
c      call pgswin(0.5,nplotx+0.5,0.5,nploty+0.5)
c want square pixels, but don't want it to skip to next panel ...
c      call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,1,-2)
      call pgwnad(0.5,nplotx+0.5,0.5,nploty+0.5)
c      call pgswin(0.5,nplotx+0.5,0.5,nploty+0.5)

c pgsitf(0) = linear, 1= log, 2= square root
      call pgsitf(2)
      call pggray(pstamp,maxpstampx,maxpstampy,
c     $     ixpmin,ixpmax,iypmin,iypmax,
     $     1,nplotx,1,nploty,
     $     fmin_im,fmax_im,trans)
      call pgswin(-yarcsecrad2,+yarcsecrad2,-yarcsecrad2,+yarcsecrad2)
c ticks out
c      call pgbox('BCINST',0.0,0,'BCIMST',0.0,0)
      call pgbox('BCNST',0.0,0,'BCMST',0.0,0)
c      call pglabel('arcsec' , '', '')

c      irow = mod(nplot,nyplot)
c      if (irow .eq. 0) irow = nyplot
c      icol = int(nplot / nyplot) +1
c      write(*,*) "icol, irow = ",icol,irow
      write(*,*) "icol, irow = ",ixpl,iypl

c if plot is on bottom, write x label
      if ((ISKIPBOTTOM .eq. 0 .and. iypl .eq. nyplot) .or.
     $     (ISKIPBOTTOM .eq. 1 .and. iypl .eq. nyplot-1)) 
     $     call pgmtxt('B',2.6,0.5,0.5,'arcsec')
c temporary for id'ing object
      if (IPRINTNAME .eq. 1) call pgmtxt('T',1.2,0.0,0.0,fname1)

      if (IFCOLTRICK .eq. 0) then
c don't finagle plot order
         call pgpage()
         if (nxplotold .ge. 0) then
c fill rows first
            if (ixpl .lt. nxplot) then
               ixpl=ixpl+1
            else
               ixpl=1
               if(iypl .lt. nyplot) then
                  iypl = iypl+1
               else
                  iypl=1
               end if
            end if
         else
c fill columns first
            if((ISKIPBOTTOM .eq. 0 .and. iypl .lt. nyplot) .or.
     $           (ISKIPBOTTOM .eq. 1 .and. iypl .lt. nyplot-1)) then
               iypl = iypl+1
            else
               iypl=1
               if(ixpl .lt. nxplot) then
                  ixpl=ixpl+1
               else
                  ixpl=1
               end if
            end if
         end if
      else
c finagling, column order
c if we are at the next-to-bottom row
         if (ISKIPBOTTOM .eq. 1 .and. iypl .eq. nyplot-1) then
            iypl=1
            if (ixpl .lt. nxplot) then
               ixpl=ixpl+1
            else
c advance to next page by going to last panel and calling pgpage
               call pgpanl(nxplot,nyplot)
               call pgpage()
               ixpl=1
            end if
            call pgpanl(ixpl,iypl)
         else
            call pgpage()
            if(iypl .lt. nyplot) then
               iypl = iypl+1
            else
               iypl = 1
               if (ixpl .lt. nxplot) then
                  ixpl = ixpl+1
               else
                  ixpl = 1
               end if
            end if
         end if
      end if

      go to 100

 666  continue

      call pgclos

      end



cccccccccccccccccccccccccccccc
c  find the index closest to a given wave in the wavelength array,
c  assumed increasing
      subroutine findwave(nw,warr,w,iw)

      real warr(nw)

c  define a buffer in Angstroms, so we can plot even if slightly off end
c  of spectrum
      wbuff = 5.0
      if (warr(1)-wbuff .gt. w) then
c iw=-1 means w not in array
         iw=-1
         return
      end if
      j=2
 100  continue
      if(warr(j) .lt. w) then
         j=j+1
c test for reaching the end
         if (j .gt. nw) then
            if (warr(nw)+wbuff .gt. w) then
c just barely off
               iw=nw
            else
               iw=-1
            end if
            return
         end if
c if not at end, try again
         go to 100
      else
c warr(j) > w, we've gone past w.  Pick the closest one
         if (abs(warr(j)-w) .lt. abs(warr(j-1)-w)) then
            iw = j
         else
            iw=j-1
         end if
         return
      end if

      end

cccccccccccccccccccccccccccccc
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

