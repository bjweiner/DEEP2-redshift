
c  read a list of images, y coord and possible redshifts
c  and plot postage stamps around locations of features
c  allow the user to pick the correct z by typing a number
c  (or hitting the cursor?)

c  This version is designed to work with the DEEP2 data
c  output by the Berkeley pipeline.  It wants two lists,
c  of B and R images

c Still todo: get fwhm, ypos for each spectrum and indicate
c   extraction window w/ horizontal ticks?
c Check why objects aren't exactly centered left/right?
c   I now believe this is because z is ~ 0.0002 too high,
c   i.e. problem is in fitz2

c 11.11.02 finally fixed the problem with narrow slits appearing in the
c          bottom half of the window: iypmin,iypmax were neither 
c          calculated nor used right.
c 10.22.02 added autosaving output file every ~5 objects
c          changed "s" behavior to set z but not output
c          changed "w" to default to whatever lineset was there before
c            fixed bug which had "w" overwriting lineset 5 rather than 6
c          substituted 3933.7 for 3889.1 in "a" lineset
c          added display of median smoothing pixel # in 1-d window
c 10.21.02 added reading in old output file
c          labeling old zchoice, qual, comments at top of 2-d window
c 10.20.02 error checking on reads so program doesn't crash
c            from typeahead in text window
c          defaults on input file names
c          indicate all linesets in 1-d window, font size bigger, ticks longer
c          label zchoice and z in 1-d window
c          "z" in 1-d window cycles through choices rather than prompting
c          added expand, replot, smoothing in 1-d window
c          anykey in 1-d indicates restwave as well as observed wavelength
c 10.18.02 added prompts for quality and comments after each object
c          made pstamp length 2*20 A
c          use g to skip-ahead rather than returning to previous undone obj
c          "i" in 1-d window does not return you to 2-d immediately
c          flux limits do not reset on new object
c          prompts for confirmation if you try to overwrite output file

      program zchoose3

c max size of 2-d image, max # of z's to plot
      parameter (maxcol=10000,maxrow=400,maxz=8)
c max size of a postage stamp in pixels
      parameter (maxpstampx=500, maxpstampy=300)
c max # of lines to plot simultaneously, # of defined sets of lines
      parameter (maxlines=15,maxlinesets=7)
c max # of 2-d images to read per object
      parameter (maximages=2)
      parameter (maxaxis=5)
c max number of spectra to examine
      parameter (maxspectra=1000)
c do statistics to determine the grayscale?
      parameter(IFSTAT=0)
c prompt for a quality level and comments?
      parameter(IFQUAL=1)

      real specdata(maxcol*maxrow)
c      real dwavedata(maxcol*maxrow)
c      real spec2d(maxcol,maxrow,maximages)
      real spec2d1(maxcol,maxrow),spec2d2(maxcol,maxrow)
      real dwave1(maxcol,maxrow),dwave2(maxcol,maxrow)
      real wave1(maxcol),wave2(maxcol)
      integer laxes1(maxaxis),laxes2(maxaxis)
      integer lwaxis1(maxaxis),lwaxis2(maxaxis)
c      real spec1d(maxcol)
c      real tmpdat(maxcol)

      real pstamp(maxpstampx,maxpstampy)
      real pstasmooth(maxpstampx,maxpstampy)

c      real wave1d(maxcol),flux1d(maxcol),ferr1d(maxcol)
      real wlopt(maxlines)
      real wlem(maxlines),wlabs(maxlines),wluv(maxlines)
      real wlstar(maxlines)
      real wlplot(maxlines,maxlinesets)
      real zcand(maxz)
c      real xpl(2),ypl(2)
c      character wlabel*10
      real trans(6)
      real cheight
c      integer isize(7)
c      integer ims(maximages)

      integer function pgopen
      integer function pgcurs

      character fname1*160,fname2*160,fnameoned*160
c      character sname1*160,sname2*160,snameoned*160
      character fname*160,cline*160
c      character answ*1,sline*160
      character toplabel*180,xlabel*10

      character sname1(maxspectra)*160
      character sname2(maxspectra)*160
      character snameoned(maxspectra)*160
      character zcandline(maxspectra)*200
      real ypos(maxspectra)
      real zout(maxspectra)
      real zqual(maxspectra)
      integer nzchoice(maxspectra)
      character zcomment(maxspectra)*80

c wavelengths of common features, taking into account
c that I only plot the first 7
      data wlopt/3728.0, 3933.7, 3968.5, 4861.3, 4959.9, 5006.8,
     $     6562.8, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0/
      data wlem/3728.0, 3868.74, 4959.9, 5006.8, 6562.8, 6583.4, 
     $     6716.5, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0/
c      data wlabs/1550.0, 2800.0, 3798.6, 3835.6, 3889.1, 3933.7,
c     $     3968.5, 4101.75, 4304.4, 4340.5, 5184.0, 5889.0,
c     $     0.0, 0.0, 0.0/
      data wlabs/3728.0, 3933.7, 4101.75, 4304.4, 4340.5, 4861.3,
     $     5183.6, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0/
c      data wluv/1215.67, 1304.37, 1334.53, 1526.71, 1548.2, 1550.77,
c     $     1608.45, 1670.79, 2344.21, 2374.46, 2382.77, 2586.65,
c     $     2600.17, 2796.35, 2803.53/
      data wluv/1215.67, 1550.0, 2344.21, 2378.0, 2586.65, 2600.17,
     $     2800.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0/
      data wlstar/4304.4, 5183.6, 5890.0, 6562.8, 8498.1, 8542.1,
     $     8662.2, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0/

      ierr=0

c max number of candidate z's to read
      nzmax = 5
      if (nzmax .gt. maxz) nzmax=maxz
      
c number of lines to actually plot
      nlplot = 7
c number of sets which are defined
      nlinesets=6
      if (nlinesets .gt. maxlinesets .or. nlplot .gt. maxlines)
     $     write(*,*) "Error in linesets, too many lines/sets!",
     $     nlinesets,maxlinesets,nlplot,maxlines
      do i=1,nlplot
         wlplot(i,1) = wlopt(i)
         wlplot(i,2) = wlabs(i)
         wlplot(i,3) = wluv(i)
         wlplot(i,4) = wlem(i)
         wlplot(i,5) = wlstar(i)
         wlplot(i,6) = wlopt(i)
      end do

c wavelength in A, and y radius in pixels to plot
c      wrad = 10.0
      wrad = 20.0
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
c font 1 appears more readable
      call pgscf(1)
      call pgsch(1.5)
      call pgask(.false.)

c      call pgqcir(ilo,ihi)
c      write(*,*) ilo,ihi
      call pgscir(16,254)

 100  continue

 110  write(*,'("List of 1-D spectra [noserendip.spec1d]: ",$)') 
      read(*,'(a)') fnameoned
      if (fnameoned(1:3) .eq. '   ') fnameoned='noserendip.spec1d'
      open(8,file=fnameoned,status='old',err=110)
 112  write(*,'("First list of images [noserendip.slit2d.B]: ",$)') 
      read(*,'(a)') fname1
      if (fname1(1:3) .eq. '   ') fname1='noserendip.slit2d.B'
      open(2,file=fname1,status='old',err=112)
 115  write(*,'("Second list of images [noserendip.slit2d.R]: ",$)') 
      read(*,'(a)') fname2
      if (fname2(1:3) .eq. '   ') fname2='noserendip.slit2d.R'
      open(7,file=fname2,status='old',err=115)
 120  write(*,'("File with y-coords [noserendip.ypos]: ",$)') 
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') fname='noserendip.ypos'
      open(3,file=fname,status='old',err=120)
 130  write
     $   (*,'(a,a,$)') "File with candidate redshifts ",
     $     "[fitz.noserendip.med9.out]: "
      read(*,'(a)') fname
c      if (fname(1:3) .eq. '   ') fname='sdssfitz.med9.out'
      if (fname(1:3) .eq. '   ') fname='fitz.noserendip.med9.out'
      open(4,file=fname,status='old',err=130)

 132  write(*,'("old output file to compare [none]: ",$)')
      read(*,'(a)') fname
      if(fname(1:3) .eq. '   ') then
         ifoldout=0
      else
         open(12,file=fname,status='old',err=132)
         ifoldout=1
      end if

 138  continue
      write(*,'("new output file [zchoose.out]: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') fname='zchoose.out'

c      open(10,file='zchoose.out',status='unknown',err=140)
      open(10,file=fname,status='new',err=140)
      go to 145
 140  continue
      write(*,'("couldnt open output file, overwrite? [y/n] ",$)')
      read(*,'(a)') cline
      if (cline(1:1) .eq. 'y') then
         open(10,file=fname,status='old',err=140)
      else
         go to 138
      end if
c      read (*,'(a)') fname
c      open(10,file=fname,status='new',err=140)
 145  continue

c read necessary information from files

      nobj=1
 160  continue
      if (nobj .gt. maxspectra) then
         write(*,*) "Max # spectra is ",maxspectra
         go to 180
      end if
      read(2,'(a)',end=180) sname1(nobj)
      read(7,'(a)',end=180) sname2(nobj)
      read(8,'(a)',end=180) snameoned(nobj)
      read(4,'(a)',end=180) zcandline(nobj)
      read(3,*,err=180,end=180) ypos(nobj)
      nobj=nobj+1
      go to 160
      
 180  continue
      nobj=nobj-1
      write(*,'("Read names and data for ",i4," objects")') nobj
      close(2)
      close(3)
      close(4)
      close(7)
      close(8)

c if we have an old output file to check against
      if (ifoldout .eq. 1) then
         do i=1,nobj
            if (IFQUAL .eq. 0) then
               read(12,'(f8.4,3x,i2,3x,i4,3x,a120)',err=190,end=190) 
     $              zout(i),nzchoice(i),idum,cline
            else
               read(12,'(f8.4,3x,f5.1,3x,i2,3x,i4,3x,a120,3x,a80)',
     $              err=190,end=190) 
     $          zout(i),zqual(i),nzchoice(i),idum,cline,zcomment(i)
            end if
         end do
         go to 192
 190     continue
         write(*,'(a,a)') "Warning - Format problems ",
     $     "or not enough lines in old output file!"
 192     continue
         close(12)
      else
c initialize output arrays
         do i=1,nobj
            zout(i) = -99.9
            nzchoice(i) = 0
            zqual(i) = 0.0
            zcomment(i) = ' '
         end do
      end if
         

c test
c      do i=1,nobj
c         write(*,'(a)') sname1(i)
c         write(*,'(a)') sname2(i)
c         write(*,'(a)') snameoned(i)
c      end do

      fmin=-12.
      fmax=25.
      write(*,*) "Default flux range: ",fmin,fmax

      call writehelp()

c main loop through objects using kobj as index
c kobjdone records how many objects an entry was made for
      kobjdone = 0
      kobj=1

 200  continue

      y = ypos(kobj)
      iy = nint(y)
c      read(4,'(a)',end=666) cline
      write(*,*) 
      write(*,*) "Object # ",kobj,"  ypos = ",iy

      nz=nzmax
c  readz resets nz
c      call readz(cline,nz,zcand)
c  just use the version that requires nz in the file for now
c      call readz2(cline,nz,zcand)
      call readz2(zcandline(kobj),nz,zcand)

c get data from images.   Use NULL to mean the image didn't exist.
c behavior may not be predictable if one of the images
c doesn't exist
      nimag= 0
      if (sname1(kobj)(1:4) .ne. 'NULL') then
         call ftgiou(im1,ierr)
         call doerr(ierr)
         call ftopen(im1,sname1(kobj),0,block,ierr)
         if (ierr .ne. 0) then
            write(*,'(a,a)') "Couldnt open image ",sname1(kobj)
         else
            nimag=nimag+1
         end if
         call doerr(ierr)
      end if
      if(sname2(kobj)(1:4) .ne. 'NULL') then
         call ftgiou(im2,ierr)
         call doerr(ierr)
         call ftopen(im2,sname2(kobj),0,block,ierr)
         if (ierr .ne. 0) then
            write(*,'(a,a)') "Couldnt open image ",sname2(kobj)
         else
            nimag=nimag+1
         end if
         call doerr(ierr)
      end if

      if (nimag .eq. 0) then
c bail
         write(*,*) "Found no images!"
         go to 400
      else if (nimag .eq. 1) then
         write(*,*) "Only got 1 image!"
c should do something         
      end if

c retrieve 2-d data from hdu 2
      ihdu = 2
      
      call readfield(im1,ihdu,'FLUX',itype1,isize1,
     $     naxis1,laxes1,specdata)
      call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,spec2d1)
      call readfield(im2,ihdu,'FLUX',itype2,isize2,
     $     naxis2,laxes2,specdata)
      call reformat(isize2,specdata,naxis2,laxes2,maxcol,maxrow,spec2d2)
c      write(*,*) "flux arrays ",isize1,isize2,naxis1,laxes1(1),laxes1(2)
c this is wasteful of memory since we don't really need the
c entire 2-d array of dwave1 and dwave2
      call readfield(im1,ihdu,'DLAMBDA',itype1,isize1,
     $     naxis1,laxes1,specdata)
      call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,dwave1)
      call readfield(im2,ihdu,'DLAMBDA',itype2,isize2,
     $     naxis2,laxes2,dwave2)
      call reformat(isize1,specdata,naxis1,laxes1,maxcol,maxrow,dwave2)
c      write(*,*) "wave arrays ",isize1,isize2,naxis1,laxes1(1),laxes1(2)
      call readfield(im1,ihdu,'LAMBDA0',itype1,iw0size1,
     $     nw0axis1,lwaxis1,wave1)
      call readfield(im2,ihdu,'LAMBDA0',itype2,iw0size2,
     $     nw0axis2,lwaxis2,wave2)

      nx1=laxes1(1)
      ny1=laxes1(2)
      nx2=laxes2(1)
      ny2=laxes2(2)
      nx = min(nx1,nx2)
      ny = min(ny1,ny2)
      if (iy .le. 0 .or. iy .gt. ny1 .or. iy .gt. ny2) then
         write(*,*) "Something's wrong: y position ",iy,
     $        " is out of bounds"
      end if
c      write(*,*) "nx, ny ",nx,ny

c compute wavelength array for the central row of extraction.
c after this we only need to use wave1 & wave2
c      do k=1,nimag
      do i=1,nx1
         wave1(i) = wave1(i) + dwave1(i,iy)
      end do
      do i=1,nx2
         wave2(i) = wave2(i) + dwave2(i,iy)
      end do
c linear approximations. dw1 and dw2 are used to determine how
c many pixels wide the windows will be
      dw1 = (wave1(nx1) - wave1(1)) / (nx1-1.)
      w01 = wave1(1) - dw1
      dw2 = (wave2(nx1) - wave2(1)) / (nx2-1.)
      w02 = wave2(1) - dw2

      write(*,'(a80,2x,f7.2,2x,f5.3)') sname1(kobj),w01,dw1
      write(*,'(a80,2x,f7.2,2x,f5.3)') sname2(kobj),w02,dw2

c dw1 and dw2 should be pretty similar so can just use dw1
c for the plot radius,  unless a wavecal is off...
      dw=dw1
      if(dw1 .lt. 0.2 .and. dw2 .gt. 0.2) then
         write(*,*) "Note, wavecal for blue 2d slit is off"
         dw=dw2
      else if(dw1 .gt. 0.2 .and. dw2 .lt. 0.2) then
         write(*,*) "Note, wavecal for red 2d slit is off"
         dw=dw1
      else if(dw1 .lt. 0.2 .and. dw2 .lt. 0.2) then
         write(*,*) "Note, wavecals for both 2d slits are likely off"
         dw=0.33
      end if

      xrad = wrad /dw
      ixrad = int(xrad)
c      xrad2 = wrad /dw2
c      ixrad2 = int(xrad2)

      iymin = max(1,iy - iyrad)
      iymax = min(ny,iy + iyrad)
      nyp = iymax-iymin+1
c      iycen = iyrad+1
      iycen = iy
c      do j=1,nyp
c         jimage=j+iymin-1
cc         call imgl2r(img,spec2d(j,1),jimage,ier)
c         call imgl2r(img,tmpdat,jimage,ier)
c         do i=1,nx
c            spec2d(i,j) = tmpdat(i)
c         end do
c      end do
c      call imclos(img,ier)
      call ftclos(im1,ierr)
      call ftclos(im2,ierr)
      call ftfiou(im1,ierr)
      call ftfiou(im2,ierr)
      call doerr(ierr)
      

c      call findends(nx,spec2d(iycen,1),imin,imax)
c      wmin = w0+imin*dw
c      wmax = w0+imax*dw

c don't reset grayscale range, do it above?
c      fmin=-7.
c      fmax=25.
c box smoothing for 2-d images: #of pixels in wave and spatial (y) dir
      nxsmoo = 1
      nysmoo = 1
c lineset to plot, 1=em+H&K, 2=optical abs, 3=UV, etc
      ilineset = 1

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
            wcen = wlplot(ipan,ilineset) * (1. + zcand(jpan))
c            ixcen = int((wcen-w0)/dw)
            call findwave(nx1,wave1,wcen,ixcen1)
            call findwave(nx2,wave2,wcen,ixcen2)
c  choose which to plot
            if (ixcen1 .gt. 0) then
               ixcen=ixcen1
               iplot=1
            else if (ixcen2 .gt. 0) then
               ixcen=ixcen2
               iplot=2
            else
               ixcen = ixcen1
               iplot=0
            end if
            rixmin=real(ixcen-ixrad)
            rixmax=real(ixcen+ixrad)
            riymin=real(iycen-iyrad)
            riymax=real(iycen+iyrad)
c  pixel coords in the pstamps array (before trimming to edges
c  of actually existing data)
            nplotx = 2*ixrad+1
            nploty = 2*iyrad+1
c square pixels with JUST (5th arg) = 1
c            call pgenv(rixmin,rixmax,riymin,riymax,1,-1)
c now we're giving the corrdinates in the pstamp array
            call pgenv(0.5,nplotx+0.5,0.5,nploty+0.5,1,-1)
c make the char height very big for labeling because the large # of 
c subpanels made it small
            call pgqch(cheight)
            call pgsch(4.0)
c label number and z over first panel
            if (ipan .eq. 1) then
               write(toplabel,'(i2,", z= ",f7.4)') jpan,zcand(jpan)
               call pgmtxt('T',1.4,0.5,0.5,toplabel)
            end if            
c label wl under each panel
            write(xlabel,'(f7.1)') wlplot(ipan,ilineset)
            call pgmtxt('B',1.4,0.5,0.5,xlabel)
c label y and image  over panel (2,1) and (last,1)
c            if (ipan .eq. 2 .and .jpan .eq. 1) then
c               write(toplabel,'("y = ",i5)') iy
c               call pgmtxt('T',1.4,0.5,0.5,toplabel)
c            end if               
            if (ipan .eq. nlplot .and .jpan .eq. 1) then
c              write(toplabel,'("object ",i4,1x,a)') kobj,sname1(kobj)
              write(toplabel,'("object ",i4,",",1x,a)') 
     $              kobj,snameoned(kobj)
c              call pgmtxt('T',1.4,1.0,1.0,sname1(kobj))
              call pgmtxt('T',1.4,1.0,1.0,toplabel)
            end if
c label your choice, z, zqual, and comments over panel 1,1
            if (ipan .eq. 1 .and. jpan .eq. 1) then
               write(toplabel,2010) nzchoice(kobj),zout(kobj),
     $              zqual(kobj),zcomment(kobj)
               call pgmtxt('T',3.0,0.0,0.0,toplabel)
 2010          format("prev. choice ",i2,", z= ",f8.4,", qual= ",
     $              f5.1," comments: ",a80)
            end if
c reset character height
            call pgsch(cheight)
c bounds checking
cc  if totally off the plot, skip pggray
c            if (rixmin .lt. real(nx) .and. rixmax .gt. 1.0 .and. 
c     $           riymin .lt. real(nyp) .and. riymax .gt. 1.0) then
            if (iplot .ne. 0) then

c            call pgpanl(ipan,jpan)
c            call pgenv(rixmin,rixmax,riymin,riymax,1,-1)

c            if (ifirst .eq. 1) then
c               write(toplabel,'(i2,"  z= ",f7.4)') jpan,zcand(jpan)
c               call pgmtxt('T',1.,0.,0.,toplabel)
c               ifirst=0
c            end if
c trim plot region to boundaries of the spec2d image
c this is the region to copy into pstamp
c  11/29/02 - try not trimming here, but in the loop that 
c    sets up pstamp
               ixmin=ixcen-ixrad
               ixmax=ixcen+ixrad
               iymin=iycen-iyrad
               iymax=iycen+iyrad
c               iymax=min(riymax,real(nyp))
c this is the region to be plotted's
c  coordinates in the pstamps image
               ixpmin=1
               ixpmax=ixpmin+ ixmax-ixmin
               iypmin=iymin+(iyrad+1-iycen)
               iypmax=iymax+(iyrad+1-iycen)
c               nplotxtrim=ixpmax-ixpmin+1
c               nplotytrim=iypmax-iypmin+1
c could determine flux limits for each plot
               if (IFSTAT .ne. 0) then
                  npix=0
                  sum=0.
                  sumsq=0.0
                  smin=1.e10
                  smax=-1.e10
                  do iy=iymin,iymax
                     do ix=ixmin,ixmax
                        npix=npix+1
                        if (iplot .eq. 1) then
                           d=spec2d1(ix,iy)
                        else if (iplot .eq. 2) then
                           d=spec2d2(ix,iy)
                        end if
                        sum = sum+d
                        sumsq=sumsq+d*d
                        smin=min(smin,d)
                        smax=max(smax,d)
                     end do
                  end do
                  write(*,*) smin,smax,sum/n,sqrt(sumsq-sum*sum/n)
               end if
c copy the region into pstamp.  we want to look from ixmin to ixmax,
c etc, but the offset in coords between spec2d1 or spec2d2 and pstamp
c is based on ixcen-ixrad,iycen-iyrad translating to 1,1
               ioffset=ixcen-ixrad-1
               joffset=iycen-iyrad-1
               do jj=iymin,iymax
                  jjnew = jj-joffset
                  do ii=ixmin,ixmax
                     if (iplot .eq. 1) then
                        if(ii .ge. 1 .and. ii .le. nx1 .and.
     $                       jj .ge. 1 .and. jj .le. ny1) then
                           pstamp(ii-ioffset,jjnew) = spec2d1(ii,jj)
                        else
                           pstamp(ii-ioffset,jjnew) = badset
                        end if
                     else if (iplot .eq. 2) then
                        if(ii .ge. 1 .and. ii .le. nx2 .and.
     $                       jj .ge. 1 .and. jj .le. ny2) then
                           pstamp(ii-ioffset,jjnew) = spec2d2(ii,jj)
                        else
                           pstamp(ii-ioffset,jjnew) = badset
                        end if
                     end if
                  end do
               end do
c smooth the pstamp
               if(nxsmoo .gt. 1 .or. nysmoo .gt. 1) then
                  call smootharr(pstamp,maxpstampx,maxpstampy,
     $                 ixpmin,ixpmax,iypmin,iypmax,
     $                 nxsmoo,nysmoo,pstasmooth)
c When using /xs it seems necessary to reverse fmax and fmin
c to get the inverted display (black on white)
c that is, the order of arguments is foreground,background
                  call pggray(pstasmooth,maxpstampx,maxpstampy,
     $                 ixpmin,ixpmax,iypmin,iypmax,
     $                 fmin,fmax,trans)
               else 
                  call pggray(pstamp,maxpstampx,maxpstampy,
     $                 ixpmin,ixpmax,iypmin,iypmax,
     $                 fmin,fmax,trans)
               end if

c              if (iplot .eq. 1) then
c              call pggray(spec2d1,maxcol,maxrow,ixmin,ixmax,iymin,iymax,
c     $              fmax,fmin,trans)
c              else if (iplot .eq. 2) then
c              call pggray(spec2d2,maxcol,maxrow,ixmin,ixmax,iymin,iymax,
c     $              fmax,fmin,trans)
c              end if
            end if
         end do
      end do
      
c  this is a total kludge, forced because pgsubp and pgpanl just 
c  do not work right. (?)
      do j=nz+1,nzmax
         do i=1,nlplot
            call pgpage()
         end do
      end do

 400  continue
      write(*,'("# of right z (0-none,f-change scale): ",$)') 
      read(*,'(a)') cline
      if (cline(1:1) .eq. 'h') then
         call writehelp()
         go to 400
      else if (cline(1:1) .eq. 'b') then
c go back one object
         kobj=kobj-1
         go to 200
      else if (cline(1:1) .eq. 'q') then
c quit out of loop
         go to 666
      else if (cline(1:1) .eq. 'g') then
c go to object N
 405     write(*,'("object # to go to: ",$)')
         read(*,*,err=405) kobj
         if (kobj .gt. nobj .or. kobj .lt. 1) then
            write(*,*) "there are only ",nobj," objects!"
            go to 405
         end if
         go to 200
      else if (cline(1:1) .eq. 'f') then
c change flux limits
 407     write(*,'("range is ",f5,1x,f5," new fmin, fmax: ",$)') 
     $        fmin,fmax
         read(*,*,err=407) fmin, fmax
         go to 300
      else if (cline(1:1) .eq. 'm') then
c median smooth
 410     write(*,'("smoothing box in x and y: ",$)')
         read(*,*,err=410) nxsmoo, nysmoo
         nxsmoo=max(1,nxsmoo)
         nysmoo=max(1,nysmoo)
         go to 300
      else if (cline(1:1) .eq. 'z') then
c z to use for last row
 412     write(*,'("enter z: ",$)')
         read(*,*,err=412) zcand(nz)
         go to 300
      else if (cline(1:1) .eq. 't') then
c tag and output
         write(*,'(a)') "Tagging with z=-10"
         zout(kobj)=-10.0
      else if (cline(1:1) .eq. 's') then
         write(*,'(a)') "star, trying z=0"
c tag, output, fall off end of loop? decided not to
c         zout(kobj)=0.0
         zcand(nz) = 0.0
         go to 300         
      else if (cline(1:1) .eq. 'o') then
         write(*,'(a)') "plotting optical features"
         ilineset=1
         go to 300
      else if (cline(1:1) .eq. 'a') then
         write(*,'(a)') "plotting absorption features"
         ilineset=2
         go to 300
      else if (cline(1:1) .eq. 'e') then
         write(*,'(a)') "plotting emission features"
         ilineset=4
         go to 300
      else if (cline(1:1) .eq. 'u') then
         write(*,'(a)') "plotting UV features"
         ilineset=3
         go to 300
      else if (cline(1:2) .eq. 'l') then
         write(*,'(a)') "plotting stellar features"
         ilineset=5
         go to 300
      else if (cline(1:1) .eq. 'w') then
c get a new lineset from input
 420     continue
         write(*,'("[Default ",7(1x,f7.1))') 
     $        (wlplot(i,nlinesets),i=1,nlplot)
         write(*,'("Input ",i2," rest wl: ",$)') nlplot
         read(*,'(a)') cline
c if line starts with blanks, don't change
         if (cline(1:5) .ne. '     ')
     $        read(cline,*,err=420,end=420) 
     $        (wlplot(i,nlinesets),i=1,nlplot)
         ilineset=nlinesets
         go to 300
      else if (cline(1:1) .eq. 'p') then
c plot 1-d spectrum
         call pgsubp(1,1)
         xpl1=-1.
         xpl2=-1.
         ypl1=0.
         ypl2=0.
         nmed1d=9
         iplz=1
         write(*,'(a,a)') "hit key in plot window, ",
     $        "i-id line, z-cycle through z choices, m-med smooth,"
         write(*,'(a)') 
     $        "e-expand, r-replot, n to continue, other-location"

c 1-d key loop
c replot and wait for new key press
 430     continue
         call show1dspec(snameoned(kobj),xpl1,xpl2,ypl1,ypl2,nmed1d,1)
c  this plots current lineset
c         call plotlineset(zcand(1),nlplot,ilineset,wlplot,
c     $        maxlines,maxlinesets)
c  this plots all defined linesets
         call plotalllines(zcand(iplz),nlplot,nlinesets,wlplot,
     $        maxlines,maxlinesets)
         write(toplabel,'("# ",i2,", z= ",f7.4)') iplz,zcand(iplz)
         call pgqch(cheight)
         call pgsch(1.1)
         call pgmtxt('T',4.0,0.0,0.0,toplabel)
         write(toplabel,'("obj# ",i4,", median smooth = ",i3)') 
     $        kobj,nmed1d
         call pgmtxt('T',4.0,1.0,1.0,toplabel)
         call pgsch(cheight)

         istat = pgcurs(xc,yc,cline)
         if (cline(1:1) .eq. 'z') then
c 432        write(*,'("# of z to plot features for: ",$)')
c            read(*,*,err=432) iplz
            iplz=iplz+1
            if (iplz .gt. nz) iplz=1
            go to 430
         else if (cline(1:1) .eq. 'i') then
 435        write(*,'("rest wave to identify as: ",$)')
            read(*,*,err=435) rwave
            zcand(nz) = xc/rwave -1.0
            iplz = nz
c if no goto 430, we fall off end and go back to 2d plot
c I decided to replot on 1-d plot
            go to 430
         else if (cline(1:1) .eq. 'm') then
 440        write(*,'("diameter for median smooth (pixels): ",$)')
            read(*,*,err=440) nmed1d
            nmed1d=max(1,nmed1d)
            go to 430
         else if (cline(1:1) .eq. 'e') then
            xpl1=xc
            ypl1=yc
            write(*,'(a)') "hit key at 2nd corner in plot window: "
            istat = pgcurs(xc,yc,cline)
c handle correctly no matter which order corners are indicated 
            xpl2=max(xpl1,xc)
            ypl2=max(ypl1,yc)
            xpl1=min(xpl1,xc)
            ypl1=min(ypl1,yc)
            go to 430
         else if (cline(1:1) .eq. 'r') then
            xpl1=-1.
            xpl2=-1.
            go to 430
         else if (cline(1:1) .eq. 's') then
c set the last z to zero, show it, and replot
            iplz=nz
            zcand(nz)=0.0
            go to 430
         else if (cline(1:1) .ne. 'n') then
            write(*,'("obs ",f8.2,"  flux ",f8.2,"  rest ",f8.2)')
     $           xc,yc,xc/(1.+zcand(iplz))
            go to 430
         end if
c end of the 1-d loop
         call pgsubp(nlplot,nzmax)
         go to 300
      else
c read number of correct z
         read(cline,*,err=400,end=400) iright
         if (iright .le. 0 .or. iright .gt. nz) then
            zout(kobj) = -99.9
            nzchoice(kobj) = 0
         else
            zout(kobj) = zcand(iright)
            nzchoice(kobj) = iright
         end if
      end if
c we get here if a number was typed (z was chosen)
c or t for tag, s for star
      write(*,'(a,i4,a,f8.4)') "Object # ",kobj,"  z = ",zout(kobj)

c ask for a quality estimate and comments
      if(IFQUAL .ne. 0) then
 500     continue
         write(*,'("z quality estimate [default: ",f5.1,"]: ",$)')
     $        zqual(kobj)
         read(*,'(a)',err=500) cline
         if (cline(1:5) .eq. '     ') then
c            zqual(kobj) = 0.0
c don't change if return is hit     
            continue
         else
            read(cline,*,err=500) zqual(kobj)
         end if
 520     continue
         write(*,'("Comments [",a60,"]: ")') zcomment(kobj)
         read(*,'(a)',err=520) cline
         if (cline(1:5) .eq. '     ') then
c don't change            
            continue
         else
            read(cline,'(a)',err=500) zcomment(kobj)
         end if
      end if


c      call pgsubp(1,1)
c      call pgpage()
c      call pgsubp(nlplot,nzmax)
c      call pgpage()

c we've done another object
      kobjdone=kobjdone+1
c auto-save every N objects, e.g. 5
      if (mod(kobjdone,5) .eq. 0) then
         write(*,*) "Saving ... "
         call writeoutput(10,IFQUAL,nobj,zout,zqual,nzchoice,
     $     snameoned,zcomment)
      end if
c increment the object# counter
      kobj=kobj+1

c if we haven't reached end, do next object
      if (kobj .le. nobj) go to 200

c loop has ended
 666  continue

      call pgclos

c write output arrays to file
      call writeoutput(10,IFQUAL,nobj,zout,zqual,nzchoice,
     $     snameoned,zcomment)
c      do i=1,nobj
c         if (IFQUAL .eq. 0) then
c            write(10,'(f8.4,3x,i2,3x,i4,3x,a120)') 
c     $           zout(i),nzchoice(i),i,snameoned(i)
c         else
c            write(10,'(f8.4,3x,f5.1,3x,i2,3x,i4,3x,a120,3x,a80)')
c     $           zout(i),zqual(i),nzchoice(i),i,
c     $           snameoned(i),zcomment(i)
c         end if
c      end do

      close(10)

      end

cccccccccccccccccccccccccccccc
c  somehow read an arbitrary number of z's out of a line
c  this may not be working
      subroutine readz(cline,nz,zarr)

      character cline*160
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

cccccccccccccccccccccccccccccc
c  read the number of z's and then read the z's
      subroutine readz2(cline,nz,zarr)
      character cline*160
      real zarr(nz)

      read(cline,*) nz
      read(cline,*,err=800,end=800) idummy,(zarr(i),i=1,nz)
 800  continue

      return
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

cccccccccccccccccccccccccccccc
      subroutine writehelp()

      write(*,'(a,a)') 
     $     "h-help, #-num of correct z, 0-no good z, ",
     $     "b-go back one obj., g-go to object #N, q-quit"
      write(*,'(a,a)') 
     $     "o-opt. lines, a-absorp, e-emiss, u-UV, l-stellar, ",
     $     "w-input a set of lines"
      write(*,'(a,a)') 
     $     "m-smooth, f-flux limits, z-enter a z, t-tag w/z=-10, ",
     $     "s-star,set z=0, p-plot 1d spec"

      return
      end

cccccccccccccccccccccccccccccc
c get and plot the one-d spectrum
c Pass in plot boundaries (if x1 and x2 are negative, choose
c automatically) and diameter to median smooth

      subroutine show1dspec(sname,x1,x2,y1,y2,ndiamed,ifploterr)
      
      parameter(maxcol=10000)
      character sname*(*)
      real wave1d(maxcol),flux1d(maxcol),ferr1d(maxcol)

c      ndiamed = 7

c  note I am using wave1d as a temporary array to store 
c  unsmoothed flux
      call get1ddeepspec(sname,n1d,wave1d,ferr1d,w01d,dw1d,ierr)
      if (ierr .ne. 0) then
         write(*,*) "couldnt open 1d spectrum ",snameoned
         ierr=0
         return
      end if

c trim out hot pix by setting to some low number,
c hopefully below the bad-flag
      do i=1,n1d
         if(wave1d(i) .gt. 5000.) wave1d(i) = -6000.
      end do
      call medsmooth(n1d,wave1d,ndiamed,flux1d)
      do i=1,n1d
         wave1d(i) = w01d + i*dw1d
      end do

c if we want to automatically choose plot boundaries      
      if (x1 .lt. 0. .and. x2 .lt. 0.) then
         if (ifploterr .ne. 0) then
            call showspecerr(n1d,wave1d,flux1d,ferr1d)
         else
            call showspec(n1d,wave1d,flux1d)
         end if
      else
         call pgenv(x1,x2,y1,y2,0,0)
         if (ifploterr .ne. 0) then
            call pgqci(ioldci)
            call pgsci(3)
            call pgline(n1d,wave1d,ferr1d)
            call pgsci(ioldci)
         end if
         call pgline(n1d,wave1d,flux1d)
      end if

      return
      end

cccccccccccccccccccccccccccccc
c  box median smooth an array
c  nxphys and nyphys are the physical dimensions
c  ix1, ix2, iy1,iy2 are the region to smooth
c  nxs, nys are the smoothing box diameters
      subroutine smootharr(arr1,nxphys,nyphys,ix1,ix2,iy1,iy2,
     $     nxs,nys,arr2)

c max number of pixels in box
      parameter(maxbox=1000)
c do average or median
      parameter(ifmed=0)

      real arr1(nxphys,nyphys),arr2(nxphys,nyphys)
      real dat(maxbox)

c bad flag
      bad = -600.

      if (ifmed .ne. 0 .and. nxs*nys .gt. maxbox) then
c override request for median
         write(*,*) "Warning, med smoothing box requested is too big"
         ifdomed = 0
      else
         ifdomed = ifmed
      end if

      ixrad1 = int(nxs/2)
      ixrad2 = nxs -ixrad1-1
      iyrad1 = int(nys/2)
      iyrad2 = nys -iyrad1-1
      
c smooth
      do j=iy1,iy2
         do i=ix1,ix2
c boundary trimming
            kx1 = max(ix1,i-ixrad1)
            kx2 = min(ix2,i+ixrad2)
            ky1 = max(iy1,j-iyrad1)
            ky2 = min(iy2,j+iyrad2)
            npts=0
            sum = 0.0
            do jj=ky1,ky2
               do ii=kx1,kx2
                  if (arr1(ii,jj) .gt. bad) then
                     npts=npts+1
                     if(ifdomed .eq. 0) then
                        sum = sum+arr1(ii,jj)
                     else
                        dat(npts) = arr1(ii,jj)
                     end if
                  end if
               end do
            end do
            if (npts .gt. 0) then
               if (ifdomed .eq. 0) then
                  arr2(i,j) = sum/npts
               else
                  itmp = int(npts/2.)
                  arr2(i,j) = select(itmp,npts,dat)
               end if
            else
               arr2(i,j) = bad
            end if
         end do
      end do

      return
      end
      
cccccccccccccccccccc
c  indicate line positions on 1-d plot using one lineset
c  wlplot is the wavelength array, nlines and nlinesets are its
c   physical dimensions, nlplot is # of lines to plot,ilineset
c   is which lineset to use, z is z
      subroutine plotlineset(z,nlplot,ilineset,wlplot,nlines,nlinesets)

      real z,wlplot(nlines,nlinesets)
      integer nlplot,ilineset
      real xpl(2),ypl(2)
      character wlabel*12

      call pgqch(oldch)
      call pgsch(1.01)
      call pgqci(ioldci)
      call pgsci(2)
      call pgqwin(xw1,xw2,yw1,yw2)
      do ii=1,nlplot
         xpl(1)=wlplot(ii,ilineset)*(1.+z)
         xpl(2)=xpl(1)
         ypl(1)=yw2-0.3*(yw2-yw1)
c         ypl(2)=yw2-0.08*(yw2-yw1)
         ypl(2)=yw2
         call pgline(2,xpl,ypl) 
         write(wlabel,'(f7.1)') wlplot(ii,ilineset)
         call pgptxt(xpl(1),ypl(2),90.,0.0,wlabel)
      end do
      call pgsch(oldch)
      call pgsci(ioldci)
      return
      end

cccccccccccccccccccc
c  indicate line positions on 1-d plot for up to N linesets
c  wlplot is the wavelength array, nlines and nlinesets are its
c   physical dimensions, nlplot is # of lines to plot,ilineset
c   is max lineset to use, z is z
      subroutine plotalllines(z,nlplot,ilineset,wlplot,nlines,nlinesets)

      real z,wlplot(nlines,nlinesets)
      integer nlplot,ilineset

      do j=1,ilineset
         call plotlineset(z,nlplot,j,wlplot,nlines,nlinesets)
      end do
      return
      end

cccccccccccccccccccccccccccccc
c  write the output arrays to a file
      subroutine writeoutput(iunit,ifqual,n,zout,zqual,nzchoice,
     $     snameoned,zcomment)
      
      real zout(n),zqual(n)
      integer nzchoice(n)
      character snameoned(n)*160
      character zcomment(n)*80
      integer iunit, ifqual

      rewind(unit=iunit)
      do i=1,n
         if (ifqual .eq. 0) then
            write(iunit,3001)
     $           zout(i),nzchoice(i),i,snameoned(i)
         else
            write(iunit,3003)
     $           zout(i),zqual(i),nzchoice(i),i,
     $           snameoned(i),zcomment(i)
         end if
      end do
      
 3001 format(f8.4,3x,i2,3x,i4,3x,a120)
 3003 format(f8.4,3x,f5.1,3x,i2,3x,i4,3x,a120,3x,a80)
      return
      end
