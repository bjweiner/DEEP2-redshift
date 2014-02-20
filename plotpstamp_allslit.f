
c Plot a 2-d spectrum postage stamp by reading the Allslits file

      program plotpstamp_allslit

      parameter (maxpstampx=512, maxpstampy=300)
      parameter (maxcol=10000,maxrow=5000)
      parameter (IFTOPLABEL=1)

c      real spec2d(maxcol,maxrow)
      real specline(maxcol),speccol(maxrow)
      real pstamp(maxpstampx,maxpstampy)
      real trans(6)
      integer isize(7)

      character iname*140,toplabel*40
      integer function pgopen

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
      call pgsch(1.6)
      call pgscir(16,254)

c      ixpl = 1
c      iypl = 1
c      call pgpanl(ixpl,iypl)

 100  continue
      write(*,'("Allslits image: ",$)')
      read(*,'(a)') iname
      if (iname(1:3) .eq. '   ') go to 666

      call imopen(iname,1,img,ier)
      if (ier .ne. 0) then
         write(*,*) "Couldnt open image "
         go to 100
      end if
      call imgsiz(img,isize,idim,itype,ier)
      nx = isize(1)
      ny = isize(2)
      call imgkwr(img,'CRPIX1',refpix,ier)
      call imgkwr(img,'CRVAL1',refval,ier)
      call imgkwr(img,'CD1_1',dwdpix,ier)
      write (*,*) "Image ",nx," by ",ny, " pixels"
      write(*,*) "Wavecal ",refpix,refval,dwdpix      
      write(*,'("slit number: ",$)')
      read (*,*) nslit
      rnslit = real(nslit)

c find slit and object by looking at pixel values in rightmost column
c the loop is a little inefficient.
      jslitmin = 20000
      jslitmax = -100
      call imgs2r(img,speccol,nx,nx,1,ny,ier)
      do j=1,ny
c         call imgl2r(img,specline,j,ier)
c         if (abs(specline(nx)-rnslit) .lt. 0.01) then
         if (abs(speccol(j)-rnslit) .lt. 0.01) then
            jslitmin = min(jslitmin,j)
            jslitmax = max(jslitmax,j)
         end if
      end do
      if (jslitmin .lt. ny+1 .and. jslitmax .gt. 0) then
         write(*,*) "Slit extent in y: ",jslitmin,jslitmax
      else
         write(*,*) "Couldnt find slit? ",nslit,jslitmin,jslitmax
         go to 100
      end if
      jslitobj = nint((jslitmin+jslitmax)/2.0)
      do j=jslitmin,jslitmax
c         call imgl2r(img,specline,j,ier)
c         if (specline(nx) .gt. 500.) jslitobj = j
         if (speccol(j) .gt. 500.) jslitobj = j
      end do
      write(*,*) "Object at y = ",jslitobj

      write(*,'("obs. wavelength, plot wave rad., y pix rad: ",$)')
      read (*,*) wobs,wrad,ypixrad
      write(*,'("grayscale flux min, max : ",$)')
      read (*,*) fmin, fmax

      yarcsecrad = ypixrad * deimospix
      xcen = (wobs-refval)/dwdpix + refpix
      xmin = xcen - wrad/dwdpix
      xmax = xcen + wrad/dwdpix
c region of Allslits to plot
      imin = int(xmin)
      imax = int(xmax)+1
      jmin = int(jslitobj-ypixrad)
      jmax = int(jslitobj+ypixrad)+1

c bounds check on region to plot
      imin = max(imin,1)
      imax = min(imax,nx)
      jmin = max(jmin,jslitmin)
      jmax = min(jmax,jslitmax)

      write(*,*) "pixel region ",imin,imax,jmin,jmax
      write(*,*) 1,imax-imin+1,jmin-jslitmin+1,jmax-jslitmin+1

c section out the part of the image to plot
c and put in pstamp

      do j=jslitmin,jslitmax
         call imgl2r(img,specline,j,ier)
         do i=imin,imax
            pstamp(i-imin+1,j-jslitmin+1) = specline(i)
         end do
      end do
      call imclos(img,ier)

c Note that it doesn't force the pixels to be square,
c would be nice to try doing that using dwdpix and deimospix.
c for 2x2, vp is 0.074 0.926 0.10 0.90
      call pgqvp(0,xvp1,xvp2,yxp1,yxp2)
      write(*,*) "viewport ",xvp1,xvp2,yxp1,yxp2
      call pgswin(wobs-wrad,wobs+wrad,-yarcsecrad,yarcsecrad)

c transformation from pstamp array coords to wavelength and
c arcsec coords.  trans(1) is wavelength of pstamp pixel i=0
c trans(4) is arcsec pos of pstamp j=0 relative to object.
      trans(1) = (imin-1.0-refpix)*dwdpix + refval
      trans(2) = dwdpix
      trans(3) = 0.0
      trans(4) = (jslitmin-1.0-jslitobj) * deimospix
      trans(5) = 0.0
      trans(6) = deimospix

      call pgsitf(0)
      call pggray(pstamp,maxpstampx,maxpstampy,
c     $     ixpmin,ixpmax,iypmin,iypmax,
     $     1,imax-imin+1,jmin-jslitmin+1,jmax-jslitmin+1,
     $     fmin,fmax,trans)
c      call pgswin(-velrad,velrad,-yarcsecrad,yarcsecrad)

      call pgbox('BCINST',0.0,0,'BCINST',0.0,0)
c      call pgmtxt('L',2.0,0.5,0.5,"arcsec")
c      call pgmtxt('B',2.6,0.5,0.5,"wavelength, A")
      call pglabel("wavelength","arcsec","")
      if (IFTOPLABEL .eq. 1) then
         write(*,'("id or top label: ",$)')
         read(*,'(a)') toplabel
         call pgmtxt('T',1.0,0.0,0.0,toplabel)
      end if
      call pgpage()
      go to 100

 666  continue
      call pgclos()

      end
