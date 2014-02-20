
c  read a bunch of spectra
c  this is a first crack and doesn't unify the B and R spectra

c      subroutine readspec(nsp,nwave)
c change this so the array of data is an argument,
c the next two are the physical dimensions, then the actual
c number of spectra is returned in nspec, and the actual
c numbers of wavelengths are returned in the nwarray array
      subroutine readspec(specarr,nsphys,nwphys,nsp,nwarray,
     $     specw0,specdw,specname)

      real specarr(nsphys,nwphys)
      real specw0(nsphys),specdw(nsphys)
      character specname(nsphys)*60
      integer nwarray(nsp)
      character fname*80,sname*80

      include 'specarray.h'

c      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
c      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname
c      character specname(NSPECMAX)*60

      real tempdat(NWAVEMAX)
      integer isize(7)

      include 'pcredshift.h'

      nsmax = NSPECMAX
      nwmax = NWAVEMAX

 100  continue
      write(*,'("File with list of spectra [quit]: ",$)')
      read (*,'(a)') fname
      if (fname(1:3) .eq. '   ') go to 666
      open(2,file=fname,status='old',err=100)

      nsp = 0
 200  continue
      read(2,'(a)', err=250,end=250) sname
      
      call imopen(sname,1,img,ier)
      call imgsiz(img,isize,idim,itype,ier)
      nx = isize(1)
      ny = isize(2)
      if (ier .eq. 0) then
         nsp=nsp+1
c         call imgl2r(img,specarr(nsp,1),nsp,ier)
         call imgl1r(img,tempdat,ier)
         do i=1,nx
            specarr(nsp,i) = tempdat(i)
         end do
c         call imgl1r(img,specarr(nsp,1),ier)
         if (ier .ne. 0) then
            nsp=nsp-1
         else
            nwarray(nsp) = nx
            specname(nsp) = sname
c  find wavelength of pixel 0, and dw/dpix
            call imgkwr(img,'CRPIX1',refpix,ier)
c            if (ier .ne. 0) refpix = badset
            if (ier .ne. 0) then
               ier=0
               refpix = 0.
            end if
            call imgkwr(img,'CRVAL1',refw,ier)
            if (ier .ne. 0) refw = badset
c            call imacck(img,'CDELT1',ier)
c            if (ier .eq. 0) then
               call imgkwr(img,'CDELT1',dw,ier)
c            end if
            if (ier .ne. 0 .or. abs(dw) .lt. 1.e-3) then
               ier=0
               call imgkwr(img,'CD1_1',dw,ier)
               if (ier .ne. 0) dw = badset
            end if
            if (refpix .gt. bad .and. refw .gt. bad .and. dw .gt. bad) 
     $           then
               specw0(nsp) = refw - refpix*dw
               specdw(nsp) = dw
            else
c  we should really do something here
               specw0(nsp) = badset
               specdw(nsp) = badset
               write(*,
     $     '("Couldnt get w0 and dw for ",a," enter w0,dw: ",$)') sname
               read(*,*) specw0(nsp),specdw(nsp)
            end if
c            write(*,'(a40,2x,f8.2,2x,f6.3)') 
c     $           sname,specw0(nsp),specdw(nsp)
         end if
      end if
      call imclos(img,ier)
      go to 200

 250  continue
c      write(*,'("Read ",i6," spectra, length ",i6)') nsp,nwave
      write(*,'("Read ",i6," spectra")') nsp

      return

 666  continue
      stop
      return
      end
