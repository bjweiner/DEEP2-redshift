
c Plot a series of DEEP2 spectra in restframe from a list of
c spec1d files and a matching list of id and z

      program plotdeepspec

c      parameter (NXSPEC=8200)
      parameter (IFMEDSMOOTH=0)
      include 'specarray.h'

      real wave(NWAVEMAX),flux(NWAVEMAX),ivar(NWAVEMAX)
      real ferr(NWAVEMAX),smooflux(NWAVEMAX)
      character sname*120,lname*120,iname*120,cline*120
      character idstring*8,clabel*20
      real xpl(2),ypl(2),wpl(5)

 100  continue
      write(*,'("File with list of z and ID: ",$)')
      read(*,'(a)') iname
      open(2,file=iname,status='old',err=100)

 110  continue
      write(*,'("File with list of spec1d files: ",$)')
      read(*,'(a)') lname
      open(3,file=lname,status='old',err=110)

      write(*,'("wave plot min, max: ",$)')
      read(*,*) wmin,wmax
      write(*,'("flux plot min, max [0 0 for auto-scale]: ",$)')
      read(*,*) fmin, fmax
      if (abs(fmin) .lt. 1.e-4 .and. abs(fmax) .lt. 1.e-4) then
         ifyauto = 1
         fmin = 0.0
         fmax = 100.
      else
         ifyauto = 0
      end if

      write(*,'("median/boxcar smooth, pixel diameter: ",$)')
      read (*,*) msmoodiam

      write(*,'("pgplot nx by ny: ",$)')
      read(*,*) nxplot,nyplot

      call pgbeg(0,'?',nxplot,nyplot)
c      call pgscf(2)
      call pgslw(2)
      call pgsch(1.8)

      iplot = 1

 300  continue
      read (3,'(a)',err=666,end=666) sname
      read (2,*,err=666,end=666) z,idnum
      write(idstring,'(i8)') idnum
      write(*,*) "read file number ",iplot," id ",idnum
cc Try to find where the ID ends.  This will fail with leading spaces
c      indx = index(cline,' ')
cc Assume the ID is in the first 8 columns
c      read(cline(1:8),'(a)') idstring
c      read(cline(9:),*) z
c Or require z to come first
c      read (2,'(a)',err=666,end=666) cline 
c      read(cline,*,err=666,end=666) z,idstring
c      read (2,*) z,idstring


c Read the spectrum, using read1dspec (no rebinning)
      ifopt=0
      call read1dspec(sname,ifopt,npts,wave,flux,ferr,iferr)
      if (iferr .ne. 0) then
         write(*,*) "Error reading ",sname
      end if

c convert to restwave
      do i=1,npts
         wave(i) = wave(i) / (1.+z)
      end do

c median smooth the spectrum.
      if (msmoodiam .gt. 1) then
         if (IFMEDSMOOTH .eq. 1) then
            call medsmooth(npts,flux,msmoodiam,smooflux)
         else
            call boxsmooth(npts,flux,msmoodiam,smooflux)
         end if
c rescale error?
      else
         do i=1,npts
            smooflux(i) = flux(i)
         end do
      end if

c auto scale if needed
      if (ifyauto .eq. 1) then
         fmin = 1.0e4
         fmax = -1.0e4
         do i=1,npts
            if (wave(i) .gt. wmin .and. wave(i) .lt. wmax .and.
     $       smooflux(i) .gt. -1.0e4 .and. smooflux(i) .lt. 1.0e5) then
               fmin = min(fmin,smooflux(i))
               fmax = max(fmax,smooflux(i))
            end if
         end do
      end if
      fmarg = 0.05*(fmax-fmin)

      ymin = fmin-fmarg
      ymax = fmax+fmarg
      call pgenv(wmin,wmax,ymin,ymax,0,0)
      call pglabel("wavelength","flux,"," ")
c      call pgline(npts,wave,smooflux)
      call pgbin(npts,wave,smooflux,.true.)
c Plot the error?
      write(clabel,'(a8," z = ",f7.4)') idstring(1:8),z
      xplot = wmin + 0.05*(wmax-wmin)
      yplot = ymax - 3.0*fmarg
      call pgptxt(xplot,yplot,0.0,0.0,clabel)

c Draw dotted lines for zero-velocity Mg II in vacuum
      nlines = 2
      wpl(1) = 2796.352
      wpl(2) = 2803.531
      wpl(3) = 2852.964
      ypl(1) = ymin
      ypl(2) = ymax
      do k=1,nlines
         xpl(1) = wpl(k)
         xpl(2) = wpl(k)
         call pgsls(4)
         call pgline(2,xpl,ypl)
         call pgsls(1)
      end do

      iplot=iplot+1
      go to 300

 666  continue
      close(2)
      close(3)
      call pgend()

      end


cccccccccccccccccccc
c boxcar smooth the spectrum.  Try to ignore bad pixels

      subroutine boxsmooth(n,x,ndiam,xsmoo)

      real x(n),xsmoo(n),s
      integer n,np

      bad = -6000.
      irad1 = int(ndiam/2)
      irad2 = ndiam-1-irad1

      do i=1,n
         j1 = max(1,i-irad1)
         j2 = min(n,i+irad2)
         s = 0.0
         np = 0
         do j=j1,j2
            if (x(j) .gt. bad) then
               s = s + x(j)
               np = np+1
            end if
         end do
         xsmoo(i) = s/np
      end do
      return
      end

