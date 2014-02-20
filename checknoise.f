
c  This is just a little test program to compare the
c  noise (RMS) calculated one way (e.g. rms along column)
c  with the RMS from photon and read noise assumptions.

      program checknoise

      include 'specarray.h'

      real sky(NWAVEMAX),skyerr(NWAVEMAX)
      real rms(NWAVEMAX),rmserr(NWAVEMAX)
      real wave(NWAVEMAX)
      integer isize(7)
      character skyfile*60,rmsfile*60
      character skyname*60,rmsname*60

      include 'pcredshift.h'

      nwmax=NWAVEMAX

c assumptions: 2 exposures, 30 min each, 7 pixel extraction window
      nexp = 2
      exptime = 30.*60.
      dpix = 7.0

      call pgbeg(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.3)

 100  continue
      write(*,'("List of sky files: ",$)')
      read(*,'(a)') skyfile
      open(2,file=skyfile,status='old',err=100)
 110  continue
      write(*,'("List of rms files: ",$)')
      read(*,'(a)') rmsfile
      open(3,file=rmsfile,status='old',err=110)
      if (skyfile(1:3) .eq. '   ' .or. rmsfile(1:3) .eq. '   ') stop

      isp=0
 200  continue
      read(2,'(a)',err=666,end=666) skyname
      read(3,'(a)',err=666,end=666) rmsname

      call imopen(skyname,1,imgs,ier)
      call imopen(rmsname,1,imgr,ier)
      if (ier .ne. 0) go to 200
      isp=isp+1
      call imgsiz(imgs,isize,idim,itype,ier)
      nws = isize(1)
      call imgsiz(imgr,isize,idim,itype,ier)
      nwr = isize(1)
      nw = max (nws,nwr)
      if (nw .gt. nwmax) 
     $   write(*,*) "Warning, maxlen ",nwmax," but files are ",nws,nwr
      call imgl1r(imgs,sky,ier)
      call imgl1r(imgr,rms,ier)
      
c  find wavelength of pixel 0, and dw/dpix
c  assume the rms and sky have same w.l. calibration
      call imgkwr(imgs,'CRPIX1',refpix,ier)
      if (ier .ne. 0) refpix = badset
      call imgkwr(imgs,'CRVAL1',refw,ier)
      if (ier .ne. 0) refw = badset
      call imgkwr(imgs,'CDELT1',dw,ier)
      if (ier .ne. 0) dw = badset
      if (refpix .gt. bad .and. refw .gt. bad .and. dw .gt. bad) 
     $     then
         w0 = refw - refpix*dw
         dw = dw
      else
         w0 = 0.0
         dw = 1.0
      end if
      do j=1,nws
         wave(j) = w0+dw*j
      end do
      call imclos(imgs,ier)
      call imclos(imgr,ier)

      call scalenoise(sky,nws,nexp,exptime,dpix,skyerr)

      do j=1,nwr
         rmserr(j) = rms(j)/sqrt(dpix)
      end do

c  auto boundaries
c      call showspec(nws,wave,skyerr)
c  find wave plot boundaries and fix y boundaries
      call findends(nws,sky,imin,imax)
      wmin = wave(imin)
      wmax = wave(imax)
      fmin = 0.
      fmax = 200.
      call pgenv(wmin,wmax,fmin,fmax,0,0)
      call pglabel("wavelength","counts"," ")

      call pgline(nws,wave,sky)
      call pgmtxt('T',2.5,0.,0.,skyname)
      call pgmtxt('T',2.5,0.6,0.,"sky")
      call pgqci(indexc)
      call pgsci(3)
      call pgline(nwr,wave,rms)
      call pgmtxt('T',1.5,0.,0.,rmsname)
      call pgmtxt('T',1.5,0.6,0.,"rms")
      call pgsci(2)
      call pgline(nws,wave,skyerr)
      call pgmtxt('T',2.5,1.0,1.,"err from sky")
      call pgsci(4)
      call pgline(nwr,wave,rmserr)
      call pgmtxt('T',1.5,1.0,1.,"err from rms")
      call pgsci(indexc)

      open(10,file='checknoise.out',status='unknown')
      do j=1,nwr
         write(10,1010) j,wave(j),sky(j),skyerr(j),rms(j),rmserr(j)
      end do
      close(10)
 1010 format(i6,5(3x,f8.2))

      go to 200

 666  continue

      call pgend()

      end
