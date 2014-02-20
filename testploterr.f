
c test various routines

      program testplot

c plot each individual spectrum?
      parameter (IFPLOTALL=1)
      parameter (IFDEBUG=1)

      include 'specarray.h'
      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname(NSPECMAX,15)
      common /specname/ specname
      character specname(NSPECMAX)*60

      real spec1(NWAVEMAX),spec2(NWAVEMAX),spec3(NWAVEMAX)
      real waverest(NLOGWMAX),specrest1(NLOGWMAX),wrestlin(NLOGWMAX)
      real sumspec(NLOGWMAX),sumsqspec(NLOGWMAX)
      real avgspec(NLOGWMAX),rmsspec(NLOGWMAX),avgsn(NLOGWMAX)
      real smoothsp(NLOGWMAX)
      integer ngoodspec(NLOGWMAX)
      real wave1(NWAVEMAX),wave2(NWAVEMAX),wave3(NWAVEMAX)
      character sname*60,xlabel*60,answ*3
      real zarray(NSPECMAX)
      integer nwarray(NSPECMAX)
      integer nwarraysky(NSPECMAX)

      include 'pcredshift.h'

      call pgbeg(0,'?',2,2)
      call pgscf(2)
      call pgsch(1.3)

      write(*,'("Full continuum fit subtraction [1,y-yes]? ",$)')
      read(*,'(a1)') answ
      if (answ(1:1) .eq. '1' .or. answ(1:1) .eq. 'y' .or.
     $     answ(1:1) .eq. 'Y') then
         ifcontsub = 1
      else
         ifcontsub = 0
      end if
      write(*,'("Do mean subtraction [1,y-yes]? ",$)')
      read(*,'(a1)') answ
      if (answ(1:1) .eq. '1' .or. answ(1:1) .eq. 'y' .or.
     $     answ(1:1) .eq. 'Y') then
         ifmeansub = 1
      else
         ifmeansub = 0
      end if
 100  write(*,
     $     '("Diameter for median smoothing, pixels [0,1=none]: ",$)')
      read(*,'(a3)') answ
      if (answ(1:3) .eq. '   ') then
         ndiamed = 0
      else
         read(answ,*,err=100) ndiamed
      end if

c      nspmax = NSPECMAX
      nspec = NSPECMAX
      nspsky = NSPECMAX

c nspec enters readspec as the dimension of the nwarray array,
c and nspec is set to the actual number of spectra by readspec
      call readspec(nspec,nwarray)

c actually we can't do this the same because of the common block...
      call readspec(nspsky,nwarraysky)

c read the redshifts
      call readz(nspec,zarray)

c      nw = nwave
c      nwlog = nwave

c wavelength zero pixel for the deshifted spectrum
c      w0rest = 1000.
      nrest = NLOGWMAX
      w0rest = 3000.
      w0rlog = log10(w0rest)
      dwrlog = 1.e-4
      do i=1,nrest
         waverest(i) = w0rlog + i*dwrlog
         wrestlin(i) = 10**waverest(i)
         ngoodspec(i) = 0
         sumspec(i) = 0.0
         sumsqspec(i) = 0.0
      end do

c for writing out the good-data-boundary for each spectrum
      open(4,file='testplot.region',status='unknown')

c loop through some spectra
      do k=1,nspec

      nw = nwarray(k)
      nwlog = NWAVEMAX
      do i=1,nw
         spec1(i) = specarr(k,i)
      end do

      sname = specname(k)
c      w0=6001.28
c      dw=1.28
c      w0log = log10(6000.)
      w0 = specw0(k)
      dw = specdw(k)
      write(*,'(a50,", w0= ",f9.3,", dw= ",f6.3)') sname,w0,dw
      w0log = log10(w0)
c      dwlog = 1.e-4
      dwlog = dwrlog
      do i=1,nw
         wave1(i) = w0+i*dw
      end do
      do i=1,nwlog
         wave2(i) = w0log+i*dwlog
      end do

c this is swiped from showspec just to print out the boundaries
c of the good data region as a test
      call findends(nw,spec1,imin,imax)
      wmin = wave1(imin)
      wmax = wave1(imax)
      write(4,'(a40,2(2x,i5),2(2x,f8.2),2(2x,f8.1))') 
     $     sname,imin,imax,wmin,wmax,spec1(imin),spec1(imax)
cccccccccccccccccccccccccccccccccc     

      if (IFPLOTALL .ne. 0) then
         call showspec(nw,wave1,spec1)
         call pglabel("wavelength (A)","counts",sname)
      end if

c median smooth the spectrum if requested
      if (ndiamed .gt. 1) then
         call medsmooth(nw,spec1,ndiamed,spec3)
      else
         do i=1,nw
            spec3(i) = spec1(i)
         end do
      end if

      if (IFDEBUG .ne. 0) then
         open(2,file='tmp1.out',status='unknown')
         do i=1,nw
            write(2,*) i,wave1(i),spec3(i)
         end do
      end if
      close(2)

      if (IFPLOTALL .ne. 0 .and. ndiamed .gt. 1) then
         call showspec(nw,wave1,spec3)
         call pglabel("wavelength (A)","median smoothed counts",sname)
      end if

      call logrebin(spec3,nw,w0,dw,nwlog,w0log,dwlog,spec2)

      call findends(nwlog,spec2,imin,imax)
      call cleanends(nwlog,spec2,imin,imax)

      if (IFDEBUG .ne. 0) then
         open(2,file='tmp2.out',status='unknown')
         do i=1,max(nw,nwlog)
            write(2,*) i,wave1(i),spec3(i),wave2(i),spec2(i)
         end do
      end if
      close(2)

      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wave2,spec2)
         call pglabel('log wavelength','counts',sname)
      end if

c      call deshiftz(spec2,nwlog,dwlog,0.43)
      if (ifcontsub .ne. 0) then
         contwrad = 100.*dwlog
         call contsubmed(spec2,nwlog,dwlog,contwrad)
      else
         call contsubconst(spec2,nwlog)
      end if

      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wave2,spec2)
         call pglabel('log wavelength','continuum-subtracted counts',
     $        sname)
      end if

      call deshiftz2(spec2,nwlog,w0log,dwlog,zarray(k),
     $     specrest1,nrest,w0rlog)
c test on ndiamed to keep 4 plots per page
      if (IFPLOTALL .ne. 0 .and. ndiamed .le. 1) then
         write(xlabel,'("log rest wavelength, z = ",f7.4)') zarray(k)
         call showspec(nrest,waverest,specrest1)
         call pglabel(xlabel,'cont-sub counts',sname)
      end if

c  accumulate the mean spectrum
      do j=1,nrest
         if (specrest1(j) .gt. bad) then
            ngoodspec(j) = ngoodspec(j) + 1
            sumspec(j) = sumspec(j) + specrest1(j)
            sumsqspec(j) = sumsqspec(j) + specrest1(j)**2
         end if
      end do

c      call pgpage()
c end the do k loop
      end do
cccccccccccccccccccc

      close(4)

c  compute mean spectrum and rms

      do j=1,nrest
         if (ngoodspec(j) .gt. 1) then
            avgspec(j) = sumspec(j) / ngoodspec(j)
            tmp = sumsqspec(j)/ngoodspec(j) - avgspec(j)**2
            rmsspec(j) = sqrt(tmp)
            if (rmsspec(j) .gt. 0.1) then
               avgsn(j) = avgspec(j) / rmsspec(j)
            else
               avgsn(j) = badset
            end if
         else if (ngoodspec(j) .eq. 1) then
            avgspec(j) = sumspec(j)
            rmsspec(j) = badset
            avgsn(j) = badset
         else
            avgspec(j) = badset
            rmsspec(j) = badset
            avgsn(j) = badset
         end if
      end do

      ismoothrad=2
      call boxcar(nrest,avgspec,ismoothrad,smoothsp)

      open(2,file='testplot.out',status='unknown')
      do j=1,nrest
         write(2,'(i5,2x,f8.5,2x,i5,4(2x,f9.3))') j,waverest(j),
     $        ngoodspec(j),avgspec(j),rmsspec(j),avgsn(j),smoothsp(j)
      end do
      close(2)

      call pgsubp(1,1)

      call showspec(nrest,waverest,avgspec)
      call pglabel("log rest wavelength","counts","average spectrum")
      call showspec(nrest,wrestlin,avgspec)
      call pglabel("rest wavelength","counts","average spectrum")

      call showspec(nrest,waverest,rmsspec)
      call pglabel("log rest wavelength","counts",
     $     "rms of average spectrum")

      call showspec(nrest,waverest,avgspec)
      call pglabel("log rest wavelength","counts","average spectrum")
      call pgsls(2)
      call pgline(nrest,waverest,rmsspec)
      call pgsls(1)

      call showspec(nrest,waverest,avgsn)
      call pglabel("log rest wavelength","sigma",
     $     "S/N - avg/rms spectrum")

      call showspec(nrest,waverest,smoothsp)
      call pglabel("log rest wavelength","counts",
     $     "average spectrum, smoothed")
      call showspec(nrest,wrestlin,smoothsp)
      call pglabel("rest wavelength","counts",
     $     "average spectrum, smoothed")

      call pgend()

      end

