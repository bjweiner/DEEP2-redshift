
c test various routines

      program testplotsky

c plot each individual spectrum?
      parameter (IFPLOTALL=1)
c write out temp files?
      parameter (IFDEBUG=0)
c do S/N rather than counts?
      parameter (IFSN=0)
c scale sky by estimating read and photon noise ?
      parameter (IFSCALENOISE=0)

      include 'specarray.h'
c      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
c      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname
      character specname(NSPECMAX)*60
      real specarr(NSPECMAX,NWAVEMAX)
      real specw0(NSPECMAX),specdw(NSPECMAX)
      character skyspecname(NSPECMAX)*60
      real skyspecarr(NSPECMAX,NWAVEMAX)
      real skyspecw0(NSPECMAX),skyspecdw(NSPECMAX)

      real spec1(NWAVEMAX),spec2(NWAVEMAX),spec3(NWAVEMAX)
      real errspec(NWAVEMAX),espec2(NWAVEMAX)
      real waverest(NLOGWMAX),specrest1(NLOGWMAX),wrestlin(NLOGWMAX)
      real sumspec(NLOGWMAX),sumsqspec(NLOGWMAX)
      real avgspec(NLOGWMAX),rmsspec(NLOGWMAX),avgsn(NLOGWMAX)
      real smoothsp(NLOGWMAX)
      integer ngoodspec(NLOGWMAX)
      real wave1(NWAVEMAX),wave2(NWAVEMAX),wave3(NWAVEMAX)
      character sname*60,xlabel*60,answ*3
      real zarray(NSPECMAX)
      integer nwarray(NSPECMAX)
      integer nskywarray(NSPECMAX)

      include 'pcredshift.h'

      call pgbeg(0,'?',2,2)
      call pgscf(2)
      call pgsch(1.3)

      write(*,'("Full continuum fit subtraction [1,y-yes]? ",$)')
      read(*,'(a1)') answ
      if (answ(1:1) .eq. '0' .or. answ(1:1) .eq. 'n' .or.
     $     answ(1:1) .eq. 'N') then
         ifcontsub = 0
      else
         ifcontsub = 1
      end if
      write(*,'("Do mean subtraction [1,y-yes]? ",$)')
      read(*,'(a1)') answ
      if (answ(1:1) .eq. '0' .or. answ(1:1) .eq. 'n' .or.
     $     answ(1:1) .eq. 'N') then
         ifmeansub = 0
      else
         ifmeansub = 1
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

c nspec enters readspec as the dimension of the nwarray array,
c and nspec is set to the actual number of spectra by readspec
c      call readspec(nspec,nwarray)
c change this so the array of data is an argument,
c the next two are the physical dimensions, then the actual
c number of spectra is returned in nspec, and the actual
c numbers of wavelengths are returned in the nwarray array
      call readspec(specarr,NSPECMAX,NWAVEMAX,nspec,nwarray,
     $     specw0,specdw,specname)
c read the sky spectra also
      call readspec(skyspecarr,NSPECMAX,NWAVEMAX,nskyspec,nskywarray,
     $     skyspecw0,skyspecdw,skyspecname)
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

c find noise spectrum, and convert to signal/noise if
c desired.  It's better to divide by noise after continuum
c subtraction.
      if (IFSCALENOISE .eq. 1) then
c convert the sky spectrum to a std.dev. spectrum
         do i=1,nw
            spec3(i) = skyspecarr(k,i)
         end do
c Placeholder assumptions: 2x30 min exposures, 5 pixel diam 
c extraction window.
         nexp = 2
         exptime = 30.*60.
c         dpix = 5.
         dpix = 1.
         call scalenoise(spec3,nw,nexp,exptime,dpix,errspec)
      else
c or if the sky spectrum is actually a rms/per pixel spectrum
c we could do nothing, or scale down by a factor of 
c sqrt(#pixels).  #pixels is most likely 7.
         do i=1,nw
c setting spec3 is vestigial
            spec3(i) = skyspecarr(k,i)
            errspec(i) = skyspecarr(k,i) / sqrt(7.0)
         end do
      end if

c  we probably would be better off dividing by the error
c  after continuum subtraction, so I've moved it to later.
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
c plot spectrum
         call showspec(nw,wave1,spec1)
         call pglabel("wavelength (A)","counts",sname)
         if (IFSN .ne. 0) then
            call pgqci(indexc)
c error spectrum
            call pgsci(3)
            call pgline(nw,wave1,errspec)
c sky spectrum
            call pgsci(4)
            call pgline(nw,wave1,spec3)
            call pgsci(indexc)
         end if
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
c  this is a total hack
      call logrebin(errspec,nw,w0,dw,nwlog,w0log,dwlog,espec2)
      do i=1,nwlog
         errspec(i) = espec2(i)
      end do

      call findends(nwlog,spec2,imin,imax)
      call cleanends(nwlog,spec2,imin,imax)

      if (IFDEBUG .ne. 0) then
         open(2,file='tmp2.out',status='unknown')
         do i=1,max(nw,nwlog)
            write(2,*) i,wave1(i),spec3(i),wave2(i),spec2(i)
         end do
      end if
      close(2)

c test on ndiamed to keep 4 plots per page
      if (IFPLOTALL .ne. 0 .and. ndiamed .le. 1) then
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

c  now, divide by the error
c  For right now we assume that the sky spectrum is wl-calibrated
c  just like the object spectrum.  We add a test for badness...
      if (IFSN .ne. 0) then
         do i=1,nw
            if (spec2(i) .gt. bad .and. errspec(i) .gt. 1.e-3) then
               spec2(i) = spec2(i) / errspec(i)
            else
               spec2(i) = badset
            end if
         end do
         if(IFPLOTALL .ne. 0) then
            call pgqci(indexc)
            call pgsci(3)
            call pgline(nwlog,wave2,spec2)
            call pgsci(indexc)
         end if
      end if


      call deshiftz2(spec2,nwlog,w0log,dwlog,zarray(k),
     $     specrest1,nrest,w0rlog)
      if (IFPLOTALL .ne. 0) then
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

