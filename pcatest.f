
c test of attempting to do a PCA

      program pcatest

c plot each individual spectrum?
c      parameter (IFPLOTALL=0)

      include 'specarray.h'
c      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
c      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname
      character specname(NSPECMAX)*60
      real specarr(NSPECMAX,NWAVEMAX)
      real specw0(NSPECMAX),specdw(NSPECMAX)

      real speclogarr(NSPECMAX,NLOGWMAX)

      real specmatrix(NLOGWMAX,NSPECMAX)
c      real eigenmatrix(NSPECMAX,NSPECMAX)
      real weights(NSPECMAX,NSPECMAX)
      real eigenvals(NSPECMAX)
      integer indeigen(NSPECMAX)

      real spec1(NWAVEMAX),spec2(NWAVEMAX),spec3(NWAVEMAX)
      real waverest(NLOGWMAX),specrest1(NLOGWMAX),wrestlin(NLOGWMAX)
      real sumspec(NLOGWMAX),sumsqspec(NLOGWMAX)
      real avgspec(NLOGWMAX),rmsspec(NLOGWMAX),avgsn(NLOGWMAX)
      real smoothsp(NLOGWMAX)
      integer ngoodspec(NLOGWMAX)
      real wave1(NWAVEMAX),wave2(NWAVEMAX),wave3(NWAVEMAX)
      character sname*60,xlabel*60,toplabel*60,answ*3
      real zarray(NSPECMAX)
      integer nwarray(NSPECMAX)

      include 'pcredshift.h'

c      call pgbeg(0,'?',2,2)
      call pgbeg(0,'?',1,1)
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

      write(*,
     $     '("Restwave spectrum region to use in PCA, min, max: ",$)')
      read(*,*) pcawmin,pcawmax
      pwminlog = log10(pcawmin)
      pwmaxlog = log10(pcawmax)

      nspmax = NSPECMAX
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
c         ngoodspec(i) = 0
c         sumspec(i) = 0.0
c         sumsqspec(i) = 0.0
      end do

c for writing out the good-data-boundary for each spectrum
c      open(4,file='testplot.region',status='unknown')

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
c      do i=1,nw
c         wave1(i) = w0+i*dw
c      end do
c      do i=1,nwlog
c         wave2(i) = w0log+i*dwlog
c      end do

cc this is swiped from showspec just to print out the boundaries
cc of the good data region as a test
c      wmin = wave1(1)
c      imin = 1
c      i = 1
c 1000 continue
c      if (spec1(i) .gt. bad) then
c         wmin = wave1(i)
c         imin = i
c      else if (i .eq. nw) then
c         go to 1010
c      else
c         i = i+1
c         go to 1000
c      end if
c 1010 continue
c      wmax = wave1(nw)
c      imax = nw
c      i=nw
c 1050 continue
c      if (spec1(i) .gt. bad) then
c         wmax = wave1(i)
c         imax = i
c      else if (i .eq. 1) then
c         go to 1060
c      else
c         i = i-1
c         go to 1050
c      end if
c 1060 continue
c      write(4,'(a40,2(2x,i6),2(2x,f8.2))') sname,imin,imax,wmin,wmax
ccccccccccccccccccccccccccccccccccc     

c median smooth the spectrum if requested
      if (ndiamed .gt. 1) then
         call medsmooth(nw,spec1,ndiamed,spec3)
      else
         do i=1,nw
            spec3(i) = spec1(i)
         end do
      end if

      call logrebin(spec3,nw,w0,dw,nwlog,w0log,dwlog,spec2)

      call findends(nwlog,spec2,imin,imax)
      call cleanends(nwlog,spec2,imin,imax)

      if (ifcontsub .ne. 0) then
         contwrad = 100.*dwlog
         call contsubmed(spec2,nwlog,dwlog,contwrad)
      else
         call contsubconst(spec2,nwlog)
      end if
      call deshiftz2(spec2,nwlog,w0log,dwlog,zarray(k),
     $     specrest1,nrest,w0rlog)
      
      do i=1,nrest
         speclogarr(k,i) = specrest1(i)
      end do

c end the do k loop
      end do

c  blank out the regions outside the region of interest
c  specified by pcawmin,pcawmax
      call blankout(speclogarr,nspec,nspmax,nrest,waverest,
     $     pwminlog,pwmaxlog)

c  find mean spectrum and plot
      call compmean(speclogarr,nspec,nspmax,nrest,avgspec,rmsspec)
      call showspec(nrest,waverest,avgspec)
      call pglabel("log rest wavelength","counts","average spectrum")
      call showspec(nrest,wrestlin,avgspec)
      call pglabel("rest wavelength","counts","average spectrum")

      call showspec(nrest,waverest,rmsspec)
      call pglabel("log rest wavelength","counts",
     $     "rms of average spectrum")
      do i=1,nrest
         if (rmsspec(i) .gt. 0.1) then
            avgsn(i) = avgspec(i)/rmsspec(i)
         else
            avgsn(i) = badset
         end if
      end do
      call showspec(nrest,waverest,avgsn)
      call pglabel("log rest wavelength","sigma",
     $     "S/N - avg/rms spectrum")
     
      ismoothrad=2
      call boxcar(nrest,avgspec,ismoothrad,smoothsp)
      call showspec(nrest,waverest,smoothsp)
      call pglabel("log rest wavelength","counts",
     $     "average spectrum, smoothed")
      call showspec(nrest,wrestlin,smoothsp)
      call pglabel("rest wavelength","counts",
     $     "average spectrum, smoothed")

      open(2,file='pcatest.out1',status='unknown')
      do j=1,nrest
         write(2,'(i5,2x,f8.5,4(2x,f9.3))') j,waverest(j),
     $        avgspec(j),rmsspec(j),avgsn(j),smoothsp(j)
c         write(2,'(i5,2x,f8.5,2x,i5,4(2x,f9.3))') j,waverest(j),
c     $        ngoodspec(j),avgspec(j),rmsspec(j),avgsn(j),smoothsp(j)
      end do
      close(2)

c subtract the mean spectrum from each
c  Will need to output the mean spectrum!!!
      if (ifmeansub .ne. 0) then
         call submean(speclogarr,nspec,nspmax,nrest,avgspec)
      end if

c normalize spectra to modulus 1
      call normalize(speclogarr,nspec,nspmax,nrest)

c zero out all the bad-flag values
      call zerobad(speclogarr,nspec,nspmax,nrest,NLOGWMAX)

c have to transpose the matrix
      do j=1,NLOGWMAX
         do i=1,NSPECMAX
            specmatrix(j,i)=speclogarr(i,j)
         end do
      end do

c diagonalize matrix.  Note the stock Numerical Recipes svdcmp
c is limited to n=500 columns (ie 500 input spectra), will need
c to recompile.
      call svdcmp(specmatrix,nrest,nspec,NLOGWMAX,NSPECMAX,
     $     eigenvals,weights)

      open(2,file='pcatest.evals',status='unknown')
      do i=1,nspec
         write(2,*) i,eigenvals(i)
      end do
      close(2)

c  Note, this step is probably unneeded because the eigenvalues
c  are already sorted when returned.
c sort an index to eigenvalues (ascending order)
c so now eigenvals(indeigen(1)) is smallest eigenvalue
c        eigenvals(indeigen(nspec)) is largest eigenvalue
      call indexx(nspec,eigenvals,indeigen)

      open(2,file='pcatest.evals.sort',status='unknown')
      do i=nspec,1,-1
         write(2,*) nspec+1-i,indeigen(i),eigenvals(indeigen(i))
      end do
      close(2)

c plot the first N eigenvectors

      call pgsubp(2,2)
      nplot=10
      do k=1,nplot
         iv = indeigen(nspec+1-k)
         call showspec(nrest,waverest,specmatrix(1,iv))
         write(toplabel,'(a,i4,2x,i4,a,1pe10.3)') "eigenvector ",k,iv,
     $        " , eigenvalue ",eigenvals(iv)
         call pglabel("log rest wavelength","normalized to modulus 1",
     $        toplabel)
      end do
c      call pgpage()
c  this is a hack to force the next series onto a new _physical_ page
      call pgsubp(1,1)
      call pgsubp(2,2)
      do k=1,nplot
         iv = indeigen(nspec+1-k)
         call showspec(nrest,wrestlin,specmatrix(1,iv))
         write(toplabel,'(a,i4,2x,i4,a,1pe10.3)') "eigenvector ",k,iv,
     $        " , eigenvalue ",eigenvals(iv)
         call pglabel("rest wavelength","normalized to modulus 1",
     $        toplabel)
      end do
      call pgsubp(1,1)

c  write out mean spectrum and first 10 eigenvectors
      open(2,file='pcatest.evects',status='unknown')
      do i=1,nrest
         write(2,'(f7.5,11(2x,1pe10.3))') waverest(i),avgspec(i),
     $        (specmatrix(i,indeigen(nspec+1-k)),k=1,10)
      end do
      close(2)

      call pgend()

      end

