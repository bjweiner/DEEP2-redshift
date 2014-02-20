
c Do the type of fitting as in fitz3.f but only at the given
c redshift, and return the coeffs of the best fit, for a
c list of objects/1d-spectra

      program fittemplates

      include 'specarray.h'

c read spectrum names from a file?
      parameter (IFFILE=1)
c plot intermediate steps?
      parameter (IFPLOTALL=0)
c plot initial spectrum (first step)?
      parameter (IFPLOTFIRST=1)
c polynomial or median continuum subtraction?
      parameter(IFPOLYSUB=1)

      common /galvectors/ koffset,z,nvmax,nwlmax,
     $     evects(NMAXVECTS,NLOGWMAX)
      real spec1(NWAVEMAX),espec1(NWAVEMAX)
      real spec2(NLOGWMAX),espec2(NLOGWMAX)
      real spec3(NLOGWMAX),espec3(NLOGWMAX)
      real tmpspec(NLOGWMAX),wavelog(NLOGWMAX),wavelin(NWAVEMAX)
      real waverest1(NLOGWMAX)
      real wrestlin(NWAVEMAX)
      integer ifit(NLOGWMAX)
      real wfit(NLOGWMAX),sfit(NLOGWMAX),efit(NLOGWMAX)
      real rifit(NLOGWMAX),fitvect(NLOGWMAX)
      character sname*160,esname*160,vname*160,oname*160
      character answ*3,tlabel*160,xlabel*160
      character namelabel*160,outfilename*80
      character fsname*160,fename*160,fzname*160,fobjname*160
c      integer isize(7)
c      character charobjno*60

c  stuff for svdfit.  note these are dimensioned to nmaxvects
c  but I also use them when fitting stars so they should really
c  be the bigger of (nmaxvects,nstarvects)
      real uu(NLOGWMAX,NMAXVECTS)
      real vv(NMAXVECTS,NMAXVECTS),ww(NMAXVECTS)
      real acoeff(NMAXVECTS),evals(NMAXVECTS)
c      real chifit(NMAXSTEP),rchifit(NMAXSTEP)
c      real chifitstar(NMAXSTEP),rchifitstar(NMAXSTEP)
      real dacoeff(NMAXVECTS),covar(NMAXVECTS,NMAXVECTS)
c      real doffit(NMAXSTEP),doffitstar(NMAXSTEP)

      external galfuncs

      include 'pcredshift.h'

      call pgbeg(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.3)
      nxplot=2
      nyplot=2

      write(*,'("Full continuum fit subtraction [1,y-yes]? ",$)')
      read(*,'(a1)') answ
      if (answ(1:1) .eq. '0' .or. answ(1:1) .eq. 'n' .or.
     $     answ(1:1) .eq. 'N') then
         ifcontsub = 0
      else
         ifcontsub = 1
      end if
 100  write(*,
     $     '("Diameter for median smoothing, pixels [0,1=none]: ",$)')
      read(*,'(a3)') answ
      if (answ(1:3) .eq. '   ') then
         ndiamed = 0
      else
         read(answ,*,err=100) ndiamed
      end if

      nvmax=NMAXVECTS
      nwlmax=NLOGWMAX

c open file and read eigenvectors. nvects is returned as
c the # of vectors to use, nw is the actual # of wavelength pixels
      call readevects(evects,nvmax,nwlmax,waverest1,nv,nw)
      nwlog = nw
      dwrlog = (waverest1(nw) - waverest1(1)) / (nwlog-1.)
      w0rlog = waverest1(1) - dwrlog
      do i=1,nwlog
         wrestlin(i) = 10**waverest1(i)
      end do

c figure out w0 and dw for the eigenvectors.  This is in log w
c There is no check yet to make sure it's evenly spaced.
      dwrlog = (waverest1(nw) - waverest1(1)) / (nwlog-1.)
      w0rlog = waverest1(1) - dwrlog
      do i=1,nwlog
         wrestlin(i) = 10**waverest1(i)
      end do
c      write(*,*) "w0rlog, dwrlog = ", w0rlog,dwrlog

c plot the eigenvectors
      call pgsubp(2,2)
      do i=1,nv
         do j=1,nwlog
            tmpspec(j) = evects(i,j)
         end do
         call showspec(nwlog,waverest1,tmpspec)
         write(tlabel,'(a,1x,i3)') "eigenvector",i-1
         call pglabel("log wavelength"," ",tlabel)
      end do
cc advance to next full page
c      ipage = mod(4-nv,4)
c      do i=1,ipage
c         call pgpage
c      end do
      call pgsubp(2,2)

      call pgsubp(nxplot,nyplot)


 200  continue
      if (IFFILE .eq. 1) then
         write(*,'("File with spectrum names: ",$)')
         read(*,'(a)') fsname
         open(10,file=fsname,status='old',err=200)

 202     write(*,'("File with actual redshifts: ",$)')
         read(*,'(a)') fzname
         open(12,file=fzname,status='old',err=202)
         ifzfile=1
      end if

 210  write(*,'("main output file [fittempl.out]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fittempl.out'
      open(2,file=oname,status='unknown',err=210)

c Here we start looping over objects
 220  continue

c open files and read spec and error, with wavelength calib.
      if (IFFILE .eq. 1) then
         read(10,'(a)',err=666,end=666) sname
         read(12,*,err=666,end=666) zreal
      end if

c label for tops of plot (objno would be better)
c use last 50 chars of spectrum name
      iend = index(sname,' ')-1
      write(namelabel,'(a)') sname(iend-49:iend)

c this linearizes the spectrum which may be less than ideal
      call get1ddeepspec(sname,nwspec,spec1,espec1,w0sp,dwsp,ier)
      if (ier .ne. 0) then
c write something to the output files?
         write(2,'("couldnt open ",a)') sname
         go to 220
      end if

      write(*,*) "nwspec, w0sp, dwsp = ", nwspec,w0sp,dwsp
c zero out anything beyond the actual spectrum
      do i=nwspec+1,NWAVEMAX
         spec1(i) = 0.0
         espec1(i) = 0.0
      end do
c try to clean out bad values
c 11.08.02 - added a check for very small error values, which
c   crept in in a test spectrum
      badmax = 10000.
      errmin = 1.0e-3
      nwspecgood=0
      do i=1,nwspec
         if(spec1(i) .lt. bad .or. spec1(i) .gt. badmax .or.
     $        espec1(i) .lt. errmin) then
            spec1(i) = badset
c         if(espec1(i) .lt. bad .or. espec1(i) .gt. badmax)
            espec1(i) = badset
         else
            nwspecgood=nwspecgood+1
         end if
      end do

      do i=1,nwspec
         wavelin(i) = w0sp+i*dwsp
c  this was wrong or not useful
c         wavelog(i) = log10(wavelin(i))
      end do
      if (IFPLOTALL .ne. 0 .or. IFPLOTFIRST .ne. 0) then
         call showspec(nwspec,wavelin,spec1)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwspec,wavelin,espec1)
         call pgsci(indexc)
         write (tlabel,'(a,", z= ",f8.4)') "spectrum and error",zreal
         call pglabel("wavelength","counts",tlabel)
      end if


      call blankoutsky(spec1,1,1,nwspec,wavelin)
      if (ndiamed .gt. 1) then
         call medsmooth(nwspec,spec1,ndiamed,spec2)
      else
         do i=1,nwspec
            spec2(i) = spec1(i)
         end do
      end if
      do i=1,nwspec
         espec2(i) = espec1(i)
      end do
      if (IFPLOTALL .ne. 0 .or. IFPLOTFIRST .ne. 0) then
         call showspec(nwspec,wavelin,spec2)
         call pglabel("wavelength","counts","median smoothed")
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwspec,wavelin,espec2)
         call pgsci(indexc)
      end if
      nwmedsmgood=0
      do i=1,nwspec
         if (spec2(i) .gt. bad) nwmedsmgood=nwmedsmgood+1
      end do

c change wavelength scale to zero redshift
      w0sp = w0sp / (1. + zreal)
      dwsp = dwsp / (1. + zreal)

cc Figure out the ref wavelength for the log rebinning so that
cc it falls on the wl scale of the eigenvectors.
c      tmp = log10(w0sp)
c      itmp = int( (tmp-w0rlog)/dwrlog )
c      w0log = w0rlog + dwrlog*itmp
cc k0offset is the offset in integer index of log wavelength 
cc from the beginning of the eigenvectors to the spectrum to be fit.
c      k0offset = itmp
c      dwlog = dwrlog
c Rather than doing the above, make it simple and rebin the
c spectrum onto the exact wavelength scale of the eigenvectors
      w0log = w0rlog
      dwlog = dwrlog
      itmp = 0
      k0offset = 0
      write(*,*) "w0log, dwlog = ", w0log,dwlog

      call logrebinerr(spec2,espec2,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,
     $     spec3,espec3)
      do i=1,nwlog
         wavelog(i) = w0log + i*dwlog
      end do
      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wavelog,spec3)
         call pglabel("log rest wavelength","counts","log rebinned")
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwlog,wavelog,espec3)
         call pgsci(indexc)
      end if
      nwlogrebgood=0
      do i=1,nwlog
         if (spec3(i) .gt. bad) nwlogrebgood=nwlogrebgood+1
      end do
      call findends(nwlog,spec3,imin,imax)
      call cleanends(nwlog,spec3,imin,imax)
      write(*,*) "imin, imax = ", imin,imax

      if (ifcontsub .ne. 0) then
         if(IFPOLYSUB .ne. 0) then
c            call contsubpoly(wavelog,spec3,espec3,nwlog,15)
c            call contsubleg(wavelog,spec3,espec3,nwlog,15)
            call contsubleg(wavelog(imin),spec3(imin),espec3(imin),
     $           imax-imin+1,15)
         else
            contwrad = 100.*dwlog
            call contsubmed(spec3,nwlog,dwlog,contwrad)
         end if
      else
         call contsubconst(spec3,nwlog)
      end if
c      write(*,*) 'got here'
      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wavelog,spec3)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwlog,wavelog,espec3)
         call pgsci(indexc)
         call pglabel("log rest wavelength","counts",
     $        "continuum subtracted")
      end if
      nwcsubgood=0
      do i=1,nwlog
         if (spec3(i) .gt. bad) nwcsubgood=nwcsubgood+1
      end do

c copy wave and spec to temp array, cleaning out bad points
      npfit=0
      do i=1,nwlog
         if(spec3(i) .gt. bad .and. espec3(i) .gt. bad) then
            npfit=npfit+1
            ifit(npfit) = i
            rifit(npfit) = real(i)
            wfit(npfit) = w0log + i*dwlog
            sfit(npfit) = spec3(i)
            efit(npfit) = espec3(i)
         end if
      end do
      write (*,*) npfit, " points to fit"
      write (*,'("npts by stage: ",6(3x,i6))') nwspec,nwspecgood,
     $     nwmedsmgood,nwlog,nwlogrebgood,nwcsubgood

      if (npfit .le. 0) then
         write(*,*) "Error! No points to fit!"
c         do i=1,nfind
c            zestarr(i) = 0.0
c         end do
c         go to 400
      end if

      call showspec(npfit,wfit,sfit)
      call pgqci(indexc)
      call pgsci(3)
      call pgline(npfit,wfit,efit)
      call pgsci(indexc)
      call pglabel("log rest wavelength","counts",
     $     "spec to fit, "//namelabel)

c The spectrum-to-fit may have gaps, while the templates should
c not.  In fitz3, rather than re-align the templates with the
c good points in the spectrum at every lag, I used an indirect
c reference to the data - the "x" over which the fit is made
c is rifit, which is just a real version of the array index;
c as opposed to using wfit, which is the actual wavelength.
c Then galfuncs uses this array index to find the template values.
c For the fittemplates application, since I'm only using a 
c zero lag, this isn't really necessary, but I'll stick with it.

c At least this means we don't have to do any complicated
c dorking around with indexes for min and max points to fit.
      iminfit = 1
      npfittmp = npfit
c koffset is the offset in index between spectrum-to-fit
c and templates, which is passed through to galfuncs in a common block.  
c Since I put the spectrum onto the log scale of the evects (zero lag),
c it is zero here.
      koffset = 0

      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
c      chifit(kindex) = chisq1
c get error bars on coefficients
      call svdvar(vv,nv,nvmax,ww,covar,nvmax)
      do i=1,nv
         dacoeff(i) = sqrt(covar(i,i))
      end do

c compute best fit spectrum
c should this be npfittmp???
      do i=1,npfit
         call galfuncs(rifit(i),evals,nv)
         fitvect(i) = 0.0
         do j=1,nv
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do

      write(2,'(1pe10.3,$)') chisq1/real(npfittmp)
      do j=1,nv
         write(2,'(2(1x,1pe10.3),$)') acoeff(j),dacoeff(j)
      end do
c      write(2,'(2x,f8.5,2x,a100)') zreal,sname
      write(2,'(2x,f8.5)') zreal

      call showspec(npfit,wfit,sfit)
      call pglabel("log rest wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      write(tlabel,'(a,", z= ",f8.4)') "spectrum and best fits",zreal
c      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgmtxt('T',1.5,0.0,0.0,tlabel)
      call pgqci(indexc)
      call pgsci(2)
      call pgline(npfit,wfit,fitvect)
c      write(tlabel,'(a,f7.4)') "chi-sq, z=",zmin
c      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)

c overplot fitted vectors

c necessary?
c      call pgsubp(nxplot,nyplot)

      go to 220

 666  continue

      close(2)
      call pgend()

      end

cccccccccccccccccccc

cccccccccccccccccccc
c  galfuncs returns the values of the nv eigenvectors
c  evaluated at position xfit, in the array evals

c  if svdfit is called with argument wfit, xfit is the
c  log wavelength.
c  if svdfit is called with argument rifit, xfit is the 
c  index in the log wavelength array, as a real. 
c  

      subroutine galfuncs(xfit,evals,nv)

      include 'specarray.h'

      real evals(nv)

      common /galvectors/ koffset,z,nvmax,nwlmax,
     $     evects(NMAXVECTS,NLOGWMAX)

c find the index corresponding to the position xfit
c and retrieve the eigenvector values

c  this is for use with rifit
      ispec = nint(xfit)
      jvect = ispec-koffset
      if(jvect .ge. 1 .and. jvect .le. nwlmax) then
         do i=1,nv
            evals(i) = evects(i,jvect)
         end do
      else
         do i=1,nv
            evals(i) = 0.0
         end do
      end if

      return
      end
