
c Do the type of fitting as in fitz3.f but only at the given
c redshift, and return the coeffs of the best fit, for a
c list of objects/1d-spectra

c rework fittemplates to fit each spectrum w/a series of 
c a single template broadened by various vel dispersions.

c fitdisp2 is like fitdisp but also does Monte Carlo realizations
c of the spectrum and refitting to get a hopefully more
c believable error bar.

      program fitdisp2

      include 'specarray.h'

c read spectrum names from a file?
      parameter (IFFILE=1)
c plot intermediate steps?
      parameter (IFPLOTALL=0)
c plot initial spectrum (first step)?
      parameter (IFPLOTFIRST=1)
c polynomial or median continuum subtraction?
      parameter(IFPOLYSUB=1)
c do cont sub/dev before blanking sky?
      parameter(IFSUBFIRST=0)
c debug level
      parameter(IFDEBUG=1)

c size of gaussian kernel for the broadening convolution
      parameter(NGAUSMAX=401)
c number of broadenings
      parameter(NBROADMAX=201)
c to conserve memory, don't allow 25 templates, just 5
c replace NMAXVECTS with NMAXTEMPL
c      parameter(NMAXTEMPL=5)
c max number of monte carlo realizations
      parameter (NMONTEMAX=201)

      common /galvectors/ koffset,z,nvmax,nwlmax,ifitvect(NMAXVECTS),
     $     evects(NMAXVECTS,NLOGWMAX)
c      real broadvects(NBROADMAX,NMAXVECTS,NLOGWMAX)
      real broadvects(NMAXVECTS,NBROADMAX,NLOGWMAX)
      real spec1(NWAVEMAX),espec1(NWAVEMAX)
      real spec2(NLOGWMAX),espec2(NLOGWMAX)
      real spec3(NLOGWMAX),espec3(NLOGWMAX)
      real specratio(NLOGWMAX),specratiofit(NLOGWMAX)
      real tmpwl(NLOGWMAX),tmpspecfit(NLOGWMAX)
      real tmpspec(NLOGWMAX),wavelog(NLOGWMAX),wavelin(NWAVEMAX)
      real gauskern(NGAUSMAX),tmpxkern(NGAUSMAX)
      real waverest1(NLOGWMAX)
      real wrestlin(NWAVEMAX),tmpspec2(NWAVEMAX)
      integer ifit(NLOGWMAX)
      real wfit(NLOGWMAX),sfit(NLOGWMAX),efit(NLOGWMAX)
      real rifit(NLOGWMAX),fitvect(NLOGWMAX)
      real smcfit(NLOGWMAX)
      real dispmonte(NMONTEMAX)
      real disparray(NBROADMAX),chidisp(NBROADMAX),rchidisp(NBROADMAX)
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
c      external gasdev

      include 'pcredshift.h'

      root2pi = sqrt(2.0*3.141593)
      root2 = sqrt(2.0)
      idum=-1
      nmonte=50

      call pgbeg(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.3)
      nxplot=1
      nyplot=1

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
      nvbroadmax=NBROADMAX

      nvbroad = 121
      dispstep = 5.

      write (*,'("Read pre-smoothed template set? 0/1 ",$)')
      read (*,*) ifpresmooth

c open file and read eigenvectors. nv is returned as
c the # of vectors to use, nw is the actual # of wavelength pixels
c so now evects is the unbroadened templates
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
      write(*,*) "w0rlog, dwrlog = ", w0rlog,dwrlog

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

 210  write(*,'("fit output file [fitdisp1.out]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fitdisp1.out'
      open(8,file=oname,status='unknown',err=210)

 215  write(*,'("disp output file [fitdisp2.out]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fitdisp2.out'
      open(9,file=oname,status='unknown',err=210)

c 217  write(*,'("chisq array output file [fitdisp3.out]: ",$)')
c      read(*,'(a)') oname
c      if (oname(1:3) .eq. '   ') oname='fitdisp3.out'
c      open(7,file=oname,status='unknown',err=210)

c old - as an expedient, construct a set of smoothed vectors
c by smoothing the first one, i.e. we are discarding all but
c the first vector

c smooth each of the input vectors in evects to fill up 
c broadvects, which has dimensions nvbroadmax,nvmax,nwlmax
c or nvmax,nvbroadmax,nwlmax

c smoothing by a dv const. w/wavelength is also smoothing by
c a d(log wavelength) const w/wave

c reading presmoothed does not yet work in the more-than-one-template
c version
      if (ifpresmooth .eq. 1) go to 180

c      dispstep = 10.
      disparray(1) = 0.0
      do k=2,nvbroad
         disparray(k) = disparray(1) + (k-1) * dispstep
      end do

c loop over the nv templates
      do itempl=1,nv
      call pgenv(waverest1(1),waverest1(nwlog),0.,1000.,0,1)
      write(tlabel,'("template ",i3)') itempl
      call pglabel("log wavelength","flux",tlabel)
      do j=1,nwlog
         tmpspec(j) = evects(itempl,j)
      end do
      call pgline(nwlog,waverest1,tmpspec)
c copy the first set in unchanged
      do j=1,nwlog
         broadvects(itempl,1,j) = evects(itempl,j)
      end do
c loop over broadenings
      do k = 2,nvbroad
         displogw = disparray(k) / 3.0e5 / 2.3026
         displogpix = displogw/dwrlog
         displogpixsq = displogpix**2
c should really move the entire gaussian convol to a subroutine
         iradkern = int(4.0*displogpix)
         iradkern = min(iradkern,(NGAUSMAX-1)/2)
         ngauskern = 2*iradkern+1
         ikernctr = 1+iradkern
         gkernsum = 0.0
         do j=1,ngauskern
c This was a direct gaussian, while we really ought to be
c using error functions.  Not sure this really makes a 
c signif difference.
            gauskern(j) = 1./root2pi/displogpix *
     $           exp(-(j-ikernctr)**2/2.0/displogpixsq)
c erf(x) is integral exp(-x^2)
c            if (j .ne. ikernctr) then
c               gauskern(j) = 1./root2pi/displogpix *
c     $              abs( erf((j+0.5-ikernctr)/displogpix/root2) -
c     $              erf((j-0.5-ikernctr)/displogpix/root2) )
c            else
c               gauskern(j) = 2. * 1./root2pi/displogpix *
c     $              abs( erf(-(j+0.5-ikernctr)/displogpix/root2) )
c            end if
c this ought to come out near 1
            gkernsum = gkernsum + gauskern(j)
         end do
c         if (IFDEBUG .ge. 1) then
c         ifplotkern=1
c         if (ifplotkern .ge. 1 .and. mod(k,20) .eq. 1) then
cc            write(*,*) k,disparray(k),displogpix,gkernsum
c            do j=1,ngauskern
c               tmpxkern(j) = (j-ikernctr)*displogpix
c            end do
c            call pgenv(tmpxkern(1),tmpxkern(ngauskern),
c     $           0.,gauskern(ikernctr)*1.05,0,1)
c            call pgline(ngauskern,tmpxkern,gauskern)
c            write(tlabel,'("disp ",f7.1," displogpix ",f5.2)') 
c     $           disparray(k),displogpix
c            call pglabel("log disp","kernel",tlabel)
c         end if
         do j=1,ngauskern
            gauskern(j) = gauskern(j) /gkernsum
         end do
c loop over wavelength, filling into the broadvects array
         do i=1,nwlog
            jmin = max(1,i-iradkern)
            jmax = min(nwlog,i+iradkern)
c            evects(k,i) = 0.0
            broadvects(itempl,k,i) = 0.0
            do j=jmin,jmax
c               evects(k,i) = evects(k,i) + 
c     $              evects(1,j)*gauskern(ikernctr+j-i)
               broadvects(itempl,k,i) = broadvects(itempl,k,i) + 
     $              evects(itempl,j)*gauskern(ikernctr+j-i)
            end do
         end do
         if (IFDEBUG .ge. 2 .and. mod(k,5) .eq. 0) then
            call pgenv(waverest1(1),waverest1(nwlog),0.,1000.,0,1)
            do j=1,nwlog
               tmpspec(j) = broadvects(itempl,k,j)
            end do
            call pgline(nwlog,waverest1,tmpspec)
         end if
      end do
c end the looping over templates
      end do
c write out smoothed vectors for next time
c this doesn't yet work in the multi-template version
      open(4,file='smoothvects.out',status='unknown')
      do j=1,nwlog
         write(4,'(f8.6,$)') waverest1(j)
         do itempl=1,nv
c            write(4,'(f8.2,$)') evects(itempl,j)
            write(4,'(f8.2,$)') broadvects(itempl,1,j)
         end do
         write(4,*)
      end do

 180  continue

c continuum subtract the template vectors
c dummy array of errors
      do i=1,nwlog
         tmpspec2(i) = 1.0
      end do
c loop over templates
      do itempl=1,nv
      do k=1,nvbroad
         do j=1,nwlog
c            tmpspec(j) = evects(k,j)
            tmpspec(j) = broadvects(itempl,k,j)
         end do
         if (ifcontsub .ne. 0) then
            contwrad = 100.*dwrlog
            call contdivpctile(tmpspec,nwlog,tmpspec2,dwrlog,contwrad)
c            if(IFPOLYSUB .ne. 0) then
cc            call contsubpoly(wavelog,spec3,espec3,nwlog,15)
cc            call contsubleg(wavelog,spec3,espec3,nwlog,15)
c               call contsubleg(waverest1,tmpspec,tmpspec2,
c     $           nwlog,15)
c            else
c               contwrad = 100.*dwrlog
c               call contsubmed(tmpspec,nwlog,dwrlog,contwrad)
c            end if
         else
c            call contsubconst(evects(k,1),nwlog)
            call contsubconst(broadvects(itempl,k,1),nwlog)
         end if
         do j=1,nwlog
c            evects(k,j) = tmpspec(j)
            broadvects(itempl,k,j) = tmpspec(j)
         end do
      end do
c end looping over templates
      end do
      
c in the multi-template version, this only writes the
c unbroadened (but continuum divided) vectors
      open(4,file='csubvects.out',status='unknown')
      do j=1,nwlog
         write(4,'(f8.6,$)') waverest1(j)
         do itempl=1,nv
            write(4,'(f8.2,$)') broadvects(itempl,1,j)
         end do
         write(4,*)
      end do

c plot some of the template vectors
      call pgsubp(2,2)
      do itempl=1,nv
      do k=1,nvbroad,20
         do j=1,nwlog
c            tmpspec(j) = evects(k,j)
            tmpspec(j) = broadvects(itempl,k,j)
         end do
         call showspec(nwlog,waverest1,tmpspec)
         write(tlabel,'(2(a,1x,i3))') "template",itempl," broadening ",k
         call pglabel("log wavelength"," ",tlabel)
      end do
      end do
cc advance to next full page
c      ipage = mod(4-nv,4)
c      do i=1,ipage
c         call pgpage
c      end do
      call pgsubp(2,2)

      call pgsubp(nxplot,nyplot)

c Here we start looping over objects
 220  continue

c open files and read spec and error, with wavelength calib.
      if (IFFILE .eq. 1) then
         read(10,'(a)',err=666,end=666) sname
         read(12,*,err=666,end=666) zreal
      end if

c label for tops of plot (objno would be better)
c use up to last 40 chars of spectrum name
      iend = index(sname,' ')-1
      istart = max(1,iend-39)
      write(namelabel,'(a)') sname(istart:iend)

c this linearizes the spectrum which may be less than ideal
c      call get1ddeepspec(sname,nwspec,spec1,espec1,w0sp,dwsp,ier)
      call get1dcoaddspec(sname,nwspec,spec1,espec1,w0sp,dwsp,ier)
      if (ier .ne. 0) then
c write something to the output files?
         write(*,'("couldnt open ",a)') sname
         write(8,'("couldnt open ",a)') sname
         write(9,'("couldnt open ",a)') sname
         go to 220
      end if

      write(*,*) "nwspec, w0sp, dwsp = ", nwspec,w0sp,dwsp
c zero out anything beyond the actual spectrum
      do i=nwspec+1,NWAVEMAX
         spec1(i) = 0.0
         espec1(i) = 0.0
      end do
c try to clean out bad values
c feb 11 2008 - added a check for very small error values, which
c   crept in in a test spectrum
c also add a computation of median SNR
      badmax = 10000.
      errmin = 1.0e-5
      nwspecgood=0
      do i=1,nwspec
         if(spec1(i) .lt. bad .or. spec1(i) .gt. badmax .or.
     $        espec1(i) .lt. errmin) then
            spec1(i) = badset
c         if(espec1(i) .lt. bad .or. espec1(i) .gt. badmax)
            espec1(i) = badset
         else
            nwspecgood=nwspecgood+1
            tmpspec2(nwspecgood) = spec1(i)/espec1(i)
         end if
      end do

c find the mean/median SNR
      snrmean = 0.0
      do i=1,nwspecgood
         snrmean = snrmean + tmpspec2(i)
      end do
      snrmean = snrmean / nwspecgood
      imed = nint(nwspecgood/2.0)
      snrmed = select(imed,nwspecgood,tmpspec2)
      write(*,'(" SNR mean, median: ",f7.3,1x,f7.3)') snrmean,snrmed

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
      if (IFPLOTALL .ne. 0 .or. IFPLOTFIRST .ne. 0 
     $     .and. ndiamed .gt. 1) then
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

c do continuum subtraction before blanking, which at the
c moment means doing it before log rebinning?  Or not
     
      if (ifcontsub .ne. 0 .and. IFSUBFIRST .eq. 1) then
c save the pre-subtractiion spectrum.  This is before log rebinning
c and blanking.  If we do this then we have to log rebin
c tmpspec2 into specratio to be consistent with later.
         do i=1,nwspec
            tmpspec2(i) = spec2(i)
         end do
         call logrebin(tmpspec2,nwspec,w0sp,dwsp,nwlog,w0rlog,dwrlog,
     $        specratio)
c note dwlog not defined yet, could use dwrlog
c         contwrad = 100.*dwlog
c         call contdivpctile(spec3,nwlog,espec3,dwlog,contwrad)
         contwrad = 150.*dwsp
         call contdivpctile(spec2,nwspec,espec2,dwsp,contwrad)
c         if(IFPOLYSUB .ne. 0) then
cc            call contsubpoly(wavelog,spec3,espec3,nwlog,15)
cc            call contsubleg(wavelog,spec3,espec3,nwlog,15)
c            call contsubleg(wavelog(imin),spec3(imin),espec3(imin),
c     $           imax-imin+1,15)
c         else
c            contwrad = 100.*dwlog
c            call contsubmed(spec3,nwlog,dwlog,contwrad)
c         end if
      else
c         call contsubconst(spec3,nwlog)
         call contsubconst(spec2,nwsp)
      end if


c easier to deal with cont sub, etc if we blank later, after
c cont subtraction etc.
c      call blankoutsky2(spec1,1,1,nwspec,wavelin)
c      call selectabsregions(spec1,1,1,nwspec,wavelin,zreal)
      call blankoutsky2(spec2,1,1,nwspec,wavelin)
c Blank out emission lines too.  This needs to know
c z so we can blank in the restframe
      call blankoutlines(spec2,1,1,nwspec,wavelin,zreal)
c      call selectabsregions(spec2,1,1,nwspec,wavelin,zreal,0)
      if (IFPLOTALL .ne. 0) then
         call showspec(nwspec,wavelin,spec2)
         call pglabel("wavelength","counts","sky blanked")
      end if

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

c cont sub if doing it later
      if (ifcontsub .ne. 0 .and. IFSUBFIRST .eq. 0) then
c save the pre-division spectrum.  This is after log rebinning
         do i=1,nwlog
            specratio(i) = spec3(i)
         end do
c contdivpctile divides or subtracts by the 90th %ile of pixels
c in some running bin.  Currently set up to have it divide.
         contwrad = 150.*dwlog
         call contdivpctile(spec3,nwlog,espec3,dwlog,contwrad)
c         contwrad = 100.*dwsp
c         call contdivpctile(spec2,nwspec,espec2,dwsp,contwrad)
c         if(IFPOLYSUB .ne. 0) then
cc            call contsubpoly(wavelog,spec3,espec3,nwlog,15)
cc            call contsubleg(wavelog,spec3,espec3,nwlog,15)
c            call contsubleg(wavelog(imin),spec3(imin),espec3(imin),
c     $           imax-imin+1,15)
c         else
c            contwrad = 100.*dwlog
c            call contsubmed(spec3,nwlog,dwlog,contwrad)
c         end if
      else
         call contsubconst(spec3,nwlog)
c         call contsubconst(spec2,nwsp)
      end if

c Set specratio to be the ratio of pre-continuum removal to after.
      do i=1,nwlog
         specratio(i) = specratio(i)/spec3(i)
      end do
         
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

c select regions with abs lines, in log wavelength.
c at this point the spectrum is already shifted to rest,
c so send in z=0.0
c      call selectabsregions(spec3,1,1,nwlog,wavelog,0.0,1)

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
c This is for plotting later
            specratiofit(npfit) = specratio(i)
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
         npfittmp = 0
         write(*,'(a,3(2x,f8.2))') "best disp and err range: ",
     $     -999,-10000,10000
         write(9,'(a40,3(2x,f8.2))') namelabel,
     $     -999,-10000,10000
         write(8,'(1pe12.5,$)') -999
         do j=1,nv
            write(8,'(2(1x,1pe10.3),$)') -999.,10000
         end do
         write(8,'(2x,f8.5)') zreal
         go to 220
      end if

c don't really need to show this if we later overplot best fit
c      call showspec(npfit,wfit,sfit)
c      call pgqci(indexc)
c      call pgsci(3)
c      call pgline(npfit,wfit,efit)
c      call pgsci(indexc)
c      call pglabel("log rest wavelength","counts",
c     $     "spec to fit, "//namelabel)

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

c      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
c     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
cc      chifit(kindex) = chisq1
cc get error bars on coefficients
c      call svdvar(vv,nv,nvmax,ww,covar,nvmax)
c      do i=1,nv
c         dacoeff(i) = sqrt(covar(i,i))
c      end do

c for fitdisp, I really just want to fit a single template at a time
c and all I really need is one coefficient.  In principle there should
c not even be a DC offset since both are continuum subtracted.
c control which vector to fit by passing in an integer array
c through galvectors common block to galfuncs; this doesn't fit
c the dc offset

c for multi-template fitdisp, loop over broadening;
c at each value, copy the appropriately broadened templates
c into the evects array, which is what galfuncs actually has
c access to through the common block.  There is probably
c a more elegant and faster way to do this.

c if we are allowing the first fit coeff to be the DC level      
      nvfit = nv+1
c      nvfit = nv

      do i=1,nvfit
         ifitvect(i) = 1
      end do
      do i=nvfit+1,nvmax
         ifitvect(i) = 0
      end do

      do k=1,nvbroad
c         do i=1,nv
c            ifitvect(i) = 0
c         end do
c         ifitvect(k) = 1
         do itempl=1,nv
            do j=1,nwlog
               evects(itempl,j) = broadvects(itempl,k,j)
            end do
         end do
         call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
c     $     npfittmp,acoeff,nvmax,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
     $     npfittmp,acoeff,nvfit,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
c      chifit(kindex) = chisq1
      chidisp(k) = chisq1
c get error bars on coefficients - no real need for error bars here
c      call svdvar(vv,nvmax,nvmax,ww,covar,nvmax)
c      call svdvar(vv,nv,nvmax,ww,covar,nvmax)
c      do i=1,nv
c         dacoeff(i) = sqrt(covar(i,i))
c      end do
c end looping over broadening
      end do

      write(*,*) "used ",npfittmp," points"
      rchimin = 1.e10
      rchimax = -1.e10
      do k=1,nvbroad
         rchidisp(k) = chidisp(k) / npfittmp
         if (rchidisp(k) .lt. rchimin) kmin = k
         rchimin = min(rchimin,rchidisp(k))
c don't count ridiculous values
         if (rchidisp(k) .lt. 50.) rchimax = max(rchimax,rchidisp(k))
c         write(9,'(1pe13.6,1x,$)') rchidisp(k)
c         write(*,'(f4,2x,1pe13.6,)') disparray(k),rchidisp(k)
      end do
c      write(9,*)
      write(*,*) "Best chi-sq at ",kmin,disparray(kmin)
c now we have an array of chi-sq vs dispersion
c need to find the minimum and error range
c look for a delta-total-chisq of 1, scaled by whatever the
c min rchisq was, as if we scaled up the error bars to get
c rchisq to come out to 1.  If rchisq was <1, don't scale down.
c  The errors this way don't seem to be big enough.
      drchisq = 1.0/npfittmp * max(1.0, rchidisp(kmin))
      call findminerr2(nvbroad,disparray,rchidisp,drchisq,
     $     dispmin,disperrlo,disperrhi)
c plot the array and min
      call pgenv(disparray(1),disparray(nvbroad),rchimin,rchimax,0,1)
      call pglabel("velocity dispersion","reduced chi-squared",
     $     namelabel)
      call pgpoint(nvbroad,disparray,rchidisp,17)
      call pgqci(indexc)
      call pgsci(2)
      call pgpt1(disparray(kmin),rchidisp(kmin),12)
c      yplot = rchidisp(kmin)*1.1
      yplot = rchimin + 0.3*(rchimax-rchimin)
      call pgpt1(dispmin,yplot,17)
      call pgerr1(1,dispmin,yplot,disperrhi-dispmin,1.0)
      call pgerr1(3,dispmin,yplot,-disperrlo+dispmin,1.0)
      call pgsci(indexc)
      write(*,'(a,3(2x,f8.2))') "best disp and err range: ",
     $     dispmin,disperrlo,disperrhi
c      write(9,'(a40,3(2x,f8.2))') namelabel,
c     $     dispmin,disperrlo,disperrhi


c To plot best fit vector:
c I think I need to figure out the min of rchi, go to that value
c of vel (really the proper index k) and rerun bigsvdfit to get acoeff
c right.
c      do i=1,nvmax
c         ifitvect(i) = 0
c      end do
c      ifitvect(kmin) = 1
         do itempl=1,nv
            do j=1,nwlog
               evects(itempl,j) = broadvects(itempl,kmin,j)
            end do
         end do
      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
c     $     npfittmp,acoeff,nvmax,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
     $     npfittmp,acoeff,nvfit,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
c get error bars on coeffs of the best fit
c write out best coeffs
      call svdvar(vv,nvfit,nvmax,ww,covar,nvmax)
      write(*,'(a,$)') "best-fit coeffs:"
      do i=1,nvfit
         dacoeff(i) = sqrt(covar(i,i))
         write(*,'(2x,1pe10.3," +-",1pe10.3,$)') acoeff(i),dacoeff(i)
      end do
      write(*,*)
c compute best fit spectrum
c should this be npfittmp???
      do i=1,npfit
c         call galfuncs(rifit(i),evals,nvmax)
         call galfuncs(rifit(i),evals,nvfit)
         fitvect(i) = 0.0
         do j=1,nvfit
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do

      write(8,'(1pe12.5,$)') chisq1/real(npfittmp)
      do j=1,nvfit
         write(8,'(2(1x,1pe10.3),$)') acoeff(j),dacoeff(j)
      end do
c      write(2,'(2x,f8.5,2x,a100)') zreal,sname
      write(8,'(2x,f8.5)') zreal

c make a lightly smoothed version of the data, for plotting
c make the smoothing constant in dw, so if we use a more finely
c sampled template we smooth more
      iwid = int(2.e-4/dwlog)
c smooth sfit and store in spec2
      call medsmooth(npfit,sfit,iwid,spec2)
c      call showspec(npfit,wfit,sfit)
      call showspec(npfit,wfit,spec2)
      call pglabel("log rest wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      write(tlabel,'(a,", z= ",f8.4)') "spectrum and best fit",zreal
c      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgmtxt('T',1.5,0.0,0.0,tlabel)
      call pgqci(indexc)
      call pgsci(2)
      call pgline(npfit,wfit,fitvect)
c      write(tlabel,'(a,f7.4)') "chi-sq, z=",zmin
c      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)

c overplot fitted vectors

c Try plotting with the continuum multiplied back in. 
c We could also try plotting on a linear wavelength scale

      do i=1,npfit
         tmpspec(i) = sfit(i) * specratiofit(i)
         tmpwl(i) = 10**wfit(i)
         tmpspecfit(i) = fitvect(i) * specratiofit(i)
      end do
c Would like to make a vector of the best-fit model that showed
c model through the blanked-out wavelengths
c      do i=1,nwlog
c         tmpwl2(i) = 10**wavelog(i)
c Can't quite do this since I would have to pass wavelengths into
c galfuncs somehow to get the evaluated model
c         tmpspecfit2(i) = 
c      end do

c smooth, reuse spec2 for temp storage of the smooth vector to plot      
      call medsmooth(npfit,tmpspec,iwid,spec2)
c plot in log wavelength
      call showspec(npfit,wfit,spec2)
      call pglabel("log rest wavelength","normalized counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      write(tlabel,'(a,", z= ",f8.4)') "spectrum and best fit",zreal
c      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgmtxt('T',1.5,0.0,0.0,tlabel)
      call pgqci(indexc)
      call pgsci(2)
      call pgline(npfit,wfit,tmpspecfit)
c      write(tlabel,'(a,f7.4)') "chi-sq, z=",zmin
c      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)

      
c do the overplot in linear wavelength
c Here I could plot the model even where the data was blanked out
c but constructing it is extra work and not done yet.
      call showspec(npfit,tmpwl,spec2)
      call pglabel("rest wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      write(tlabel,'(a,", z= ",f8.4)') "spectrum and best fit",zreal
c      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgmtxt('T',1.5,0.0,0.0,tlabel)
      call pgqci(indexc)
      call pgsci(2)
c      call pgline(npfit,tmpwl,tmpspecfit)
      call pgline(npfit,tmpwl,tmpspecfit)
c      write(tlabel,'(a,f7.4)') "chi-sq, z=",zmin
c      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)


c Try to do a Monte Carlo to find the errors on dispersion.
c Do this by taking the best fit model, rewiggling each point 
c by the error bar to make a new spectrum, and refitting that.
c      nmonte = 50
      do ii=1,nmonte
         do i=1,npfit
            smcfit(i) = fitvect(i) + gasdev(idum)*efit(i)
         end do
         do k=1,nvbroad
            do itempl=1,nv
               do j=1,nwlog
                evects(itempl,j) = broadvects(itempl,k,j)
               end do
            end do
            call bigsvdfit(rifit(iminfit),smcfit(iminfit),efit(iminfit),
c     $     npfittmp,acoeff,nvmax,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
     $      npfittmp,acoeff,nvfit,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
c      chifit(kindex) = chisq1
            chidisp(k) = chisq1
         end do
         rchimin = 1.e10
         rchimax = -1.e10
         do k=1,nvbroad
            rchidisp(k) = chidisp(k) / npfittmp
            if (rchidisp(k) .lt. rchimin) kmin = k
            rchimin = min(rchimin,rchidisp(k))
c don't count ridiculous values
            if (rchidisp(k) .lt. 50.) rchimax = max(rchimax,rchidisp(k))
c         write(9,'(1pe13.6,1x,$)') rchidisp(k)
c         write(*,'(f4,2x,1pe13.6,)') disparray(k),rchidisp(k)
         end do
         dispmonte(ii) = disparray(kmin)
      end do
      sumd = 0.0
      sumdsq = 0.0
      do ii=1,nmonte
         sumd = sumd + dispmonte(ii)
         sumdsq = sumdsq + dispmonte(ii)**2
      end do
      avgdispmonte = sumd/nmonte
      rmsdispmonte = sqrt(sumdsq/nmonte - avgdispmonte**2)
      write(9,'(a40,6(2x,f8.2))') namelabel,
     $     dispmin,disperrlo,disperrhi,avgdispmonte,rmsdispmonte,snrmed
      write(*,'(a40,6(2x,f8.2))') namelabel,
     $     dispmin,disperrlo,disperrhi,avgdispmonte,rmsdispmonte,snrmed

c necessary?
c      call pgsubp(nxplot,nyplot)

      go to 220

 666  continue

      close(8)
      close(9)
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

c if we are allowing a DC component the first func is a constant

      subroutine galfuncs(xfit,evals,nterms)

      include 'specarray.h'

      real evals(nterms)

      common /galvectors/ koffset,z,nvmax,nwlmax,ifitvect(NMAXVECTS),
     $     evects(NMAXVECTS,NLOGWMAX)

c find the index corresponding to the position xfit
c and retrieve the eigenvector values

c  this is for use with rifit
      ispec = nint(xfit)
      jvect = ispec-koffset
      if(jvect .ge. 1 .and. jvect .le. nwlmax) then
         evals(1) = 1.0
         do i=2,nterms
c allow ifitvect to control what vectors are actually fit
            if (ifitvect(i) .ne. 0) then
               evals(i) = evects(i-1,jvect)
            else
               evals(i) = 0.0
            end if
         end do
      else
         do i=1,nterms
            evals(i) = 0.0
         end do
      end if

      return
      end

c given an array of y and a delta-y, find the x where y is minimum
c and x error range

      subroutine findminerr(n,x,y,deltay,xmin,xerrlo,xerrhi)

      real x(n),y(n)

c first just find the min point
      ymin = 1.0e20
      do i=1,n
         if (y(i) .lt. ymin) then
            ymin = y(i)
            xmin = x(i)
            imin = i
         end if
      end do
      ythresh = ymin + deltay

c now step away to find where we go above ythresh.
c This assumes that the array is fairly smooth.  could
c also look for the min and max points above ythresh
c which is what findminerr2 does
      ilo=imin-1
 2010 continue
      if (y(ilo) .lt. ythresh .and. ilo .gt. 1) then
         ilo = ilo-1
         go to 2010
      end if
      if (ilo .eq. 1) then
         xerrlo = x(1)
      else
         frac = (ythresh-y(ilo)) / (y(ilo+1)-y(ilo))
         xerrlo = frac*x(ilo+1) + (1.-frac)*x(ilo)
      end if
c now step up
      ihi = imin+1
 2020 continue
      if (y(ihi) .lt. ythresh .and. ihi .lt. n) then
         ihi = ihi+1
         go to 2020
      end if
      if (ihi .eq. n) then
         xerrhi = x(n)
      else
         frac = (ythresh-y(ihi-1)) / (y(ihi)-y(ihi-1))
         xerrhi = frac*x(ihi) + (1.-frac)*x(ihi-1)
      end if

c it would be better to do the whole thing on a splined version     
      
      return
      end

c given an array of y and a delta-y, find the x where y is minimum
c and the x error range

      subroutine findminerr2(n,x,y,deltay,xmin,xerrlo,xerrhi)

      real x(n),y(n)

c first just find the min point
      ymin = 1.0e20
      do i=1,n
         if (y(i) .lt. ymin) then
            ymin = y(i)
            xmin = x(i)
            imin = i
         end if
      end do
      ythresh = ymin + deltay

c now step away to find where we go above ythresh.
c This assumes that the array is fairly smooth. 
c in this version, we go all the way to the end of teh
c array, in case there are secondary minima.
c Start with the xerrlo and xerrhi being at least 1 grid
c point away (I didn''t bother interpolating)
      xerrlo = x(imin-1)
      do i=imin-1,1,-1
         if (y(i) .lt. ythresh) xerrlo = x(i)
      end do
      xerrhi = x(imin+1)
      do i=imin+1,n
         if (y(i) .lt. ythresh) xerrhi = x(i)
      end do

      return
      end


