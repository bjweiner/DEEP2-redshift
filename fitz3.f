
c read a few eigenvectors, a spectrum, and try to fit z

c fitz2 - this version is intended to read DEEP format 1-d spectra
c fitz3 - this version is to use polynomial continuum sub,
c   go over a fixed range of z, intended to use long templates
c   padded with zeros, and eventually incorporate a separate
c   looop with stellar templates.

      program fitz3

      include 'specarray.h'

c this is now moved to specarray.h
c      parameter(NMAXVECTS=12,NSTARVECTS=25)
      parameter(NMAXSTEP=2*NLOGWMAX)
c scale sky by e2stimating read and photon noise ?
      parameter (IFSCALENOISE=0)
c plot intermediate steps?
      parameter (IFPLOTALL=0)
c write intermediate spectra?
      parameter (IFWRITEALL=0)
c read spectrum names from a file?
      parameter (IFFILE=1)
c overplot the actual z, when known?
      parameter (IFREALZ=1)
c ask for slit and object numbers?
      parameter (IFOBJNO=1)
c max number of minima to find and zestimates to return
      parameter (NMAXEST=12)
c polynomial or median continuum subtraction?
      parameter(IFPOLYSUB=1)
c number of minima to find
      parameter(NMINFIND=5)

c      real evects(NMAXVECTS,NLOGWMAX)
      integer koffset,nvmax,nwlmax
      integer koffsetstar,nvmaxstar,nwlmaxstar
      real z,zstar
      common /galvectors/ koffset,z,nvmax,nwlmax,
     $     evects(NMAXVECTS,NLOGWMAX)
      common /starvectors/ koffsetstar,zstar,nvmaxstar,nwlmaxstar,
     $     starvects(NSTARVECTS,NLOGWMAX)

      real spec1(NWAVEMAX),espec1(NWAVEMAX)
      real spec2(NLOGWMAX),espec2(NLOGWMAX)
      real spec3(NLOGWMAX),espec3(NLOGWMAX)
      real tmpspec(NLOGWMAX),wavelog(NLOGWMAX),wavelin(NWAVEMAX)
      integer ifit(NLOGWMAX)
      real wfit(NLOGWMAX),sfit(NLOGWMAX),efit(NLOGWMAX)
      real rifit(NLOGWMAX),fitvect(NLOGWMAX)
      real waverest1(NLOGWMAX)
      real wrestlin(NWAVEMAX)
      real zp1log(NMAXSTEP),zarr(NMAXSTEP)
      real rnparr(NMAXSTEP)
      integer indkarr(NMAXEST)
      real diffarr(NMAXEST),zestarr(NMAXEST)
c      real xpl(2),ypl(2)
c      real zout(NLOGWMAX),zchisq(NLOGWMAX)
      character sname*160,esname*160,vname*160,oname*160
      character answ*3,tlabel*160,xlabel*160
      character namelabel*160,outfilename*80
      character fsname*160,fename*160,fzname*160,fobjname*160
      integer isize(7)
      character charobjno*60

      real waverest1star(NLOGWMAX)
      real zarrstar(NMAXSTEP)

c  stuff for svdfit.  note these are dimensioned to nmaxvects
c  but I also use them when fitting stars so they should really
c  be the bigger of (nmaxvects,nstarvects)
      real uu(NLOGWMAX,NMAXVECTS)
      real vv(NMAXVECTS,NMAXVECTS),ww(NMAXVECTS)
      real acoeff(NMAXVECTS),evals(NMAXVECTS)
      real chifit(NMAXSTEP),rchifit(NMAXSTEP)
      real chifitstar(NMAXSTEP),rchifitstar(NMAXSTEP)
      real doffit(NMAXSTEP),doffitstar(NMAXSTEP)

      external galfuncs
      external starfuncs
      
      include 'pcredshift.h'

c      call pgbeg(0,'?',2,2)
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

c      write(*,
c     $     '("Restwave spectrum region to use in PCA, min, max: ",$)')
c      read(*,*) pcawmin,pcawmax
c      pwminlog = log10(pcawmin)
c      pwmaxlog = log10(pcawmax)

      nvmax=NMAXVECTS
      nvmaxstar=NSTARVECTS
      nwlmax=NLOGWMAX
      nwlmaxstar=NLOGWMAX
      nfind=NMINFIND

c open file and read eigenvectors. nvects is returned as
c the # of vectors to use, nw is the actual # of wavelength pixels
      call readevects(evects,nvmax,nwlmax,waverest1,nv,nw)
c      write(*,*) "nv, nw = ",nv,nw
      nwlog = nw
c read star eigenvectors
      call readevects(starvects,nvmaxstar,nwlmaxstar,waverest1star,
     $     nvstar,nwstar)
      if (waverest1star(1) .ne. waverest1(1) .or. nwstar .ne. nw) then
c I should really have it deal with this problem rather than
c complaining.
         write(*,'(a,a,i5,2x,i5,2(2x,f11.4)) ')
     $   "Warning, star and galaxy templates are not ",
     $     "on same wavelength ",nw,nwstar,waverest1(1),waverest1star(1)
      end if


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
c plot the star eigenvectors
      do i=1,nvstar
         do j=1,nwlog
            tmpspec(j) = starvects(i,j)
         end do
         call showspec(nwlog,waverest1,tmpspec)
         write(tlabel,'(a,1x,i3)') "star eigenvector",i-1
         call pglabel("log wavelength"," ",tlabel)
      end do

      call pgsubp(nxplot,nyplot)

 200  continue
      if (IFFILE .eq. 1) then
         write(*,'("File with spectrum names: ",$)')
         read(*,'(a)') fsname
         open(10,file=fsname,status='old',err=200)
c         write(*,'("File with sky/rms names: ",$)')
c         read(*,'(a)') fename
c         open(11,file=fename,status='old',err=200)
         if (IFREALZ .eq. 1) then
            write(*,'("File with actual zs for plot [none]: ",$)')
            read(*,'(a)') fzname
            if (fzname(1:3) .ne. '   ') then
               open(12,file=fzname,status='old',err=202)
               ifzfile=1
            else
c  no actual z file
               ifzfile=0
            end if
 202        continue
         else
            ifzfile=0
         end if
         if (IFOBJNO .ne. 0) then
            write(*,'("File with slit.object numbers: ",$)')
            read(*,'(a)') fobjname
            if (fobjname(1:3) .ne. '   ') then
               open(15,file=fobjname,status='old',err=204)
               ifobjnofile=1
            else
c  no actual z file
               ifobjnofile=0
            end if
 204        continue
         else
            ifobjnofile=0
         end if            
      end if

 210  write(*,'("main output file [fitz.out]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fitz.out'
      open(2,file=oname,status='unknown',err=210)
 212  write(*,'("aux output file 1 [fitz.out1]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fitz.out1'
      open(3,file=oname,status='unknown',err=212)
 214  write(*,'("aux output file 2 [fitz.out2]: ",$)')
      read(*,'(a)') oname
      if (oname(1:3) .eq. '   ') oname='fitz.out2'
      open(4,file=oname,status='unknown',err=214)
c 215  write(*,'("aux output file 3 [fitz.out3]: ",$)')
c      read(*,'(a)') oname
c      if (oname(1:3) .eq. '   ') oname='fitz.out3'
      oname='fitz.out3'
      open(7,file=oname,status='unknown',err=214)

cccccccccccccccccccc
c Here we start looping over objects
      
 220  continue
c open files and read spec and error, with wavelength calib.
      if (IFFILE .eq. 1) then
         read(10,'(a)',err=666,end=666) sname
c         read(11,'(a)',err=666,end=666) esname
         zreal=1.0e6
         if (IFREALZ .eq. 1 .and. ifzfile .eq. 1) then
            read(12,*,err=222,end=222) zreal
         end if
c label for tops of plot.  the file we read should have one line
c for each object like 001.22024832
c 
         if (ifobjnofile .ne. 0) then
            read(15,'(a)',err=222,end=222) charobjno
            icharobjlen = index(charobjno,' ') - 1
            write(namelabel,'(a)') charobjno(1:icharobjlen)
            write(*,'("doing slit, object: ",a)') namelabel
         else
            write(namelabel,'(a)') sname
         end if
 222     continue     
      else
         write(*,'("fits file with spectrum [quit]: ",$)')
         read(*,'(a)') sname
         if (sname(1:3) .eq. '   ') go to 666
c         write(*,'("fits file with sky/rms [quit]: ",$)')
c         read(*,'(a)') esname
c         if (esname(1:3) .eq. '   ') go to 666
         if (IFREALZ .eq. 1) then
            write(*,'("actual z for plotting [none]: ",$)')
            read(*,'(a)') fzname
            zreal=1.0e6
            if (fzname(1:3) .ne. '   ') then
               read(fzname,*,err=232) zreal
            end if
 232        continue
         end if
      end if

      call get1ddeepspec(sname,nwspec,spec1,espec1,w0sp,dwsp,ier)
      if (ier .ne. 0) then
c write something to the output files?
         write(2,'("couldnt open ",a)') sname
         go to 220
      end if
c trim a buffer at end to get rid of the region where Marc encoded
c the slit number (in Allslits).  This may not be necessary for 
c those retrieved from the spec1d files.
      ibuff2 = 25
      nwspec = nwspec - ibuff2
      write(*,*) "nwspec, w0sp, dwsp = ", nwspec,w0sp,dwsp

c zero out anything beyond the actual spectrum
c      do i=nwspec+1,nw
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

c Figure out the ref wavelength for the log rebinning so that
c it falls on the wl scale of the eigenvectors.
      tmp = log10(w0sp)
      itmp = int( (tmp-w0rlog)/dwrlog )
      w0log = w0rlog + dwrlog*itmp
c k0offset is the offset in integer index of log wavelength 
c from the beginning of the eigenvectors to the spectrum to be fit.
      k0offset = itmp
      dwlog = dwrlog
      write(*,*) "w0log, dwlog = ", w0log,dwlog
      

c convert sky spectrum into std.dev spectrum.  This might be better
c done after smoothing and log rebinning?
c Placeholder assumptions: 2x30 min exposures, 5 pixel diam 
c extraction window.
      if (IFSCALENOISE .eq. 1) then
         nexp = 2
         exptime = 30.*60.
c         dpix = 5.
         dpix = 1.
         call scalenoise(espec1,nwspec,nexp,exptime,dpix,espec2)
      else
c or if the sky spectrum is actually a rms/per pixel spectrum
c we could do nothing, or scale down by a factor of 
c sqrt(#pixels).  #pixels is most likely 7.
         do i=1,nwspec
c            sqrtnpix=sqrt(7.0)
            sqrtnpix = 1.0
            espec2(i) = espec1(i) / sqrtnpix
         end do
      end if

c      do i=1,nw
      do i=1,nwspec
         wavelin(i) = w0sp+i*dwsp
c  this was wrong or not useful
c         wavelog(i) = log10(wavelin(i))
      end do
      if (IFPLOTALL .ne. 0) then
         call showspec(nwspec,wavelin,spec1)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwspec,wavelin,espec2)
         call pgsci(indexc)
         call pglabel("wavelength","counts","spectrum and error")
      end if
      if (IFWRITEALL .ne. 0) then
         if (ifobjnofile .ne. 0) then
            outfilename = charobjno(1:icharobjlen) // '.spec'
         else
            outfilename = 'tmp.out1'
         end if
         open(13,file=outfilename,status='unknown')
         do i=1,nwspec
c            write(13,*) wavelog(i),wavelin(i),spec1(i),espec2(i)
            write(13,*) wavelin(i),spec1(i),espec2(i)
         end do
         close(13)
      end if

c blank out?  blankoutsky is supposed to run on a 2-d array

      call blankoutsky(spec1,1,1,nwspec,wavelin)

cc find region of good data
c      call findends(nwspec,spec1,imin,imax)
c      imingood = imin
c      nwspecgood = imax-imin+1

c median smooth
      if (ndiamed .gt. 1) then
         call medsmooth(nwspec,spec1,ndiamed,spec2)
c         call medsmooth(nwspecgood,spec1(imingood),ndiamed,
c     $        spec2(imingood))
c attempt to compensate the errors.  I divided ndiamed by 2 as
c a hack since adjacent pixels are not independent (assume the #
c of indep measurements is about pixels/2 if well sampled)
c   It's not even clear I should do this since smoothed pixels
c   are correlated ...
c         tmp = sqrt(max(ndiamed/2.,1.))
c         do i=1,nwspec
c            espec2(i) = espec2(i) / tmp
c         end do
      else
         do i=1,nwspec
            spec2(i) = spec1(i)
            espec2(i) = espec1(i)
         end do
      end if

      if (IFPLOTALL .ne. 0) then
         call showspec(nwspec,wavelin,spec2)
         call pglabel("wavelength","counts","median smoothed")
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwspec,wavelin,espec2)
         call pgsci(indexc)
      end if
      if (IFWRITEALL .ne. 0) then
         if (ifobjnofile .ne. 0) then
            outfilename = charobjno(1:icharobjlen) // '.medsm'
         else
            outfilename = 'tmp.out2'
         end if
         open(13,file=outfilename,status='unknown')
         do i=1,nwspec
            write(13,*) wavelin(i),spec2(i),espec2(i)
         end do
         close(13)
      end if
      nwmedsmgood=0
      do i=1,nwspec
         if (spec2(i) .gt. bad) nwmedsmgood=nwmedsmgood+1
      end do

c 10.29.02 - This was old and wrong!
cc log rebin both.  log rebinning the std dev the same way is a
cc total hack, thus it is better done on the variance
c      call logrebin(spec2,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,spec3)
cc      call logrebin(spec2(imingood),nwspecgood,w0sp,dwsp,
cc     $     nwlog,w0log,dwlog,spec3)
cc convert std dev to variance
c      do i=1,nwspec
c         espec3(i) = espec2(i)**2
c      end do
c      call logrebin(espec3,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,espec2)
cc      call logrebin(espec3(imingood),nwspecgood,w0sp,dwsp,
cc     $     nwlog,w0log,dwlog,espec2)
cc convert back
c      do i=1,nwlog
c         espec3(i) = sqrt(max(espec2(i),1.e-2))
c      end do

c 10.29.02 This is new and hopefully correct
      call logrebinerr(spec2,espec2,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,
     $     spec3,espec3)
      do i=1,nwlog
         wavelog(i) = w0log + i*dwlog
      end do

      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wavelog,spec3)
         call pglabel("log wavelength","counts","log rebinned")
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwlog,wavelog,espec3)
         call pgsci(indexc)
      end if
      if (IFWRITEALL .ne. 0) then
         if (ifobjnofile .ne. 0) then
            outfilename = charobjno(1:icharobjlen) // '.logreb'
         else
            outfilename = 'tmp.out3'
         end if
         open(13,file=outfilename,status='unknown')
         do i=1,nwlog
            write(13,*) wavelog(i),spec3(i),espec3(i)
         end do
         close(13)
      end if
      nwlogrebgood=0
      do i=1,nwlog
         if (spec3(i) .gt. bad) nwlogrebgood=nwlogrebgood+1
      end do

c find&clean ends
c  11.07.02 - this may be wrong because findends excludes zeros
c     yet the spectrum might legitimately be zero or small.
      call findends(nwlog,spec3,imin,imax)
      call cleanends(nwlog,spec3,imin,imax)
      write(*,*) "imin, imax = ", imin,imax

c continuum subtract
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

      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wavelog,spec3)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwlog,wavelog,espec3)
         call pgsci(indexc)
         call pglabel("log wavelength","counts","continuum subtracted")
      end if
      if (IFWRITEALL .ne. 0) then
         if (ifobjnofile .ne. 0) then
            outfilename = charobjno(1:icharobjlen) // '.csub'
         else
            outfilename = 'tmp.out4'
         end if
         open(13,file=outfilename,status='unknown')
         do i=1,nwlog
            write(13,*) wavelog(i),spec3(i),espec3(i)
         end do
         close(13)
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
         do i=1,nfind
            zestarr(i) = 0.0
         end do
         go to 400
      end if

      call showspec(npfit,wfit,sfit)
      call pgqci(indexc)
      call pgsci(3)
      call pgline(npfit,wfit,efit)
      call pgsci(indexc)
      call pglabel("log wavelength","counts",
     $     "spectrum and error to fit, "//namelabel)

c What are the good data regions?  The spectrum was
c good from imin to imax (before cleanends).  Eigenvector 1
c (actually the mean) is good from jmin to jmax.
c  11.07.02 - this causes a problem because findends 
c  excludes zero pixels.  I substituted findbadends.
c  Of course one would like the eigenvectors not to have
c  bad data, anyway.
      do j=1,nwlog
         spec2(j) = evects(1,j)
      end do
c      call findends(nwlog,spec2,jmin,jmax)
      call findbadends(nwlog,spec2,jmin,jmax)
      write(*,*) 'imin,imax,jmin,jmax, k0off: ',imin,imax,jmin,jmax,
     $     k0offset

c This is old and comes from determining the range
c automatically rather than using a fixed range:

c loop over steps in log w.  the spectrum, indexed by i,
c is offset to the red of the eigenvectors, indexed by j,
c by k steps: j+k = i, k = i-j
c  and note that the two scales are offset by k0offset, so when
c  z=0, i=j-k0offset.  Typically the spectrum starts redder
c  than the eigenv. so k0offset>0.  If we started at 
c  z=0, k would start at -koffset.
c log w = log wrest + log(1+z), so  k*dwlog = log(1+z)
c Start with an offset of z=-0.01
c      kmin = int(log10(1-0.01)/dwlog) - k0offset
c For pos. redshift, i+k0offset>j (i is position in spectrum, j in evect)
c Go to the point where only 10 pts overlap, which doesn't
c depend on k0offset.
c      kmax = imax-10-jmin

c try going over fixed z range rather than determining it
c from the templates.  Probably want to run this with
c templates that are long, or padded with zeros.
      zfitmin=-0.01
      zfitmax=1.501
      kmin = int(log10(1+zfitmin)/dwlog) - k0offset
      kmax = int(log10(1+zfitmax)/dwlog) - k0offset

      nk=kmax-kmin+1
c      write(*,*) "kmin, kmax: ",kmin,kmax
      tmp1 = (kmin+k0offset)*dwlog
      tmp2 = (kmax+k0offset)*dwlog
c      write(*,*) "log(1+z): ",tmp1,tmp2,"  z: ",
c     $     10**tmp1-1.0,10**tmp2-1.0

c      open(2,file='fitz.out1',status='unknown')

 300  continue

      do k=kmin,kmax
c  kindex is 1 when k=kmin
         kindex = k-kmin+1
c  for passing k to the funcs subroutine
         koffset=k
         zp1log(kindex) = (k+k0offset)*dwlog
         zarr(kindex) = 10**zp1log(kindex) - 1.0
         z = zarr(kindex)

c  we should only use the data that overlaps the eigenvectors,
c  the eigenvectors extend from jmin to jmax, and i=j+k.
c  But, since we cleaned out bad points, the spectrum to fit
c  is no longer evenly spaced in dwlog, so ifitmax is not
c  trivial to calculate.  Bummer!  
c  that is, what I've been calling "i" is the index of spec3(i),
c  and the contents of ifit() and rifit(), but it is not the
c  index of ifit
         iminorig = jmin+k
         imaxorig = jmax+k
c now find the corresponding index of ifit - i.e. 
c we want iminfit where ifit(iminfit) >= iminorig.  Confused yet?
         call findindex(npfit,ifit,iminorig,imaxorig,iminfit,imaxfit)

         npfittmp = imaxfit-iminfit+1
         rnparr(kindex) = real(npfittmp)
c         write(*,*) "iminfit,imaxfit: ",iminfit,imaxfit,npfittmp
c I should test on npfittmp so that we bail if there are
c less than some minimum number.  It's not clear what to do
c for chisq, though.  If we use long or padded templates,
c this shouldn't happen.
         doffit(kindex) = real(npfittmp-nv)
         if (doffit(kindex) .gt. 0.99) then
c now we gotta only use the matching points in the eigenvectors,
c which is a PITA.  Actually we've solved this through the way
c that funcs() works.

c do fit using svdfit or something like it.
c Each of the functions will be an eigenvector

c            call bigsvdfit(wfit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $        nwlmax,nvmax,chisq1,funcs)
c            call bigsvdfit(rifit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $        nwlmax,nvmax,chisq1,funcs)
            call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $        npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
            chifit(kindex) = chisq1
c to get rid of division by zero when there were no points
c to fit.
c            rchifit(kindex) = chisq1/max(real(npfittmp),1.e-4)
            rchifit(kindex) = chisq1/max(doffit(kindex),1.)
c            write(2,*) k,kindex,npfittmp,zp1log(kindex),zarr(kindex),
c     $        chifit(kindex),rchifit(kindex)
         else
            if (kindex .gt. 1) then
               chifit(kindex) = chifit(kindex-1)
               rchifit(kindex) = rchifit(kindex-1)
            else
               chifit(kindex) = 1.e5
               rchifit(kindex) = 1.e5
            end if
         end if

      end do
c      close(2)

c--------------------
c do a smaller loop in z for stars
      zfitmin=-0.01
      zfitmax= 0.051
      kminstar = int(log10(1+zfitmin)/dwlog) - k0offset
      kmaxstar = int(log10(1+zfitmax)/dwlog) - k0offset
      nkstar=kmaxstar-kminstar+1
      do k=kminstar,kmaxstar
         kindex = k-kminstar+1
         koffsetstar=k
c         zp1log(kindex) = (k+k0offset)*dwlog
         zarrstar(kindex) = 10**zp1log(kindex) - 1.0
         zstar = zarrstar(kindex)
         iminorig = jmin+k
         imaxorig = jmax+k
         call findindex(npfit,ifit,iminorig,imaxorig,iminfit,imaxfit)
         npfittmp = imaxfit-iminfit+1
         rnparr(kindex) = real(npfittmp)
         doffitstar(kindex) = real(npfittmp-nvstar)
         if (doffitstar(kindex) .gt. 0.99) then
            call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $        npfittmp,acoeff,nvstar,uu,vv,ww,nwlmaxstar,nvmaxstar,
     $           chisq1,starfuncs)
            chifitstar(kindex) = chisq1
            rchifitstar(kindex) = chisq1/max(doffitstar(kindex),1.)
         else
            if (kindex .gt. 1) then
               chifitstar(kindex) = chifitstar(kindex-1)
               rchifitstar(kindex) = rchifitstar(kindex-1)
            else
               chifitstar(kindex) = 1.e5
               rchifitstar(kindex) = 1.e5
            end if
         end if
      end do

c--------------------
      if (ifobjnofile .ne. 0) then
c write chisq arrays to files
         outfilename = charobjno(1:icharobjlen) // '.chisq'
         open(14,file=outfilename,status='unknown')
         do k=1,nk
            write(14,*) k,doffit(k),chifit(k),rchifit(k)
         end do
         close(14)
         outfilename = charobjno(1:icharobjlen) // '.chisq.star'
         open(14,file=outfilename,status='unknown')
         do k=1,nkstar
            write(14,*) k,doffitstar(k),chifitstar(k),rchifitstar(k)
         end do
         close(14)
      end if      

c find the absolute minimum in chi-squared.  there are probably
c better ways to do this in the long run
c      call findabsmin(nk,chifit,indkmin,chimin)
c      call findabsmin(nk,rchifit,indkrmin,rchimin)
c  try finding the point with the biggest drop below "continuum"
c  i.e. deepest local minimum
      call findlocmin(nk,chifit,indkmin,chimin,chiminavg)
      call findlocmin(nk,rchifit,indkrmin,rchimin,rchiminavg)

c find the n deepest local minima in chi squared
      call findmultmin(nk,chifit,nfind,indkarr,diffarr)
      
      zmin=zarr(indkmin)
      zrmin=zarr(indkrmin)
      do i=1,nfind
c this is just the index of the closest minimum
c         zestarr(i) = zarr(indkarr(i))
c this finds the centroid of the minimum
         zestarr(i) = centroidmin(nk,zarr,chifit,indkarr(i))
c want to write out the coeffs(eigenvalues) of the best fits
c but this requires storing them ,or recalculating the fit
c         do j=1,nv
c            coeffmin(i,j) =
c         end do
      end do

      write(*,*) "chisq min at ",indkmin,zmin,chimin,chiminavg
      write(*,*) "reduced chisq min at ",indkrmin,zrmin,rchimin,
     $     rchiminavg

c find the best star fit
      call findlocmin(nkstar,chifitstar,indkminstar,chiminstar,
     $     chiminavgstar)
      call findlocmin(nkstar,rchifitstar,indkrminstar,rchiminstar,
     $     rchiminavgstar)
      zstarrmin=zarrstar(indkrminstar)
c if the best star fit is 0.01 in rchisq better than the best gal
c fit, and it is at reasonably low z, promote it to best fit.
c note that there can be a wide minimum in rchisq for stars
c so the z test is fairly relaxed.  it would be better to test
c on deltarchisq*degrees of freedom, perhaps, but too complex -
c the 0.01 is a rule of thumb I got from looking at mask 3206.
      if(rchiminstar .lt. (rchimin-0.01)
     $     .and. abs(zstarrmin) .lt. 0.005) then
         zrmin = zstarrmin
         indkmin = indkminstar
         indkrmin = indkrminstar
         chimin = chiminstar
         chiminavg = chiminavgstar
         rchimin = rchiminstar
         rchiminavg = rchiminavgstar
         do i=2,nfind
            zestarr(i) = zestarr(i-1)
            indkarr(i) = indkarr(i-1)
            diffarr(i) = diffarr(i-1)
         end do
         zestarr(1) = zrmin
         indkarr(1) = indkmin
         diffarr(1) = chiminavg - chimin
      end if

 400  continue
c write stuff to output files.
c 2 - fitz.out gets # of z-estimates, z-estimates, spec name
c 3 - fitz.out1 gets specname(truncated), z(deepest chisq min) ,
c          z(deepest rchisq min), chi-min, chi-"continuum",
c          rchi-min, rchi-"continuum"
c 4 - fitz.out2 gets for each z-est: k-index, z-est, depth in chisq
c 7 - fitz.out3 is for storing coeffs/eigenvalues of best fits

c      write(3,'(a,f7.4,3x,f7.4)') sname,zmin,zrmin
      write(3,'(a25,2(3x,f7.4),3x,4(1x,1pe10.3))') 
     $     sname,zmin,zrmin,chimin,chiminavg,rchimin,rchiminavg
      write(2,'(i2,$)') nfind
      do i=1,nfind
         write(2,'(2x,f7.4$)') zestarr(i)
         write(4,'(i5,2x,f7.4,2x,f9.1,$)') indkarr(i),zestarr(i),
     $        diffarr(i)
c         write(7,'(i2,$)') i
c         do j=1,nv
c            write(7,'(1x,1pe10.3,$)') coeffmin(i,j)
c         end do
      end do
      write(2,'(2x,a100)') sname
      write(4,'(2x,a100)') sname
c      write(7,'(2x,a100)') sname

c plot
c      call showspec(nk,zp1log,chifit)
c      call plotz(log10(1.+zreal),log10(1.+zmin),log10(1+zrmin))
c      write(xlabel,1010) "log(1+z)",zreal,zmin,zrmin
c      call pglabel(xlabel,"chi-squared",sname)
      call showspec(nk,zarr,chifit)
      call plotz(zreal,zmin,zrmin)
      write(xlabel,1010) "z",zreal,zmin,zrmin
      call pglabel(xlabel,"chi-squared",namelabel)
      call showspec(nk,zarr,rchifit)
      call plotz(zreal,zmin,zrmin)
      call pglabel(xlabel,"reduced chi-squared",namelabel)
 1010 format(a,",  z=",f7.4,"  zest1=",f7.4,"  zest2=",f7.4)
c      call showspec(nk,zarr,rnparr)
c      call pglabel("z","number of points in fit",sname)

c  calculate the fit at chimin.  Redo the fit to get the coeffs
c  for the best-fit spectrum
      koffset=indkmin+kmin-1
      call findindex(npfit,ifit,jmin+koffset,jmax+koffset,
     $     iminfit,imaxfit)
      npfittmp=imaxfit-iminfit+1
c      write(*,*) "iminfit,imaxfit: ",iminfit,imaxfit,npfittmp
c      call bigsvdfit(rifit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $     nwlmax,nvmax,chisq1,funcs)
      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
      write(*,*) "coeffs ",(acoeff(i),i=1,nv)
c should this be npfittmp???
      do i=1,npfit
         call galfuncs(rifit(i),evals,nv)
         fitvect(i) = 0.0
         do j=1,nv
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do
c HERE we can write out the best fit coeffs to a file
      write(7,'(i2,$)') 1
      do j=1,nv
            write(7,'(1x,1pe10.3,$)') acoeff(j)
         end do
      write(7,'(2x,a100)') sname

      call showspec(npfit,wfit,sfit)
      call pglabel("log wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgqci(indexc)
      call pgsci(3)
      call pgline(npfit,wfit,fitvect)
      write(tlabel,'(a,f7.4)') "chi-sq, z=",zmin
      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)

c  calculate the fit at rchimin
      koffset=indkrmin+kmin-1
      call findindex(npfit,ifit,jmin+koffset,jmax+koffset,
     $     iminfit,imaxfit)
      npfittmp=imaxfit-iminfit+1
c      call bigsvdfit(rifit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $     nwlmax,nvmax,chisq1,funcs)
      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,galfuncs)
c      write(*,*) "coeffs ",(acoeff(i),i=1,nv)
c should this be npfittmp???
      do i=1,npfit
         call galfuncs(rifit(i),evals,nv)
         fitvect(i) = 0.0
         do j=1,nv
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do
      call pgqci(indexc)
      call pgsci(2)
      call pgline(npfit,wfit,fitvect)
      write(tlabel,'(a,f7.4)') "red.chisq, z=",zrmin
      call pgmtxt('T',1.5,1.0,1.0,tlabel)
      call pgsci(indexc)

c new page 
      call pgsubp(nxplot,nyplot)

c do the star fits
      call findlocmin(nkstar,chifitstar,indkmin,chimin,chiminavg)
      call findlocmin(nkstar,rchifitstar,indkrmin,rchimin,rchiminavg)
c      nfind=5
c      call findmultmin(nkstar,chifitstar,nfind,indkarr,diffarr)
      zstarmin=zarrstar(indkmin)
      zstarrmin=zarrstar(indkrmin)
      write(*,*) "star chisq min at ",indkmin,zstarmin,chimin,chiminavg
      write(*,*) "star rchisq min at ",indkrmin,
     $     zstarrmin,rchimin,rchiminavg

c do the star plots
      call showspec(nkstar,zarrstar,chifitstar)
      call plotz(zreal,zstarmin,zstarrmin)
      write(xlabel,1010) "z",zreal,zstarmin,zstarrmin
      call pglabel(xlabel,"chi-squared",namelabel)
      call showspec(nkstar,zarrstar,rchifitstar)
      call plotz(zreal,zstarmin,zstarrmin)
      write(xlabel,1010) "z",zreal,zmin,zrmin
      call pglabel(xlabel,"reduced chi-squared",namelabel)
c  calculate the fit at chimin.  Redo the fit to get the coeffs
c  for the best-fit spectrum
      koffset=indkmin+kminstar-1
      call findindex(npfit,ifit,jmin+koffset,jmax+koffset,
     $     iminfit,imaxfit)
      npfittmp=imaxfit-iminfit+1
      call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $     npfittmp,acoeff,nvstar,uu,vv,ww,nwlmaxstar,nvmaxstar,chisq1,
     $     starfuncs)
      do i=1,npfit
         call starfuncs(rifit(i),evals,nvstar)
         fitvect(i) = 0.0
         do j=1,nvstar
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do
      call showspec(npfit,wfit,sfit)
      call pglabel("log wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,namelabel)
      call pgmtxt('T',1.5,0.0,0.0,"spectrum and best fits")
      call pgqci(indexc)
      call pgsci(3)
      call pgline(npfit,wfit,fitvect)
      write(tlabel,'(a,f7.4)') "chi-sq, z=",zminstar
      call pgmtxt('T',2.5,1.0,1.0,tlabel)
      call pgsci(indexc)

c new page
      call pgsubp(nxplot,nyplot)

      go to 220

 666  continue

      close(2)
      close(3)
      close(4)
      close(7)
      close(10)
      close(11)
      if (ifzfile .ne. 0) close(12)
      if (ifobjnofile .ne. 0) then
         close(15)
      end if

      call pgend()

      end

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

cccccccccccccccccccccccccccccccccccccccc
c the only difference between galfuncs and starfuncs
c is which set of vectors is getting passed in sideways
c through the common block
      subroutine starfuncs(xfit,evals,nv)

      include 'specarray.h'

      real evals(nv)

      common /starvectors/ koffset,z,nvmax,nwlmax,
     $     starvects(NSTARVECTS,NLOGWMAX)

c find the index corresponding to the position xfit
c and retrieve the eigenvector values

c  this is for use with rifit
      ispec = nint(xfit)
      jvect = ispec-koffset
      if(jvect .ge. 1 .and. jvect .le. nwlmax) then
         do i=1,nv
            evals(i) = starvects(i,jvect)
         end do
      else
         do i=1,nv
            evals(i) = 0.0
         end do
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccc
c  Draw vertical lines for the actual z and z-estimates
c  on the chi-squared plots.

      subroutine plotz(zreal,zest1,zest2)
    
      real xpl(2),ypl(2)
      integer indexc,indexls

      ypl(1) = -100.
      ypl(2) = 1.0e10
      call pgqci(indexc)
      call pgqls(indexls)
c      call pgsls(1)
      xpl(1) = zreal
      xpl(2) = zreal
      call pgsci(4)
      call pgsls(5)
      call pgline(2,xpl,ypl)
      xpl(1) = zest1
      xpl(2) = zest1
      call pgsci(3)
      call pgsls(4)
      call pgline(2,xpl,ypl)
      xpl(1) = zest2
      xpl(2) = zest2
      call pgsci(2)
      call pgsls(3)
      call pgline(2,xpl,ypl)
      call pgsci(indexc)
      call pgsls(indexls)

      return
      end

      
