
c read a few eigenvectors, a spectrum, and try to fit z

      program fitz

      include 'specarray.h'

      parameter(NMAXVECTS=12)
      parameter(NMAXSTEP=2*NLOGWMAX)
c scale sky by estimating read and photon noise ?
      parameter (IFSCALENOISE=1)
c plot intermediate steps?
      parameter (IFPLOTALL=0)
c read spectrum names from a file?
      parameter (IFFILE=1)
c overplot the actual z, when known?
      parameter (IFREALZ=1)
c max number of minima to find and zestimates to return
      parameter (NMAXEST=12)

c      real evects(NMAXVECTS,NLOGWMAX)
      integer koffset,nvmax,nwlmax
      real z
      common /vectorfit/ koffset,z,nvmax,nwlmax,
     $     evects(NMAXVECTS,NLOGWMAX)

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
      character sname*60,esname*60,vname*60
      character answ*3,tlabel*60,xlabel*60
      character fsname*60,fename*60,fzname*60
      integer isize(7)

c  stuff for svdfit
      real uu(NLOGWMAX,NMAXVECTS)
      real vv(NMAXVECTS,NMAXVECTS),ww(NMAXVECTS)
      real acoeff(NMAXVECTS),evals(NMAXVECTS)
      real chifit(NMAXSTEP),rchifit(NMAXSTEP)

      external funcs
      
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

c      write(*,
c     $     '("Restwave spectrum region to use in PCA, min, max: ",$)')
c      read(*,*) pcawmin,pcawmax
c      pwminlog = log10(pcawmin)
c      pwmaxlog = log10(pcawmax)

      nvmax=NMAXVECTS
      nwlmax=NLOGWMAX

c open file and read eigenvectors. nvects is returned as
c the # of vectors to use, nw is the actual # of wavelength pixels
      call readevects(evects,nvmax,nwlmax,waverest1,nv,nw)
      nwlog = nw

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
c      call pgsubp(1,1)
      call pgsubp(2,2)

 200  continue
      if (IFFILE .eq. 1) then
         write(*,'("File with spectrum names: ",$)')
         read(*,'(a)') fsname
         open(10,file=fsname,status='old',err=200)
         write(*,'("File with sky/rms names: ",$)')
         read(*,'(a)') fename
         open(11,file=fename,status='old',err=200)
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
         end if
      end if

      open(2,file='fitz.out',status='unknown')
      open(3,file='fitz.out1',status='unknown')
      open(4,file='fitz.out2',status='unknown')
      
 220  continue
c open files and read spec and error, with wavelength calib.
      if (IFFILE .eq. 1) then
         read(10,'(a)',err=666,end=666) sname
         read(11,'(a)',err=666,end=666) esname
         zreal=1.0e6
         if (IFREALZ .eq. 1 .and. ifzfile .eq. 1) then
            read(12,*,err=222,end=222) zreal
         end if
 222     continue            
      else
         write(*,'("fits file with spectrum [quit]: ",$)')
         read(*,'(a)') sname
         if (sname(1:3) .eq. '   ') go to 666
         write(*,'("fits file with sky/rms [quit]: ",$)')
         read(*,'(a)') esname
         if (esname(1:3) .eq. '   ') go to 666
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

      call imopen(sname,1,imgs,ier)
      if (ier .ne. 0) go to 220
      call imgsiz(imgs,isize,idim,itype,ier)
c trim a buffer at end to get rid of the region where Marc encoded
c the slit number
      ibuff2 = 25
      nwspec = isize(1) - ibuff2
      call imopen(esname,1,imge,ier)
c      call imgsiz(imge,isize,idim,itype,ier)
      if (ier .ne. 0) go to 220
      
      call imgl1r(imgs,spec1,ier)
      call imgl1r(imge,espec1,ier)

c  find wavelength of pixel 0, and dw/dpix
      call imgkwr(imgs,'CRPIX1',refpix,ier)
c      if (ier .ne. 0) refpix = badset
      if (ier .ne. 0) then
         ier=0
         refpix = 0.
      end if
      call imgkwr(imgs,'CRVAL1',refw,ier)
      if (ier .ne. 0) refw = badset
      call imgkwr(imgs,'CDELT1',dwsp,ier)
      if (ier .ne. 0 .or. abs(dwsp) .lt. 1.e-3) then
         ier=0
         call imgkwr(imgs,'CD1_1',dwsp,ier)
         if (ier .ne. 0) dwsp = badset
      end if
      if (refpix .gt. bad .and. refw .gt. bad .and. dwsp .gt. bad) 
     $     then
         w0sp = refw - refpix*dwsp
c         dws = dwsp
      else
c  we should really do something here
c         w0s = badset
c         dws = badset
         write(*,
     $     '("Couldnt get w0 and dw for ",a," enter w0,dw: ",$)') sname
         read(*,*) w0sp,dwsp
      end if
      write(*,*) "w0sp, dwsp = ", w0sp,dwsp

      call imclos(imgs,ier)
      call imclos(imge,ier)

c zero out anything beyond the actual spectrum
      do i=nwspec+1,nw
         spec1(i) = 0.0
         espec1(i) = 0.0
      end do

c try to clean out bad values
      badmax = 5000.
      do i=1,nw
         if(spec1(i) .lt. bad .or. spec1(i) .gt. badmax)
     $        spec1(i) = badset
         if(espec1(i) .lt. bad .or. espec1(i) .gt. badmax)
     $        espec1(i) = badset
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
c      write(*,*) "w0log, dwlog = ", w0log,dwlog

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
            espec2(i) = espec1(i) / sqrt(7.0)
         end do
      end if

      do i=1,nw
         wavelin(i) = w0sp+i*dwsp
         wavelog(i) = log10(wavelin(i))
      end do
      if (IFPLOTALL .ne. 0) then
         call showspec(nwspec,wavelin,spec1)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwspec,wavelin,espec2)
         call pgsci(indexc)
         call pglabel("wavelength","counts","spectrum and error")
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
C of indep measurements is about pixels/2 if well sampled)
         tmp = sqrt(max(ndiamed/2.,1.))
         do i=1,nwspec
            espec2(i) = espec2(i) / tmp
         end do
      else
         do i=1,nwspec
            spec2(i) = spec1(i)
            espec2(i) = espec2(i)
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

c log rebin both.  log rebinning the std dev the same way is a
c total hack, thus it is better done on the variance

      call logrebin(spec2,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,spec3)
c      call logrebin(spec2(imingood),nwspecgood,w0sp,dwsp,
c     $     nwlog,w0log,dwlog,spec3)
c convert std dev to variance
      do i=1,nwlog
         espec3(i) = espec2(i)**2
      end do
      call logrebin(espec3,nwspec,w0sp,dwsp,nwlog,w0log,dwlog,espec2)
c      call logrebin(espec3(imingood),nwspecgood,w0sp,dwsp,
c     $     nwlog,w0log,dwlog,espec2)
c convert back
      do i=1,nwlog
         espec3(i) = sqrt(max(espec2(i),1.e-2))
      end do

      if (IFPLOTALL .ne. 0) then
         call showspec(nwlog,wavelog,spec3)
         call pglabel("log wavelength","counts","log rebinned")
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nwlog,wavelog,espec3)
         call pgsci(indexc)
      end if

c find&clean ends
      call findends(nwlog,spec3,imin,imax)
      call cleanends(nwlog,spec3,imin,imax)
      write(*,*) "imin, imax = ", imin,imax

c continuum subtract
      if (ifcontsub .ne. 0) then
         contwrad = 100.*dwlog
         call contsubmed(spec3,nwlog,dwlog,contwrad)
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

c      if (IFPLOTALL .ne. 0) then
         call showspec(npfit,wfit,sfit)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(npfit,wfit,efit)
         call pgsci(indexc)
         call pglabel("log wavelength","counts",
     $        "spectrum and error to fit, "//sname)
c      end if


c What are the good data regions?  The spectrum was
c good from imin to imax (before cleanends).  Eigenvector 1
c (actually the mean) is good from jmin to jmax.
      do j=1,nwlog
         spec2(j) = evects(1,j)
      end do
      call findends(nwlog,spec2,jmin,jmax)
c      write(*,*) 'imin,imax,jmin,jmax, k0off: ',imin,imax,jmin,jmax,
c     $     k0offset

c loop over steps in log w.  the spectrum, indexed by i,
c is offset to the red of the eigenvectors, indexed by j,
c by k steps: j+k = i, k = i-j
c  and note that the two scales are offset by k0offset, so when
c  z=0, i=j-k0offset.  Typically the spectrum starts redder
c  than the eigenv. so k0offset>0.  If we started at 
c  z=0, k would start at -koffset.
c log w = log wrest + log(1+z), so  k*dwlog = log(1+z)
c Start with an offset of z=-0.01
      kmin = int(log10(1-0.01)/dwlog) - k0offset
c For pos. redshift, i+k0offset>j (i is position in spectrum, j in evect)
c Go to the point where only 10 pts overlap, which doesn't
c depend on k0offset.
      kmax = imax-10-jmin
      nk=kmax-kmin+1
c      write(*,*) "kmin, kmax: ",kmin,kmax
      tmp1 = (kmin+k0offset)*dwlog
      tmp2 = (kmax+k0offset)*dwlog
c      write(*,*) "log(1+z): ",tmp1,tmp2,"  z: ",
c     $     10**tmp1-1.0,10**tmp2-1.0

c      open(2,file='fitz.out1',status='unknown')

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
            
c now we gotta only use the matching points in the eigenvectors,
c which is a PITA.  Actually we've solved this through the way
c that funcs() works.

c do fit using svdfit or something like it.
c Each of the functions will be an eigenvector

c         call bigsvdfit(wfit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $        nwlmax,nvmax,chisq1,funcs)
c         call bigsvdfit(rifit,sfit,efit,npfit,acoeff,nv,uu,vv,ww,
c     $        nwlmax,nvmax,chisq1,funcs)
         call bigsvdfit(rifit(iminfit),sfit(iminfit),efit(iminfit),
     $        npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,funcs)
         chifit(kindex) = chisq1
c to get rid of division by zero when there were no points
c to fit.
         rchifit(kindex) = chisq1/max(real(npfittmp),1.e-4)
c         write(2,*) k,kindex,npfittmp,zp1log(kindex),zarr(kindex),
c     $        chifit(kindex),rchifit(kindex)

      end do
c      close(2)

c find the absolute minimum in chi-squared.  there are probably
c better ways to do this in the long run
c      call findabsmin(nk,chifit,indkmin,chimin)
c      call findabsmin(nk,rchifit,indkrmin,rchimin)
c  try finding the point with the biggest drop below "continuum"
c  i.e. deepest local minimum
      call findlocmin(nk,chifit,indkmin,chimin,chiminavg)
      call findlocmin(nk,rchifit,indkrmin,rchimin,rchiminavg)

c find the n deepest local minima in chi squared
      nfind=5
      call findmultmin(nk,chifit,nfind,indkarr,diffarr)
      
      zmin=zarr(indkmin)
      zrmin=zarr(indkrmin)
      do i=1,nfind
         zestarr(i) = zarr(indkarr(i))
      end do

      write(*,*) "chisq min at ",indkmin,zmin,chimin,chiminavg
      write(*,*) "reduced chisq min at ",indkrmin,zrmin,rchimin,
     $     rchiminavg

c write stuff to output files.
c 2 - fitz.out gets # of z-estimates, z-estimates, spec name
c 3 - fitz.out1 gets specname(truncated), z(deepest chisq min) ,
c          z(deepest rchisq min), chi-min, chi-"continuum",
c          rchi-min, rchi-"continuum"
c 4 - fitz.out2 gets for each z-est: k-index, z-est, depth in chisq

c      write(3,'(a,f7.4,3x,f7.4)') sname,zmin,zrmin
      write(3,'(a25,2(3x,f7.4),3x,4(1x,1pe10.3))') 
     $     sname,zmin,zrmin,chimin,chiminavg,rchimin,rchiminavg
      write(2,'(i2,$)') nfind
      do i=1,nfind
         write(2,'(2x,f7.4$)') zestarr(i)
         write(4,'(i5,2x,f7.4,2x,f9.1,$)') indkarr(i),zestarr(i),
     $        diffarr(i)
      end do
      write(2,'(2x,a25)') sname
      write(4,'(2x,a25)') sname

c plot
c      call showspec(nk,zp1log,chifit)
c      call plotz(log10(1.+zreal),log10(1.+zmin),log10(1+zrmin))
c      write(xlabel,1010) "log(1+z)",zreal,zmin,zrmin
c      call pglabel(xlabel,"chi-squared",sname)
      call showspec(nk,zarr,chifit)
      call plotz(zreal,zmin,zrmin)
      write(xlabel,1010) "z",zreal,zmin,zrmin
      call pglabel(xlabel,"chi-squared",sname)
      call showspec(nk,zarr,rchifit)
      call plotz(zreal,zmin,zrmin)
      write(xlabel,1010) "z",zreal,zmin,zrmin
      call pglabel(xlabel,"reduced chi-squared",sname)
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
     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,funcs)
      write(*,*) "coeffs ",(acoeff(i),i=1,nv)
      do i=1,npfit
         call funcs(rifit(i),evals,nv)
         fitvect(i) = 0.0
         do j=1,nv
            fitvect(i) = fitvect(i) + acoeff(j)*evals(j)
         end do
      end do

      call showspec(npfit,wfit,sfit)
      call pglabel("log wavelength","counts"," ")
      call pgmtxt('T',2.5,0.0,0.0,sname)
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
     $     npfittmp,acoeff,nv,uu,vv,ww,nwlmax,nvmax,chisq1,funcs)
c      write(*,*) "coeffs ",(acoeff(i),i=1,nv)
      do i=1,npfit
         call funcs(rifit(i),evals,nv)
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


      go to 220

 666  continue

      close(2)
      close(3)
      close(4)
      close(10)
      close(11)
      if (ifzfile .ne. 0) close(12)

      call pgend()

      end

cccccccccccccccccccc
c  funcs returns the values of the nv eigenvectors
c  evaluated at position xfit, in the array evals

c  if svdfit is called with argument wfit, xfit is the
c  log wavelength.
c  if svdfit is called with argument rifit, xfit is the 
c  index in the log wavelength array, as a real. 
c  

      subroutine funcs(xfit,evals,nv)

      include 'specarray.h'
c having NMAXVECTS in here as well as fitz is a bad kludge
      parameter(NMAXVECTS=12)

      real evals(nv)

      common /vectorfit/ koffset,z,nvmax,nwlmax,
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

      
