
c simple program to read a set of linefits - like an NNN.concat file -
c and calculate sigma according to some scheme.  provide several
c schemes, and make it easy to add schemes so we can fool
c around with them.  this replaces compsigma.sh and compsigma_wave.sh,
c which worked fine but are a pain in the butt to implement
c anything more than mean or weighted mean.


      program compsigma

c max number of lines in a .concat file      
      parameter (MAXL=128)
c max number of wavelengths of lines
      parameter (MAXW=20)
c setup to restrict to only good lines?
      parameter (IFGLINESETUP=1)
cc intensity cut
c      parameter (IFINTCUT=0)

      character lname*80,fname*80,oname*80
      character objid*10
      real pars(6)
      real z(MAXL),wrest(MAXL),totint(MAXL),etotint(MAXL)
      real cont(MAXL),econt(MAXL)
      real wcen(MAXL),ewcen(MAXL),disp(MAXL),edisp(MAXL)
      real contbiw(MAXL),contbierr(MAXL),ewrest(MAXL),ewresterr(MAXL)
      real sigcorrsq(MAXL),esigcorrsq(MAXL),sigkms(MAXL),esigkms(MAXL)
      real ratioint(MAXL),eratioint(MAXL)
      integer ifuse(MAXL)
      real wgoodline(MAXW)
      real wobs(MAXL),sigsqkms(MAXL),esigsqkms(MAXL)
      real x(MAXL),ex(MAXL)
      real xout(2),xrms(2)

 100  write(*,'("file with list of .concat files: ",$)')
      read(*,'(a)') lname
      open(3,file=lname,status='old',err=100)

      open(11,file='compsigma.out',status='unknown')

      write(*,'("Take mean in linear sigma (1), squared (2): ",$)')
      read (*,*) imethod

c set up goodlist.  This is inelegant
      if (IFGLINESETUP .eq. 1) then
         ngoodline = 4
         wgoodline(1) = 3728.0
         wgoodline(2) = 4861.3
         wgoodline(3) = 5006.9
         wgoodline(4) = 6562.8
c         write(*,'("goodlines: ",ngoodline(2x,f7.1))') 
         write(*,*)
     $        (wgoodline(i),i=1,ngoodline)
         write(*,'("Use only goodlines above (1-yes): ",$)')
         read(*,*) ifgoodlines
      else
         ngood = 0
         ifgoodlines = 0
      end if

c intensity cut
      write (*,'("Signif cut on intensity, nsigma (<=0: none): ",$)')
      read (*,*) nsigintcut
      if (nsigintcut .gt. 1.e-3) then
         ifintcut = 1
c         nsigintcut = 4.0
      else
         ifintcut = 0
         nsigintcut = -1.e6
      end if

      nnolines = 0
      nobj = 0
ccc start file loop

 200  continue
      read(3,'(a)',err=666,end=666) fname
      open(10,file=fname,status='old',err=666)
      nobj=nobj+1
      
      k=0
 300  continue
c format cribbed from linefit.f should work
c sigkms is in kms
c sigcorrsq is (W_obs^2 - W_inst^2) in A
      k=k+1
      read(10,1030,err=400,end=400) objid,z(k),wrest(k),
     $     cont(k),econt(k),totint(k),etotint(k),
     $     wcen(k),ewcen(k),disp(k),edisp(k),
     $     contbiw(k),contbierr(k),ewrest(k),ewresterr(k),
     $     sigcorrsq(k),esigcorrsq(k),
     $     sigkms(k),esigkms(k),ratioint(k),eratioint(k)
 1030 format(a8,2x,f9.5,2x,f8.2,4(2x,1pe9.2),2(2x,1pe13.6),
     $        12(2x,1pe9.2))
      
      wobs(k) = wrest(k) * (1.0 + z(k))
      kmsang = 3.0e5/wobs(k)
      sigsqkms(k) = sigcorrsq(k) * kmsang**2
      esigsqkms(k) = esigcorrsq(k) * kmsang**2

      go to 300

 400  continue
      close(10)
      nlines=k-1

      do i=1,nlines
         ifuse(i) = 1
      end do

      if (ifgoodlines .eq. 1 .or. ngood .ne. 0) 
     $     call trimrestwave(nlines,wrest,ngoodline,wgoodline,ifuse)

      if (ifintcut .eq. 1)
     $     call intsigcut(nlines,totint,etotint,nsigintcut,ifuse)

c trim to use only good points, decide whether to take
c mean in sigma or sigma^2
      nfit = 0
      do i=1,nlines
         if (ifuse(i) .eq. 1) then
            nfit = nfit + 1
            if (imethod .eq. 1) then
               x(nfit) = sigkms(i)
               ex(nfit) = esigkms(i)
            else if (imethod .eq. 2) then
               x(nfit) = sigsqkms(i)
               ex(nfit) = esigsqkms(i)
            end if
         end if
      end do

c      write(*,*) 'Object ',objid,' read ',nlines,' lines, fitting ',nfit
      if (nfit .eq. 0) then
         nnolines = nnolines + 1
         write(*,*) 'Object ',objid,' read ',nlines,' lines, 0 to fit'
      end if

c two-step for mean and wtmean
c convert from sigma^2 to sigma if necessary
      do k=1,2
      if (nfit .ge. 1) then
         if (k .eq. 1) then
            call calcmean(nfit,x,ex,xout(k),xrms(k))
         else
            call calcwtmean(nfit,x,ex,xout(k),xrms(k),xchisq)
         end if
c if we only have one measurement, use its error as the error
c better to do this in the mean/rms routines?
         if (nfit .eq. 1) then
            xrms(k) = ex(1)
         end if
c convert squared to linear as necessary
         if (imethod .eq. 1) then
            continue
         else if (imethod .eq. 2) then
            if (xout(k) .gt. 0) then
               xout(k) = sqrt(xout(k))
               xrms(k) = xrms(k) * 0.5 /xout(k)
            else 
               xout(k) = -sqrt(-xout(k))
               xrms(k) = xrms(k) * 0.5 /abs(xout(k))
            end if
         end if
      else
         xout(k) = -999.
         xrms(k) = 1.e4
      end if
      end do

      write(11,1060) objid,nlines,xout(1),xrms(1),xout(2),xrms(2),xchisq
 1060 format(a10,2x,i3,5(2x,f8.2))
      go to 200

 666  continue
      write(*,*) nobj,' objects ',nnolines,' with no lines'
      close(11)
      close(3)

      end
      
cccccccccccccccccccccccccccccc
c Check to see if line is on the goodlist.  wtol allows for minor
c differences in line naming (even air-to-vac diff)
      subroutine trimrestwave(n,wave,ngood,wgood,ifuse)

      real wave(n),wgood(ngood)
      integer ifuse(n)

      wtol = 4.0
      do i= 1,n
         if (ifuse(i) .eq. 1) then
            j=1
 1100       continue
            if (abs(wave(i)-wgood(j)) .gt. wtol) then
               if (j .lt. ngood) then
                  j=j+1
                  go to 1100
               else
c got to end w/o match                  
                  ifuse(i) = 0
               end if
            end if
         end if
      end do
      return
      end
               
c check how many sigma lines was detected by
      subroutine intsigcut(n,rint,erint,nsig,ifuse)

      real rint(n),erint(n)
      integer ifuse(n)
      
      do i=1,n
         if (abs(rint(i)/erint(i)) .lt. nsig)
     $        ifuse(i) = 0
      end do
      return
      end

cccccccccccccccccccc
c  compute mean and rms
      subroutine calcmean(n,x,ex,xmean,xrms)
      real x(n),ex(n)

      sum = 0.0
      sumsq = 0.0
      do i=1,n
         sum = sum + x(i)
         sumsq = sumsq + x(i)**2
      end do
      xmean = sum/n
c      xrms = sqrt(sumsq/n - xmean**2)
c proper n-1 way?  no, that's for std error
      if (n .gt. 1) then
c         xrms = sqrt(sumsq/(n-1.) - xmean**2)
         xrms = sqrt(sumsq/n - xmean**2)
      else
         xrms = ex(1)
      end if
      return
      end

c  compute weighted mean and rms
      subroutine calcwtmean(n,x,ex,xwtmean,xwtrms,chisq)
      real x(n),ex(n)

      sumwtx = 0.0
      wtsum = 0.0
      do i=1,n
         wt = 1./ex(i)**2
         sumwtx = sumwtx + wt*x(i)
         wtsum = wtsum + wt
      end do
      xwtmean = sumwtx / wtsum
      sumwtx = 0.0
c      wtsum = 0.0
      chisq = 0.0
      do i=1,n
         wt = 1./ex(i)**2
         devsq = (x(i)-xwtmean)**2
         sumwtx = sumwtx + wt*devsq
c         wtsum = wtsum + wt
         chisq = chisq + wt*devsq
      end do
c      xwtrms = sqrt( sumwtx / wtsum)
c proper n-1 way?
      if (n .gt. 1) then
c         xwtrms = sqrt(sumsq/(n-1.) - xmean**2)
         xwtrms = sqrt( sumwtx / wtsum)
      else
         xrms = ex(1)
      end if
      chisq = chisq / n
      return
      end 

