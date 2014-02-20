c version of linefit to run on DEEP1 spectra and deal with 
c the fact that there are both B and R spectra which may
c overlap, and have different instrumental resolutions.

c given lists of object numbers, spec1d files, and spectral features,
c fit lines.  output to a file

      program linefit_deep1

      include 'specarray.h'

      parameter (MAXWFIT=3000)
      parameter (NMAXLINES=50)
      parameter(NPARAMS=6)
      parameter(IFPLOT=1)

      real wave(NWAVEMAX),spec(NWAVEMAX),espec(NWAVEMAX)
      real wdata(MAXWFIT),sdata(MAXWFIT),edata(MAXWFIT)
      real yfit(MAXWFIT)
      real wrest(NMAXLINES)
      integer ifit(NPARAMS)
      real pars(NPARAMS),epars(NPARAMS)

      character fname*160,sname*200,cline*160,objid*20
      character xlabel*60,toplabel*60,charspec*1

      fitrad = 20.
      wmaxoffset = 10.
      contrad = 100.
      ifoptextract=0
      ndphys = MAXWFIT
      npars = NPARAMS

 90   continue
      write(*,'("linelist: ",$)')
      read(*,'(a)') fname
      open(10,file=fname,status='old',err=90)
 100  continue
      write(*,'("list of object numbers and zs: ",$)')
      read(*,'(a)') fname
      open(2,file=fname,status='old',err=100)
c 110  continue
c      write(*,'("list of redshifts: ",$)')
c      read(*,'(a)') fname
c      open(3,file=fname,status='old',err=110)
 130  continue
      write(*,'("list of 1-d spectra: ",$)')
      read(*,'(a)') fname
      open(4,file=fname,status='old',err=130)

c ideally we would have a list of instrumental sigma for each object
      write(*,'("instrumental sigmas in A, for B and R: ",$)')
      read(*,*) sigib,sigir
      sigibsq = sigib**2
      sigirsq = sigir**2

c read lines from linelist

      nlines=1
 200  continue
      read(10,*,err=210,end=210) wrest(nlines)
      nlines=nlines+1
      go to 200
 210  continue
      nlines=nlines-1
      write(*,*) "Read ", nlines," wavelengths"
      close(10)

      open(11,file='linefit.out',status='unknown')

      nobj=1
 300  continue
c loop through objects
c      read(2,*,err=666,end=666) iobjno
c      read(3,*,err=666,end=666) z
c      read(2,*,err=666,end=666) iobjno,z
      read(2,'(a)',err=666,end=666) cline
      read(cline(1:8),'(a)') objid
      read(cline(9:),*) z
      read(4,'(a)',err=666,end=666) sname

      if (IFPLOT .eq. 1) then
c         write(cline,"('\"',i8,'.lfit.gif\"/gif')") iobjno
         write(cline,"('\"',a8,'.lfit.gif\"/gif')") objid
c         write(*,'(a)') cline
         call pgbeg(0,cline,3,3)
         call pgscf(1)
         call pgsch(1.8)
c make sure we get black on white
         call pgscr(0,1.0,1.0,1.0)
         call pgscr(1,0.0,0.0,0.0)
      end if

c      write(cline,'(i8,".linedat")') iobjno
      write(cline,'(a8,".linedat")') objid
      open(10,file=cline,status='unknown')

      do ispec=1,2
         if (ispec .eq. 1) then
            charspec = 'B'
            sigisq = sigibsq
         else
            charspec = 'R'
            sigisq = sigirsq
         end if
c      call read1dspec(sname,ifoptextract,np,wave,spec,espec,ifer)
      call read1dspec_deep1(sname,charspec,
     $        ifoptextract,np,wave,spec,espec,ifer)

      if (ifer .ne. 0) then
         write(*,*) "Didnt get spectrum for ",objid," ",charspec
c         go to 300
         go to 400
      end if

      do iline=1,nlines
         wfit = wrest(iline) * (1.+z)
         if (wfit .lt. wave(1) .or. wfit .gt. wave(np)) then
            ifoutofrange = 1
            nd = 0
         else
            ifoutofrange = 0
c get data in a window around the line
            call getfitdata(np,wave,spec,espec,wfit,fitrad,
     $           ndphys,nd,wdata,sdata,edata)
         end if

c         open(12,file='linefit.tmpdat',status='unknown')
c         do i=1,nd
c            write(12,*) i,wdata(i),sdata(i),edata(i)
c         end do
c         close(12)

c         write(*,*) nobj,iline,wrest(iline),wfit,nd
         if (IFPLOT .eq. 1) then
            fmin = 1.e10
            fmax=-1.e10
            do j=1,nd
               fmin = min(fmin,sdata(j))
               fmax = max(fmax,sdata(j))
            end do
c            call pgenv(wdata(1),wdata(nd),fmin,fmax,0,0)
            if (nd .ge. 1) then
               fmarg = (fmax-fmin)*0.075
               plotmin = min(fmin-marg,0.0)
               call pgenv(wfit-fitrad,wfit+fitrad,
     $              plotmin,fmax+fmarg,0,0)
               call pgpt(nd,wdata,sdata,17)
               call pgerrb(6,nd,wdata,sdata,edata,0.0)
               write(xlabel,'(f7.1," at ",f7.4)') wrest(iline),z
c            write(toplabel,'(i8,1x,a1)') iobjno,charspec
               write(toplabel,'(a8,1x,a1)') objid,charspec
c            call pglabel(xlabel," ",toplabel)
               call pgmtxt('B',2.6,0.5,0.5,xlabel)
               call pgmtxt('T',2.3,0.0,0.0,toplabel)
            end if
         end if
 
c check if we have enough data, should we test if nd=0?
         if (nd .gt. 9) then
            ifbad = 0
c need to test for 3727 and fit a constrained doublet.

            if (wrest(iline) .gt. 3725.0 .and. 
     $           wrest(iline) .lt. 3729.5) then
               nfitlines=2
               do j=1,5
                  ifit(j) = 1
                  pars(j) = 0.0
               end do
               ifit(6) = 0
               pars(6) = 3726.03 /3728.82
               iguess=2
               pars(3) = wfit
               call fitdoublet(wdata,sdata,edata,nd,pars,epars,
     $           ifit,iguess,chisq,ier)
            else
               nfitlines=1
c parameters are cont, intensity, mean, sigma
               do j=1,4
                  ifit(j)=1
                  pars(j)=0.0
               end do
               ifit(5) =0
               ifit(6)=0
               pars(5)=0.0
c let the routine make initial guesses & guess em or abs line.
c pass in the guess for the wavelength
               iguess=2
               pars(3) = wfit
c fit emission line
c               iguess=1
               call fitgaus2(wdata,sdata,edata,nd,pars,epars,
     $              ifit,iguess,chisq,ier)
            end if

c test to see if fit was bad?
c for example, if dispersion is very large,
c or cwl is out of range of the fit data, or cwl
c is too far from initial guess
            if (ier .eq. 1 .or. pars(4) .gt. 100.0 .or. 
     $           pars(3) .lt. wdata(1) .or. pars(3) .gt. wdata(nd) .or. 
     $           abs(pars(3)-wfit) .gt. wmaxoffset) then
               ifbad =1
            end if
         else
            ifbad=1
         end if

         if (ifbad .eq. 1) then
            do i=1,5
               pars(i) = -999.
               epars(i) = 1.0e3
            end do
            totint = -999.
            etotint = 1.0e3
         end if

c plot fit?
         if (IFPLOT .eq. 1 .and. nd .gt. 1) then
c            fmin = 1.e10
c            fmax=-1.e10
c            do j=1,nd
c               fmin = min(fmin,sdata(j))
c               fmax = max(fmax,sdata(j))
c            end do
cc            call pgenv(wdata(1),wdata(nd),fmin,fmax,0,0)
c            call pgenv(wfit-fitrad,wfit+fitrad,fmin,fmax,0,0)
c            if (nd .ge. 1) then
c               call pgpt(nd,wdata,sdata,17)
c               call pgerrb(6,nd,wdata,sdata,edata,0.0)
c            end if
            if (ifbad .eq. 0) then
               if (nfitlines .eq. 1) then
                  call evalgaus(nd,wdata,npars,pars,yfit)
               else
                  call evaldoublet(nd,wdata,npars,pars,yfit)
               end if
               call pgsci(4)
               call pgline(nd,wdata,yfit)
               call pgsci(1)
            end if
         end if
            
c get data in the continuum window excluding the line
         call getcontdata(np,wave,spec,espec,wfit,contrad,fitrad,
     $        ndphys,nd,wdata,sdata,edata)
         if (nd .ge. 1) then
            call biwgt(sdata,nd,contbiw,contbisc)
            contbierr = contbisc / sqrt(real(nd))
         else
            contbiw = 0.
            contbierr = 1.e3
         end if

c output fit and continuum result and EW, check ifoutofrange and ifbad
c to output some kind of flag?

c note there is a difference between observed and restframe EW
c a 1+z factor.

         if (ifbad .eq. 0) then
            if (nfitlines .eq. 1) then
               totint = pars(2)
               etotint = epars(2)
            else
c for doublet, add intensities
               totint = pars(2) + pars(5)
c               etotint = epars(2)+epars(5)
               etotint = sqrt(epars(2)**2 + epars(5)**2)
            end if
            ew = totint / contbiw
            ewrest = ew / (1.+z)
            ewerr = sqrt((etotint/contbiw)**2 + 
     $           (contbierr*totint/contbiw**2)**2)
c     $           (totint/contbierr**2)**2)
            ewresterr = ewerr / (1.+z)
            sigcorrsq = pars(4)**2 - sigisq
            esigcorrsq = 2.*pars(4)*epars(4)
            if (sigcorrsq .gt. 0.) then
               sigkms = sqrt(sigcorrsq) *3.0e5/pars(3)
               esigkms = 0.5/sqrt(sigcorrsq)*esigcorrsq *3.0e5/pars(3)
            else
c this happens if the measured sigma is less than instrumental
c so the vel disp is formally imaginary
               sigkms = -sqrt(-sigcorrsq) *3.0e5/pars(3)
c this is clearly wrong for the error estimate.  what to do?
c the error is formally infinite.
c               esigkms = 0.5/sqrt(-sigcorrsq)*esigcorrsq *3.0e5/pars(3)
               esigkms = 999.0
            end if
            if (IFPLOT .eq. 1) then
               write(toplabel,'("EW_0 ",f7.1)') ewrest
               call pgmtxt('T',1.1,0.0,0.0,toplabel)
               write(toplabel,'("Disp  ",f5.2," +- ",f5.2," A")') 
     $              pars(4),epars(4)
               call pgmtxt('T',2.3,0.4,0.0,toplabel)
               write(toplabel,'("Sigma ",f5.1," +- ",f5.1)') 
     $              sigkms,esigkms
               call pgmtxt('T',1.1,0.4,0.0,toplabel)
            end if               
         else
            ew = -999.
            ewerr = 1.0e4
            ewrest = -999.
            ewresterr = 1.0e4
            sigcorrsq = -999.
            sigkms = -999.
            esigcorrsq = 1.0e4
            esigkms = 1.0e4
c            if (IFPLOT .eq. 1 .and. nd  .ge. 1) then
c               write(toplabel,'(a)') "EW = ***"
c               call pgmtxt('T',1.1,0.0,0.0,toplabel)
c            end if
         end if
         
c         write (10,1020) iobjno,z,wrest(iline),(pars(i),epars(i),i=1,4),
c     $    contbiw,contbierr,ew,ewerr,sigcorrsq,esigcorrsq,sigkms,esigkms
c         write (11,1020) iobjno,z,wrest(iline),(pars(i),epars(i),i=1,4),
c     $    contbiw,contbierr,ew,ewerr,sigcorrsq,esigcorrsq,sigkms,esigkms
c 1020    format(i8,2x,f9.5,2x,f8.2,4(2x,1pe9.2),2(2x,1pe13.6),
c     $        10(2x,1pe9.2))

c now writing out EW restframe, not EW observed
         write (10,1030) objid,z,wrest(iline),pars(1),epars(1),
     $        totint,etotint,(pars(i),epars(i),i=3,4),contbiw,contbierr,
     $        ewrest,ewresterr,sigcorrsq,esigcorrsq,sigkms,esigkms,
     $        charspec
         write (11,1030) objid,z,wrest(iline),pars(1),epars(1),
     $        totint,etotint,(pars(i),epars(i),i=3,4),contbiw,contbierr,
     $        ewrest,ewresterr,sigcorrsq,esigcorrsq,sigkms,esigkms,
     $        charspec
 1030    format(a8,2x,f9.5,2x,f8.2,4(2x,1pe9.2),2(2x,1pe13.6),
     $        10(2x,1pe9.2),2x,a1)
      end do

 400  continue
c end the ispec loop
      end do

      close(10)
      if (IFPLOT .eq. 1) then
         call pgend()
      end if

      nobj=nobj+1
      go to 300

 666  continue
      close(11)
      close(2)
c      close(3)
      close(4)
      nobj=nobj-1
      write(*,*) "fit lines for ",nobj," spectra"

c      call pgend()

      end

cccccccccccccccccccccccccccccc
c return data within xrad of xcen
c assume input array is sorted ascending
      subroutine getfitdata(n,x,y,yerr,xcen,xrad,nftot,
     $     nfit,xfit,yfit,eyfit)

      real x(n),y(n),yerr(n)
      real xfit(nftot),yfit(nftot),eyfit(nftot)

      include 'pcredshift.h'

      xmin = xcen-xrad
      xmax = xcen+xrad

      call locate(x,n,xmin,jmin)
      call locate(x,n,xmax,jmax)
c      nfit = jmax - jmin
      nfit = 0
      do j=jmin+1,jmax
c         i=j-jmin
         if(y(j) .gt. bad .and. yerr(j) .gt. 1.0e-2) then
            nfit = nfit+1
            xfit(nfit) = x(j)
            yfit(nfit) = y(j)
            eyfit(nfit) = yerr(j)
         end if
      end do
      return
      end
         
cccccccccccccccccccccccccccccc
c return data within xrad of xcen excluding a window xexclude in radius
c assume input array is sorted ascending
      subroutine getcontdata(n,x,y,yerr,xcen,xrad,xexclude,nftot,
     $     nfit,xfit,yfit,eyfit)

      real x(n),y(n),yerr(n)
      real xfit(nftot),yfit(nftot),eyfit(nftot)

      include 'pcredshift.h'

      xmin = xcen-xrad
      xmax = xcen+xrad
      exclmin = xcen - xexclude
      exclmax = xcen + xexclude

      call locate(x,n,xmin,jmin)
      call locate(x,n,xmax,jmax)
      call locate(x,n,exclmin,jexmin)
      call locate(x,n,exclmax,jexmax)

      nfit = 0
      do j=jmin+1,jexmin
         if (y(j) .gt. bad .and. yerr(j) .gt. 0.01) then
            nfit = nfit+1
            xfit(nfit) = x(j)
            yfit(nfit) = y(j)
            eyfit(nfit) = yerr(j)
         end if
      end do
      do j=jexmax+1,jmax
         if (y(j) .gt. bad .and. yerr(j) .gt. 0.01) then
            nfit = nfit+1
            xfit(nfit) = x(j)
            yfit(nfit) = y(j)
            eyfit(nfit) = yerr(j)
         end if
      end do
      return
      end
         
cccccccccccccccccccc
c from Numerical Recipes
      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END

cccccccccc
c evaluate gaussian on points x(i) from parameters pars(npars)
      subroutine evalgaus(n,x,npars,pars,y)
      
      real x(n),y(n)
      real pars(npars)

      root2pi = sqrt(2.0*3.1415927)
      r2pisig = root2pi*abs(pars(4))

      do i=1,n
         y(i) = pars(1) + pars(2)/r2pisig * 
     $        exp( -(x(i)-pars(3))**2 / 2. / pars(4)**2)
      end do
      return
      end

cccccccccc
c evaluate doublet on points x(i) from parameters pars(npars)
      subroutine evaldoublet(n,x,npars,pars,y)
      
      real x(n),y(n)
      real pars(npars)

      root2pi = sqrt(2.0*3.1415927)
      r2pisig = root2pi*abs(pars(4))

      do i=1,n
         y(i) = pars(1) + pars(2)/r2pisig * 
     $        exp( -(x(i)-pars(3))**2 / 2. / pars(4)**2) +
     $        pars(5)/r2pisig * 
     $        exp(-(x(i)-pars(3)*pars(6))**2 /2./pars(4)**2)
      end do
      return
      end

