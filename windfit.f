

c this program is to fit the Mg II absorption doublet
c with a two-component model

      program windfit

c      include 'specarray.h'

      parameter (NWMAX=30000)
      parameter (MAXWFIT=4000)
      parameter(NPARAMS=4)
      parameter(NEMPARAMS=6)

c if we are doing the Mg I singlet, otherwise assume Mg II
      parameter (IFSINGLE=0)
c if we are using the "fitgaus2special" routine that puts the
c variables in a different order
      parameter(IFITSPECIAL=1)
c also fit a gaussian to the 3727 emission doublet?
      parameter(IFITEM=1)

      real wave2(NWMAX),spec2(NWMAX),espec2(NWMAX)
      real wave(MAXWFIT),spec(MAXWFIT),espec(MAXWFIT)
      real tmp(MAXWFIT),contin(MAXWFIT)
c      real wdata(MAXWFIT)
      real sdata(MAXWFIT),edata(MAXWFIT),vdata(MAXWFIT),snorm(MAXWFIT)
      real sdat2(MAXWFIT),edat2(MAXWFIT),vdat2(MAXWFIT),snorm2(MAXWFIT)
      real vtrim(MAXWFIT),ctrim(MAXWFIT)
      real xpl(2),ypl(2)

      real scomp1(MAXWFIT),scomp2(MAXWFIT),scomptot(MAXWFIT)
      real bcomp(MAXWFIT),scompplot(MAXWFIT),bcompplot(MAXWFIT)
      real emplot(MAXWFIT)
      real cumulopac(MAXWFIT),ecumulopac(MAXWFIT)
      real opacthresh(10),vthresh(10),othresh(10),errcumul(10)
      real evthresh1(10),evthresh2(10),evthresh(10)
      
      integer ifit(NPARAMS),iemfit(NEMPARAMS)
      real pars(NPARAMS),epars(NPARAMS)
      real empars(NEMPARAMS),eempars(NEMPARAMS)
      character fname*160,sname*160
      character xlabel*60,toplabel*60

      call pgbegin (0,"?",1,1)
      call pgsch(1.4)
      call pgslw(2)
      call pgscf(2)
      call pgask (.true.)

      root2pi = sqrt(2.0*3.1416)
      clight = 3.0e5
c air
c      wline1 = 2795.528
c      wline2 = 2802.705
c vacuum
      if (IFSINGLE .eq. 1) then
         wline1 = 2840.
         wline2 = 2852.964
      else
         wline1 = 2796.352
         wline2 = 2803.531
      end if
c velocity of the bluer line in the restframe of the redder one
      vline1 = wavetovel(wline2,wline1)
      write(*,*) vline1
c intensity ratio line1/line2 for MgII
      ratioint12 = 1.0

c spectrograph dispersion in A
      siginsta = 0.56
c typical z, need this to convert siginsta
      ztypical = 1.4

c continuum windows in velocity relative to line 2
      cv1min = -3000.
      cv1max = -2000.
      cv2min = 800.
      cv2max = 1800.
      cw1min = veltowave(wline2,cv1min)
      cw1max = veltowave(wline2,cv1max)
      cw2min = veltowave(wline2,cv2min)
      cw2max = veltowave(wline2,cv2max)


c is this going to work?
 100  continue
      write(*,'("name of spec1d file: ",$)')
      read(*,'(a)') sname
c read spectrum.  This would work on an original DEEP2 spec1d file,
c but not on the coadded spectrum data structure, which is simpler.
c      ifoptextract=0
c      call read1dspec(sname,ifoptextract,np2,wave2,spec2,espec2,ifer)
c readcoadd should just read the SPEC,LAMBDA,IVAR,TOTN from a 
c idl-written bintable
      call readcoadd(sname,np2,wave2,spec2,espec2,ifer)
      if (ifer .ne. 0) then
         write(*,*) "error reading file"
         go to 100
      end if

c if spectrum not restframe
c      do i=1,np2
c         wave2(i) = wave2(i) / (1.0 + z)
c      end do

c fit emission?
      if (IFITEM .eq. 1) then
c convert to vacuum
         wemline = 3728.82 * 1.00029
         ifdoublet = 1
         wtrim1 = wemline - 20.
         wtrim2 = wemline + 20.
         call locate(wave2,np2,wtrim1,ntrim1)
         call locate(wave2,np2,wtrim2,ntrim2)
         npemfit = ntrim2-ntrim1+1
c         write(*,*) "em: ",npemfit,wtrim1,wtrim2,ntrim1,ntrim2
         do i=1,npemfit
            iold = i+ntrim1-1
            wave(i) = wave2(iold)
            spec(i) = spec2(iold)
            espec(i) = espec2(iold)
         end do
         do i=1,5
            iemfit(i) = 1
         end do
         iemfit(6) = 0
         empars(6) = 3726.03 / 3728.82 
         empars(3) = wemline
         call fitdoublet(wave,spec,espec,
     $        npemfit,empars,eempars,iemfit,1,emchisq,ier)
c         write(*,*) "emfit : ",(empars(j),j=1,6)
c         write(*,*) "emerr : ",(eempars(j),j=1,6)
c         write(*,*) "chisq, ier: ",chisq,ier
c for doublet, add intensities
         emchisq = emchisq / npemfit
         emcont = empars(1)
         eminty = empars(2) + empars(5)
         eeminty = sqrt(eempars(2)**2 + eempars(5)**2)
c make a robust continuum est or just use the fit continuum?
         emew = eminty / empars(1)
         eemew = sqrt((eeminty/empars(1))**2 + 
     $        (eminty*eempars(1)/empars(1)**2)**2)
c we're in rest w.l. so don't need to do a 1+z correction.
c will need to correct for siginst though.
         emvdisp = empars(4) * 3.0e5 / wemline
         eemvdisp = eempars(4) * 3.0e5 / wemline
         sigsq = empars(4)**2
         do i=1,npemfit
            x = wave(i)
            emplot(i) = empars(1) + empars(2)/root2pi/empars(4) *
     $           exp( -(x-empars(3))**2/2.0/sigsq ) +
     $           empars(5)/root2pi/empars(4) * 
     $           exp( -(x-empars(3)*empars(6))**2/2.0/sigsq )
         end do
         yplmax = empars(1) + 1.2 * empars(2)/root2pi/empars(4) 
         call pgenv(wave(1),wave(npemfit),0.0,yplmax,0,0)
         call pglabel("wavelength","DN/pixel","emission line fit")
c         call pgbin(npemfit,wave,spec,.true.)
         call pgpoint(npemfit,wave,spec,17)
         call pgerrb(6,npemfit,wave,spec,espec,0.0)
c         call pgqci(indexc)
c         call pgsci(2)
         call pgline(npemfit,wave,emplot)
c         call pgsci(indexc)
      end if

c trim the spectrum to some region around the lines
      wtrim1 = wline2 - 60.
      wtrim2 = wline2 + 40.

      call locate(wave2,np2,wtrim1,ntrim1)
      call locate(wave2,np2,wtrim2,ntrim2)

      np = ntrim2-ntrim1+1
      if (np .gt. MAXWFIT) then 
         write(*,*) "too many points in fit trimming region ",np,MAXWFIT
      end if
c      write(*,*) np," points in trim region ",ntrim1,ntrim2
      do i=ntrim1,ntrim2
         inew = i-ntrim1+1
         wave(inew) = wave2(i)
         spec(inew) = spec2(i)
         espec(inew) = espec2(i)
      end do
         
      call setcontvec(np,wave,spec,espec,cw1min,cw1max,cw2min,cw2max,
     $     contin)

c      do i=ntrim1,ntrim2
c         inew = i-ntrim1+1
c         vtrim(inew) = wavetovel(wline2,wave2(i))
c         ctrim(inew) = contin(i)
c      end do
      do i=1,np
         vtrim(i) = wavetovel(wline2,wave(i))
         ctrim(i) = contin(i)
      end do

c      write(*,'("Intensity ratio for I_2796/I_2803: ",$)')
c      read(*,*) ratioint12

c get the data for the >0 vel side of the redder line

      vfitmin = 0.
      vfitmax = 2000.
      wfitmin = veltowave(wline2,vfitmin)
      wfitmax = veltowave(wline2,vfitmax)
      call locate(wave,np,wfitmin,ifitmin)
      ifitmin = ifitmin+1
      call locate(wave,np,wfitmax,ifitmax)
      call locate(wave,np,wline2,iline2center)
c      write(*,*) ifitmin,ifitmax,iline2center

      npfit = ifitmax-ifitmin+1
      do i=ifitmin,ifitmax
c check for bad data here?
         inew = i-ifitmin+1
c         wdata(inew) = wave(i)
         vdata(inew) = wavetovel(wline2,wave(i))
c         cdata(inew) = contin(i)
c         sdata(inew) = spec(i)
         snorm(inew) = spec(i) / contin(i)
         edata(inew) = espec(i) / contin(i)
      end do
c symmetrize the data before fitting
      do i=1,npfit
         jpos = npfit+i
         jneg = npfit+1-i
         vdat2(jpos) = vdata(i)
         vdat2(jneg) = -vdata(i)
         snorm2(jpos) = snorm(i)
         snorm2(jneg) = snorm(i)
         edat2(jpos) = edata(i)
         edat2(jneg) = edata(i)
      end do
      npfit2 = 2*npfit
c      write(*,*) "npfit2 = ",npfit2

c      open(4,file='test1.out',status='unknown')
c      do i=1,npfit2
c         write(4,*) i,vdat2(i),snorm2(i),edat2(i)
c      end do
c      close(4)

c fit a gaussian to it
c hold velocity fixed at 0 (or let it float if data symmetrized)
c maybe i should force the continuum to be a value obtained from
c setcontvec.  see if we can get it to fit right with fixed
c vel and iguess=0, passing initial guesses  in with a negative
c intensity
      if (IFITSPECIAL .eq. 0) then
         ifit(1) = 0
         ifit(2) = 1
         ifit(3) = 0
         ifit(4) = 1
         pars(1) = 1.0
         pars(2) = -100.
         pars(3) = 0.0
         pars(4) = 50.
      else
         ifit(3) = 0
         ifit(1) = 1
         ifit(4) = 0
         ifit(2) = 1
         pars(3) = 1.0
         pars(1) = -100.
         pars(4) = 0.0
         pars(2) = 50.
      end if
c      iguess=-1
      iguess = 0
      if (IFITSPECIAL .eq. 0) then
c         call fitgaus2(vdat2,snorm2,edat2,npfit2,pars,epars,
c     $        ifit,iguess,chisq,ier)
      else
         call fitgaus2special(vdat2,snorm2,edat2,npfit2,pars,epars,
     $     ifit,iguess,chisq,ier)
      end if
      contcenter = contin(iline2center)
      write(*,*) "contin value at line center: ",contcenter
      write(*,*) "fit pars:  ",(pars(i),i=1,4)
      write(*,*) "fit errs:  ",(epars(i),i=1,4)
      write(*,*) "n*2,chisq: ",npfit2,chisq
c convert intensity of opacity from km/s to EW in A
      symfitew = -pars(1) * wline2/3.0e5
      esymfitew = epars(1) * wline2/3.0e5
c we should be correcting the vdisp for siginst
      vdispsym = pars(2)
      evdispsym = epars(2)
      write(*,'("EW,err: ",f6.2,1x,f6.2,"  vdisp,err: ",f7.2,1x,f6.2)')
     $     symfitew,esymfitew,vdispsym,evdispsym

c compute the gaussian on the velocity scale vdata
c     duplicate the gaussian for the bluer line and multiply to
c get symmetric component 

c we are adding the fitted continuum of the opacity, which
c is often off 1.0 by 1-2%.  this is annoying.  I could force
c the continuum to be 1.0
      if (IFITSPECIAL .eq. 0) then
         parcont = pars(1)
         parinty = pars(2)
         parvel  = pars(3)
         parsig  = pars(4)
      else
         parcont = pars(3)
         parinty = pars(1)
         parvel  = pars(4)
         parsig  = pars(2)
      end if
      fitsigsq = parsig**2
      do i=1,np
         v = vtrim(i)
         scomp2(i) = parcont + parinty/root2pi/parsig * 
     $           exp(-(v-parvel)**2/2.0/fitsigsq)
c            scomp2(i) = pars(3) + pars(1)/root2pi/pars(2) * 
c     $           exp(-(v-pars(4))**2/2.0/fitsigsq)
c compute velocity relative to vline1 where the blue line is predicted to be
c in principle, should also transform the sigma, but the change in sigma
c is really tiny.  we are assuming that the ratio of the two line intensities
c ratioint12 is 1.0, which may be wrong.
         if (IFSINGLE .eq. 1) then
            scomp1(i) = 1.0
         else
          scomp1(i) = parcont + ratioint12*parinty/root2pi/parsig * 
     $        exp(-(v-(vline1+parvel))**2/2.0/fitsigsq)
c         scomp1(i) = pars(3) + ratioint12*pars(1)/root2pi/pars(2) * 
c     $        exp(-(v-vline1)**2/2.0/fitsigsq)
         end if
         scomptot(i) = scomp1(i) * scomp2(i)
c divide to get the residual blueshifted component
         bcomp(i) = spec(i)/ctrim(i) / scomptot(i)
c multiply back by the continuum for plotting
         scompplot(i) = scomptot(i) * ctrim(i)
         bcompplot(i) = bcomp(i) * ctrim(i)
      end do

c Plot the fit to the symmetric component

      call pgenv(vdat2(1),vdat2(npfit2),0.,2.,0,0)
      write (xlabel,'("velocity relative to ",f6.1)') wline2
      call pglabel(xlabel,"DN/pixel/continuum",
     $     "symmetrized component fit")
      call pgpoint(npfit2,vdat2,snorm2,17)
c      call pgerry(npfit2,vdat2,snorm2,edat2)
      call pgerrb(6,npfit2,vdat2,snorm2,edat2,0.0)
      call pgline(np,vtrim,scomp2)
c      write (*,*) "got here 1"
         
c make plot
      xbd1 = -2500.
      xbd2 = 1000.
      ybd1 = 0.0
      ybd2 = 1.5*contcenter
      call pgenv(xbd1,xbd2,ybd1,ybd2,0,0)
c      write (xlabel,'(a,f6.1)') "velocity relative to ",wline2
      write (toplabel,'(a)') sname
      call pglabel(xlabel,"DN/pixel",toplabel)
c vertical lines for zero velocity
      xpl(1) = 0.0
      xpl(2) = 0.0
      ypl(1) = ybd1
      ypl(2) = ybd2
      call pgline(2,xpl,ypl)
      if (IFSINGLE .eq. 0) then
         xpl(1) = vline1
         xpl(2) = vline1
         call pgline(2,xpl,ypl)
      end if
c plot the spectra as histograms 
      call pgqci(icolor)
      call pgbin(np,vtrim,spec,.true.)
      call pgbin(np,vtrim,espec,.true.)
      call pgsci(3)
      call pgbin(np,vtrim,ctrim,.true.)
      call pgsci(2)
      call pgbin(np,vtrim,scompplot,.true.)
      call pgsci(4)
      call pgbin(np,vtrim,bcompplot,.true.)
      call pgsci(icolor)

      
c print out components
c the fields are wavelength, velocity, obsflux, error, continuum,
c sym_component/continuum, wind_component/continuum (i.e. normalized to 1)
c      write(*,*) np
      open(3,file='windfit.out',status='unknown')
      do i=1,np
         write(3,1100) wave(i),vtrim(i),spec(i),espec(i),
     $        ctrim(i),scomptot(i),bcomp(i)
      end do
 1100 format(f7.2,2x,f12.2,3(2x,f8.3),2(2x,f8.4))
      close(3)

c try to cumulate the opacity in the bluer line
c now try cumulating up from the min velocity, so
c the error at large negative velocity is more meaningful.
      if (IFSINGLE .eq. 0) then
         vopacmin = -1200.0
      else
         vopacmin = -800.0
      end if
      vopacmax = 0.0
      opacthresh(1) = 0.5
      opacthresh(2) = 0.75
      opacthresh(3) = 0.85
      nthresh = 3
c bcomp can be >1, especially if there is "P Cygni" emission,
c leading to negative opacity
c Don't count >1 pixels above some velocity to exclude the Pcyg.
c can count the ones at large negative velocity that are due to noise.
c Set vopacignore to a large negative number to also ignore
c >1 pixels at large neg vel.
      vopacignore = -300.0
      if (IFSINGLE .eq. 0) then
         voffset = vline1
      else 
         voffset = 0.0
      endif
      vopacmin = vopacmin + voffset
      vopacmax = vopacmax + voffset
      vopacignore = vopacignore + voffset
      call locate(vtrim,np,vopacmin,imin)
      call locate(vtrim,np,vopacmax,imax)
c      cumulopac(imax+1) = 0.0
c      ecumulopac(imax+1) = 0.0
      cumulopac(imin) = 0.0
      ecumulopac(imin) = 0.0
      do i = imin+1,imax
         opac = 1.0-bcomp(i)
c bcomp is = spec/continuum/scomptot, and is 1 when no opacity
c         erropac = espec(i)/spec(i) * bcomp(i)
         erropac = espec(i)/ctrim(i)/scomptot(i)
c positive pixels at e.g. -300<v<0 might be excess emission
         if (vtrim(i) .gt. vopacignore .and. opac .lt. 0.0) then
            opac = 0.0
c            erropac = espec(i) / spec(i)
            erropac = espec(i) / ctrim(i)
         end if
c         cumulopac(i) = cumulopac(i+1) + opac
c         ecumulopac(i) = sqrt(ecumulopac(i+1)**2 + erropac**2)
         cumulopac(i) = cumulopac(i-1) + opac
         ecumulopac(i) = sqrt(ecumulopac(i-1)**2 + erropac**2)
      end do
c multiply by pixel width in km/s
      dv = (vtrim(imax)-vtrim(imin)) / real(imax-imin)
c      totopac = cumulopac(imin) * dv
c      etotopac = ecumulopac(imin) * dv
      totopac = cumulopac(imax) * dv
      etotopac = ecumulopac(imax) * dv
      if (IFSINGLE .eq. 1) then
         wlinewind = wline2
      else
c we used the bluer line        
         wlinewind = wline1
      end if
c convert to EW in A by lambda/c
      totew = totopac * wlinewind / 3.0e5
      etotew = etotopac * wlinewind / 3.0e5

c find the vel thresholds at which the opacity reached X fraction
c of total, or 1-X if we're counting up from the min value
      npcumul = imax-imin+1
      do j=1,nthresh
c find vthresh(j). normalize 1-opacthresh to total opacity
         othresh(j) = (1.-opacthresh(j)) * totopac/dv
         call myinterp(npcumul,cumulopac(imin),vtrim(imin),
     $        othresh(j),ithresh,vthresh(j))
c having found the threshold, use error on the cumul. opacity
c at the threshold to estimate the error on velocity.
c first, get the error on opac at the location of the threshold
         call myinterp(npcumul,cumulopac(imin),ecumulopac(imin),
     $        othresh(j),iout,errcumul(j))
         eothresh1 = othresh(j) - errcumul(j)
         eothresh2 = othresh(j) + errcumul(j)
c now find the velocities where it crosses the error thresholds.
c Basically we are multiplying err(cumulopac) by dv/d(cumulopac).
c However these values can be bogus since the cumulative function
c does skip around.
c A problem was excessive error ranges caused by skipping.
c For ex suppose we are looking for the place where it hits 0.15,
c and the formal error is 0.06 there.  then the error range is
c given by where it hits 0.09 and 0.21.  however, it might hit
c 0.09 more than once, e.g. as velocity increases toward 0, it
c could go up to 0.11, back down, then through 0.09 and up to 0.15.
c taking the first occurrence of 0.09 gives an unreasonably large vel range.
c         call myinterp(npcumul,cumulopac(imin),vtrim(imin),
c     $        eothresh1,iout,evthresh1(j))
c         call myinterp(npcumul,cumulopac(imin),vtrim(imin),
c     $        eothresh2,iout,evthresh2(j))
         call finderrrange(npcumul,vtrim(imin),cumulopac(imin),
     $      vthresh(j),eothresh1,eothresh2,evthresh1(j),evthresh2(j))
         evthresh(j) = (evthresh2(j) - evthresh1(j))/2.0
c transform into -+ errors?
         evthresh1(j) = vthresh(j) - evthresh1(j)
         evthresh2(j) = evthresh2(j) - vthresh(j)
c for debugging
c         write(*,1170) "op.frac ",opacthresh(j),othresh(j),errcumul(j),
c     $        vthresh(j),evthresh1(j),evthresh2(j)
c 1170    format(a,2x,f5.2,2(2x,f6.2),3(2x,f7.1))
      end do

c plot velocities
      ypl(1) = contcenter*0.8
      ypl(2) = contcenter*1.1
      xpl(1) = vopacmin
      xpl(2) = vopacmin
      call pgline(2,xpl,ypl)
      xpl(1) = vopacmax
      xpl(2) = vopacmax
      call pgline(2,xpl,ypl)
      write(*,'("wind absorption and err in km/s: ",f8.2,2x,f8.2)') 
     $     totopac,etotopac
      write(*,'("wind absorption and err in EW,A: ",f8.4,2x,f8.4)')
     $     totew,etotew
      do j=1,nthresh
         xpl(1) = vthresh(j)
         xpl(2) = vthresh(j)
         call pgline(2,xpl,ypl)
c if doublet, put this vel relative to bluer line
         vthresh(j) = vthresh(j) - voffset
         write(*,1190) opacthresh(j),vthresh(j),
     $        evthresh1(j),evthresh2(j)
 1190    format("Velocity containing ",f6.3," of the opacity: ",
     $        f8.2," - ",f7.2," + ",f7.2)
      end do
      open(4,file='windfit.cumulopac',status='unknown')
      do i=imin,imax
         write(4,'(f8.2,2x,f8.2,2x,f8.3,2x,f7.4,2(2x,f7.4))') vtrim(i),
     $        vtrim(i)-voffset,cumulopac(i)*dv,ecumulopac(i)*dv,
     $        cumulopac(i)*dv/totopac,ecumulopac(i)*dv/totopac
      end do
      close(4)

c write some data out neatly for pasting into data files
      write(*,1200) wline2,contcenter,pars(1),epars(1),
     $     symfitew,esymfitew,vdispsym,evdispsym
 1200 format("sym    ",f6.1,1x,f4.1,3x,f6.1,1x,f5.1,3x,f5.2,1x,f5.2,
     $     3x,f6.1,1x,f6.1)
      write(*,1210) wlinewind,contcenter,totopac,etotopac,totew,etotew
 1210 format("wind   ",f6.1,1x,f4.1,3x,f6.1,1x,f5.1,3x,f5.2,1x,f5.2,
     $     3x,$)
      do j=1,nthresh
         write(*,1220) vthresh(j),evthresh1(j),evthresh2(j)
 1220    format(1x,f7.1,2(1x,f5.1),$)
      end do
      write(*,*)
      if (IFITEM .eq. 1) 
     $     write(*,1230) wemline,emcont,eminty,eeminty,emew,eemew,
     $     emvdisp,eemvdisp,emchisq
 1230 format("emfit  ",f6.1,2x,f4.1,2x,f7.1,1x,f6.1,3x,f6.2,1x,f6.2,
     $     3x,f6.1,1x,f6.1,3x,f6.3)


      end

cccccccccc
      function veltowave(wrest,dv)
      real wrest,dv
      clight = 3.0e5
      veltowave = wrest * (1.0 + dv/clight)
      return
      end

cccccccccc
      function wavetovel(wrest,wobs)
      real wrest,wobs
      clight = 3.0e5
      wavetovel = clight * (wobs/wrest - 1.0)
      return
      end

cccccccccc
c read a few fields out of an IDL bintable
c this follows somewhat read1dspec.f but is simplified
      subroutine readcoadd(fname,n,wave,spec,espec,ifer)
      parameter (MAXNWAVE=30000)
      character fname*160
      real wave(MAXNWAVE),spec(MAXNWAVE),espec(MAXNWAVE)
c      integer ntot(MAXNWAVE)

      ierr=0
      ifer=0
      call ftgiou(im,ierr)
      call doerr(ierr)
      call ftopen(im,fname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         ifer = 1
         write(*,'("Error opening ",a)') fname
         go to 666
      end if

      indhdu = 2
      call readfield(im,indhdu,'SPEC',itype1,isize1,naxis,laxes,spec)
      call readfield(im,indhdu,'LAMBDA',itype2,isize2,naxis,laxes,wave)
      call readfield(im,indhdu,'IVAR',itype3,isize3,naxis,laxes,espec)
c      call readfield(im,indhdu,'TOTN',itype4,isize4,naxis,laxes,ntot)
      call ftclos(im,ierr)
      call doerr(ierr)
      call ftfiou(im,ierr)
      call doerr(ierr)

      n = isize1
c      write(*,*) "read: ",n,wave(1),wave(n)
      do i=1,n
         if (espec(i) .gt. 1.0e-10) then
            espec(i) = 1.0/sqrt(espec(i))
         else
            spec(i) = -800.0
            espec(i) = 1.e4
         end if
      end do
 666  continue
      return
      end     


cccccccccc
c compute continuum
      subroutine setcontvec(n,wave,spec,espec,cw1a,cw1b,cw2a,cw2b,cvec)
      real wave(n),spec(n),espec(n),cvec(n)
      real stmp(200)
      
      bad = -600.0
      call locate(wave,n,cw1a,n1a)
      call locate(wave,n,cw1b,n1b)
      call locate(wave,n,cw2a,n2a)
      call locate(wave,n,cw2b,n2b)
c do a straight sum, not weighted
      s1=0.0
      ngood1 = 0.0
      wavetot1 = 0.0
      do i=n1a,n1b
         if (spec(i) .gt. bad) then
            ngood1 = ngood1+1
            wavetot1 = wavetot1 + wave(i)
            s1 = s1 + spec(i)
         end if
      end do
      if (ngood1 .le. 0) 
     $     write(*,*) "error: no good points for continuum window 1"
      avgs1 = s1 / ngood1
      avgw1 = wavetot1 / ngood1
c same for second window
      s2=0.0
      ngood2 = 0.0
      wavetot2 = 0.0
      do i=n2a,n2b
         if (spec(i) .gt. bad) then
            ngood2 = ngood2+1
            wavetot2 = wavetot2 + wave(i)
            s2 = s2 + spec(i)
         end if
      end do
      if (ngood2 .le. 0) 
     $     write(*,*) "error: no good points for continuum window 2"
      avgs2 = s2 / ngood2
      avgw2 = wavetot2 / ngood2
c make a linear continuum vector
      do i=1,n
         cvec(i) = avgs1 + 
     $        (wave(i) - avgw1) * (avgs2-avgs1)/(avgw2-avgw1)
      end do
      return
      end

c find e.g. the velocity y at which the opacity x goes
c above some threshold - like interp but don't use locate
c x is assumed increasing or mostly so
      subroutine myinterp(n,x,y,xthresh,iout,yout)
      real x(n),y(n)

      i=1
 2010 continue
      if (x(i) .lt. xthresh) then
         if (i .lt. n) then
            i = i+1
            go to 2010
         else
            iout = n
            yout = y(n)
            return
         end if
      else
c found where it crosses the threshold
         if (i .eq. 1) then
            iout = 1
            yout = y(1) 
            return
         else
            frac = (xthresh-x(i-1)) / (x(i)-x(i-1))
            iout = i-1
            yout = (1.0-frac) * y(i) + frac * y(i-1)
         end if
      end if
      return
      end

cccccccccccccccccccc
c given a y center value and yr1,yr2 which are some low and high error
c ranges about it, e.g. 0.15 and 0.09, 0.21, find the values of xr1 and
c xr2 that correspond.  we are assuming that the x array is increasing
c and that y is generally increasing, but always - for ex it 
c might go up to 0.11, drop back, then up through 0.09 to 0.15,
c and we want the occurence that is closer to 0.15.

      subroutine finderrrange(n,x,y,xcen,yr1,yr2,xr1,xr2)
      real x(n),y(n)

c find the y value at xcen; not really necessary to find ycen,
c but need the index icen where xcen occurs - xcen is between
c icen and icen+1
      call myinterp(n,x,y,xcen,icen,ycen)
      if (yr1 .gt. ycen .or. yr2 .lt. ycen) then
         write(*,*) "finderrrange: warning: y ranges out of order ",
     $        yr1,ycen,yr2
      end if
c count down from icen+1 to find the lower error bound
      j=icen+1
 2200 continue
      if (y(j) .gt. yr1) then
         j=j-1
         if (j .ge. 1) then
            go to 2200
         else
c hit bottom of array
            xr1 = x(1)
         end if
      else
c now y(j) is below yr1
         frac = (yr1-y(j)) / (y(j+1)-y(j))
         xr1 = frac * x(j+1) + (1.0-frac) * x(j)
      end if

c count up from icen to find the upper error bound
      j=icen
 2210 continue
      if (y(j) .lt. yr2) then
         j=j+1
         if (j .le. n-1) then
            go to 2210
         else
c hit top of array
            xr2 = x(n)
         end if
      else
c now y(j) is above yr2
         frac = (yr2-y(j-1)) / (y(j)-y(j-1))
         xr2 = frac * x(j) + (1.0-frac) * x(j-1)
      end if
      return
      end


