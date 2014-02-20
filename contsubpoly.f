
c do continuum subtraction by subtracting a low order polynomial
c fit to the good data points. 
c I added a sigma-clipping routine before the fit, which may not
c work great since the continuum varies, avsigclip or something
c would be better.

c do the subtraction in place, or not
c      subroutine contsubmed(spec,nsp,dwlog,wrad,specout)
      subroutine contsubpoly(wave,spec,serr,nsp,iord)

c      parameter (IORDER=15)

      real wave(nsp),spec(nsp),serr(nsp)
      real wscale(nsp),specclip(nsp)
      real wfit(nsp),ffit(nsp),efit(nsp)
      real a(iord),uu(nsp,iord),vv(iord,iord),ww(iord)
      real avals(iord)

      external polyfuncs

      include 'pcredshift.h'

      ncliprad=100
      clipsig=2.0
      call localsigclip(nsp,spec,specclip,ncliprad,clipsig)

      wrange = wave(nsp) - wave(1)
      do i=1,nsp
c rescale wfit to run between 0 and 1 for improved polynomial fitting
         wscale(i) = (wave(i) - wave(1)) / wrange
      end do

      npoly=iord
      nfit=0
      do i=1,nsp
         if(specclip(i) .gt. bad .and. serr(i) .gt. bad 
     $        .and. serr(i) .lt. 1.0e6) then
            nfit=nfit+1
            wfit(nfit)=wscale(i)
            ffit(nfit)=specclip(i)
            efit(nfit)=serr(i)
         end if
      end do
      if (nfit .le. npoly) return

c do fit
      call bigsvdfit(wfit,ffit,efit,nfit,a,npoly,uu,vv,ww,nsp,iord,
     $     chisq,polyfuncs)

c evaluate and subtract fit
      do i=1,nsp
c         call polyfuncs(wave(i),avals,npoly)
         call polyfuncs(wscale(i),avals,npoly)
         fitval = 0.
         do j=1,npoly
            fitval = fitval + a(j)*avals(j)
         end do
         if (spec(i) .gt. bad) then
            spec(i) = spec(i) - fitval
         end if
      end do

      return
      end

cccccccccccccccccccccccccccccc
      subroutine polyfuncs(x,a,nord)
c IPOLY: 1 for simple polynomial: 1,x,x^2,...
      parameter(IPOLY=1)

      real a(nord)

      if (IPOLY .eq. 1) then
         a(1) = 1.0
         do i=2,nord
            a(i) = a(i-1)*x
         end do
      end if

      return
      end
