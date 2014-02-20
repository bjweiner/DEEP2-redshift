
c do continuum subtraction by subtracting a low order Legendre polynomial
c fit to the good data points. 
c I added a sigma-clipping routine before the fit, which may not
c work great since the continuum varies, avsigclip or something
c would be better.

c do the subtraction in place, or not
      subroutine contsubleg(wave,spec,serr,nsp,iord)

c      parameter (IORDER=15)

      real wave(nsp),spec(nsp),serr(nsp)
      real wscale(nsp),specclip(nsp)
      real wfit(nsp),ffit(nsp),efit(nsp)
      real a(iord),uu(nsp,iord),vv(iord,iord),ww(iord)
      real avals(iord)

      external fleg

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
     $     chisq,fleg)

c evaluate and subtract fit
      do i=1,nsp
c         call fleg(wave(i),avals,npoly)
         call fleg(wscale(i),avals,npoly)
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
c from Numerical Recipes, generates Legendre polynomials
c via recurrence relation

      SUBROUTINE fleg(x,pl,nl)
      INTEGER nl
      REAL x,pl(nl)
      INTEGER j
      REAL d,f1,f2,twox
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      END
