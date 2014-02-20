
c rebin a spectrum from linear into log w space
c  This is a derivation from logrebin.  It does the
c  same thing but also readjusts the noise (std dev)
c  spectrum at the same time.  I just square the noise
c  to get variance and do the same additions. 
c  This means adjacent bins are no longer independent,
c  but I'm not sure it's a real problem when the dominant
c  noise sources are read noise and Poisson.

c note that we ASSUME that there is a good noise value
c for every pixel where there's a good flux value.
c this could be modified but that might also cause problems.

      subroutine logrebin2(spec,nsp,espec,w0,dw,nslog,w0log,dwlog,
     $     slog,eslog)

      real spec(nsp), slog(nslog)
      real espec(nsp), eslog(nslog)
      real varspec(nsp)
      real w0,dw,w0log,dwlog
      integer nsp,nslog

      include 'pcredshift.h'

      do i=1,nspec
c         if (espec(i) .gt. bad) then
         if (spec(i) .gt. bad) then
            varspec(i) = espec(i)*espec(i)
         else
            varspec(i) = espec(i)
         end if
      end do

c      wmin = w0 + 0.5*dw
c      wmax = w0 + (nsp+0.5)*dw

c  assume that the bins of both linear and log spectra
c  are centered on the appropriate wavelength values.
c  this means throwing around a lot of 0.5's to get the
c  bin boundary values.

      do i=1,nslog
         fsum = 0.0
         wtgood = 0.0
         wtbad = 0.0
         efsum = 0.0
c  boundaries of log bin
         w1 = 10**(w0log + (i-0.5)*dwlog)
         w2 = 10**(w0log + (i+0.5)*dwlog)
c compute what lin bin the log boundaries fall in
         j1 = int( (w1-w0)/dw +0.5 )
         j2 = int( (w2-w0)/dw +0.5 )
         j1 = min(max(j1,1),nspec)
         j2 = min(max(j2,1),nspec)
         frac1 = (w1 - (w0+(j1-0.5)*dw)) / dw
         frac2 = (w2 - (w0+(j2-0.5)*dw)) / dw
c         if (j1 .ge. 1 .and. j1 .le. nsp .and. spec(j1) .gt. bad) then
         if (spec(j1) .gt. bad) then
            wtgood = wtgood + (1.-frac1)
            fsum = fsum + (1.-frac1) * spec(j1)
            efsum = efsum + (1.-frac1) * varspec(j1)
         else
            wtbad = wtbad + (1.-frac1)
         end if
         do j=j1+1,j2-1
c            if (j .ge. 1 .and. j .le. nsp .and. spec(j) .gt. bad) then
            if (spec(j) .gt. bad) then
               wtgood = wtgood + 1.
               fsum  = fsum + spec(j)
               efsum = efsum + varspec(j)
            else
               wtbad = wtbad + 1.
            end if
         end do
c         if (j2 .ge. 1 .and. j2 .le. nsp .and. spec(j2) .gt. bad) then
         if (spec(j2) .gt. bad) then
            wtgood = wtgood + frac2
            fsum = fsum + frac2 * spec(j1)
            efsum = efsum + frac2 * varspec(j1)
         else
            wtbad = wtbad + frac2
         end if
c  conserve counts (photons) so don't use the weight to normalize
c  this approach has a problem if, say, one pixel is bad
c         slog(i) = fsum
c  so rather, we estimate the flux that belonged in this bin by scaling
c  total/good
         if (wtgood .gt. 1.e-3) then
            slog(i) = fsum * (wtgood + wtbad)/wtgood
            eslog(i) = sqrt(efsum * (wtgood + wtbad)/wtgood)
         else
            slog(i) = bad
         end if

      end do

      return
      end

