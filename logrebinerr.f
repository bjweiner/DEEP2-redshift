
c rebin a spectrum from linear into log w space
c  This is a derivation from logrebin.  It does the
c  same thing but also readjusts the noise (std dev)
c  spectrum at the same time.  If  F = sum(a_i f_i),
c   then  var(F) = sum(a_i^2 var(f_i))
c         err(F) = sqrt(sum(a_i^2 err_i^2))
c  This means adjacent bins are no longer independent,
c  but I'm not sure it's a real problem when the dominant
c  noise sources are read noise and Poisson.

      subroutine logrebinerr(spec,espec,nsp,w0,dw,nslog,w0log,dwlog,
     $     slog,eslog)

c controls whether we rebin the variance "correctly," sum (a_i^2 var_i) 
c  or not, sum (a_i var_i) which seems to be more robust and possibly
c  more correct if you think in terms of Poisson noise (IFVARCORR=0).
      parameter (IFVARCORR=0)

      real spec(nsp), slog(nslog)
      real espec(nsp), eslog(nslog)
      real varspec(nsp)
      real w0,dw,w0log,dwlog
      integer nsp,nslog

      include 'pcredshift.h'

c      wmin = w0 + 0.5*dw
c      wmax = w0 + (nsp+0.5)*dw

c note that we ASSUME that there is a good noise value
c for every pixel where there's a good flux value.
c this could be modified but that might also cause problems.
      do i=1,nsp
c         if (espec(i) .gt. bad) then
c         if (spec(i) .gt. bad) then
            varspec(i) = espec(i)*espec(i)
c         else
c            varspec(i) = espec(i)
c         end if
      end do

c  assume that the bins of both linear and log spectra
c  are centered on the appropriate wavelength values.
c  this means throwing around a lot of 0.5's to get the
c  bin boundary values.

      do i=1,nslog
         fsum = 0.0
         varsum = 0.0
         wtgood = 0.0
         wtbad = 0.0
c  boundaries of log bin
         w1 = 10**(w0log + (i-0.5)*dwlog)
         w2 = 10**(w0log + (i+0.5)*dwlog)
c  corresponding pixel coords in the linear spectrum
         r1 = (w1-w0)/dw
         r2 = (w2-w0)/dw
c compute what lin bin the log boundaries fall in
c         j1 = int( (w1-w0)/dw +0.5 )
c         j2 = int( (w2-w0)/dw +0.5 )
         j1 = int(r1+0.5)
         j2 = int(r2+0.5)
c if we're completely out of the region where there is data in
c the linear spectrum, bail
         if (j2 .lt. 1 .or. j1 .gt. nsp) go to 200
            
c just trimming to the min/max is not ideal but leave it for now
         j1 = min(max(j1,1),nsp)
         j2 = min(max(j2,1),nsp)
c         frac1 = (w1 - (w0+(j1-0.5)*dw)) / dw
c         frac2 = (w2 - (w0+(j2-0.5)*dw)) / dw
         frac1 = r1 - (j1-0.5)
         frac2 = r2 - (j2-0.5)
c I can't see how this would happen
c         if (frac1 .gt. 1.0 .or. frac2 .gt. 1.0) then
c            write(*,*) "warning in logrebin: ",
c     $           frac1,frac2,i,r1,r2,j1,j2
c         end if

c first handle the special case if the boundaries are
c in the same bin
         if (j1 .eq. j2) then
            if (spec(j1) .gt. bad) then
               wtgood = frac2 - frac1
               fsum = (frac2-frac1) * spec(j1)
               if (IFVARCORR .eq. 0) then
                  varsum = (frac2-frac1) * varspec(j1)
               else
                  varsum = (frac2-frac1)**2 * varspec(j1)
               end if
            else
               wtbad = frac2-frac1
            end if
            go to 200
         end if

c Do the first, fractional pixel
c         if (j1 .ge. 1 .and. j1 .le. nsp .and. spec(j1) .gt. bad) then
         if (spec(j1) .gt. bad) then
            wtgood = wtgood + (1.-frac1)
            fsum = fsum + (1.-frac1) * spec(j1)
            if (IFVARCORR .eq. 0) then
               varsum = varsum + (1.-frac1) * varspec(j1)
            else
               varsum = varsum + (1.-frac1)**2 * varspec(j1)
            end if
         else
            wtbad = wtbad + (1.-frac1)
         end if
c Loop through the whole pixels
         do j=j1+1,j2-1
c            if (j .ge. 1 .and. j .le. nsp .and. spec(j) .gt. bad) then
            if (spec(j) .gt. bad) then
               wtgood = wtgood + 1.
               fsum  = fsum + spec(j)
               if (IFVARCORR .eq. 0) then
                  varsum  = varsum + varspec(j)
               else
                  varsum  = varsum + varspec(j)
               end if
            else
               wtbad = wtbad + 1.
            end if
         end do
c Do the last, fractional pixel
c         if (j2 .ge. 1 .and. j2 .le. nsp .and. spec(j2) .gt. bad) then
         if (spec(j2) .gt. bad) then
            wtgood = wtgood + frac2
            fsum = fsum + frac2 * spec(j1)
            if (IFVARCORR .eq. 0) then
               varsum = varsum + frac2 * varspec(j1)
            else
               varsum = varsum + frac2**2 * varspec(j1)
            end if
         else
            wtbad = wtbad + frac2
         end if

 200     continue

c  conserve counts (photons) so don't use the weight to normalize
c  this approach has a problem if, say, one pixel is bad
c         slog(i) = fsum
c  so rather, we estimate the flux that belonged in this bin by scaling
c  total/good
         if (wtgood .gt. 1.e-3) then
            slog(i) = fsum * (wtgood + wtbad)/wtgood
c varsum is the variance of fsum, so we have to scale
c the _error_ by the same ratio of weights as we scaled fsum
            if (IFVARCORR .eq. 0) then
               eslog(i) = sqrt(varsum * (wtgood + wtbad)/wtgood )
            else
               eslog(i) = sqrt(varsum) * (wtgood + wtbad)/wtgood
            end if
         else
            slog(i) = badset
            eslog(i) = 1.0e10
         end if

      end do

      return
      end

