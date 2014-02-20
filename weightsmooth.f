
c  smooth a spectrum with a window of diameter = iwidth pixels
c  using inverse variance weighting
c  return the smoothed spectrum and new rms
c  the original spectrum is preserved.

c Weighted mean generally:
c   xmean = sum(w_i * x_i) / sum(w_i)
c   sigma(xmean) = sqrt(sum(w_i^2 * sigma_i^2)) / sum(w_i)
c here the weight is w_i = 1/sigma_i^2
c   so sigma(xmean) = (sum(1/sigma_i^2))^-0.5 = (sum(w_i))^-0.5

c Note that the resulting pixels are correlated, so this might
c not be ideal.  I had some problems in the linrebinerr and logrebinerr
c routines.  One can also rebin the variance just like the flux
c (as if Poisson noise), and if IFVARCORR=0 this happens.
c  then var(xmean) = sum(w_i * var_i) / sum(w_i)

      subroutine weightsmooth(n,spin,rmsin,iwidth,spout,rmsout)

      parameter(IFVARCORR=0)

      real spin(n),spout(n)
      real rmsin(n),rmsout(n)
      integer iwidth
c      real temparr(n)
c      external select

      include 'pcredshift.h'

c high cut on input rms, i.e. anything with very 
c large rms is ignored to avoid underflows, very low
c rms is assumed to be specious
      rmscuthi = 1.e10
      rmscutlo = 1.e-3
      
c  if smoothing is 1 or 0 pixels, array is unchanged
      if (iwidth .le. 1) then
         do i=1,n
            spout(i) = spin(i)
            rmsout(i) = rmsin(i)
         end do
         return
      end if

c compute radii for smoothing window
c if iwidth is odd we get a symmetric window, if it's even,
c truncate low side down, to
c make the high side the long side (no particular reason)
      irlow = int((iwidth-1)/2)
      irhigh = iwidth-irlow-1
      
      do i=1,n
         i1=max(1,i-irlow)
         i2=min(n,i+irhigh)
         np=0
         do j=i1,i2
            sum = 0.0
            varsum = 0.0
            wtsum = 0.0
            if(spin(j) .gt. bad .and. rmsin(j) .gt. rmscutlo
     $           .and. rmsin(j) .lt. rmscuthi) then
               np = np+1
               var = rmsin(j)**2
               weight = 1. / var
               sum = sum + spin(j)*weight
               wtsum = wtsum + weight
               if (IFVARCORR .eq. 0) then
                  varsum = varsum + var*weight
c               else
c                  varsum = varsum + var* weight**2
               end if
            end if
         end do
         if (np .eq. 0) then
            spout(i) = badset
            rmsout(i) = rmsin(i)
         else 
            spout(i) = sum / wtsum
            if (IFVARCORR .eq. 0) then
               rmsout(i) = sqrt(varsum / wtsum)
            else
c               rmsout(i) = sqrt(varsum) / wtsum
               rmsout(i) = 1. / sqrt(wtsum)
            end if
         end if
      end do

      return
      end
