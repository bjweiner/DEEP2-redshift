c 12 Dec 2002 - fixed a problem which occurred at the inter-CCD gap,
c    where the wavelength of the input pixels takes a jump, but there
c    is really no data.  Now tests for this by effectively checking
c    for an adjacent-pixel spacing that is very different on left and
c    right.

c  this is a rebinning/pixelizing routine.
c  given n points and x(i), y(i), rebin into
c  nout pixels that have x(0)=x0, and linear spacing dx.
c  The input arrays x and y are assumed to be sorted in 
c  increasing x.

c linrebinerr rebins the error at the same time.
c   If  F = sum(a_i f_i),
c   then  var(F) = sum(a_i^2 var(f_i))
c         err(F) = sqrt(sum(a_i^2 err_i^2))

c Because we are looping over old pixels rather than
c new pixels, accumulating the flux, and especially the
c variance, is kind of complicated.

c linrebinerr expects to get the error (rms) as input

      subroutine linrebinerr(n,x,y,yerr,nout,x0,dx,yout,youterr)

c controls whether we rebin the variance "correctly," sum (a_i^2 var_i) 
c  or not, sum (a_i var_i) which seems to be more robust and possibly
c  correct if you think in terms of Poission noise.
      parameter (IFVARCORR=0)

      integer n,nout
      real x(n),y(n),yout(nout)
      real yerr(n), youterr(nout)
      real yvar(n), youtvar(nout)
      real x0,dx
      real weight(nout)

c      include 'pcredshift.h'
      bad = -1.e4
      badset = -1.e6

      do i=1,n
         yvar(i) = yerr(i)**2
      end do

c for the moment we just calculate x0 and dx from the 
c extreme points.  This might or might not be a bad idea.
c If the x are nearly linear it should be ok. 

      xmin=x(1)
      xmax=x(n)
      dx = (xmax-xmin)/(nout-1.)
      x0 = xmin-dx

      do j=1,nout
         yout(j) = 0.0
         youtvar(j) = 0.0
         weight(j) = 0.0
      end do

c loop through OLD pixels and figure out how much
c flux to put in each new pixel and how much of the 
c new pixel is covered.  Hopefully this has a lot in
c common with logrebin, but here I'm looping over the
c OLD bins while in logrebin I looped over the NEW bins.

c I'm trying to conserve total counts (photons or whatever)

      do i=1,n
c if the cell value is bad skip the whole thing
         if(y(i) .lt. bad) go to 200
c find the x boundaries of this cell.  Treat first and last cells
c  i=1, i=n, as special cases
         if (i .eq. 1) then
            x1 = x(1) - (x(2) - x(1))/2.
            x2 = (x(1) + x(2)) /2.            
         else if (i .eq. n) then
            x1 = (x(n-1) + x(n)) /2.
            x2 = x(n) + (x(n) - x(n-1))/2.
         else
            delx1 = (x(i) - x(i-1)) / 2.
            delx2 = (x(i+1) - x(i)) / 2.
c There is going to be a problem because there are no
c "old" pixels in the inter-CCD gap.  Check this by testing for
c delx in one direction significantly larger than the other.
c If so, we have an "end" pixel and just make it symmetric.
c this would be a problem if the variation from linearity was
c really bad.
            if (delx1 .gt. 2.*delx2) then
               delx1 = delx2
            else if (delx2 .gt. 2.*delx1) then
               delx2 = delx1
            end if
            x1 = x(i) - delx1
            x2 = x(i) + delx2
         end if
c boundaries in pixel coords
         r1 = (x1-x0)/dx
         r2 = (x2-x0)/dx
c compute what lin bin the cell boundaries fall in
         j1 = int(r1+0.5)
         j2 = int(r2+0.5)
c just trimming to the min/max is not ideal but leave it for now
         j1 = min(max(j1,1),nout)
         j2 = min(max(j2,1),nout)
         frac1 = r1 - (j1-0.5)
         frac2 = r2 - (j2-0.5)
c for testing
c         write(*,'(i4,2(2x,f8.3),2(2x,i4),2(2x,f7.4))') i,r1,r2,
c     $        j1,j2,frac1,frac2
c We have to keep track of both how much of the old pixel to
c put in each new bin, e.g. (1-frac1)/(r2-r1) for the first new pixel;
c and how much of each new pixel is being covered, this is weight(j)
c
c first handle the special case if the boundaries are
c in the same bin
         if (j1 .eq. j2) then
            yout(j1) = yout(j1) + y(i)
            if (IFVARCORR .eq. 0) then
               youtvar(j1) = youtvar(j1) + yvar(i)
            else
               youtvar(j1) = youtvar(j1) + yvar(i)
            end if
            weight(j1) = weight(j1) + (frac2-frac1)
            go to 200
         end if
c Do the first, fractional pixel
         yout(j1) = yout(j1) + y(i) * (1.-frac1) / (r2-r1)
         if (IFVARCORR .eq. 0) then
            youtvar(j1) = youtvar(j1) + yvar(i) * (1.-frac1) / (r2-r1)
         else            
            youtvar(j1) = youtvar(j1) + yvar(i) * 
     $           ((1.-frac1) / (r2-r1))**2
         end if
         weight(j1) = weight(j1) + (1.-frac1) 
c Loop through the whole pixels
         do j=j1+1,j2-1
            yout(j)  = yout(j) + y(i) / (r2-r1)
            if (IFVARCORR .eq. 0) then
               youtvar(j)  = youtvar(j) + yvar(i) / (r2-r1)
            else
               youtvar(j)  = youtvar(j) + yvar(i) / ((r2-r1)**2)
            end if
            weight(j) = weight(j) + 1.
         end do
c Do the last, fractional pixel
         yout(j2) = yout(j2) + y(i) * frac2 / (r2-r1)
         if (IFVARCORR .eq. 0) then
            youtvar(j2) = youtvar(j2) + yvar(i) * frac2 / (r2-r1)
         else
            youtvar(j2) = youtvar(j2) + yvar(i) * (frac2 / (r2-r1))**2
         end if
         weight(j2) = weight(j2) + frac2

 200     continue
      end do

c  conserve counts (photons) so don't use the weight to normalize

c  Hopefully most of the time we will have covered all of the new
c  pixel so weight(j) will be 1.0 (within roundoff).  However if
c  a bad old pixel was mapped onto part of the new pixel, weight(j)
c  will be less than 1 and yout(j) will be lower than it ought to be.

c youtvar is the variance of the accumulated yout, so
c the error youterr is scaled by the same factor as yout

      do j=2,nout-1
         if(weight(j) .gt. 1.001) then
c this shouldn't happen, but does
c            write(*,*) "Weight problem > 1 in linrebin ",
c     $           j,x0+j*dx,weight(j),yout(j)
c don't adjust yout(j), and youterr does not need to be rescaled
            if (youtvar(j) .gt. 1.0e-8) then
               youterr(j) = sqrt(youtvar(j))
            else
               youterr(j) = 1.0e10
            end if
c 0.05 is arbitrary - only use an output pixel if it has at least 5% real data
         else if (weight(j) .gt. 0.05) then
            yout(j) = yout(j) * 1.0 / weight(j)
            if (youtvar(j) .gt. 1.0e-8) then
            if (IFVARCORR .eq. 0) then
               youterr(j) = sqrt( youtvar(j) * 1.0 / weight(j) )
            else
               youterr(j) = sqrt(youtvar(j)) * 1.0 / weight(j)
            end if
            else
               youterr(j) = 1.0e10
            end if
         else
c this may be happening occasionally 
c and should happen in the gap
c            write(*,*) "Weight problem < 0.05 in linrebin"
c     $           ,j,x0+j*dx,weight(j),yout(j)
            yout(j) = badset
            youterr(j) = 1.0e10
         end if
      end do

      return
      end
