
c  this is a rebinning/pixelizing routine.
c  given n points and x(i), y(i), rebin into
c  nout pixels that have x(0)=x0, and linear spacing dx.
c  The input arrays x and y are assumed to be sorted in 
c  increasing x.

      subroutine linrebin(n,x,y,nout,x0,dx,yout)

      integer n,nout
      real x(n),y(n),yout(nout)
      real x0,dx
      real weight(nout)

c      include 'pcredshift.h'
      bad = -1.e4
      badset = -1.e6

c for the moment we just calculate x0 and dx from the 
c extreme points.  This might or might not be a bad idea.
c If the x are nearly linear it should be ok. 

      xmin=x(1)
      xmax=x(n)
      dx = (xmax-xmin)/(nout-1.)
      x0 = xmin-dx

      do j=1,nout
         yout(j) = 0.0
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
c boundaries of this cell
         if (i .eq. 1) then
            x1 = x(i)
         else
            x1 = (x(i-1) + x(i)) /2.
         end if
         if (i .eq. n) then
            x2 = x(i)
         else
            x2 = (x(i) + x(i+1)) /2.
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
c We have to keep track of both how much of the old pixel to
c put in each new bin, e.g. (1-frac1)/(r2-r1) for the first new pixel;
c and how much of each new pixel is being covered, this is weight(j)
c
c first handle the special case if the boundaries are
c in the same bin
         if (j1 .eq. j2) then
            yout(j1) = yout(j1) + y(i)
            weight(j1) = weight(j1) + (frac2-frac1)
            go to 200
         end if
c Do the first, fractional pixel
         yout(j1) = yout(j1) + y(i) * (1.-frac1) / (r2-r1)
         weight(j1) = weight(j1) + (1.-frac1) 
c Loop through the whole pixels
         do j=j1+1,j2-1
            yout(j)  = yout(j) + y(i) / (r2-r1)
            weight(j) = weight(j) + 1.
         end do
c Do the last, fractional pixel
         yout(j2) = yout(j2) + y(i) * frac2 / (r2-r1)
         weight(j2) = weight(j2) + frac2

 200     continue
      end do

c  conserve counts (photons) so don't use the weight to normalize

c  Hopefully most of the time we will have covered all of the new
c  pixel so weight(j) will be 1.0 (within roundoff).  However if
c  a bad old pixel was mapped onto part of the new pixel, weight(j)
c  will be less than 1 and yout(j) will be lower than it ought to be.

      do j=2,nout-1
         if(weight(j) .gt. 1.001) then
c            write(*,*) "Weight problem ",j,x0+j*dx,weight(j),yout(j)
         else if (weight(j) .gt. 0.1) then
            yout(j) = yout(j) * 1.0 / weight(j)
         else
c            write(*,*) "Weight problem ",j,x0+j*dx,weight(j),yout(j)
            yout(j) = badset
         end if
      end do

      return
      end
