c find the centroid of a minimum in chisq given its
c approx location.  indx is the index of the minimum
c we found with findmultmin or whatever

      function centroidmin(n,x,y,indx)

      integer indx
      real x(n),y(n)
c      real centroidmin
      
c look in a window this size around the min
      nwtmin = 5
c buffer
      nbuff = 3
c window for "continuum"
      ncontwin = 10
      
      csum=0.0
      ncont=0
      i1 = max(1,indx-nwtmin-nbuff-ncontwin)
      i2 = max(1,indx-nwtmin-nbuff-1)
      do i=i1,i2
         csum=csum+y(i)
         ncont=ncont+1
      end do
      i1= min(n,indx+nwtmin+nbuff+1)
      i2= min(n,indx+nwtmin+nbuff+ncontwin)
      do i=i1,i2
         csum = csum+y(i)
         ncont=ncont+1
      end do
      cont=csum/ncont

c this is like finding a center of mass:
c take the y-weighted mean of the x_i
c (but make the weights y minus continuum)

      xysum = 0.0
      wtsum = 0.0
      i1 = min(n,indx-nwtmin)
      i2 = min(n,indx+nwtmin)
      do i=i1,i2
         wt = y(i)-cont
         xysum = xysum+wt*x(i)
         wtsum = wtsum+wt
      end do
      
      centroidmin = xysum/wtsum
      return
      end

