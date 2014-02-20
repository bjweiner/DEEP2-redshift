
c  find the biggest local minimum value in an array
c  by looking for the difference between the value and
c  a moving average around it.  Return the index of the 
c  minimum, min value, and avg value at the minimum.

c  It would probably be a good idea to look for the
c  minimum with the greatest area rather than just 
c  the deepest pixel value, but that's more involved.

      subroutine findlocmin(n,arr,minbest,amin,avgbest)

      real arr(n)

c radius of window in which to take average
      iwin = 10
      dwin = 2.0*iwin
c radius of buffer around point in question
      ibuff = 5
      irad = iwin+ibuff

      jmin=1+irad
      jmax=n-irad

c calculate the difference value at jmin
      sum1=0.0
      sum2=0.0
c      np = 0
      do i=1,jmin-ibuff-1
c I could put a test for >bad here
         sum1=sum1+arr(i)
      end do
      do i=jmin+ibuff+1,jmin+irad
         sum2=sum2+arr(i)
      end do
      avg = (sum1+sum2) / dwin
      diffmin = arr(jmin) - avg
      minbest = jmin

      do j=jmin+1,jmax
c  four points move in/out of the windows
         sum1 = sum1 - arr(j-irad-1)
         sum1 = sum1 + arr(j-ibuff-1)
         sum2 = sum2 - arr(j+ibuff)
         sum2 = sum2 + arr(j+irad)
         avg = (sum1 + sum2) / dwin
         diff = arr(j) - avg
         if (diff .lt. diffmin) then
c check to make sure it's actually lower then each of the windows,
c rather than being the point at the bottom of a "cliff"
            if (arr(j) .lt. sum1/iwin .and.
     $           arr(j) .lt. sum2/iwin) then
               diffmin = diff
               minbest = j
               avgbest = avg
            end if
         end if
      end do

      amin = arr(minbest)

      return
      end

