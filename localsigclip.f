
c do a clipping based on the local mean and rms value
c and return an array with clipped values set to a bad flag
c the output array is suitable for eg continuum fitting

c averaging window radius irad, number of sigma to clip is clipsig
c note this has window in pixels, not in x

c      subroutine localsigclip(n,x,y,yout,irad,clipsig)
      subroutine localsigclip(n,y,yout,irad,clipsig)

      real y(n),yout(n)
      integer n,irad
      real clipsig

      include 'pcredshift.h'

c j is index of point under consideration
      j=1
c accumulate sums for the first point
      sum = 0.0
      sumsq = 0.0
      np = 0
      do i=1,irad+1
         if (y(i) .gt. bad) then
            np=np+1
            sum = sum+y(i)
            sumsq=sumsq+y(i)*y(i)
         end if
      end do

 100  continue
c calculate rms and do clipping
      if (np .gt. 1) then
         smean = sum/np
         srms = sqrt((sumsq - sum*sum/np)/np)
         clip1 = smean - clipsig*srms
         clip2 = smean + clipsig*srms
         if (y(j) .gt. clip1 .and. y(j) .lt. clip2) then
            yout(j) = y(j)
         else
            yout(j) = badset
         end if
      else
         yout(j) = badset
      end if
      
      j=j+1
c test if we're done
      if (j .gt. n) go to 200
c if not, remove the last point from sums
      jlast = j-irad
      if (jlast .ge. 1 .and. y(jlast) .gt. bad) then
         sum = sum - y(jlast)
         sumsq = sumsq - y(jlast)*y(jlast)
         np = np-1
      end if
c add new point to sums
      jnew = j+irad
      if(jnew .le. n .and. y(jnew) .gt. bad) then
         sum = sum + y(jnew)
         sumsq = sumsq + y(jnew)*y(jnew)
         np = np+1
      end if
      go to 100

c we've reached the end
 200  continue
      return
      end
