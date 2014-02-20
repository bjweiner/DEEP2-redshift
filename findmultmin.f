
c  derived from findlocmin - same idea, but save the
c  N deepest minima

c  find the deepest N local minimum values in an array
c  by looking for the difference between the value and
c  a moving average around it.  Return the index of the 
c  minimum, min value, and avg value at the minimum.

c  It would probably be a good idea to look for the
c  minimum with the greatest area rather than just 
c  the deepest pixel value, but that's more involved.

      subroutine findmultmin(n,arr,nfind,mini,diffarr)

      real arr(n)
      integer mini(nfind)
      real diffarr(nfind)

c tolerance - don't want to return several values that are
c right next to each other so force them to be this far apart:
      itol = 8

c initialize the return array
      do i=1,nfind
         mini(i) = 1
         diffarr(i) = 0.0
      end do      
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
c      diffmin = arr(jmin) - avg
      diffmin = arr(jmin) - avg
      mini(1) = jmin
      diffarr(1) = diffmin

c  now loop through the data

      do j=jmin+1,jmax
c  four points move in/out of the windows
         sum1 = sum1 - arr(j-irad-1)
         sum1 = sum1 + arr(j-ibuff-1)
         sum2 = sum2 - arr(j+ibuff)
         sum2 = sum2 + arr(j+irad)
         avg = (sum1 + sum2) / dwin
         diff = arr(j) - avg

c  if it's less than the nfind'th deepest minimum
         if (diff .lt. diffarr(nfind)) then
c check to make sure it's actually lower then each of the windows,
c rather than being the point at the bottom of a "cliff"
            if (arr(j) .lt. sum1/iwin .and.
     $           arr(j) .lt. sum2/iwin) then

c check if it's within itol of one of the values.
               k=1
 110           continue
c  test if too close to an old min
               if (abs(j-mini(k)) .lt. itol) then
c  if the new min is lower, replace the old
                  if (diff .lt. diffarr(k)) then
                     mini(k)=j
                     diffarr(k) = diff
                  end if
c  now bail out of the if clauses, going to the end
c  of the do loop.  This is probably causing a problem 
c  because we haven't resorted the minima, so I add
c  a sort here and at the end.
                  call piksr2int(nfind,diffarr,mini)
                  go to 300
               else
                  k=k+1
                  if (k .le. nfind) go to 110
               end if

c so now we have a new min that is not within itol               
c find where it belongs among the nfind values
               k=1
 150           continue
c if we haven't fallen off the end
               if (k .le. nfind) then
                  if(diff .gt. diffarr(k)) then
c advance to next index
                     k=k+1
                     go to 150
                  else
c move the old values up one index                     
                     do l=nfind,k+1,-1
                        mini(l) = mini(l-1)
                        diffarr(l) = diffarr(l-1)
                     end do
c insert new value
                     mini(k) = j
                     diffarr(k) = diff
                  end if
               end if

            end if
         end if
c end the if clauses for finding a new min

 300     continue
      end do

c      amin = arr(minbest)

c  before returning, make sure the minima really are sorted
c      call piksr2int(nfind,diffarr,mini)

      return
      end

c  Two-array insertion sort from Numerical Recipes
c  modified so the second array is integer.  Insertion sort
c  should not be used for N>20 or so, which should not be a
c  problem here.

      SUBROUTINE piksr2int(n,arr,brr)
      INTEGER n
      REAL arr(n)
      integer brr(n)
      INTEGER i,j
      REAL a
      integer b
      do 12 j=2,n
        a=arr(j)
        b=brr(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
11      continue
        i=0
10      arr(i+1)=a
        brr(i+1)=b
12    continue
      return
      END
