
c given y1(x1) interpolate onto the x values in x2, returning y2
c x1 and x2 are assumed sorted increasing

      subroutine lininterp(n1,x1,y1,n2,x2,y2)

      real x1(n1),y1(n1),x2(n2),y2(n2)
      external xinterp
      real xinterp

      i1 = 1
      j2 = 1

 100  continue
      if (x2(j2) .lt. x1(i1)) then
         y2(j2) = y1(i1)
         j2 = j2+1
         go to 100
      end if

 200  continue
      if (x2(j2) .gt. x1(i1)) then
         if (x2(j2) .le. x1(i1+1)) then
            y2(j2) = xinterp(x1(i1),x1(i1+1),y1(i1),y1(i1+1),x2(j2))
            j2 = j2 + 1
            if (j2 .gt. n2) go to 666
            go to 200
         else
            i1 = i1 + 1
            if (i1 .gt. n1) go to 300
            go to 200
         end if
      else
c probably shouldn't happen
         go to 100
      end if


 300  continue
c get here if we're off the end of the first array
      do j=j2,n2
         y2(j) = y1(n1)
      end do

 666  continue
      return
      end

cccccccccccccccccccccccccccccc

      function xinterp(x1,x2,y1,y2,xval)

      real x1,x2,y1,y2,xval

      if (xval .le. x1) then
         xinterp = y1
      else if (xval .ge. x2) then
         xinterp = y2
      else if (x2 .eq. x1) then
c shouldn't happen
         xinterp = (y1+y2)/2.
      else
         frac = (xval-x1) / (x2-x1)
         xinterp = (1.-frac)*y1 + frac*y2
      end if

      return
      end

         
