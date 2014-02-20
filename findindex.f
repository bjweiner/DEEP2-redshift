
c  return jmin and jmax which are the indices at 
c  which the array values in iarr equal imin and imax

      subroutine findindex(n,iarr,imin,imax,jmin,jmax)

      integer iarr(n)

      jmin=1
      jmax=n
c count up to first point
      ii=1
 100  continue
      if(iarr(ii) .ge. imin) then
c success
         jmin = ii
      else if(ii .eq. n) then
c hopefully this won't happen
         jmin = n
      else
         ii=ii+1
         go to 100
      end if
 110  continue

c count down to last point
      ii=n
 200  continue
      if(iarr(ii) .le. imax) then
         jmax = ii
      else if (ii .eq. 1) then
         jmax=1
      else
         ii=ii-1
         go to 200
      end if
 210  continue

      return
      end

         
         

