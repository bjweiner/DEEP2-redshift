
c  median smooth a spectrum with a window of diameter = iwidth pixels
c  the original spectrum is preserved.

c  11.08.02 - fixed an error, it was setting sout(j) rather than sout(i)
c     this increased flattopping of emission lines, etc, and was worse
c     with larger smoothing

      subroutine medsmooth(n,spin,iwidth,sout)

      real spin(n),sout(n)
      integer iwidth
      real temparr(n)
      external select

      include 'pcredshift.h'

c  if smoothing is 1 or 0 pixels, array is unchanged
      if (iwidth .le. 1) then
         do i=1,n
            sout(i) = spin(i)
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
            if(spin(j) .gt. bad) then
               np = np+1
               temparr(np) = spin(j)
            end if
         end do
         if (np .eq. 0) then
            sout(i) = badset
         else if (np .eq. 1) then
            sout(i) = temparr(1)
         else if (np .eq. 2) then
            sout(i) = (temparr(1) + temparr(2))/2.0
         else
c do the median right for once
            itmp = np/2
            if (mod(np,2) .ne. 0) then
               sout(i) = select(itmp+1,np,temparr)
            else 
               sout(i) = (select(itmp,np,temparr) + 
     $              select(itmp+1,np,temparr)) / 2.0
            end if
         end if
      end do

      return
      end
