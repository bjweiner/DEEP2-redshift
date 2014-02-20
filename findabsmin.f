
c  find the absolute minimum value in an array

      subroutine findabsmin(n,arr,imin,amin)

      real arr(n)

      amin = 1.0e30
      do i=1,n
         if(arr(i) .lt. amin) then
            imin=i
            amin=arr(i)
         end if
      end do

      return
      end

