
      subroutine readz(n,zarr)

      real zarr(n)

      character*60 fname

 100  continue
      write(*,'("File with redshifts [use 0]: ",$)')
      read(*,'(a)',err=100) fname
      if (fname(1:3) .eq. '   ') then
         do i=1,n
            zarr(i) = 0.0
         end do
         return
      end if

      open(3,file=fname,status='old',err=100)

      i=1
 200  continue
      read(3,*,end=300) zarr(i)
      i=i+1
      go to 200

 300  continue
      if (i .le. n) then
         write (*,*) "Warning - didnt get enough redshifts?",i,n
         do j=i,n
            zarr(j) = 0.0
         end do
      end if
      close(3)

      return
      end
