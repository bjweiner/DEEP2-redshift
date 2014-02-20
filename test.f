
      program test

      parameter (maxz=8)
      character fname*80,sname*80,cline*80,answ*1
      real zcand(maxz)

      nzmax = 5
      if (nzmax .gt. maxz) nzmax=maxz

 130  write(*,'("File with candidate redshifts: ",$)') 
      read(*,'(a)') fname
c      if (fname(1:3) .eq. '   ') stop
      open(4,file=fname,status='old',err=130)

 200  continue
      read(4,'(a)',end=666) cline

      write(*,'(a)') cline
      nz=nzmax
c  readz resets nz
c      call readz(cline,nz,zcand)
      call readz2(cline,nz,zcand)

      write(*,*) nz
      write(*,*) (zcand(i),i=1,nz)

      go to 200

 666  continue
      close(4)
      end

c  somehow read an arbitrary number of z's out of a line
      subroutine readz(cline,nz,zarr)

      character cline*80
      real zarr(nz)

      do i=1,nz
         zarr(i)=0.0
      end do
      read(cline,*,err=700,end=700) (zarr(i),i=1,nz)
 700  continue
      i=1
 720  continue
      if (zarr(i) .ne. 0.) then
         go to 720
      else
         nz=i-1
      end if

      return
      end

c  read the number of z's and then read the z's
      subroutine readz2(cline,nz,zarr)
      character cline*80
      real zarr(nz)

      read(cline,*) nz
      read(cline,*,err=800,end=800) idummy,(zarr(i),i=1,nz)
 800  continue

      return
      end
