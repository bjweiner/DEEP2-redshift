
c retrieve a 1-d Berkeley format BINTABLE spectrum
c and write linearized version to ascii file.

      program writedeepspec

      include 'specarray.h'

      real wave(NWAVEMAX),spec(NWAVEMAX),serr(NWAVEMAX)
      character fname*160,oname*160

      call pgbeg(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.2)

 100  continue
      write(*,'("1-d spec fits file [quit]: ",$)')
      read(*,'(a)') fname
      if (fname(1:5) .eq. '     ') go to 666
         
 200  continue
      write(*,'("output ascii file [writedeepspec.out]: ",$)')
      read(*,'(a)') oname
      if (oname(1:5) .eq. '     ') then
         oname='writedeepspec.out'
      end if
      open(2,file=oname,status='unknown',err=200)
      
      call get1ddeepspec(fname,nsp,spec,serr,w0,dw,ifer)
      if (ifer .ne. 0) then
         write(*,*) "Error reading spectrum!"
      end if

      do i=1,nsp
         wave(i) = w0 + i*dw
         write(*,'(i6,3x,f9.3,2(3x,f10.4))') i,wave(i),spec(i),serr(i)
      end do
      close(2)

      call showspecerr(nsp,wave,spec,serr)
      call pglabel("linearized wavelength","counts and error",fname)
      go to 100

 666  continue
      call pgend()

      end


      
