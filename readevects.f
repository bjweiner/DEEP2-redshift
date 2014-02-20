c read the wavelengths, mean vector, and eigenvectors from a file
c such as that made by pcatest

      subroutine readevects(evects,nvphys,nwphys,wlrest,nv,nw)

      real evects(nvphys,nwphys)
      real wlrest(nwphys)
      character fname*200

      include 'specarray.h'

      include 'pcredshift.h'

 100  continue
      write(*,'("File with w.l. and eigenvectors [quit]: ",$)')
      read (*,'(a)') fname
      if (fname(1:3) .eq. '   ') go to 666
      open(2,file=fname,status='old',err=100)

 120  continue
      write(*,'("# of vectors to read (includes mean): ",$)')
      read (*,*) nv

      nw=1
     
 200  continue
      read(2,*,err=666,end=300) wlrest(nw),(evects(i,nw),i=1,nv)
      nw=nw+1
      go to 200

 300  continue
      close(2)
      nw = nw-1
      write(*,*) "Read ",nw," wavelength pixels"

      return
      
 666  continue
      close(2)
      write(*,*) "Error reading eigenvectors - not enough columns?"
      stop

      end
