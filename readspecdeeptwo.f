
c  read a bunch of spectra
c  this version reads from the Berkeley pipeline DEEP 2 spectra
c  fits tables

c      subroutine readspec(nsp,nwave)
c change this so the array of data is an argument,
c the next two are the physical dimensions, then the actual
c number of spectra is returned in nspec, and the actual
c numbers of wavelengths are returned in the nwarray array
      subroutine readspecdeeptwo(specarr,skyspecarr,nsphys,nwphys,
     $     nsp,nwarray,
     $     specw0,specdw,specname)

      real specarr(nsphys,nwphys)
      real skyspecarr(nsphys,nwphys)
      real specw0(nsphys),specdw(nsphys)
      character specname(nsphys)*60
      integer nwarray(nsp)
      character fname*80,sname*80

      include 'specarray.h'

c      common /specarray/ specarr(NSPECMAX,NWAVEMAX)
c      common /specdata/ specw0(NSPECMAX),specdw(NSPECMAX)
c      common /specname/ specname
c      character specname(NSPECMAX)*60

      real tempdat(NWAVEMAX),temperr(NWAVEMAX)
      integer isize(7)

      include 'pcredshift.h'

      nsmax = NSPECMAX
      nwmax = NWAVEMAX

 100  continue
      write(*,'("File with list of spectra [quit]: ",$)')
      read (*,'(a)') fname
      if (fname(1:3) .eq. '   ') go to 666
      open(2,file=fname,status='old',err=100)

      nsp = 0
 200  continue
      read(2,'(a)', err=250,end=250) sname
      
      call get1ddeepspec(sname,nx,tempdat,temperr,w0,dw,ifer)
      if (ifer .ne. 0) go to 200
      nsp=nsp+1
      nwarray(nsp)=nx
      specname(nsp)=sname
      specw0(nsp)=w0
      specdw(nsp)=dw
      
      do i=1,nx
         specarr(nsp,i) = tempdat(i)
         skyspecarr(nsp,i) = temperr(i)
      end do

      go to 200

 250  continue
c      write(*,'("Read ",i6," spectra, length ",i6)') nsp,nwave
      write(*,'("Read ",i6," spectra")') nsp

      return

 666  continue
      stop
      return
      end
