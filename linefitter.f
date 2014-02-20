
c given lists of object numbers, spec1d files, and spectral features,
c fit lines.

      program linefitter

      include 'specarray.h'

      parameter (NMAXLINES=50)

      real wave(NWAVEMAX),spec(NWAVEMAX),espec(NWAVEMAX)

      real wrest(NMAXLINES)
      character fname*160,sname*200

      fitrad = 30.
      wmaxoffset = 10.
      contrad = 120.
      ifoptextract=0


 90  continue
      write(*,'("linelist: ",$)')
      read(*,'(a)') fname
      open(10,file=fname,status='old',err=100)
 100  continue
      write(*,'("list of object numbers: ",$)')
      read(*,'(a)') fname
      open(2,file=fname,status='old',err=100)
 110  continue
      write(*,'("list of 1-d spectra: ",$)')
      read(*,'(a)') fname
      open(3,file=fname,status='old',err=110)

c read lines from linelist

      nlines=1
 200  continue
      read(10,*,err=210,end=210) wrest(nlines)
      nlines=nlines+1
      go to 200
 210  continue
      nlines=nlines-1
      write(*,*) "Read ", nlines," wavelengths"
      close(10)

      nobj=1
 300  continue
c loop through objects
      read(2,*,err=666,end=666) iobjno
      read(2,'(a)',err=666,end=666) sname

      call read1dspec(sname,ifoptextract,np,wave,spec,espec,ifer)
      if (ifer .ne. 0) then
         write(*,*) "Error getting spectrum for ",iobjno
         go to 300
      end if

      do i=1,nlines
         


      end do

