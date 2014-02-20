
c  given a list of spec1d files, retrieve the
c  y position from which the spectrum was extracted
c  and write a file

      program getypos

      parameter (maxaxis=10)
      integer laxes(maxaxis)
      real ypos,fwhm

      character fname*200,iname*200

      ierr=0
      
 100  continue
 110  write(*,'("List of files: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') stop
      open(2,file=fname,status='old',err=110)
      
      open(3,file='getypos.out',status='unknown')
      open(4,file='getypos.out2',status='unknown')
      ihdu = 2
      i=1
      ifail=0

 200  continue
      read(2,'(a)',err=500,end=500) iname
      
      call ftgiou(im,ierr)
      call doerr(ierr)
      call ftopen(im,iname,0,block,ierr)
      call doerr(ierr)
           if (ierr .ne. 0) then
         call doerr(ierr)
         ifer = 1
         write(*,'("Error opening ",a)') iname
         write(3,'(f7.2)') 0
         ifail=ifail+1
         go to 300
      end if

      call readfield(im,ihdu,'OBJPOS',itype,isize,naxis,laxes,ypos)
      call readfield(im,ihdu,'FWHM',itype,isize,naxis,laxes,fwhm)
      write (3,'(f7.2)') ypos
      write (4,'(f7.2,3x,f7.2)') ypos,fwhm
      call ftclos(im,ierr)
      call doerr(ierr)
      call ftfiou(im,ierr)
      call doerr(ierr)

      i=i+1

 300  continue
      go to 200

 500  continue
      i=i-1
      write(*,*) "Read ",i," files OK, ",ifail," failures"
      close(3)
      close(4)

      end
