
c get and print skysigmas for a list of DEEP 2 1-d spectra

      program getskysigma

      parameter(maxhdu=12)
      parameter (IFOPTIMAL=0)

      integer indhdu(maxhdu)
      real sig(maxhdu)
      character cobjno*80

      character fname*160,sname*200
      
      ierr=0
      ifopt=IFOPTIMAL

c HDU's to get spectra from.  Blue and red chips are separate
c In typical file, HDU 1 is common, HDUs 2 and 3 are boxcar
c  4 and 5 are optimal
      nhdu = 2
c ifopt switches boxcar/optimal
      if (ifopt .ne. 0) then
         indhdu(1) = 4
         indhdu(2) = 5
      else
         indhdu(1) = 2
         indhdu(2) = 3
      end if

 100  continue
      write(*,'("File with list of 1-d spectra: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') stop
      open(2,file=fname,status='old',err=100)
      open(3,file='getskysigma.out',status='unknown')

 200  continue
      read(2,'(a)',err=666,end=666) sname
      
      call ftgiou(im,ierr)
      call doerr(ierr)
      call ftopen(im,sname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         write(*,'("Error opening ",a)') sname
         do i=1,nhdu
            sig(i)=0.0
         end do
         go to 300
      end if
      
      call readobjno(im,indhdu(1),cobjno)
      do ihdu=1,nhdu
         call readskysigma(im,indhdu(ihdu),sig(ihdu))
      end do
      call ftclos(im,ierr)
      call doerr(ierr)
      call ftfiou(im,ierr)
      call doerr(ierr)

 300  continue
c      write(3,*) (sig(i),i=1,nhdu),sname
      write(3,1000) (sig(i),i=1,2),cobjno,sname
 1000 format(2(f7.3,2x),a16,a120)
      
      go to 200

 666  continue
      close(2)
      close(3)

      end

