
c test getting fields by name out of a binary table

      program testtable

      parameter (maxpix=10000,maxrow=100)
      parameter (maxdim=maxrow*maxpix)
      parameter (maxfield=128)

      real fdata(maxdim)
      integer laxes(maxfield)

      character tname*120,oname*120

      ierr=0

 100  continue
      write(*,'("Name of table file [quit]: ",$)')
      read(*,'(a)') tname
      if (tname(1:3) .eq. '   ') go to 666
      call ftgiou(imt,ierr)
      call doerr(ierr)
      call ftopen(imt,tname,0,block,ierr)
      if (ierr .ne. 0 ) then
         call doerr(ierr)
         go to 100
      end if

 120  write(*,'("New image file: ",$)')
      read (*,'(a)') oname
      if (oname(1:3) .eq. '   ') go to 666
      call ftgiou(imo,ierr)
      call doerr(ierr)
      call ftinit(imo,oname,block,ierr)
      if (ierr .ne. 0 ) then
         call doerr(ierr)
         go to 120
      end if

c copy the primary header (CHDU)
c      call imhcpy(imt,imo,ierr)
c      call doerr(ierr)
      call copyhead(imt,imo,ierr)
      write(*,'("which HDU: ",$)')
      read (*,*) ihdu

c      call readfield(imt,ihdu,'flux',itype,isize,naxis,laxes,fdata)
      call readfield(imt,ihdu,'spec',itype,isize,naxis,laxes,fdata)
c      write (*,*) "did readfield"
      call writeimage(imo,itype,isize,naxis,laxes,fdata)
c make extension
      call readfield(imt,ihdu,'ivar',itype,isize,naxis,laxes,fdata)
      call writeimage(imo,itype,isize,naxis,laxes,fdata)
      call readfield(imt,ihdu,'lambda',itype,isize,naxis,laxes,fdata)
      call writeimage(imo,itype,isize,naxis,laxes,fdata)
c      call readfield(imt,ihdu,'lambda0',itype,isize,naxis,laxes,fdata)
c      call writeimage(imo,itype,isize,naxis,laxes,fdata)

      call ftclos(imt,ierr)
      call ftclos(imo,ierr)
      call doerr(ierr)
      call ftfiou(imt,ierr)
      call ftfiou(imo,ierr)
      call doerr(ierr)

c  loop back for another
      go to 100

 666  continue
      end 


      subroutine writeimage(imo,itype,isize,naxis,laxes,fdata)

      real fdata(isize)
      integer laxes(naxis)

c      write (*,*) "writeimage 1"

      if(itype .eq. 21) then
         call ftiimg(imo,16,naxis,laxes,ierr)
         call doerr(ierr)
         call ftppri(imo,0,1,isize,fdata,ierr)
      else if(itype .eq. 31 .or. itype .eq. 41) then
         call ftiimg(imo,32,naxis,laxes,ierr)
         call doerr(ierr)
         call ftpprj(imo,0,1,isize,fdata,ierr)
      else if (itype .eq. 82) then
         call ftiimg(imo,-64,naxis,laxes,ierr)
         call doerr(ierr)
         call ftpprd(imo,0,1,isize,fdata,ierr)
      else if(itype .eq. 42) then
         call ftiimg(imo,-32,naxis,laxes,ierr)
         call doerr(ierr)
         call ftppre(imo,0,1,isize,fdata,ierr)
      else
         call ftiimg(imo,-32,naxis,laxes,ierr)
         call doerr(ierr)
         call ftppre(imo,0,1,isize,fdata,ierr)
      end if
c      write (*,*) "writeimage 2"
      call doerr(ierr)
c      write (*,*) "writeimage 3"

      return
      end

cccccccccccccccccccccccccccccccccccccccc

      subroutine imhcpy( im1, im2, ierr)

      parameter(ifskipbasic=0)
      integer im1, im2
      character*80 card, comment
      integer nexist,nadd,ntmp

      ierr = 0
      call ftgkyj(im2, 'naxis', naxis, comment, ierr)
      call doerr( ierr)

      if ( naxis.eq.1) js=5
      if ( naxis.eq.2) js=6

      ierr = 0
      do j=js,js+4
      call ftdrec(im2, 6, ierr)
      call doerr( ierr)
      end do
  50      continue

c      write(6,*) 'copying header'

c  instead of looping to end, cleaner to use
c  ftghsp to find # of existing keywords
      call ftghsp( im1,nexist,ntmp,ierr)
      write(*,*) "nexist1, nadd1: ",nexist,ntmp
      call doerr(ierr)
      call ftghsp( im2,ntmp,nadd,ierr)
      write(*,*) "nexist2, nadd2: ",ntmp,nadd
      call doerr(ierr)
      if (nadd .ne. -1 .and. nexist .gt. ntmp+nadd) then
         write(*,*) 'fitsio: not enough space to copy header? '
c  do something about it?
      end if
      do i=1,nexist
         call ftgrec( im1, i, card, ierr)
         call doerr(ierr)
         if ( ierr.ne.0) then
c  assume we've reached end of header
c            write (*,*) '  fitsio/imhcpy: error getting record',i,card(1:8)
            goto 10
         end if
         if (ifskipbasic .ne. 0) then
          if ( card(1:8).eq.'SIMPLE  ') goto 5
          if ( card(1:8).eq.'BITPIX  ') goto 5
          if ( card(1:8).eq.'NAXIS   ') goto 5
          if ( card(1:8).eq.'NAXIS1  ') goto 5
          if ( card(1:8).eq.'NAXIS2  ') goto 5
          if ( card(1:8).eq.'BSCALE  ') goto 5
          if ( card(1:8).eq.'BZERO   ') goto 5
          if ( card(1:8).eq.'EXTEND  ') goto 5
          if ( card(1:8).eq.'UNSIGN  ') goto 5
         end if
c  test
         if (ichar(card(80:80)).eq.10) card(80:80)=' '
c         write(6,'(i2,i3,2x,a72)') ierr,i,card(1:72)
         call ftprec( im2, card, ierr)
c  Here it might be possible to run out of space in the header
c  of the second image.         
c         if ( ierr.ne.0) then
c        write (*,*) '  fitsio/imhcpy: error putting record',i,card(1:8)
c         end if
         call doerr( ierr)
 5       continue
      end do

 10   continue
      ierr = 0
      call doerr( ierr)

      return
      end


      subroutine doerr( status)

      integer status
      character*80 error

      if ( status.ne.0) then
      call ftgerr( status, error)
      write(6,25) error
c  25      format(20x,a30)
  25      format(2x,'fitsio: ',a30)
      end if

      status = 0

      call ftcmsg()

      return
      end

cccccccccccc
      subroutine copyhead(im1,im2,ierr)

      integer im1, im2, ierr
      character card*80, comment*80, keyw*8

      call ftghsp(im1,nexist1,nadd1,ierr)
      call ftghsp(im2,nexist2,nadd2,ierr)

      do i=1,nexist1
         call ftgrec(im1,i,card,ierr)
         call doerr(ierr)
         if (ierr .ne. 0) then
c  assume we've reached end of header
            go to 100
         end if
         keyw = card(1:8)
c test to avoid overwriting, SIMPLE,BITPIX, etc?
         
         call ftprec(im2,card,ierr)
         call doerr(ierr)
      end do

 100  continue
      return
      end

      
         


