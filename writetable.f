
c read columns from a fits bintable and write as ascii

      program writetable

c max number of columns to read/write
      parameter (maxcol=15)
c max length of column
      parameter (maxlen=200000)
      parameter (maxaxis=8)

      character bname*80,oname*80,cline*80,colname*80
      character comment*80
      real datarr(maxlen,maxcol)
      real rdata(maxlen)
      integer idata(maxlen)
      integer laxes(maxaxis)
      integer icoltype(maxcol),nelem(maxcol)
      
c I despise equivalence but this may be the right way to do it?
c      equivalence(idata(1),rdata(1))

      nmax=maxlen
      ierr=0

 100  write(*,'("bintable file: ",$)')
      read(*,'(a)') bname
      
      call ftgiou(im,ierr)
      call doerr(ierr)
      call ftopen(im,bname,0,block,ierr)
      if (ierr .ne. 0) then
         call doerr(ierr)
         write(*,'("Error opening ",a)') bname
         go to 100
      end if

      write(*,'("output file: ",$)')
      read(*,'(a)') oname
      open(2,file=oname,status='unknown')

      ihdu=2

      nfield=1

 200  continue
      write(*,'("field name to get [end loop]: ",$)')
      read(*,'(a)') cline
      if (cline(1:3) .eq. '   ') go to 300

c move to requested HDU
      call ftmahd(im,ihdu,ihdutype,ierr)
      call doerr(ierr)
      if(ihdutype .ne. 2) then
         write(*,*) "Warning - HDU ",ihdu," is not a binary table"
         go to 200
      else
         write(*,*) "HDU ",ihdu," type ",ihdutype
      end if

c check table size with usual header keywords
      call ftgkyj(im,'NAXIS1',laxis1,comment,ierr)
      call ftgkyj(im,'NAXIS2',laxis2,comment,ierr)
      write(*,*) "table is ",laxis1," bytes, ",laxis2," entries"

c find number of column      
c     call ftgcnn(im,.false.,cline,colname,ifield,ierr)
c find datatype of column
c      call ftgtcl(im,ifield,itype,nrepeat,nwidth,ierr)

c don't check datatype yet
c      call readfield(im,ihdu,cline,itype,isize,naxis,laxes,rdata)
      call readfield(im,ihdu,cline,itype,isize,naxis,laxes,
     $     datarr(1,nfield))
      if (naxis .gt. 1) then
         write(*,*) "Warning, naxis >  1, naxis=",naxis
      end if
c      nelem(nfield) = laxes(1)
      nelem(nfield) = isize
      icoltype(nfield) = itype
      write(*,*) nfield, icoltype(nfield),nelem(nfield)
c      do i=1,nelem(nfield)
c         datarr(i,nfield) = rdata(i)
c      end do

      nfield=nfield+1
      go to 200

 300  continue
      nfield=nfield-1

      nmaxel = 0
      write(*,'("nelements: ",$)')
      do i=1,nfield
         nmaxel=max(nmaxel,nelem(i))
         write(*,'(1x,i7,$)') nelem(i)
      end do
      write(*,*)

c this is necessary if the table is screwed up (no TDIM)
      if (nmaxel .lt. laxis2) then
         nmaxel = laxis2
         write(*,*) "setting #elem = ",laxis2
      end if

      do iel=1,nmaxel
         do jcol=1,nfield
            if (icoltype(jcol) .eq. 41) then
c 32 bit integer
               write(2,'(1x,i11,$)') datarr(iel,jcol)
            else if (icoltype(jcol) .eq. 16) then
c ascii - this won't work right
               write(2,'(1x,a11,$)') datarr(iel,jcol)
            else if (icoltype(jcol) .eq. 21) then
c 16 bit int - this may not work right
               write(2,'(1x,i11,$)') datarr(iel,jcol)
c            else if (icoltype(jcol) .eq. 42) then
c 32 bit real
            else
c default to real
               write(2,'(1x,1pe11.4,$)') datarr(iel,jcol)
            end if
         end do
         write(2,*)
      end do

      close(2)
      call ftclos(im,ierr)
      call ftfiou(im,ierr)
      call doerr(ierr)
      
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
