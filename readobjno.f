
c read object number out of the header of a given HDU of a
c DEEP 2 fits table 1-d spectrum

      subroutine readobjno(im,ihdu,cobjno)

      integer im,ihdu
      character cobjno*80,comment*72

      ierr=0
      call ftmahd(im,ihdu,ihdutype,ierr)
      call doerr(ierr)
      call ftgkys(im,'OBJNO',cobjno,comment,ierr)
      call doerr(ierr)

      return
      end
