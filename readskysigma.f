
c read a sky sigma out of the header of a given HDU of a
c DEEP 2 fits table 1-d spectrum

      subroutine readskysigma(im,ihdu,sigma)

      integer im,ihdu
      real sigma
      character comment*72

      ierr=0
      call ftmahd(im,ihdu,ihdutype,ierr)
      call doerr(ierr)
      call ftgkye(im,'SKYSIGMA',sigma,comment,ierr)
      call doerr(ierr)

      return
      end
