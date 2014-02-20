cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The problem with this routine is that it clears the error
c  flag rather than returning it ...

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

