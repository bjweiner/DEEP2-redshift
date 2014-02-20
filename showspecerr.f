
c Plot a spectrum using sensible limits.  Try to trim out
c regions with zero or bad data.

      subroutine showspecerr(n,wave,spec,serr)

      real wave(n),spec(n),serr(n)

      include 'pcredshift.h'

      small = 1.e-4

      call findends(n,spec,imin,imax)
      wmin = wave(imin)
      wmax = wave(imax)
c avoid problems with invalid limits
      if (imin .ge. imax .or. wmin .ge. wmax) then
         wmin = wave(1)
         wmax = wave(n)
      end if

c      write(*,*) imin,imax,wmin,wmax

c find the flux max/min within the good window, but ignore a
c buffer of ibuff pixels because there is sometimes high/low garbage
c in the last few pixels

      smin = 1.e10
      smax= -1.e10
c      do i=1,n
      ibuff = 20
      do i=imin+ibuff,imax-ibuff
         if (spec(i) .gt. bad) then
            smin = min(smin,spec(i))
            smax = max(smax,spec(i))
         end if
      end do
c trim plot to -120,120 to plot for now
      if (smin .gt. 9.e9) smin = badset
      smin=max(smin,-120.)
c      smax=max(min(smax,120.),-120.)
c make sure plot includes 0
      smin=min(smin,0.)
      if (smin .gt. smax-0.001) then
         smin=-120.
         smax=120.
      end if

c      write(*,'("showspec limits:",2(1x,i5),4(1x,f9.2))') imin,imax,
c     $      wmin,wmax,smin,smax      
      call pgenv(wmin,wmax,smin,smax,0,0)
      call pgqci(ioldci)
      call pgsci(3)
      call pgline(n,wave,serr)
      call pgsci(ioldci)
      call pgline(n,wave,spec)

c      open(10,file='showspec.out',status='unknown')
c      do i=1,n
c         write(10,'(i6,2x,f9.2,2x,f9.2)') i,wave(i),spec(i)
c      end do
c      close(10)

      return
      end
