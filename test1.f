
c test various routines

      program test1

      include 'specarray.h'
      common /specarray/ specarr(NSPECMAX,NWAVEMAX)

      real spec1(NWAVEMAX),spec2(NWAVEMAX),spec3(NWAVEMAX)
      real wave1(NWAVEMAX),wave2(NWAVEMAX),wave3(NWAVEMAX)

      include 'pcredshift.h'

      call pgbeg(0,'?',2,2)
      call pgscf(2)
      call pgsch(1.5)

      call readspec(nspec,nwave)

      do i=1,nwave
         spec1(i) = specarr(1,i)
      end do

      nw = nwave
      nwlog = nwave
      
      w0=6001.28
      dw=1.28
      w0log = log10(6000.)
      dwlog = 1.e-4
      do i=1,nw
         wave1(i) = w0+i*dw
      end do
      do i=1,nwlog
         wave2(i) = w0log+i*dwlog
      end do

c do a plot
      call showspec(nw,wave1,spec1)

      call logrebin(spec1,nw,w0,dw,nwlog,w0log,dwlog,spec2)

      call showspec(nwlog,wave2,spec2)

c      call deshiftz(spec2,nwlog,dwlog,0.43)
      contwrad = 100.*dwlog
      call contsubmed(spec2,nwlog,dwlog,contwrad)

      call showspec(nwlog,wave2,spec2)

      end

