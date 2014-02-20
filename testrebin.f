
c test rebinning subroutines

      program testrebin

      parameter (IFCONST=0)
      parameter (nmax=1000)

      real w(nmax),warb(nmax),wlog(nmax)
      real wlin(nmax),wlin2(nmax)
      real flux1(nmax),fluxlog(nmax),fluxlin(nmax)
      real flux2(nmax),flux2log(nmax),flux2lin(nmax)
      real ferr1(nmax),ferrlog(nmax),ferrlin(nmax)
      character toplabel*100

      integer ilines(10)
      real rlinekern(11)

      n=300
      nlog=300
      nlin=400
c      n=100
c      nlog=100
c      nlin=100
      narb=n
      w0=6500
      dw=0.5

c make the log scale cover slightly more than the linear scale
      w0log= log10(w0-10.)
      wlogmax = log10(w0+dw*n+10.)
      dwlog = (wlogmax-w0log)/nlog

      write(*,'(a,4(2x,f11.5))') "w0,dw,w0log,dwlog",w0,dw,w0log,dwlog

      nlines=4
      do j=1,nlines
         ilines(j) = 61*j
      end do
c construct gaussian kernel for lines
      iradkern=4
      nlinekern=2*iradkern+1
      gsigma = 1.1
      do j=1,nlinekern
         rlinekern(j) = 10.* exp(-(iradkern+1-j)**2/2./gsigma**2)
      end do

      do i=1,n
c standard linear wl scale
         w(i) = w0+i*dw
c wl not quite linear
         warb(i) = w(i) + 0.02*(i-1)
      end do
      do i=1,nlog
         wlog(i) = w0log + dwlog*i
      end do

      if (IFCONST .eq. 1) then
c make it all flat
         do i=1,n/2
            flux1(i) =15.
         end do
         do i=n/2+1,n
            flux1(i) =25.
         end do
      else
c constant part
         do i=1,50
            flux1(i) = 15.
         end do
c quadratic increase
         do i=51,150
            flux1(i) = flux1(50) + 0.001*(i-50)**2
         end do
c lin decrease
         do i=151,250
            flux1(i) = flux1(150) - 0.1*(i-150)
         end do
c zero patch
         do i=251,260
            flux1(i)=0.
         end do
c constant w/offset
         do i=261,270
            flux1(i) = flux1(250) + 5.
         end do
c lin increase
         do i=271,300
            flux1(i) = flux1(270) + 0.1*(i-270)
         end do
c add a few lines
         do j=1,nlines
            do ikern=1,nlinekern
               i=ilines(j)+ikern-iradkern-1
               flux1(i) = flux1(i) + rlinekern(ikern)
            end do
         end do
      end if

      rnsq=16.0
      do i=1,n
         ferr1(i) = sqrt(flux1(i) + rnsq)
      end do

      call pgbeg(0,'?',2,2)
      call pgscf(2)
      call pgsch(1.5)

      call showspecerr(n,w,flux1,ferr1)
      call pglabel('linear wavelength','counts','linear original')

      call logrebinerr(flux1,ferr1,n,w0,dw,nlog,w0log,dwlog,
     $     fluxlog,ferrlog)
      call logrebin(flux1,n,w0,dw,nlog,w0log,dwlog,
     $     flux2log)

      call showspecerr(nlog,wlog,fluxlog,ferrlog)
      write(toplabel,700) n,nlog,dw,dwlog
 700  format("log rebin, nlin= ",i4,", nlog= ",i4,
     $     ", dw= ",f4.2,", dwlog= ",f7.5)
      call pglabel('log wavelength','rebinned counts',toplabel)
      call pgsci(2)
      call pgsls(2)
      call pgline(nlog,wlog,flux2log)
      call pgsci(1)
      call pgsls(1)
      
      call pgpage()
      call pgpage()

      call showspecerr(narb,warb,flux1,ferr1)
      call pglabel('non-linear wavelength','counts',
     $     'non-linear original')

      call linrebinerr(narb,warb,flux1,ferr1,nlin,w0lin,dwlin,
     $     fluxlin,ferrlin)
      call linrebin(narb,warb,flux1,nlin,w0lin2,dwlin2,flux2lin)

      do i=1,nlin
         wlin(i) = w0lin +i*dwlin
         wlin2(i) = w0lin2 + i*dwlin2
      end do
      call showspecerr(nlin,wlin,fluxlin,ferrlin)
      call pglabel('linearized wavelength','rebinned counts',
     $     'linearized rebin of non-linear original')
      call pgsci(2)
      call pgsls(2)
      call pgline(nlin,wlin2,flux2lin)
      call pgsci(1)
      call pgsls(1)

      call logrebinerr(fluxlin,ferrlin,nlin,w0lin,dwlin,nlog,w0log,
     $     dwlog,fluxlog,ferrlog)
      call logrebin(fluxlin,nlin,w0lin2,dwlin2,nlog,w0log,dwlog,
     $     flux2log)

      call showspecerr(nlog,wlog,fluxlog,ferrlog)
      write(toplabel,710) nlin,nlog,dwlin,dwlog
 710  format("log rebin, nlin= ",i4,", nlog= ",i4,
     $     ", dw= ",f4.2,", dwlog= ",f7.5)
      call pglabel('log wavelength','rebinned counts',toplabel)
      call pgsci(2)
      call pgsls(2)
      call pgline(nlog,wlog,flux2log)
      call pgsci(1)
      call pgsls(1)
     
      call pgend()

      end
