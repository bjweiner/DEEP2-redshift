
c use my continuum subtraction routine on an ascii file,
c such as the one I made from Sloan eignevectors

      program asciicontsub

      parameter(maxwave=10000, maxvect=50)

c      parameter(IFPOLYSUB=1)

      real wlog(maxwave)
      real vects(maxwave,maxvect)
      real serr(maxwave),sspec(maxwave),scont(maxwave)
      character fname*80,oname*80

      include 'pcredshift.h'

      call pgbeg(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.2)      

      write(*,'("Fit method: 0-median, 1-poly, 2-Legendre: ",$)')
      read(*,*) isubtype
      
 100  write(*,'("File with vectors: ",$)')
      read(*,'(a)') fname
      open(2,file=fname,status='old',err=100)

 120  write(*,'("Output file: ",$)')
      read(*,'(a)') oname
      open(3,file=oname,status='new',err=120)
      
 130  write(*,'("Number of vectors: ",$)')
      read(*,*) nv
      if(nv .gt. maxvect) then
         write(*,*) "Maximum # of vectors is ",maxvect
         go to 130
      end if

      nwave = 1
 200  continue
      read(2,*,err=300,end=300) wlog(nwave),(vects(nwave,j),j=1,nv)
      nwave=nwave+1
      if (nwave .gt. maxwave) then
         write(*,*) "File length exceeded max # ",maxwave
         stop
      end if
      go to 200

 300  continue
      close(2)
      nwave=nwave-1
      write(*,'("Read ",i6," lines")') nwave
      
      dwlog = (wlog(nwave)-wlog(1)) / (nwave-1.0)
      w0log = wlog(1)-dwlog
      write(*,'("for evenly spaced scale, w0, dw = ",f11.5,2x,f9.5)')
     $     w0log,dwlog
      if (isubtype .eq. 0) then
         write(*,'("radius for median subtraction in pixels: ",$)')
         read(*,*) radpix
         contrad = radpix*dwlog
      else
         write(*,'("Polynomial order: ",$)') 
         read(*,*) iorder
      end if

      do i=1,nwave
         serr(i) = 1.0
      end do

      do j=1,nv
c plot boundaries
         wlmin=wlog(1)
         wlmax=wlog(nwave)
         ngood=0
         sum=0.0
         sumsq=0.0
         do i=1,nwave
            if (vects(i,j) .gt. bad) then
               ngood=ngood+1
               sspec(i) = vects(i,j)
               sum=sum+vects(i,j)
               sumsq=sumsq+vects(i,j)**2
            end if
         end do
         smean = sum/ngood
         srms = sqrt((sumsq - sum*sum/ngood)/ngood)
c subtracted spec is about zero, so make minimum of -rms
         s1 = min(-rms,smean-3*srms)
         s2 = max(rms,smean+3*srms)
         call pgenv(wlmin,wlmax,s1,s2,0,0)
         call pglabel("wavelength","intensity",
     $        "spectrum, continuum, subtracted")
         call pgline(nwave,wlog,sspec)

         if (isubtype .eq. 0) then
            call contsubmed(vects(1,j),nwave,dwlog,contrad)
         else if (isubtype .eq. 1) then
            call contsubpoly(wlog,vects(1,j),serr,nwave,iorder)
         else
            call contsubleg(wlog,vects(1,j),serr,nwave,iorder)
         end if
c reverse engineer continuum fit
         do i=1,nwave
            scont(i) = sspec(i) - vects(i,j)
         end do
         call pgsci(2)
         call pgline(nwave,wlog,vects(1,j))
         call pgsci(4)
         call pgline(nwave,wlog,scont)
         call pgsci(1)
      end do

      do i=1,nwave
         write(3,1000) wlog(i)
c this do loop is needed since we don't know how many vectors to
c put in the format statement
         do j=1,nv
            write(3,1010) vects(i,j)
         end do
c write the newline
         write(3,*)
      end do
 1000 format(f7.5,$)
 1010 format(2x,1pe10.3,$)

      close(3)
      call pgend()

      end

      
