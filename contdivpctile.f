      
c do continuum subtraction or division by taking a high percentile
c value in a (large) moving window
c for example, take the 90%ile value.  we don't want to take the
c maximum because some pixels will be bad or emission lines.
c this was inspired by Karl's method of taking the biweight of
c the top 1/3rd of the data - if you take the median of the
c top 1/3rd, it's the 83%ile by def.

c can be more efficient, as Karl did, by only taking the 
c %ile every Nsamp pixels and interpolating, rather than
c recomputing at every pixel.

c wrad is the radius of the subtraction window, in wavelength
c or log(1+z), dw is the x-spacing
c do the subtraction in place, or not
c      subroutine contsubmed(spec,nsp,dwlog,wrad,specout)
      subroutine contdivpctile(spec,nsp,serr,dw,wrad)

c divide or subtract?
      parameter(IFDIVIDE=1)
c how many pixels apart the places to compute the %ile are at
      parameter(NSAMP=200)
c how many pixels to use to compute median.  This is redundant
c with having wrad as an input param.
c      parameter(NWIND=201)
c nvals is the max number of places to compute the %ile at, so
c nmaxvals*nsamp must be greater than the dimension of the input array, nsp
      parameter(NMAXVALS=250)
      parameter(IFDEBUG=0)

      real spec(nsp)
      real serr(nsp)
c      real specout(nsp)
      real swork(nsp)
      real tplot(nsp)
      real cvals(NMAXVALS),rivals(NMAXVALS)
      integer ivals(NMAXVALS)

c      external function select
c      external function selip
      external select
      external selip

      include 'pcredshift.h'

      pctfrac = 0.9

      nrad = int(wrad/dw)
      if (nrad .lt. 1) nrad=1
      if (IFDEBUG .ge. 2) write(*,*) "contdivpct: nrad = ",nrad

c indexes of central locations of windows
      ivals(1) = nrad+1
      nvals = int( (nsp-ivals(1))/real(NSAMP) ) + 1
      do j=2,nvals
         ivals(j) = ivals(1) + (j-1)*NSAMP
         if (IFDEBUG .ge. 2) write(*,*) "cdp: ",j,ivals(j)
      end do

c loop over locations; load pixels into window using swork, compute the
c 90th or whatever %ile with select, save in cvals
      do j=1,nvals
         imin = max(ivals(j)-nrad,1)
         imax = min(ivals(j)+nrad,nsp)
c         ntodo = imax-imin+1
         np = 0
         do i=imin,imax
            if (spec(i) .gt. bad) then
               np = np+1
               swork(np) = spec(i)
            end if
         end do
         if (np .gt. 0) then
            itmp = int(pctfrac*np)
            cvals(j) = select(itmp,np,swork)
            if (IFDEBUG .ge. 2) then
               write(*,*) "cdp: ",j,ivals(j),np,itmp,cvals(j)
            end if
         else
            cvals(j) = badset
         end if
      end do
      if (IFDEBUG .ge. 2) then
         do j=1,nvals
            write(*,*) "contdivpctile: ",j,ivals(j),cvals(j)
         end do
      end if

c go through locations and fix any that are equal to badset by
c interpolating between good values.  This shouldn't really happen.
      if (cvals(1) .lt. bad) cvals(1) = 0.0
      if (cvals(nvals) .lt. bad) cvals(nvals) = 0.0
      do j=2,nvals-1
         if (cvals(j) .lt. bad) then
            jprev = j-1
            jnext = j+1
c loop up to find next good value
 110        continue
            if (cvals(jnext) .lt. bad) then
               jnext = jnext+1
               go to 110
            end if
            if (IFDEBUG .ge. 2) then
               write(*,*) "contdivpctile fix: ",j,ivals(j),cvals(j),
     $              jprev,jnext,cvals(jprev),cvals(jnext)
            end if
            cvals(j) = cvals(jprev) +
     $        (j-jprev)/real(jnext-jprev) * (cvals(jnext)-cvals(jprev))
         end if
      end do

c now use swork at the full length of the spec array; interpolate the 
c cvals points onto it to create a continuum array.  Points lower than
c the first or higher than the last are special cases
      do i=1,ivals(1)-1
         swork(i) = cvals(1) + 
     $        (i-ivals(1))/real(ivals(2)-ivals(1)) * (cvals(2)-cvals(1))
      end do
      do j=1,nvals-1
         swork(ivals(j)) = cvals(j)
         do i=ivals(j),ivals(j+1)-1
            swork(i) = cvals(j) + 
     $    (i-ivals(j))/real(ivals(j+1)-ivals(j)) * (cvals(j+1)-cvals(j))
         end do
      end do
      do i=ivals(nvals),nsp
         swork(i) = cvals(nvals) + 
     $        (i-ivals(nvals-1))/real(ivals(nvals)-ivals(nvals-1)) * 
     $        (cvals(nvals)-cvals(nvals-1))
      end do

      if (IFDEBUG .ge. 1) then
         smax = -1.e6
         smin = 1.e6
         do i=1,nsp
            tplot(i) = real(i)
            if (spec(i) .gt. bad) then
               smax = max(smax,spec(i))
               smin = min(smin,spec(i))
            end if
         end do
         if (smin .gt. 0.0) smin=0.0
         call pgenv(0.,real(nsp),smin,smax,0,1)
         call pglabel("index","flux","spectrum and continuum")
         call pgline(nsp,tplot,spec)
         call pgqci(indexc)
         call pgsci(3)
         call pgline(nsp,tplot,swork)
         call pgsci(2)
         do j=1,nvals
            rivals(j) = real(ivals(j))
         end do
         call pgpt(nvals,rivals,cvals,17)
         call pgsci(indexc)
c         call showspec(nsp,tplot,spec)
c         call showspec(nsp,tplot,swork)
      end if

c now we either divide or subtract the continuum in place
c if we do division, also have to divide the error.
      if (IFDIVIDE .eq. 0) then
         do i=1,nsp
            if (spec(i) .gt. bad) spec(i) = spec(i) - swork(i)
         end do
      else
         do i=1,nsp
            if (spec(i) .gt. bad) then
               if (swork(i) .gt. 1.e-4) then
                  spec(i) = spec(i) / swork(i)
                  serr(i) = serr(i) / swork(i)
               else
                  spec(i) = spec(i) / 1.e-4
                  serr(i) = serr(i) / 1.e-4
               end if
            end if
         end do
      end if

      return
      end

         
