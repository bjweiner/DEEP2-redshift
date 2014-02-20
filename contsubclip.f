
c do continuum subtraction with a moving average after sigma clipping

      subroutine contsubavgclip(spec,nsp,dw,wrad)
c spec is the array, nsp the dimension,
c dw the array spacing in wavelength or log w.l.
c wrad the radius of the window in wl or log wl space

      parameter(NSIGCLIP=3.0)
      parameter(IALGOR=1)

      real spec(nsp)
      real specout(nsp)
      real swork(nsp)

      include 'pcredshift.h'

      nrad = int(wrad/dw)
      if (nrad .lt. 1) nrad=1

      ngood1 = 0
      sum = 0.0
      sumsq=0.0
      do i=1,nsp
         if(spec(nsp) .gt. bad) then
            ngood1=ngood1+1
            sum = sum + spec(nsp)
            sumsq = sumsq + spec(nsp)*spec(nsp)
         end if
      end do
      if (ngood1 .lt. 1) then
         return
      end if
      avg = sum/ngood1
      rms = sqrt(sumsq/ngood1 - avg*avg)
      clipmin = avg - NSIGCLIP*rms
      clipmax = avg + NSIGCLIP*rms

c choose the way of doing it
      if (IALGOR .eq. 1) then
         go to 100
c      else if (IALGOR .eq. 2) then
      else
         go to 200
      end if

c  try 1 this is the direct but inefficient way

 100  continue

      do i = 1,nsp
         jmin = max(1,min(nsp,i-nrad))
         jmax = max(1,min(nsp,i+nrad))
         np = 0
         sum = 0.0
         do j=jmin,jmax
            if (spec(j) .gt. bad .and. spec(j) .gt. clipmin
     $           .and. spec(j) .lt. clipmax) then
               np = np+1
               sum = sum+spec(j)
            end if
         end do
         specout(i) = sum/np
      end do

      go to 600

c try 2

 200  continue

c i is the central point in the spectral arrray
c jmin,jmax are the indexes of the window ends in the spectral array
c npmin, npmax are the window ends in the work array
      i=1
      jmin = 1
      jmax = max(1,min(nsp,i+nrad))
      npmin = 1
      npmax = 0
c compute the mean in the window around the first point
      sum = 0.0
      np = 0
      do j=jmin,jmax
         if(spec(j) .gt. bad .and. spec(j) .gt. clipmin
     $        .and. spec(j) .lt. clipmax) then
            npmax=npmax+1
c            np = np + 1
            sum = sum + spec(j)
         end if
      end do

 220  continue
      np = npmax-npmin+1
      specout(i) = sum/np
      
c increment the point under consideration
      i=i+1
c if we've reached the end, bail
      if(i .gt. nsp) go to 230
c move the window
      jminold = jmin
      jmaxold = jmax
      jmin = max(1,min(nsp,i-nrad)
      jmax = max(1,min(nsp,i+nrad)
c  if we are eliminating a point at the low end
      if (jmin .gt. jminold .and. spec(jminold) .gt. bad
     $     .and. spec(jminold) .gt. clipmin
     $        .and. spec(jminold) .lt. clipmax) then
         npmin=npmin+1
c         np = np-1
         sum = sum - spec(jminold)
      end if
c  if we are adding a point at the high end
      if (jmax .gt. jmaxold .and. spec(jmax) .gt. bad
     $     .and. spec(jmax) .gt. clipmin 
     $     .and. spec(jmax) .lt. clipmax) then
         npmax = npmax+1
c         np = np+1
         sum = sum + spec(jmax)
      end if

c  go back to compute the mean
      go to 220

 230  continue
      go to 600

 600  continue

c  do actual subtraction
      do i=1,nsp
         if (spec(i) .gt. bad) spec(i) = spec(i) - specout(i)
      end do

      return
      end


