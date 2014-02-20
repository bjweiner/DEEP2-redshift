
c  This is the same as contsubmed.f but it does the whole array
c  at once

c do continuum subtraction by median in a (large) moving window

c The clever way to this is to have the wavelength be the
c first (most rapidly varying) index in specarray, i.e.
c specarray(nwave,spec).  Then you can avoid copying it 
c in and out of the spec temp array.  Unfortunately I didn't
c set specarray up that way.  Maybe next time.

c wrad is the radius of the subtraction window, in wavelength
c or log(1+z), dw is the x-spacing
c do the subtraction in place, or not
c      subroutine contsubmed(spec,nsp,dwlog,wrad,specout)
c      subroutine contsubmed(spec,nsp,dw,wrad)
      subroutine contsubmedarray(specarray,nspec,nwave,dw,wrad)

      parameter(IALGOR=2)

      real specarray(nspec,nwave)
      real spec(nwave),specout(nwave)
      real swork(nwave)

c      external function select
c      external function selip
      external select
      external selip

      include 'pcredshift.h'

      nrad = int(wrad/dw)
      if (nrad .lt. 1) nrad=1

      do ii = 1,nspec
      do j=1,nwave
         spec(j) = specarray(ii,j)
      end do

c choose the way of doing it
      if (IALGOR .eq. 1) then
         go to 100
c      else if (IALGOR .eq. 2) then
      else
         go to 200
      end if

cTry 1
c this must be an incredibly inefficient way of convolving
c with a median filter
 100  continue
      do i=1,nwave
         jmin = max(1,min(nwave,i-nrad))
         jmax = max(1,min(nwave,i+nrad))
         np = 0
         do j=jmin,jmax
            if (spec(j) .gt. bad) then
               np = np+1
               swork(np) = spec(j)
            end if
         end do
c         itmp = int(np/2)
         itmp = np/2
         if (np .gt. 0) then
c sloppy median
            rmed = select(itmp+1,np,swork)
         else
            rmed = badset
         end if
         specout(i) = rmed
      end do

      go to 600

c Try 2

 200  continue

c i is the central point in the spectral arrray
c jmin,jmax are the indexes of the window ends in the spectral array
c npmin, npmax are the window ends in the work array
      i=1
      jmin = 1
      jmax = max(1,min(nwave,i+nrad))
      npmin = 1
      npmax = 0
c load the array with the window around the first point
      do j=jmin,jmax
         if(spec(j) .gt. bad) then
            npmax=npmax+1
            swork(npmax) = spec(j)
         end if
      end do

 220  continue
c sloppy median. np is the number of points in the work array to median.
      np = npmax-npmin+1
      k = np/2+1
      if (np .lt. 1) then
         specout(i) = badset
      else
         specout(i) = selip(np/2+1,np,swork(npmin))
      end if

c increment the point under consideration
      i=i+1
c if we've reached the end, bail
      if(i .gt. nwave) go to 230
c move the window
      jminold = jmin
      jmaxold = jmax
      jmin = max(1,min(nwave,i-nrad))
      jmax = max(1,min(nwave,i+nrad))
c  if we are eliminating a point at the low end
      if (jmin .gt. jminold .and. spec(jminold) .gt. bad) then
         npmin=npmin+1
      end if
c  if we are adding a point at the high end
      if (jmax .gt. jmaxold .and. spec(jmax) .gt. bad) then
         npmax = npmax+1
         swork(npmax) = spec(jmax)
      end if
c go back to compute the median
      go to 220

c get here when done
 230  continue
      go to 600
      
c here there's space for more algorithms.         

 300  continue

 400  continue

 600  continue

cc  do actual subtraction
c      do i=1,nwave
c         if (spec(i) .gt. bad) spec(i) = spec(i) - specout(i)
c      end do
      do i=1,nwave
         if (specarray(ii,i) .gt. bad) then
            specarray(ii,i) = specarray(ii,i) - specout(i)
         end if
      end do

c end the do ii loop
      end do

      return
      end

