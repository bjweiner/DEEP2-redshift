
c Given arrays of flux and error, calculate the median SNR/pixel.
c Could add blanking values that are way off for whatever reason.


      subroutine calcsnrmed(np,spec,espec,snrmed)

      real spec(np),espec(np)
      real snr(np)

      external select

      do i=1,np
         snr(i) = spec(i)/espec(i)
      end do

      itmp = np/2
      snrmed = select(itmp,np,snr)

      return
      end

