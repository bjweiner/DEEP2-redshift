
c test the linrebin subroutine

      program testlinear

      parameter(maxp=12000)
      real wave(maxp),flux(maxp),fnew(maxp)
      character fname*80,oname*80

 100  continue
      write(*,'("file with wave and flux [quit]: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') stop
      open(2,file=fname,status='old',err=100)

 120  continue
      write(*,'("output file: ",$)')
      read(*,'(a)') oname
      if (fname(1:3) .eq. '   ') then
         open(3,file='testlin.out',status='unknown',err=120)
      else
         open(3,file=oname,status='new',err=120)
      end if

      np=1
 200  continue
      read(2,*,err=250,end=250) wave(np),flux(np)
      np=np+1
      go to 200
 250  continue
      np = np-1
      close(2)
      write(*,*) "read ",np," points"
      write(*,*) "wave range ",wave(1),wave(np)

      call linrebin(np,wave,flux,np,w0,dw,fnew)

      write(*,*) "new w0, dw ",w0,dw

      do i=1,np
         wi = w0 + i*dw
         write(3,'(f15.6,3x,f15.6)') wi,fnew(i)
      end do
      close(3)

      end


         

      
