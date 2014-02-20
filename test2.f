
      program test2

      character cline*80,fname*80,chfield*80
      character concat1*60,concat2*60

 100  continue
      write(*,'("file with fields [quit]: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') stop
      open(2,file=fname,status='old',err=100)

      i=0
 200  continue
      read(2,'(a)',err=250,end=250) cline
      read(cline,'(a)') chfield
      i=i+1
      ilen = index(chfield,' ')-1
      write(*,'(i4,2x,i3,2x,a)') i,ilen,chfield
      concat1 = chfield // ".sufx1"
      concat2 = chfield(1:ilen) // ".sufx2"
      write(*,'(a)') concat1
      write(*,'(a)') concat2

      go to 200
      
 250  close(2)
      end
