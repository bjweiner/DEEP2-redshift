
c  this routine is supposed to retrieve a field/column
c  from a FITS binary table, using fitsio

c  input
c   im     pointer to image
c   ihdu   number of hdu to access
c   fname  name of column wanted
c  output
c   itype  datatype
c   isize  total length of column
c   naxis  number of axes
c   laxes  array of axis lengths
c   fdata  data array

      subroutine readfield(im,ihdu,fname,itype,isize,naxis,laxes,fdata)

c maxsize is max size of the column to return, obtained from maxpix
c and maxrow, but it doesn't have to be 2-d.  maxaxis is the
c maximum number of dimensions
      parameter(maxpix=4100,maxrow=300,maxaxis=8)
      parameter(maxsize=maxpix*maxrow)
c max number of columns in the binary table
      parameter(maxcolumn=128)
      parameter(IFDEBUG=0)

      integer im,ihdu,itype,isize,naxis
      character fname*(*)
      integer laxes(maxaxis)
      real fdata(maxsize)
      double precision ddata(maxsize)
      integer idata(maxsize)

      integer nfields
      character ttypes(maxcolumn)*80,tforms(maxcolumn)*80
      character tunits(maxcolumn)*80,extenname*80
      character comment*32,keyname*12
      character colname*32,colform*32,coldim*32
      character cnum*3

c  this stuff is for ftbtcl which probably is unnecessary
      character dtypestr*1
      double precision tscal,tzero

      ierr=0
c max number of dimensions allowed in a data array
      nmaxaxis = maxaxis

c check number of HDUs in file
      call ftthdu(im,nhdutot,ierr)
      call doerr(ierr)
      if(ihdu .gt. nhdutot) then
         write(*,*) "Error, requested hdu ",ihdu," out of ",nhdutot
         itype=-1
         isize=0
         go to 666
      end if
      if(IFDEBUG .eq. 1) then
         write(*,*) "readfield: hdu ",ihdu," of ",nhdutot
      end if
      
c move to the requested HDU.  Note 1=primary, 2=first extension, etc.
      call ftmahd(im,ihdu,ihdutype,ierr)
      if (ierr .ne. 0) then
         write(*,*) "Error, couldn't move to HDU ",ihdu
      end if
      call doerr(ierr)
      if(ihdutype .ne. 2) then
         write(*,*) "Warning - HDU ",ihdu," is not a binary table"
         itype =-1
         isize=0
         go to 666
      end if

c get number of fields
c      call ftgkyj(im,'TFIELDS',nfields,comment,ierr)
      call ftgncl(im,nfields,ierr)
      call doerr(ierr)
      if(nfields .gt. maxcolumn) then
         write(*,*) "Warning - too many fields"
      else
c         write(*,*) "found ",nfields," fields"
      end if

cc This is unnecessary since ftgcnn exists
cc loop through fields looking for the field with the right name
c      ifield=1
c 100  continue
cc this is a hack to left-justify the field number i.e. TTYPEn
c      write(cnum,'(i3)') ifield
c      write(keyname,'(a,a)') 'TTYPE',cnum
c      call ftgkys(im,keyname,colname,comment,ierr)
c      call doerr(ierr)
cc is it not what we want?
c      if (colname .ne. fname) then
c         ifield = ifield +1
c         if (ifield .gt. nfields) then
cc we've run out of fields, bail
c            write(*,'(a,a)') "Couldnt find field ",fname
cc should return an error here, for the moment set itype = -1
c            itype = -1
c            isize = 0
c            go to 666
c         else
c            go to 100
c         end if
c      end if

c find the field with requested name, not case-sensitive
      call ftgcnn(im,.false.,fname,colname,ifield,ierr)
      if (ierr .eq. 237) then
         write(*,'(a,a)') "Warning - column name not unique: ",fname
      else if (ierr .eq. 219) then
         write(*,'(a,a)') "Error - column was not found: ",fname
c hack to return an error 
         itype = -1
         isize = 0
         call doerr(ierr)
         go to 666
      else
c         write(*,'(a,i4,1x,a)') "Found column: ",ifield,fname
      end if
      call doerr(ierr)

c Now we've supposedly found what we want

cc This should only be needed as a check
c      write(keyname,'(a,a)') 'TFORM',cnum
c      call ftgkys(imt,keyname,colform,comment,ierr)
c      call doerr(ierr)
c      write(keyname,'(a,a)') 'TDIM',cnum
c      call ftgkys(imt,keyname,coldim,comment,ierr)
c      call doerr(ierr)
c      write(*,'(a,1x,a,1x,a)') colname(1:12),colform(1:12),coldim(1:12)

c get dimensions of the column from TDIMn keyword
      call ftgtdm(im,ifield,nmaxaxis,naxis,laxes,ierr)
      call doerr(ierr)

c get datatype for column. 1=bit X, 11=byte B, 14=logical L, 16=ascii A, 
c 21=short int(TINT) I, 41=int(TLONG) J, 42=real(TFLOAT) E,82=TDOUBLE D
c nrepeat is the number of elements in the column and 
c width is the length of each element (only relevant for strings?)
      call ftgtcl(im,ifield,itype,nrepeat,nwidth,ierr)
c      write(*,'(a,i4,1x,a,1x,i3,1x,i5)') 
c     $     "#, col, type, size ",ifield,fname,itype,nrepeat

cc this is for binary tables only, returns type as a 1-char string in dtypestr
c      call ftbtcl(im,ifield,ttype,tunit,dtypestr,nrepeat,
c     $     tscal,tzero,tnull,tdisp,ierr
c      call doerr(ierr)

c ought to do some size checking here
      nelem = nrepeat
c      isize = nrepeat
      if (itype .eq. 1) then
         nbits = 1
      else if (itype .eq. 11 .or. itype .eq. 14 .or. itype .eq. 16) then
         nbits = 8
      else if (itype .eq. 21) then
         nbits = 16
      else if (itype .eq. 41 .or. itype .eq. 42) then
         nbits = 32
      else if (itype .eq. 82) then
         nbits = 64
      else
         nbits=32
         write (*,*) "Warning, unknown data type ",itype
      end if
      n32bits = int(real(nrepeat*nbits)/32)
      if (n32bits .gt. maxsize) then
         write(*,*) "Error, column is too long ",n32bits,maxsize
      end if
      isize = nrepeat
      if (nbits .lt. 32)
     $     write(*,*) "Warning, readfield found a <4-byte datatype, ",
     $        "the data may be garbled."


c now even though we declared the data array as real, we'll
c dump the data in as the correct type

c Revised - try to read at least integer, real and double correctly

c  I may not have the null-flag type right for bits, bytes, short ints

c      write(*,*) "got here readfield 1 ",ifield,nelem,ianynull
      if (itype .eq. 1) then
c bits???
         call ftgcvb(im,ifield,1,1,nelem,0,fdata,ianynull,ierr)
      else if (itype .eq. 11 .or. itype .eq. 14) then
c byte or logical
         call ftgcvb(im,ifield,1,1,nelem,0,fdata,ianynull,ierr)
      else if(itype .eq. 16) then
c ascii
         call ftgcvs(im,ifield,1,1,nelem,' ',fdata,ianynull,ierr)
      else if(itype .eq. 21) then
c short int
         call ftgcvi(im,ifield,1,1,nelem,0,fdata,ianynull,ierr)
      else if(itype .eq. 41) then
c int
         call ftgcvj(im,ifield,1,1,nelem,0,idata,ianynull,ierr)
      else if(itype .eq. 42) then
c real
         call ftgcve(im,ifield,1,1,nelem,0.0,fdata,ianynull,ierr)
      else if(itype .eq. 82) then
c double
         call ftgcvd(im,ifield,1,1,nelem,0.d0,ddata,ianynull,ierr)
      else
c default to real
         call ftgcve(im,ifield,1,1,nelem,0.0,fdata,ianynull,ierr)
      end if
      call doerr(ierr)
c      write(*,*) "got here readfield 2 ",ifield,nelem,ianynull

c Always return an array of reals because the calling
c routine will expect that; convert integers and doubles.
      if (itype .eq. 41) then
         do i=1,nelem
            fdata(i) = real(idata(i))
            itype = 42
         end do
      else if (itype .eq. 82) then
         do i=1,nelem
            fdata(i) = real(ddata(i))
            itype = 42
         end do
      end if

 666  continue
      return

      end
