MODULE SUBS

!       Subroutines in this module:

!       xsort  -> Sorts multicolumn arrays by direct insertion.
!       qsort  -> Sorts multicolumn arrays by useing the QUICKSORT routine
!                 taken from Numerical Recipies p.236
!       xsize  -> Determines the number of columns and rows in an ASCII file.
!       readit -> Reads an ASCII data file into an array.
!       writeit-> Writes an array into a file in ASCII format. 

!       Module compiles with: f90 -c subs.f90

        implicit none

        contains

!----------------------------------------------------------------------------

subroutine qsort(cols, rows, sortcol, data)

!       This subroutine will sort a multicolumn array into ascending order.
!       The column to be sorted is specified by the sortcol variable and the
!       number of columns and rows by the columns and points variables respectively.
!       Calculations are done in double precission. Routine taken from 
!       Numerical Recipies p.323 (8.2 Quicksorting )and adapted for a 
!       multicolumn array.

!       Variable definitions.
!       *********************************************************************
        integer :: n, m, nstack, i, ir, j, jstack, k, l
        integer :: cols, rows, sortcol

        double precision :: testval

        double precision :: data(cols,rows), repl(cols,1), temp(cols,1)

        integer, allocatable, dimension(:) :: istack
!       *********************************************************************

!       Set the constants.
        jstack = 0
        l = 1
        ir = rows
        m = 7
        nstack = 50
        
        allocate(istack(nstack))

1       if (ir-l<m) then
           do j=l+1, ir
              testval = data(sortcol,j)
              temp(:,1) = data(:,j)
              do i=j-1, l, -1
                 if (data(sortcol,i)<=testval) goto 2
                 data(:,i+1) = data(:,i)
              end do
              i = l - 1
2             data(:,i+1) = temp(:,1)
           end do
           if (jstack==0) return
           ir = istack(jstack)
           l = istack(jstack - 1)
           jstack = jstack - 2
        else
           k=(l+ir)/2
           temp(:,1) = data(:,k)
           data(:,k) = data(:,l+1)
           data(:,l+1) = temp(:,1)
           if (data(sortcol,l)>data(sortcol,ir)) then
              temp(:,1)=data(:,l)
              data(:,l) = data(:,ir)
              data(:,ir) = temp(:,1)
           end if
           if (data(sortcol,l+1)>data(sortcol,ir)) then
              temp(:,1) = data(:,l+1)
              data(:,l+1) = data(:,ir)
              data(:,ir) = temp(:,1)
           end if
           if (data(sortcol,l)>data(sortcol,l+1)) then
              temp(:,1) = data(:,l)
              data(:,l) = data(:,l+1)
              data(:,l+1) = temp(:,1)
           end if
           i=l+1
           j=ir
           testval = data(sortcol,l+1)
           repl(:,1) = data(:,l+1)
3          continue
           i=i+1
           if (data(sortcol,i)<testval) goto 3
4          continue
           j=j-1
           if (data(sortcol,j)>testval) goto 4
           if (j<i) goto 5
           temp(:,1) = data(:,i)
           data(:,i) = data(:,j)
           data(:,j) = temp(:,1)
           goto 3
5          data(:,l+1) = data(:,j)
           data(:,j) = repl(:,1)
           jstack = jstack + 2
           if (jstack>nstack) pause 'nstack too small!'
           if (ir-i+1>=j-l) then
              istack(jstack) = ir
              istack(jstack-1) = i
              ir = j-1
           else
              istack(jstack) = j-1
              istack(jstack-1) = l
              l = i
           end if
        end if
        goto 1

end subroutine qsort

!----------------------------------------------------------------------------

subroutine xsize(unit, cols, rows, blines)

!       This subroutine finds the number of columns and rows in an ASCII file.
!       The program can handle up to 3 blank lines at the top of the file and
!       the file line length must be 80 characters.

!       Variable definitions.
!       *********************************************************************

        character(len=80) :: line

        integer :: i, j, nline, unit, eof1, cols, rows, flag, blines 

        double precision :: dummy

!       *********************************************************************

        cols = 0
        rows = 0
        nline = 4

        if (blines/=0) then
           do j=1, blines
              read(unit,*)
           end do
           nline = 1
        end if

        blanks: do j=1, nline

                   read(unit,'(a)')line

                   columns: do i=1, 80
 
                            if (line(i:i)==' ') then
                               flag = 0
                            else
                               flag = 1
                            end if

                            if (i==1.and.flag==1) then
                               cols = cols + 1
                               cycle columns 
                            end if
                            if (flag==1.and.line(i-1:i-1)==' ') then
                               cols = cols + 1
                            end if

                   end do columns

                   if (cols==0) then
                      blines = blines + 1
                      cycle blanks
                   else
                      exit blanks
                   end if

        end do blanks

        if (cols==0) then
           print*
           print*, '*** WARNING: Failed to find any columns, please check file!'
           print*
           stop
        end if

!       Now measure the number of rows in the file.
        rewind(unit)
        eof1 = 0

        if (blines/=0) then
           do i=1, blines
              read(unit,*)
           end do
        end if

        do
          read (unit,*,iostat=eof1)dummy
          if (eof1/=0) exit
          rows = rows + 1
        end do

end subroutine xsize

!----------------------------------------------------------------------------

subroutine readit(unit, cols, rows, blines, data)

!       This subroutine reads in files into a data array.
!       The program can handle arrays with up to 5 columns.


!       Variable definitions.
!       *********************************************************************

        integer :: i, unit, eof1, cols, rows, blines

        double precision :: data(cols,rows)   

!       *********************************************************************

        eof1 = 0

!       Go past any blank lines that might exist at the top of the file.
        if (blines/=0) then
           do i=1, blines
              read(unit,*)
           end do
        end if

        do i=1, rows
           if (eof1/=0) exit
           if (cols==1) then
              read(unit,*,iostat=eof1)data(1,i)
           else if (cols==2) then
              read(unit,*,iostat=eof1)data(1,i), data(2,i)
           else if (cols==3) then
              read(unit,*,iostat=eof1)data(1,i), data(2,i), data(3,i)
           else if (cols==4) then
              read(unit,*,iostat=eof1)data(1,i), data(2,i), data(3,i), data(4,i)
           else if (cols==5) then
              read(unit,*,iostat=eof1)data(1,i), data(2,i), data(3,i), data(4,i), data(5,i)
           else
              print*
              print*, '*** WARNING: File contains more than 5 columns,'
              print*, '             I can not handle that!! Exiting...'
              print*
              stop
           end if
        end do

end subroutine readit

!----------------------------------------------------------------------------
subroutine writeit(unit, cols, rows, blines, data)

!       This subroutine writes an array into a file in real format.
!       The program can handle arrays with up to 5 columns.

!       Variable definitions.
!       *********************************************************************

        integer :: i, unit, cols, rows, blines

        double precision :: data(cols,rows)   

!       *********************************************************************

!       Go past any blank lines that might exist at the top of the file.
        if (blines/=0) then
           do i=1, blines
              write(unit,*)
           end do
        end if

        do i=1, rows
           if (cols==1) then
              write(unit,*)real(data(1,i))
           else if (cols==2) then
              write(unit,*)real(data(1,i)), real(data(2,i))
           else if (cols==3) then
              write(unit,*)real(data(1,i)), real(data(2,i)), real(data(3,i))
           else if (cols==4) then
              write(unit,*)real(data(1,i)), real(data(2,i)), real(data(3,i)), real(data(4,i))
           else if (cols==5) then
              write(unit,*)real(data(1,i)), real(data(2,i)), real(data(3,i)), real(data(4,i)),&
                           real(data(5,i))
           else
              print*
              print*, '*** WARNING: Array contains more than 5 columns,'
              print*, '             I can not handle that!! Exiting...'
              print*
              stop
           end if
        end do

end subroutine writeit

!----------------------------------------------------------------------------

END MODULE SUBS
