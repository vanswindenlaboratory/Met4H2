    !*************************************************************************************
    !				TREND Version 5.0
    !		   Thermodynamic Reference & Engineering Data
    !
    !- software for the calculation of thermodynamic and other properties -
    !
    !Copyright (C) 2020,  Prof. Dr.-Ing. R.Span
    !                     Lehrstuhl fuer Thermodynamik
    !                     Ruhr-Universitaet Bochum
    !                     Universitaetsstr. 150
    !                     D-44892 Bochum
    !
    !Cite as: Span, R.; Beckmüller, R.; Hielscher, S.; Jäger, A.; Mickoleit, E.; 
	!          Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020): 	
    !          TREND. Thermodynamic Reference and Engineering Data 5.0. 
    !          Lehrstuhl für Thermodynamik, Ruhr-Universität Bochum.

    !
    !This program is free software: you can redistribute it and/or modify
    !it under the terms of the GNU General Public License as published by
    !the Free Software Foundation, either version 3 of the License, or
    !(at your option) any later version.
    !
    !This program is distributed in the hope that it will be useful,
    !but WITHOUT ANY WARRANTY; without even the implied warranty of
    !MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !GNU General Public License for more details.
    !
    !You should have received a copy of the GNU General Public License
    !along with this program.  If not, see < http://www.gnu.org/licenses/ >.
    !
    !*************************************************************************************

    ! module for file utility.f90
    module utility_module
    
    contains

    subroutine uppertolower_char(string, n)

    implicit none

    !Variable for changing the fluid file names to lower case letters
    integer, intent(in) :: n
    character*n :: string
    character (26):: lower,upper

    integer:: i, j, k                                           ! count variables for the moles vector
    integer :: f                                              ! amount of fluid in vector

    !---------------------------------------------------------------------------
    !-----------------------------   LOWER CASE STRING   -----------------------

    !For the SRK the fluid names have to be lower case
    !TRANSFORM ALL FLUID NAMES TO LOWER CASE LETTERS
    !Andreas March 2013
    lower = 'abcdefghijklmnopqrstuvwxyz'
    upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    !write upper cases in lower cases for the components vector
    do i = 1,len(trim(string))
        do j = 1,26
            if (string(i:i) == upper(j:j)) then
                string(i:i) = lower(j:j)
                exit
            end if
        end do
    end do


    end subroutine uppertolower_char



    function is_infinity (value) result(check)
    logical :: check
    double precision :: value

    if (dabs(value) == -dlog(0.d0)) then
        check = .true.
    else
        check = .false.
    endif
    end function

    !Andreas, May 2014
    !Interface routine for TREND implementation into MATLAB 64 bit
    !Based on earlier work of J. Gernert

    !*******************************************************
    ! subroutine converts a string passed by a C application
    ! to a string readable by a fortran function
    !*******************************************************
    subroutine cstringtof(stringf, nout, stringc)
    ! Convert fluid from C to fortran srting
    ! /!\ Watch out if nout is less than 2...
    USE, intrinsic:: ISO_C_BINDING

    implicit none

    integer::nout,i,leng
    character(kind=c_char), dimension(*):: stringc


    character(nout)::stringf

    leng=0
    do
        if (stringc(leng+1) == C_NULL_CHAR) exit
        leng = leng + 1
    end do
    write(stringf,*)(stringc(i),i=1,leng)
    stringf=stringf(2:nout)


    end subroutine cstringtof

    function countsubstring(s1, s2) result(c)
    !DEC$ ATTRIBUTES DLLEXPORT :: countsubstring
    character(*), intent(in) :: s1, s2
    integer :: c, p, posn

    c = 0
    if(len(s2) == 0) return
    p = 1
    do
        posn = index(s1(p:), s2)
        if(posn == 0) return
        c = c + 1
        p = p + posn + len(s2)
    end do
    end function

    recursive subroutine split_char(char_arg, vector, index_arg, separate)
    !DEC$ ATTRIBUTES DLLEXPORT :: split_char
    implicit none
    character (*) :: char_arg
    character (*), dimension(*) :: vector
    character (1) :: separate
    integer :: n, index_arg
    n = countsubstring(char_arg,separate)
    if (n /= 0) then
        read(char_arg(1:index(char_arg,separate)-1),*)vector(index_arg)
        index_arg = index_arg + 1
        char_arg = char_arg(index(char_arg,separate)+1:)
        call split_char(char_arg,vector,index_arg, separate)
    endif

    read(char_arg,*)vector(index_arg)

    end subroutine

    recursive subroutine split_dbl(char_arg, vector, index_arg, separate)
    !DEC$ ATTRIBUTES DLLEXPORT :: split_dbl
    implicit none
    character (*) :: char_arg
    double precision, dimension(*) :: vector
    character (1) :: separate
    integer :: n, index_arg
    n = countsubstring(char_arg,separate)
    if (n /= 0) then
        read(char_arg(1:index(char_arg,separate)-1),*)vector(index_arg)
        index_arg = index_arg + 1
        char_arg = char_arg(index(char_arg,separate)+1:)
        call split_dbl(char_arg,vector,index_arg, separate)
    endif

    read(char_arg,*)vector(index_arg)

    end subroutine


    end module utility_module

    ! -----------------------------------------------
MODULE String_Functions  ! by David Frank  dave_frank@hotmail.com
IMPLICIT NONE            ! 

! Copy (generic) char array to string or string to char array
! Clen           returns same as LEN      unless last non-blank char = null
! Clen_trim      returns same as LEN_TRIM    "              "
! Ctrim          returns same as TRIM        "              "
! Count_Items    in string that are blank or comma separated
! Reduce_Blanks  in string to 1 blank between items, last char not blank
! Replace_Text   in all occurances in string with replacement string
! Spack          pack string's chars == extract string's chars
! Tally          occurances in string of text arg
! Translate      text arg via indexed code table
! Upper/Lower    case the text arg

INTERFACE Copy    ! generic
   MODULE PROCEDURE copy_a2s, copy_s2a
END INTERFACE Copy

CONTAINS
! ------------------------
PURE FUNCTION Copy_a2s(a)  RESULT (s)    ! copy char array to string
CHARACTER(*),INTENT(IN) :: a(:)
CHARACTER(SIZE(a,1)) :: s
INTEGER :: i
i=1
s = ''

do while(len_trim(a(i)) /= 0)
    s = trim(s) // ' ' // trim(a(i))
    i = i + 1
end do

END FUNCTION Copy_a2s

! ------------------------
PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
CHARACTER(*),INTENT(IN) :: s
CHARACTER :: a(LEN(s))
INTEGER :: i
DO i = 1,LEN(s)
   a(i) = s(i:i)
END DO
END FUNCTION Copy_s2a

! ------------------------
PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
INTEGER :: i
Clen = LEN(s)
i = LEN_TRIM(s)
IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
END FUNCTION Clen

! ------------------------
PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
INTEGER :: i                       ! then len of C string is returned, note:
                                   ! Ctrim is only user of this function
i = LEN_TRIM(s) ; Clen_trim = i
IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
END FUNCTION Clen_trim

! ----------------
FUNCTION Ctrim(s1)  RESULT(s2)     ! returns same result as TRIM unless:
CHARACTER(*),INTENT(IN)  :: s1     ! last non-blank char is null in which
CHARACTER(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
s2 = s1                            ! are output
END FUNCTION Ctrim

! --------------------
INTEGER FUNCTION Count_Items(s1)  ! in string or C string that are blank or comma separated
CHARACTER(*) :: s1
CHARACTER(Clen(s1)) :: s
INTEGER :: i, k

s = s1                            ! remove possible last char null
k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
DO i = 1,LEN_TRIM(s)-1
   IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
                    .AND.s(i+1:i+1) == ' '.OR.s(i+1:i+1) == ',') k = k+1
END DO
Count_Items = k
END FUNCTION Count_Items

! --------------------
FUNCTION Reduce_Blanks(s)  RESULT (outs)
CHARACTER(*)      :: s
CHARACTER(LEN_TRIM(s)) :: outs
INTEGER           :: i, k, n

n = 0  ; k = LEN_TRIM(s)          ! k=index last non-blank (may be null)
DO i = 1,k-1                      ! dont process last char yet
   n = n+1 ; outs(n:n) = s(i:i)
   IF (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
END DO
n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
IF (n < k) outs(n+1:) = ' '       ! pad trailing blanks
END FUNCTION Reduce_Blanks

! ------------------
FUNCTION Replace_Text (s,text,rep)
CHARACTER(*)        :: s,text,rep
CHARACTER(len(s)) :: Replace_Text
integer:: pos
Replace_Text = s
pos = index(Replace_Text,text)
if(pos == 0) return
do while(index(Replace_Text,text) /= 0) 
    pos = index(Replace_Text,text)
    Replace_Text(pos:pos) = trim(rep)
end do
 

END FUNCTION Replace_Text

! ---------------------------------
FUNCTION Spack (s,ex)  RESULT (outs)
CHARACTER(*) :: s,ex
CHARACTER(LEN(s)) :: outs
CHARACTER :: aex(LEN(ex))   ! array of ex chars to extract
INTEGER   :: i, n

n = 0  ;  aex = Copy(ex)
DO i = 1,LEN(s)
   IF (.NOT.ANY(s(i:i) == aex)) CYCLE   ! dont pack char
   n = n+1 ; outs(n:n) = s(i:i)
END DO
outs(n+1:) = ' '     ! pad with trailing blanks
END FUNCTION Spack

! --------------------
INTEGER FUNCTION Tally (s,text)
CHARACTER(*) :: s, text
INTEGER :: i, nt

Tally = 0 ; nt = LEN_TRIM(text)
DO i = 1,LEN(s)-nt+1
   IF (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
END DO
END FUNCTION Tally

! ---------------------------------
FUNCTION Translate(s1,codes)  RESULT (s2)
CHARACTER(*)       :: s1, codes(2)
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER            :: i, j

DO i = 1,LEN(s1)
   ch = s1(i:i)
   j = INDEX(codes(1),ch) ; IF (j > 0) ch = codes(2)(j:j)
   s2(i:i) = ch
END DO
END FUNCTION Translate

! ---------------------------------
FUNCTION Upper(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
   s2(i:i) = ch
END DO
END FUNCTION Upper

! ---------------------------------
FUNCTION Lower(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
   s2(i:i) = ch
END DO
END FUNCTION Lower

END MODULE String_Functions