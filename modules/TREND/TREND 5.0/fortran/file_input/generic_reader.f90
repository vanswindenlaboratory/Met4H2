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

    !*************************************************************************************
    ! Written by Sven Pohl, Lehrstuhl fuer Thermodynamik RUB
    !*************************************************************************************
    !****************************************************************************************************************************
    ! Module for reading a file by path
    ! Single input: path to the file
    ! Returns: structure with the file content
    !##################################################################################################
    module file_by_path
    !************************************************************
    implicit none

    !************************************************************
    type ptr
        character(256),dimension(:), allocatable ::p
    end type
    !************************************************************
    type f_t
        type(ptr) :: act_file_ptr
    end type
    !************************************************************

    integer, parameter:: size_content = 4000

    type list_char
        character(256),allocatable,dimension(:) :: c
        integer :: nr
        !contains
        !procedure extend_int
    end type


    type list_int
        integer,allocatable,dimension(:) :: c
        integer :: nr
        !contains
        !procedure extend_int
    end type

    !interface list_char
    !module procedure extend_char
    !end interface
    !interface list_int
    !module procedure extend_int
    !end interface
    !************************************************************
    !interface for write_content
    interface write_content
    module procedure write_content_char, write_content_double, write_content_integer
    end interface
    !************************************************************
    !interface for standard constructor
    interface f_t
    module procedure ini
    end interface
    !************************************************************

    contains

    !subroutine append(this,element)
    !class(list) :: this
    !class(*), intent(in) :: element
    !class(*), allocatable :: tmp(:)
    !integer :: l
    !
    !if (.not.allocated(this%c)) then
    !    !select type(element)
    !    !type is (charcater)
    !    allocate(this%c(1),source=element)
    !    !end select
    !    this%c(1) = element
    !    this%num_nodes = 1
    !    return
    !else
    !    l = size(this%c,1)
    !    allocate(tmp(l+1),source=element)
    !    tmp(1:l) = this%c(1:l)
    !    tmp(l+1) = element
    !    deallocate(this%c)
    !    call move_alloc(tmp,this%c)
    !    this%num_nodes = this%num_nodes + 1
    !endif
    !end subroutine



    subroutine write_content_char(input,output)
    character(*), intent(in) :: input(:)
    character(*), intent(out) :: output(:)
    output = input
    end subroutine

    subroutine write_content_integer(input,output)
    integer,dimension(:), intent(in) :: input(:)
    integer,dimension(:), intent(out) :: output(:)
    output = input
    end subroutine

    subroutine write_content_double(input,output)
    double precision,dimension(:), intent(in) :: input
    double precision, dimension(:),intent(out) :: output
    output = input
    end subroutine

    !##################################################################################################
    type(f_t) function ini(path)
    character(*), intent(in) :: path
    type(list_char) :: lines
    character(256), allocatable :: output(:)
    character(256):: line
    integer:: err,cnt,i,file_unit

    !Init:
    cnt = 1

    allocate(lines%c(size_content))
    lines%c = ''
    lines%nr = 1

    !Open the correspondig file: defined by path
    open(newunit = file_unit,file=path,status='old',action='read',iostat=err)

    !read the lines from the file into a linked list ------------
    DO
        IF (err > 0)  THEN
            !ini%act_file_ptr%p => NULL()
            if (allocated(ini%act_file_ptr%p)) deallocate(ini%act_file_ptr%p)
            return
        ELSE IF (err < 0) THEN
            exit
        ELSE
            read(file_unit,'(A)',iostat=err) line
            if(err < 0) exit

            lines%c(lines%nr) = line

            if(trim(line) .eq. '@END') then
                exit ! only take content until @END
            end if

            lines%nr = lines%nr + 1

        END IF

    end do

    !Get the content and write it to the output structure
    allocate(output(lines%nr)) !ggf -1
    call write_content(lines%c,output)
    !select type (lines%c)
    !type is(character(len=*))
    !    write(*,*)len(lines%c)
    !end select
    !do i=1,lines%num_nodes
    !    call get_item(lines,i,output(i))
    !end do

    !pointer to output
    !Return value of the function
    allocate(ini%act_file_ptr%p,mold=output)
    ini%act_file_ptr%p = output

    close(file_unit)


    end function ini
    !##################################################################################################
    end module
    !****************************************************************************************************************************

    !###########################################################
    !Generic reader for TREND 4.0, Sven Pohl & Monika Thol
    !###########################################################
    !This generic file strcuture is for the reading of thermodynamic model files
    ! Methods are:
    ! get_content(path,name,folder)
    ! get_sub_string_lines()
    ! get_empty_string_lines()
    ! get_content_until_empty()
    ! get_ref_from_sub_content()
    ! delete_lines()
    !*********************************************************************************************************
    module file
    use hard_coded_files
    use file_by_path
    implicit none

    type file_t
        integer:: nr_lines,act_line
        character(256),dimension(:), allocatable::content
        logical :: is_hc
        integer:: errorflag
    contains
    !procedure:: get_content! =>get_content
    !procedure:: get_sub_string_lines !=> get_sub_string_lines
    !procedure:: get_lines !=> get_lines
    !procedure:: get_empty_string_lines !=> get_empty_string_lines
    !procedure:: get_content_until_empty_line!=>get_content_until_empty_line
    !procedure:: get_ref_from_sub_content!=>get_ref_from_sub_content
    !procedure:: delete_lines!=> delete_lines
    end type

    interface file_t
    module procedure init_file_t
    end interface

    type return_pos_t
        character(256),allocatable,dimension(:)::p
        integer,allocatable,dimension(:) :: pos_id
    end type

    !************************************************************
    contains

    !constructor
    type(file_t) function init_file_t(is_hc_in)
    logical :: is_hc_in
    init_file_t%is_hc = is_hc_in
    init_file_t%nr_lines = 0
    init_file_t%act_line = 0
    init_file_t%errorflag = 0
    end function
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the head of a fluid file
    ! Input: path<character>: absolut path to the model
    !        name<character array>: name of the fluid
    !        folder<character>: name of the folder
    !--------------------------------------------------------------------------------------------------
    subroutine get_content(this,path,name,folder)
    type(file_t) :: this
    character(*), intent(in) :: path,folder
    character(*),dimension(:), intent(in):: name
    character(256),target, dimension(:),allocatable:: err_msg
    type(hc_t):: hc
    type(f_t):: ff
    integer:: last_slash,i,tab_pos, index_end
    !Decide wheter a hardcoded file or a "real" file is readed
    last_slash = index(path,'\',BACK=.true.)
    !HC is present
    if(this%is_hc) then
        hc = hc_t(name,folder) !create struct and get the file content
        if (associated(hc%act_file_ptr%p)) then
            this%content = hc%act_file_ptr%p
            index_end = findloc(this%content(:)(1:4),'@END',dim=1) 
            if (index_end .ne. 0) this%content = this%content(1:index_end)
            this%nr_lines = size(this%content)
            this%act_line = 1
        end if
        !Standard input file
    else
        ff = f_t(path)
        if (allocated(ff%act_file_ptr%p)) then
            this%content = ff%act_file_ptr%p
            this%nr_lines = size(this%content)
            this%act_line = 1
        endif
    end if

    if(.not.allocated(this%content)) then
        allocate(err_msg(1))
        err_msg  = 'ERROR: The content could not be read'
        this%content = err_msg
        this%errorflag = -7878
    else
        this%errorflag = 0
        do i = 1,size(this%content,1)
            tab_pos = index(this%content(i),achar(9))
            do while(tab_pos /= 0)
                if(tab_pos /= 0) then
                    this%content(i)(tab_pos:tab_pos) = ''
                end if
                tab_pos = index(this%content(i),achar(9))
            end do
        enddo
        this%content = ADJUSTL(this%content)
    end if
    end subroutine
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Get a  sub-content of a file
    ! Input: start<integer>: line nr where to start
    !        ende<integer>: line nr where to end
    !--------------------------------------------------------------------------------------------------
    function get_lines(this,start,ende)
    type(file_t):: this
    integer, intent(in):: start,ende
    character(256),target, dimension(:),allocatable:: err_msg
    character(256),dimension(:),allocatable::get_lines
    if(allocated(this%content)) then
        get_lines = this%content(start:ende)
    else
        err_msg = 'No content available'
        get_lines = err_msg
    end if
    end function
    ! ********** GET LINES OF CONTENT **************************
    !****************************************************************************************************************************************************************
    ! ********** GET ALL LINES CONTAINING SUBSTRING ************
    function get_sub_string_lines(this,substrings,errorflag)
    type(file_t):: this
    character(*),dimension(:):: substrings
    character(256),target, dimension(:),allocatable:: err_msg
    type(return_pos_t) :: get_sub_string_lines
    integer:: loop_file,loop_sub_strings,lines_cntr,i,errorflag
    character(256),dimension(:),allocatable:: output
    integer,dimension(:),allocatable:: positions
    character(256) :: line
    type(list_char):: lines
    type(list_int):: pos

    allocate(lines%c(size_content))
    lines%c = ''
    allocate(pos%c(size_content))
    pos%c =  0
    lines%nr = 1
    pos%nr = 1

    !Init -----------------------------------------------------------
    errorflag = 0
    lines_cntr = 0
    if(allocated(this%content)) then
        do loop_sub_strings=1,size(substrings)
            do loop_file=1,this%nr_lines
                line = this%content(loop_file)
                if(INDEX( line, trim(substrings(loop_sub_strings)) ) /= 0) then
                    !call append(lines,line)
                    !call append(pos,loop_file)
                    lines%c(lines%nr) = line
                    lines%nr = lines%nr + 1
                    pos%c(pos%nr) = loop_file
                    pos%nr = pos%nr + 1

                    lines_cntr = lines_cntr + 1
                end if
            end do
        end do
        if(lines_cntr == 0) then
            errorflag = -45257 !none of the substrings was found
            return
        end if
        !get the content --------------------------------------------
        allocate(output(lines_cntr),positions(lines_cntr),get_sub_string_lines%pos_id(lines_cntr))
        call write_content (pos%c,positions )
        call write_content (lines%c,output )
        get_sub_string_lines%pos_id = positions
        get_sub_string_lines%p = output(:)
    else
        err_msg = 'No content available'
        get_sub_string_lines%p = err_msg
    end if
    end function
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Get all empty character lines of file ( with line index)
    !--------------------------------------------------------------------------------------------------
    function get_empty_string_lines(this)
    type(file_t):: this
    character(256),target, dimension(:),allocatable:: err_msg
    type(return_pos_t) :: get_empty_string_lines
    integer:: loop_file,lines_cntr,i
    character(256),dimension(:),allocatable:: output
    integer,dimension(:),allocatable:: positions
    character(256):: line
    type(list_char):: lines
    type(list_int):: pos

    !Init -----------------------------------------------------------
    allocate(lines%c(size_content))
    lines%c = ''
    allocate(pos%c(size_content))
    pos%c = 0
    lines%nr = 1
    pos%nr = 1

    lines_cntr = 0
    if(allocated(this%content)) then

        do loop_file=1,this%nr_lines
            line = this%content(loop_file)
            if( len_trim(line) .eq. 0) then
                !call append(lines,line)
                !call append(pos,loop_file)
                lines%c(lines%nr) = line
                lines%nr = lines%nr + 1
                pos%c(pos%nr) = loop_file
                pos%nr = pos%nr + 1
                lines_cntr = lines_cntr + 1
            end if
        end do
        if (lines_cntr > 0) then
            !get the content --------------------------------------------
            allocate(output(lines_cntr),positions(lines_cntr),get_empty_string_lines%pos_id(lines_cntr))
            output = ''
            call write_content (pos%c,positions )
            get_empty_string_lines%pos_id = positions
            get_empty_string_lines%p = output(:)
        end if
    else
        err_msg = 'No content available'
        get_empty_string_lines%p = err_msg
    end if
    end function get_empty_string_lines
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Get a  sub-content of a file froms line nr start to the next empty line
    ! Input: start<integer>: line nr where to start
    !--------------------------------------------------------------------------------------------------
    function get_content_until_empty_line(this,start)
    type(file_t):: this
    integer:: loop_file,lines_cntr,i,start
    character(256),dimension(:),allocatable::get_content_until_empty_line
    character(256):: line
    type(list_char):: lines

    !allocate(lines%c(this%content))
    allocate(lines%c(size_content))
    lines%c = ''
    lines%nr = 1


    if (start == 0) return !return when given model is not present for current substance
    !Init -----------------------------------------------------------
    lines_cntr = 0

    do loop_file=start,size(this%content)
        line = this%content(loop_file)
        if( len_trim(line) .eq. 0 ) then
            get_content_until_empty_line = this%content(start:loop_file)
            return
        elseif (line(1:len_trim(line)) == '@END' .or. loop_file == size(this%content)) then
            get_content_until_empty_line = this%content(start:loop_file)
            return
        end if
    end do

    end function
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Get a  reference block from a file
    ! Input: start<integer>: line nr where to start
    !        size_in<integer>: size of the character array
    !--------------------------------------------------------------------------------------------------
    function get_ref_from_sub_content(this,start,size_in)
    type(file_t) :: this
    integer:: size_in
    character(size_in)::get_ref_from_sub_content
    character(256),dimension(:),allocatable:: help_array
    integer:: loop_file,lines_cntr,i,start,lit_lines_found,start_line
    character(256):: line

    !x1    write(*,*) start
    if (start .eq. 0) then
        !return !return when given model is not present for current substance
    else
        !Init -----------------------------------------------------------
        lines_cntr = start
        lit_lines_found = 0
        start_line= 0
        get_ref_from_sub_content = ''
        do loop_file=start,size(this%content)
            line = this%content(loop_file)
            lines_cntr = lines_cntr + 1
            if( line(1:5) .eq. "?````" .or.  trim(line) .eq. '?LITERATURE REFERENCE \'.or.  trim(line) .eq. '?LITERATURE REFERENCE\') then
                start_line = lines_cntr
                do while(trim(line) .ne. '?' .and. trim(line) .ne. '?\' .and. trim(line) .ne. '!end of info section')
                    lit_lines_found = lit_lines_found + 1
                    line = this%content(lines_cntr)
                    lines_cntr = lines_cntr + 1
                end do
                !get the content --------------------------------------------
                help_array = this%content(start_line:start_line+lit_lines_found-2)

                do i=1,size(help_array)
                    if (index(help_array(i),'?') /= 0 ) then
                        help_array(i) = help_array(i)(index(line,'?')+1:)
                    elseif (index(help_array(i),'!') /= 0) then
                        help_array(i) = help_array(i)(index(line,'!')+1:)
                    end if
                    get_ref_from_sub_content = trim(get_ref_from_sub_content)  // trim(help_array(i))
                end do
                return
            end if
        end do
    end if
    end function
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Delete all lines that start with: "?" or "!"  or ":"
    !--------------------------------------------------------------------------------------------------
    subroutine delete_lines(this,sub_content)
    type(file_t):: this
    character(256),dimension(:),allocatable :: sub_content
    character(256):: line
    integer:: loop_file,i
    type(list_char):: lines
    character(256),dimension(:),allocatable:: output

    allocate(lines%c(size_content))
    lines%c = ''
    lines%nr = 0

    if (.not. allocated(sub_content)) return !return when given model is not present for current substance
    !Init -----------------------------------------------------------
    do loop_file=1,size(sub_content)
        line = sub_content(loop_file)
        if(line(1:1) .ne. '?' .and. line(1:1) .ne. '!' .and. line(1:1) .ne. ':') then
            !call append(lines,line)
            lines%nr = lines%nr + 1
            lines%c(lines%nr) = line
        end if
    end do


    allocate(output(lines%nr))
    call write_content(lines%c,output)
    sub_content = ''
    sub_content(1:size(output)) = output(:)

    end subroutine
    !##########################################################################################################
    end module file