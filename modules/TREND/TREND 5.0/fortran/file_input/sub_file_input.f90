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

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Written by Sven Pohl & Monika Thol, Lehrstuhl fuer Thermodynamik RUB
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Module for collecting infos given by the path and model (eos_type and mix_indicator)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    module file_info

    ! Load modules
    use module_all_types
    use file

    implicit none

    ! Type definitions
    type file_info_t
        character(255) :: user_path
        character(255) :: hc_pure,hc_mix
        character(255) :: pathforsingles
        integer :: path_len
        logical :: is_hc
        ! Possible HC Folders
        character(50) , dimension(5) :: hc_path =   ['hc','hc_eos-cg-2016','hc_eos-cg-2019','hc_gerg-2008','hc_eos-lng']
        character(50) , dimension(5) :: hc_pure_all =   ['fluids_pure','eos-cg-2016_pure','eos-cg-2019_pure','gerg-2008_pure','eos-lng_pure']
        character(50) , dimension(5) :: hc_mix_all  =   ['binary_mix_files','eos-cg-2016_mix','eos-cg-2019_mix','gerg-2008_mix','eos-lng_mix']

        character(20) , dimension(7) :: hc_pure_models  = ['srk_pure','pr_pure','lkp_pure','gen_eq','saft_pure','rkm','costald']

        ! Lists for the paths for pure models
        character(255) , allocatable :: hc_list_pure(:,:)
        character(255) , allocatable :: list_pure(:,:)

        ! Lists for the paths for mixture models
        character(:) , allocatable:: hc_mix_path,mix_path

    contains

    ! Methods of the type file_info_t
    procedure :: set_paths => set_paths
    procedure :: set_hc_pure_folders => set_hc_pure_folders
    procedure :: set_eos_model_from_path_pure => set_eos_model_from_path_pure
    procedure :: set_hc_mix_folders => set_hc_mix_folders
    procedure :: set_eos_model_from_path_mix => set_eos_model_from_path_mix
    procedure :: get_helmholtz_path_pure_hc => get_helmholtz_path_pure_hc
    procedure :: get_helmholtz_path_pure => get_helmholtz_path_pure
    procedure :: get_helmholtz_path_mix_hc => get_helmholtz_path_mix_hc

    end type

    ! Interface for constructor
    interface file_info_t
    module procedure ini_f
    end interface

    contains

    ! #####################################################################################################################################
    ! Constructor for file_info_t
    ! Checks if a path exists / or HC models exists
    !-------------------------------------------------------------------------------------
    type(file_info_t) function ini_f(gl,path_in,errorflag)
    type(type_gl) :: gl
    character(*) :: path_in
    integer :: errorflag
    character(1), dimension(2) :: enders = ['/','\']
    integer :: hc_index
    logical :: DIR_EXISTS

    ! save user path
    ini_f%user_path =  path_in

    ! check if the path is existent
    inquire(DIRECTORY=ini_f%user_path,EXIST=DIR_EXISTS)

    if(.not. DIR_EXISTS) then
        hc_index = index(ini_f%user_path,'hc')
        if( hc_index == 0) then
            ini_f%is_hc =  .false.
        else
            ini_f%is_hc =  .true.
        end if
    end if

    ! If the dir does not exists, check fo "hc", if no hc return with error
    if(DIR_EXISTS .or. ini_f%is_hc ) then

        ! set error to 0
        errorflag = 0

        ! Allocate arrays of paths
        if(.not.allocated(ini_f%list_pure)) then
            allocate(ini_f%list_pure(gl%ncomp,2))
        end if

        if(.not.allocated(ini_f%hc_list_pure)) then
            allocate(ini_f%hc_list_pure(gl%ncomp,2))
        end if

        ! Get the length of the path
        ini_f%path_len = len_trim(ini_f%user_path)

        ! Check if at the end of the path some Slash or Backslash is present
        if( any( enders == ini_f%user_path(ini_f%path_len:ini_f%path_len) ) ) then
            ini_f%user_path(ini_f%path_len:ini_f%path_len) = ' '
        end if

        ! Set Slash at the end of path (consistent paths!)
        ini_f%user_path(1:ini_f%path_len) = trim(ini_f%user_path) // '\'

        
        ! init hc paths
        ini_f%hc_pure = ''
        ini_f%hc_mix = ''
        
        if(DIR_EXISTS) then ! normal path
            ini_f%is_hc =  .false.
        else ! hc path
            ini_f%is_hc =  .true.
        end if

        ! Return with error
    else
        errorflag = -7878
    end if

    end function
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Set the path for the models
    ! The paths are setted for the pure and mixture models
    !-------------------------------------------------------------------------------------
    subroutine set_paths(this,gl,errorflag)

    class(file_info_t) :: this
    type(type_gl) :: gl
    integer :: nrsubst,errorflag
    logical :: FILE_EXISTS

    errorflag = 0

    ! Check if any comp name is a cas number, if yes, replace with correct fluid name given in cas list
    do nrsubst =1,gl%ncomp

        ! HC CASE:
        if(this%is_hc) then

            this%mix_path = ''
            call this%set_hc_pure_folders(gl%eq_type(nrsubst),nrsubst,errorflag)

            ! CASE with normal file path
        else

            this%hc_mix_path = ' '
            call  this%set_eos_model_from_path_pure(gl,gl%eq_type(nrsubst),nrsubst,errorflag)

        end if

    end do   ! end of nrsubst loop

    ! Check if a mixture is present
    ! If Yes: create the mixing paths for hc or normal path
    if(gl%ncomp > 1 .and. gl%eq_type(nrsubst) .ne. 1) then
        if(this%is_hc) then
            this%mix_path = 'hc/'
            call  this%set_hc_mix_folders(gl%mix_type,nrsubst,errorflag)
        else
            call  this%set_eos_model_from_path_mix(gl%mix_type)
        end if
    end if

    end subroutine
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Set the paths for HC -> pure fluids
    !-------------------------------------------------------------------------------------
    subroutine set_hc_pure_folders(this,eqtype,nrsubst,errorflag)

    class(file_info_t) :: this
    integer :: eqtype,nrsubst,errorflag,hc_index

    hc_index = findloc(this%hc_path,this%user_path(1:this%path_len-1),1)


    if (eqtype .eq. 1) then
        this%hc_list_pure(nrsubst,1) = this%get_helmholtz_path_pure_hc(errorflag)
    else
        if(hc_index .eq. 1) then
            if((eqtype .eq. 2) .or.  (eqtype .eq. 7))    then
                this%hc_list_pure(nrsubst,1) = 'srk_pure'
            elseif(    eqtype .eq. 3)     then
                this%hc_list_pure(nrsubst,1) = 'pr_pure'
            elseif(    eqtype .eq. 4)   then
                this%hc_list_pure(nrsubst,1) = 'lkp_pure'
            else if(    eqtype .eq. 6)    then
                this%hc_list_pure(nrsubst,1) = 'saft_pure'
            elseif(    eqtype .eq. 9)  then
                this%hc_list_pure(nrsubst,1) = 'srk_pure'
                this%hc_list_pure(nrsubst,2) = 'costald'
            elseif( (eqtype .eq. 51) .or. (eqtype .eq. 52) .or. (eqtype .eq. 53) ) then
                this%hc_list_pure(nrsubst,1) = 'gen_eq'
            elseif(    (eqtype .eq. 8) .or. (eqtype .eq. 81))    then
                this%hc_list_pure(nrsubst,1) = 'srk_pure'
                this%hc_list_pure(nrsubst,2) = 'rkm'
            end if
        else
            errorflag =  -7878
        end if
    end if

    end subroutine
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Set the paths from directory -> pure fluids
    !-------------------------------------------------------------------------------------
    subroutine set_eos_model_from_path_pure(this,gl,eqtype,nrsubst,errorflag)

    class(file_info_t) :: this
    type(type_gl) :: gl
    integer :: eqtype,nrsubst,errorflag

    if( eqtype .eq. 1) then
        this%list_pure(nrsubst,1) =   this%get_helmholtz_path_pure(gl,nrsubst,errorflag)
    elseif((eqtype .eq. 2) .or.  (eqtype .eq. 7))    then
        this%list_pure(nrsubst,1) =  'srk/srk.fld'
    elseif(    eqtype .eq. 3)     then
        this%list_pure(nrsubst,1) = 'pr/pr.fld'
    elseif(    eqtype .eq. 4)   then
        this%list_pure(nrsubst,1) = 'lkp/lkp.fld'
    else if(    eqtype .eq. 6)    then
        this%list_pure(nrsubst,1) = 'saft/saft.fld'
    elseif( eqtype .eq. 9)  then
        this%list_pure(nrsubst,1) = 'srk/srk.fld'
        this%list_pure(nrsubst,2) = 'costald/costald.fld'
    elseif( (eqtype .eq. 51) .or. (eqtype .eq. 52) .or. (eqtype .eq. 53) ) then
        this%list_pure(nrsubst,1) = 'gen_eq/gen_eq.fld'
    elseif( (eqtype .eq. 8) .or. (eqtype .eq. 81))    then
        this%list_pure(nrsubst,1) = 'srk/srk.fld'
        this%list_pure(nrsubst,2) = 'rkm'
    end if

    ! add the directory to the path
    if( eqtype .ne. 1) then
        this%list_pure(nrsubst,1) = trim(this%user_path) //  trim(this%list_pure(nrsubst,1))
        this%list_pure(nrsubst,2) = trim(this%user_path) //  trim(this%list_pure(nrsubst,2))
    end if

    end subroutine
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Set the paths for HC -> mixtures
    !-------------------------------------------------------------------------------------
    subroutine set_hc_mix_folders(this,mixtype,nrsubst,errorflag)

    class(file_info_t) :: this
    integer :: mixtype,nrsubst,errorflag

    if( (mixtype .eq. 1) .or. (mixtype .eq. 110) .or. (mixtype .eq. 111) .or. (mixtype .eq. 120) .or. (mixtype .eq. 121)) then
        this%hc_mix_path = this%get_helmholtz_path_mix_hc(errorflag)
    elseif(    mixtype .eq. 2)     then
        this%hc_mix_path = 'srk_mix'
    elseif(    mixtype .eq. 3)   then
        this%hc_mix_path = 'pr_mix'
    elseif(    mixtype .eq. 4 ) then
        this%hc_mix_path = 'lkp_mix'
    else if(    mixtype .eq. 6)    then
        this%hc_mix_path = 'binary_saft'
    elseif( (mixtype .eq. 12) .or. (mixtype .eq.  22))    then
        this%hc_mix_path = 'helm_ge'
    elseif(    mixtype .eq. 13)    then
        this%hc_mix_path = ''
    end if

    end subroutine
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Set the paths from directory -> pure fluids
    !-------------------------------------------------------------------------------------
    subroutine set_eos_model_from_path_mix(this,mixtype)

    class(file_info_t) :: this
    integer :: mixtype
    !PSRK mixing rules for parameter "a". For PSRK the UNIFAC parameters are needed which are stored in the files "Interac.par" and "rk_qk.par"
    !Andreas Jäger, January 2017
    !Use same routine for reading Helmholtz+gE parameters (Mix_type == 12). Andreas Jäger, August 2017
    if( (mixtype .eq. 1) .or. (mixtype .eq. 110) .or. (mixtype .eq. 111) .or. (mixtype .eq. 120) .or. (mixtype .eq. 121)) then
        this%mix_path =  'binary_mix_files/'
    elseif(    mixtype .eq. 2)     then
        this%mix_path = 'srk/srk.mix'
    elseif(    mixtype .eq. 3)   then
        this%mix_path = 'pr/pr.mix'
    elseif(    mixtype .eq. 4 ) then
        this%mix_path = 'lkp/lkp.mix'
    else if(    mixtype .eq. 6)    then
        this%mix_path = 'saft/saft.mix'
    else
        this%mix_path = ''
    end if

    ! add the directory to the path
    this%mix_path = trim(this%user_path) //  trim(this%mix_path)

    end subroutine
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Get the paths for HC for all pure fluids for Helmholtz EOS models
    !-------------------------------------------------------------------------------------
    function get_helmholtz_path_pure_hc(this,errorflag)

    class(file_info_t) :: this
    character(:),allocatable :: get_helmholtz_path_pure_hc
    integer :: hc_index,errorflag

    ! get index of hc model
    hc_index = findloc(this%hc_path,this%user_path(1:this%path_len-1),1)
    if(hc_index > 0) then
        this%hc_pure = this%hc_pure_all(hc_index)
    else
        errorflag = -7878
    end if

    get_helmholtz_path_pure_hc = this%hc_pure

    end function

    function get_helmholtz_path_mix_hc(this,errorflag)

    class(file_info_t) :: this
    character(:),allocatable :: get_helmholtz_path_mix_hc
    integer :: hc_index,errorflag

    ! get index of hc model
    hc_index = findloc(this%hc_path,this%user_path(1:this%path_len-1),1)
    if(hc_index > 0) then
        this%hc_mix = this%hc_mix_all(hc_index)
    else
        errorflag = -7878
    end if

    get_helmholtz_path_mix_hc = this%hc_mix

    end function
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Get the paths from directory for all pure fluids for Helmholtz EOS models
    !-------------------------------------------------------------------------------------
    function get_helmholtz_path_pure(this,gl,nrsubst,errorflag)

    character(:),allocatable :: get_helmholtz_path_pure
    class(file_info_t) :: this
    type(type_gl) :: gl
    integer :: errorflag,nrsubst
    logical ::  FILE_EXISTS

    this%pathforsingles =   trim(this%user_path) // 'fluids/'

    ! Check if the fluid name is CAS-NR
    ! If yes: replace CAS-NR with correct fluid name
    if(index(gl%components(nrsubst),'-') > 1) then
        gl%components(nrsubst)  = cas_to_fld(gl%components(nrsubst) ,nrsubst,this%pathforsingles,this%is_hc,this%hc_pure,errorflag)
    end if

    get_helmholtz_path_pure = trim(this%pathforsingles) // trim(gl%components(nrsubst)) // '.fld'

    ! check if the file exists
    inquire(FILE=get_helmholtz_path_pure,EXIST=FILE_EXISTS)

    if(.not. FILE_EXISTS) then
        get_helmholtz_path_pure = trim(this%pathforsingles) // trim(gl%components(nrsubst)) // '.ppf'

        ! check again if the ppf file exists
        inquire(FILE=get_helmholtz_path_pure,EXIST=FILE_EXISTS)
        if(.not. FILE_EXISTS) then
            errorflag = -7878
        end if
    end if

    end function
    ! #####################################################################################################################################

    ! #####################################################################################################################################
    ! Check if the fluid name is a CAS number
    ! Input: gl<type_gl>: global variable structure
    !        nrsubst<integer>: Number of actual component
    !        pathforsingles<character>: absolut path to the model
    !        hc_exists<logical>: pointer for hc
    !        hc_folder<character>: name of the hardcoded folder
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    character(30)  function cas_to_fld(cas_nr,nrsubst,pathforsingles,hc_exists,hc_folder,errorflag)
    ! Variable declaration
    implicit none
    integer:: errorflag,nrsubst
    character(255):: CAS_list,pathforsingles,hc_folder
    logical:: CAS_list_exists,hc_exists
    type(file_t) :: my_file
    character(20),dimension(1):: tmp_name,cas_nr
    type(return_pos_t):: sub_part
    character(50):: dummy
    !---------------------------------------------------------------------------------------------------

    tmp_name = ''
    CAS_list_exists = .false.

    !If Cas_list_exists is false: get cas.txt by file path
    if(.not.hc_exists) then
        CAS_list = trim (pathforsingles) // 'cas.txt'
        inquire (file = CAS_list, exist = CAS_list_exists)
    else
        !here HC Case
        tmp_name = 'cas'
        CAS_list = pathforsingles
        CAS_list_exists = .true.
    end if

    !If a file could be found or HC is present
    if(CAS_list_exists) then
        call get_content(my_file,CAS_list,tmp_name,hc_folder)
        !cas_nr = gl%components(nrsubst)
        sub_part = get_sub_string_lines(my_file,cas_nr,errorflag)
        if(errorflag /= 0) return
        read(sub_part%p(1),*) dummy,cas_to_fld
    else
        !Write out error
        errorflag = -7884
        return
    end if

    end function cas_to_fld
    ! #####################################################################################################################################

    end module
    ! - end of module collecting all the file infos
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Written by Sven Pohl & Monika Thol, Lehrstuhl fuer Thermodynamik RUB
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Module for reading thermodynamic models
    ! The models are identified by a equation indicator (eos_type) and if a mixture is present by an mix_type
    ! EOS_TYPES: 1 | 7        -> Helmholtz models
    !            2 | 3 | 4    -> Cubic models
    !            6            -> PC-SAFT model
    !            51 | 52 | 53 -> gen. Helmholtz models
    !            8            -> RKM model
    !            9            -> Costald model
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    module sub_file_input_module

    ! Load modules
    use module_all_types
    use sub_file_input_helper
    use helmholtz_reader
    use other_reader
    use ge_reader
    use file
    use file_info

    implicit none

    contains

    ! #####################################################################################################################################
    ! Main read function for reading pure fluid model parameters and info
    ! Input: gl<type_gl>: global variable structure
    !        path<character>: absolut path to the model
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    ! -7878 - Fluid file not found / path wrong
    !-------------------------------------------------------------------------------------
    subroutine read_from_fluidfiles (gl,path, errorflag)
    !-------------------------------------------------------------------------------------
    type(type_gl) :: gl
    character(*), intent(in):: path
    integer:: errorflag
    !-------------------------------------------------------------------------------------
    type(file_info_t) :: f_info
    integer:: nrsubst, j
    !-------------------------------------------------------------------------------------
    ! end of declaration section

    ! 1) Collect all file infos
    f_info = file_info_t(gl,path,errorflag)

    if(errorflag .eq. 0) then

        ! 2) Set the paths for the models
        call f_info%set_paths(gl,errorflag)

    end if


    ! 3) Check if an error occured in info collection
    ! Read the files for each fluid / model
    if(errorflag .eq. 0) then

        ! -------------------------------------------------------------------------------------------------------------
        ! Pure fluid section ------------------------------------------------------------------------------------------
        ! -------------------------------------------------------------------------------------------------------------
        do nrsubst=1,gl%ncomp

            ! Helmholtz equations of state
            if((gl%Eq_type(nrsubst) == 1) ) then

                call read_helmholtz(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)

                ! Cubic equations of state
            elseif(gl%Eq_type(nrsubst) == 2 .or. gl%Eq_type(nrsubst) == 3 .or. gl%Eq_type(nrsubst) == 4 .or. (gl%Eq_type(nrsubst) == 7)) then

                call read_simples(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)

                ! Generalized Helmholtz equation of state
            elseif(gl%Eq_type(nrsubst) == 51 .or. gl%Eq_type(nrsubst) == 52 .or. gl%Eq_type(nrsubst) == 53) then

                call read_gen_eq(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)

                ! PC-SAFT equation of state
            elseif(gl%Eq_type(nrsubst) == 6) then

                call read_saft(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)

                ! RKM equation of state
            elseif (gl%Eq_type(nrsubst) == 8 .or. gl%Eq_type(nrsubst) == 81) then

                call read_simples(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)
                call read_rkm(gl,f_info%list_pure(nrsubst,2),f_info%hc_list_pure(nrsubst,2),nrsubst,f_info%is_hc,errorflag)

                ! Costald equation of state
            elseif(gl%Eq_type(nrsubst) == 9) then

                call read_simples(gl,f_info%list_pure(nrsubst,1),f_info%hc_list_pure(nrsubst,1),nrsubst,f_info%is_hc,errorflag)
                call read_costald(gl,f_info%list_pure(nrsubst,2),f_info%hc_list_pure(nrsubst,2),nrsubst,f_info%is_hc,errorflag)

            end if
        end do


        ! -------------------------------------------------------------------------------------------------------------
        ! Mix fluid section -------------------------------------------------------------------------------------------
        ! -------------------------------------------------------------------------------------------------------------
        If (gl%ncomp > 1 .and. errorflag .eq. 0) then !read mixingcoefficients if more than one component:

            !Loop over all components
            do nrsubst=1,gl%ncomp

                ! Check if no error is present
                if(errorflag .eq. 0) then

                    if (gl%components(nrsubst+1) /= '') then !'' means the current substance is also the last one

                        do j = nrsubst+1,gl%ncomp !always compare current substance to next ones
                            
                            if(errorflag .eq. 0) then
                                
                                if(gl%Mix_type == 1  .or. gl%Mix_type == 110 .or. &
                                    & gl%Mix_type == 111 .or. gl%Mix_type == 120 .or. gl%Mix_type == 121) then

                                    if(gl%Mix_type == 110 .or. gl%Mix_type == 120 ) then

                                        ! set default mixing rule 
                                        call set_default_mixing_rules(gl,nrsubst,j)

                                    else

                                        !check if binary mix files are available
                                        call read_helmholtz_mix(gl,nrsubst,j,f_info%mix_path, f_info%hc_mix_path, f_info%is_hc, errorflag)

                                        ! set default mixing rule if error occured
                                        if (errorflag /= 0) then
                                            call set_default_mixing_rules(gl,nrsubst,j)
                                        endif


                                    endif

                                elseif(gl%Mix_type == 12 .or. gl%Mix_type == 22) then

                                    call read_unifac_mix(gl, nrsubst, j, f_info%mix_path, f_info%hc_mix_path, f_info%is_hc, errorflag)

                                elseif (gl%Mix_type == 13) then

                                    call read_cosmo_mix(gl, nrsubst, j, f_info%mix_path, f_info%hc_mix_path, f_info%is_hc, errorflag)

                                elseif (gl%Mix_type == 2 .or. gl%Mix_type == 3 .or. gl%Mix_type == 4 .or. gl%Mix_type == 6) then

                                    call read_simples_mix(gl, nrsubst, j, f_info%mix_path, f_info%hc_mix_path,f_info%is_hc, errorflag)

                                endif
                            endif
                        end do
                    endif
                end if ! end checking error
            end do  ! end of looping gl%ncomp
            ! END OF MIXING LOOP
            if ((gl%Mix_type == 110).or.(gl%Mix_type == 111).or.(gl%Mix_type == 120).or.(gl%Mix_type == 121)) then
                gl%mix_type = 1
            endif
        endif
    end if

    end subroutine read_from_fluidfiles
    ! #####################################################################################################################################


    ! #####################################################################################################################################
    ! default parameters for mixing rules (helmholtz energy equations)
    !-------------------------------------------------------------------------------------
    subroutine set_default_mixing_rules(gl,nrsubst,j)

    type(type_gl) :: gl
    integer :: nrsubst,j

    !110, 120: special cases => all components are mixed according to linear or LB mixing rule
    if (gl%Mix_type == 110 .or. gl%Mix_type == 111) then
        !All components are mixed according to the linear mixing rule
        gl%rfbetat(nrsubst,j) = 1.d0
        gl%rfbetarho(nrsubst,j) = 1.d0
        gl%rfbetat(j,nrsubst) = 1.d0
        gl%rfbetarho(j,nrsubst) = 1.d0
        gl%rfgammat(nrsubst,j) = 0.5d0*(gl%tc(nrsubst)+gl%tc(j))/dsqrt(gl%tc(nrsubst)*gl%tc(j))
        gl%rfgammat(j,nrsubst) = gl%rfgammat(nrsubst,j)
        gl%rfgammarho(nrsubst,j) = 4.d0*(1.d0/gl%rhoc(nrsubst)+1.d0/gl%rhoc(j))/(gl%rhoc(nrsubst)**(-1.d0/3.d0)+gl%rhoc(j)**(-1.d0/3.d0))**3
        gl%rfgammarho(j,nrsubst) = gl%rfgammarho(nrsubst,j)

        !All components are mixed according to the Lorentz-Berthelot mixing rule.
    elseif (gl%Mix_type == 120 .or. gl%Mix_type == 121) then
        gl%rfbetat(nrsubst,j) = 1.d0
        gl%rfbetarho(nrsubst,j) = 1.d0
        gl%rfgammat(nrsubst,j) = 1.d0
        gl%rfgammarho(nrsubst,j) = 1.d0
        gl%rfbetat(j,nrsubst) = 1.d0
        gl%rfbetarho(j,nrsubst) = 1.d0
        gl%rfgammat(j,nrsubst) = 1.d0
        gl%rfgammarho(j,nrsubst) = 1.d0
    end if

    end subroutine
    ! #####################################################################################################################################

    end  module sub_file_input_module
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


