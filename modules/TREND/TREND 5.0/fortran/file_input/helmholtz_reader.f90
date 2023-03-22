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
    
    module helmholtz_reader
    use module_all_types
    use sub_file_input_helper
    use calc_functions
    use file
    use initialize_module

    implicit none
    contains
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the head of a fluid file
    ! Input: gl<type_gl>: global variable structure
    !        path<character>: absolut path to the model
    !        folder<character>: name of the folder
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helmholtz(gl,path,folder,nrsubst,is_hc,errorflag)
    implicit none
    character(*) :: path,folder
    type(type_gl) :: gl
    type(file_t) :: my_file
    type(return_pos_t):: sub_part,sub_part_end
    integer:: nrsubst,errorflag,loop_aux
    character(256),dimension(:),allocatable::head,EOS
    character(10), dimension(2):: EOS_name
    character(:),dimension(:),allocatable:: dummy
    logical :: is_hc ! info if file is hc or normal path file

    errorflag = 0

    ! Init the file variable
    my_file = file_t(is_hc)

    !get the content of the fluid file
    call get_content(my_file,path,(/gl%components(nrsubst)/),folder)

    if(my_file%errorflag .ne. 0 ) then
        errorflag = my_file%errorflag
        return
    end if

    !       ############    EOS INFO SECTION  #####################
    !--------------------------------------------------------------
    !Get the head of the file and read it in gl
    sub_part = get_empty_string_lines(my_file)
    head = get_lines(my_file,1,sub_part%pos_id(1)-1)
    call read_helm_head(gl,head,nrsubst)
    !--------------------------------------------------------------


    !       ############    EOS PARAMETERS  #######################
    !--------------------------------------------------------------
    !read equation of state
    EOS_name = '#EOS'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)

    if(errorflag == 0) then

        if(sub_part%pos_id(1) .ne. sub_part%pos_id(2)) then
            errorflag = -8000
            return
        end if

        if(.not.allocated(gl%litref)) allocate(gl%litref)
        gl%litref%lit_ref_res(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_res(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_EOS(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------

    !       ############    IDEAL AUX ############################
    !--------------------------------------------------------------
    if(.not.allocated(gl%litref)) allocate(gl%litref)

    EOS_name(1) = '#AUX'
    EOS_name(2) = '@AUX'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)

    if(errorflag == 0) then
        !cases for the different ideal models
        do loop_aux=1,size(sub_part%pos_id)
            if(gl%eos_coeff%cppcheck(nrsubst)(1:3) .eq.  my_file%content(sub_part%pos_id(loop_aux)+1)(1:3)) then
                gl%litref%lit_ref_ideal(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(loop_aux),len(gl%litref%lit_ref_ideal(nrsubst)))
                EOS = get_content_until_empty_line(my_file,sub_part%pos_id(loop_aux))
                !delete all comment lines
                call delete_lines(my_file,EOS)
                call read_helm_cp0(gl,EOS,nrsubst,errorflag)
                exit
            end if
        end do
    end if
    !--------------------------------------------------------------
    !       ############    VISCOSITY  ############################
    !--------------------------------------------------------------
    EOS_name = '#ETA'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)

    if(errorflag == 0) then

        if(sub_part%pos_id(1) .ne. sub_part%pos_id(2)) then
            errorflag = -8000
            return
        end if

        gl%litref%lit_ref_eta(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_eta(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_eta(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    VISCOSITY AUX ############################
    !--------------------------------------------------------------
    if (errorflag == 0) then
        EOS_name(1) = '#AUX'
        EOS_name(2) = '@AUX'
        sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
        if(errorflag == 0) then
            do loop_aux=1,size(sub_part%pos_id)
                if(gl%visco%pointereta(nrsubst)(1:3) .eq.  my_file%content(sub_part%pos_id(loop_aux)+1)(1:3)) then
                    gl%litref%lit_ref_eta_aux(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(loop_aux),len(gl%litref%lit_ref_ideal(nrsubst)))
                    EOS = get_content_until_empty_line(my_file,sub_part%pos_id(loop_aux))
                    !delete all comment lines
                    call delete_lines(my_file,EOS)
                    call read_helm_eta_aux(gl,EOS,nrsubst,errorflag)
                    exit
                end if
            end do
        end if
    end if


    !--------------------------------------------------------------
    !       ############    THERMAL CONDUCTIVITY  ############################
    !--------------------------------------------------------------
    EOS_name = '#TCX'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)

    if(errorflag == 0) then

        if(sub_part%pos_id(1) .ne. sub_part%pos_id(2)) then
            errorflag = -8000
            return
        end if

        gl%litref%lit_ref_tcx(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_eta(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_tcx(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    THERMAL AUX CONDUCTIVITY  ############################
    !--------------------------------------------------------------
    if (errorflag == 0) then
        EOS_name(1) = '#AUX'
        EOS_name(2) = '@AUX'
        sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
        if(errorflag == 0) then
            do loop_aux=1,size(sub_part%pos_id)
                if(gl%tcx%tkmodel(nrsubst)(1:3) .eq.  my_file%content(sub_part%pos_id(loop_aux)+1)(1:3)) then
                    gl%litref%lit_ref_tcx_aux(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(loop_aux),len(gl%litref%lit_ref_ideal(nrsubst)))
                    EOS = get_content_until_empty_line(my_file,sub_part%pos_id(loop_aux))
                    !delete all comment lines
                    call delete_lines(my_file,EOS)
                    call read_helm_tcx_aux(gl,EOS,nrsubst,errorflag)
                    exit
                end if
            end do
        end if
    end if


    !--------------------------------------------------------------
    !       ############    PMELT AUX ############################
    !--------------------------------------------------------------
    EOS_name = '#MLT'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_mlt(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_mlt(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_mlt(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    SUBL AUX ############################
    !--------------------------------------------------------------
    EOS_name = '#SBL'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        if(.not.allocated(sub_part%pos_id)) then
            write(*,*) 'error', errorflag
            pause
        end if
        gl%litref%lit_ref_sbl(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_sbl(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_sbl(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    VAP PRESS AUX   #######################
    !--------------------------------------------------------------
    EOS_name = '#PS'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_pv(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_pv(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_pv(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    DL AUX   #######################
    !--------------------------------------------------------------
    EOS_name = '#DL'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_dl(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_dl(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_dl(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    DV AUX   #######################
    !--------------------------------------------------------------
    EOS_name = '#DV'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_dl(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_dv(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_dv(gl,EOS,nrsubst,errorflag)
    end if
    !--------------------------------------------------------------
    !       ############    SURFACE TENSION AUX  ##################
    !--------------------------------------------------------------
    EOS_name = '#STN'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_stn(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_stn(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_stn(gl,EOS,nrsubst,errorflag)
    end if

    !--------------------------------------------------------------
    !       ############    DIELEC. AUX          ##################
    !--------------------------------------------------------------
    EOS_name = '#DE'
    sub_part = get_sub_string_lines(my_file,EOS_name,errorflag)
    if(errorflag == 0) then
        gl%litref%lit_ref_de(nrsubst) = get_ref_from_sub_content(my_file,sub_part%pos_id(1),len(gl%litref%lit_ref_de(nrsubst)))
        EOS = get_content_until_empty_line(my_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(my_file,EOS)
        call read_helm_de(gl,EOS,nrsubst,errorflag)
    end if

    if(errorflag /= -7878) then
        errorflag = 0
    end if

    if(allocated(head)) deallocate(head)
    if(allocated(head)) deallocate(EOS)
    if(allocated(head)) deallocate(dummy)

    end subroutine read_helmholtz
    !#################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the head of a fluid file
    ! Input: gl<type_gl>: global variable structure
    !        head<character array>: information block of fluid file
    !        nrsubst<integer>: Number of actual component
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_head(gl,head,nrsubst)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::head
    double precision , dimension(gl%ncomp):: wmcur,ttpcur,tnbpcur,tccur,pccur,rhoccur,accencur
    integer:: nrsubst
    ! Body of read_info_section
    !allocate divers fluid names
    !DEC$ IF DEFINED(HCgen)
    !DEC$ ELSE
    !if (.not.allocated(gl%substshortname)) then
    !    allocate(gl%substshortname(gl%ncomp))
    !    allocate(gl%substfullname(gl%ncomp))
    !    allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
    !    allocate(gl%substcasnr,mold=gl%substshortname)
    !else
    !    if (size(gl%substshortname,1) /= gl%ncomp) then
    !        deallocate(gl%substshortname,gl%substfullname,gl%substchemname,gl%substsynonym,gl%substcasnr)
    !        allocate(gl%substshortname(gl%ncomp))
    !        allocate(gl%substfullname(gl%ncomp))
    !        allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
    !        allocate(gl%substcasnr,mold=gl%substshortname)
    !    endif
    !end if
    !DEC$ ENDIF
    !---------------------------------------------

    !first paragraph till UN-number:
    !1. Name of Substance:
    read(head(1),*,end = 9999) gl%substshortname(nrsubst)
    read(head(2),*,end = 9999) gl%substcasnr(nrsubst)
    read(head(3),*,end = 9999) gl%substfullname(nrsubst)
    read(head(4),*,end = 9999) gl%substchemname(nrsubst)
    read(head(5),*,end = 9999) gl%substsynonym(nrsubst)
    !2. Current Properties of Substance:
    read(head(6),*,end = 9999)  wmcur(nrsubst)             !current molecular weight [g/mol]
    read(head(7),*,end = 9999)  ttpcur(nrsubst)            !current temperature at triple point [K]
    read(head(8),*,end = 9999)  tnbpcur(nrsubst)           !current temperature at normal boiling point [K]
    read(head(9),*,end = 9999)  tccur(nrsubst)             !current temperature at critical point [K]
    read(head(10),*,end = 9999) pccur(nrsubst)             !current pressure at critical point [kPa]
    !pccur(nrsubst)*factor_inv   ![MPa]
    read(head(11),*,end = 9999) rhoccur(nrsubst)           !current density at critical point [mol/L]
    !rhoccur(nrsubst)=rhoccur(nrsubst)*gl%Factor ![mol/m³]
    read(head(12),*,end = 9999) accencur(nrsubst)          !current acentric factor [-]
    read(head(13),*,end = 9999) gl%dipole(nrsubst)            !dipole moment [debye]
    !3. Reference State:
    read(head(14),*) gl%refstate(nrsubst)           !default reference state
    if (gl%refstate(nrsubst)(1:2) == 'OT') then !if 'OTher' reference state then read in reference T,p,h,s
        read(head(15),*)  gl%tref(nrsubst),gl%pref(nrsubst),gl%href(nrsubst),gl%sref(nrsubst)
        gl%pref(nrsubst) = gl%pref(nrsubst)/gl%factor        ![MPa]
    elseif (gl%refstate(nrsubst)(1:3) == 'NBP') then
        gl%tref(nrsubst) = 0.d0 !later in this subroutine set to normal boiling point temperature
        gl%pref(nrsubst) = 0.101325d0 !always 1 atm = 0.101325 MPa
        gl%href(nrsubst) = 0.d0 ![J/mol]
        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
    elseif (gl%refstate(nrsubst)(1:3) == 'ASH') then
        gl%tref(nrsubst) = 233.15d0
        gl%pref(nrsubst) = -1.d0 !"-1" shows that pref must be calculated later in another subroutine
        gl%href(nrsubst) = 0.d0 ![J/mol]
        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
    elseif (gl%refstate(nrsubst)(1:3) == 'IIR') then !Pay attention when using Lennard Jones Units
        gl%tref(nrsubst) = 273.15d0
        gl%pref(nrsubst) = -1.d0  !"-1" shows that pref must be calculated later in another subroutine
        gl%href(nrsubst) = 200.d3 ![J/kg]     -later changed with wm into [J/mol]
        gl%sref(nrsubst) = 1.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
    else !set to default values
        gl%tref(nrsubst) = 0.d0
        gl%pref(nrsubst) = 0.d0
        gl%href(nrsubst) = 0.d0 ![J/mol]
        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
    endif

    !                    !In the fluids files there is additional information given (e.g., version number, family, etc.).
    !                    !These information is not needed here so that it is not read in.

9999 return

    end subroutine read_helm_head
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of the equation of state
    ! Input: gl<type_gl>: global variable structure
    !        EOS<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_EOS(gl,EOS,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::EOS
    character(10):: eos_type
    integer:: errorflag
    integer:: nrsubst,err_read,loop,ibwr,lbwr
    integer:: j,m,dummycounter
    errorflag = 0

    !read the equation type
    read(EOS(2),*,end= 9999) eos_type


    !Large decision block for different
    !check for equation types
    if ((eos_type(1:2) == 'FE')) then
        if (eos_type(1:3) == 'FEH') gl%eos_coeff%hard_sphere(nrsubst) = .true.
        !data paragraph
        read(EOS(3),*,end= 9999) gl%tminfluid(nrsubst)          !lower temperature limit [K]
        read(EOS(4),*,end= 9999) gl%tmaxfluid(nrsubst)          !upper temperature limit [K]
        read(EOS(5),*,end= 9999) gl%pmaxfluid(nrsubst)          !upper pressure limit [kPa]
        gl%pmaxfluid(nrsubst) = gl%pmaxfluid(nrsubst)/gl%Factor    ![MPa]
        read(EOS(6),*,end= 9999) gl%rhomaxfluid(nrsubst)        !upper density limit [mol/L]
        gl%rhomaxfluid(nrsubst) = gl%rhomaxfluid(nrsubst)*gl%Factor ![mol/m³]
        read(EOS(7),*,end= 9999) gl%eos_coeff%cppcheck(nrsubst)        !pointer to Cp0 model
        read(EOS(8),*,end= 9999) gl%wm(nrsubst)                 !molecular weight [g/mol]
        gl%wm(nrsubst) =  gl%wm(nrsubst)/gl%Factor

        if (gl%refstate(nrsubst)(1:3) == 'IIR') then
            gl%href(nrsubst) = gl%href(nrsubst)*gl%wm(nrsubst)![j/mol]
            gl%sref(nrsubst) = gl%sref(nrsubst)*gl%wm(nrsubst)![j/(mol*k)]
        endif

        read(EOS(9),*) gl%ttp(nrsubst)          !triple point temperature

        !     (if "&" occurs in the input variable)
        if ((.not. gl%hold_limits).and.(.not. gl%check_solid)) then
            gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the vle curves
        end if

        read(EOS(10),*) gl%ptp(nrsubst)          !pressure at triple point
        gl%ptp(nrsubst) = gl%ptp(nrsubst)/gl%Factor  ![mpa]
        read(EOS(11),*) gl%rhotp(nrsubst)        !density at triple point
        gl%rhotp(nrsubst) = gl%rhotp(nrsubst)*gl%factor ![mol/m³]
        read(EOS(12),*) gl%tnbp(nrsubst)         !temperature at normal boiling point

        if (gl%refstate(nrsubst)(1:3) == 'NBP') then
            gl%tref(nrsubst) = gl%tnbp(nrsubst) !could not be set to normal boiling point temperature before
        endif

        read(EOS(13),*) gl%accen(nrsubst)        !acentric factor
        read(EOS(14),*) gl%tc(nrsubst),gl%pc(nrsubst),gl%rhoc(nrsubst)   !critical parameters
        gl%pc(nrsubst) = gl%pc(nrsubst)/gl%Factor       ![mpa]
        gl%rhoc(nrsubst) = gl%rhoc(nrsubst)*gl%factor   ![mol/m³]
        read(EOS(15),*) gl%tred(nrsubst),gl%rhored(nrsubst)  !reducing parameters
        gl%rhored(nrsubst) = gl%rhored(nrsubst)*gl%factor ![mol/m³]
        read(EOS(16),*) gl%req(nrsubst)          !gas constant used in fit


        !Read the parameters for the EOS
        !data columns
        err_read = 0
        read(EOS(17),*,iostat=err_read) gl%eos_coeff%nreg(nrsubst),gl%eos_coeff%nlreg(nrsubst), gl%eos_coeff%ncrt(nrsubst),gl%eos_coeff%nlcrt(nrsubst), gl%eos_coeff%nna(nrsubst), dummycounter, gl%eos_coeff%nlna(nrsubst)
        if (err_read /= 0) then
            read(EOS(17),*,iostat=err_read) gl%eos_coeff%nreg(nrsubst),gl%eos_coeff%nlreg(nrsubst), gl%eos_coeff%ncrt(nrsubst),gl%eos_coeff%nlcrt(nrsubst), gl%eos_coeff%nna(nrsubst), gl%eos_coeff%nlna(nrsubst)
        endif

        !  the gamma term is a multiplier for the (rho or del) in only the exponential
        !  terms; it is needed for e.g. the Bender EOS; set to 1.0 if not present

        do j=1,gl%eos_coeff%nreg(nrsubst) !for normal terms
            !p_i_double=0d0
            if (gl%eos_coeff%nlreg(nrsubst) == 5) then ! if 5 columns are stated for regular terms, then gama = 1 (except for BWR equations)
                read(EOS(17+j),*) gl%eos_coeff%ni(j,nrsubst),gl%eos_coeff%ti(j,nrsubst),gl%eos_coeff%di(j,nrsubst),gl%eos_coeff%p_i(j,nrsubst),gl%eos_coeff%gama(j,nrsubst)
                !gl%eos_coeff%p_i(j,nrsubst) = int(p_i_double)
            else
                read(EOS(17+j),*) gl%eos_coeff%ni(j,nrsubst),gl%eos_coeff%ti(j,nrsubst),gl%eos_coeff%di(j,nrsubst),gl%eos_coeff%p_i(j,nrsubst)
                !gl%eos_coeff%p_i(j,nrsubst) = int(p_i_double)
                if (gl%eos_coeff%p_i(j,nrsubst) /= 0) gl%eos_coeff%gama(j,nrsubst) = 1.d0
            endif
        end do



        if (gl%eos_coeff%ncrt(nrsubst) > 0) then
            do j=1,gl%eos_coeff%ncrt(nrsubst) !for critical terms
                read(EOS(17+gl%eos_coeff%nreg(nrsubst)+j),*) gl%eos_coeff%ni(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
                    gl%eos_coeff%ti(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
                    gl%eos_coeff%di(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
                    gl%eos_coeff%pli(j,nrsubst), &
                    gl%eos_coeff%tli(j,nrsubst), &
                    gl%eos_coeff%eta(j,nrsubst), &
                    gl%eos_coeff%beta(j,nrsubst), &
                    gl%eos_coeff%gam(j,nrsubst), &
                    gl%eos_coeff%eps(j,nrsubst), &
                    gl%eos_coeff%etana(j,nrsubst), &!!
                    gl%eos_coeff%eidna(j,nrsubst), &!!
                    gl%eos_coeff%eitna(j,nrsubst)!!
            enddo
        endif

        if (gl%eos_coeff%nna(nrsubst) > 0) then
            do j=1,gl%eos_coeff%nna(nrsubst) !for special gaussian terms
                read(EOS(17+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst)+j),*) gl%eos_coeff%ni(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%ti(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%di(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%pli(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%tli(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%eta(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%beta(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%gam(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%eps(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%etana(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%eidna(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
                    gl%eos_coeff%eitna(j+gl%eos_coeff%ncrt(nrsubst),nrsubst)
            enddo
        endif

        !                                tc_est_inp = gl%tc(nrsubst)
        !                                rhoc_est_inp = gl%rhoc(nrsubst)
        !                                gl%rhoredmix = 1.d0
        !                                ! find the number of gaussian bell-shaped terms in the equation
        gl%eos_coeff%I_GBS(nrsubst) = 0
        do m = 1, gl%eos_coeff%ncrt(nrsubst)
            if (gl%eos_coeff%etana(m,nrsubst) == 0.d0) gl%eos_coeff%I_GBS(nrsubst) = gl%eos_coeff%I_GBS(nrsubst) + 1
        end do

        !packing fraction parameters for hard sphere terms
        if (gl%eos_coeff%hard_sphere(nrsubst)) then
            read (EOS(17+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst)+gl%eos_coeff%nna(nrsubst)+1),*) gl%eos_coeff%h1(nrsubst), gl%eos_coeff%h2(nrsubst), gl%eos_coeff%h3(nrsubst), gl%eos_coeff%h4(nrsubst), gl%eos_coeff%dof2(nrsubst)
        end if
        !END OF STANDARD HELMHOLTZ EOS *******************************************************************************************************

    elseif (eos_type(1:3) == 'ECS') then
        errorflag = -7999	!chosen file has wrong format
        return

    elseif ((eos_type(1:3) == 'BWR') .or. (eos_type(1:3) == 'BEN') .or. (eos_type(1:3) == 'STA'))  then
        if (eos_type(1:3) == 'BWR') then
            gl%bwr_read = .true.
            loop = 11
        elseif (eos_type(1:3) == 'BEN') then
            gl%ben_read = .true.
            loop = 7
        elseif (eos_type(1:3) == 'STA') then
            gl%sta_read = .true.
            loop = 5
        endif

        !reading parameters
        !data paragraph
        read(EOS(3),*) gl%tminfluid(nrsubst)                            !lower temperature limit [k]
        read(EOS(4),*) gl%tmaxfluid(nrsubst)                            !upper temperature limit [k]
        read(EOS(5),*) gl%pmaxfluid(nrsubst)                            !upper pressure limit [kpa]
        gl%pmaxfluid(nrsubst) = gl%pmaxfluid(nrsubst)/gl%Factor      ![mpa]
        read(EOS(6),*) gl%rhomaxfluid(nrsubst)                          !upper density limit [mol/l]
        gl%rhomaxfluid(nrsubst) = gl%rhomaxfluid(nrsubst)*gl%factor  ![mol/m³]
        read(EOS(7),*) gl%eos_coeff%cppcheck(nrsubst)                             !pointer to cp0 model
        read(EOS(8),*) gl%wm(nrsubst)                                   !molecular weight [g/mol]
        gl%wm(nrsubst) = gl%wm(nrsubst)/gl%factor                       !!molecular weight [kg/mol]
        read(EOS(9),*) gl%ttp(nrsubst)                                  !triple point temperature
        !     (if "&" occurs in the input variable)
        if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the vle curves
        read(EOS(10),*) gl%ptp(nrsubst)                                  !pressure at triple point
        gl%ptp(nrsubst) = gl%ptp(nrsubst)/gl%Factor                  ![mpa]
        read(EOS(11),*) gl%rhotp(nrsubst)                                !density at triple point
        gl%rhotp(nrsubst) = gl%rhotp(nrsubst)*gl%factor              ![mol/m³]
        read(EOS(12),*) gl%tnbp(nrsubst)                                 !temperature at normal boiling point
        read(EOS(13),*) gl%accen(nrsubst)                                !acentric factor
        read(EOS(14),*) gl%tc(nrsubst),gl%pc(nrsubst),gl%rhoc(nrsubst)         !critical parameters
        gl%pc(nrsubst) = gl%pc(nrsubst)/gl%Factor                    ![mpa]
        gl%rhoc(nrsubst) = gl%rhoc(nrsubst)*gl%factor                ![mol/m³]
        read(EOS(15),*) gl%tred(nrsubst),gl%rhored(nrsubst)                 !reducing parameters
        gl%rhored(nrsubst) = gl%rhored(nrsubst)*gl%factor            ![mol/m³]
        read(EOS(16),*) gl%gama_bwr(nrsubst)
        read(EOS(17),*) gl%req(nrsubst)                                  !gas constant used in fit [l-bar/mol-k]
        gl%req(nrsubst) = gl%req(nrsubst)*gl%factorrbwr                  ![j/mol-k]
        read(EOS(18),*) gl%eos_coeff%nreg(nrsubst)                                 !number of terms in bwr equation

        lbwr = 0
        do ibwr=1, loop
            if (lbwr+ibwr < (gl%eos_coeff%nreg(nrsubst)-1.d0)) then
                read(EOS(18+ibwr),*) gl%ncoeff(nrsubst,ibwr+lbwr), gl%ncoeff(nrsubst, (ibwr+lbwr+1)), gl%ncoeff(nrsubst, (ibwr+lbwr+2))
            else
                if (gl%bwr_read) then
                    read(EOS(18+ibwr),*)  gl%ncoeff(nrsubst,ibwr+lbwr), gl%ncoeff(nrsubst, (ibwr+lbwr+1))
                elseif (gl%ben_read) then
                    read(EOS(18+ibwr),*)  gl%ncoeff(nrsubst,ibwr+lbwr)
                endif
            endif
            lbwr = lbwr+2
        enddo

        call bwr_recalc(gl,nrsubst)

    endif

    !gl%rhoredmix = gl%rhored(nrsubst) ! needed, because the module variable hasn't been set, yet
    !gl%Tredmix = gl%Tred(nrsubst)     ! needed, because the module variable hasn't been set, yet
    !gl%pc(nrsubst) = P_CALC(gl, gl%tc(nrsubst), gl%rhoc(nrsubst),nrsubst)  !critical pressure calculated from EOS [MPa]


9999 return

    end subroutine read_helm_EOS
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of cp0 model
    ! Input: gl<type_gl>: global variable structure
    !        IDEAL_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_cp0(gl,IDEAL_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::IDEAL_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    character(4)::cpdummy
    errorflag = 0

    read(IDEAL_AUX(2),*) eos_type

    !allocate ideal gas parameters
    if (.not.allocated(gl%cp0coeff)) then
        allocate(gl%cp0coeff(20,gl%ncomp))
        allocate(gl%cp0exp,gl%cp0hyp1,gl%cp0hyp2,gl%cp0hyp3,mold=gl%cp0coeff)
    else
        if (size(gl%cp0coeff,2) /= gl%ncomp) then
            deallocate(gl%cp0coeff,gl%cp0exp,gl%cp0hyp1,gl%cp0hyp2,gl%cp0hyp3)
            allocate(gl%cp0coeff(20,gl%ncomp))
            allocate(gl%cp0exp,gl%cp0hyp1,gl%cp0hyp2,gl%cp0hyp3,mold=gl%cp0coeff)
        endif

    end if
    !allocate divers fluid names
    !if (.not.allocated(gl%substshortname)) then
    !    allocate(gl%substshortname(gl%ncomp))
    !    allocate(gl%substfullname(gl%ncomp))
    !    allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
    !    allocate(gl%substcasnr,mold=gl%substshortname)
    !else
    !    if (size(gl%substshortname,1) /= gl%ncomp) then
    !        deallocate(gl%substshortname,gl%substfullname,gl%substchemname,gl%substsynonym,gl%substcasnr)
    !        allocate(gl%substshortname(gl%ncomp))
    !        allocate(gl%substfullname(gl%ncomp))
    !        allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
    !        allocate(gl%substcasnr,mold=gl%substshortname)
    !    endif
    !end if


    !Large decision block for different
    !check for equation types
    if (eos_type(1:2) == 'CJ' .and. (gl%eos_coeff%cppcheck(nrsubst)(1:2) == 'CJ')) then
        gl%cpmodel(nrsubst) = .true.

        !read reducing parameters for cp0:
        read(IDEAL_AUX(7),*) gl%tcp0red(nrsubst), gl%cp0red(nrsubst)
        !read coeffs ande exps:
        read(IDEAL_AUX(8),*) gl%ncp0poly(nrsubst),gl%ncp0pl(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
        do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
            read(IDEAL_AUX(8+j),*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
        enddo

        cpdummy=gl%eos_coeff%cppcheck(nrsubst)
        cpdummy(1:2) = 'CP'
        gl%cp0red(nrsubst) = 8.3144621d0
        gl%cp0coeff(1,nrsubst) = (gl%cp0coeff(1,nrsubst)-37.93d0)/gl%cp0red(nrsubst)
        gl%cp0coeff(2,nrsubst) = (gl%cp0coeff(2,nrsubst)+0.21d0)/gl%cp0red(nrsubst)
        gl%cp0coeff(3,nrsubst) = (gl%cp0coeff(3,nrsubst)-3.91d-4)/gl%cp0red(nrsubst)
        gl%cp0coeff(4,nrsubst) = (gl%cp0coeff(4,nrsubst)+2.06d-7)/gl%cp0red(nrsubst)

    elseif(eos_type(1:2) == 'CP' .and. (gl%eos_coeff%cppcheck(nrsubst)(1:2) == 'CP')) then

        gl%cpmodel(nrsubst) = .true.

        !read reducing parameters for cp0:
        read(IDEAL_AUX(7),*) gl%tcp0red(nrsubst), gl%cp0red(nrsubst)

        !read coeffs ande exps:
        read(IDEAL_AUX(8),*) gl%ncp0poly(nrsubst),gl%ncp0pl(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
        do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
            read(IDEAL_AUX(8+j),*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
        enddo

        if ((gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) > 0) then
            do j=1,(gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) !for hyperbolic
                read(IDEAL_AUX(8+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))+j),*) gl%cp0coeff(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
                    gl%cp0exp(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
                    gl%cp0hyp1(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
                    gl%cp0hyp2(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
                    gl%cp0hyp3(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst)
            enddo

        endif

    elseif ( ((eos_type(1:3) == 'PHK') .and. (gl%eos_coeff%cppcheck(nrsubst)(1:3) == 'PHK'))  ) then !modified by andreas july 2011 , modified by judith august 2011
        gl%phkmodel(nrsubst) = .true.
        gl%cpmodel(nrsubst) = .false.

        !set reference state and reducing parameters for cp0:
        gl%refstate(nrsubst)='OT0'
        gl%tref(nrsubst) = 298.15d0
        gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
        gl%href(nrsubst) = 0.d0 ![j/kg]     -later changed with wm into [j/mol]
        gl%sref(nrsubst) = 0.d3   ![j/(kg k)] -later changed with wm into [j/(mol k)]
        gl%tcp0red(nrsubst) = gl%tc(nrsubst)
        gl%cp0red(nrsubst) = 8.31451d0
        !read coeffs and exps:
        read(IDEAL_AUX(7),*) gl%ncp0log(nrsubst),gl%ncp0(nrsubst),gl%ncp0logexp(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
        !read polynomials and constants ci and cii
        do j=1,(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst))
            read(IDEAL_AUX(7+j),*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
        end do
        if ((gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) > 0) then
            do j=1,(gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) !for hyperbolic terms
                read(IDEAL_AUX(7+(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst))+j),*) gl%cp0coeff(j+(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst)),nrsubst), &
                    gl%cp0exp(j+(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst)),nrsubst)
            end do
        end if

    endif

    !change coefficients of the ideal gas if cp0red = 1.d0 (instead of the molar gas constant)
    if((gl%cp0red(nrsubst)-1.d0) < 1.d-8) then
        do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
            gl%cp0coeff(j,nrsubst)=gl%cp0coeff(j,nrsubst)/gl%req(nrsubst)
        end do
        gl%cp0red(nrsubst)=gl%req(nrsubst)
    end if


    end subroutine
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of melting line model
    ! Input: gl<type_gl>: global variable structure
    !        MLT_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_mlt(gl,MLT_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::MLT_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0

    if(.not.(allocated(MLT_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%pmeltcoeff)) then
        allocate(gl%pmeltcoeff(20,gl%ncomp))
        allocate(gl%pmeltexp, mold=gl%pmeltcoeff)
    else
        if (size(gl%pmeltcoeff,2) /= gl%ncomp) then
            deallocate(gl%pmeltcoeff,gl%pmeltexp,gl%substchemname,gl%substsynonym,gl%substcasnr)
            allocate(gl%pmeltcoeff(20,gl%ncomp))
            allocate(gl%pmeltexp, mold=gl%pmeltcoeff)
        endif
    end if

    read(MLT_AUX(2),*) eos_type
    gl%pmelttype(nrsubst) = ichar(eos_type(3:3)) - 48 !get the number at 3. place
    !paragraph with comments can be ignored (every line starts with "?")
    if (gl%pmelttype(nrsubst) == 39) then

        !read reducing parameters for melting pressure:
        read(MLT_AUX(3),*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
        gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)/gl%Factor ![mpa]

        !read coeffs and exponents:
        read(MLT_AUX(4),*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst)
        do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst))
            read(MLT_AUX(4+j),*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
        enddo
    elseif (gl%pmelttype(nrsubst) == 24) then                                        !fluid is water

        !paragraph with doubled parameters can be ignored (always 4 lines)
        read(MLT_AUX(3),*) gl%pmeltmintemp(nrsubst)
        read(MLT_AUX(4),*) gl%pmeltmaxtemp(nrsubst)

        !read reducing parameters for melting pressure:
        read(MLT_AUX(7),*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
        gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)/gl%Factor ![mpa]

        !read coeffs ande exps:
        read(MLT_AUX(8),*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst),&
            gl%npmeltcoeff4(nrsubst), gl%npmeltcoeff5(nrsubst)
        do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst)&
            +gl%npmeltcoeff4(nrsubst)+ gl%npmeltcoeff5(nrsubst))
            read(MLT_AUX(8+j),*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
        enddo
    elseif (gl%pmelttype(nrsubst) == 20) then                                        !fluid is heavy water
        !paragraph with doubled parameters can be ignored (always 4 lines)
        read(MLT_AUX(3),*) gl%pmeltmintemp(nrsubst)
        read(MLT_AUX(4),*) gl%pmeltmaxtemp(nrsubst)
        read(MLT_AUX(7),*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)

    else
        !paragraph with doubled parameters can be ignored (always 4 lines)
        read(MLT_AUX(3),*) gl%pmeltmintemp(nrsubst)
        read(MLT_AUX(4),*) gl%pmeltmaxtemp(nrsubst)
        read(MLT_AUX(7),*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
        !andreas, september 2014, i guess the 1.d-30 cant be right? replaced by 1.d-3
        !pmeltred(nrsubst) = pmeltred(nrsubst)*1.d-30 ![mpa]
        gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)/gl%Factor ![mpa]
        read(MLT_AUX(8),*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst),&
            gl%npmeltcoeff4(nrsubst), gl%npmeltcoeff5(nrsubst)
        do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst)&
            +gl%npmeltcoeff4(nrsubst)+ gl%npmeltcoeff5(nrsubst))
            read(MLT_AUX(8+j),*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
        enddo
    end if



    end subroutine


    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of viscosity aux model
    ! Input: gl<type_gl>: global variable structure
    !        ETA<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_eta_aux(gl,ETA_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::ETA_AUX
    character(20):: eos_type
    integer:: nrsubst, errorflag,k
    double precision, dimension(30) :: qc_inv, qd_inv   !input parameters viscosity
    character(255):: dummy
    errorflag = 0

    if(.not.(allocated(ETA_AUX))) return !return when given model is not present for current substance
    read(ETA_AUX(2),*)   eos_type


    if ((eos_type(1:3) == 'H2O') .AND. (gl%visco%pointer_lambda(nrsubst) == 'H2O')) then

        read(ETA_AUX(3),*) gl%visco%lowbo_temp_mue2(nrsubst)       !lower boundary [K] of region around the crit. Temp.
        read(ETA_AUX(4),*) gl%visco%uppbo_temp_mue2(nrsubst)       !upper boundary [K] of region around the crit. Temp.
        read(ETA_AUX(5),*) gl%visco%lowbo_dens_mue2(nrsubst)       !lower boundary	[kg/m^3] of region around the crit. Dens.
        read(ETA_AUX(6),*) gl%visco%uppbo_dens_mue2(nrsubst)       !upper boundary	[kg/m^3] of region around the crit. Dens.
        read(ETA_AUX(7),*) gl%visco%corr_ref_mue2(nrsubst)         !reference correlation length
        read(ETA_AUX(8),*) gl%visco%X_mue2(nrsubst)                !x_mue
        read(ETA_AUX(9),*) qc_inv(nrsubst)               !q_c^-1  [1/nm]
        gl%visco%qc_mue2(nrsubst)=1.d0/(qc_inv(nrsubst)*1.d-9)    ![m]
        read(ETA_AUX(10),*) qd_inv(nrsubst)               !q_d^-1  [1/nm]
        gl%visco%qd_mue2(nrsubst)=1.d0/(qd_inv(nrsubst)*1.d-9)    ![m]
        read(ETA_AUX(11),*) gl%visco%nue_mue2(nrsubst)              !nue
        read(ETA_AUX(12),*) gl%visco%gamma_mue2(nrsubst)            !gamma
        read(ETA_AUX(13),*) gl%visco%Xi0_mue2(nrsubst)              !Xi [nm]
        gl%visco%Xi0_mue2(nrsubst)=gl%visco%Xi0_mue2(nrsubst)*1.d-9  ![m]
        read(ETA_AUX(14),*) gl%visco%gamma0_mue2(nrsubst)           !GAMMA_0
        read(ETA_AUX(15),*) gl%visco%T_mue2(nrsubst)                !T_quer_R
        return


        !elseif ((eos_type(1:3) == 'I08') .AND. (gl%visco%pointerceaux(nrsubst) == 'I08')) then
    elseif (eos_type(1:3) == 'I08') then


        read(ETA_AUX(3),*) gl%visco%lowbo_temp_mue2(nrsubst)       !lower boundary [K] of region around the crit. Temp.
        read(ETA_AUX(4),*) gl%visco%uppbo_temp_mue2(nrsubst)       !upper boundary [K] of region around the crit. Temp.
        read(ETA_AUX(5),*) gl%visco%lowbo_dens_mue2(nrsubst)       !lower boundary	[kg/m^3] of region around the crit. Dens.
        read(ETA_AUX(6),*) gl%visco%uppbo_dens_mue2(nrsubst)       !upper boundary	[kg/m^3] of region around the crit. Dens.
        read(ETA_AUX(7),*) gl%visco%corr_ref_mue2(nrsubst)         !reference correlation length
        read(ETA_AUX(8),*) gl%visco%X_mue2(nrsubst)                !x_mue
        read(ETA_AUX(9),*) gl%visco%qc_mue2(nrsubst)               !q_c^-1
        read(ETA_AUX(10),*) gl%visco%qd_mue2(nrsubst)               !q_d^-1
        read(ETA_AUX(11),*) gl%visco%nue_mue2(nrsubst)              !nue
        read(ETA_AUX(12),*) gl%visco%gamma_mue2(nrsubst)            !gamma
        read(ETA_AUX(13),*) gl%visco%Xi0_mue2(nrsubst)              !Xi
        read(ETA_AUX(14),*) gl%visco%gamma0_mue2(nrsubst)           !GAMMA_0
        read(ETA_AUX(15),*) gl%visco%T_mue2(nrsubst)                !T_quer_R
        gl%visco%lowbo_dens_mue2(nrsubst)=gl%visco%lowbo_dens_mue2(nrsubst)*gl%Factor/gl%wm(nrsubst)
        gl%visco%uppbo_dens_mue2(nrsubst)=gl%visco%uppbo_dens_mue2(nrsubst)*gl%Factor/gl%wm(nrsubst)
        gl%visco%h2o_read= .true.
        return

    elseif ((eos_type(1:3) == 'CI1') .AND. (gl%visco%pointereta(nrsubst) == 'CI1')) then

        read(ETA_AUX(3),*) gl%visco%low_tempcoll(nrsubst)               !lower temperature limit [K]
        read(ETA_AUX(4),*) gl%visco%upp_tempcoll(nrsubst)               !upper temperature limit [K]
        read(ETA_AUX(5),*) gl%visco%low_prescoll(nrsubst)               !(dummy) upper pressure limit
        read(ETA_AUX(6),*) gl%visco%upp_denscoll(nrsubst)               !(dummy) maximum density
        read(ETA_AUX(7),*) gl%visco%colltermnr(nrsubst)               !number of terms
        do k = 1, gl%visco%colltermnr(nrsubst)          !coeff, power of Tstar
            read(ETA_AUX(7+k),*) gl%visco%coeffcoll(nrsubst,k),gl%visco%powercoll(nrsubst,k)
        end do
        gl%visco%coll_read=.true.
        return

    elseif ((dummy(1:3) == 'CI2') .AND. (gl%visco%pointereta(nrsubst) == 'CI2')) then

        read(ETA_AUX(3),*) gl%visco%low_tempcoll(nrsubst)               !lower temperature limit [K]
        read(ETA_AUX(4),*) gl%visco%upp_tempcoll(nrsubst)               !upper temperature limit [K]
        read(ETA_AUX(5),*) gl%visco%low_prescoll(nrsubst)               !(dummy) upper pressure limit
        read(ETA_AUX(6),*) gl%visco%upp_denscoll(nrsubst)               !(dummy) maximum density
        read(ETA_AUX(7),*) gl%visco%colltermnr(nrsubst)                 !number of terms
        do k = 1, gl%visco%colltermnr(nrsubst)                 !coeff, power of Tstar
            read(ETA_AUX(8+k),*) gl%visco%coeffcoll(nrsubst,k)
        end do
        gl%visco%coll_read=.true.
        return
    end if

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of viscosity model
    ! Input: gl<type_gl>: global variable structure
    !        ETA<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_eta(gl,ETA,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::ETA
    integer:: errorflag,nrsubst,j,last_sum,i,ieta,jeta,keta,leta,meta,neta,err_read
    double precision:: dummy1,dummy2
    character(20):: eos_type
    character(255) :: dummy
    errorflag = 0

    if(.not.(allocated(ETA))) then
        errorflag = -1 !model does not exist
        return !return when given model is not present for current substance
    endif
    if(.not.allocated(gl%visco)) then
        allocate(gl%visco)
        call initialize_viscosity(gl)
    end if
    gl%visco%eta_read = .true.
    last_sum = 0

    read(ETA(2),*) gl%visco%etamodel(nrsubst)

    !reading parameters for VS0 model specification
    if (gl%visco%etamodel(nrsubst) == 'VS0') then
        read(ETA(3),*) gl%visco%tmineta(nrsubst)    ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)    ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)   ![kPa]
        gl%visco%pmaxeta(nrsubst)=gl%visco%pmaxeta(nrsubst)/gl%Factor
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)    ![mol/L]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%factor
        read(ETA(7),*) gl%visco%pointer_lambda(nrsubst)     !Pointer to hardcoded model

        if (trim(gl%visco%pointer_lambda(nrsubst))  == 'LJF' .or. trim(gl%visco%pointer_lambda(nrsubst))  == 'LJ1' ) then
            return
        end if

        if((trim(gl%visco%pointer_lambda(nrsubst)) == 'HE') .or.(trim(gl%visco%pointer_lambda(nrsubst)) == 'ETY').or.(trim(gl%visco%pointer_lambda(nrsubst)) == 'NEO')) then
            return
        end if


        read(ETA(8),*) gl%visco%nterm1_eta(nrsubst), gl%visco%nterm2_eta(nrsubst), gl%visco%nterm3_eta(nrsubst), gl%visco%nterm4_eta(nrsubst), &
            & gl%visco%nterm5_eta(nrsubst), gl%visco%nterm6_eta(nrsubst), gl%visco%nterm7_eta(nrsubst), gl%visco%nterm8_eta(nrsubst)

        if (trim(gl%visco%pointer_lambda(nrsubst)) == 'R23') then
            read(ETA(9),*) gl%visco%pointereta(nrsubst)
            read(ETA(10),*) gl%visco%lejo_sigma(nrsubst)
            read(ETA(11),*) gl%visco%lejo_epka(nrsubst)
            read(ETA(12),*) gl%visco%red_temp_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)
            read(ETA(13),*) gl%visco%Chapman_Enskog1(nrsubst), gl%visco%Chapman_Enskog2(nrsubst)

            read(ETA(14),*) gl%visco%red_temp_eta(nrsubst) , gl%visco%red_dens_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)    !reducing parameters
            gl%visco%red_dens_eta(nrsubst)=gl%visco%red_dens_eta(nrsubst)*gl%Factor
            read(ETA(15),*) gl%visco%rhol_R23
            gl%visco%rhol_R23=gl%visco%rhol_R23*gl%factor
            read(ETA(16),*) gl%visco%c1_R23
            read(ETA(17),*) gl%visco%c2_R23
            read(ETA(18),*) gl%visco%delg_R23
            read(ETA(19),*) gl%visco%etamax_R23
            read(ETA(20),*) gl%visco%pointercrit_eta(nrsubst)
            return
        end if
        !end of reading R23 vis model

        read(ETA(9),*) gl%visco%red_temp_eta(nrsubst) , gl%visco%red_dens_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)    !reducing parameters
        gl%visco%red_dens_eta(nrsubst)=gl%visco%red_dens_eta(nrsubst)*gl%Factor

        do i = 1,gl%visco%nterm1_eta(nrsubst)   ! Loop for reduced viscosity at high density according to the Enskog theory
            read(ETA(9+last_sum),*) gl%visco%coeff_ei(nrsubst,i),gl%visco%exp_ei(nrsubst,i)
            last_sum = last_sum + 1
        end do

        do i = 1,gl%visco%nterm2_eta(nrsubst)   ! Loop for dilute-gas viscosity Mue0
            read(ETA(9+last_sum),*) gl%visco%coeff_dil(nrsubst,i),gl%visco%exp_dil(nrsubst,i)
            last_sum = last_sum + 1
        end do
        !Moni
        last_sum = 1
        do i = 1,(gl%visco%nterm3_eta(nrsubst)+gl%visco%nterm4_eta(nrsubst)) ! Loop for water and ethylene
            read(ETA(9+last_sum),*) gl%visco%coeff_hieta(nrsubst,i),dummy1,dummy2, gl%visco%exp1eta(nrsubst,i), gl%visco%exp2eta(nrsubst,i)
            gl%visco%iexpeta(nrsubst,i) = int(dummy1)
            gl%visco%jexpeta(nrsubst,i) = int(dummy2)
            !if (trim(gl%components(nrsubst)) == 'water') then
            if ((trim(gl%substcasnr(nrsubst)) == '7732-18-5') .or. (trim(gl%substcasnr(nrsubst)) == '7789-20-0')) then   !water and heavy water
                gl%visco%hieta(nrsubst, gl%visco%iexpeta(nrsubst,i)+1, gl%visco%jexpeta(nrsubst,i)+1)=gl%visco%coeff_hieta(nrsubst,i)
            end if
            last_sum = last_sum + 1
        end do

        do i = 1,gl%visco%nterm5_eta(nrsubst)   ! Loop for viscosity in the zero density limit
            read(ETA(9+last_sum),*) gl%visco%coeff_ai(nrsubst,i),gl%visco%exp_ai(nrsubst,i)
            last_sum = last_sum + 1
        end do

        do i = 1,gl%visco%nterm6_eta(nrsubst)   ! Loop for second viscosity virial coefficient
            read(ETA(9+last_sum),*) gl%visco%coeff_bi(nrsubst,i),gl%visco%exp_bi(nrsubst,i)
            last_sum = last_sum + 1
        end do

        do i = 1,gl%visco%nterm7_eta(nrsubst)   ! Loop for third viscosity virial coefficient
            read(ETA(9+last_sum),*) gl%visco%coeff_ci(nrsubst,i),gl%visco%exp_ci(nrsubst,i)
            last_sum = last_sum + 1
        end do

        do i = 1,gl%visco%nterm8_eta(nrsubst)   ! Loop for reduced viscosity at high density according to the Enskog theory
            read(ETA(9+last_sum),*) gl%visco%coeff_di(nrsubst,i),gl%visco%exp_di(nrsubst,i)
            last_sum = last_sum + 1
        end do

        read(ETA(9+last_sum),*) gl%visco%pointercrit_eta(nrsubst)


        !reading parameters for VS1 model specification
    elseif (gl%visco%etamodel(nrsubst) == 'VS1') then
        read(ETA(3),*) gl%visco%tmineta(nrsubst)    ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)    ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)   ![kPa]
        gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)/gl%factor
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
        read(ETA(7),*) gl%visco%dilgas_term(nrsubst)     !number of terms allocated with dilute-gas function
        read(ETA(8),*) gl%visco%pointereta(nrsubst)     !Points out which collision integral model is used
        read(ETA(9),*) gl%visco%slj(nrsubst)      !Lennard-Jones coefficient sigma [nm]
        read(ETA(10),*) gl%visco%eklj(nrsubst)       !Lennard-Jones coefficient epsilon/kappa [K]
        read(ETA(11),*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
        read(ETA(12),*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)   !ce (=0.0266958*sqrt(M)), cep (=Power of T)
        if((gl%substcasnr(nrsubst) == '64-17-5').or.(gl%substcasnr(nrsubst) == '306-83-2')) then !ethanol and R123
            read(ETA(13),*) gl%visco%a0_vs1(nrsubst)
            read(ETA(14),*) gl%visco%a1_vs1(nrsubst)
            read(ETA(15),*) gl%visco%a2_vs1(nrsubst)
            last_sum = 3
        end if
        read(ETA(13+last_sum),*) gl%visco%denstermnr(nrsubst)

        if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN
            last_sum = last_sum + 1
            read(ETA(13+last_sum),*) gl%visco%treddens(nrsubst), gl%visco%etaB2(nrsubst) !reducing parameters for T (= eps/k), etaB2 (= 0.6022137*sigma**3)
            DO jeta=1,gl%visco%denstermnr(nrsubst),1
                read(ETA(14+last_sum),*) gl%visco%coeffbetastar(nrsubst,jeta), gl%visco%powerbetastar(nrsubst,jeta)
                last_sum = last_sum + 1
            END DO
        end if
        read(ETA(14+last_sum),*) gl%visco%term_num1_eta(nrsubst), gl%visco%term_num2_eta(nrsubst), gl%visco%term_num3_eta(nrsubst), &
            & gl%visco%term_num4_eta(nrsubst), gl%visco%term_num5_eta(nrsubst), gl%visco%term_num6_eta(nrsubst)

        read(ETA(15+last_sum),*) gl%visco%tredeta(nrsubst) , gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst)    !reducing parameters
        gl%visco%rhoredeta(nrsubst)=gl%visco%rhoredeta(nrsubst)*gl%Factor

        if (gl%visco%term_num1_eta(nrsubst) /= 0) then
            do ieta = 1, abs(gl%visco%term_num1_eta(nrsubst)), 1
                read(ETA(16+last_sum),*) gl%visco%a_vs1(nrsubst,ieta), gl%visco%powertau(nrsubst,ieta)
                last_sum = last_sum + 1
            end do
        end if

        if (gl%visco%term_num2_eta(nrsubst) /= 0) then
            do jeta = 1,gl%visco%term_num2_eta(nrsubst),1
                read(ETA(16+last_sum),*) gl%visco%b_vs1(nrsubst,jeta), gl%visco%ptausp(nrsubst,jeta),gl%visco%pdelsp(nrsubst,jeta),gl%visco%pdel0sp(nrsubst,jeta),gl%visco%expdel(nrsubst,jeta)
                last_sum = last_sum + 1
            enddo
        end if
        if (gl%visco%term_num3_eta(nrsubst) /= 0) then
            do keta = gl%visco%term_num2_eta(nrsubst)+1,gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst),1
                read(ETA(16+last_sum),*) gl%visco%b_vs1(nrsubst,keta), gl%visco%ptausp(nrsubst,keta),gl%visco%pdelsp(nrsubst,keta),gl%visco%pdel0sp(nrsubst,keta),gl%visco%expdel(nrsubst,keta)
                last_sum = last_sum + 1
            enddo
        end if
        if (gl%visco%term_num4_eta(nrsubst) /= 0) then
            do leta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+1,&
                & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst),1
                read(ETA(16+last_sum),*) gl%visco%b_vs1(nrsubst,leta), gl%visco%ptausp(nrsubst,leta),gl%visco%pdelsp(nrsubst,leta),gl%visco%pdel0sp(nrsubst,leta),gl%visco%expdel(nrsubst,leta)
                last_sum = last_sum + 1
            enddo
        end if
        if (gl%visco%term_num5_eta(nrsubst) /= 0) then
            do meta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+1, &
                & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst),1
                read(ETA(16+last_sum),*) gl%visco%b_vs1(nrsubst,meta), gl%visco%ptausp(nrsubst,meta),gl%visco%pdelsp(nrsubst,meta),gl%visco%pdel0sp(nrsubst,meta),gl%visco%expdel(nrsubst,meta)
                last_sum = last_sum + 1
            enddo
        end if
        if (gl%visco%term_num6_eta(nrsubst) /= 0) then
            do neta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst)+1,&
                & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst)+gl%visco%term_num6_eta(nrsubst),1
                read(ETA(16+last_sum),*) gl%visco%b_vs1(nrsubst,neta), gl%visco%ptausp(nrsubst,neta),gl%visco%pdelsp(nrsubst,neta),gl%visco%pdel0sp(nrsubst,neta),gl%visco%expdel(nrsubst,neta)
                last_sum = last_sum + 1
            enddo
        end if

        !reading parameters for VS2 model specification
        !B. A. Younglove, J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982
    elseif (gl%visco%etamodel(nrsubst) == 'VS2') then
        read(ETA(3),*) gl%visco%tmineta(nrsubst)          ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)          ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)          ![kPa]
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)        ![mol/l]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
        read(ETA(7),*) gl%visco%pointereta(nrsubst)       !Points out which collision integral model is used
        read(ETA(8),*) gl%visco%slj(nrsubst)              !Lennard-Jones coefficient sigma [nm]
        read(ETA(9),*) gl%visco%eklj(nrsubst)             !Lennard-Jones coefficient epsilon/kappa [K]
        read(ETA(10),*) gl%visco%const19_VS2(nrsubst)      !const in Eq 19
        read(ETA(11),*) gl%visco%exp19_VS2(nrsubst)        !exponent in Eq 19 for T
        read(ETA(12),*) gl%visco%Fv_VS2(nrsubst,1)         !coefficients for initial density dependence of viscosity for eq. (21)

        do i=2,4
            read(ETA(12+i-1),*) gl%visco%Fv_VS2(nrsubst,i)
        end do
        read(ETA(17),*) gl%visco%Ev_VS2(nrsubst,1)         !coefficients for residual viscosity for eqs. (22) - (25)
        do i=2,8
            read(ETA(17+i-1),*) gl%visco%Ev_VS2(nrsubst,i)
        end do
        read(ETA(26),*) gl%visco%pointercrit_eta(nrsubst)  !pointer to critical enhancement auxiliary function (none used)
        !! unter #AUX CI2 könnten die Koeffizienten für Omega eingelesen werden

        !reading parameters for VS4 model specification
    elseif (gl%visco%etamodel(nrsubst) == 'VS4') then
        read(ETA(3),*) gl%visco%tmineta(nrsubst)    ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)    ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)   ![kPa]
        gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)/gl%factor
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
        read(ETA(7),*) gl%visco%dilgas_term(nrsubst)     !number of terms allocated with dilute-gas function
        read(ETA(8),*) gl%visco%pointereta(nrsubst)     !Points out which collision integral model is used
        read(ETA(9),*) gl%visco%slj(nrsubst)      !Lennard-Jones coefficient sigma [nm]
        read(ETA(10),*) gl%visco%eklj(nrsubst)       !Lennard-Jones coefficient epsilon/kappa [K]
        read(ETA(11),*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
        read(ETA(12),*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)   !ce (=0.0266958*sqrt(M)), cep (=Power of T)
        do i = 1, gl%visco%dilgas_term(nrsubst)-1
            read(ETA(12+i),*) gl%visco%d0_vs4(nrsubst,i), gl%visco%d0exp(nrsubst,i)
        end do
        last_sum = i
        read(ETA(12+last_sum),*) gl%visco%denstermnr(nrsubst)
        if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN !?????
            read(ETA(13+last_sum),*) gl%visco%treddens(nrsubst), gl%visco%etaB2(nrsubst) !reducing parameters for T (= eps/k), etaB2 (= 0.6022137*sigma**3)
            last_sum = last_sum + 1
            DO jeta=1,gl%visco%denstermnr(nrsubst),1
                read(ETA(13+last_sum),*) gl%visco%coeffbetastar(nrsubst,jeta), gl%visco%powerbetastar(nrsubst,jeta)
                last_sum = last_sum + 1
            END DO
        end if
        read(ETA(13+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%a_0(nrsubst), gl%visco%a_1(nrsubst), gl%visco%a_2(nrsubst)
        read(ETA(14+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%b0_vs4(nrsubst), gl%visco%b1_vs4(nrsubst), gl%visco%b2_vs4(nrsubst)
        read(ETA(15+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%c0(nrsubst), gl%visco%c_1(nrsubst), gl%visco%c_2(nrsubst)
        read(ETA(16+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%aa0(nrsubst), gl%visco%aa1(nrsubst), gl%visco%aa2(nrsubst)
        read(ETA(17+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%bb0(nrsubst), gl%visco%bb1(nrsubst), gl%visco%bb2(nrsubst)
        read(ETA(18+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%cc0(nrsubst), gl%visco%cc1(nrsubst), gl%visco%cc2(nrsubst)
        read(ETA(19+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%dd_0(nrsubst), gl%visco%dd_1(nrsubst), gl%visco%dd_2(nrsubst)
        read(ETA(20+last_sum),'(A255)')dummy
        if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%e_0(nrsubst), gl%visco%e_1(nrsubst), gl%visco%e_2(nrsubst)

        !reading parameters for VS5 model specification
        !T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.
    elseif (gl%visco%etamodel(nrsubst) == 'VS5') then
        read(ETA(3),*) gl%visco%tmineta(nrsubst)                     ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)                     ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)                     ![kPa]
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)                   ![mol/l]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
        read(ETA(7),*) gl%visco%dilgas_term(nrsubst)                 !number of terms allocated with dilute-gas function
        read(ETA(8),*) gl%visco%pointereta(nrsubst)                  !Points out which collision integral model is used
        read(ETA(9),*) gl%visco%slj(nrsubst)                         !Lennard-Jones coefficient sigma [nm]
        gl%visco%slj_VS5(nrsubst) = gl%visco%slj(nrsubst)*gl%factor_VS5slj  !unit conversion to Angstroem
        read(ETA(10),*) gl%visco%eklj(nrsubst)                        !Lennard-Jones coefficient epsilon/kappa [K] =Tc/1.2593
        gl%visco%tc_VS5(nrsubst) = gl%visco%eklj(nrsubst)*1.2593d0       !Chung uses a different tc value
        read(ETA(11),*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
        read(ETA(12),*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)            !ce (=0.021357*SQRT(MW)), cep (=Power of T)
        read(ETA(13),*) gl%visco%denstermnr(nrsubst)
        if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN
            !not yet used for VS5
        end if
        read(ETA(14),*) gl%visco%accen_VS5(nrsubst), &
            & gl%visco%dipolered_VS5(nrsubst), gl%visco%kappa_VS5(nrsubst)       !w, mur, kappa for Chung, fit
        read(ETA(15),*) gl%visco%addchung_VS5(nrsubst)                !additional parameters for Chung used if residual parameters are not the same as dilute gas
        read(ETA(16),*) gl%visco%pointercrit_eta(nrsubst)             !pointer to critical enhancement auxiliary function (none used)

        !reading parameters for VS9 model specification
        !cf. Vogel, Span, Herrmann, J. Phys. Chem. Ref. Data 44, 043101 (2015)
    elseif (gl%visco%etamodel(nrsubst) == 'VS7') then
        gl%visco%eta_read = .false.
        errorflag = -5246
    elseif (gl%visco%etamodel(nrsubst) == 'VS9') then
        !Info section
        read(ETA(3),*) gl%visco%tmineta(nrsubst)    ![K]
        read(ETA(4),*) gl%visco%tmaxeta(nrsubst)    ![K]
        read(ETA(5),*) gl%visco%pmaxeta(nrsubst)   ![kPa]
        gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)/gl%factor
        read(ETA(6),*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
        gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
        read(ETA(7),*) gl%visco%omega_model(nrsubst)
        read(ETA(8),*) gl%visco%dilgas_term(nrsubst)
        read(ETA(9),*) gl%visco%tredetadg(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredetadg(nrsubst) !reducing parameters for T, rho, eta

        !-------------------------------------------------------------------------------------------------------------
        !read diluted gas terms
        do ieta=1, gl%visco%dilgas_term(nrsubst)
            read(ETA(9+ieta),*,iostat = err_read) gl%visco%coeff_hi(nrsubst,ieta)
        end do

        last_sum = ieta

        read(ETA(9+ieta),*) gl%visco%ndens_initial(nrsubst)

        !-------------------------------------------------------------------------------------------------------------
        !read initial density terms
        if (gl%visco%ndens_initial(nrsubst) .ne. 0) then

            read(ETA(10+ieta),*) gl%visco%tinit_red(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%etainit_red(nrsubst)  !reducing parameters for T, rho, eta
            gl%visco%rhoredeta(nrsubst)=gl%visco%rhoredeta(nrsubst)*gl%Factor
            do ieta = 1, gl%visco%ndens_initial(nrsubst)
                read(ETA(10+last_sum+ieta),*,iostat = err_read) gl%visco%cieta(nrsubst,ieta), gl%visco%tieta(nrsubst,ieta)
            end do
            last_sum = last_sum + ieta
        end if

        read(ETA(10+last_sum),*) gl%visco%term_num1_eta(nrsubst)

        read(ETA(11+last_sum),*) gl%visco%tredeta(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst)  !reducing parameters for T, rho, eta
        gl%visco%rhoredeta(nrsubst) = gl%visco%rhoredeta(nrsubst) * gl%factor

        !-------------------------------------------------------------------------------------------------------------
        !read polynomial terms
        if (gl%visco%term_num1_eta(nrsubst) .ne. 0) then
            do ieta = 1,gl%visco%term_num1_eta(nrsubst)
                read(ETA(11+last_sum+ieta),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta)
            end do
            last_sum = last_sum + ieta
        end if

        !-------------------------------------------------------------------------------------------------------------
        if (gl%substcasnr(nrsubst)=="106-97-8") then    !--------------------------BUTANE-----------------------------

            !read complex exponential terms
            read(ETA(11+last_sum),*) gl%visco%term_com_expo_eta(nrsubst)

            do ieta =gl%visco%term_num1_eta(nrsubst)+1,gl%visco%term_num1_eta(nrsubst)+gl%visco%term_com_expo_eta(nrsubst)
                read(ETA(12+last_sum),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta)
                last_sum = last_sum + 1
            end do
            gl%visco%term_num1_eta(nrsubst) = gl%visco%term_num1_eta(nrsubst)+gl%visco%term_com_expo_eta(nrsubst)  !number of simple terms (polynomial + exponential terms)

        elseif((gl%substcasnr(nrsubst)=="74-84-0") .or. (gl%substcasnr(nrsubst)=="74-98-6")) then  !--------------------------ETHANE and PROPANE------------------------------
            !read simple exponential terms
            read(ETA(11+last_sum),*,iostat = err_read) gl%visco%term_expo_eta(nrsubst)
            if (err_read .eq. 0) then
                read(ETA(12+last_sum),*) gl%visco%tredeta(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst)  !reducing parameters for T, rho, eta
                gl%visco%rhoredeta(nrsubst) = gl%visco%rhoredeta(nrsubst) * gl%factor
                do
                    read(ETA(13+last_sum),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta)
                    if (err_read .ne. 0) exit
                    ieta = ieta + 1
                    last_sum = last_sum  + 1
                end do
            end if
            gl%visco%term_num1_eta(nrsubst) = gl%visco%term_num1_eta(nrsubst) + gl%visco%term_expo_eta(nrsubst)  !number of simple terms (polynomial + exponential terms)

            !read complex terms
            read(ETA(13+last_sum),*,iostat = err_read) gl%visco%term_num2_eta(nrsubst)
            if (err_read .eq. 0) then
                do ieta = gl%visco%term_num1_eta(nrsubst)+1,gl%visco%term_num1_eta(nrsubst)+gl%visco%term_num2_eta(nrsubst)
                    read(ETA(14+last_sum),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta), gl%visco%eta_VS9(nrsubst,ieta), gl%visco%beta_VS9(nrsubst,ieta), gl%visco%par1(nrsubst,ieta), gl%visco%par2(nrsubst,ieta)
                    last_sum = last_sum + 1
                end do
            end if
        else        !---------------ALL OTHER VS9 FLUIDS------------------------------------------------------------------
            !-------------------------------------------------------------------------------------------------------------
            !read simple exponential terms
            read(ETA(11+last_sum),*,iostat = err_read) gl%visco%term_expo_eta(nrsubst)
            if (err_read .eq. 0) then
                read(ETA(12+last_sum),*) gl%visco%tredeta(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst)  !reducing parameters for T, rho, eta
                gl%visco%rhoredeta(nrsubst) = gl%visco%rhoredeta(nrsubst) * gl%factor
                do
                    read(ETA(13+last_sum),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta)
                    if (err_read .ne. 0) exit
                    ieta = ieta + 1
                    last_sum = last_sum  + 1
                end do
            end if
            gl%visco%term_num1_eta(nrsubst) = gl%visco%term_num1_eta(nrsubst) + gl%visco%term_expo_eta(nrsubst)  !number of simple terms (polynomial + exponential terms)

            !-------------------------------------------------------------------------------------------------------------
            !read complex terms
            read(ETA(12+last_sum),*,iostat = err_read) gl%visco%term_num2_eta(nrsubst)
            if (err_read .eq. 0) then
                do ieta = gl%visco%term_num1_eta(nrsubst)+1,gl%visco%term_num1_eta(nrsubst)+gl%visco%term_num2_eta(nrsubst)
                    read(ETA(13+last_sum),*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS9(nrsubst,ieta), gl%visco%li_VS9(nrsubst,ieta), gl%visco%eta_VS9(nrsubst,ieta), gl%visco%beta_VS9(nrsubst,ieta), gl%visco%par1(nrsubst,ieta), gl%visco%par2(nrsubst,ieta)
                    last_sum = last_sum + 1
                end do
            end if
        end if
    else
        gl%visco%eta_read = .false.
        errorflag = -5243

    END IF




    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of thermal conductivity model
    ! Input: gl<type_gl>: global variable structure
    !        TCX_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_tcx(gl,TCX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::TCX
    integer:: nrsubst,errorflag,last_term,i
    character(20):: eos_type

    errorflag = 0
    if(.not.(allocated(TCX))) then
        errorflag = -1
        return !return when given model is not present for current substance
    endif
    last_term = 0

    if(.not.allocated(gl%tcx)) then
        allocate(gl%tcx)
        if (.not. allocated(gl%tcx%accen_TC5)) then
            allocate(gl%tcx%accen_TC5(gl%ncomp),gl%tcx%tc_TC5(gl%ncomp),gl%tcx%dipolered_TC5(gl%ncomp),gl%tcx%kappa_TC5(gl%ncomp),gl%tcx%addchung_TC5(gl%ncomp))
        endif

        call initialize_thermal_conductivity(gl)
    end if

    read(TCX(2),*) gl%tcx%tcmodel(nrsubst)
    gl%tcx%tcx_read = .true.


    if (gl%tcx%tcmodel(nrsubst) == 'TC1') then
        read(TCX(3),*) gl%tcx%tmin_tc(nrsubst)                                  ! K
        read(TCX(4),*) gl%tcx%tmax_tc(nrsubst)                                  ! K
        read(TCX(5),*) gl%tcx%pmax_tc(nrsubst)                                  ! kPa
        gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)/gl%factor             !MPa
        read(TCX(6),*) gl%tcx%rhomax_tc(nrsubst)                                ! mol/L
        gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor          ! mol/m?
        read(TCX(7),*) gl%tcx%num_dil_tc(nrsubst), gl%tcx%den_dil_tc(nrsubst)
        read(TCX(8),*) gl%tcx%tred_dil(nrsubst), gl%tcx%tcx_dil(nrsubst)
        do i=1, gl%tcx%num_dil_tc(nrsubst)
            read(TCX(9+last_term),*) gl%tcx%num_tc_coeff_dil(nrsubst,i), gl%tcx%num_tc_powerT_dil(nrsubst,i)
            last_term = last_term + 1
        end do
        do i=1, gl%tcx%den_dil_tc(nrsubst)
            read(TCX(9+last_term),*) gl%tcx%den_tc_coeff_dil(nrsubst,i), gl%tcx%den_tc_exp_dil(nrsubst,i)
            last_term = last_term + 1
        enddo
        read(TCX(9+last_term),*) gl%tcx%num_bgrd_tc(nrsubst), gl%tcx%den_bgrd_tc(nrsubst)
        read(TCX(10+last_term),*) gl%tcx%tred_bgrd(nrsubst), gl%tcx%rhored_bgrd(nrsubst), gl%tcx%tcx_bgrd(nrsubst)
        do i=1, gl%tcx%num_bgrd_tc(nrsubst)
            read(TCX(11+last_term),*) gl%tcx%num_tc_coeff_bgrd(nrsubst,i), gl%tcx%num_tc_powerT_bgrd(nrsubst,i), gl%tcx%num_tc_rho_bgrd(nrsubst,i), gl%tcx%rho_exp_bgrd(nrsubst,i)
            last_term = last_term + 1
        end do
        read(TCX(11+last_term),*) gl%tcx%tkmodel(nrsubst)

    elseif (gl%tcx%tcmodel(nrsubst) == 'TC0') then
        read(TCX(3),*) gl%tcx%tmin_tc(nrsubst)                                  ! K
        read(TCX(4),*) gl%tcx%tmax_tc(nrsubst)                                  ! K
        read(TCX(5),*) gl%tcx%pmax_tc(nrsubst)                                  ! kPa
        gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)/gl%factor             ! MPa
        read(TCX(6),*) gl%tcx%rhomax_tc(nrsubst)                                ! mol/L
        gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor
        read(TCX(7),*) gl%tcx%pointer_hardtc(nrsubst)
        read(TCX(8),*) gl%tcx%tcxlam0_num(nrsubst), gl%tcx%tcxlam0_den(nrsubst), gl%tcx%dellam(nrsubst), gl%tcx%dellamcr(nrsubst)

        If(gl%tcx%pointer_hardtc(nrsubst) == 'R23') then
            read(TCX(9),*) gl%tcx%tred_tc(nrsubst), gl%tcx%rhored_tc(nrsubst), gl%tcx%etared_tc(nrsubst)
            read(TCX(10),*) gl%tcx%rholl_r23
            gl%tcx%rholl_r23=gl%tcx%rholl_r23*gl%factor
            read(TCX(11),*) gl%tcx%b1_r23
            read(TCX(12),*) gl%tcx%b2_r23
            read(TCX(13),*) gl%tcx%c1l_r23
            read(TCX(14),*) gl%tcx%c2l_r23
            read(TCX(15),*) gl%tcx%delgl_r23
            read(TCX(16),*) gl%tcx%lmax_r23
            return
        end if
        last_term = last_term  + 1
        if (gl%tcx%pointer_hardtc(nrsubst) /= 'ETY') then
            read(TCX(8+last_term),*) gl%tcx%tred_tc(nrsubst), gl%tcx%rhored_tc(nrsubst), gl%tcx%etared_tc(nrsubst)      ! use tred_bgrd/rhored_tc instead of tred_tc/rhored_tc?!
            gl%tcx%rhored_tc(nrsubst) = gl%tcx%rhored_tc(nrsubst) * gl%factor   !mol/l  =>  mol/m³
            do i=1, gl%tcx%tcxlam0_num(nrsubst)
                read(TCX(9+last_term),*) gl%tcx%lam0_coeff(nrsubst,i), gl%tcx%lam0_exp(nrsubst,i)
                last_term = last_term  + 1
            enddo
            ! here lamo0_den_coeff...
            do i=1, gl%tcx%dellam(nrsubst)
                read(TCX(9+last_term),*) gl%tcx%dellam_coeff(nrsubst,i), gl%tcx%dellam_exp(nrsubst,i), gl%tcx%dellam_exp2(nrsubst,i), gl%tcx%dellam_ln(nrsubst,i)
                last_term = last_term  + 1
            enddo
            do i=1, gl%tcx%dellamcr(nrsubst)
                read(TCX(9+last_term),*) gl%tcx%dellamcr_coeff(nrsubst,i), gl%tcx%dellamcr_exp(nrsubst,i)
                last_term = last_term  + 1
            enddo
        endif

        read(TCX(9+last_term),*) gl%tcx%tkmodel(nrsubst)

    elseif (gl%tcx%tcmodel(nrsubst) == 'TC3') then                 ! xenon
        read(TCX(3),*) gl%tcx%tmin_tc(nrsubst)                          ! K
        read(TCX(4),*) gl%tcx%tmax_tc(nrsubst)                          ! K
        read(TCX(5),*) gl%tcx%pmax_tc(nrsubst)                          ! kPa
        gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)/gl%factor      ! MPa
        read(TCX(6),*) gl%tcx%rhomax_tc(nrsubst)                        ! mol/L
        gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor  ! mol/m?
        read(TCX(7),*) gl%tcx%lj_sigma(nrsubst)                         ! nm
        read(TCX(8),*) gl%tcx%lj_epskap(nrsubst)                        ! K
        read(TCX(9),*) gl%tcx%const_eq20(nrsubst)
        read(TCX(10),*) gl%tcx%Texp_eq20(nrsubst)
        do i=1, 9   ! manually
            read(TCX(11+last_term),*) gl%tcx%eta0(nrsubst,i)
            last_term = last_term + 1
        enddo
        do i=1, 4   ! manually
            read(TCX(11+last_term),*) gl%tcx%Fvi(nrsubst,i)
            last_term = last_term + 1
        enddo
        do i=1, 8   ! manually
            read(TCX(11+last_term),*) gl%tcx%Evi(nrsubst,i)
            last_term = last_term + 1
        enddo
        read(TCX(11+last_term),*) gl%tcx%F_tc3(nrsubst)
        read(TCX(12+last_term),*) gl%tcx%rm(nrsubst)
        read(TCX(13+last_term),*) gl%tcx%tkmodel(nrsubst)
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC5') then
        read(TCX(3),*) gl%tcx%tmin_tc(nrsubst)                          ! K
        read(TCX(4),*) gl%tcx%tmax_tc(nrsubst)                          ! K
        read(TCX(5),*) gl%tcx%pmax_tc(nrsubst)                          ! kPa
        gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)/gl%factor      ! MPa
        read(TCX(6),*) gl%tcx%rhomax_tc(nrsubst)                        ! mol/L
        gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor  ! mol/m?
        read(TCX(7),*) gl%tcx%lj_sigma(nrsubst)                         ! nm
        read(TCX(8),*) gl%tcx%lj_epskap(nrsubst)                        ! K
        gl%tcx%tc_TC5(nrsubst) = gl%tcx%lj_epskap(nrsubst)*1.2593d0       !Chung uses a different tc value
        read(TCX(9),*) gl%tcx%accen_TC5(nrsubst), gl%tcx%dipolered_TC5(nrsubst), gl%tcx%kappa_TC5(nrsubst)       !w, mur, kappa for Chung, fit
        read(TCX(10),*) gl%tcx%addchung_TC5(nrsubst)                !additional parameters for Chung used if residual parameters are not the same as dilute gas
        read(TCX(11),*) gl%tcx%tkmodel(nrsubst)

    else
        gl%tcx%tcx_read = .false.
    endif

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of thermal conductivity model
    ! Input: gl<type_gl>: global variable structure
    !        TCX_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_tcx_aux(gl,TCX_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::TCX_AUX
    character(20):: eos_type
    integer:: nrsubst,errorflag,last_sum,i
    character(255) :: dummy
    errorflag = 0

    if(.not.(allocated(TCX_AUX))) return !return when given model is not present for current substance
    last_sum = 0


    read(TCX_AUX(2),*) eos_type

    if ((eos_type(1:3) == 'TK1') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK1')) then

        read(TCX_AUX(3),*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
        read(TCX_AUX(4),*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
        read(TCX_AUX(5),*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
        read(TCX_AUX(6),*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
        read(TCX_AUX(7),*) gl%tcx%npnum_tk1(nrsubst),gl%tcx%npdenom_tk1(nrsubst),gl%tcx%nexp_tk1(nrsubst),gl%tcx%nspare_tk1(nrsubst)        !number of terms
        read(TCX_AUX(8),*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
        do i=1,gl%tcx%npnum_tk1(nrsubst)
            read(TCX_AUX(8+last_sum),*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
            last_sum = last_sum +1
        end do

        do i=gl%tcx%npnum_tk1(nrsubst)+1,gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)
            read(TCX_AUX(8+last_sum),*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
            last_sum = last_sum +1
        end do

        read(TCX_AUX(8+last_sum),*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)

        last_sum = last_sum +1
        gl%tcx%rhored_tk(nrsubst)=gl%tcx%rhored_tk(nrsubst)*gl%factor

        do i=gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)+1,gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)+gl%tcx%nexp_tk1(nrsubst)
            read(TCX_AUX(8+last_sum),*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
            last_sum = last_sum +1
        end do
        gl%tcx%tcx_crit_read=.true.
        return

    elseif ((eos_type(1:3) == 'TK3') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK3')) then

        read(TCX_AUX(3),*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
        read(TCX_AUX(4),*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
        read(TCX_AUX(5),*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
        read(TCX_AUX(6),*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
        read(TCX_AUX(7),*) gl%tcx%term_tk(nrsubst)                 !number of terms
        read(TCX_AUX(8),*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
        read(TCX_AUX(9),*) gl%tcx%gnu_tk(nrsubst)
        read(TCX_AUX(10),*) gl%tcx%gamma_tk(nrsubst)
        read(TCX_AUX(11),*) gl%tcx%R0_tk(nrsubst)
        read(TCX_AUX(12),*) gl%tcx%z_visonly_tk(nrsubst)
        read(TCX_AUX(13),*) gl%tcx%c_tk(nrsubst)
        read(TCX_AUX(14),*) gl%tcx%xi0_tk(nrsubst)
        read(TCX_AUX(15),*) gl%tcx%gam0_tk(nrsubst)
        read(TCX_AUX(16),*) gl%tcx%qd_inverse_tk(nrsubst)
        read(TCX_AUX(17),*) gl%tcx%tref_tk(nrsubst)
        gl%tcx%tcx_crit_read=.true.
        return

    elseif ((eos_type(1:3) == 'TK6') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK6')) then

        read(TCX_AUX(3),*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
        read(TCX_AUX(4),*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
        read(TCX_AUX(5),*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
        read(TCX_AUX(6),*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
        read(TCX_AUX(7),*) gl%tcx%term_tk(nrsubst)                 !number of terms
        read(TCX_AUX(8),*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
        read(TCX_AUX(9),*) gl%tcx%gnu_tk(nrsubst)
        read(TCX_AUX(10),*) gl%tcx%gamma_tk(nrsubst)
        read(TCX_AUX(11),*) gl%tcx%R0_tk(nrsubst)
        read(TCX_AUX(12),*) gl%tcx%z_visonly_tk(nrsubst)
        read(TCX_AUX(13),*) gl%tcx%c_tk(nrsubst)
        read(TCX_AUX(14),*) gl%tcx%xi0_tk(nrsubst)
        read(TCX_AUX(15),*) gl%tcx%gam0_tk(nrsubst)
        read(TCX_AUX(16),*) gl%tcx%qd_inverse_tk(nrsubst)
        read(TCX_AUX(17),*) gl%tcx%tref_tk(nrsubst)
        gl%tcx%tcx_crit_read=.true.
        return
    end if

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of sublimation line model
    ! Input: gl<type_gl>: global variable structure
    !        SBL_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_sbl(gl,SBL_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::SBL_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0

    if(.not.(allocated(SBL_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%psubcoeff)) then
        allocate(gl%psubcoeff(20,gl%ncomp))
        allocate(gl%psubexp, mold=gl%psubcoeff)
    else
        if (size(gl%psubcoeff,2) /= gl%ncomp) then
            deallocate(gl%psubcoeff,gl%psubexp)
            allocate(gl%psubcoeff(20,gl%ncomp))
            allocate(gl%psubexp, mold=gl%psubcoeff)
        endif
    end if

    read(SBL_AUX(2),*) eos_type
    gl%psubtype(nrsubst) = ichar(eos_type(3:3)) - 48 !get the number at 3. place
    !read reducing parameters for sublimation pressure:
    read(SBL_AUX(7),*) gl%tpsubred(nrsubst), gl%psubred(nrsubst)
    gl%psubred(nrsubst) = gl%psubred(nrsubst)/gl%Factor ![MPa]
    !read coeffs ande exps:
    read(SBL_AUX(8),*) gl%npsubcoeff1(nrsubst), gl%npsubcoeff2(nrsubst), gl%npsubcoeff3(nrsubst)
    do j=1,(gl%npsubcoeff1(nrsubst) + gl%npsubcoeff2(nrsubst) + gl%npsubcoeff3(nrsubst))
        read(SBL_AUX(8+j),*) gl%psubcoeff(j,nrsubst),gl%psubexp(j,nrsubst)
    enddo


    end subroutine
    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of vapour pressure model
    ! Input: gl<type_gl>: global variable structure
    !        PV_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_pv(gl,PV_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::PV_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0


    if(.not.(allocated(PV_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%vpcoeff)) then
        allocate(gl%vpcoeff(20,gl%ncomp))
        allocate(gl%vpexp, mold=gl%vpcoeff)
    else
        if (size(gl%vpcoeff,2) /= gl%ncomp) then
            deallocate(gl%vpcoeff,gl%vpexp)
            allocate(gl%vpcoeff(20,gl%ncomp))
            allocate(gl%vpexp, mold=gl%vpcoeff)
        endif
    end if

    read(PV_AUX(2),*) eos_type
    gl%vptype(nrsubst) = ichar(eos_type(3:3)) - 48 !get the number at 3. place

    !read reducing parameters for vapour pressure:
    read(PV_AUX(7),*) gl%tvpred(nrsubst), gl%vpred(nrsubst)
    gl%vpred(nrsubst) = gl%vpred(nrsubst)/gl%factor ![MPa]
    !read coeffs ande exps:
    read(PV_AUX(8),*) gl%nvpcoeff(nrsubst)

    do j=1,gl%nvpcoeff(nrsubst)
        read(PV_AUX(8+j),*) gl%vpcoeff(j,nrsubst),gl%vpexp(j,nrsubst)
    enddo
    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of liquid saturation density model
    ! Input: gl<type_gl>: global variable structure
    !        DL_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_dl(gl,DL_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::DL_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0

    if(.not.(allocated(DL_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%dlcoeff)) then
        allocate(gl%dlcoeff(20,gl%ncomp))
        allocate(gl%dlexp,mold=gl%dlcoeff)
    else
        if (size(gl%dlcoeff,2) /= gl%ncomp) then
            deallocate(gl%dlcoeff,gl%dlexp)
            allocate(gl%dlcoeff(20,gl%ncomp))
            allocate(gl%dlexp,mold=gl%dlcoeff)
        endif
    end if

    read(DL_AUX(2),*) eos_type
    gl%dltype(nrsubst) = ichar(eos_type(3:3)) - 48 !get the number at 3. place

    !read reducing parameters for saturated liquid density:
    read(DL_AUX(7),*) gl%tdlred(nrsubst), gl%dlred(nrsubst)
    gl%dlred(nrsubst) = gl%dlred(nrsubst)*gl%factor ![mol/m³]
    !read coeffs ande exps:
    read(DL_AUX(8),*) gl%ndlcoeff(nrsubst)
    do j=1,gl%ndlcoeff(nrsubst)
        read(DL_AUX(8+j),*) gl%dlcoeff(j,nrsubst),gl%dlexp(j,nrsubst)
    enddo
    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of vapour saturation density model
    ! Input: gl<type_gl>: global variable structure
    !        DV_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_dv(gl,DV_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::DV_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0

    if(.not.(allocated(DV_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%dvcoeff)) then
        allocate(gl%dvcoeff(20,gl%ncomp))
        allocate(gl%dvexp,mold=gl%dvcoeff)
    else
        if (size(gl%dvcoeff,2) /= gl%ncomp) then
            deallocate(gl%dvcoeff,gl%dvexp)
            allocate(gl%dvcoeff(20,gl%ncomp))
            allocate(gl%dvexp,mold=gl%dvcoeff)
        endif
    end if

    read(DV_AUX(2),*) eos_type
    gl%dvtype(nrsubst) = ichar(eos_type(3:3)) - 48 !get the number at 3. place

    !read reducing parameters for saturated vapour density:
    read(DV_AUX(7),*) gl%tdvred(nrsubst), gl%dvred(nrsubst)
    gl%dvred(nrsubst) = gl%dvred(nrsubst)*gl%factor ![mol/m³]
    !read coeffs ande exps:
    read(DV_AUX(8),*) gl%ndvcoeff(nrsubst)
    do j=1,gl%ndvcoeff(nrsubst)
        read(DV_AUX(8+j),*) gl%dvcoeff(j,nrsubst),gl%dvexp(j,nrsubst)
    enddo

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of surface tension model
    ! Input: gl<type_gl>: global variable structure
    !        STN_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_stn(gl,STN_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::STN_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0
    gl%stn_read = .true.

    if(.not.(allocated(STN_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%stn)) then
        allocate(surface_tension_t::gl%stn(gl%ncomp))
    else
        if (size(gl%stn,1) /= gl%ncomp) then
            deallocate(gl%stn)
            allocate(surface_tension_t::gl%stn(gl%ncomp))
        endif
    endif
    read(STN_AUX(2),*) gl%stn(nrsubst)%stmodel

    !reading Parameters for ST1
    if (gl%stn(nrsubst)%stmodel == 'ST1') then
        read(STN_AUX(3),*) gl%stn(nrsubst)%low_temp_st
        read(STN_AUX(4),*) gl%stn(nrsubst)%upp_temp_st
        read(STN_AUX(5),*) gl%stn(nrsubst)%upp_press_st
        read(STN_AUX(6),*) gl%stn(nrsubst)%max_dens_st
        read(STN_AUX(7),*) gl%stn(nrsubst)%term_num_st
        read(STN_AUX(8),*) gl%stn(nrsubst)%Temp_crit_st

        do j=1,(gl%stn(nrsubst)%term_num_st)
            read(STN_AUX(8+j),*) gl%stn(nrsubst)%sigma0_st(j),gl%stn(nrsubst)%n_exp_st(j)
        end do
    else
        gl%stn_read = .false.
    endif

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the parameter of dielectric model
    ! Input: gl<type_gl>: global variable structure
    !        DE_AUX<character array>: text block in fluid file
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helm_de(gl,DE_AUX,nrsubst,errorflag)
    ! Variable declaration
    implicit none
    type(type_gl) :: gl
    character(256),dimension(:),allocatable::DE_AUX
    integer:: errorflag,nrsubst,j
    character(20):: eos_type
    errorflag = 0

    if(.not.(allocated(DE_AUX))) return !return when given model is not present for current substance
    if(.not.allocated(gl%de)) allocate(gl%de)

    read(DE_AUX(2),*) gl%de%demodel(nrsubst)

    !reading Parameters for DE3 and DE4 together, spereated from DE2
    if ((gl%de%demodel(nrsubst) == 'DE3') .OR. (gl%de%demodel(nrsubst) == 'DE4'))  then
        read(DE_AUX(3),*) gl%de%tmin_de(nrsubst)
        read(DE_AUX(4),*) gl%de%tmax_de(nrsubst)
        read(DE_AUX(5),*) gl%de%pmax_de(nrsubst)
        read(DE_AUX(6),*) gl%de%rhomax_de(nrsubst)
        read(DE_AUX(7),*) gl%de%ref_temp_de(nrsubst)
        read(DE_AUX(8),*) gl%de%term_num1_de(nrsubst), gl%de%term_num2_de(nrsubst), gl%de%term_num3_de(nrsubst), &
            gl%de%term_num4_de(nrsubst), gl%de%term_num5_de(nrsubst), gl%de%term_num6_de(nrsubst)

        !reading coefficients
        if ((gl%de%term_num1_de(nrsubst) /= 0).or.(gl%de%term_num2_de(nrsubst) /= 0).or.(gl%de%term_num3_de(nrsubst) /= 0)) then
            do j=1,(gl%de%term_num1_de(nrsubst)+gl%de%term_num2_de(nrsubst)+gl%de%term_num3_de(nrsubst))
                read(DE_AUX(8+j),*) gl%de%coeffde(nrsubst,j),gl%de%texpde(nrsubst,j),gl%de%dexpde(nrsubst,j)
            end do
        end if
        !reading Parameters for DE2
    elseif (gl%de%demodel(nrsubst) == 'DE2') then
        read(DE_AUX(3),*) gl%de%tmin_de(nrsubst)
        read(DE_AUX(4),*) gl%de%tmax_de(nrsubst)
        read(DE_AUX(5),*) gl%de%pmax_de(nrsubst)
        read(DE_AUX(6),*) gl%de%rhomax_de(nrsubst)
        read(DE_AUX(7),*) gl%de%ref_temp_de(nrsubst) , gl%de%ref_dens_de(nrsubst), gl%de%ref_press_de(nrsubst)
        read(DE_AUX(8),*) gl%de%term_num1_de(nrsubst), gl%de%term_num2_de(nrsubst), gl%de%term_num3_de(nrsubst), &
            gl%de%term_num4_de(nrsubst), gl%de%term_num5_de(nrsubst), gl%de%term_num6_de(nrsubst)

        if ((gl%de%term_num1_de(nrsubst) /= 0).or.(gl%de%term_num2_de(nrsubst) /= 0).or.(gl%de%term_num3_de(nrsubst) /= 0)) then
            do j=1,(gl%de%term_num1_de(nrsubst)+gl%de%term_num2_de(nrsubst)+gl%de%term_num3_de(nrsubst))
                read(DE_AUX(8+j),*) gl%de%coeffde(nrsubst,j),gl%de%texpde(nrsubst,j),gl%de%dexpde(nrsubst,j),gl%de%pexpde(nrsubst,j)
            end do
        end if

    else

        gl%de%de_read =.false.

    endif

    end subroutine

    !##################################################################################################
    !****************************************************************************************************************************************************************
    !##################################################################################################
    ! Read the head of a fluid file
    ! Input: gl<type_gl>: global variable structure
    !        path<character>: absolut path to the model
    !        folder<character>: name of the folder
    !        nrsubst<integer>: Number of actual component
    !        errorflag<integer>: Error information ( 0 =  NO Error )
    !--------------------------------------------------------------------------------------------------
    subroutine read_helmholtz_mix(gl, substance1, substance2 , path, folder, is_hc, errorflag)
    implicit none
    character(*) :: path, folder
    character(255) :: path1, path2, casnr1, casnr2, dummy
    type(type_gl) :: gl
    type(file_t) :: mix_file
    type(return_pos_t):: sub_part
    integer:: substance1,substance2,errorflag,j
    integer :: firstsubst, secondsubst
    integer :: dfnreg, dfnlreg, dfnspec, dfnlspec, dfngauss, dfnlgauss,dfpli,dftli,dummynr1
    logical :: is_hc ! info if file is hc or normal path file

    character(256),dimension(:),allocatable::head,MIX
    character(10), dimension(2):: mix_name
    !character(:),dimension(:),allocatable:: dummy
    character(len=3) :: modeltype
    double precision :: betat,gammat,betav,gammav,fij_dat
    double precision :: dfni,dfti,dfdi,dfli,dfetai,dfepsi,dfbetai,dfgammai,dfgetai,dfgepsi,dfgbetai,dfggami,dfgdummy1,dfgdummy2,dfgdummy3, dfpi
    integer, dimension(2) :: bin_combo

    errorflag = 0
    firstsubst = 0
    secondsubst = 0
    dfnreg = 0
    dfnlreg = 0
    dfnspec = 0
    dfnlspec = 0
    dummynr1 = 0
    betat=0.d0
    gammat=0.d0
    betav=0.d0
    gammav=0.d0
    fij_dat=0.d0
    dfni=0.d0
    dfti=0.d0
    dfdi=0.d0
    dfli=0.D0
    dfetai=0.d0
    dfepsi=0.d0
    dfbetai=0.d0
    dfgammai=0.d0
    dfpi = 0.d0


    path1 = trim(path) // trim(gl%components(substance1)) //'-'// trim(gl%components(substance2)) //'.mix'
    path2 = trim(path) // trim(gl%components(substance2)) //'-'// trim(gl%components(substance1)) //'.mix'

    bin_combo(1) = substance1
    bin_combo(2) = substance2

    ! Init the file variable
    mix_file = file_t(is_hc)

    !get the content of the fluid file
    call get_content(mix_file,path1,gl%components(bin_combo),folder)

    if(mix_file%errorflag == 0 ) then
        firstsubst = substance1
        secondsubst = substance2
    else
        deallocate(mix_file%content)
        bin_combo(1) = substance2
        bin_combo(2) = substance1
        call get_content(mix_file,path2,gl%components(bin_combo),folder)
        firstsubst = substance2
        secondsubst = substance1
    end if

    if(mix_file%errorflag .ne. 0 ) then
        errorflag = mix_file%errorflag
        return
    end if

    if(.not.allocated(gl%litref)) allocate(gl%litref)
    gl%litref%lit_ref_mix = get_ref_from_sub_content(mix_file,1,len(gl%litref%lit_ref_mix))
    !sub_part%pos_id(1) = 1
    !gl%litref%lit_ref_mix = mix_file%content(2)
    !gl%litref%lit_ref_mix = trim(gl%litref%lit_ref_mix) // ' ' // get_ref_from_sub_content(mix_file,sub_part%pos_id(1),len(gl%litref%lit_ref_mix))

    read(mix_file%content(3)(index(mix_file%content(3),' '):),*) modeltype,betat,gammat,betav,gammav,fij_dat,dummy !(3)

    If ((modeltype(1:2) == 'KW').or.(modeltype(1:1) == "B")) then
        gl%rfbetat(firstsubst,secondsubst) = betat
        gl%rfbetarho(firstsubst,secondsubst) = betav !gammat
        gl%rfgammat(firstsubst,secondsubst) = gammat !betav
        gl%rfgammarho(firstsubst,secondsubst) = gammav
        gl%rfbetat(secondsubst,firstsubst) = 1/betat
        gl%rfbetarho(secondsubst,firstsubst) = 1/betav !gammat
        gl%rfgammat(secondsubst,firstsubst) = gammat !1/betav
        gl%rfgammarho(secondsubst,firstsubst) = gammav
    else
        !if not modeltype not "KW" then set default values to 1:
        gl%rfbetat(firstsubst,secondsubst) = 1.d0
        gl%rfbetarho(firstsubst,secondsubst) = 1.d0
        gl%rfgammat(firstsubst,secondsubst) = 1.d0
        gl%rfgammarho(firstsubst,secondsubst) = 1.d0
        gl%rfbetat(secondsubst,firstsubst) = 1.d0
        gl%rfbetarho(secondsubst,firstsubst) = 1.d0
        gl%rfgammat(secondsubst,firstsubst) = 1.d0
        gl%rfgammarho(secondsubst,firstsubst) = 1.d0
    endif

    if (fij_dat == 0.0d0) return

    MIX_name = modeltype
    sub_part = get_sub_string_lines(mix_file,MIX_name,errorflag)
    if (errorflag == 0) then
        MIX = get_content_until_empty_line(mix_file,sub_part%pos_id(1))
        !delete all comment lines
        call delete_lines(mix_file,MIX)
        MIX_name = 'BetaT'
        sub_part = get_sub_string_lines(mix_file,MIX_name,errorflag)
        if ((modeltype(1:2) /= 'KW').and.(modeltype(1:1) /= "B")) then !use the default values:

            read(mix_file%content(sub_part%POS_ID(1)+1),*) betat,gammat,betav,gammav,gl%fij,dummy
            gl%rfbetat(firstsubst,secondsubst) = betat
            gl%rfbetarho(firstsubst,secondsubst) = gammat
            gl%rfgammat(firstsubst,secondsubst) = betav
            gl%rfgammarho(firstsubst,secondsubst) = gammav
            gl%rfbetat(secondsubst,firstsubst) = 1/betat
            gl%rfbetarho(secondsubst,firstsubst) = gammat
            gl%rfgammat(secondsubst,firstsubst) = 1/betav
            gl%rfgammarho(secondsubst,firstsubst) = gammav
        else
            read(mix_file%content(sub_part%POS_ID(1)+2),*) dfnreg, dfnlreg, dummynr1, dfnspec, dfnlspec, dfngauss, dfnlgauss
            gl%dfpol(firstsubst,secondsubst) = dfnreg
            gl%dfpol(secondsubst,firstsubst) = dfnreg
            gl%dfexp(firstsubst,secondsubst) = dfnspec
            gl%dfexp(secondsubst,firstsubst) = dfnspec
            gl%dfgau(firstsubst,secondsubst) = dfngauss
            gl%dfgau(secondsubst,firstsubst) = dfngauss
            gl%fij(firstsubst,secondsubst) = fij_dat
            gl%fij(secondsubst,firstsubst) = fij_dat
        endif

        do j=1,dfnreg !for normal terms
            if (dfnlreg == 4) then  !checks, if a value for dfli exists for this mixture - J.Gernert, Sept. 2010
                read(mix_file%content(j+sub_part%POS_ID(1)+2),*) dfni,dfti,dfdi, dfli
                if (dabs(dfli) > 1d-12) dfpi = 1d0
            elseif (dfnlreg == 5) then
                read(mix_file%content(j+sub_part%POS_ID(1)+2),*) dfni,dfti,dfdi, dfli, dfpi
            else
                read(mix_file%content(j+sub_part%POS_ID(1)+2),*) dfni,dfti,dfdi
            end if
            !---------------------------------------
            !Andreas Aug 2010
            gl%dfn(j,firstsubst,secondsubst) = dfni
            gl%dfn(j,secondsubst,firstsubst) = dfni
            gl%dfd(j,firstsubst,secondsubst) = dfdi
            gl%dfd(j,secondsubst,firstsubst) = dfdi
            gl%dft(j,firstsubst,secondsubst) = dfti
            gl%dft(j,secondsubst,firstsubst) = dfti
            gl%dfl(j,firstsubst,secondsubst) = dfli
            gl%dfl(j,secondsubst,firstsubst) = dfli
            gl%dfp(j,firstsubst,secondsubst) = dfpi
            gl%dfp(j,secondsubst,firstsubst) = dfpi

            select case(gl%dfd_coeff_structure(firstsubst,secondsubst))

            case (coeff_structure_all_zeros)
                if (dfdi /= 0d0) then
                    if (DBLE(INT(dfdi)) /= dfdi) then
                        gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
                    else
                        gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_int
                    endif
                endif

            case (coeff_structure_int)
                if (DBLE(INT(dfdi)) /= dfdi) then
                    gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
                endif

            end select


            select case(gl%dfl_coeff_structure(firstsubst,secondsubst))

            case (coeff_structure_all_zeros)
                if (dfli /= 0d0) then
                    if (DBLE(INT(dfli)) /= dfli) then
                        gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
                    else
                        gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_int
                    endif
                endif

            case (coeff_structure_int)
                if (DBLE(INT(dfli)) /= dfli) then
                    gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
                endif

            end select

            !---------------------------------------
        enddo

        gl%dfd_coeff_structure(secondsubst,firstsubst) = gl%dfd_coeff_structure(firstsubst,secondsubst)
        gl%dfl_coeff_structure(secondsubst,firstsubst) = gl%dfl_coeff_structure(firstsubst,secondsubst)

        if (dfnspec > 0) then !for special terms
            do j=1,dfnspec
                read(mix_file%content(j+sub_part%POS_ID(1)+2+dfnreg),*) dfni,dfti,dfdi,dfetai,dfepsi,dfbetai,dfgammai

                !---------------------------------------
                !Andreas Aug 2010
                gl%dfn(j+dfnreg,firstsubst,secondsubst) = dfni
                gl%dfn(j+dfnreg,secondsubst,firstsubst) = dfni
                gl%dfd(j+dfnreg,firstsubst,secondsubst) = dfdi
                gl%dfd(j+dfnreg,secondsubst,firstsubst) = dfdi
                gl%dft(j+dfnreg,firstsubst,secondsubst) = dfti
                gl%dft(j+dfnreg,secondsubst,firstsubst) = dfti

                gl%dfeta(j,firstsubst,secondsubst) = dfetai
                gl%dfeta(j,secondsubst,firstsubst) = dfetai
                gl%dfeps(j,firstsubst,secondsubst) = dfepsi
                gl%dfeps(j,secondsubst,firstsubst) = dfepsi
                gl%dfbeta(j,firstsubst,secondsubst) = dfbetai
                gl%dfbeta(j,secondsubst,firstsubst) = dfbetai
                gl%dfgamma(j,firstsubst,secondsubst) = dfgammai
                gl%dfgamma(j,secondsubst,firstsubst) = dfgammai
                !----------------------------------------

            enddo
        endif

        if (dfngauss > 0) then !for gaussian terms
            do j=1,dfngauss
                read(mix_file%content(j+sub_part%POS_ID(1)+2+dfnreg+dfnspec),*) dfni,dfti,dfdi,dfpli,dftli,dfgetai,dfgbetai,dfggami,dfgepsi,dfgdummy1,dfgdummy2,dfgdummy3

                !---------------------------------------
                !Andreas Aug 2010
                gl%dfn(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfni
                gl%dfn(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfni
                gl%dfd(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfdi
                gl%dfd(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfdi
                gl%dft(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfti
                gl%dft(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfti

                gl%dfgeta(j,firstsubst,secondsubst) = dfgetai
                gl%dfgeta(j,secondsubst,firstsubst) = dfgetai
                gl%dfgeps(j,firstsubst,secondsubst) = dfgepsi
                gl%dfgeps(j,secondsubst,firstsubst) = dfgepsi
                gl%dfgbeta(j,firstsubst,secondsubst) = dfgbetai
                gl%dfgbeta(j,secondsubst,firstsubst) = dfgbetai
                gl%dfggam(j,firstsubst,secondsubst) = dfggami
                gl%dfggam(j,secondsubst,firstsubst) = dfggami
                !----------------------------------------

            enddo
        endif
    endif

    end subroutine read_helmholtz_mix


    end module