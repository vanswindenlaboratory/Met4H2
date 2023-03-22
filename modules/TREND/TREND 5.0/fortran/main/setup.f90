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
    !          TREND. Thermodynamic Reference and Engineering Data 5.0.s
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

    ! module for file setup.f90
    submodule (setup_module) impl
    !global use inclusion
    use module_all_types
    use utility_module
    use ancillary_equations_mix_module
    use initialize_module
    use seawater_module
    !use fluids_gerg_module
    !use fluids_module
    !use fluids_cubic_module
    !use fluids_eoscg_module
    !use fluids_lkp_module
    !use fluids_pr_module
    use sub_file_input_module
    use flash_pure_module
    use fniderivs_module
    use hdrt_properties_modules_module
    use reduced_parameters_calc_module
    use mixtures_AGA8_module
    use electrolytes
    use unit_convertion
    contains




    module subroutine setup (gl,input, prop1, prop2, fluids, moles, path, eqtype, mixtype, errorflag)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "setup" :: setup
    
    implicit none

    type(type_gl) :: gl


    !define input variables
    character (12), intent(inout) :: input                        ! TD, DT, TP, PT, PH, HP, PS, SP
    double precision, intent(inout) :: prop1, prop2               ! input values
    !character (30), dimension(30), intent(inout) :: fluidl                      ! fluid list
    !character (255), intent(inout) :: molesl                      ! fluid list
    double precision, dimension (30), intent(inout) :: moles      ! molefraction vector
    double precision, dimension (30) :: x_spec
    character (255), intent(inout) :: path                                     ! path where to find fluid and mix files
    !character (255) :: path_lower
    ! In the first version EOS_indicator was an array of kind integer with 31 elements.
    ! Analog to the fluidlist and the list of moles the EOStype is passed to the interface routine as character(255) now. Andreas March 2013
    ! integer, dimension(31):: EOS_indicator
    !integer, dimension(30), intent(inout) :: EOS_indicator
    integer :: errorflag                        ! information about warnings or errors

    character (30), dimension (30), intent(inout) :: fluids                    ! fluid vector
    integer:: errorfld, i, nf, nm, ne
    integer, dimension(30), intent(inout):: eqtype
    integer, intent(inout) :: mixtype
    integer:: convert, converttype, modelflag, comp_save, n_zero
    logical :: resorted, electrolytes, el_present_Save                                        !flag, whether fluid list needed to be resorted due to hydrate formers


    ! A new parameter was added to the setup routine. EOS_indicator specifies, which equations should be used
    ! Andreas, Oct 2012
    ! EOS_indicator 1 - 30 -- >  Choose which EoS is used for which component :        1: Helmholtz Eos
    !                                                                                 2: SRK
    ! EOS_indicator(31) =>  Which mixing rules should be used?                       1: GERG-mixing rules (ALWAYS POSSIBLE, SINCE CUBICS ARE INTEGRATED TO HELMHOLTZ)
    !                                                                                2: SRK-Mixing rules (ONLY SRK EoS FOR ALL COMPONENTS!!!)
    !crossover model is loaded
    !call init_crs(gl,.true.,path)

    errorfld = 0
    gl%same_components = .false.
    x_spec = 0.d0


    ! Catching wrong inputs
    call setpropinput (gl,input, prop1, prop2, errorflag)
    if (errorflag /= 0) return

    !block before el_merge
    !! set given fluids
    !call setfluid (gl, fluids, nf, errorflag)
    !if (errorflag /= 0) return
    !! set given eqtypes
    !call seteqtype (gl, eqtype, mixtype, ne, errorflag)
    !if (errorflag /= 0) return
    !end block before el_merge

    ! set given eqtypes, order changed for easier handling of electrolytes
    call seteqtype (gl, eqtype, mixtype, ne, errorflag)
    if (errorflag /= 0) return


    ! set given fluids
    call setfluid (gl, fluids, nf, errorflag)
    if (errorflag /= 0) return

    !DEC$ IF DEFINED(MESSER)
    continue
    !DEC$ ELSE
    if(gl%el_present) then
        if(gl%el%n_salts .ge. 1) then
            nf = gl%el%n_fluids
            gl%ncomp = nf
            gl%components = gl%el%fluids_wo_salts
            gl%eq_type = gl%el%eq_type_el
            eqtype = gl%eq_type
        end if
    end if
    if(gl%seawater) then
        ne = ne-1
    elseif (nf /= ne) then !number of eqtypes does not match given number of fluids
        errorflag = -9953
        return
    end if
    !DEC$ END IF

    call setpath (gl,path)


    if ((gl%eq_type(1) == 6) .and. (gl%ncomp > 1))  then !TE vorerst keine Gemische mit SAFT
        !   errorflag = -2224
        !   return
    end if
    !At the moment, always read in the parameters again when mixture model 12 (or 13) (Helm+gE) or 22 (PSRK) is used. (Problems with px diagram routine occured). Andreas Jäger, September 2017
    if ((gl%mix_type == 12) .or. (gl%mix_type == 13) .or. (gl%mix_type == 22)) then
        gl%same_components = .false.
        resorted = .false.
    end if

    call allocate_subtypes(gl,errorflag)

    if (.not. gl%same_components) then
        !allocation of eos parameter-arrays to the necessary size
        ! Resetting module variables
        el_present_Save = gl%el_present
        modelflag = gl%modelflag
        call initialize(gl)
        gl%el_present = el_present_Save
        gl%modelflag = modelflag
        ! module variables that have already been evaluated need to be reset
        !gl%components = fluids
        !gl%Eq_type = eqtype
        !gl%mix_type = mixtype
        !gl%ncomp = nf
    elseif (gl%same_components) then
        !if mix_types 110, 111, 120, 121 are used, the same mixing rules are used as in the case of mix_type = 1.
        !In order to avoid adding the four special mix_types in every single routine, it is changed to mix_type = 1 here.
        if ((gl%mix_type .eq. 110) .or. (gl%mix_type .eq. 111) .or. (gl%mix_type .eq. 120) .or. (gl%mix_type .eq. 121)) gl%mix_type = 1
        call initialize_vle(gl)
        call initialize_flash(gl)
    end if


    ! set given composition
    call setmoles (gl, moles, nm, errorflag)

    if (errorflag /= 0) return
    !if (nf /= nm) then !number of compositions does not match given number of fluids
    !    errorflag = -9953
    !    return
    !end if

    if((gl%seawater)) then
        call setseawater(gl)
        gl%seawatercalled = .true.

        if ((trim(input)=='tp') .or. (trim(input)=='tp+')) then
            gl%sea%seap = prop2
            errorflag = 0
            if(gl%hold_limits) then
                call check_seawater_limits(gl, prop1, prop2, errorflag)
                if(errorflag .ne. 0)then
                    return
                end if
            end if

        else
            errorflag = -12800
        end if
    elseif(gl%el_present) then      !kann man das nicht in einer funkiton machen
        gl%gepress = prop2
        gl%el%press = prop2
        gl%el%x_salt = gl%el%molality(1) / (1.d0/18.015268d-3 + gl%el%molality(1))
        gl%el%n2 = gl%el%molality(1) * 18.015268d-3
        gl%el%mol_in = gl%el%molality(1)
        !call read_el_params(gl)
        call check_el_limits(gl, prop1, prop2, input, errorflag)
        if(errorflag .lt. 0) return
    end if

    if((gl%check_solid) .and. (gl%seawater) .and. (nf .ge. 2)) then
        errorflag = -12900
        return
    end if


    call setsolid (gl,fluids, nf, moles, resorted)


    if (.not. gl%same_components .OR. resorted) then
        call uppertolower_char(path, len(path))

        !if (trim(path_lower(1:2)) == 'hc') then
        !
        !    if(.not.allocated(gl%vpcoeff)) then
        !        allocate(gl%vpcoeff(20,gl%ncomp))
        !        allocate(gl%vpexp,gl%dlcoeff,gl%dlexp,gl%dvcoeff,gl%dvexp,gl%pmeltcoeff,gl%pmeltexp,gl%psubcoeff,gl%psubexp, mold=gl%vpcoeff)
        !    end if
        !
        !    !allocate the stn type to the size of components
        !    if(.not.allocated(gl%stn)) allocate(gl%stn(gl%ncomp))
        !
        !    !using hardcoded fluids
        !    do i=1,size(gl%Eq_type)
        !        if ((gl%Eq_type(i) /= 1).and.(gl%Eq_type(i) /= 2).and.(gl%Eq_type(i) /= 3).and.(gl%Eq_type(i) /= 4).and.(gl%Eq_type(i) /= 0)) then
        !            errorflag = -7892 !only Helmholtz, SRK, PR, LKP EOS are hardcoded at the moment
        !            exit
        !        elseif (gl%Eq_type(i) == 0) then
        !            exit
        !        endif
        !    enddo
        !    if (gl%ncomp > 1) then
        !        if ((gl%mix_type /= 1).and.(gl%mix_type /= 110).and.(gl%mix_type /= 111).and.(gl%mix_type /= 120).and.(gl%mix_type /= 121).and.(gl%mix_type /= 2).and.(gl%mix_type /= 3) .and. (gl%mix_type /= 4)) errorflag = -7892 !only Eqtype 1(HelmEOS),2(SRK),4(LKP) are hardcoded for mixtures
        !    endif
        !
        !    if (errorflag == 0) then
        !if ((gl%Eq_type(1) == 1) .and. ((gl%mix_type == 1) .or. (gl%mix_type == 110) .or. (gl%mix_type == 111) .or. (gl%mix_type == 120) .or. (gl%mix_type == 121) .or. (gl%mix_type == 1))) then
        !    if ((trim(path_lower) == 'hc\').or.(trim(path_lower) == 'hc/')) then
        !        call fluids_hardcoded(gl,errorflag)
        !    elseif ((trim(path_lower) == 'hc-gerg\').or.(trim(path_lower) == 'hc-gerg/')) then
        !        call fluids_hardcoded_gerg(gl,errorflag)
        !    elseif ((trim(path_lower) == 'hc-eoscg\').or.(trim(path_lower) == 'hc-eoscg/')) then
        !        call fluids_hardcoded_eoscg(gl,errorflag)
        !    endif
        !    gl%ref_set = .true.
        !elseif ((gl%Eq_type(1) == 2) .and. (gl%mix_type == 2)) then
        !    call fluids_srk_hardcoded(gl,errorflag)
        !elseif ((gl%Eq_type(1) == 3) .and. (gl%mix_type == 3)) then
        !    call fluids_pr_hardcoded(gl,errorflag)
        !elseif ((gl%Eq_type(1) == 4) .and. (gl%mix_type == 4)) then
        !    call fluids_lkp_hardcoded(gl,errorflag)
        !endif

        !    if ((((gl%Eq_type(1) == 2) .and. (gl%mix_type ==2)) .or. ((gl%Eq_type(1) == 3) .and. (gl%mix_type == 3)) .or. ((gl%Eq_type(1) == 4) .and. (gl%mix_type == 4))) .and. (errorflag == 0)) then
        !        if (gl%calc_ref) then
        !            call ref_state(gl,errorflag)   !reference point needed
        !        end if
        !    end if
        !
        !endif
        !if (errorflag /= 0) return

        !else

        if(gl%mix_type.eq.7) then !aga8
            call compare_fluids(gl, errorflag)
            !call check_limits(gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)
            if (errorflag /= 0) then
                return
            end if
        end if

        ! Check if it is necessary to re-read  fluid files
        if (.not.gl%same_components) then
            call read_from_fluidfiles (gl,path, errorflag)
        end if

        ! Calculate fluid specific properties (pc,tpp etc.)
        call calc_fld_specific_props(gl)


        !if error in read_from_fluidfiles: return
        if (errorflag /= 0) return

        !endif



        !New setup routine for one-fluid model, Andreas Jäger, July 2016
        if (gl%mix_type == 19) then
            call setup_one_fluid_model(gl)
        end if

        !New setup routine for gE models, Andreas, Erik, Dezember 2018
        if ((gl%mix_type == 12) .or. (gl%mix_type == 13)) then
            call setup_HelmgE(gl, errorflag)
            if (errorflag /= 0) then
                return
            end if
        end if

        ! save origianl (argument) path
        gl%path_arg = path

    endif

    if (gl%unitin==2) then
        converttype = 1
        call convert_fractions(gl, converttype, 1234.d0,  x_spec)           !if inout==1, then convert specific to molar; if inout ==2: convert molar to sp  !tine has to be defnied !!
        moles = x_spec
        gl%molfractions = x_spec
        gl%moleslist_hydrate = gl%molfractions
        if (input == 'td') then
            call convert_inputprop(gl, 1, prop2)
            !d = prop2
        elseif(input == 'ps') then
            call convert_inputprop(gl, 2, prop2)
            !s = prop2
        elseif(input == 'ph') then
            call convert_inputprop(gl, 2, prop2)
            !h = prop2
        end if
        if(gl%zero_comp) then  !!converting zero vector here!
            converttype = 1
            !comp_save = gl%ncomp
            !gl%ncomp = maxval(gl%comp_map)
            x_spec = gl%moles_zero_orig
            n_zero = maxval(gl%comp_map) - gl%ncomp
            call convert_fractions(gl, converttype, 1234.d0,  x_spec, n_zero)           !if inout==1, then convert specific to molar; if inout ==2: convert molar to sp  !tine has to be defnied !!
            gl%moles_zero_orig_spec = x_spec
            gl%zero_spec = .true.
            !gl%ncomp = comp_save
        end if
    end if




    if (errorflag /= 0) return
    ! "R"-Correction of the ideal and residual parts of the equations
    ! The residual part and the ideal part are corrected by multiplying
    ! with the fluid specific R and reducing with a universal R
    ! (see module_fluid_parameters for more detail)
    ! Andreas Nov. 2013
    ! NOT YET IMPLEMENTED!!!!!

    call reduced_parameters_calc(gl,300.D0) !Dummy temperature 300 K for the SRK

    ! same components: reference state is needed although it has not been set
    if ((gl%calc_ref) .and. (.not.gl%ref_set)) then
        call ref_state(gl,errorflag)   !reference point needed
    end if
    !range of validity of transport equations may be different to the EOS limits
    if (gl%transport /= 0) call set_limits_trans(gl)

    !if(gl%unitin == 2) moles = moles * gl%wm

    if(gl%seawater) then
        gl%sea%wm_sea = molar_sea(gl)
    end if

    end subroutine setup



    module subroutine setpropinput (gl,input, prop1, prop2, errorflag)






    implicit none

    type(type_gl) :: gl



    !define input variables
    character (12), intent(inout) :: input                      ! TD, DT, TP, PT, PH, HP, PS, SP
    double precision, intent(inout) :: prop1, prop2             ! input values
    integer, intent(out) :: errorflag                           ! information about warnings or errors


    integer:: i, j                                              ! count variables for the moles vector
    !define supportive variable
    double precision :: a, wm_mix                                       ! needed to switch property position (prop1, prop2)
    integer :: noc,m
    !character (26):: lower,upper
    logical :: LJunit                                           ! LJunit == .true.: reduced units (Lennard-Jones fluids) can be calculated
    integer :: ilim, iLJunit


    errorflag = 0
    iLJunit=0
    ilim=0


    !lower = 'abdhilmpqsutv'
    !upper = 'ABDHILMPQSUTV'

    gl%inptorig=input
    !Andreas, Feb. 2015
    gl%check_solid = .false.
    LJunit = .false.
    gl%hold_limits = .true.

    call uppertolower_char(input,len(input))

    noc = 0     !number of characters in input
    do m =1,12
        if (input(m:m) /= ' ') then
            if(input(m:m) == '&') then
                gl%hold_limits = .false.
                ilim=m
            end if
            !if(input(m:m) == '#') then
            !    LJunit = .true.
            !    iLJunit = m
            !end if
            !Andreas Feb 2014. Solids are indicated by "+", e.g. "tp+" indicates that also solid formation is checked at the specified T and p
            if(input(m:m) == '+') then
                gl%hold_limits = .false.             !Do NOT hold limits
                ilim = m                  !Save at which position the "+" stands
                gl%check_solid = .true.
            end if
            noc = noc + 1
        else
            exit
        end if
    end do

    if (gl%unitin==3) then !for calculating in reduced properties
        gl%Factor = 1.D0
        gl%factorpress=1.D6
        gl%factortrans=1.D0
        gl%factorrbwr=1.D0
        gl%factor_VS5eta=1.D0
        gl%factor_VS5slj=1.D0


    else                         !for calculating in unreduced properties
        gl%Factor = 1.D3
        gl%factorpress=1.D0
        gl%factortrans=1.D6
        gl%factorrbwr=1.D2
        gl%factor_VS5eta=1.D5
        gl%factor_VS5slj=1.D1

    end if

    if(gl%hold_limits .eqv. .false.) then
        input=input(1:ilim-1)//input(ilim+1:noc)
        if (iLJunit > ilim) iLJunit=iLJunit-1
        noc=noc-1
    end if

    !if(LJunit) then !for calculating with Lennard-Jones units
    !    gl%Factor = 1.D0
    !    gl%factorpress=1.D6
    !    gl%factortrans=1.D0
    !    gl%factorrbwr=1.D0
    !    gl%factor_VS5eta=1.D0
    !    gl%factor_VS5slj=1.D0
    !    input=input(1:iLJunit-1)//input(iLJunit+1:noc)
    !    noc=noc-1
    !else                         !for calculating real fluids
    !    gl%Factor = 1.D3
    !    gl%factorpress=1.D0
    !    gl%factortrans=1.D6
    !    gl%factorrbwr=1.D2
    !    gl%factor_VS5eta=1.D5
    !    gl%factor_VS5slj=1.D1
    !
    !end if


    !write lower cases in upper cases for the input variables
    !do i = 1,noc
    !    do j = 1,26
    !        if (input(i:i) == lower(j:j)) then
    !            input(i:i) = upper(j:j)
    !            exit
    !        end if
    !    end do
    !end do

    !Catch inputs TSUBV, TSUBS, PSUBV, PSUBS, TMLTL, TMLTS, PMLTL, PMTLS
    !No "+" necessary for these inputs, it is clear that a solid quantity is asked for
    !Andreas March 2014
    if ((input == 'tsubv') .or. (input == 'tsubs') &
        .or. (input == 'psubv') .or. (input == 'psubs') &
        .or. (input == 'tmltl') .or. (input == 'tmlts') &
        .or. (input == 'pmltl') .or. (input == 'pmlts')) then
        gl%hold_limits = .false.
        gl%check_solid = .true.
    end if

    !reduce input variable combinations
    if ((input == 'td') .or. (input == 'tp') .or. (input == 'ph') .or. (input == 'ps') &
        .or. (input == 'tliq') .or. (input == 'tvap') &
        .or. (input == 'pliq') .or. (input == 'pvap') &
        .or. (input == 'tliqd') .or. (input == 'tvapd') &
        .or. (input == 'pliqd') .or. (input == 'pvapd') &
        .or. (input == 'tq') .or. (input == 'pq') &
        .or. (input == 'tsubv') .or. (input == 'tsubs') &
        .or. (input == 'psubv') .or. (input == 'psubs') &
        .or. (input == 'tmltl') .or. (input == 'tmlts') &
        .or. (input == 'pmltl') .or. (input == 'pmlts')) then
        continue
        !switch DT
    elseif (input == 'dt') then
        input = 'td'
        a = prop1
        prop1 = prop2
        prop2 = a
        !switch PT
    elseif (input == 'pt') then
        input = 'tp'
        a = prop1
        prop1 = prop2
        prop2 = a
        !switch HP
    elseif (input == 'hp') then
        input = 'ph'
        a = prop1
        prop1 = prop2
        prop2 = a
        !switch SP
    elseif (input == 'sp') then
        input = 'ps'
        a = prop1
        prop1 = prop2
        prop2 = a
        !switch QT
    elseif (input == 'qt') then
        input = 'tq'
        a = prop1
        prop1 = prop2
        prop2 = a
        !switch QP
    elseif (input == 'qp') then
        input = 'pq'
        a = prop1
        prop1 = prop2
        prop2 = a
    else
        errorflag = -9955
        !write(*,*) 'wrong input combination'
    end if

    if ((input == 'ph').or.(input == 'ps')) then
        gl%calc_ref = .true.
    endif

    !catch T, D, P <= 0
    !T
    if ((index(input, 't') == 1) .and. (prop1 <= 0.d0)) then
        errorflag = -9911
        return
    end if

    !D
    if ((index(input, 'd') == 2) .and. (prop2 <= 0.d0)) then
        errorflag = -9922
        return
    end if

    !P
    if (((index(input, 'p') == 1) .and. (prop1 <= 0.d0)) .or. &
        ((index(input, 'p') == 2) .and. (prop2 <= 0.d0))) then
        errorflag = -9933
        return
    end if


    !converting input props not possible as molar weigth gets defined while reading fluid file! Benedikt, May 2018

    !save information of input combination and state point
    if ((gl%s_p%input /= input).or.(trim(input) /= "tp").or.(dabs(gl%s_p%prop(1) - prop1) > 1e-14).or.(dabs(gl%s_p%prop(2) - prop2) > 1e-14)) gl%s_p%fluid_present = .false.
    gl%s_p%input = input
    gl%s_p%prop(:) = (/prop1,prop2/)
    gl%inputprops(1) = prop1
    gl%inputprops(2) = prop2

    end subroutine setpropinput




    module subroutine setfluid (gl, fluids, nf, errorflag)



    implicit none

    type(type_gl) :: gl



    !define input variables
    !character (255), intent(inout) :: fluidl                ! fluid list
    integer, intent(out) :: errorflag                       ! information about warnings or errors
    integer:: j,t                                             ! count variable for the fluid vector
    character (30), dimension (30) :: fluids                ! fluid vector
    !character (255) :: fluidlwc                             ! fluid list working copy
    !define supportive variable
    !integer :: e,t,l                                        ! needed to transform fluid string to vector
    character(30) :: tf                                     ! needed to check for doubled fluids in fluidlist
    integer :: nf                                           ! amount of fluid in vector

    errorflag = 0
    nf = 0
    nf = count(fluids(:).ne.' ')
    !do j = 1, 30
    !    if (fluids(j) /= '') then
    !        nf = nf + 1
    !    elseif (fluids(j) == '') then
    !        exit
    !    endif
    !enddo
    !fluids = ''
    !
    !
    !fluidlwc = fluidl       ! write fluidlist in working copy
    !nf = 0                  ! initializing amount of fluids
    !
    !do j = 1, 30            ! maximum of 30 fluids allowed
    !    l = 0               ! amount of spaces in front of next fluids
    !    nf = nf+1           ! amount of fluids in list
    !
    !    e = index(fluidlwc, ';')    ! searches for position of separator ; in fluid list
    !
    !    if (e /= 0) then            ! a ; could be found in the list
    !
    !        ! check for amount of spaces in front of next fluid, max. 4 allowed
    !        do t = 1, 5
    !            if (fluidlwc(e+1:) == " ") then ! last fluid in list
    !                exit
    !            end if
    !            if (fluidlwc(e+l+1:e+l+1) == " ") then
    !                l=l+1   ! counting the spaces
    !            else
    !                exit    ! if there are none continue
    !            end if
    !            if (t == 5) then    ! maximum number of spaces reached, blank fluid - >  error in setup
    !                errorflag = -9960
    !                return
    !            end if
    !        end do
    !
    !        ! saving fluid in fluid vector, which is used for all further functions
    !            !Monika 04/2018
    !            !if (index(fluidlwc(1:e-1)," ") == 0) then       ! there are no blanks allowed in the fluidname
    !                fluids(j) = fluidlwc(1:e-1)
    !                fluidlwc = fluidlwc(e+1+l:)                 ! cutting current fluid from list
    !            !else
    !            !    errorflag = -9960                           ! fluidname contains blanks -- >  error
    !            !    return
    !            !end if
    !    else if ((e == 0) .and. (fluidlwc /= '')) then  ! string not empty, but no ; found -- >  last fluid in list
    !        fluids(j) = fluidlwc(1:26)
    !        fluidlwc = fluidlwc(e+1:)
    !        exit
    !    else        ! string empty and no ; found -- >  previous fluid was last fluid in list
    !        nf = nf-1
    !        exit
    !    end if
    !end do
    !
    !

    !
    do j = 1, nf
        call uppertolower_char(fluids(j),len(fluids(j)))
    enddo


    if(.not. gl%seawater) then
        gl%seawatercalled = .false.
        do j = 1, nf
            tf = trim(fluids(j))  ! each fluid is compared to the others
            if ((tf .eq. ('seawater'))) then
                if(.not. allocated(gl%sea)) allocate(gl%sea)
                if(j .ne. 1) then
                    errorflag = -9960
                    return
                end if
                gl%seawater = .true.
                gl%seacalc = .false.    !only set true, wenn seawater prop is needed !!
                gl%gecarrier = .false.
                gl%modelflag = 1
                gl%sea%wm_salt = 0.03140382d0 !set as parameter in module ?? S. Pohl ??
                !gl%check_solid = .true.
                gl%seapos = j
                fluids(j) = 'water'

            end if
        end do

        if((gl%seawater) .and. (.not. gl%seawatercalled)) then            !to make calls with gl_handle possible
            do j = 1, nf
                tf = fluids(j)  ! each fluid is compared to the others

                if (trim(tf) == 'salinity') then
                    gl%salpos = j
                    ! if(fluids(j) .ge. 1d-15)then    !checking if salinity is geater zero, else only water used!
                    gl%sal = .true.
                    fluids (j) = fluids(j+1)        !salinity not needed in futher steps, value stored in gl%sea%salinity in setmoles
                    !end if
                end if
                if((gl%sal) .and. (j .lt. 30) .and. (j .ge. gl%salpos) ) then
                    fluids(j) = fluids(j+1)
                elseif((gl%sal) .and. (j .eq. 30) .and. (j .ge. gl%salpos)) then
                    fluids(j) = ' '
                    !else
                    !   fluids(j) = fluids(j)
                end if
            end do
            !end if

            !if (gl%seawater) then               ! reducing number of fluids for seawater / salintiy
            if(.not. gl%sal) then
                errorflag = -9960
                return
            end if
            nf = nf-1
        end if
        !gl%seawatercalled = .true.
        !else
        !gl%seawatercalled = .true.
    end if
    ! check whether each fluid is only set once

    !sorting for electrolyte model
    !Benedikt September 2019
    if(any(gl%eq_type(:) == 15)) then
        call electrolytes_setfluids(gl,fluids, nf)
        fluids = fluids
        nf = nf
    end if

    do j = 1, nf
        tf = fluids(j)  ! each fluid is compared to the others
        do t = j+1, nf
            if (tf == fluids(t)) then
                errorflag = -9960
                return
            end if
        end do

    end do
    !
    !Andreas July 2014
    ! set number of given components
    gl%ncomp = nf
    gl%components = fluids

    if ((gl%eq_type(1) .eq. 8) .or. (gl%eq_type(1) .eq. 81)) then
        call check_fld_RKM(gl,errorflag)
    end if

    end subroutine setfluid


    !------------------------------------------------------------
    subroutine check_fld_RKM(gl,errorflag)

    implicit none

    type(type_gl) :: gl
    character*40,dimension(20):: fld_list=(/'methane','methan','ethane','ethan','propane','propane','butane','n-butane','butan','isobutane','ibutane','isobutan','nitrogen','stickstoff','pentane','n-pentane','pentan','isopentane','isopentan','ipentane'/)
    integer:: wrong_fld, i, errorflag

    !check if the specified eq_type is valid
    !list of valid eq_types is defined in gl%m_info%eq_type_list

    do i = 1,gl%ncomp
        wrong_fld = 0
        wrong_fld = all(.not.gl%components(i) .eq. fld_list)
        if (wrong_fld .ne. 0) then
            errorflag = -8002
            return
        end if
    end do

    end subroutine


    module subroutine setmoles (gl, moles, nm, errorflag)






    implicit none

    type(type_gl) :: gl



    !define input variables
    !character (255), intent(inout) :: molesl                      ! fluid list
    double precision, dimension (30), intent(inout) :: moles      ! molefraction vector
    integer, intent(out) :: errorflag                        ! information about warnings or errors


    integer:: nm, j, k, l, nf                                           ! count variables for the moles vector
    !character (255) :: moleslwc                                 ! comp list working copy
    double precision:: sum                                      ! sum of all mole fractions
    !define supportive variable
    !integer :: e,t,l                                          ! needed to transform fluid string to vector
    !character(12) :: dd                                         ! needed to add "." if it's missing
    !integer :: b
    !integer :: c                                              ! amount of fluid in vector
    integer :: error

    errorflag = 0
    nm = 0
    nm = count(moles(1:gl%ncomp) >=0.d0 )

    !! check if the given mole fractions are reasonable
    !moleslwc = molesl       ! overwrite complist to working copy
    !do j = 1, 30            ! maximum of 30 fluids allowed
    !    l = 0               ! amount of spaces in front of next fluids
    !    nm = nm+1           ! amount of fluids in list
    !
    !    e = index(moleslwc, ';')    ! searches for position of separator ; in comp list
    !
    !    if (e /= 0) then            ! a ; could be found in the list
    !
    !        ! check for amount of spaces in front of next comp, max. 4 allowed
    !        do t = 1, 5
    !            if (moleslwc(e+1:) == " ") then ! last comp in list
    !                exit
    !            end if
    !            if (moleslwc(e+l+1:e+l+1) == " ") then
    !                l=l+1   ! counting the spaces
    !            else
    !                exit    ! if there are none continue
    !            end if
    !            if (t == 5) then    ! maximum number of spaces reached, blank fluid - >  error in setup
    !                errorflag = -9954
    !            end if
    !        end do
    !
    !        ! saving comp in composition vector, which is used for all further functions
    !        if (index(moleslwc(1:e-1)," ") == 0) then       ! there are no blanks allowed in the composition
    !            ! this is basically not needed here, since it could only be used for entry 1 and 0, and 0 is forbidden, 1 is treated below
    !            if (index(moleslwc(1:e-1),".") /= 0) then   ! a "." is needed in order to read number as double
    !                read(moleslwc(1:e-1),*, iostat=error) moles(j)
    !            else
    !                dd = moleslwc(1:e-1)//"."
    !                read(dd(1:e),*, iostat=error) moles(j)
    !            end if
    !            moleslwc = moleslwc(e+1+l:)                 ! cutting current comp from list
    !        else
    !            errorflag = -9954                           ! composition contains blanks -- >  error
    !            return
    !        end if
    !
    !        if (error /= 0) then
    !            errorflag = -9954
    !            return
    !        end if
    !
    !    else if ((e == 0) .and. (moleslwc /= '')) then  ! string not empty, but no ; found -- >  only composition in list
    !        b = index(moleslwc," ")
    !        if (index(moleslwc(1:b),".") /= 0) then   ! a "." is needed in order to read number as double
    !            read(moleslwc(1:b),*, iostat=error) moles(j)
    !        else
    !            dd = moleslwc(1:b)//"."
    !            read(dd(1:b+1),*, iostat=error) moles(j)
    !        end if
    !        if (error /= 0) then
    !            errorflag = -9954
    !            return
    !        end if
    !        moleslwc = moleslwc(e+1:)
    !        exit
    !    else        ! string empty and no ; found -- >  previous composition was last composition in list
    !        nm = nm-1
    !        exit
    !    end if
    !end do

    ! check if all entries are within 0 < x_i < 1 and fill the vector to the size of 30 with zeros
    sum = 0.d0


    if((gl%sal) .and. (.not. gl%seawatercalled)) then!.and. (.not. gl%seawatercalled)) then
        gl%sea%salinity = moles(gl%salpos)
        if (moles(gl%salpos) .gt. 0.12d0) then
            errorflag = -12800
        end if

        if(gl%sea%salinity .le. 1.d-9) then     !
            gl%seawater = .false.           ! only water
            !elseif(gl%sea%salinity .eq. 99.d0) then
            !    gl%seawater = .true.                !for consistency check only!!
            !    gl%sea%salinity = 0.d0
        end if

        moles(gl%salpos) = 0.d0

        do l = (gl%salpos), nm
            if((gl%sal) .and. (l .lt. 30)) then
                moles(l) = moles(l+1)
            elseif((gl%sal) .and. (l .eq. 30)) then
                moles(l) = 0.d0
            end if
        end do
        !nm = nm - 1
        !gl%seawatercalled = .true.
    end if


    if(gl%el_present) then
        call  electrolytes_adjustmoles(gl, moles, nm)
        !gl%mix_type=15
    end if


    do k = 1, nm
        !These lines caused trouble, because the criterion that the mole fraction must be greater than 10^-11 are below double precision which means we cannot specify mixtures with smaller amount of impurities.
        !This might be needed for some calculations, so the criterion was changed here to 1.d-15. Andreas Jäger, May 2017
        !if (((1.d-11 >= moles(k)) .OR. (moles(k) > 1.d0)) .AND. (k == 1)) then    !check if 0 < xi <= 1
        if (((1.d-15 >= moles(k)) .OR. (moles(k) > 1.d0)) .AND. (k == 1)) then    !check if 0 < xi <= 1
            errorflag = -9951   ! there is a valid molefraction given in the first entry
            return
        end if
        !if (((1.d-11 > moles(k)) .OR. (moles(k) >= 1.d0)) .AND. (k > 1)) then    !check if 0 < xi < 1
        if (((moles(k) < 1.d-15) .OR. (moles(k) >= 1.d0)) .AND. (k > 1)) then    !check if 0 < xi < 1
            errorflag = -9951
            return
        end if

        sum = sum + moles(k)

    end do





    if (abs(sum - 1.d0) > 1.d-10) then  !check if sum /= 1
        errorflag = -9952
        return
    end if

    gl%molfractions = moles
    gl%moles_read = moles

    end subroutine setmoles


    !**************
    module subroutine setseawater(gl)




    implicit none

    type(type_gl) :: gl

    !MBialdyga, May 2018
    !Handling of Seawater:
    !Checking Seawater Flag
    if (gl%seawater) then
        call seawater_limits(gl)
        call set_coef_Seawater(gl)
        gl%sea%wm_sea = molar_sea(gl)
        gl%modelflag = 1
        gl%gemix = .false.
    endif

    end subroutine setseawater


    module subroutine setsolid(gl,fluids, nf, moles, resorted)






    implicit none

    type(type_gl) :: gl


    integer :: nf                                               ! amount of fluid in vector
    integer :: j, t
    character (30), dimension (30) :: fluids                    ! fluid vector
    double precision, dimension (30), intent(out) :: moles      ! molefraction vector
    character(30) :: tf                                         ! needed to check for doubled fluids in fluidlist
    integer, dimension(30) :: eqtype                            ! needed to resort Eq_Type
    logical :: resorted                                         ! flag whether fluid list has been resorted  Theresa
    integer :: hydrate_structure_flag


    resorted = .false.

    !Andreas, May 2012
    !Handling of solids:
    !Go through the list of components and check whether solid equations are available for the components selected
    if (gl%check_solid) then

        gl%nrofhydrateformers = 0
        gl%nrofsolids = 0

        !Check the existance and position of all components for which hydrate formation is possible AND equations are available (including water)
        do j = 1, count(gl%Hydrate_list(:) /= "")
            tf = gl%Hydrate_list(j)  ! each fluid is compared to the others
            if(index(tf," ") .eq. 1) exit
            if(gl%nrofhydrateformers .eq. nf) exit !Theresa
            do t = 1, nf
                if (tf .eq. fluids(t)) then
                    gl%nrofhydrateformers = gl%nrofhydrateformers + 1
                    gl%Hydrate_formers(gl%nrofhydrateformers) = fluids(t)
                    gl%moleslist_hydrate(gl%nrofhydrateformers) = moles(t)
                    !Hydrate_pos(nrofhydrateformers) = t
                    gl%mapping(gl%nrofhydrateformers) = t
                end if
            end do
            !if((j .eq. 1) .and. (nrofhydrateformers .eq. 0)) exit !No water found, no hydrate formation possible
        end do
        !For easier calculations of hydrates, the fluid vector and composition vector is mapped
        !If water and other hydrate forming substances are present in the mixture, internally the fluid vector has a certain structure as follows:
        !Fluidlist_hydrate = "water"                    -- .gt.  Always position 1, counts as hydrate former 1
        !                    "hydrate former 2"
        !                    "Hydrate former 3"
        !                    .
        !                    "Hydrate former m"
        !                    "Additional component 1"
        !                    "Additional component 2"
        !                    .
        !                    "Additional component n-m"
        gl%Fluidlist_hydrate = gl%Hydrate_formers

        do j = 1, nf
            tf = fluids(j)  ! each fluid is compared to the others
            t = 1
            do while (t .lt. nf+1)
                if (index(gl%Fluidlist_hydrate(t)," ") .eq. 1) then
                    gl%Fluidlist_hydrate(t) = tf
                    gl%moleslist_hydrate(t) = moles(j)
                    gl%mapping(t) = j
                    exit
                else
                    if (tf .eq. gl%Fluidlist_hydrate(t)) then
                        exit
                    else
                        t = t + 1
                    end if
                end if
            end do
        end do

        !At least water and one hydrate forming component is needed in order to form hydrates (nrofhydrateformers .gt. 1)
        if ((gl%nrofhydrateformers .gt. 1) .and. (gl%Fluidlist_hydrate(1) .eq. "water")) then !Hydrate formation possible
            gl%solidtype(2) = 1
        else
            gl%nrofhydrateformers = 0
        End if


        !Check the existence and position of all components for which solid formation equations are available
        do j = 1, count(gl%solid_list(:) /= "")
            tf = gl%solid_list(j)  ! each fluid is compared to the others
            if(index(tf," ") .eq. 1) exit
            do t = 1, nf
                if (tf .eq. gl%Fluidlist_hydrate(t)) then
                    gl%nrofsolids = gl%nrofsolids + 1
                    gl%solid_indicator(j) = t
                    exit
                elseif (trim(gl%Fluidlist_hydrate(j)) ==  '7732-18-5')then
                    gl%nrofsolids = gl%nrofsolids + 1
                    gl%solid_indicator(j) = t
                    gl%Fluidlist_hydrate(j) = 'water'
                    exit
                elseif (trim(gl%Fluidlist_hydrate(j)) ==  '124-38-9') then
                    gl%nrofsolids = gl%nrofsolids + 1
                    gl%solid_indicator(j) = t
                    gl%Fluidlist_hydrate(j) = 'co2'
                    exit
                end if
            end do
        end do

        !Check, whether fluidlist needs to be resorted
        do j = 1, nf
            if (gl%mapping(j) .ne. j) then
                resorted = .TRUE.
                call initialize_ideal_gas_coefficients(gl)
            end if
        end do

        if (resorted) then
            !Change the fluid vector according to sorted list for solids
            gl%components = gl%Fluidlist_hydrate
            gl%molfractions = gl%moleslist_hydrate
            fluids = gl%Fluidlist_hydrate
            moles = gl%moleslist_hydrate
            ! resort Eq_Type
            eqtype = gl%Eq_type
            do j = 1, nf
                gl%Eq_type(j) = eqtype(gl%mapping(j))
            end do
        end if

        !Set fixed solid former if a pure substance is calculated
        if (nf .eq. 1) then
            if (gl%Fluidlist_hydrate(1) .eq. "water") then
                gl%solidtype(1) = 1
                gl%solidtype_akt_phase = 1
            elseif (gl%Fluidlist_hydrate(1) .eq. "co2") then
                gl%solidtype(1) = 2
                gl%solidtype_akt_phase = 2
            end if
            gl%solid_pos = 1
            gl%solidpos_akt_phase = 1
        end if

        !Store the quadruple point information for binary mixtures forming hydrate
        if (gl%nrofhydrateformers .eq. 2) then
            gl%p_Q_Hyd_2C = 0.D0
            gl%T_Q_Hyd_2C = 0.D0
            select case (trim(gl%Fluidlist_hydrate(2)))
            case ("co2")
                gl%p_Q_Hyd_2C(1) = 0.970952D0      !VLwHIw
                gl%p_Q_Hyd_2C(2) = 4.480149D0      !VLwLcH
                gl%p_Q_Hyd_2C(3) = 0.517772D0      !VLcHIc
                gl%p_Q_Hyd_2C(4) = 486.743236D0    !LwLcHIc
                gl%T_Q_Hyd_2C(1) = 271.1885D0      !VLwHIw
                gl%T_Q_Hyd_2C(2) = 283.0787D0      !VLwLcH
                gl%T_Q_Hyd_2C(3) = 216.5869D0      !VLcHIc
                gl%T_Q_Hyd_2C(4) = 292.9047D0      !LwLcHIc
                hydrate_structure_flag = 1
            case ("methane")
                gl%p_Q_Hyd_2C(1) = 2.549980D0      !VLwHIw
                gl%T_Q_Hyd_2C(1) = 272.9695D0      !VLwHIw
                hydrate_structure_flag = 1
            case ("ethane")
                gl%p_Q_Hyd_2C(1) = 0.461073D0      !VLwHIw
                gl%p_Q_Hyd_2C(2) = 3.244851D0      !VLwLeH
                gl%T_Q_Hyd_2C(1) = 272.9278D0      !VLwHIw
                gl%T_Q_Hyd_2C(2) = 287.3880D0      !VLwLeH
                hydrate_structure_flag = 1
            case ("co")
                gl%p_Q_Hyd_2C(1) = 11.559560D0     !VLwHIw
                gl%T_Q_Hyd_2C(1) = 272.1444D0      !VLwHIw
                hydrate_structure_flag = 1
            case ("nitrogen")
                gl%p_Q_Hyd_2C(1) = 14.159368D0     !VLwHIw
                gl%T_Q_Hyd_2C(1) = 271.8518D0      !VLwHIw
                hydrate_structure_flag = 2
            case ("oxygen")
                gl%p_Q_Hyd_2C(1) = 10.882913D0     !VLwHIw
                gl%T_Q_Hyd_2C(1) = 272.0910D0     !VLwHIw
                hydrate_structure_flag = 2
            case ("propane")
                gl%p_Q_Hyd_2C(1) = 0.167532D0      !VLwHIw
                gl%p_Q_Hyd_2C(2) = 0.557926D0      !VLwLpH
                gl%T_Q_Hyd_2C(1) = 273.1476D0      !VLwHIw
                gl%T_Q_Hyd_2C(2) = 278.5100D0      !VLwLpH
                hydrate_structure_flag = 2
            case ("argon")
                gl%p_Q_Hyd_2C(1) = 8.583812D0      !VLwHIw
                gl%T_Q_Hyd_2C(1) = 272.3696D0      !VLwHIw
                hydrate_structure_flag = 2
                case default

            end select
        elseif (gl%nrofhydrateformers > 2) then
            if ((trim(gl%Fluidlist_hydrate(2)) == 'co2').or.(trim(gl%Fluidlist_hydrate(2)) == 'methane').or.(trim(gl%Fluidlist_hydrate(2)) == 'ethane').or.(trim(gl%Fluidlist_hydrate(2)) == 'co') ) then
                hydrate_structure_flag = 1 !set arbitrary value
            elseif ((trim(gl%Fluidlist_hydrate(2)) == 'nitrogen').or.(trim(gl%Fluidlist_hydrate(2)) == 'oxygen').or.(trim(gl%Fluidlist_hydrate(2)) == 'propane').or.(trim(gl%Fluidlist_hydrate(2)) == 'argon') ) then
                hydrate_structure_flag = 2 !set arbitrary value
            else
                hydrate_structure_flag = 3
            endif
        else !only water or co2
            hydrate_structure_flag = 0
        end if

        call hdrt_structure_definition(gl,hydrate_structure_flag)

        !If solids are considered then always calculate the reference point
        gl%calc_ref = .true.

    end if
    end subroutine setsolid

    module subroutine seteqtype (gl, eqtype, mixtype, ne, errorflag)




    implicit none

    type(type_gl) :: gl



    !define input variables
    !character (255), intent(inout) :: EOS_indicator
    integer, intent(out) :: errorflag                        ! information about warnings or errors


    integer:: j, i                                           ! count variables for the moles vector
    !character (255) :: eos_indicatorwc                           ! eos indicator working copy
    !define supportive variable
    integer :: e,t,l                                          ! needed to transform fluid string to vector
    integer :: b
    integer :: ne                                              ! amount of fluid in vector
    integer :: etype
    integer :: error
    integer, dimension(30), intent(inout) :: eqtype
    integer, intent(inout):: mixtype
    integer:: wrong_eq


    errorflag = 0
    ne = 0
    ne = count(eqtype(:).ne.0)

    gl%mix_type = mixtype
    gl%Eq_type = eqtype

    !check if the specified eq_type is valid
    !list of valid eq_types is defined in gl%m_info%eq_type_list
    do i = 1,ne
        wrong_eq = 0
        wrong_eq = all(.not.gl%Eq_type(i) .eq. gl%m_info%eq_type_list)
        if (wrong_eq .ne. 0) then
            errorflag = -9957
            return
        end if
    end do

    !j=0
    !ne=0
    !e=0
    !t=0
    !l=0
    !b=0
    !etype=0
    !error=0
    !eqtype=0
    !mixtype=0
    !
    !! Andreas, March 2013
    !!---------------------------------------------------------------------------
    !!-----------------------------   EOSTYPE  ----------------------------------
    !! check if the given mole fractions are reasonable
    !eos_indicatorwc = eos_indicator       ! overwrite eos indicator to working copy
    !do j = 1, 21    ! maximum of 20 fluids allowed plus mixing rules = 21 entries
    !    l = 0       ! amount of spaces in front of next element
    !    etype = etype+1     ! count up amount of entries in list
    !    e = index(eos_indicatorwc, ';')    ! searches for position of separator ; in eos indicator list
    !    if (e /= 0) then    ! a ; could be found in the list
    !
    !        ! check for amount of spaces in front of next element, max. 4 allowed
    !        do t = 1, 5
    !            if (eos_indicatorwc(e+1:) == " ") then ! last comp in list
    !                exit
    !            end if
    !            if (eos_indicatorwc(e+l+1:e+l+1) == " ") then
    !                l=l+1   ! counting the spaces
    !            else
    !                exit    ! if there are none continue
    !            end if
    !            if (t == 5) then    ! maximum number of spaces reached, blank fluid - >  error in setup
    !                errorflag = -9957
    !            end if
    !        end do
    !        ! saving EOStype in Eqtype vector
    !        if (index(eos_indicatorwc(1:e-1)," ") == 0) then       ! there are no blanks allowed in the Eqtype (read between two ";")
    !            ! since the eos type is an integer, no "." or "," is allowed
    !            if ((index(eos_indicatorwc(1:e-1),".") == 0) .and. (index(eos_indicatorwc(1:e-1),",") == 0)) then
    !                read(eos_indicatorwc(1:e-1),*, iostat=error) gl%Eq_type(j)
    !                if (error /= 0) then
    !                    errorflag = -9959
    !                    return
    !                end if
    !            else
    !                errorflag = -9957                           ! Eostype contains blanks (like 1;2;1 2;1) -- >  error
    !                return
    !            end if
    !            eos_indicatorwc = eos_indicatorwc(e+1+l:)                 ! cutting current entry from list
    !        else
    !            errorflag = -9957                           ! Eostype contains blanks (like 1;2;1 2;1) -- >  error
    !            return
    !        end if
    !
    !    else if ((e == 0) .and. (eos_indicatorwc /= '')) then  ! string not empty, but no ; found -- >  last entry in list
    !        b = index(eos_indicatorwc," ")
    !        ! since the eos type is an integer, no "." or "," is allowed
    !        if ((index(eos_indicatorwc(1:b),".") == 0) .and. (index(eos_indicatorwc(1:b),",") == 0)) then
    !            read(eos_indicatorwc(1:b),*, iostat=error) gl%Mix_type ! The last entry in the vector indicates the mixture model type
    !            if (error /= 0) then
    !                errorflag = -9959
    !                return
    !            end if
    !            if (etype == 1) then    ! Theresa, January 2015: save number of given eqtypes
    !                gl%Eq_type(1) = gl%Mix_type ! If the eostype has only one entry and the number of fluids = 1, Eqtype(1) = Mix_type
    !                ne = 1
    !                exit
    !            end if
    !            ne = j-1 ! save number of given eqtypes
    !        else
    !            errorflag = -9957                           ! Eostype contains blanks (like 1;2;1 2;1) -- >  error
    !            return
    !        end if
    !        eos_indicatorwc = eos_indicatorwc(e+1:)
    !        exit
    !    else        ! string empty and no ; found -- >  previous fluid was last fluid in list
    !        etype = etype-1
    !        gl%Mix_type = gl%Eq_type(etype)
    !        exit
    !    end if
    !end do
    !For pure fluids, set Mix_type according to the equation used
    !if (ne == 1) then
    !    if (etype == 0) then
    !        errorflag = -9957
    !        return
    !    end if
    !    gl%Mix_type = gl%Eq_type(1)
    !end if

    !!error catching: PC-SAFT does not exist
    !if ((eq_type(1) == 6) .or. (Mix_type == 6)) then
    !    errorflag = -9957
    !end if

    if ((gl%Mix_type == 51) .or. (gl%Mix_type == 52) .or. (gl%Mix_type == 53)) gl%mix_type = 120  !Lorentz-Berthelot combining rule used for generalized EOS

    !! check if the eostype vector is of the same length as the compositions vector and if the mixture model is given (not needed for pure fluids)
    !if ((f > 1) .and. (etype < (f+1))) then
    !    errorflag = -9957
    !    return
    !end if

    ! STEFAN Feb 2014
    ! Check valid mixture model. Helmholtz EOS (1) can only be mixed with Helmholtz mixing rules (1). SRK(gl,2) and PR (3) mixing rules can only be used with SRK(gl,2) or PR (3) EOS, respectively. However, cubic and Helmholtz EOS may be mixed with Helmholtz mixing rules since the cubic EOS are integrated to the Helmholtz energy
    if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21)  .or. (gl%Mix_type == 22)) then !Andreas Feb 2016: Mix_type = 21: Quadratic mixing rules for b (SRK), Mix_type = 22: PSRK
        do j = 1, ne
            if (gl%Eq_type(j) /= 2) then
                errorflag = -9958
                return
            end if
        end do
    elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then !Andreas Feb 2016: Mix_type = 31: Quadratic mixing rules for b (PR)
        do j = 1, ne
            if (gl%Eq_type(j) /= 3) then
                errorflag = -9958
                return
            end if
        end do
    elseif (gl%Mix_type == 4) then
        do j = 1, ne
            if (gl%Eq_type(j) /= 4) then
                errorflag = -9958
                return
            end if
        end do
    elseif ((gl%Mix_type == 51) .or. (gl%Mix_type == 52) .or. (gl%Mix_type == 53)) then
        do j = 1, ne
            if ((gl%Eq_type(j) /= 51) .and. (gl%Eq_type(j) /= 52) .and. (gl%Eq_type(j) /= 53)) then
                errorflag = -9958
                return
            end if
        end do
    end if
    !Erik, July 2019, allow mixing of EOS with multi fluid mixture model with standard mixing rules
    !elseif ((gl%Mix_type == 110) .or. (gl%Mix_type == 111) .or. (gl%Mix_type == 121)) then  !fixed mixing rules if binary mix file not available
    !    do j = 1, ne
    !        if (gl%Eq_type(j) /= 1) then
    !            errorflag = -9958
    !            return
    !        end if
    !    end do
    !elseif (gl%Mix_type == 120) then  !fixed mixing rule (Lorentz-Berthelot) if binary mix file not available
    !    do j = 1, ne
    !        if ((gl%Eq_type(j) /= 1) .and. (gl%Eq_type(j) /= 51) .and. (gl%Eq_type(j) /= 52) .and. (gl%Eq_type(j) /= 53)) then
    !            errorflag = -9958
    !            return
    !        end if
    !    end do
    !end if
    !Check if the Eq_type exists. (At the moment: 1: Helmholtz, 2: SRK, 3: PR, 4:LKP, 5:generalized EOS, 6: PC-SAFT, 7: AGA8, 9:Costald) All other inputs cause errors
    if (gl%Mix_type < 1) then
        errorflag = -9957
        return

    elseif ((gl%Mix_type > 4) .and. (gl%Mix_type /= 6) .and. (gl%Mix_type /= 7) .and. (gl%Mix_type /= 8) .and. (gl%Mix_type /= 81) .and. (gl%Mix_type /= 9) .and. (gl%Mix_type /= 51) .and. (gl%Mix_type /= 52) .and. (gl%Mix_type /= 53)) then

        !Andreas, September 2015, Feb 2016; Moni June 2016
        if ((gl%Mix_type /= 11) .and. (gl%Mix_type /= 12) .and. (gl%Mix_type /= 13) .and. (gl%Mix_type /= 21) .and. (gl%Mix_type /= 22) .and. (gl%Mix_type /= 31) .and. &                !Mix_type = 11,12,13: Different mixing rules for Helmholtz mixtures, Mix_type = 21: Quadratic mixing rules for b (SRK), Mix_type = 31: Quadratic mixing rules for b (PR) , Mix_type = 22: PSRK
            & (gl%Mix_type /= 110) .and. (gl%Mix_type /= 111) .and. (gl%Mix_type /= 120) .and. (gl%Mix_type /= 121)) then     !fixed mixing rules if binary mix file not available
            errorflag = -9957
            return
        end if

        !AGA 8 - check if mixtype and eqtype fit
    elseif(gl%mix_type.eq.7) then
        !errorflag = -9957
        if(any(gl%eq_type(1:ne).ne.7)) then
            errorflag = -9957
            return
        end if
    end if

    do j = 1, ne
        if (gl%Eq_type(j) < 1) then
            errorflag = -9957
            return

        elseif ((gl%Mix_type > 4) .and. (gl%Mix_type /= 6) .and. (gl%Mix_type /= 7) .and. (gl%Mix_type /= 8) .and. (gl%Mix_type /= 81) .and. (gl%Mix_type /= 9) .and. (gl%Mix_type /= 51) .and. (gl%Mix_type /= 52) .and. (gl%Mix_type /= 53)) then

            !Andreas, September 2015, Feb 2016
            if ((gl%Mix_type /= 11) .and. (gl%Mix_type /= 12) .and. (gl%Mix_type /= 13) .and. (gl%Mix_type /= 21) .and. (gl%Mix_type /= 22) .and. (gl%Mix_type /= 31) .and. &     !Mix_type = 11,12,13: Different mixing rules for Helmholtz mixtures, Mix_type = 21: Quadratic mixing rules for b (SRK), Mix_type = 31: Quadratic mixing rules for b (PR) , Mix_type = 22: PSRK
                & (gl%Mix_type /= 110) .and. (gl%Mix_type /= 111) .and. (gl%Mix_type /= 120) .and. (gl%Mix_type /= 121)) then     !fixed mixing rules if binary mix file not available
                errorflag = -9957
                return
            end if
        end if
    end do

    !IN THE ACTUAL VERSION (TREND 1.1) MIXING OF CUBIC AND HELMHOLTZ EOS IS NOT ALLOWED
    !ANDREAS Nov 2013
    !if (Mix_type == 1) then
    !    do j = 1, ne
    !        if (Eq_type(j) /= 1) then
    !            errorflag = -9958
    !            return
    !        end if
    !    end do
    !end if

    !Solids can only be checked if (certain!!!) Helmholtz equations are used
    !ANDREAS Feb 2014
    if (gl%check_solid) then
        do j = 1, ne
            if ((gl%components(j) == "co2") .or. (gl%components(j) == "water")) then    !Andreas, Erik Juni 2019, at the moment only equations for solid water and solid co2 implemented
                if (gl%Eq_type(j) /= 1) then
                    if(gl%eq_type(j) /= 15) then    !exception for electrolytes
                        errorflag = -19900
                        return
                    end if
                end if
            end if
        end do
        if ((gl%Mix_type == 1).or.(gl%Mix_type == 110).or.(gl%Mix_type == 111).or.(gl%Mix_type == 120).or.(gl%Mix_type == 121).or.(gl%Mix_type == 12).or.(gl%Mix_type == 13)) then
            continue
        else
            errorflag = -19900
            return
        end if
    end if

    if(.not. gl%seawatercalled) then
        if(any(gl%components == 'seawater')) gl%sal = .true.
        if((gl%sal) ) then !.and. (.not. gl%seawatercalled)
            ne = ne - 1

        end if
    end if

    if(any(eqtype(:) == 15)) then
        if(.not. allocated(gl%el)) allocate(gl%el)
        call initialize_electrolytes(gl)
        call electrolytes_seteq(gl,eqtype, ne, errorflag)
        !call electrolytes_setfluids(gl, gl%components, gl%ncomp)
        ne = ne
    end if

    !eqtype = gl%Eq_type
    !mixtype = gl%mix_type
    end subroutine seteqtype


    module subroutine setpath (gl,path)


    implicit none

    type(type_gl) :: gl


    !define input variables
    character (255), intent(inout) :: path                                     ! path where to find fluid and mix files
    integer :: lenpath

    lenpath = len(trim(path))
    if (lenpath == 0) lenpath = lenpath + 1
    if (lenpath == 255) lenpath = lenpath - 1
    !changed for compatibility with gfortran
    !if (path(lenpath:lenpath) /= "\") path(lenpath+1:lenpath+1) = "\"
    if ( (path(lenpath:lenpath) /= "\") .and. (path(lenpath:lenpath) /= "/")) path(lenpath+1:lenpath+1) = "/"
    !---------------------------------------------------------------------------

    end subroutine setpath



    !  ---------------  -----------------------------------
    !  Sub for calculating the Molmass of Mixtures
    !  J. Moeller, Denmark, 09-2009
    !  --------------------------------------------------


    module subroutine wm_mix_calc(gl,wm_mix)


    implicit none

    type(type_gl) :: gl


    double precision :: wm_mix !Output of Subroutine

    integer :: j
    j=1

    wm_mix = 0.d0

    do j=1,gl%ncomp
        wm_mix = wm_mix + gl%molfractions(j)*gl%wm(j)
    enddo

    end subroutine wm_mix_calc
    !-------------------------------------------------------------------------

    !    !  --------------------------------------------------
    !  Sub for calculating the Molmass of Reactive Mixtures
    !  B. S. 10/2018
    !  --------------------------------------------------
    module subroutine wm_mix_reac_calc(gl, wm_comps, molfrac, wm_reac_mix)


    implicit none

    type(type_gl) :: gl

    double precision :: wm_reac_mix
    double precision, dimension(30), intent(in) :: wm_comps, molfrac
    integer :: j

    wm_reac_mix = 0.d0

    do j=1, gl%ncomp
        wm_reac_mix = wm_reac_mix + molfrac(j) * wm_comps(j)
    end do

    end subroutine wm_mix_reac_calc
    !-------------------------------------------------------------------


    ! --------------------------------------------------
    ! Suroutine for the calculation of the Universal Gas
    ! Constant R for mixtures. This is neccessary because
    ! different pure fluid equations use different values
    ! for R
    ! J. Gernert, Aug. 2010
    ! --------------------------------------------------
    module subroutine R_mix_calc(gl,R_mix)



    implicit none

    type(type_gl) :: gl


    double precision :: R_mix !Output of Subroutine

    integer :: j
    j=1

    R_mix = 0.d0

    do j=1,gl%ncomp
        R_mix = R_mix + gl%molfractions(j)*gl%REQ(j)
    enddo

    end subroutine R_mix_calc


    module subroutine ref_state(gl,errorflag)







    implicit none

    type(type_gl) :: gl

    integer :: errorflag                        ! information about warnings or errors
    integer :: nrsubst, ref_error

    !calculate rhoref (and pref in some cases)
    flloop: do nrsubst = 1, gl%ncomp

        call ref_calc(gl,nrsubst, ref_error)
        if (ref_error == 0) then

            !calculate c1 and c2:
            IF (gl%cpmodel(nrsubst) .OR. (gl%EQ_TYPE(nrsubst) == 2) .OR. (gl%EQ_TYPE(nrsubst) == 3) .OR. (gl%EQ_TYPE(nrsubst) == 4) .OR. (gl%EQ_TYPE(nrsubst) == 6)) THEN
                call idealconsts(gl,nrsubst)
            elseif (allocated(gl%eos_coeff)) then
                if (gl%eos_coeff%cppcheck(nrsubst) == 'PHK') then
                    call idealconsts(gl,nrsubst)
                end if
            endif

            gl%ref_set = .true.
        else
            errorflag = ref_error
            return
        end if

    enddo flloop

    end subroutine ref_state




    !subroutine R_correct()
    !



    !
    !n = 0.D0    !Linear coefficients of the residual part
    !dfn = 0.D0  !Linear coefficients of the departure function
    !cp0coeff    !Linear coefficients for the ideal gas part
    !
    !return
    !end subroutine



    !  ----------------------------------------------------------------------------------
    !  this subroutine checks whether the range of validity of the EOS is maintained
    !  T. Wiens, April 2012
    !  ----------------------------------------------------------------------------------


    module subroutine check_limits (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)





    implicit none

    type(type_gl) :: gl


    character(12) :: input
    double precision :: t, d, dliq, dvap, p
    double precision :: x_Phase(30,5)               ! vector containing the compositions of all phases
    integer :: nrofphases
    !Indicate which phases are present
    integer :: phasetype(5)                          !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- =>  nrofphases = 2
    !-- =>  phasetype(1) = 2 (light liquid)
    !-- =>  phasetype(2) = 3 (heavy liquid)
    integer :: ILIMITS                           ! states whether the limits of the EOS have been hold, ONE VALUE BACK

    double precision :: tau, del, taumin, taumax, delmax, tmelt, pmelt
    integer :: i
    integer, dimension(3) :: check                            ! states whether each boundary has been hold
    !-- =>  check(1) : TMIN
    !-- =>  check(2) : TMAX
    !-- =>  check(3) : DMAX

    double precision :: rhoredmix_orig, tredmix_orig, rhomix_crit, rholim_ag8
    double precision, dimension(:), allocatable :: x, rho_crit
    check = 0
    ILIMITS = 0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! AGA8 limits
    if(.not.allocated(x)) then
        x = pack(gl%molfractions, gl%molfractions /= 0)         ! mole fraction vector without zeros
    end if

    if(.not.allocated(rho_crit)) then
        rho_crit = pack(gl%rhoc, gl%rhoc /= 0)                  ! critical density vector of pure fluids without zeros
    end if

    if (gl%Eq_type(1) == 7) then
        ! calculate critical density of mixture
        rhomix_crit = 1.d0/sum(x(1:size(x))/rho_crit(1:size(x)))
        rholim_ag8 = 0.5d0*rhomix_crit

        do i = 1, 5
            if ((nrofphases /= 1) .or. (phasetype(i) > 1)) then
                ILIMITS = -9994                                     ! only valid in gas phase
                return
            end if
        end do

        if (d >= rholim_ag8) then
            ILIMITS = -9997                                     !valid only up to 50% of critical density of mixture
            return
        end if

        ! AGA8-DC92 can be applied for temperatures between 143.15K and 453.15K and pressures up to 280MPa
        if ((t < 143.15) .or. (t > 453.15)) then
            ILIMITS = -9998
            return
        end if

        if (p > 280) then
            ILIMITS = -9999
            return
        end if
    end if


    ! Check whether all components are calculated with Helmholz
    ! If not then no limits have to be checked => exit the routine
    Do i = 1, gl%ncomp
        if (gl%Eq_type(i) /= 1) return
    end do

    !Andreas December 2013
    !Check limits deactivated for mixtures
    if (gl%ncomp > 1) return

    if ((index(input,"LIQ") /= 0) .or. (index(input,"VAP") /= 0)) return
    ! possible entries are errvals or values from trip to crit


    ! check singlephase results
    if(nrofphases == 1) then

        ! this procedure is the same for pure components and mixtures!
        tau = gl%tredmix / t
        del = d / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tminfluid(i)
            taumax = gl%tredmix / gl%tmaxfluid(i)
            delmax = gl%rhomaxfluid(i) / gl%rhoredmix
            if (tau > taumin) check(1) = 1
            if (tau < taumax) check(2) = 1
            if (del > delmax) check(3) = 1
            if (p > gl%pmaxfluid(i)) ILIMITS = -9932
        end do

        !pure fluid: if lower temperature limit is violated, check if melting curve is available
        ! and if the corresponding temperature smaller than the lower temperature limit
        ! that might be possible for water and heavy water
        ! Monika 2019/02/05
        do i = 1, gl%NCOMP
            if ((gl%ncomp == 1) .and. (check(1) == 1)) then
                tmelt = tmelt_eq(gl,p, i)
                if (tmelt > 1.d-5) then
                    taumin = gl%tredmix / tmelt
                    if (tau < taumin) check(1) = 0
                end if
            end if
        end do

        if (nrofphases == 6) then !p >= pmelt!
            if (gl%hold_limits) ILIMITS = -9915
        end if

        ! check twophase results
    else

        ! get reducing parameters for the vapor phase of the mixture
        if (gl%ncomp > 1) then
            gl%molfractions = x_phase(:,phasetype(1))
            call reduced_parameters_calc(gl,t)
        end if

        ! again the same for pure components and mixtures
        tau = gl%tredmix / t
        del = dvap / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tminfluid(i)
            taumax = gl%tredmix / gl%tmaxfluid(i)
            delmax = gl%rhomaxfluid(i) / gl%rhoredmix
            if (tau > taumin) check(1) = 1
            if (tau < taumax) check(2) = 1
            if (del > delmax) check(3) = 1
        end do

        ! get reducing parameters for the liquid phase of the mixture
        if (gl%ncomp > 1) then
            gl%molfractions = x_phase(:,phasetype(2))
            call reduced_parameters_calc(gl,t)
        end if

        tau = gl%tredmix / t
        del = dliq / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tminfluid(i)
            taumax = gl%tredmix / gl%tmaxfluid(i)
            delmax = gl%rhomaxfluid(i) / gl%rhoredmix
            if (tau > taumin) check(1) = 1   !lower temperature limit violated
            if (tau < taumax) check(2) = 1   !upper temperature limit violated
            if (del > delmax) check(3) = 1
        end do


        !Set the module variables back
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        !gl%molfractions = z
    end if


    !give warnings or throw errors, according to hold_limits
    do i = 1, 3
        if (check(i) == 1) then
            if (gl%hold_limits) ILIMITS = -9911 - i
        end if
    end do


    end subroutine check_limits


    module subroutine  set_limits_trans(gl)       ! will not be called with new version of interface
    !Monika June 2016









    implicit none

    type(type_gl) :: gl


    integer:: nrsubst

    if (gl%transport == 1) then !viscosity
        do nrsubst=1,gl%ncomp
            gl%tmintrans(nrsubst) = gl%visco%tmineta(nrsubst)          !lower temperature limit [K]
            gl%tmaxtrans(nrsubst) = gl%visco%tmaxeta(nrsubst)          !upper temperature limit [K]
            gl%pmaxtrans(nrsubst) = gl%visco%pmaxeta(nrsubst)          !upper pressure limit [MPa]
            gl%rhomaxtrans(nrsubst) = gl%visco%rhomaxeta(nrsubst)      !upper density limit [mol/m³]
        end do
    elseif (gl%transport == 2) then
        do nrsubst=1,gl%ncomp
            gl%tmintrans(nrsubst) = gl%tcx%tmin_tc(nrsubst)          !lower temperature limit [K]
            gl%tmaxtrans(nrsubst) = gl%tcx%tmax_tc(nrsubst)          !upper temperature limit [K]
            gl%pmaxtrans(nrsubst) = gl%tcx%pmax_tc(nrsubst)          !upper pressure limit [MPa]
            gl%rhomaxtrans(nrsubst) = gl%tcx%rhomax_tc(nrsubst)      !upper density limit [mol/m³]
        end do
    elseif (gl%transport == 3) then
        do nrsubst=1,gl%ncomp
            gl%tmintrans(nrsubst) = gl%de%tmin_de(nrsubst)          !lower temperature limit [K]
            gl%tmaxtrans(nrsubst) = gl%de%tmax_de(nrsubst)          !upper temperature limit [K]
            gl%pmaxtrans(nrsubst) = gl%de%pmax_de(nrsubst)          !upper pressure limit [MPa]
            gl%rhomaxtrans(nrsubst) = gl%de%rhomax_de(nrsubst)      !upper density limit [mol/m³]
            if (gl%rhomaxtrans(nrsubst) == 0.d0) gl%rhomaxtrans(nrsubst) = 1.d12   !dummy
        end do
    endif

    end subroutine set_limits_trans

    module subroutine check_limits_trans (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, prop, ILIMITS)
    !Monika June 2016


    implicit none

    type(type_gl) :: gl


    character(12) :: input
    double precision :: t, p, d, dliq, dvap
    double precision :: x_Phase(30,5)               ! vector containing the compositions of all phases
    integer :: nrofphases
    !Indicate which phases are present
    integer :: phasetype(5)                          !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- =>  nrofphases = 2
    !-- =>  phasetype(1) = 2 (light liquid)
    !-- =>  phasetype(2) = 3 (heavy liquid)
    integer :: ILIMITS                           ! states whether the limits of the EOS have been hold, ONE VALUE BACK

    double precision :: tau, del, taumin, taumax, delmax, pmax
    integer :: i, prop
    integer, dimension(4) :: check                            ! states whether each boundary has been hold
    !-- =>  check(1) : TMIN
    !-- =>  check(2) : TMAX
    !-- =>  check(3) : DMAX

    double precision :: rhoredmix_orig, tredmix_orig

    check = 0
    ILIMITS = 0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    !Check limits deactivated for mixtures
    if (gl%ncomp > 1) return

    if ((index(input,"LIQ") /= 0) .or. (index(input,"VAP") /= 0)) return
    ! possible entries are errvals or values from trip to crit

    ! check singlephase results
    if(nrofphases == 1) then

        ! this procedure is the same for pure components and mixtures!
        tau = gl%tredmix / t
        del = d / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tmintrans(i)
            taumax = gl%tredmix / gl%tmaxtrans(i)
            delmax = gl%rhomaxtrans(i) / gl%rhoredmix
            pmax = gl%visco%pmaxeta(i)
            if (tau > taumin) check(1) = 1
            if (tau < taumax) check(2) = 1
            if (del > delmax) check(3) = 1
            if (p   > pmax)   check(4) = 1
        end do

        if (nrofphases == 6) then !p >= pmelt!
            if (gl%hold_limits) ILIMITS = -9915
        end if

        ! check twophase results
    else

        ! get reducing parameters for the vapor phase of the mixture
        if (gl%ncomp > 1) then
            gl%molfractions = x_phase(:,phasetype(1))
            call reduced_parameters_calc(gl,t)
        end if

        ! again the same for pure components and mixtures
        tau = gl%tredmix / t
        del = dvap / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tmintrans(i)
            taumax = gl%tredmix / gl%tmaxtrans(i)
            delmax = gl%rhomaxtrans(i) / gl%rhoredmix
            pmax = gl%visco%pmaxeta(i)
            if (tau > taumin) check(1) = 1
            if (tau < taumax) check(2) = 1
            if (del > delmax) check(3) = 1
            if (p   > pmax)   check(4) = 1
        end do

        ! get reducing parameters for the liquid phase of the mixture
        if (gl%ncomp > 1) then
            gl%molfractions = x_phase(:,phasetype(2))
            call reduced_parameters_calc(gl,t)
        end if

        tau = gl%tredmix / t
        del = dliq / gl%rhoredmix

        do i = 1, gl%NCOMP
            taumin = gl%tredmix / gl%tmintrans(i)
            taumax = gl%tredmix / gl%tmaxtrans(i)
            delmax = gl%rhomaxtrans(i) / gl%rhoredmix
            pmax = gl%visco%pmaxeta(i)
            if (tau > taumin) check(1) = 1
            if (tau < taumax) check(2) = 1
            if (del > delmax) check(3) = 1
            if (p   > pmax)   check(4) = 1
        end do

        !Set the module variables back
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        !molfractions = z
    end if

    !give warnings or throw errors, according to hold_limits
    do i = 1, 4
        if (check(i) == 1) then
            if (gl%hold_limits) ILIMITS = -9917 - i
        end if
    end do

    end subroutine check_limits_trans




    !setup routine for one-fluid mixture model, Andreas Jäger, July 2016
    module subroutine setup_one_fluid_model(gl)






    implicit none

    type(type_gl) :: gl


    integer:: i, j, k, m
    logical:: newterm
    integer:: max_nreg,nr_of_equal_terms
    integer, dimension(maxval(gl%mix_nreg)):: equal_terms

    !Setup the "departure function" for new mixing rules (mix_type = 19)
    gl%mix_nreg = 0
    nr_of_equal_terms = 0
    newterm = .true.

    max_nreg=maxval(gl%mix_nreg)
    if(.not.allocated(gl%mix_nij)) allocate(gl%mix_nij(gl%ncomp,gl%ncomp,max_nreg,2))
    if(.not.allocated(gl%mix_nij)) then
        allocate(gl%mix_tij(gl%ncomp,gl%ncomp,max_nreg))
        allocate(gl%mix_p_ij(gl%ncomp,gl%ncomp,max_nreg))
        allocate(gl%mix_dij,gl%mix_gama,mold=gl%mix_tij)
        allocate(gl%mix_nreg(gl%ncomp,gl%ncomp))
    end if



    do i = 1, gl%ncomp
        do j = i+1, gl%ncomp
            equal_terms = 0
            nr_of_equal_terms = 0
            do k = 1, gl%eos_coeff%nreg(i)
                !Add the k-th regular term of fluid i to the mixture function
                gl%mix_nreg(i,j) = gl%mix_nreg(i,j) + 1
                gl%mix_tij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%ti(k,i)    !temp exponents
                gl%mix_dij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%di(k,i)           !density exponents
                gl%mix_p_ij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%p_i(k,i)         !density exponents in exponential term
                gl%mix_gama(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%gama(k,i)        !factor in exponential term, mostly 1
                gl%mix_nij(i,j,gl%mix_nreg(i,j),1) = gl%eos_coeff%ni(k,i)         !Coefficient of fluid i of term k
                gl%mix_nij(i,j,gl%mix_nreg(i,j),2) = 0.D0            !Coefficient of fluid j of term k (Set to 0 and maybe overwrite later)

                !Search for terms m of fluid j that are equal to the k-th regular term of fluid i
                do m = 1, gl%eos_coeff%nreg(j)

                    if ((dabs(gl%eos_coeff%ti(k,i) - gl%eos_coeff%ti(m,j)) < 1.D-10) .and. (dabs(gl%eos_coeff%di(k,i) - gl%eos_coeff%di(m,j)) < 1.D-10) .and. &
                        &  (gl%eos_coeff%p_i(k,i) == gl%eos_coeff%p_i(m,j)) .and. (dabs(gl%eos_coeff%gama(k,i) - gl%eos_coeff%gama(m,j)) < 1.D-10)) then

                        nr_of_equal_terms = nr_of_equal_terms + 1
                        equal_terms(nr_of_equal_terms) = m
                        gl%mix_nij(i,j,gl%mix_nreg(i,j),2) = gl%eos_coeff%ni(m,j)            !Coefficient of fluid j of term m

                    end if

                end do
            end do

            !Go through the terms of fluid j and add terms which are not the same as those of fluid i to the mixing function
            do k = 1, gl%eos_coeff%nreg(j)
                !Check if this term is already present in fluid i
                newterm = .true.
                do m = 1, nr_of_equal_terms
                    if (equal_terms(m) == k) then
                        newterm = .false.
                    end if
                end do
                if (newterm) then
                    gl%mix_nreg(i,j) = gl%mix_nreg(i,j) + 1
                    gl%mix_tij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%ti(k,j)            !temp exponents
                    gl%mix_dij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%di(k,j)            !density exponents
                    gl%mix_p_ij(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%p_i(k,j)          !density exponents in exponential term
                    gl%mix_gama(i,j,gl%mix_nreg(i,j)) = gl%eos_coeff%gama(k,j)         !factor in exponential term, mostly 1
                    gl%mix_nij(i,j,gl%mix_nreg(i,j),2) = gl%eos_coeff%ni(k,j)         !Coefficient of fluid j of term k
                    gl%mix_nij(i,j,gl%mix_nreg(i,j),1) = 0.D0            !Coefficient of fluid i of term k
                end if
            end do

        end do
    end do

    end subroutine

    !setup routine for for combination of gE model with multi fluid mixture model, Andreas, Erik, Dezember 2018
    module subroutine setup_HelmgE(gl, errval)






    implicit none

    type(type_gl) :: gl


    integer :: i
    double precision :: rhovap_est, rholiq_est   ![mol/m³]
    integer :: iFlash
    double precision :: press, Temp              ![MPa, K]
    integer :: errval, iter

    rhovap_est = 0.D0
    rholiq_est = 0.D0
    iFlash = 2
    !Dummy Temperature
    Temp = 0.D0

    !to calculate density of saturated liquid at 1 atm

    do i = 1, gl%ncomp
        press = 0.101325D0
        !test whether press is below triple point pressure of fluid
        if (gl%ptp(i) > press) then
            press = gl%ptp(i) + 0.0001D0
        end if

        call Flash_Pure_PhaseBoundary (gl, press, Temp, rhovap_est, rholiq_est, iFlash, errval, iter, i)
        if (errval /= 0) then
            return
        end if

        gl%rho_i_ref(i) = rholiq_est

        gl%bi_HelmgE(i) = 1.D0 / gl%rho_i_ref(i) / gl%u_pack

        rhovap_est = 0.D0
        rholiq_est = 0.D0
        Temp = 0.D0
    end do


    end subroutine

    ! Add comment
    subroutine calc_fld_specific_props(gl)
    implicit none
    type(type_gl) ::gl
    integer :: nrsubst,errorflag



    do nrsubst =1,gl%ncomp
        if( (gl%eq_type(nrsubst) .eq. 1) .or. (gl%eq_type(nrsubst) .eq. 51) .or. (gl%eq_type(nrsubst) .eq. 52) .or. (gl%eq_type(nrsubst) .eq. 53)) then
            gl%ref = .true.
            gl%rhoredmix = gl%rhored(nrsubst)
            gl%Tredmix = gl%Tred(nrsubst)
            gl%pc(nrsubst) = P_CALC(gl, gl%tc(nrsubst), gl%rhoc(nrsubst),nrsubst)  !critical pressure calculated from EOS [MPa]
            gl%ref = .false.
        elseif(gl%eq_type(nrsubst) .eq. 7) then
            gl%eq_type(nrsubst) = 2
            gl%rhored(nrsubst) = rho_SRK (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)
            gl%rhoc(nrsubst)  = gl%rhored(nrsubst)
            gl%eq_type(nrsubst) = 7
        end if
    end do



    end subroutine

    !**********************************************************************************************************************************
    !subroutine for specefic input handle
    !Benedikt April 2018

    !subroutine set_molarspecific(gl, convert,input, fluids, moles, prop1, prop2, x_spec)
    !

    !
    !implicit none
    !
    !    type(type_gl) :: gl
    !
    !
    !double precision :: mw_mix, molar_weight, wm_mix, sum_weight
    !double precision :: prop1, prop2
    !double precision, dimension (30) :: x_specific, x_molar, mw_comp, moles, x_spec
    !character(12) :: eos_plus_mix, converttype, input
    !character(12), dimension(30) :: fluids
    !integer :: i, convert                           ! convert=1 --> specific to molar, convert=2 --> molar to specific
    !convert=1
    !!gl%unitin = 1          !1=molar 2=specific
    !!gl%unitout = 1          !1=molar 2=specific
    !        !1=molar 2=specific
    !if ((convert == 1) .and. (gl%ncomp > 1)) then
    !    convert =1
    !elseif ((convert ==1) .and. (gl%ncomp==1)) then
    !    convert =12
    !else
    !    convert = 2
    !end if
    !
    !
    !if((convert == 1))then!.and. (gl%ncomp > 2)) then !.and. (gl%unitin /= gl%unitstat)) then
    !
    !    sum_weight = 0
    !
    !    do i=1,gl%ncomp                                                                     !this loop will calculate the the moles from specific input
    !        sum_weight = sum_weight + gl%molfractions(i) * (gl%wm(i)**(-1.d0))
    !    end do
    !
    !    do i=1,gl%ncomp
    !        x_molar(i) = (gl%molfractions(i) * (gl%wm(i)**(-1.d0)))  / (sum_weight)
    !        gl%molfractions(i) = x_molar(i)
    !        x_spec(i) = x_molar(i)                                  !for giving back the changed composition array
    !    end do
    !
    !        call wm_mix_calc(gl, wm_mix)
    !    ! mal schauen, ob man dass noch anders machen kann, damit der aufruf für die Molbruchumrechnung noch besser geht(mol_strtodbl)
    !    select case (trim(input))
    !        case('tp')
    !            prop1=prop1
    !            prop2=prop2
    !       case('td')
    !            prop1 = prop1
    !            prop2 = prop2 / wm_mix !* 1000
    !        case('ps')
    !            prop1 = prop1
    !            prop2 = prop2 * wm_mix!/ 1000
    !        case('ph')
    !            prop1 = prop1
    !            prop2 = prop2 * wm_mix!/ 1000
    !    end select
    !
    !
    !
    !        gl%unitstat = 1
    !
    !
    !elseif(convert==12) then
    !    call wm_mix_calc(gl, wm_mix)
    !
    !    select case (trim(input))
    !        case('tp')
    !            prop1=prop1
    !            prop2=prop2
    !       case('td')
    !            prop1 = prop1
    !            prop2 = prop2 / wm_mix !* 1000
    !        case('ps')
    !            prop1 = prop1
    !            prop2 = prop2 * wm_mix!/ 1000
    !        case('ph')
    !            prop1 = prop1
    !            prop2 = prop2 * wm_mix!/ 1000
    !    end select
    !end if
    !
    ! !end if
    !
    !!!elseif(conerttype == 2) then
    !!!
    !!!    do i=0,gl%ncomp
    !!!        !eos_plus_mix = TRIM(gl%eostype(i))//TRIM(";")//TRIM(gl%mixtype)
    !!!        x_spec(i) = moles(i) * gl%wm(i)
    !!!        !gl%molfractions(i) = x_mass(i)      !oder besser moles?!
    !!!
    !!!    end do
    !!!
    !!!    gl%unitstat = 2
    !!!
    !!!
    !!!
    !!!    !mw_mix = MW(fluids, moles, EOSType, path)                                           !converting d,h,s input to molar with mixture molecular mass mw_mix
    !!!    !gl%prop2 = gl%prop2 / mw_mix                                                        !be aware of cases: h,s is different from rho!!!!!
    !
    !!!elseif(convert == 2) then
    !!!
    !!!    mw_mix = MW(fluids, moles, gl%eostype, path)                                           !converting output to from molar to specific
    !!!
    !!!    if(gl%unitstat .neq. gl%unitin) then
    !!!    Select Case (gl%calctype_interanl)
    !!!
    !!!    Case(2)
    !!!    TREND_EOS = TREND_EOS * mw_mix
    !!!    end select
    !!!    end if
    !
    !!elseif(converttype .eq. "subout") then
    !
    !
    !!!else
    !!!
    !!!    errval = -123456        !error in specific in output
    !!!
    !
    !
    !!convertion :
    !
    !
    ! !!   Select Case Trim(input)
    ! !!   Case ("td")
    ! !!       !T[K]
    ! !!       !DENSITY [kg/m3] --> [mol/m3]
    ! !!       prop2 = prop2 / molar_weight / 1000.d0
    ! !!   Case ("dt")
    ! !!       !DENSITY [kg/m3] --> [mol/m3]
    ! !!       !T[K]
    ! !!       prop1 = prop1 / molar_weight / 1000.d0
    ! !!   Case ("ph")
    ! !!       ! P[MPa]
    ! !!       ! ENTHALPY [kJ/kg] --> [J/mol]
    ! !!       prop2 = prop2 * molar_weight * 1000
    ! !!   Case ("hp")
    ! !!       ! ENTHALPY [kJ/kg] --> [J/mol]
    ! !!       ! P[MPa]
    ! !!       prop1 = prop1 * molar_weight * 1000
    ! !!   Case ("ps")
    ! !!       ! P[MPa]
    ! !!       ! ENTROPY [kJ/kg/K] --> [J/mol/K]
    ! !!       prop2 = prop2 * molar_weight * 1000
    ! !!   Case ("sp")
    ! !!       ! ENTROPY [kJ/kg/K] --> [J/mol/K]
    ! !!       ! P[MPa]
    ! !!       prop1 = prop1 * molar_weight * 1000
    ! !!
    ! !!       end select
    !
    !end subroutine

    !*****************************************************************************************
    !subroutine wm_phase(gl, nrofphases, x_phase, wm_phase)
    !!calcutaes the molecular weight for each phase in flash or whatever routine for molar <--> specific handling
    !! Benedikt, April 2018
    !!
    !!

    !implicit none
    !
    !    type(type_gl) :: gl
    !
    !    integer :: i, j, k, nrofphases
    !    double precision, dimension(30,5) :: x_phase
    !    double precision, dimension(5), intent(out):: wm_phase
    !
    !    do i=1,nrofphases
    !        do j=1,gl%ncomp
    !        wm_phase(i) = x_phase(j,i) * gl%wm(j)
    !        end do
    !    end do
    !end subroutine wm_phase
    !*******************************************************************************************
    !
    !
    !!*******************************************************************************************
    !subroutine convert_spec_to_mol(gl, molesv, x_spec)
    !!converting specific input to molfractions for molar TREND routines
    !!Benedikt May 2018
    !!
    !!

    !implicit none
    !
    !type(type_gl) :: gl
    !
    !    integer :: i
    !    double precision :: sum_weight
    !    double precision, dimension (30) :: x_spec, molesv, x_molar
    !    double precision, dimension (30,5) :: prop_phase
    !
    !
    !    sum_weight = 0
    !
    !    do i=1,gl%ncomp                                                                     !this loop will calculate the the moles from specific input
    !        sum_weight = sum_weight + gl%molfractions(i) * (gl%wm(i)**(-1.d0))
    !    end do
    !
    !    do i=1,gl%ncomp
    !        x_molar(i) = (gl%molfractions(i) * (gl%wm(i)**(-1.d0)))  / (sum_weight)
    !        gl%molfractions(i) = x_molar(i)
    !        x_spec(i) = x_molar(i)                                  !for giving back the changed composition array
    !    end do
    !
    !end subroutine
    !*******************************************************************************************
    !
    !
    !
    !module subroutine convert_flash_spec(gl, prop_overall,wm_mix)
    !!converting flash properties for output handle
    !!Benedikt May 2018
    !!
    !!
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !integer :: i, j, nrofphases, IS_ERROR
    !double precision :: sum_weight
    !double precision :: wm_mix, wm_water
    !double precision, dimension (30), intent(inout) :: prop_overall
    !double precision, dimension (30,5) :: prop_phase, x_phase
    !character(255) :: errorcodes
    !
    !if(gl%seawater) then
    !    wm_water = gl%wm(1)
    !    gl%wm(1) = gl%sea%wm_sea
    !    call wm_mix_calc(gl, wm_mix)
    !elseif(gl%el_present) then
    !    wm_water = gl%wm(1)
    !    gl%wm(1) = wm_brine(gl)
    !    call wm_mix_calc(gl, wm_mix)
    !end if
    !
    !do i=3, 10
    !    !if (IS_ERROR(int(prop_overall(i))) == 1).and.((prop_overall(i) < -1000.d0) .and. (dabs(prop_overall(i) - int(prop_overall(i))) < 1.d-14))) then
    !    !    prop_overall(i) = prop_overall(i)
    !    if(i == 3) then
    !        prop_overall(i) = prop_overall(i) * wm_mix
    !    elseif (( i .gt. 3 ) .and. (i .lt. 11)) then
    !        prop_overall(i) = prop_overall(i) / wm_mix
    !    end if
    !end do
    !
    !if(gl%seawater .or. gl%el_present) then
    !    gl%wm(1) = wm_water
    !end if
    !
    !
    !end subroutine
    !!***************************************************************************************************
    !!
    !!
    !!***************************************************************************************************
    !module subroutine convert_fractions(gl, converttype, wm_phase, x_spec)!, wm_mix)
    !!converting phase composition form molar to specific or from specific to molar
    !!
    !!
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !integer :: i,j, converttype
    !double precision :: wm_mix, sum_weight, wm_phase      !Mixed molar weight for (phase) composition
    !! double precision, dimension (5) :: wm_phase
    !double precision, dimension (30) :: x_molar
    !double precision, dimension (30), intent(inout) :: x_spec
    !! double precision, dimension (30,5) :: x_phase
    !
    !!call wm_mix_calc(gl, wm_mix)
    !sum_weight = 0.d0
    !x_molar =0.d0
    !
    !
    !if(converttype == 1) then    !convert specific input to moles
    !    do i=1,gl%ncomp                                                                     !this loop will calculate the the moles from specific input
    !        sum_weight = sum_weight + gl%molfractions(i) * (gl%wm(i)**(-1.d0))
    !    end do
    !
    !    do i=1,gl%ncomp
    !        x_molar(i) = (gl%molfractions(i) * (gl%wm(i)**(-1.d0)))  / (sum_weight)
    !        gl%molfractions(i) = x_molar(i)
    !        x_spec(i) = x_molar(i)
    !    end do
    !
    !
    !elseif (converttype == 2) then  !convert calculated molar compositopns from TREND to specific output
    !    !call wm_mix_calc(gl, wm_mix)
    !    do i=1, gl%ncomp
    !        x_molar(i) = x_spec(i) * (gl%wm(i) / wm_phase)
    !    end do
    !    x_spec = x_molar
    !end if
    !
    !end subroutine
    !!**************************<**********************************************************************************
    !!
    !!
    !!
    !module subroutine convert_inputprop(gl, inputflag, prop)
    !!converting specific input to molar props for TREND calculations
    !!Benedikt May 2018
    !!
    !!
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !integer :: inputflag
    !double precision :: wm_mix
    !double precision, intent(inout) :: prop
    !
    !if(gl%already_converted == 1) then
    !    prop = prop
    !else
    !
    !    call wm_mix_calc(gl, wm_mix)
    !    gl%already_converted = 1
    !    if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
    !        prop = prop / wm_mix
    !    elseif(inputflag == 2) then     ! entropy, enthalpy,....
    !        prop = prop * wm_mix
    !    end if
    !end if
    !!futher routines for estimated densities in flash expert?!
    !
    !end subroutine
    !!**************************************************************************************************
    !!**************************************************************************************************
    !
    !module subroutine convert_single_prop(gl, inputflag, prop, nrsubst)
    !!converting single prop molar<--> specific for single component
    !!Benedikt June 2018
    !
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !integer :: inputflag, nrsubst
    !double precision :: wm, wm_sea
    !double precision, intent(inout) :: prop
    !
    !
    !if(nrsubst .eq. 9999)then
    !    !wm_sea = molar_sea(gl)
    !    wm_sea = gl%wm(1)
    !    !wm_sea = 0.0314038218d0
    !    ! wm_Sea = 1.d0
    !    if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
    !        prop = prop / wm_sea
    !    elseif(inputflag == 2) then     ! entropy, enthalpy,....
    !        prop = prop * wm_sea
    !    end if
    !else
    !    if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
    !        prop = prop / gl%wm(nrsubst)
    !    elseif(inputflag == 2) then     ! entropy, enthalpy,....
    !        prop = prop * gl%wm(nrsubst)
    !    end if
    !end if
    !end subroutine




    !This subroutine controls the allocation of the vectors and matrices in the sub-types
    module subroutine allocate_subtypes(gl, errorflag)



    implicit none

    type(type_gl), intent(inout) :: gl

    integer:: size1 !variable for first dimension allocation
    integer:: size2 !variable for second dimension allocation
    integer:: errorflag

    errorflag = 0

    if ((.not. allocated(gl%ge)).and.((gl%mix_type == 12).or.(gl%mix_type == 13).or.(gl%mix_type == 22))) then
        allocate(gl%ge)

        ! intalize ge variables / moved from init
        gl%ge%Temp_prev = 0.D0
        gl%ge%GETDER_prev = 0
        !Dens_prev = 0.d0
        gl%ge%molfrac_prev = 0.D0
        gl%ge%Eq_type_prev = 0
        gl%ge%mixtype_prev = 0
        gl%ge%C_or_R_prev = -1
        gl%ge%gE_C_prev = 0.D0
        gl%ge%gE_R_prev = 0.D0
        gl%ge%ln_gamma_C_prev = 0.D0
        gl%ge%ln_gamma_R_prev = 0.D0

        gl%ge%Temp_dxa_prev = 0.D0
        gl%ge%GETDER_dxa_prev = 0
        !Dens_dxa_prev = 0.d0
        gl%ge%molfrac_dxa_prev = 0.D0
        gl%ge%Eq_type_dxa_prev = 0
        gl%ge%mixtype_dxa_prev = 0
        gl%ge%C_or_R_dxa_prev = -1
        gl%ge%gE_C_dxa_prev = 0.D0
        gl%ge%gE_R_dxa_prev = 0.D0
        gl%ge%ln_gamma_C_dxa_prev = 0.D0
        gl%ge%ln_gamma_R_dxa_prev = 0.D0

        gl%ge%Temp_dxadxb_prev = 0.D0
        gl%ge%GETDER_dxadxb_prev = 0
        !Dens_dxadxb_prev = 0.d0
        gl%ge%molfrac_dxadxb_prev = 0.D0
        gl%ge%Eq_type_dxadxb_prev = 0
        gl%ge%mixtype_dxadxb_prev = 0
        gl%ge%C_or_R_dxadxb_prev = -1
        gl%ge%gE_C_dxadxb_prev = 0.D0
        gl%ge%gE_R_dxadxb_prev = 0.D0


        gl%ge%ln_gamma_C_dxadxb_prev = 0.D0
        gl%ge%ln_gamma_R_dxadxb_prev = 0.D0
    end if

    if ((.not. allocated(gl%cosmo)) .and. (gl%mix_type == 13)) allocate(gl%cosmo)
    if (.not. allocated(gl%visco) .and. ((gl%transport == -1) .or. (gl%transport == 1))) then
        allocate (gl%visco)
        call initialize_viscosity(gl)
        call initialize_transport_properties(gl)
        if (gl%same_components) call read_from_fluidfiles (gl,gl%path_arg, errorflag)
    end if

    if (.not. allocated(gl%tcx) .and. ((gl%transport == -1) .or. (gl%transport == 2))) then
        allocate (gl%tcx)
        if (.not. allocated(gl%visco)) then
            allocate (gl%visco)
            call initialize_viscosity(gl)
        end if
        if (.not. allocated(gl%tcx%accen_TC5)) then
            allocate(gl%tcx%accen_TC5(gl%ncomp),gl%tcx%tc_TC5(gl%ncomp),gl%tcx%dipolered_TC5(gl%ncomp),gl%tcx%kappa_TC5(gl%ncomp),gl%tcx%addchung_TC5(gl%ncomp))
        endif
        call initialize_thermal_conductivity(gl)
        call initialize_transport_properties(gl)
        if (gl%same_components) call read_from_fluidfiles (gl,gl%path_arg, errorflag)
    end if

    if (.not. allocated(gl%de) .and. ((gl%transport == -1) .or. (gl%transport == 3))) then
        allocate (gl%de)
        call initialize_transport_properties(gl)
        if (gl%same_components) call read_from_fluidfiles (gl,gl%path_arg, errorflag)
    end if

    if (.not. allocated(gl%stn)) then
        allocate (gl%stn(gl%ncomp))
        if (gl%same_components) call read_from_fluidfiles (gl,gl%path_arg, errorflag)
    end if

    if ((.not. allocated(gl%eos_coeff)) .and. (any((gl%Eq_Type(:)==1) .or. (gl%Eq_Type(:)==7) .or. (gl%Eq_Type(:)==51) .or. (gl%Eq_Type(:)==52) .or. (gl%Eq_Type(:)==53)))) allocate(gl%eos_coeff)


    !check if a sub-type is allocated:
    ! IF Yes: allocate the members which are allocatable to the necessary size
    ! IF  NO: skip the sub type


    size1 = 100
    size2 = gl%ncomp !second dimension
    !Sub type for the eos-coefficients (EOS-TYPE = 1) (aga8 eqtype = 7)
    if(allocated(gl%eos_coeff) .and. any((gl%Eq_type(1:size2)==1) .or. (gl%Eq_type(1:size2)==7) .or. (gl%Eq_type(1:size2)==51) .or. (gl%Eq_type(1:size2)==52) .or. (gl%Eq_type(1:size2)==53))) then

        !Info section variables allocation: ################################################################
        if (.not.allocated(gl%substshortname)) then
            allocate(gl%substshortname(gl%ncomp))
            allocate(gl%substfullname(gl%ncomp))
            allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
            allocate(gl%substcasnr,mold=gl%substshortname)
        else
            if (size(gl%substshortname,1) /= gl%ncomp) then
                deallocate(gl%substshortname,gl%substfullname,gl%substchemname,gl%substsynonym,gl%substcasnr)
                allocate(gl%substshortname(gl%ncomp))
                allocate(gl%substfullname(gl%ncomp))
                allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
                allocate(gl%substcasnr,mold=gl%substshortname)
            endif
        end if
        ! end info section allocation #######################################################################

        !size1 =  maxval(gl%eos_coeff%nreg(1:size2))+maxval(gl%eos_coeff%ncrt(1:size2))+maxval(gl%eos_coeff%nna(1:size2))

        if(.not.allocated(gl%eos_coeff%ni)) then
            size1 = 100
            allocate(gl%eos_coeff%ni(size1,size2))     !coefficients
            allocate(gl%eos_coeff%ti(size1,size2))      !temp exponents
            allocate(gl%eos_coeff%di(size1,size2))      !density exponents
            size1 = 70
            allocate(gl%eos_coeff%p_i(size1,size2))      !density exponents in exponential term
            allocate(gl%eos_coeff%gama(size1,size2))     !factor in exponential term, mostly 1




        endif
        !critical + nonanalytic (why always both even if one does not exist? - see read from fluid files)
        if(.not.allocated(gl%eos_coeff%eps)) then
            !#######################################################
            !uncomment this part when the revised hardcode generation routines are finished
            !if(any(gl%eos_coeff%ncrt(1:size2)/=0).or.any(gl%eos_coeff%nna(1:size2)/=0)) then
            size1 = 20
            allocate(gl%eos_coeff%eps(size1,size2))      !epsilon in gaussian term
            allocate(gl%eos_coeff%beta(size1,size2))     !beta in gaussian term
            allocate(gl%eos_coeff%gam(size1,size2))      !gamma in gaussian term
            allocate(gl%eos_coeff%eta(size1,size2))      !eta in gaussian term
            allocate(gl%eos_coeff%tli(size1,size2))      !exponent always 2 (not needed?)
            allocate(gl%eos_coeff%pli(size1,size2))      !exponent always 2 (not needed?)

            allocate(gl%eos_coeff%etana(size1,size2))    !nonanalyticterm 1
            allocate(gl%eos_coeff%eidna(size1,size2))    !nonanalyticterm 2
            allocate(gl%eos_coeff%eitna(size1,size2))    !nonanalyticterm 3
            !end if
            !#######################################################
        end if

        !hardsphere terms
        !allocated every time because info about hard sphere in subfile-input
        if (.not.allocated(gl%eos_coeff%h1)) then
            ! if(any(gl%eos_coeff%hard_sphere(1:size2))) then
            allocate(gl%eos_coeff%h1(size2))
            allocate(gl%eos_coeff%h2(size2))
            allocate(gl%eos_coeff%h3(size2))
            allocate(gl%eos_coeff%h4(size2))
            allocate(gl%eos_coeff%dof2(size2))
            ! end if
        end if

        !Sub type for eos-coefficietns (EOS_TYPE = 51) !calculation based on Alexandrov
        !elseif(allocated(gl%eos_coeff).and.((gl%Eq_type(1)==51).or.(gl%Eq_type(1)==52).or.(gl%Eq_type(1)==53))) then


        !if(gl%Eq_type(1)==51) then
        !    size1 =  gl%anz_term_igor
        !elseif(gl%Eq_type(1)==52) then
        !    size1 = gl%anz_term_span
        !elseif(gl%Eq_type(1)==53) then
        !    size1 = gl%anz_term_sun
        !end if
        !size1 = 100
        !if(.not.allocated(gl%eos_coeff%ni))   allocate(gl%eos_coeff%ni(size1,size2))     !coefficients
        !if(.not.allocated(gl%eos_coeff%ti))   allocate(gl%eos_coeff%ti(size1,size2))      !temp exponents
        !if(.not.allocated(gl%eos_coeff%di))   allocate(gl%eos_coeff%di(size1,size2))      !density exponents
        !if(.not.allocated(gl%eos_coeff%p_i))  allocate(gl%eos_coeff%p_i(size1,size2))      !density exponents in exponential term
        !if(.not.allocated(gl%eos_coeff%gama)) allocate(gl%eos_coeff%gama(size1,size2))     !factor in exponential term, mostly 1

    elseif (allocated(gl%ge) .and.((gl%mix_type == 12).or.(gl%mix_type == 13).or.(gl%mix_type == 22))) then
        if(.not.allocated(gl%ge%ln_gamma_C_dxa))       allocate(gl%ge%ln_gamma_C_dxa(nderivs,30,30))
        if(.not.allocated(gl%ge%ln_gamma_R_dxa))       allocate(gl%ge%ln_gamma_R_dxa(nderivs,30,30))
        if(.not.allocated(gl%ge%gE_C_dxadxb))          allocate(gl%ge%gE_C_dxadxb(nderivs,30,30))
        if(.not.allocated(gl%ge%gE_R_dxadxb))          allocate(gl%ge%gE_R_dxadxb(nderivs,30,30))
        if(.not.allocated(gl%ge%ln_gamma_C_dxadxb))    allocate(gl%ge%ln_gamma_C_dxadxb(nderivs,30,30,30))
        if(.not.allocated(gl%ge%ln_gamma_R_dxadxb))    allocate(gl%ge%ln_gamma_R_dxadxb(nderivs,30,30,30))
    end if
    size1 = 0
    if (gl%check_solid) then
        if(.not.allocated(gl%hydrate_list)) then
            allocate(gl%hydrate_list(9))
            size1 = size1 + 1
        endif
        if(.not.allocated(gl%solid_list)) then
            allocate(gl%solid_list(2))
            size1 = size1 + 1
        endif
        if (size1 > 0) call initialize_solids(gl)

    endif

    ! ###################################################################################################
    ! Check if one of the eq_types is  a cubic model
    ! allocation of cubic models parameters
    if(any(gl%Eq_type(:) == 2) .or. any(gl%Eq_type(:) == 3) .or. any(gl%Eq_type(:) == 4)) then
        if(.not.allocated(gl%kij_SRK)) then
            allocate(gl%kij_SRK(gl%ncomp,gl%ncomp))
            allocate(gl%aij_SRK,gl%kij_PR,gl%aij_PR,gl%lij_SRK,gl%bij_SRK,gl%lij_PR,gl%bij_PR,mold=gl%kij_SRK)
            gl%aij_SRK = 0d0
            gl%bij_SRK = 0d0
            gl%kij_SRK = 0d0
            gl%lij_SRK = 0d0
            gl%kij_PR = 0d0
            gl%aij_PR = 0d0
            gl%lij_PR = 0d0
            gl%bij_PR = 0d0
        end if
    end if
    ! ###################################################################################################

    ! check if the mixtype fits one of the models (1,7,110,111,120,121)
    if(gl%Mix_type == 1 .or. gl%Mix_type == 7 .or. gl%Mix_type == 110 .or. gl%Mix_type == 111 .or. gl%Mix_type == 120 .or. gl%Mix_type == 121) then
        !allocation of variables used in mixtures_i
        if(.not.allocated(gl%dfn)) then
            allocate(gl%dfn(25,gl%ncomp,gl%ncomp))
            allocate(gl%dfd,gl%dft,gl%dfl,gl%dfp,gl%dfeta,gl%dfeps,gl%dfbeta,gl%dfgamma,gl%dfgeta,gl%dfgbeta,gl%dfggam,gl%dfgeps,mold=gl%dfn)
            gl%dfn     = 0d0
            gl%dfd     = 0d0
            gl%dft     = 0d0
            gl%dfl     = 0d0
            gl%dfp     = 0d0
            gl%dfeta   = 0d0
            gl%dfeps   = 0d0
            gl%dfbeta  = 0d0
            gl%dfgamma = 0d0
            gl%dfgeta  = 0d0
            gl%dfgbeta = 0d0
            gl%dfggam  = 0d0
            gl%dfgeps  = 0d0

            allocate(gl%dfd_coeff_structure(gl%ncomp,gl%ncomp))
            allocate(gl%dfl_coeff_structure(gl%ncomp,gl%ncomp))
            gl%dfd_coeff_structure = coeff_structure_all_zeros
            gl%dfl_coeff_structure = coeff_structure_all_zeros
        end if

    end if

    end subroutine allocate_subtypes



    end submodule impl
