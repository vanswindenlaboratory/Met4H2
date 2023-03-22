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


    !DEC$ IF DEFINED(HDRT)
    submodule (hdrt_data) impl

    use interface_support_module
    use module_all_types
    use module_hdrt_phaselines
    use ptflash_2c_3p_module
    use ptflash_solids_mix_module
    use setup_module
    use hdrt_main_module
    use hdrt_properties_modules_module
    use hdrt_chem_pot_module
    use hdrt_phaselines_calc


    !type hdrt_experimental_data
    !    character(255) :: author
    !    character(1000) :: readr_in
    !    double precision :: weightIn
    !    double precision :: TempIn, pressIn, VcellIn, VliqIn, Tempfill, pressfill
    !    double precision, dimension(30) :: xIn, x_vap_givenIn, x_fluid1In, x_fluid2In, x_fluid3In, x_fluid4In, x_solIn, x_hydIn, n_vapIn
    !    double precision :: x_h2oIn
    !
    !    integer :: ScenarioIn, iPhase_sol
    !    character(10) :: x_vap_givenIn_Equil
    !    character(3) :: hdrt_structure_exp
    !    character(20) :: EquiltypeIn, EquiltypeGen        !Given Type of Equilibrium, which should be calculated
    !
    !    double precision :: tempcalc, presscalc
    !    double precision, dimension(2) :: ttest
    !    !variables for properties of different phases
    !    !##############################################
    !    double precision, dimension(30) :: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, n_vap
    !    double precision, dimension(30,5) :: x_Phase
    !    double precision, dimension(4) :: rho_ph
    !    double precision, dimension(5) :: rho
    !    double precision, dimension(30) :: Chempot, lnfi, hyd_nrIn, hyd_nr
    !
    !    double precision, dimension(3):: rho_sol
    !
    !    ! Mixed hydrates 2015
    !    double precision, dimension(3,30):: occup, occup_single, occup_double, C_iJ
    !    double precision, dimension(3,30):: occupIn
    !    double precision, dimension(:), allocatable :: Phasefrac
    !    !##############################################
    !end type hdrt_experimental_data

    contains

    module subroutine hdrt_init_expdata(expdata)
    implicit none
    type(hdrt_experimental_data) :: expdata
    expdata%author = ""
    expdata%x_vap_givenIn_Equil = ""
    expdata%hdrt_structure_exp = ""
    expdata%EquiltypeIn = ""
    expdata%EquiltypeGen = ""
    expdata%ScenarioIn = 0
    expdata%iPhase_sol = 0

    expdata%weightIn        = 0.d0
    expdata%TempIn          = 0.d0
    expdata%pressIn         = 0.d0
    expdata%VcellIn         = 0.d0
    expdata%VliqIn          = 0.d0
    expdata%Tempfill        = 0.d0
    expdata%pressfill       = 0.d0
    expdata%xIn             = 0.d0
    expdata%x_vap_givenIn   = 0.d0
    expdata%x_h2oIn         = 0.d0
    expdata%tempcalc        = 0.d0
    expdata%presscalc       = 0.d0
    expdata%ttest           = 0.d0
    expdata%x_fluid1        = 0.d0
    expdata%x_fluid2        = 0.d0
    expdata%x_fluid3        = 0.d0
    expdata%x_fluid4        = 0.d0
    expdata%x_sol           = 0.d0
    expdata%x_hyd           = 0.d0
    expdata%n_vap           = 0.d0
    expdata%x_fluid1In      = 0.d0
    expdata%x_fluid2In      = 0.d0
    expdata%x_fluid3In      = 0.d0
    expdata%x_fluid4In      = 0.d0
    expdata%x_solIn         = 0.d0
    expdata%x_hydIn         = 0.d0
    expdata%n_vapIn         = 0.d0
    expdata%x_Phase         = 0.d0
    expdata%rho_ph          = 0.d0
    expdata%rho             = 0.d0
    expdata%Chempot         = 0.d0
    expdata%lnfi            = 0.d0
    expdata%hyd_nr = 0.d0
    expdata%rho_sol         = 0.d0
    expdata%occupIn         = 0.d0
    expdata%occup           = 0.d0
    expdata%occup_single    = 0.d0
    expdata%occup_double    = 0.d0
    expdata%C_ij            = 0.d0
    !    expdata%Phasefrac       = 0.d0
    end subroutine hdrt_init_expdata


    !******************************************************************************
    module subroutine hydratedatadev(gl,pel,DatenfilepathIn, pathIn, fluidspath, error)!, Q_point, Qexist, Temp_Q, press_Q, &
    !& x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    !******************************************************************************


    implicit none

    type(type_gl), dimension(:) :: gl
    type(type_hdrt_pheq) :: pel
    type(hdrt_experimental_data), dimension(:), allocatable :: expdata

    character(255), intent(inout) :: pathIn, fluidspath
    character(255), intent(in) :: DatenfilepathIn
    integer, intent(inout) :: error

    integer :: Datenfileunit
    integer :: Resultfileunit, Resulttwophasefileunit, Resulthydcompfileunit


    character(12) :: input
    character(255) :: devfilepath, twophasefilepath, hydratecompositionpath
    character(255) :: dummy1, dummy2
    character(30), dimension(30) :: components_orig

    !Variables at/for the three-phase equilibrium
    double precision :: deviation, deviation_abs_cumm, deviation_cumm, rho_feed, n_water
    double precision,dimension(4):: rhoph_est

    integer, dimension(30) :: EOS_Indicator
    integer :: MIX_indicator, bisect_mode
    integer :: errorin, errorout, i, ii, i_old, j, pt, pts_err, iTp, column, skip_data
    integer :: iFlash2p, iFlash3p, iFlash4p, iter, Phasefrac_0
    integer :: pointsexp_total, ncomp_orig
    integer, dimension(2) :: solidtype_copy
    integer :: hdrt_structure_flag, gl_indicator

    logical :: show_once = .false., next_author

    !Variables at the Quadruple point
    character(len=6), dimension(4)  :: Q_point
    double precision, dimension (4) :: press_Q, Temp_Q
    double precision, dimension(:,:), allocatable :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist

    double precision, dimension(30) :: composition, chem_pot, fug_g


    !//////////////////////////////////////////////////////////////////

    devfilepath = trim(pathIn) // trim(pel%CompsHdrtsAll) // '\three_phase_data\deviation.csv'
    twophasefilepath = trim(pathIn) // trim(pel%CompsHdrtsAll) // '\three_phase_data\deviation-two-phase.csv'
    hydratecompositionpath = trim(pathIn) // trim(pel%CompsHdrtsAll) // '\three_phase_data\datadev-' // trim(pel%CompsHdrtsAll) // '-hydrate-composition.csv'
    if (.not.allocated(x_Q_ph1))allocate(x_Q_ph1(4,30))
    if (.not.allocated(x_Q_ph2))allocate(x_Q_ph2(4,30))
    if (.not.allocated(x_Q_ph3))allocate(x_Q_ph3(4,30))
    if (.not.allocated(x_Q_ph4))allocate(x_Q_ph4(4,30))
    error = 0
    pel%phl3%eqtypestart = 1
    do i = 1,size(pel%phl4)
        pel%phl4(i)%Q_map = 1
    enddo
    !components_orig = gl%components
    !ncomp_orig = gl%ncomp
    input = trim('tp+')
    next_author = .false.

    iTp = 2         ! p is given
    !Setup for the three phase flash - T or p solved
    if (iTp == 1) then
        iFlash2p = 3        ! p, T, x given x' is determined
        iFlash3p = 4        ! Calculate Freezing point, T given
        iFlash4p = 2        ! Calculate Freezing point, T given
    elseif (iTp == 2) then
        iFlash2p = 3        ! p, T, x given x' is determined
        iFlash3p = 2        ! Calculate Freezing point, p given
        iFlash4p = 1        ! Calculate Freezing point, T given
    end if

    call hydratedatadev_readfile(gl, pel, expdata, Datenfilepathin, fluidspath)
    open(newunit = Resultfileunit, file=trim(devfilepath), status='unknown', action='write', iostat=errorout)
    if (errorout == 0) open(newunit = Resulttwophasefileunit, file=trim(twophasefilepath), status='unknown', action='write', iostat=errorout)
    if (errorout == 0) open(newunit = Resulthydcompfileunit, file=trim(hydratecompositionpath), status='unknown', action='write', iostat=errorout)
    if (errorout /= 0) then
        call messages(gl(1), pel,201,trim(''),trim(''),0,0,0)
        return
    endif

    pts_err = 0
    deviation_abs_cumm = 0.d0
    deviation_cumm = 0.d0
    if (pel%print_head .eqv..true.) then
        ii = len(trim('author'))
        j = len(trim('EqType'))
        write(Resultfileunit,3022)'author','EqType','pressExp', 'PropExp', 'PropCalc', 'Exp - Calc'
        ii = len(trim('author'))
        j = len(trim('EqType'))
        write(Resulttwophasefileunit,3022)'author','EqType','pressExp', 'TempExp', 'TempCalc', 'Exp - Calc'
        ii = len(trim('author'))
        j = len(trim('EqType'))
        write(Resulthydcompfileunit,3021)'author','EqType','pressExp', 'TempExp', 'nhexp', 'nhcalc','occupexpS','occupcalcS','occupexpL','occupcalcL'
    endif
    !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
    write(*,*)
    !DEC$ END IF ! WO_WRITE_TO_CONSOLE
    deviationloop: do pt = 1,size(expdata)
        gl_indicator = 1    !if N_hdrts = 2 -> pure hydrates always
        !                   2: pure hydrates with guest 1
        !                   n: pure hydrates with guest n-1
        call expdata_fill(gl, gl_indicator, pel, expdata, pt, fluidspath)
        if ((expdata(pt)%readr_in(1:1) == "!").or.((expdata(pt)%weightin == 0.d0).and.(expdata(pt)%pressin == 0.d0).and.(expdata(pt)%tempin == 0.d0))) cycle

        !force calculation of specified hydrate structure
        call hdrt_ref_a_gw0_hw0_COC(gl(gl_indicator))
        if ((trim(expdata(pt)%EquiltypeGen) /= 'VLxLw').and.(trim(expdata(pt)%EquiltypeGen) /= 'VLwHIw').and.(trim(expdata(pt)%EquiltypeGen) /= 'VLwLxH')&
            & .and.(trim(expdata(pt)%EquiltypeGen) /= 'VH').and.(trim(expdata(pt)%EquiltypeGen) /= 'LxH').and.(trim(expdata(pt)%EquiltypeGen) /= 'LwH')) then
            if (.not.(allocated(expdata(pt)%phasefrac)))allocate(expdata(pt)%phasefrac(3))
            if (expdata(pt)%HDRT_STRUCTURE_EXP == 'non') then
                do i = 1,2
                    call hdrt_structure_definition(gl(gl_indicator),i)
                    expdata(pt)%ttest(i) = expdata(pt)%tempcalc
                    call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%ttest(i), expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                        & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                enddo
                if (expdata(pt)%ttest(1) > expdata(pt)%ttest(2)) then
                    expdata(pt)%HDRT_STRUCTURE_EXP = 's1'
                    i = 1
                else
                    expdata(pt)%HDRT_STRUCTURE_EXP = 's2'
                    i = 2
                endif
                call hdrt_structure_definition(gl(gl_indicator),i)
                call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                    & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
            else
                call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                    & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
            endif
            if (gl(gl_indicator)%n_guests == 2) then
                if (((expdata(pt)%ScenarioIn == 1).and.(gl_indicator == 1)).and.(error == 0)) then !only for mixed hydrates
                    call threephase_equil_iter(gl(gl_indicator), pel, expdata(pt)%EquiltypeIn, 1, iTp, expdata(pt)%tempcalc, expdata(pt)%pressIn, expdata(pt)%xIn, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_hyd, expdata(pt)%rho_ph(1), expdata(pt)%rho_ph(2), expdata(pt)%rho_ph(3),  &
                        & expdata(pt)%Chempot, expdata(pt)%lnfi, expdata(pt)%occup, error, expdata(pt)%x_vap_givenIn, 0)!, Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4, expdata(pt)%x_vap_givenIn)
                    !MB Maximilian This If-Statement gives the oppurtunity for a second calculation in the routine "threephase_equil_iter" with another mixrule of x_upper and x_lower (detection with a errorcode) and determination with "bisect_mode" flag
                    if (error /= 0) then
                        do bisect_mode = 1, 2
                            call hydratedatadev_processdata(gl, expdata, pt, pt)
                            call expdata_fill(gl, gl_indicator, pel, expdata, pt, fluidspath)
                            !force calculation of specified hydrate structure
                            call hdrt_ref_a_gw0_hw0_COC(gl(gl_indicator))
                            if (trim(expdata(pt)%EquiltypeGen) /= 'VLxLw') then
                                call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                                    & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                                if (gl(gl_indicator)%n_guests == 2) then
                                    if (((expdata(pt)%ScenarioIn == 1).and.(gl_indicator == 1)).and.(error == 0)) then !only for mixed hydrates
                                        call threephase_equil_iter(gl(gl_indicator), pel, expdata(pt)%EquiltypeIn, 1, iTp, expdata(pt)%tempcalc, expdata(pt)%pressIn, expdata(pt)%xIn, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_hyd, expdata(pt)%rho_ph(1), expdata(pt)%rho_ph(2), expdata(pt)%rho_ph(3),  &
                                            & expdata(pt)%Chempot, expdata(pt)%lnfi, expdata(pt)%occup, error, expdata(pt)%x_vap_givenIn, bisect_mode)!, Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4, expdata(pt)%x_vap_givenIn)
                                    endif
                                endif
                            endif
                            if (error == 0) then
                                if (expdata(pt)%HDRT_STRUCTURE_EXP == 'non') then
                                    do i = 1,2
                                        call hdrt_structure_definition(gl(gl_indicator),i)
                                        expdata(pt)%ttest(i) = expdata(pt)%tempcalc
                                        call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%ttest(i), expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                                            & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                                    enddo
                                    if (expdata(pt)%ttest(1) > expdata(pt)%ttest(2)) then
                                        expdata(pt)%HDRT_STRUCTURE_EXP = 's1'
                                        i = 1
                                    else
                                        expdata(pt)%HDRT_STRUCTURE_EXP = 's2'
                                        i = 2
                                    endif
                                    call hdrt_structure_definition(gl(gl_indicator),i)
                                    call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                                        & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                                else
                                    call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                                        & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                                endif
                                exit
                            endif
                        enddo
                    else
                        continue
                    endif
                endif
            endif
        elseif ((trim(expdata(pt)%EquiltypeGen) == 'VLwHIw').or.(trim(expdata(pt)%EquiltypeGen) == 'VLwLxH')) then
            Phasefrac_0 = 3
            if (.not.allocated(expdata(pt)%phasefrac)) allocate(expdata(pt)%phasefrac(4))
            if (expdata(pt)%HDRT_STRUCTURE_EXP == 'non') then
                do i = 1,2
                    call hdrt_structure_definition(gl(gl_indicator),i)
                    expdata(pt)%ttest(i) = expdata(pt)%tempcalc
                    call ptflash_solid_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%ttest(i), expdata(pt)%xIn, expdata(pt)%rho_sol, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                        & rhoph_est(2), expdata(pt)%phasefrac, iFlash3p, expdata(pt)%iPhase_sol, iter, error)
                enddo
                if (expdata(pt)%ttest(1) > expdata(pt)%ttest(2)) then
                    expdata(pt)%HDRT_STRUCTURE_EXP = 's1'
                    i = 1
                else
                    expdata(pt)%HDRT_STRUCTURE_EXP = 's2'
                    i = 2
                endif
                call hdrt_structure_definition(gl(gl_indicator),i)
                call ptflash_solid_NC_4P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_fluid3, expdata(pt)%x_fluid4, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                    & rhoph_est(2), rhoph_est(3), rhoph_est(4), expdata(pt)%phasefrac, iFlash4p, expdata(pt)%iPhase_sol, Phasefrac_0, error, iter)
            else
                call ptflash_solid_NC_4P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_fluid3, expdata(pt)%x_fluid4, expdata(pt)%x_sol, expdata(pt)%x_hyd, rhoph_est(1), &
                    & rhoph_est(2), rhoph_est(3), rhoph_est(4), expdata(pt)%phasefrac, iFlash4p, expdata(pt)%iPhase_sol, Phasefrac_0, error, iter)
            endif
        elseif ((trim(expdata(pt)%EquiltypeGen) == 'VH').or.(trim(expdata(pt)%EquiltypeGen) == 'LxH').or.(trim(expdata(pt)%EquiltypeGen) == 'LwH')) then
            if (.not.allocated(expdata(pt)%phasefrac)) allocate(expdata(pt)%phasefrac(2))
            call ptflash_solid_NC_2P (gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempIn, expdata(pt)%rho_sol, expdata(pt)%xIn, expdata(pt)%x_sol, expdata(pt)%x_hyd, expdata(pt)%x_fluid1, rhoph_est(1), expdata(pt)%phasefrac(1), iFlash2p, expdata(pt)%iPhase_sol, error, iter)
            expdata(pt)%phasefrac(2) = 1.d0 - expdata(pt)%phasefrac(1)
        else
            if (.not.allocated(expdata(pt)%phasefrac)) allocate(expdata(pt)%phasefrac(3))
            call ptflash_NC_3P(gl(gl_indicator),expdata(pt)%pressIn, expdata(pt)%tempcalc, expdata(pt)%xIn, expdata(pt)%rho, expdata(pt)%x_fluid1, expdata(pt)%x_fluid2, expdata(pt)%x_fluid3, expdata(pt)%rho_ph(1), expdata(pt)%rho_ph(2), expdata(pt)%rho_ph(3), expdata(pt)%phasefrac, iFlash3p, iter, error)
        endif
        if ((error == 0) .and. (expdata(pt)%EquiltypeIn(1:1) == 'h')) then
            call hdrt_ancillary_hs(gl(gl_indicator),expdata(pt)%tempcalc, expdata(pt)%rho_sol(1), expdata(pt)%pressIn *1.d6, expdata(pt)%x_fluid1, chem_pot, fug_g, error)
            call hdrt_mole_fract(gl(gl_indicator),expdata(pt)%tempcalc,expdata(pt)%pressIn *1.d6,fug_g,expdata(pt)%occup,expdata(pt)%C_iJ,expdata(pt)%x_hyd,expdata(pt)%occup_single, expdata(pt)%occup_double)
            do J = 1,gl(gl_indicator)%N_guests
                do ii = 1,gl(gl_indicator)%N_cavi
                    expdata(pt)%hyd_nr(J) = expdata(pt)%hyd_nr(J) + gl(gl_indicator)%v_cavi(ii) * expdata(pt)%occup(ii,J)
                end do
                expdata(pt)%hyd_nr(J) = 1.d0/expdata(pt)%hyd_nr(J)
            end do
        endif
        if (error == 0) then
            if (gl(gl_indicator)%N_guests == 1) then
                if ((trim(gl(gl_indicator)%components(2)) == 'co2')) then
                    continue
                elseif ((trim(gl(gl_indicator)%components(2)) == 'methane').or.(trim(gl(gl_indicator)%components(2)) == 'ethane').or.(trim(gl(gl_indicator)%components(2)) == 'propane').or. &
                    & (trim(gl(gl_indicator)%components(2)) == 'oxygen')) then
                    if (expdata(pt)%pressIn > 100.d0) then
                        expdata(pt)%author = trim(expdata(pt)%author) // ' p > 100 MPa'
                        expdata(pt)%EquiltypeIn = trim(expdata(pt)%EquiltypeIn) // ' p > 100 MPa'
                    endif
                elseif ((trim(gl(gl_indicator)%components(2)) == 'argon')) then
                    if (expdata(pt)%pressIn > 720.d0) then
                        expdata(pt)%author = trim(expdata(pt)%author) // ' p > 720 MPa'
                        expdata(pt)%EquiltypeIn = trim(expdata(pt)%EquiltypeIn) // ' p > 720 MPa'
                    endif
                elseif ((trim(gl(gl_indicator)%components(2)) == 'nitrogen')) then
                    if (expdata(pt)%pressIn > 650.d0) then
                        expdata(pt)%author = trim(expdata(pt)%author) // ' p > 650 MPa'
                        expdata(pt)%EquiltypeIn = trim(expdata(pt)%EquiltypeIn) // ' p > 650 MPa'
                    endif
                endif
            endif
            ii = len(trim(expdata(pt)%author))
            j = len(trim(expdata(pt)%EquiltypeIn))
            if (trim(expdata(pt)%EquiltypeIn(1:1)) == 'h') then
                if (size(expdata(pt)%phasefrac,1) == 2) then
                    write(Resulthydcompfileunit,3027)expdata(pt)%author,expdata(pt)%EquiltypeIn(2:),expdata(pt)%pressIn, expdata(pt)%TempIn, expdata(pt)%hyd_nrIn(1),  expdata(pt)%hyd_nr(1),  expdata(pt)%occupIn(1,1),  expdata(pt)%occup(1,1), &
                        & expdata(pt)%occupIn(2,1),  expdata(pt)%occup(2,1)
                elseif (size(expdata(pt)%phasefrac,1) == 3) then
                    write(Resulthydcompfileunit,3027)expdata(pt)%author,expdata(pt)%EquiltypeIn(2:),expdata(pt)%pressIn, expdata(pt)%Tempcalc, expdata(pt)%hyd_nrIn(1),  expdata(pt)%hyd_nr(1),  expdata(pt)%occupIn(1,1),  expdata(pt)%occup(1,1), &
                        & expdata(pt)%occupIn(2,1),  expdata(pt)%occup(2,1)
                endif
            else
                if (size(expdata(pt)%phasefrac,1) == 2) then
                    write(Resultfileunit,3026)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, expdata(pt)%x_fluid1(2), expdata(pt)%x_fluid1In(2) - expdata(pt)%x_fluid1(2)
                elseif (size(expdata(pt)%phasefrac,1) == 3) then
                    deviation_abs_cumm = deviation_abs_cumm + dabs(expdata(pt)%TempIn - expdata(pt)%tempcalc)
                    deviation_cumm = deviation_cumm + expdata(pt)%TempIn - expdata(pt)%tempcalc
                    write(Resultfileunit,3023)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, expdata(pt)%tempcalc, expdata(pt)%TempIn - expdata(pt)%tempcalc, expdata(pt)%HDRT_STRUCTURE_EXP
                elseif (size(expdata(pt)%phasefrac,1) == 4) then
                    if (gl(gl_indicator)%N_guests == 1) then
                        write(Resultfileunit,3025)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, expdata(pt)%tempcalc, -9123456
                    else
                        write(Resultfileunit,3023)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, expdata(pt)%tempcalc, expdata(pt)%TempIn - expdata(pt)%tempcalc, expdata(pt)%HDRT_STRUCTURE_EXP
                    endif
                else
                    write(*,*)'Warning. Deviation of current experimental data point could not be calculated. AD and AAD will be wrong!'
                endif
            endif
        else
            if (trim(expdata(pt)%EquiltypeIn(1:1)) == 'h') then
                write(Resultfileunit,3028)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, Error, Error, Error, Error, Error, Error
            else
                write(Resultfileunit,3024)expdata(pt)%author,expdata(pt)%EquiltypeIn,expdata(pt)%pressIn, expdata(pt)%TempIn, Error, Error
            endif
            pts_err = pts_err + 1
        endif
    enddo deviationloop

500 call messages(gl(gl_indicator), pel,205,trim(''),trim(''),0,0,0)

    write(Resultfileunit,'(";",";",A,";",f10.6,";",A)')'AD = ',deviation_cumm/(size(expdata) - pts_err),' K'
    write(Resultfileunit,'(";",";",A,";",f10.6,";",A)')'AAD = ',deviation_abs_cumm/(size(expdata) - pts_err),' K'
    write(Resultfileunit,'(";",";",A,";",";",I6)')'Total exp data points: ',size(expdata)
    write(Resultfileunit,'(";",";",A,";",";",I6)')'Failed exp data points: ',pts_err

    close(Resultfileunit)
    close(Resulttwophasefileunit)
    close(Resulthydcompfileunit)


3021 format (a<ii>, ";", a<j>, ";", a10, ";", a10, ";", a10, ";", a10, ";", a10, ";", a10, ";", a10, ";", a10)
3022 format (a<ii>, ";", a<j>, ";", a10, ";", a10, ";", a10, ";", a10)
3023 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", f10.4, ";", f10.6, ";",a3)
3024 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", I10, ";", I10)
3025 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", f10.4, ";", I)
3026 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", f12.8, ";", f12.8)
3027 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", f10.6, ";", f10.6, ";", f10.8, ";", f10.8, ";", f10.8, ";", f10.8)
3028 format (a<ii>, ";", a<j>, ";", f12.6, ";", f10.4, ";", I10, ";", I10, ";", I10, ";", I10, ";", I10, ";", I10)


    end subroutine hydratedatadev
    !******************************************************************************
    !******************************************************************************

    module subroutine hydratedatadev_readfile(gl, pel, expdata, Datenfilepathin, fluidspath)



    implicit none
    type(type_gl), dimension(:) :: gl
    type(type_hdrt_pheq) :: pel
    type(hdrt_experimental_data), dimension(:), allocatable :: expdata

    integer :: Datenfileunit
    character(255) :: Datenfilepathin, dummy1

    integer :: errorin, errorout, head, pt, pts_total, i, i_old, column
    logical :: next_author

    !TREND setup variables
    !#####################
    character(12) :: input
    integer, dimension(30) :: EOS_indicator
    integer :: MIX_indicator
    double precision, dimension(30) :: moles
    character(255) :: fluidspath
    integer :: error
    !#####################

    !if (gl(1)%n_guests == 1) then
    !    dummy1 = DatenfilepathIn(1:len(trim(DatenfilepathIn))-4) // '_un.hdt'
    !    open(newunit = Datenfileunit, file=trim(dummy1), status='unknown', action='read', iostat=errorin)
    !    if (errorin /= 0) open(newunit = Datenfileunit, file=trim(DatenfilepathIn), status='unknown', action='read', iostat=errorin)
    !else
    open(newunit = Datenfileunit, file=trim(DatenfilepathIn), status='unknown', action='read', iostat=errorin)
    !endif

    if (errorin /= 0) then
        !write(*,*)'Input file does not exist or is used by antoher program.'
        call messages(gl(1), pel,200,trim(''),trim(''),0,0,0)
        return
    endif
    head = 0
    if (gl(1)%N_hdrts == 2) then
        do
            read(Datenfileunit,*, end = 400)dummy1
            if (dummy1(1:1) /= "!") exit
            head = head + 1
        enddo
    endif
    pt = -1
    do
        read(Datenfileunit,*, end = 400)dummy1
        pt = pt + 1
    enddo
400 rewind(Datenfileunit)
    pts_total = pt
    if (.not. allocated(expdata))allocate (expdata(pt))

    if (errorin /= 0) then
        !write(*,*)'Input file does not exist or is used by antoher program.'
        call messages(gl(1), pel,200,trim(''),trim(''),0,0,0)
        return
    end if
    if (gl(1)%N_hdrts == 2) then
        do i = 1,head
            read(Datenfileunit,*, end = 400)dummy1
        enddo
        next_author = .true.
    else
        read(Datenfileunit,'(A)')dummy1
    endif
    read(Datenfileunit,'(A)')expdata%readr_in!head
    close(Datenfileunit)
    !readloop: do pt = 1, pts_total
    !    call hydratedatadev_processdata(gl,expdata)
    !enddo readloop

    pt = 1
    call hydratedatadev_processdata(gl,expdata, pt, pts_total)

500 call messages(gl(1), pel,204,trim(''),trim(''),0,0,0)
    end subroutine hydratedatadev_readfile

    module subroutine hydratedatadev_processdata(gl, expdata, pt, pts_total)



    implicit none
    type(type_gl), dimension(:) :: gl
    type(type_hdrt_pheq) :: pel
    type(hdrt_experimental_data), dimension(:), allocatable :: expdata

    integer :: Datenfileunit
    character(255) :: Datenfilepathin, dummy1

    integer :: errorin, errorout, head, pt, pts_total, i, i_old, column, pt_orig, ii
    logical :: next_author

    !TREND setup variables
    !#####################
    character(12) :: input
    integer, dimension(30) :: EOS_indicator
    integer :: MIX_indicator
    double precision, dimension(30) :: moles
    character(255) :: fluidspath
    integer :: error
    !#####################

    !type(hdrt_experimental_data), dimension(:), allocatable :: expdata
    !integer :: pt, column

    !type(type_hdrt_pheq) :: pel
    !type(hdrt_experimental_data), dimension(:), allocatable :: expdata


    !############################################

    pt_orig = pt

    readloop: do pt = pt_orig, pts_total
        call hdrt_init_expdata(expdata(pt))

        if (gl(1)%N_hdrts == 2) then   ! deviation to pure hydrates => NEVER mixed hydrates during calculation
            if (expdata(pt)%readr_in(1:1) == '!') then
                next_author = .true.
                cycle
            else
                if ((next_author .eqv. .true.).or.(pt == 1)) then
                    do while (index(expdata(pt)%readr_in,achar(9)) /= 0)
                        expdata(pt)%readr_in(index(expdata(pt)%readr_in,achar(9)):index(expdata(pt)%readr_in,achar(9))) = " "
                    enddo
                    read(expdata(pt)%readr_in,'(A)')expdata(pt)%author
                    next_author = .false.
                    cycle
                endif
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn
                expdata(pt)%weightIn = 1.d0
                expdata(pt)%author = expdata(pt-1)%author
                expdata(pt)%xIn(1) = 0.5d0
                expdata(pt)%xIn(2) = 0.5d0
                if (errorout /= 0) then
                    continue
                    cycle
                elseif  (dabs(expdata(pt)%weightIn) <= 1.0d-12) then
                    continue
                    cycle
                endif
            endif



        elseif (gl(1)%N_hdrts > 2) then    ! mixed hydrates in general, some data points might be pure hydrates with one of the components
            i_old = 1
            column = 1
            if (expdata(pt)%readr_in(1:1) /= '!') then
                lineloop: do i = 1,len(trim(expdata(pt)%readr_in))
                    if (expdata(pt)%readr_in(i:i) == ";") then
                        if (i > i_old) then
                            select case (column)
                            case (1)
                                expdata(pt)%author = expdata(pt)%readr_in(i_old:i-1)
                                if (expdata(pt)%author(1:1) == '!') then
                                endif
                            case (2)
                                expdata(pt)%EquiltypeIn = expdata(pt)%readr_in(i_old:i-1)
                            case (3)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%weightIn
                                if (errorin /= 0) exit
                            case (4)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%TempIn
                                if (errorin /= 0) exit
                            case (5)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%pressIn
                                if (errorin /= 0) exit
                            case (6)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%ScenarioIn
                                if (errorin /= 0) exit
                            case (7)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%VcellIn
                                if (errorin /= 0) exit
                            case (8)
                                read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%VliqIn
                                if (errorin /= 0) exit
                            case (9)
                                if (expdata(pt)%ScenarioIn == 0) then
                                    dummy1 = expdata(pt)%readr_in(i_old:i-1)
                                    ii = 1
                                    call split_dbl(dummy1, expdata(pt)%xIn, ii, ",")
                                    if (gl(1)%N_hdrts == 3) then
                                        read(dummy1,*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3)
                                        !read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3)
                                        if (errorin /= 0) exit
                                    elseif (gl(1)%N_hdrts == 4) then
                                        read(dummy1,*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3),expdata(pt)%xIn(4)
                                        !read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3),expdata(pt)%xIn(4)
                                        if (errorin /= 0) exit
                                    elseif (gl(1)%N_hdrts == 5) then
                                        read(dummy1,*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3),expdata(pt)%xIn(4),expdata(pt)%xIn(5)
                                        !read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%xIn(1),expdata(pt)%xIn(2),expdata(pt)%xIn(3),expdata(pt)%xIn(4)
                                        if (errorin /= 0) exit
                                    endif
                                endif
                            case (10)
                                if ((expdata(pt)%ScenarioIn == 1).or.(expdata(pt)%ScenarioIn == 2).or.(expdata(pt)%ScenarioIn == 3)) then
                                    if (gl(1)%N_hdrts == 3) then
                                        read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%x_vap_givenIn_Equil,expdata(pt)%x_vap_givenIn(2),expdata(pt)%x_vap_givenIn(3)
                                        if (errorin /= 0) exit
                                    elseif (gl(1)%N_hdrts == 4) then
                                        read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%x_vap_givenIn_Equil,expdata(pt)%x_vap_givenIn(2),expdata(pt)%x_vap_givenIn(3),expdata(pt)%x_vap_givenIn(4)
                                        if (errorin /= 0) exit
                                    elseif (gl(1)%N_hdrts == 5) then
                                        read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%x_vap_givenIn_Equil,expdata(pt)%x_vap_givenIn(2),expdata(pt)%x_vap_givenIn(3),expdata(pt)%x_vap_givenIn(4),expdata(pt)%x_vap_givenIn(5)
                                        if (errorin /= 0) exit
                                    endif
                                endif
                            case (11)
                                if (expdata(pt)%ScenarioIn == 3) then
                                    read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%x_h2oIn
                                    if (errorin /= 0) exit
                                endif
                            case (14)
                                if (expdata(pt)%readr_in(i_old:i-1) /= "") then
                                    read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%hdrt_structure_exp
                                    if (errorin /= 0) exit
                                endif
                            case (15)
                                if (expdata(pt)%readr_in(i_old:i-1) /= "") then
                                    read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%Tempfill
                                    if (errorin /= 0) exit
                                endif
                            case (16)
                                if (expdata(pt)%readr_in(i_old:i-1) /= "") then
                                    read(expdata(pt)%readr_in(i_old:i-1),*,iostat=errorin)expdata(pt)%pressfill
                                    if (errorin /= 0) exit
                                endif
                            end select
                        endif
                        i_old = i + 1
                        column = column + 1
                    endif
                enddo lineloop
            endif
        endif
        expdata(pt)%tempcalc = expdata(pt)%TempIn

    enddo readloop

    pt = pt_orig

    end subroutine hydratedatadev_processdata

    module subroutine expdata_fill(gl, gl_indicator, pel, expdata, pt, fluidspath)



    implicit none

    type(type_gl), dimension(:) :: gl
    type(type_gl), allocatable :: gl_fluid
    type(type_hdrt_pheq) :: pel
    type(hdrt_experimental_data), dimension(:), allocatable :: expdata

    character(255), intent(inout) :: fluidspath !,pathin
    character(255) :: path_local
    !integer, intent(inout) :: error  !MB It was with the intent(inout) but then we have a compiling error....
    integer:: error

    integer :: Datenfileunit
    integer :: Resultfileunit


    character(12) :: input
    !character(255) :: devfilepath
    character(30), dimension(30) :: components_orig

    !Variables at/for the three-phase equilibrium
    double precision :: rho_feed, n_water
    double precision:: rhoph1_est, rhoph2_est

    integer, dimension(30) :: EOS_Indicator
    integer :: MIX_indicator
    integer :: errorin, errorout, i, i_old, pt, pts_err, iTp, column, skip_data, pt_orig
    integer :: iFlash, iter
    integer :: pointsexp_total, ncomp_orig
    integer, dimension(2) :: solidtype_copy
    integer :: hdrt_structure_flag, gl_indicator

    logical :: show_once = .false., next_author

    !Variables at the Quadruple point
    !character(len=6), dimension(4)  :: Q_point
    !double precision, dimension (4) :: press_Q, Temp_Q
    !double precision, dimension(4,30) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    !logical, dimension (4):: Qexist

    double precision, dimension(30) :: composition

    error = 0


    gl_indicator = 1 !if N_hdrts = 2 -> pure hydrates always
    !if N_hdrts > 2 -> 1: mixed hydrates
    !                  2: pure hydrates with guest 1
    !                  n: pure hydrates with guest n-1
    iter = 0
    skip_data = 0
    error = 0
    !gl%N_hdrts = ncomp_orig

    if (pel%show_progress .eqv. .true.) then
        !write(*,'(''+'',A,I5,A4,I5)')'Calculating Deviation of Point ',pts_total,' of ',pointsexp_total
        call messages(gl(gl_indicator), pel,202,trim(''),trim(''),pt-count(expdata(1:pt)%weightin<1.d-14),count(expdata(:)%weightin>1.d-14),0)
    endif

    if ((expdata(pt)%readr_in(1:1) == "!").or.(expdata(pt)%weightin == 0.d0)) return
    !if fill state not given -> guess
    if (expdata(pt)%Tempfill < 1.d-14) expdata(pt)%Tempfill = 293.15d0
    !if (expdata(pt)%pressfill < 0.d0) expdata(pt)%pressfill = 3.d0
    if (expdata(pt)%pressfill < 1.d-14) expdata(pt)%pressfill = expdata(pt)%pressIn

    if (gl(1)%N_hdrts == 2) then
        if ((trim(expdata(pt)%EquiltypeIn) == 'LwLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLpH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLeH')) then
            expdata(pt)%EquiltypeGen = 'LwLxH'
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.9999d0
            expdata(pt)%x_fluid2(2) = 0.0001d0
            expdata(pt)%iPhase_sol = 1 !1 - Assume liquid / liquid equilibrium --> with solidtype = LLH equilibrium, 2: 2 - Assume vapor / liquid equilibrium --> with solidtype = VLH equilibrium
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VLwH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVLwH')) then
            if (trim(expdata(pt)%EquiltypeIn) == 'VLwH') then
                expdata(pt)%EquiltypeGen = 'VLwH'
            elseif (trim(expdata(pt)%EquiltypeIn) == 'hVLwH') then
                expdata(pt)%EquiltypeGen = 'hVLwH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_fluid2In(1), expdata(pt)%x_fluid2In(2), expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.9999d0
            expdata(pt)%x_fluid2(2) = 0.0001d0
            expdata(pt)%iPhase_sol = 2 !1 - Assume liquid / liquid equilibrium --> with solidtype = LLH equilibrium, 2: 2 - Assume vapor / liquid equilibrium --> with solidtype = VLH equilibrium
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLpH').or.&
            & (trim(expdata(pt)%EquiltypeIn) == 'hVLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVLpH')) then
            if ((trim(expdata(pt)%EquiltypeIn) == 'VLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLpH')) then
                expdata(pt)%EquiltypeGen = 'VLxH'
            elseif ((trim(expdata(pt)%EquiltypeIn) == 'hVLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVLpH')) then
                expdata(pt)%EquiltypeGen = 'VLxH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_fluid2In(1), expdata(pt)%x_fluid2In(2), expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.01d0
            expdata(pt)%x_fluid2(2) = 0.99d0
            expdata(pt)%iPhase_sol = 2 !1 - Assume liquid / liquid equilibrium --> with solidtype = LLH equilibrium, 2: 2 - Assume vapor / liquid equilibrium --> with solidtype = VLH equilibrium
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VHIw').or.(trim(expdata(pt)%EquiltypeIn) == 'hVHIw')) then
            if (trim(expdata(pt)%EquiltypeIn) == 'VHIw') then
                expdata(pt)%EquiltypeGen = 'VHIw'
            elseif (trim(expdata(pt)%EquiltypeIn) == 'hVHIw') then
                expdata(pt)%EquiltypeGen = 'VHIw'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%x_solIn(1), expdata(pt)%x_solIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.01d0
            expdata(pt)%x_fluid1(2) = 0.99d0
            expdata(pt)%x_sol = 0.d0   !initialize
            expdata(pt)%x_sol(1) = 1.d0 !Pure solid water
            gl(gl_indicator)%solidtype(1) = 1                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VLcLw').or.(trim(expdata(pt)%EquiltypeIn) == 'VLpLw').or.(trim(expdata(pt)%EquiltypeIn) == 'VLeLw')) then
            expdata(pt)%EquiltypeGen = 'VLxLw'
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.01d0
            expdata(pt)%x_fluid2(2) = 0.99d0
            expdata(pt)%x_fluid3(1) = 0.9999d0
            expdata(pt)%x_fluid3(2) = 0.0001d0
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 0                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif (trim(expdata(pt)%EquiltypeIn) == 'LwHIc') then
            expdata(pt)%EquiltypeGen = 'LwHIc'
            gl(gl_indicator)%solid_pos = 2
            expdata(pt)%x_fluid1(1) = 0.99d0
            expdata(pt)%x_fluid1(2) = 0.01d0
            expdata(pt)%x_sol = 0.d0   !initialize
            expdata(pt)%x_sol(2) = 1.d0 !Pure solid co2
            gl(gl_indicator)%solidtype(1) = 2                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif (trim(expdata(pt)%EquiltypeIn) == 'VLwIw') then
            expdata(pt)%EquiltypeGen = 'VLwIw'
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.01d0
            expdata(pt)%x_fluid1(2) = 0.99d0
            expdata(pt)%x_fluid2(1) = 0.999d0
            expdata(pt)%x_fluid2(2) = 0.001d0
            expdata(pt)%x_sol = 0.d0   !initialize
            expdata(pt)%x_sol(1) = 1.d0 !Pure solid water
            expdata(pt)%iPhase_sol = 2
            gl(gl_indicator)%solidtype(1) = 1                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 0                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif (trim(expdata(pt)%EquiltypeIn) == 'LwHIw') then
            expdata(pt)%EquiltypeGen = 'LwHIw'
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.999d0
            expdata(pt)%x_fluid1(2) = 0.001d0
            expdata(pt)%x_sol = 0.d0   !initialize
            expdata(pt)%x_sol(1) = 1.d0 !Pure solid water
            gl(gl_indicator)%solidtype(1) = 1                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif (trim(expdata(pt)%EquiltypeIn) == 'VLwHIw') then
            expdata(pt)%EquiltypeGen = 'VLwHIw'
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.9999d0
            expdata(pt)%x_fluid2(2) = 0.0001d0
            expdata(pt)%x_sol = 0.d0   !initialize
            expdata(pt)%x_sol(1) = 1.d0 !Pure solid water
            expdata(pt)%iPhase_sol = 2
            gl(gl_indicator)%solidtype(1) = 1                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VLwLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLwLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLwLpH')) then
            expdata(pt)%EquiltypeGen = 'VLwLxH'
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%x_fluid2(1) = 0.9999d0
            expdata(pt)%x_fluid2(2) = 0.0001d0
            expdata(pt)%x_fluid3(1) = 0.0001d0
            expdata(pt)%x_fluid3(2) = 0.9999d0
            expdata(pt)%iPhase_sol = 2
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1                                       !0: no hydrate phase exists, 1: hydrate phase exists
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'VH').or.(trim(expdata(pt)%EquiltypeIn) == 'hVH')) then
            if (trim(expdata(pt)%EquiltypeIn) == 'VH') then
                expdata(pt)%EquiltypeGen = 'VH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2)
            elseif (trim(expdata(pt)%EquiltypeIn) == 'hVH') then
                expdata(pt)%EquiltypeGen = 'VH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.001d0
            expdata(pt)%x_fluid1(2) = 0.999d0
            expdata(pt)%iPhase_sol = 2                                             !0 - Let the density solver choose based on the gibbs energy, 1 - Assume liquid, 2 - Assume vapor
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'LcH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'LeH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'LpH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLpH')) then
            if ((trim(expdata(pt)%EquiltypeIn) == 'LcH').or.(trim(expdata(pt)%EquiltypeIn) == 'LeH').or.(trim(expdata(pt)%EquiltypeIn) == 'LpH')) then
                expdata(pt)%EquiltypeGen = 'LxH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2)
            elseif ((trim(expdata(pt)%EquiltypeIn) == 'hLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLpH')) then
                expdata(pt)%EquiltypeGen = 'LxH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            expdata(pt)%EquiltypeGen = 'LxH'
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.01d0
            expdata(pt)%x_fluid1(2) = 0.99d0
            expdata(pt)%iPhase_sol = 1                                             !0 - Let the density solver choose based on the gibbs energy, 1 - Assume liquid, 2 - Assume vapor
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1
        elseif ((trim(expdata(pt)%EquiltypeIn) == 'LwH').or.(trim(expdata(pt)%EquiltypeIn) == 'hLwH')) then
            if (trim(expdata(pt)%EquiltypeIn) == 'LwH') then
                expdata(pt)%EquiltypeGen = 'LwH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2)
            elseif (trim(expdata(pt)%EquiltypeIn) == 'hLwH') then
                expdata(pt)%EquiltypeGen = 'LwH'
                read(expdata(pt)%readr_in,*,iostat=errorout)expdata(pt)%EquiltypeIn,expdata(pt)%weightIn,expdata(pt)%TempIn,expdata(pt)%pressIn, expdata(pt)%xIn(1), expdata(pt)%xIn(2), expdata(pt)%x_fluid1In(1), expdata(pt)%x_fluid1In(2), &
                    & expdata(pt)%x_hydIn(1), expdata(pt)%x_hydIn(2), expdata(pt)%hyd_nrIn(1), expdata(pt)%occupIn(1,1), expdata(pt)%occupIn(2,1)
            endif
            gl(gl_indicator)%solid_pos = 1
            expdata(pt)%x_fluid1(1) = 0.99d0
            expdata(pt)%x_fluid1(2) = 0.01d0
            expdata(pt)%iPhase_sol = 1                                             !0 - Let the density solver choose based on the gibbs energy, 1 - Assume liquid, 2 - Assume vapor
            gl(gl_indicator)%solidtype(1) = 0                                       !0: no solid phase exists, 1: solid phase is Iw, 2: solid phase is Ic
            gl(gl_indicator)%solidtype(2) = 1
        else
            return
        endif
        if (maxval(expdata(pt)%xIn) < 1.d-14) then
            expdata(pt)%xIn(1) = 0.5d0
            expdata(pt)%xIn(2) = 0.5d0
        endif
    elseif (gl(1)%N_hdrts > 2) then
        if ((expdata(pt)%ScenarioIn == 0).or.(expdata(pt)%ScenarioIn == 3)) then !overall composition is given
            if (expdata(pt)%ScenarioIn == 3) then
                expdata(pt)%xIn(1) = expdata(pt)%x_h2oIn
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%xIn(i) = (1.d0 - expdata(pt)%xIn(1))*expdata(pt)%x_vap_givenIn(i)
                enddo
            endif
            expdata(pt)%x_fluid1(1) = 1.d-3    !vap  water - small amount
            do i = 2,gl(gl_indicator)%N_hdrts
                expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%xIn(i)/(sum(expdata(pt)%xIn(2:gl(gl_indicator)%N_hdrts))) ! vap  guest1 with ratio of guest1/guest2 in overall comp
            enddo
            if (trim(expdata(pt)%EquiltypeIn) == 'VLwH') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'VLwH'
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
            elseif (trim(expdata(pt)%EquiltypeIn) == 'VHIw') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'VHIw'
                expdata(pt)%x_sol = 0.d0   !initialize
                expdata(pt)%x_sol(1) = 1.d0 !Pure solid water
                gl(gl_indicator)%solid_pos = 1
            elseif (trim(expdata(pt)%EquiltypeIn) == 'VLwLxH') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(4))
                expdata(pt)%EquiltypeGen = 'VLwLxH'
                expdata(pt)%x_fluid1(1) = 1.d-5    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%xin(i)/(1.d0-expdata(pt)%xin(1)) ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
                expdata(pt)%x_fluid3(1) = 0.001d0    !liq3 water
                expdata(pt)%x_fluid3(2) = 0.97d0
                expdata(pt)%x_fluid3(3) = 1.d0 - sum(expdata(pt)%x_fluid3(1:gl(gl_indicator)%N_hdrts))
            endif
            expdata(pt)%x_hyd = 0.D0
            expdata(pt)%phasefrac(1) = 1.d0 - expdata(pt)%xIn(1) !Phasefrac of vapor phase
            expdata(pt)%phasefrac(3) = expdata(pt)%xIn(1) !Phasefrac of liquidphase
        elseif (expdata(pt)%ScenarioIn == 1) then ! vapor composition at equilibrium is given
            !VLwH
            !Estimate expdata(pt)%xIn from given vapor composition
            expdata(pt)%xIn(1) = 0.5d0 !water
            do i = 2,gl(gl_indicator)%N_hdrts
                expdata(pt)%xIn(i) = (1.d0 - expdata(pt)%xIn(1))*expdata(pt)%x_vap_givenIn(i)
            enddo
            if ((trim(expdata(pt)%EquiltypeIn) == 'VH').or.(trim(expdata(pt)%EquiltypeIn) == 'HV')) then
                expdata(pt)%EquiltypeIn = 'VH'
                expdata(pt)%EquiltypeGen = 'VH'
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))

            elseif (trim(expdata(pt)%EquiltypeIn) == 'VLwH') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'VLwH'
                expdata(pt)%x_fluid1(1) = 1.d-3    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i) ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
                !pure hydrates?
                if (expdata(pt)%x_vap_givenIn(2) < 1.d-14) then
                    gl_indicator = 3
                    expdata(pt)%xIn(2) = expdata(pt)%xIn(3)
                    expdata(pt)%xIn(3) = 0.d0
                    expdata(pt)%x_fluid1(2) = expdata(pt)%x_fluid1(3)
                    expdata(pt)%x_fluid1(3) = 0.d0
                    expdata(pt)%x_fluid2(2) = 1.d-3
                    expdata(pt)%x_fluid2(3) = 0.d0
                elseif  (expdata(pt)%x_vap_givenIn(3) < 1.d-14) then
                    gl_indicator = 2
                    expdata(pt)%xIn(3) = 0.d0
                    expdata(pt)%x_fluid2(2) = 1.d-3
                    expdata(pt)%x_fluid2(3) = 0.d0
                endif
            elseif ((trim(expdata(pt)%EquiltypeIn) == 'LwLpH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLcH')) then          !MB Copy&paste block from VLwH above because LwLpH had no allocation for fluid in expdata
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'LwLxH'
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq1 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid1(1) = 1.d-3    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i)
                enddo
            elseif (trim(expdata(pt)%EquiltypeIn) == 'VHIw') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'VHIw'
                expdata(pt)%x_fluid1(1) = 1.d-5    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i) ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_sol = 0.d0
                expdata(pt)%x_sol(1) = 1.d0
                if (expdata(pt)%x_vap_givenIn(2) < 1.d-14) then
                    gl_indicator = 3
                    expdata(pt)%xIn(2) = expdata(pt)%xIn(3)
                    expdata(pt)%xIn(3) = 0.d0
                    expdata(pt)%x_fluid1(2) = expdata(pt)%x_fluid1(3)
                    expdata(pt)%x_fluid1(3) = 0.d0
                elseif  (expdata(pt)%x_vap_givenIn(3) < 1.d-14) then
                    gl_indicator = 2
                    expdata(pt)%xIn(3) = 0.d0
                endif
            elseif (trim(expdata(pt)%EquiltypeIn) == 'VLwHIw') then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(4))
                expdata(pt)%EquiltypeGen = 'VLwHIw'
                expdata(pt)%x_fluid1(1) = 1.d-5    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i) ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
                expdata(pt)%x_sol = 0.d0
                expdata(pt)%x_sol(1) = 1.d0
                if (expdata(pt)%x_vap_givenIn(2) < 1.d-14) then
                    gl_indicator = 3
                    expdata(pt)%xIn(2) = expdata(pt)%xIn(3)
                    expdata(pt)%xIn(3) = 0.d0
                    expdata(pt)%x_fluid1(2) = expdata(pt)%x_fluid1(3)
                    expdata(pt)%x_fluid1(3) = 0.d0
                elseif  (expdata(pt)%x_vap_givenIn(3) < 1.d-14) then
                    gl_indicator = 2
                    expdata(pt)%xIn(3) = 0.d0
                endif
            elseif ((trim(expdata(pt)%EquiltypeIn) == 'VLwLcH').or.(trim(expdata(pt)%EquiltypeIn) == 'VLwLxH')) then
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(4))
                expdata(pt)%EquiltypeGen = 'VLwLxH'
                expdata(pt)%x_fluid1(1) = 1.d-5    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i) ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
                expdata(pt)%x_fluid3(1) = 0.001d0    !liq3 water
                expdata(pt)%x_fluid3(2) = 0.93d0
                expdata(pt)%x_fluid3(3) = 1.d0 - sum(expdata(pt)%x_fluid3(1:gl(gl_indicator)%N_hdrts))
                !do i = 2,gl(gl_indicator)%N_hdrts
                !    expdata(pt)%x_fluid3(i) = (1.d0 - expdata(pt)%x_fluid3(1))/gl(gl_indicator)%N_guests
                !enddo
                if (expdata(pt)%x_vap_givenIn(2) < 1.d-14) then
                    gl_indicator = 3
                    expdata(pt)%xIn(2) = expdata(pt)%xIn(3)
                    expdata(pt)%xIn(3) = 0.d0
                    expdata(pt)%x_fluid1(2) = expdata(pt)%x_fluid1(3)
                    expdata(pt)%x_fluid1(3) = 0.d0
                elseif  (expdata(pt)%x_vap_givenIn(3) < 1.d-14) then
                    gl_indicator = 2
                    expdata(pt)%xIn(3) = 0.d0
                endif
            endif

            expdata(pt)%x_hyd = 0.D0
            expdata(pt)%phasefrac(1) = 1.d0 - expdata(pt)%xIn(1) !Phasefrac of vapor phase
            expdata(pt)%phasefrac(3) = expdata(pt)%xIn(1) !Phasefrac of liquidphase

        elseif (expdata(pt)%ScenarioIn == 2) then! vapor composition of feed is given
            if (.not.allocated(gl_fluid)) allocate(gl_fluid)

            input = 'tp'
            EOS_indicator = 0

            EOS_indicator(1:gl(gl_indicator)%N_guests) = 1
            MIX_indicator = 1
            gl_fluid = gl(gl_indicator)
            components_orig = ''
            components_orig(1:gl_fluid%N_guests) = gl_fluid%components(2:gl_fluid%N_guests+1)
            composition = 0.d0
            composition(1:29) = expdata(pt)%x_vap_givenIn(2:30)
            rho_feed = prop_eos(gl_fluid,2,input, expdata(pt)%Tempfill, expdata(pt)%pressfill, components_orig, composition, EOS_indicator, MIX_indicator, fluidspath, error)
            expdata(pt)%n_vap = 0.d0
            do i = 2,gl(gl_indicator)%N_hdrts
                expdata(pt)%n_vap(i) = rho_feed * expdata(pt)%VcellIn * 1.d-3 * expdata(pt)%x_vap_givenIn(i)
            enddo
            !fluidl = 'water'
            gl_fluid%components = ''
            gl_fluid%components(1) = 'water'
            expdata(pt)%xIn = 0.d0
            expdata(pt)%xIn(1) =1.d0
            !molesl = '1'
            EOS_indicator = 0
            EOS_indicator(1) = 1
            MIX_indicator = 1
            if (expdata(pt)%VliqIn < 1.d-14) then
                expdata(pt)%xIn(1) = 0.5
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%xIn(i) = expdata(pt)%x_vap_givenIn(i)*(1.d0-expdata(pt)%xIn(1))
                enddo
            else

                n_water = prop_eos(gl_fluid,2,input, expdata(pt)%Tempfill, expdata(pt)%pressfill, gl_fluid%components, expdata(pt)%xIn, EOS_indicator, MIX_indicator, fluidspath, error) * expdata(pt)%VliqIn * 1.d-3
                expdata(pt)%xIn(1) = n_water/(n_water + sum(expdata(pt)%n_vap))
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%xIn(i) = expdata(pt)%n_vap(i)/(n_water + sum(expdata(pt)%n_vap))
                enddo
            endif
            expdata(pt)%x_fluid1(1) = 1.d-3    !vap  water - small amount
            do i = 2,gl(gl_indicator)%N_hdrts
                expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%xIn(i)/(sum(expdata(pt)%xIn(2:gl(gl_indicator)%N_hdrts))) ! vap  guest1 with ratio of guest1/guest2 in overall comp
            enddo
            if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
            if (trim(expdata(pt)%EquiltypeIn) == 'VLwH') then
                expdata(pt)%EquiltypeGen = 'VLwH'
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq2 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests
                enddo
            elseif (trim(expdata(pt)%EquiltypeIn) == 'VHIw') then
                expdata(pt)%EquiltypeGen = 'VHIw'
                expdata(pt)%x_sol = 0.d0
                expdata(pt)%x_sol(1) = 1.d0
            elseif ((trim(expdata(pt)%EquiltypeIn) == 'LwLpH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLeH').or.(trim(expdata(pt)%EquiltypeIn) == 'LwLcH')) then          !MB Copy&paste block from VLwH above because LwLpH had no allocation for fluid in expdata
                if (.not. (allocated(expdata(pt)%phasefrac))) allocate(expdata(pt)%phasefrac(3))
                expdata(pt)%EquiltypeGen = 'LwLxH'
                expdata(pt)%x_fluid2(1) = 0.999d0    !liq1 water
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid2(i) = (1.d0 - expdata(pt)%x_fluid2(1))/gl(gl_indicator)%N_guests ! vap  guest1 with ratio of guest1/guest2 in overall comp
                enddo
                expdata(pt)%x_fluid1(1) = 1.d-3    !vap  water - small amount
                do i = 2,gl(gl_indicator)%N_hdrts
                    expdata(pt)%x_fluid1(i) = (1.d0-expdata(pt)%x_fluid1(1))*expdata(pt)%x_vap_givenIn(i)
                enddo
            endif
            expdata(pt)%x_hyd = 0.D0
            expdata(pt)%phasefrac(1) = 1.d0 - expdata(pt)%xIn(1) !Phasefrac of vapor phase
            expdata(pt)%phasefrac(3) = expdata(pt)%xIn(1) !Phasefrac of liquidphase
        endif
    endif

    if (expdata(pt)%hdrt_structure_exp /= "") then
        call uppertolower_char(expdata(pt)%hdrt_structure_exp,len(trim(expdata(pt)%hdrt_structure_exp)))
        if ((trim(expdata(pt)%hdrt_structure_exp) == "s1") .or. (trim(expdata(pt)%hdrt_structure_exp) == "si")) then
            hdrt_structure_flag = 1
        elseif ((trim(expdata(pt)%hdrt_structure_exp) == "s2") .or. (trim(expdata(pt)%hdrt_structure_exp) == "sii")) then
            hdrt_structure_flag = 2
        else
            read(expdata(pt)%hdrt_structure_exp,*)hdrt_structure_flag
        endif
        call hdrt_structure_definition(gl(gl_indicator),hdrt_structure_flag)
    else
        expdata(pt)%hdrt_structure_exp = 'non'
    endif

    if ((trim(expdata(pt)%EquiltypeGen) == 'VLwH').or.(trim(expdata(pt)%EquiltypeGen) == 'VLxH')) then
        gl(gl_indicator)%solidtype(1) = 0  !no solid phase
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 2 !2 - Assume vapor / liquid equilibrium --> with solidtype = VLH equilibrium
    elseif (trim(expdata(pt)%EquiltypeGen) == 'LwLxH') then
        gl(gl_indicator)%solidtype(1) = 0  !no solid phase
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 1 !1 - Assume liquid / liquid equilibrium --> with solidtype = LLH equilibrium
    elseif (trim(expdata(pt)%EquiltypeGen) == 'VHIw') then
        gl(gl_indicator)%solidtype(1) = 1  !1: solid water, 2: solid co2
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 2    !Liqhter phase is vapor
        gl(gl_indicator)%solid_pos = 1
    elseif (trim(expdata(pt)%EquiltypeGen) == 'LwHIc') then
        gl(gl_indicator)%solidtype(1) = 2  !1: solid water, 2: solid co2
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 1    !Single Liquid phase
        gl(gl_indicator)%solid_pos = 2
    elseif (trim(expdata(pt)%EquiltypeGen) == 'LwHIw') then
        gl(gl_indicator)%solidtype(1) = 1
        gl(gl_indicator)%solidtype(2) = 1
        expdata(pt)%iPhase_sol = 1
    elseif (trim(expdata(pt)%EquiltypeGen) == 'VLwIw') then
        gl(gl_indicator)%solidtype(1) = 1
        expdata(pt)%iPhase_sol = 2
    elseif (trim(expdata(pt)%EquiltypeGen) == 'VLwHIw') then
        gl(gl_indicator)%solidtype(1) = 1  !1: solid water, 2: solid co2
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 2    !Liqhter phase is vapor
        gl(gl_indicator)%solid_pos = 1
    elseif (trim(expdata(pt)%EquiltypeGen) == 'VLwLxH') then
        gl(gl_indicator)%solidtype(1) = 0  !1: solid water, 2: solid co2
        gl(gl_indicator)%solidtype(2) = 1  !hydrate phase is emerging phase
        expdata(pt)%iPhase_sol = 2    !Liqhter phase is vapor
        gl(gl_indicator)%solid_pos = 1
    endif
    if (show_once .eqv. .false.) then
        do i = 1,gl(gl_indicator)%ncomp
            if ((expdata(pt)%TempIn < gl(gl_indicator)%tminfluid(i)) .or. (expdata(pt)%TempIn > gl(gl_indicator)%tmaxfluid(i)) .or. (expdata(pt)%pressIn > gl(gl_indicator)%pmaxfluid(i))) then
                call messages(gl(gl_indicator), pel,203,trim(''),trim(''),0,0,0)
                show_once = .true.
            endif
        enddo
    endif



    end subroutine expdata_fill

    end submodule impl
    !DEC$ END IF