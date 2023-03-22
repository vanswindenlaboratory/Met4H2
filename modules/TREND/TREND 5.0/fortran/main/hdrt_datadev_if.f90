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
    module hdrt_data

    use module_all_types
    use module_hdrt_phaselines


    type hdrt_experimental_data
        character(255) :: author
        character(1000) :: readr_in
        double precision :: weightIn
        double precision :: TempIn, pressIn, VcellIn, VliqIn, Tempfill, pressfill
        double precision, dimension(30) :: xIn, x_vap_givenIn, x_fluid1In, x_fluid2In, x_fluid3In, x_fluid4In, x_solIn, x_hydIn, n_vapIn
        double precision :: x_h2oIn

        integer :: ScenarioIn, iPhase_sol
        character(10) :: x_vap_givenIn_Equil
        character(3) :: hdrt_structure_exp
        character(20) :: EquiltypeIn, EquiltypeGen        !Given Type of Equilibrium, which should be calculated

        double precision :: tempcalc, presscalc
        double precision, dimension(2) :: ttest
        !variables for properties of different phases
        !##############################################
        double precision, dimension(30) :: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, n_vap
        double precision, dimension(30,5) :: x_Phase
        double precision, dimension(4) :: rho_ph
        double precision, dimension(5) :: rho
        double precision, dimension(30) :: Chempot, lnfi, hyd_nrIn, hyd_nr

        double precision, dimension(3):: rho_sol

        ! Mixed hydrates 2015
        double precision, dimension(3,30):: occup, occup_single, occup_double, C_iJ
        double precision, dimension(3,30):: occupIn
        double precision, dimension(:), allocatable :: Phasefrac
        !##############################################
    end type hdrt_experimental_data

    interface

    module subroutine hdrt_init_expdata(expdata)
    implicit none
    type(hdrt_experimental_data) :: expdata

    end subroutine hdrt_init_expdata


    !******************************************************************************
    module subroutine hydratedatadev(gl,pel,DatenfilepathIn, pathIn, fluidspath, error)
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
    double precision :: deviation, deviation_abs_cumm, deviation_cumm, prop_eos, rho_feed, n_water
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
    double precision :: prop_eos, rho_feed, n_water
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



    end subroutine expdata_fill

    end interface

    end module hdrt_data
    !DEC$ ELSE
    module hdrt_data
    end module hdrt_data
    !DEC$ END IF