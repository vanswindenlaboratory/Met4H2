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

    ! module for file phasedet_mix.f90
    module phasedet_mix_module
    !global use inclusion
    use module_all_types
    use module_regula_falsi_support

    interface

    module subroutine PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x
    double precision, dimension(5),intent(out):: rho
    double precision, dimension(30,5),intent(out):: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, intent(out):: vapfrac
    integer, intent(out):: errval, nrofphases
    integer, dimension(5), intent(out):: phasetype
    end subroutine PhaseDet


    module subroutine PhaseDet2(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x
    double precision, dimension(5):: rho
    double precision, dimension(30,5):: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision:: vapfrac
    integer:: errval, nrofphases
    integer, dimension(5) :: phasetype
    end subroutine PhaseDet2


    double precision module function tpd(gl,press, temp, x, rho_orig, x_trial, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, rho_orig
    double precision, dimension (30):: x, x_trial
    integer :: errval
    End function tpd

    module subroutine d_tpd_dxj(gl,press, temp, x, x_trial, dtpddxj, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x, x_trial
    integer:: errval
    double precision, dimension(30):: dtpddxj
    End subroutine d_tpd_dxj

    module subroutine Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_trial, tpd_min, iphase, Stable, NrofIter, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x
    double precision:: rho_orig
    double precision, dimension(30):: x_trial
    integer :: NrofIter
    integer :: errval
    double precision :: tpd_min
    logical :: Stable
    integer :: iphase
    End subroutine Succ_Sub_tpd

    module subroutine PhaseDet_td(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, rho_spec, nrofphases, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, rho_spec
    double precision, dimension(30) :: x
    double precision, dimension(5):: rho
    double precision, dimension(30,5) :: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision:: vapfrac
    double precision:: dens_residuum
    integer:: errval, nrofphases
    integer, dimension(5) :: phasetype
    end

    Double Precision module Function Density_dif(gl,press, Parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    Double Precision :: press             ! Inputvariables
    End Function


    module subroutine PhaseDet_ps(gl,press, Temp, x, rho, x_phase, phasetype, vapfrac, s_spec, nrofphases, errval)
    implicit none
    type(type_gl) :: gl
    double precision :: press, s_spec, vapfrac, temp
    double precision, dimension(30) :: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_phase
    double precision, dimension(30):: x_vap, x_liq2
    integer :: errval, nrofphases
    integer, dimension(5) :: phasetype
    end subroutine PhaseDet_ps

    Double Precision module Function Entropy_dif(gl,Temp, Parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    Double Precision :: Temp, press, s_spec  ! Inputvariables
    End Function Entropy_dif


    module subroutine PhaseDet_ph(gl,press, Temp, x, rho, x_phase, phasetype, vapfrac, h_spec, nrofphases, errval)
    implicit none
    type(type_gl) :: gl
    double precision :: press, h_spec, vapfrac, temp
    double precision, dimension(30) :: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_phase
    double precision, dimension(30):: x_vap, x_liq2
    integer :: errval, nrofphases
    integer, dimension(5) :: phasetype
    end subroutine PhaseDet_ph

    Double Precision module Function Enthalpy_dif(gl,Temp, Parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    Double Precision :: Temp, press, h_spec  ! Inputvariables
    End Function Enthalpy_dif

    end interface

    contains

    end module phasedet_mix_module
