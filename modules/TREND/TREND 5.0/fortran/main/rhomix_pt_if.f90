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

    ! module for file rhomix_pt.f90
    module rhomix_pt_module
    !global use inclusion
    use module_all_types
    use module_regula_falsi_support

    ! interface for routines with circular dependencies (implemented in submodule in file <this-file-name>_impl)
    interface

    double precision module function rhomix_calc(gl,TEMPERATURE, PRESSURE, RHO_EST_GIVEN, IPHASE, nrsubst)
    !DEC$ ATTRIBUTES DLLEXPORT :: rhomix_calc
    implicit none
    type(type_gl):: gl
    double precision:: TEMPERATURE, PRESSURE, RHO_EST_GIVEN
    integer:: IPHASE,nrsubst
    end function rhomix_calc

    Double Precision module Function dp_drho_zero(gl,rho, parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    Double Precision :: rho
    End Function dp_drho_zero

    module subroutine rho_test (gl,T, p, rho, IPhase, IFound, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, rho
    integer:: IPhase, nrsubst
    integer:: IFound
    endsubroutine rho_test

    Double Precision module Function rhomix_iter(gl,T, p, rho_min, rho_max, rho_min_allowed, rho_max_allowed,nrsubst)
    implicit none
    type(type_gl) :: gl
    Double Precision :: T, p, rho_min, rho_max, rho_min_allowed, rho_max_allowed! Inputvariablen
    integer :: nrsubst
    End Function rhomix_iter

    Double Precision module Function pressure_dif(gl,rho, Parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    Double Precision :: rho               ! Inputvariablen
    End Function pressure_dif
    !**************************************************************************

    end interface

    contains

    end module rhomix_pt_module
