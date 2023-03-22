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

    ! module for file reduced_parameters_calc.f90
    module reduced_parameters_calc_module
    !global use inclusion
    use module_all_types

    interface

    module subroutine update_parameters(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: Temp, T
    integer :: nrsubst
    end subroutine update_parameters



    module subroutine reduced_parameters_calc(gl,Temp)
    !DEC$ ATTRIBUTES DLLEXPORT :: reduced_parameters_calc
    implicit none
    type(type_gl) :: gl
    double precision :: Temp
    end subroutine reduced_parameters_calc


    module subroutine ndYr_dni(gl,ndTred_dni, ndrhored_dni)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30) ::ndTred_dni, ndrhored_dni   ! Return vectors, maximum 30 components
    end subroutine ndYr_dni


    module subroutine dYr_dxi(gl,dTred_dxi, drhored_dxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::dTred_dxi, drhored_dxi
    end subroutine dYr_dxi


    module subroutine d2Yr_dxi2(gl,d2Tred_dxi2, d2rhored_dxi2)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::d2Tred_dxi2, d2rhored_dxi2   ! Return vectors, maximum 30 components
    end subroutine d2Yr_dxi2


    module subroutine d2Yr_dxidxj(gl,d2Tred_dxidxj, d2rhored_dxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30,30) :: d2Tred_dxidxj, d2rhored_dxidxj
    end subroutine d2Yr_dxidxj


    module subroutine dndYr_dnidxj(gl,dndTr_dnidxj, dndrhor_dnidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30)::dndTr_dnidxj, dndrhor_dnidxj   ! Return matrices, maximum 30 components
    end subroutine dndYr_dnidxj


    module subroutine d3Yr_dxidxjdxk(gl,d3Tred_dxidxjdxk, d3rhored_dxidxjdxk, d2Tred_dxidxj, d2rhored_dxidxj, dTred_dxi, drhored_dxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30,30,30) :: d3Tred_dxidxjdxk, d3rhored_dxidxjdxk
    double precision, dimension(30):: dTred_dxi
    double precision, dimension(30,30)::d2Tred_dxidxj
    double precision, dimension(30):: drhored_dxi
    double precision, dimension(30,30)::d2rhored_dxidxj
    end subroutine d3Yr_dxidxjdxk


    module subroutine d2ndYr_dnidxjdxk(gl,d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30, 30)::d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk   ! Return matrices, maximum 30 components
    end subroutine d2ndYr_dnidxjdxk


    module subroutine Help_PSI_Derivs(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine Help_PSI_Derivs


    module subroutine ndYr_dni_old(gl,ndTred_dni, ndrhored_dni)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::ndTred_dni, ndrhored_dni   ! Return vectors, maximum 30 components
    end subroutine ndYr_dni_old


    module subroutine dYr_dxi_old(gl,dTred_dxi, drhored_dxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::dTred_dxi, drhored_dxi
    end subroutine dYr_dxi_old


    module subroutine d2Yr_dxi2_old(gl,d2Tred_dxi2, d2rhored_dxi2)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::d2Tred_dxi2, d2rhored_dxi2   ! Return vectors, maximum 30 components
    end subroutine d2Yr_dxi2_old


    module subroutine d2Yr_dxidxj_old(gl,d2Tred_dxidxj, d2rhored_dxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30,30):: d2Tred_dxidxj, d2rhored_dxidxj
    end subroutine d2Yr_dxidxj_old


    module subroutine dndYr_dnidxj_old(gl,dndTr_dnidxj, dndrhor_dnidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30)::dndTr_dnidxj, dndrhor_dnidxj   ! Return matrices, maximum 30 components
    end subroutine dndYr_dnidxj_old


    end interface

    contains

    end module reduced_parameters_calc_module
