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

    ! module for file pc-saft_ancillary_routines.f90
    module pc_saft_ancillary_routines_module
    !global use inclusion
    use module_all_types
    use variables_transformation_module

    
! interface for routines with circular dependencies (implemented in submodule in file <this-file-name>_impl)    
    interface

    module subroutine calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_Trho

    module subroutine calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (inout) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_QQ_Trho

    module subroutine calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x1_QQ

    module subroutine calculate_PCSAFT_functionparts_x2_QQ(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x2_QQ

    module subroutine calculate_PCSAFT_functionparts_x2_DD(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x2_DD

    module subroutine calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x1_DD

    module subroutine calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (inout) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_DD_Trho

    module subroutine calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x1

    module subroutine calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x2

    module subroutine calculate_PCSAFT_functionparts_x3(gl,T,D,getprevious)
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    end subroutine calculate_PCSAFT_functionparts_x3

    


!    module subroutine add_T_conversion_derivs(gl,get)
!    implicit none
!    type(type_gl) :: gl
!    integer, dimension(nderivs), intent (inout) :: get
!    end subroutine add_T_conversion_derivs
!    
!    
!    module subroutine convert_T_derivs_Trho(gl,get,result_vector_T,result_vector_tau)
!    implicit none
!    type(type_gl) :: gl
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nderivs), intent (in) :: result_vector_T
!    double precision, dimension(nderivs), intent (out) :: result_vector_tau
!    end subroutine convert_T_derivs_Trho
!    
!    
!    module subroutine convert_T_derivs_x1(gl,get,result_vector_T,result_vector_tau)
!    implicit none
!    type(type_gl) :: gl
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx1derivs,gl%ncomp), intent (in) :: result_vector_T
!    double precision, dimension(nx1derivs,gl%ncomp), intent (out) :: result_vector_tau
!    end subroutine convert_T_derivs_x1
!    
!
!    module subroutine convert_T_derivs_x2(gl,get,result_vector_T,result_vector_tau)
!    implicit none
!    type(type_gl) :: gl
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
!    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
!    end subroutine convert_T_derivs_x2
!    
!    
!    module subroutine convert_T_derivs_x3(gl,get,result_vector_T,result_vector_tau)
!    implicit none
!    type(type_gl) :: gl
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
!    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
!    end subroutine convert_T_derivs_x3
    
    
    module subroutine add_previous_derivs(gl,getvector,getprevious)
    implicit none
    type(type_gl) :: gl
    integer, dimension(nderivs), intent (in) :: getvector
    integer, dimension(nderivs), intent (out) :: getprevious
    end subroutine
    
    
    module subroutine ABPARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABPARAMETERS
    
    
    module subroutine ABCQPARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABCQPARAMETERS


    module subroutine ABCDPARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABCDPARAMETERS

    
    module subroutine ABX1PARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABX1PARAMETERS
    
    
    module subroutine ABX2PARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABX2PARAMETERS
    

    module subroutine ABX3PARAMETERS(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine ABX3PARAMETERS

    
    module subroutine allocate_arrays_PCSAFT(gl,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer, intent (inout) :: nrsubst
    end subroutine allocate_arrays_PCSAFT


    module subroutine init_pure_mix(gl,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer :: nrsubst
    end subroutine init_pure_mix


    module subroutine init_mixture_PCSAFT(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine init_mixture_PCSAFT
  

    module subroutine init_derivs_PCSAFT(gl,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer, intent (inout) :: nrsubst
    end subroutine init_derivs_PCSAFT
    
    module subroutine check_calculated_PCSAFT(gl,T,D,nrsubst,calculated_derivs,getvector,already_calculated)
    implicit none
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: calculated_derivs
    integer, intent (in) :: nrsubst
    integer, dimension (nderivs), intent (inout) :: getvector
    logical, intent (out) :: already_calculated
    end subroutine check_calculated_PCSAFT

    
    module subroutine remove_redundant_getvector_entries(gl,calculated_derivs,getvector)
    implicit none
    type(type_gl) :: gl
    integer, dimension (nderivs), intent (in) :: calculated_derivs
    integer, dimension (nderivs), intent (inout) :: getvector
    end subroutine remove_redundant_getvector_entries

    
    module subroutine update_modulevariables_PCSAFT(gl,T,D,nrsubst,calculated_derivs,get)
    implicit none
    type(type_gl) :: gl
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: get
    integer, intent (in) :: nrsubst
    integer, dimension (nderivs), intent (inout) :: calculated_derivs
    end subroutine update_modulevariables_PCSAFT


    
    end interface
    
    
    contains



    end module pc_saft_ancillary_routines_module
