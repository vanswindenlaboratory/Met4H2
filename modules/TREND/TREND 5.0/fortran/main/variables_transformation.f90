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
!          Jäger, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (201
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
    
module variables_transformation_module

use module_all_types

contains

subroutine add_T_conversion_derivs(gl,get)
! Henning Markgraf, July 2016
!DEC$ ATTRIBUTES REFERENCE :: add_T_conversion_derivs
    ! Adds the additional T derivatives to the getvector, which are
    ! required for the conversion to the tau derivatives
    
    
implicit none



    type(type_gl) :: gl

    
    ! I. Declarations
    ! input and output:
    integer, dimension(nderivs), intent (inout) :: get

    
    integer:: errorfld

    
    ! II. Add the additional T derivatives required for the conversion
    if (get(5) .EQ. 1) then 
        get(4) = 1
    end if
    
    if (get(7) .EQ. 1) then 
        get(6) = 1
    end if
    
    if (get(9) .EQ. 1) then 
        get(5) = 1
        get(4) = 1
    end if

    if (get(13) .EQ. 1) then 
        get(10) = 1
    end if

    if (get(14) .EQ. 1) then 
        get(7) = 1
        get(6) = 1
    end if

    if (get(15) .EQ. 1) then 
        get(9) = 1
        get(5) = 1
        get(4) = 1
    end if
    
end subroutine add_T_conversion_derivs
    
    
subroutine convert_T_derivs_Trho(gl,get,result_vector_T,result_vector_tau)
!DEC$ ATTRIBUTES REFERENCE :: convert_T_derivs_Trho
!INCOMPLETE
! Henning Markgraf, July 2016

    ! Converts derivatives of a/(RT) with respect to T to derivatives
    ! of alpha with respect to tau, according to Span.2000, page 36
    
    ! Input:
    !       get                 A getvector with dimensions (nderivs) to specify
    !                           which T derivatives have to be converted
    !       result_vector_T     A result vector with derivatives of dimensions (nderivs)
    ! Output:
    !       result_vector_tau   The same result vector with T derivatives converted
    !                           to tau derivatives
    
    
implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input:
    integer, dimension(nderivs), intent (in) :: get
    double precision, dimension(nderivs), intent (in) :: result_vector_T
    ! output:
    double precision, dimension(nderivs), intent (out) :: result_vector_tau

    integer:: errorfld

    
    ! Initializations
    result_vector_tau = 0.d0
    
    ! Copy result_vector_T into result_vector_tau
    result_vector_tau = result_vector_T
    
    ! Convert T derivatives:
    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
    if (get(4) .EQ. 1) then
        result_vector_tau(4) = - result_vector_T(4)
    end if
    
    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
    if (get(5) .EQ. 1) then
        result_vector_tau(5) = result_vector_T(5) + 2.d0 * result_vector_T(4)
    end if
    
    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
    if (get(6) .EQ. 1) then
        result_vector_tau(6) = - result_vector_T(6)
    end if
    
    ! 7: tau^2*delta*d^3(alpha)/(dtau^2*ddelta) = T^2*rho*d^3(a/RT)/(dT^2*drho) + T*rho*d^2(a/RT)/(dT*drho)
    if (get(7) .EQ. 1) then
        result_vector_tau(7) = result_vector_T(7) + 2.d0*result_vector_T(6)
    end if
    
    ! 9:
    if (get(9) .EQ. 1) then
        result_vector_tau(9) = - result_vector_T(9) - 6.d0*(result_vector_T(5)+result_vector_T(4))
    end if
    
    ! 10:
    if (get(10) .EQ. 1) then
        result_vector_tau(10) = - result_vector_T(10)
    end if
    
    ! 12:
    if (get(12) .EQ. 1) then
        result_vector_tau(12) = - result_vector_T(12)
    end if
    
    ! 13:
    if (get(13) .EQ. 1) then
        result_vector_tau(13) = result_vector_T(13) + 2.d0*result_vector_T(10)
    end if

    ! 14:
    if (get(14) .EQ. 1) then
        result_vector_tau(14) = - result_vector_T(14) - 6.d0*(result_vector_T(7)+result_vector_T(6))
    end if
    
    ! 15:
    if (get(15) .EQ. 1) then
        result_vector_tau(15) = result_vector_T(15) + 12.d0*result_vector_T(9) + 36.d0*result_vector_T(5) + 24.d0*result_vector_T(4)
    end if
    
    
    end subroutine convert_T_derivs_Trho
    
subroutine convert_T_derivs_x1(gl,get,result_vector_T,result_vector_tau)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "convert_T_derivs_x1" :: convert_T_derivs_x1
!INCOMPLETE
! Henning Markgraf, July 2016

    ! Converts derivatives of a/(RT) with respect to T to derivatives
    ! of alpha with respect to tau, according to Span.2000, page 36
    
    ! Input:
    !       get                 A getvector with dimensions (nderivs) to specify
    !                           which T derivatives have to be converted
    !       result_vector_T     A result vector with derivatives of dimensions (nx1derivs,ncomp)
    ! Output:
    !       result_vector_tau   The same result vector with T derivatives converted
    !                           to tau derivatives
    
    
implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input:
    integer, dimension(nderivs), intent (in) :: get
    double precision, dimension(nx1derivs,gl%ncomp), intent (in) :: result_vector_T
    ! output:
    double precision, dimension(nx1derivs,gl%ncomp), intent (out) :: result_vector_tau
    
    integer:: errorfld

    
    ! Initializations
    result_vector_tau = 0.d0
    
    ! Copy result_vector_T into result_vector_tau
    result_vector_tau = result_vector_T
    
    ! Convert T derivatives:
    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
    if (get(4) .EQ. 1) then
        result_vector_tau(4,1:gl%ncomp) = - result_vector_T(4,1:gl%ncomp)
    end if
    
    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
    if (get(5) .EQ. 1) then
        result_vector_tau(5,1:gl%ncomp) = result_vector_T(5,1:gl%ncomp) + 2.d0 * result_vector_T(4,1:gl%ncomp)
    end if
    
    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
    if (get(6) .EQ. 1) then
        result_vector_tau(6,1:gl%ncomp) = - result_vector_T(6,1:gl%ncomp)
    end if
    
    ! 7: tau^2*delta*d^3(alpha)/(dtau^2*ddelta) = T^2*rho*d^3(a/RT)/(dT^2*drho) + T*rho*d^2(a/RT)/(dT*drho)
    if (get(7) .EQ. 1) then
        result_vector_tau(7,1:gl%ncomp) = result_vector_T(7,1:gl%ncomp) + 2.d0*result_vector_T(6,1:gl%ncomp)
    end if
    
    ! 9:
    if (get(9) .EQ. 1) then
        result_vector_tau(9,1:gl%ncomp) = - result_vector_T(9,1:gl%ncomp) - 6.d0*(result_vector_T(5,1:gl%ncomp)+result_vector_T(4,1:gl%ncomp))
    end if
    
    ! 10:
    if (get(10) .EQ. 1) then
        result_vector_tau(10,1:gl%ncomp) = - result_vector_T(10,1:gl%ncomp)
    end if
    
    
    end subroutine convert_T_derivs_x1
    
subroutine convert_T_derivs_x2(gl,get,result_vector_T,result_vector_tau)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "convert_T_derivs_x2" :: convert_T_derivs_x2
! Henning Markgraf, July 2016

    ! Converts derivatives of a/(RT) with respect to T to derivatives
    ! of alpha with respect to tau, according to Span.2000, page 36
    
    ! Input:
    !       get                 A getvector with dimensions (nderivs) to specify
    !                           which T derivatives have to be converted
    !       result_vector_T     A result vector with derivatives of dimensions (nx2derivs,ncomp)
    ! Output:
    !       result_vector_tau   The same result vector with T derivatives converted
    !                           to tau derivatives
    
    
implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input:
    integer, dimension(nderivs), intent (in) :: get
    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
    ! output:
    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
    
    integer:: errorfld

    
    ! Initializations
    result_vector_tau = 0.d0
    
    ! Copy result_vector_T into result_vector_tau
    result_vector_tau = result_vector_T
    
    ! Convert T derivatives:
    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
    if (get(4) .EQ. 1) then
        result_vector_tau(4,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(4,1:gl%ncomp,1:gl%ncomp)
    end if
    
    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
    if (get(5) .EQ. 1) then
        result_vector_tau(5,1:gl%ncomp,1:gl%ncomp) = result_vector_T(5,1:gl%ncomp,1:gl%ncomp) + 2.d0 * result_vector_T(4,1:gl%ncomp,1:gl%ncomp)
    end if
    
    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
    if (get(6) .EQ. 1) then
        result_vector_tau(6,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(6,1:gl%ncomp,1:gl%ncomp)
    end if

    end subroutine convert_T_derivs_x2
    
    
subroutine convert_T_derivs_x3(gl,get,result_vector_T,result_vector_tau)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "convert_T_derivs_x3" :: convert_T_derivs_x3
! Henning Markgraf, July 2016

    ! Converts derivatives of a/(RT) with respect to T to derivatives
    ! of alpha with respect to tau, according to Span.2000, page 36
    
    ! Input:
    !       get                 A getvector with dimensions (nderivs) to specify
    !                           which T derivatives have to be converted
    !       result_vector_T     A result vector with derivatives of dimensions (nx3derivs,ncomp)
    ! Output:
    !       result_vector_tau   The same result vector with T derivatives converted
    !                           to tau derivatives
    
    
implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input:
    integer, dimension(nderivs), intent (in) :: get
    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
    ! output:
    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
    
    integer:: errorfld

    
    ! Initializations
    result_vector_tau = 0.d0
    
    ! Copy result_vector_T into result_vector_tau
    result_vector_tau = result_vector_T
    
    ! Convert T derivatives:
    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
    if (get(4) .EQ. 1) then
        result_vector_tau(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) ! careful again with the different numbering in the third comp. result vector, nr. 4 is position 3!
    end if

end subroutine convert_T_derivs_x3

end module variables_transformation_module