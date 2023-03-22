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
    SUBMODULE (pc_saft_ancillary_routines_module) impl
    !global use inclusion
    use module_all_types
    use pc_saft_module
    use pc_saft_A2X_derivs_module
    use pc_saft_A3X_derivs_module
    use pc_saft_ADISPX_derivs_module
    use pc_saft_AHCX_derivs_module
    use pc_saft_AHSX_derivs_module
    use pc_saft_ANCX_derivs_module
    use pc_saft_ARX_derivs_module
    use pc_saft_CX1_derivs_module
    use pc_saft_CX2_derivs_module
    use pc_saft_CX3_derivs_module
    use pc_saft_GIIX_derivs_module
    use pc_saft_IX_derivs_module
    use pc_saft_JX_derivs_module
    use pc_saft_MEOX_derivs_module


    contains




module subroutine calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)

! Henning Markgraf, June 2016

    ! calculates the PC SAFT function parts needed for the calculation of a_hc and a_disp, except for a_hs
    ! i.e.: d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1

    ! calculates the PC SAFT function parts needed for the calculation of a_qq 
    ! i.e.: J_2, J_3, zeta_n, a_ij, b_ij, c_ij 



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    ! output: module variables: di_PCSAFT, z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, mmean_PCSAFT,
    !                           ab_PCSAFT, gii_PCSAFT, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, C_PCSAFT

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Calculate the function parts d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1
    ! calculate d_i
    call DIDERIVS(gl,T,getprevious)
    ! calculate zeta_n
    call ZETADERIVS(gl,D,getprevious,0)
    call ZETADERIVS(gl,D,getprevious,1)
    call ZETADERIVS(gl,D,getprevious,2)
    call ZETADERIVS(gl,D,getprevious,3)
        
    ! mmean and a_ii, b_ii which only depend on the composition
    ! are already calculated in init_derivs_PCSAFT before
    
    ! VIII. calculate g_ii, I_1, I_2, meo1, meo2, C and J_2, J_3 (for quardrupolar and polar molecules accordingly) which depend on the previously calculated functions
    call GIIDERIVS(gl,getprevious)
    call IDERIVS(gl,getprevious,1)
    call IDERIVS(gl,getprevious,2)
    call MEODERIVS(gl,T,getprevious,1)
    call MEODERIVS(gl,T,getprevious,2)
    call CDERIVS(gl,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_Trho


module subroutine calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (inout) :: getprevious
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    call ABCQPARAMETERS(gl)
    
    call J2DERIVS(gl,T,D,getprevious)
    call J3DERIVS(gl,T,D,getprevious)
    call A2DERIVS(gl,T,D,getprevious)
    call A3DERIVS(gl,T,D,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_QQ_Trho

module subroutine calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call J2X1DERIVS(gl,T,D,getprevious)
    call J3X1DERIVS(gl,T,D,getprevious)
    call A2X1DERIVS(gl,T,D,getprevious)
    call A3X1DERIVS(gl,T,D,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_x1_QQ



module subroutine calculate_PCSAFT_functionparts_x2_QQ(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call J2X2DERIVS(gl,T,D,getprevious)
    call J3X2DERIVS(gl,T,D,getprevious)
    call A2X2DERIVS(gl,T,D,getprevious)
    call A3X2DERIVS(gl,T,D,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_x2_QQ

module subroutine calculate_PCSAFT_functionparts_x2_DD(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call J2X2DERIVS_D(gl,T,D,getprevious)
    call J3X2DERIVS_D(gl,T,D,getprevious)
    call A2X2DERIVS_D(gl,T,D,getprevious)
    call A3X2DERIVS_D(gl,T,D,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_x2_DD

module subroutine calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call J2X1DERIVS_D(gl,T,D,getprevious)
    call J3X1DERIVS_D(gl,T,D,getprevious)
    call A2X1DERIVS_D(gl,T,D,getprevious)
    call A3X1DERIVS_D(gl,T,D,getprevious)
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_x1_DD


module subroutine calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)



    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (inout) :: getprevious

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call ABCDPARAMETERS(gl)

    call J2DERIVS_D(gl,T,D,getprevious)
    call J3DERIVS_D(gl,T,D,getprevious)
    call A2DERIVS_D(gl,T,D,getprevious)
    call A3DERIVS_D(gl,T,D,getprevious)    
    
!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_DD_Trho

    
module subroutine calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)

! Henning Markgraf, June 2016

    ! calculates the PC SAFT function parts needed for the calculation of ahc_x and adisp_x, except for ahs_x
    ! i.e.: di_x, zetan_x, ab_x, gii_x, I_x, meo_x, C_x




    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    ! output: module variables: z0x1_PCSAFT, z1x1_PCSAFT, z2x1_PCSAFT, z3x1_PCSAFT, abx1_PCSAFT, giix1_PCSAFT, 
    !                           i1x1_PCSAFT, i2x1_PCSAFT, meo1x1_PCSAFT, meo2x1_PCSAFT, cx1_PCSAFT

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! Calculation of zeta_n_x, ahs_x, I_x, meo_x, C_x
    ! calculate zeta_n_x
    call ZETAX1DERIVS(gl,D,getprevious,0)
    call ZETAX1DERIVS(gl,D,getprevious,1)
    call ZETAX1DERIVS(gl,D,getprevious,2)
    call ZETAX1DERIVS(gl,D,getprevious,3)
    
    ! mmean and a_ii, b_ii which only depend on the composition
    ! are already calculated in init_derivs_PCSAFT before
    
    ! calculate ahs_x, gii_x, I1_x, I2_x, meo1_x, meo2_x, C_x
    call GIIX1DERIVS(gl,getprevious)
    call IX1DERIVS(gl,getprevious,1)
    call IX1DERIVS(gl,getprevious,2)
    call MEOX1DERIVS(gl,T,getprevious,1)
    call MEOX1DERIVS(gl,T,getprevious,2)
    call CX1DERIVS(gl,getprevious)
    
!DEC$ END IF
    
    end subroutine calculate_PCSAFT_functionparts_x1

module subroutine calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)

! Henning Markgraf, June 2016

    ! calculates the PC SAFT function parts needed for the calculation of ahc_xx and adisp_xx, except for ahs_xx
    ! i.e.: ab_xx, gii_xx, I_xx, meo_xx, C_xx




    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    ! output: module variables: abx2_PCSAFT, giix2_PCSAFT, 
    !                           i1x2_PCSAFT, i2x2_PCSAFT, meo1x2_PCSAFT, meo2x2_PCSAFT, cx2_PCSAFT
    ! mmean and a_ii, b_ii which only depend on the composition
    ! are already calculated in init_derivs_PCSAFT before

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call GIIX2DERIVS(gl,getprevious)
    call IX2DERIVS(gl,getprevious,1)
    call IX2DERIVS(gl,getprevious,2)
    call MEOX2DERIVS(gl,T,getprevious,1)
    call MEOX2DERIVS(gl,T,getprevious,2)
    call CX2DERIVS(gl,getprevious)

!DEC$ END IF
    
    end subroutine calculate_PCSAFT_functionparts_x2
    
module subroutine calculate_PCSAFT_functionparts_x3(gl,T,D,getprevious)

! Henning Markgraf, June 2016

    ! calculates the PC SAFT function parts needed for the calculation of ahc_xxx and adisp_xxx, except for ahs_xxx
    ! i.e.: ab_xxx, gii_xxx, I_xxx, C_xxx




    

implicit none

    type(type_gl) :: gl

    
    !Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: getprevious
    ! output: module variables: abx3_PCSAFT, giix3_PCSAFT, 
    !                           i1x3_PCSAFT, i2x3_PCSAFT, cx3_PCSAFT
    ! mmean and a_ii, b_ii which only depend on the composition
    ! are already calculated in init_derivs_PCSAFT before

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    call GIIX3DERIVS(gl,getprevious)
    call IX3DERIVS(gl,getprevious,1)
    call IX3DERIVS(gl,getprevious,2)
    call CX3DERIVS(gl,getprevious)
    !call J2X3DERIVS(T,D,getprevious)
    !call J3X3DERIVS(D,getprevious)
    !call J2X3DERIVS_D(T,D,getprevious)
    !call J3X3DERIVS_D(D,getprevious)

!DEC$ END IF
    
end subroutine calculate_PCSAFT_functionparts_x3



!module subroutine add_T_conversion_derivs(gl,get)
!! Henning Markgraf, July 2016
!
!    ! Adds the additional T derivatives to the getvector, which are
!    ! required for the conversion to the tau derivatives
!    
!    
!implicit none
!
!
!
!    type(type_gl) :: gl
!
!    
!    ! I. Declarations
!    ! input and output:
!    integer, dimension(nderivs), intent (inout) :: get
!
!    
!    integer:: errorfld
!
!!DEC$ IF DEFINED(WO_PCSAFT)
!    errorfld = -7878
!!DEC$ ELSE
!    
!    ! II. Add the additional T derivatives required for the conversion
!    if (get(5) .EQ. 1) then 
!        get(4) = 1
!    end if
!    
!    if (get(7) .EQ. 1) then 
!        get(6) = 1
!    end if
!    
!    if (get(9) .EQ. 1) then 
!        get(5) = 1
!        get(4) = 1
!    end if
!
!    if (get(13) .EQ. 1) then 
!        get(10) = 1
!    end if
!
!    if (get(14) .EQ. 1) then 
!        get(7) = 1
!        get(6) = 1
!    end if
!
!    if (get(15) .EQ. 1) then 
!        get(9) = 1
!        get(5) = 1
!        get(4) = 1
!    end if
!!DEC$ END IF
!    
!end subroutine add_T_conversion_derivs
!    
!    
!module subroutine convert_T_derivs_Trho(gl,get,result_vector_T,result_vector_tau)
!!INCOMPLETE
!! Henning Markgraf, July 2016
!
!    ! Converts derivatives of a/(RT) with respect to T to derivatives
!    ! of alpha with respect to tau, according to Span.2000, page 36
!    
!    ! Input:
!    !       get                 A getvector with dimensions (nderivs) to specify
!    !                           which T derivatives have to be converted
!    !       result_vector_T     A result vector with derivatives of dimensions (nderivs)
!    ! Output:
!    !       result_vector_tau   The same result vector with T derivatives converted
!    !                           to tau derivatives
!    
!    
!implicit none
!
!    type(type_gl) :: gl
!
!    
!    ! I. Declarations
!    ! input:
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nderivs), intent (in) :: result_vector_T
!    ! output:
!    double precision, dimension(nderivs), intent (out) :: result_vector_tau
!
!    integer:: errorfld
!
!!DEC$ IF DEFINED(WO_PCSAFT)
!    errorfld = -7878
!!DEC$ ELSE    
!    
!    ! Initializations
!    result_vector_tau = 0.d0
!    
!    ! Copy result_vector_T into result_vector_tau
!    result_vector_tau = result_vector_T
!    
!    ! Convert T derivatives:
!    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
!    if (get(4) .EQ. 1) then
!        result_vector_tau(4) = - result_vector_T(4)
!    end if
!    
!    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
!    if (get(5) .EQ. 1) then
!        result_vector_tau(5) = result_vector_T(5) + 2.d0 * result_vector_T(4)
!    end if
!    
!    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
!    if (get(6) .EQ. 1) then
!        result_vector_tau(6) = - result_vector_T(6)
!    end if
!    
!    ! 7: tau^2*delta*d^3(alpha)/(dtau^2*ddelta) = T^2*rho*d^3(a/RT)/(dT^2*drho) + T*rho*d^2(a/RT)/(dT*drho)
!    if (get(7) .EQ. 1) then
!        result_vector_tau(7) = result_vector_T(7) + 2.d0*result_vector_T(6)
!    end if
!    
!    ! 9:
!    if (get(9) .EQ. 1) then
!        result_vector_tau(9) = - result_vector_T(9) - 6.d0*(result_vector_T(5)+result_vector_T(4))
!    end if
!    
!    ! 10:
!    if (get(10) .EQ. 1) then
!        result_vector_tau(10) = - result_vector_T(10)
!    end if
!    
!    ! 12:
!    if (get(12) .EQ. 1) then
!        result_vector_tau(12) = - result_vector_T(12)
!    end if
!    
!    ! 13:
!    if (get(13) .EQ. 1) then
!        result_vector_tau(13) = result_vector_T(13) + 2.d0*result_vector_T(10)
!    end if
!
!    ! 14:
!    if (get(14) .EQ. 1) then
!        result_vector_tau(14) = - result_vector_T(14) - 6.d0*(result_vector_T(7)+result_vector_T(6))
!    end if
!    
!    ! 15:
!    if (get(15) .EQ. 1) then
!        result_vector_tau(15) = result_vector_T(15) + 12.d0*result_vector_T(9) + 36.d0*result_vector_T(5) + 24.d0*result_vector_T(4)
!    end if
!    
!!DEC$ END IF
!    
!    end subroutine convert_T_derivs_Trho
!    
!module subroutine convert_T_derivs_x1(gl,get,result_vector_T,result_vector_tau)
!!INCOMPLETE
!! Henning Markgraf, July 2016
!
!    ! Converts derivatives of a/(RT) with respect to T to derivatives
!    ! of alpha with respect to tau, according to Span.2000, page 36
!    
!    ! Input:
!    !       get                 A getvector with dimensions (nderivs) to specify
!    !                           which T derivatives have to be converted
!    !       result_vector_T     A result vector with derivatives of dimensions (nx1derivs,ncomp)
!    ! Output:
!    !       result_vector_tau   The same result vector with T derivatives converted
!    !                           to tau derivatives
!    
!    
!implicit none
!
!    type(type_gl) :: gl
!
!    
!    ! I. Declarations
!    ! input:
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx1derivs,gl%ncomp), intent (in) :: result_vector_T
!    ! output:
!    double precision, dimension(nx1derivs,gl%ncomp), intent (out) :: result_vector_tau
!    
!    integer:: errorfld
!
!!DEC$ IF DEFINED(WO_PCSAFT)
!    errorfld = -7878
!!DEC$ ELSE
!    
!    ! Initializations
!    result_vector_tau = 0.d0
!    
!    ! Copy result_vector_T into result_vector_tau
!    result_vector_tau = result_vector_T
!    
!    ! Convert T derivatives:
!    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
!    if (get(4) .EQ. 1) then
!        result_vector_tau(4,1:gl%ncomp) = - result_vector_T(4,1:gl%ncomp)
!    end if
!    
!    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
!    if (get(5) .EQ. 1) then
!        result_vector_tau(5,1:gl%ncomp) = result_vector_T(5,1:gl%ncomp) + 2.d0 * result_vector_T(4,1:gl%ncomp)
!    end if
!    
!    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
!    if (get(6) .EQ. 1) then
!        result_vector_tau(6,1:gl%ncomp) = - result_vector_T(6,1:gl%ncomp)
!    end if
!    
!    ! 7: tau^2*delta*d^3(alpha)/(dtau^2*ddelta) = T^2*rho*d^3(a/RT)/(dT^2*drho) + T*rho*d^2(a/RT)/(dT*drho)
!    if (get(7) .EQ. 1) then
!        result_vector_tau(7,1:gl%ncomp) = result_vector_T(7,1:gl%ncomp) + 2.d0*result_vector_T(6,1:gl%ncomp)
!    end if
!    
!    ! 9:
!    if (get(9) .EQ. 1) then
!        result_vector_tau(9,1:gl%ncomp) = - result_vector_T(9,1:gl%ncomp) - 6.d0*(result_vector_T(5,1:gl%ncomp)+result_vector_T(4,1:gl%ncomp))
!    end if
!    
!    ! 10:
!    if (get(10) .EQ. 1) then
!        result_vector_tau(10,1:gl%ncomp) = - result_vector_T(10,1:gl%ncomp)
!    end if
!    
!!DEC$ END IF
!    
!    end subroutine convert_T_derivs_x1
!    
!module subroutine convert_T_derivs_x2(gl,get,result_vector_T,result_vector_tau)
!
!! Henning Markgraf, July 2016
!
!    ! Converts derivatives of a/(RT) with respect to T to derivatives
!    ! of alpha with respect to tau, according to Span.2000, page 36
!    
!    ! Input:
!    !       get                 A getvector with dimensions (nderivs) to specify
!    !                           which T derivatives have to be converted
!    !       result_vector_T     A result vector with derivatives of dimensions (nx2derivs,ncomp)
!    ! Output:
!    !       result_vector_tau   The same result vector with T derivatives converted
!    !                           to tau derivatives
!    
!    
!implicit none
!
!    type(type_gl) :: gl
!
!    
!    ! I. Declarations
!    ! input:
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
!    ! output:
!    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
!    
!    integer:: errorfld
!
!!DEC$ IF DEFINED(WO_PCSAFT)
!    errorfld = -7878
!!DEC$ ELSE
!    
!    ! Initializations
!    result_vector_tau = 0.d0
!    
!    ! Copy result_vector_T into result_vector_tau
!    result_vector_tau = result_vector_T
!    
!    ! Convert T derivatives:
!    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
!    if (get(4) .EQ. 1) then
!        result_vector_tau(4,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(4,1:gl%ncomp,1:gl%ncomp)
!    end if
!    
!    ! 5: tau^2*d^2(alpha)/dtau^2 = T^2*d^2(a/RT)/dT^2 + 2 * T*d(a/RT)/dT
!    if (get(5) .EQ. 1) then
!        result_vector_tau(5,1:gl%ncomp,1:gl%ncomp) = result_vector_T(5,1:gl%ncomp,1:gl%ncomp) + 2.d0 * result_vector_T(4,1:gl%ncomp,1:gl%ncomp)
!    end if
!    
!    ! 6: tau*delta*d^2(alpha)/(dtau*ddelta) = - T*rho*d^2(a/RT)/(dT*drho)
!    if (get(6) .EQ. 1) then
!        result_vector_tau(6,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(6,1:gl%ncomp,1:gl%ncomp)
!    end if
!
!!DEC$ END IF
!    
!    end subroutine convert_T_derivs_x2
!    
!    
!module subroutine convert_T_derivs_x3(gl,get,result_vector_T,result_vector_tau)
!
!! Henning Markgraf, July 2016
!
!    ! Converts derivatives of a/(RT) with respect to T to derivatives
!    ! of alpha with respect to tau, according to Span.2000, page 36
!    
!    ! Input:
!    !       get                 A getvector with dimensions (nderivs) to specify
!    !                           which T derivatives have to be converted
!    !       result_vector_T     A result vector with derivatives of dimensions (nx3derivs,ncomp)
!    ! Output:
!    !       result_vector_tau   The same result vector with T derivatives converted
!    !                           to tau derivatives
!    
!    
!implicit none
!
!    type(type_gl) :: gl
!
!    
!    ! I. Declarations
!    ! input:
!    integer, dimension(nderivs), intent (in) :: get
!    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (in) :: result_vector_T
!    ! output:
!    double precision, dimension(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (out) :: result_vector_tau
!    
!    integer:: errorfld
!
!!DEC$ IF DEFINED(WO_PCSAFT)
!    errorfld = -7878
!!DEC$ ELSE
!    
!    ! Initializations
!    result_vector_tau = 0.d0
!    
!    ! Copy result_vector_T into result_vector_tau
!    result_vector_tau = result_vector_T
!    
!    ! Convert T derivatives:
!    ! 4: tau*dalpha/dtau = - T*d(a/RT)/dT
!    if (get(4) .EQ. 1) then
!        result_vector_tau(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) = - result_vector_T(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) ! careful again with the different numbering in the third comp. result vector, nr. 4 is position 3!
!    end if
!
!!DEC$ END IF
!    
!end subroutine convert_T_derivs_x3
    
    
module subroutine add_previous_derivs(gl,getvector,getprevious)
    ! returns a new getvector which additionally has a 1 at the spots of the previous (lower order) derivatives
    ! ( e.g. previous derivative of (d^2f)/(dTdD) are (df)/(dT) and (df)/(dD) )
    
    ! Input:    getvector           - array of dimension (nderivs) filled with 0 or 1
    ! Output:   getprevious   - array of dimension (nderivs) filled with 0 or 1
    !                                   new getvector array, but derivatives which are prior derivatives to the one in the input getvector
    !                                   also have a 1 to be calculated
    !                                   Example: getvector only has a "1" at element 5, then the getprevious will have a "1" at element
    !                                              4 and 5 because, nr. 5 represents the derivative wrt tau^2 and nr. 4 is the one wrt tau
    
    
implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    integer, dimension(nderivs), intent (in) :: getvector
    integer, dimension(nderivs), intent (out) :: getprevious
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Initialization
    getprevious = getvector
    
    ! III. determine which previous derivatives are required
    getprevious(1) = 1 ! the basic equation is always needed
    ! check delta derivatives:
    if (getprevious(11) .EQ. 1) then ! 4 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
        getprevious(8) = 1
        getprevious(11) = 1
    elseif (getprevious(8) .EQ. 1) then ! 3 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
        getprevious(8) = 1
    elseif (getprevious(3) .EQ. 1) then ! 2 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
    end if
    ! pure tau derivatives:
    if (getprevious(15) .EQ. 1) then ! 4 tau derivatives
        getprevious(4) = 1
        getprevious(5) = 1
        getprevious(9) = 1
        getprevious(15) = 1
    elseif (getprevious(9) .EQ. 1) then ! 3 tau derivatives
        getprevious(4) = 1
        getprevious(5) = 1
        getprevious(9) = 1
    elseif (getprevious(5) .EQ. 1) then ! 2 tau derivatives
        getprevious(4) = 1
        getprevious(5) = 1
    end if
    ! mixed tau, delta derivatives
    if (getprevious(14) .EQ. 1) then ! 3 tau, 1 delta derivatives
        getprevious(2) = 1
        getprevious(4) = 1
        getprevious(5) = 1
        getprevious(6) = 1
        getprevious(7) = 1
        getprevious(9) = 1
        getprevious(14) = 1
    end if
    if (getprevious(13) .EQ. 1) then ! 2 tau, 2 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
        getprevious(4) = 1
        getprevious(5) = 1
        getprevious(6) = 1
        getprevious(7) = 1
        getprevious(10) = 1
        getprevious(13) = 1
    elseif (getprevious(7) .EQ. 1) then ! 2 tau, 1 delta derivatives
        getprevious(2) = 1
        getprevious(4) = 1
        getprevious(5) = 1
        getprevious(6) = 1
        getprevious(7) = 1
    elseif (getprevious(6) .EQ. 1) then ! 1 tau, 1 delta derivatives
        getprevious(2) = 1
        getprevious(4) = 1
        getprevious(6) = 1
    end if
    if (getprevious(12) .EQ. 1) then ! 1 tau, 3 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
        getprevious(4) = 1
        getprevious(6) = 1
        getprevious(8) = 1
        getprevious(10) = 1
        getprevious(12) = 1
    elseif (getprevious(10) .EQ. 1) then ! 1 tau, 2 delta derivatives
        getprevious(2) = 1
        getprevious(3) = 1
        getprevious(4) = 1
        getprevious(6) = 1
        getprevious(10) = 1
    end if

!DEC$ END IF
    
    end subroutine
    
module subroutine ABPARAMETERS(gl)

! Henning Markgraf, June 2016

    ! a_ii and b_ii: parameters in the PC SAFT equation
    ! defined by eq. A.18 and A.19 in Gross, Sadowski 2001:
    ! dependent on mmean

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    ! output: ab_PCSAFT (gl,module variable)
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. calculate a_ii
    do i = 0, 6
        gl%ab_PCSAFT(1,i) = a0PCSAFT(i) + (gl%mmean_PCSAFT-1.d0)/gl%mmean_PCSAFT*a1PCSAFT(i) + (gl%mmean_PCSAFT-1.d0)/gl%mmean_PCSAFT*(gl%mmean_PCSAFT-2.d0)/gl%mmean_PCSAFT*a2PCSAFT(i)
    end do
    
    ! III. calculate b_ii
    do i = 0, 6
        gl%ab_PCSAFT(2,i) = b0PCSAFT(i) + (gl%mmean_PCSAFT-1.d0)/gl%mmean_PCSAFT*b1PCSAFT(i) + (gl%mmean_PCSAFT-1.d0)/gl%mmean_PCSAFT*(gl%mmean_PCSAFT-2.d0)/gl%mmean_PCSAFT*b2PCSAFT(i)
    end do
    
!DEC$ END IF
    
    end subroutine ABPARAMETERS
    
module subroutine ABCQPARAMETERS(gl)

! T. Eckermann, 02/2017

    ! a_ij, b_ij, and c_ijk: parameters in the PC SAFT equation
    ! defined by equation 15, 16 and 17 in Gross 2005
    ! a_ij and b_ij dependent on mmeanQij_PCSAFT, c_ijk dependent on mmeanQij_PCSAFT and mmeanQijk_PCSAFT

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    ! output: abc_PCSAFTQ (gl,module variable)
    integer :: n,i,j,k
    double precision, dimension(30) :: mi
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! T. Eckermann, 02/2017
    gl%mmeanQij_PCSAFT = 0.d0
    gl%mmeanQijk_PCSAFT = 0.d0
    mi = 0.d0
    
    do i = gl%n_start, gl%n_end
        if (gl%mPCSAFT(i) > 2.d0) then
            mi(i) = 2.d0
        else
            mi(i) = gl%mPCSAFT(i)
        end if
    end do
    
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            gl%mmeanQij_PCSAFT(i,j) = (mi(i)*mi(j))**0.5d0
            do k = gl%n_start, gl%n_end
                gl%mmeanQijk_PCSAFT(i,j,k) = (mi(i)*mi(j)*mi(k))**(1.d0/3.d0)
            end do
        end do
    end do
    
    
    ! II. calculate a_ij ( eq. 15, p. 2560)
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            do n = 0, 4
                gl%ab_PCSAFTQ(1,n,i,j) = a0PCSAFTQ(n) + (gl%mmeanQij_PCSAFT(i,j)-1.d0)/gl%mmeanQij_PCSAFT(i,j)*a1PCSAFTQ(n) &
                    & + (gl%mmeanQij_PCSAFT(i,j)-1.d0)/gl%mmeanQij_PCSAFT(i,j)*(gl%mmeanQij_PCSAFT(i,j)-2.d0)/gl%mmeanQij_PCSAFT(i,j)*a2PCSAFTQ(n)
            end do
        end do
    end do
    
    ! III. calculate b_ij ( eq. 16, p. 2560)
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            do n = 0, 4
                gl%ab_PCSAFTQ(2,n,i,j) = b0PCSAFTQ(n) + (gl%mmeanQij_PCSAFT(i,j)-1.d0)/gl%mmeanQij_PCSAFT(i,j)*b1PCSAFTQ(n) &
                    & + (gl%mmeanQij_PCSAFT(i,j)-1.d0)/gl%mmeanQij_PCSAFT(i,j)*(gl%mmeanQij_PCSAFT(i,j)-2.d0)/gl%mmeanQij_PCSAFT(i,j)*b2PCSAFTQ(n)
            end do
        end do
    end do
    
    ! III. calculate c_ijk ( eq. 17, p. 2560)
    do k = gl%n_start, gl%n_end
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do n = 0, 4
                    gl%c_PCSAFTQ(n,i,j,k) = c0PCSAFTQ(n) + (gl%mmeanQijk_PCSAFT(i,j,k)-1.d0)/gl%mmeanQijk_PCSAFT(i,j,k)*c1PCSAFTQ(n) &
                        & + (gl%mmeanQijk_PCSAFT(i,j,k)-1.d0)/gl%mmeanQijk_PCSAFT(i,j,k)*(gl%mmeanQijk_PCSAFT(i,j,k)-2.d0)/gl%mmeanQijk_PCSAFT(i,j,k)*c2PCSAFTQ(n)
                end do
            end do
        end do
    end do
    
!DEC$ END IF
    
end subroutine ABCQPARAMETERS


module subroutine ABCDPARAMETERS(gl)

! I. Schuelling 07/2017

    ! a_ij, b_ij, and c_ijk: parameters in the PC SAFT equation
    ! defined by equation 12, 13 and 14 in Gross 2006
    ! a_ij and b_ij dependent on mmeanDij_PCSAFT, c_ijk dependent on mmeanDij_PCSAFT and mmeanDijk_PCSAFT

   
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    ! output: abc_PCSAFTD (gl,module variable)
    integer :: n,i,j,k
    
    integer:: errorfld
    double precision, dimension(30) :: mi

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! I. Schuelling 07/2017
    ! IMPORTANT:  the restriction mij<=2 and mijk<=2 is still missing here
    gl%mmeanDij_PCSAFT = 0.d0
    gl%mmeanDijk_PCSAFT = 0.d0
    mi = 0.d0
    
    do i = gl%n_start, gl%n_end
        if (gl%mPCSAFT(i) > 2.d0) then
            mi(i) = 2.d0
        else
            mi(i) = gl%mPCSAFT(i)
        end if
    end do
    
    
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
        
            gl%mmeanDij_PCSAFT(i,j) = (mi(i)*mi(j))**0.5d0
           
            do k = gl%n_start, gl%n_end
                gl%mmeanDijk_PCSAFT(i,j,k) = (mi(i)*mi(j)*mi(k))**(1.d0/3.d0) 
             !   if (gl%mmeanDijk_PCSAFT(i,j,k) > 2.d0) gl%mmeanDijk_PCSAFT(i,j,k) = 2.d0
            end do
        end do
    end do
    
    
    ! II. calculate a_ij ( eq. 12, p. 1196)
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            do n = 0, 4
                gl%ab_PCSAFTD(1,n,i,j) = a0PCSAFTD(n) + (gl%mmeanDij_PCSAFT(i,j)-1.d0)/gl%mmeanDij_PCSAFT(i,j)*a1PCSAFTD(n) &
                    & + (gl%mmeanDij_PCSAFT(i,j)-1.d0)/gl%mmeanDij_PCSAFT(i,j)*(gl%mmeanDij_PCSAFT(i,j)-2.d0)/gl%mmeanDij_PCSAFT(i,j)*a2PCSAFTD(n)
            end do
        end do
    end do
    
    ! III. calculate b_ij ( eq. 13, p. 1196)
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            do n = 0, 4
                gl%ab_PCSAFTD(2,n,i,j) = b0PCSAFTD(n) + (gl%mmeanDij_PCSAFT(i,j)-1.d0)/gl%mmeanDij_PCSAFT(i,j)*b1PCSAFTD(n) &
                    & + (gl%mmeanDij_PCSAFT(i,j)-1.d0)/gl%mmeanDij_PCSAFT(i,j)*(gl%mmeanDij_PCSAFT(i,j)-2.d0)/gl%mmeanDij_PCSAFT(i,j)*b2PCSAFTD(n)
            end do
        end do
    end do
    
    ! III. calculate c_ijk ( eq. 14, p. 1196)
    do k = gl%n_start, gl%n_end
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do n = 0, 4
                    gl%c_PCSAFTD(n,i,j,k) = c0PCSAFTD(n) + (gl%mmeanDijk_PCSAFT(i,j,k)-1.d0)/gl%mmeanDijk_PCSAFT(i,j,k)*c1PCSAFTD(n) &
                        & + (gl%mmeanDijk_PCSAFT(i,j,k)-1.d0)/gl%mmeanDijk_PCSAFT(i,j,k)*(gl%mmeanDijk_PCSAFT(i,j,k)-2.d0)/gl%mmeanDijk_PCSAFT(i,j,k)*c2PCSAFTD(n)
                end do
            end do
        end do
    end do

!DEC$ END IF
    
end subroutine ABCDPARAMETERS

    
module subroutine ABX1PARAMETERS(gl)

! Henning Markgraf, June 2016

    ! a_ii and b_ii: parameters in the PC SAFT equation
    ! defined by eq. A.18 and A.19 in Gross, Sadowski 2001:
    ! dependent on mmean

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !output: abx1_PCSAFT (gl,module variable)
    !working variables
    double precision :: factor_x
    integer :: i, xi
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. calculate a_ii
    do i = 0, 6
	    ! factor to save computation time in the loops
	    factor_x = (a1PCSAFT(i)/gl%mmean_PCSAFT**2 + 3.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**2 - 4.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**3)
        do xi = gl%n_start, gl%n_end
            gl%abx1_PCSAFT(1,i,xi) = gl%mPCSAFT(xi)*factor_x
        end do
    end do
    
    ! III. calculate b_ii
    do i = 0, 6
        ! factor to save computation time in the loops
	    factor_x = (b1PCSAFT(i)/gl%mmean_PCSAFT**2 + 3.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**2 - 4.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**3)
        do xi = gl%n_start, gl%n_end
            gl%abx1_PCSAFT(2,i,xi) = gl%mPCSAFT(xi)*factor_x
        end do
    end do
    
!DEC$ END IF
    
    end subroutine ABX1PARAMETERS
    
    
module subroutine ABX2PARAMETERS(gl)

! Henning Markgraf, June 2016

    ! a_ii and b_ii: parameters in the PC SAFT equation
    ! defined by eq. A.18 and A.19 in Gross, Sadowski 2001:
    ! dependent on mmean

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !output: abx2_PCSAFT
    !working variables
    double precision :: factor_xx
    integer :: i, xi, xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. calculate a_ii
    do i = 0, 6
	    ! factor to save computation time in the loops
	    factor_xx = (-2.d0*a1PCSAFT(i)/gl%mmean_PCSAFT**3 - 6.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**3 + 12.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**4)
        do xj = gl%n_start, gl%n_end
            do xi = gl%n_start, gl%n_end!xj
                gl%abx2_PCSAFT(1,i,xi,xj) = gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*factor_xx
		    end do
        end do
    end do
    
    ! III. calculate b_ii
    do i = 0, 6
        ! factor to save computation time in the loops
	    factor_xx = (-2.d0*b1PCSAFT(i)/gl%mmean_PCSAFT**3 - 6.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**3 + 12.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**4)
        do xj = gl%n_start, gl%n_end
            do xi = gl%n_start, gl%n_end
		        gl%abx2_PCSAFT(2,i,xi,xj) = gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*factor_xx
		    end do
        end do
    end do
    
!DEC$ END IF
    
    end subroutine ABX2PARAMETERS
    
module subroutine ABX3PARAMETERS(gl)

! Henning Markgraf, June 2016

    ! a_ii and b_ii: parameters in the PC SAFT equation
    ! defined by eq. A.18 and A.19 in Gross, Sadowski 2001:
    ! dependent on mmean

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !output: abx3_PCSAFT
    !working variables
    double precision :: factor_xxx
    integer :: i, xi, xj, xk
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. calculate a_ii
    do i = 0, 6
	    ! factor to save computation time in the loops
	    factor_xxx = (6.d0*a1PCSAFT(i)/gl%mmean_PCSAFT**4 + 18.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**4 - 48.d0*a2PCSAFT(i)/gl%mmean_PCSAFT**5)
        do xk = gl%n_start, gl%n_end
            do xj = gl%n_start, xk 
            	do xi = gl%n_start, xj
			        gl%abx3_PCSAFT(1,i,xi,xj,xk) = gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*factor_xxx
		        end do
		    end do
        end do
    end do
    
    ! III. calculate b_ii
    do i = 0, 6
	    ! factor to save computation time in the loops
	    factor_xxx = (6.d0*b1PCSAFT(i)/gl%mmean_PCSAFT**4 + 18.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**4 - 48.d0*b2PCSAFT(i)/gl%mmean_PCSAFT**5)
        do xk = gl%n_start, gl%n_end
            do xj = gl%n_start, xk 
            	do xi = gl%n_start, xj
			        gl%abx3_PCSAFT(2,i,xi,xj,xk) = gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*factor_xxx
		        end do
		    end do
        end do
    end do
    
!DEC$ END IF
    
    end subroutine ABX3PARAMETERS

module subroutine allocate_arrays_PCSAFT(gl,nrsubst)
! Henning Markgraf, June 2016

    ! allocates all the module variable arrays for avoiding the double calculations in module_PCSAFT
    ! if they are not already allocated (the result is that the allocation will only be executed on the first call)

    
implicit none

    type(type_gl) :: gl

    integer, intent (inout) :: nrsubst
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    if (.not. gl%already_allocated_PCSAFT) then
        
        ! Allocate all arrays with their dimensions:
        ! T. Eckermann, 02/2017
        allocate(gl%mmeanQij_PCSAFT(gl%ncomp,gl%ncomp))
        allocate(gl%mmeanQijk_PCSAFT(gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%ab_PCSAFTQ(2,0:4,gl%ncomp,gl%ncomp))
        allocate(gl%c_PCSAFTQ(0:4,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%mmeanDij_PCSAFT(gl%ncomp,gl%ncomp))
        allocate(gl%mmeanDijk_PCSAFT(gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%ab_PCSAFTD(2,0:4,gl%ncomp,gl%ncomp))
        allocate(gl%c_PCSAFTD(0:4,gl%ncomp,gl%ncomp,gl%ncomp))
        

        
        !-----------------------------------------------------------------
        ! arrays that store which derivative numbers have been calculated
        
        ! the T, rho derivatives:
        allocate(gl%ar_calculated_PCSAFT(nderivs))
        
        ! the dxi (x1) derivatives:
        allocate(gl%arx1_calculated_PCSAFT(nderivs))
        
        ! the dxi dxj (x2) derivatives
        allocate(gl%arx2_calculated_PCSAFT(nderivs))

        ! the dxi dxj dxk (x3) derivatives
        allocate(gl%arx3_calculated_PCSAFT(nderivs))

        
        !--------------------------------------------------------------
        ! arrays that store the results of the calculated derivatives:
        ! T,rho-derivative variables
        allocate(gl%ar_PCSAFT(nderivs)) ! note: is the option to catch errors allocate(ar_PCSAFT(nderivs),stat=allocation_status) necessary?
        allocate(gl%adisp_PCSAFT(nderivs))
        allocate(gl%c_PCSAFT(nderivs))
        allocate(gl%i1_PCSAFT(nderivs))
        allocate(gl%i2_PCSAFT(nderivs))
        allocate(gl%meo1_PCSAFT(nderivs))
        allocate(gl%meo2_PCSAFT(nderivs))
        allocate(gl%ahc_PCSAFT(nderivs))
        allocate(gl%ahs_PCSAFT(nderivs))
        allocate(gl%z0_PCSAFT(nderivs))
        allocate(gl%z1_PCSAFT(nderivs))
        allocate(gl%z2_PCSAFT(nderivs))
        allocate(gl%z3_PCSAFT(nderivs))
        allocate(gl%gii_PCSAFT(gl%ncomp,nderivs))
        allocate(gl%di_PCSAFT(gl%ncomp,nderivs))
        allocate(gl%gij_PCSAFT(gl%ncomp,nderivs))
        
        ! dxi (x1) derivative variables
        allocate(gl%arx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%adispx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%cx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%i1x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%i2x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%meo1x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%meo2x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%ahcx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%ahsx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%z0x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%z1x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%z2x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%z3x1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%giix1_PCSAFT(gl%ncomp,nx1derivs,gl%ncomp))
        allocate(gl%abx1_PCSAFT(2,0:6,gl%ncomp))
        
        ! dxi dxj (x2) derivative variables
        allocate(gl%arx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%adispx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%cx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%i1x2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%i2x2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%meo1x2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%meo2x2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%ahcx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%ahsx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%giix2_PCSAFT(gl%ncomp,nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%abx2_PCSAFT(2,0:6,gl%ncomp,gl%ncomp))
        
        ! dxi dxj dxk (x3) derivative variables
        allocate(gl%arx3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%adispx3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%cx3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%i1x3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%i2x3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%ahcx3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%ahsx3_PCSAFT(nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%giix3_PCSAFT(gl%ncomp,nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%abx3_PCSAFT(2,0:6,gl%ncomp,gl%ncomp,gl%ncomp))
        
        ! T,rho-derivative variables
        allocate(gl%J2_PCSAFTQ(nderivs,gl%ncomp,gl%ncomp))
        allocate(gl%J3_PCSAFTQ(nderivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2_PCSAFTQ(nderivs))
        allocate(gl%A3_PCSAFTQ(nderivs))
        allocate(gl%AQQ_PCSAFTQ(nderivs))
        allocate(gl%J2_PCSAFTD(nderivs,gl%ncomp,gl%ncomp))
        allocate(gl%J3_PCSAFTD(nderivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2_PCSAFTD(nderivs))
        allocate(gl%A3_PCSAFTD(nderivs))
        allocate(gl%ADD_PCSAFTD(nderivs))
        allocate(gl%AASSOC_PCSAFT(nderivs))
        allocate(gl%xA_PCSAFT(4,2)) ! binary mixtures, maximum of 4 sites
        !allocate(delta_AB(4,2,4,2,nderivs))
        allocate(gl%delta_AB(gl%ncomp,gl%ncomp,4,4))
        
        ! dxi (x1) derivative variables
        allocate(gl%AASSOCx1_PCSAFT(nx1derivs,gl%ncomp))
        allocate(gl%AASSOCx2_PCSAFT(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%J2X1_PCSAFTQ(nx1derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%J3X1_PCSAFTQ(nx1derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2X1_PCSAFTQ(nx1derivs,gl%ncomp))
        allocate(gl%A3X1_PCSAFTQ(nx1derivs,gl%ncomp))
        allocate(gl%AQQX1_PCSAFTQ(nx1derivs,gl%ncomp))
        allocate(gl%J2X1_PCSAFTD(nx1derivs,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%J3X1_PCSAFTD(nx1derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2X1_PCSAFTD(nx1derivs,gl%ncomp))
        allocate(gl%A3X1_PCSAFTD(nx1derivs,gl%ncomp))
        allocate(gl%ADDX1_PCSAFTD(nx1derivs,gl%ncomp))

        ! dxi dxj (x2) derivative variables
        allocate(gl%J2X2_PCSAFTQ(nx2derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%J3X2_PCSAFTQ(nx2derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2X2_PCSAFTQ(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%A3X2_PCSAFTQ(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%AQQX2_PCSAFTQ(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%J2X2_PCSAFTD(nx2derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%J3X2_PCSAFTD(nx2derivs,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp,gl%ncomp))
        allocate(gl%A2X2_PCSAFTD(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%A3X2_PCSAFTD(nx2derivs,gl%ncomp,gl%ncomp))
        allocate(gl%ADDX2_PCSAFTD(nx2derivs,gl%ncomp,gl%ncomp))
        
        
        !-----------------------------------------------------
        
        call init_pure_mix(gl,nrsubst)
        call init_derivs_PCSAFT(gl,nrsubst)        
        gl%already_allocated_PCSAFT = .TRUE.
    end if

!DEC$ END IF
    
end subroutine allocate_arrays_PCSAFT


module subroutine init_pure_mix(gl,nrsubst)

    
implicit none

    type(type_gl) :: gl

    
    integer :: nrsubst
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    if (nrsubst /= 0) then
        if (.not. gl%mol_save == .true.) then !not yet overwritten
            gl%mol_save = .true.
            gl%n_start = nrsubst
            gl%n_end = nrsubst
            gl%molfractions_save = gl%molfractions(nrsubst)
            gl%molfractions(nrsubst) = 1.d0
        end if
    else 
        gl%n_start = 1
        gl%n_end = gl%ncomp
    end if
    
!DEC$ END IF
    
end subroutine init_pure_mix


module subroutine init_mixture_PCSAFT(gl)
! T. Eckermann, 18/03

    
implicit none

    type(type_gl) :: gl

    
    integer :: i

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! Important to calculate mmean and ab here already for the first program run because in ARDERIVS they are only calculated again if x_PCSAFT != molfractions!
    ! mmean
    gl%mmean_PCSAFT = 0.d0
    do i = gl%n_start, gl%n_end
        gl%mmean_PCSAFT = gl%mmean_PCSAFT + gl%molfractions(i)*gl%mPCSAFT(i)
    end do
    
    call ABPARAMETERS(gl)
    
!DEC$ END IF
    
end subroutine init_mixture_PCSAFT
  

module subroutine init_derivs_PCSAFT(gl,nrsubst)
! Henning Markgraf, June 2016

    ! initializes all the module variables for avoiding the double calculations in module_PCSAFT

    
implicit none

    type(type_gl) :: gl

    
    ! Declarations
    integer, intent (inout) :: nrsubst
    integer :: i, j, k
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! Initialize all:
    
    !----------------------------------------------------
    ! vectors that save which derivatives are calculated
    
    ! the T, rho derivatives:
    gl%ar_calculated_PCSAFT = 0
        
    ! the dxi (x1) derivatives:
    gl%arx1_calculated_PCSAFT = 0
        
    ! the dxi dxj (x2) derivatives
    gl%arx2_calculated_PCSAFT = 0

    ! the dxi dxj dxk (x3) derivatives
    gl%arx3_calculated_PCSAFT = 0
     
    !----------------------------------------
    
    ! depending variables:
    gl%T_PCSAFT = 0.d0
    gl%D_PCSAFT = 0.d0
    gl%x_PCSAFT = 0.d0
    gl%nrsubst_PCSAFT = 0
    
    !----------------------------------------
    
    ! Important to calculate mmean and ab here already for the first program run because in ARDERIVS they are only calculated again if x_PCSAFT != molfractions!
    ! mmean
    gl%mmean_PCSAFT = 0.d0
    do i = gl%n_start, gl%n_end
        gl%mmean_PCSAFT = gl%mmean_PCSAFT + gl%molfractions(i)*gl%mPCSAFT(i)
    end do
    
    

    
    
    ! ab parameters
    gl%ab_PCSAFT = 0.d0
    gl%ab_PCSAFTQ = 0.d0
    gl%c_PCSAFTQ = 0.d0
    gl%ab_PCSAFTD = 0.d0
    gl%c_PCSAFTD = 0.d0
    gl%abx1_PCSAFT = 0.d0
    gl%abx2_PCSAFT = 0.d0
    gl%abx3_PCSAFT = 0.d0
    call ABPARAMETERS(gl)
    call ABX1PARAMETERS(gl)
    call ABX2PARAMETERS(gl)
    call ABX3PARAMETERS(gl)
    
    ! T,rho-derivative variables
    gl%AQQ_PCSAFTQ = 0.d0
    gl%A2_PCSAFTQ = 0.d0
    gl%A3_PCSAFTQ = 0.d0
    gl%J2_PCSAFTQ = 0.d0
    gl%J3_PCSAFTQ = 0.d0
   

    gl%ADD_PCSAFTD = 0.d0
    gl%A2_PCSAFTD = 0.d0
    gl%A3_PCSAFTD = 0.d0
    gl%J2_PCSAFTD = 0.d0
    gl%J3_PCSAFTD = 0.d0
    
    gl%AASSOC_PCSAFT = 0.d0
    gl%AASSOCx1_PCSAFT = 0.d0
    gl%AASSOCx2_PCSAFT = 0.d0
    
    ! dxi (gl,x1) derivative variables
    gl%J2X1_PCSAFTQ = 0.d0
    gl%J3X1_PCSAFTQ = 0.d0
    gl%A2X1_PCSAFTQ = 0.d0
    gl%A3X1_PCSAFTQ = 0.d0
    gl%AQQX1_PCSAFTQ = 0.d0
    
    gl%J2X1_PCSAFTD = 0.d0
    gl%J3X1_PCSAFTD = 0.d0
    gl%A2X1_PCSAFTD = 0.d0
    gl%A3X1_PCSAFTD = 0.d0
    gl%ADDX1_PCSAFTD = 0.d0
    
    ! dxj (x2) derivative variables
    gl%J2X2_PCSAFTQ = 0.d0
    gl%J3X2_PCSAFTQ = 0.d0
    gl%A2X2_PCSAFTQ = 0.d0
    gl%A3X2_PCSAFTQ = 0.d0
    gl%AQQX2_PCSAFTQ = 0.d0
    
    gl%J2X2_PCSAFTD = 0.d0
    gl%J3X2_PCSAFTD = 0.d0
    gl%A2X2_PCSAFTD = 0.d0
    gl%A3X2_PCSAFTD = 0.d0
    gl%ADDX2_PCSAFTD = 0.d0
    !---------------------------------------
    ! arrays that store the results of the calculated derivatives:
    
    ! T,rho-derivative variables
    gl%ar_PCSAFT = 0.d0
    gl%adisp_PCSAFT = 0.d0
    gl%c_PCSAFT = 0.d0
    gl%i1_PCSAFT = 0.d0
    gl%i2_PCSAFT = 0.d0
    gl%meo1_PCSAFT = 0.d0
    gl%meo2_PCSAFT = 0.d0
    gl%ahc_PCSAFT = 0.d0
    gl%ahs_PCSAFT = 0.d0
    gl%z0_PCSAFT = 0.d0
    gl%z1_PCSAFT = 0.d0
    gl%z2_PCSAFT = 0.d0
    gl%z3_PCSAFT = 0.d0
    gl%gii_PCSAFT = 0.d0
    gl%di_PCSAFT = 0.d0
    gl%gij_PCSAFT = 0.d0
    
    ! dxi (x1) derivative variables
    gl%arx1_PCSAFT = 0.d0
    gl%adispx1_PCSAFT = 0.d0
    gl%cx1_PCSAFT = 0.d0
    gl%i1x1_PCSAFT = 0.d0
    gl%i2x1_PCSAFT = 0.d0
    gl%meo1x1_PCSAFT = 0.d0
    gl%meo2x1_PCSAFT = 0.d0
    gl%ahcx1_PCSAFT = 0.d0
    gl%ahsx1_PCSAFT = 0.d0
    gl%z0x1_PCSAFT = 0.d0
    gl%z1x1_PCSAFT = 0.d0
    gl%z2x1_PCSAFT = 0.d0
    gl%z3x1_PCSAFT = 0.d0
    gl%giix1_PCSAFT = 0.d0
    
    ! dxi dxj (x2) derivative variables
    gl%arx2_PCSAFT = 0.d0
    gl%adispx2_PCSAFT = 0.d0
    gl%cx2_PCSAFT = 0.d0
    gl%i1x2_PCSAFT = 0.d0
    gl%i2x2_PCSAFT = 0.d0
    gl%meo1x2_PCSAFT = 0.d0
    gl%meo2x2_PCSAFT = 0.d0
    gl%ahcx2_PCSAFT = 0.d0
    gl%ahsx2_PCSAFT = 0.d0
    gl%giix2_PCSAFT = 0.d0
    
    ! dxi dxj dxk (x3) derivative variables
    gl%arx3_PCSAFT = 0.d0
    gl%adispx3_PCSAFT = 0.d0
    gl%cx3_PCSAFT = 0.d0
    gl%i1x3_PCSAFT = 0.d0
    gl%i2x3_PCSAFT = 0.d0
    gl%ahcx3_PCSAFT = 0.d0
    gl%ahsx3_PCSAFT = 0.d0
    gl%giix3_PCSAFT = 0.d0
    
    gl%delta_AB = 0.d0
   
!DEC$ END IF
    
end subroutine init_derivs_PCSAFT
    
module subroutine check_calculated_PCSAFT(gl,T,D,nrsubst,calculated_derivs,getvector,already_calculated)

! Henning Markgraf, June 2016

    ! checks if PC SAFT has already been calculated with the same T, D, molfractions and getvector
    
    ! 3 possibilities:
    ! 1. All needed PC SAFT derivatives have already been saved in the module variables:
    !           --> Nothing shall be calculated again
    !           --> getvector set to 0
    ! 2. T, D and molfractions have already been calculated, but not all the derivatives:
    !           --> Only the missing derivatives shall be calculated
    !           --> all positions in the getvector which are already saved in the module variables are set to 0
    ! 3. this T, D and molfractions combination is not saved in the module variables:
    !           --> Everything has to be calculated
    !           --> all module variables are initialized

    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: calculated_derivs
    integer, intent (in) :: nrsubst
    ! input and output at the same time
    integer, dimension (nderivs), intent (inout) :: getvector
    ! output
    logical, intent (out) :: already_calculated
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! Initializations
    already_calculated = .FALSE.
    
    ! III. check if the values stored in the module variables are the needed values already
    !Andreas, Erik, October 2019, issues regarding saved values; that is why everything is commented out
    !if (gl%T_PCSAFT .EQ. T .AND. gl%D_PCSAFT .EQ. D .AND. gl%nrsubst_PCSAFT .EQ. nrsubst .AND. all(gl%molfractions .EQ. gl%x_PCSAFT) .AND. nrsubst .EQ. gl%nrsubst_PCSAFT) then ! if T, D and molfracs the same
    !    already_calculated = .TRUE.
	   ! call remove_redundant_getvector_entries(gl,calculated_derivs,getvector) ! only calculate the missing ones
    !    ! The getvector can be 0 after this subroutine, if all required derivs have already been calculated
    !else ! this T, D and molfracs combination has NOT been saved in the module variables already
    !    already_calculated = .FALSE.
    !end if
    
!DEC$ END IF
    
    end subroutine check_calculated_PCSAFT

module subroutine remove_redundant_getvector_entries(gl,calculated_derivs,getvector)

! Henning Markgraf, June 2016

    ! compares the two given arrays and sets all positions of getvector to 0, where both vectors have a 1
    
    ! Explanation:
    ! Input:           calculated_derivs        must be one of the module_PCSAFT arrays Trho_calculated_PCSAFT,
    !                                           x1_calculated_PCSAFT, x2_calculated_PCSAFT or x3_calculated_PCSAFT
    !                                           which save which derivatives have already been calculated
    !
    ! Input/ Output:    getvector               the array which specifies which derivatives are required
    !
    ! The derivatives with a 1 in the getvector shall not be calculated again, if they are
    ! already saved in calculated_derivs. Positions in the getvector are therefore set to 0
    ! if they are a 1 in both arrays
    
    
implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: calculated_derivs
    ! input and output at the same time
    integer, dimension (nderivs), intent (inout) :: getvector
    !working variables
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. remove redundant entries
    do i = 1, nderivs
        if (getvector(i) .EQ. 1 .AND. calculated_derivs(i) .EQ. 1) then
            getvector(i) = 0
        end if
    end do
    
!DEC$ END IF
    
    end subroutine remove_redundant_getvector_entries
    
module subroutine update_modulevariables_PCSAFT(gl,T,D,nrsubst,calculated_derivs,get)

! Henning Markgraf, June 2016

    ! save the depending variables at which the PC SAFT routines are executed


implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input
    double precision, intent (in) :: T, D
    integer, dimension (nderivs), intent (in) :: get
    integer, intent (in) :: nrsubst
    ! input and output into the module variable passed to the subroutine:
    integer, dimension (nderivs), intent (inout) :: calculated_derivs
    ! working variable
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%T_PCSAFT = T
    gl%D_PCSAFT = D
    gl%x_PCSAFT(1:gl%ncomp) = gl%molfractions(1:gl%ncomp)
    gl%nrsubst_PCSAFT = nrsubst
    do i = 1, nderivs
        if (get(i) .EQ. 1) then
            calculated_derivs(i) = 1
        end if
    end do
    
!DEC$ END IF
    
end subroutine update_modulevariables_PCSAFT



    end submodule impl
