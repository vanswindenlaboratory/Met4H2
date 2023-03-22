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

    ! module for file pc-saft_AHSX_derivs.f90
    module pc_saft_AHSX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains

    
    
subroutine AHSX1DERIVS(gl,GETDERAHS)

! Henning Markgraf, June 2016

    ! a_hs: Helmholtz free energy of the hard sphere fluid
    ! defined by eq. A.6 in Gross, Sadowski 2001:
    ! a_hs = 1/zeta_0 * ( 3*zeta_1*zeta_2/(1-zeta_3) + zeta_2**3/(zeta_3*(1-zeta_3)**2)
    !           + (zeta_2**3/zeta_3**2 - zeta_0)*ln(1-zeta_3) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERAHS
    !output: ahsx1_PCSAFT (module variable)
    !working variables
    double precision :: a_hs
    integer :: i, xi
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    
    ! II. X1 DERIVATIVES of a_hs
    ! The basic equation a_hs is a recurring factor in the derivatives
    a_hs = (3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)** &
	      & 2) + (-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)
    
    !calculate the derivatives of a_hs
    ! 1: a_hs
    if (GETDERAHS(1) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(1,xi) = (-a_hs*gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2* &
			& gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z2_PCSAFT(1) &
			& **2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0) &
			& **2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)) - gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*( &
			& gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) - (gl%z3_PCSAFT(1) - 1.d0)**3*( &
			& gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0))/(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3*( &
			& gl%z3_PCSAFT(1) - 1.d0)**3)
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, d_i, d_j
    if (GETDERAHS(2) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(2,xi) = (a_hs*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z0_PCSAFT(1)*(gl%z3_PCSAFT(1) &
			& - 1.d0)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*( &
			& gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + 3.d0* &
			& gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1) &
			& ) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1) &
			& **3) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
			& )) + gl%z0_PCSAFT(1)*(6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 4.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) &
			& *(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& *gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)) - 6.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
			& gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0) &
			& **2*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(6.d0*gl%z2_PCSAFT(1)** &
			& 2*gl%z2x1_PCSAFT(1,xi) + gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - 2.d0* &
			& gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) &
			& + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(-gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) - (gl%z3_PCSAFT(1) - 1.d0)**4*( &
			& gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*( &
			& gl%z3_PCSAFT(1) - 1.d0)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 2.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
			& gl%z2_PCSAFT(1)**3)*dlog(-gl%z3_PCSAFT(1) + 1.d0)))/(gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - &
			& 1.d0)**4)
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(3) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(3,xi) = (gl%z0x1_PCSAFT(1,xi)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)** &
			& 2 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 - 5.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(1)**2 + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) + 3.d0*gl%z2_PCSAFT(1)**3) + gl%z1x1_PCSAFT(1,xi)*(-6.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)) + gl%z2x1_PCSAFT(1,xi)*(-6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 3.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2) + gl%z3x1_PCSAFT(1,xi)*(-2.d0*gl%z0_PCSAFT(1)**2 &
			& *gl%z3_PCSAFT(1)**3 + 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - &
			& 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) + 2.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 10.d0*gl%z0_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) - 16.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3))/(gl%z0_PCSAFT(1)**2*(gl%z3_PCSAFT(1)**5 - &
			& 5.d0*gl%z3_PCSAFT(1)**4 + 10.d0*gl%z3_PCSAFT(1)**3 - 10.d0*gl%z3_PCSAFT(1)**2 + 5.d0*gl%z3_PCSAFT(1) - 1.d0))
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(4) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(4,xi) = (gl%z0_PCSAFT(1)*(-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - &
			& 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2 - 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*( &
			& gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1) &
			& *gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)) - gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
			& -6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)* &
			& gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0* &
			& gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)))*dlog( &
			& -gl%z3_PCSAFT(1) + 1.d0) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1) + &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& *gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)) - gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& **2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(4,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
			& gl%z2_PCSAFT(1)**3)) - gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*(2.d0* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
			& + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + gl%z2_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3* &
			& (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)))/(gl%z0_PCSAFT(1) &
			& **2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4)
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(5) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(5,xi) = (18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) &
			& - 1.d0) - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - &
			& 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) &
			& **2 - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - &
			& 1.d0)**3 + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) - 2.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& **2*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 2.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
			& 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))) - 2.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + &
			& 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0 &
			& *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0* &
			& gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0* &
			& gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi) + 2.d0* &
			& gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(5,xi)*gl%z2_PCSAFT(1) &
			& + 2.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(5)) + 3.d0*gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 2.d0*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(4,xi) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1x1_PCSAFT(1,xi)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) - 6.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z1_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
			& + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4)**2 + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2) + 3.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1) &
			& **2*gl%z2x1_PCSAFT(5,xi) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)) - gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& **3*gl%z3x1_PCSAFT(5,xi)*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) - &
			& gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi) + 6.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
			& + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
			& gl%z2_PCSAFT(1)**3) - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3 &
			& )) - 2.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(4,xi) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 &
			& + gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + gl%z0_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) - gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi))) - gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(6.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)** &
			& 2)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) - gl%z0_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0 &
			& )**5*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0 &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 2.d0* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4) &
			& **2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*( &
			& gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1) &
			& *gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)))*dlog(-gl%z3_PCSAFT(1) + &
			& 1.d0) + gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)** &
			& 4*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
			& 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*( &
			& gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 2.d0*gl%z2_PCSAFT(1)** &
			& 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4) &
			& *gl%z3_PCSAFT(1)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + &
			& 2.d0*gl%z2_PCSAFT(4)**2) - gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 &
			& - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0* &
			& gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2))*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
			& ) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(4) + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*( &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%z1_PCSAFT(4) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)**2*(-gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 + gl%z2_PCSAFT(1)**3))))/(gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**5)
        end do
    end if
    
    
    ! 6: 2ND MIXED DERIVATIVE OF a_hs WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(6) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(6,xi) = (gl%z0x1_PCSAFT(1,xi)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**3 + 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) - 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 9.d0* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0*gl%z1_PCSAFT(4) &
			& *gl%z2_PCSAFT(1) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 10.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 8.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 &
			& + 15.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 21.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
			& + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)) + gl%z1x1_PCSAFT(4,xi)*(3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 &
			& - 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)) + gl%z1x1_PCSAFT(1,xi)*(-6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 9.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)) + gl%z2x1_PCSAFT(4,xi)*(3.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3 - &
			& 15.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 21.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 9.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2) + gl%z2x1_PCSAFT(1,xi)*(-6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1) &
			& **2 + 9.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 30.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 30.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 42.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)) + &
			& gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 3.d0* &
			& gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1) - gl%z0_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + &
			& 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) - 2.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 + 10.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) - 8.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3) + gl%z3x1_PCSAFT(1,xi)*(-2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
			& 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
			& - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 30.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3 &
			& *gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 30.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)))/( &
			& gl%z0_PCSAFT(1)**2*(gl%z3_PCSAFT(1)**5 - 5.d0*gl%z3_PCSAFT(1)**4 + 10.d0*gl%z3_PCSAFT(1)**3 - 10.d0*gl%z3_PCSAFT(1)**2 &
			& + 5.d0*gl%z3_PCSAFT(1) - 1.d0))
        end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hs WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(7) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(7,xi) = (gl%z0x1_PCSAFT(1,xi)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 18.d0*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(5) + 36.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 6.d0*gl%z1_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 36.d0* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 36.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) - 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) - 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**4 + 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**2 + 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) - 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) + 12.d0* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 36.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
			& *gl%z3_PCSAFT(4) + 36.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 + 24.d0*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 36.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 24.d0*gl%z1_PCSAFT(4) &
			& *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4) - 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**4 + 12.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 + 12.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1) + 2.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 &
			& - 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 36.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
			& **2 + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 30.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
			& 8.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
			& - 72.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 108.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 48.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) - 3.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 &
			& - 36.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 + 30.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**4 + &
			& 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**2 &
			& + 60.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) - 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2) + &
			& gl%z1x1_PCSAFT(5,xi)*(3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& **3 + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)) + gl%z1x1_PCSAFT(4,xi)*(-12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)) + &
			& gl%z1x1_PCSAFT(1,xi)*(-6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(5) - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z0_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(5)) + &
			& gl%z2x1_PCSAFT(5,xi)*(3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& **3 + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)** &
			& 2*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 30.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2 &
			& *gl%z3_PCSAFT(1) + 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2) + gl%z2x1_PCSAFT(4,xi)*(-12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(4) - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 72.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 108.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 48.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 - 72.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 144.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 120.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1) + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)) + gl%z2x1_PCSAFT(1,xi)*(-6.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
			& **2 + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1) &
			& *gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4) - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(5)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(5)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) &
			& **3*gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 108.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 54.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 90.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 24.d0 &
			& *gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 144.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 216.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 96.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 + 72.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) &
			& **4 - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)**2* &
			& gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(4)**2 &
			& ) + gl%z3x1_PCSAFT(5,xi)*(gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4 - 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3 + 6.d0 &
			& *gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1) + gl%z0_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) + 8.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3) + gl%z3x1_PCSAFT(4,xi)*(-4.d0*gl%z0_PCSAFT(1)** &
			& 2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 12.d0* &
			& gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 72.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
			& *gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) - 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4) &
			& *gl%z2_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 72.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 60.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 108.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 48.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)) + gl%z3x1_PCSAFT(1,xi)*(-2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z0_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 6.d0*gl%z0_PCSAFT(1) &
			& **2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 2.d0*gl%z0_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 72.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 72.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
			& **2*gl%z3_PCSAFT(4) - 72.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) + 36.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 72.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
			& **2 - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + &
			& 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(5) - 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 36.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 144.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + &
			& 30.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 216.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4) + 180.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
			& gl%z3_PCSAFT(1)**2 - 54.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 24.d0*gl%z0_PCSAFT(1)* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**3 + 72.d0 &
			& *gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**2 - 108.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) &
			& **2*gl%z3_PCSAFT(1) + 48.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2))/(gl%z0_PCSAFT(1)**2*(gl%z3_PCSAFT(1)**6 &
			& - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0*gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1)**2 - 6.d0* &
			& gl%z3_PCSAFT(1) + 1.d0))
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(8) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(8,xi) = 2.d0*(gl%z0x1_PCSAFT(1,xi)*(-9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3 - 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4 + 6.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 8.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
			& ) + gl%z1x1_PCSAFT(1,xi)*(9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2) + gl%z2x1_PCSAFT(1,xi)*(9.d0*gl%z0_PCSAFT(1)* &
			& gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**2 + 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 24.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xi)*(3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4 - 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1) &
			& **3 + 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + &
			& 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 + &
			& 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 + 42.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) + 8.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3))/(gl%z0_PCSAFT(1)**2*(gl%z3_PCSAFT(1)**6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0* &
			& gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z3_PCSAFT(1) + 1.d0))
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF g_ij WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT
    !           z0_PCSAFT_TTT, z1_PCSAFT_TTT, z2_PCSAFT_TTT, z3_PCSAFT_TTT
    if (GETDERAHS(9) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(9,xi) = (-72.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) &
			& - 1.d0) + 120.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi) + 96.d0 &
			& *gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + &
			& 72.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) &
			& **2 + 48.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*( &
			& gl%z3_PCSAFT(1) - 1.d0)**3 - 72.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - &
			& 1.d0)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
			& *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) - 6.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**3*(6.d0*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))) - 3.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi)*(gl%z3_PCSAFT(1) - 1.d0)**5*(2.d0*gl%z2_PCSAFT(1) &
			& *gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0) &
			& **2*(3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 + 9.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
			& + 18.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4)**2) + gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z2_PCSAFT(1) &
			& **2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 2.d0 &
			& *gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(1) &
			& *gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)** &
			& 2 + 18.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(4,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 3.d0*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 36.d0 &
			& *gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4)**2) + 18.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2) - 3.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0) &
			& **5*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(9,xi) + 3.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(5,xi) + 3.d0*gl%z1_PCSAFT(5) &
			& *gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(9)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(9,xi)*gl%z2_PCSAFT(1) + 3.d0* &
			& gl%z1x1_PCSAFT(5,xi)*gl%z2_PCSAFT(4) + 3.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(5) + gl%z1x1_PCSAFT(1,xi)* &
			& gl%z2_PCSAFT(9)) + 3.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3x1_PCSAFT(9,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(5) + gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 6.d0* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + &
			& 6.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) &
			& + 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(9)*gl%z2_PCSAFT(1)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1x1_PCSAFT(5,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%z1x1_PCSAFT(4,xi)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 3.d0*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 3.d0* &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)) - 6.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0) &
			& **3*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
			& *gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
			& gl%z3_PCSAFT(5) + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 3.d0*gl%z1_PCSAFT(5)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 3.d0 &
			& *gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 3.d0*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(4)**2) + 3.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**2* &
			& gl%z2x1_PCSAFT(9,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z2x1_PCSAFT(1,xi) + 6.d0* &
			& gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)) - 2.d0* &
			& gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(9,xi) + 9.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(5) + 3.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2* &
			& gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z2_PCSAFT(1) &
			& *gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%z2_PCSAFT(4) &
			& **3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(9.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 27.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 &
			& - gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - gl%z0_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)**4*gl%z3x1_PCSAFT(9,xi)*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
			& gl%z2_PCSAFT(1)**3) - gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(9,xi) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 9.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 9.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)* &
			& gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 18.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
			& + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + &
			& 6.d0*gl%z2_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 3.d0 &
			& *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) - 3.d0* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) - gl%z3_PCSAFT(9) &
			& *gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - 2.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**4* &
			& (gl%z3_PCSAFT(1) - 1.d0)**3*(-6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) - 6.d0*gl%z2_PCSAFT(1)**3 &
			& *gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) - 36.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) &
			& - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - &
			& 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 3.d0*gl%z3_PCSAFT(4)**2* &
			& gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + 3.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
			& *gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - 2.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) - gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2*( &
			& gl%z3_PCSAFT(1) - 1.d0)**5*(3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi))) - 3.d0*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*(6.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)** &
			& 2)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) - 3.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**4*(-8.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0* &
			& gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi))) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(6.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) - &
			& gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) - gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
			& gl%z3_PCSAFT(1) - 1.d0)**5*(3.d0*gl%z3_PCSAFT(4)*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - &
			& 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2* &
			& gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
			& *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) &
			& *gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi) + 4.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4) &
			& **2*gl%z2x1_PCSAFT(1,xi))) + gl%z3x1_PCSAFT(1,xi)*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 18.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
			& ) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) &
			& ) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**3))) - gl%z0_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**6*(-120.d0* &
			& gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi) + 72.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*( &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) - 6.d0*gl%z2_PCSAFT(1) &
			& *gl%z3_PCSAFT(1)**2*(3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
			& gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 9.d0 &
			& *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 + 9.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
			& gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(9,xi) + &
			& 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(4,xi) + &
			& 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(4,xi) + 6.d0* &
			& gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)) + 2.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(9,xi) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 9.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 9.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)* &
			& gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 18.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
			& + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + &
			& 6.d0*gl%z2_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)))* &
			& dlog(-gl%z3_PCSAFT(1) + 1.d0) + gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(-18.d0*gl%z1_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 1.d0) + 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
			& **4*gl%z3_PCSAFT(4)**3 + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 1.d0 &
			& ) + 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1) &
			& **3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 1.d0)**3 - 18.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4* &
			& gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + &
			& 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**3*(-2.d0*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(5) - 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)*(-2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + &
			& 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 1.d0 &
			& )**4*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 &
			& *(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) &
			& ) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + 18.d0* &
			& gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)** &
			& 5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5) + 3.d0* &
			& gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4) + gl%z1_PCSAFT(9)*gl%z2_PCSAFT(1)) - 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0 &
			& )**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 3.d0* &
			& gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0* &
			& gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) - 3.d0* &
			& gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**3) + 2.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*( &
			& gl%z3_PCSAFT(1) - 1.d0)**2*(-6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)** &
			& 3*gl%z3_PCSAFT(9)*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + &
			& gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 9.d0*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) - 3.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
			& gl%z2_PCSAFT(1)**3)) + (gl%z3_PCSAFT(1) - 1.d0)**5*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 18.d0* &
			& gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
			& ) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(5) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) &
			& ) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**3))*dlog(-gl%z3_PCSAFT(1) + 1.d0)))/(gl%z0_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**6)
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(10) .eq. 1) then
        do xi = 1, gl%ncomp
            gl%ahsx1_PCSAFT(10,xi) = (-2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*gl%z3x1_PCSAFT(4,xi) + 4.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) - 6.d0* &
			& gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2* &
			& gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z0_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(4,xi) - 36.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 36.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& *gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - &
			& 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**4 + 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 18.d0* &
			& gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 6.d0* &
			& gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(4,xi)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 6.d0* &
			& gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
			& gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)* &
			& gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - &
			& 12.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
			& gl%z3x1_PCSAFT(4,xi) + 90.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 16.d0*gl%z0_PCSAFT(1) &
			& *gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3* &
			& gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - &
			& 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 48.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& **2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)** &
			& 4 + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1) &
			& **2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) &
			& + 9.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2* &
			& gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 48.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
			& gl%z3_PCSAFT(1)**3 - 36.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 12.d0* &
			& gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z2_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) - 12.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3* &
			& gl%z3_PCSAFT(4) + 18.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 6.d0* &
			& gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 18.d0* &
			& gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 18.d0* &
			& gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)* &
			& gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 2.d0* &
			& gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 12.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3 &
			& *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 16.d0 &
			& *gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 3.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
			& gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 18.d0* &
			& gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
			& gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 9.d0*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4))/(gl%z0_PCSAFT(1)**2*( &
			& gl%z3_PCSAFT(1)**6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0*gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1) &
			& **2 - 6.d0*gl%z3_PCSAFT(1) + 1.d0))
        end do
    end if
    
    !DEC$ END IF
end subroutine AHSX1DERIVS
    
   
    
subroutine AHSX2DERIVS(gl,GETDERAHS)

! Henning Markgraf, June 2016

    ! a_hs: Helmholtz free energy of the hard sphere fluid
    ! defined by eq. A.6 in Gross, Sadowski 2001:
    ! a_hs = 1/zeta_0 * ( 3*zeta_1*zeta_2/(1-zeta_3) + zeta_2**3/(zeta_3*(1-zeta_3)**2)
    !           + (zeta_2**3/zeta_3**2 - zeta_0)*ln(1-zeta_3) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERAHS
    !output: ahsx2_PCSAFT (module variable)
    !working variables
    integer :: xi, xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! recurring factors in the derivatives
    ! note: to do

    !calculate the derivatives of a_hs
    ! 1: a_hs
    if (GETDERAHS(1) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
              !  if (xi .GE. xj) then
		        gl%ahsx2_PCSAFT(1,xi,xj) = (-6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*( &
             & gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)* &
             & gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
             & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) &
             & *gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*( &
             & gl%z3_PCSAFT(1) - 1.d0)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + 6.d0* &
             & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2 &
             & + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
             & gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)* &
             & gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + &
             & 1.d0) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1x1_PCSAFT(1,xj)* &
             & gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*( &
             & gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
             & gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
             & gl%z3x1_PCSAFT(1,xj)) + gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-3.d0*gl%z2_PCSAFT(1)**2* &
             & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
             & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - gl%z0_PCSAFT(1)**2 &
             & *gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0* &
             & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + &
             & gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0* &
             & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*( &
             & gl%z0x1_PCSAFT(1,xj)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + &
             & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*( &
             & gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + 3.d0* &
             & gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1) &
             & ) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1) &
             & **3) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
             & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
             & )) + gl%z0x1_PCSAFT(1,xi)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - &
             & 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
             & gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*( &
             & gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) + &
             & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)* &
             & gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + &
             & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))*dlog( &
             & -gl%z3_PCSAFT(1) + 1.d0))) - 2.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - &
             & 1.d0)*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
             & (-gl%z3_PCSAFT(1) + 1.d0) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)* &
             & dlog(-gl%z3_PCSAFT(1) + 1.d0)))/(gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4)
             !   end if
            end do
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, d_i, d_j
    if (GETDERAHS(2) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		gl%ahsx2_PCSAFT(2,xi,xj) = (-2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 2.d0*gl%z0_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0 &
     & *gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)** &
     & 2 + 9.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0 &
     & *gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)** &
     & 2 + 9.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 30.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 30.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2 &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj) + 30.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 30.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 42.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) &
     & + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1) &
     & + 2.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 10.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 15.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - &
     & 21.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 9.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) &
     & - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) &
     & - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1) &
     & + 2.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 10.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 8.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 15.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - &
     & 21.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 9.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 + 18.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 6.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) + 2.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 - 10.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)** &
     & 2 + 14.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) - 6.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3)/(gl%z0_PCSAFT(1)**3*(gl%z3_PCSAFT(1)**5 - 5.d0*gl%z3_PCSAFT(1)**4 + 10.d0* &
     & gl%z3_PCSAFT(1)**3 - 10.d0*gl%z3_PCSAFT(1)**2 + 5.d0*gl%z3_PCSAFT(1) - 1.d0))
                end if
	        end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(3) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		gl%ahsx2_PCSAFT(3,xi,xj) = (4.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z0_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - &
     & 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 36.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) &
     & - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **2 + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 + 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 90.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)** &
     & 2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) &
     & + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 6.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 16.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
     & *gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi) - 12.d0 &
     & *gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1) &
     & *gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - 6.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 16.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 3.d0 &
     & *gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 - 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj) - 12.d0 &
     & *gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 36.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1) - 2.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4 + 12.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 4.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(1) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3)/(gl%z0_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & **6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0*gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1)**2 - &
     & 6.d0*gl%z3_PCSAFT(1) + 1.d0))
                end if
	        end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(4) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		gl%ahsx2_PCSAFT(4,xi,xj) = (18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)** &
     & 2 - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3 + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0 &
     & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 2.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*( &
     & gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) &
     & ) + gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)** &
     & 2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2) - &
     & gl%z3x1_PCSAFT(1,xj)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi))) - gl%z3x1_PCSAFT(1,xi)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)))) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(4,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - &
     & 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) + gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0 &
     & )**3*(-3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
     & - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z2_PCSAFT(1)** &
     & 3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0 &
     & *gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) - gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3x1_PCSAFT(4,xj)*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(4,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + gl%z0_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + 6.d0* &
     & gl%z0_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**5*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)))*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0) + gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z0x1_PCSAFT(1,xj)*(6.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)) + gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0 &
     & *gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0* &
     & gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)))*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xi)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)) + gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1))))) - 2.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*( &
     & -3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + gl%z2_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0 &
     & )**3*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
     & ) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%z1_PCSAFT(4)* &
     & gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)))/(gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**5)
                end if
	        end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(5) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		gl%ahsx2_PCSAFT(5,xi,xj) = (-72.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 120.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 96.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & **4*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 72.d0*gl%z0_PCSAFT(1)**2 &
     & *gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - &
     & 1.d0)**2 + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**3 - 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*( &
     & gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xj)) - 4.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(6.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 18.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3_PCSAFT(1))) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) &
     & **2*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(4)**2) + 2.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
     & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3_PCSAFT(1)) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(-gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) - gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2) + 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2) - 3.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z1x1_PCSAFT(5,xj)*gl%z2x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(5,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(4,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(5,xi) + gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(5,xj)) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) &
     & + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(5,xj) &
     & + 2.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(5,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z1x1_PCSAFT(5,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 2.d0*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(4,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xj) + 2.d0*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + &
     & 2.d0*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 2.d0*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 2.d0 &
     & *gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)**2) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5 &
     & *(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(5,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) + gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) + gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(5,xj) + 4.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) &
     & + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + 4.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 4.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + &
     & 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + gl%z0_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4*(-3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) - 3.d0 &
     & *gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(5,xj) - 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 12.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) - 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1) &
     & *gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) - &
     & 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) - 6.d0*gl%z2_PCSAFT(4)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 12.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) &
     & *(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1) &
     & *gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + 2.d0*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 12.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 &
     & - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3 &
     & ) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1) &
     & **3) - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1) &
     & **3)) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(-6.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - 9.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - 2.d0*gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) - gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) &
     & **5*(2.d0*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi))) + 2.d0*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*(-6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj))) + gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(5,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi)*(gl%z0x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1))) + gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**5*(12.d0*gl%z3_PCSAFT(4)*(-4.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0 &
     & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
     & *(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0 &
     & *gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi))) - gl%z3x1_PCSAFT(1,xj)*(24.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*( &
     & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1) &
     & *gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi) + 4.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi))) - gl%z3x1_PCSAFT(1,xi)*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(5,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xj) + 4.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(4) &
     & **2*gl%z2x1_PCSAFT(1,xj)))) + gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(24.d0*gl%z2_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi))) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(-6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj))) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & (gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*(gl%z0x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) - 6.d0*gl%z0_PCSAFT(1)**2*( &
     & gl%z3_PCSAFT(1) - 1.d0)**6*(-20.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & 4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1) &
     & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) - gl%z3_PCSAFT(1)**4*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi)*gl%z2x1_PCSAFT(1,xj) + &
     & 2.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(5)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) + &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) + gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(5,xj) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) &
     & + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) &
     & + 4.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)*(gl%z0x1_PCSAFT(1,xj)*(-18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z3_PCSAFT(1) - 1.d0) + 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3 - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) - 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0* &
     & gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + 3.d0* &
     & gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi) + 2.d0*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(5,xi)*gl%z2_PCSAFT(1) + 2.d0* &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(5)) - 3.d0*gl%z3_PCSAFT(1)**5*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z1x1_PCSAFT(4,xi) &
     & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 6.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) &
     & **2 + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)**2) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(5,xi) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xi)) + gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi)* &
     & (gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) &
     & - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - 2.d0*gl%z2_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)** &
     & 2*(gl%z3_PCSAFT(1) - 1.d0)**3*(-gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) &
     & + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi))) - gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1) &
     & **2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) + (gl%z3_PCSAFT(1) - 1.d0)**5*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xi) + 4.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4) &
     & **2*gl%z2x1_PCSAFT(1,xi)))*dlog(-gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xi)*(-18.d0*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) + 24.d0*gl%z2_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) + 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**3 - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*( &
     & gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*(gl%z3_PCSAFT(1) - 1.d0)**4* &
     & (2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) - 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4) &
     & *gl%z3_PCSAFT(1))) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(5,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(5,xj) + 2.d0* &
     & gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj) + gl%z1_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(5,xj)*gl%z2_PCSAFT(1) &
     & + 2.d0*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(4) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(5)) - 3.d0*gl%z3_PCSAFT(1)**5*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xj) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z1x1_PCSAFT(4,xj) &
     & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 6.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) &
     & **2 + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)**2) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(5,xj) + 4.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(4)**2*gl%z2x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xj)* &
     & (gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5) + 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) &
     & - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) - 2.d0*gl%z2_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) - 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)** &
     & 2*(gl%z3_PCSAFT(1) - 1.d0)**3*(-gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) &
     & + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj))) - gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1) &
     & **2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*( &
     & gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + (gl%z3_PCSAFT(1) - 1.d0)**5*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(5,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) &
     & *gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(5,xj) + 4.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(4) &
     & **2*gl%z2x1_PCSAFT(1,xj)))*dlog(-gl%z3_PCSAFT(1) + 1.d0))) - 2.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 &
     & *(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 4.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & **2*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) - 3.d0*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2) - &
     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & *(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4) + &
     & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
     & - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 6.d0 &
     & *gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)**2*(-gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + &
     & gl%z2_PCSAFT(1)**3))))/(gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**6)
                end if
	        end do
        end do
    end if
    
    
    ! 6: 2ND MIXED DERIVATIVE OF a_hs WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(6) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		gl%ahsx2_PCSAFT(6,xi,xj) = (-2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 2.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 12.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z0_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) - 72.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 72.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 36.d0 &
     & *gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) &
     & + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xj) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 &
     & *gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
     & 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) &
     & + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **3 + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) &
     & + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) &
     & **3 + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(4,xi)* &
     & gl%z2x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) &
     & + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) - 36.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
     & 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 &
     & *gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + &
     & 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi) &
     & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)** &
     & 4 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)**2 &
     & *gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) &
     & *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) &
     & - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 144.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 30.d0*gl%z0_PCSAFT(1)**2 &
     & *gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 30.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)** &
     & 2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 108.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 90.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) - 54.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) &
     & + 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1) &
     & **2*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 54.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 &
     & *gl%z3x1_PCSAFT(4,xi) - 108.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xi) - 54.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) + 90.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) + 18.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 36.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) - 108.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 54.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 90.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)** &
     & 2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 24.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 72.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 108.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) + 72.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj) - 108.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj) + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
     & 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 &
     & - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 72.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 72.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 108.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 48.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 &
     & + 72.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 60.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 6.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(4,xi)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1) &
     & *gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)* &
     & gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4) + 2.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi) - 30.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) &
     & - 8.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj) &
     & *gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 54.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**3 - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2 + 30.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi) + &
     & 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 36.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 54.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 24.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1) &
     & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1) &
     & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 3.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 6.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 18.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - &
     & 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - &
     & 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)* &
     & gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4) + 2.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) - 30.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - 8.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj) + 6.d0*gl%z0_PCSAFT(1) &
     & *gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 54.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**4 + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**3 - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)**2 + 30.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj) + &
     & 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 36.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 54.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 24.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) - 12.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
     & **4 - 24.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 36.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **4 - 24.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 36.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1) - 4.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) + 24.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 36.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 16.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 6.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 18.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4))/(gl%z0_PCSAFT(1)**3*(gl%z3_PCSAFT(1)**6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0* &
     & gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z3_PCSAFT(1) + 1.d0))
                end if
	        end do
        end do
    end if

    !DEC$ END IF
end subroutine AHSX2DERIVS
    

    
subroutine AHSX3DERIVS(gl,GETDERAHS)

! Henning, March 2016

    ! a_hs: Helmholtz free energy of the hard sphere fluid
    ! defined by eq. A.6 in Gross, Sadowski 2001:
    ! a_hs = 1/zeta_0 * ( 3*zeta_1*zeta_2/(1-zeta_3) + zeta_2**3/(zeta_3*(1-zeta_3)**2)
    !           + (zeta_2**3/zeta_3**2 - zeta_0)*ln(1-zeta_3) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERAHS
    !output: ahsx3_PCSAFT (module variable)
    !working variables
    integer :: i, xi, xj, xk

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! note: possible to enhance calculation time, when calculating the parts which are not derivatives wrt to composition before the loops
    ! recurring factors in the derivatives
    ! note: to do
    
    ! III. X3 DERIVATIVES of a_hs
    ! 1: a_hs
    if (GETDERAHS(1) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	    gl%ahsx3_PCSAFT(1,xi,xj,xk) = (18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - 24.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - 12.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) &
     & **2 - 6.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0)**3 + 18.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*( &
     & gl%z3_PCSAFT(1) - 1.d0)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 12.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 &
     & *(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) - 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)** &
     & 4*(gl%z3x1_PCSAFT(1,xj)*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) + gl%z3x1_PCSAFT(1,xi)*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) + gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**2)) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) &
     & **4*(gl%z3_PCSAFT(1) - 1.d0)**3 + 3.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) &
     & - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi) &
     & *gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)) + 2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) - 6.d0*gl%z0_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**5*(4.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) - &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3)*dlog(-gl%z3_PCSAFT(1) + 1.d0) &
     & - gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z0x1_PCSAFT(1,xj)*(-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**4*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & **3*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)) + 3.d0* &
     & gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk) &
     & *gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - gl%z3_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xi)*(-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)) + 3.d0* &
     & gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk) &
     & *gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - gl%z3_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xk)*(-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)) + 3.d0* &
     & gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi) &
     & *gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - gl%z3_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))))) - 2.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*( &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*( &
     & gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1))*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z2_PCSAFT(1) &
     & **2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0) &
     & **2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + (gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & (-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) &
     & - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + ( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0 &
     & *gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))*dlog(-gl%z3_PCSAFT(1) + 1.d0))) + 6.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*( &
     & -gl%z3_PCSAFT(1) + 1.d0) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)*dlog &
     & (-gl%z3_PCSAFT(1) + 1.d0)))/(gl%z0_PCSAFT(1)**4*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**5)
                    end if
        	    end do
            end do 
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, d_i, d_j
    if (GETDERAHS(2) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	    gl%ahsx3_PCSAFT(2,xi,xj,xk) = (6.d0*gl%z0_PCSAFT(1)**4*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 12.d0* &
     & gl%z0_PCSAFT(1)**4*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**4* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 72.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 72.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) &
     & *gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) &
     & *gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
     & - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) + 18.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1) &
     & **3*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) + 18.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1) &
     & **3*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1) &
     & **3*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 144.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 90.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 90.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 90.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) + 72.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 48.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 12.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) &
     & *gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 72.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) &
     & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 48.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) &
     & *gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 72.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) &
     & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 108.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 48.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) &
     & - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 6.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **4 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 12.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 30.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)** &
     & 2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + 54.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 24.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 54.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 24.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi) &
     & *gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 36.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1) &
     & **2*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) &
     & **3 - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 + &
     & 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)** &
     & 2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj) &
     & - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) - 30.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)** &
     & 3*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) + 54.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1) &
     & **2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 54.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
     & **4 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - &
     & 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 12.d0* &
     & gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj) &
     & - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - 30.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)** &
     & 3*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 54.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1) &
     & **2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) + 54.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - 24.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj) &
     & *gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1) - 4.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) + 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + 16.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk) - &
     & 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 6.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 24.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 36.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 24.d0* &
     & gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1) - 4.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 16.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi) - 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi) &
     & *gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)* &
     & gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)* &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1) - 4.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) + 24.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj) - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 16.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**4 - 36.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 72.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**2 - 60.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj) - &
     & 18.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 + 72.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 108.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 72.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 18.d0* &
     & gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) - 6.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4 + 36.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 + 60.d0*gl%z0x1_PCSAFT(1,xj)* &
     & gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) - 18.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi) &
     & *gl%z0x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)**3)/(gl%z0_PCSAFT(1)**4*(gl%z3_PCSAFT(1)**6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0* &
     & gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1)**2 - 6.d0*gl%z3_PCSAFT(1) + 1.d0))
                    end if
        	    end do
            end do 
        end do
    end if

    ! 4: 1ST DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(4) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	    gl%ahsx3_PCSAFT(3,xi,xj,xk) = (-72.d0*gl%z0_PCSAFT(1)**3*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) + 120.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 96.d0*gl%z0_PCSAFT(1)** &
     & 3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z3_PCSAFT(1) - 1.d0) + 72.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 48.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3 - 24.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) &
     & *gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 12.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) &
     & *gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) &
     & *gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 18.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1) &
     & *gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z3x1_PCSAFT(4,xj)*(gl%z2_PCSAFT(1)** &
     & 2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) + &
     & gl%z3x1_PCSAFT(4,xi)*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) + gl%z3x1_PCSAFT(4,xk)*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2)) + gl%z0_PCSAFT(1)**3* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z2_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) - 6.d0* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) - 6.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk)*(gl%z2_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2) + &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk))) + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(-6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi))) + gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)))) + 3.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk) + gl%z1x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(4,xi) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)) - 6.d0*gl%z0_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) &
     & *gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) &
     & *gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z1x1_PCSAFT(4,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)) + 18.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + &
     & 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z2x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(4,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) - 12.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) &
     & + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) - 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) + 2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - &
     & 1.d0)**3*(6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)** &
     & 2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
     & + 12.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 12.d0*gl%z2_PCSAFT(1) &
     & *gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk) &
     & *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj) - gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
     & - gl%z2_PCSAFT(1)**3) - gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)** &
     & 2 - gl%z2_PCSAFT(1)**3) - gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3)) + 6.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(-3.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
     & *(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - 2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*( &
     & gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + &
     & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + 6.d0* &
     & gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**5*(-gl%z3_PCSAFT(4)*(4.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*( &
     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) - gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3) + gl%z3x1_PCSAFT(1,xj)*(-4.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) &
     & + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk))) + gl%z3x1_PCSAFT(1,xi)*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk))) + gl%z3x1_PCSAFT(1,xk)*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)))) + gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(24.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) - 6.d0*gl%z0_PCSAFT(1)**3*( &
     & gl%z3_PCSAFT(1) - 1.d0)**6*(-20.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)* &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z3_PCSAFT(1)**4*( &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk) + gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) + 2.d0*gl%z3_PCSAFT(1)**3*( &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) &
     & + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)))*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0) - gl%z0_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z0x1_PCSAFT(1,xj)*(18.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 18.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - &
     & 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0 &
     & )**2 - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3 + 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - &
     & 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) - 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
     & + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(4)) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)) + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) - &
     & gl%z3x1_PCSAFT(1,xi)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xk))) - gl%z3x1_PCSAFT(1,xk)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)))) - 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1x1_PCSAFT(4,xi) &
     & *gl%z2x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xk) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(4,xi)) + 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) &
     & - 6.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(4)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + 6.d0*gl%z3_PCSAFT(1)**4 &
     & *(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)) + gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**3*(-3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(4,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) - 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) + &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - &
     & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3x1_PCSAFT(4,xi)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + &
     & gl%z3x1_PCSAFT(4,xk)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
     & ) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) &
     & *(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1))) + 6.d0*(gl%z3_PCSAFT(1) - 1.d0)**5*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)))*dlog(-gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xi)*(18.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 18.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - &
     & 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0 &
     & )**2 - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3 + 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - &
     & 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)* &
     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
     & + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(4)) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**2) - &
     & gl%z3x1_PCSAFT(1,xj)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xk))) - gl%z3x1_PCSAFT(1,xk)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)))) - 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1x1_PCSAFT(4,xj) &
     & *gl%z2x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xk) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(4,xj)) + 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + &
     & gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) &
     & - 6.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(4)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z3_PCSAFT(1)**4 &
     & *(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)) + gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**3*(-3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(4,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) + &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + &
     & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xk) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - &
     & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3x1_PCSAFT(4,xj)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + &
     & gl%z3x1_PCSAFT(4,xk)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
     & ) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) &
     & *(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + 6.d0*(gl%z3_PCSAFT(1) - 1.d0)**5*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)))*dlog(-gl%z3_PCSAFT(1) + 1.d0)) + gl%z0x1_PCSAFT(1,xk)*(18.d0* &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 18.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - &
     & 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0 &
     & )**2 - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3 + 6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - &
     & 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - 6.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + &
     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
     & + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) + gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)**2* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)) + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2) - &
     & gl%z3x1_PCSAFT(1,xj)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi))) - gl%z3x1_PCSAFT(1,xi)*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)))) - 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1x1_PCSAFT(4,xj) &
     & *gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(4,xi) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(4,xj)) + 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) &
     & - 6.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(4)* &
     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) + 6.d0*gl%z3_PCSAFT(1)**4 &
     & *(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)) + gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**3*(-3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3x1_PCSAFT(4,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) - 6.d0*gl%z2_PCSAFT(1)* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3x1_PCSAFT(1,xj) - 6.d0*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + &
     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + &
     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) - &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) - &
     & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3x1_PCSAFT(4,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + &
     & gl%z3x1_PCSAFT(4,xi)*(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1) &
     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
     & ) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z0x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) &
     & *(gl%z0x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1))) + 6.d0*(gl%z3_PCSAFT(1) - 1.d0)**5*(-4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
     & gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(4)* &
     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)) - gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) &
     & + 2.d0*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + gl%z2_PCSAFT(4) &
     & *gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)))*dlog(-gl%z3_PCSAFT(1) + 1.d0))) - 2.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) &
     & *(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk)*(gl%z3_PCSAFT(1) - &
     & 1.d0)**2 + 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + &
     & 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)) &
     & + gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + 2.d0* &
     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk)))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xk) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(4,xk)* &
     & gl%z2_PCSAFT(1) + gl%z1x1_PCSAFT(1,xk)*gl%z2_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*( &
     & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xk) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + gl%z1_PCSAFT(1)* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xk) + gl%z1x1_PCSAFT(1,xk)* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xk)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xk) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) + 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xk)*gl%z3_PCSAFT(4) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xk)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xk) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xk)* &
     & gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xk)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xi) - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)) + gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2x1_PCSAFT(1,xi) &
     & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xi)))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(4,xi)*gl%z2_PCSAFT(1) + &
     & gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(4,xi) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%z1x1_PCSAFT(1,xi)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xi)* &
     & gl%z3_PCSAFT(1)))) + gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xj) - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - &
     & 1.d0) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0* &
     & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4) &
     & *gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)) - 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & - 1.d0)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)) + gl%z2_PCSAFT(1)*( &
     & gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2x1_PCSAFT(1,xj) &
     & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + 2.d0*gl%z2_PCSAFT(4)* &
     & gl%z2x1_PCSAFT(1,xj)))*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
     & gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj) + gl%z1_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(4,xj)*gl%z2_PCSAFT(1) + &
     & gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) &
     & *gl%z3x1_PCSAFT(4,xj) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%z1_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(4) + gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%z1x1_PCSAFT(1,xj)*gl%z2_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)) + gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z2_PCSAFT(1)**3* &
     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z2_PCSAFT(1)**2* &
     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) - gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
     & gl%z2_PCSAFT(1)**3)) + gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*(2.d0* &
     & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1)**3 + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xj) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2x1_PCSAFT(1,xj)* &
     & gl%z3_PCSAFT(1))))) + 6.d0*gl%z0x1_PCSAFT(1,xj)*gl%z0x1_PCSAFT(1,xi)*gl%z0x1_PCSAFT(1,xk)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) &
     & - 1.d0)**3*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0) + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) &
     & - 1.d0) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) + gl%z2_PCSAFT(1)**2*( &
     & gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) + &
     & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1) &
     & **2 - gl%z2_PCSAFT(1)**3)))/(gl%z0_PCSAFT(1)**4*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**6)
                    end if
        	    end do
            end do 
        end do
    end if
    

    !DEC$ END IF
end subroutine AHSX3DERIVS
    




    end module pc_saft_AHSX_derivs_module
