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

    ! module for file pc-saft_GIIX_derivs.f90
    module pc_saft_GIIX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains


 

subroutine GIIX1DERIVS(gl,GETDERGII)

! Henning Markgraf, June 2016

    ! g_ij_hs: site-site radial distribution function of hard sphere fluid
    ! defined by eq. A.7 in Gross, Sadowski 2001:
    ! g_ij_hs = 1/(1-zeta_3) + (d_i*d_j/(d_i+d_j)) * 3*zeta_2/(1-zeta_3)**2
    !           + (d_i*d_j/(d_i+d_j))**2 * 2*zeta_2**2/(1-zeta_3)**3
    ! dependent on D and T

!--------------------------------------------------------------------------------------------------
! All first composition derivatives of g_ii_hs are calculated in this subroutine
!--------------------------------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERGII
    !output: giix1_PCSAFT (module variable)
    !working variables
    integer :: nrsubst, xi

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. X1 DERIVATIVES
    ! calculate gii_x for all components nrsubst
    do nrsubst = 1, gl%ncomp
        !1: g_ii_hs
        if (GETDERGII(1) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,1,xi) = (3.d0/2.d0)*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)/(-gl%z3_PCSAFT(1) + 1.d0)**4 + gl%di_PCSAFT(nrsubst,1)** &
		         & 2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)/(-gl%z3_PCSAFT(1) + 1.d0)**3 + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)/( &
		         & -gl%z3_PCSAFT(1) + 1.d0)**3 + (3.d0/2.d0)*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + &
		         & gl%z3x1_PCSAFT(1,xi)/(-gl%z3_PCSAFT(1) + 1.d0)**2
            end do
        end if
    
        ! 2: 1ST DERIVATIVE OF g_ii_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
        ! requires di_PCSAFT, d_j, z2_PCSAFT(1), z3_PCSAFT(1)
        if (GETDERGII(2) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,2,xi) = (1.d0/2.d0)*(-12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(1) + 6.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi) &
		         & + 2.d0*gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + &
		         & 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**5
            end do
        end if
    
        ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j
        if (GETDERGII(3) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,3,xi) = (30.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
		         & gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 3.d0* &
		         & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) &
		         & + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + &
		         & 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) &
		         & ))/(gl%z3_PCSAFT(1) - 1.d0)**6
            end do
        end if
    
        ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), di_PCSAFT(nrsubst,4), d_j_T
        if (GETDERGII(4) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,4,xi) = (1.d0/2.d0)*(-12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0 &
		         & )**2*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 2.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**5
            end do
        end if
    
        ! 5: 2ND DERIVATIVE OF g_ij WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), di_PCSAFT(nrsubst,4), d_j_T, z3_PCSAFT(5), z2_PCSAFT(5), di_PCSAFT(nrsubst,5), d_j_TT
        if (GETDERGII(5) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,5,xi) = (1.d0/2.d0)*(60.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) - 12.d0*gl%di_PCSAFT(nrsubst,1) &
		         & *gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) &
		         & + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**4*(3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(5,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2x1_PCSAFT(1,xi) + 2.d0* &
		         & gl%z3x1_PCSAFT(5,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi) + 2.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 4.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + &
		         & 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 6.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(nrsubst,4)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4) &
		         & *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)) + 3.d0 &
		         & *(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) &
		         & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 4.d0* &
		         & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1) &
		         & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 2.d0* &
		         & gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**6
            end do
        end if
    
        ! 6: 2ND MIXED DERIVATIVE OF g_ij WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), di_PCSAFT(nrsubst,4), d_j_T
        if (GETDERGII(6) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,6,xi) = (1.d0/2.d0)*(60.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 12.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**4*(3.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)) - 2.d0*(gl%z3_PCSAFT(1) &
		         & - 1.d0)**3*(2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)* &
		         & gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
		         & 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 4.d0*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + &
		         & 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) &
		         & *gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) &
		         & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) &
		         & + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)** &
		         & 6
            end do
        end if
    
        ! 7: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO rho, T, AND T, MULTIPLIED BY T*T*rho
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), z3_PCSAFT(5), z2_PCSAFT(5), di_PCSAFT(nrsubst,4), d_j_T, di_PCSAFT(nrsubst,5), d_j_TT
        if (GETDERGII(7) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,7,xi) = (1.d0/2.d0)*(-240.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) &
		         & + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 25.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 12.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2* &
		         & gl%z3x1_PCSAFT(1,xi)) + 2.d0*gl%z3_PCSAFT(1)*(-60.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2* &
		         & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
		         & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)) - (gl%z3_PCSAFT(1) &
		         & - 1.d0)**4*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(5,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,5)* &
		         & gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(5,xi)) + 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(5,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)* &
		         & gl%z2x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)) - 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2 &
		         & *gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + &
		         & 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0* &
		         & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1) &
		         & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4) &
		         & *gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 2.d0*gl%di_PCSAFT(nrsubst,4)**2* &
		         & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 4.d0* &
		         & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi))) + (gl%z3_PCSAFT(1) - 1.d0)**5*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(5,xi) + &
		         & 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(5,xi)) - 4.d0*( &
		         & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(4,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)* &
		         & gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)) + ( &
		         & gl%z3_PCSAFT(1) - 1.d0)**3*(9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 36.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
		         & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(1) + 36.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 36.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0* &
		         & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 36.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2* &
		         & gl%z3x1_PCSAFT(4,xi) + 72.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 72.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5) &
		         & *gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 108.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 54.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 108.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + &
		         & 54.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) &
		         & + 4.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(4,xi) + 108.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
		         & + 6.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%z3_PCSAFT(4)**2* &
		         & gl%z3x1_PCSAFT(1,xi)) - 6.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(5,xi) + 16.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)** &
		         & 2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
		         & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 32.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1) &
		         & **2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 16.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(4)**2 + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2* &
		         & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 32.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + &
		         & 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) &
		         & *gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 48.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2* &
		         & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 2.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 4.d0* &
		         & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**7
            end do
        end if
    
        ! 8: 3RD DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j
        if (GETDERGII(8) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,8,xi) = 3.d0*(-60.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(1,xi) + 10.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)*(9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi)) - 4.d0*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) &
		         & **2*(9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(1) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1) &
		         & **2*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**7
            end do
        end if
    
        ! 9: 3RD DERIVATIVE OF g_ij WITH RESPECT TO T, MULTIPLIED BY T^3
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), z3_PCSAFT(5), z2_PCSAFT(5),
        ! z3_PCSAFT(9), z2_PCSAFT(9), di_PCSAFT(nrsubst,4), d_j_T, di_PCSAFT(nrsubst,5), d_j_TT, di_PCSAFT(nrsubst,9), d_j_TTT
        if (GETDERGII(9) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,9,xi) = (1.d0/2.d0)*(-360.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi) + 60.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
		         & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) &
		         & - 1.d0)**5*(3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(9,xi) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(5,xi) + 9.d0*gl%di_PCSAFT(nrsubst,5)* &
		         & gl%z2x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,9)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(9,xi)) - 2.d0*(gl%z3_PCSAFT(1) &
		         & - 1.d0)**4*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(9,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(5,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(9) &
		         & *gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 2.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,9)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(9,xi) + 9.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 9.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(5) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 6.d0*gl%di_PCSAFT(nrsubst,4)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4) &
		         & *gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 9.d0*gl%di_PCSAFT(nrsubst,5)* &
		         & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,5)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,9)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) &
		         & ) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(9,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
		         & gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
		         & gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(9) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
		         & + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0 &
		         & *gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 24.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2* &
		         & gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,9)*gl%z2_PCSAFT(1)**2* &
		         & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 36.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)* &
		         & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
		         & 6.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 36.d0*gl%di_PCSAFT(nrsubst,4)* &
		         & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + &
		         & 18.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(4,xi) &
		         & + 12.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)) - 12.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)** &
		         & 2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)* &
		         & gl%z3x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
		         & gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)**2* &
		         & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + &
		         & 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1) &
		         & **2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2* &
		         & gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)** &
		         & 3 + 6.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0) &
		         & **7
            end do
        end if
    
        ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
        ! requires z2_PCSAFT(1), z3_PCSAFT(1), di_PCSAFT, d_j, z3_PCSAFT(4), z2_PCSAFT(4), di_PCSAFT(nrsubst,4), d_j_T
        if (GETDERGII(10) .eq. 1) then
            do xi = 1, gl%ncomp
                gl%giix1_PCSAFT(nrsubst,10,xi) = (-180.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 30.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1) &
		         & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + &
		         & 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
		         & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi) &
		         & *gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
		         & + 2.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0) &
		         & **3*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) &
		         & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 6.d0* &
		         & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		         & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		         & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2 + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 2.d0 &
		         & *gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) + 8.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) - 12.d0*( &
		         & gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
		         & gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		         & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
		         & gl%z3_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
		         & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) &
		         & + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
		         & gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 3.d0* &
		         & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
		         & gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**7
            end do
        end if
    end do
    
!DEC$ END IF
end subroutine GIIX1DERIVS
         
    
subroutine GIIX2DERIVS(gl,GETDERGII)

! Henning Markgraf, June 2016

    ! g_ij_hs: site-site radial distribution function of hard sphere fluid
    ! defined by eq. A.7 in Gross, Sadowski 2001:
    ! g_ij_hs = 1/(1-zeta_3) + (d_i*d_j/(d_i+d_j)) * 3*zeta_2/(1-zeta_3)**2
    !           + (d_i*d_j/(d_i+d_j))**2 * 2*zeta_2**2/(1-zeta_3)**3
    ! dependent on D and T

!--------------------------------------------------------------------------------------------------
! All second composition derivatives of g_ii_hs are calculated in this subroutine
!--------------------------------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERGII
    !output: giix2_PCSAFT
    !working variables
    integer :: nrsubst, xi, xj
    
    double precision :: help1,help2,help3,help4,help5,help6,help7,help8,help9
    double precision :: z31,z312,z313,z314,z315,z316,z317
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. X2 DERIVATIVES of g_ii_hs for all components i
    ! calculate gii_x for all components i
    do nrsubst = 1, gl%ncomp
	    ! 1: d^2(g_ii_hs)/(dxi dxj)
	    if (GETDERGII(1) .eq. 1) then
			do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
				gl%giix2_PCSAFT(nrsubst,1,xi,xj) = (-6.d0*gl%di_PCSAFT(nrsubst,1)**2 *gl%z2_PCSAFT(1)**2 *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) & 
                    & + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                    & + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) &
                    & - (gl%z3_PCSAFT(1) - 1.d0)**2 *(gl%di_PCSAFT(nrsubst,1)**2 *gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                    & + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**5
                    end if
				end do
			end do
	    end if
	    
	    ! 2: 1ST DERIVATIVE OF g_ii_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
	    ! d^3(g_ii_hs)/(dxi dxj drho)
	    if (GETDERGII(2) .eq. 1) then
            z31 = (gl%z3_PCSAFT(1) - 1.d0)
            z312 = z31*z31
            z313 = z312*z31
            z316 = z313*z313
            help4 = gl%di_PCSAFT(nrsubst,1)**2 * gl%z2_PCSAFT(1)
            help8 = gl%di_PCSAFT(nrsubst,1) * gl%z2_PCSAFT(1)
		    do xj = 1, gl%ncomp
		        do xi = 1, gl%ncomp
                    if (xi .GE. xj) then
                    help1 = gl%di_PCSAFT(nrsubst,1)**2 * gl%z2x1_PCSAFT(1,xj) * gl%z2x1_PCSAFT(1,xi)
                    help2 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)
                    help3 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)
                    help5 = gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)
                    help6 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xj) * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xi)
                    help7 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xi) * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xj)
                
			    gl%giix2_PCSAFT(nrsubst,2,xi,xj) = (18.d0 * help4 * gl%z2_PCSAFT(1) * gl%z3_PCSAFT(1) * help5 &
                    & - 6.d0 * help8 * z31 * (4.d0 * help8 * help5 + help6 + help7 + 3.d0 * gl%z3_PCSAFT(1) * help5) &
                    & + 2.d0 * gl%z3_PCSAFT(1) * (6.d0 * help4 * gl%z2_PCSAFT(1) * help5 - 3.d0 * help8 * z31 * (help2 + help3 + 3.d0 * help5) + z312 * (help1 + 3.d0 * help2 + 3.d0 * help3 + 2.d0 * help5)) &
                    & - 2.d0 * z313 * (help1 + 3.d0 * help2 + 3.d0 * help3 + 2.d0 * help5) &
                    & + z312 * (9.d0 * help4 * gl%z2x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi) &
                    & + 9.d0 * help4 * gl%z2x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) + help1 * gl%z3_PCSAFT(1) + 27.d0*help8*help5 + 3.d0*help6 + 3.d0*help7 + 2.d0*gl%z3_PCSAFT(1) * help5))/z316
                    end if
			    end do
		    end do
	    end if
	    
	    ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
	    ! d^4(g_ii_hs)/(dxi dxj drho^2)
	    if (GETDERGII(3) .eq. 1) then
        
            help1 = gl%di_PCSAFT(nrsubst,1)**2 * gl%z2_PCSAFT(1)
            z312 = gl%z3_PCSAFT(1)*gl%z3_PCSAFT(1)
            z313 = z312*gl%z3_PCSAFT(1)
            z314 = z313*gl%z3_PCSAFT(1)
            z315 = z314*gl%z3_PCSAFT(1)
            z316 = z315*gl%z3_PCSAFT(1)
            z317 = z316*gl%z3_PCSAFT(1)
            
		    do xj = 1, gl%ncomp
		        do xi = 1, gl%ncomp
                    if (xi .GE. xj) then
                    
                    help2 = help1 * gl%z2_PCSAFT(1) * gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi) !di_PCSAFT(nrsubst,1)**2 * z2_PCSAFT(1)**2 * z3x1_PCSAFT(1,xj) * z3x1_PCSAFT(1,xi)
                    help3 = help1 * gl%z2x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)    !di_PCSAFT(nrsubst,1)**2 * z2_PCSAFT(1) * z2x1_PCSAFT(1,xj) * z3x1_PCSAFT(1,xi)
                    help4 = help1 * gl%z2x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)  !di_PCSAFT(nrsubst,1)**2 * z2_PCSAFT(1) * z2x1_PCSAFT(1,xi) * z3x1_PCSAFT(1,xj)
                    help5 = gl%di_PCSAFT(nrsubst,1)**2 * gl%z2x1_PCSAFT(1,xj) * gl%z2x1_PCSAFT(1,xi)
                    help6 = gl%di_PCSAFT(nrsubst,1) * gl%z2_PCSAFT(1) * gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)
                    help7 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)
                    help8 = gl%di_PCSAFT(nrsubst,1) * gl%z2x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)
                    help9 = gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi)
                    
			    gl%giix2_PCSAFT(nrsubst,3,xi,xj) = ( &
                    & - 12.d0 * help2 * z312 - 96.d0 * help2 * gl%z3_PCSAFT(1) - 72.d0 * help2 &
                    &  + 6.d0 * help3 * z313 + 30.d0 * help3 * z312 - 18.d0 * help3 * gl%z3_PCSAFT(1) - 18.d0 * help3 &
                    &  + 6.d0 * help4 * z313 + 30.d0 * help4 * z312 - 18.d0 * help4 * gl%z3_PCSAFT(1) - 18.d0 * help4 &
                    &  - 2.d0 * help5 * z314 - 4.d0 * help5 * z313 + 12.d0 * help5 * z312 - 4.d0 * help5 * gl%z3_PCSAFT(1) - 2.d0 * help5 &
                    & + 18.d0 *help6 * z313 + 90.d0 *help6 * z312 - 54.d0 *help6 * gl%z3_PCSAFT(1) - 54.d0 *help6 &
                    &  - 6.d0 * help7 * z314 - 12.d0 * help7 * z313 + 36.d0 * help7 * z312 - 12.d0 * help7 * gl%z3_PCSAFT(1) - 6.d0 * help7 &
                    &  - 6.d0 * help8 * z314 - 12.d0 * help8 * z313 + 36.d0 * help8 * z312 - 12.d0 * help8 * gl%z3_PCSAFT(1) - 6.d0 * help8 &
                    &  - 4.d0 * z314 * help9 - 8.d0 * z313 * help9 + 24.d0 * z312 * help9 - 8.d0 * gl%z3_PCSAFT(1) * help9 - 4.d0 * help9) &
                    & / (z317 - 7.d0 * z316 + 21.d0 * z315 - 35.d0 * z314 + 35.d0 * z313 - 21.d0 * z312 + 7.d0 * gl%z3_PCSAFT(1) - 1.d0)
                    end if
		        end do
		    end do
	    end if
	    
	    ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
	    ! d^3(g_ii_hs)/(dxi dxj dT)
	    if (GETDERGII(4) .eq. 1) then
		do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
			gl%giix2_PCSAFT(nrsubst,4,xi,xj) = (18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1) &
	     & *(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 2.d0*gl%z3_PCSAFT(4)*(6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj) &
	     & *gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & ) + (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) - (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xj)* &
	     & gl%z2x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)* &
	     & gl%z3x1_PCSAFT(1,xj)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4) &
	     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**6
                end if
			end do
		end do
	    end if
	    
	    ! 5: 2ND DERIVATIVE OF g_ij WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
	    ! d^4(g_ii_hs)/(dxi dxj dT^2)
	    if (GETDERGII(5) .eq. 1) then
		do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
			gl%giix2_PCSAFT(nrsubst,5,xi,xj) = (-72.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
	     & gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + 6.d0* &
	     & gl%z3_PCSAFT(4)**2*(-6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) + 4.d0*gl%z3_PCSAFT(4)*( &
	     & -18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi) &
	     & *gl%z3x1_PCSAFT(1,xj)) - (gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4) &
	     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) + 2.d0*gl%z3_PCSAFT(5)*( &
	     & gl%z3_PCSAFT(1) - 1.d0)*(6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) - (gl%z3_PCSAFT(1) - 1.d0) &
	     & **4*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(5,xj)*gl%z2x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(5,xi)* &
	     & gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(4,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + &
	     & 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(5,xj) + 2.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + &
	     & 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + &
	     & 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + &
	     & 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,5)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
	     & 2.d0*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 4.d0* &
	     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) &
	     & + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi) &
	     & + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(5,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
	     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) &
	     & + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(5)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xj) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) + &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 4.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(4,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)* &
	     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,4)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)**2 &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 18.d0*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,5)* &
	     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 4.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & + 4.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)** &
	     & 2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0* &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(5)* &
	     & gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) &
	     & **2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4) &
	     & *gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2 + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2 &
	     & *gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
	     & 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,5)*gl%z2_PCSAFT(1)**2* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)**2* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj) + 6.d0*gl%di_PCSAFT(nrsubst,4)**2 &
	     & *gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))/( &
	     & gl%z3_PCSAFT(1) - 1.d0)**7
                end if
			end do
		end do
	    end if
	    
	    ! 6: 2ND MIXED DERIVATIVE OF g_ij WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
	    ! d^4(g_ii_hs)/(dxi dxj dT drho)
	    if (GETDERGII(6) .eq. 1) then
		do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
			gl%giix2_PCSAFT(nrsubst,6,xi,xj) = (-72.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & + gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 5.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) &
	     & *gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) &
	     & - 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)** &
	     & 2*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) + 2.d0*gl%z3_PCSAFT(1) &
	     & *(-18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi) &
	     & *gl%z3x1_PCSAFT(1,xj)) - (gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4) &
	     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) + 2.d0*gl%z3_PCSAFT(4)*( &
	     & gl%z3_PCSAFT(1) - 1.d0)*(6.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi))) - 2.d0*gl%z3_PCSAFT(4)*( &
	     & 18.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*( &
	     & gl%z3_PCSAFT(1) - 1.d0)*(4.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)* &
	     & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)* &
	     & gl%z2x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(9.d0*gl%di_PCSAFT(nrsubst,1)**2 &
	     & *gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + 27.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
	     & )) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi) + gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1) &
	     & *gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 3.d0*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0* &
	     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)) + (gl%z3_PCSAFT(1) - 1.d0) &
	     & **3*(9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + &
	     & 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(4,xi)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4) + 18.d0* &
	     & gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 18.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1) &
	     & + 27.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 27.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 27.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
	     & 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 27.d0*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 2.d0*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 6.d0* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(12.d0*gl%di_PCSAFT(nrsubst,1)**2* &
	     & gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2* &
	     & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj) + 24.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
	     & gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0* &
	     & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(4,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
	     & gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xj) + 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + 3.d0 &
	     & *gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) &
	     & *gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%di_PCSAFT(nrsubst,1)**2*gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
	     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + &
	     & 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
	     & gl%z2_PCSAFT(1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
	     & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xi)* &
	     & gl%z3x1_PCSAFT(1,xj) + 36.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 9.d0*gl%di_PCSAFT(nrsubst,1) &
	     & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xj)*gl%z3_PCSAFT(1)* &
	     & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2x1_PCSAFT(1,xi)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) &
	     & + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
	     & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))/(gl%z3_PCSAFT(1) - 1.d0)**7
                end if
			end do
		end do
	    end if
    end do
    
    !DEC$ END IF
end subroutine GIIX2DERIVS


subroutine GIIX3DERIVS(gl,GETDERGII)

! Henning Markgraf, June 2016

    ! g_ij_hs: site-site radial distribution function of hard sphere fluid
    ! defined by eq. A.7 in Gross, Sadowski 2001:
    ! g_ij_hs = 1/(1-zeta_3) + (d_i*d_j/(d_i+d_j)) * 3*zeta_2/(1-zeta_3)**2
    !           + (d_i*d_j/(d_i+d_j))**2 * 2*zeta_2**2/(1-zeta_3)**3
    ! dependent on D and T

!--------------------------------------------------------------------------------------------------
! All third composition derivatives of g_ii_hs are calculated in this subroutine
!--------------------------------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERGII
    !output: giix3_PCSAFT (module variable)
    !working variables
    integer :: nrsubst, xi, xj, xk
    double precision :: part1, part2, part3, part4, part5, part6, part7, part8
    double precision :: part9, part10, part11, part12, part13, part14, part15, part16
    double precision :: part17, part18, part19, part20, part21, part22, part23, part24
    double precision :: part25, part26, part27, part28, part29
  
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. X3 DERIVATIVES of g_ii for all components i
	do nrsubst = 1, gl%ncomp
		!1: d^3(g_ii_hs)/(dxi dxj dxk)
		if (GETDERGII(1) .eq. 1) then
		part1 = 3.d0/(gl%z3_PCSAFT(1) - 1.d0)**6
		part2 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 - 2.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1) + gl%di_PCSAFT(nrsubst,1)**2)
		part3 = (gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 - 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) + gl%di_PCSAFT(nrsubst,1)**2)
		part4 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%z3_PCSAFT(1)**2 - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1))
		part5 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 - 2.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) + gl%di_PCSAFT(nrsubst,1)**2)
		part6 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 4.d0*gl%di_PCSAFT(nrsubst,1) &
		 & **2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1))
		part7 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) &
		 & + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1))
		part8 = (10.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**2 - 4.d0*gl%z3_PCSAFT(1) + 2.d0)
		    do xk = 1 , gl%ncomp
                do xj = 1 , gl%ncomp
                    do xi = 1 , gl%ncomp
                        if (xi .GE. xj .AND. xj .GE. xk) then
		        	    gl%giix3_PCSAFT(nrsubst,1,xi,xj,xk) = part1*(gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part2 + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)*part3 + gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
		 & gl%z3x1_PCSAFT(1,xk)*part4 + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)* &
		 & gl%z3x1_PCSAFT(1,xj)*part5 + &
		 & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part6 + gl%z2x1_PCSAFT(1,xk) &
		 & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part7 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
		 & *gl%z3x1_PCSAFT(1,xk)*part8)
                        end if
		    	    end do
		        end do 
		    end do
		end if
		
		! 2: 1ST DERIVATIVE OF g_ij_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
		! d^4(g_ii_hs)/(dxi dxj dxk drho)
		if (GETDERGII(2) .eq. 1) then
		part1 = -3.d0/(gl%z3_PCSAFT(1) - 1.d0)**7
		part2 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 + gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 - 5.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2)
		part3 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 + gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 - 5.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) &
		 & + 3.d0*gl%di_PCSAFT(nrsubst,1)**2)
		part4 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		 & gl%z3_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 16.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%z3_PCSAFT(1)**3 + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 - 15.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 9.d0*gl%di_PCSAFT(nrsubst,1))
		part5 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 + gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 - 5.d0* &
		 & gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)**2)
		part6 = (-4.d0* &
		 & gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 16.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 - 15.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 9.d0* &
		 & gl%di_PCSAFT(nrsubst,1))
		part7 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
		 & **2 - 12.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 16.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)** &
		 & 3 + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 - 15.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) + 9.d0*gl%di_PCSAFT(nrsubst,1))
		part8 = (10.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 50.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 36.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 48.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**3 + 2.d0*gl%z3_PCSAFT(1)**2 - 10.d0*gl%z3_PCSAFT(1) + 6.d0)
			do xk = 1 , gl%ncomp
                do xj = 1 , gl%ncomp
                    do xi = 1 , gl%ncomp
                        if (xi .GE. xj .AND. xj .GE. xk) then 
		        	    gl%giix3_PCSAFT(nrsubst,2,xi,xj,xk) = part1*(gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part2 + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)* &
		 & gl%z3x1_PCSAFT(1,xi)*part3 + gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part4 + gl%z2x1_PCSAFT(1,xi) &
		 & *gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*part5 + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part6 + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part7 + gl%z3x1_PCSAFT(1,xj)* &
		 & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part8)
                        end if
		    	    end do
		        end do 
		    end do
		end if
		
		! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
		! d^4(g_ii_hs)/(dxi dxj dxk dT)
		if (GETDERGII(4) .eq. 1) then
		part1 = 3.d0/(gl%z3_PCSAFT(1) - 1.d0)**7
		part2 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part3 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part4 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part5 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part6 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part7 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part8 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part9 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part10 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part11 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part12 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1) - 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4))
		part13 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 + 3.d0 &
		 & *gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part14 = (-4.d0* &
		 & gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1) - 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4))
		part15 = ( &
		 & -4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1) &
		 & )
		part16 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
		 & + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part17 = (20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
		 & - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 16.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
		 & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
		 & 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 9.d0* &
		 & gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,4))
		part18 = (gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**3 - 3.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2 + 3.d0 &
		 & *gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1) - gl%di_PCSAFT(nrsubst,1)**2)
		part19 = (-4.d0* &
		 & gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 6.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(nrsubst,1)* &
		 & gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1) - 2.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4))
		part20 = ( &
		 & -4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1) &
		 & )
		part21 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
		 & + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part22 = (20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
		 & - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 16.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4) &
		 & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
		 & 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 9.d0* &
		 & gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,4))
		part23 = (-4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) &
		 & **2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1))
		part24 = ( &
		 & -4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,1) &
		 & )
		part25 = (20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)** &
		 & 2 + 8.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 4.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(4) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		 & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 16.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 8.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		 & gl%z2_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0* &
		 & gl%di_PCSAFT(nrsubst,1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%di_PCSAFT(nrsubst,4)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(nrsubst,4)* &
		 & gl%z3_PCSAFT(1) - 3.d0*gl%di_PCSAFT(nrsubst,4))
		part26 = (10.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 10.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
		 & + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**3 - 6.d0*gl%z3_PCSAFT(1)** &
		 & 2 + 6.d0*gl%z3_PCSAFT(1) - 2.d0)
		part27 = (10.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 10.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
		 & + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**3 - 6.d0*gl%z3_PCSAFT(1)** &
		 & 2 + 6.d0*gl%z3_PCSAFT(1) - 2.d0)
		part28 = (10.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 10.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
		 & + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**3 - 6.d0*gl%z3_PCSAFT(1)** &
		 & 2 + 6.d0*gl%z3_PCSAFT(1) - 2.d0)
		part29 = (-60.d0*gl%di_PCSAFT(nrsubst,1)**2* &
		 & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 20.d0*gl%di_PCSAFT(nrsubst,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 20.d0*gl%di_PCSAFT(nrsubst,1)**2 &
		 & *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 20.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 20.d0*gl%di_PCSAFT(nrsubst,1)*gl%di_PCSAFT(nrsubst,4)* &
		 & gl%z2_PCSAFT(1)**2 + 60.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 60.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 24.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - &
		 & 12.d0*gl%di_PCSAFT(nrsubst,1)*gl%z2_PCSAFT(4) - 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 24.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1)* &
		 & gl%z3_PCSAFT(1) - 12.d0*gl%di_PCSAFT(nrsubst,4)*gl%z2_PCSAFT(1) - 8.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 16.d0*gl%z3_PCSAFT(1)* &
		 & gl%z3_PCSAFT(4) - 8.d0*gl%z3_PCSAFT(4))
		    do xk = 1 , gl%ncomp
                do xj = 1 , gl%ncomp
                    do xi = 1 , gl%ncomp
                        if (xi .GE. xj .AND. xj .GE. xk) then 
		        	    gl%giix3_PCSAFT(nrsubst,3,xi,xj,xk) = part1*(gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part2 + gl%z2x1_PCSAFT(4,xj)*gl%z2x1_PCSAFT(1,xk)* &
		 & gl%z3x1_PCSAFT(1,xi)*part3 + gl%z2x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part4 + &
		 & gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part5 + gl%z2x1_PCSAFT(4,xi)*gl%z2x1_PCSAFT(1,xk)* &
		 & gl%z3x1_PCSAFT(1,xj)*part6 + gl%z2x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part7 + &
		 & gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part8 + gl%z2x1_PCSAFT(4,xk)*gl%z2x1_PCSAFT(1,xi)* &
		 & gl%z3x1_PCSAFT(1,xj)*part9 + gl%z2x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part10 + &
		 & gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk)*part11 + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xi)* &
		 & gl%z3x1_PCSAFT(1,xk)*part12 + gl%z2x1_PCSAFT(1,xj)* &
		 & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)*part13 + gl%z2x1_PCSAFT(1,xj)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xi)*part14 + &
		 & gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk)*part15 + gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xi)*part16 + gl%z2x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
		 & gl%z3x1_PCSAFT(1,xk)*part17 + gl%z2x1_PCSAFT(1,xi)* &
		 & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xj)*part18 + gl%z2x1_PCSAFT(1,xi)*gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*part19 + &
		 & gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk)*part20 + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*part21 + gl%z2x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xj)* &
		 & gl%z3x1_PCSAFT(1,xk)*part22 + gl%z2x1_PCSAFT(1,xk)* &
		 & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*part23 + gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*part24 + &
		 & gl%z2x1_PCSAFT(1,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part25 + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part26 + &
		 & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part27 + gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part28 + &
		 & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part29)
                        end if
		    	    end do
		        end do 
		    end do
		end if
	end do

    !DEC$ END IF
end subroutine GIIX3DERIVS



    end module pc_saft_GIIX_derivs_module
