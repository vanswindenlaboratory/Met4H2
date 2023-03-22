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

    ! module for file pc-saft_AHCX_derivs.f90
    module pc_saft_AHCX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains


    
subroutine AHCX1DERIVS(gl,GETDERAHC)

! Henning Markgraf, June 2016

    ! a_hc: Hard chain reference contribution to the Helmholtz free energy
    ! defined by eq. A.4 in Gross, Sadowski 2001:
    ! a_hc = mmean*a_hs - Sum( x_i*(m_i-1)*ln(g_ii_hs) )
    ! dependent on D and T

!--------------------------------------------------------------------------------------------------
! All first composition derivatives of g_ij_hs are calculated in this subroutine
!--------------------------------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAHC
    !output: ahcx1_PCSAFT
    !working variables
    double precision :: sum
    integer :: i, xi
    integer:: errorfld
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. X1 DERIVATIVES
    !calculate the derivatives of a_hc
    ! 1: a_hc
    if (GETDERAHC(1) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + gl%giix1_PCSAFT(i,1,xi)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
            end do
            gl%ahcx1_PCSAFT(1,xi) = gl%ahs_PCSAFT(1)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(1,xi)*gl%mmean_PCSAFT - (gl%mPCSAFT(xi) - 1.d0)*dlog(gl%gii_PCSAFT(xi,1)) - sum
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires ...
    if (GETDERAHC(2) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%gii_PCSAFT(i,1)*gl%giix1_PCSAFT(i,2,xi) - gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,1,xi))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
            end do
            gl%ahcx1_PCSAFT(2,xi) = gl%ahs_PCSAFT(2)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(2,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,2)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) - sum
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_D, gii_PCSAFT_hs_DD, ahs_PCSAFT_DD
    if (GETDERAHC(3) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%giix1_PCSAFT(i,3,xi) - 2.d0*gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,2,xi)/gl%gii_PCSAFT(i,1) - &
				& gl%gii_PCSAFT(i,3)*gl%giix1_PCSAFT(i,1,xi)/gl%gii_PCSAFT(i,1) + 2.d0*gl%gii_PCSAFT(i,2)**2*gl%giix1_PCSAFT(i,1,xi) &
				& /gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
            end do
            gl%ahcx1_PCSAFT(3,xi) = gl%ahs_PCSAFT(3)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(3,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,3)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) &
			& + gl%gii_PCSAFT(xi,2)**2*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - sum
        end do
    end if

    ! 4: 1ST DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires ...
    if (GETDERAHC(4) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%gii_PCSAFT(i,1)*gl%giix1_PCSAFT(i,4,xi) - gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,1,xi))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
            end do
            gl%ahcx1_PCSAFT(4,xi) = gl%ahs_PCSAFT(4)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,4)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) - sum
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_T, gii_PCSAFT_hs_TT, ahs_PCSAFT_TT
    if (GETDERAHC(5) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%giix1_PCSAFT(i,5,xi) - 2.d0 &
			     & *gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,4,xi)/gl%gii_PCSAFT(i,1) - gl%gii_PCSAFT(i,5)*gl%giix1_PCSAFT(i,1,xi)/gl%gii_PCSAFT(i,1) &
			     & + 2.d0*gl%gii_PCSAFT(i,4)**2*gl%giix1_PCSAFT(i,1,xi)/gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
            end do
            gl%ahcx1_PCSAFT(5,xi) = gl%ahs_PCSAFT(5)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(5,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,5)*(gl%mPCSAFT(xi) - 1.d0) &
			& /gl%gii_PCSAFT(xi,1) + gl%gii_PCSAFT(xi,4)**2*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - sum
        end do
    end if


    ! 6: 2ND MIXED DERIVATIVE OF a_hc WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_T, gii_PCSAFT_hs_D, gii_PCSAFT_hs_TD, ahs_PCSAFT_TD
    if (GETDERAHC(6) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + gl%giix1_PCSAFT(i,6,xi)*(gl%mPCSAFT(i) - 1.d0 &
		     & )*gl%molfractions(i)/gl%gii_PCSAFT(i,1) - gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,2,xi)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1) &
			& **2 - gl%gii_PCSAFT(i,6)*gl%giix1_PCSAFT(i,1,xi)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2 - &
			& gl%giix1_PCSAFT(i,4,xi)*gl%gii_PCSAFT(i,2)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2 + 2.d0*gl%gii_PCSAFT(i,4) &
			& *gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,1,xi)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3
            end do
            gl%ahcx1_PCSAFT(6,xi) = gl%ahs_PCSAFT(6)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(6,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,6)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) &
			& + gl%gii_PCSAFT(xi,4)*gl%gii_PCSAFT(xi,2)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - sum
        end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires ahs_PCSAFT(7) and gii_PCSAFT_hs(7,6,5,4,2,1)  (derivative numbers)
    if (GETDERAHC(7) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)**3 &
			     & *gl%giix1_PCSAFT(i,7,xi) - gl%gii_PCSAFT(i,1)**2*(2.d0*gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,6,xi) + &
			     & gl%gii_PCSAFT(i,5)*gl%giix1_PCSAFT(i,2,xi) + gl%gii_PCSAFT(i,7)*gl%giix1_PCSAFT(i,1,xi) + 2.d0* &
			     & gl%gii_PCSAFT(i,6)*gl%giix1_PCSAFT(i,4,xi)) + gl%gii_PCSAFT(i,1)*(2.d0*gl%gii_PCSAFT(i,4)**2* &
			     & gl%giix1_PCSAFT(i,2,xi) + 4.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,6)*gl%giix1_PCSAFT(i,1,xi) + 2.d0* &
			     & gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,4,xi)*gl%gii_PCSAFT(i,2) + gl%gii_PCSAFT(i,5)*gl%gii_PCSAFT(i,2)* &
			     & gl%giix1_PCSAFT(i,1,xi)) - 4.d0*gl%gii_PCSAFT(i,4)**2*gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,1,xi) - &
			     & gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,1)**2*gl%giix1_PCSAFT(i,5,xi) - gl%gii_PCSAFT(i,1)*(2.d0*gl%gii_PCSAFT(i,4)* &
			     & gl%giix1_PCSAFT(i,4,xi) + gl%gii_PCSAFT(i,5)*gl%giix1_PCSAFT(i,1,xi)) + 2.d0*gl%gii_PCSAFT(i,4)**2* &
			     & gl%giix1_PCSAFT(i,1,xi)))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
            end do
            gl%ahcx1_PCSAFT(7,xi) = gl%ahs_PCSAFT(7)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(7,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,7)*( &
			     & gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + 2.d0*gl%gii_PCSAFT(xi,4)*gl%gii_PCSAFT(xi,6)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 &
			     & + gl%gii_PCSAFT(xi,5)*gl%gii_PCSAFT(xi,2)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0*gl%gii_PCSAFT(xi,4)**2* &
			     & gl%gii_PCSAFT(xi,2)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**3 - sum
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_D, gii_PCSAFT_hs_DD, gii_PCSAFT_hs_DDD, ahs_PCSAFT_DDD
    if (GETDERAHC(8) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*( &
				& -gl%gii_PCSAFT(i,1)**3*gl%giix1_PCSAFT(i,8,xi) + gl%gii_PCSAFT(i,1)**2*(3.d0*gl%gii_PCSAFT(i,2)* &
				& gl%giix1_PCSAFT(i,3,xi) + 3.d0*gl%gii_PCSAFT(i,3)*gl%giix1_PCSAFT(i,2,xi) + &
				& gl%gii_PCSAFT(i,8)*gl%giix1_PCSAFT(i,1,xi)) - 6.d0*gl%gii_PCSAFT(i,1)*gl%gii_PCSAFT(i,2)*( &
				& gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,2,xi) + gl%gii_PCSAFT(i,3)*gl%giix1_PCSAFT(i,1,xi)) + 6.d0* &
				& gl%gii_PCSAFT(i,2)**3*gl%giix1_PCSAFT(i,1,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
            end do
            gl%ahcx1_PCSAFT(8,xi) = gl%ahs_PCSAFT(8)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(8,xi)*gl%mmean_PCSAFT - &
			& gl%gii_PCSAFT(xi,8)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + 3.d0*gl%gii_PCSAFT(xi,2)* &
			& gl%gii_PCSAFT(xi,3)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0*gl%gii_PCSAFT(xi,2)**3*(gl%mPCSAFT(xi) - 1.d0)/ &
			& gl%gii_PCSAFT(xi,1)**3 + sum
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires ahs_PCSAFT(9) and gii_PCSAFT_hs(9,5,4,1)  (derivative numbers)
    if (GETDERAHC(9) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(-gl%gii_PCSAFT(i,1)**3* &
			& gl%giix1_PCSAFT(i,9,xi) + gl%gii_PCSAFT(i,1)**2*(3.d0*gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,5,xi) + 3.d0* &
			& gl%gii_PCSAFT(i,5)*gl%giix1_PCSAFT(i,4,xi) + gl%gii_PCSAFT(i,9)*gl%giix1_PCSAFT(i,1,xi)) - 6.d0*gl%gii_PCSAFT(i,1)* &
			& gl%gii_PCSAFT(i,4)*(gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,4,xi) + gl%gii_PCSAFT(i,5)*gl%giix1_PCSAFT(i,1,xi)) + 6.d0* &
			& gl%gii_PCSAFT(i,4)**3*gl%giix1_PCSAFT(i,1,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
            end do
            gl%ahcx1_PCSAFT(9,xi) = gl%ahs_PCSAFT(9)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(9,xi)*gl%mmean_PCSAFT - &
			& gl%gii_PCSAFT(xi,9)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + 3.d0*gl%gii_PCSAFT(xi,4)*gl%gii_PCSAFT(xi,5)*(gl%mPCSAFT(xi) - 1.d0)/ &
			& gl%gii_PCSAFT(xi,1)**2 - 2.d0*gl%gii_PCSAFT(xi,4)**3*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**3 + sum
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires ahs_PCSAFT(10) and gii_PCSAFT_hs(10,6,4,3,2,1)  (derivative numbers)
    if (GETDERAHC(10) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 1, gl%ncomp
                sum = sum + (gl%mPCSAFT(i) - 1.d0)*( &
			& -gl%gii_PCSAFT(i,1)**3*gl%giix1_PCSAFT(i,10,xi) + gl%gii_PCSAFT(i,1)**2*(gl%gii_PCSAFT(i,4)* &
			& gl%giix1_PCSAFT(i,3,xi) + 2.d0*gl%gii_PCSAFT(i,6)*gl%giix1_PCSAFT(i,2,xi) + &
			& gl%gii_PCSAFT(i,10)*gl%giix1_PCSAFT(i,1,xi) + 2.d0*gl%giix1_PCSAFT(i,6,xi)*gl%gii_PCSAFT(i,2) + &
			& gl%giix1_PCSAFT(i,4,xi)*gl%gii_PCSAFT(i,3)) - 2.d0*gl%gii_PCSAFT(i,1)*(2.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2) &
			& *gl%giix1_PCSAFT(i,2,xi) + gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,3)*gl%giix1_PCSAFT(i,1,xi) + 2.d0* &
			& gl%gii_PCSAFT(i,6)*gl%gii_PCSAFT(i,2)*gl%giix1_PCSAFT(i,1,xi) + gl%giix1_PCSAFT(i,4,xi)*gl%gii_PCSAFT(i,2)** &
			& 2) + 6.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2)**2*gl%giix1_PCSAFT(i,1,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
            end do
            gl%ahcx1_PCSAFT(10,xi) = gl%ahs_PCSAFT(10)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(10,xi)*gl%mmean_PCSAFT - gl%gii_PCSAFT(xi,10)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + gl%gii_PCSAFT(xi,4)* &
			& gl%gii_PCSAFT(xi,3)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 + 2.d0*gl%gii_PCSAFT(xi,6)*gl%gii_PCSAFT(xi,2)* &
			& (gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0*gl%gii_PCSAFT(xi,4)*gl%gii_PCSAFT(xi,2)**2*(gl%mPCSAFT(xi) - 1.d0)/ &
			& gl%gii_PCSAFT(xi,1)**3 + sum
        end do
    end if
    
    !DEC$ END IF
end subroutine AHCX1DERIVS
    

    
subroutine AHCX2DERIVS(gl,GETDERAHC)

! Henning Markgraf, June 2016

    ! a_hc: Hard chain reference contribution to the Helmholtz free energy
    ! defined by eq. A.4 in Gross, Sadowski 2001:
    ! a_hc = mmean*a_hs - Sum( x_i*(m_i-1)*ln(g_ii_hs) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAHC
    !output: ahcx2_PCSAFT
    !working variables
    double precision :: sum
    integer :: i, xi, xj
    integer:: errorfld
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Initializations
    sum = 0.d0
    
	! III. X2 Derivatives
    !calculate the derivatives of a_hc
    ! 1: a_hc
    if (GETDERAHC(1) .eq. 1) then
	    do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum + (gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)) &
			            & *(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
		            end do
		            gl%ahcx2_PCSAFT(1,xi,xj) = gl%ahsx1_PCSAFT(1,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(1,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) - &
                 & gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) - sum
                end if
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires ...
    if (GETDERAHC(2) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum -(gl%gii_PCSAFT(i,1)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,2,xi,xj) + gl%gii_PCSAFT(i,2)* &
                 & gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,2,xi) - gl%giix1_PCSAFT(i,2,xj)* &
                 & gl%giix1_PCSAFT(i,1,xi)) - 2.d0*gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3
		            end do
		            gl%ahcx2_PCSAFT(2,xi,xj) = gl%ahsx1_PCSAFT(2,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(2,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(2,xi,xj)*gl%mmean_PCSAFT - gl%giix1_PCSAFT(xi,2,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + gl%gii_PCSAFT(xi,2)* &
                 & gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - gl%giix1_PCSAFT(xj,2,xi)*(gl%mPCSAFT(xj) - 1.d0)/ &
                 & gl%gii_PCSAFT(xj,1) + gl%gii_PCSAFT(xj,2)*gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 + sum
                end if
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_D, gii_PCSAFT_hs_DD, a_hs_DD
    if (GETDERAHC(3) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,3,xi,xj) + 2.d0* &
                 & gl%gii_PCSAFT(i,2)*gl%giix2_PCSAFT(i,2,xi,xj) + gl%gii_PCSAFT(i,3)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,3,xi) - 2.d0*gl%giix1_PCSAFT(i,2,xj)* &
                 & gl%giix1_PCSAFT(i,2,xi) - gl%giix1_PCSAFT(i,3,xj)*gl%giix1_PCSAFT(i,1,xi) - 4.d0*gl%gii_PCSAFT(i,2)* &
                 & (gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,2,xi,xj) + gl%gii_PCSAFT(i,2)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,2,xi) - gl%giix1_PCSAFT(i,2,xj)*gl%giix1_PCSAFT(i,1,xi))/ &
                 & gl%gii_PCSAFT(i,1) - 2.d0*gl%gii_PCSAFT(i,3)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix1_PCSAFT(i,1,xi))/gl%gii_PCSAFT(i,1) + 6.d0*gl%gii_PCSAFT(i,2)**2*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))/gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
		            end do
		            gl%ahcx2_PCSAFT(3,xi,xj) = gl%ahsx1_PCSAFT(3,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(3,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(3,xi,xj)* &
                 & gl%mmean_PCSAFT - gl%giix1_PCSAFT(xi,3,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + 2.d0*gl%gii_PCSAFT(xi,2)* &
                 & gl%giix1_PCSAFT(xi,2,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 + gl%gii_PCSAFT(xi,3)*gl%giix1_PCSAFT(xi,1,xj) &
                 & *(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0*gl%gii_PCSAFT(xi,2)**2*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/ &
                 & gl%gii_PCSAFT(xi,1)**3 - gl%giix1_PCSAFT(xj,3,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) + 2.d0*gl%gii_PCSAFT(xj,2) &
                 & *gl%giix1_PCSAFT(xj,2,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 + gl%gii_PCSAFT(xj,3)* &
                 & gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 - 2.d0*gl%gii_PCSAFT(xj,2)**2*gl%giix1_PCSAFT(xj,1,xi)*( &
                 & gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**3 - sum
                end if
            end do
        end do
    end if

    ! 4: 1ST DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires ...
    if (GETDERAHC(4) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum + -(gl%gii_PCSAFT(i,1)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,4,xi,xj) + gl%gii_PCSAFT(i,4)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,4,xi) - gl%giix1_PCSAFT(i,4,xj)*gl%giix1_PCSAFT(i,1,xi)) - 2.d0* &
                 & gl%gii_PCSAFT(i,4)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3
		            end do
		            gl%ahcx2_PCSAFT(4,xi,xj) = gl%ahsx1_PCSAFT(4,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(4,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%giix1_PCSAFT(xi,4,xj)*(gl%mPCSAFT(xi) - 1.d0 &
                 & )/gl%gii_PCSAFT(xi,1) + gl%gii_PCSAFT(xi,4)*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - &
                 & gl%giix1_PCSAFT(xj,4,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) + gl%gii_PCSAFT(xj,4)*gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0) &
                 & /gl%gii_PCSAFT(xj,1)**2 + sum
                end if
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_T, gii_PCSAFT_hs_TT, a_hs_TT
    if (GETDERAHC(5) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,5,xi,xj) + 2.d0*gl%gii_PCSAFT(i,4)*gl%giix2_PCSAFT(i,4,xi,xj) &
                 & + gl%gii_PCSAFT(i,5)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,5,xi) - 2.d0* &
                 & gl%giix1_PCSAFT(i,4,xj)*gl%giix1_PCSAFT(i,4,xi) - gl%giix1_PCSAFT(i,5,xj)*gl%giix1_PCSAFT(i,1,xi) - 4.d0* &
                 & gl%gii_PCSAFT(i,4)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,4,xi,xj) + gl%gii_PCSAFT(i,4)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,4,xi) - gl%giix1_PCSAFT(i,4,xj)*gl%giix1_PCSAFT(i,1,xi))/gl%gii_PCSAFT(i,1) - &
                 & 2.d0*gl%gii_PCSAFT(i,5)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))/ &
                 & gl%gii_PCSAFT(i,1) + 6.d0*gl%gii_PCSAFT(i,4)**2*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix1_PCSAFT(i,1,xi))/gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
		            end do
		            gl%ahcx2_PCSAFT(5,xi,xj) = gl%ahsx1_PCSAFT(5,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(5,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(5,xi,xj)*gl%mmean_PCSAFT - &
                 & gl%giix1_PCSAFT(xi,5,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + 2.d0*gl%gii_PCSAFT(xi,4)*gl%giix1_PCSAFT(xi,4,xj)*(gl%mPCSAFT(xi) &
                 & - 1.d0)/gl%gii_PCSAFT(xi,1)**2 + gl%gii_PCSAFT(xi,5)*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0 &
                 & *gl%gii_PCSAFT(xi,4)**2*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**3 - gl%giix1_PCSAFT(xj,5,xi)*( &
                 & gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) + 2.d0*gl%gii_PCSAFT(xj,4)*gl%giix1_PCSAFT(xj,4,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 &
                 & + gl%gii_PCSAFT(xj,5)*gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 - 2.d0*gl%gii_PCSAFT(xj,4)**2* &
                 & gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**3 - sum
                end if
            end do
        end do
    end if

    ! 6: 2ND MIXED DERIVATIVE OF a_hc WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires gii_PCSAFT_hs, gii_PCSAFT_hs_T, gii_PCSAFT_hs_D, gii_PCSAFT_hs_TD, a_hs_TD
    if (GETDERAHC(6) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 1, gl%ncomp
		                sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,6,xi,xj) + gl%gii_PCSAFT(i,4)* &
                 & gl%giix2_PCSAFT(i,2,xi,xj) + gl%gii_PCSAFT(i,6)*gl%giix2_PCSAFT(i,1,xi,xj) + gl%gii_PCSAFT(i,2)* &
                 & gl%giix2_PCSAFT(i,4,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,6,xi) - gl%giix1_PCSAFT(i,4,xj)* &
                 & gl%giix1_PCSAFT(i,2,xi) - gl%giix1_PCSAFT(i,6,xj)*gl%giix1_PCSAFT(i,1,xi) - gl%giix1_PCSAFT(i,2,xj)* &
                 & gl%giix1_PCSAFT(i,4,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2 - 2.d0*gl%gii_PCSAFT(i,4)*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)* &
                 & gl%giix2_PCSAFT(i,2,xi,xj) + gl%gii_PCSAFT(i,2)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix1_PCSAFT(i,2,xi) - gl%giix1_PCSAFT(i,2,xj)*gl%giix1_PCSAFT(i,1,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3 - &
                 & 2.d0*gl%gii_PCSAFT(i,6)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)) &
                 & *(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3 - 2.d0*gl%gii_PCSAFT(i,2)*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)* &
                 & gl%giix2_PCSAFT(i,4,xi,xj) + gl%gii_PCSAFT(i,4)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix1_PCSAFT(i,4,xi) - gl%giix1_PCSAFT(i,4,xj)*gl%giix1_PCSAFT(i,1,xi))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3 + 6.d0* &
                 & gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix1_PCSAFT(i,1,xi))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
		            end do
		            gl%ahcx2_PCSAFT(6,xi,xj) = gl%ahsx1_PCSAFT(6,xj)*gl%mPCSAFT(xi) + gl%ahsx1_PCSAFT(6,xi)*gl%mPCSAFT(xj) + gl%ahsx2_PCSAFT(6,xi,xj)*gl%mmean_PCSAFT - &
                 & gl%giix1_PCSAFT(xi,6,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) + gl%gii_PCSAFT(xi,4)*gl%giix1_PCSAFT(xi,2,xj)*( &
                 & gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 + gl%gii_PCSAFT(xi,6)*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1) &
                 & **2 + gl%gii_PCSAFT(xi,2)*gl%giix1_PCSAFT(xi,4,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**2 - 2.d0* &
                 & gl%gii_PCSAFT(xi,4)*gl%gii_PCSAFT(xi,2)*gl%giix1_PCSAFT(xi,1,xj)*(gl%mPCSAFT(xi) - 1.d0)/gl%gii_PCSAFT(xi,1)**3 - &
                 & gl%giix1_PCSAFT(xj,6,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) + gl%gii_PCSAFT(xj,4)*gl%giix1_PCSAFT(xj,2,xi)*( &
                 & gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 + gl%gii_PCSAFT(xj,6)*gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1) &
                 & **2 + gl%gii_PCSAFT(xj,2)*gl%giix1_PCSAFT(xj,4,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**2 - 2.d0* &
                 & gl%gii_PCSAFT(xj,4)*gl%gii_PCSAFT(xj,2)*gl%giix1_PCSAFT(xj,1,xi)*(gl%mPCSAFT(xj) - 1.d0)/gl%gii_PCSAFT(xj,1)**3 - sum
                end if
            end do
        end do
    end if

    !DEC$ END IF
end subroutine AHCX2DERIVS
    

    
subroutine AHCX3DERIVS(gl,GETDERAHC)

! Henning Markgraf, June 2016

    ! a_hc: Hard chain reference contribution to the Helmholtz free energy
    ! defined by eq. A.4 in Gross, Sadowski 2001:
    ! a_hc = mmean*a_hs - Sum( x_i*(m_i-1)*ln(g_ii_hs) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAHC
    !output: ahcx3_PCSAFT
    !working variables
    double precision :: sum
    integer :: i, xi, xj, xk
    integer:: errorfld
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X3 DERIVATIVES of a_hc
    ! 1: d^3(a_hc)/(dxi dxj dxk)
    if (GETDERAHC(1) .eq. 1) then
	    do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 1, gl%ncomp
		                    sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)*gl%giix3_PCSAFT(i,1,xi,xj,xk) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix2_PCSAFT(i,1,xi,xk) - gl%giix2_PCSAFT(i,1,xj,xk)*gl%giix1_PCSAFT(i,1,xi) + &
                 & gl%giix2_PCSAFT(i,1,xi,xj)*gl%giix1_PCSAFT(i,1,xk) - 2.d0*gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)* &
                 & gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))/gl%gii_PCSAFT(i,1))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)** &
                 & 2
		                end do
		                gl%ahcx3_PCSAFT(1,xi,xj,xk) = gl%ahsx2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) + gl%ahsx2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) + gl%ahsx3_PCSAFT(1,xi,xj,xk)*gl%mmean_PCSAFT + &
                 & gl%ahsx2_PCSAFT(1,xi,xk)*gl%mPCSAFT(xj) - (gl%mPCSAFT(xk) - 1.d0)*(gl%gii_PCSAFT(xk,1)*gl%giix2_PCSAFT(xk,1,xi,xj) - &
                 & gl%giix1_PCSAFT(xk,1,xj)*gl%giix1_PCSAFT(xk,1,xi))/gl%gii_PCSAFT(xk,1)**2 - (gl%mPCSAFT(xi) - 1.d0)*(gl%gii_PCSAFT(xi,1)* &
                 & gl%giix2_PCSAFT(xi,1,xj,xk) - gl%giix1_PCSAFT(xi,1,xj)*gl%giix1_PCSAFT(xi,1,xk))/gl%gii_PCSAFT(xi,1)**2 - (gl%mPCSAFT(xj) - 1.d0)* &
                 & (gl%gii_PCSAFT(xj,1)*gl%giix2_PCSAFT(xj,1,xi,xk) - gl%giix1_PCSAFT(xj,1,xi)*gl%giix1_PCSAFT(xj,1,xk))/gl%gii_PCSAFT(xj,1)**2 - sum
                    end if
		        end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires ...
    if (GETDERAHC(2) .eq. 1) then
	    do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 1, gl%ncomp
		                    sum = sum -(gl%mPCSAFT(i) - 1.d0)*(-gl%gii_PCSAFT(i,1)**2*( &
                 & -gl%gii_PCSAFT(i,1)*gl%giix3_PCSAFT(i,2,xi,xj,xk) - gl%gii_PCSAFT(i,2)*gl%giix3_PCSAFT(i,1,xi,xj,xk) + &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix2_PCSAFT(i,2,xi,xk) + gl%giix1_PCSAFT(i,2,xj)*gl%giix2_PCSAFT(i,1,xi,xk) + &
                 & gl%giix2_PCSAFT(i,1,xj,xk)*gl%giix1_PCSAFT(i,2,xi) + gl%giix2_PCSAFT(i,2,xj,xk)*gl%giix1_PCSAFT(i,1,xi) - &
                 & gl%giix2_PCSAFT(i,1,xi,xj)*gl%giix1_PCSAFT(i,2,xk) - gl%giix2_PCSAFT(i,2,xi,xj)*gl%giix1_PCSAFT(i,1,xk)) &
                 & - 2.d0*gl%gii_PCSAFT(i,1)*(gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,2,xi,xj) + &
                 & gl%gii_PCSAFT(i,2)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,2,xi) - &
                 & gl%giix1_PCSAFT(i,2,xj)*gl%giix1_PCSAFT(i,1,xi)) + gl%giix1_PCSAFT(i,2,xk)*(gl%gii_PCSAFT(i,1)* &
                 & gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))) + 2.d0*gl%gii_PCSAFT(i,2)* &
                 & gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)) + &
                 & 2.d0*gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,1)*(-gl%gii_PCSAFT(i,1)*gl%giix3_PCSAFT(i,1,xi,xj,xk) + gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix2_PCSAFT(i,1,xi,xk) + gl%giix2_PCSAFT(i,1,xj,xk)*gl%giix1_PCSAFT(i,1,xi) - gl%giix2_PCSAFT(i,1,xi,xj)* &
                 & gl%giix1_PCSAFT(i,1,xk)) + 2.d0*gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))))*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**4
		                end do
		                gl%ahcx3_PCSAFT(2,xi,xj,xk) = gl%ahsx2_PCSAFT(2,xj,xk)*gl%mPCSAFT(xi) + gl%ahsx2_PCSAFT(2,xi,xj)*gl%mPCSAFT(xk) + gl%ahsx3_PCSAFT(2,xi,xj,xk)* &
                 & gl%mmean_PCSAFT + gl%ahsx2_PCSAFT(2,xi,xk)*gl%mPCSAFT(xj) - ( &
                 & gl%mPCSAFT(xk) - 1.d0)*(gl%gii_PCSAFT(xk,1)*gl%giix2_PCSAFT(xk,2,xi,xj) + gl%gii_PCSAFT(xk,2)*gl%giix2_PCSAFT(xk,1,xi,xj) &
                 & - gl%giix1_PCSAFT(xk,1,xj)*gl%giix1_PCSAFT(xk,2,xi) - gl%giix1_PCSAFT(xk,2,xj)*gl%giix1_PCSAFT(xk,1,xi))/ &
                 & gl%gii_PCSAFT(xk,1)**2 + 2.d0*gl%gii_PCSAFT(xk,2)*(gl%mPCSAFT(xk) - 1.d0)*(gl%gii_PCSAFT(xk,1)*gl%giix2_PCSAFT(xk,1,xi,xj) - &
                 & gl%giix1_PCSAFT(xk,1,xj)*gl%giix1_PCSAFT(xk,1,xi))/gl%gii_PCSAFT(xk,1)**3 - (gl%mPCSAFT(xi) - 1.d0)*(gl%gii_PCSAFT(xi,1)* &
                 & gl%giix2_PCSAFT(xi,2,xj,xk) + gl%gii_PCSAFT(xi,2)*gl%giix2_PCSAFT(xi,1,xj,xk) - gl%giix1_PCSAFT(xi,1,xj)* &
                 & gl%giix1_PCSAFT(xi,2,xk) - gl%giix1_PCSAFT(xi,2,xj)*gl%giix1_PCSAFT(xi,1,xk))/gl%gii_PCSAFT(xi,1)**2 + 2.d0* &
                 & gl%gii_PCSAFT(xi,2)*(gl%mPCSAFT(xi) - 1.d0)*(gl%gii_PCSAFT(xi,1)*gl%giix2_PCSAFT(xi,1,xj,xk) - gl%giix1_PCSAFT(xi,1,xj)* &
                 & gl%giix1_PCSAFT(xi,1,xk))/gl%gii_PCSAFT(xi,1)**3 - (gl%mPCSAFT(xj) - 1.d0)*(gl%gii_PCSAFT(xj,1)*gl%giix2_PCSAFT(xj,2,xi,xk) + &
                 & gl%gii_PCSAFT(xj,2)*gl%giix2_PCSAFT(xj,1,xi,xk) - gl%giix1_PCSAFT(xj,1,xi)*gl%giix1_PCSAFT(xj,2,xk) - &
                 & gl%giix1_PCSAFT(xj,2,xi)*gl%giix1_PCSAFT(xj,1,xk))/gl%gii_PCSAFT(xj,1)**2 + 2.d0*gl%gii_PCSAFT(xj,2)*(gl%mPCSAFT(xj) - 1.d0) &
                 & *(gl%gii_PCSAFT(xj,1)*gl%giix2_PCSAFT(xj,1,xi,xk) - gl%giix1_PCSAFT(xj,1,xi)*gl%giix1_PCSAFT(xj,1,xk))/gl%gii_PCSAFT(xj,1)**3 + sum
                    end if
		        end do
            end do
        end do
    end if

    ! 4: 1ST DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires ...
    if (GETDERAHC(4) .eq. 1) then
	    do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 1, gl%ncomp
		                    sum = sum -(gl%mPCSAFT(i) - 1.d0)*(-gl%gii_PCSAFT(i,1)**2*(-gl%gii_PCSAFT(i,1)* &
                 & gl%giix3_PCSAFT(i,3,xi,xj,xk) - gl%gii_PCSAFT(i,4)*gl%giix3_PCSAFT(i,1,xi,xj,xk) + gl%giix1_PCSAFT(i,1,xj)* &
                 & gl%giix2_PCSAFT(i,4,xi,xk) + gl%giix1_PCSAFT(i,4,xj)*gl%giix2_PCSAFT(i,1,xi,xk) + gl%giix2_PCSAFT(i,1,xj,xk)* &
                 & gl%giix1_PCSAFT(i,4,xi) + gl%giix2_PCSAFT(i,4,xj,xk)*gl%giix1_PCSAFT(i,1,xi) - gl%giix2_PCSAFT(i,1,xi,xj)* &
                 & gl%giix1_PCSAFT(i,4,xk) - gl%giix2_PCSAFT(i,4,xi,xj)*gl%giix1_PCSAFT(i,1,xk)) - 2.d0*gl%gii_PCSAFT(i,1)*( &
                 & gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,4,xi,xj) + gl%gii_PCSAFT(i,4)*gl%giix2_PCSAFT(i,1,xi,xj) &
                 & - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,4,xi) - gl%giix1_PCSAFT(i,4,xj)*gl%giix1_PCSAFT(i,1,xi)) + &
                 & gl%giix1_PCSAFT(i,4,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi) &
                 & )) + 2.d0*gl%gii_PCSAFT(i,4)*gl%giix1_PCSAFT(i,1,xk)*(gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - &
                 & gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi)) + 2.d0*gl%gii_PCSAFT(i,4)*(gl%gii_PCSAFT(i,1)*(-gl%gii_PCSAFT(i,1)* &
                 & gl%giix3_PCSAFT(i,1,xi,xj,xk) + gl%giix1_PCSAFT(i,1,xj)*gl%giix2_PCSAFT(i,1,xi,xk) + gl%giix2_PCSAFT(i,1,xj,xk)* &
                 & gl%giix1_PCSAFT(i,1,xi) - gl%giix2_PCSAFT(i,1,xi,xj)*gl%giix1_PCSAFT(i,1,xk)) + 2.d0*gl%giix1_PCSAFT(i,1,xk)*( &
                 & gl%gii_PCSAFT(i,1)*gl%giix2_PCSAFT(i,1,xi,xj) - gl%giix1_PCSAFT(i,1,xj)*gl%giix1_PCSAFT(i,1,xi))))*gl%molfractions(i)/gl%gii_PCSAFT(i,1) &
                 & **4
		                end do
		                gl%ahcx3_PCSAFT(3,xi,xj,xk) = gl%ahsx2_PCSAFT(4,xj,xk)*gl%mPCSAFT(xi) + gl%ahsx2_PCSAFT(4,xi,xj)*gl%mPCSAFT(xk) + gl%ahsx3_PCSAFT(3,xi,xj,xk)*gl%mmean_PCSAFT + &
                 & gl%ahsx2_PCSAFT(4,xi,xk)*gl%mPCSAFT(xj) - (gl%mPCSAFT(xk) - 1.d0)*(gl%gii_PCSAFT(xk,1)*gl%giix2_PCSAFT(xk,4,xi,xj) + &
                 & gl%gii_PCSAFT(xk,4)*gl%giix2_PCSAFT(xk,1,xi,xj) - gl%giix1_PCSAFT(xk,1,xj)*gl%giix1_PCSAFT(xk,4,xi) - &
                 & gl%giix1_PCSAFT(xk,4,xj)*gl%giix1_PCSAFT(xk,1,xi))/gl%gii_PCSAFT(xk,1)**2 + 2.d0*gl%gii_PCSAFT(xk,4)*(gl%mPCSAFT(xk) - 1.d0)*( &
                 & gl%gii_PCSAFT(xk,1)*gl%giix2_PCSAFT(xk,1,xi,xj) - gl%giix1_PCSAFT(xk,1,xj)*gl%giix1_PCSAFT(xk,1,xi))/gl%gii_PCSAFT(xk,1)**3 - ( &
                 & gl%mPCSAFT(xi) - 1.d0)*(gl%gii_PCSAFT(xi,1)*gl%giix2_PCSAFT(xi,4,xj,xk) + gl%gii_PCSAFT(xi,4)*gl%giix2_PCSAFT(xi,1,xj,xk) - &
                 & gl%giix1_PCSAFT(xi,1,xj)*gl%giix1_PCSAFT(xi,4,xk) - gl%giix1_PCSAFT(xi,4,xj)*gl%giix1_PCSAFT(xi,1,xk))/gl%gii_PCSAFT(xi,1)**2 &
                 & + 2.d0*gl%gii_PCSAFT(xi,4)*(gl%mPCSAFT(xi) - 1.d0)*(gl%gii_PCSAFT(xi,1)*gl%giix2_PCSAFT(xi,1,xj,xk) - gl%giix1_PCSAFT(xi,1,xj)* &
                 & gl%giix1_PCSAFT(xi,1,xk))/gl%gii_PCSAFT(xi,1)**3 - (gl%mPCSAFT(xj) - 1.d0)*(gl%gii_PCSAFT(xj,1)*gl%giix2_PCSAFT(xj,4,xi,xk) + &
                 & gl%gii_PCSAFT(xj,4)*gl%giix2_PCSAFT(xj,1,xi,xk) - gl%giix1_PCSAFT(xj,1,xi)*gl%giix1_PCSAFT(xj,4,xk) - &
                 & gl%giix1_PCSAFT(xj,4,xi)*gl%giix1_PCSAFT(xj,1,xk))/gl%gii_PCSAFT(xj,1)**2 + 2.d0*gl%gii_PCSAFT(xj,4)*(gl%mPCSAFT(xj) - 1.d0)*( &
                 & gl%gii_PCSAFT(xj,1)*gl%giix2_PCSAFT(xj,1,xi,xk) - gl%giix1_PCSAFT(xj,1,xi)*gl%giix1_PCSAFT(xj,1,xk))/gl%gii_PCSAFT(xj,1)**3 + sum
                    end if
		        end do
            end do
        end do
    end if

    !DEC$ END IF
end subroutine AHCX3DERIVS
    



    end module pc_saft_AHCX_derivs_module
