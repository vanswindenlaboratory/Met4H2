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

    ! module for file pc-saft_ANCX_derivs.f90
    module pc_saft_ANCX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module
    use pc_saft_ancillary_routines_module
    use variables_transformation_module


    contains


    
subroutine ZETAX1DERIVS(gl,D, GETDERZETA, n_zeta)

! Henning Markgraf, June 2016 

    ! zeta: abbreviation defined by eq. 9 in Gross, Sadowski 2001:
    ! zeta_n = pi/6*rho*Sum(x_i*m_i*d_i**n) , n element of {0, 1, 2, 3}
    ! dependent on x_i, D and T

!--------------------------------------------------------------------------------------------------
! All first composition derivatives of zeta_n are calculated in this subroutine
! (all second and third composition derivatives equal zero)
!--------------------------------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, intent (in) :: n_zeta
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERZETA
    !output: z0x1_PCSAFT, z1x1_PCSAFT, z2x1_PCSAFT, z3x1_PCSAFT (module variables)
    !working variables
    double precision :: derivative, sum
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Initializations
    sum = 0.d0
    
    ! III. X1 DERIVATIVES
    ! calculate the component derivatives of zeta_0, zeta_1, zeta_2 and zeta_3
    ! 1: zeta_n
    ! 2: 1ST DERIVATIVE OF zeta_n WITH RESPECT TO x_i and D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERZETA(1) .eq. 1 .OR. GETDERZETA(2) .eq. 1) then
        do i = 1, gl%ncomp
            derivative = piPCSAFT/6.d0*D*gl%mPCSAFT(i)*gl%di_PCSAFT(i,1)**n_zeta
            select case (n_zeta)
                case (0)
                    gl%z0x1_PCSAFT(1,i) = derivative
                    gl%z0x1_PCSAFT(2,i) = derivative
                case (1)
                    gl%z1x1_PCSAFT(1,i) = derivative
                    gl%z1x1_PCSAFT(2,i) = derivative
                case (2)
                    gl%z2x1_PCSAFT(1,i) = derivative
                    gl%z2x1_PCSAFT(2,i) = derivative
                case (3)
                    gl%z3x1_PCSAFT(1,i) = derivative
                    gl%z3x1_PCSAFT(2,i) = derivative
            end select
        end do
    end if    
    
    ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! DERZETA(3) = 0.d0 (do nothing, already initialized as 0.d0)
    
    ! possible to use derivative number one as the basis for derivative 4 (would it save significant calculation time)
    ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERZETA(4) .eq. 1 .OR. GETDERZETA(6) .EQ. 1) then
        do i = 1, gl%ncomp
            derivative = piPCSAFT/6.d0*D*gl%mPCSAFT(i)*n_zeta*gl%di_PCSAFT(i,1)**(n_zeta-1)*gl%di_PCSAFT(i,4)
            select case (n_zeta)
                case (0)
                    gl%z0x1_PCSAFT(4,i) = derivative
                    gl%z0x1_PCSAFT(6,i) = derivative
                case (1)
                    gl%z1x1_PCSAFT(4,i) = derivative
                    gl%z1x1_PCSAFT(6,i) = derivative
                case (2)
                    gl%z2x1_PCSAFT(4,i) = derivative
                    gl%z2x1_PCSAFT(6,i) = derivative
                case (3)
                    gl%z3x1_PCSAFT(4,i) = derivative
                    gl%z3x1_PCSAFT(6,i) = derivative
            end select
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERZETA(5) .eq. 1 .OR. GETDERZETA(7) .EQ. 1) then
        do i = 1, gl%ncomp
            derivative = (1.d0/6.d0)*piPCSAFT*gl%di_PCSAFT(i,1)**(n_zeta-2)*n_zeta*D* &
			& (gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5) + gl%di_PCSAFT(i,4)**2*n_zeta - gl%di_PCSAFT(i,4)**2)*gl%mPCSAFT(i)
            select case (n_zeta)
                case (0)
                    gl%z0x1_PCSAFT(5,i) = derivative
                    gl%z0x1_PCSAFT(7,i) = derivative
                case (1)
                    gl%z1x1_PCSAFT(5,i) = derivative
                    gl%z1x1_PCSAFT(7,i) = derivative
                case (2)
                    gl%z2x1_PCSAFT(5,i) = derivative
                    gl%z2x1_PCSAFT(7,i) = derivative
                case (3)
                    gl%z3x1_PCSAFT(5,i) = derivative
                    gl%z3x1_PCSAFT(7,i) = derivative
            end select
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! DERZETA(6) = DERZETA(4)
    
    ! 7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO rho, T, AND T, MULTIPLIED BY T*T*rho
    ! DERZETA(7) = DERZETA(5)
    
    ! 8: 3RD DERIVATIVE OF F WITH RESPECT TO rho, MULTIPLIED BY rho^3
    ! DERZETA(8) = 0.d0 (do nothing, already initialized as 0.d0)
    
    ! 9: 3RD DERIVATIVE OF F WITH RESPECT TO T, MULTIPLIED BY T^3
    if (GETDERZETA(9) .EQ. 1) then
        do i = 1, gl%ncomp
            derivative = (1.d0/6.d0)*piPCSAFT*gl%di_PCSAFT(i,1)**(n_zeta - 3)*n_zeta*D*(gl%di_PCSAFT(i,1)**2*gl%di_PCSAFT(i,9) + 3.d0* &
			     & gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*(n_zeta - 1) + gl%di_PCSAFT(i,4)**3*(n_zeta**2 - 3*n_zeta + &
			     & 2))*gl%mPCSAFT(i)
            select case (n_zeta)
                case (0)
                    gl%z0x1_PCSAFT(9,i) = derivative
                case (1)
                    gl%z1x1_PCSAFT(9,i) = derivative
                case (2)
                    gl%z2x1_PCSAFT(9,i) = derivative
                case (3)
                    gl%z3x1_PCSAFT(9,i) = derivative
            end select
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO T, rho, AND rho, MULTIPLIED BY T*rho*rho
    ! DERZETA(10) = 0.d0 (do nothing, already initialized as 0.d0)
    
!DEC$ END IF
end subroutine ZETAX1DERIVS

    
subroutine AASSOCX1(gl,TEMP, DENS, GETDERAASSOC1)

    ! a_assoc: association contribution to the Helmholtz free energy
    ! defined by eq. 21 in Chapman et al. 1990:
    ! dependent on T and D






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (inout) :: TEMP, DENS
    integer, dimension (nx1derivs), intent (inout) :: GETDERAASSOC1
    integer, dimension (nderivs) :: GETDERAASSOC
    !output: AASSOC_PCSAFT (module variable)
    integer :: xi,A,i,j
    double precision :: sum,DELTA,DP,DM,TP,TM,DP2,DM2,TP2,TM2,DP3,DM3,TP3,TM3
    double precision, dimension(gl%ncomp) :: AASSOC_P,AASSOC_M, AASSOC_0, AASSOC_P2,AASSOC_M2,AASSOC_PP,AASSOC_MM,AASSOC_PM,AASSOC_MP, AASSOC_P3,AASSOC_M3
    double precision, dimension(gl%ncomp) :: AASSOC_21,AASSOC_01,AASSOC_M21,AASSOC_2M1,AASSOC_0M1,AASSOC_M2M1
    double precision, dimension(30) :: molfractions_orig
    double precision, dimension(nx1derivs,gl%ncomp) ::ASSO_P, ASSO_M,ASSO_P2, ASSO_M2
    
    integer:: errorfld
    
    molfractions_orig = 0.d0

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    GETDERAASSOC = 0
    GETDERAASSOC(1:nx1derivs) = GETDERAASSOC1
    DELTA = 5.d-4
    
    molfractions_orig = gl%molfractions
    
    do j = 1, nx1derivs
        if (GETDERAASSOC1(j) == 1) then
            do i = 1, gl%ncomp
                !increase xi 1
                gl%molfractions(i) = gl%molfractions(i) * (1.d0 + DELTA)
                call init_mixture_PCSAFT(gl)
                call AASSOC(gl,TEMP, DENS, GETDERAASSOC)
                ASSO_P(:,i) = gl%AASSOC_PCSAFT(1:nx1derivs)
                !increase xi 2
                gl%molfractions(i) = molfractions_orig(i) * (1.d0 + 2.d0*DELTA)
                call init_mixture_PCSAFT(gl)
                call AASSOC(gl,TEMP, DENS, GETDERAASSOC)
                ASSO_P2(:,i) = gl%AASSOC_PCSAFT(1:nx1derivs)
                !decrease xi 1
                gl%molfractions(i) = molfractions_orig(i) * (1.d0 - DELTA)
                call init_mixture_PCSAFT(gl)
                call AASSOC(gl,TEMP, DENS, GETDERAASSOC)
                ASSO_M(:,i) = gl%AASSOC_PCSAFT(1:nx1derivs)
                !decrease xi 2
                gl%molfractions(i) = molfractions_orig(i) * (1.d0 - 2.d0*DELTA)
                call init_mixture_PCSAFT(gl)
                call AASSOC(gl,TEMP, DENS, GETDERAASSOC)
                ASSO_M2(:,i) = gl%AASSOC_PCSAFT(1:nx1derivs)
        
                gl%AASSOCx1_PCSAFT(j,i) = (8.d0*ASSO_P(j,i) - 8.d0*ASSO_M(j,i) - ASSO_P2(j,i) + ASSO_M2(j,i)) / (12.d0 * DELTA * gl%molfractions(i))
             end do
        end if
    end do
    gl%molfractions = molfractions_orig
  
    
!DEC$ END IF
end subroutine AASSOCX1

subroutine AASSOCX2(gl,TEMP, DENS, GETDERAASSOC)

    ! a_assoc: association contribution to the Helmholtz free energy
    ! defined by eq. 21 in Chapman et al. 1990:
    ! dependent on T and D






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (inout) :: TEMP, DENS
    integer, dimension (nx2derivs), intent (inout) :: GETDERAASSOC
    integer, dimension (nx1derivs) :: GETDERAASSOC1
    !output: AASSOC_PCSAFT (module variable)
    integer :: xi,A,i,j
    double precision :: sum,DELTA,DP,DM,TP,TM,DP2,DM2,TP2,TM2,DP3,DM3,TP3,TM3
    double precision, dimension(30) :: molfractions_orig
    double precision, dimension(nx2derivs,gl%ncomp,gl%ncomp) ::ASSO_P, ASSO_M,ASSO_P2, ASSO_M2

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
   
    GETDERAASSOC1 = 0
    GETDERAASSOC1(1:nx2derivs) = GETDERAASSOC
    DELTA = 1.d-3
    
    molfractions_orig = gl%molfractions
    
    do j = 1, nx2derivs
        if (GETDERAASSOC(j) == 1) then
            do xi = 1, gl%ncomp
                do i = 1, gl%ncomp
                    !increase xi 1
                    gl%molfractions(i) = molfractions_orig(i) * (1.d0 + DELTA)
                    call init_mixture_PCSAFT(gl)
                    call AASSOCX1(gl,TEMP, DENS, GETDERAASSOC1)
                    ASSO_P(:,xi,i) = gl%AASSOCx1_PCSAFT(1:nx2derivs,xi)
                    !increase xi 2
                    gl%molfractions(i) = molfractions_orig(i) * (1.d0 + 2.d0*DELTA)
                    call init_mixture_PCSAFT(gl)
                    call AASSOCX1(gl,TEMP, DENS, GETDERAASSOC1)
                    ASSO_P2(:,xi,i) = gl%AASSOCx1_PCSAFT(1:nx2derivs,xi)
                    !decrease xi 1
                    gl%molfractions(i) = molfractions_orig(i) * (1.d0 - DELTA)
                    call init_mixture_PCSAFT(gl)
                    call AASSOCX1(gl,TEMP, DENS, GETDERAASSOC1)
                    ASSO_M(:,xi,i) = gl%AASSOCx1_PCSAFT(1:nx2derivs,xi)
                    !decrease xi 2
                    gl%molfractions(i) = molfractions_orig(i) * (1.d0 - 2.d0*DELTA)
                    call init_mixture_PCSAFT(gl)
                    call AASSOCX1(gl,TEMP, DENS, GETDERAASSOC1)
                    ASSO_M2(:,xi,i) = gl%AASSOCx1_PCSAFT(1:nx2derivs,xi)
        
                    gl%AASSOCx2_PCSAFT(j,xi,i) = (8.d0*ASSO_P(j,xi,i) - 8.d0*ASSO_M(j,xi,i) - ASSO_P2(j,xi,i) + ASSO_M2(j,xi,i)) / (12.d0 * DELTA * gl%molfractions(i))
                end do
            end do
        end if
    end do
    gl%molfractions = molfractions_orig
    
!DEC$ END IF

end subroutine AASSOCX2

subroutine ADDX1DERIVS(gl,T,DENS,GETDERADD)
    
    ! a_dd: contribution to the residual helmholz engergy from quadrupol-quadrupol interactions 
    ! defined by eq. A.8 in Gross, Sadowski 2005:
    ! a_dd = a_2dd*a_2dd/(a_2dd-a_3dd)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    integer :: k
    !output: add_PCSAFT (module variable)
    !working variables

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_dd
    if (GETDERADD(1) .eq. 1) then
        do k = 1, gl%ncomp
        gl%ADDX1_PCSAFTD(1,k) = (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)**2) / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF a_dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_2dd, a_3dd
    if (GETDERADD(2) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(2,k) = (-2.d0 * gl%A2_PCSAFTD(2) + 2.d0 * gl%A3_PCSAFTD(2)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))&
                & * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k))*  gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 + (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))*gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if
    ! I. Schuelling 05/17
    ! 3: 2ND DERIVATIVE OF a_dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires a_2dd, a_3dd
    if (GETDERADD(3) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(3,k) = (-2.0 * gl%A2_PCSAFTD(2)**2 * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 6.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))**2 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1))*gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 2.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(3,k) - gl%A3X1_PCSAFTD(3,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
     end if
     
    ! 4: 1ST DERIVATIVE OF a_dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_2dd, a_3dd
    if (GETDERADD(4) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(4,k) = (-2.d0 * gl%A2_PCSAFTD(4)+ 2.d0 * gl%A3_PCSAFTD(4)) &
                & * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k))*gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires a_2dd, a_3dd
     if (GETDERADD(5) .eq. 1) then
        do k = 1, gl%ncomp
            !ADDX1_PCSAFTD(5,k) = (-2.d0 * A2_PCSAFTD(4)**2 * (A2X1_PCSAFTD(1,k) - A3X1_PCSAFTD(1,k)) &
            !    & + 4.d0 * A2_PCSAFTD(4) * A2X1_PCSAFTD(4,k) * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) &
            !    & + 4.d0 * A2_PCSAFTD(4) * (A2_PCSAFTD(4) - A3_PCSAFTD(4)) * A2X1_PCSAFTD(1,k) &
            !    & - 4.d0 * A2_PCSAFTD(4) * (A2X1_PCSAFTD(4,k) - A3X1_PCSAFTD(4,k)) * A2_PCSAFTD(1) &
            !    & + 2.d0 * A2_PCSAFTD(5) * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2X1_PCSAFTD(1,k) &
            !    & - 2.d0 * A2_PCSAFTD(5) * (A2X1_PCSAFTD(1,k) - A3X1_PCSAFTD(1,k)) * A2_PCSAFTD(1) &
            !    & + 4.d0 * A2X1_PCSAFTD(4,k) * (A2_PCSAFTD(4) - A3_PCSAFTD(4)) * A2_PCSAFTD(1) &
            !    & + 2.d0 * A2X1_PCSAFTD(5,k) * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2_PCSAFTD(1) &
            !    & + 6.d0 * (A2_PCSAFTD(4) - A3_PCSAFTD(4))**2 * (2.d0 * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2X1_PCSAFTD(1,k) &
            !    & - (A2X1_PCSAFTD(1,k) - A3X1_PCSAFTD(1,k)) * A2_PCSAFTD(1)) * A2_PCSAFTD(1) &
            !    & / (A2_PCSAFTD(1) - A3_PCSAFTD(1))**2 &
            !    & - 4.d0 * (A2_PCSAFTD(4) - A3_PCSAFTD(4)) * (2.d0 * A2_PCSAFTD(4) * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2X1_PCSAFTD(1,k) &
            !    & - 2.d0 * A2_PCSAFTD(4) * (A2X1_PCSAFTD(1,k) - A3X1_PCSAFTD(1,k)) * A2_PCSAFTD(1) &
            !    & + 2.d0 * A2X1_PCSAFTD(4,k) * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2_PCSAFTD(1) &
            !    & + 2.d0 * (A2_PCSAFTD(4) - A3_PCSAFTD(4)) * A2_PCSAFTD(1) * A2X1_PCSAFTD(1,k) &
            !    & - (A2X1_PCSAFTD(4,k) - A3X1_PCSAFTD(4,k)) * A2_PCSAFTD(1)**2) &
            !    & / (A2_PCSAFTD(1) - A3_PCSAFTD(1)) &
            !    & - 2.d0 * (A2_PCSAFTD(5) - A3_PCSAFTD(5)) * (2.d0 * (A2_PCSAFTD(1) - A3_PCSAFTD(1)) * A2X1_PCSAFTD(1,k) &
            !    & - (A2X1_PCSAFTD(1,k) - A3X1_PCSAFTD(1,k)) * A2_PCSAFTD(1)) * A2_PCSAFTD(1) &
            !    & /(A2_PCSAFTD(1) - A3_PCSAFTD(1)) &
            !    & + 2.d0 * (A2_PCSAFTD(5) - A3_PCSAFTD(5)) * A2_PCSAFTD(1) * A2X1_PCSAFTD(1,k) &
            !    & - (A2X1_PCSAFTD(5,k) - A3X1_PCSAFTD(5,k)) * A2_PCSAFTD(1)**2) &
            !    & / (A2_PCSAFTD(1) - A3_PCSAFTD(1))**2
        
            gl%ADDX1_PCSAFTD(5,k) = (6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (2.d0 *(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))* &
                & gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k))*gl%A2_PCSAFTD(1))*gl%A2_PCSAFTD(1) - 2.d0 *(2.d0*(gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))*(2.d0 *gl%A2_PCSAFTD(4)*(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))*gl%A2X1_PCSAFTD(1,k) - 2.d0 *gl%A2_PCSAFTD(4)*(gl%A2X1_PCSAFTD(1,k) &
                & - gl%A3X1_PCSAFTD(1,k))*gl%A2_PCSAFTD(1) + 2.d0*gl%A2X1_PCSAFTD(4,k)*(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))*gl%A2_PCSAFTD(1) + 2.d0*(gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))*gl%A2_PCSAFTD(1)* &
                & gl%A2X1_PCSAFTD(1,k) + (-gl%A2X1_PCSAFTD(4,k) + gl%A3X1_PCSAFTD(4,k))*gl%A2_PCSAFTD(1)**2) + (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5))*(2.d0*(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & *gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k))*gl%A2_PCSAFTD(1))*gl%A2_PCSAFTD(1))*(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) + (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2* &
                & (2.d0 *gl%A2_PCSAFTD(4)**2*(-gl%A2X1_PCSAFTD(1,k) + gl%A3X1_PCSAFTD(1,k)) + 4.d0*gl%A2_PCSAFTD(4)*gl%A2X1_PCSAFTD(4,k)*(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) + 4.d0*gl%A2_PCSAFTD(4)*(gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))*gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0*gl%A2_PCSAFTD(4)*( gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k))*gl%A2_PCSAFTD(1) + 2.d0*gl%A2_PCSAFTD(5)*(gl%A2_PCSAFTD(1) &
                &  - gl%A3_PCSAFTD(1))*gl%A2X1_PCSAFTD(1,k) - 2.d0*gl%A2_PCSAFTD(5)*(gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k))*gl%A2_PCSAFTD(1) + 4.d0*gl%A2X1_PCSAFTD(4,k)*(gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))*gl%A2_PCSAFTD(1) + 2.d0*gl%A2X1_PCSAFTD(5,k)*(gl%A2_PCSAFTD(1) - &
                & gl%A3_PCSAFTD(1))*gl%A2_PCSAFTD(1) + 2.d0*(gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5))*gl%A2_PCSAFTD(1)*gl%A2X1_PCSAFTD(1,k) + (-gl%A2X1_PCSAFTD(5,k) + gl%A3X1_PCSAFTD(5,k))*gl%A2_PCSAFTD(1)**2))/(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**4
        
        end do
     end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_dd WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_2dd, a_3dd
    if (GETDERADD(6) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(6,k) = (-2.d0 * gl%A2_PCSAFTD(4) + 2.d0 * gl%A3_PCSAFTD(4)) * (-3.d0 * gl%A2_PCSAFTD(2) + 3.d0 * gl%A3_PCSAFTD(2)) &
                & * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**4 &
                & + (-2.d0 * gl%A2_PCSAFTD(4) + 2.d0 * gl%A3_PCSAFTD(4)) &
                & * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
                & /(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTD(6) + 2.d0 * gl%A3_PCSAFTD(6)) &
                & * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTD(2) + 2.d0 * gl%A3_PCSAFTD(2)) &
                & * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))*gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1)- gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(6,k) - gl%A3X1_PCSAFTD(6,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
            end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires a_2dd, a_3dd
    
    if (GETDERADD(7) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(7,k) = (-2.d0 * gl%A2_PCSAFTD(2) + 2.d0 * gl%A3_PCSAFTD(2)) * (-2.d0 * gl%A2_PCSAFTD(4)**2 * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) - 2.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(5,k) - gl%A3X1_PCSAFTD(5,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTD(4)**2 * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) &
                & - 4.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & - 4.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(6,k) - gl%A3X1_PCSAFTD(6,k)) * gl%A2_PCSAFTD(1) &
                & - 2.d0 * gl%A2_PCSAFTD(5) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTD(5) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(7) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(7) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2_PCSAFTD(6) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(6) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 6.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(5,k) - gl%A3X1_PCSAFTD(5,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(7,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2_PCSAFTD(1) &
                & + 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (-2.d0 * gl%A2_PCSAFTD(2) + 2.d0 * gl%A3_PCSAFTD(2)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (-gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & + 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (2.d0 * gl%A2_PCSAFTD(6) - 2.d0 * gl%A3_PCSAFTD(6)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (-gl%A2_PCSAFTD(2) + gl%A3_PCSAFTD(2)) * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 - 4.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (-2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(6,k) - gl%A3X1_PCSAFTD(6,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 2.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (-gl%A2_PCSAFTD(2) + gl%A3_PCSAFTD(2)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 2.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (-gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 2.d0 * (gl%A2_PCSAFTD(7) - gl%A3_PCSAFTD(7)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(7) - gl%A3_PCSAFTD(7)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - (gl%A2X1_PCSAFTD(7,k) - gl%A3X1_PCSAFTD(7,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if
    
    ! I. Schuelling 05/17
    ! 8: 3ND DERIVATIVE OF a_dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires a_2dd, a_3dd
    if (GETDERADD(8) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(8,k) = (-6.d0 * gl%A2_PCSAFTD(2)**2 * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) &
            & - 6.d0 * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
            & + 12.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) &
            & + 6.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & + 6.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2X1_PCSAFTD(1,k) &
            & - 6.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(3,k) - gl%A3X1_PCSAFTD(3,k)) * gl%A2_PCSAFTD(1) &
            & + 6.d0 * gl%A2_PCSAFTD(3) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & + 6.d0 * gl%A2_PCSAFTD(3) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
            & - 6.d0 * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2_PCSAFTD(8) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - 2.d0 * gl%A2_PCSAFTD(8) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
            & + 6.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2_PCSAFTD(1) &
            & + 6.d0 * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2X1_PCSAFTD(8,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & - 24.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))**3 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
            & + 18.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))**2 * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
            & + 18.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
            & - 6.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (-2.d0 * gl%A2_PCSAFTD(2)**2 * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
            & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & + 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
            & - 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
            & + 4.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(3,k) - gl%A3X1_PCSAFTD(3,k)) * gl%A2_PCSAFTD(1)**2) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & - 6.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & - 2.d0 * (gl%A2_PCSAFTD(8) - gl%A3_PCSAFTD(8)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & + 2.d0 * (gl%A2_PCSAFTD(8) - gl%A3_PCSAFTD(8)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
            & - (gl%A2X1_PCSAFTD(8,k) - gl%A3X1_PCSAFTD(8,k)) * gl%A2_PCSAFTD(1)**2) &
            & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if   
    
    ! 9: 3ND DERIVATIVE OF a_dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3
    ! requires a_2dd, a_3dd
    if (GETDERADD(9) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(9,k) =  (-6.d0 * gl%A2_PCSAFTD(4)**2 * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) &
                & - 6.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(5) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 12.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) &
                & + 6.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 6.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2X1_PCSAFTD(1,k) &
                & - 6.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(5,k) - gl%A3X1_PCSAFTD(5,k)) * gl%A2_PCSAFTD(1) &
                & + 6.d0 * gl%A2_PCSAFTD(5) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 6.d0 * gl%A2_PCSAFTD(5) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 6.d0 * gl%A2_PCSAFTD(5) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(9) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(9) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 6.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2_PCSAFTD(1) &
                & + 6.0d0 * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(9,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & - 24.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**3 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + 18.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))**2 * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & + 18.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (-2.d0 * gl%A2_PCSAFTD(4)**2 * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(5) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(5,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(5,k) - gl%A3X1_PCSAFTD(5,k)) * gl%A2_PCSAFTD(1)**2) & 
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 6.d0 * (gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5)) * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 2.d0 * (gl%A2_PCSAFTD(9) - gl%A3_PCSAFTD(9)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(9) - gl%A3_PCSAFTD(9)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(9,k) - gl%A3X1_PCSAFTD(9,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2    
        end do
    end if   
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires a_2dd, a_3dd
    if (GETDERADD(10) .eq. 1) then
        do k = 1, gl%ncomp
            gl%ADDX1_PCSAFTD(10,k) = (-4.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(3,k) - gl%A3X1_PCSAFTD(3,k)) * gl%A2_PCSAFTD(1) &
                & - 4.d0 * gl%A2_PCSAFTD(6) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(6) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(6) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(10) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(10) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & - 2.d0 * gl%A2_PCSAFTD(2)**2 * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) & 
                & + 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(6,k) - gl%A3X1_PCSAFTD(6,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(3) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(10,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) &
                & - 24.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))**2 * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) & 
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 &
                & + 12.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & + 6.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * (-2.d0 * gl%A2_PCSAFTD(2)**2 * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 4.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(3) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 4.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(3,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(3,k) - gl%A3X1_PCSAFTD(3,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 12.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * (2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - 2.d0 * (gl%A2_PCSAFTD(10) - gl%A3_PCSAFTD(10)) * (2.d0 * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * (gl%A2_PCSAFTD(10) - gl%A3_PCSAFTD(10)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & + 6.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))**2 * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * (-2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(2,k) - gl%A3X1_PCSAFTD(2,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(6) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & + 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(2) * (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(6,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(2,k) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))* gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(6) - gl%A3_PCSAFTD(6)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(6,k) - gl%A3X1_PCSAFTD(6,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) - 2.d0 * (gl%A2_PCSAFTD(3) - gl%A3_PCSAFTD(3)) * (2.d0 * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X1_PCSAFTD(1,k) &
                & - 2.d0 * gl%A2_PCSAFTD(4) * (gl%A2X1_PCSAFTD(1,k) - gl%A3X1_PCSAFTD(1,k)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * gl%A2X1_PCSAFTD(4,k) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2_PCSAFTD(1) &
                & + 2.d0 * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) * gl%A2_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,k) &
                & - (gl%A2X1_PCSAFTD(4,k) - gl%A3X1_PCSAFTD(4,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
                & - (gl%A2X1_PCSAFTD(10,k) - gl%A3X1_PCSAFTD(10,k)) * gl%A2_PCSAFTD(1)**2) &
                & / (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**2
        end do
    end if

!DEC$ END IF
end subroutine ADDX1DERIVS

subroutine ADDX2DERIVS(gl,T,DENS,GETDERADD)
    
    ! a_dd: contribution to the residual helmholz engergy from dipolar interactions 
    ! defined by eq. A.7 in Gross, Sadowski 2006:
    ! a_dd = a_2dd*a_2dd/(a_2dd-a_3dd)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    integer :: o, u
    !output: add_PCSAFT (module variable)
    !variable for the numerical derivative
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision, dimension(gl%ncomp,gl%ncomp) :: ADDX2, ADDX2DP, ADDX2DP2, ADDX2DM, ADDX2DM2, ADDX2TP, ADDX2TP2, ADDX2TM, ADDX2TM2, ADDX2PP, ADDX2PM, ADDX2MP, ADDX2MM

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    
    DELTA = 1.d-5
    
    DP = DENS * (1.d0 + DELTA)
    DP2 = DENS * (1.d0 + 2.d0*DELTA)
    DM = DENS * (1.d0 - DELTA)
    DM2 = DENS * (1.d0 - 2.d0*DELTA)
    TP = T * (1.d0 + DELTA)
    TP2 = T * (1.d0 + 2.d0*DELTA)
    TM = T * (1.d0 - DELTA)
    TM2 = T * (1.d0 - 2.d0*DELTA)
    
    DEL = 5.d-4 !gut für 2. Temperauturableitung
    
    DPx = DENS * (1.d0 + DEL)
    DP2x = DENS * (1.d0 + 2.d0*DEL)
    DMx = DENS * (1.d0 - DEL)
    DM2x = DENS * (1.d0 - 2.d0*DEL)
    TPx = T * (1.d0 + DEL)
    TP2x = T * (1.d0 + 2.d0*DEL)
    TMx = T * (1.d0 - DEL)
    TM2x = T * (1.d0 - 2.d0*DEL)
    
    DELT = 1.d-4 
    
    DPxx = DENS * (1.d0 + DELT)
    DP2xx = DENS * (1.d0 + 2.d0*DELT)
    DMxx = DENS * (1.d0 - DELT)
    DM2xx = DENS * (1.d0 - 2.d0*DELT)
    TPxx = T * (1.d0 + DELT)
    TP2xx = T * (1.d0 + 2.d0*DELT)
    TMxx = T * (1.d0 - DELT)
    TM2xx = T * (1.d0 - 2.d0*DELT)
    
    !calculate the derivatives of a_dd
    if (GETDERADD(1) .eq. 1) then
        
        call base_ADDX2(gl,T,DENS,ADDX2)
        
        gl%ADDX2_PCSAFTD(1,:,:) = ADDX2
        
    end if
    
     !  2: 1ST DERIVATIVE OF a_2dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_dd
    if (GETDERADD(2) .eq. 1) then   
               
        call base_ADDX2(gl,T,DP,ADDX2DP)
        call base_ADDX2(gl,T,DP2,ADDX2DP2)
        call base_ADDX2(gl,T,DM,ADDX2DM)
        call base_ADDX2(gl,T,DM2,ADDX2DM2)
        
        gl%ADDX2_PCSAFTD(2,:,:) = (8.d0*ADDX2DP - 8.d0*ADDX2DM - ADDX2DP2 + ADDX2DM2) / (12.d0 * DELTA) 
    

    end if
    

   ! 3: 2ND DERIVATIVE OF a_2dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires a_dd
    if (GETDERADD(3) .eq. 1) then
        
        call base_ADDX2(gl,T,DENS,ADDX2)    
        call base_ADDX2(gl,T,DPx,ADDX2DP)
        call base_ADDX2(gl,T,DP2x,ADDX2DP2)
        call base_ADDX2(gl,T,DMx,ADDX2DM)
        call base_ADDX2(gl,T,DM2x,ADDX2DM2)
        
        gl%ADDX2_PCSAFTD(3,:,:) = (-ADDX2DM2 + 16.d0*ADDX2DM - 30.d0*ADDX2 + 16.d0*ADDX2DP - ADDX2DP2) / (12.d0*DEL**2) 
    
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_dd
    if (GETDERADD(4) .eq. 1) then
        call base_ADDX2(gl,TP,DENS,ADDX2TP)
        call base_ADDX2(gl,TP2,DENS,ADDX2TP2)
        call base_ADDX2(gl,TM,DENS,ADDX2TM)
        call base_ADDX2(gl,TM2,DENS,ADDX2TM2)
        
        gl%ADDX2_PCSAFTD(4,:,:) = (8.d0*ADDX2TP - 8.d0*ADDX2TM - ADDX2TP2 + ADDX2TM2) / (12.d0 * DELTA) 
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires a_dd
    if (GETDERADD(5) .eq. 1) then
    
        call base_ADDX2(gl,T,DENS,ADDX2)    
        call base_ADDX2(gl,TPx,DENS,ADDX2TP)
        call base_ADDX2(gl,TP2x,DENS,ADDX2TP2)
        call base_ADDX2(gl,TMx,DENS,ADDX2TM)
        call base_ADDX2(gl,TM2x,DENS,ADDX2TM2)
        
        gl%ADDX2_PCSAFTD(5,:,:) = (-ADDX2TM2 + 16.d0*ADDX2TM - 30.d0*ADDX2 + 16.d0*ADDX2TP - ADDX2TP2) / (12.d0*DEL**2) !https://en.wikipedia.org/wiki/Five-point_stencil
    
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2dd WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_dd
    if (GETDERADD(6) .eq. 1) then
    
        call base_ADDX2(gl,TPxx,DPxx,ADDX2PP)
        call base_ADDX2(gl,TMxx,DPxx,ADDX2MP)
        call base_ADDX2(gl,TMxx,DMxx,ADDX2MM)
        call base_ADDX2(gl,TPxx,DMxx,ADDX2PM)
        
        gl%ADDX2_PCSAFTD(6,:,:) = (ADDX2PP - ADDX2PM - ADDX2MP + ADDX2MM) / (4.d0*DELT**2) 
    end if
    
!DEC$ END IF
end subroutine ADDX2DERIVS

subroutine base_ADDX2(gl,T,D,ADDX2)
    
    ! a_dd: contribution to the residual helmholz engergy from dipolar interactions 
    ! defined by eq. A.7 in Gross, Sadowski 2006:
    ! a_dd = a_2dd*a_2dd/(a_2dd-a_3dd)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    double precision :: T, D, A2_PCSAFTD_, A3_PCSAFTD_
    double precision, dimension(gl%ncomp) :: A2X1_PCSAFTD_, A3X1_PCSAFTD_
    integer :: o, u
    double precision, dimension(gl%ncomp,gl%ncomp) :: ADDX2
    integer, dimension(nderivs) :: getprevious
    integer :: nrsubst
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !output: add_PCSAFT (module variable)
    !1.866820186682000D-005
    getprevious = 1
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)
    A2_PCSAFTD_ = gl%A2_PCSAFTD(1)
    A3_PCSAFTD_ = gl%A3_PCSAFTD(1)
    A2X1_PCSAFTD_ = gl%A2X1_PCSAFTD(1,:)
    A3X1_PCSAFTD_ = gl%A3X1_PCSAFTD(1,:)
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x2_DD(gl,T,D,getprevious)
    
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)    
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
    
    !calculate the base function of the second composition derivative of a_dd
        do o = 1, gl%ncomp
			do u = 1, gl%ncomp  
                if (u == o) then !-1.843408237705531D-003
                !    ADDX2(o,u) = (A2_PCSAFTD_ * ((A2_PCSAFTD_ - 2.d0 * A3_PCSAFTD_) * (A2_PCSAFTD_ - A3_PCSAFTD_) * A2X2_PCSAFTD(1,o,u) + A2_PCSAFTD_ &
                !                            & * (A2_PCSAFTD_ * A3X2_PCSAFTD(1,o,u) - A3_PCSAFTD_ * A3X2_PCSAFTD(1,o,u) + 2.d0 * A3X1_PCSAFTD_(o)**2)) &
                !                            & - 4.d0 * A2_PCSAFTD_ * A3_PCSAFTD_ * A2X1_PCSAFTD_(o) * A3X1_PCSAFTD_(o) + 2.d0 * A3_PCSAFTD_**2 * A2X1_PCSAFTD_(o)**2) &
                !                            & /(A2_PCSAFTD_ - A3_PCSAFTD_)**3 
                !else                                                                       
                !    ADDX2(o,u) =-(2.d0 * A3_PCSAFTD_ * A2X1_PCSAFTD_(u) * (A2_PCSAFTD_ * A3X1_PCSAFTD_(o) - A3_PCSAFTD_ * A2X1_PCSAFTD_(o)) + A2_PCSAFTD_ * (- A2_PCSAFTD_ &
                !                            & * (2.d0 * A3X1_PCSAFTD_(o) * A3X1_PCSAFTD_(u) + A3_PCSAFTD_ * (-A3X2_PCSAFTD(1,o,u) - 3.d0 * A2X2_PCSAFTD(1,o,u))) &
                !                            & + A2_PCSAFTD_ * A2_PCSAFTD_ * (- (A3X2_PCSAFTD(1,o,u) + A2X2_PCSAFTD(1,o,u))) + 2.d0 * A3_PCSAFTD_ * (A2X1_PCSAFTD_(o) * A3X1_PCSAFTD_(u) - A3_PCSAFTD_ * A2X2_PCSAFTD(1,o,u)))) &
                !                            & / (A2_PCSAFTD(1)-A3_PCSAFTD(1))**3 
                    ADDX2(o,u) = (gl%A2_PCSAFTD(1) * ((gl%A2_PCSAFTD(1) - 2.d0 * gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * gl%A2X2_PCSAFTD(1,o,u) + gl%A2_PCSAFTD(1) &
                                            & * (gl%A2_PCSAFTD(1) * gl%A3X2_PCSAFTD(1,o,u) - gl%A3_PCSAFTD(1) * gl%A3X2_PCSAFTD(1,o,u) + 2.d0 * gl%A3X1_PCSAFTD(1,o)**2)) &
                                            & - 4.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,o) * gl%A3X1_PCSAFTD(1,o) + 2.d0 * gl%A3_PCSAFTD(1)**2 * gl%A2X1_PCSAFTD(1,o)**2) &
                                            & /(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))**3 
                else                                                                       
                    ADDX2(o,u) =-(2.d0 * gl%A3_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,u) * (gl%A2_PCSAFTD(1) * gl%A3X1_PCSAFTD(1,o) - gl%A3_PCSAFTD(1) * gl%A2X1_PCSAFTD(1,o)) + gl%A2_PCSAFTD(1) * (- gl%A2_PCSAFTD(1) &
                                            & * (2.d0 * gl%A3X1_PCSAFTD(1,o) * gl%A3X1_PCSAFTD(1,u) + gl%A3_PCSAFTD(1) * (-gl%A3X2_PCSAFTD(1,o,u) - 3.d0 * gl%A2X2_PCSAFTD(1,o,u))) &
                                            & + gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * (- (gl%A3X2_PCSAFTD(1,o,u) + gl%A2X2_PCSAFTD(1,o,u))) + 2.d0 * gl%A3_PCSAFTD(1) * (gl%A2X1_PCSAFTD(1,o) * gl%A3X1_PCSAFTD(1,u) - gl%A3_PCSAFTD(1) * gl%A2X2_PCSAFTD(1,o,u)))) &
                                            & / (gl%A2_PCSAFTD(1)-gl%A3_PCSAFTD(1))**3 
                end if
			end do
		end do

!DEC$ END IF
end subroutine base_ADDX2


subroutine AQQX1DERIVS(gl,T,DENS,GETDERAQQ)
    
    ! a_qq: contribution to the residual helmholz engergy from quadrupol-quadrupol interactions 
    ! defined by eq. A.8 in Gross, Sadowski 2005:
    ! a_qq = a_2qq*a_2qq/(a_2qq-a_3qq)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    integer :: k
    !output: aqq_PCSAFT (module variable)
    !working variables

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    
    !calculate the derivatives of a_qq
    if (GETDERAQQ(1) .eq. 1) then
        do k = 1, gl%ncomp
        gl%AQQX1_PCSAFTQ(1,k) = (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)**2) / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(2) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(2,k) = (-2.d0 * gl%A2_PCSAFTQ(2) + 2.d0 * gl%A3_PCSAFTQ(2)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))&
                & * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k))*  gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 + (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))*gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if
    ! I. Schuelling 05/17
    ! 3: 2ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(3) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(3,k) = (-2.0 * gl%A2_PCSAFTQ(2)**2 * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 6.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2))**2 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1))*gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 2.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(3,k) - gl%A3X1_PCSAFTQ(3,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
     end if
     
    ! 4: 1ST DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(4) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(4,k) = (-2.d0 * gl%A2_PCSAFTQ(4)+ 2.d0 * gl%A3_PCSAFTQ(4)) &
                & * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k))*gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires a_2qq, a_3qq
     if (GETDERAQQ(5) .eq. 1) then
        do k = 1, gl%ncomp
            !AQQX1_PCSAFTQ(5,k) = (-2.d0 * A2_PCSAFTQ(4)**2 * (A2X1_PCSAFTQ(1,k) - A3X1_PCSAFTQ(1,k)) &
            !    & + 4.d0 * A2_PCSAFTQ(4) * A2X1_PCSAFTQ(4,k) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) &
            !    & + 4.d0 * A2_PCSAFTQ(4) * (A2_PCSAFTQ(4) - A3_PCSAFTQ(4)) * A2X1_PCSAFTQ(1,k) &
            !    & - 4.d0 * A2_PCSAFTQ(4) * (A2X1_PCSAFTQ(4,k) - A3X1_PCSAFTQ(4,k)) * A2_PCSAFTQ(1) &
            !    & + 2.d0 * A2_PCSAFTQ(5) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2X1_PCSAFTQ(1,k) &
            !    & - 2.d0 * A2_PCSAFTQ(5) * (A2X1_PCSAFTQ(1,k) - A3X1_PCSAFTQ(1,k)) * A2_PCSAFTQ(1) &
            !    & + 4.d0 * A2X1_PCSAFTQ(4,k) * (A2_PCSAFTQ(4) - A3_PCSAFTQ(4)) * A2_PCSAFTQ(1) &
            !    & + 2.d0 * A2X1_PCSAFTQ(5,k) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2_PCSAFTQ(1) &
            !    & + 6.d0 * (A2_PCSAFTQ(4) - A3_PCSAFTQ(4))**2 * (2.d0 * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2X1_PCSAFTQ(1,k) &
            !    & - (A2X1_PCSAFTQ(1,k) - A3X1_PCSAFTQ(1,k)) * A2_PCSAFTQ(1)) * A2_PCSAFTQ(1) &
            !    & / (A2_PCSAFTQ(1) - A3_PCSAFTQ(1))**2 &
            !    & - 4.d0 * (A2_PCSAFTQ(4) - A3_PCSAFTQ(4)) * (2.d0 * A2_PCSAFTQ(4) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2X1_PCSAFTQ(1,k) &
            !    & - 2.d0 * A2_PCSAFTQ(4) * (A2X1_PCSAFTQ(1,k) - A3X1_PCSAFTQ(1,k)) * A2_PCSAFTQ(1) &
            !    & + 2.d0 * A2X1_PCSAFTQ(4,k) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2_PCSAFTQ(1) &
            !    & + 2.d0 * (A2_PCSAFTQ(4) - A3_PCSAFTQ(4)) * A2_PCSAFTQ(1) * A2X1_PCSAFTQ(1,k) &
            !    & - (A2X1_PCSAFTQ(4,k) - A3X1_PCSAFTQ(4,k)) * A2_PCSAFTQ(1)**2) &
            !    & / (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) &
            !    & - 2.d0 * (A2_PCSAFTQ(5) - A3_PCSAFTQ(5)) * (2.d0 * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2X1_PCSAFTQ(1,k) &
            !    & - (A2X1_PCSAFTQ(1,k) - A3X1_PCSAFTQ(1,k)) * A2_PCSAFTQ(1)) * A2_PCSAFTQ(1) &
            !    & /(A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) &
            !    & + 2.d0 * (A2_PCSAFTQ(5) - A3_PCSAFTQ(5)) * A2_PCSAFTQ(1) * A2X1_PCSAFTQ(1,k) &
            !    & - (A2X1_PCSAFTQ(5,k) - A3X1_PCSAFTQ(5,k)) * A2_PCSAFTQ(1)**2) &
            !    & / (A2_PCSAFTQ(1) - A3_PCSAFTQ(1))**2
        
            gl%AQQX1_PCSAFTQ(5,k) = (6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (2.d0 *(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))* &
                & gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k))*gl%A2_PCSAFTQ(1))*gl%A2_PCSAFTQ(1) - 2.d0 *(2.d0*(gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))*(2.d0 *gl%A2_PCSAFTQ(4)*(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))*gl%A2X1_PCSAFTQ(1,k) - 2.d0 *gl%A2_PCSAFTQ(4)*(gl%A2X1_PCSAFTQ(1,k) &
                & - gl%A3X1_PCSAFTQ(1,k))*gl%A2_PCSAFTQ(1) + 2.d0*gl%A2X1_PCSAFTQ(4,k)*(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))*gl%A2_PCSAFTQ(1) + 2.d0*(gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))*gl%A2_PCSAFTQ(1)* &
                & gl%A2X1_PCSAFTQ(1,k) + (-gl%A2X1_PCSAFTQ(4,k) + gl%A3X1_PCSAFTQ(4,k))*gl%A2_PCSAFTQ(1)**2) + (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5))*(2.d0*(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & *gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k))*gl%A2_PCSAFTQ(1))*gl%A2_PCSAFTQ(1))*(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) + (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2* &
                & (2.d0 *gl%A2_PCSAFTQ(4)**2*(-gl%A2X1_PCSAFTQ(1,k) + gl%A3X1_PCSAFTQ(1,k)) + 4.d0*gl%A2_PCSAFTQ(4)*gl%A2X1_PCSAFTQ(4,k)*(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) + 4.d0*gl%A2_PCSAFTQ(4)*(gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))*gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0*gl%A2_PCSAFTQ(4)*( gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k))*gl%A2_PCSAFTQ(1) + 2.d0*gl%A2_PCSAFTQ(5)*(gl%A2_PCSAFTQ(1) &
                &  - gl%A3_PCSAFTQ(1))*gl%A2X1_PCSAFTQ(1,k) - 2.d0*gl%A2_PCSAFTQ(5)*(gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k))*gl%A2_PCSAFTQ(1) + 4.d0*gl%A2X1_PCSAFTQ(4,k)*(gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))*gl%A2_PCSAFTQ(1) + 2.d0*gl%A2X1_PCSAFTQ(5,k)*(gl%A2_PCSAFTQ(1) - &
                & gl%A3_PCSAFTQ(1))*gl%A2_PCSAFTQ(1) + 2.d0*(gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5))*gl%A2_PCSAFTQ(1)*gl%A2X1_PCSAFTQ(1,k) + (-gl%A2X1_PCSAFTQ(5,k) + gl%A3X1_PCSAFTQ(5,k))*gl%A2_PCSAFTQ(1)**2))/(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**4
        
        end do
     end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(6) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(6,k) = (-2.d0 * gl%A2_PCSAFTQ(4) + 2.d0 * gl%A3_PCSAFTQ(4)) * (-3.d0 * gl%A2_PCSAFTQ(2) + 3.d0 * gl%A3_PCSAFTQ(2)) &
                & * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**4 &
                & + (-2.d0 * gl%A2_PCSAFTQ(4) + 2.d0 * gl%A3_PCSAFTQ(4)) &
                & * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
                & /(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTQ(6) + 2.d0 * gl%A3_PCSAFTQ(6)) &
                & * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTQ(2) + 2.d0 * gl%A3_PCSAFTQ(2)) &
                & * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))*gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1)- gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(6,k) - gl%A3X1_PCSAFTQ(6,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
            end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires a_2qq, a_3qq
    
    if (GETDERAQQ(7) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(7,k) = (-2.d0 * gl%A2_PCSAFTQ(2) + 2.d0 * gl%A3_PCSAFTQ(2)) * (-2.d0 * gl%A2_PCSAFTQ(4)**2 * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) - 2.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(5,k) - gl%A3X1_PCSAFTQ(5,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + (-2.d0 * gl%A2_PCSAFTQ(4)**2 * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) &
                & - 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & - 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(6,k) - gl%A3X1_PCSAFTQ(6,k)) * gl%A2_PCSAFTQ(1) &
                & - 2.d0 * gl%A2_PCSAFTQ(5) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTQ(5) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(7) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(7) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2_PCSAFTQ(6) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(6) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 6.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(5,k) - gl%A3X1_PCSAFTQ(5,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(7,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2_PCSAFTQ(1) &
                & + 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (-2.d0 * gl%A2_PCSAFTQ(2) + 2.d0 * gl%A3_PCSAFTQ(2)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (-gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & + 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (2.d0 * gl%A2_PCSAFTQ(6) - 2.d0 * gl%A3_PCSAFTQ(6)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (-gl%A2_PCSAFTQ(2) + gl%A3_PCSAFTQ(2)) * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 - 4.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (-2.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(6,k) - gl%A3X1_PCSAFTQ(6,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 2.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (-gl%A2_PCSAFTQ(2) + gl%A3_PCSAFTQ(2)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 2.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (-gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 2.d0 * (gl%A2_PCSAFTQ(7) - gl%A3_PCSAFTQ(7)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(7) - gl%A3_PCSAFTQ(7)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - (gl%A2X1_PCSAFTQ(7,k) - gl%A3X1_PCSAFTQ(7,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if
    
    ! I. Schuelling 05/17
    ! 8: 3ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(8) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(8,k) = (-6.d0 * gl%A2_PCSAFTQ(2)**2 * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) &
            & - 6.d0 * gl%A2_PCSAFTQ(2) * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
            & + 12.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) &
            & + 6.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & + 6.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 6.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(3,k) - gl%A3X1_PCSAFTQ(3,k)) * gl%A2_PCSAFTQ(1) &
            & + 6.d0 * gl%A2_PCSAFTQ(3) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & + 6.d0 * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 6.d0 * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2_PCSAFTQ(8) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 2.d0 * gl%A2_PCSAFTQ(8) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
            & + 6.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2_PCSAFTQ(1) &
            & + 6.d0 * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2X1_PCSAFTQ(8,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & - 24.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2))**3 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
            & + 18.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2))**2 * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
            & + 18.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
            & - 6.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (-2.d0 * gl%A2_PCSAFTQ(2)**2 * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
            & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & + 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
            & + 4.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(3,k) - gl%A3X1_PCSAFTQ(3,k)) * gl%A2_PCSAFTQ(1)**2) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & - 6.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & - 2.d0 * (gl%A2_PCSAFTQ(8) - gl%A3_PCSAFTQ(8)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
            & + 2.d0 * (gl%A2_PCSAFTQ(8) - gl%A3_PCSAFTQ(8)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
            & - (gl%A2X1_PCSAFTQ(8,k) - gl%A3X1_PCSAFTQ(8,k)) * gl%A2_PCSAFTQ(1)**2) &
            & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if   
    
    ! 9: 3ND DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(9) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(9,k) =  (-6.d0 * gl%A2_PCSAFTQ(4)**2 * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) &
                & - 6.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(5) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 12.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) &
                & + 6.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 6.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 6.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(5,k) - gl%A3X1_PCSAFTQ(5,k)) * gl%A2_PCSAFTQ(1) &
                & + 6.d0 * gl%A2_PCSAFTQ(5) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 6.d0 * gl%A2_PCSAFTQ(5) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 6.d0 * gl%A2_PCSAFTQ(5) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(9) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(9) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 6.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2_PCSAFTQ(1) &
                & + 6.0d0 * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(9,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & - 24.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**3 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + 18.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))**2 * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & + 18.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (-2.d0 * gl%A2_PCSAFTQ(4)**2 * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(5) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(5,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(5,k) - gl%A3X1_PCSAFTQ(5,k)) * gl%A2_PCSAFTQ(1)**2) & 
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 6.d0 * (gl%A2_PCSAFTQ(5) - gl%A3_PCSAFTQ(5)) * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 2.d0 * (gl%A2_PCSAFTQ(9) - gl%A3_PCSAFTQ(9)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(9) - gl%A3_PCSAFTQ(9)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(9,k) - gl%A3X1_PCSAFTQ(9,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2    
        end do
    end if   
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(10) .eq. 1) then
        do k = 1, gl%ncomp
            gl%AQQX1_PCSAFTQ(10,k) = (-4.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(3,k) - gl%A3X1_PCSAFTQ(3,k)) * gl%A2_PCSAFTQ(1) &
                & - 4.d0 * gl%A2_PCSAFTQ(6) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(6) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(6) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(10) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(10) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & - 2.d0 * gl%A2_PCSAFTQ(2)**2 * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) & 
                & + 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(6,k) - gl%A3X1_PCSAFTQ(6,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(3) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(10,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) &
                & - 24.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2))**2 * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) & 
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 &
                & + 12.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & + 6.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * (-2.d0 * gl%A2_PCSAFTQ(2)**2 * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 4.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(3) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 4.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(3,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(3,k) - gl%A3X1_PCSAFTQ(3,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 12.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * (2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - 2.d0 * (gl%A2_PCSAFTQ(10) - gl%A3_PCSAFTQ(10)) * (2.d0 * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * (gl%A2_PCSAFTQ(10) - gl%A3_PCSAFTQ(10)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & + 6.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2))**2 * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2 &
                & - 4.d0 * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * (-2.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(2,k) - gl%A3X1_PCSAFTQ(2,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(6) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & + 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(2) * (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(2) - gl%A3_PCSAFTQ(2)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(6,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(2,k) * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4))* gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(6) - gl%A3_PCSAFTQ(6)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(6,k) - gl%A3X1_PCSAFTQ(6,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) - 2.d0 * (gl%A2_PCSAFTQ(3) - gl%A3_PCSAFTQ(3)) * (2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X1_PCSAFTQ(1,k) &
                & - 2.d0 * gl%A2_PCSAFTQ(4) * (gl%A2X1_PCSAFTQ(1,k) - gl%A3X1_PCSAFTQ(1,k)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * gl%A2X1_PCSAFTQ(4,k) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2_PCSAFTQ(1) &
                & + 2.d0 * (gl%A2_PCSAFTQ(4) - gl%A3_PCSAFTQ(4)) * gl%A2_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,k) &
                & - (gl%A2X1_PCSAFTQ(4,k) - gl%A3X1_PCSAFTQ(4,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) &
                & - (gl%A2X1_PCSAFTQ(10,k) - gl%A3X1_PCSAFTQ(10,k)) * gl%A2_PCSAFTQ(1)**2) &
                & / (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**2
        end do
    end if

!DEC$ END IF
end subroutine AQQX1DERIVS

subroutine AQQX2DERIVS(gl,T,DENS,GETDERAQQ)
    
    ! a_qq: contribution to the residual helmholz engergy from quadrupol-quadrupol interactions 
    ! defined by eq. A.8 in Gross, Sadowski 2005:
    ! a_qq = a_2qq*a_2qq/(a_2qq-a_3qq)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    integer :: o, u
    !output: aqq_PCSAFT (module variable)
    !variable for the numerical derivative
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision, dimension(gl%ncomp,gl%ncomp) :: AQQX2, AQQX2DP, AQQX2DP2, AQQX2DM, AQQX2DM2, AQQX2TP, AQQX2TP2, AQQX2TM, AQQX2TM2, AQQX2PP, AQQX2PM, AQQX2MP, AQQX2MM

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    DELTA = 5.d-7
    DELTA = 1.d-5
    
    DP = DENS * (1.d0 + DELTA)
    DP2 = DENS * (1.d0 + 2.d0*DELTA)
    DM = DENS * (1.d0 - DELTA)
    DM2 = DENS * (1.d0 - 2.d0*DELTA)
    TP = T * (1.d0 + DELTA)
    TP2 = T * (1.d0 + 2.d0*DELTA)
    TM = T * (1.d0 - DELTA)
    TM2 = T * (1.d0 - 2.d0*DELTA)
    
    DEL = 5.d-5 !gut für 2. Temperauturableitung
    
    DPx = DENS * (1.d0 + DEL)
    DP2x = DENS * (1.d0 + 2.d0*DEL)
    DMx = DENS * (1.d0 - DEL)
    DM2x = DENS * (1.d0 - 2.d0*DEL)
    TPx = T * (1.d0 + DEL)
    TP2x = T * (1.d0 + 2.d0*DEL)
    TMx = T * (1.d0 - DEL)
    TM2x = T * (1.d0 - 2.d0*DEL)
    
    DELT = 1.d-4
    
    DPxx = DENS * (1.d0 + DELT)
    DP2xx = DENS * (1.d0 + 2.d0*DELT)
    DMxx = DENS * (1.d0 - DELT)
    DM2xx = DENS * (1.d0 - 2.d0*DELT)
    TPxx = T * (1.d0 + DELT)
    TP2xx = T * (1.d0 + 2.d0*DELT)
    TMxx = T * (1.d0 - DELT)
    TM2xx = T * (1.d0 - 2.d0*DELT)
    
    !calculate the derivatives of a_qq
    if (GETDERAQQ(1) .eq. 1) then

        call base_AQQX2(gl,T,DENS,AQQX2)
        
        gl%AQQX2_PCSAFTQ(1,:,:) = AQQX2
        
    end if
    
     !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_qq
    if (GETDERAQQ(2) .eq. 1) then   
               
        call base_AQQX2(gl,T,DP,AQQX2DP)
        call base_AQQX2(gl,T,DP2,AQQX2DP2)
        call base_AQQX2(gl,T,DM,AQQX2DM)
        call base_AQQX2(gl,T,DM2,AQQX2DM2)
        
        gl%AQQX2_PCSAFTQ(2,:,:) = (8.d0*AQQX2DP - 8.d0*AQQX2DM - AQQX2DP2 + AQQX2DM2) / (12.d0 * DELTA) 
    
    end if
    

   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires a_qq
    if (GETDERAQQ(3) .eq. 1) then
        
        call base_AQQX2(gl,T,DENS,AQQX2)    
        call base_AQQX2(gl,T,DPx,AQQX2DP)
        call base_AQQX2(gl,T,DP2x,AQQX2DP2)
        call base_AQQX2(gl,T,DMx,AQQX2DM)
        call base_AQQX2(gl,T,DM2x,AQQX2DM2)
        
        gl%AQQX2_PCSAFTQ(3,:,:) = (-AQQX2DM2 + 16.d0*AQQX2DM - 30.d0*AQQX2 + 16.d0*AQQX2DP - AQQX2DP2) / (12.d0*DEL**2) 
    
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_qq
    if (GETDERAQQ(4) .eq. 1) then
        call base_AQQX2(gl,TP,DENS,AQQX2TP)
        call base_AQQX2(gl,TP2,DENS,AQQX2TP2)
        call base_AQQX2(gl,TM,DENS,AQQX2TM)
        call base_AQQX2(gl,TM2,DENS,AQQX2TM2)
        
        gl%AQQX2_PCSAFTQ(4,:,:) = (8.d0*AQQX2TP - 8.d0*AQQX2TM - AQQX2TP2 + AQQX2TM2) / (12.d0 * DELTA) 
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires a_qq
    if (GETDERAQQ(5) .eq. 1) then
    
        call base_AQQX2(gl,T,DENS,AQQX2)    
        call base_AQQX2(gl,TPx,DENS,AQQX2TP)
        call base_AQQX2(gl,TP2x,DENS,AQQX2TP2)
        call base_AQQX2(gl,TMx,DENS,AQQX2TM)
        call base_AQQX2(gl,TM2x,DENS,AQQX2TM2)
        
        gl%AQQX2_PCSAFTQ(5,:,:) = (-AQQX2TM2 + 16.d0*AQQX2TM - 30.d0*AQQX2 + 16.d0*AQQX2TP - AQQX2TP2) / (12.d0*DEL**2) !https://en.wikipedia.org/wiki/Five-point_stencil
    
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_qq
    if (GETDERAQQ(6) .eq. 1) then
    
        call base_AQQX2(gl,TPxx,DP,AQQX2PP)
        call base_AQQX2(gl,TMxx,DP,AQQX2MP)
        call base_AQQX2(gl,TMxx,DM,AQQX2MM)
        call base_AQQX2(gl,TPxx,DM,AQQX2PM)
        
        gl%AQQX2_PCSAFTQ(6,:,:) = (AQQX2PP - AQQX2PM - AQQX2MP + AQQX2MM) / (4.d0*DELT*DELTA) 
    end if
    
!DEC$ END IF
end subroutine AQQX2DERIVS

subroutine base_AQQX2(gl,T,D,AQQX2)
    
    ! a_qq: contribution to the residual helmholz engergy from quadrupol-quadrupol interactions 
    ! defined by eq. A.8 in Gross, Sadowski 2005:
    ! a_qq = a_2qq*a_2qq/(a_2qq-a_3qq)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    double precision :: T, D, A2_PCSAFTQ_, A3_PCSAFTQ_
    integer :: o, u
    double precision, dimension(gl%ncomp,gl%ncomp) :: AQQX2
    double precision, dimension(gl%ncomp) :: A2X1_PCSAFTQ_, A3X1_PCSAFTQ_
    integer, dimension(nderivs) :: getprevious
    integer :: nrsubst
    !output: aqq_PCSAFT (module variable)
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    getprevious = 1
    getprevious = 1
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)
    A2_PCSAFTQ_ = gl%A2_PCSAFTQ(1)
    A3_PCSAFTQ_ = gl%A3_PCSAFTQ(1)
    A2X1_PCSAFTQ_ = gl%A2X1_PCSAFTQ(1,:)
    A3X1_PCSAFTQ_ = gl%A3X1_PCSAFTQ(1,:)
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x2_QQ(gl,T,D,getprevious)
    
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)    
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)

    
    !calculate the base function of the second composition derivative of a_qq
        do o = 1, gl%ncomp
			do u = 1, gl%ncomp  
                if (u == o) then
                !    AQQX2(o,u) = (A2_PCSAFTQ(1) * ((A2_PCSAFTQ(1) - 2.d0 * A3_PCSAFTQ(1)) * (A2_PCSAFTQ(1) - A3_PCSAFTQ(1)) * A2X2_PCSAFTQ(1,o,u) + A2_PCSAFTQ(1) &
                !                            & * (A2_PCSAFTQ(1) * A3X2_PCSAFTQ(1,o,u) - A3_PCSAFTQ(1) * A3X2_PCSAFTQ(1,o,u) + 2.d0 * A3X1_PCSAFTQ(1,o)**2)) &
                !                            & - 4.d0 * A2_PCSAFTQ(1) * A3_PCSAFTQ(1) * A2X1_PCSAFTQ(1,o) * A3X1_PCSAFTQ(1,o) + 2.d0 * A3_PCSAFTQ(1)**2 * A2X1_PCSAFTQ(1,o)**2) &
                !                            & /(A2_PCSAFTQ(1) - A3_PCSAFTQ(1))**3 
                !else                                                                       
                !    AQQX2(o,u) =-(2.d0 * A3_PCSAFTQ(1) * A2X1_PCSAFTQ(1,u) * (A2_PCSAFTQ(1) * A3X1_PCSAFTQ(1,o) - A3_PCSAFTQ(1) * A2X1_PCSAFTQ(1,o)) + A2_PCSAFTQ(1) * (- A2_PCSAFTQ(1) &
                !                            & * (2.d0 * A3X1_PCSAFTQ(1,o) * A3X1_PCSAFTQ(1,u) + A3_PCSAFTQ(1) * (-A3X2_PCSAFTQ(1,o,u) - 3.d0 * A2X2_PCSAFTQ(1,o,u))) &
                !                            & + A2_PCSAFTQ(1) * A2_PCSAFTQ(1) * (- (A3X2_PCSAFTQ(1,o,u) + A2X2_PCSAFTQ(1,o,u))) + 2.d0 * A3_PCSAFTQ(1) * (A2X1_PCSAFTQ(1,o) * A3X1_PCSAFTQ(1,u) - A3_PCSAFTQ(1) * A2X2_PCSAFTQ(1,o,u)))) &
                !                            & / (A2_PCSAFTQ(1)-A3_PCSAFTQ(1))**3 
                    AQQX2(o,u) = (gl%A2_PCSAFTQ(1) * ((gl%A2_PCSAFTQ(1) - 2.d0 * gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)) * gl%A2X2_PCSAFTQ(1,o,u) + gl%A2_PCSAFTQ(1) &
                                            & * (gl%A2_PCSAFTQ(1) * gl%A3X2_PCSAFTQ(1,o,u) - gl%A3_PCSAFTQ(1) * gl%A3X2_PCSAFTQ(1,o,u) + 2.d0 * gl%A3X1_PCSAFTQ(1,o)**2)) &
                                            & - 4.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,o) * gl%A3X1_PCSAFTQ(1,o) + 2.d0 * gl%A3_PCSAFTQ(1)**2 * gl%A2X1_PCSAFTQ(1,o)**2) &
                                            & /(gl%A2_PCSAFTQ(1) - gl%A3_PCSAFTQ(1))**3 
                else                                                                       
                    AQQX2(o,u) =-(2.d0 * gl%A3_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,u) * (gl%A2_PCSAFTQ(1) * gl%A3X1_PCSAFTQ(1,o) - gl%A3_PCSAFTQ(1) * gl%A2X1_PCSAFTQ(1,o)) + gl%A2_PCSAFTQ(1) * (- gl%A2_PCSAFTQ(1) &
                                            & * (2.d0 * gl%A3X1_PCSAFTQ(1,o) * gl%A3X1_PCSAFTQ(1,u) + gl%A3_PCSAFTQ(1) * (-gl%A3X2_PCSAFTQ(1,o,u) - 3.d0 * gl%A2X2_PCSAFTQ(1,o,u))) &
                                            & + gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(1) * (- (gl%A3X2_PCSAFTQ(1,o,u) + gl%A2X2_PCSAFTQ(1,o,u))) + 2.d0 * gl%A3_PCSAFTQ(1) * (gl%A2X1_PCSAFTQ(1,o) * gl%A3X1_PCSAFTQ(1,u) - gl%A3_PCSAFTQ(1) * gl%A2X2_PCSAFTQ(1,o,u)))) &
                                            & / (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))**3 
                end if
			end do
		end do

!DEC$ END IF
end subroutine base_AQQX2






    end module pc_saft_ANCX_derivs_module
