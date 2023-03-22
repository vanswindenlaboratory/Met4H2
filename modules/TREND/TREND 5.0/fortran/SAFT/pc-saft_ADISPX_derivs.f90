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

    ! module for file pc-saft_ADISPX_derivs.f90
    module pc_saft_ADISPX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains



subroutine ADISPX1DERIVS(gl,D, GETDERADISP)

! Henning Markgraf, June 2016
    
    ! a_disp: dispersion contribution to the Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2001:
    ! a_disp = -2*pi*rho*I_1*meo1-pi*rho*mmean*C_1*I_2*meo2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERADISP ! input
    !output: adispx1_PCSAFT 
    !working variable
    double precision :: part1, part2, part3, part4, part5, part6, part7, part8, part9, part10, part11, part12
    double precision :: part13, part14, part15, part16
    integer :: i, xi

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_disp
    ! 1: a_disp
    if (GETDERADISP(1) .eq. 1) then
        part1 = -piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*D
        part2 = - piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*D
        part3 = - piPCSAFT*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT*D
        part4 = - piPCSAFT*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT*D
        part5 = - 2.d0*piPCSAFT*gl%i1_PCSAFT(1)*D
        part6 = - 2.d0*piPCSAFT*gl%meo1_PCSAFT(1)*D
        do xi = 1, gl%ncomp
            gl%adispx1_PCSAFT(1,xi) = part1*gl%mPCSAFT(xi) + part2*gl%meo2x1_PCSAFT(1,xi) + part3*gl%i2x1_PCSAFT(1,xi) + part4*gl%cx1_PCSAFT(1,xi) &
                        & + part5*gl%meo1x1_PCSAFT(1,xi) + part6*gl%i1x1_PCSAFT(1,xi)
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2)
    if (GETDERADISP(2) .eq. 1) then
        part1 = -gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part2 = - gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part3 = -gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part4 = - 2.d0*gl%meo1_PCSAFT(1)
        part5 = -gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part6 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)
        part7 = -2.d0*gl%i1_PCSAFT(1) - 2.d0*gl%i1_PCSAFT(2)
        part8 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT
        do xi = 1, gl%ncomp
            gl%adispx1_PCSAFT(2,xi) = piPCSAFT*D*(part1*gl%i2x1_PCSAFT(2,xi) + gl%cx1_PCSAFT(2,xi)*part2 + gl%cx1_PCSAFT(1,xi)*part3 + &
                        & part4*(gl%i1x1_PCSAFT(2,xi) + gl%i1x1_PCSAFT(1,xi)) + gl%i2x1_PCSAFT(1,xi)*part5 + gl%mPCSAFT(xi)*part6 + &
                        & gl%meo1x1_PCSAFT(1,xi)*part7 + gl%meo2x1_PCSAFT(1,xi)*part8)
        end do
    end if

    
    ! 3: 2ND DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires c_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(3), i1_PCSAFT(3), i2_PCSAFT(3)
    if (GETDERADISP(3) .eq. 1) then
        part1 = -gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part2 = gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part3 = (-2.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part4 = (-2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part5 = 2.d0*gl%meo1_PCSAFT(1)
        part6 = 4.d0*gl%meo1_PCSAFT(1)
        part7 = (-2.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part8 = (-2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part9 = (-2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - &
     & gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1))
        part10 = (-4.d0*gl%i1_PCSAFT(2) - 2.d0*gl%i1_PCSAFT(3))
        part11 = (-2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2) &
     & *gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT)
        do xi = 1, gl%ncomp
            gl%adispx1_PCSAFT(3,xi) = piPCSAFT*D*(part1*gl%i2x1_PCSAFT(3,xi) - gl%cx1_PCSAFT(3,xi)*part2 + gl%cx1_PCSAFT(2,xi)*part3 + &
     & gl%cx1_PCSAFT(1,xi)*part4 - part5*gl%i1x1_PCSAFT(3,xi) - part6*gl%i1x1_PCSAFT(2,xi) + gl%i2x1_PCSAFT(2,xi)*part7 + &
     & gl%i2x1_PCSAFT(1,xi)*part8 + gl%mPCSAFT(xi)*part9 + gl%meo1x1_PCSAFT(1,xi)*part10 + &
     & gl%meo2x1_PCSAFT(1,xi)*part11)
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4)
    if (GETDERADISP(4) .eq. 1) then
        part1 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT
        part2 = - gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part3 = - gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part4 = -gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part5 = - 2.d0*gl%i1_PCSAFT(1)
        part6 = - 2.d0*gl%i1_PCSAFT(4)
        part7 = - 2.d0*gl%meo1_PCSAFT(1)
        part8 = - 2.d0*gl%meo1_PCSAFT(4)
        part9 = -gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part10 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)
        part11 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT
        do xi = 1, gl%ncomp
            gl%adispx1_PCSAFT(4,xi) = piPCSAFT*D*(part1*gl%meo2x1_PCSAFT(4,xi) + part2*gl%i2x1_PCSAFT(4,xi) + gl%cx1_PCSAFT(4,xi)*part3 + gl%cx1_PCSAFT(1,xi)*part4 + &
                        & part5*gl%meo1x1_PCSAFT(4,xi) + part6*gl%meo1x1_PCSAFT(1,xi) + part7*gl%i1x1_PCSAFT(4,xi) + part8*gl%i1x1_PCSAFT(1,xi) + &
                        & gl%i2x1_PCSAFT(1,xi)*part9 + gl%mPCSAFT(xi)*part10 + gl%meo2x1_PCSAFT(1,xi)*part11)
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(5), meo2_PCSAFT(5), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5)
    if (GETDERADISP(5) .eq. 1) then
        part1 = -gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT
        part2 =  gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part3 = gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT
        part4 = (-2.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part5 = (-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)* &
     & gl%mmean_PCSAFT - gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part6 = 2.d0*gl%i1_PCSAFT(1)
        part7 = 4.d0*gl%i1_PCSAFT(4)
        part8 = 2.d0*gl%i1_PCSAFT(5)
        part9 = 2.d0*gl%meo1_PCSAFT(1)
        part10 = 4.d0*gl%meo1_PCSAFT(4)
        part11 = 2.d0*gl%meo1_PCSAFT(5)
        part12 = (-2.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT)
        part13 = (-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT)
        part14 = (-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4) - &
     & gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - &
     & gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1))
        part15 = (-2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)* &
     & gl%mmean_PCSAFT)
        part16 = (-gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT)
        do xi = 1, gl%ncomp
            gl%adispx1_PCSAFT(5,xi) = piPCSAFT*D*(part1*gl%meo2x1_PCSAFT(5,xi) -part2*gl%i2x1_PCSAFT(5,xi) - &
     & gl%cx1_PCSAFT(5,xi)*part3 + gl%cx1_PCSAFT(4,xi)*part4 + gl%cx1_PCSAFT(1,xi)*part5 - &
     & part6*gl%meo1x1_PCSAFT(5,xi) - part7*gl%meo1x1_PCSAFT(4,xi) - part8*gl%meo1x1_PCSAFT(1,xi) &
     & - part9*gl%i1x1_PCSAFT(5,xi) - part10*gl%i1x1_PCSAFT(4,xi) - part11*gl%i1x1_PCSAFT(1,xi) + &
     & gl%i2x1_PCSAFT(4,xi)*part12 + gl%i2x1_PCSAFT(1,xi)*part13 + gl%mPCSAFT(xi)*part14 + &
     & gl%meo2x1_PCSAFT(4,xi)*part15 + gl%meo2x1_PCSAFT(1,xi)*part16)
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_disp WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), c_PCSAFT(6), i1_PCSAFT(6), i2_PCSAFT(6)
    if (GETDERADISP(6) .eq. 1) then
        gl%adispx1_PCSAFT(6,:) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(6,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(6,:)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
     & gl%cx1_PCSAFT(4,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(2,:)*( &
     & -gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(1,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)* &
     & gl%mmean_PCSAFT - gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)* &
     & gl%mmean_PCSAFT) - 2.d0*gl%i1x1_PCSAFT(6,:)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(4,:)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(2,:)* &
     & gl%meo1_PCSAFT(4) - 2.d0*gl%i1x1_PCSAFT(1,:)*gl%meo1_PCSAFT(4) + gl%i2x1_PCSAFT(4,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(2,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%i2x1_PCSAFT(1,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(6)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%mPCSAFT(:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)* &
     & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - &
     & gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - &
     & gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)) + gl%meo1x1_PCSAFT(4,:)*(-2.d0*gl%i1_PCSAFT(1) - 2.d0*gl%i1_PCSAFT(2)) + gl%meo1x1_PCSAFT(1,:)*( &
     & -2.d0*gl%i1_PCSAFT(4) - 2.d0*gl%i1_PCSAFT(6)) + gl%meo2x1_PCSAFT(4,:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)* &
     & gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%meo2x1_PCSAFT(1,:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(6)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(6) &
     & *gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT)) !note: does it work as fast this way?
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_disp WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(4), meo2_PCSAFT(4),
    !          c_PCSAFT(6), i1_PCSAFT(6), i2_PCSAFT(6), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5), meo1_PCSAFT(5), meo2_PCSAFT(5),
    !          i1_PCSAFT(7), i2_PCSAFT(7), c_PCSAFT(7)
    ! c_PCSAFT(1,2,4,5,6), I(1,2,4,5,6,7), meo(1,4,5)
    if (GETDERADISP(7) .eq. 1) then
        gl%adispx1_PCSAFT(7,:) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(7,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(6,:)*gl%meo2_PCSAFT(4)* &
     & gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(6,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(7,:)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT + gl%cx1_PCSAFT(5,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) - 2.d0* &
     & gl%cx1_PCSAFT(6,:)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(6,:)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
     & gl%cx1_PCSAFT(4,:)*(-2.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(6)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(2,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)* &
     & gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(1,:)*(-gl%i2_PCSAFT(1) &
     & *gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
     & gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(5)* &
     & gl%mmean_PCSAFT) - 2.d0*gl%i1x1_PCSAFT(7,:)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(5,:)*gl%meo1_PCSAFT(1) - 4.d0*gl%i1x1_PCSAFT(6,:)* &
     & gl%meo1_PCSAFT(4) - 4.d0*gl%i1x1_PCSAFT(4,:)*gl%meo1_PCSAFT(4) - 2.d0*gl%i1x1_PCSAFT(2,:)*gl%meo1_PCSAFT(5) - 2.d0*gl%i1x1_PCSAFT(1,:)* &
     & gl%meo1_PCSAFT(5) + gl%i2x1_PCSAFT(5,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%i2x1_PCSAFT(4,:)*(-2.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(2,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(5)* &
     & gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(1,:)*(-gl%c_PCSAFT(1) &
     & *gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%meo2_PCSAFT(5)* &
     & gl%mmean_PCSAFT) + gl%mPCSAFT(:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)* &
     & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)* &
     & gl%meo2_PCSAFT(5) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(4)* &
     & gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(5) &
     & *gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(7)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - 2.d0* &
     & gl%c_PCSAFT(6)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)* &
     & gl%meo2_PCSAFT(4) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)) + gl%meo1x1_PCSAFT(5,:)*(-2.d0*gl%i1_PCSAFT(1) - 2.d0*gl%i1_PCSAFT(2)) + &
     & gl%meo1x1_PCSAFT(4,:)*(-4.d0*gl%i1_PCSAFT(4) - 4.d0*gl%i1_PCSAFT(6)) + gl%meo1x1_PCSAFT(1,:)*(-2.d0*gl%i1_PCSAFT(5) - 2.d0* &
     & gl%i1_PCSAFT(7)) + gl%meo2x1_PCSAFT(5,:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%meo2x1_PCSAFT(4,:)*(-2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(6)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - 2.d0* &
     & gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%meo2x1_PCSAFT(1,:)*(-gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(5)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(7)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0* &
     & gl%c_PCSAFT(4)*gl%i2_PCSAFT(6)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(7)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(5)* &
     & gl%mmean_PCSAFT))
    end if
    
    ! 8: 3RD DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires c_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(3), i1_PCSAFT(3), i2_PCSAFT(3), i1_PCSAFT(8), i2_PCSAFT(8), c_PCSAFT(8)
    if (GETDERADISP(8) .eq. 1) then
        gl%adispx1_PCSAFT(8,:) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(8,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(8,:)*gl%i2_PCSAFT(1)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%cx1_PCSAFT(3,:)*(-3.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 3.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT) + gl%cx1_PCSAFT(2,:)*(-6.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 3.d0*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT) + gl%cx1_PCSAFT(1,:)*(-3.d0*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT) - 2.d0*gl%i1x1_PCSAFT(8,:)*gl%meo1_PCSAFT(1) - 6.d0*gl%i1x1_PCSAFT(3,:)*gl%meo1_PCSAFT(1) + &
     & gl%i2x1_PCSAFT(3,:)*(-3.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%i2x1_PCSAFT(2,:)*(-6.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%i2x1_PCSAFT(1,:)*(-3.d0*gl%c_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(8)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%mPCSAFT(:)*(-3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(1) - 6.d0*gl%c_PCSAFT(2)* &
     & gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1) - 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) &
     & - 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(8)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)) + gl%meo1x1_PCSAFT(1,:)* &
     & (-6.d0*gl%i1_PCSAFT(3) - 2.d0*gl%i1_PCSAFT(8)) + gl%meo2x1_PCSAFT(1,:)*(-3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)* &
     & gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(8)*gl%mmean_PCSAFT - 6.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - 3.d0* &
     & gl%c_PCSAFT(2)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(3)* &
     & gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(8)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT))
    end if
    
    ! 9: 3RD DERIVATIVE OF a_disp WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(5), meo2_PCSAFT(5), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5)
    !          c_PCSAFT(9), i1_PCSAFT(9), i2_PCSAFT(9), meo1_PCSAFT(9), meo2_PCSAFT(9)
    ! c_PCSAFT(1,4,5,9), I(1,4,5,9), meo(1,4,5,9)
    if (GETDERADISP(9) .eq. 1) then
        gl%adispx1_PCSAFT(9,:) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(9,:)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(9,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
     & gl%cx1_PCSAFT(9,:)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%cx1_PCSAFT(5,:)*(-3.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0* &
     & gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(4,:)*(-3.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 6.d0*gl%i2_PCSAFT(4)* &
     & gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(1,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9)* &
     & gl%mmean_PCSAFT - 3.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 3.d0*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(9)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) - 2.d0*gl%i1_PCSAFT(1)*gl%meo1x1_PCSAFT(9,:) - 6.d0*gl%i1_PCSAFT(4)*gl%meo1x1_PCSAFT(5,:) - 6.d0*gl%i1_PCSAFT(5)* &
     & gl%meo1x1_PCSAFT(4,:) - 2.d0*gl%i1_PCSAFT(9)*gl%meo1x1_PCSAFT(1,:) - 2.d0*gl%i1x1_PCSAFT(9,:)*gl%meo1_PCSAFT(1) - 6.d0*gl%i1x1_PCSAFT(5,:)* &
     & gl%meo1_PCSAFT(4) - 6.d0*gl%i1x1_PCSAFT(4,:)*gl%meo1_PCSAFT(5) - 2.d0*gl%i1x1_PCSAFT(1,:)*gl%meo1_PCSAFT(9) + gl%i2x1_PCSAFT(5,:)*(-3.d0* &
     & gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(4,:)*(-3.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(5) &
     & *gl%mmean_PCSAFT - 6.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(1,:)*( &
     & -gl%c_PCSAFT(1)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(4)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(5)*gl%meo2_PCSAFT(4)* &
     & gl%mmean_PCSAFT - gl%c_PCSAFT(9)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%mPCSAFT(:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9) - 3.d0*gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(5) - 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(9)*gl%meo2_PCSAFT(1) - 3.d0*gl%c_PCSAFT(4)* &
     & gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5) - 6.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4) - 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1) - 3.d0* &
     & gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(9)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)) + &
     & gl%meo2x1_PCSAFT(5,:)*(-3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%meo2x1_PCSAFT(4,:)*( &
     & -3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%mmean_PCSAFT - 6.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%meo2x1_PCSAFT(1,:)*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(9)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(5)*gl%mmean_PCSAFT - 3.d0*gl%c_PCSAFT(5)* &
     & gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(9)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT))
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF a_disp WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! c_PCSAFT(1,2,3,4,6,10), I(1,2,3,4,6,10), meo(1,4)
    if (GETDERADISP(10) .eq. 1) then
        gl%adispx1_PCSAFT(10,:) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(10,:)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(10,:)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT + gl%cx1_PCSAFT(6,:)*(-2.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%cx1_PCSAFT(4,:)*(-2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + &
     & gl%cx1_PCSAFT(3,:)*(-gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(2,:) &
     & *(-2.d0*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%cx1_PCSAFT(1,:)*(-2.d0*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
     & - gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%i2_PCSAFT(3)* &
     & gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) - 2.d0*gl%i1x1_PCSAFT(10,:)*gl%meo1_PCSAFT(1) - 4.d0*gl%i1x1_PCSAFT(6,:)*gl%meo1_PCSAFT(1) - 2.d0* &
     & gl%i1x1_PCSAFT(3,:)*gl%meo1_PCSAFT(4) - 4.d0*gl%i1x1_PCSAFT(2,:)*gl%meo1_PCSAFT(4) + gl%i2x1_PCSAFT(6,:)*(-2.d0*gl%c_PCSAFT(1)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(4,:)*(-2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(1)* &
     & gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(3,:)*(-gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT &
     & - gl%c_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%i2x1_PCSAFT(2,:)*(-2.d0*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)* &
     & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + &
     & gl%i2x1_PCSAFT(1,:)*(-2.d0*gl%c_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0* &
     & gl%c_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT) + gl%mPCSAFT(:)*(-2.d0*gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1) &
     & *gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1) &
     & - 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(10)*gl%i2_PCSAFT(1) &
     & *gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)* &
     & gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4) &
     & - gl%c_PCSAFT(3)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)) + gl%meo1x1_PCSAFT(4,:)*(-4.d0*gl%i1_PCSAFT(2) - 2.d0*gl%i1_PCSAFT(3)) &
     & + gl%meo1x1_PCSAFT(1,:)*(-4.d0*gl%i1_PCSAFT(6) - 2.d0*gl%i1_PCSAFT(10)) + gl%meo2x1_PCSAFT(4,:)*(-2.d0*gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0* &
     & gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT) + gl%meo2x1_PCSAFT(1,:)*(-2.d0*gl%c_PCSAFT(1)* &
     & gl%i2_PCSAFT(6)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(10)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - &
     & gl%c_PCSAFT(4)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(6)* &
     & gl%i2_PCSAFT(2)*gl%mmean_PCSAFT - gl%c_PCSAFT(10)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0 &
     & *gl%c_PCSAFT(2)*gl%i2_PCSAFT(6)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT))
    end if
    !DEC$ END IF
end subroutine ADISPX1DERIVS

    

subroutine ADISPX2DERIVS(gl,D, GETDERADISP)
    
! Henning Markgraf, June 2016
    
    ! a_disp: dispersion contribution to the Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2001:
    ! a_disp = -2*pi*rho*I_1*meo1-pi*rho*mmean*C_1*I_2*meo2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERADISP
    !output: adispx2_PCSAFT (module variable)
    !working variable
    double precision :: part1, part2, part3, part4, part5, part6, part7, part8, part9, part10, part11, part12
    integer :: i, xi, xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Calculate the derivatives of a_disp
    ! 1: d(a_disp)/(dxi dxj)
    ! note: possible to put the non-x-dependent parts before the loops and save computation time
    if (GETDERADISP(1) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
              !  if (xi .GE. xj) then
         gl%adispx2_PCSAFT(1,xi,xj) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) + gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) + gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)* &
     & gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) + gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)* &
     & gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) + gl%c_PCSAFT(1)* &
     & gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT + gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) + gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)* &
     & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT + gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) &
     & *gl%mmean_PCSAFT + gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) + gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT + &
     & gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%i1_PCSAFT(1)*gl%meo1x2_PCSAFT(1,xi,xj) + 2.d0*gl%i1x1_PCSAFT(1,xj)* &
     & gl%meo1x1_PCSAFT(1,xi) + 2.d0*gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x1_PCSAFT(1,xj))
             !   end if
            end do
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2)
    if (GETDERADISP(2) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%adispx2_PCSAFT(2,xi,xj) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)* &
         & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2_PCSAFT(1)* &
         & gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
         & gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)* &
         & gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)* &
         & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)* &
         & gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(2,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)* &
         & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)* &
         & gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(2,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0* &
         & gl%i1_PCSAFT(1)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1_PCSAFT(2)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1x1_PCSAFT(2,xj)*gl%meo1x1_PCSAFT(1,xi) - &
         & 2.d0*gl%i1x2_PCSAFT(2,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(2,xi)*gl%meo1x1_PCSAFT(1,xj) - 2.d0*gl%i1x1_PCSAFT(1,xj)*gl%meo1x1_PCSAFT(1,xi) &
         & - 2.d0*gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x1_PCSAFT(1,xj))
                end if
            end do
        end do
    end if

    
    ! 3: 2ND DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires c_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(3), i1_PCSAFT(3), i2_PCSAFT(3)
    if (GETDERADISP(3) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%adispx2_PCSAFT(3,xi,xj) = piPCSAFT*D*(-2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - &
         & 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(3,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(3,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(3,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(3,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - &
         & gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(3,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0* &
         & gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0* &
         & gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - 2.d0*gl%c_PCSAFT(2) &
         & *gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - 2.d0* &
         & gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xi,xj)* &
         & gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(2,xj)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0* &
         & gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0* &
         & gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0* &
         & gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(3)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(3)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(3)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(3)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - &
         & gl%c_PCSAFT(3)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(3,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - &
         & gl%cx1_PCSAFT(3,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(3,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)* &
         & gl%mmean_PCSAFT - gl%cx2_PCSAFT(3,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(3,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)* &
         & gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(3,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(3,xi)*gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(1)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(2,xj)* &
         & gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx2_PCSAFT(2,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & 2.d0*gl%cx2_PCSAFT(2,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - &
         & 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - &
         & 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(2,xj)* &
         & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)* &
         & gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(3)* &
         & gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(3)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2x1_PCSAFT(3,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
         & - 2.d0*gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(3)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(3)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(3,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(2,xj)*gl%meo2_PCSAFT(1)* &
         & gl%mmean_PCSAFT - 4.d0*gl%i1_PCSAFT(2)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1_PCSAFT(3)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0* &
         & gl%i1x1_PCSAFT(3,xj)*gl%meo1x1_PCSAFT(1,xi) - 2.d0*gl%i1x2_PCSAFT(3,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0* &
         & gl%i1x1_PCSAFT(3,xi)*gl%meo1x1_PCSAFT(1,xj) - 4.d0*gl%i1x1_PCSAFT(2,xj)*gl%meo1x1_PCSAFT(1,xi) - 4.d0*gl%i1x2_PCSAFT(2,xi,xj)* &
         & gl%meo1_PCSAFT(1) - 4.d0*gl%i1x1_PCSAFT(2,xi)*gl%meo1x1_PCSAFT(1,xj))
                end if
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4)
    if (GETDERADISP(4) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%adispx2_PCSAFT(4,xi,xj) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)* &
         & gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - &
         & gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(4,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
         & gl%i2x1_PCSAFT(4,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi) &
         & *gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(4)* &
         & gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)* &
         & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - &
         & gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(4,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(1) &
         & *gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)* &
         & gl%meo2_PCSAFT(4) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(4,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(4)* &
         & gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(4,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%i1_PCSAFT(1)*gl%meo1x2_PCSAFT(4,xi,xj) - 2.d0*gl%i1_PCSAFT(4)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0* &
         & gl%i1x1_PCSAFT(4,xj)*gl%meo1x1_PCSAFT(1,xi) - 2.d0*gl%i1x2_PCSAFT(4,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(4,xi)*gl%meo1x1_PCSAFT(1,xj) - 2.d0* &
         & gl%i1x1_PCSAFT(1,xj)*gl%meo1x1_PCSAFT(4,xi) - 2.d0*gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1_PCSAFT(4) - 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x1_PCSAFT(4,xj))
                end if
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(5), meo2_PCSAFT(5), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5)
    if (GETDERADISP(5) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%adispx2_PCSAFT(5,xi,xj) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(5,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(5,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1) &
         & *gl%meo2x2_PCSAFT(5,xi,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)* &
         & gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%mPCSAFT(xj)* &
         & gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT &
         & - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(5,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(5,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
         & gl%i2x2_PCSAFT(5,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(5,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(5,xi)* &
         & gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj)* &
         & gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(4,xi,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xi) &
         & *gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)* &
         & gl%meo2_PCSAFT(5) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(5,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(5)* &
         & gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(5) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(5,xj)*gl%mmean_PCSAFT - &
         & 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2x1_PCSAFT(4,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(4,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2x2_PCSAFT(4,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(4,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2x1_PCSAFT(4,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(4)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - 2.d0*gl%c_PCSAFT(4)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0* &
         & gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - 2.d0*gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - &
         & gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(5)*gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(5)*gl%i2x1_PCSAFT(1,xi)* &
         & gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(5)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(5,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) &
         & - gl%cx1_PCSAFT(5,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(5,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(5,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(5,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(5,xi)* &
         & gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(5,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(4,xj)* &
         & gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - 2.d0*gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(4,xj)* &
         & gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(4,xj)* &
         & gl%i2x1_PCSAFT(4,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(4,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx2_PCSAFT(4,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%cx2_PCSAFT(4,xi,xj)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0 &
         & *gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - 2.d0*gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - 2.d0*gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx1_PCSAFT(4,xi)*gl%i2x1_PCSAFT(4,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(4,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(5) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(5,xi)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - 2.d0*gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(5)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(5)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj) &
         & *gl%i2x1_PCSAFT(5,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(4,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj) &
         & *gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(5) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(5,xj)*gl%mmean_PCSAFT - 2.d0* &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - 2.d0*gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(5)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(5)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi) &
         & *gl%i2x1_PCSAFT(5,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - 2.d0*gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(4,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi) &
         & *gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT - 2.d0*gl%i1_PCSAFT(1)*gl%meo1x2_PCSAFT(5,xi,xj) - 4.d0*gl%i1_PCSAFT(4)* &
         & gl%meo1x2_PCSAFT(4,xi,xj) - 2.d0*gl%i1_PCSAFT(5)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1x1_PCSAFT(5,xj)*gl%meo1x1_PCSAFT(1,xi) - 2.d0* &
         & gl%i1x2_PCSAFT(5,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(5,xi)*gl%meo1x1_PCSAFT(1,xj) - 4.d0*gl%i1x1_PCSAFT(4,xj)*gl%meo1x1_PCSAFT(4,xi) - &
         & 4.d0*gl%i1x2_PCSAFT(4,xi,xj)*gl%meo1_PCSAFT(4) - 4.d0*gl%i1x1_PCSAFT(4,xi)*gl%meo1x1_PCSAFT(4,xj) - 2.d0*gl%i1x1_PCSAFT(1,xj)* &
         & gl%meo1x1_PCSAFT(5,xi) - 2.d0*gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1_PCSAFT(5) - 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x1_PCSAFT(5,xj))
                end if
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_disp WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires c_PCSAFT(1), i1_PCSAFT(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), c_PCSAFT(6), i1_PCSAFT(6), i2_PCSAFT(6)
    if (GETDERADISP(6) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%adispx2_PCSAFT(6,xi,xj) = piPCSAFT*D*(-gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)* &
         & gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%mPCSAFT(xj)* &
         & gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%meo2x2_PCSAFT(1,xi,xj)* &
         & gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(6,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(6,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT &
         & - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(6,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(6,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)* &
         & gl%i2x1_PCSAFT(6,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xj) &
         & *gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(4,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xi)*gl%mPCSAFT(xj)* &
         & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(4,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - &
         & gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
         & gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xj)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)* &
         & gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)* &
         & gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xj) &
         & *gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)* &
         & gl%meo2_PCSAFT(4) - gl%c_PCSAFT(1)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)* &
         & gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(4)*gl%i2_PCSAFT(2) &
         & *gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(2,xj) &
         & *gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)* &
         & gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)* &
         & gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(4)*gl%i2x1_PCSAFT(1,xi)* &
         & gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(6)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)* &
         & gl%meo2_PCSAFT(1) - gl%c_PCSAFT(6)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(6)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(1) &
         & *gl%mmean_PCSAFT - gl%c_PCSAFT(6)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(6)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj)* &
         & gl%mmean_PCSAFT - gl%cx1_PCSAFT(6,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(6,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(6,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(6,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(6,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(6,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(6,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(4,xj) &
         & *gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(4,xj)* &
         & gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(4,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(4,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(4,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(4,xi) &
         & *gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(4,xi)* &
         & gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(4,xi)*gl%i2x1_PCSAFT(2,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(4,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(4,xi) - gl%c_PCSAFT(2) &
         & *gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(4,xj) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2x2_PCSAFT(4,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)* &
         & gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(4,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(4,xj)* &
         & gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(4,xi,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(4,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(4,xi)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)* &
         & gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - gl%cx1_PCSAFT(2,xj)* &
         & gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xj)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(2,xj)* &
         & gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(4,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(2,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT &
         & - gl%cx2_PCSAFT(2,xi,xj)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - &
         & gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - &
         & gl%cx1_PCSAFT(2,xi)*gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(4,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(6)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2_PCSAFT(6)*gl%meo2x1_PCSAFT(1,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(6,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(4,xi)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2_PCSAFT(4) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(4,xi)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(1)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(4)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2_PCSAFT(4)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(6)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(1) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2_PCSAFT(6)*gl%meo2x1_PCSAFT(1,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(6,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(4,xj)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2_PCSAFT(4) - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(4,xj)*gl%mmean_PCSAFT - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(2,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT &
         & - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT - 2.d0*gl%i1_PCSAFT(1)*gl%meo1x2_PCSAFT(4,xi,xj) - 2.d0*gl%i1_PCSAFT(4)* &
         & gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1_PCSAFT(6)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1x1_PCSAFT(6,xj)*gl%meo1x1_PCSAFT(1,xi) - 2.d0 &
         & *gl%i1x2_PCSAFT(6,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(6,xi)*gl%meo1x1_PCSAFT(1,xj) - 2.d0*gl%i1x1_PCSAFT(4,xj)* &
         & gl%meo1x1_PCSAFT(1,xi) - 2.d0*gl%i1x2_PCSAFT(4,xi,xj)*gl%meo1_PCSAFT(1) - 2.d0*gl%i1x1_PCSAFT(4,xi)*gl%meo1x1_PCSAFT(1,xj) - 2.d0*gl%i1_PCSAFT(2)* &
         & gl%meo1x2_PCSAFT(4,xi,xj) - 2.d0*gl%i1x1_PCSAFT(2,xj)*gl%meo1x1_PCSAFT(4,xi) - 2.d0*gl%i1x2_PCSAFT(2,xi,xj)*gl%meo1_PCSAFT(4) - &
         & 2.d0*gl%i1x1_PCSAFT(2,xi)*gl%meo1x1_PCSAFT(4,xj) - 2.d0*gl%i1x1_PCSAFT(1,xj)*gl%meo1x1_PCSAFT(4,xi) - 2.d0*gl%i1x2_PCSAFT(1,xi,xj)* &
         & gl%meo1_PCSAFT(4) - 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x1_PCSAFT(4,xj))
                end if
            end do
        end do
    end if

    !DEC$ END IF
end subroutine ADISPX2DERIVS

    

subroutine ADISPX3DERIVS(gl,D, GETDERADISP)
    
! Henning Markgraf, June 2016

    ! a_disp: dispersion contribution to the Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2001:
    ! a_disp = -2*pi*rho*I_1*meo1-pi*rho*mmean*C_1*I_2*meo2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERADISP
    !output: adispx3_PCSAFT
    !working variable
    double precision :: part1, part2, part3, part4, part5, part6, part7, part8, part9, part10, part11, part12
    integer :: i, xi, xj, xk
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X3 DERIVATIVES of a_disp
    ! note: still possible to separate non-x-dependent parts and calculate them before the loops to save calculation time?
    !calculate the derivatives of a_disp
    ! 1: a_disp
    if (GETDERADISP(1) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then 
            	        gl%adispx3_PCSAFT(1,xi,xj,xk) = -piPCSAFT*D*(gl%c_PCSAFT(1)*(gl%i2_PCSAFT(1)*(gl%mPCSAFT(xj)*gl%meo2x2_PCSAFT(1,xi,xk) + gl%mPCSAFT(xi)*gl%meo2x2_PCSAFT(1,xj,xk) + gl%mPCSAFT(xk)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)) + gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) + gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) + &
         & gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) + gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) + gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj)* &
         & gl%meo2x1_PCSAFT(1,xi) + gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) + gl%meo2_PCSAFT(1)*(gl%i2x2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) + &
         & gl%i2x2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) + gl%i2x2_PCSAFT(1,xi,xk)*gl%mPCSAFT(xj))) + 2.d0*gl%i1x1_PCSAFT(1,xj)*gl%meo1x2_PCSAFT(1,xi,xk) + 2.d0* &
         & gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1x1_PCSAFT(1,xk) + 2.d0*gl%i1x3_PCSAFT(1,xi,xj,xk)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1x2_PCSAFT(1,xj,xk)* &
         & gl%meo1x1_PCSAFT(1,xi) + 2.d0*gl%i1x1_PCSAFT(1,xi)*gl%meo1x2_PCSAFT(1,xj,xk) + 2.d0*gl%i1x2_PCSAFT(1,xi,xk)*gl%meo1x1_PCSAFT(1,xj) + 2.d0* &
         & gl%i1x1_PCSAFT(1,xk)*gl%meo1x2_PCSAFT(1,xi,xj) + gl%i2_PCSAFT(1)*(gl%cx1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) + gl%cx1_PCSAFT(1,xj)*gl%mPCSAFT(xk)* &
         & gl%meo2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) + gl%cx1_PCSAFT(1,xi)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) + gl%cx1_PCSAFT(1,xk)* &
         & gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) + gl%meo2_PCSAFT(1)*(gl%cx2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) + &
         & gl%cx2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) + gl%cx2_PCSAFT(1,xi,xk)*gl%mPCSAFT(xj))) + gl%meo2_PCSAFT(1)*(gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xk) + &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xi) + gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk) + gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj) + &
         & gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi) + gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)) + gl%mmean_PCSAFT*(gl%c_PCSAFT(1)*(gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2x2_PCSAFT(1,xi,xk) + gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(1,xk) + gl%i2x2_PCSAFT(1,xj,xk)*gl%meo2x1_PCSAFT(1,xi) + gl%i2x1_PCSAFT(1,xi) &
         & *gl%meo2x2_PCSAFT(1,xj,xk) + gl%i2x2_PCSAFT(1,xi,xk)*gl%meo2x1_PCSAFT(1,xj) + gl%i2x1_PCSAFT(1,xk)*gl%meo2x2_PCSAFT(1,xi,xj)) + &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xk) + gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xi)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xk) + gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xj) + gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xj)* &
         & gl%meo2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj) + gl%i2_PCSAFT(1)*(gl%cx1_PCSAFT(1,xj)*gl%meo2x2_PCSAFT(1,xi,xk) + &
         & gl%cx2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(1,xk) + gl%cx2_PCSAFT(1,xj,xk)*gl%meo2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xi)*gl%meo2x2_PCSAFT(1,xj,xk) + &
         & gl%cx2_PCSAFT(1,xi,xk)*gl%meo2x1_PCSAFT(1,xj) + gl%cx1_PCSAFT(1,xk)*gl%meo2x2_PCSAFT(1,xi,xj)) + gl%meo2_PCSAFT(1)*(gl%c_PCSAFT(1)* &
         & gl%i2x3_PCSAFT(1,xi,xj,xk) + gl%cx1_PCSAFT(1,xj)*gl%i2x2_PCSAFT(1,xi,xk) + gl%cx2_PCSAFT(1,xi,xj)*gl%i2x1_PCSAFT(1,xk) + &
         & gl%cx3_PCSAFT(1,xi,xj,xk)*gl%i2_PCSAFT(1) + gl%cx2_PCSAFT(1,xj,xk)*gl%i2x1_PCSAFT(1,xi) + gl%cx1_PCSAFT(1,xi)*gl%i2x2_PCSAFT(1,xj,xk) + &
         & gl%cx2_PCSAFT(1,xi,xk)*gl%i2x1_PCSAFT(1,xj) + gl%cx1_PCSAFT(1,xk)*gl%i2x2_PCSAFT(1,xi,xj))))
                    end if
        	    end do
            end do 
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires c_PCSAFT(1), I_1(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), c_PCSAFT(2), I_1(2), i2_PCSAFT(2)
    if (GETDERADISP(2) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then 
            	        gl%adispx3_PCSAFT(2,xi,xj,xk) = piPCSAFT*D*(gl%c_PCSAFT(1)*(gl%i2_PCSAFT(1)*(-gl%mPCSAFT(xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%mPCSAFT(xi)*gl%meo2x2_PCSAFT(1,xj,xk) - gl%mPCSAFT(xk)* &
         & gl%meo2x2_PCSAFT(1,xi,xj)) - gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x2_PCSAFT(1,xj,xk) &
         & - gl%i2_PCSAFT(2)*gl%mPCSAFT(xk)*gl%meo2x2_PCSAFT(1,xi,xj) - gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) - gl%i2x1_PCSAFT(2,xj)* &
         & gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) - gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) &
         & - gl%i2x1_PCSAFT(2,xk)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%i2x1_PCSAFT(2,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi) &
         & *gl%meo2x1_PCSAFT(1,xk) - gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) - gl%i2x1_PCSAFT(1,xi) &
         & *gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) - gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) + &
         & gl%meo2_PCSAFT(1)*(-gl%i2x2_PCSAFT(2,xi,xj)*gl%mPCSAFT(xk) - gl%i2x2_PCSAFT(2,xj,xk)*gl%mPCSAFT(xi) - gl%i2x2_PCSAFT(2,xi,xk)*gl%mPCSAFT(xj) &
         & - gl%i2x2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) - gl%i2x2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) - gl%i2x2_PCSAFT(1,xi,xk)*gl%mPCSAFT(xj))) - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(2)* &
         & gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)* &
         & gl%i2_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xk)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) - 2.d0* &
         & gl%i1x1_PCSAFT(2,xj)*gl%meo1x2_PCSAFT(1,xi,xk) - 2.d0*gl%i1x2_PCSAFT(2,xi,xj)*gl%meo1x1_PCSAFT(1,xk) - 2.d0* &
         & gl%i1x2_PCSAFT(2,xj,xk)*gl%meo1x1_PCSAFT(1,xi) - 2.d0*gl%i1x1_PCSAFT(2,xi)*gl%meo1x2_PCSAFT(1,xj,xk) - 2.d0* &
         & gl%i1x2_PCSAFT(2,xi,xk)*gl%meo1x1_PCSAFT(1,xj) - 2.d0*gl%i1x1_PCSAFT(2,xk)*gl%meo1x2_PCSAFT(1,xi,xj) - 2.d0*gl%i1x1_PCSAFT(1,xj)* &
         & gl%meo1x2_PCSAFT(1,xi,xk) - 2.d0*gl%i1x2_PCSAFT(1,xi,xj)*gl%meo1x1_PCSAFT(1,xk) - 2.d0*gl%i1x2_PCSAFT(1,xj,xk)*gl%meo1x1_PCSAFT(1,xi) - 2.d0* &
         & gl%i1x1_PCSAFT(1,xi)*gl%meo1x2_PCSAFT(1,xj,xk) - 2.d0*gl%i1x2_PCSAFT(1,xi,xk)*gl%meo1x1_PCSAFT(1,xj) - 2.d0*gl%i1x1_PCSAFT(1,xk)* &
         & gl%meo1x2_PCSAFT(1,xi,xj) + gl%i2_PCSAFT(1)*(-gl%c_PCSAFT(2)*gl%mPCSAFT(xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%c_PCSAFT(2)*gl%mPCSAFT(xi)* &
         & gl%meo2x2_PCSAFT(1,xj,xk) - gl%c_PCSAFT(2)*gl%mPCSAFT(xk)*gl%meo2x2_PCSAFT(1,xi,xj) - gl%cx1_PCSAFT(2,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) - &
         & gl%cx1_PCSAFT(2,xj)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(2,xi)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(2,xi)* &
         & gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(2,xk)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(2,xk)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) &
         & - gl%cx1_PCSAFT(1,xj)*gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xj)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xi)*gl%mPCSAFT(xj)* &
         & gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xi)*gl%mPCSAFT(xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)*gl%mPCSAFT(xj)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xk)* &
         & gl%mPCSAFT(xi)*gl%meo2x1_PCSAFT(1,xj) + gl%meo2_PCSAFT(1)*(-gl%cx2_PCSAFT(2,xi,xj)*gl%mPCSAFT(xk) - gl%cx2_PCSAFT(2,xj,xk)*gl%mPCSAFT(xi) - &
         & gl%cx2_PCSAFT(2,xi,xk)*gl%mPCSAFT(xj) - gl%cx2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) - gl%cx2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) - gl%cx2_PCSAFT(1,xi,xk)* &
         & gl%mPCSAFT(xj))) + gl%meo1_PCSAFT(1)*(-2.d0*gl%i1x3_PCSAFT(2,xi,xj,xk) - 2.d0*gl%i1x3_PCSAFT(1,xi,xj,xk)) + gl%meo2_PCSAFT(1)*( &
         & -gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xi,xj)*gl%mPCSAFT(xk) - gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xj,xk)*gl%mPCSAFT(xi) - gl%c_PCSAFT(2)* &
         & gl%i2x2_PCSAFT(1,xi,xk)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xk) - gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xi) &
         & - gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk) - gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(2,xk)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi) - gl%cx1_PCSAFT(2,xk)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xk) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xk)*gl%mPCSAFT(xi) - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xk) - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xk)* &
         & gl%mPCSAFT(xi) - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xk) - gl%cx2_PCSAFT(1,xj,xk)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xi) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2x1_PCSAFT(2,xj)*gl%mPCSAFT(xk) - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(2,xk)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xk) - &
         & gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xk)*gl%mPCSAFT(xj) - gl%cx2_PCSAFT(1,xi,xk)*gl%i2_PCSAFT(2)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(2,xj) &
         & *gl%mPCSAFT(xi) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(2,xi)*gl%mPCSAFT(xj) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xj)*gl%mPCSAFT(xi) - gl%cx1_PCSAFT(1,xk)* &
         & gl%i2x1_PCSAFT(1,xi)*gl%mPCSAFT(xj)) + gl%mmean_PCSAFT*(gl%c_PCSAFT(1)*(-gl%i2x1_PCSAFT(2,xj)*gl%meo2x2_PCSAFT(1,xi,xk) - &
         & gl%i2x2_PCSAFT(2,xi,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%i2x2_PCSAFT(2,xj,xk)*gl%meo2x1_PCSAFT(1,xi) - gl%i2x1_PCSAFT(2,xi)* &
         & gl%meo2x2_PCSAFT(1,xj,xk) - gl%i2x2_PCSAFT(2,xi,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%i2x1_PCSAFT(2,xk)*gl%meo2x2_PCSAFT(1,xi,xj) - &
         & gl%i2x1_PCSAFT(1,xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%i2x2_PCSAFT(1,xj,xk)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%i2x1_PCSAFT(1,xi)*gl%meo2x2_PCSAFT(1,xj,xk) - gl%i2x2_PCSAFT(1,xi,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%i2x1_PCSAFT(1,xk)*gl%meo2x2_PCSAFT(1,xi,xj)) &
         & - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(1,xk) - &
         & gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xj,xk)*gl%meo2x1_PCSAFT(1,xi) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x2_PCSAFT(1,xj,xk) - &
         & gl%c_PCSAFT(2)*gl%i2x2_PCSAFT(1,xi,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%c_PCSAFT(2)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x2_PCSAFT(1,xi,xj) - &
         & gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(2,xj)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(2,xi)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xj) - &
         & gl%cx1_PCSAFT(2,xk)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(2,xk)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xk) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(2,xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xj)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xj)* &
         & gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx2_PCSAFT(1,xi,xj)*gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xk) - gl%cx2_PCSAFT(1,xj,xk)* &
         & gl%i2_PCSAFT(2)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xi)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xj,xk) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2x1_PCSAFT(2,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(2,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xi)* &
         & gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx1_PCSAFT(1,xi)*gl%i2x1_PCSAFT(1,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx2_PCSAFT(1,xi,xk)*gl%i2_PCSAFT(2)* &
         & gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)*gl%i2_PCSAFT(2)*gl%meo2x2_PCSAFT(1,xi,xj) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(2,xj)* &
         & gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(2,xi)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xj)*gl%meo2x1_PCSAFT(1,xi) &
         & - gl%cx1_PCSAFT(1,xk)*gl%i2x1_PCSAFT(1,xi)*gl%meo2x1_PCSAFT(1,xj) + gl%i2_PCSAFT(1)*(-gl%cx1_PCSAFT(2,xj)*gl%meo2x2_PCSAFT(1,xi,xk) - &
         & gl%cx2_PCSAFT(2,xi,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx2_PCSAFT(2,xj,xk)*gl%meo2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(2,xi)* &
         & gl%meo2x2_PCSAFT(1,xj,xk) - gl%cx2_PCSAFT(2,xi,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(2,xk)*gl%meo2x2_PCSAFT(1,xi,xj) - &
         & gl%cx1_PCSAFT(1,xj)*gl%meo2x2_PCSAFT(1,xi,xk) - gl%cx2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(1,xk) - gl%cx2_PCSAFT(1,xj,xk)*gl%meo2x1_PCSAFT(1,xi) - &
         & gl%cx1_PCSAFT(1,xi)*gl%meo2x2_PCSAFT(1,xj,xk) - gl%cx2_PCSAFT(1,xi,xk)*gl%meo2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)*gl%meo2x2_PCSAFT(1,xi,xj)) &
         & + gl%meo2_PCSAFT(1)*(-gl%c_PCSAFT(1)*gl%i2x3_PCSAFT(2,xi,xj,xk) - gl%c_PCSAFT(1)*gl%i2x3_PCSAFT(1,xi,xj,xk) - gl%c_PCSAFT(2)* &
         & gl%i2x3_PCSAFT(1,xi,xj,xk) - gl%cx1_PCSAFT(2,xj)*gl%i2x2_PCSAFT(1,xi,xk) - gl%cx2_PCSAFT(2,xi,xj)*gl%i2x1_PCSAFT(1,xk) - &
         & gl%cx3_PCSAFT(2,xi,xj,xk)*gl%i2_PCSAFT(1) - gl%cx2_PCSAFT(2,xj,xk)*gl%i2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(2,xi)* &
         & gl%i2x2_PCSAFT(1,xj,xk) - gl%cx2_PCSAFT(2,xi,xk)*gl%i2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(2,xk)*gl%i2x2_PCSAFT(1,xi,xj) - &
         & gl%cx1_PCSAFT(1,xj)*gl%i2x2_PCSAFT(2,xi,xk) - gl%cx1_PCSAFT(1,xj)*gl%i2x2_PCSAFT(1,xi,xk) - gl%cx2_PCSAFT(1,xi,xj)* &
         & gl%i2x1_PCSAFT(2,xk) - gl%cx2_PCSAFT(1,xi,xj)*gl%i2x1_PCSAFT(1,xk) - gl%cx3_PCSAFT(1,xi,xj,xk)*gl%i2_PCSAFT(1) - &
         & gl%cx3_PCSAFT(1,xi,xj,xk)*gl%i2_PCSAFT(2) - gl%cx2_PCSAFT(1,xj,xk)*gl%i2x1_PCSAFT(2,xi) - gl%cx2_PCSAFT(1,xj,xk)* &
         & gl%i2x1_PCSAFT(1,xi) - gl%cx1_PCSAFT(1,xi)*gl%i2x2_PCSAFT(2,xj,xk) - gl%cx1_PCSAFT(1,xi)*gl%i2x2_PCSAFT(1,xj,xk) - gl%cx2_PCSAFT(1,xi,xk) &
         & *gl%i2x1_PCSAFT(2,xj) - gl%cx2_PCSAFT(1,xi,xk)*gl%i2x1_PCSAFT(1,xj) - gl%cx1_PCSAFT(1,xk)*gl%i2x2_PCSAFT(2,xi,xj) - &
         & gl%cx1_PCSAFT(1,xk)*gl%i2x2_PCSAFT(1,xi,xj))))
                    end if
        	    end do
            end do 
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires c_PCSAFT(1), I_1(1), i2_PCSAFT(1), meo1_PCSAFT(1), meo2_PCSAFT(1), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), I_1(4), i2_PCSAFT(4)
    if (GETDERADISP(4) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	        gl%adispx3_PCSAFT(3,xi,xj,xk) = piPCSAFT*D*( &
                            & gl%c_PCSAFT(1) * (gl%i2_PCSAFT(1) * (&
                            & - gl%mPCSAFT(xj) * gl%meo2x2_PCSAFT(4,xi,xk) &
                            & - gl%mPCSAFT(xi) * gl%meo2x2_PCSAFT(4,xj,xk) &
                            & - gl%mPCSAFT(xk) * gl%meo2x2_PCSAFT(4,xi,xj)) &
                            & - gl%i2_PCSAFT(4) * gl%mPCSAFT(xj) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%i2_PCSAFT(4) * gl%mPCSAFT(xi) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%i2_PCSAFT(4) * gl%mPCSAFT(xk) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%i2x1_PCSAFT(4,xj) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%i2x1_PCSAFT(4,xj) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%i2x1_PCSAFT(4,xi) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%i2x1_PCSAFT(4,xi) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%i2x1_PCSAFT(4,xk) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%i2x1_PCSAFT(4,xk) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%i2x2_PCSAFT(1,xi,xj) * gl%mPCSAFT(xk) * gl%meo2_PCSAFT(4) &
                            & - gl%i2x2_PCSAFT(1,xj,xk) * gl%mPCSAFT(xi) * gl%meo2_PCSAFT(4) &
                            & - gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(4,xj) &
                            & - gl%i2x2_PCSAFT(1,xi,xk) * gl%mPCSAFT(xj) * gl%meo2_PCSAFT(4) &
                            & - gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(4,xj) &
                            & + gl%meo2_PCSAFT(1) * &
                            & (-gl%i2x2_PCSAFT(4,xi,xj) * gl%mPCSAFT(xk) &
                            & - gl%i2x2_PCSAFT(4,xj,xk) * gl%mPCSAFT(xi) &
                            & - gl%i2x2_PCSAFT(4,xi,xk) * gl%mPCSAFT(xj))) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xi) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xj) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xi) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xj) * gl%meo2_PCSAFT(4) &
                            & - 2.d0 * gl%i1x1_PCSAFT(4,xj) *gl%meo1x2_PCSAFT(1,xi,xk) &
                            & - 2.d0 * gl%i1x2_PCSAFT(4,xi,xj) * gl%meo1x1_PCSAFT(1,xk) &
                            & - 2.d0 * gl%i1x3_PCSAFT(3,xi,xj,xk) * gl%meo1_PCSAFT(1) &
                            & - 2.d0 * gl%i1x2_PCSAFT(4,xj,xk) * gl%meo1x1_PCSAFT(1,xi) &
                            & - 2.d0 * gl%i1x1_PCSAFT(4,xi) * gl%meo1x2_PCSAFT(1,xj,xk) &
                            & - 2.d0 * gl%i1x2_PCSAFT(4,xi,xk) * gl%meo1x1_PCSAFT(1,xj) &
                            & - 2.d0 * gl%i1x1_PCSAFT(4,xk) * gl%meo1x2_PCSAFT(1,xi,xj) &
                            & - 2.d0 * gl%i1x1_PCSAFT(1,xj) * gl%meo1x2_PCSAFT(4,xi,xk) &
                            & - 2.d0 * gl%i1x2_PCSAFT(1,xi,xj) * gl%meo1x1_PCSAFT(4,xk) &
                            & - 2.d0 * gl%i1x3_PCSAFT(1,xi,xj,xk) * gl%meo1_PCSAFT(4) &
                            & - 2.d0 * gl%i1x2_PCSAFT(1,xj,xk) * gl%meo1x1_PCSAFT(4,xi) &
                            & - 2.d0 * gl%i1x1_PCSAFT(1,xi) * gl%meo1x2_PCSAFT(4,xj,xk) &
                            & - 2.d0 * gl%i1x2_PCSAFT(1,xi,xk) * gl%meo1x1_PCSAFT(4,xj) &
                            & - 2.d0 * gl%i1x1_PCSAFT(1,xk) * gl%meo1x2_PCSAFT(4,xi,xj) &
                            & + gl%i2_PCSAFT(1) * ( &
                            & - gl%c_PCSAFT(4) * gl%mPCSAFT(xj) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%c_PCSAFT(4) * gl%mPCSAFT(xi) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%c_PCSAFT(4) * gl%mPCSAFT(xk) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%mPCSAFT(xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%mPCSAFT(xi) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%mPCSAFT(xk) * gl%meo2x1_PCSAFT(4,xj) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%mPCSAFT(xj) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%mPCSAFT(xj) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%mPCSAFT(xi) * gl%meo2x1_PCSAFT(4,xj) &
                            & + gl%meo2_PCSAFT(1) * ( &
                            & - gl%cx2_PCSAFT(4,xi,xj) * gl%mPCSAFT(xk) &
                            & - gl%cx2_PCSAFT(4,xj,xk) * gl%mPCSAFT(xi) &
                            & - gl%cx2_PCSAFT(4,xi,xk) * gl%mPCSAFT(xj))) &
                            & + gl%meo2_PCSAFT(1) * ( &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xi,xj) * gl%mPCSAFT(xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xj,xk) * gl%mPCSAFT(xi) &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xi,xk) * gl%mPCSAFT(xj) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xk) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xi) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xk) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%i2x1_PCSAFT(1,xk) * gl%mPCSAFT(xj) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%i2x1_PCSAFT(1,xj) * gl%mPCSAFT(xi) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%i2x1_PCSAFT(1,xi) * gl%mPCSAFT(xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(4,xi) * gl%mPCSAFT(xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(4,xk) * gl%mPCSAFT(xi) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xk) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xi) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(4,xj) * gl%mPCSAFT(xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(4,xk) * gl%mPCSAFT(xj) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%i2_PCSAFT(4) * gl%mPCSAFT(xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(4,xj) * gl%mPCSAFT(xi) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(4,xi) * gl%mPCSAFT(xj)) &
                            & + gl%mmean_PCSAFT * (gl%c_PCSAFT(1) * ( &
                            & - gl%i2x1_PCSAFT(4,xj) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%i2x2_PCSAFT(4,xi,xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%i2x2_PCSAFT(4,xj,xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%i2x1_PCSAFT(4,xi) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%i2x2_PCSAFT(4,xi,xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%i2x1_PCSAFT(4,xk) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%i2x1_PCSAFT(1,xj) * gl%meo2x2_PCSAFT(4,xi,xk) &
                            & - gl%i2x2_PCSAFT(1,xi,xj)*gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%i2x3_PCSAFT(1,xi,xj,xk) * gl%meo2_PCSAFT(4) &
                            & - gl%i2x2_PCSAFT(1,xj,xk) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%i2x1_PCSAFT(1,xi) * gl%meo2x2_PCSAFT(4,xj,xk) &
                            & - gl%i2x2_PCSAFT(1,xi,xk) * gl%meo2x1_PCSAFT(4,xj) &
                            & - gl%i2x1_PCSAFT(1,xk) * gl%meo2x2_PCSAFT(4,xi,xj)) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xj) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xi,xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xj,xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xi) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x2_PCSAFT(1,xi,xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%c_PCSAFT(4) * gl%i2x1_PCSAFT(1,xk) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%i2x1_PCSAFT(1,xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%i2x1_PCSAFT(1,xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%i2x1_PCSAFT(1,xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%i2x1_PCSAFT(1,xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%i2x1_PCSAFT(1,xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%i2x1_PCSAFT(1,xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2_PCSAFT(4) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(4,xi) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(4,xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(1,xi) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x2_PCSAFT(1,xi,xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x1_PCSAFT(1,xk) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%i2_PCSAFT(4) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%i2x1_PCSAFT(1,xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%i2_PCSAFT(4) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%i2x1_PCSAFT(1,xi) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2_PCSAFT(4) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(4,xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(4,xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(1,xj) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x2_PCSAFT(1,xj,xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x1_PCSAFT(1,xk) * gl%meo2x1_PCSAFT(4,xj) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%i2_PCSAFT(4) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%i2x1_PCSAFT(1,xj) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2_PCSAFT(4) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(4,xj) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(4,xi) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(1,xj) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x2_PCSAFT(1,xi,xj) * gl%meo2_PCSAFT(4) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x1_PCSAFT(1,xi) * gl%meo2x1_PCSAFT(4,xj) &
                            & + gl%i2_PCSAFT(1) * ( &
                            & - gl%cx1_PCSAFT(4,xj) * gl%meo2x2_PCSAFT(1,xi,xk) &
                            & - gl%cx2_PCSAFT(4,xi,xj) * gl%meo2x1_PCSAFT(1,xk) &
                            & - gl%cx2_PCSAFT(4,xj,xk) * gl%meo2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%meo2x2_PCSAFT(1,xj,xk) &
                            & - gl%cx2_PCSAFT(4,xi,xk) * gl%meo2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%meo2x2_PCSAFT(1,xi,xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%meo2x2_PCSAFT(4,xi,xk) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%meo2x1_PCSAFT(4,xk) &
                            & - gl%cx3_PCSAFT(1,xi,xj,xk) * gl%meo2_PCSAFT(4) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%meo2x1_PCSAFT(4,xi) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%meo2x2_PCSAFT(4,xj,xk) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%meo2x1_PCSAFT(4,xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%meo2x2_PCSAFT(4,xi,xj)) &
                            & + gl%meo2_PCSAFT(1) * ( &
                            & - gl%c_PCSAFT(1) * gl%i2x3_PCSAFT(3,xi,xj,xk) &
                            & - gl%c_PCSAFT(4) * gl%i2x3_PCSAFT(1,xi,xj,xk) &
                            & - gl%cx1_PCSAFT(4,xj) * gl%i2x2_PCSAFT(1,xi,xk) &
                            & - gl%cx2_PCSAFT(4,xi,xj) * gl%i2x1_PCSAFT(1,xk) &
                            & - gl%cx3_PCSAFT(3,xi,xj,xk) * gl%i2_PCSAFT(1) &
                            & - gl%cx2_PCSAFT(4,xj,xk) * gl%i2x1_PCSAFT(1,xi) &
                            & - gl%cx1_PCSAFT(4,xi) * gl%i2x2_PCSAFT(1,xj,xk) &
                            & - gl%cx2_PCSAFT(4,xi,xk) * gl%i2x1_PCSAFT(1,xj) &
                            & - gl%cx1_PCSAFT(4,xk) * gl%i2x2_PCSAFT(1,xi,xj) &
                            & - gl%cx1_PCSAFT(1,xj) * gl%i2x2_PCSAFT(4,xi,xk) &
                            & - gl%cx2_PCSAFT(1,xi,xj) * gl%i2x1_PCSAFT(4,xk) &
                            & - gl%cx3_PCSAFT(1,xi,xj,xk) * gl%i2_PCSAFT(4) &
                            & - gl%cx2_PCSAFT(1,xj,xk) * gl%i2x1_PCSAFT(4,xi) &
                            & - gl%cx1_PCSAFT(1,xi) * gl%i2x2_PCSAFT(4,xj,xk) &
                            & - gl%cx2_PCSAFT(1,xi,xk) * gl%i2x1_PCSAFT(4,xj) &
                            & - gl%cx1_PCSAFT(1,xk) * gl%i2x2_PCSAFT(4,xi,xj))))
                    end if
        	    end do
            end do 
        end do
    end if

    
!DEC$ END IF
end subroutine ADISPX3DERIVS





    end module pc_saft_ADISPX_derivs_module
