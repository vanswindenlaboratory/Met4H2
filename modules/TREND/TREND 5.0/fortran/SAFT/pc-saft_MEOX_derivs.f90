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

    ! module for file pc-saft_MEOX_derivs.f90
    module pc_saft_MEOX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains


    
subroutine MEOX1DERIVS(gl,T, GETDERMEO, index)

! Henning Markgraf, June 2016

    ! meo1, meo2: abbreviations
    ! defined by eq. A.12 and A.13 in Gross, Sadowski 2001:
    ! meo1 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))*sigma_ij**3))     if index==1
    ! meo2 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))**2*sigma_ij**3))  if index==2
    ! combining rules for sigma_ij and epsilon_ij
    ! sigma_ij = 1/2*(sigma_i+sigma_j)
    ! epsilon_ij = sqrt(epsilon_i*epsilon_j)*(1-kij)
    ! dependent only on T

    ! note: kij, module variable in module_PCSAFT
    !                   can be used to correlate binary mixtures better






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    double precision, intent (in) :: T
    integer, dimension (nderivs), intent (in) :: GETDERMEO
    integer, intent (in) :: index
    !output: meo1x1_PCSAFT or meo2x1_PCSAFT (module variables)
    !working variables
    double precision :: sigma_ij, epsk_ij
    double precision, dimension (gl%ncomp) :: basiceqn, derivative
    integer :: xi, j
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Derivative of the basic equation (meo1 or meo2 depending on index) wrt x_i
    ! All T derivatives of meo1 or meo2 are just this term times a factor
    basiceqn = 0.d0
    do xi = 1, gl%ncomp
        do j = 1, gl%ncomp
            sigma_ij = 1.d0/2.d0*(gl%sigPCSAFT(xi)+gl%sigPCSAFT(j))
            epsk_ij = dsqrt(gl%epskPCSAFT(xi)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(xi,j)) ! note: where does kij get its values?
            basiceqn(xi) = basiceqn(xi) + sigma_ij**3*(epsk_ij/T)**index*gl%mPCSAFT(j)*gl%molfractions(j)
        end do
	    basiceqn(xi) = 2.d0*gl%mPCSAFT(xi)*basiceqn(xi)
    end do

    ! IV. X1 DERIVATIVES of meo
    ! calculate the derivatives of meo1, meo2 depending on index
    ! 1: meo1 or meo2
    if (GETDERMEO(1) .eq. 1) then
        derivative = basiceqn
	select case (index)
            case (1)
                gl%meo1x1_PCSAFT(1,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x1_PCSAFT(1,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 2, 3: derivatives with respect to D equal zero
    
    ! 4: 1ST DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! T*d(meo_index)/dT = -index*meo_index
    if (GETDERMEO(4) .eq. 1) then
        derivative = -index*basiceqn
	select case (index)
            case (1)
                gl%meo1x1_PCSAFT(4,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x1_PCSAFT(4,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 5: 2ND DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    if (GETDERMEO(5) .eq. 1) then
        derivative = index*(index+1)*basiceqn
	select case (index)
            case (1)
                gl%meo1x1_PCSAFT(5,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x1_PCSAFT(5,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 6, 7, 8: derivatives with respect to D equal zero
    
    ! 9: 3RD DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3 
    if (GETDERMEO(9) .eq. 1) then
        derivative = -index*(index**2+3*index+2)*basiceqn
	select case (index)
            case (1)
                gl%meo1x1_PCSAFT(9,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x1_PCSAFT(9,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 10: derivatives with respect to D equal zero
    
    !DEC$ END IF
end subroutine MEOX1DERIVS


        
subroutine MEOX2DERIVS(gl,T, GETDERMEO, index)

! Henning Markgraf, June 2016
    
    ! meo1, meo2: abbreviations
    ! defined by eq. A.12 and A.13 in Gross, Sadowski 2001:
    ! meo1 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))*sigma_ij**3))     if index==1
    ! meo2 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))**2*sigma_ij**3))  if index==2
    ! combining rules for sigma_ij and epsilon_ij
    ! sigma_ij = 1/2*(sigma_i+sigma_j)
    ! epsilon_ij = sqrt(epsilon_i*epsilon_j)*(1-kij)
    ! dependent only on T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    double precision, intent (in) :: T
    integer, dimension (nderivs), intent (in) :: GETDERMEO
    integer, intent (in) :: index
    !output: meo1x2_PCSAFT or meo2x2_PCSAFT (module variables) depending on input index
    !working variables
    double precision :: sigma_ij, epsk_ij
    double precision, dimension (gl%ncomp,gl%ncomp) :: basiceqn, derivative
    integer :: xi, xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Derivative of the basic equation wrt x
    ! All derivatives of meo1 or meo2 wrt T are just this basic equation times a factor
    basiceqn = 0.d0
    do xj = 1, gl%ncomp
        do xi = 1, gl%ncomp
  !          if (xi .GE. xj) then    
                sigma_ij = 1.d0/2.d0*(gl%sigPCSAFT(xi)+gl%sigPCSAFT(xj))
                epsk_ij = dsqrt(gl%epskPCSAFT(xi)*gl%epskPCSAFT(xj))*(1.d0-gl%kij_PCSAFT(xi,xj)) ! note: where does kij get its value?
                !basiceqn(xi,xj) = mPCSAFT(xi)*mPCSAFT(xj)*sigma_ij**3*(epsk_ij/T)**index
                basiceqn(xi,xj) = 2.d0*gl%mPCSAFT(xi)*gl%mPCSAFT(xj)*sigma_ij**3*(epsk_ij/T)**index
!            end if
        end do
    end do

    ! IV. X2 DERIVATIVES of meo
    !calculate the derivatives of meo1, meo2
    ! 1: meo1 or meo2
    if (GETDERMEO(1) .eq. 1) then
        derivative = basiceqn
		select case (index)
            case (1)
                gl%meo1x2_PCSAFT(1,1:gl%ncomp,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x2_PCSAFT(1,1:gl%ncomp,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 2, 3: derivatives with respect to D equal zero
    
    ! 4: 1ST DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! T*d(meo_index)/dT = -index*meo_index
    if (GETDERMEO(4) .eq. 1) then
        derivative = -index*basiceqn
		select case (index)
            case (1)
                gl%meo1x2_PCSAFT(4,1:gl%ncomp,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x2_PCSAFT(4,1:gl%ncomp,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 5: 2ND DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    if (GETDERMEO(5) .eq. 1) then
        derivative = index*(index+1)*basiceqn
		select case (index)
            case (1)
                gl%meo1x2_PCSAFT(5,1:gl%ncomp,1:gl%ncomp) = derivative
            case (2)
                gl%meo2x2_PCSAFT(5,1:gl%ncomp,1:gl%ncomp) = derivative
        end select
    end if
    
    ! 6: derivatives with respect to D equal zero

    !DEC$ END IF
end subroutine MEOX2DERIVS




    end module pc_saft_MEOX_derivs_module
