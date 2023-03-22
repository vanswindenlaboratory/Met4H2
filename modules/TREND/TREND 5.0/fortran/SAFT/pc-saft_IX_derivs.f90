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

    ! module for file pc-saft_IX_derivs.f90
    module pc_saft_IX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains


    

subroutine IX1DERIVS(gl,GETDERI, index)

! Henning Markgraf, June 2016

    ! I_1, I_2: integrals of the perturbation theory
    ! defined by eq. A.14 and A.15 in Gross, Sadowski 2001:
    ! I_1 = Sum_0_6(a_i*eta**i)     if index==1
    ! I_2 = Sum_0_6(b_i*eta**i)     if index==2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    integer, intent (in) :: index
    !output: i1x1_PCSAFT or i2x1_PCSAFT (module variables)
    !working variables
    integer :: i, xi
    double precision :: sum
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X1 DERIVATIVES of I
    !calculate the derivatives of I_1, I_2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)**(i-1)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**i
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(1,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(1,xi) = sum
            end select
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*(gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)**(i-1)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**i)
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(2,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(2,xi) = sum
            end select
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i-1)*(gl%ab_PCSAFT(index,i)*i**2*gl%z3x1_PCSAFT(1,xi) - gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)* &
			& i*gl%z3_PCSAFT(1) - gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(3,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(3,xi) = sum
            end select
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i-1)*(gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(4) + (gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(4)* &
                                        & gl%z3x1_PCSAFT(1,xi)*(i-1)/gl%z3_PCSAFT(1))*gl%ab_PCSAFT(index,i))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(4,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(4,xi) = sum
            end select
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i - 3)*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(i**2 - 3*i + 2 &
			& ) + gl%z3_PCSAFT(1)**2*(gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(5,xi) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(5)) + gl%z3_PCSAFT(1)*(2.d0* &
			& gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) - 2.d0*gl%ab_PCSAFT(index,i)* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i* &
			& gl%z3_PCSAFT(4)**2 - gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(4)**2))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(5,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(5,xi) = sum
            end select
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i - 2)*(gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(i - 1) + gl%z3_PCSAFT(1)* &
			& (gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(4,xi) + gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(6,xi) - gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(4,xi) + &
			& gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(4)))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(6,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(6,xi) = sum
            end select
        end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i - 3)*(gl%ab_PCSAFT(index,i)*i**3*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%ab_PCSAFT(index,i)*i** &
			& 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)**2* &
			& gl%z3x1_PCSAFT(5,xi) + 2.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi) - 4.d0*gl%ab_PCSAFT(index,i)*i* &
			& gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + &
			& 2.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi) + gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(7,xi) - &
			& gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(5,xi) - 2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi) &
			& + 2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%abx1_PCSAFT(index,i,xi)*i**2*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)**2 + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1)* &
			& gl%z3_PCSAFT(4)**2)
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(7,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(7,xi) = sum
            end select
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i - 1)*(gl%ab_PCSAFT(index,i)*i**3*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%ab_PCSAFT(index,i)*i**2*gl%z3x1_PCSAFT(1,xi) &
			& + 2.d0*gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i**2*gl%z3_PCSAFT(1) - 3.d0*gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1) + &
			& 2.d0*gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(8,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(8,xi) = sum
            end select
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF I_1, I_2 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i*gl%z3_PCSAFT(1)**(i - 4)*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(i**3 - 6*i**2 + &
			& 11*i - 6) + gl%z3_PCSAFT(1)**3*(gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(9,xi) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(9)) + &
			& gl%z3_PCSAFT(1)**2*(3.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(5)* &
			& gl%z3x1_PCSAFT(4,xi) + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)* &
			& gl%z3x1_PCSAFT(5,xi) - 3.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(9)* &
			& gl%z3x1_PCSAFT(1,xi) + 3.d0*gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3.d0*gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(4)* &
			& gl%z3_PCSAFT(5)) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(3.d0*gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 3.d0 &
			& *gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) - 9.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - 9.d0 &
			& *gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + 6.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 6.d0*gl%ab_PCSAFT(index,i)* &
			& gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i**2*gl%z3_PCSAFT(4)**2 - 3.d0*gl%abx1_PCSAFT(index,i,xi)*i* &
			& gl%z3_PCSAFT(4)**2 + 2.d0*gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(4)**2))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(9,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(9,xi) = sum
            end select
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
        do xi = 1, gl%ncomp
            sum = 0.d0
            do i = 0, 6
                sum = sum + i**2*gl%z3_PCSAFT(1)**(i - 2)*(gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) - 2.d0*gl%ab_PCSAFT(index,i)*i* &
			& gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%z3_PCSAFT(1)*(gl%ab_PCSAFT(index,i)*i* &
			& gl%z3x1_PCSAFT(4,xi) - gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(4,xi) + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(4) - gl%abx1_PCSAFT(index,i,xi)* &
			& gl%z3_PCSAFT(4)))
            end do
            select case (index)
	        case (1)
		    gl%i1x1_PCSAFT(10,xi) = sum
	        case (2)
		    gl%i2x1_PCSAFT(10,xi) = sum
            end select
        end do
    end if
    
    !DEC$ END IF
end subroutine IX1DERIVS
    
    

subroutine IX2DERIVS(gl,GETDERI, index)

    ! I_1, I_2: integrals of the perturbation theory
    ! defined by eq. A.14 and A.15 in Gross, Sadowski 2001:
    ! I_1 = Sum_0_6(a_i*eta**i)     if index==1
    ! I_2 = Sum_0_6(b_i*eta**i)     if index==2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    integer, intent (in) :: index
    !output: i1x2_PCSAFT, i2x2_PCSAFT (module variables)
    !working variables
    integer :: i, xi, xj
    double precision :: sum
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X2 DERIVATIVES of I
    !calculate the derivatives of I_1, I_2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
!                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + gl%ab_PCSAFT(index,i)*i*(i - 1)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                                & + gl%abx1_PCSAFT(index,i,xj)*i*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(1,xi) &
                                & + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(1,xj) &
                                & + gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**i
		            end do
				    select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(1,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(1,xi,xj) = sum
                    end select
 !               end if
	        end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                !if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + i*gl%z3_PCSAFT(1)**(i - 2)*(gl%ab_PCSAFT(index,i)*i**2*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%ab_PCSAFT(index,i)*i* &
                 & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xj)*i*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i* &
                 & gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**2)
		            end do
		            select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(2,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(2,xi,xj) = sum
                    end select
                !end if
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + i*gl%z3_PCSAFT(1)**(i - 2)*(gl%ab_PCSAFT(index,i)*i**3*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 2.d0*gl%ab_PCSAFT(index,i)*i**2 &
                 & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xj)*i &
                 & **2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xj)*i*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i** &
                 & 2*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) - gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj) + gl%abx2_PCSAFT(index,i,xi,xj)*i* &
                 & gl%z3_PCSAFT(1)**2 - gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**2)
		            end do
		            select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(3,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(3,xi,xj) = sum
                    end select
                end if
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + i*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 2)* &
                 & (i - 1) + gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*(i - 1) + &
                 & gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*(i - 1) + gl%abx2_PCSAFT(index,i,xi,xj)* &
                 & gl%z3_PCSAFT(1)**i*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(1)*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(1,xj)* &
                 & gl%z3x1_PCSAFT(4,xi)*(i - 1) + gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(4,xj)* &
                 & gl%z3x1_PCSAFT(1,xi)*(i - 1) + gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(4,xi) + &
                 & gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(4,xj)))/gl%z3_PCSAFT(1)
		            end do
		            select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(4,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(4,xi,xj) = sum
                    end select
                end if
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + i*(gl%z3_PCSAFT(1)**2*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(5,xi)*(i &
                 & - 1) + 2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi)*(i - 1) + &
                 & gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 1) + gl%abx1_PCSAFT(index,i,xj)* &
                 & gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(5,xi) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)* &
                 & gl%z3x1_PCSAFT(5,xj)) + gl%z3_PCSAFT(1)*(2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3_PCSAFT(4)* &
                 & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi)*(i - 2)*(i - 1) + 2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)* &
                 & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 2)*(i - 1) + gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i &
                 & - 2)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 2)*(i - 1) + 2.d0*gl%abx1_PCSAFT(index,i,xj)* &
                 & gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*(i - 1) + gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**( &
                 & i - 1)*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi)*(i - 1) + 2.d0*gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)* &
                 & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xj)*(i - 1) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(5) &
                 & *gl%z3x1_PCSAFT(1,xj)*(i - 1) + gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**i*gl%z3_PCSAFT(5)) + gl%z3_PCSAFT(4) &
                 & **2*(gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 2)**2*(i - 1 &
                 & ) - gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**(i - 2)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*(i - 2)*(i - 1) + &
                 & gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(1,xi)*(i - 1)**2 - gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**( &
                 & i - 1)*gl%z3x1_PCSAFT(1,xi)*(i - 1) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(1,xj)*(i &
                 & - 1)**2 - gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**(i - 1)*gl%z3x1_PCSAFT(1,xj)*(i - 1) + gl%abx2_PCSAFT(index,i,xi,xj) &
                 & *i*gl%z3_PCSAFT(1)**i - gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**i))/gl%z3_PCSAFT(1)**2
		            end do
		            select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(5,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(5,xi,xj) = sum
                    end select
                end if
	        end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
		            sum = 0.d0
		            do i = 0, 6
		                sum = sum + i**2*gl%z3_PCSAFT(1)**(i - 3)*(gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                    & + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)* &
                    & gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - 3.d0*gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                    & - gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(4,xi) - gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(4,xj)* &
                    & gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xj)*i* &
                    & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xi) - &
                    & gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
                    & gl%z3x1_PCSAFT(1,xj) + gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)**2*gl%z3x1_PCSAFT(4,xj) - gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(1)* &
                    & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4))
		            end do
		            select case (index)
			    	    case (1)
						    gl%i1x2_PCSAFT(6,xi,xj) = sum
			    	    case (2)
						    gl%i2x2_PCSAFT(6,xi,xj) = sum
                    end select
                end if
            end do
        end do
    end if

    !DEC$ END IF
end subroutine IX2DERIVS
    


subroutine IX3DERIVS(gl,GETDERI, index)
    
! Henning Markgraf, June 2016

    ! I_1, I_2: integrals of the perturbation theory
    ! defined by eq. A.14 and A.15 in Gross, Sadowski 2001:
    ! I_1 = Sum_0_6(a_i*eta**i)     if index==1
    ! I_2 = Sum_0_6(b_i*eta**i)     if index==2
    ! dependent on D and T






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    integer, intent (in) :: index
    !output: i1x3_PCSAFT, i2x3_PCSAFT (module variables) depending on input index
    !working variables
    integer :: i, xi, xj, xk
    double precision :: sum
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. X3 DERIVATIVES of I
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 0, 6
		                    sum = sum + gl%z3_PCSAFT(1)**(i - 3)*(gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(i**2 &
                 & - 3*i + 2) + gl%abx3_PCSAFT(index,i,xi,xj,xk)*gl%z3_PCSAFT(1)**3 + i*gl%z3_PCSAFT(1)**2*(gl%abx2_PCSAFT(index,i,xi,xj)* &
                 & gl%z3x1_PCSAFT(1,xk) + gl%abx2_PCSAFT(index,i,xj,xk)*gl%z3x1_PCSAFT(1,xi) + gl%abx2_PCSAFT(index,i,xi,xk)*gl%z3x1_PCSAFT(1,xj)) + i* &
                 & gl%z3_PCSAFT(1)*(gl%abx1_PCSAFT(index,i,xj)*i*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xj)*gl%z3x1_PCSAFT(1,xi)* &
                 & gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xi)*gl%z3x1_PCSAFT(1,xj) &
                 & *gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xk)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xk)* &
                 & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))
		                end do
					    select case (index)
						    case (1)
							    gl%i1x3_PCSAFT(1,xi,xj,xk) = sum
						    case (2)
							    gl%i2x3_PCSAFT(1,xi,xj,xk) = sum
                        end select
                    end if
		        end do
	        end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 0, 6
		                    sum = sum + i*gl%z3_PCSAFT(1)**(i - 3)*(gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(i** &
                 & 2 - 3*i + 2) + gl%abx3_PCSAFT(index,i,xi,xj,xk)*gl%z3_PCSAFT(1)**3 + i*gl%z3_PCSAFT(1)**2*(gl%abx2_PCSAFT(index,i,xi,xj)* &
                 & gl%z3x1_PCSAFT(1,xk) + gl%abx2_PCSAFT(index,i,xj,xk)*gl%z3x1_PCSAFT(1,xi) + gl%abx2_PCSAFT(index,i,xi,xk)*gl%z3x1_PCSAFT(1,xj)) + i* &
                 & gl%z3_PCSAFT(1)*(gl%abx1_PCSAFT(index,i,xj)*i*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xj)*gl%z3x1_PCSAFT(1,xi)* &
                 & gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xi)*gl%z3x1_PCSAFT(1,xj) &
                 & *gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xk)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xk)* &
                 & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)))
		                end do
		                select case (index)
						    case (1)
							    gl%i1x3_PCSAFT(2,xi,xj,xk) = sum
						    case (2)
							    gl%i2x3_PCSAFT(2,xi,xj,xk) = sum
                        end select
                    end if
		        end do
	        end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
		                sum = 0.d0
		                do i = 0, 6
		                    sum = sum + i*gl%z3_PCSAFT(1)**(i - 4)*(-3*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
                 & gl%z3x1_PCSAFT(1,xk)*(i**2 - 3*i + 2) + gl%z3_PCSAFT(1)**3*(gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3x1_PCSAFT(4,xk) + &
                 & gl%abx2_PCSAFT(index,i,xj,xk)*gl%z3x1_PCSAFT(4,xi) + gl%abx2_PCSAFT(index,i,xi,xk)*gl%z3x1_PCSAFT(4,xj)) + gl%z3_PCSAFT(1)**2*( &
                 & gl%abx1_PCSAFT(index,i,xj)*i*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xj)*i*gl%z3x1_PCSAFT(4,xk)* &
                 & gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xj)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xj)*gl%z3x1_PCSAFT(4,xk) &
                 & *gl%z3x1_PCSAFT(1,xi) - gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xk) - gl%abx2_PCSAFT(index,i,xj,xk)* &
                 & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xi)*i*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) + gl%abx1_PCSAFT(index,i,xi)* &
                 & i*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) - gl%abx1_PCSAFT(index,i,xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk) - &
                 & gl%abx1_PCSAFT(index,i,xi)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj) - gl%abx2_PCSAFT(index,i,xi,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj) + &
                 & gl%abx1_PCSAFT(index,i,xk)*i*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) + gl%abx1_PCSAFT(index,i,xk)*i*gl%z3x1_PCSAFT(4,xi)* &
                 & gl%z3x1_PCSAFT(1,xj) - gl%abx1_PCSAFT(index,i,xk)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xk)*gl%z3x1_PCSAFT(4,xi) &
                 & *gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(1)*(gl%ab_PCSAFT(index,i)*i**2*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
                 & + gl%ab_PCSAFT(index,i)*i**2*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + gl%ab_PCSAFT(index,i)*i**2* &
                 & gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 3*gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(4,xj)* &
                 & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - 3*gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) &
                 & - 3*gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) + 2*gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(4,xj)* &
                 & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2*gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + &
                 & 2*gl%ab_PCSAFT(index,i)*gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - 2*gl%abx1_PCSAFT(index,i,xj)*i*gl%z3_PCSAFT(4)* &
                 & gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + 2*gl%abx1_PCSAFT(index,i,xj)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) &
                 & - 2*gl%abx1_PCSAFT(index,i,xi)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + 2*gl%abx1_PCSAFT(index,i,xi)*gl%z3_PCSAFT(4) &
                 & *gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - 2*gl%abx1_PCSAFT(index,i,xk)*i*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)* &
                 & gl%z3x1_PCSAFT(1,xi) + 2*gl%abx1_PCSAFT(index,i,xk)*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)) + gl%z3_PCSAFT(4) &
                 & *(gl%ab_PCSAFT(index,i)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*(i**2 - 3*i + 2) + &
                 & gl%abx3_PCSAFT(index,i,xi,xj,xk)*gl%z3_PCSAFT(1)**3 + i*gl%z3_PCSAFT(1)**2*(gl%abx2_PCSAFT(index,i,xi,xj)*gl%z3x1_PCSAFT(1,xk) + &
                 & gl%abx2_PCSAFT(index,i,xj,xk)*gl%z3x1_PCSAFT(1,xi) + gl%abx2_PCSAFT(index,i,xi,xk)*gl%z3x1_PCSAFT(1,xj)) + i*gl%z3_PCSAFT(1)*( &
                 & gl%abx1_PCSAFT(index,i,xj)*i*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk) + &
                 & gl%abx1_PCSAFT(index,i,xi)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) - gl%abx1_PCSAFT(index,i,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk) + &
                 & gl%abx1_PCSAFT(index,i,xk)*i*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) - gl%abx1_PCSAFT(index,i,xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi) &
                 & )))
		                end do
		                select case (index)
						    case (1)
							    gl%i1x3_PCSAFT(3,xi,xj,xk) = sum
						    case (2)
							    gl%i2x3_PCSAFT(3,xi,xj,xk) = sum
                        end select
                    end if
		        end do
	        end do
        end do
    end if
        
    !DEC$ END IF
end subroutine IX3DERIVS
    



    end module pc_saft_IX_derivs_module
