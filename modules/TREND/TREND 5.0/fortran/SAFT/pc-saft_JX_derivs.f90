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

    ! module for file pc-saft_JX_derivs.f90
    module pc_saft_JX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains



subroutine J2X1DERIVS(gl,T,DENS,GETDERI)
    
    ! Ina Schülling, September 2017

    ! J2: integral of the perturbation theory
    ! defined by eq. eq 12 in Gross 2005:
    ! J2 = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2x1_PCSAFT (module variables)
    integer :: i,j,n,xi
    double precision :: eps_ij, sum1, sum2, sum3, sum4
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    eps_ij = 0.d0
    
    gl%J2X1_PCSAFTQ = 0.d0
    
    ! III. X1 DERIVATIVES of I
    !calculate the derivatives of J2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
    do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTQ(1,i,j,xi) = gl%J2X1_PCSAFTQ(1,i,j,xi) +  n*gl%z3_PCSAFT(1)**(n - 1.d0)*gl%z3x1_PCSAFT(1,xi)*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T) 
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
      do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTQ(2,i,j,xi) = gl%J2X1_PCSAFTQ(2,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 2.d0)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                            & *(gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(2,xi) + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0))/T
                        !J2X1_PCSAFTQ(2,i,j,xi) = J2X1_PCSAFTQ(2,i,j,xi) + n**2 *z3_PCSAFT(1)**(n - 1)* z3x1_PCSAFT(1,xi)*(ab_PCSAFTQ(1,n,i,j) + ab_PCSAFTQ(2,n,i,j)*eps_ij/T)

                    end do
                end do
            end do
        end do
    end if  
   
    
    ! 3: 2ND DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
      do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTQ(3,i,j,xi) = gl%J2X1_PCSAFTQ(3,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * gl%z3_PCSAFT(2) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij)  &
                             & * (2.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n - 1.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0)**2) / T  
                    end do
                end do
            end do
        end do  
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTQ(4,i,j,xi) = gl%J2X1_PCSAFTQ(4,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 2) * (T * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(4,xi) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij)  &
                             & + T * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * (n - 1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                             & - T * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xi)) / T**2
                    end do
                end do
            end do
        end do  
       
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
       do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTQ(5,i,j,xi) = J2X1_PCSAFTQ(5,i,j,xi) + n*(ab_PCSAFTQ(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(5,xi)&
                        !&+ 2*ab_PCSAFTQ(2,n,i,j)*eps_ij*z3_PCSAFT(1)**(n - 3)*z3_PCSAFT(4)**2*z3x1_PCSAFT(1,xi)/T)
                    
                        gl%J2X1_PCSAFTQ(5,i,j,xi) = gl%J2X1_PCSAFTQ(5,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (T * (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) * (T * gl%ab_PCSAFTQ(1,n,i,j) + &
                             & gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) * (2.d0 * n * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(4,xi) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) + 2.d0 * &
                             & T * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) - 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * &
                             & eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(4,xi) * T + 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(1,xi) * T) + gl%z3x1_PCSAFT(1,xi) * (T**2 * n * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * &
                             & eps_ij) + T**2 * n * gl%z3_PCSAFT(4)**2 * (n - 1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) - 2.d0 * T * &
                             & gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * n * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * T + 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * T**2 )) / T**3

                    end do
                end do
            end do
        end do  
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTQ(6,i,j,xi) = gl%J2X1_PCSAFTQ(6,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(6,xi) * (T * gl%ab_PCSAFTQ(1,n,i,j) + &
                             & gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) * (n * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(2,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + n * gl%z3x1_PCSAFT(4,xi) * &
                             & gl%z3_PCSAFT(2) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - &
                             & gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) + T * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (T * &
                             & gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) * (n**2 - 3.d0 * n + 2.d0) - gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * &
                             & gl%z3x1_PCSAFT(2,xi) * T + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * T * (-n + 1.d0))/ &
                             & T**2
                        !J2X1_PCSAFTQ(6,i,j,xi) = J2X1_PCSAFTQ(6,i,j,xi) + n*(ab_PCSAFTQ(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(6,xi)- ab_PCSAFTQ(2,n,i,j)*eps_ij*n*z3_PCSAFT(1)**(n - 2)*z3_PCSAFT(4)*z3x1_PCSAFT(1,xi)/T)

                    end do
                end do
            end do
        end do  
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTQ(7,i,j,xi) = J2X1_PCSAFTQ(7,i,j,xi) + n*(ab_PCSAFTQ(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(7,xi)+ 2*ab_PCSAFTQ(2,n,i,j)*eps_ij*n*z3_PCSAFT(1)**(n - 3)*z3_PCSAFT(4)**2*z3x1_PCSAFT(1,xi)/T)
                        
                        gl%J2X1_PCSAFTQ(7,i,j,xi) = gl%J2X1_PCSAFTQ(7,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (T * n * (gl%z3_PCSAFT(2) *&
                             & (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                             & + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                             & * (2 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) &
                             & + 2.d0 * T * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij)  &
                             & - 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(4,xi) * T + 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * T)  &
                             
                             & + gl%z3x1_PCSAFT(1,xi) * (T * gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij)  &
                             & + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) * (2.d0 * n * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) - gl%z3_PCSAFT(5) * gl%z3_PCSAFT(2)) &
                             & + 2.d0 * T * gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                             & - 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(6) * T + 2.d0 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * T)) &
                             
                             & + T*(T*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(7,xi)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)  &
                             & + T*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(2.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi) &
                             & + 2.d0*n*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(4,xi) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi)  &
                             & - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(2,xi) - gl%z3_PCSAFT(7)*gl%z3x1_PCSAFT(1,xi) - gl%z3x1_PCSAFT(5,xi)*gl%z3_PCSAFT(2)  &
                             & - 2.d0*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(4,xi)) + 2.d0*T*gl%z3_PCSAFT(1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                             & * (-n*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(2,xi) - 2.d0*n*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) - 2.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2)  &
                             & + gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(2,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2)+ gl%z3_PCSAFT(5)*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi))   &
                             
                             & + 6.d0*T*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                             & - 2.d0*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(6,xi)*T + 2.d0*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2 &
                             & * (gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) + gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi)+ gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2))*T  &
                             & - 4.d0*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*T) &
                             & + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(2,xi)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                             & + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1.d0)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) - 2.d0*T*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*T &
                             & + 2.d0*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*T**2) &
                             & + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTQ(1,n,i,j)+ gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)  &
                             & + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1.d0)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                             &  - 2.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2)) &
                             & /T**3

                    end do
                end do
            end do
        end do  
    end if
    
    ! 8: 3RD DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTQ(8,i,j,xi) = J2X1_PCSAFTQ(8,i,j,xi) + (n**2*z3_PCSAFT(1)**n*z3x1_PCSAFT(1,xi)*(ab_PCSAFTQ(1,n,i,j)+ ab_PCSAFTQ(2,n,i,j)*eps_ij/T)*(n**2*z3_PCSAFT(1) - 3*n*z3_PCSAFT(1) + 2*z3_PCSAFT(1))/z3_PCSAFT(1)**2 &
                        !    & + n*z3_PCSAFT(1)**n*(ab_PCSAFTQ(1,n,i,j)+ ab_PCSAFTQ(2,n,i,j)*eps_ij/T)*(n**2*z3x1_PCSAFT(1,xi) - 3*n*z3x1_PCSAFT(1,xi) + 2*z3x1_PCSAFT(1,xi))/z3_PCSAFT(1) &
                        !    & - n*z3_PCSAFT(1)**n*z3x1_PCSAFT(1,xi)*(ab_PCSAFTQ(1,n,i,j)+ ab_PCSAFTQ(2,n,i,j)*eps_ij/T)*(n**2*z3_PCSAFT(1) - 3*n*z3_PCSAFT(1) + 2*z3_PCSAFT(1))/z3_PCSAFT(1)**2)
                        gl%J2X1_PCSAFTQ(8,i,j,xi) = gl%J2X1_PCSAFTQ(8,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * gl%z3_PCSAFT(2)**2 * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                        & * (3.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n**2 - 3.d0 * n + 2.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n**2 - 3.d0 * n + 2.d0) + 2.d0 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (-n**2 + 3.d0 * n - 2.d0)) / T
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J2 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
       do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTQ(9,i,j,xi) = J2X1_PCSAFTQ(9,i,j,xi) + n*(ab_PCSAFTQ(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(9,xi)&
                        !    & +6*ab_PCSAFTQ(2,n,i,j)*eps_ij*z3_PCSAFT(1)**(n - 3)*z3_PCSAFT(4)*z3_PCSAFT(5)*z3x1_PCSAFT(1,xi)/T)
                    
                        gl%J2X1_PCSAFTQ(9,i,j,xi) = gl%J2X1_PCSAFTQ(9,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 4)*(T*(T**2*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(9,xi)*(T* &
                             & gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) + T**2*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(3.d0*n* &
                             & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*n*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - 3.d0*gl%z3_PCSAFT(4)* &
                             & gl%z3x1_PCSAFT(5,xi) - 3.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi))  &
                             
                             & + 3.d0*T**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(n**2*gl%z3_PCSAFT(4)* &
                             & gl%z3x1_PCSAFT(4,xi) - 3.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - 2.d0*n*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
                             & + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi))  &
                             
                             & + 3.d0*T**2*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(-n**2 + 3.d0*n - 2.d0) &
                             & - 3.d0 *T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi) + 3.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*(-2.d0 &
                             & *n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(5)* &
                             & gl%z3x1_PCSAFT(1,xi)) &
                             & + 6.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0 &
                             & ) + 6.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) - 6.d0*T**2*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + gl%z3x1_PCSAFT(1,xi)*(T**3*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9)*(T &
                             & *gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) + 3.d0*T**3*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1.d0)*( &
                             & T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) + T**3*n*gl%z3_PCSAFT(4)**3*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(n &
                             & **2 - 3.d0*n + 2) - 3.d0*T**3*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 3.d0*T**3* &
                             & gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(-n + 1.d0) + 6.d0*T**3*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1) &
                             & **2*gl%z3_PCSAFT(4) - 6.d0*T**3*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3))/T**4

                    end do
                end do
            end do
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
      do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTQ(10,i,j,xi) = J2X1_PCSAFTQ(10,i,j,xi) + (ab_PCSAFTQ(1,n,i,j)*n* z3_PCSAFT(1)**(n - 1)* z3x1_PCSAFT(10,xi) &
                        !    & +2*ab_PCSAFTQ(2,n,i,j)* eps_ij *n**2* z3_PCSAFT(1)**(n - 3)* z3_PCSAFT(4)**2* z3x1_PCSAFT(1,xi)/T)/(T*Dens**2)
                    gl%J2X1_PCSAFTQ(10,i,j,xi) = gl%J2X1_PCSAFTQ(10,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 4.d0)*(2.d0*T*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                         & *(n*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(2,xi) + n*gl%z3x1_PCSAFT(6,xi)*gl%z3_PCSAFT(2) - &
                         & gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(2,xi) - gl%z3x1_PCSAFT(6,xi)*gl%z3_PCSAFT(2)) + T*gl%z3_PCSAFT(1) &
                         & *gl%z3_PCSAFT(2)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(2.d0*n**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) + &
                         & n**2*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2) - 6.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) - 2.d0*n* &
                         & gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) - 3.d0*n*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2) + 4.d0*gl%z3_PCSAFT(4)* &
                         & gl%z3x1_PCSAFT(2,xi) + 2.d0*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)* &
                         & gl%z3_PCSAFT(2)) + 2.d0*T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)**2*gl%z3x1_PCSAFT(1,xi)*(T*gl%ab_PCSAFTQ(1,n,i,j) + &
                         & gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(-n**2 + 3.d0*n - 2.d0) + 2.d0*T*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(2)* &
                         & gl%z3x1_PCSAFT(2,xi)*(-n + 1.d0) + T*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(2)**2* &
                         & gl%z3x1_PCSAFT(1,xi)*(n - 1.d0) + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(2.d0*T*gl%z3_PCSAFT(1)* &
                         & gl%z3_PCSAFT(6)*(n - 1.d0)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) + T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)* &
                         & (T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(n**2 - 3.d0*n + 2.d0) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)* &
                         & gl%z3_PCSAFT(2)*T*(-n + 1.d0)))/T**2

                    end do
                end do
            end do
        end do
    end if
    
    !DEC$ END IF
end subroutine J2X1DERIVS
      
subroutine J2X2DERIVS(gl,T,DENS,GETDERI)
    
    ! Ina Schülling, September 2017

    ! J2: integral of the perturbation theory
    ! defined by eq. eq 12 in Gross 2005:
    ! J2 = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2x1_PCSAFT (module variables)
    integer :: i,j,n,xi,xj
    double precision :: eps_ij, sum1, sum2, sum3, sum4
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    eps_ij = 0.d0
    
    gl%J2X2_PCSAFTQ = 0.d0
    
    ! III. X1 DERIVATIVES of I
    !calculate the derivatives of J2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        do xj = 1, gl%ncomp
            do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(1,i,j,xi,xj) = gl%J2X2_PCSAFTQ(1,i,j,xi,xj) + n * (n - 1.d0) * gl%z3_PCSAFT(1)**(n - 2.d0) * gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi) * (gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij / T) 
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
    do xj = 1, gl%ncomp
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(2,i,j,xi,xj) = gl%J2X2_PCSAFTQ(2,i,j,xi,xj) +  n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) * (gl%z3_PCSAFT(1)*(gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) &
                                & + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0)) / T
                        end do
                    end do
                end do
            end do
        end do 
    end if  
   
    
    ! 3: 2ND DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
    do xj = 1, gl%ncomp  
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(3,i,j,xi,xj) = gl%J2X2_PCSAFTQ(3,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (2.d0 * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(2,xj) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (n - 2.d0) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0)) / T
                        end do
                    end do
                end do
            end do 
        end do 
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
    do xj = 1, gl%ncomp   
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(4,i,j,xi,xj) = gl%J2X2_PCSAFTQ(4,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) + T * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0) &
                                & * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) - T * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) / T**2
                        end do
                    end do
                end do
            end do  
        end do 
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
    do xj = 1, gl%ncomp  
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(5,i,j,xi,xj) = gl%J2X2_PCSAFTQ(5,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T**2 * gl%z3_PCSAFT(1)**2 * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(5,xj) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(5,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T**2 * gl%z3_PCSAFT(1) * (n - 2.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T**2 * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & - 2.d0 * T**2 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & - 2.d0 * T**2 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0) &
                                & + 2.d0 * T**2 * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) / T**3
                        end do
                    end do
                end do
            end do
        end do 
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
    do xj = 1, gl%ncomp   
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTQ(6,i,j,xi,xj) = gl%J2X2_PCSAFTQ(6,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T * gl%z3_PCSAFT(1)**2 * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(6,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(6,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(4,xj)) &
                                & + T * gl%z3_PCSAFT(1) * (n - 2.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j) * eps_ij) &
                                & * (gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0) * (T * gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                                & - T * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & - T * gl%ab_PCSAFTQ(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0)) / T**2
                        end do
                    end do
                end do
            end do  
        end do 
    end if
  
    
!DEC$ END IF
end subroutine J2X2DERIVS    

subroutine J2X2DERIVS_D(gl,T,DENS,GETDERI)
    
    ! Ina Schülling, September 2017

    ! J2: integral of the perturbation theory
    ! defined by eq. eq 12 in Gross 2005:
    ! J2 = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2x1_PCSAFT (module variables)
    integer :: i,j,n,xi,xj
    double precision :: eps_ij, sum1, sum2, sum3, sum4
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    eps_ij = 0.d0
    gl%J2X2_PCSAFTD = 0.d0
    
    ! III. X1 DERIVATIVES of I
    !calculate the derivatives of J2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        do xj = 1, gl%ncomp
            do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(1,i,j,xi,xj) = gl%J2X2_PCSAFTD(1,i,j,xi,xj) + n * (n - 1.d0) * gl%z3_PCSAFT(1)**(n - 2.d0) * gl%z3x1_PCSAFT(1,xj) * gl%z3x1_PCSAFT(1,xi) * (gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij / T) 
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
    do xj = 1, gl%ncomp
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(2,i,j,xi,xj) = gl%J2X2_PCSAFTD(2,i,j,xi,xj) +  n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) * (gl%z3_PCSAFT(1)*(gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) &
                                & + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0)) / T
                        end do
                    end do
                end do
            end do
        end do 
    end if  
   
    
    ! 3: 2ND DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
    do xj = 1, gl%ncomp  
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(3,i,j,xi,xj) = gl%J2X2_PCSAFTD(3,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (2.d0 * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(2,xj) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (n - 2.d0) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0)) / T
                        end do
                    end do
                end do
            end do 
        end do 
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
    do xj = 1, gl%ncomp   
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(4,i,j,xi,xj) = gl%J2X2_PCSAFTD(4,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) + T * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0) &
                                & * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) - T * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) / T**2
                        end do
                    end do
                end do
            end do  
        end do 
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
    do xj = 1, gl%ncomp  
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(5,i,j,xi,xj) = gl%J2X2_PCSAFTD(5,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T**2 * gl%z3_PCSAFT(1)**2 * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(5,xj) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(5,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T**2 * gl%z3_PCSAFT(1) * (n - 2.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T**2 * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & - 2.d0 * T**2 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & - 2.d0 * T**2 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0) &
                                & + 2.d0 * T**2 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) / T**3
                        end do
                    end do
                end do
            end do
        end do 
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
    do xj = 1, gl%ncomp   
        do xi = 1, gl%ncomp   
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                        do n = 0, 4
                            gl%J2X2_PCSAFTD(6,i,j,xi,xj) = gl%J2X2_PCSAFTD(6,i,j,xi,xj) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (T * gl%z3_PCSAFT(1)**2 * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(6,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(6,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(4,xj)) &
                                & + T * gl%z3_PCSAFT(1) * (n - 2.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                                & * (gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & + T * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                                & - T * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                & - T * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0)) / T**2
                        end do
                    end do
                end do
            end do  
        end do 
    end if
  
    
!DEC$ END IF
end subroutine J2X2DERIVS_D    

subroutine J3X2DERIVS(gl,T,DENS,GETDERI)
! Ina Schülling, September 2017

    ! J3: integral of the perturbation theory
    ! defined by eq. eq 13 in Gross 2005:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS, T
    !output: j3x1_PCSAFT (module variables)
    integer :: i,j,k,n,xi,xj

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X1 DERIVATIVES of J3
    !calculate the derivatives of J3
    gl%J3X2_PCSAFTQ = 0.d0

    if (GETDERI(1) .eq. 1) then
    do xj = 1, gl%ncomp    
        do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(1,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(1,i,j,k,xi,xj)  + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 2.d0) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 1.d0)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xj = 1, gl%ncomp
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(2,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(2,i,j,k,xi,xj) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (gl%z3_PCSAFT(1) &
                                    & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xj = 1, gl%ncomp    
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(3,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(3,i,j,k,xi,xj) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0)  &
                                    & * (2.d0 * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(2,xj) &
                                    & + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (n - 2.d0) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do   
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xj = 1, gl%ncomp      
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(4,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(4,i,j,k,xi,xj) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0)   &
                                    & * (gl%z3_PCSAFT(1) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xj = 1, gl%ncomp       
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(5,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(5,i,j,k,xi,xj) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) &
                                    & * (gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(5,xj) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(5,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(1) * (n - 2.d0) * (2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj) &
                                    & + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xj = 1, gl%ncomp       
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTQ(6,i,j,k,xi,xj) = gl%J3X2_PCSAFTQ(6,i,j,k,xi,xj) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (gl%z3_PCSAFT(1)**2 &
                                    & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(6,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(6,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(4,xj)) &
                                    & + gl%z3_PCSAFT(1) * (n - 2.d0) * (gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj) &
                                    & + gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
  
    !DEC$ END IF
end subroutine J3X2DERIVS
   
subroutine J3X2DERIVS_D(gl,T,DENS,GETDERI)
! Ina Schülling, September 2017

    ! J3: integral of the perturbation theory
    ! defined by eq. eq 13 in Gross 2005:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS, T
    !output: j3x1_PCSAFT (module variables)
    integer :: i,j,k,n,xi,xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%J3X2_PCSAFTD = 0.d0
   
    ! III. X1 DERIVATIVES of J3
    !calculate the derivatives of J3

    if (GETDERI(1) .eq. 1) then
    do xj = 1, gl%ncomp    
        do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(1,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(1,i,j,k,xi,xj)  + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 2.d0) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 1.d0)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xj = 1, gl%ncomp
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(2,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(2,i,j,k,xi,xj) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0) * (gl%z3_PCSAFT(1) &
                                    & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xj = 1, gl%ncomp    
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(3,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(3,i,j,k,xi,xj) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0)  &
                                    & * (2.d0 * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(2,xj) &
                                    & + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (n - 2.d0) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do   
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xj = 1, gl%ncomp      
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(4,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(4,i,j,k,xi,xj) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3.d0) * (n - 1.d0)   &
                                    & * (gl%z3_PCSAFT(1) * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xj = 1, gl%ncomp       
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(5,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(5,i,j,k,xi,xj) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) &
                                    & * (gl%z3_PCSAFT(1)**2 * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(5,xj) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3x1_PCSAFT(5,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(1) * (n - 2.d0) * (2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj) &
                                    & + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj)) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xj = 1, gl%ncomp       
            do xi = 1, gl%ncomp    
                do k = gl%n_start, gl%n_end
                    do j = gl%n_start, gl%n_end
                        do i = gl%n_start, gl%n_end
                            do n = 0, 4
                                gl%J3X2_PCSAFTD(6,i,j,k,xi,xj) = gl%J3X2_PCSAFTD(6,i,j,k,xi,xj) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4.d0) * (n - 1.d0) * (gl%z3_PCSAFT(1)**2 &
                                    & * (gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(6,xj) + gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3x1_PCSAFT(6,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(4,xj)) &
                                    & + gl%z3_PCSAFT(1) * (n - 2.d0) * (gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(2,xj) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) * gl%z3x1_PCSAFT(1,xj) &
                                    & + gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(4,xj) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(4,xi) * gl%z3x1_PCSAFT(1,xj)) &
                                    & + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * gl%z3x1_PCSAFT(1,xj) * (n - 3.d0) * (n - 2.d0))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    
  
!DEC$ END IF
end subroutine J3X2DERIVS_D
    
subroutine J3X1DERIVS(gl,T,DENS,GETDERI)
! Ina Schülling, September 2017

    ! J3: integral of the perturbation theory
    ! defined by eq. eq 13 in Gross 2005:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS, T
    !output: j3x1_PCSAFT (module variables)
    integer :: i,j,k,n,xi

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. X1 DERIVATIVES of J3
    !calculate the derivatives of J3
    
    gl%J3X1_PCSAFTQ = 0.d0

    if (GETDERI(1) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(1,i,j,k,xi) = gl%J3X1_PCSAFTQ(1,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**n * gl%z3x1_PCSAFT(1,xi) / gl%z3_PCSAFT(1)
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(2,i,j,k,xi) = gl%J3X1_PCSAFTQ(2,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n**2 * gl%z3_PCSAFT(1)**(n - 1) * gl%z3x1_PCSAFT(1,xi)
                        end do
                    end do
                end do
            end do
        end do    
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(3,i,j,k,xi) = gl%J3X1_PCSAFTQ(3,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*gl%z3_PCSAFT(2)*(2.d0*gl%z3_PCSAFT(1)* &
                             &  gl%z3x1_PCSAFT(2,xi)*(n - 1.d0) + gl%z3_PCSAFT(2)* gl%z3x1_PCSAFT(1,xi)*(-n + 1.d0) + &
                             & gl%z3_PCSAFT(2)* gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)**2)
                        end do
                    end do
                end do
            end do
        end do    
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(4,i,j,k,xi) = gl%J3X1_PCSAFTQ(4,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 2) * (gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(5,i,j,k,xi) = gl%J3X1_PCSAFTQ(5,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) + gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(6,i,j,k,xi) = gl%J3X1_PCSAFTQ(6,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(6,xi) &
                                & + gl%z3_PCSAFT(1) * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) &
                                & - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n**2 - 3.d0 * n + 2.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(7,i,j,k,xi) = gl%J3X1_PCSAFTQ(7,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) &
                                & * (n * gl%z3_PCSAFT(2) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0)) &
                                & + n * gl%z3x1_PCSAFT(1,xi) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0))  &
                                & + gl%z3_PCSAFT(1)**3 * gl%z3x1_PCSAFT(7,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(6,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(6,xi) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(4,xi))  &
                                & + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * (- n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - 2.d0 * n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - 2.d0 * n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + 2.d0 * gl%z3_PCSAFT(6) *  gl%z3x1_PCSAFT(1,xi) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) &
                                & + gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n - 1.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2) &
                                & + 2.d0 * gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n**2 - 3.d0 * n + 2.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2)  &
                                & - gl%z3_PCSAFT(2) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0))&
                                & - gl%z3x1_PCSAFT(1,xi) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(8,i,j,k,xi) = gl%J3X1_PCSAFTQ(8,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * gl%z3_PCSAFT(2)**2 * (3.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n**2 - 3.d0 * n + 2.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n**2 - 3.d0 * n + 2.d0) + 2.d0 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J3 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(9,i,j,k,xi) = gl%J3X1_PCSAFTQ(9,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * (gl%z3_PCSAFT(1)**3 * gl%z3x1_PCSAFT(9,xi)  &
                                & + 3.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(5,xi) + n * gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(5,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(4,xi))  &
                                & + 3.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * (n**2 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - 3.d0 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) -  n * gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) + 2.d0 * gl%z3_PCSAFT(4)**3 * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0) &
                                & + gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(9) + 3.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(5) * (n - 1.d0) + gl%z3_PCSAFT(4)**3 * (n**2 - 3.d0 * n + 2.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTQ(10,i,j,k,xi) = gl%J3X1_PCSAFTQ(10,i,j,k,xi) + gl%c_PCSAFTQ(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * (2.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(2,xi) &
                                & + n * gl%z3x1_PCSAFT(6,xi) * gl%z3_PCSAFT(1) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(2,xi) - gl%z3x1_PCSAFT(6,xi) * gl%z3_PCSAFT(2)) &
                                & + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (2.d0 * n**2 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + n**2 * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) &
                                & - 6.d0 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - 2.d0 * n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) &
                                & - 3.d0 * n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) + 4.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + 2.d0 * gl%z3x1_PCSAFT(4,xi) *gl%z3_PCSAFT(2)) &
                                & + 2.d0 * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(6) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * (n**2 - 3.d0 * n + 2.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
!DEC$ END IF
end subroutine J3X1DERIVS

subroutine J2X1DERIVS_D(gl,T, DENS, GETDERI)

! I.Schuelling 07/17
    
    ! J2: integral of the perturbation theory
    ! defined by eq. eq 10 in Gross 2006:
    ! J2_D = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2_PCSAFT (module variables)
    integer :: i,j,n,xi
    double precision :: eps_ij, sum1, sum2, sum3, sum4
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    eps_ij = 0.d0
    gl%J2X1_PCSAFTD = 0.d0
    
    ! III. X1 DERIVATIVES of I
    !calculate the derivatives of J2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
    do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTD(1,i,j,xi) = gl%J2X1_PCSAFTD(1,i,j,xi) +  n*gl%z3_PCSAFT(1)**(n - 1.d0)*gl%z3x1_PCSAFT(1,xi)*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T) 
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
      do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTD(2,i,j,xi) = gl%J2X1_PCSAFTD(2,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 2.d0)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                            & *(gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(2,xi) + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0))/T
                        !J2X1_PCSAFTD(2,i,j,xi) = J2X1_PCSAFTD(2,i,j,xi) + n**2 *z3_PCSAFT(1)**(n - 1)* z3x1_PCSAFT(1,xi)*(ab_PCSAFTD(1,n,i,j) + ab_PCSAFTD(2,n,i,j)*eps_ij/T)

                    end do
                end do
            end do
        end do
    end if  
   
    
    ! 3: 2ND DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
      do xi = 1, gl%ncomp   
        do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTD(3,i,j,xi) = gl%J2X1_PCSAFTD(3,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * gl%z3_PCSAFT(2) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij)  &
                             & * (2.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n - 1.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0)**2) / T  
                    end do
                end do
            end do
        end do  
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTD(4,i,j,xi) = gl%J2X1_PCSAFTD(4,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 2) * (T * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(4,xi) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij)  &
                             & + T * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * (n - 1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                             & - T * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(1,xi)) / T**2
                    end do
                end do
            end do
        end do  
       
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
       do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                   
                        gl%J2X1_PCSAFTD(5,i,j,xi) = gl%J2X1_PCSAFTD(5,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (T * (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) * (T * gl%ab_PCSAFTD(1,n,i,j) + &
                             & gl%ab_PCSAFTD(2,n,i,j) * eps_ij) + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) * (2.d0 * n * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(4,xi) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) + 2.d0 * &
                             & T * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) - 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * &
                             & eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(4,xi) * T + 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(1,xi) * T) + gl%z3x1_PCSAFT(1,xi) * (T**2 * n * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * &
                             & eps_ij) + T**2 * n * gl%z3_PCSAFT(4)**2 * (n - 1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) - 2.d0 * T * &
                             & gl%ab_PCSAFTD(2,n,i,j) * eps_ij * n * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * T + 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * T**2 )) / T**3

                    end do
                end do
            end do
        end do  
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        gl%J2X1_PCSAFTD(6,i,j,xi) = gl%J2X1_PCSAFTD(6,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 3.d0) * (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(6,xi) * (T * gl%ab_PCSAFTD(1,n,i,j) + &
                             & gl%ab_PCSAFTD(2,n,i,j) * eps_ij) + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) * (n * gl%z3_PCSAFT(4) * &
                             & gl%z3x1_PCSAFT(2,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + n * gl%z3x1_PCSAFT(4,xi) * &
                             & gl%z3_PCSAFT(2) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - &
                             & gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) + T * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (T * &
                             & gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) * (n**2 - 3.d0 * n + 2.d0) - gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * &
                             & gl%z3x1_PCSAFT(2,xi) * T + gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * T * (-n + 1.d0))/ &
                             & T**2
                        !J2X1_PCSAFTD(6,i,j,xi) = J2X1_PCSAFTD(6,i,j,xi) + n*(ab_PCSAFTD(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(6,xi)- ab_PCSAFTD(2,n,i,j)*eps_ij*n*z3_PCSAFT(1)**(n - 2)*z3_PCSAFT(4)*z3x1_PCSAFT(1,xi)/T)

                    end do
                end do
            end do
        end do  
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTD(7,i,j,xi) = J2X1_PCSAFTD(7,i,j,xi) + n*(ab_PCSAFTD(1,n,i,j)*z3_PCSAFT(1)**(n - 1)*z3x1_PCSAFT(7,xi)+ 2*ab_PCSAFTD(2,n,i,j)*eps_ij*n*z3_PCSAFT(1)**(n - 3)*z3_PCSAFT(4)**2*z3x1_PCSAFT(1,xi)/T)
                        
                        gl%J2X1_PCSAFTD(7,i,j,xi) = gl%J2X1_PCSAFTD(7,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * (T * n * (gl%z3_PCSAFT(2) *&
                             & (T * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                             & + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                             & * (2 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) &
                             & + 2.d0 * T * gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij)  &
                             & - 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(4,xi) * T + 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * T)  &
                             
                             & + gl%z3x1_PCSAFT(1,xi) * (T * gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij)  &
                             & + T * gl%z3_PCSAFT(1) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) * (2.d0 * n * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) - 2.d0 * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) - gl%z3_PCSAFT(5) * gl%z3_PCSAFT(2)) &
                             & + 2.d0 * T * gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0) * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                             & - 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(6) * T + 2.d0 * gl%ab_PCSAFTD(2,n,i,j) * eps_ij * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * T)) &
                             
                             & + T*(T*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(7,xi)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)  &
                             & + T*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(2.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi) &
                             & + 2.d0*n*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(4,xi) - 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(6,xi)  &
                             & - gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(2,xi) - gl%z3_PCSAFT(7)*gl%z3x1_PCSAFT(1,xi) - gl%z3x1_PCSAFT(5,xi)*gl%z3_PCSAFT(2)  &
                             & - 2.d0*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(4,xi)) + 2.d0*T*gl%z3_PCSAFT(1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                             & * (-n*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(2,xi) - 2.d0*n*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) - 2.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2)  &
                             & + gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(2,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2)+ gl%z3_PCSAFT(5)*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi))   &
                             
                             & + 6.d0*T*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                             & - 2.d0*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(6,xi)*T + 2.d0*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2 &
                             & * (gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) + gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi)+ gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2))*T  &
                             & - 4.d0*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*T) &
                             & + gl%z3_PCSAFT(1)*gl%z3x1_PCSAFT(2,xi)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                             & + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1.d0)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) - 2.d0*T*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*T &
                             & + 2.d0*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*T**2) &
                             & + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTD(1,n,i,j)+ gl%ab_PCSAFTD(2,n,i,j)*eps_ij)  &
                             & + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1.d0)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                             &  - 2.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2)) &
                             & /T**3

                    end do
                end do
            end do
        end do  
    end if
    
    ! 8: 3RD DERIVATIVE OF J2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                        !J2X1_PCSAFTD(8,i,j,xi) = J2X1_PCSAFTD(8,i,j,xi) + (n**2*z3_PCSAFT(1)**n*z3x1_PCSAFT(1,xi)*(ab_PCSAFTD(1,n,i,j)+ ab_PCSAFTD(2,n,i,j)*eps_ij/T)*(n**2*z3_PCSAFT(1) - 3*n*z3_PCSAFT(1) + 2*z3_PCSAFT(1))/z3_PCSAFT(1)**2 &
                        !    & + n*z3_PCSAFT(1)**n*(ab_PCSAFTD(1,n,i,j)+ ab_PCSAFTD(2,n,i,j)*eps_ij/T)*(n**2*z3x1_PCSAFT(1,xi) - 3*n*z3x1_PCSAFT(1,xi) + 2*z3x1_PCSAFT(1,xi))/z3_PCSAFT(1) &
                        !    & - n*z3_PCSAFT(1)**n*z3x1_PCSAFT(1,xi)*(ab_PCSAFTD(1,n,i,j)+ ab_PCSAFTD(2,n,i,j)*eps_ij/T)*(n**2*z3_PCSAFT(1) - 3*n*z3_PCSAFT(1) + 2*z3_PCSAFT(1))/z3_PCSAFT(1)**2)
                        gl%J2X1_PCSAFTD(8,i,j,xi) = gl%J2X1_PCSAFTD(8,i,j,xi) + n * gl%z3_PCSAFT(1)**(n - 4.d0) * gl%z3_PCSAFT(2)**2 * (T * gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) * eps_ij) &
                        & * (3.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n**2 - 3.d0 * n + 2.d0) + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n**2 - 3.d0 * n + 2.d0) + 2.d0 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (-n**2 + 3.d0 * n - 2.d0)) / T
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J2 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
       do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                    
                        gl%J2X1_PCSAFTD(9,i,j,xi) = gl%J2X1_PCSAFTD(9,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 4)*(T*(T**2*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(9,xi)*(T* &
                             & gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + T**2*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(3.d0*n* &
                             & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(5,xi) + 3.d0*n*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - 3.d0*gl%z3_PCSAFT(4)* &
                             & gl%z3x1_PCSAFT(5,xi) - 3.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(9)*gl%z3x1_PCSAFT(1,xi))  &
                             
                             & + 3.d0*T**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(n**2*gl%z3_PCSAFT(4)* &
                             & gl%z3x1_PCSAFT(4,xi) - 3.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) - 2.d0*n*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi) &
                             & + 2.d0*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2.d0*gl%z3_PCSAFT(5)*gl%z3x1_PCSAFT(1,xi))  &
                             
                             & + 3.d0*T**2*gl%z3_PCSAFT(4)**3*gl%z3x1_PCSAFT(1,xi)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(-n**2 + 3.d0*n - 2.d0) &
                             & - 3.d0 *T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(5,xi) + 3.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*(-2.d0 &
                             & *n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + 2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(5)* &
                             & gl%z3x1_PCSAFT(1,xi)) &
                             & + 6.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0 &
                             & ) + 6.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3*gl%z3x1_PCSAFT(4,xi) - 6.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(1,xi)) + gl%z3x1_PCSAFT(1,xi)*(T**3*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9)*(T &
                             & *gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + 3.d0*T**3*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1.d0)*( &
                             & T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + T**3*n*gl%z3_PCSAFT(4)**3*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(n &
                             & **2 - 3.d0*n + 2) - 3.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 3.d0*T**3* &
                             & gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(-n + 1.d0) + 6.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1) &
                             & **2*gl%z3_PCSAFT(4) - 6.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3))/T**4

                    end do
                end do
            end do
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
      do xi = 1, gl%ncomp   
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                    do n = 0, 4
                    gl%J2X1_PCSAFTD(10,i,j,xi) = gl%J2X1_PCSAFTD(10,i,j,xi) + n*gl%z3_PCSAFT(1)**(n - 4.d0)*(2.d0*T*gl%z3_PCSAFT(1)**2*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                         & *(n*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(2,xi) + n*gl%z3x1_PCSAFT(6,xi)*gl%z3_PCSAFT(2) - &
                         & gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(2,xi) - gl%z3x1_PCSAFT(6,xi)*gl%z3_PCSAFT(2)) + T*gl%z3_PCSAFT(1) &
                         & *gl%z3_PCSAFT(2)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(2.d0*n**2*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) + &
                         & n**2*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2) - 6.d0*n*gl%z3_PCSAFT(4)*gl%z3x1_PCSAFT(2,xi) - 2.d0*n* &
                         & gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) - 3.d0*n*gl%z3x1_PCSAFT(4,xi)*gl%z3_PCSAFT(2) + 4.d0*gl%z3_PCSAFT(4)* &
                         & gl%z3x1_PCSAFT(2,xi) + 2.d0*gl%z3_PCSAFT(6)*gl%z3x1_PCSAFT(1,xi) + 2.d0*gl%z3x1_PCSAFT(4,xi)* &
                         & gl%z3_PCSAFT(2)) + 2.d0*T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)**2*gl%z3x1_PCSAFT(1,xi)*(T*gl%ab_PCSAFTD(1,n,i,j) + &
                         & gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(-n**2 + 3.d0*n - 2.d0) + 2.d0*T*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(2)* &
                         & gl%z3x1_PCSAFT(2,xi)*(-n + 1.d0) + T*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(2)**2* &
                         & gl%z3x1_PCSAFT(1,xi)*(n - 1.d0) + gl%z3_PCSAFT(2)*gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)*(2.d0*T*gl%z3_PCSAFT(1)* &
                         & gl%z3_PCSAFT(6)*(n - 1.d0)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)* &
                         & (T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(n**2 - 3.d0*n + 2.d0) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)* &
                         & gl%z3_PCSAFT(2)*T*(-n + 1.d0)))/T**2

                    end do
                end do
            end do
        end do
    end if
    
!DEC$ END IF
end subroutine J2X1DERIVS_D

subroutine J3X1DERIVS_D(gl,T,DENS,GETDERI)

! I. Schuelling 07/2017
    
    ! J3_D: integral over the reference-fluid pair-correlation function and over three-body correlation functions
    ! defined by eq. eq 11 in Gross 2006:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T
    





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS, T
    !output: j3x1_PCSAFT (module variables)
    integer :: i,j,k,n,xi
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%J3X1_PCSAFTD = 0.d0
    
    
    ! III. X1 DERIVATIVES of J3
    !calculate the derivatives of J3

    if (GETDERI(1) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(1,i,j,k,xi) = gl%J3X1_PCSAFTD(1,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**n * gl%z3x1_PCSAFT(1,xi) / gl%z3_PCSAFT(1)
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(2,i,j,k,xi) = gl%J3X1_PCSAFTD(2,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n**2 * gl%z3_PCSAFT(1)**(n - 1) * gl%z3x1_PCSAFT(1,xi)
                        end do
                    end do
                end do
            end do
        end do    
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(3,i,j,k,xi) = gl%J3X1_PCSAFTD(3,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*gl%z3_PCSAFT(2)*(2.d0*gl%z3_PCSAFT(1)* &
                             &  gl%z3x1_PCSAFT(2,xi)*(n - 1.d0) + gl%z3_PCSAFT(2)* gl%z3x1_PCSAFT(1,xi)*(-n + 1.d0) + &
                             & gl%z3_PCSAFT(2)* gl%z3x1_PCSAFT(1,xi)*(n - 1.d0)**2)
                        end do
                    end do
                end do
            end do
        end do    
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(4,i,j,k,xi) = gl%J3X1_PCSAFTD(4,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 2) * (gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(5,i,j,k,xi) = gl%J3X1_PCSAFTD(5,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0) + gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(6,i,j,k,xi) = gl%J3X1_PCSAFTD(6,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 3) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(6,xi) &
                                & + gl%z3_PCSAFT(1) * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) &
                                & - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n**2 - 3.d0 * n + 2.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(7,i,j,k,xi) = gl%J3X1_PCSAFTD(7,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) &
                                & * (n * gl%z3_PCSAFT(2) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0)) &
                                & + n * gl%z3x1_PCSAFT(1,xi) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0))  &
                                & + gl%z3_PCSAFT(1)**3 * gl%z3x1_PCSAFT(7,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(6,xi) + n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(6,xi) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(4,xi))  &
                                & + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * (- n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - 2.d0 * n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) - 2.d0 * n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) + gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + 2.d0 * gl%z3_PCSAFT(6) *  gl%z3x1_PCSAFT(1,xi) + 2.d0 * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2)) &
                                & + gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n - 1.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2) &
                                & + 2.d0 * gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n**2 - 3.d0 * n + 2.d0) * (n * gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2)  &
                                & - gl%z3_PCSAFT(2) * (gl%z3_PCSAFT(1)**2 * gl%z3x1_PCSAFT(5,xi) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3x1_PCSAFT(1,xi) * (-n + 1.d0))&
                                & - gl%z3x1_PCSAFT(1,xi) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(7) + 2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(6) * (n - 1.d0) + gl%z3_PCSAFT(4)**2 * gl%z3_PCSAFT(2) * (-n + 1.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(8,i,j,k,xi) = gl%J3X1_PCSAFTD(8,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * gl%z3_PCSAFT(2)**2 * (3.d0 * gl%z3_PCSAFT(1) * gl%z3x1_PCSAFT(2,xi) * (n**2 - 3.d0 * n + 2.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (n**2 - 3.d0 * n + 2.d0) + 2.d0 * gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J3 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(9,i,j,k,xi) = gl%J3X1_PCSAFTD(9,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * (gl%z3_PCSAFT(1)**3 * gl%z3x1_PCSAFT(9,xi)  &
                                & + 3.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(5,xi) + n * gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(4,xi) - gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(5,xi) - gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(4,xi))  &
                                & + 3.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * (n**2 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) - 3.d0 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) -  n * gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(4,xi) + gl%z3_PCSAFT(5) * gl%z3x1_PCSAFT(1,xi)) + 2.d0 * gl%z3_PCSAFT(4)**3 * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0) &
                                & + gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (gl%z3_PCSAFT(1)**2 * gl%z3_PCSAFT(9) + 3.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(5) * (n - 1.d0) + gl%z3_PCSAFT(4)**3 * (n**2 - 3.d0 * n + 2.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
        do xi = 1, gl%ncomp    
            do k = gl%n_start, gl%n_end
                do j = gl%n_start, gl%n_end
                    do i = gl%n_start, gl%n_end
                        do n = 0, 4
                            gl%J3X1_PCSAFTD(10,i,j,k,xi) = gl%J3X1_PCSAFTD(10,i,j,k,xi) + gl%c_PCSAFTD(n,i,j,k) * n * gl%z3_PCSAFT(1)**(n - 4) * (2.d0 * gl%z3_PCSAFT(1)**2 * (n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(2,xi) &
                                & + n * gl%z3x1_PCSAFT(6,xi) * gl%z3_PCSAFT(1) - gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(2,xi) - gl%z3x1_PCSAFT(6,xi) * gl%z3_PCSAFT(2)) &
                                & + gl%z3_PCSAFT(1) * gl%z3_PCSAFT(2) * (2.d0 * n**2 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) + n**2 * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) &
                                & - 6.d0 * n * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) - 2.d0 * n * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) &
                                & - 3.d0 * n * gl%z3x1_PCSAFT(4,xi) * gl%z3_PCSAFT(2) + 4.d0 * gl%z3_PCSAFT(4) * gl%z3x1_PCSAFT(2,xi) &
                                & + 2.d0 * gl%z3_PCSAFT(6) * gl%z3x1_PCSAFT(1,xi) + 2.d0 * gl%z3x1_PCSAFT(4,xi) *gl%z3_PCSAFT(2)) &
                                & + 2.d0 * gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2)**2 * gl%z3x1_PCSAFT(1,xi) * (- n**2 + 3.d0 * n - 2.d0) &
                                & + gl%z3_PCSAFT(2) * gl%z3x1_PCSAFT(1,xi) * (n - 1.d0) * (2.d0 * gl%z3_PCSAFT(1) * gl%z3_PCSAFT(6) * (n - 1.d0) &
                                & + gl%z3_PCSAFT(4) * gl%z3_PCSAFT(2) * (n**2 - 3.d0 * n + 2.d0)))
                        end do
                    end do
                end do
            end do
        end do
    end if
    
!DEC$ END IF
end subroutine J3X1DERIVS_D




    end module pc_saft_JX_derivs_module
