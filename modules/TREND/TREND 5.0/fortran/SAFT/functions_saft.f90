    !*************************************************************************************
    !				TREND Version 4.0
    !		   Thermodynamic Reference & Engineering Data
    !
    !- software for the calculation of thermodynamic and other properties -
    !
    !Copyright (C) 2019,  Prof. Dr.-Ing. R.Span
    !                     Lehrstuhl fuer Thermodynamik
    !                     Ruhr-Universitaet Bochum
    !                     Universitaetsstr. 150
    !                     D-44892 Bochum
    !
    !Cite as: Span, R.; Beckmüller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.; 
	!          Jäger, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2019): 	
    !          TREND. Thermodynamic Reference and Engineering Data 4.0. 
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

subroutine test_ar_x3(T,D,delta,dx3)
    ! The alternative routine doesn't change the last molfraction to 1-Sum(other molfractions)
    
    use module_nderivs
    use module_fluid_parameters
    
    implicit none
    
    double precision :: T, D
    double precision :: delta, sum_xi
    double precision, dimension(30) :: molfractions_original
    double precision, dimension(ncomp,ncomp,ncomp) :: numerical_derivative, analytical_derivative, func_plus_delta, func_minus_delta
    double precision, dimension(nx2derivs,ncomp,ncomp) :: x2
    double precision, dimension(nx3derivs,ncomp,ncomp,ncomp) :: x3
    double precision, dimension(3,30,30,30) :: dx3
    integer, dimension(nderivs) :: getvector
    integer :: i, xi, j, xj, xk
    
    !only works when ncomp > 1
    if (ncomp .EQ. 1) then
        RETURN
    end if
    
    !save the original molfractions
    molfractions_original = molfractions
    
    ! The derivatives to xi
    do i = 1, 3
        
        !re-initialize the getvector in every loop
        getvector = 0
        getvector(i) = 1
        
        ! in x3 derivatives, the number 3 is the number 4 of the x2 derivatives
        j = i ! for 1 and 2 it is the same
        if (i .EQ. 3) then
            getvector(i) = 0
            getvector(4) = 1
            j = 4
        end if
        
        !Analytical Derivative number i
        x3 = 0.d0
        call ARX3DERIVS(T, D, getvector, x3)
        analytical_derivative = x3(i,1:ncomp,1:ncomp,1:ncomp)
        
        !Numerical Derivative number i
        !f(x_i + delta_x_i)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                do xk = 1, ncomp
                    !add delta_x to the xi
                    molfractions(xk) = molfractions(xk) + delta
            
                    ! calculate the function with the changed xj (2nd component which is derived to)
                    x2 = 0.d0
                    call ARX2DERIVS(T,D,getvector,x2)
                    func_plus_delta(xi,xj,xk) = x2(j,xi,xj) ! the number xi and xj gets derived by xk
            
                    molfractions = molfractions_original
                end do
            end do
        end do
        
        !f(x_i - delta_x_i)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                do xk = 1, ncomp
                    !subtract delta_x from the xi
                    molfractions(xk) = molfractions(xk) - delta
            
                    ! calculate the function with the changed xi
                    x2 = 0.d0
                    call ARX2DERIVS(T,D,getvector,x2)
                    func_minus_delta(xi,xj,xk) = x2(j,xi,xj)
            
                    molfractions = molfractions_original
                end do
            end do
        end do
        
        !Numerical derivative = (func_plus_delta - func_minus_delta) / (2*delta)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                do xk = 1, ncomp
                    numerical_derivative(xi,xj,xk) = (func_plus_delta(xi,xj,xk) - func_minus_delta(xi,xj,xk))/(2.d0*delta)
                end do
            end do
        end do
        
        !calculate the difference dx1 between numerical and analytical derivative
        ! assignment with dx3(i,:,:,:) = ... not possible! because the array numbering gets mixed up as dx3 has dimensions (i,30,30,30) and the others have (ncomp,ncomp,ncomp)
        do xi = 1, ncomp
            do xj = 1, ncomp
                do xk = 1, ncomp
                    dx3(i,xi,xj,xk) = analytical_derivative(xi,xj,xk) - numerical_derivative(xi,xj,xk)
                end do
            end do
        end do
        
    end do
    
end subroutine test_ar_x3
    
subroutine test_ar_x2(T,D,delta,dx2)
    ! The alternative routine doesn't change the last molfraction to 1-Sum(other molfractions)
    
    use module_nderivs
    use module_fluid_parameters
    
    implicit none
    
    double precision :: T, D
    double precision :: delta, sum_xi
    double precision, dimension(30) :: molfractions_original
    double precision, dimension(ncomp,ncomp) :: numerical_derivative, analytical_derivative, func_plus_delta, func_minus_delta
    double precision, dimension(nx1derivs,ncomp) :: x1
    double precision, dimension(nx2derivs,ncomp,ncomp) :: x2
    double precision, dimension(6,30,30) :: dx2
    integer, dimension(nderivs) :: getvector
    double precision, dimension(nderivs) ::  resultvector
    integer :: i, xi, j, xj
    
    !only works when ncomp > 1
    if (ncomp .EQ. 1) then
        RETURN
    end if
    
    !save the original molfractions
    molfractions_original = molfractions
    
    ! The derivatives to xi
    do i = 1, 6
        
        !re-initialize the getvector in every loop
        getvector = 0
        getvector(i) = 1
        
        !Analytical Derivative number i
        x2 = 0.d0
        call ARX2DERIVS(T, D, getvector, x2)
        analytical_derivative = x2(i,1:ncomp,1:ncomp)
        
        !Numerical Derivative number i
        !f(x_i + delta_x_i)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                !add delta_x to the xi
                molfractions(xj) = molfractions(xj) + delta
            
                ! calculate the function with the changed xj (second component which is derived to)
                x1 = 0.d0
                call ARX1DERIVS(T,D,getvector,x1)
                func_plus_delta(xi,xj) = x1(i,xi) ! xi in x1 here as the first component xi derived to xj again
            
                molfractions = molfractions_original
            end do
        end do
        
        !f(x_i - delta_x_i)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                !subtract delta_x from the xi
                molfractions(xj) = molfractions(xj) - delta
            
                ! calculate the function with the changed xi
                x1 = 0.d0
                call ARX1DERIVS(T,D,getvector,x1)
                func_minus_delta(xi,xj) = x1(i,xi)! not sure if this is correct: xi in x1 here as the first component xi derived to xj again
            
                molfractions = molfractions_original
            end do
        end do
        
        !Numerical derivative = (func_plus_delta - func_minus_delta) / (2*delta)
        do xi = 1, ncomp
            do xj = 1, ncomp ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                numerical_derivative(xi,xj) = (func_plus_delta(xi,xj) - func_minus_delta(xi,xj))/(2.d0*delta)
            end do
        end do
        
        !calculate the difference dx1 between numerical and analytical derivative
        ! assignment with dx2(i,:,:) = ... not possible! because the array numbering gets mixed up as dx2 has dimensions (i,30,30) and the others have (ncomp,ncomp)
        do xi = 1, ncomp
            do xj = 1, ncomp
                dx2(i,xi,xj) = analytical_derivative(xi,xj) - numerical_derivative(xi,xj)
            end do
        end do
        
    end do
    
    end subroutine test_ar_x2
    
subroutine test_adisp_x2(T,D_in,delta,dx2)
    ! The alternative routine doesn't change the last molfraction to 1-Sum(other molfractions)
    
    use module_nderivs
    use module_fluid_parameters
    use module_general_eos_parameters !for N_A
    use module_PCSAFT
    
    implicit none
    
    double precision :: T, D_in, D
    double precision :: delta, sum_xi
    double precision, dimension(30) :: molfractions_original
    double precision, dimension(ncomp,ncomp) :: numerical_derivative, analytical_derivative, func_plus_delta, func_minus_delta
    double precision, dimension(nx1derivs,ncomp) :: x1
    double precision, dimension(nx2derivs,ncomp,ncomp) :: x2
    double precision, dimension(6,30,30) :: dx2
    integer, dimension(nderivs) :: getvector, getprevious
    double precision, dimension(nderivs) ::  resultvector
    integer :: i, xi, j, xj, index
    
    !only works when ncomp > 1
    if (ncomp .EQ. 1) then
        RETURN
    end if
    
    !save the original molfractions
    molfractions_original = molfractions
    
    ! Change the units of the density from mol/m^3 to 1/Angstrom
    D = D_in*N_A/1.d30
    
    ! The derivatives to xi
    do i = 1, 6
        
        !re-initialize the getvector in every loop
        getvector = 0
        getvector(i) = 1
        
        !Analytical Derivative number i
        x2 = 0.d0
        
        ! for ahcx2
        ! get mmean, a_hs, g_ii, ahs_x, gii_x, ahs_xx, gii_xx
        ! get d_i, zeta_0, zeta_1, zeta_2, zeta_3, zeta_0_x, zeta_1_x, zeta_2_x, zeta_3_x
        call add_previous_derivs(getvector,getprevious)
        call calculate_PCSAFT_functionparts_Trho(T,D,getprevious)
        call AHSDERIVS(getvector)
        ! X derivs
        call calculate_PCSAFT_functionparts_x1(T,D,getprevious)
        call AHSX1DERIVS(getvector)
        ! XX derivs
        call calculate_PCSAFT_functionparts_x2(T,D,getprevious)
        call AHSX2DERIVS(getvector)
        
        call ADISPX2DERIVS(D, getvector)
        analytical_derivative = adispx2_PCSAFT(i,1:ncomp,1:ncomp)
        
        !Numerical Derivative number i
        !f(x_i + delta_x_i)
        do xj = 1, ncomp
            do xi = 1, xj ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                !add delta_x to the xi
                molfractions(xj) = molfractions(xj) + delta
            
                ! calculate the function with the changed xj (second component which is derived to)
                x1 = 0.d0
                
                !-----
                ! calculate the needed derivs again with new molfracs
                call calculate_PCSAFT_functionparts_Trho(T,D,getprevious)
                call AHSDERIVS(getvector)
                ! X derivs
                call calculate_PCSAFT_functionparts_x1(T,D,getprevious)
                call AHSX1DERIVS(getvector)
                !-----
                
                call ADISPX1DERIVS(D,getvector)
                func_plus_delta(xi,xj) = adispx1_PCSAFT(i,xi) ! xi in x1 here as the first component xi derived to xj again
            
                molfractions = molfractions_original
            end do
        end do
        
        !f(x_i - delta_x_i)
        do xj = 1, ncomp
            do xi = 1, xj ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                !subtract delta_x from the xi
                molfractions(xj) = molfractions(xj) - delta
            
                ! calculate the function with the changed xi
                x1 = 0.d0
                
                !-----
                ! calculate the needed derivs again with new molfracs
                call calculate_PCSAFT_functionparts_Trho(T,D,getprevious)
                call AHSDERIVS(getvector)
                ! X derivs
                call calculate_PCSAFT_functionparts_x1(T,D,getprevious)
                call AHSX1DERIVS(getvector)
                !-----
                
                call ADISPX1DERIVS(D,getvector)
                func_minus_delta(xi,xj) = adispx1_PCSAFT(i,xi)! not sure if this is correct: xi in x1 here as the first component xi derived to xj again
            
                molfractions = molfractions_original
            end do
        end do
        
        !Numerical derivative = (func_plus_delta - func_minus_delta) / (2*delta)
        do xj = 1, ncomp
            do xi = 1, xj ! calculating the whole symmetrical matrix (double calculation of the same results not important here in the test routine)
                numerical_derivative(xi,xj) = (func_plus_delta(xi,xj) - func_minus_delta(xi,xj))/(2.d0*delta)
            end do
        end do
        
        !calculate the difference dx1 between numerical and analytical derivative
        ! assignment with dx2(i,:,:) = ... not possible! because the array numbering gets mixed up as dx2 has dimensions (i,30,30) and the others have (ncomp,ncomp)
        do xj = 1, ncomp
            do xi = 1, xj
                dx2(i,xi,xj) = analytical_derivative(xi,xj) - numerical_derivative(xi,xj)
            end do
        end do
        
    end do
    
    end subroutine test_adisp_x2
    
subroutine test_ar_x1(T,D,delta,dx1)
    ! The alternative routine doesn't change the last molfraction to 1-Sum(other molfractions)
    
    use module_nderivs
    use module_fluid_parameters
    
    implicit none
    
    double precision :: T, D
    double precision :: delta, sum_xi
    double precision, dimension(30) :: molfractions_original
    double precision, dimension(ncomp) :: numerical_derivative, analytical_derivative, func_plus_delta, func_minus_delta
    double precision, dimension(10,30) :: dx1
    double precision, dimension(nx1derivs,ncomp) :: x1
    double precision, dimension(nx2derivs,ncomp,ncomp) :: x2
    double precision, dimension(nx3derivs,ncomp,ncomp,ncomp) :: x3
    integer, dimension(nderivs) :: getvector
    double precision, dimension(nderivs) ::  resultvector
    integer :: i, xi, j
    
    !only works when ncomp > 1
    if (ncomp .EQ. 1) then
        RETURN
    end if
    
    !save the original molfractions
    molfractions_original = molfractions
    
    ! The derivatives to xi
    do i = 1, 10
        
        !re-initialize the getvector in every loop
        getvector = 0
        getvector(i) = 1
        
        !Analytical Derivative number i
        x1 = 0.d0
        x2 = 0.d0
        x3 = 0.d0
        call ARX1DERIVS(T, D, getvector, x1)
        analytical_derivative = x1(i,1:ncomp)
        
        !Numerical Derivative number i
        !f(x_i + delta_x_i)
        do xi = 1, ncomp
            !add delta_x to the xi
            molfractions(xi) = molfractions(xi) + delta
            
            ! calculate the function with the changed xi
            call ARDERIVS(T,D,getvector,resultvector,xi)
            func_plus_delta(xi) = resultvector(i)
            
            molfractions = molfractions_original
        end do
        
        !f(x_i - delta_x_i)
        do xi = 1, ncomp
            !subtract delta_x from the xi
            molfractions(xi) = molfractions(xi) - delta
            
            ! calculate the function with the changed xi
            call ARDERIVS(T,D,getvector,resultvector,xi)
            func_minus_delta(xi) = resultvector(i)
            
            molfractions = molfractions_original
        end do
        
        !Numerical derivative = (func_plus_delta - func_minus_delta) / (2*delta)
        do xi = 1, ncomp
            numerical_derivative(xi) = (func_plus_delta(xi) - func_minus_delta(xi))/(2*delta)
        end do
        
        !calculate the difference dx1 between numerical and analytical derivative
        dx1(i,1:ncomp) = analytical_derivative(1:ncomp) - numerical_derivative(1:ncomp)
        
    end do
    
    end subroutine test_ar_x1
    
subroutine test_ar(T_in, D_in, delta, difference_numerical_analytical)

! Henning Markgraf, April 2016

    ! Subroutine differentiates numerically with respect to T and rho and takes the difference to the analytical subroutine

    ! WARNING:
    ! July, 2016, Henning Markgraf
    ! Subroutine DOES NOT WORK WITH THE new T conversions in the program!
    ! To be able to use this subroutine, the last part in subroutine ARDERIVS
    ! needs to be changed like this:
    ! Add the following two lines.
    ! DERAR = 0.d0
    ! DERAR(1:nderivs) = ar_PCSAFT(1:nderivs)
    ! Make the following line a comment.
    ! call convert_T_derivs_Trho(get,ar_PCSAFT,DERAR)
    ! Didn't have time to make this routine differentiate numerically wrt tau :(

    use module_nderivs
    
    implicit none
    
    double precision :: T, D, T_in, D_in
    integer :: nrderiv, nrbasiceqn, nD, nT, i
    logical :: changeD
    double precision :: delta, deltaT, deltaD, divisor, factor
    double precision :: numerical_derivative, analytical_derivative, func_plus_delta, func_minus_delta
    integer, dimension(nderivs) :: getvector
    double precision, dimension(nderivs) ::  resultvector, difference_numerical_analytical

    difference_numerical_analytical = 0.d0
 
    do i = 2, nderivs
        
        deltaT = 0.d0
        deltaD = 0.d0
        T = T_in
        D = D_in
    
        call basic_data_of_deriv(i, nrbasiceqn, nD, nT, T, D, delta, deltaD, deltaT, changeD)
    
        !Analytical Derivative number i
        resultvector = 0.d0
        getvector = 0
        getvector(i) = 1
        call ARDERIVS(T, D, getvector, resultvector)
        analytical_derivative = resultvector(i)
    
        !numerical derivative of the Basic Equation
        !f(x+deltaX)
        D = D + deltaD
        T = T + deltaT
        getvector(nrbasiceqn) = 1
        getvector(i) = 0
        call ARDERIVS(T, D, getvector, resultvector)
        if (changeD) then
            divisor = D**(nD-1)
        else
            divisor = T**(nT-1)
        end if
        func_plus_delta = resultvector(nrbasiceqn)/divisor
        !f(x-deltaX)
        D = D - 2.d0*deltaD
        T = T - 2.d0*deltaT
        call ARDERIVS(T, D, getvector, resultvector)
        if (changeD) then
            divisor = D**(nD-1)
        else
            divisor = T**(nT-1)
        end if
        func_minus_delta = resultvector(nrbasiceqn)/divisor

        !((f(x+deltaX)-f(x-deltaX))/(2*deltaX))*X^nX
        T = T + deltaT
        D = D + deltaD
        if (changeD) then
            factor = D**nD
            numerical_derivative = factor*(func_plus_delta-func_minus_delta)/(2.d0*deltaD)
        else
            factor = T**nT
            numerical_derivative = factor*(func_plus_delta-func_minus_delta)/(2.d0*deltaT)
        end if
    
        !difference analytical and numerical (should be zero)
        difference_numerical_analytical(i) = analytical_derivative - numerical_derivative
    end do
    
    end subroutine test_ar
    
subroutine basic_data_of_deriv(nrderiv, nrbasiceqn, nD, nT, T, D, delta, deltaD, deltaT, changeD)

    implicit none
    
    double precision :: T, D, delta, deltaT, deltaD
    integer :: nrderiv, nrbasiceqn, nD, nT
    logical :: changeD

    !nrbasiceqn important!
    !Damit es klappt muss nrbasiceqn gleich eine Ableitung niedriger sein, und mit der gleichen Anzahl
    !Ableitungen zur anderen Variable, z.B. Ableitung nach T^2 und D^2 ist die basic equation dann
    !nach T und D^2 oder nach T^2 und D

    !get order of derivatives to T and D
    if (nrderiv .EQ. 1) then
        !error message
        !basic equation can't be tested
    else if (nrderiv .EQ. 2) then
        nrbasiceqn = 1
        nD = 1
        nT = 0
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 3) then
        nrbasiceqn = 2
        nD = 2
        nT = 0
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 4) then
        nrbasiceqn = 1
        nD = 0
        nT = 1
        deltaT = delta*T
        deltaD = 0.d0
        changeD = .FALSE.
    else if (nrderiv .EQ. 5) then
        nrbasiceqn = 4
        nD = 0
        nT = 2
        deltaT = delta*T
        deltaD = 0.d0
        changeD = .FALSE.
    else if (nrderiv .EQ. 6) then
        nrbasiceqn = 2 ! or 4 if changeD is true
        nD = 1
        nT = 1
        deltaT = delta*T
        deltaD = 0.d0
        changeD = .FALSE.
    else if (nrderiv .EQ. 7) then
        nrbasiceqn = 6
        nD = 1
        nT = 2
        deltaT = delta*T
        deltaD = 0.d0
        changeD = .FALSE.
    else if (nrderiv .EQ. 8) then
        nrbasiceqn = 3
        nD = 3
        nT = 0
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 9) then
        nrbasiceqn = 5
        nD = 0
        nT = 3
        deltaT = delta*T
        deltaD = 0.d0
        changeD = .FALSE.
    else if (nrderiv .EQ. 10) then
        nrbasiceqn = 6
        nD = 2
        nT = 1
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 11) then
        nrbasiceqn = 8
        nD = 4
        nT = 0
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 12) then
        nrbasiceqn = 8
        nD = 3
        nT = 1
        deltaD = 0.d0
        deltaT = delta*T
        changeD = .FALSE.
    else if (nrderiv .EQ. 13) then
        nrbasiceqn = 10
        nD = 2
        nT = 2
        deltaD = 0.d0
        deltaT = delta*T
        changeD = .FALSE.
    else if (nrderiv .EQ. 14) then
        nrbasiceqn = 9
        nD = 1
        nT = 3
        deltaD = delta*D
        deltaT = 0.d0
        changeD = .TRUE.
    else if (nrderiv .EQ. 15) then
        nrbasiceqn = 9
        nD = 0
        nT = 4
        deltaD = 0.d0
        deltaT = delta*T
        changeD = .FALSE.
    end if

end subroutine basic_data_of_deriv
