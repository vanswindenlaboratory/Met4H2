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

    ! module for file ptflash_2c_3p.f90
    module ptflash_2c_3p_module
    !global use inclusion
    use module_all_types
    use rhomix_pt_module
    use flash_module    
    use setup_module
    use vle_derivs_module
    use reduced_parameters_calc_module

    contains




!**************************************************************************
subroutine ptflash_NC_3P(gl,press, Temp, x_known, rho, x_vap, x_liq1, x_liq2, rhovap_est, &
                               & rholiq1_est, rholiq2_est, Phasefrac, iFlash, iter, errval)
!**************************************************************************
! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THREE-PHASE-EQUILIBRIA FOR MULTI COMPONENTS MIXTURES. FOLLOWING THE GIBBS' PHASE RULE 
! THE THREE PHASE REGIONS ARE LINES IN THE T-P-DIAGRAM FOR BINARY MIXTURES AND MAY BE AREAS FOR THREE AND MORE COMPONENT MIXTURES 
!
!--------------------------------------------------------------------------
! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
!   - THREE PHASE LINE ("Bubble" Point) : p, x Given beta_V = 0     --  iFlash = 1
!   - THREE PHASE LINE ("Dew" Point)    : p, x Given beta_L1 = 0    --  iFlash = 2
!   - THREE PHASE LINE                  : p, x Given beta_L2 = 0    --  iFlash = 3
!   - THREE PHASE LINE ("Bubble" Point) : T, x Given beta_V = 0     --  iFlash = 4
!   - THREE PHASE LINE ("Dew" Point)    : T, x Given beta_L1 = 0    --  iFlash = 5
!   - THREE PHASE LINE                  : T, x Given beta_L2 = 0    --  iFlash = 6
!   - THREE PHASE Flash: p, T and x VECTOR GIVEN                    --  iFlash = 7
!--------------------------------------------------------------------------
! Variables:
! INPUT:
!   press           - Pressure [MPa]
!   Temp            - Temperature [K]
!   x_vap           - Composition of the vapor phase
!   x_liq1          - Composition of the (first) fluid phase (lighter)
!   x_liq2          - Composition of the (second) fluid phase (heavier)
!   rhovap_est      - Estimated vapor phase density
!   rholiq1_est     - Estimated (first) liquid phase density 
!   rholiq2_est     - Estimated (second) liquid phase density
!   iFlash          - Flash type
!   Phasefrac       - Phasefractions according to:  Phasefrac(1)    :   nV / n      - Molar amount of vapor phase divided by overall molar amount of the mixture 
!                                                   Phasefrac(2)    :   nL1 / n     - Molar amount of (lighter) liquid phase 1 divided by overall molar amount of the mixture
!                                                   Phasefrac(3)    :   nL2 / n     - Molar amount of (heavier) liquid phase 2 divided by overall molar amount of the mixture


! OUTPUT:
!   errval          - Error value
!   iter            - number of iterations carried out
!   rho             - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, (solid, hydrate -- >  Seperate routine for solid substances) )
!   x_vap           - Composition of the vapor phase
!   x_liq1          - Composition of the (first) fluid phase (lighter)
!   x_liq2          - Composition of the (second) fluid phase (heavier)
!   Phasefrac       - Phasefractions


!------------------0--------------------------------------------------------
! A. Jäger August 2012 






implicit none

    type(type_gl) :: gl


double precision, dimension(30)::  x_known
double precision :: rhovap_est, rholiq1_est, rholiq2_est
integer :: iFlash

double precision :: press, temp
double precision, dimension(30) :: x_vap, x_liq1, x_liq2
double precision, dimension(3) :: Phasefrac

double precision, dimension(5) :: rho  
integer :: errval, iter

double precision, dimension(3):: Phasefrac_new
double precision, dimension(30):: x_vap_new, x_liq1_new, x_liq2_new
double precision, dimension(60, 60):: JacMatrix
double precision, dimension(60):: GibbsEQN, Delta_X, Var_X

double precision:: eps_Gibbs, eps_del, max_del, stepsize, eps_Gibbs_alt
double precision:: sum_vap, sum_liq1, sum_liq2, Temp_new, press_new
integer:: i, j, k, eqn

character(255):: herr

errval = 0
Delta_X = 1.D0
x_vap_new = 0.D0
x_liq1_new = 0.D0
x_liq2_new = 0.D0
phasefrac_new = 0.D0
Temp_new = 0.D0
press_new = 0.D0    
stepsize = 1.D0 

!Catch NaN
if ((press /= press) .or. (Temp /= Temp)) then
    errval = -4321
    return
endif

!If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
eps_Gibbs = 1.d-7
!If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
eps_del = 1.d-12
!If the criterion for the relative change of the unknowns is fulfilled, nevertheless the gibbs criterion has to be checked for convergence
!In such a case the gibbs criterion is lowered a little, to avoid slow convergence to a lower residuum when in principle converged
eps_Gibbs_alt = 1.d-5

!Iteration counter is set back to 1
iter = 1
!The GibbsEQN get a high starting value
GibbsEQN = 1.D10
!Starting values for max_del and Var_X
!max_del = maximum relative difference for the set of unknowns between two iterations 
max_del = 0.D0
!Current value for all variables
Var_X = 1.D0

!Check for starting values
Do i = 1, gl%ncomp
    if ((x_vap(i) <= 0.D0) .or. (x_liq1(i) <= 0.D0) .or. (x_liq2(i) <= 0.D0)) then
        errval = -1111   !Wrong (missing) inputs to a routine
        return
    End if
End Do

If ((Temp <= 0.D0) .or. (press <= 0.D0)) then
    errval = -1111   !Wrong (missing) inputs to a routine
    return
End if

if ((iFlash == 7) .and. (gl%ncomp < 3)) then
    errval = -9903   !Calculation not possible for binary mixtures (or pure substances)
    return    
end if

if ((iFlash == 1) .or. (iFlash == 4)) phasefrac(1) = 0   !Bubble point -- >  beta_V = 0
if ((iFlash == 2) .or. (iFlash == 5)) phasefrac(2) = 0   !Dew point -- >  beta_L1 = 0
if ((iFlash == 3) .or. (iFlash == 6)) phasefrac(3) = 0   !Dew point -- >  beta_L2 = 0

do i = 1, 60
                 
    call SysOfEqs_NC_3P(gl,press, Temp, x_known, x_vap, x_liq1, x_liq2, rhovap_est, rholiq1_est, rholiq2_est, &
                            & iFlash, Phasefrac, GibbsEQN, errval)
                           
    rho(1) = gl%rho_vap
    rho(2) = gl%rho_liq1
    rho(3) = gl%rho_liq2
    
    !Break Condition
    ! This is the first exit condition for the VLE iteration!
    if (maxval(dabs(GibbsEQN)) < eps_Gibbs) return
    
    if (errval /= 0) return
    call Jacobi_NC_3P(gl,press, temp, x_vap, x_liq1, x_liq2, iFlash,Phasefrac, JacMatrix, errval)
    Delta_X = - GibbsEQn
    
    eqn = 3 * gl%ncomp         
    
    call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)   
    
    ! Initialize
    sum_vap = 0.D0
    sum_liq1 = 0.D0
    sum_liq2 = 0.D0

    iter = i

    select case (iFlash)
! --------------------------------------------------------------------------------    
        case(1) !Bubble point (beta_V = 0)
            j = 1
            press_new = press
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    temp_new = temp + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(4) = temp_new
                    !If (ncomp > 2) then
                    phasefrac_new(2) = phasefrac(2) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(2)
                    phasefrac_new(3) = phasefrac(3) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(3)
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0) & 
                    !& .OR. (phasefrac_new(3) < 0.d0) .OR. (phasefrac_new(3) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
        case(2) !Dew point (beta_L1 = 0)
            j = 1
            press_new = press
            phasefrac_new(2) = phasefrac(2)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    temp_new = temp + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(4) = temp_new
                    !If (ncomp > 2) then
                    phasefrac_new(1) = phasefrac(1) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(1)
                    phasefrac_new(3) = phasefrac(3) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(3)
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) &
                    !& .OR. (phasefrac_new(3) < 0.d0) .OR. (phasefrac_new(3) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
        case(3) !beta_L2 = 0
            j = 1
            press_new = press
            phasefrac_new(3) = phasefrac(3)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    temp_new = temp + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(4) = temp_new
                    !If (ncomp > 2) then
                    phasefrac_new(1) = phasefrac(1) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(1)
                    phasefrac_new(2) = phasefrac(2) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(2)
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) &
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
         case(4) !Bubble point (beta_V = 0)
            j = 1
            Temp_new = Temp
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    press_new = press + Delta_X(j+2*gl%NCOMP-2) * 1.D-6
                    Var_X(4) = press_new
                    !If (ncomp > 2) then
                    phasefrac_new(2) = phasefrac(2) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(2)
                    phasefrac_new(3) = phasefrac(3) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(3)                    
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0) &
                    !& .OR. (phasefrac_new(3) < 0.d0) .OR. (phasefrac_new(3) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
        case(5) !Dew point (beta_L1 = 0)
            j = 1
            Temp_new = Temp
            phasefrac_new(2) = phasefrac(2)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    press_new = press + Delta_X(j+2*gl%NCOMP-2) * 1.D-6
                    Var_X(4) = press_new
                    !If (ncomp > 2) then
                    phasefrac_new(1) = phasefrac(1) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(1)
                    phasefrac_new(3) = phasefrac(3) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(3)
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) &
                    !& .OR. (phasefrac_new(3) < 0.d0) .OR. (phasefrac_new(3) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
        case(6) !beta_L2 = 0
            j = 1
            Temp_new = Temp
            phasefrac_new(3) = phasefrac(3)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    press_new = press + Delta_X(j+2*gl%NCOMP-2) * 1.D-6
                    Var_X(4) = press_new
                    !If (ncomp > 2) then
                    phasefrac_new(1) = phasefrac(1) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(1)
                    phasefrac_new(2) = phasefrac(2) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(2)
                    !end if
                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) ) then
                    !Commented, Andreas March 2016
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) &
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
        case(7) !T,p flash 
            j = 1
            Temp_new = Temp
            press_new = press
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    x_liq1_new(j) = x_liq1(j) + Delta_X(j+gl%NCOMP-1)
                    x_liq2_new(j) = x_liq2(j) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(1) = x_vap_new(j)
                    Var_X(2) = x_liq1_new(j)
                    Var_X(3) = x_liq2_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    x_liq1_new(gl%ncomp) = 1.D0 - sum_liq1
                    x_liq2_new(gl%ncomp) = 1.D0 - sum_liq2
                    phasefrac_new(1) = phasefrac(1) + Delta_X(j+2*gl%NCOMP-2)
                    Var_X(4) = phasefrac_new(1)
                    phasefrac_new(2) = phasefrac(2) + Delta_X(j+2*gl%NCOMP-1)
                    Var_X(5) = phasefrac_new(2)
                    phasefrac_new(3) = phasefrac(3) + Delta_X(j+2*gl%NCOMP)
                    Var_X(6) = phasefrac_new(3)

                End if
                sum_vap = sum_vap + x_vap_new(j)
                sum_liq1 = sum_liq1 + x_liq1_new(j)
                sum_liq2 = sum_liq2 + x_liq2_new(j)
                ! check if the stepsize is too large. 
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) < 0.D0) .OR. (x_liq1_new(j) < 0.D0) .OR. (x_liq2_new(j) < 0.D0) &
                    & .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq1_new(j) > 1.D0) .OR. (x_liq2_new(j) > 1.D0) )then
                    ! Andreas, March 2014
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) .OR. (phasefrac_new(2) < 0.d0) &
                    !& .OR. (phasefrac_new(2) > 1.d0) .OR. (phasefrac_new(3) > 1.d0) .OR. (phasefrac_new(3) < 0.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0 
                    sum_vap = 0.D0
                    sum_liq1 = 0.D0
                    sum_liq2 = 0.D0
                    phasefrac_new = phasefrac
                end if    
                j = j + 1            
            end do
! --------------------------------------------------------------------------------
    end select
    
    ! write the new values to the variables
    x_vap = x_vap_new
    x_liq1 = x_liq1_new
    x_liq2 = x_liq2_new
    Temp = Temp_new
    press = press_new  
    phasefrac = phasefrac_new
    !phasefrac(3) = 1.D0 - phasefrac(1) - phasefrac(2)
      
    !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
    max_del = 0.D0
    Do k = 1, j - 1     
        if(abs(delta_X(k) / Var_X(k)) > max_del) then
            max_del = abs(delta_X(k) / Var_X(k))
        end if
    end do
    
    if (((max_del) < eps_del)  .and. (maxval(dabs(GibbsEQN)) < eps_Gibbs)) return

    !Catch the case, that the phase compositions become the same
    !Andreas, October 2014: Added density check. Phase compositions can be quite similar, so a density check is required to be sure of the phases are the same and a trivial solution was found
    if (       ( (maxval(dabs(x_vap - x_liq1)) < 0.0000042d0) .AND. (dabs(gl%rho_vap - gl%rho_liq1) < 0.0000042d0) ) &
        & .OR. ( (maxval(dabs(x_vap - x_liq2)) < 0.0000042d0) .AND. (dabs(gl%rho_vap - gl%rho_liq2) < 0.0000042d0) ) &
        & .OR. ( (maxval(dabs(x_liq1 - x_liq2)) < 0.0000042d0) .AND. (dabs(gl%rho_liq1 - gl%rho_liq2) < 0.0000042d0) ) ) then
             errval = -4323
             return
    end if
 
      
end do

! Iteration failed!
if (i > 60) then
    errval = -2222
else
    rho(1) = gl%rho_vap
    rho(2) = gl%rho_liq1
    rho(3) = gl%rho_liq2
End if

end subroutine ptflash_NC_3P


!**************************************************************************
subroutine SysOfEqs_NC_3P(gl,P, T, x_known, x_vap, x_liq1, x_liq2, rhovap_est, rholiq1_est, rholiq2_est, iFlash, &
                        & Phasefrac, GibbsEQN, errval)
!**************************************************************************
! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
! EQUILIBRIUM CALCULATIONS.
! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
!--------------------------------------------------------------------------
!   - THREE PHASE LINE ("Bubble" Point) : p, x Given beta_V = 0     --  iFlash = 1
!   - THREE PHASE LINE ("Dew" Point)    : p, x Given beta_L1 = 0    --  iFlash = 2
!   - THREE PHASE LINE                  : p, x Given beta_L2 = 0    --  iFlash = 3
!   - THREE PHASE LINE ("Bubble" Point) : T, x Given beta_V = 0     --  iFlash = 4
!   - THREE PHASE LINE ("Dew" Point)    : T, x Given beta_L1 = 0    --  iFlash = 5
!   - THREE PHASE LINE                  : T, x Given beta_L2 = 0    --  iFlash = 6
!   - THREE PHASE Flash: p, T and x VECTOR GIVEN                    --  iFlash = 7
!--------------------------------------------------------------------------
! Variables:
! INPUT:
!   P               - Pressure
!   T               - Temperature
!   x_vap           - Composition of the vapor phase
!   x_liq1          - Composition of the (first) fluid phase
!   x_liq2          - Composition of the (second) fluid phase
!   rhovap_est      - Estimated vapor phase density
!   rholiq1_est     - Estimated (first) liquid phase density 
!   rholiq2_est     - Estimated (second) liquid phase density
!   iFlash          - Flash mode
!   Phasefrac       - Phasefractions according to:  Phasefrac(1)    :   nV / n      - Molar amount of vapor phase divided by overall molar amount of the mixture 
!                                                   Phasefrac(2)    :   nL1 / n     - Molar amount of (lighter) liquid phase 1 divided by overall molar amount of the mixture
!                                                   Phasefrac(3)    :   nL2 / n     - Molar amount of (heavier) liquid phase 2 divided by overall molar amount of the mixture

! OUTPUT:
!   errval          - Error value
!   GibbsEQN        - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
!--------------------------------------------------------------------------
! A. Jäger,  August 2012
 





implicit none

    type(type_gl) :: gl


double precision :: P, T
double precision, dimension(30):: x_known, x_vap, x_liq1, x_liq2
double precision:: rhovap_est, rholiq1_est, rholiq2_est
double precision, dimension(3):: Phasefrac
integer:: iFlash

double precision, dimension(60):: GibbsEQN
integer:: errval

double precision, dimension(30):: lnfi_vap, lnfi_liq1, lnfi_liq2      
double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq1, d_liq2
double precision:: Rmix
integer:: errorflag, Iphase, i

errval = 0
GibbsEQN = 0.D0

!Just in case, maybe not necessary
!------------------------------
gl%molfractions = x_known
call reduced_parameters_calc(gl,T)
!------------------------------ 
rhoredmix_orig = gl%rhoredmix
tredmix_orig = gl%tredmix

! write the vapor and liquid phase 1 density
! into the module variables for the new set of T, p, x_vap, and x_liq1:
Iphase = 5 !VL
call newvars (gl,P, T, x_vap, x_liq1, rhovap_est, rholiq1_est, Iphase, errorflag) 
if (errorflag /= 0) then
    errval = errorflag
    return
end if
gl%rho_liq1 = gl%rho_liq

! Get the density of liquid phase 2
gl%molfractions = x_liq2
call reduced_parameters_calc(gl,T)
iPhase = 1  !Liquid
! calculate the liquid phase 2 density
gl%rho_liq2 = rhomix_calc(gl,T, P, rholiq2_est, iPhase,0)

if (gl%rho_liq2 == 0.D0)then
    errval = -8888
    return
end if

!----------------------------------
! get the gas phase properties
gl%molfractions = x_vap
call reduced_parameters_calc(gl,T)
! get the gas phase density from the module
d_vap = gl%rho_vap
! Calculate the chemical potentials for the gas phase
call R_mix_calc(gl,Rmix)
!call dna_dni(T, d_vap, ChemPot_vap, 0)
!ChemPot_vap = ChemPot_vap !* Rmix * T    !Chemical Potentials of the vapor phase [J / mol]
call lnf_mix(gl,T, d_vap, p, lnfi_vap)

!----------------------------------
! get the liquid phase 1 properties
gl%molfractions = x_liq1
call reduced_parameters_calc(gl,T)
! get the gas phase density from the module
d_liq1 = gl%rho_liq1
! Calculate the chemical potentials for the liquid phase
call R_mix_calc(gl,Rmix)
!call dna_dni(T, d_liq1, ChemPot_liq1, 0)
!ChemPot_liq1 = ChemPot_liq1 !* Rmix * T    !Chemical Potentials of the liquid 1 phase [J / mol]
call lnf_mix(gl,T, d_liq1, p, lnfi_liq1)

!----------------------------------
! get the liquid phase 2 properties
gl%molfractions = x_liq2
call reduced_parameters_calc(gl,T)
! get the gas phase density from the module
d_liq2 = gl%rho_liq2
! Calculate the chemical potentials for the liquid phase
call R_mix_calc(gl,Rmix)
!call dna_dni(T, d_liq2, ChemPot_liq2, 0)
!ChemPot_liq2 = ChemPot_liq2 !* Rmix * T    !Chemical Potentials of the liquid 2 phase [J / mol]
call lnf_mix(gl,T, d_liq2, p, lnfi_liq2)

!----------------------------------
! Setting up the system of equations 
! for the minimization of the Gibbs free energy.

!Equality of fugacities. These equations are needed for each type of flash calculations and all numbers of components
Do i = 1, gl%ncomp
    GibbsEQN(i) = lnfi_vap(i) - lnfi_liq1(i)  
    GibbsEQN(gl%ncomp+i) = lnfi_vap(i) - lnfi_liq2(i)
End Do

!Additional mass balance equations. These equations are necessary in the following cases:
!Former version, changes made to enable setting beta_L1 to 0 and to calculate phasefractions in case of two components and phaseboundaries (IFlash 1 - 6)
!!1) 3 components and Flash type 1, 2, 3, or 4 ("Dew" or "Bubble" Line, Transformation of L1 to V or vice versa. L2 should stay unchanged)
!!2) more than 3 components in the mixture
!If ((ncomp > 3) .or. ((ncomp == 3) .and. (iFlash < 5))) then
!    Do i = 1, ncomp - 1
!       !Mass balance according to: beta_V * xi_V + beta_L1 * xi_L1 + beta_L2 * xi_L2 = zi  with beta_v + beta_L1 + beta_L2 = 1 
!       GibbsEQN(2*ncomp+i) = Phasefrac(1) * x_vap(i) + Phasefrac(2) * x_liq1(i) & 
!                           & + (1.D0 - Phasefrac(1) - Phasefrac(2)) * x_liq2(i) - x_known(i)
!    End do
!End if

!New System of equations, Andreas June 2013
!No limitations necessary, in every case the mass balance has to be solved
Do i = 1, gl%ncomp - 1
    !Mass balance according to: beta_V * xi_V + beta_L1 * xi_L1 + beta_L2 * xi_L2 = zi  with beta_v + beta_L1 + beta_L2 = 1 
    GibbsEQN(2*gl%ncomp+i) = Phasefrac(1) * x_vap(i) + Phasefrac(2) * x_liq1(i) + Phasefrac(3) * x_liq2(i) - x_known(i)
End do
!Additional equation, ensuring that beta_V + beta_L1 + beta_L2 = 1
GibbsEQN(3*gl%ncomp) = Phasefrac(1) + Phasefrac(2) + Phasefrac(3) - 1.D0 

! set the module variables back to original values
gl%rhoredmix = rhoredmix_orig
gl%tredmix = tredmix_orig
gl%molfractions = x_known

end subroutine SysOfEqs_NC_3P
!**************************************************************************

!**************************************************************************
subroutine Jacobi_NC_3P (gl,P, T, x_vap, x_liq1, x_liq2, iFlash, Phasefrac, JacMatrix, errval)
!**************************************************************************
! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM
! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
!--------------------------------------------------------------------------
!   - THREE PHASE LINE ("Bubble" Point) : p, x Given beta_V = 0     --  iFlash = 1
!   - THREE PHASE LINE ("Dew" Point)    : p, x Given beta_L1 = 0    --  iFlash = 2
!   - THREE PHASE LINE                  : p, x Given beta_L2 = 0    --  iFlash = 3
!   - THREE PHASE LINE ("Bubble" Point) : T, x Given beta_V = 0     --  iFlash = 4
!   - THREE PHASE LINE ("Dew" Point)    : T, x Given beta_L1 = 0    --  iFlash = 5
!   - THREE PHASE LINE                  : T, x Given beta_L2 = 0    --  iFlash = 6
!   - THREE PHASE Flash: p, T and x VECTOR GIVEN                    --  iFlash = 7
!--------------------------------------------------------------------------
! Variables:
! INPUT:
!   P               - Pressure
!   T               - Temperature
!   x_vap           - Composition of the vapor phase
!   x_liq1          - Composition of the (first) fluid phase
!   x_liq2          - Composition of the (second) fluid phase
!   iFlash          - Flash mode
!   Phasefrac       - Phasefractions according to:  Phasefrac(1)    :   nV / n      - Molar amount of vapor phase divided by overall molar amount of the mixture 
!                                                   Phasefrac(2)    :   nL1 / n     - Molar amount of (lighter) liquid phase 1 divided by overall molar amount of the mixture
!                                                   Phasefrac(3)    :   nL2 / n     - Molar amount of (heavier) liquid phase 2 divided by overall molar amount of the mixture
! OUTPUT:
!   errval          - Error value
!   JacMatrix       - 60 x 60 matrix containing the derivatives of all Gibbs-equations
!                     F_i with respect to all independent variables X_i
!--------------------------------------------------------------------------
! A. Jäger, August 2012 






implicit none

    type(type_gl) :: gl


double precision :: T, p
double precision, dimension(30) :: x_vap, x_liq1, x_liq2
double precision, dimension(3) :: Phasefrac
integer :: iFlash
double precision, dimension(60, 60) :: JacMatrix
integer :: errval

integer:: i, j
double precision:: d_vap, d_liq1, d_liq2
double precision, dimension(30)::  z

double precision:: rhoredmix_orig, tredmix_orig, Rmix
double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq1, dlnfidXj_liq2
double precision, dimension(30):: dlnphiidT_vap, dlnphiidP_vap, dlnphiidT_liq1,  dlnphiidP_liq1
double precision, dimension(30):: dlnphiidT_liq2, dlnphiidP_liq2


JacMatrix = 0.D0
z = 0.D0
errval = 0

rhoredmix_orig = gl%rhoredmix
tredmix_orig = gl%tredmix
z = gl%molfractions

! get the vapor phase properties
d_vap = gl%rho_vap
gl%molfractions = x_vap
call reduced_parameters_calc(gl,T)
call R_mix_calc(gl,Rmix) 
if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then  
    !Calculate the derivative of the fugacity with respect to T at constant
    !X and p for the vapor phase   
    call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
End if    
if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then   
    !Calculate the derivative of the fugacity with respect to p at constant
    !X and T for the vapor phase  
    call dlnphii_dP(gl,T, d_vap, dlnphiidP_vap)
End if       
!Calculate the derivative of the fugacity with respect to xj at constant
!p and T for the vapor phase
call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)

! get the liquid phase 1 properties
d_liq1 = gl%rho_liq1
gl%molfractions = x_liq1
call reduced_parameters_calc(gl,T)    
call R_mix_calc(gl,Rmix)
if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then  
    !Calculate the derivative of the fugacity with respect to T at constant
    !X and p for the vapor phase   
    call dlnphii_dT(gl,T, d_liq1, dlnphiidT_liq1)
End if    
if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then   
    !Calculate the derivative of the fugacity with respect to p at constant
    !X and T for the vapor phase  
    call dlnphii_dP(gl,T, d_liq1, dlnphiidP_liq1)
End if      
!Calculate the derivative of the fugacity with respect to xj at constant
!p and T for the (lighter) liquid phase
call dlnfi_dxj_TP (gl,T, d_liq1, dlnfidXj_liq1)      

! get the liquid phase 2 properties
d_liq2 = gl%rho_liq2
gl%molfractions = x_liq2
call reduced_parameters_calc(gl,T)    
call R_mix_calc(gl,Rmix)
!Calculate the derivative of the fugacity with respect to T at constant
!X and p for the (heavier) liquid phase 
if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then  
    !Calculate the derivative of the fugacity with respect to T at constant
    !X and p for the vapor phase   
    call dlnphii_dT(gl,T, d_liq2, dlnphiidT_liq2)
End if    
if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then   
    !Calculate the derivative of the fugacity with respect to p at constant
    !X and T for the vapor phase  
    call dlnphii_dP(gl,T, d_liq2, dlnphiidP_liq2)
End if            
!Calculate the derivative of the fugacity with respect to xj at constant
!p and T for the (heavier) liquid phase
call dlnfi_dxj_TP (gl,T, d_liq2, dlnfidXj_liq2)


!The Matrix is split into different section (compare to notes)
!A,B,C ..... M
!Furthermore the special case of bubble or dew points for three component mixtures must be considered


!All Sections that are present for ALL FLASH TYPES, independent of the number of components, i.e., the fugacities w.r.t. composition of all phases
Do j = 1, gl%ncomp -1
    Do i = 1, gl%ncomp
        !A
        JacMatrix(j,i) = dlnfidXj_vap(j, i)
        JacMatrix(j,gl%ncomp+i) = dlnfidXj_vap(j, i)
        !B
        JacMatrix(gl%ncomp-1+j, i) = -dlnfidXj_liq1(j, i)
        !C
        JacMatrix(2*gl%ncomp-2+j, i) = 0.D0
        !D
        JacMatrix(gl%ncomp-1+j, gl%ncomp+i) = 0.D0
        !E
        JacMatrix(2*gl%ncomp-2+j, gl%ncomp+i) = -dlnfidXj_liq2(j, i)
    End do
End do

!In case of "Dew" or "Bubble" point calculations the derivatives of the fugacities wrt T or p are needed
if(iFlash < 7) then 
    Do i = 1, gl%ncomp        
        if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then  
            !G
            JacMatrix(3*gl%ncomp-2,i) = dlnphiidT_vap(i) - dlnphiidT_liq1(i)
            !H
            JacMatrix(3*gl%ncomp-2,gl%ncomp+i) = dlnphiidT_vap(i) - dlnphiidT_liq2(i)
        end if
        if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
            !G
            JacMatrix(3*gl%ncomp-2,i) = dlnphiidp_vap(i) - dlnphiidp_liq1(i)
            !H
            JacMatrix(3*gl%ncomp-2,gl%ncomp+i) = dlnphiidp_vap(i) - dlnphiidp_liq2(i)
        end if       
    End do
End if

!If ncomp > 3 or ncomp == 3 and iflash < 5 the mass balance needs to be solved
!If ((ncomp > 3) .or. ((ncomp == 3) .and. (iFlash < 5))) then - not longer needed
Do j = 1, gl%ncomp - 1
    Do i = 1, gl%ncomp-1
        if (iFlash == 7) then
            !G  
            JacMatrix(3*gl%ncomp-2,i) = 0.D0
            !H
            JacMatrix(3*gl%ncomp-2,gl%ncomp+i) = 0.D0
        End if
        !F
            !UNSCHÖN!!!!
            JacMatrix(3*gl%ncomp-1,i) = 0.D0
            JacMatrix(3*gl%ncomp-1,gl%ncomp+i) = 0.D0   
        !I
            if (i == j) then 
                JacMatrix(j,2*gl%ncomp+i) = Phasefrac(1)   !In case of a bubble point, this is 0
            else
                JacMatrix(j,2*gl%ncomp+i) = 0.D0
            end if
        !J
            if (i == j) then
                JacMatrix(j+gl%ncomp-1,2*gl%ncomp+i) = Phasefrac(2)   !In case of a dew point, this is 0
            else
                JacMatrix(j+gl%ncomp-1,2*gl%ncomp+i) = 0.D0
            end if
        !K
            if (i == j) then
                JacMatrix(j+2*gl%ncomp-2,2*gl%ncomp+i) = Phasefrac(3) !1.D0 - Phasefrac(1) - Phasefrac(2)   
            else
                JacMatrix(j+2*gl%ncomp-2,2*gl%ncomp+i) = 0.D0
            end if
        !L
            if (iFlash < 7) then
                JacMatrix(3*gl%ncomp-2,2*gl%ncomp+i) = 0.D0    
            else
                JacMatrix(3*gl%ncomp-2,2*gl%ncomp+i) = x_vap(i)!x_vap(i) - x_liq2(i)
            end if
        !M
            if ((iFlash == 1) .or. (iFlash == 4).or. (iFlash == 7)) then  !"Bubble point" phasefrac(1) (beta_V) = 0 (Same equation needed for iFlash = 5)
                JacMatrix(3*gl%ncomp-1,2*gl%ncomp+i) = x_liq1(i)!x_liq1(i) - x_liq2(i)   
            end if
            if ((iFlash == 2) .or. (iFlash == 5)) then  !"Dew point" phasefrac(2) (beta_L1) = 0
                JacMatrix(3*gl%ncomp-1,2*gl%ncomp+i) = x_vap(i)!x_vap(i) - x_liq2(i)
            end if
            if ((iFlash == 3) .or. (iFlash == 6)) then  !phasefrac(3) (beta_L2) = 0
                JacMatrix(3*gl%ncomp-1,2*gl%ncomp+i) = x_liq1(i)!x_vap(i) - x_liq2(i)
            end if
        !N
            if ((iFlash == 1) .or. (iFlash == 4).or. (iFlash == 7)) then  !"Bubble point" phasefrac(1) (beta_V) = 0 (Same equation needed for iFlash = 5)
                JacMatrix(3*gl%ncomp,2*gl%ncomp+i) = x_liq2(i)  
            end if
            if ((iFlash == 2) .or. (iFlash == 5)) then  !"Dew point" phasefrac(2) (beta_L1) = 0
                JacMatrix(3*gl%ncomp,2*gl%ncomp+i) = x_liq2(i)
            end if
            if ((iFlash == 3) .or. (iFlash == 6)) then  !phasefrac(3) (beta_L2) = 0
                JacMatrix(3*gl%ncomp,2*gl%ncomp+i) = x_liq1(i)
            end if
        !O
            JacMatrix(3*gl%ncomp,i) = 0.D0
            JacMatrix(3*gl%ncomp,gl%ncomp+i) = 0.D0   
        !P  
            JacMatrix(j,3*gl%ncomp) = 0.D0
            JacMatrix(j+gl%ncomp-1,3*gl%ncomp) = 0.D0
            JacMatrix(j+2*gl%ncomp-2,3*gl%ncomp) = 0.D0
    End do
End do
!End if
if (iflash == 7) then
    JacMatrix(3*gl%ncomp-2,3*gl%ncomp) = 1.D0
else
    JacMatrix(3*gl%ncomp-2,3*gl%ncomp) = 0.D0
end if
JacMatrix(3*gl%ncomp-1,3*gl%ncomp) = 1.D0
JacMatrix(3*gl%ncomp,3*gl%ncomp) = 1.D0
      
! set the module variables back to original values
gl%rhoredmix = rhoredmix_orig
gl%tredmix = tredmix_orig
gl%molfractions = z


end subroutine Jacobi_NC_3P
!**************************************************************************



    end module ptflash_2c_3p_module
