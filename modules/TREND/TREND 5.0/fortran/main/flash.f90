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

    ! module for file flash.f90
    submodule (flash_module) impl
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use rhomix_pt_module
    use setup_module
    use phasenv_vbased_module
    use vle_derivs_module
    use reduced_parameters_calc_module
    use gibbsderivs_module

    contains






    !**************************************************************************
    module subroutine Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
        &  Nr_x_given, errval, iter)
    !DEC$ ATTRIBUTES DLLEXPORT :: Flash_PhaseBoundary
    !**************************************************************************
    !THIS SUBROUTINE CALLS THE FLASH CALCULATION ROUTINE WITH T,P,x INPUT AND DETERMINES WHICH TYPE OF PHASE EQUILIBRIUM IS PRESENT
    !A. Jäger August 2011, renamed and restructured November 2011
    ! Modified by J.Gernert, 01.2012
    ! Modified the call of the phasenv routine - J.Gernert, A.Jäger, 08.2012
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: errval, iter

    integer:: iFlash, Nr_x_given, IPhase_try
    double precision:: press_given, temp_given, p_spec, T_spec
    double precision:: x_vap_given(30), x_liq_given(30), x_vap_start(30), x_liq_start(30)
    double precision:: pT_calc(6), rho_calc(6,2)
    double precision:: x_calc(30,6), T_start, p_start
    integer:: points_found, n, ptsfnd, critpt, i
    logical:: converged

    !Save the given start values (if there are any)
    press_given = press
    temp_given = temp
    x_vap_given = x_vap
    x_liq_given = x_liq
    IPhase_try = 0
    !---------------------------------------------------------------------------
    !1) Assume, that 2 phases are present
    !---------------------------------------------------------------------------
    !Generate start values for the phase split
    call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
    if (errval == 0) then   !Andreas Januar 2014
        x_vap_start = x_vap
        x_liq_start = x_liq
        T_start = Temp
        p_start = press

        !!---------------------------------------------------------------------------
        !!2)  do n steps of successive substitution
        !!    and check whether the Gibbs-Energy of the split system is lower than the
        !!    Gibbs-Energy of the initial system
        !!---------------------------------------------------------------------------
        n = 5
        converged = .false.
        ! use the given (experimental) composition values instead of the calculated ones as start values
        if (x_vap_given(1) /= 0.d0) x_vap = x_vap_given
        if (x_liq_given(1) /= 0.d0) x_liq = x_liq_given
        call Succ_Sub(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, n, converged)
        !for pvap = 0.101325 nitrogen;water '0.9993;0.0007' the density solver fails, due to a too small estimation of the saturation temperature
        !try max 3 times to increase the estimated temperature
        if ((errval /= 0).and. ((iFlash == 1).or.(iFlash == 2))) then
            do i = 1,3
                Temp = Temp * 1.05d0
                call Succ_Sub(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, n, converged)
                if (errval == 0) exit
            enddo
        endif
        !check if the successive substitution already converged
        !This line was taken out, since for some mixture the succ_sub converged to two phases with the same composition giving a wrong point -- Andreas July 2012
        !if ((converged) .and. (abs(1.d0-rho_vap/rho_liq) > 0.05)) return

        !!---------------------------------------------------------------------------
        ! 3)  Try to calculate the phase equilibrium with Newton-Raphson,
        !     using the results of the successive substitution as start values.
        !!---------------------------------------------------------------------------
        !Andreas August 2012: Do only use the start values of Succ_Sub if the succ Sub did not converge to the trivial solution!
        if (maxval(dabs(x_vap - x_liq)) < 1.D-3) then
            press = p_start
            Temp = T_start
            x_vap = x_vap_start
            x_liq = x_liq_start
        End if

        ! first, assume a VLE situation for the density solver - added by J.G. 11.2012
        IPhase_try = 5
        call Flash_PhaseBoundary_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, &
            & vapfrac, iFlash, IPhase_try, Nr_x_given, errval, iter)

        ! if that fails, try to let the density solver find the correct density by means of the lower Gibbs energy
        if (errval /= 0) then
            press = p_start
            Temp = T_start
            x_vap = x_vap_start
            x_liq = x_liq_start
            IPhase_try = 0
            call Flash_PhaseBoundary_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, &
                & vapfrac, iFlash, IPhase_try, Nr_x_given, errval, iter)
        end if

        !!---------------------------------------------------------------------------
        ! 4)  If the phase dinsities merged, try again with the start values from
        !     the K-Value method (or given values)
        !!---------------------------------------------------------------------------
        if (errval == 0) then   ! catch the case that the iteration above was cancelled du to density calculation errors - J.G., 05. 2012
            if ((abs(1.d0-x_vap(1)/x_liq(1)) < 0.01) .AND. (abs(1.d0-gl%rho_vap/gl%rho_liq) < 0.01)) then
                if (x_vap_given(1) /= 0.d0) then
                    x_vap = x_vap_given
                else
                    x_vap = x_vap_start
                end if
                if (x_liq_given(1) /= 0.d0) then
                    x_liq = x_liq_given
                else
                    x_liq = x_liq_start
                end if
                if (temp_given /= 0.d0) then
                    temp = temp_given
                else
                    temp = t_start
                end if
                if (press_given /= 0.d0) then
                    press = press_given
                else
                    press = p_start
                end if
                call Flash_PhaseBoundary_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, &
                    & vapfrac, iFlash, IPhase_try, Nr_x_given, errval, iter)
            end if
        end if
    end if

    !!---------------------------------------------------------------------------
    ! 5)  If the flash calculation fails, try calculating the point from the
    !     phase envelope routine (very slow!)
    !!---------------------------------------------------------------------------
    !if ((errval /= 0) .OR. ((abs(1.d0-x_vap(1)/x_liq(1)) < 0.01) .AND. (abs(1.d0-rho_vap/rho_liq) < 0.01))) then
    !Changed that criterion, since for CO2 / O2 mixtures at 13 MPa a wrong solution was found -- Andreas September 2012
    ! Additionally changed, that in some cases the algorithm might find a bubble point, when a dew point is required and vice versa, in this case calculate the envelope
    if ((errval /= 0) .or. (((abs(1.d0-x_vap(1)/x_liq(1)) < 0.05d0) .AND. (abs(1.d0-gl%rho_vap/gl%rho_liq) < 0.05d0)))) then
        !Robin, Moni 08.03.19:
        ! For the input Tliq, it might be possible that gl%rho_vap > gl%rho_vap if an LLE occurs.
        ! Therefore, this test is deleted.
        !    if ((errval /= 0) .OR. ((abs(1.d0-x_vap(1)/x_liq(1)) < 0.05d0) .AND. (abs(1.d0-gl%rho_vap/gl%rho_liq) < 0.05d0)) .OR. &
        !& (gl%rho_vap > gl%rho_liq)) then
        if (iFlash <= 2) then   ! calculate phasenv up to specified pressure (for p,x'- or p,x"-flash)
            p_spec = press_given
            T_spec = 0.d0
        else if (iFlash <= 4) then   ! calculate phasenv up to specified temperature (for T,x'- or T,x"-flash)
            p_spec = 0.d0
            T_spec = temp_given
        end if
        call phasenv(gl,x_known, p_spec, T_spec, pT_calc, rho_calc, x_calc, points_found, iFlash, errval)
        ! if the bubble line calculation failed or the dew line proceeds to infitine pressures
        ! and the specified T/p is above the cricondenbar, try again with iFlash = 0 (dew line calc. with jump over crit point)
        if ((points_found == 0) .and. (((errval /= 0) .and. ((iFlash == 1) .OR. (iFlash == 3))) &
            & .OR. ((iFlash  == 2) .or. (iFlash == 4)))) then
            errval = 0
            call phasenv(gl,x_known, p_spec, T_spec, pT_calc, rho_calc, x_calc, points_found, 0, errval)
            errval = -1
            ! check if specified point might be on the line expected to be the bubble line (LLE?) BS 10/2020
            if((points_found == 0) .and. (((errval /= 0) .and. (iFlash == 4)))) then
                errval = 0
                call phasenv(gl,x_known, p_spec, T_spec, pT_calc, rho_calc, x_calc, points_found, -1, errval)
                errval = -2
            end if
        end if
        select case (iFlash)
        case (1)    ! p,x' flash
            if ((points_found > 0) .and. (errval /= 0)) then
                ! search in the phase boundary points for a specified pressure on the bubble line
                points_found = 0; critpt = 0; ptsfnd = 0
                do i = 1, gl%phasenv_pts   !search for the first point left from the critical point
                    if (gl%pointID(i) == 1) critpt = i
                    if (gl%pointID(i) == 5) ptsfnd = ptsfnd + 1
                    if ((gl%pointID(i) == 5) .and. (i > critpt) .and. (ptsfnd > 1)) then
                        points_found = ptsfnd
                        exit
                    end if
                end do
            end if
            if (points_found == 0) then
                errval = -4401  !error: no bubble point existing
                return
            end if
            temp = pT_calc(points_found)
            x_vap = x_calc(:,points_found)
            gl%rho_vap = rho_calc(points_found,1)
            gl%rho_liq = rho_calc(points_found,2)
            errval = 0
        case (2)    !! p,x" flash
            if (points_found == 0) then
                errval = -4402  !error: no dew point existing
                return
            end if
            ! this may lead to pressures on the bubble line in some cases (retrograde bubble line, p > pcrit)!
            temp = pT_calc(1)
            x_liq = x_calc(:,1)
            gl%rho_vap = rho_calc(1,2)
            gl%rho_liq = rho_calc(1,1)
            errval = 0  !Andreas April 2012 -- >  If one point was found the errval is set to 0 since it is assumed that a correct dew point was found and the routine failed later

        case (3) ! T,x' flash
            critpt = 1000; ptsfnd = 0 !SH 11/2016: changed init value for critpt from 0 to 1000
            if ((points_found > 0) .and. (errval /= 0)) then
                ! search in the phase boundary points for a specified temperature on the bubble line
                points_found = 0
                do i = 1, gl%phasenv_pts   !search for the first point left from the critical point
                    if (gl%pointID(i) == 1) critpt = i
                    if (gl%pointID(i) == 4) ptsfnd = ptsfnd + 1
                    if ((gl%pointID(i) == 4) .and. (i > critpt)) then
                        points_found = ptsfnd
                        exit
                    end if
                end do
            end if
            if (points_found == 0) then
                errval = -4403  !error: no bubble point existing
                return
            end if
            press = pT_calc(points_found)
            x_vap = x_calc(:,points_found)
            gl%rho_vap = rho_calc(points_found,1)
            gl%rho_liq = rho_calc(points_found,2)
            errval = 0

        case (4) ! T,x" flash
            ! this does not cover some cases with an open phase envelope and T < Tc, and LLE!
            if ((points_found > 0) .and. (errval /= 0)) then
                points_found = 0
                do i = 1, gl%phasenv_pts   !search for the first point right from the critical point
                    if (gl%pointID(i) == 1) exit
                    if (gl%pointID(i) == 4) then
                        points_found = 1
                        exit
                    end if
                end do
            end if

            if (points_found == 0) then
                errval = -4404  !error: no dew point existing
                return
            end if
            if(errval == -2) then!change liq <-> vap props here, dew poitn was found with "bubble" line
                press = pT_calc(1)
                x_liq = x_calc(:,1)
                gl%rho_vap = rho_calc(1,1)
                gl%rho_liq = rho_calc(1,2)
                errval = 0
            else
                press = pT_calc(1)
                x_liq = x_calc(:,1)
                gl%rho_vap = rho_calc(1,2)
                gl%rho_liq = rho_calc(1,1)
                errval = 0
            end if
        end select
    end if

    if ((errval == 0).and.(gl%savebounds_p .eqv. .true.)) then
        gl%pmin_old = press * gl%min_factor
        gl%pmax_old = press
    end if

    End subroutine Flash_PhaseBoundary
    !**************************************************************************


    !**************************************************************************
    module subroutine Flash_PhaseBoundary_calc(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
        & iPhase_try, Nr_x_given, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac     - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !---------------------------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines






    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(:, :), allocatable:: JacMatrix
    double precision, dimension(:), allocatable:: GibbsEQN, Delta_X, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn, max_iter
    character(255):: herr
    if (.not.allocated(JacMatrix)) allocate(JacMatrix(60,60))
    if (.not.allocated(GibbsEQN)) allocate(GibbsEQN(60))
    if (.not.allocated(Delta_X)) allocate(Delta_X(60))
    if (.not.allocated(Var_X)) allocate(Var_X(60))

    !z = x_known
    errval = 0
    Delta_X = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-12
    !maximum iterations
    max_iter = 100!50!25!30
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible -- >  exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------
    do i = 1, max_iter
        call SysOfEqs_PhaseBoundary(gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, iPhase_try, GibbsEQN, errval)

        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) return

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_PhaseBoundary(gl,press, temp, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JacMatrix, errval)
        Delta_X = - GibbsEQn

        !For all flash types in this routine the number of equations to be solved equals the number of components (equality of fugacities for each component)
        eqn = gl%ncomp

        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + Delta_X(j)
                    Var_X(j) = Temp_new
                end if
                sum_vap = sum_vap + x_vap_new(j)

                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(j)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + Delta_X(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + Delta_X(j)
                    Var_X(j) = Temp_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(j)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    press_new = press + Delta_X(j)*1.D-6
                    Var_X(j) = press_new*1.D6
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: press_new > 0
                ! the 4th condition: abs(Dp) < 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (press_new < 0.d0) .OR. (dabs(Delta_X(j)) > 1.d6)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + Delta_X(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    press_new = press + Delta_X(j)*1.D-6
                    Var_X(j) = press_new*1.D6
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: press_new < 0
                ! the 4th condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (press_new < 0.d0) .OR. (abs(Delta_X(j)) > 1.d6)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + Delta_X(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + Delta_X(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + Delta_X(gl%NCOMP - 1)
                        press_new = press + Delta_x(gl%NCOMP) *1.D-6
                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = press_new*1.D6
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: press_new < 0
                ! the 4th condition: Dp > 1 MPa
                ! The 5th condition: Temp_new < 0
                ! The 6th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (press_new < 0.d0) .OR. (abs(Delta_X(gl%NCOMP)) > 1.d6) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + Delta_X(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + Delta_X(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + Delta_X(gl%NCOMP - 1)
                        press_new = press + Delta_x(gl%NCOMP) *1.D-6
                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = press_new*1.D6
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: press_new < 0
                ! the 4th condition: Dp > 1 MPa
                ! The 5th condition: Temp_new < 0
                ! The 6th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (press_new < 0.d0) .OR. (Delta_X(gl%NCOMP) > 1.d6) &
                    & .OR. (Temp_new < 0.d0) .OR. (Delta_X(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        ! write the new values to the variables
        x_vap = x_vap_new!/sum_vap
        x_liq = x_liq_new!/sum_liq
        Temp = Temp_new
        press = press_new

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
        max_del = 0.D0
        Do k = 1, j - 1
            if(abs(delta_X(k) / Var_X(k)) > max_del) then
                max_del = abs(delta_X(k) / Var_X(k))
            end if
        end do

        if ((max_del) < eps_del) return

        !Catch unphysical temperatures and pressures
        if ((temp < 0.D0) .or. temp > 1000.D0) then
            errval = -4321
            exit
        end if
        if ((press < 0.D0) .or. press > 1000.D0) then
            errval = -4322
            exit
        end if
    end do

    ! If convergence was not reached after 30 iterations -- >  algorithm failed!
    if ((i > max_iter) .AND. (maxval(dabs(GibbsEQN)) > 1.D-6)) then
        errval = -2222
    End if

    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_calc
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_PhaseBoundary(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, iPhase_try, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS ON THE PHASE BOUNDARY.
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac     - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines







    implicit none

    type(type_gl) :: gl


    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag

    errval = 0
    GibbsEQN = 0.D0
    !If generation of new values fails, molfractions module variable will be set to zero -> no second try of calculation, TE and SH 02/17
    !z = 0.D0
    z = gl%molfractions
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:
    call newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    if (errorflag /= 0) then
        errval = errorflag
        ! set the module variables back to original values, Andreas Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if
    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = gl%rho_vap
    call lnf_mix(gl,T, d_vap, p, lnfi_vap)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
        ! set the module variables back to original values, Andreas Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = gl%rho_liq
    call lnf_mix(gl,T, d_liq, p, lnfi_liq)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        ! set the module variables back to original values, Andreas Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !!NEW VLE FOR DROPLETS
    !if (droplet) then
    !    !Pl > Pv, from module variable
    !    d_liq = rhomix_calc(gl,T, pl, 0.d0, 1, 1)
    !    call lnf_mix(T, d_liq, pl, lnfi_liq)
    !end if
    !!NEW VLE FOR DROPLETS

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp-1
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do
    GibbsEQN(gl%ncomp) = lnfi_vap(gl%ncomp) - lnfi_liq(gl%ncomp)

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_PhaseBoundary
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_PhaseBoundary (gl,P, T, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JacMatrix, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM. THIS ROUTINE GENERATES
    ! THE MATRIX NEEDED FOR PHASE BOUNDARY CALCULATIONS
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode (explanation, see routine "SysOfEqs")
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 6 or 7
    ! OUTPUT:
    !   errval  - Error value
    !   JacMatrix  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (iFlash 6 & 7)






    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: iFlash, Nr_x_given
    integer:: i, j
    double precision:: d_vap, d_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z, dlnphiidP_vap, dlnphiidP_liq
    double precision:: rhoredmix_orig, tredmix_orig

    JacMatrix = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. T and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        d_vap = gl%rho_vap
        z = gl%molfractions
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
        call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JacMatrix(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp, i) = dlnphiidT_vap(i) -  dlnphiidT_liq(i)
        end do

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. T and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
        call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JacMatrix(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp, i) = dlnphiidT_vap(i) - dlnphiidT_liq(i)
        end do
    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given. p and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        d_vap = gl%rho_vap
        z = gl%molfractions
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
        call dlnphii_dP(gl,T, d_vap, dlnphiidP_vap)
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnphii_dP(gl,T, d_liq, dlnphiidP_liq)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JacMatrix(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp, i) = dlnphiidP_vap(i) - dlnphiidP_liq(i)
        end do

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given. p and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
        call dlnphii_dP(gl,T, d_liq, dlnphiidP_liq)
        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnphii_dP(gl,T, d_vap, dlnphiidP_vap)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JacMatrix(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp, i) = dlnphiidP_vap(i) - dlnphiidP_liq(i)
        end do

    case (5)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, xi" and x' vector are given. T, p and
        ! xk" need to be calculated
        !--------------------------------------------------------------------------
        d_vap = gl%rho_vap
        z = gl%molfractions
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
        call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
        call dlnphii_dP(gl,T, d_vap, dlnphiidP_vap)
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
        call dlnphii_dP(gl,T, d_liq, dlnphiidP_liq)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JacMatrix(j, i) = dlnfidXj_vap(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JacMatrix(j, i) = dlnfidXj_vap(j+1, i)
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp-1, i) = dlnphiidT_vap(i) - dlnphiidT_liq(i)
            JacMatrix(gl%ncomp, i) = dlnphiidP_vap(i) - dlnphiidP_liq(i)
        end do
    case (6)
        !--------------------------------------------------------------------------
        ! Dew point calculation, xi' and x" vector are given. T, p and
        ! xk' need to be calculated
        !--------------------------------------------------------------------------
        d_vap = gl%rho_vap
        z = gl%molfractions
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
        call dlnphii_dP(gl,T, d_vap, dlnphiidP_vap)
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
        call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
        call dlnphii_dP(gl,T, d_liq, dlnphiidP_liq)
        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JacMatrix(j, i) = -dlnfidXj_liq(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JacMatrix(j, i) = -dlnfidXj_liq(j+1, i)
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JacMatrix(gl%ncomp-1, i) = dlnphiidT_vap(i) - dlnphiidT_liq(i)
            JacMatrix(gl%ncomp, i) = dlnphiidP_vap(i) - dlnphiidP_liq(i)
        end do
        !--------------------------------------------------------------------------
        case default
        errval = -1111
        !--------------------------------------------------------------------------
    end select

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_PhaseBoundary
    !**************************************************************************

    !**************************************************************************
    module subroutine PTX_startvals_PhaseBoundary (gl,P, T, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF ESTIMATED START VALUES FOR TEMPERATURE,
    ! PRESSURE AND COMPOSITION FOR A PHASE BOUNDARY CALCULATION.
    ! START VALUES CAN E CALCULATED FOR ONE OF THE FOLLOWING FLASH CALCULATIONS:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !-------------------------------------------------------------------------
    ! Variables
    !
    !   P       - Pressure - can be INPUT or OUTPUT, depending on iFlash
    !   T       - Temperature - can be INPUT or OUTPUT, depending on iFlash
    !   x_known - Composition vector, INPUT - can be overall composition x, x' or x"
    !   x_vap   - Vapor phase composition, OUTPUT
    !   x_liq   - Liquid phase composition, OUTPUT
    !   vapfrac    - Molar vapor fraction, OUTPUT; the given vapfrac works as an initial estimate (INPUT) for the rachford rice solver
    !   iFlash  - Flash mode, INPUT
    !   errval  - Error value, OUTPUT
    !--------------------------------------------------------------------------
    ! NOTE:
    ! All input or input/output variables are mandatory, depending on the flash
    ! situation. E.g. for iFlash = 1 the temperature is unknown and cannot be given.
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       ALL INPUT OR INPUT/OUTPUT VARIABLES THAT ARE UNKNOWN HAVE TO BE
    !       SET TO ZERO WHEN CALLING PTX_startvals
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !--------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: P, T, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: iFlash
    integer:: errval

    double precision, dimension(30):: psat
    double precision:: sum_liq, sum_vap, help
    integer:: i

    x_vap = 0.D0
    x_liq = 0.D0
    errval = 0
    psat = 0.D0
    sum_liq = 1.D0
    sum_vap = 1.D0
    help = 0.d0

    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. Start values for T and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        x_liq = x_known
        !Calculate the Temperature:
        if (T == 0.d0) then
            T = Tsat_iter (gl,p, x_known, 1)
        end if
        if (T <= 0.d0) then
            errval = -4321
            return
        end if
        ! vapor fraction
        vapfrac = 0.D0
        ! calculate the vapor phase composition
        sum_vap = 0.D0
        do i = 1, gl%ncomp
            !check if the component is supercritical at the given temperature
            if (gl%Tc(i) > T) then
                ! calculate the saturation pressure for each pure component from a generalized vapor pressure equation
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/T))

            else if (gl%Tc(i) < 50.d0) then
                !exception for hydrogen and helium, extrapolaption of the generalized vapor pressure equation, Beckmüller
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/(0.8d0*gl%Tc(i))))
                psat(i) = gl%pc(i)+((gl%pc(i)-psat(i))/(gl%Tc(i)-0.8d0*gl%Tc(i)))*(T-gl%Tc(i))
                if ((((p/gl%pc(i))>8.d0).and.(t/gl%Tc(i))>2.d0).or.(p/gl%pc(i))>20.d0) then !"highly scientific" factor based on decades of experience
                    psat(i) = psat(i)/(p/gl%pc(i))
                end if

            else
                ! here a reasonable estimate for a virtual vapor pressure must be set
                psat(i) = gl%pc(i) ! does that work ??? not for hydrogen mixtures
            end if
            !Andreas January 2014
            !x_liq(i) = x_known(i)
            x_vap(i) = x_liq(i)*psat(i)/p
            if (x_vap(i) < 1.D-14) x_vap(i) = 1.D-14
            sum_vap = sum_vap + x_vap(i)
        end do

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. Start values for T and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        !Calculate the Temperature if no start value is given:
        x_vap = x_known
        if (T == 0.d0) then
            T = Tsat_iter (gl,p, x_known, 2)
        end if
        if (T <= 0.d0) then
            errval = -4321
            return
        end if
        ! vapor fraction
        vapfrac = 1.D0
        ! calculate the liquid phase composition
        sum_liq = 0.D0
        do i = 1, gl%ncomp
            !check if the component is supercritical at the given temperature
            if (gl%Tc(i) > T) then
                ! calculate the saturation pressure for each pure component from a generalized vapor pressure equation
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/T))

            else if (gl%Tc(i) < 50.d0) then
                !exception for hydrogen and helium, extrapolaption of the generalized vapor pressure equation, Beckmüller
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/(0.8d0*gl%Tc(i))))
                psat(i) = gl%pc(i)+((gl%pc(i)-psat(i))/(gl%Tc(i)-0.8d0*gl%Tc(i)))*(T-gl%Tc(i))
                if ((((p/gl%pc(i))>8.d0).and.(t/gl%Tc(i))>2.d0).or.(p/gl%pc(i))>20.d0) then !"highly scientific" factor based on decades of experience
                    psat(i) = psat(i)/(p/gl%pc(i))
                end if

            else
                ! here a reasonable estimate for a virtual vapor pressure must be set
                psat(i) = gl%pc(i) ! does that work ???
            end if
            !Andreas January 2014
            !x_vap(i) = x_known(i)
            x_liq(i) = x_vap(i)*p/psat(i)




            if (x_liq(i) < 1.D-14) x_liq(i) = 1.D-14
            sum_liq = sum_liq + x_liq(i)
        end do

    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given. Start values for p and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        x_liq = x_known
        do i = 1, gl%ncomp
            !check if the component is supercritical at the given temperature
            if (gl%Tc(i) > T) then
                ! calculate the saturation pressure for each pure component from a generalized vapor pressure equation
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/T))

            else if (gl%Tc(i) < 50.d0) then
                !exception for hydrogen and helium, extrapolaption of the generalized vapor pressure equation, Beckmüller
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/(0.8d0*gl%Tc(i))))
                psat(i) = gl%pc(i)+((gl%pc(i)-psat(i))/(gl%Tc(i)-0.8d0*gl%Tc(i)))*(T-gl%Tc(i))
                if ((((p/gl%pc(i))>8.d0).and.(t/gl%Tc(i))>2.d0).or.(p/gl%pc(i))>20.d0) then !"highly scientific" factor based on decades of experience
                    psat(i) = psat(i)/(p/gl%pc(i))
                end if

            else
                ! here a reasonable estimate for a virtual vapor pressure must be set
                psat(i) = gl%pc(i) ! does that work ???
            end if
            ! calculate the extimated bubble point pressure from Raoult's law and Dalton's law
            ! This equation uses the following relations:
            !       p*x(i)" = psat(i)*x(i)'     (1) (Dalton'law + Raoult's law)
            !       Sum(x(i)") = 1              (2)
            ! Eqn. (1) is put into Eqn (2), then solved for p
            help = help + x_known(i)*psat(i)
        end do
        ! Calculate the pressure if no start values are given:
        if (p == 0.d0) p = help
        ! vapor fraction
        vapfrac = 0.D0
        ! Calculate the vapor phase composition from Dalton's law
        sum_vap = 0.D0
        do i = 1, gl%ncomp
            !Andreas January 2014
            !x_liq(i) = x_known(i)
            x_vap(i) = x_liq(i)*psat(i)/p
            if (x_vap(i) < 1.D-14) x_vap(i) = 1.D-14
            sum_vap = sum_vap + x_vap(i)
        end do

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given. Start values for p and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        x_vap = x_known
        do i = 1, gl%ncomp
            !check if the component is supercritical at the given temperature
            if (gl%Tc(i) > T) then
                ! calculate the saturation pressure for each pure component from a generalized vapor pressure equation
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/T))

            else if (gl%Tc(i) < 50.d0) then
                !exception for hydrogen and helium, extrapolaption of the generalized vapor pressure equation, Beckmüller
                psat(i) = gl%pc(i)*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/(0.8d0*gl%Tc(i))))
                psat(i) = gl%pc(i)+((gl%pc(i)-psat(i))/(gl%Tc(i)-0.8d0*gl%Tc(i)))*(T-gl%Tc(i))
                if ((((p/gl%pc(i))>8.d0).and.(t/gl%Tc(i))>2.d0).or.(p/gl%pc(i))>25.d0) then !"highly scientific" factor based on decades of experience
                    psat(i) = psat(i)/(p/gl%pc(i))
                end if

            else
                ! here a reasonable estimate for a virtual vapor pressure must be set
                psat(i) = gl%pc(i) ! does that work ???
            end if
            ! calculate the extimated bubble point pressure from Raoult's law and Dalton's law
            ! This equation uses the following relations:
            !       p*x(i)" = psat(i)*x(i)'     (1) (Dalton's law + Raoult's law)
            !       Sum(x(i)') = 1              (2)
            ! Eqn. (1) is put into Eqn (2), then solved for p
            help = help + x_known(i)/psat(i)
        end do
        ! Calculate the pressure if no start values are given:
        if (p == 0.d0) p = 1.d0/help
        ! vapor fraction
        vapfrac = 1.D0
        ! Calculate the liquid phase composition from Dalton's law
        sum_liq = 0.D0
        do i = 1, gl%ncomp
            !Andreas January 2014
            !x_vap(i) = x_known(i)
            x_liq(i) = x_vap(i)*p/psat(i)
            if (x_liq(i) < 1.D-14) x_liq(i) = 1.D-14
            sum_liq = sum_liq + x_liq(i)
        end do
        !--------------------------------------------------------------------------
        case default
        errval = -1111
    end select


    ! Check for Sum(x_i) = 1
    do i = 1, gl%ncomp
        x_vap(i) = x_vap(i)/sum_vap
        x_liq(i) = x_liq(i)/sum_liq
    end do

    end subroutine PTX_startvals_PhaseBoundary
    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_pT(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, errval, iter)
    !**************************************************************************
    !THIS SUBROUTINE CALLS THE FLASH CALCULATION ROUTINE WITH T,P,x INPUT AND DETERMINES WHICH TYPE OF PHASE EQUILIBRIUM IS PRESENT
    !A. Jäger August 2011, renamed and restructured November 2011








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: errval, iter
    double precision:: GIBBS_VapLiq, GIBBS_LiqLiq, G_vap, G_liq !, G_CALC

    double precision, dimension(30):: x_vap_vapliq, x_liq_vapliq, x_vap_old, x_liq_old
    double precision::beta_vapliq, rho_vap_vapliq, rho_liq_vapliq, rho_liq_liqliq, rho_vap_liqliq
    integer:: iPhase_try

    Gibbs_Vapliq = 1.D10
    Gibbs_LiqLiq = 1.D10
    x_vap_old = x_vap
    x_liq_old = x_liq
    ! Try to calculate the phase equilibrium with iphase 0 first. This means, that the density solver chooses the density with the lower Gibbs-Energy for each phase
    ! This procedure might fail during iteration such that the density switches between two values. In that case the iteration will fail.
    iPhase_try = 5
    call Flash_pT_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
    if (errval /= 0) then
        x_vap = x_vap_old
        x_liq = x_liq_old
        iPhase_try = 0
        call Flash_pT_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
    end if
    ! try again with improved start values from a few steps of successive substitution
    !    if (errval /= 0) then
    !        x_vap = x_vap_old
    !        x_liq = x_liq_old
    !        call Succ_Sub(press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, 5, errval, 5, converged)
    !        !check if the successive substitution lead to identical phase compositions
    !        do k = 1, ncomp
    !            if (abs(x_vap(k) - x_liq(k))/x_vap(k) < 1.d-3) then
    !                x_vap = 0.d0    ! delete the wrong phase compositions
    !                x_liq = 0.d0
    !                exit
    !            end if
    !        end do
    !        call Flash_pT_calc(press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
    !    end if

    !If the iteration fails, first an vapor liquid equilibrium is assumed and calculated and then an liquid liquid equilibrium. For both calculations the overall Gibbs-Energy is
    !calculated. The solution with the lower overall Gibbs-Energy is the correct one
    if (errval /= 0) then
        x_vap = x_vap_old
        x_liq = x_liq_old
        iPhase_try = 2  !try calculation by fixing the vapor density first
        call Flash_pT_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
        if (errval /= 0) then
            x_vap = x_vap_old
            x_liq = x_liq_old
            iPhase_try = 3 !try calculation by fixing the liquid density next
            call Flash_pT_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
        end if
        if (errval == 0) then
            !Save the solution of this calculation
            x_vap_vapliq = x_vap
            x_liq_vapliq = x_liq
            beta_vapliq = vapfrac
            rho_vap_vapliq = gl%rho_vap
            rho_liq_vapliq = gl%rho_liq
            !CALCULATE THE OVERALL GIBBS-ENERGY OF THE VAPOR LIQUID EQUILIBRIUM SOLUTION:
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)
            G_vap = G_CALC(gl,Temp, rho_vap_vapliq, 0)
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            G_liq = G_CALC(gl,Temp, rho_liq_vapliq, 0)
            !This is the overall Gibbs-Energy of the vapor liquid equlibrium
            GIBBS_VapLiq = vapfrac * G_vap + (1.D0 - vapfrac) * G_liq
        End if

        iPhase_try = 1
        x_vap = x_vap_old
        x_liq = x_liq_old
        call Flash_pT_calc(gl,press, temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)

        if (errval == 0) then
            !CALCULATE THE OVERALL GIBBS-ENERGY OF THE VAPOR LIQUID EQUILIBRIUM SOLUTION:
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)
            rho_vap_liqliq = gl%rho_vap
            rho_liq_liqliq = gl%rho_liq
            G_vap = G_CALC(gl,Temp, rho_vap_liqliq, 0)
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            G_liq = G_CALC(gl,Temp, rho_liq_liqliq, 0)
            !This is the overall Gibbs-Energy of the liquid liquid equlibrium
            GIBBS_LiqLiq = vapfrac * G_vap + (1.D0 - vapfrac) * G_liq
        End if

        !Error catching, if both tries failed the density iteration failed at one point:
        If ((Gibbs_vapliq == 1.D10) .And. (Gibbs_liqliq == 1.D10)) then
            return!errval = -8888
        End if

        !Choose the solution with the lower overall Gibbs-Energy:
        If (Gibbs_vapliq < Gibbs_liqliq) then
            x_vap = x_vap_vapliq
            x_liq = x_liq_vapliq
            vapfrac = beta_vapliq
            gl%rho_vap = rho_vap_vapliq
            gl%rho_liq = rho_liq_vapliq
        End if

    End if

    End subroutine Flash_pT
    !**************************************************************************


    !**************************************************************************
    module subroutine Flash_pT_calc(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF TWO PHASE EQUILIBRIA
    ! WITH THE INPUT VARIABLES T,P AND OVERALL COMPOSITION x
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !-------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATION CAN BE PERFORMED:
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 6 or 7
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger March 2011 (flash 6 & 7)
    ! A. Jäger Oct 2011 (flash 8)
    ! A. Jäger November 2011 Splittet in separate routine






    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: z, x_vap_new, x_liq_new, x_save
    double precision:: rholiq_est, rhovap_est, Delta_X_min
    integer:: iPhase_try
    integer:: errval, iter

    double precision, dimension(:, :), allocatable:: JacMatrix
    double precision, dimension(:), allocatable:: GibbsEQN, Delta_X, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new, T, p
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn, x_given, itermax
    character(255):: herr
    if(.not.(allocated(JacMatrix)))allocate(JacMatrix(60,60))
    if(.not.(allocated(GibbsEQN)))allocate(GibbsEQN(60))
    if(.not.(allocated(Delta_X)))allocate(Delta_X(60))
    if(.not.(allocated(Var_X)))allocate(Var_X(60))


    !double precision, dimension(50,4) :: x_compstore

    z = x_known
    errval = 0
    Delta_X = 1.D0
    x_save = 0.D0
    x_given = 0


    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    if (gl%seawater .or. gl%el_present) then
        eps_Gibbs = 1.d-6
        itermax = 200
        Delta_X_min = 1.d-9
    else
        eps_Gibbs = 1.d-8
        itermax = 30
        Delta_X_min = 1.d-10
    endif
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-12

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if only one phase composition is given, save it!
    if ((x_vap(1) == 0.D0) .AND. (x_liq(1) > 0.D0)) then
        x_save = x_liq
        x_given = 1
    else if ((x_vap(1) > 0.D0) .AND. (x_liq(1) == 0.D0)) then
        x_save = x_vap
        x_given = 2
    end if

    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        p = press
        T = Temp
        call PTX_startvals_pT (gl,p, T, x_known, x_vap, x_liq, vapfrac, errval)
    end if

    if (errval /= 0) then
        return
    end if

    ! use the given values for x_vap or x_liq (if there are any!) as start values instead of the calculated values from K-Factors
    if (x_given == 1) x_liq = x_save
    if (x_given == 2) x_vap = x_save
    if (Temp == 0.D0) Temp = T
    if (press == 0.D0) press = p

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !-----------------------------------------------------------

    do i = 1, itermax !30

        call SysOfEqs_pT(gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, GibbsEQN, errval)
        ! This is the exit condition for the VLE iteration!
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            ! calculate vapfrac as mean value from all components
            vapfrac = 0.d0
            do j = 1, gl%ncomp
                vapfrac = vapfrac + (z(j) - x_liq(j))/(x_vap(j) - x_liq(j))
            end do
            vapfrac = vapfrac / gl%ncomp
            if ((vapfrac .lt. 0.d0) .or. (vapfrac .gt. 1.d0)) then
                errval = -5566
            end if
            return
        end if


        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_pT(gl,press, temp, x_vap, x_liq, vapfrac, JacMatrix, errval)
        Delta_X = - GibbsEQn


        !When T and p are given, the phase compositions are unknown. So 2 * ncomp - 2 equations have to be solved
        eqn = 2*(gl%ncomp - 1)

        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
        if (errval == 1) then
            errval = -4444
            !write(*,*) 'Error in LUdecomp: ', herr
            return
        end if

        ! Additional break critertion in case GibbsEQN gets stuck at about 1.d-8   --- J.G. 09.2012
        if (maxval(dabs(Delta_X)) < Delta_X_min) then
            ! calculate vapfrac as mean value from all components
            vapfrac = 0.d0
            do j = 1, gl%ncomp
                vapfrac = vapfrac + (z(j) - x_liq(j))/(x_vap(j) - x_liq(j))
            end do
            vapfrac = vapfrac / gl%ncomp
            if ((vapfrac .lt. 0.d0) .or. (vapfrac .gt. 1.d0)) then
                errval = -5566
            end if
            return
        end if

        ! Reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        !Update the variables
        j = 1
        Temp_new = Temp
        press_new = press
        do while (j < gl%ncomp + 1)
            if (j < gl%ncomp) then
                x_vap_new(j) = x_vap(j) + Delta_X(j)
                x_liq_new(j) = x_liq(j) + Delta_X(gl%ncomp -1 + j)!*1.05d0
                Var_X(j) = x_vap_new(j)
                Var_X(gl%ncomp-1+j) = x_liq_new(j)
            else
                x_liq_new(gl%ncomp) = 1.D0 - sum_liq !- x_liq(1) * 0.02d0
                x_vap_new(gl%ncomp) = 1.D0 - sum_vap
            End if
            sum_liq = sum_liq + x_liq_new(j)
            sum_vap = sum_vap + x_vap_new(j)
            ! check if the stepsize is too large.
            ! The 1st condition: (x(j)_vap - z(j))*(x(j)_liq - z(j)) < 0  ( z is inbetween xliq and xvap)
            ! The 2nd condition: x(j)_vap AND x(j)_liq > 0
            ! have to be fulfilled for all components
            ! Otherwise the stepsize is reduced
            !if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0)) then

            !------------------------------------------------------------------------
            !UNCOMMENTED FOR THE CALCULATION OF VLc EQUILIBRIA!!! Andreas July 2011
            !                if (abs(x_vap_new(j) - x_liq_new(j)) < 1.D-4) then
            !                    x_vap_new(j) = 1 - x_liq_new(j)
            !                end if
            !                if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0) &
            !                    &   .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq_new(j) > 1.D0)) then
            if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0) &
                & .OR. (x_vap_new(j) >= 1.D0) .OR. (x_liq_new(j) >= 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                & .OR. (x_liq_new(j) <= 0.D0)) then          ! addtional check included for x <= 0  J.G., 09.2012
                !------------------------------------------------------------------------
                stepsize = 2.D0
                Delta_X = Delta_X/stepsize
                if (maxval(dabs(Delta_X)) < 1.D-10) then
                    errval = -3333
                    return
                end if
                stepsize = 1.D0
                j = 0
                sum_liq = 0.D0
                sum_vap = 0.D0
            end if


            j = j + 1
        end do

        ! write the new values to the variables
        x_vap = x_vap_new!/sum_vap
        x_liq = x_liq_new!/sum_liq
        Temp = Temp_new
        press = press_new

        ! calculate vapfrac as mean value from all components
        vapfrac = 0.d0
        do j = 1, gl%ncomp
            vapfrac = vapfrac + (z(j) - x_liq(j))/(x_vap(j) - x_liq(j))
        end do
        vapfrac = vapfrac / gl%ncomp

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
        !    max_del = 0.D0
        !    Do k = 1, j - 1
        !        if(abs(delta_X(k) / Var_X(k)) > max_del) then
        !            max_del = abs(delta_X(k) / Var_X(k))
        !        end if
        !    end do
        !
        !    if ((max_del) < eps_del) return

        !Catch the case, that the phase compositions and densities become the same (trivial solution or close to the crit. point) - J.G. 11.2012
        do k = 1, gl%ncomp
            if ((abs(x_vap(k) - x_liq(k)) < 0.0000042d0) .AND. (abs(gl%rho_vap - gl%rho_liq) < 1.d-3)) then
                errval = -4323
                return
            end if
        end do
    end do

    ! Iteration failed!
    if((i > itermax)) then
        errval = -2222      !30
        ! else if (i > 30) then
        ! errval = -2222
        !else
        !errval = -2222
    end if

    end subroutine Flash_pT_calc
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_pT(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS.
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !-------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATION CAN BE PERFORMED:
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (Modified)
    ! A. Jäger November 2011 Splittet in separate routine
    ! B. Semrau January 2019 Seawater integration







    implicit none

    type(type_gl) :: gl


    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap, chempot_liq, chempot_vap, chempot_orig, diff, chempot
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq, dxk_dxj_vap, dxk_dxj_liq, rmix, g, p_vap, p_liq
    integer:: i, j, errorflag

    !gl%seawater = .true.

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    lnfi_vap = 0.d0
    lnfi_liq = 0.d0
    chempot_vap = 0.d0
    chempot_liq = 0.d0
    errorflag = 0

    if(gl%seawater) then
        gl%seawaterflash = .true.
    else
        gl%seawaterflash = .false.
    end if

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:

    call newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    !


    if (errorflag /= 0) then
        errval = errorflag
        ! set the module variables back to original values  !Andreas, Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if
    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = gl%rho_vap
    call R_mix_calc(gl, Rmix)

    if (gl%seawater .or. gl%el_present) then
        if(gl%seawater) then
            p = gl%sea%seap
        else
            p = gl%gepress
        end if

        call dna_dni(gl,T, d_vap, lnfi_vap, 0)          !chempot / RT written on lnfi_vap for more simple handling
        lnfi_vap = lnfi_vap * Rmix * t

    else
        call lnf_mix(gl,T, d_vap, p, lnfi_vap)
    end if

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if ((lnfi_vap(1) == 0.D0) .and. (chempot_vap(1) .eq. 0.d0)) then
        ! set the module variables back to original values  !Andreas, Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        errval = -7777
        return
    end if

    !----------------------------------
    ! get the liquid phase properties
    !g =g_Sea_calc(gl, 280.d0, 0.1d0, 55500d0)
    x_liq(1) = x_liq(1)! - 0.02d0 * x_liq(1)  !new try
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = gl%rho_liq
    call R_mix_calc(gl, Rmix)

    if(gl%seawater .or. gl%el_present) then
        if(gl%seawater) then
            p = gl%sea%seap
        else
            p = gl%gepress
        end if
        if(gl%seawater) then
            call chempot_num_reac(gl, t, d_liq, lnfi_liq,1)
        elseif(gl%el_present) then
            call chempot_brine(gl, t, d_liq, lnfi_liq)
        end if
        p_liq = p_calc(gl, t, d_liq, 0)
        gl%seacalc = .false.
    else
        call lnf_mix(gl,T, d_liq, p, lnfi_liq)

    end if

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if ((lnfi_liq(1) == 0.D0) .and. (chempot_liq(1) .eq. 0.d0))then
        errval = -7778
        ! set the module variables back to original values  !Andreas, Jan 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    !if(gl%seawater) then
    !    lnfi_vap = chempot_vap
    !    lnfi_liq = chempot_liq
    !end if

    do i = 1, gl%ncomp-1
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        !if((gl%seawater)) then
        !    GibbsEQN(i) = GibbsEQN(i) / lnfi_vap(i) * 100.d0
        !end if
        ! the nth equation
    end do
    GibbsEQN(gl%ncomp) = lnfi_vap(gl%ncomp) - lnfi_liq(gl%ncomp)



    ! for the case of a p,T-flash another n-2 set of equations
    ! is needed (comes from the mass balance):
    ! Error handling required: if one of the comositions gets too close to the overall composition the next step would be a division through 0. Andreas August 2012
    ! Increased to 1.D-10, since this causes unnecessary errors for e.g. VLc equilibria, where V and Lc are very rich on CO2, Andreas March 2014 (Originally 1.D-5)
    !Commented Andy/Moni 2018/11/16
    !if ((maxval(abs(x_vap - z)) < 1.D-10) .or. (maxval(abs(x_liq - z)) < 1.D-10)) then
    !errval = -5566
    !! set the module variables back to original values  !Andreas, Jan 2013
    !gl%rhoredmix = rhoredmix_orig
    !gl%tredmix = tredmix_orig
    !gl%molfractions = z
    !return
    !end if
    do j = 1, gl%ncomp-2
        !calculate the derivative d(xk)/(d(xj) for vapor and liquid phase, k = ncomp-1
        dxk_dxj_vap = (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/(x_vap(j) - z(j))
        dxk_dxj_liq = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/(x_liq(j) - z(j))
        GibbsEQN(gl%ncomp + j) = dxk_dxj_vap - dxk_dxj_liq
    end do

    if(gl%seawater .or. gl%el_present) then
        GibbsEQN(gl%ncomp + gl%ncomp -1) = gl%gepress - p_calc(gl, t, d_liq,0)
    end if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_pT
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_pT (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac        - Molar vapor fraction
    ! OUTPUT:
    !   errval  - Error value
    !   JacMatrix  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (iFlash 6 & 7)






    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: i, j
    double precision:: d_vap, d_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq, dchempotidxj_liq, dchempotidxj_vap, dchempotdx
    double precision, dimension(30):: Z
    double precision:: rhoredmix_orig, tredmix_orig, Rmix

    JacMatrix = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    !--------------------------------------------------------------------------
    ! p-T flash calculation, p and T and the overall x vector are given.
    ! x' and x" need to be calculated
    !--------------------------------------------------------------------------
    ! get the liquid phase properties


    d_liq = gl%rho_liq
    z = gl%molfractions
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    call R_mix_calc(gl, Rmix)

    if (gl%seawater .or. gl%el_present) then
        !gl%seacalc = .true.
        call d2na_dnidxj_PT(gl, T, d_liq,dlnfidXj_liq, 0)
        dlnfidXj_liq = dlnfidXj_liq * Rmix * T
        !gl%seacalc = .false.
    else
        call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
        do i =1,gl%ncomp
            dlnfidXj_liq(:,i) = dlnfidXj_liq(:,i) !* Rmix * T! * gl%wm(i)
        end do
    end if

    ! get the vapor phase properties
    d_vap = gl%rho_vap
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    call R_mix_calc(gl, Rmix)

    if (gl%seawater .or. gl%el_present) then
        call d2na_dnidxj_PT(gl, T, d_vap, dchempotdx, 0)
        do i =1,gl%ncomp
            dlnfidXj_vap(:,i) = dchempotdx(:,i) * Rmix * T !* gl%wm(i)
        end do
    else
        call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
    end if



    do j = 1, (gl%ncomp-1)   ! goes through the 2*(n-1) columns, which are the derivativs to x_j
        do i = 1, gl%ncomp   ! goes through the first n-1 rows, which are the functions F_i
            JacMatrix(j, i) = dlnfidXj_vap(j, i)
            ! goes through the columns j = n to 2*(n-1) and the rows i = 1 to n-1, which are the derivs d(Fi)/d(x_liq(j)
            JacMatrix(gl%ncomp-1 + j, i) = - dlnfidXj_liq(j, i)
        end do
    end do
    do j = 1, gl%ncomp-2
        JacMatrix(j, gl%ncomp+j) = - (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/((x_vap(j) - z(j))**2)   !Ableitung von F_n+j nach x_j_vap
        JacMatrix(gl%ncomp-1, gl%ncomp+j) = 1.D0/(x_vap(j) - z(j))    !Ableitung von f_n+j nach x_n-1_vap
        JacMatrix(gl%ncomp-1+j, gl%ncomp+j) = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/((x_liq(j) - z(j))**2)    !Ableitung von F_n+j nach x_j_liq
        JacMatrix(2*(gl%ncomp-1), gl%ncomp+j) = - 1.D0/(x_liq(j) - z(j))  !Ableitung von f_n+j nach x_n-1_liq
    end do
    !--------------------------------------------------------------------------

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_pT
    !**************************************************************************

    !**************************************************************************
    module subroutine PTX_startvals_pT (gl,P, T, x_known, x_vap, x_liq, vapfrac, errval)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF ESTIMATED START VALUES FOR TEMPERATURE,
    ! PRESSURE AND COMPOSITION FOR A P-T FLASH CALCULATION.
    ! STARTVALUES FOR THE FOLLOWING CASE WILL BE GENERATED:
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! THE PROCEDURE IS TAKEN FROM
    !Michels & Mollerup 2004
    !Thermodynamic Models: Fundamentals & Computational Aspects
    !pp.220 - 222
    !--------------------------------------------------------------------------
    ! Variables
    !
    !   P       - Pressure - can be INPUT or OUTPUT, dependingon iFlash
    !   T       - Temperature - can be INPUT or OUTPUT, dependingon iFlash
    !   x_known - Composition vector, INPUT - can be overall composition x, x' or x"
    !   x_vap   - Vapor phase composition, OUTPUT
    !   x_liq   - Liquid phase composition, OUTPUT
    !   vapfrac    - Molar vapor fraction, OUTPUT; the given vapfrac works as an initial estimate (INPUT) for the rachford rice solver
    !   errval  - Error value, OUTPUT
    !--------------------------------------------------------------------------
    ! NOTE:
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       ALL INPUT OR INPUT/OUTPUT VARIABLES THAT ARE UNKNOWN HAVE TO BE
    !       SET TO ZERO WHEN CALLING PTX_startvals
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !--------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: P, T, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: errval


    double precision, dimension(30):: K_Val

    integer:: i

    x_vap = 0.D0
    x_liq = 0.D0
    errval = 0


    !--------------------------------------------------------------------------
    ! p-T flash, p, T and overall composition vector x are given.
    ! Start values for x" and x' need to be calculated
    !--------------------------------------------------------------------------
    !Calculate the K-Values from a generalized Wilson-Equation (Gerg 2004, page 62, Eq. (5.61))
    !and then calculate the molefractions for the vapor and liquid phase with an initial estimate
    !for the vaporfraction ß

    !*********************************************************************
    !Johannes & Andreas Aug 2011: Changed the starting values according to
    !Michels & Mollerup 2004
    !Thermodynamic Models: Fundamentals & Computational Aspects
    !pp.220 - 222
    !*********************************************************************

    K_val = 0.D0

    !Gernert et al. (2014), Eq. (10)
    Do i = 1, gl%ncomp
        !Wilson-Equation- Calculate the K-values
        !K_Val(i) = dexp(dlog(gl%pc(i)/p) + 5.373D0*(1.D0+gl%accen(i))*(1.D0-gl%Tc(i) /T))
        !Changed the mathematical formtulation of Wilsons equation for the calculation of negative pressures (sometimes necessary for critical points)
        !Andreas Jäger, August 2016.
        K_Val(i) = gl%pc(i)/p * dexp(5.373D0*(1.D0+gl%accen(i))*(1.D0-gl%Tc(i) /T))
    End do

    call RachRice (gl,K_Val, x_known, x_vap, x_liq, vapfrac, errval)
    !    Do i = 1, ncomp
    !        !Vapor fraction
    !        x_vap(i) = K_Val(i) * x_known(i) / (1 - vapfrac + vapfrac * K_Val(i))
    !        sum_vap = sum_vap + x_vap(i)
    !        !Liquid fraction
    !        x_liq(i) = x_known(i) / (1 - vapfrac + vapfrac * K_Val(i))
    !        sum_liq = sum_liq + x_liq(i)
    !    End Do

    ! Check for Sum(x_i) = 1
    !do i = 1, ncomp
    !    x_vap(i) = x_vap(i)/sum_vap
    !    x_liq(i) = x_liq(i)/sum_liq
    !end do

    end subroutine PTX_startvals_pT
    !**************************************************************************


    !**************************************************************************
    module subroutine Flash_ph(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, h_spec, iPhase_try, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF A TWO PHASE EQUILIBRIUM
    ! FOR THE INPUT VARIABLES: p, h and x
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PH-FLASH:     p, h AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   h_spec      - Specified (given) enthalpy value
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger Oct 2011 (flash 8)





    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est, h_spec
    integer:: iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(60):: GibbsEQN, Delta_X, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    character(255):: herr

    errval = 0
    Delta_X = 1.D0

    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-12
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-8
    !If the difference between the calculated and specified enthalpy is smaller than eps_h convergence is reached
    !eps_h = 1.d-6

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    ! The startvalue for T must be given from outside elsewise -- >  error
    ! So the same routine for startvalues as in the pt-flash is used (p and x given, estimated T)
    if (Temp == 0.D0) then
        errval = -1111
        return
    end if
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        call PTX_startvals_pT (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, errval)
    end if

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------
    do i = 1, 100
        call SysOfEqs_ph(gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, h_spec, iPhase_try, GibbsEQN, errval)

        ! This is the exit condition for the VLE iteration!
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) return
        ! usually the equation for the enthalpy (last equation) has the biggest residuum. So a separate convergence criterion for the enthalpy is needed (iteration to enthalpy differences of 10^-8 are senseless)
        !if (dabs(GibbsEQN(2*NCOMP - 1)) < eps_h) return

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_ph(gl,press, temp, x_vap, x_liq, vapfrac, JacMatrix, errval)
        Delta_X = - GibbsEQn

        !In case of a ph flash 2 * ncomp - 1 equations have to be solved. the unknowns are the compositions of both phases and the temperature
        eqn = 2*gl%ncomp - 1

        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
        if (errval == 1) then
            errval = -4444
            !write(*,*) 'Error in LUdecomp: ', herr
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        !Update the variables
        j = 1
        press_new = press
        do while (j < gl%ncomp + 1)
            if (j < gl%ncomp) then
                x_vap_new(j) = x_vap(j) + Delta_X(j)
                x_liq_new(j) = x_liq(j) + Delta_X(gl%ncomp -1 + j)
                Var_X(j) = x_vap_new(j)
                Var_X(gl%ncomp-1+j) = x_liq_new(j)
            else
                x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                Temp_new = Temp + Delta_X(2* gl%ncomp - 1)
                Var_X(2* gl%ncomp - 1) = Temp_new
            End if
            sum_liq = sum_liq + x_liq_new(j)
            sum_vap = sum_vap + x_vap_new(j)
            ! check if the stepsize is too large.
            ! The 1st condition: (x(j)_vap - z(j))*(x(j)_liq - z(j)) < 0  ( z is inbetween xliq and xvap)
            ! The 2nd condition: x(j)_vap AND x(j)_liq > 0
            ! have to be fulfilled for all components
            ! Otherwise the stepsize is reduced
            !if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0)) then

            !------------------------------------------------------------------------
            !UNCOMMENTED FOR THE CALCULATION OF VLc EQUILIBRIA!!! Andreas July 2011
            !                if (abs(x_vap_new(j) - x_liq_new(j)) < 1.D-4) then
            !                    x_vap_new(j) = 1 - x_liq_new(j)
            !                end if
            !                if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0) &
            !                    &   .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq_new(j) > 1.D0)) then
            if ((x_vap_new(j) < 0.D0) .OR. (x_liq_new(j) < 0.D0).OR. (x_vap_new(j) > 1.D0) .OR. (x_liq_new(j) > 1.D0) &
                & .Or. ((x_vap_new(j) - x_known(j)) * (x_liq_new(j) - x_known(j)) > 0.D0)) then
                !------------------------------------------------------------------------
                stepsize = 2.D0
                Delta_X = Delta_X/stepsize
                if (maxval(dabs(Delta_X)) < 1.D-10) then
                    errval = -3333
                    return
                end if
                stepsize = 1.D0
                j = 0
                sum_liq = 0.D0
                sum_vap = 0.D0
            end if


            j = j + 1
        end do
        ! --------------------------------------------------------------------------------


        ! write the new values to the variables
        x_vap = x_vap_new!/sum_vap
        x_liq = x_liq_new!/sum_liq
        Temp = Temp_new
        press = press_new

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del
        max_del = 0.D0
        Do k = 1, j - 1
            if(abs(delta_X(k) / Var_X(k)) > max_del) then
                max_del = abs(delta_X(k) / Var_X(k))
            end if
        end do

        if ((max_del) < eps_del) return

        !Catch unphysical temperatures
        if ((temp < 0.D0) .or. temp > 1000.D0) then
            errval = -4321
            exit
        end if
    end do

    ! Iteration failed!
    if (i > 100) errval = -2222

    end subroutine Flash_ph
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_ph(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, h_spec, iPhase_try, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS.
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PH-FLASH:     p, h AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   h_spec      - Specified (given) enthalpy value
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger October 2011








    implicit none

    type(type_gl) :: gl


    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est, h_spec
    integer::  iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq, dxk_dxj_vap, dxk_dxj_liq,  h_liq, h_vap, beta_h !H_CALC
    integer:: i, j, errorflag, reasonable_beta

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:
    call newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    if (errorflag /= 0) then
        errval = errorflag
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if
    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = gl%rho_vap
    call lnf_mix(gl,T, d_vap, p, lnfi_vap)
    h_vap = H_CALC(gl,T,d_vap, 0)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = gl%rho_liq
    call lnf_mix(gl,T, d_liq, p, lnfi_liq)
    h_liq = H_CALC(gl,T,d_liq, 0)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp-1
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do
    GibbsEQN(gl%ncomp) = lnfi_vap(gl%ncomp) - lnfi_liq(gl%ncomp)

    ! for the case of a p,h-flash another n-2 set of equations
    ! is needed (comes from the mass balance):
    do j = 1, gl%ncomp-2
        !calculate the derivative d(xk)/(d(xj) for vapor and liquid phase, k = ncomp-1
        dxk_dxj_vap = (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/(x_vap(j) - z(j))
        dxk_dxj_liq = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/(x_liq(j) - z(j))
        GibbsEQN(gl%ncomp + j) = dxk_dxj_vap - dxk_dxj_liq
    end do

    ! In case of a ph-flash, an additional equation is needed: h(T,p,x) - h_specified = 0
    beta_h = 0.D0
    vapfrac = 0.D0
    reasonable_beta = 0
    do j = 1, gl%ncomp-1
        beta_h = (z(i) - x_liq(i))/(x_vap(i) - x_liq(i))
        !if ((beta_h >= 0.D0) .AND. (beta_h <= 1.D0)) then
        vapfrac = vapfrac + (z(i) - x_liq(i))/(x_vap(i) - x_liq(i))
        reasonable_beta = reasonable_beta + 1
        !end if
    end do
    if (reasonable_beta == 0) then
        vapfrac = 0.65D0
    else
        vapfrac = vapfrac / reasonable_beta
    End if

    GibbsEQN(2*gl%ncomp - 1) = (h_vap * vapfrac + (1.D0 - vapfrac) * h_liq) - h_spec

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_ph
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_ph (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PH-FLASH:     p, h AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac        - Molar vapor fraction
    ! OUTPUT:
    !   errval  - Error value
    !   JacMatrix  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! A. Jäger October 2011








    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: i, j
    double precision:: d_vap, d_liq, dh_dT_vap, dh_dT_liq,dbeta_dxj_vap, dbeta_dxj_liq, h_vap, h_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z
    double precision, dimension(30):: dh_dxj_vap, dh_dxj_liq
    double precision:: rhoredmix_orig, tredmix_orig

    JacMatrix = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    h_vap = 0.D0
    h_liq = 0.D0
    dh_dT_liq = 0.D0
    dh_dT_vap = 0.D0
    dbeta_dxj_vap = 0.D0
    dbeta_dxj_liq = 0.D0
    dh_dxj_vap = 0.D0
    dh_dxj_liq = 0.D0

    !--------------------------------------------------------------------------
    ! p-h flash calculation, p and h and the overall x vector are given.
    ! x',x" and T need to be calculated
    !--------------------------------------------------------------------------
    ! get the liquid phase properties
    d_liq = gl%rho_liq
    z = gl%molfractions
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
    call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
    call dh_dx_TP (gl,T, d_liq, dh_dxj_liq)
    dh_dT_liq = dh_dT_px (gl,T, d_liq)
    h_liq = H_CALC(gl,T, d_liq, 0)
    ! get the vapor phase properties
    d_vap = gl%rho_vap
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
    call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
    call dh_dx_TP (gl,T, d_vap, dh_dxj_vap)
    dh_dT_vap = dh_dT_px (gl,T, d_vap)
    h_vap = H_CALC(gl,T, d_vap, 0)
    do j = 1, (gl%ncomp-1)   ! goes through the 2*(n-1) columns, which are the derivativs to x_j
        do i = 1, gl%ncomp   ! goes through the first n-1 rows, which are the functions F_i
            JacMatrix(j, i) = dlnfidXj_vap(j, i)
            ! goes through the columns j = n to 2*(n-1) and the rows i = 1 to n-1, which are the derivs d(Fi)/d(x_liq(j)
            JacMatrix(gl%ncomp-1 + j, i) = - dlnfidXj_liq(j, i)
        end do
        dbeta_dxj_liq = (z(j) - x_vap(j)) / (x_vap(j) - x_liq(j))**2 / (gl%ncomp-1) !Added the /ncomp Andreas Jun 2013
        dbeta_dxj_vap = (x_liq(j) - z(j)) / (x_vap(j) - x_liq(j))**2 / (gl%ncomp-1) !Added the /ncomp Andreas Jun 2013
        JacMatrix(j,2*gl%ncomp - 1) = dbeta_dxj_vap * (h_vap - h_liq) + vapfrac * dh_dxj_vap(j)             !Goes through the columns 1 to n-1 in the 2n-1th row, which are the derivatives of h w.r.t. x in the vapor phase
        JacMatrix(gl%ncomp - 1 +j,2*gl%ncomp - 1) = dbeta_dxj_liq * (h_vap - h_liq) + (1.D0 - vapfrac) * dh_dxj_liq(j)  !Goes through the columns n to n-2 in the 2n-1th row, which are the derivatives of h w.r.t. x in the liquid phase
    end do
    do j = 1, gl%ncomp-2
        JacMatrix(j, gl%ncomp+j) = - (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/((x_vap(j) - z(j))**2)   !Ableitung von F_n+j nach x_j_vap
        JacMatrix(gl%ncomp-1, gl%ncomp+j) = 1.D0/(x_vap(j) - z(j))    !Ableitung von f_n+j nach x_n-1_vap
        JacMatrix(gl%ncomp-1+j, gl%ncomp+j) = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/((x_liq(j) - z(j))**2)    !Ableitung von F_n+j nach x_j_liq
        JacMatrix(2*(gl%ncomp-1), gl%ncomp+j) = - 1.D0/(x_liq(j) - z(j))  !Ableitung von f_n+j nach x_n-1_liq
    end do
    do i = 1, gl%ncomp
        JacMatrix(2*gl%ncomp - 1, i) = dlnphiidT_vap(i) - dlnphiidT_liq(i) !Derivatives of ln(fi_vap) - ln(fi_vap) with respect to T at const. p and x
    end do
    if (gl%ncomp > 2) then
        do i = gl%ncomp + 1, 2*gl%ncomp - 2
            JacMatrix(2*gl%ncomp - 1, i) = 0.D0 !Derivatives of Fn+j w.r.t. T
        end do
    end if
    JacMatrix(2*gl%ncomp -1, 2*gl%ncomp -1) = vapfrac * dh_dT_vap + (1.D0-vapfrac) * dh_dT_liq
    !--------------------------------------------------------------------------

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_ph
    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_ps(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, s_spec, iPhase_try, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF A TWO PHASE EQUILIBRIUM
    ! FOR THE INPUT VARIABLES: p, s and x
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PS-FLASH:     p, s AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   s_spec      - Specified (given) entropy value
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger Oct 2011 (flash 8)





    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est, s_spec
    integer:: iPhase_try
    integer:: errval, iter


    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(60):: GibbsEQN, Delta_X, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    character(255):: herr

    errval = 0
    Delta_X = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-8
    !If the difference between the calculated and specified entropy is smaller than eps_s convergence is reached
    !eps_s = 1.d-6

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    ! The startvalue for T must be given from outside elsewise -- >  error
    ! So the same routine for startvalues as in the pt-flash is used (p and x given, estimated T)
    if (Temp == 0.D0) then
        errval = -1111
        return
    end if
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        call PTX_startvals_pT (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, errval)
    end if

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------
    do i = 1, 100
        call SysOfEqs_ps(gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, s_spec, iPhase_try, GibbsEQN, errval)

        ! This is the exit condition for the VLE iteration!
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) return
        ! usually the equation for the enthalpy (last equation) has the biggest residuum. So a separate convergence criterion for the enthalpy is needed (iteration to enthalpy differences of 10^-8 are senseless)
        !if (dabs(GibbsEQN(2*NCOMP - 1)) < eps_s) return

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_ps(gl,press, temp, x_vap, x_liq, vapfrac, JacMatrix, errval)
        Delta_X = - GibbsEQn

        !In case of a ps flash 2 * ncomp - 1 equations have to be solved. the unknowns are the compositions of both phases and the temperature
        eqn = 2*gl%ncomp - 1

        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
        if (errval == 1) then
            errval = -4444
            !write(*,*) 'Error in LUdecomp: ', herr
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        !Update the variables
        j = 1
        press_new = press
        do while (j < gl%ncomp + 1)
            if (j < gl%ncomp) then
                x_vap_new(j) = x_vap(j) + Delta_X(j)
                x_liq_new(j) = x_liq(j) + Delta_X(gl%ncomp -1 + j)
                Var_X(j) = x_vap_new(j)
                Var_X(gl%ncomp-1+j) = x_liq_new(j)
            else
                x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                Temp_new = Temp + Delta_X(2* gl%ncomp - 1)
                Var_X(2* gl%ncomp - 1) = Temp_new
            End if
            sum_liq = sum_liq + x_liq_new(j)
            sum_vap = sum_vap + x_vap_new(j)
            ! check if the stepsize is too large.
            ! The 1st condition: (x(j)_vap - z(j))*(x(j)_liq - z(j)) < 0  ( z is inbetween xliq and xvap)
            ! The 2nd condition: x(j)_vap AND x(j)_liq > 0
            ! have to be fulfilled for all components
            ! Otherwise the stepsize is reduced
            !if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0)) then

            !------------------------------------------------------------------------
            !UNCOMMENTED FOR THE CALCULATION OF VLc EQUILIBRIA!!! Andreas July 2011
            !                if (abs(x_vap_new(j) - x_liq_new(j)) < 1.D-4) then
            !                    x_vap_new(j) = 1 - x_liq_new(j)
            !                end if
            !                if (((x_vap_new(j) - z(j))*(x_liq_new(j) - z(j)) > 0.D0) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0) &
            !                    &   .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq_new(j) > 1.D0)) then
            !        if ((x_liq_new(j)*x_vap_new(j) < 0.D0) .OR. (x_vap_new(j) > 1.D0) .OR. (x_liq_new(j) > 1.D0)) then
            ! took the criteria from flash_pt_calc!  - J.G. 2012-10
            if (((x_vap_new(j) - x_known(j))*(x_liq_new(j) - x_known(j)) > -1.D-10) .OR. (x_liq_new(j)*x_vap_new(j) < 0.D0) &
                & .OR. (x_vap_new(j) >= 1.D0) .OR. (x_liq_new(j) >= 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                & .OR. (x_liq_new(j) <= 0.D0)) then          ! addtional check included for x <= 0  J.G., 09.2012
                !------------------------------------------------------------------------
                stepsize = 2.D0
                Delta_X = Delta_X/stepsize
                if (maxval(dabs(Delta_X)) < 1.D-10) then
                    errval = -3333
                    return
                end if
                stepsize = 1.D0
                j = 0
                sum_liq = 0.D0
                sum_vap = 0.D0
            end if


            j = j + 1
        end do
        ! --------------------------------------------------------------------------------


        ! write the new values to the variables
        x_vap = x_vap_new!/sum_vap
        x_liq = x_liq_new!/sum_liq
        Temp = Temp_new
        press = press_new

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del
        max_del = 0.D0
        Do k = 1, j - 1
            if(abs(delta_X(k) / Var_X(k)) > max_del) then
                max_del = abs(delta_X(k) / Var_X(k))
            end if
        end do

        if ((max_del) < eps_del) return

        !Catch unphysical temperatures
        if ((temp < 0.D0) .or. temp > 1000.D0) then
            errval = -4321
            exit
        end if
    end do

    ! Iteration failed!
    if (i > 30) errval = -2222

    end subroutine Flash_ps
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_ps(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, s_spec, iPhase_try, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS.
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PS-FLASH:     p, s AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   s_spec      - Specified (given) entropy value
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger October 2011








    implicit none

    type(type_gl) :: gl


    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est, s_spec
    integer::  iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq, dxk_dxj_vap, dxk_dxj_liq,  s_liq, s_vap, beta_s !S_CALC
    integer:: i, j, errorflag, reasonable_beta

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:
    call newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    if (errorflag /= 0) then
        errval = errorflag
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if
    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = gl%rho_vap
    call lnf_mix(gl,T, d_vap, p, lnfi_vap)
    s_vap = S_CALC(gl,T,d_vap, 0)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = gl%rho_liq
    call lnf_mix(gl,T, d_liq, p, lnfi_liq)
    s_liq = S_CALC(gl,T,d_liq, 0)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        ! set the module variables back to original values
        ! Andreas Jun 2013
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        gl%molfractions = z
        return
    end if

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp-1
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do
    GibbsEQN(gl%ncomp) = lnfi_vap(gl%ncomp) - lnfi_liq(gl%ncomp)

    ! for the case of a p,h-flash another n-2 set of equations
    ! is needed (comes from the mass balance):
    do j = 1, gl%ncomp-2
        !calculate the derivative d(xk)/(d(xj) for vapor and liquid phase, k = ncomp-1
        dxk_dxj_vap = (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/(x_vap(j) - z(j))
        dxk_dxj_liq = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/(x_liq(j) - z(j))
        GibbsEQN(gl%ncomp + j) = dxk_dxj_vap - dxk_dxj_liq
    end do

    ! In case of a ph-flash, an additional equation is needed: h(T,p,x) - h_specified = 0
    beta_s = 0.D0
    vapfrac = 0.D0
    reasonable_beta = 0
    do j = 1, gl%ncomp-1
        beta_s = (z(i) - x_liq(i))/(x_vap(i) - x_liq(i))
        !if ((beta_s > 0.D0) .AND. (beta_s  < 1.D0)) then
        vapfrac = vapfrac + (z(i) - x_liq(i))/(x_vap(i) - x_liq(i))
        reasonable_beta = reasonable_beta + 1
        !end if
    end do
    if (reasonable_beta == 0) then
        vapfrac = 0.65D0
    else
        vapfrac = vapfrac / reasonable_beta
    End if

    GibbsEQN(2*gl%ncomp - 1) = (s_vap * vapfrac + (1.D0 - vapfrac) * s_liq) - s_spec

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_ps
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_ps (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - PS-FLASH:     p, s AND x VECTOR GIVEN
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac        - Molar vapor fraction
    ! OUTPUT:
    !   errval  - Error value
    !   JacMatrix  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! A. Jäger October 2011








    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: i, j
    double precision:: d_vap, d_liq, ds_dT_vap, ds_dT_liq,dbeta_dxj_vap, dbeta_dxj_liq, s_vap, s_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z
    double precision, dimension(30):: ds_dxj_vap, ds_dxj_liq
    double precision:: rhoredmix_orig, tredmix_orig

    JacMatrix = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    s_vap = 0.D0
    s_liq = 0.D0
    ds_dT_liq = 0.D0
    ds_dT_vap = 0.D0
    dbeta_dxj_vap = 0.D0
    dbeta_dxj_liq = 0.D0
    ds_dxj_vap = 0.D0
    ds_dxj_liq = 0.D0

    !--------------------------------------------------------------------------
    ! p-h flash calculation, p and h and the overall x vector are given.
    ! x',x" and T need to be calculated
    !--------------------------------------------------------------------------
    ! get the liquid phase properties
    d_liq = gl%rho_liq
    z = gl%molfractions
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    call dlnfi_dxj_TP (gl,T, d_liq, dlnfidXj_liq)
    call dlnphii_dT(gl,T, d_liq, dlnphiidT_liq)
    call ds_dx_TP (gl,T, d_liq, ds_dxj_liq)
    ds_dT_liq = ds_dT_px (gl,T, d_liq)
    s_liq = S_CALC(gl,T, d_liq, 0)
    ! get the vapor phase properties
    d_vap = gl%rho_vap
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    call dlnfi_dxj_TP (gl,T, d_vap, dlnfidXj_vap)
    call dlnphii_dT(gl,T, d_vap, dlnphiidT_vap)
    call ds_dx_TP (gl,T, d_vap, ds_dxj_vap)
    ds_dT_vap = ds_dT_px (gl,T, d_vap)
    s_vap = S_CALC(gl,T, d_vap, 0)
    do j = 1, (gl%ncomp-1)   ! goes through the 2*(n-1) columns, which are the derivativs to x_j
        do i = 1, gl%ncomp   ! goes through the first n-1 rows, which are the functions F_i
            JacMatrix(j, i) = dlnfidXj_vap(j, i)
            ! goes through the columns j = n to 2*(n-1) and the rows i = 1 to n-1, which are the derivs d(Fi)/d(x_liq(j)
            JacMatrix(gl%ncomp-1 + j, i) = - dlnfidXj_liq(j, i)
        end do
        dbeta_dxj_liq = (z(j) - x_vap(j)) / (x_vap(j) - x_liq(j))**2 / (gl%ncomp-1)    !Added the /ncomp Andreas Jun 2013
        dbeta_dxj_vap = (x_liq(j) - z(j)) / (x_vap(j) - x_liq(j))**2 / (gl%ncomp-1)    !Added the /ncomp Andreas Jun 2013
        JacMatrix(j,2*gl%ncomp - 1) = dbeta_dxj_vap * (s_vap - s_liq) + vapfrac * ds_dxj_vap(j)             !Goes through the columns 1 to n-1 in the 2n-1th row, which are the derivatives of h w.r.t. x in the vapor phase
        JacMatrix(gl%ncomp - 1 +j,2*gl%ncomp - 1) = dbeta_dxj_liq * (s_vap - s_liq) + (1.D0 - vapfrac) * ds_dxj_liq(j)  !Goes through the columns n to n-2 in the 2n-1th row, which are the derivatives of h w.r.t. x in the liquid phase
    end do
    do j = 1, gl%ncomp-2
        JacMatrix(j, gl%ncomp+j) = - (x_vap(gl%ncomp-1) - z(gl%ncomp-1))/((x_vap(j) - z(j))**2)   !Ableitung von F_n+j nach x_j_vap
        JacMatrix(gl%ncomp-1, gl%ncomp+j) = 1.D0/(x_vap(j) - z(j))    !Ableitung von f_n+j nach x_n-1_vap
        JacMatrix(gl%ncomp-1+j, gl%ncomp+j) = (x_liq(gl%ncomp-1) - z(gl%ncomp-1))/((x_liq(j) - z(j))**2)    !Ableitung von F_n+j nach x_j_liq
        JacMatrix(2*(gl%ncomp-1), gl%ncomp+j) = - 1.D0/(x_liq(j) - z(j))  !Ableitung von f_n+j nach x_n-1_liq
    end do
    do i = 1, gl%ncomp
        JacMatrix(2*gl%ncomp - 1, i) = dlnphiidT_vap(i) - dlnphiidT_liq(i) !Derivatives of ln(fi_vap) - ln(fi_vap) with respect to T at const. p and x
    end do
    if (gl%ncomp > 2) then
        do i = gl%ncomp + 1, 2*gl%ncomp - 2
            JacMatrix(2*gl%ncomp - 1, i) = 0.D0 !Derivatives of Fn+j w.r.t. T
        end do
    end if
    JacMatrix(2*gl%ncomp -1, 2*gl%ncomp -1) = vapfrac * ds_dT_vap + (1.D0-vapfrac) * ds_dT_liq
    !--------------------------------------------------------------------------

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_ps
    !**************************************************************************

    !native of branch, lead to errors, copied from trunk, insersted below
    !!**************************************************************************
    !Subroutine RachRice (gl,K_Val, x_known, x_vap, x_liq, vapfrac, errorflag)
    !!**************************************************************************
    !! THIS SUBROUTINE SOLVES THE RACHFORD-RICE EQUATION Sum(xi'' - xi') = 0 WHERE
    !! xi'' = Ki * xi / (1 - vapfrac + vapfrac * Ki) AND  xi' =  xi / (1 - vapfrac + vapfrac * Ki)
    !! FOR the vapor fraction vapfrac
    !!--------------------------------------------------------------------------
    !! Variables
    !!
    !!   K_val       - Vector containing the K_values
    !!   x_known     - The given overall composition vector
    !!   x_vap       - Return the composition of the vapor phase
    !!   x_liq       - Return the composition of the liquid phase
    !!   vapfrac        - vapor fraction (Input is startvalue, output the result of this routine)
    !!
    !!--------------------------------------------------------------------------
    !! Johannes and Andreas Aug 2011

    !

    !implicit none
    !
    !type(type_gl) :: gl
    !
    !
    !double precision, dimension(30) :: x_known, K_val
    !double precision, dimension(30) :: x_vap, x_liq
    !double precision :: vapfrac
    !
    !double precision:: sumx_bubble, sumx_dew, sumx_2phase
    !double precision:: frac_min, frac_max, frac_min_allowed, frac_max_allowed, frac
    !double precision:: Delta_allowed
    !double precision, dimension(65)::Parameters
    !integer:: i, frac_type, Max_iterations, Iterations, errorflag
    !
    !x_vap = 0.d0
    !x_liq = 0.d0
    !sumx_bubble = 0.D0
    !sumx_dew = 0.D0
    !sumx_2phase = 0.D0
    !
    !!Equations: Gernert et al. (2014)
    !Do i = 1, gl%Ncomp
    !    !sumx_bubble is the Rachford Rice equation with the assumption vapfrac = 0
    !    sumx_bubble = sumx_bubble + x_known(i) * K_val(i)   !Eq. (15)
    !    !sumx_dew is the Rachford Rice equation with the assumption vapfrac = 1
    !    sumx_dew = sumx_dew + x_known(i) / K_val(i)         !Eq. (17)
    !    !sumx_2phase is the Rachford Rice equation with the assumption vapfrac = 0.5
    !    sumx_2phase = sumx_2phase + x_known(i) * (K_val(i) - 1.D0) / (K_val(i) + 1.D0)  !Eq. (14) with beta = 0.5d0 and some weird things
    !End do
    !
    !!Gernert et al. (2014), Eq. (16)
    !If (sumx_bubble < 1.D0) then !The mixture is assumed to be at the bubble point or below
    !    vapfrac = 0
    !    Do i = 1, gl%Ncomp
    !        x_vap(i) =   x_known(i) * K_Val(i) / sumx_bubble
    !    end do
    !    x_liq = x_known
    !    return
    !end if
    !
    !!Gernert et al. (2014), Eq. (18)
    !If (sumx_dew < 1.D0) then !The mixture is assumed to be at the dew point or above
    !    vapfrac = 0
    !    Do i = 1, gl%Ncomp
    !        x_liq(i) =   x_known(i) / K_Val(i) / sumx_dew
    !    end do
    !    x_vap = x_known
    !    return
    !end if
    !
    !!The mixture is in the 2 phase region with 0.5 < vapfrac < 1.
    !!Instead of vapfrac the liquid fraction alpha is used, since vapfrac near unity might cause round off errors (see Michelsen & Mollerup)
    !If (sumx_2phase > 0.D0) then
    !    !The liquid fraction is chosen as independent variable
    !    frac_type = 1
    !else
    !    frac_type = 0
    !End if
    !Delta_allowed = 1.D-8
    !Parameters = 0.D0
    !Parameters(1:30) = x_known
    !Parameters(31:60) = K_val
    !Parameters(61) = frac_type
    !frac_min = 0.D0
    !frac_min_allowed = 0.D0
    !frac_max = 0.5D0
    !frac_max_allowed = 1.D0
    !Max_iterations = 1.d3 !50
    !frac = 0.D0
    !errorflag = 0
    !
    !!solve Gernert et al. (2014), Eq. (10) with Regular Falsi
    !call Regula_Falsi(gl,rac_func, frac, frac_min, frac_max, Delta_allowed, frac_min_allowed, frac_max_allowed, &
    !    &   Max_iterations, Iterations, Errorflag, Parameters)
    !
    !!If liquid fraction was used, change to vapor fraction
    !if(frac_type == 1) then
    !    vapfrac = 1.D0 - frac
    !else
    !    vapfrac = frac
    !end if
    !
    !!Gernert et al. (2014), Eq. (13)
    !sumx_dew = 0.d0
    !sumx_bubble = 0.d0
    !Do i = 1, gl%Ncomp
    !    x_vap(i) =   x_known(i) * K_Val(i) / (1.D0 - vapfrac + vapfrac * K_val(i))
    !    sumx_dew = sumx_dew + x_vap(i)
    !    x_liq(i) =   x_known(i) / (1.D0 - vapfrac + vapfrac * K_val(i))
    !    sumx_bubble = sumx_bubble + x_liq(i)
    !end do
    !x_vap = x_vap/sumx_dew
    !x_liq = x_liq/sumx_bubble
    !
    !End Subroutine RachRice
    !!**************************************************************************


    !**************************************************************************
    module subroutine RachRice (gl,K_Val, x_known, x_vap, x_liq, vapfrac, errorflag)
    !**************************************************************************
    ! THIS SUBROUTINE SOLVES THE RACHFORD-RICE EQUATION Sum(xi'' - xi') = 0 WHERE
    ! xi'' = Ki * xi / (1 - vapfrac + vapfrac * Ki) AND  xi' =  xi / (1 - vapfrac + vapfrac * Ki)
    ! FOR the vapor fraction vapfrac
    !--------------------------------------------------------------------------
    ! Variables
    !
    !   K_val       - Vector containing the K_values
    !   x_known     - The given overall composition vector
    !   x_vap       - Return the composition of the vapor phase
    !   x_liq       - Return the composition of the liquid phase
    !   vapfrac        - vapor fraction (Input is startvalue, output the result of this routine)
    !
    !--------------------------------------------------------------------------
    ! Johannes and Andreas Aug 2011



    implicit none

    type(type_gl) :: gl


    double precision, dimension(30) :: x_known, K_val
    double precision, dimension(30) :: x_vap, x_liq
    double precision :: vapfrac

    double precision:: sumx_bubble, sumx_dew, sumx_2phase
    double precision:: frac_min, frac_max, frac_min_allowed, frac_max_allowed, frac
    double precision:: Delta_allowed
    type(type_additional_parameters) :: parameters
    integer:: i, frac_type, Max_iterations, Iterations, errorflag

    x_vap = 0.d0
    x_liq = 0.d0
    sumx_bubble = 0.D0
    sumx_dew = 0.D0
    sumx_2phase = 0.D0

    !Equations: Gernert et al. (2014)
    Do i = 1, gl%Ncomp
        !sumx_bubble is the Rachford Rice equation with the assumption vapfrac = 0
        sumx_bubble = sumx_bubble + x_known(i) * K_val(i)   !Eq. (15)
        !sumx_dew is the Rachford Rice equation with the assumption vapfrac = 1
        sumx_dew = sumx_dew + x_known(i) / K_val(i)         !Eq. (17)
        !sumx_2phase is the Rachford Rice equation with the assumption vapfrac = 0.5
        sumx_2phase = sumx_2phase + x_known(i) * (K_val(i) - 1.D0) / (K_val(i) + 1.D0)  !Eq. (14) with beta = 0.5d0 and some weird things
    End do

    !Gernert et al. (2014), Eq. (16)
    If (sumx_bubble < 1.D0 + 1d-15) then !The mixture is assumed to be at the bubble point or below
        vapfrac = 0
        Do i = 1, gl%Ncomp
            x_vap(i) =   x_known(i) * K_Val(i) / sumx_bubble
        end do
        x_liq = x_known
        return
    end if

    !Gernert et al. (2014), Eq. (18)
    If (sumx_dew < 1.D0 + 1d-15) then !The mixture is assumed to be at the dew point or above
        vapfrac = 0
        Do i = 1, gl%Ncomp
            x_liq(i) =   x_known(i) / K_Val(i) / sumx_dew
        end do
        x_vap = x_known
        return
    end if

    !The mixture is in the 2 phase region with 0.5 < vapfrac < 1.
    !Instead of vapfrac the liquid fraction alpha is used, since vapfrac near unity might cause round off errors (see Michelsen & Mollerup)
    If (sumx_2phase > 0.D0) then
        !The liquid fraction is chosen as independent variable
        frac_type = 1
    else
        frac_type = 0
    End if
    Delta_allowed = 1.D-8
    !Parameters = 0.D0
    parameters%a_p(1:30) = x_known
    parameters%a_p(31:60) = K_val
    parameters%a_p(61) = frac_type
    frac_min = -1.D-16 !SH 09/18 changed the boundaries from 0.d0 and 1.d0 to -1.d-16 and 1.d0 + 1.d-16 because in release mode the old boundaries produce errors
    frac_min_allowed = -1.D-16
    frac_max = 0.5D0
    frac_max_allowed = 1.D0 + 1.d-16
    Max_iterations = 50
    frac = 0.D0
    errorflag = 0

    !solve Gernert et al. (2014), Eq. (10) with Regular Falsi
    call Regula_Falsi(gl,rac_func, frac, frac_min, frac_max, Delta_allowed, frac_min_allowed, frac_max_allowed, &
        &   Max_iterations, Iterations, Errorflag, Parameters)
    if (errorflag /= 0) then
        errorflag = -2111
    endif
    !If liquid fraction was used, change to vapor fraction
    if(frac_type == 1) then
        vapfrac = 1.D0 - frac
    else
        vapfrac = frac
    end if

    !Gernert et al. (2014), Eq. (13)
    sumx_dew = 0.d0
    sumx_bubble = 0.d0
    Do i = 1, gl%Ncomp
        x_vap(i) =   x_known(i) * K_Val(i) / (1.D0 - vapfrac + vapfrac * K_val(i))
        sumx_dew = sumx_dew + x_vap(i)
        x_liq(i) =   x_known(i) / (1.D0 - vapfrac + vapfrac * K_val(i))
        sumx_bubble = sumx_bubble + x_liq(i)
    end do
    x_vap = x_vap/sumx_dew
    x_liq = x_liq/sumx_bubble

    End Subroutine RachRice
    !**************************************************************************

    !old one of branch
    !!**************************************************************************
    !double precision function rac_func(gl,frac, parameters)
    !!**************************************************************************
    !! THIS SUBROUTINE IS!! THE RACHFORD-RICE EQUATION Sum(xi'' - xi') = 0 WHERE
    !! xi'' = Ki * xi / (1 - vapfrac + vapfrac * Ki) AND  xi' =  xi / (1 - vapfrac + vapfrac * Ki)
    !! FOR the vapor fraction vapfrac or the liquid fraction alpha
    !!--------------------------------------------------------------------------
    !! Variables
    !!
    !!   K_val       - Vector containing the K_values
    !!   x_known     - The given overall composition vector
    !!   frac        - The vapor or liquid fraction
    !!   frac_type   - Indicates which fraction is given (0 -- >  vapfrac, 1 -- >  alpha)
    !!
    !!--------------------------------------------------------------------------
    !! Johannes and Andreas Aug 2011
    !
    !!Gernert et al. (2014), Eq. (14)
    !

    !

    !implicit none
    !
    !type(type_gl) :: gl
    !
    !
    !double precision, dimension(30):: x_known, K_val
    !double precision:: frac
    !integer:: frac_type, i
    !
    !Double Precision, dimension(65):: Parameters            ! Inputvariables
    !
    !x_known = Parameters(1:30)
    !K_val = Parameters(31:60)
    !frac_type = int(Parameters(61))
    !
    !rac_func = 0.D0
    !
    !if(frac_type == 0) then
    !    Do i = 1, gl%ncomp
    !        rac_func = rac_func + x_known(i) * (K_val(i) - 1.D0) / (1.D0 - frac + frac * K_val(i))
    !    End do
    !else
    !    Do i = 1, gl%ncomp
    !        rac_func = rac_func + x_known(i) * (K_val(i) - 1.D0) / (frac + (1.D0 - frac) * K_val(i))
    !    End do
    !end if
    !
    !End function rac_func
    !!**************************************************************************


    !**************************************************************************
    double precision module function rac_func(gl,frac, parameters)
    !**************************************************************************
    ! THIS SUBROUTINE IS!! THE RACHFORD-RICE EQUATION Sum(xi'' - xi') = 0 WHERE
    ! xi'' = Ki * xi / (1 - vapfrac + vapfrac * Ki) AND  xi' =  xi / (1 - vapfrac + vapfrac * Ki)
    ! FOR the vapor fraction vapfrac or the liquid fraction alpha
    !--------------------------------------------------------------------------
    ! Variables
    !
    !   K_val       - Vector containing the K_values
    !   x_known     - The given overall composition vector
    !   frac        - The vapor or liquid fraction
    !   frac_type   - Indicates which fraction is given (0 -- >  vapfrac, 1 -- >  alpha)
    !
    !--------------------------------------------------------------------------
    ! Johannes and Andreas Aug 2011

    !Gernert et al. (2014), Eq. (14)




    implicit none

    type(type_gl) :: gl


    double precision, dimension(30):: x_known, K_val
    double precision:: frac
    integer:: frac_type, i

    type(type_additional_parameters) :: Parameters            ! Inputvariables

    x_known = parameters%a_p(1:30)
    K_val = parameters%a_p(31:60)
    frac_type = int(parameters%a_p(61))

    rac_func = 0.D0

    if(frac_type == 0) then
        Do i = 1, gl%ncomp
            rac_func = rac_func + x_known(i) * (K_val(i) - 1.D0) / (1.D0 - frac + frac * K_val(i))
        End do
    else
        Do i = 1, gl%ncomp
            rac_func = rac_func + x_known(i) * (K_val(i) - 1.D0) / (frac + (1.D0 - frac) * K_val(i))
        End do
    end if

    End function rac_func
    !**************************************************************************



    !**************************************************************************
    double precision module function Tsat_iter (gl,P, x, IPhase)
    !**************************************************************************
    ! FUNCTION FOR THE ITERATIVE ESTIMATION OF THE BUBBLE POINT OR DEW POINT
    ! TEMPERATURE USING THE WILSON K-FACTOR APPROXIMATION AND THE NEWTON METHOD
    !--------------------------------------------------------------------------
    !   j. Gernert, Jan. 2011





    implicit none

    type(type_gl) :: gl


    double precision :: T_old, Zerofunc, Deriv, delta
    double precision :: P, T_new
    double precision:: Tmin, Tmax, Tmin_allowed, Tmax_allowed, Delta_allowed
    type(type_additional_parameters) :: parameters
    double precision, dimension(30):: x
    double precision, dimension(30):: K
    integer:: IPhase, i, Max_iterations, errorflag, Iterations
    integer:: count

    double precision:: T1, T2, racrice1, racrice2, racrice_new

    ! Set start values
    T_new = 0
    K = 0.D0
    Zerofunc = -1.D0
    Deriv = 0.D0
    delta = 1
    count = 0

    ! create a start value for T
    T_old = 0.d0
    do i = 1, gl%ncomp
        T_old = T_old + (gl%ttp(i) + gl%tc(i))/2.d0*gl%molfractions(i)
    end do

    !iteration loop
    do while (abs(delta) > 0.001)
        do i = 1, gl%ncomp
            ! calculate K values
            K(i) = gl%pc(i)/p*dexp(5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/T_old))
            if (IPhase == 1) then
                ! The function F:= Sum(x(i)") - 1 = 0 has to be fulfilled
                Zerofunc = Zerofunc + K(i)*x(i)
                ! Calculate dF/dT
                Deriv = Deriv + K(i)*x(i)*5.373D0*(1.D0 + gl%accen(i))*gl%Tc(i)/(T_old**2)
            else
                !Error handling required if K gets 0 or negative -- >  Iteration failed
                !Andreas March. 2013
                if (K(i) < 1.d-12) then
                    T_new = -1.D0 !Dummy value
                    exit
                end if
                ! The function F:= Sum(x(i)") - 1 = 0 has to be fulfilled
                Zerofunc = Zerofunc + x(i)/K(i)
                ! Calculate dF/dT
                Deriv = Deriv + x(i)/K(i)*(-5.373D0)*(1.D0 + gl%accen(i))*gl%Tc(i)/(T_old**2)
            end if
        end do
        !Error handling required if K gets 0 or negative -- >  Iteration failed
        !Andreas March. 2013
        if (T_new < 0.D0) exit
        !Newton-Raphson Method
        T_new = T_old - Zerofunc/Deriv
        delta = Zerofunc
        T_old = T_new
        if (T_new < -1.d-1) exit
        Zerofunc = -1.D0
        Deriv = 0.D0
        count = count + 1
        !Andreas August 2012: increased the maximum number of iterations to 60, since for a hydrogen rich mixture (fluidl = 'methane;ethane;propane;ammonia;hydrogen' molesl = '0.2;0.04;0.01;0.0012;0.7488' at 1Mpa) 33 iteration are necessary to converge
        if (count > 60) exit
    end do

    if ((T_new < 0.d0) .OR. (count > 60)) then
        !parameters = 0.d0
        Delta_allowed = 1.d-6
        Max_Iterations = 30
        parameters%a_p(1) = p
        parameters%a_p(2) = IPhase
        Tmin = 1.d6
        do i = 1, gl%ncomp
            if (Tmin > gl%ttp(i)) Tmin = gl%ttp(i)
        end do
        ! Andreas August 2012:
        ! The regula falsi routine found a wrong solution for the mixture: fluidl = 'methane;ethane;propane;ammonia;hydrogen' molesl = '0.2;0.04;0.01;0.0012;0.7488' at 1Mpa
        ! Problem yet unsolved
        Tmax = maxval(gl%tc)*1.2d0
        Tmin = Tmin*0.8d0
        Tmin_allowed = Tmin*0.6
        Tmax_allowed = Tmax*1.5
        call Regula_Falsi(gl,RacRice_div, T_new, Tmin, Tmax, Delta_allowed, &
            Tmin_allowed, Tmax_allowed, Max_Iterations, Iterations, errorflag, parameters)


        !Check whether the Regula_falsi converged
        Racrice1 = RacRice_div(gl,T_new, Parameters)
        if ((errorflag /= 0) .or. (dabs(racrice1) > 1.D-6)) then
            ! Andreas Jul 2013
            ! For the mixture nitrogen and toluol with xN2 = 0.999934 and xTol = 0.000066 the Regula_falsi claims to have found a solution at the upper interval limit.
            ! This is not the correct solution but occurs from RacRice_div giving a very steep slope.
            ! Thus another try was implemented here: If the Regula_falsi fails, try to decrease the interval "manually" by a few steps of the bisection method
            ! 1) Take original starting values
            T2 = maxval(gl%tc)*1.2d0
            T1 = Tmin*0.8d0
            ! 2) calculate Rac Rice at upper and lower limit
            parameters%a_p(1) = p
            parameters%a_p(2) = IPhase
            Racrice1 = RacRice_div(gl,T1, Parameters)
            Racrice2 = RacRice_div(gl,T2, Parameters)
            if (Racrice2 * Racrice1 > 0.D0) then
                !Root not in the interval, Quit with error
                errorflag = -4321
                Tsat_iter = -1.D0
                return
            end if

            ! 3) Start bisection
            Do i = 1, 30
                T_new= (T1 + T2) / 2.D0
                Racrice_new = RacRice_div(gl,T_new, Parameters)
                if (racrice_new * racrice1 < 0.D0) then
                    T2 = T_new
                else
                    T1 = T_new
                    racrice1 = racrice_new
                end if
                if(dabs(racrice_new) < 1.D-1) exit
            End do

            Tmax = T2
            Tmin = T1
            Tmin_allowed = Tmin*0.6
            Tmax_allowed = Tmax*1.5
            call Regula_Falsi(gl,RacRice_div, T_new, Tmin, Tmax, Delta_allowed, &
                Tmin_allowed, Tmax_allowed, Max_Iterations, Iterations, errorflag, parameters)

            !Check if Regula Falsi converged
            Racrice1 = RacRice_div(gl,T_new, Parameters)
            if ((errorflag /= 0) .or. (dabs(Racrice1) > 1.D-6)) then
                T_new = -1.D0
            end if

        end if

    end if

    Tsat_iter = T_new

    end function Tsat_iter
    !**************************************************************************

    !**************************************************************************
    Double Precision module Function RacRice_div(gl,temp, Parameters)
    !**************************************************************************
    ! Modified to check, if the exponent does not create double precision overflow - J. Gernert, A.Jäger, 08.2012


    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: temp, p              ! Inputvariablen
    Double Precision :: sum_calc, exponent
    integer::IPhase, i
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    IPhase = int(parameters%a_p(2))

    sum_calc = 0.d0
    if (IPhase == 1) then
        do i = 1, gl%ncomp
            !Changed by Andreas Feb. 2012
            exponent = 5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/temp)
            if (exponent > 50.d0) exponent = 50.d0
            sum_calc = sum_calc + gl%molfractions(i)*(gl%pc(i)/p*dexp(exponent))
            !sum_calc = sum_calc + molfractions(i)*dlog(pc(i)/p*dexp(5.373D0*(1.D0 + accen(i))*(1.D0 - Tc(i)/temp)))

        end do
    else
        do i = 1, gl%ncomp
            !sum_calc = sum_calc + molfractions(i) / K_val(i)
            !Changed by Andreas Feb. 2012
            exponent = 5.373D0*(1.D0 + gl%accen(i))*(1.D0 - gl%Tc(i)/temp)
            if (exponent < -50.d0) exponent = -50.d0
            sum_calc = sum_calc + gl%molfractions(i)/(gl%pc(i)/p*dexp(exponent))
            !sum_calc = sum_calc + molfractions(i)/dlog(pc(i)/p*dexp(5.373D0*(1.D0 + accen(i))*(1.D0 - Tc(i)/temp)))

        end do
    end if

    RacRice_div = (1.d0 - sum_calc)

    End Function RacRice_div
    !**************************************************************************

    !**************************************************************************
    module subroutine newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, iPhase_try, errorflag)
    !**************************************************************************
    ! CALCULATES DENSITIES FOR THE PURE FLUIDS AT GIVEN TEMPERATURE AND PRESSURE
    ! AND MIXTURE DENSITIES FOR THE VAPOR AND LIQUID PHASE
    ! THE RESULTS ARE STORED IN A MODULE
    ! INPUTS:
    !   P               -   Pressure in MPa
    !   T               -   Temperature in K
    !   x_vap           -   Composition of vapor phase (THIS MIGHT BE A SECOND LIQUID COMPOSITION IN SOME CASES)
    !   x_liq           -   Composition of the liquid phase
    !   rho_vap_est     -   Estimate for vapor phase density
    !   rho_liq_est     -   Estimate for the liquid phase density
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !
    ! OUTPUTS:
    !   errorflag       -   Indicates if an error occured during the calculations
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan 2011
    ! A. Jäger, Modified August 2011






    implicit none

    type(type_gl) :: gl


    double precision:: T, P, rhovap_est, rholiq_est
    integer:: iPhase_try
    double precision, dimension(30):: x_vap, x_liq, x_liq_reac, x_vap_reac
    integer:: errorflag, iPhaseVap, iPhaseLiq
    double precision:: rhoredmix_orig, tredmix_orig
    double precision, dimension(30):: z
    type(type_additional_parameters) :: parameters
    logical :: react



    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    errorflag = 0
    react = .false.
    !parameters = 0.d0
    !--------------------------------------------------------------------------
    ! calculate the liquid and vapor phase density for the mixture and save it in the
    ! module variables "rho_vap" and "rho_liq" in the module "module_VLE"
    !--------------------------------------------------------------------------
    z = 0.D0

    select case (Iphase_try)
    case(0)     ! phase distribution unknown
        iPhaseVap = 0
        iPhaseLiq = 0
    case(1)     ! LLE assumed, light phase fixed
        iPhaseVap = 1
        iPhaseLiq = 0
    case(2)     ! VLE assumed, vapor phase fixed
        iPhaseVap = 2
        iPhaseLiq = 0
    case(3)     ! LLE of VLE, Liquid phase fixed
        iPhaseVap = 0
        iPhaseLiq = 1
    case(4)     ! LLE assumed, both phases fixed
        iPhaseVap = 1
        iPhaseLiq = 1
    case(5)     ! VLE assumed, both phases fixed
        iPhaseVap = 2
        iPhaseLiq = 1
    case(6)     ! VLE assumed, phases switched
        iPhaseVap = 1
        iPhaseLiq = 2
    end select

    ! Since the Phase is not known (LLE possible!) both liquid and vapor densities have to be checked!

    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! calculate the gas phase density
    ! Same call in both branches of if-statement. Not needed --> commented, Andreas May 2016
    !if (mix_type == 1) then !Helmholtz equations
    gl%rho_vap = rhomix_calc(gl,T, P, rhovap_est, iPhaseVap, 0)
    if (gl%rho_vap < 1.d-14) then
        gl%rho_vap = rhomix_calc(gl,T, P, rhovap_est, 0, 0)
    endif
    !else               !SRK
    !    rho_vap = rhomix_calc(gl,T, P, rhovap_est, iPhaseVap, 0)
    !!    rho_vap = rhomix_calc(gl,T, P, rhovap_est, 2,0)
    !End if


    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! calculate the liquid phase density
    ! Same call in both branches of if-statement. Not needed --> commented, Andreas May 2016
    !if (mix_type == 1) then !Helmholtz equations
    if(gl%seawaterflash) gl%seacalc = .true.
    gl%rho_liq = rhomix_calc(gl,T, P, rholiq_est, iPhaseLiq, 0)
    gl%seacalc = .false.
    if (gl%rho_liq < 1.d-14) then
        if(gl%seawater) gl%seacalc = .true.
        gl%rho_liq = rhomix_calc(gl,T, P, rholiq_est, 0, 0)
        gl%seacalc = .false.
    endif
    !else               !SRK
    !    rho_liq = rhomix_calc(gl,T, P, rholiq_est, iPhaseLiq, 0)
    !    rho_liq = rhomix_calc(gl,T, P, rholiq_est, 1,0)
    !End if


    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    !Replaced "==" "<" comparison Andreas, August 2016
    !if ((rho_liq == 0.D0) .OR. (rho_vap == 0.D0)) then
    !when no root is found, rho_liq is 0, thus it IS an error, but not catched
    !if ((rho_liq < 0.D0) .OR. (rho_vap < 0.D0)) then
    if ((gl%rho_liq <= 1.D-16) .OR. (gl%rho_vap <= 1.D-16)) then
        if (p > 0.D0) then   !Only for positive pressures throw error here, for negative pressure negative densities are allowed for checking the stability of critical points at negative pressures
            errorflag = -8888
        end if
        !write (*,*) 'error ', errorflag, ' -- mixture density iteration failed in subroutine "newvars"'
        return
    end if

    !If the two densities found are very close recalculate the densities with iphase = 1 and iphase = 2
    !if ((abs(rho_liq - rho_vap) / rho_liq) < 0.01D0) then
    !    IphaseLiq = 1
    !    ! calculate the liquid phase density
    !    rho_liq = rhomix_calc(gl,T, P, 0.d0, IphaseLiq,0)
    !
    !    IPhaseVap = 2
    !    z = molfractions
    !    molfractions = x_vap
    !    call reduced_parameters_calc(gl,T)
    !    ! calculate the gas phase density
    !    rho_vap = rhomix_calc(gl,T, P, 0.d0, IPhaseVap,0)
    !end if

    !This is checked two times -> not necessary here, Andreas Jäger, August 2016
    !if ((rho_liq == 0.D0) .OR. (rho_vap == 0.D0)) then
    !    errorflag = -8888
    !    !write (*,*) 'error ', errorflag, ' -- mixture density iteration failed in subroutine "newvars"'
    !    return
    !end if

    !What is this if-statement good for? Commented, Andreas Jäger, August 2016
    !if ((rho_liq > 1.D5) .OR. (rho_vap > 1.D5)) then
    !    errorflag = -8888
    !    !write (*,*) 'error ', errorflag, ' -- mixture density iteration failed in subroutine "newvars"'
    !    return
    !end if

    ! set the module variables back to original values
    errorflag = 0
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine newvars
    !**************************************************************************


    !**************************************************************************
    module subroutine LUdecomp (gl,n,aMatrix,cMatrix,ierr,herr)
    !**************************************************************************
    ! subroutine to solve a linear system of quations
    !                               X*A = C
    ! for the unknown vector X
    ! written by E.Lemmon
    ! adapted by J. Gernert, Jan. 2011
    !--------------------------------------------------------------------------
    ! Variables:
    !   n       - Number of components in the mixture
    !   aMatrix - the matrix A
    !   cMatrix - the vector C
    !--------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    !max number of components in mixture
    double precision, dimension(60, 60):: amatrix
    double precision, dimension(60):: cmatrix,  ctemp, sdecomp
    integer, dimension(60)::iord
    character(255):: herr
    integer:: ierr, i, j, k, n
    double precision:: sum
    !dimension i,j,k,sum
    ierr=0
    do i=1,n
        iord(i)=i
        sdecomp(i)=Abs(aMatrix(1,i))
        do j=2,n
            If (Abs(aMatrix(j,i)) > sdecomp(i)) sdecomp(i)=Abs(aMatrix(j,i))
        end do
        !Singular matrix
        If (sdecomp(i) == 0) Then
            ierr=1
            herr='Singular matrix'
            RETURN
        end if
    end do

    j=1
    Call Pivot(gl,n,j,iord,aMatrix,sdecomp)
    If (aMatrix(1,iord(1)) == 0) Then
        ierr=1
        RETURN
    end if

    do j=2,n
        aMatrix(j,iord(1))=aMatrix(j,iord(1))/aMatrix(1,iord(1))
    end do
    do j=2,n-1
        do i=j,n
            sum=0
            do k=1,j-1
                sum=sum+aMatrix(k,iord(i))*aMatrix(j,iord(k))
            end do
            aMatrix(j,iord(i))=aMatrix(j,iord(i))-sum
        end do
        Call Pivot(gl,n,j,iord,aMatrix,sdecomp)
        do k=j+1,n
            sum=0
            do i=1,j-1
                sum=sum+aMatrix(i,iord(j))*aMatrix(k,iord(i))
            end do
            If (aMatrix(j,iord(j)) == 0) aMatrix(j,iord(j))=1E+20
            aMatrix(k,iord(j))=(aMatrix(k,iord(j))-sum)/aMatrix(j,iord(j))
        end do
    end do
    sum=0
    do k=1,n-1
        sum=sum+aMatrix(k,iord(n))*aMatrix(n,iord(k))
    end do
    aMatrix(n,iord(n))=aMatrix(n,iord(n))-sum
    If (aMatrix(n,iord(n)) == 0) aMatrix(n,iord(n))=1E+20

    cMatrix(iord(1))=cMatrix(iord(1))/aMatrix(1,iord(1))
    do i=2,n
        sum=0
        do j=1,i-1
            sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        end do
        cMatrix(iord(i))=(cMatrix(iord(i))-sum)/aMatrix(i,iord(i))
    end do

    do i=n-1,1,-1
        sum=0
        do j=i+1,n
            sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        end do
        cMatrix(iord(i))=cMatrix(iord(i))-sum
    end do
    do i=1,n
        ctemp(i)=cMatrix(iord(i))
    end do
    do i=1,n
        cMatrix(i)=ctemp(i)
    end do

    END subroutine LUdecomp
    !**************************************************************************

    !**************************************************************************
    module subroutine Pivot (gl,n,j,iord,aMatrix,sdecomp)
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60)::amatrix
    integer, dimension(60):: iord
    double precision, dimension(60):: sdecomp

    integer:: ipivt, j, i, n, idummy
    double precision:: big, dummy


    ipivt=j
    big=Abs(aMatrix(j,iord(j))/sdecomp(iord(j)))
    do i=j+1,n
        dummy=Abs(aMatrix(j,iord(i))/sdecomp(iord(i)))
        If (dummy > big) Then
            big=dummy
            ipivt=i
        end if
    end do
    idummy=iord(ipivt)
    iord(ipivt)=iord(j)
    iord(j)=idummy

    END SUBROUTINE Pivot
    !**************************************************************************

    !**************************************************************************
    module subroutine Succ_Sub(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, Nr_of_iter, converged)
    !**************************************************************************
    ! Solution method "Successive Substitution" as published by
    !           Prausnitz and Chueh
    !           "Computer Calculations for high-pressure vapor-liquid equilibria."
    !           Prentice-Hall, Englewood Cliffs, New Jersey, 1968
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN --  iFlash = 5
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   Nr_of_iter  - How many calculation loops should be carried out?
    ! OUTPUT:
    !   errval      - Error value
    !   Depending on the iFlash chosen:
    !       iFlash = 1:     T and x_vap
    !       iFlash = 2:     T and x_liq
    !       iFlash = 3:     p and x_vap
    !       iFlash = 4:     p and x_liq
    !       iFlash = 5:     x_liq and x_vap
    !--------------------------------------------------------------------------
    ! A. Jäger Feb. 2011
    ! Modified the exit criteria, so that the module variables rhoredmix etc. are written back to their original values - J.Gernert, A.Jäger, 08.2012








    implicit none

    type(type_gl) :: gl


    double precision:: Temp, press, vapfrac
    double precision:: rhovap_est, rholiq_est
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: iFlash, Nr_of_iter
    integer:: errval
    logical:: converged

    double precision :: d_vap, d_liq, sum_vap, sum_liq, press_new
    double precision :: del_sum, temp_orig, int_fac
    double precision :: rhoredmix_orig, tredmix_orig, eps_conv, Diff_fug
    double precision, dimension(30):: lnfi_liq, lnfi_vap, fugcoef_vap, fugcoef_liq, z, K_val
    double precision, dimension(3):: temp1, sum_vap1, sum_liq1
    integer:: i, j, k, count, errorflag
    integer :: Iphase_try  !Indicates how many phases are present

    errval = 0
    !Save the original reducing parameters
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    z = x_known
    converged = .false.
    Iphase_try = 0
    lnfi_liq = 0.d0
    lnfi_vap = 0.d0

    !convergence criterion
    eps_conv = 1.D-8

    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        if (iFlash == 5) then
            call PTX_startvals_pT (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, errval)
        else
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        end if
    else
        if ((iFlash == 1) .OR. (iFlash == 3)) then   !1,3 bubble point
            x_liq = x_known
        else if (iFlash < 5) then   !2,4 dew point Andy was too lazy!!!!!
            x_vap = x_known
        end if
    end if

    !SH 03/2017: errorhandling from startvalue generation was missing here
    if (errval /= 0) return

    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. T and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        temp_orig = temp
        int_fac = 0.01D0
        i = 1
        vapfrac = 0.d0
        Do while (i < Nr_of_iter)
            sum_vap = 0.D0
            del_sum = 1.D0

            !The temperature iteration is carried out by using the bisection method
            !The algorithm searches for two temperatures where the sum of the vapor
            !fractions is < 1 for one and > 1 for the other temperature.
            if (i == 1) then
                temp1(1) = Temp_orig * (1.0D0 - int_fac)
                if (temp1(1) <= 0.d0) then
                    do while (temp1(1) <= 0.d0)
                        int_fac = int_fac*0.9d0
                        temp1(1) = Temp_orig * (1.0D0 - int_fac)
                    end do
                end if
                temp = temp1(1)
            elseif(i == 2) then
                temp1(2) = temp_orig * (1.0D0 + int_fac)
                temp = temp1(2)
            else
                Temp1(3) = (Temp1(1) + Temp1(2))/2
                temp = temp1(3)
            end if

            ! write the vapor and liquid phase density and the pure fluid densities
            ! into the module variables for the new set of T, p, x_vap, and x_liq:
            call newvars (gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try,  errorflag)
            ! If an error occurs at the density iteration leave the program
            if (errorflag /= 0) then
                errval = errorflag
                exit
            end if

            ! get the liquid phase properties (rho_liq is a module variable, calculated in newvars)
            d_liq = gl%rho_liq
            z = gl%molfractions
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacities of the liquid phase:
            call lnf_mix(gl,Temp, d_liq, press, lnfi_liq)

            !Recalculate the vapor properties until the sum of the vapor fractions attains a constant value
            !Here a break criterion is needed, Andreas September 2012
            count = 0
            Do while (abs(del_sum) > 1.D-6)
                ! get the vapor phase properties
                gl%molfractions = x_vap
                call reduced_parameters_calc(gl,Temp)

                d_vap = rhomix_calc(gl,Temp, press, rhovap_est, 2,0)

                if (d_vap == 0.D0) then
                    errorflag = -8888
                    !write (*,*) 'error ', errorflag, ' -- mixture density iteration failed in subroutine "newvars"'
                    errval = errorflag
                    gl%rhoredmix = rhoredmix_orig
                    gl%tredmix = tredmix_orig
                    gl%molfractions = z
                    return
                end if

                !Calculate the fugacity coefficient of the vapor phase
                call FUGCO_CALC_MIX(gl,Temp, d_vap, fugcoef_vap)

                !The old sum over the vaporfractions is saved
                del_sum = sum_vap

                sum_vap = 0.D0
                Do j = 1, gl%NCOMP
                    x_vap(j) = dexp(lnfi_liq(j)) / press / fugcoef_vap(j)
                    sum_vap = sum_vap + x_vap(j)
                end do

                !del_sum is the new sum over the vaporfractions minus the old one
                del_sum = del_sum - sum_vap

                !Normalize
                x_vap = x_vap / sum_vap

                !Here a break criterion is needed, Andreas September 2012
                count = count + 1
                if (count > 60) then
                    errval = -2222
                    !Set module variables back, Andreas, Jan 2013
                    gl%rhoredmix = rhoredmix_orig
                    gl%tredmix = tredmix_orig
                    gl%molfractions = z
                    return
                end if
            End do

            if (i == 1) then
                sum_vap1(1) = sum_vap
            elseif(i == 2) then
                sum_vap1(2) = sum_vap
            else
                sum_vap1(3) = sum_vap
                If ((1.D0 - sum_vap1(3)) * (1.D0 - sum_vap1(2)) < 0) then
                    Temp1(1) = Temp1(3)
                else
                    Temp1(2) = Temp1(3)
                End if
            end if

            if (dabs(sum_vap - 1.D0) < 1.D-6) then
                converged = .true.
                exit
            end if

            !Check if the root is in the interval.
            If(i == 2) then
                If ((1.D0 - sum_vap1(2)) * (1.D0 - sum_vap1(1)) > 0) then
                    !The root is still outside the interval, so the interval has to be adjusted
                    int_fac = int_fac * 2.D0
                    !An error handle is needed here, since int_fac must never exeed 1!! Andreas September 2012
                    if (int_fac > 0.9D0) then
                        errval = -2222
                        !Set module variables back, Andreas, Jan 2013
                        gl%rhoredmix = rhoredmix_orig
                        gl%tredmix = tredmix_orig
                        gl%molfractions = z
                        return
                    end if
                    i = 0
                End if
            End if
            i = i + 1
        End do

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. T and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        temp_orig = temp
        int_fac = 0.01D0
        i = 1
        vapfrac = 1.d0
        Do while (i < Nr_of_iter)
            sum_liq = 0.D0
            del_sum = 1.D0

            !The temperature iteration is carried out by using the bisection method
            !The algorithm searches for two temperatures where the sum of the vapor
            !fractions is < 1 for one and > 1 for the other temperature.
            if (i == 1) then
                temp1(1) = Temp_orig * (1.0D0 - int_fac)
                temp = temp1(1)
            elseif(i == 2) then
                temp1(2) = temp_orig * (1.0D0 + int_fac)
                temp = temp1(2)
            else
                Temp1(3) = (Temp1(1) + Temp1(2))/2
                temp = temp1(3)
            end if

            ! write the vapor and liquid phase density and the pure fluid densities
            ! into the module variables for the new set of T, p, x_vap, and x_liq:
            call newvars (gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
            ! If an error occurs at the density iteration leave the program
            if (errorflag /= 0) then
                errval = errorflag
                exit
            end if

            ! get the vapor phase properties (rho_vap is a module variable, calculated in newvars)
            d_vap = gl%rho_vap
            z = gl%molfractions
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacities of the vapor phase:
            call lnf_mix(gl,Temp, d_vap, press, lnfi_vap)
            call FUGCO_CALC_MIX(gl,Temp, d_vap, fugcoef_vap)
            !Recalculate the liquid properties until the sum of the liquid fractions attains a constant value
            !Here a break criterion is needed, Andreas September 2012
            count = 0
            Do while (abs(del_sum) > 1.D-6)
                ! get the liquid phase properties
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)

                d_liq = rhomix_calc(gl,Temp, Press,rholiq_est , 1,0)

                if (d_liq == 0.D0) then
                    errorflag = -8888
                    !write (*,*) 'error ', errorflag, ' -- mixture density iteration failed in subroutine "newvars"'
                    errval = errorflag
                    gl%rhoredmix = rhoredmix_orig
                    gl%tredmix = tredmix_orig
                    gl%molfractions = z
                    return
                end if

                !Calculate the fugacity coefficient of the liquid phase
                call FUGCO_CALC_MIX(gl,Temp, d_liq, fugcoef_liq)

                !The old sum over the vaporfractions is saved
                del_sum = sum_liq

                sum_liq = 0.D0
                Do j = 1, gl%NCOMP
                    x_liq(j) = dexp(lnfi_vap(j)) / press / fugcoef_liq(j)
                    sum_liq = sum_liq + x_liq(j)
                end do

                !del_sum is the new sum over the vaporfractions minus the old one
                del_sum = del_sum - sum_liq

                !Normalize
                x_liq = x_liq / sum_liq

                !Here a break criterion is needed, Andreas September 2012
                count = count + 1
                if (count > 60) then
                    errval = -2222
                    !Set module variables back, Andreas, Jan 2013
                    gl%rhoredmix = rhoredmix_orig
                    gl%tredmix = tredmix_orig
                    gl%molfractions = z
                    return
                end if
            End do

            if (i == 1) then
                sum_liq1(1) = sum_liq
            elseif(i == 2) then
                sum_liq1(2) = sum_liq
            else
                sum_liq1(3) = sum_liq
                If ((1.D0 - sum_liq1(3)) * (1.D0 - sum_liq1(2)) < 0) then
                    Temp1(1) = Temp1(3)
                else
                    Temp1(2) = Temp1(3)
                End if
            end if

            if (dabs(sum_liq - 1.D0) < 1.D-8) then
                converged = .true.
                exit
            end if

            !Check if the root is in the interval.
            If(i == 2) then
                If ((1.D0 - sum_liq1(2)) * (1.D0 - sum_liq1(1)) > 0) then
                    !The root is still outside the interval, so the interval has to be adjusted
                    int_fac = int_fac * 2.D0
                    !An error handle is needed here, since int_fac must never exeed 1!! Andreas September 2012
                    if (int_fac > 0.9D0) then
                        errval = -2222
                        !Set module variables back, Andreas, Jan 2013
                        gl%rhoredmix = rhoredmix_orig
                        gl%tredmix = tredmix_orig
                        gl%molfractions = z
                        return
                    end if
                    i = 0
                End if
            End if
            i = i + 1
        End do

    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given. p and
        ! x" need to be calculated
        !--------------------------------------------------------------------------
        vapfrac = 0.d0
        Do i = 1, Nr_of_iter
            sum_vap = 0.D0
            press_new = 0.D0

            ! write the vapor and liquid phase density and the pure fluid densities
            ! into the module variables for the new set of T, p, x_vap, and x_liq:
            call newvars (gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
            ! If an error occurs at the density iteration leave the program
            if (errorflag /= 0) then
                errval = errorflag
                exit
            end if

            ! get the liquid phase properties (rho_liq is a module variable, calculated in newvars)
            d_liq = gl%rho_liq
            z = gl%molfractions
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacities of the liquid phase:
            call lnf_mix(gl,Temp, d_liq, press, lnfi_liq)

            ! get the vapor phase properties (rho_vap is a module variable, calculated in newvars)
            d_vap = gl%rho_vap
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacity coefficient of the vapor phase
            call FUGCO_CALC_MIX(gl,Temp, d_vap, fugcoef_vap)

            Do j = 1, gl%NCOMP
                x_vap(j) = dexp(lnfi_liq(j)) / press / fugcoef_vap(j)
                sum_vap = sum_vap + x_vap(j)
                press_new = press_new + dexp(lnfi_liq(j)) / fugcoef_vap(j)
            end do

            if ((dabs(press_new - press) < 1.D-8).AND.(dabs(sum_vap - 1.D0) < 1.D-8)) then
                converged = .true.
                press = press_new
                x_vap = x_vap / sum_vap
                exit
            end if

            !Take the new pressure
            press = press_new
            !Normieren
            x_vap = x_vap / sum_vap

        End do

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given. p and
        ! x' need to be calculated
        !--------------------------------------------------------------------------
        ! Calculation of bubble point pressure and liquid phase composition
        ! according to the algorithm described in
        ! Ding Y.P., "Accelerated Successive Substitution Schemes for Bubble-Point and Dew-Point calculations"
        ! CAN J CHEM ENG, 69(4), 1991
        vapfrac = 1.d0
        Do i = 1, Nr_of_iter
            sum_liq = 0.D0
            press_new = 0.D0
            ! write the vapor and liquid phase density and the pure fluid densities
            ! into the module variables for the new set of T, p, x_vap, and x_liq:
            call newvars (gl,press, Temp, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
            ! If an error occurs at the density iteration leave the program
            if (errorflag /= 0) then
                errval = errorflag
                exit
            end if

            ! get the liquid phase properties (rho_liq is a module variable, calculated in newvars)
            d_liq = gl%rho_liq
            z = gl%molfractions
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacities of the liquid phase:
            call lnf_mix(gl,Temp, d_liq, press, lnfi_liq)
            call FUGCO_CALC_MIX(gl,Temp, d_liq, fugcoef_liq)

            ! get the vapor phase properties (rho_vap is a module variable, calculated in newvars)
            d_vap = gl%rho_vap
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)

            !Calculate the fugacity coefficient of the vapor phase
            call lnf_mix(gl,Temp, d_vap, press, lnfi_vap)
            call FUGCO_CALC_MIX(gl,Temp, d_vap, fugcoef_vap)

            Do j = 1, gl%NCOMP
                x_liq(j) = dexp(lnfi_vap(j)) / press / fugcoef_liq(j)
                sum_liq = sum_liq + x_liq(j)
                press_new = press_new + dexp(lnfi_liq(j)) / fugcoef_vap(j)
            end do

            if ((dabs(press_new - press) < 1.D-8).AND.(dabs(sum_liq - 1.D0) < 1.D-8)) then
                converged = .true.
                press = press_new
                x_liq = x_liq / sum_liq
                exit
            end if

            !Take the new pressure
            press = press_new

            !Normieren
            x_liq = x_liq / sum_liq
        End do

    case (5)
        !--------------------------------------------------------------------------
        ! p-T flash calculation, p and T and the overall x vector are given.
        ! x' and x" need to be calculated
        ! Algorithm according to
        !   Michelsen: Thermodynamic Models: Fundamentals & Computational Aspects
        !--------------------------------------------------------------------------
        Diff_fug = 0.D0
        K_val = 0.D0

        !0) Generate initial estimates for the phase fractions and K-values
        !   Assumption: 2 phases exist

        !Calculate K values
        do i = 1, gl%ncomp
            K_val(i) = x_vap(i) / x_liq(i)
        end do

        !Start iteration process
        Do k = 1, Nr_of_iter

            !1) Calculate new K-values
            ! Recalculate the fugacity coefficients for each phase
            ! write the vapor and liquid phase density and the pure fluid densities
            ! into the module variables for the new set of T, p, x_vap, and x_liq:
            ! Assume VLE here, Andreas Feb 2015
            !call newvars (Press, Temp, x_vap, x_liq, 0.D0, 0.D0, 0, errval)
            call newvars (gl,Press, Temp, x_vap, x_liq, 0.D0, 0.D0, 5, errval)

            if (errval /= 0) then
                exit
            end if
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)
            !dens = rhomix_calc(gl,Temp, press, 0.D0, 0,0)
            call FUGCO_CALC_MIX(gl,Temp,gl%rho_vap, fugcoef_vap)

            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            !dens = rhomix_calc(gl,Temp, press, 0.D0, 0,0)
            call FUGCO_CALC_MIX(gl,Temp,gl%rho_liq, fugcoef_liq)


            !Convergence criterion: If the maximum difference of the fugacities is smaller than Gibbs_
            do i = 1, gl%ncomp
                if ((fugcoef_vap(i) <= 0.d0) .OR. (fugcoef_liq(i) <= 0.d0)) then
                    errval = -2222
                    !Set module variables back, Andreas, Jan 2013
                    gl%rhoredmix = rhoredmix_orig
                    gl%tredmix = tredmix_orig
                    gl%molfractions = z
                    return
                end if
                !logarithm because of numerical stability
                lnfi_vap(i) = dlog(fugcoef_vap(i) * x_vap(i) * press)
                lnfi_liq(i) = dlog(fugcoef_liq(i) * x_liq(i) * press)
            end do
            Diff_fug = maxval(dabs(lnfi_vap - lnfi_liq))
            if (Diff_fug < eps_conv) then
                converged = .true.
                exit
            end if
            Diff_fug = 0.D0

            !Calculate K values
            do i = 1, gl%ncomp
                K_val(i) = fugcoef_liq(i) / fugcoef_vap(i)
            end do


            !2) Solve Rachford Rice system of equations for the phase fractions vapfrac
            !   using Newton Raphson method
            call RachRice (gl,K_Val, x_known, x_vap, x_liq, vapfrac, errval)

            if (errval /= 0) then
                !Set module variables back, Andreas, Jan 2013
                gl%rhoredmix = rhoredmix_orig
                gl%tredmix = tredmix_orig
                gl%molfractions = z
                return
            end if
        end do

        !--------------------------------------------------------------------------
        case default
        errval = -1111
    end select

    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z



    end subroutine Succ_Sub
    !**************************************************************************


    !**************************************************************************
    module subroutine Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errorflag)
    !Gauß-algorithm to bring a matrix into upper diagonal form
    !This algorithm can be used to solve a linear system of equation of the
    !form A*x = b, where A is the matrix, x is the solution vector and b is
    !the right side vector of the system of equations.
    !Furthermore, this algorithm can be used to calculate the determinant of
    !a matrix
    !INPUT:
    !MatrixA: A quadratic matrix containing a maximum of 60 rows and columns
    !vectorb: The right hand side vector, containing up to 60 entries
    !RankA: How many columns and rows does the matrix have? (means: how many
    !       equations need to be solved
    !OUTPUT:
    !vectorx: Vector containing the solutions

    !Andreas, February 2016
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision, dimension(60):: vectorb
    double precision, dimension(60):: vectorx
    double precision:: det_A
    integer:: errorflag

    integer:: i,j,k,l
    double precision:: eps_zero
    double precision:: factor_lc
    double precision, dimension(60,60):: MatA_Gauss

    eps_zero = 1.D-16   !Tolerance for values close to 0
    det_A = 0.D0


    !A rank smaller than 0 does not exist. Rank > 60 is not allowed because the matrix is restricted to 60 entries -> quit with error in both cases
    if ((rankA <= 0) .or. (rankA > 60)) then
        errorflag = -3322
        return
    end if

    !If the rank is 1, the solution is simple and treated separately
    if (rankA == 1) then
        Det_A = MatrixA(1,1)
        if (dabs(MatrixA(1,1)) > eps_zero) then
            vectorx(1) = vectorb(1) / MatrixA(1,1)
            return
        end if
    end if

    !Copy Matrix A to another matrix for Gauss transformations
    do i= 1, rankA
        do j = 1, rankA
            MatA_Gauss(i,j) = MatrixA(i,j)
        end do
    end do


    !If the rank is larger than 1:
    !Bring the matrix in form of an upper triangular matrix
    do j = 1, rankA-1     !Loop over the colums of matrix A
        !When starting to look at a new column, the element on the diagonal must have a value which is not equal to 0
        !This is ensured here:
        !1) Check whether the entry is 0, if yes, go through other rows k > j and subtract a row which has not a value of 0 in column j . If the entry is not equal to 0, nothing needs to be done
        If (dabs(MatA_Gauss(j,j)) < eps_zero) then
            Do k = j+1, rankA
                if (dabs(MatA_Gauss(k,j)) > eps_zero) then
                    do l = j, rankA
                        MatA_Gauss(j,l) = MatA_Gauss(j,l) - MatA_Gauss(k,l)
                    end do
                    vectorb(j) = vectorb(j) - vectorb(k)
                    exit
                end if
                if (k == rankA) then
                    !No row was found that has an entry in the needed column --> quit with error
                    det_A = 0.D0
                    errorflag = -3322
                    return
                end if
            end do
        end if

        Do i = j+1, rankA     !Loop over the rows of matrix A
            !For all rows i>j, make entry i,j 0 by subtracting row j
            if (dabs(MatA_Gauss(i,j)) > eps_zero) then
                factor_lc = MatA_Gauss(i,j)/MatA_Gauss(j,j)
                do l = j, rankA
                    MatA_Gauss(i,l) = MatA_Gauss(i,l) - factor_lc * MatA_Gauss(j,l)
                end do
                vectorb(i) = vectorb(i) - factor_lc * vectorb(j)
            end if
        end do
    end do


    !Calculate the determinant
    Det_A = 1.D0
    Do i = rankA, 1, -1
        !Calculate determinant
        Det_A = Det_A * MatA_Gauss(i,i)
    end do

    !Calculate the solution of the system of linear equations, if the determinant is not equal (close to) 0
    if (abs(Det_A) > eps_zero) then
        Do i = rankA, 1, -1
            !Calculate the solution x(i)
            vectorx(i) = vectorb(i)/MatA_Gauss(i,i)
            do k= i+1, rankA
                vectorx(i) = vectorx(i) - MatA_Gauss(i,k) * vectorx(k) / MatA_Gauss(i,i)
            end do
        end do
    end if

    end subroutine Gauss_algorithm
    !**************************************************************************


    !**************************************************************************
    module subroutine Adjugate(gl,MatrixA, rankA, adj_A, errorflag)
    !Algorithm to calculate the adjugate of a matrix A
    !INPUT:
    !MatrixA: A quadratic matrix containing a maximum of 60 rows and columns
    !rankA: Number of rows and columns of matrix A
    !OUTPUT:
    !adj_A: the adjugate of matrix A
    !errorflag: indicates if an error occured

    !Andreas, February 2016
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision, dimension(60,60):: adj_A
    integer:: errorflag

    double precision, dimension(60,60):: Minor_A        !The minor of A by leaving out a line and a row
    double precision, dimension(60,60):: cof_A_T        !Transposed of the cofactor matrix of A
    double precision:: det_minor                        !The determinant of the minor of A
    integer:: i,j,k,l

    !Variables needed for Gauß algorithm
    double precision, dimension(60) :: vector_b, vector_x

    vector_b = 0.D0
    vector_x = 0.D0

    !The adjugate of A is calculated according to the following formula:
    !adj(A) = cof(A)^T
    !where cof(A)^T is the transpose of the cofactor matrix of A
    !An entry cij of the cofactor matrix of A is obtained by the following equation:
    !cij = (-1)^(i+j) * Mij, where i is the row and j the column of matrix A and Mij is the determinant of the matrix that is created from A by leaving out row i and column j
    !
    if (rankA == 1) then
        adj_A(1,1) = 1.D0
    else
        do i = 1, rankA
            do j = 1, rankA

                !create the minor of A by taking out row i and column j
                do k = 1, rankA-1
                    do l = 1, rankA-1
                        if ((k < i) .and. (l < j)) then
                            Minor_A(k,l) = MatrixA(k,l)
                        end if
                        if ((k < i) .and. (l >= j)) then
                            Minor_A(k,l) = MatrixA(k,l+1)
                        end if
                        if ((k >= i) .and. (l < j)) then
                            Minor_A(k,l) = MatrixA(k+1,l)
                        end if
                        if ((k >= i) .and. (l >= j)) then
                            Minor_A(k,l) = MatrixA(k+1,l+1)
                        end if
                    end do
                end do
                call Gauss_algorithm(gl,Minor_A, rankA-1, vector_b, vector_x, det_minor, errorflag)
                !if (errorflag /= 0) then
                !    return
                !end if
                cof_A_T(i,j) = (-1.D0)**(i+j)* det_minor
                adj_A(j,i) = cof_A_T(i,j)
            end do
        end do
    end if

    end subroutine Adjugate
    !**************************************************************************

    !**************************************************************************
    module subroutine Trace(gl,MatrixA, rankA, trace_A, errorflag)
    !Algorithm to calculate the trace of a matrix A
    !INPUT:
    !MatrixA: A quadratic matrix containing a maximum of 60 rows and columns
    !rankA: Number of rows and columns of matrix A
    !OUTPUT:
    !trace_A: the trace of A
    !errorflag: indicates if an error occured

    !Andreas, February 2016
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision:: trace_A
    integer:: errorflag

    integer:: i

    trace_A = 0.D0

    !The trace is simply the sum of the entries along the diagonal of A
    do i = 1, rankA
        trace_A = trace_A + MatrixA(i,i)
    end do

    end subroutine Trace
    !**************************************************************************



    !**************************************************************************
    module subroutine Mat_mult(gl,MatrixA, rowA, colA, MatrixB, rowB, colB, MatrixC, errorflag)
    !Algorithm to multipy two matrices A * B = C
    !INPUT:
    !MatrixA:   Left matrix for multiplication
    !rowA:      Number of rows of matrix A
    !colA:      Number of columns of matrix A
    !MatrixB:   Right matrix for multiplication
    !rowB:      Number of rows of matrix B
    !colB:      Number of columns of matrix B
    !OUTPUT:
    !Matrix C   The product of the matrix multiplication with dimension rowA X colB
    !errorflag: indicates if an error occured

    !Andreas, February 2016
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60):: MatrixA
    double precision, dimension(60,60):: MatrixB
    double precision, dimension(60,60):: MatrixC
    integer:: rowA, rowB, colA, colB
    integer:: errorflag

    integer:: i, j, k

    MatrixC = 0.D0

    if (colA /= rowB) then
        errorflag = -3322
        return
    end if

    !Multiply the matrices
    do i = 1, rowA
        do j = 1, colB
            do k = 1, rowB
                MatrixC(i,j) = MatrixC(i,j) + MatrixA(i,k) * MatrixB(k,j)
            end do
        end do
    end do

    end subroutine Mat_mult
    !**************************************************************************



    module subroutine crit_pt(gl,z_given, T_crit, p_crit, d_crit, nr_crit_pts, errorflag)
    !**************************************************************************
    !Routine to calculate the critical point(s) of a mixture
    !For an explanation of the algorithm, see
    !Bell and Jäger, XXX
    !
    !INPUT:
    !   z_given: overall composition of the mixture
    !OUTPUT:
    !   T_crit: Vector with critical temperatures
    !   p_crit: vector with critical pressures
    !   d_crit: vector with critical densities
    !   nr_crit_pts: Number of critical points found
    !   errorflag: indicates an error during calculations
    !**************************************************************************






    implicit none

    type(type_gl) :: gl


    double precision, dimension(30):: z_given
    double precision, dimension(10):: T_crit
    double precision, dimension(10):: p_crit
    double precision, dimension(10):: d_crit
    integer:: nr_crit_pts
    integer:: errorflag

    integer:: i, j, k, l, m, n

    Double precision::Temp, dens, press, delta, tau

    double precision:: rhoredmix_orig, tredmix_orig
    double precision, dimension(30):: z_orig

    double precision, dimension(:,:),allocatable:: Matrix_L       !dimension(60,60) needed for Gauß algorithm
    double precision, dimension(:,:),allocatable:: Matrix_M       !dimension(60,60) needed for Gauß algorithm
    integer:: rankLM                                    !rank of matrix L and M

    double precision, dimension(:,:),allocatable:: Adj_L              !The adjugate of L
    double precision, dimension(:,:),allocatable:: n2dL_dni           !Derivative of matrix L with respect to mole numbers ni, multiplied by n²
    double precision, dimension(:,:),allocatable:: MultLdL            !Result of matrix multiplication Adj_L * n2dL_dni

    double precision, dimension(:,:),allocatable:: nd2nar_dnidnj      !second derivative of alphar at constant T and V with respect to ni and nj multiplied with n
    double precision, dimension(:,:),allocatable:: n_dlnfi_dnj        !derivative of fi with respect to nj at constant T and V multiplied with n

    !Double precision, dimension(30,30,30):: n2_dlnfi_dnjdnk  !Second derivative of the logarithm of the fugacity of component i with respect to nj and nk at constant T and V, multiplied with n²
    double precision, dimension(:,:,:),allocatable:: n_dndlnfi_dnjdnk  !derivative of the derivative of the logarithm of the fugacity of component i with respect to nj multiplied by n with respect to nk. The whole term is multiplied by n

    double precision, dimension(60):: n3_detL_dni       !Derivative of the determinant of L with respect to mole number ni, multiplied by n³

    double precision:: det_L, det_M
    double precision, dimension(60):: vector_b, vector_x               !Dummy arrays needed for the gauß algorithm

    !Variables needed for Newton-Raphson
    double precision, dimension(2,2):: Jac_inv
    double precision, dimension(2):: delta_X
    double precision:: Det_Jac
    double precision:: tau_new, del_new
    double precision:: eps_det, residuum
    double precision:: ddetL_dtau, ddetL_ddel, ddetM_dtau, ddetM_ddel           !Derivatives of the matrices L and M with respect to tau and delta
    double precision, dimension(:,:),allocatable:: d_n_dlnfi_dnj_dtau, d_n_dlnfi_dnj_ddel
    double precision, dimension(:,:,:),allocatable:: d_n2_dlnfi_dnjdnk_dtau, d_n2_dlnfi_dnjdnk_ddel
    double precision, dimension(:,:),allocatable :: dL_dtau, dL_ddel, dM_dtau ,dM_ddel     !Derivatives of matrices L and M with respect to tau and delta
    double precision, dimension(:,:),allocatable :: Adj_M                                  !The adjugate of M
    double precision, dimension(:,:),allocatable :: MultLdLdtau                            !Result of matrix multiplication Adj_L * dL_dtau
    double precision, dimension(:,:),allocatable :: MultLdLddel                            !Result of matrix multiplication Adj_L * dL_ddel
    double precision, dimension(:,:),allocatable :: MultMdMdtau                            !Result of matrix multiplication Adj_M * dM_dtau
    double precision, dimension(:,:),allocatable :: MultMdMddel                            !Result of matrix multiplication Adj_M * dM_ddel
    double precision, dimension(:,:),allocatable :: Mult_L2dni_adjLdtau, Mult_L2dni_adjLddel       !Result of matrix multiplication n2dLdni * d(adj(L))_dtau AND n2dLdni * d(adj(L))_ddel
    double precision, dimension(:,:),allocatable :: Mult_adjL_dn2Ldnidtau, Mult_adjL_dn2Ldniddel   !Result of matrix multiplication adj(L) * d(n^2 dL / dni)/dtau AND adj(L) * d(n^2 dL / dni)/ddel
    Double precision, dimension(:,:),allocatable :: helpM_tau, helpM_del
    double precision, dimension(:,:),allocatable :: d_adjL_dtau, d_adjL_ddel
    double precision, dimension(:,:),allocatable :: d_n2dL_dni_dtau, d_n2dL_dni_ddel
    double precision:: dMNi_dtau, dMNi_ddel

    !double precision:: P_CALC

    !Help variables for cubic EOS
    double precision:: rhoc_est, Tc_est

    !Variables for bisection of first point on detL=0 contour
    double precision:: Temp_high, Temp_low

    !Variables for contour tracer
    double precision:: theta_detL, Rtau, Rdel, dL_dtheta, delta_theta
    double precision, dimension(1000):: T_pts_contourL, D_pts_contourL, detL_pts_contourL, detM_pts_contourL, theta_detL_contourL

    if(.not.(allocated(dL_dtau))) then
        allocate(dL_dtau(60,60))
        allocate( dL_ddel, dM_dtau ,dM_ddel,Adj_M,MultLdLdtau,MultLdLddel, MultMdMdtau,MultMdMddel,Mult_L2dni_adjLdtau, Mult_L2dni_adjLddel,     &
            Mult_adjL_dn2Ldnidtau, Mult_adjL_dn2Ldniddel, &
            helpM_tau, helpM_del, &
            d_adjL_dtau, d_adjL_ddel, &
            d_n2dL_dni_dtau, d_n2dL_dni_ddel, &
            Matrix_L, Matrix_M,Adj_L    ,n2dL_dni ,MultLdL  ,mold=dL_dtau)
        allocate(d_n2_dlnfi_dnjdnk_dtau(30,30,30), d_n2_dlnfi_dnjdnk_ddel(30,30,30),n_dndlnfi_dnjdnk(30,30,30))
        allocate(d_n_dlnfi_dnj_dtau(30,30))
        allocate(d_n_dlnfi_dnj_ddel,nd2nar_dnidnj , n_dlnfi_dnj   ,mold=d_n_dlnfi_dnj_dtau)
        allocate(gl%n2d2lnfidnjdnk(30,30,30),&
            gl%ndndlnfidnjdnk(30,30,30),&
            gl%dndndlnfidnjdnkddel(30,30,30),&
            gl%dn2d2lnfidnjdnkddel(30,30,30),&
            gl%dn2d2lnfidnjdnkdtau(30,30,30),&
            gl%dndndlnfidnjdnkdtau(30,30,30),&
            gl%d_n_dlnfi_dnj_dxm(30,30,30),&
            gl%d2_n_dlnfi_dnj_dxmdtau(30,30,30),&
            gl%d2_n_dlnfi_dnj_dxmddel(30,30,30),&
            gl%dndnardnidnjALL(30,30,30),&
            gl%dndlnfidnjdxk(30,30,30),&
            gl%ndardnidxi_all(30,30,30),&
            gl%dndndardnidnjdxk_ALL(30,30,30,30),&
            gl%d2ndlnfidnjdxkddel(30,30,30),&
            gl%d2ndlnfidnjdxkdtau(30,30,30),&
            gl%d2PSI_Y_dxjdxk(30,30,30,30),&
            gl%MIXDERIVFNR_dxidxj(30,30,30),&
            gl%dPSI_Y_dxj(30,30,30),&
            gl%MIXDERIVFNR_dxidxjdxk(30,30,30,30),&
            gl%PSI_Y(30,30),&
            gl%ndndardnidnj_ALL(30,30,30),&
            gl%dndardnidxi_ALL(30,30,30),&
            gl%d2ndardnidxidxj_ALL(30,30,30,30))
    end if

    n3_detL_dni = 0.D0
    vector_b = 0.D0
    vector_x = 0.D0

    !Save original overall composition and reducing temperature and density
    z_orig = gl%molfractions
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    gl%molfractions = z_given
    !Temp = 350.322D0       !Dummy value
    !tau = Tredmix / temp
    !Dens = 5327.48D0       !Dummy value
    !delta = dens / rhoredmix
    call reduced_parameters_calc(gl,Temp)


    rankLM = gl%ncomp          !Rank of matrices L and M
    eps_det = 1.D-8         !Exit criterion for Gauß algorithm
    nr_crit_pts = 0



    !----------------------------------------------------------------------------
    !Find critical points by tracing along the det(L) = 0 contour
    !----------------------------------------------------------------------------
    !Step 1:    Find the first point on the contour at high temperature
    !           Set T = 1.5 Tredmix
    !           Set rho = 0.5 rhoredmix
    !
    If ((gl%mix_type == 1) .or. (gl%mix_type == 11) .or. (gl%mix_type == 12)  .or. (gl%mix_type == 13) .or. (gl%mix_type == 19)) then    !Helmholtz models
        rhoc_est = gl%rhoredmix
        Tc_est = gl%tredmix
    elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .or. (gl%mix_type == 3) .or. (gl%mix_type == 31)) then  !Cubic EOS
        rhoc_est = 0.D0
        Tc_est = 0.D0
        do i = 1,gl%ncomp
            rhoc_est = rhoc_est + gl%molfractions(i)/gl%rhoc(i)    !Estimate a critical density by linear combination of the pure fluid densities
            Tc_est = Tc_est + gl%molfractions(i) * gl%tc(i)        !Estimate a critical temperature by linear cobination of the pure fluid temperatures
        end do
        rhoc_est = 1.D0/rhoc_est
    else
        !wrong mixtype or model not yet implemented
        errorflag = -4408
        return
    end if

    Temp = Tc_est * 1.5D0
    Dens = rhoc_est * 0.5D0
    call reduced_parameters_calc(gl,Temp)


    !++
    !The assumption is that det_L > 0 for very high temperatures
    !Decrease the Temperature slowly until detL < 0
    do k = 1, 1000
        !Calculate the determinant of L at the first point
        call ndlnfi_dnj(gl,Temp, Dens, n_dlnfi_dnj)
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                Matrix_L(i,j) = n_dlnfi_dnj(i,j)
            end do
        end do
        !Calculate the determinant of matrix L
        call Gauss_algorithm(gl,Matrix_L, rankLM, vector_b, vector_x, det_L, errorflag)
        if (errorflag /= 0) then
            return
        end if
        if (det_L > 0.D0) then
            Temp_high = Temp
            Temp = Temp - 5.0D0
            call reduced_parameters_calc(gl,Temp)
        else
            Temp_low = Temp
            exit
        end if
    end do
    !++


    !++
    Do k = 1, 100
        !Do slow but robust bisection for finding the correct interval  --> THIS COULD BE IMPROVED ACCORDING TO THE PAPER!!!
        Temp = 0.5D0 * (Temp_high + Temp_low)
        call reduced_parameters_calc(gl,Temp)
        !Calculate the determinant of L at the first point
        call ndlnfi_dnj(gl,Temp, Dens, n_dlnfi_dnj)
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                Matrix_L(i,j) = n_dlnfi_dnj(i,j)
            end do
        end do
        !Calculate the determinant of matrix L
        call Gauss_algorithm(gl,Matrix_L, rankLM, vector_b, vector_x, det_L, errorflag)
        if (errorflag /= 0) then
            return
        end if

        If (dabs(det_L) < 1.D-6) then
            exit
        end if

        if (det_L > 0.D0) then
            Temp_high = Temp
        else
            Temp_low = Temp
        end if

    End do
    !++

    !Step 2: Trace along the contour and stop when there is a change of signs in detM
    if ((gl%mix_type == 1)  .or. (gl%mix_type == 11) .or. (gl%mix_type == 12)  .or. (gl%mix_type == 13) .or. (gl%mix_type == 19)) then   !Search more carefully in case of multifluid models
        Rtau = 0.1D0
        Rdel = 0.025D0
    elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .or. (gl%mix_type == 3) .or. (gl%mix_type == 31)) then   !Cubic EOS, scale the radii
        Rtau = 0.1D0 / Tc_est
        Rdel = 0.025D0 * rhoc_est
    end if

    do n = 1, 1000

        !Calculate detM
        call ndlnfi_dnj(gl,Temp, Dens, n_dlnfi_dnj)
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                Matrix_M(i,j) = n_dlnfi_dnj(i,j)
            end do
        end do
        call Adjugate(gl,Matrix_L, rankLM, Adj_L, errorflag)
        if (errorflag /= 0) then
            return
        end if
        call n2d2lnfi_dnjdnk(gl,Temp, Dens)
        do k = 1, gl%ncomp
            !Calculate the derivative of n_dlnfi_dnj with respect to nk
            Do i = 1, gl%ncomp
                do j = 1, gl%ncomp
                    n2dL_dni(i,j) = gl%n2d2lnfidnjdnk(i,j,k)
                end do
            end do
            !Matrix multiplication of the adjugate of L with dL_dni
            call Mat_mult(gl,Adj_L, rankLM, rankLM, n2dL_dni, rankLM, rankLM, MultLdL, errorflag)
            if (errorflag /= 0) then
                return
            end if
            call Trace(gl,MultLdL, rankLM, n3_detL_dni(k), errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Add last row to Matrix_M
            Matrix_M(gl%ncomp,k) = n3_detL_dni(k)
        end do
        !Calculate the determinant of matrix M
        call Gauss_algorithm(gl,Matrix_M, rankLM, vector_b, vector_x, det_M, errorflag)
        if (errorflag /= 0) then
            return
        end if

        T_pts_contourL(n) = Temp
        D_pts_contourL(n) = Dens
        detL_pts_contourL(n) = det_L
        detM_pts_contourL(n) = det_M
        tau = gl%tredmix/temp
        delta = dens/gl%rhoredmix

        !Break criterion for iteration
        !If the density exceeds 5.0 * rhoc_est, stop the search
        if (dens > 10.0D0 * rhoc_est) then
            return
        end if


        if (n>1) then
            if ((detM_pts_contourL(n-1) * det_M) < 0.D0) then       !Do NEWTON-RAPHSON TO CALCULATE THE CRITICAL POINT

                !Loop for Newton-Raphson to find T and rho at critical point
                do m = 1, 40

                    if (m == 40) then
                        errorflag = -4408
                        return
                    end if

                    !-----------------------------------------------------------------------------------------
                    !First step: Calculate the matrices L* and M* and the determinants of the matrices
                    !-----------------------------------------------------------------------------------------
                    !Calculate the entries of the Matrix L* and M*
                    !The last row of Matrix M* will be calculated later
                    call ndlnfi_dnj(gl,Temp, Dens, n_dlnfi_dnj)
                    do i = 1, gl%ncomp
                        do j = 1, gl%ncomp
                            Matrix_L(i,j) = n_dlnfi_dnj(i,j)
                            Matrix_M(i,j) = n_dlnfi_dnj(i,j)
                        end do
                    end do
                    !Last row of matrix M is the derivative of the determinant of matrix L* with respect to the mole number n1...nN
                    !Therefore, first calculate the derivatives of the determinant of matrix L* with respect to mole numbers n1...nN
                    !According to Jacobis Formula, the derivative of a determinant of the matrix L* can be obtained according to
                    !----
                    ! ddet(L)_dni = tr(adj(L) * dL_dni) --> See Eq. 17 and Eq. 18 in Bell and Jäger (2016)
                    !----

                    !Get the adjugate of L
                    call Adjugate(gl,Matrix_L, rankLM, Adj_L, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    !Calculate the derivative needed for the last row of matrix M
                    !Get n^2 d2 (lnfi) / dnidnk at constant T,V,nm
                    !call n2d2lnfi_dnjdnk_TV_num (Temp, Dens, gl%n2d2lnfidnjdnk)  !- NUMERICALLY
                    call n2d2lnfi_dnjdnk(gl,Temp, Dens) !- ANALYTICALLY
                    do k = 1, gl%ncomp
                        !Calculate the derivative of n_dlnfi_dnj with respect to nk
                        Do i = 1, gl%ncomp
                            do j = 1, gl%ncomp
                                n2dL_dni(i,j) = gl%n2d2lnfidnjdnk(i,j,k)
                            end do
                        end do
                        !Matrix multiplication of the adjugate of L with dL_dni
                        call Mat_mult(gl,Adj_L, rankLM, rankLM, n2dL_dni, rankLM, rankLM, MultLdL, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if
                        call Trace(gl,MultLdL, rankLM, n3_detL_dni(k), errorflag)
                        if (errorflag /= 0) then
                            return
                        end if
                        !Add last row to Matrix_M
                        Matrix_M(gl%ncomp,k) = n3_detL_dni(k)
                    end do
                    !End of calculation of last row of M


                    !Calculate the determinant of matrix L
                    call Gauss_algorithm(gl,Matrix_L, rankLM, vector_b, vector_x, det_L, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Calculate the determinant of matrix M
                    call Gauss_algorithm(gl,Matrix_M, rankLM, vector_b, vector_x, det_M, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    residuum = (det_L**2 + det_M**2)**0.5
                    if (residuum < eps_det) then
                        nr_crit_pts = nr_crit_pts + 1
                        T_crit(nr_crit_pts) = Temp
                        d_crit(nr_crit_pts) = dens
                        p_crit(nr_crit_pts) = P_CALC(gl,Temp, Dens, 0)
                        Temp = T_pts_contourL(n)
                        Dens = D_pts_contourL(n)
                        det_L = detL_pts_contourL(n)
                        det_M = detM_pts_contourL(n)
                        tau = gl%tredmix / temp
                        delta = dens / gl%rhoredmix
                        call reduced_parameters_calc(gl,Temp)
                        exit
                    end if
                    !-----------------------------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------
                    !Second step:   Calculate the tau and delta derivatives of the determinants of matrices
                    !               L and M for the Jacobian Matrix of Newton_Raphson.
                    !               Then, calculate the step in delta and tau
                    !-----------------------------------------------------------------------------------------
                    !Get tau-derivative of n*d(lnfi)/dnj
                    call dndlnfi_dnjdtau(gl,Temp, Dens, d_n_dlnfi_dnj_dtau)
                    d_n_dlnfi_dnj_dtau = d_n_dlnfi_dnj_dtau / tau
                    !Get del-derivative of n*d(lnfi)/dnj
                    call dndlnfi_dnjddel(gl,Temp, Dens, d_n_dlnfi_dnj_ddel)
                    d_n_dlnfi_dnj_ddel = d_n_dlnfi_dnj_ddel / delta
                    do i = 1, gl%ncomp
                        do j = 1, gl%ncomp
                            dL_dtau(i,j) = d_n_dlnfi_dnj_dtau(i,j)
                            dL_ddel(i,j) = d_n_dlnfi_dnj_ddel(i,j)
                            dM_dtau(i,j) = d_n_dlnfi_dnj_dtau(i,j)    !LAST ROW OF MATRIX M IS STILL INCORRECT!!
                            dM_ddel(i,j) = d_n_dlnfi_dnj_ddel(i,j)    !LAST ROW OF MATRIX M IS STILL INCORRECT!!
                        end do
                    end do

                    !Matrix multiplication of the adjugate of L with dL_dtau
                    call Mat_mult(gl,Adj_L, rankLM, rankLM, dL_dtau, rankLM, rankLM, MultLdLdtau, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Matrix multiplication of the adjugate of L with dL_ddel
                    call Mat_mult(gl,Adj_L, rankLM, rankLM, dL_ddel, rankLM, rankLM, MultLdLddel, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Get the derivative of the determinant of L with respect to tau
                    call Trace(gl,MultLdLdtau, rankLM, ddetL_dtau, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Get the derivative of the determinant of L with respect to del
                    call Trace(gl,MultLdLddel, rankLM, ddetL_ddel, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    !Calculate the derivatives needed for the last row of matrices dM_dtau and dM_ddel
                    !Get d (n^2 d2 (lnfi) / dnidnk) / dtau at constant del and x
                    call dn2d2lnfi_dnjdnkdtau(gl,Temp, Dens)
                    d_n2_dlnfi_dnjdnk_dtau = d_n2_dlnfi_dnjdnk_dtau / tau
                    !Get d (n^2 d2 (lnfi) / dnidnk) / ddel at constant tau and x
                    call dn2d2lnfi_dnjdnkddel(gl,Temp, Dens)
                    d_n2_dlnfi_dnjdnk_ddel = d_n2_dlnfi_dnjdnk_ddel / delta

                    !Calculate the derivative of the adjugate of L with respect to tau and del
                    call dAdjL_dX(gl,Matrix_L, dL_dtau, dL_ddel, rankLM, d_adjL_dtau, d_adjL_ddel, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    do k = 1, gl%ncomp

                        !Calculate the derivative of n_dlnfi_dnj with respect to nk
                        Do i = 1, gl%ncomp
                            do j = 1, gl%ncomp
                                n2dL_dni(i,j) = gl%n2d2lnfidnjdnk(i,j,k)
                            end do
                        end do

                        !Matrix multiplication of the adjugate of n2dL/dni with d(adj(L))/dtau
                        call Mat_mult(gl,n2dL_dni, rankLM, rankLM, d_adjL_dtau, rankLM, rankLM, Mult_L2dni_adjLdtau, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if
                        !Matrix multiplication of the adjugate of n2dL/dni with d(adj(L))/ddel
                        call Mat_mult(gl,n2dL_dni, rankLM, rankLM, d_adjL_ddel, rankLM, rankLM, Mult_L2dni_adjLddel, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if

                        !Calculate the derivative of n2dL/dni with respect to tau and del
                        Do i = 1, gl%ncomp
                            do j = 1, gl%ncomp
                                d_n2dL_dni_dtau(i,j) = d_n2_dlnfi_dnjdnk_dtau(i,j,k)
                                d_n2dL_dni_ddel(i,j) = d_n2_dlnfi_dnjdnk_ddel(i,j,k)
                            end do
                        end do

                        !Matrix multiplication of the adjugate of L with d(n2dL/dni)/dtau
                        call Mat_mult(gl,Adj_L, rankLM, rankLM, d_n2dL_dni_dtau, rankLM, rankLM, Mult_adjL_dn2Ldnidtau, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if
                        do i = 1, gl%ncomp
                            do j = 1, gl%ncomp
                                helpM_tau(i,j) = Mult_L2dni_adjLdtau(i,j) + Mult_adjL_dn2Ldnidtau(i,j)
                            end do
                        end do
                        call Trace(gl,helpM_tau, rankLM, dMNi_dtau, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if

                        !Matrix multiplication of the adjugate of L with d(n2dL/dni)/ddel
                        call Mat_mult(gl,Adj_L, rankLM, rankLM, d_n2dL_dni_ddel, rankLM, rankLM, Mult_adjL_dn2Ldniddel, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if
                        do i = 1, gl%ncomp
                            do j = 1, gl%ncomp
                                helpM_del(i,j) = Mult_L2dni_adjLddel(i,j) + Mult_adjL_dn2Ldniddel(i,j)
                            end do
                        end do
                        call Trace(gl,helpM_del, rankLM, dMNi_ddel, errorflag)
                        if (errorflag /= 0) then
                            return
                        end if

                        !Add last row to derivative of Matrix_M with respect to tau and delta
                        dM_dtau(gl%ncomp,k) = dMNi_dtau
                        dM_ddel(gl%ncomp,k) = dMNi_ddel
                    end do

                    !Compute the adjugate of Matrix M
                    call Adjugate(gl,Matrix_M, rankLM, Adj_M, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    !Matrix multiplication of the adjugate of M with dM_dtau
                    call Mat_mult(gl,Adj_M, rankLM, rankLM, dM_dtau, rankLM, rankLM, MultMdMdtau, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Matrix multiplication of the adjugate of M with dM_ddel
                    call Mat_mult(gl,Adj_M, rankLM, rankLM, dM_ddel, rankLM, rankLM, MultMdMddel, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    !Get the derivative of the determinant of M with respect to tau
                    call Trace(gl,MultMdMdtau, rankLM, ddetM_dtau, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if
                    !Get the derivative of the determinant of M with respect to del
                    call Trace(gl,MultMdMddel, rankLM, ddetM_ddel, errorflag)
                    if (errorflag /= 0) then
                        return
                    end if

                    !Calculate inverse Jacobian
                    Det_Jac = ddetL_dtau * ddetM_ddel - ddetL_ddel * ddetM_dtau
                    Jac_inv(1,1) = ddetM_ddel / Det_Jac
                    Jac_inv(1,2) = -ddetL_ddel / Det_Jac
                    Jac_inv(2,1) = -ddetM_dtau / Det_Jac
                    Jac_inv(2,2) = ddetL_dtau / Det_Jac

                    !Calculate the step in delta and tau
                    Delta_X(1) = -(Jac_inv(1,1) * det_L + Jac_inv(1,2) * det_M)
                    Delta_X(2) = -(Jac_inv(2,1) * det_L + Jac_inv(2,2) * det_M)
                    !-----------------------------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------
                    !Third step:   Update temperature and density and do next step
                    !-----------------------------------------------------------------------------------------
                    tau_new = tau + Delta_X(1)
                    del_new = delta + Delta_X(2)
                    tau = tau_new
                    delta = del_new
                    Temp = gl%tredmix / tau
                    dens = gl%rhoredmix * delta
                    call reduced_parameters_calc(gl,Temp)
                    !-----------------------------------------------------------------------------------------

                end do

            end if
        end if

        !Search the next point on the detL=0 contour
        if (n == 1) then
            theta_detL = 3.1415D0/2.D0  !pi/2 as startvalue for first point
        else
            theta_detL = theta_detL_contourL(n-1) !Take last theta as start value
        end if
        Do m = 1, 30
            tau_new = tau + Rtau * dcos(theta_detL)
            del_new = delta + Rdel * dsin(theta_detL)
            Temp = gl%tredmix / tau_new
            dens = gl%rhoredmix * del_new
            call reduced_parameters_calc(gl,Temp)
            !Calculate the new value of det_L
            call ndlnfi_dnj(gl,Temp, Dens, n_dlnfi_dnj)
            do i = 1, gl%ncomp
                do j = 1, gl%ncomp
                    Matrix_L(i,j) = n_dlnfi_dnj(i,j)
                end do
            end do
            !Calculate the determinant of matrix L
            call Gauss_algorithm(gl,Matrix_L, rankLM, vector_b, vector_x, det_L, errorflag)
            if (errorflag /= 0) then
                return
            end if

            If (dabs(det_L) < 1.D-6) then
                if (n>1) then
                    if (dabs(theta_detL - theta_detL_contourL(n-1)) > 1.57D0) then        !Sharp bench, stop calculations
                        return
                    end if
                end if
                theta_detL_contourL(n) = theta_detL
                exit
            end if

            !Calculate derivative of detL with respect to theta
            !Get tau-derivative of n*d(lnfi)/dnj
            call dndlnfi_dnjdtau(gl,Temp, Dens, d_n_dlnfi_dnj_dtau)
            d_n_dlnfi_dnj_dtau = d_n_dlnfi_dnj_dtau / tau_new
            !Get del-derivative of n*d(lnfi)/dnj
            call dndlnfi_dnjddel(gl,Temp, Dens, d_n_dlnfi_dnj_ddel)
            d_n_dlnfi_dnj_ddel = d_n_dlnfi_dnj_ddel / del_new
            do i = 1, gl%ncomp
                do j = 1, gl%ncomp
                    dL_dtau(i,j) = d_n_dlnfi_dnj_dtau(i,j)
                    dL_ddel(i,j) = d_n_dlnfi_dnj_ddel(i,j)
                end do
            end do
            !Matrix multiplication of the adjugate of L with dL_dtau
            call Mat_mult(gl,Adj_L, rankLM, rankLM, dL_dtau, rankLM, rankLM, MultLdLdtau, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Matrix multiplication of the adjugate of L with dL_ddel
            call Mat_mult(gl,Adj_L, rankLM, rankLM, dL_ddel, rankLM, rankLM, MultLdLddel, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Get the derivative of the determinant of L with respect to tau
            call Trace(gl,MultLdLdtau, rankLM, ddetL_dtau, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Get the derivative of the determinant of L with respect to del
            call Trace(gl,MultLdLddel, rankLM, ddetL_ddel, errorflag)
            if (errorflag /= 0) then
                return
            end if
            dL_dtheta = -ddetL_dtau * Rtau * dsin(theta_detL) + ddetL_ddel * Rdel * dcos(theta_detL)

            !Update theta
            delta_theta = -det_L/dL_dtheta
            theta_detL = theta_detL + delta_theta

        end do

    end do
    !----------------------------------------------------------------------------

    !Set back overall composition and reducing temperature and density
    gl%molfractions = z_orig
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig

    end subroutine crit_pt
    !**************************************************************************


    !**************************************************************************
    module subroutine dAdjL_dX(gl,MatrixL, dL_dtau, dL_ddel, rankL, dadjL_dtau, dadjL_ddel, errorflag)
    !Algorithm to calculate the derivative of the adjugate of matrix L with
    !respect to tau and with respect to delta
    !
    !INPUT:
    !MatrixL: A quadratic matrix containing a maximum of 60 rows and columns
    !dL_dtau: The derivatives of matrix L with respect to tau
    !dL_ddel: The derivatives of matrix L with respect to del
    !rankL: Number of rows and columns of matrix L
    !OUTPUT:
    !dadjL_dtau: the derivative of the adjugate of matrix L wrt tau
    !dadjL_ddel: the derivative of the adjugate of matrix L wrt del
    !errorflag: indicates if an error occured

    !Andreas, March 2016
    !**************************************************************************


    implicit none

    type(type_gl) :: gl


    double precision, dimension(60,60):: MatrixL
    double precision, dimension(60,60):: dL_dtau
    double precision, dimension(60,60):: dL_ddel
    integer:: rankL
    double precision, dimension(60,60):: dadjL_dtau
    double precision, dimension(60,60):: dadjL_ddel
    integer:: errorflag

    double precision, dimension(60,60):: Minor_L        !The minor of L by leaving out a line and a row
    double precision, dimension(60,60):: dMinorL_dtau   !The derivative of the minor of L with respect to tau
    double precision, dimension(60,60):: dMinorL_ddel   !The derivative of the minor of L with respect to del
    double precision, dimension(60,60):: adj_M          !The adjugate of the minor
    double precision, dimension(60,60):: Mult_adjM_dMdtau, Mult_adjM_dMddel
    double precision:: tr_adjM_dMdtau, tr_adjM_dMddel

    integer:: i,j,k,l

    !Explantion of the routine given in article of Bell and Jäger (2016), page 6, Eq. 30-32 XXX
    if (rankL == 1) then
        dadjL_dtau(1,1) = 1.D0
        dadjL_ddel(1,1) = 1.D0
    else
        do i = 1, rankL
            do j = 1, rankL

                !create the minor of L by taking out row i and column j
                do k = 1, rankL-1
                    do l = 1, rankL-1
                        if ((k < i) .and. (l < j)) then
                            Minor_L(k,l) = MatrixL(k,l)
                            dMinorL_dtau(k,l) = dL_dtau(k,l)
                            dMinorL_ddel(k,l) = dL_ddel(k,l)
                        end if
                        if ((k < i) .and. (l >= j)) then
                            Minor_L(k,l) = MatrixL(k,l+1)
                            dMinorL_dtau(k,l) = dL_dtau(k,l+1)
                            dMinorL_ddel(k,l) = dL_ddel(k,l+1)
                        end if
                        if ((k >= i) .and. (l < j)) then
                            Minor_L(k,l) = MatrixL(k+1,l)
                            dMinorL_dtau(k,l) = dL_dtau(k+1,l)
                            dMinorL_ddel(k,l) = dL_ddel(k+1,l)
                        end if
                        if ((k >= i) .and. (l >= j)) then
                            Minor_L(k,l) = MatrixL(k+1,l+1)
                            dMinorL_dtau(k,l) = dL_dtau(k+1,l+1)
                            dMinorL_ddel(k,l) = dL_ddel(k+1,l+1)
                        end if
                    end do
                end do

                !Calculate the adjugate of the minor of L
                call Adjugate(gl,Minor_L, rankL-1, adj_M, errorflag)
                if (errorflag /= 0) then
                    return
                end if

                !Compute the tau derivative of the adjugate of L
                call Mat_mult(gl,adj_M, rankL-1, rankL-1, dMinorL_dtau, rankL-1, rankL-1, Mult_adjM_dMdtau, errorflag)
                if (errorflag /= 0) then
                    return
                end if
                call trace(gl,Mult_adjM_dMdtau, rankL-1, tr_adjM_dMdtau, errorflag)
                if (errorflag /= 0) then
                    return
                end if
                dadjL_dtau(j,i) = (-1.D0)**(i+j)*tr_adjM_dMdtau

                !Compute the del derivative of the adjugate of L
                call Mat_mult(gl,adj_M, rankL-1, rankL-1, dMinorL_ddel, rankL-1, rankL-1, Mult_adjM_dMddel, errorflag)
                if (errorflag /= 0) then
                    return
                end if
                call trace(gl,Mult_adjM_dMddel, rankL-1, tr_adjM_dMddel, errorflag)
                if (errorflag /= 0) then
                    return
                end if
                dadjL_ddel(j,i) = (-1.D0)**(i+j)*tr_adjM_dMddel

            end do
        end do
    end if

    end subroutine dAdjL_dX
    !**************************************************************************


    !**************************************************************************
    !OLD SUCCESSIVE SUBSTITUTION (iflash = 5) FOR UP TO THREE PHASES


    !   !0) Generate initial estimates for the phase fractions and K-values
    !    !   Assumption: 3 phases exist
    !
    !    K_ij = 0.D0
    !    x_j = 0.D0
    !    F_beta = 0.D0
    !    J_beta = 0.D0
    !    Delta_beta = 0.D0
    !    x_old = 0.D0
    !
    !    NrofPhases = 2
    !    beta_j(1) = 0.5D0
    !    beta_j(2) = 0.5D0
    !    x_j(1,:) = x_vap
    !    x_j(2,:) = x_liq
    !    !Calculate K values
    !    do i = 1, ncomp
    !        K_ij(1,i) = x_j(1,i) / x_j(2,i)
    !    end do
    !
    !!----------------------------------------------------------
    !!    !Initial guess: three phases are in equilibrium
    !!    NrofPhases = 3
    !!
    !!    !Generate start values for x1 and x3 (Assuming a 2 phase equilibrium)
    !!    call PTX_startvals (press, Temp, x_known, x_j(3,:), x_j(1,:), beta_j(1), iFlash, errval)
    !!
    !!    beta_j(1) = 0.4D0
    !!    beta_j(2) = 0.3D0
    !!    beta_j(3) = 1.D0 - beta_j(1) - beta_j(2)
    !!
    !!    !Calculate K values
    !!    do i = 1, ncomp
    !!        K_ij(1,i) = x_j(1,i) / x_j(3,i)
    !!    end do
    !!
    !!    !Generate third phase composition and K-values
    !!    do i = 1, ncomp
    !!        x_j(2,i) = (x_known(i) - beta_j(1) *x_j(1,i) - beta_j(3) *x_j(3,i)) / beta_j(2)
    !!        K_ij(2,i) = x_j(2,i) / x_j(3,i)
    !!    end do
    !!------------------------------------------------------------
    !
    !    !Start iteration process
    !    Do k = 1, Nr_of_iter
    !
    !        !1) Solve Rachford Rice system of equations for the phase fractions vapfrac
    !        !   using Newton Raphson method
    !        Do o = 1, 30
    !
    !            F_beta = 0.D0
    !            J_beta = 0.D0
    !            !Set up the system of equations
    !            Do j = 1, NrofPhases - 1
    !                Do i = 1, ncomp
    !                    denominator = 1.d0
    !                    Do n = 1, NrofPhases - 1
    !                        denominator = denominator + beta_j(n) * (K_ij(n,i) - 1.D0)
    !                    End do
    !                    F_beta(j) = F_beta(j) + x_known(i) * (K_ij(j,i) -1.D0)  / denominator
    !                end do
    !            end do
    !
    !!            F_beta_test(1) = x_known(1)*(K_ij(1,1) - 1.D0) &
    !!                            & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0)) &
    !!                            & + x_known(2)*(K_ij(1,2) - 1.D0) &
    !!                            & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))
    !!            F_beta_test(2) = x_known(1)*(K_ij(2,1) - 1.D0) &
    !!                            & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0)) &
    !!                            & + x_known(2)*(K_ij(2,2) - 1.D0) &
    !!                            & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))
    !            !Exit criterion for vapfrac
    !            if (maxval(dabs(F_beta)) < 1.D-8) exit
    !
    !            !Calculate the Jacobi matrix
    !            Do m = 1, NrofPhases - 1
    !                Do j = 1, NrofPhases - 1
    !                    Do i = 1, ncomp
    !                        denominator = 1.d0
    !                        Do n = 1, NrofPhases - 1
    !                            denominator = denominator + beta_j(n) * (K_ij(n,i) - 1.D0)
    !                        End do
    !                        J_beta(j,m) = J_beta(j,m) - x_known(i)* (K_ij(j,i) - 1.D0) * (K_ij(m,i) - 1.D0) / &
    !                                        & (denominator * denominator)
    !                    end do
    !                end do
    !            end do
    !!            J_beta_test(1,1) = - x_known(1)*(K_ij(1,1) - 1.D0) *(K_ij(1,1)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0))**2 &
    !!                             & - x_known(2)*(K_ij(1,2) - 1.D0) *(K_ij(1,2)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))**2
    !!            J_beta_test(1,2) = - x_known(1)*(K_ij(1,1) - 1.D0) *(K_ij(2,1)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0))**2 &
    !!                             & - x_known(2)*(K_ij(1,2) - 1.D0) *(K_ij(2,2)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))**2
    !!            J_beta_test(2,1) = - x_known(1)*(K_ij(2,1) - 1.D0) *(K_ij(1,1)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0))**2 &
    !!                             & - x_known(2)*(K_ij(2,2) - 1.D0) *(K_ij(1,2)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))**2
    !!            J_beta_test(2,2) = - x_known(1)*(K_ij(2,1) - 1.D0) *(K_ij(2,1)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,1) - 1.D0) + beta_j(2) * (K_ij(2,1) - 1.D0))**2 &
    !!                             & - x_known(2)*(K_ij(2,2) - 1.D0) *(K_ij(2,2)-1.D0) &
    !!                             & / (1.D0 + beta_j(1) * (K_ij(1,2) - 1.D0) + beta_j(2) * (K_ij(2,2) - 1.D0))**2
    !            !Solve the system of equations
    !            Delta_beta = -F_beta
    !            call LUdecomp (NrofPhases-1,J_beta,Delta_beta,ierr,herr)
    !
    !            !Change the variables
    !            beta_j(NrofPhases) = 1.D0
    !            do j = 1, NrofPhases - 1
    !                beta_j(j) = beta_j(j) + Delta_beta(j)
    !                beta_j(NrofPhases) = beta_j(NrofPhases) - beta_j(j)
    !            end do
    !
    !        end do
    !
    !        ! If the phase fraction of a trial phase exceeds the limits of 0 < vapfrac < 1 then the phase is assumed to be unstable
    !        ! Delete that phase from the calculations and carry on with the other phases
    !        ! Check the reference phase (last phase) first
    !        if (beta_j(NrofPhases) < 0.D0) then
    !            if (NrofPhases == 2) then
    !               x_vap = x_j(1,:)
    !               x_liq = x_j(2,:)
    !               vapfrac = beta_j(1)
    !               return    !If only one phase is left this calculation does not make sense any more
    !            End if
    !            beta_j(NrofPhases) = 0.D0
    !            x_j(NrofPhases,:) = 0.D0
    !            NrofPhases = NrofPhases - 1
    !            beta_j = 1.D0 / Nrofphases
    !        else
    !            do j = 1, NrofPhases-1
    !                if (beta_j(j) < 0.D0) then
    !                     if (NrofPhases == 2) then
    !                         x_vap = x_j(1,:)
    !                         x_liq = x_j(2,:)
    !                         vapfrac = beta_j(1)
    !                         return    !If only one phase is left this calculation does not make sense any more
    !                    End if
    !                    !Rearrange the array and matrices used: The jth phase is unstable so it will no longer be considered and therefore erased
    !                    Do m = j, NrofPhases - 1
    !                        beta_j(m) = beta_j(m + 1)
    !                        x_j(m,:) = x_j(m+1,:)
    !                    End do
    !                    NrofPhases = NrofPhases - 1
    !                    beta_j = 1.D0 / Nrofphases
    !               end if
    !            end do
    !        end if
    !
    !
    !        !2) Calculate new compositions
    !        do i = 1, ncomp
    !            denominator = 1.d0
    !            Do n = 1, NrofPhases - 1
    !                denominator = denominator + beta_j(n) * (K_ij(n,i) - 1.D0)
    !            End do
    !            x_j(NrofPhases,i) = x_known(i) / denominator
    !            do m = 1, NrofPhases - 1
    !                x_j(m,i) = (x_known(i) * K_ij(m,i)) / denominator
    !            end do
    !        end do
    !
    !        sum = 0.D0
    !        !Normieren
    !        Do j = 1, NrofPhases
    !            Do i = 1, ncomp
    !                sum(j) = sum(j) + x_j(j,i)
    !            end do
    !            x_j(j,:) = x_j(j,:) / sum(j)
    !        end do
    !
    !        !Convergence criterion: If the molfractions don't change more than 10^-6 % the calculations converged
    !        max_dev = 0.D0
    !        Do j = 1, NrofPhases
    !            Do i = 1, ncomp
    !                if (abs((x_j(j,i) - x_old(j,i))/x_j(j,i)) > max_dev) then
    !                    max_dev = abs((x_j(j,i) - x_old(j,i))/x_j(j,i))
    !                end if
    !            end do
    !        End do
    !        x_old = x_j
    !        if (max_dev < 1.D-6) then
    !            x_vap = x_j(1,:)
    !            x_liq = x_j(2,:)
    !            return
    !        end if
    !
    !        !3) Calculate new K-values
    !        ! Recalculate the fugacity coefficients for each phase
    !        Do j = 1, NrofPhases
    !            molfractions = x_j(j,:)
    !            call reduced_parameters_calc
    !            Dens =  rhomix_calc(gl,TEMP, PRESS, 0.D0, 0,0)
    !            call FUGCO_CALC_MIX(Temp,Dens, fug_coeff(j,:))
    !        end do
    !
    !        !Calculate K values
    !        do j = 1, NrofPhases - 1
    !            do i = 1, ncomp
    !                K_ij(j,i) = fug_coeff(NrofPhases,i) / fug_coeff(j,i)
    !            end do
    !        end do
    !
    !   end do
    !   x_vap = x_j(1,:)
    !   x_liq = x_j(2,:)
    !   vapfrac = beta_j(1)
    !
    !!--------------------------------------------------------------------------
    !**************************************************************************

    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN     --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN     --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN     --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN     --  iFlash = 4
    !   - BUBBLE POINT: xk" AND x' VECTOR GIVEN   --  iFlash = 5
    !   - DEW POINT:    xk' AND x" VECTOR GIVEN   --  iFlash = 6
    !   - BUBBLE POINT: rho' AND x' VECTOR GIVEN  --  iFlash = 7
    !   - DEW POINT:    rho" AND x" VECTOR GIVEN  --  iFLash = 8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new, GibbsEQN_b!, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JT
    double precision, dimension(60):: GibbsEQN, Delta_X, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del, eps_gibbs2
    integer:: i, j, k, eqn
    ! warnings (Theresa) integer, dimension(1):: maxID
    character(255):: herr

    !z = x_known
    press_new = 0.d0
    errval = 0
    Delta_X = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished if the maximum difference of fugacities is lower than eps_gibbs2
    eps_del = 1.d-12
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs2 = 1.d-6
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible --> exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------

    if (iFlash < 3) then
        eqn = gl%ncomp+2
    else
        eqn = gl%ncomp+1
    end if


    do i = 1, 100
        call SysOfEqs_PhaseBoundary_Vbased(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)
        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if
        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_PhaseBoundary_Vbased(gl,press, temp, rhovap, rholiq, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JT, errval)
        Delta_X = - GibbsEQn

        !For all flash types in this routine the number of equations to be solved equals the number of components (equality of fugacities for each component)
        if (iFlash < 3) then
            eqn = gl%ncomp+2
        else
            eqn = gl%ncomp+1
        end if

        call LUdecomp(gl,eqn,JT,Delta_X,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + Delta_X(j)
                    rhovap_new = rhovap + Delta_X(gl%ncomp+1)
                    rholiq_new = rholiq + Delta_X(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)

                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(j)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + Delta_X(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + Delta_X(j)
                    rhovap_new = rhovap + Delta_X(gl%ncomp+1)
                    rholiq_new = rholiq + Delta_X(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(j)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    rhovap_new = rhovap + Delta_X(gl%ncomp)
                    rholiq_new = rholiq + Delta_X(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + Delta_X(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    rhovap_new = rhovap + Delta_X(gl%ncomp)
                    rholiq_new = rholiq + Delta_X(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + Delta_X(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + Delta_X(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + Delta_X(gl%ncomp - 1)
                        rhovap_new = rhovap + Delta_x(gl%ncomp)
                        rholiq_new = rholiq + Delta_x(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + Delta_X(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + Delta_X(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + Delta_X(gl%ncomp - 1)
                        rhovap_new = rhovap + Delta_x(gl%ncomp)
                        rholiq_new = rholiq + Delta_x(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! the 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (Delta_X(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(7)
            x_liq_new = x_liq
            Temp_new = Temp
            rholiq_new = rholiq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + Delta_X(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + Delta_X(gl%ncomp)
                    rhovap_new = rhovap + Delta_X(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(8)
            x_vap_new = x_vap
            Temp_new = Temp
            rhovap_new = rhovap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + Delta_X(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + Delta_X(gl%ncomp)
                    rholiq_new = rholiq + Delta_X(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(Delta_X(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        ! write the new values to the variables
        x_vap = x_vap_new!/sum_vap
        x_liq = x_liq_new!/sum_liq
        Temp = Temp_new
        press = press_new
        rhovap = rhovap_new
        rholiq = rholiq_new

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
        max_del = 0.D0

        if (iFlash < 3)then
            j = gl%ncomp+2
        else
            j = gl%ncomp+1
        end if

        Do k = 1, j
            if(abs(delta_X(k) / Var_X(k)) > max_del) then
                max_del = abs(delta_X(k) / Var_X(k))
            end if
        end do


        if ((max_del) < eps_del)then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if

        !if (((max_del) < eps_del) .and. (maxval(dabs(GibbsEQN)) < eps_Gibbs2)) then
        !    if (press_new == 0.D0) then
        !        molfractions = x_liq
        !        call reduced_parameters_calc(gl,Temp)
        !        press = P_CALC(gl,Temp, rholiq, 0)
        !        molfractions = z
        !        call reduced_parameters_calc(gl,Temp)
        !    end if
        !    return
        !end if

        !Catch unphysical temperatures and pressures
        if ((temp < 0.D0) .or. temp > 1000.D0) then
            errval = -4321
            exit
        end if
        if ((press < 0.D0) .or. press > 1000.D0) then
            errval = -4322
            exit
        end if
    end do

    ! If convergence was not reached after 30 iterations --> algorithm failed!
    if ((i > 100) .AND. (maxval(dabs(GibbsEQN)) > 1.D-6)) then
        errval = -2222
    End if

    if ((i > 100) .AND. (GibbsEQN_b < 1.D-6)) then
        if (press_new == 0.D0) then
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            press = P_CALC(gl,Temp, rholiq, 0)
            gl%molfractions = z
            call reduced_parameters_calc(gl,Temp)
        end if
    end if

    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_Vbased_calc
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_PhaseBoundary_Vbased(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS ON THE PHASE BOUNDARY.
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines









    implicit none

    type(type_gl) :: gl


    double precision:: P, T,  p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:

    !call newvars (P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    !In newvars wird der Dichtesolver aufgerufen mit den aktualisierten Zusammensetzungen / Drücken / Temperaturen --> Für v-based Berechnungen nicht nötig!!!

    !if (errorflag /= 0) then
    !    errval = errorflag
    !    return
    !end if

    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = rhovap
    p_vap = P_CALC(gl,T, d_vap, 0)
    call lnf_mix(gl,T, d_vap, p_vap, lnfi_vap)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = rholiq
    p_liq = P_CALC(gl,T, d_liq, 0)
    call lnf_mix(gl,T, d_liq, p_liq, lnfi_liq)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        return
    end if

    !!NEW VLE FOR DROPLETS
    !if (droplet) then
    !    !Pl > Pv, from module variable
    !    d_liq = rhomix_calc(gl,T, pl, 0.d0, 1, 1)
    !    call lnf_mix(T, d_liq, pl, lnfi_liq)
    !end if
    !!NEW VLE FOR DROPLETS

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do

    if (iFlash < 3) then
        GibbsEQN(gl%ncomp+1) = (p - p_vap) !* 1.D6
        GibbsEQN(gl%ncomp+2) = (p - p_liq) !* 1.D6
    else
        GibbsEQN(gl%ncomp+1) = p_vap - p_liq
    end if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_PhaseBoundary_Vbased
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_PhaseBoundary_Vbased (gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM. THIS ROUTINE GENERATES
    ! THE MATRIX NEEDED FOR PHASE BOUNDARY CALCULATIONS
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac     - Molar vapor fraction
    !   iFlash      - Flash mode (explanation, see routine "SysOfEqs")
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 6 or 7
    ! OUTPUT:
    !   errval  - Error value
    !   JT  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (iFlash 6 & 7)






    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    integer:: i, j
    double precision:: d_vap, d_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z, dlnphiidP_vap, dlnphiidP_liq, dlnfidT_liq, dlnfidT_vap, dlnfidrho_liq, dlnfidrho_vap, dPdXj_liq, dPdXj_vap
    double precision:: rhoredmix_orig, tredmix_orig

    JT = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. T,
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_vap(j)
            JT(J,gl%ncomp+2) = 0
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. T,
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = 0
            JT(J,gl%ncomp+2) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given.
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given.
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (5)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, xk" and x' vector are given.
        ! T, xi"(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = dPdXj_vap(j)
            else
                JT(j, gl%ncomp+1) = dPdXj_vap(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (6)
        !--------------------------------------------------------------------------
        ! Dew point calculation, xk' and x" vector are given.
        ! T, xi'(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = -dPdXj_liq(j)
            else
                JT(j, gl%ncomp+1) = -dPdXj_liq(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (7)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, rho' and x' vector are given.
        ! T, x", rho" need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = dlnfidXj_vap(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = dPdrho_vap

    case (8)
        !--------------------------------------------------------------------------
        ! Dew point calculation, rho" and x" vector are given.
        ! T, x', rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq = dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = -dlnfidXj_liq(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

        !--------------------------------------------------------------------------
        case default
        errval = -1111
        !--------------------------------------------------------------------------
    end select

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_PhaseBoundary_Vbased
    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN     --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN     --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN     --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN     --  iFlash = 4
    !   - BUBBLE POINT: xk" AND x' VECTOR GIVEN   --  iFlash = 5
    !   - DEW POINT:    xk' AND x" VECTOR GIVEN   --  iFlash = 6
    !   - BUBBLE POINT: rho' AND x' VECTOR GIVEN  --  iFlash = 7
    !   - DEW POINT:    rho" AND x" VECTOR GIVEN  --  iFLash = 8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, P_CALC
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JT, A, D, DT, J_matrix
    double precision, dimension(60):: GibbsEQN, GibbsEQN_new, h_lm, Var_X, g
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    ! warnings (Theresa) integer, dimension(1):: maxID
    character(255):: herr

    double precision:: mu, nu, GibbsEQN_b,GibbsEQN_b_new,rho,tol,h_lm_b, eps_Gibbs2
    double precision, dimension(2):: p_vec

    !z = x_known
    errval = 0
    h_lm = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished if the maximum difference of fugacities is lower than eps_gibbs2
    eps_del = 1.d-12
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs2 = 1.d-6
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible --> exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    press_new = 0.d0
    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------

    tol = 1.D-12
    mu = 1.D0
    nu = 2.D0

    !For all flash types in this routine the number of equations to be solved equals the number of components (equality of fugacities for each component)
    if (iFlash < 3) then
        eqn = gl%ncomp+2
    else
        eqn = gl%ncomp+1
    end if

    call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

    do i = 1, 100

        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return

        call Jacobi_PhaseBoundary_Vbased_LevMar(gl,press, temp, rhovap, rholiq, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JT, errval)

        J_matrix = transpose(JT)

        A = matmul(JT,J_matrix)
        call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
        g = matmul(JT,GibbsEQN)

        call D_rechner(gl,eqn,mu,A,D)

        DT = transpose(D)

        h_lm = - g

        call LUdecomp(gl,eqn,DT,h_lm,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        call h_lm_betrag(gl,eqn,h_lm,h_lm_b)

        !if (h_lm_b < eps_del) then
        !    return
        !end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_lm(j)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)
                    rholiq_new = rholiq + h_lm(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_lm(j)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)
                    rholiq_new = rholiq + h_lm(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    rhovap_new = rhovap + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    rhovap_new = rhovap + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + h_lm(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + h_lm(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + h_lm(gl%ncomp - 1)
                        rhovap_new = rhovap + h_lm(gl%ncomp)
                        rholiq_new = rholiq + h_lm(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + h_lm(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + h_lm(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + h_lm(gl%ncomp - 1)
                        rhovap_new = rhovap + h_lm(gl%ncomp)
                        rholiq_new = rholiq + h_lm(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! the 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (h_lm(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(7)
            x_liq_new = x_liq
            Temp_new = Temp
            rholiq_new = rholiq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_lm(gl%ncomp)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(8)
            x_vap_new = x_vap
            Temp_new = Temp
            rhovap_new = rhovap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press_new, Temp_new, x_vap_new, x_liq_new, rhovap_new, rholiq_new, vapfrac, iFlash, GibbsEQN_new, errval)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN_new,GibbsEQN_b_new)

        call rho_rechner(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,mu,h_lm,g,rho)

        if (rho > 0) then
            !______________________________________________________________________________________________
            !x = x_new:

            ! write the new values to the variables
            x_vap = x_vap_new!/sum_vap
            x_liq = x_liq_new!/sum_liq
            Temp = Temp_new
            press = press_new
            rhovap = rhovap_new
            rholiq = rholiq_new

            !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
            max_del = 0.D0

            if (iFlash < 3)then
                j = gl%ncomp+2
            else
                j = gl%ncomp+1
            end if

            Do k = 1, j
                if(abs(h_lm(k) / Var_X(k)) > max_del) then
                    max_del = abs(h_lm(k) / Var_X(k))
                end if
            end do

            !if ((max_del) < eps_del .and. (maxval(dabs(GibbsEQN)) < eps_Gibbs2)) then
            !    if (press_new == 0.D0) then
            !        molfractions = x_liq
            !        call reduced_parameters_calc(gl,Temp)
            !        press = P_CALC(gl,Temp, rholiq, 0)
            !        molfractions = z
            !        call reduced_parameters_calc(gl,Temp)
            !    end if
            !    return
            !end if

            if ((max_del) < eps_del) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if

            !!Catch unphysical temperatures and pressures
            !if ((temp < 0.D0) .or. temp > 1000.D0) then
            ! errval = -4321
            ! exit
            !end if
            !if ((press < 0.D0) .or. press > 1000.D0) then
            ! errval = -4322
            ! exit
            !end if
            !_______________________________________________________________________________________________________

            call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
            call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

            if(GibbsEQN_b <= tol) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if

            p_vec = 0.D0

            p_vec(1) = 1.D0/3.D0
            p_vec(2) = 1.D0 - (2.D0*rho - 1.D0)**3

            mu = mu * maxval(p_vec)
            nu = 2
        else
            mu = mu*nu
            nu = 2*nu
        end if
    end do

    ! If convergence was not reached after 30 iterations --> algorithm failed!
    if ((i > 100) .AND. (maxval(dabs(GibbsEQN)) > 1.D-6)) then
        errval = -2222
    End if

    if ((i > 100) .AND. (GibbsEQN_b < 1.D-6)) then
        if (press_new == 0.D0) then
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            press = P_CALC(gl,Temp, rholiq, 0)
            gl%molfractions = z
            call reduced_parameters_calc(gl,Temp)
        end if
    end if
    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_Vbased_calc_LevMar
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS ON THE PHASE BOUNDARY.
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines









    implicit none

    type(type_gl) :: gl


    double precision:: P, T, p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:

    !call newvars (P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    !In newvars wird der Dichtesolver aufgerufen mit den aktualisierten Zusammensetzungen / Drücken / Temperaturen --> Für v-based Berechnungen nicht nötig!!!

    !if (errorflag /= 0) then
    !    errval = errorflag
    !    return
    !end if

    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = rhovap
    p_vap = P_CALC(gl,T, d_vap, 0)
    call lnf_mix(gl,T, d_vap, p_vap, lnfi_vap)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = rholiq
    p_liq = P_CALC(gl,T, d_liq, 0)
    call lnf_mix(gl,T, d_liq, p_liq, lnfi_liq)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        return
    end if

    !!NEW VLE FOR DROPLETS
    !if (droplet) then
    !    !Pl > Pv, from module variable
    !    d_liq = rhomix_calc(gl,T, pl, 0.d0, 1, 1)
    !    call lnf_mix(T, d_liq, pl, lnfi_liq)
    !end if
    !!NEW VLE FOR DROPLETS

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do

    if (iFlash < 3) then
        GibbsEQN(gl%ncomp+1) = (p - p_vap) !* 1.D6
        GibbsEQN(gl%ncomp+2) = (p - p_liq) !* 1.D6
    else
        GibbsEQN(gl%ncomp+1) = p_vap - p_liq
    end if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_PhaseBoundary_Vbased_LevMar
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_PhaseBoundary_Vbased_LevMar(gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM. THIS ROUTINE GENERATES
    ! THE MATRIX NEEDED FOR PHASE BOUNDARY CALCULATIONS
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac     - Molar vapor fraction
    !   iFlash      - Flash mode (explanation, see routine "SysOfEqs")
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 6 or 7
    ! OUTPUT:
    !   errval  - Error value
    !   JT  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (iFlash 6 & 7)






    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    integer:: i, j
    double precision:: d_vap, d_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z, dlnphiidP_vap, dlnphiidP_liq, dlnfidT_liq, dlnfidT_vap, dlnfidrho_liq, dlnfidrho_vap, dPdXj_liq, dPdXj_vap
    double precision:: rhoredmix_orig, tredmix_orig

    JT = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. T,
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_vap(j)
            JT(J,gl%ncomp+2) = 0
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. T,
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = 0
            JT(J,gl%ncomp+2) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given.
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given.
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (5)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, xk" and x' vector are given.
        ! T, xi"(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = dPdXj_vap(j)
            else
                JT(j, gl%ncomp+1) = dPdXj_vap(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (6)
        !--------------------------------------------------------------------------
        ! Dew point calculation, xk' and x" vector are given.
        ! T, xi'(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = -dPdXj_liq(j)
            else
                JT(j, gl%ncomp+1) = -dPdXj_liq(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (7)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, rho' and x' vector are given.
        ! T, x", rho" need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = dlnfidXj_vap(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = dPdrho_vap

    case (8)
        !--------------------------------------------------------------------------
        ! Dew point calculation, rho" and x" vector are given.
        ! T, x', rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq = dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = -dlnfidXj_liq(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

        !--------------------------------------------------------------------------
        case default
        errval = -1111
        !--------------------------------------------------------------------------
    end select

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_PhaseBoundary_Vbased_LevMar
    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_PhaseBoundary_Vbased_calc_LevMar_mix(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter, GibbsEQN_b)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN     --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN     --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN     --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN     --  iFlash = 4
    !   - BUBBLE POINT: xk" AND x' VECTOR GIVEN   --  iFlash = 5
    !   - DEW POINT:    xk' AND x" VECTOR GIVEN   --  iFlash = 6
    !   - BUBBLE POINT: rho' AND x' VECTOR GIVEN  --  iFlash = 7
    !   - DEW POINT:    rho" AND x" VECTOR GIVEN  --  iFLash = 8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new!, P_CALC
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JT, A, D, DT, J_matrix
    double precision, dimension(60):: GibbsEQN, GibbsEQN_new, h_lm, Var_X, g
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    ! warnings (Theresa) integer, dimension(1):: maxID
    character(255):: herr

    double precision:: mu, nu, GibbsEQN_b,GibbsEQN_b_new,rho,tol,h_lm_b, eps_Gibbs2
    double precision, dimension(2):: p_vec

    !z = x_known
    errval = 0
    h_lm = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished if the maximum difference of fugacities is lower than eps_gibbs2
    eps_del = 1.d-12
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs2 = 1.d-6
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible --> exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    press_new = 0.d0
    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------

    tol = 1.D-12
    mu = 1.D-10
    nu = 2.D0

    !For all flash types in this routine the number of equations to be solved equals the number of components (equality of fugacities for each component)
    if (iFlash < 3) then
        eqn = gl%ncomp+2
    else
        eqn = gl%ncomp+1
    end if

    call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

    do i = 1, 100

        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return

        call Jacobi_PhaseBoundary_Vbased_LevMar(gl,press, temp, rhovap, rholiq, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JT, errval)

        J_matrix = transpose(JT)

        A = matmul(JT,J_matrix)
        call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
        g = matmul(JT,GibbsEQN)

        call D_rechner(gl,eqn,mu,A,D)

        DT = transpose(D)

        h_lm = - g

        call LUdecomp(gl,eqn,DT,h_lm,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        call h_lm_betrag(gl,eqn,h_lm,h_lm_b)

        !if (h_lm_b < eps_del) then
        !    return
        !end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_lm(j)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)
                    rholiq_new = rholiq + h_lm(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_lm(j)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)
                    rholiq_new = rholiq + h_lm(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    rhovap_new = rhovap + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    rhovap_new = rhovap + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + h_lm(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + h_lm(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + h_lm(gl%ncomp - 1)
                        rhovap_new = rhovap + h_lm(gl%ncomp)
                        rholiq_new = rholiq + h_lm(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + h_lm(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + h_lm(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + h_lm(gl%ncomp - 1)
                        rhovap_new = rhovap + h_lm(gl%ncomp)
                        rholiq_new = rholiq + h_lm(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! the 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (h_lm(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(7)
            x_liq_new = x_liq
            Temp_new = Temp
            rholiq_new = rholiq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_lm(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_lm(gl%ncomp)
                    rhovap_new = rhovap + h_lm(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(8)
            x_vap_new = x_vap
            Temp_new = Temp
            rhovap_new = rhovap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_lm(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_lm(gl%ncomp)
                    rholiq_new = rholiq + h_lm(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_lm(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_lm = h_lm/stepsize
                    if (maxval(dabs(h_lm)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press_new, Temp_new, x_vap_new, x_liq_new, rhovap_new, rholiq_new, vapfrac, iFlash, GibbsEQN_new, errval)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN_new,GibbsEQN_b_new)

        call rho_rechner(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,mu,h_lm,g,rho)

        if (rho > 0) then
            !______________________________________________________________________________________________
            !x = x_new:

            ! write the new values to the variables
            x_vap = x_vap_new!/sum_vap
            x_liq = x_liq_new!/sum_liq
            Temp = Temp_new
            press = press_new
            rhovap = rhovap_new
            rholiq = rholiq_new

            !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
            max_del = 0.D0

            if (iFlash < 3)then
                j = gl%ncomp+2
            else
                j = gl%ncomp+1
            end if

            Do k = 1, j
                if(abs(h_lm(k) / Var_X(k)) > max_del) then
                    max_del = abs(h_lm(k) / Var_X(k))
                end if
            end do

            if ((max_del) < eps_del .and. (maxval(dabs(GibbsEQN)) < eps_Gibbs2)) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if

            !!Catch unphysical temperatures and pressures
            !if ((temp < 0.D0) .or. temp > 1000.D0) then
            ! errval = -4321
            ! exit
            !end if
            !if ((press < 0.D0) .or. press > 1000.D0) then
            ! errval = -4322
            ! exit
            !end if
            !_______________________________________________________________________________________________________

            call SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
            call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

            if(GibbsEQN_b <= tol) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if

            p_vec = 0.D0

            p_vec(1) = 1.D0/3.D0
            p_vec(2) = 1.D0 - (2.D0*rho - 1.D0)**3

            mu = mu * maxval(p_vec)
            nu = 2
        else
            mu = mu*nu
            nu = 2*nu
        end if

        !Wenn gewisse Anzahl an Iterationen erreicht (hier 20) verlasse den Algorithmus und übergebe die Werte an Gauß-Newton
        !Oder wenn gewisse Genauigkeit erreicht
        if ((i > 20) .OR. (GibbsEQN_b < 1.D-4)) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
        end if

    end do

    if ((i > 100) .AND. (GibbsEQN_b < 1.D-6)) then
        if (press_new == 0.D0) then
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            press = P_CALC(gl,Temp, rholiq, 0)
            gl%molfractions = z
            call reduced_parameters_calc(gl,Temp)
        end if
    end if
    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_Vbased_calc_LevMar_mix
    !**************************************************************************



    !**************************************************************************
    module subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN     --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN     --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN     --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN     --  iFlash = 4
    !   - BUBBLE POINT: xk" AND x' VECTOR GIVEN   --  iFlash = 5
    !   - DEW POINT:    xk' AND x" VECTOR GIVEN   --  iFlash = 6
    !   - BUBBLE POINT: rho' AND x' VECTOR GIVEN  --  iFlash = 7
    !   - DEW POINT:    rho" AND x" VECTOR GIVEN  --  iFLash = 8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JT
    double precision, dimension(60):: GibbsEQN, h_dl, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    ! warnings (Theresa) integer, dimension(1):: maxID
    character(255):: herr

    double precision:: tol,delta,alpha,h_gn_b,h_sd_b,g_b,beta_lc,h_dl_b,GibbsEQN_b,GibbsEQN_b_new,rho
    double precision, dimension(2):: q
    double precision, dimension(60)::g,h_sd,h_gn,GibbsEQN_new
    double precision, dimension(60,60)::J_matrix,A,AT

    !z = x_known
    errval = 0
    h_dl = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-12
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible --> exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    press_new = 0.d0
    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------

    if (iFlash < 3) then
        eqn = gl%ncomp+2
    else
        eqn = gl%ncomp+1
    end if

    delta = 10000.D0
    tol = 1.D-12

    call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

    do i = 1, 100

        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_PhaseBoundary_Vbased_DogLeg(gl,press, temp, rhovap, rholiq, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JT, errval)

        J_matrix = transpose(JT)

        g = matmul(JT,GibbsEQN)

        call alpha_rechner(gl,eqn,g,J_matrix,alpha)

        h_sd = -alpha * g

        !solve g_n:

        A = matmul(JT,J_matrix)
        AT = transpose(A)
        h_gn = - g

        call LUdecomp(gl,eqn,AT,h_gn,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        !______________________ compute h_dl ______________________ -->

        call h_gn_betrag(gl,eqn,h_gn,h_gn_b)

        call h_sd_betrag(gl,eqn,h_sd,h_sd_b)

        if(h_gn_b <= delta)then
            h_dl = h_gn

        else if(h_sd_b >= delta) then
            call ohnenamen(gl,eqn,delta,g,g_b,h_dl)

        else
            call beta_rechner(gl,eqn,h_sd,h_gn,delta,beta_lc)

            h_dl = h_sd + beta_lc * (h_gn - h_sd)

        end if

        !______________________ compute h_dl ende______________________

        call h_dl_betrag(gl,eqn,h_dl,h_dl_b)

        if (h_dl_b < eps_del) then
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_dl(j)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)
                    rholiq_new = rholiq + h_dl(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_dl(j)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)
                    rholiq_new = rholiq + h_dl(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    rhovap_new = rhovap + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0


            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    rhovap_new = rhovap + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0

            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + h_dl(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + h_dl(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + h_dl(gl%ncomp - 1)
                        rhovap_new = rhovap + h_dl(gl%ncomp)
                        rholiq_new = rholiq + h_dl(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + h_dl(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + h_dl(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + h_dl(gl%ncomp - 1)
                        rhovap_new = rhovap + h_dl(gl%ncomp)
                        rholiq_new = rholiq + h_dl(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! the 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (h_dl(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0

            ! --------------------------------------------------------------------------------
        case(7)
            x_liq_new = x_liq
            Temp_new = Temp
            rholiq_new = rholiq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_dl(gl%ncomp)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(8)
            x_vap_new = x_vap
            Temp_new = Temp
            rhovap_new = rhovap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press_new, Temp_new, x_vap_new, x_liq_new, rhovap_new, rholiq_new, vapfrac, iFlash, GibbsEQN_new, errval)
        call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN_new,GibbsEQN_b_new)

        call rho_rechner_DogLeg(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,h_dl,g,rho)

        if (rho > 0)then

            ! write the new values to the variables
            x_vap = x_vap_new!/sum_vap
            x_liq = x_liq_new!/sum_liq
            Temp = Temp_new
            press = press_new
            rhovap = rhovap_new
            rholiq = rholiq_new

            !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
            max_del = 0.D0

            if (iFlash < 3)then
                j = gl%ncomp+2
            else
                j = gl%ncomp+1
            end if

            Do k = 1, j
                if(abs(h_dl(k) / Var_X(k)) > max_del) then
                    max_del = abs(h_dl(k) / Var_X(k))
                end if
            end do

            if ((max_del) < eps_del) return

            !Catch unphysical temperatures and pressures
            !if ((temp < 0.D0) .or. temp > 1000.D0) then
            ! errval = -4321
            ! exit
            !end if
            !if ((press < 0.D0) .or. press > 1000.D0) then
            ! errval = -4322
            ! exit
            !end if
            call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
            call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

            if(GibbsEQN_b <= tol) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if
        end if

        if (rho > 0.75) then
            call q_rechner(gl,delta, h_dl_b, q)
            delta = maxval(q)
        else if (rho < 0.25) then
            delta = delta / 2
        end if

    end do

    ! If convergence was not reached after 30 iterations --> algorithm failed!
    if ((i > 100) .AND. (maxval(dabs(GibbsEQN)) > 1.D-6)) then
        errval = -2222
    End if

    if ((i > 100) .AND. (GibbsEQN_b < 1.D-6)) then
        if (press_new == 0.D0) then
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            press = P_CALC(gl,Temp, rholiq, 0)
            gl%molfractions = z
            call reduced_parameters_calc(gl,Temp)
        end if
    end if

    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg
    !**************************************************************************

    !**************************************************************************
    module subroutine SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS ON THE PHASE BOUNDARY.
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !   - BUBBLE POINT: xi" AND x' VECTOR GIVEN --  iFlash = 5
    !   - DEW POINT:    xi' and x" VECTOR GIVEN --  iFlash = 6
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   rhovap_est  - Estimated vapor phase density
    !   rholiq_est  - Estimated liquid phase density
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try      -   Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                       1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                       2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    ! OUTPUT:
    !   errval      - Error value
    !   GibbsEQN    - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines









    implicit none

    type(type_gl) :: gl


    double precision:: P, T, p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag

    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! write the vapor and liquid phase density and the pure fluid densities
    ! into the module variables for the new set of T, p, x_vap, and x_liq:

    !call newvars (P, T, x_vap, x_liq, rhovap_est, rholiq_est, Iphase_try, errorflag)
    !In newvars wird der Dichtesolver aufgerufen mit den aktualisierten Zusammensetzungen / Drücken / Temperaturen --> Für v-based Berechnungen nicht nötig!!!

    !if (errorflag /= 0) then
    !    errval = errorflag
    !    return
    !end if

    !----------------------------------
    ! get the gas phase properties
    z = gl%molfractions
    gl%molfractions = x_vap
    call reduced_parameters_calc(gl,T)
    ! get the gas phase density from the module
    d_vap = rhovap
    p_vap = P_CALC(gl,T, d_vap, 0)
    call lnf_mix(gl,T, d_vap, p_vap, lnfi_vap)

    !Errorhandling: Calculation of the fugacities of the vapor phase failed
    if (lnfi_vap(1) == 0.D0) then
        errval = -7777
    end if

    !----------------------------------
    ! get the liquid phase properties
    gl%molfractions = x_liq
    call reduced_parameters_calc(gl,T)
    ! get the liquid phase density from the module
    d_liq = rholiq
    p_liq = P_CALC(gl,T, d_liq, 0)
    call lnf_mix(gl,T, d_liq, p_liq, lnfi_liq)

    !Errorhandling: Calculation of the fugacities of the liquid phase failed
    if (lnfi_liq(1) == 0.D0) then
        errval = -7778
        return
    end if

    !!NEW VLE FOR DROPLETS
    !if (droplet) then
    !    !Pl > Pv, from module variable
    !    d_liq = rhomix_calc(gl,T, pl, 0.d0, 1, 1)
    !    call lnf_mix(T, d_liq, pl, lnfi_liq)
    !end if
    !!NEW VLE FOR DROPLETS

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    do i = 1, gl%ncomp
        ! The first n equations:
        GibbsEQN(i) = lnfi_vap(i) - lnfi_liq(i)
        ! the nth equation
    end do

    if (iFlash < 3) then
        GibbsEQN(gl%ncomp+1) = (p - p_vap) !* 1.D6
        GibbsEQN(gl%ncomp+2) = (p - p_liq) !* 1.D6
    else
        GibbsEQN(gl%ncomp+1) = p_vap - p_liq
    end if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine SysOfEqs_PhaseBoundary_Vbased_DogLeg
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_PhaseBoundary_Vbased_DogLeg(gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE JACOBI MATRIX OF THE SYSTEM OF EQUATIONS
    ! FOR THE GIBBS FREE ENERGY MINIMIZATION ALGORITHM. THIS ROUTINE GENERATES
    ! THE MATRIX NEEDED FOR PHASE BOUNDARY CALCULATIONS
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 20075
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P           - Pressure
    !   T           - Temperature
    !   x_vap       - Vapor phase composition
    !   x_liq       - Liquid phase composition
    !   vapfrac     - Molar vapor fraction
    !   iFlash      - Flash mode (explanation, see routine "SysOfEqs")
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 6 or 7
    ! OUTPUT:
    !   errval  - Error value
    !   JT  - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !             F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! J. Gernert, Jan. 2011
    ! A. Jäger, March 2011 (iFlash 6 & 7)






    implicit none

    type(type_gl) :: gl


    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    integer:: i, j
    double precision:: d_vap, d_liq
    double precision, dimension(30, 30):: dlnfidXj_vap, dlnfidXj_liq
    double precision, dimension(30):: dlnphiidT_liq, dlnphiidT_vap, z, dlnphiidP_vap, dlnphiidP_liq, dlnfidT_liq, dlnfidT_vap, dlnfidrho_liq, dlnfidrho_vap, dPdXj_liq, dPdXj_vap
    double precision:: rhoredmix_orig, tredmix_orig

    JT = 0.D0
    z = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, p and x' vector are given. T,
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_vap(j)
            JT(J,gl%ncomp+2) = 0
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (2)
        !--------------------------------------------------------------------------
        ! Dew point calculation, p and x" vector are given. T,
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+2,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = 0
            JT(J,gl%ncomp+2) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = -dPdT_vap
        JT(gl%ncomp,gl%ncomp+2) = -dPdT_liq

        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+2) = 0

        JT(gl%ncomp+2,gl%ncomp+1) = 0
        JT(gl%ncomp+2,gl%ncomp+2) = -dPdrho_liq

    case (3)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, T and x' vector are given.
        ! x", rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = dlnfidXj_vap(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (4)
        !--------------------------------------------------------------------------
        ! Dew point calculation, T and x" vector are given.
        ! x', rho", rho' need to be calculated
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-1   ! goes through the n-1 rows, which are the derivativs to x_j
            do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                JT(j, i) = -dlnfidXj_liq(j, i)
            end do
        end do
        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (5)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, xk" and x' vector are given.
        ! T, xi"(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------
        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = dlnfidXj_vap(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = dPdXj_vap(j)
            else
                JT(j, gl%ncomp+1) = dPdXj_vap(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (6)
        !--------------------------------------------------------------------------
        ! Dew point calculation, xk' and x" vector are given.
        ! T, xi'(i=/k), rho", rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = gl%rho_liq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq =  dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = gl%rho_vap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap )
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        ! calculate the derivatives of the first n equations w.r.t the composition
        do j = 1, gl%ncomp-2   ! goes through the n-1 rows, which are the derivativs to x_j
            if (j < Nr_x_given) then
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j, i)
                end do
            else
                do i = 1, gl%ncomp   ! goes through the first n-1 columns, which are the functions F_i
                    JT(j, i) = -dlnfidXj_liq(j+1, i) !lässt quasi die Komponente aus, die bekannt ist bei den Ableitungen
                end do
            End if
        end do
        do i = 1, gl%ncomp   ! goes through the first n-1 F_i and calcs their derivative w.r.t. temperature
            JT(gl%ncomp-1, i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp, i) = dlnfidrho_vap(i)
            JT(gl%ncomp+1, i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-2
            if (j < Nr_x_given) then
                JT(j, gl%ncomp+1) = -dPdXj_liq(j)
            else
                JT(j, gl%ncomp+1) = -dPdXj_liq(j+1)
            end if
        end do

        JT(gl%ncomp-1,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp,gl%ncomp+1) = dPdrho_vap
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

    case (7)
        !--------------------------------------------------------------------------
        ! Bubble point calculation, rho' and x' vector are given.
        ! T, x", rho" need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_vap, dlnfidXj_vap)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dlnfi_drho_TX(gl,T, d_vap, dlnfidrho_vap)
        call dP_dXi_Trho(gl,T, d_vap, dPdXj_vap)
        dPdXj_vap = dPdXj_vap / 1.D6
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6
        call dP_drho (gl,T, d_vap, dPdrho_vap)
        dPdrho_vap = dPdrho_vap / 1.D6

        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = dlnfidXj_vap(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = dlnfidrho_vap(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = dPdXj_vap(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = dPdrho_vap

    case (8)
        !--------------------------------------------------------------------------
        ! Dew point calculation, rho" and x" vector are given.
        ! T, x', rho' need to be calculated.
        !--------------------------------------------------------------------------

        ! get the liquid phase properties
        d_liq = rholiq
        z = gl%molfractions
        gl%molfractions = x_liq
        call reduced_parameters_calc(gl,T)
        call dlnfi_dXi_TrhoXj(gl,T, d_liq, dlnfidXj_liq)
        call dlnfi_drho_TX(gl,T, d_liq, dlnfidrho_liq)
        call dlnfi_dT_rhoX(gl,T, d_liq, dlnfidT_liq)
        call dP_dT (gl,T, d_liq, dPdT_liq)
        dPdT_liq = dPdT_liq / 1.D6
        call dP_dXi_Trho(gl,T, d_liq, dPdXj_liq)
        dPdXj_liq = dPdXj_liq / 1.D6
        call dP_drho (gl,T, d_liq, dPdrho_liq)
        dPdrho_liq = dPdrho_liq / 1.D6

        ! get the vapor phase properties
        d_vap = rhovap
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,T)
        call dlnfi_dT_rhoX(gl,T, d_vap, dlnfidT_vap)
        call dP_dT (gl,T, d_vap, dPdT_vap)
        dPdT_vap = dPdT_vap / 1.D6


        do j = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                JT(j,i) = -dlnfidXj_liq(j,i)
            end do
        end do

        do i = 1, gl%ncomp
            JT(gl%ncomp,i) = dlnfidT_vap(i) - dlnfidT_liq(i)
            JT(gl%ncomp+1,i) = -dlnfidrho_liq(i)
        end do

        do j = 1, gl%ncomp-1
            JT(j,gl%ncomp+1) = -dPdXj_liq(j)
        end do

        JT(gl%ncomp,gl%ncomp+1) = dPdT_vap - dPdT_liq
        JT(gl%ncomp+1,gl%ncomp+1) = -dPdrho_liq

        !--------------------------------------------------------------------------
        case default
        errval = -1111
        !--------------------------------------------------------------------------
    end select

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine Jacobi_PhaseBoundary_Vbased_DogLeg
    !**************************************************************************

    !**************************************************************************
    module subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg_mix(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE PHASE BOUNDARY EQUILIBRIA
    ! (DEW POINT FLASH ROUTINES AND BUBBLE POINT FLASH ROUTINES)
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATIONS:
    !--------------------------------------------------------------------------
    !           Iglesias-Silva at al.,
    !           Fluid Phase Equilibria 210 (2003), 229-245
    !           Michelson, M.L. ; Mollerup, J.M.
    !           "Thermodynamic Models: Fundamentals & Computational Aspects"
    !           Tie-Line Publications, Denmark 2004
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - BUBBLE POINT: P AND x' VEXTOR GIVEN     --  iFlash = 1
    !   - DEW POINT:    P AND x" VECTOR GIVEN     --  iFlash = 2
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN     --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN     --  iFlash = 4
    !   - BUBBLE POINT: xk" AND x' VECTOR GIVEN   --  iFlash = 5
    !   - DEW POINT:    xk' AND x" VECTOR GIVEN   --  iFlash = 6
    !   - BUBBLE POINT: rho' AND x' VECTOR GIVEN  --  iFlash = 7
    !   - DEW POINT:    rho" AND x" VECTOR GIVEN  --  iFLash = 8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x_known     - Overall composion
    !   x_vap       - Vapor phase composition (if empty startvalues will be generated)
    !   x_liq       - Liquid phase composition (if empty startvalues will be generated)
    !   rhovap_est  - Estimated vapor phase density (if not given, the PSRK will be used for startvalues)
    !   rholiq_est  - Estimated liquid phase density (if not given, the PSRK will be used for startvalues)
    !   vapfrac        - Molar vapor fraction
    !   iFlash      - Flash mode
    !   iPhase_try  - Necessary to distinguish between vapor / liquid and liquid / liquid equilibria  0 : Try iPhase 0 for both phases (Let the density solver try to find the correct equilibrium)
    !                                                                                                 1 : Try iPhase 1 for vapor phase  (Assume liquid / liquid equilibrium)
    !                                                                                                 2 : Try iphase 2 for vapor phase  (Assume vapor / liquid equilibrium)
    !   Nr_x_given  - Nr.(i) of the component given in case of iflash 5 or 6
    ! OUTPUT:
    !   errval  - Error value
    !   GibbsEQN- 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !------------------0--------------------------------------------------------
    ! A. Jäger November 2011: Old routine containing all phase equilibrium calculations was split into several routines








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter

    double precision, dimension(60, 60):: JT
    double precision, dimension(60):: GibbsEQN, h_dl, Var_X
    double precision:: sum_vap, sum_liq, stepsize, Temp_new, press_new
    double precision:: eps_Gibbs, eps_del, max_del
    integer:: i, j, k, eqn
    ! warnings (Theresa) integer, dimension(1):: maxID
    character(255):: herr

    double precision:: tol,delta,alpha,h_gn_b,h_sd_b,g_b,beta_lc,h_dl_b,GibbsEQN_b,GibbsEQN_b_new,rho
    double precision, dimension(2):: q
    double precision, dimension(60)::g,h_sd,h_gn,GibbsEQN_new
    double precision, dimension(60,60)::J_matrix,A,AT

    !z = x_known
    errval = 0
    h_dl = 1.D0
    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
    eps_del = 1.d-12
    !write the overall composition into the module variable.
    z = gl%molfractions
    gl%molfractions = x_known

    ! -------------------------------------------------------------
    !  GENERATION OF START VALUES FOR THE GIBBS MINIMIZATION
    ! -------------------------------------------------------------
    ! if no start values are given for one of the phases calculate start values from Wilson K-Factors
    if ((x_vap(1) == 0.D0) .OR. (x_liq(1) == 0.D0)) then
        !In case of iFlash = 5 or 6, estimation of start values is not possible --> exit with errorcode
        if (iFlash < 5) then
            call PTX_startvals_PhaseBoundary (gl,press, Temp, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
        else
            errval = -1111
            return
        end if
    end if
    !check if start value generation failed
    if (errval < 0) return

    press_new = 0.d0
    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0
    !------------------------------------------------------------
    !   GIBBS FREE ENERGY MINIMIZATION ALGORITHM
    !------------------------------------------------------------

    if (iFlash < 3) then
        eqn = gl%ncomp+2
    else
        eqn = gl%ncomp+1
    end if

    delta = 1.D0
    tol = 1.D-12

    call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

    do i = 1, 100

        ! This is the exit condition for the VLE iteration!
        ! If the residua of the equations to be solved are lower than eps_Gibbs, the algorithm converged
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (press_new == 0.D0) then
                gl%molfractions = x_liq
                call reduced_parameters_calc(gl,Temp)
                press = P_CALC(gl,Temp, rholiq, 0)
                gl%molfractions = z
                call reduced_parameters_calc(gl,Temp)
            end if
            return
        end if

        !In case an error occurs in SysOfEqs, exit the routine
        if (errval /= 0) return
        call Jacobi_PhaseBoundary_Vbased_DogLeg(gl,press, temp, rhovap, rholiq, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JT, errval)

        J_matrix = transpose(JT)

        g = matmul(JT,GibbsEQN)

        call alpha_rechner(gl,eqn,g,J_matrix,alpha)

        h_sd = -alpha * g

        !solve g_n:

        A = matmul(JT,J_matrix)
        AT = transpose(A)
        h_gn = - g

        call LUdecomp(gl,eqn,AT,h_gn,errval,herr)
        if (errval == 1) then
            !write(*,*) 'Error in LUdecomp: ', herr
            errval = -4444
            return
        end if

        !______________________ compute h_dl ______________________ -->

        call h_gn_betrag(gl,eqn,h_gn,h_gn_b)

        call h_sd_betrag(gl,eqn,h_sd,h_sd_b)

        if(h_gn_b <= delta)then
            h_dl = h_gn

        else if(h_sd_b >= delta) then
            call ohnenamen(gl,eqn,delta,g,g_b,h_dl)

        else
            call beta_rechner(gl,eqn,h_sd,h_gn,delta,beta_lc)

            h_dl = h_sd + beta_lc * (h_gn - h_sd)

        end if

        !______________________ compute h_dl ende______________________

        call h_dl_betrag(gl,eqn,h_dl,h_dl_b)

        if (h_dl_b < eps_del) then
            return
        end if

        ! reset the new set of variables
        sum_vap = 0.D0
        sum_liq = 0.D0
        x_vap_new = 0.D0
        x_liq_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0

        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(1)
            x_liq_new = x_liq
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_dl(j)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)
                    rholiq_new = rholiq + h_dl(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(2)
            x_vap_new = x_vap
            press_new = press
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_dl(j)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)
                    rholiq_new = rholiq + h_dl(gl%ncomp+2)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                    Var_X(gl%ncomp+2) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(j)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    Temp_new = Temp
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        case(3)
            x_liq_new = x_liq
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    rhovap_new = rhovap + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0


            ! --------------------------------------------------------------------------------
        case(4)
            x_vap_new = x_vap
            Temp_new = Temp
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    rhovap_new = rhovap + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = rhovap_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Dp > 1 MPa
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0

            ! --------------------------------------------------------------------------------
        case(5)
            x_liq_new = x_liq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_vap_new(j) = x_vap(j) + h_dl(j)
                            Var_X(j) = x_vap_new(j)
                        else
                            x_vap_new(j) = x_vap(j) + h_dl(j-1)
                            Var_X(j-1) = x_vap_new(j)
                        end if

                    else
                        x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                        Temp_new = Temp + h_dl(gl%ncomp - 1)
                        rhovap_new = rhovap + h_dl(gl%ncomp)
                        rholiq_new = rholiq + h_dl(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_vap_new(j) = x_vap(j)
                end if

                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP - 1)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0

            ! --------------------------------------------------------------------------------
        case(6)
            x_vap_new = x_vap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j /= Nr_x_given) then
                    if (j < gl%ncomp) then
                        if (j < Nr_x_given) then
                            x_liq_new(j) = x_liq(j) + h_dl(j)
                            Var_X(j) = x_liq_new(j)
                        else
                            x_liq_new(j) = x_liq(j) + h_dl(j-1)
                            Var_X(j-1) = x_liq_new(j)
                        end if
                    else
                        x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                        Temp_new = Temp + h_dl(gl%ncomp - 1)
                        rhovap_new = rhovap + h_dl(gl%ncomp)
                        rholiq_new = rholiq + h_dl(gl%ncomp+1)

                        Var_X(gl%ncomp-1) = Temp_new
                        Var_X(gl%ncomp) = rhovap_new
                        Var_X(gl%ncomp+1) = rholiq_new
                    end if
                else    !The j-th mole fraction stays the same
                    x_liq_new(j) = x_liq(j)
                end if

                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! the 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (h_dl(gl%NCOMP - 1) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0

            ! --------------------------------------------------------------------------------
        case(7)
            x_liq_new = x_liq
            Temp_new = Temp
            rholiq_new = rholiq
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_vap_new(j) = x_vap(j) + h_dl(j)
                    Var_X(j) = x_vap_new(j)
                else
                    x_vap_new(gl%ncomp) = 1.D0 - sum_vap
                    Temp_new = Temp + h_dl(gl%ncomp)
                    rhovap_new = rhovap + h_dl(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rhovap_new
                end if
                sum_vap = sum_vap + x_vap_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_liq < 1
                ! The 2nd condition: x(j)_liq > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_vap_new(j) > 1.D0) .OR. (x_vap_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_vap = 0.D0
                end if
                j = j + 1
            end do
            sum_liq = 1.D0
            ! --------------------------------------------------------------------------------
        case(8)
            x_vap_new = x_vap
            Temp_new = Temp
            rhovap_new = rhovap
            j = 1
            do while (j < gl%ncomp + 1)
                if (j < gl%ncomp) then
                    x_liq_new(j) = x_liq(j) + h_dl(j)
                    Var_X(j) = x_liq_new(j)
                else
                    x_liq_new(gl%ncomp) = 1.D0 - sum_liq
                    Temp_new = Temp + h_dl(gl%ncomp)
                    rholiq_new = rholiq + h_dl(gl%ncomp+1)

                    Var_X(j) = Temp_new
                    Var_X(gl%ncomp+1) = rholiq_new
                end if
                sum_liq = sum_liq + x_liq_new(j)
                ! check if the stepsize is too large.
                ! The 1st condition: x(j)_vap < 1
                ! The 2nd condition: x(j)_vap > 0
                ! The 3rd condition: Temp_new < 0
                ! The 4th condition: DT > 10 K
                ! have to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_liq_new(j) > 1.D0) .OR. (x_liq_new(j) <= 0.D0) &
                    & .OR. (Temp_new < 0.d0) .OR. (abs(h_dl(gl%NCOMP)) > 10.d0)) then
                    stepsize = 2.D0
                    h_dl = h_dl/stepsize
                    if (maxval(dabs(h_dl)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    press_new = press
                    j = 0
                    sum_liq = 0.D0
                end if
                j = j + 1
            end do
            sum_vap = 1.D0
            ! --------------------------------------------------------------------------------
        end select

        call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press_new, Temp_new, x_vap_new, x_liq_new, rhovap_new, rholiq_new, vapfrac, iFlash, GibbsEQN_new, errval)
        call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)
        call GibbsEQN_betrag(gl,eqn,GibbsEQN_new,GibbsEQN_b_new)

        call rho_rechner_DogLeg(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,h_dl,g,rho)

        if (rho > 0)then

            ! write the new values to the variables
            x_vap = x_vap_new!/sum_vap
            x_liq = x_liq_new!/sum_liq
            Temp = Temp_new
            press = press_new
            rhovap = rhovap_new
            rholiq = rholiq_new

            !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
            max_del = 0.D0

            if (iFlash < 3)then
                j = gl%ncomp+2
            else
                j = gl%ncomp+1
            end if

            Do k = 1, j
                if(abs(h_dl(k) / Var_X(k)) > max_del) then
                    max_del = abs(h_dl(k) / Var_X(k))
                end if
            end do

            if ((max_del) < eps_del) return

            !Catch unphysical temperatures and pressures
            !if ((temp < 0.D0) .or. temp > 1000.D0) then
            ! errval = -4321
            ! exit
            !end if
            !if ((press < 0.D0) .or. press > 1000.D0) then
            ! errval = -4322
            ! exit
            !end if
            call SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,press, Temp, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
            call GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)

            if(GibbsEQN_b <= tol) then
                if (press_new == 0.D0) then
                    gl%molfractions = x_liq
                    call reduced_parameters_calc(gl,Temp)
                    press = P_CALC(gl,Temp, rholiq, 0)
                    gl%molfractions = z
                    call reduced_parameters_calc(gl,Temp)
                end if
                return
            end if
        end if

        if (rho > 0.75) then
            call q_rechner(gl,delta, h_dl_b, q)
            delta = maxval(q)
        else if (rho < 0.25) then
            delta = delta / 2
        end if

    end do

    ! If convergence was not reached after 30 iterations --> algorithm failed!
    if ((i > 100) .AND. (maxval(dabs(GibbsEQN)) > 1.D-6)) then
        errval = -2222
    End if

    if ((i > 100) .AND. (GibbsEQN_b < 1.D-6)) then
        if (press_new == 0.D0) then
            gl%molfractions = x_liq
            call reduced_parameters_calc(gl,Temp)
            press = P_CALC(gl,Temp, rholiq, 0)
            gl%molfractions = z
            call reduced_parameters_calc(gl,Temp)
        end if
    end if

    gl%molfractions = z


    end subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg_mix
    !**************************************************************************







    end submodule impl
