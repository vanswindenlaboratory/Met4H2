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
    !Cite as: Span, R.; Beckmuller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.;
    !          Jager, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020):
    !          TREND. Thermodynamic Reference and Engineering Data 5.0.
    !          Lehrstuhl fur Thermodynamik, Ruhr-Universitat Bochum.

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

    ! module for file ptflash_solids_mix.f90
    submodule (ptflash_solids_mix_module) impl
    !global use inclusion
    use module_all_types
    use calc_functions
    use rhomix_pt_module
    use flash_module
    use waterice_module
    use dryice_module
    use hdrt_chem_pot_module
    use spline_module
    use setup_module
    use vle_derivs_module
    use utility_module
    use reduced_parameters_calc_module
    use gibbsderivs_module

    contains


    !NEW ALGORITHM TO CALCULATE THE EQUILIBRIUM OF FOUR PHASES FOR N COMPONENTS
    !**************************************************************************
    module subroutine ptflash_solid_NC_4P(gl,press, Temp, x_known, rho, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, errval, iter)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF FOUR PHASE EQUILIBRIA (QUADRUPLE POINTS, LINES, AND AREAS) FOR N COMPONENT MIXTURES
    ! ACCORDING TO THE PHASE RULE OF GIBBS:
    !   - FOR TWO COMPONENT MIXTURES (N=2) FOUR PHASE EQUILIBRIA ARE POINTS IN THE P-T-DIAGRAM
    !   - FOR THREE COMPONENT MIXTURES (N=3) FOUR PHASE EQUILIBRIA ARE LINES IN THE P-T-DIAGRAM
    !   - FOR MORE THAN THREE COMPONENTS (N>3) FOUR PHASE EQUILIBRIA ARE AREAS IN THE P-T-DIAGRAM
    !---------------------------------------------------------------------------------------------
    ! THE FOLLOWING CALCULATIONS MIGHT BE PERFORMED
    ! Case 1:   a) VLwHIw      Indicators: solidtype(1) = 1     -> Pure solid is solid water (Attention: set solid_pos to the position of water in the fluid vector)
    !                                      solidtype(2) = 1     -> Hydrate phase exists
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LxLwHIw     Indicators: solidtype(1) = 1     -> Pure solid is solid water
    !                                      solidtype(2) = 1     -> Hydrate phase exists
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 2:   a) VLwHIc      Indicators: solidtype(1) = 2     -> Pure solid is solid co2
    !                                      solidtype(2) = 1     -> Hydrate phase exists
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LxLwHIc     Indicators: solidtype(1) = 2     -> Pure solid is solid co2
    !                                      solidtype(2) = 1     -> Hydrate phase exists
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 3:   a) VLxLwH      Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 1     -> Hydrate forms
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LyLxLwH     Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 1     -> Hydrate forms
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 4:   a) VLxLwIw     Indicators: solidtype(1) = 1     -> Pure water ice forms
    !                                      solidtype(2) = 0     -> No Hydrate forms
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LyLxLwIw    Indicators: solidtype(1) = 1     -> Pure water ice forms
    !                                      solidtype(2) = 0     -> No Hydrate forms
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 5:   a) VLxLwIc     Indicators: solidtype(1) = 2     -> Dry ice forms
    !                                      solidtype(2) = 0     -> No Hydrate forms
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LyLxLwIc    Indicators: solidtype(1) = 2     -> Dry ice forms
    !                                      solidtype(2) = 0     -> No Hydrate forms
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 6:   a) VLxLyLz     Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 0     -> No Hydrate forms
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LwLxLyLz    Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 0     -> Hydrate forms
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    !MAYBE NEEDED IN THE FUTURE --> Quadruple points with different hydrate structures HsI and HsII
    ! Case 7:   a) VLwHsIHsII  Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 1     -> Hydrate forms --> ADDITIONAL INDICATOR REQUIRED, BECAUSE TWO HYDRATE STRUCTURES FORM!!!
    !                                      iphase = 2           -> "Lighter" fluid phase is vapor (take vapor density solution in density solver)
    !           b) LxLwHsIHsII Indicators: solidtype(1) = 0     -> No pure solid
    !                                      solidtype(2) = 1     -> Hydrate forms --> ADDITIONAL INDICATOR REQUIRED, BECAUSE TWO HYDRATE STRUCTURES FORM!!!
    !                                      iphase = 1           -> "Lighter" fluid phase is liquid (take liquid density solution in density solver)
    ! Case 8:   a) VIwHsIHsII  Indicators: solidtype(1) = 1     -> Solid water forms
    !                                      solidtype(2) = 1     -> Hydrate forms --> ADDITIONAL INDICATOR REQUIRED, BECAUSE TWO HYDRATE STRUCTURES FORM!!!
    !                                      iphase = 2           -> Fluid phase is vapor (take vapor density solution in density solver)
    !           b) LxIwHsIHsII Indicators: solidtype(1) = 1     -> Solid water forms
    !                                      solidtype(2) = 1     -> Hydrate forms --> ADDITIONAL INDICATOR REQUIRED, BECAUSE TWO HYDRATE STRUCTURES FORM!!!
    !                                      iphase = 1           -> Fluid phase is any liquid (take liquid density solution in density solver)
    !
    !--------------------------------------------------------------------------
    ! DEPENDING ON THE NUMBER OF COMPONENTS IN THE MIXTURE DIFFERENT CALCULATIONS CAN BE PERFORMED WHICH ARE INDICATED BY
    ! THE VARIABLE IFLASH. FOR THE CALCULATION OF PHASE BOUNDARIES, THE PHASE FRACTION (BETA) WHICH SHOULD BE 0 MUST BE GIVEN!!
    ! NOTE THAT THE BETAs ARE GIVEN IN THE ORDER AS WRITTEN ABOVE, E.G.:
    !
    ! QUADRUPLE POINT: VLwHIw  --> beta_V  = Phasefrac(1), beta_Lw = Phasefrac(2), beta_H  = Phasefrac(3), beta_Iw = Phasefrac(4)
    ! QUADRUPLE POINT: VLxLwH  --> beta_V  = Phasefrac(1), beta_Lx = Phasefrac(2), beta_lw = Phasefrac(3), beta_H  = Phasefrac(4)
    ! QUADRUPLE POINT: VLxHIc  --> beta_V  = Phasefrac(1), beta_Lx = Phasefrac(2), beta_H  = Phasefrac(3), beta_Ic = Phasefrac(4)
    ! QUADRUPLE POINT: LxLwHIc --> beta_Lx = Phasefrac(1), beta_Lw = Phasefrac(2), beta_H  = Phasefrac(3), beta_Ic = Phasefrac(4)
    !
    ! The variable Phasefrac_0 then indicates which of the betas should be 0
    !
    !
    !---------------------------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - FOUR PHASE LINE BOUNDARY        : p, x Given. one phase fraction = 0 (indicated by Phasefrac_0)             --  iFlash = 1
    !   - FOUR PHASE LINE BOUNDARY        : T, x Given. one phase fraction = 0 (indicated by Phasefrac_0)             --  iFlash = 2
    !   - FOUR PHASE FLASH: p, T and x VECTOR GIVEN                                                                   --  iFlash = 3
    !---------------------------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------------------------
    ! Variables:
    ! INPUT (and partly output):
    !   press           - Pressure [MPa]
    !   Temp            - Temperature [K]
    !   x_known         - Composition of the overall mixture
    !   rho             - Vector with the density of all phases
    !   x_sol           - Composition of the solid phase (0 except for the position solidnr)
    !   x_fluid1        - Composition of the first fluid phase
    !   x_fluid2        - Composition of the second fluid phase
    !   x_fluid3        - Composition of the third fluid phase
    !   x_fluid4        - Composition of the fourth fluid phase
    !   rhofluid1_est   - Estimated density of fluid phase 1
    !   rhofluid2_est   - Estimated density of fluid phase 2
    !   rhofluid3_est   - Estimated density of fluid phase 3
    !   rhofluid4_est   - Estimated density of fluid phase 4
    !   Phasefrac       - Vector with the phase fractions of the four phases
    !   iFlash          - Indicates if a four phase boundary at given pressure (iFlash = 1) or at given temperature (iFlash = 2) or a flash calculation at given T and p (iflash = 3) shall be performed
    !   iphase          - Control variable for the density solver:
    !                   - iphase = 1 --> Calculate liquid density
    !                   - iphase = 2 --> Calculate vapor density
    !                   - iphase = 0 --> Let density solver choose density
    !   Phasefrac_0     - Indicates which phasefraction should be 0, if a four phase boundary is calculated
    ! OUTPUT:
    !   errval          - Error value
    !   iter            - number of iterations

    !---------------------------------------------------------------------------
    ! A. J輍er Feb. 2016

    !**************************************************************************






    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x_known
    double precision, dimension(4):: rho
    double precision, dimension(30):: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_hyd, x_sol
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    double precision, dimension(4) :: Phasefrac
    integer:: errval
    integer:: iter
    integer:: iFlash
    integer:: iphase
    integer:: Phasefrac_0

    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(60):: GibbsEQN, Delta_X
    double precision:: stepsize, Temp_new, press_new
    double precision, dimension(30):: x_fluid1_new, x_fluid2_new, x_fluid3_new, x_fluid4_new
    double precision, dimension(4):: Phasefrac_new
    integer:: i, j, k, eqn
    !integer, dimension(1):: maxID
    character(255):: herr
    double precision:: sum_fluid1, sum_fluid2, sum_fluid3, sum_fluid4, sum_test

    !Variables for the break criterion
    double precision:: eps_Gibbs, eps_del, max_del, eps_Gibbs_min
    double precision, dimension(60):: Var_X

    !Variables which store how many fluid phases are present
    integer:: nr_of_fluidphases

    !Mapping variable for phasefraction.
    integer, dimension(3):: mapp_index     !This variable helps figuring out which phasefractions need to be updated. The index indicated by Phasefrac_0 disappears in the variable
    !Example: If Phasefrac_0 = 2 --> mapp_index(1) = 1, mapp_index(2) = 3, mapp_index(3) = 4
    logical:: unphysical_x      !Variable that indicates if physically unreasonable values for the mole fractions occur during iteration

    double precision, dimension(30):: fug_guest     !Fugacity of all guest components in the hydrate phase


    errval = 0
    Delta_X = 1.D0

    iter = 1
    GibbsEQN = 1.D10
    unphysical_x = .false.

    !Check input values for plausiblity
    if ((Temp .lt. 1.D-16) .or. (press .lt. 1.D-16)) then
        errval = -15571
        return
    end if
    if ((iphase .lt. 0) .or. (iphase .gt. 2)) then
        errval = -15571
        return
    end if
    if ((iFlash .lt. 1) .or. (iFlash .gt. 3)) then
        errval = -15571
        return
    end if
    !For more than 3 components in the mixture, the overall compositon is always needed in order to calculate the four-phase equilibrium
    if (gl%ncomp .gt. 3) then
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_known(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_known(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
    end if

    !If the maximum residuum of the system of equations is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del and the maximum value of Gibbs_EQN is smaller than eps_Gibbs_min, the calculation is finished
    eps_del = 1.d-4
    eps_Gibbs_min = 1.D-4
    if ((gl%solidtype(1) == 2) .and. (gl%solidtype(2) == 1)) eps_Gibbs_min = 1.d-2

    !FOR BINARY MIXTURES, USE DIFFERENT CONVERGENCE CRITERIA IN ORDER TO GET THE SAME RESULTS AS BEFORE!!
    !THIS SHOULD BE DELETED AFTER ALL TESTS ARE FINISHED!!
    If (gl%ncomp .eq. 2) then
        !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
        eps_Gibbs = 1.d-8
        !If the relative change of the unknowns (T,p or x) is below eps_del, the calculation is finished
        eps_del = 1.d-12
        eps_Gibbs_min = 1000.D0
    end if


    if ((iFlash .eq. 1) .or. (iFlash .eq. 2)) then
        if ((Phasefrac_0 < 1) .or. (Phasefrac_0 > 4)) then
            errval = -15571
            return
        end if
        Select case (Phasefrac_0)
        Case(1) !Phase fraction of phase 1 is 0
            mapp_index(1) = 2
            mapp_index(2) = 3
            mapp_index(3) = 4
        Case(2) !Phase fraction of phase 2 is 0
            mapp_index(1) = 1
            mapp_index(2) = 3
            mapp_index(3) = 4
        Case(3) !Phase fraction of phase 3 is 0
            mapp_index(1) = 1
            mapp_index(2) = 2
            mapp_index(3) = 4
        Case(4) !Phase fraction of phase 4 is 0
            mapp_index(1) = 1
            mapp_index(2) = 2
            mapp_index(3) = 3
        end select
    end if


    eqn = 0
    !For the different quadruple points a different amount of equations needs to be solved
    !---------------------------------------------------------------------------
    !Case: VLwHIw, LxLwHIw, VLwHIc, LxLwHIc
    if ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then
        !Check if the start values for the phase compositions are ok
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid1(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid1(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid2(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid2(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        nr_of_fluidphases = 2
        if (gl%ncomp .eq. 2) then  !binary mixture
            eqn = 4
        else
            eqn = 2*gl%ncomp + 2
        end if
    end if
    !Case: VLxLwH, LyLxLwH, VLxLwIw, LyLxLwIw, VLxLwIc, LyLxLwIc
    if ( ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1)) &
        & .or. ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 0))) then
        !Check if the start values for the phase compositions are ok
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid1(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid1(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid2(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid2(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid3(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid3(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        nr_of_fluidphases = 3
        if (gl%ncomp .eq. 2) then  !binary mixture
            eqn = 5
        else
            eqn = 3*gl%ncomp + 1
        end if
    end if
    !Case: VLxLyLz, LwLxLyLz
    if ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 0)) then
        !Check if the start values for the phase compositions are ok
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid1(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid1(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid2(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid2(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid3(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid3(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        sum_test = 0.D0
        Do k = 1, gl%ncomp
            if (x_fluid4(k) < 0.D0) then
                errval = -15571
                return
            end if
            sum_test = sum_test + x_fluid4(k)
        end do
        if (dabs(sum_test - 1.D0) > 1.D-14) then
            errval = -15571
            return
        end if
        !---
        nr_of_fluidphases = 4
        if (gl%ncomp .eq. 2) then  !binary mixture
            !Impossible for binary mixture
            errval = -15571
            return
        else
            eqn = 4*gl%ncomp
        end if
    end if
    !---------------------------------------------------------------------------
    if (eqn .eq. 0) then
        errval = -15571
        return
    end if

    do i = 1, 200

        call SysOfEqs_solid_NC_4P(gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
            & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iphase, nr_of_fluidphases, &
            & fug_guest, rho, GibbsEQN, errval)

        ! Break Condition
        ! This is the exit condition for the VLE iteration!
        if (maxval(dabs(GibbsEQN)) .lt. eps_Gibbs) return
        if (errval /= 0) return

        call Jacobi_solid_NC_4P (gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, &
            & Phasefrac, iFlash, nr_of_fluidphases, fug_guest, mapp_index, JacMatrix, errval)

        Delta_X = - GibbsEQn

        !Solve the system of equations
        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
        if (errval /= 0) then
            errval = -2222
            return
        endif
        x_fluid1_new = 0.D0
        x_fluid2_new = 0.D0
        x_fluid3_new = 0.D0
        x_fluid4_new = 0.D0
        Temp_new = 0.D0
        press_new = 0.D0
        stepsize = 0.D0
        iter = i
        Phasefrac_new = 0.D0

        !Update the variables
        !ONLY FOR BINARY MIXTURES
        !********************************************************************************************************************************
        If (gl%ncomp == 2) then
            ! --------------------------------------------------------------------------------
            if ( ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1)) &
                & .or. ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 0))) then !VLxLwH, LyLxLwH, VLxLwIw, LyLxLwIw, VLxLwIc, LyLxLwIc
                j = 1
                Do while(j < 2)
                    x_fluid1_new(1) = x_fluid1(1) + Delta_x(1)
                    x_fluid1_new(2) = 1.D0 - x_fluid1_new(1)
                    x_fluid2_new(1) = x_fluid2(1) + Delta_x(2)
                    x_fluid2_new(2) = 1.D0 - x_fluid2_new(1)
                    x_fluid3_new(1) = x_fluid3(1) + Delta_x(3)
                    x_fluid3_new(2) = 1.D0 - x_fluid3_new(1)
                    temp_new = temp + Delta_X(4)
                    press_new = press + Delta_X(5) /1.D6
                    Var_X(1) = x_fluid1_new(1)
                    Var_X(2) = x_fluid2_new(1)
                    Var_X(3) = x_fluid3_new(1)
                    Var_X(4) = temp_new
                    Var_X(5) = press_new
                    ! check if the stepsize is too large.
                    ! The condition: x(j)_vap AND x(j)_liq > 0
                    ! hase to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    if ((x_fluid1_new(1) .lt. 0.D0) .or. (x_fluid1_new(1) .gt. 1.D0) .or. &
                        & (x_fluid2_new(1) .lt. 0.D0) .or. (x_fluid2_new(1) .gt. 1.D0) .or. &
                        & (x_fluid3_new(1) .lt. 0.D0) .or. (x_fluid3_new(1) .gt. 1.D0) .or. &
                        & (Temp_new < 0.D0) .or. (press_new < 0.D0) ) then
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) .lt. 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                    end if
                    j = j + 1
                end do
                ! --------------------------------------------------------------------------------
            elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then  !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc
                j = 1
                Do while(j < 2)
                    x_fluid1_new(1) = x_fluid1(1) + Delta_x(1)
                    x_fluid1_new(2) = 1.D0 - x_fluid1_new(1)
                    x_fluid2_new(1) = x_fluid2(1) + Delta_x(2)
                    x_fluid2_new(2) = 1.D0 - x_fluid2_new(1)
                    temp_new = temp + Delta_X(3)
                    press_new = press + Delta_X(4) /1.D6
                    Var_X(1) = x_fluid1_new(1)
                    Var_X(2) = x_fluid2_new(1)
                    Var_X(3) = temp_new
                    Var_X(4) = press_new
                    ! check if the stepsize is too large.
                    ! The condition: x(j)_vap AND x(j)_liq > 0
                    ! hase to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    if ((x_fluid1_new(1) .lt. 0.D0) .or. (x_fluid1_new(1) .gt. 1.D0) .or. &
                        & (x_fluid2_new(1) .lt. 0.D0) .or. (x_fluid2_new(1) .gt. 1.D0) .or. &
                        & (Temp_new < 0.D0) .or. (press_new < 0.D0) ) then
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) .lt. 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                    end if
                    j = j + 1
                end do
                ! --------------------------------------------------------------------------------
            end if
        End if
        !********************************************************************************************************************************

        !********************************************************************************************************************************
        !UPDATE OF VARIABLES FOR MORE THAN 2 COMPONENTS
        If (gl%ncomp .gt. 2) then
            press_new = press
            temp_new = temp
            ! --------------------------------------------------------------------------------
            if ( ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1)) &
                & .or. ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 0))) then !VLxLwH, LyLxLwH, VLxLwIw, LyLxLwIw, VLxLwIc, LyLxLwIc
                j = 1
                Do while(j < 2)
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    sum_fluid3 = 0.D0
                    sum_fluid4 = 0.D0
                    Do k=1, gl%ncomp-1
                        x_fluid1_new(k) = x_fluid1(k) + Delta_x(k)
                        Var_X(k) = x_fluid1_new(k)
                        sum_fluid1 = sum_fluid1 + x_fluid1_new(k)

                        x_fluid2_new(k) = x_fluid2(k) + Delta_x(gl%ncomp-1+k)
                        Var_X(gl%ncomp-1+k) = x_fluid2_new(k)
                        sum_fluid2 = sum_fluid2 + x_fluid2_new(k)

                        x_fluid3_new(k) = x_fluid3(k) + Delta_x(2*gl%ncomp-2+k)
                        Var_X(2*gl%ncomp-2+k) = x_fluid3_new(k)
                        sum_fluid3 = sum_fluid3 + x_fluid3_new(k)
                    end do
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                    x_fluid3_new(gl%ncomp) = 1.D0 - sum_fluid3
                    if (iFlash == 1) then
                        temp_new = temp + Delta_X(3*gl%ncomp-3+1)
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(3*gl%ncomp-3+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(3*gl%ncomp-3+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(3*gl%ncomp-3+4)
                        Var_X(3*gl%ncomp-3+1) = temp_new
                        Var_X(3*gl%ncomp-3+2) = Phasefrac_new(mapp_index(1))
                        Var_X(3*gl%ncomp-3+3) = Phasefrac_new(mapp_index(2))
                        Var_X(3*gl%ncomp-3+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 2) then
                        press_new = press + Delta_X(3*gl%ncomp-3+1) /1.D6
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(3*gl%ncomp-3+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(3*gl%ncomp-3+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(3*gl%ncomp-3+4)
                        Var_X(3*gl%ncomp-3+1) = press_new
                        Var_X(3*gl%ncomp-3+2) = Phasefrac_new(mapp_index(1))
                        Var_X(3*gl%ncomp-3+3) = Phasefrac_new(mapp_index(2))
                        Var_X(3*gl%ncomp-3+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 3) then
                        Phasefrac_new(1) = Phasefrac(1) + Delta_X(3*gl%ncomp-3+1)
                        Phasefrac_new(2) = Phasefrac(2) + Delta_X(3*gl%ncomp-3+2)
                        Phasefrac_new(3) = Phasefrac(3) + Delta_X(3*gl%ncomp-3+3)
                        Phasefrac_new(4) = Phasefrac(4) + Delta_X(3*gl%ncomp-3+4)
                        Var_X(3*gl%ncomp-3+1) = Phasefrac_new(1)
                        Var_X(3*gl%ncomp-3+2) = Phasefrac_new(2)
                        Var_X(3*gl%ncomp-3+3) = Phasefrac_new(3)
                        Var_X(3*gl%ncomp-3+4) = Phasefrac_new(4)
                    end if


                    ! check if the stepsize is too large.
                    ! The condition: x(j)_vap AND x(j)_liq > 0
                    ! hase to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    Do k=1, gl%ncomp
                        if ((x_fluid1_new(k) .lt. 0.D0) .or. (x_fluid1_new(k) .gt. 1.D0) .or. &
                            & (x_fluid2_new(k) .lt. 0.D0) .or. (x_fluid2_new(k) .gt. 1.D0)  .or. &
                            & (x_fluid3_new(k) .lt. 0.D0) .or. (x_fluid3_new(k) .gt. 1.D0)) then

                            unphysical_x = .true. !One of the mole factions is out of physically reasonable bounds

                        end if
                    end do
                    if ((temp_new .lt. 0.D0) .or. (press_new .lt. 0.D0) .or. unphysical_x) then
                        unphysical_x = .false.
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) .lt. 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                    end if
                    j = j + 1
                end do
                ! --------------------------------------------------------------------------------
            elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then  !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc
                j = 1
                Do while(j < 2)
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    sum_fluid3 = 0.D0
                    sum_fluid4 = 0.D0
                    Do k=1, gl%ncomp-1
                        x_fluid1_new(k) = x_fluid1(k) + Delta_x(k)
                        Var_X(k) = x_fluid1_new(k)
                        sum_fluid1 = sum_fluid1 + x_fluid1_new(k)

                        x_fluid2_new(k) = x_fluid2(k) + Delta_x(gl%ncomp-1+k)
                        Var_X(gl%ncomp-1+k) = x_fluid2_new(k)
                        sum_fluid2 = sum_fluid2 + x_fluid2_new(k)
                    end do
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                    if (iFlash == 1) then
                        temp_new = temp + Delta_X(2*gl%ncomp-2+1)
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(2*gl%ncomp-2+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(2*gl%ncomp-2+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(2*gl%ncomp-2+4)
                        Var_X(2*gl%ncomp-2+1) = temp_new
                        Var_X(2*gl%ncomp-2+2) = Phasefrac_new(mapp_index(1))
                        Var_X(2*gl%ncomp-2+3) = Phasefrac_new(mapp_index(2))
                        Var_X(2*gl%ncomp-2+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 2) then
                        press_new = press + Delta_X(2*gl%ncomp-2+1) /1.D6
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(2*gl%ncomp-2+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(2*gl%ncomp-2+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(2*gl%ncomp-2+4)
                        Var_X(2*gl%ncomp-2+1) = press_new
                        Var_X(2*gl%ncomp-2+2) = Phasefrac_new(mapp_index(1))
                        Var_X(2*gl%ncomp-2+3) = Phasefrac_new(mapp_index(2))
                        Var_X(2*gl%ncomp-2+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 3) then
                        Phasefrac_new(1) = Phasefrac(1) + Delta_X(2*gl%ncomp-2+1)
                        Phasefrac_new(2) = Phasefrac(2) + Delta_X(2*gl%ncomp-2+2)
                        Phasefrac_new(3) = Phasefrac(3) + Delta_X(2*gl%ncomp-2+3)
                        Phasefrac_new(4) = Phasefrac(4) + Delta_X(2*gl%ncomp-2+4)
                        Var_X(2*gl%ncomp-2+1) = Phasefrac_new(1)
                        Var_X(2*gl%ncomp-2+2) = Phasefrac_new(2)
                        Var_X(2*gl%ncomp-2+3) = Phasefrac_new(3)
                        Var_X(2*gl%ncomp-2+4) = Phasefrac_new(4)
                    end if

                    ! check if the stepsize is too large.
                    ! The condition: x(j)_vap AND x(j)_liq > 0
                    ! hase to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    Do k=1, gl%ncomp
                        if ((x_fluid1_new(k) .lt. 0.D0) .or. (x_fluid1_new(k) .gt. 1.D0) .or. &
                            & (x_fluid2_new(k) .lt. 0.D0) .or. (x_fluid2_new(k) .gt. 1.D0)) then

                            unphysical_x = .true. !One of the mole factions is out of physically reasonable bounds

                        end if
                    end do
                    if ((temp_new .lt. 0.D0) .or. (press_new .lt. 0.D0) .or. unphysical_x) then
                        unphysical_x = .false.
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) .lt. 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                    end if
                    j = j + 1
                end do
                ! --------------------------------------------------------------------------------
            elseif ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 0)) then    !Case: VLxLyLz, LwLxLyLz

                j = 1
                Do while(j < 2)
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    sum_fluid3 = 0.D0
                    sum_fluid4 = 0.D0
                    Do k=1, gl%ncomp-1
                        x_fluid1_new(k) = x_fluid1(k) + Delta_x(k)
                        Var_X(k) = x_fluid1_new(k)
                        sum_fluid1 = sum_fluid1 + x_fluid1_new(k)

                        x_fluid2_new(k) = x_fluid2(k) + Delta_x(gl%ncomp-1+k)
                        Var_X(gl%ncomp-1+k) = x_fluid2_new(k)
                        sum_fluid2 = sum_fluid2 + x_fluid2_new(k)

                        x_fluid3_new(k) = x_fluid3(k) + Delta_x(2*gl%ncomp-2+k)
                        Var_X(2*gl%ncomp-2+k) = x_fluid3_new(k)
                        sum_fluid3 = sum_fluid3 + x_fluid3_new(k)

                        x_fluid4_new(k) = x_fluid4(k) + Delta_x(3*gl%ncomp-3+k)
                        Var_X(3*gl%ncomp-3+k) = x_fluid4_new(k)
                        sum_fluid4 = sum_fluid4 + x_fluid4_new(k)
                    end do
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                    x_fluid3_new(gl%ncomp) = 1.D0 - sum_fluid3
                    x_fluid4_new(gl%ncomp) = 1.D0 - sum_fluid4
                    if (iFlash == 1) then
                        temp_new = temp + Delta_X(4*gl%ncomp-4+1)
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(4*gl%ncomp-4+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(4*gl%ncomp-4+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(4*gl%ncomp-4+4)
                        Var_X(4*gl%ncomp-4+1) = temp_new
                        Var_X(4*gl%ncomp-4+2) = Phasefrac_new(mapp_index(1))
                        Var_X(4*gl%ncomp-4+3) = Phasefrac_new(mapp_index(2))
                        Var_X(4*gl%ncomp-4+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 2) then
                        press_new = press + Delta_X(4*gl%ncomp-4+1) /1.D6
                        Phasefrac_new(mapp_index(1)) = Phasefrac(mapp_index(1)) + Delta_X(4*gl%ncomp-4+2)
                        Phasefrac_new(mapp_index(2)) = Phasefrac(mapp_index(2)) + Delta_X(4*gl%ncomp-4+3)
                        Phasefrac_new(mapp_index(3)) = Phasefrac(mapp_index(3)) + Delta_X(4*gl%ncomp-4+4)
                        Var_X(4*gl%ncomp-4+1) = press_new
                        Var_X(4*gl%ncomp-4+2) = Phasefrac_new(mapp_index(1))
                        Var_X(4*gl%ncomp-4+3) = Phasefrac_new(mapp_index(2))
                        Var_X(4*gl%ncomp-4+4) = Phasefrac_new(mapp_index(3))
                    end if
                    if (iFlash == 3) then
                        Phasefrac_new(1) = Phasefrac(1) + Delta_X(4*gl%ncomp-4+1)
                        Phasefrac_new(2) = Phasefrac(2) + Delta_X(4*gl%ncomp-4+2)
                        Phasefrac_new(3) = Phasefrac(3) + Delta_X(4*gl%ncomp-4+3)
                        Phasefrac_new(4) = Phasefrac(4) + Delta_X(4*gl%ncomp-4+4)
                        Var_X(4*gl%ncomp-4+1) = Phasefrac_new(1)
                        Var_X(4*gl%ncomp-4+2) = Phasefrac_new(2)
                        Var_X(4*gl%ncomp-4+3) = Phasefrac_new(3)
                        Var_X(4*gl%ncomp-4+4) = Phasefrac_new(4)
                    end if


                    ! check if the stepsize is too large.
                    ! The condition: x(j)_vap AND x(j)_liq > 0
                    ! hase to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    Do k=1, gl%ncomp
                        if ((x_fluid1_new(k) .lt. 0.D0) .or. (x_fluid1_new(k) .gt. 1.D0) .or. &
                            & (x_fluid2_new(k) .lt. 0.D0) .or. (x_fluid2_new(k) .gt. 1.D0) .or. &
                            & (x_fluid3_new(k) .lt. 0.D0) .or. (x_fluid3_new(k) .gt. 1.D0) .or. &
                            & (x_fluid4_new(k) .lt. 0.D0) .or. (x_fluid4_new(k) .gt. 1.D0)) then

                            unphysical_x = .true. !One of the mole factions is out of physically reasonable bounds

                        end if
                    end do
                    if ((temp_new .lt. 0.D0) .or. (press_new .lt. 0.D0) .or. unphysical_x) then
                        unphysical_x = .false.
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) .lt. 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                    end if
                    j = j + 1
                end do

                ! --------------------------------------------------------------------------------
            end if
        end if
        !********************************************************************************************************************************


        ! write the new values to the variables
        x_fluid1 = x_fluid1_new
        x_fluid2 = x_fluid2_new
        x_fluid3 = x_fluid3_new
        x_fluid4 = x_fluid4_new
        Temp = Temp_new
        press = press_new
        Phasefrac = Phasefrac_new

        !Second exit criterion: If the maximum relative change of the variables is lower than eps_del
        max_del = 0.D0
        Do k = 1, j - 1
            if(abs(delta_X(k) / Var_X(k)) .gt. max_del) then
                max_del = abs(delta_X(k) / Var_X(k))
            end if
        end do

        if (((max_del) .lt. eps_del) .and. (maxval(dabs(GibbsEQN)) .lt. eps_Gibbs_min)) then

            return

        end if


    end do

    ! Iteration failed!
    if (i > 200) errval = -2222

    end subroutine ptflash_solid_NC_4P



    !**************************************************************************
    module subroutine SysOfEqs_solid_NC_4P (gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iphase_given, nr_of_fluidphases, &
        & fug_guest, rho, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS FOR FOUR PHASES IN EQUILIBRIUM.
    ! FOR VARIABLE DEFINITIONS AND OTHER EXPLANATIONS, SEE ROUTINE ptflash_sol_NC_4P
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !A. J輍er, February 2016








    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4
    double precision, dimension(30):: x_sol, x_hyd
    double precision, dimension(60):: GibbsEQN
    double precision, dimension(30):: fug_guest
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    double precision, dimension(4):: Phasefrac, rho
    integer:: errval, iphase_given, nr_of_fluidphases

    double precision, dimension(30):: ChemPot_fluid1, ChemPot_fluid2, ChemPot_fluid3, ChemPot_fluid4
    double precision:: ChemPot_sol
    double precision:: rhoredmix_orig, tredmix_orig
    double precision:: Rmix
    double precision, dimension(30):: lnf!, occup_ls, occup_ld, occup_sms, occup_smd

    integer:: k

    !Variables needed for hydrate model
    double precision:: chpw, Chempot_hyd
    double precision, dimension(3,30):: occup, CiJ, occup_single, occup_double
    double precision:: rho_H

    integer:: iphase

    errval = 0
    GibbsEQN = 0.D0

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    !Fluid Phase 1
    !**********************************************************
    gl%molfractions = x_fluid1
    call reduced_parameters_calc(gl,Temp)
    ! calculate the density of fluid phase 1
    gl%rho_vap = rhomix_calc(gl,Temp, Press, rhofluid1_est, iphase_given, 0)

    if (gl%rho_vap < 1.D-12) then
        errval = -8888
        return
    end if
    rho(1) = gl%rho_vap
    !----------------------------------
    ! Calculate the chemical potentials of fluid phase 1
    call R_mix_calc(gl,Rmix)
    call dna_dni(gl,Temp, gl%rho_vap, ChemPot_fluid1, 0)
    ChemPot_fluid1 = ChemPot_fluid1 * Rmix * Temp    !Chemical Potentials of the first fluid phase [J / mol]
    !Get the fugacities of the first fluid phase (THIS IS NEEDED FOR THE HYDRATE MODEL)
    if (gl%solidtype(2) .eq. 1) then
        call lnf_mix(gl,Temp, gl%rho_vap, press, lnf)
    end if
    !----------------------------------
    !**********************************************************

    if (nr_of_fluidphases > 1) then
        !Fluid Phase 2
        !**********************************************************
        iPhase = 1
        gl%molfractions = x_fluid2
        call reduced_parameters_calc(gl,Temp)
        ! calculate the density of phase 2
        gl%rho_liq1 = rhomix_calc(gl,Temp, press, rhofluid2_est, iPhase,0)

        if (gl%rho_liq1 < 1.D-12) then
            errval = -8888
            return
        end if
        rho(2) = gl%rho_liq1
        ! Calculate the chemical potentials of the second fluid phase
        call R_mix_calc(gl,Rmix)
        call dna_dni(gl,Temp, gl%rho_liq1, ChemPot_fluid2, 0)
        ChemPot_fluid2 = ChemPot_fluid2 * Rmix * Temp    !Chemical Potentials of fluid phase 2 [J / mol]
    End if

    if (nr_of_fluidphases > 2) then
        !Fluid Phase 3
        !**********************************************************
        iPhase = 1
        gl%molfractions = x_fluid3
        call reduced_parameters_calc(gl,Temp)
        ! calculate the gas phase density
        gl%rho_liq2 = rhomix_calc(gl,Temp, press, rhofluid3_est, iPhase,0)

        if (gl%rho_liq2 < 1.D-12) then
            errval = -8888
            return
        end if
        rho(3) = gl%rho_liq2
        !----------------------------------
        ! Calculate the chemical potentials of the (second) liquid phase
        call R_mix_calc(gl,Rmix)
        call dna_dni(gl,Temp, gl%rho_liq2, ChemPot_fluid3, 0)
        ChemPot_fluid3 = ChemPot_fluid3 * Rmix * Temp    !Chemical Potentials of the (second) liquid phase [J / mol]
        !----------------------------------
    End if

    if (nr_of_fluidphases > 3) then
        !Fluid Phase 4
        !**********************************************************
        iPhase = 1
        gl%molfractions = x_fluid4
        call reduced_parameters_calc(gl,Temp)
        ! calculate the gas phase density
        gl%rho_liq3 = rhomix_calc(gl,Temp, press, rhofluid4_est, iPhase,0)

        if (gl%rho_liq3 < 1.D-12) then
            errval = -8888
            return
        end if
        rho(4) = gl%rho_liq3
        !----------------------------------
        ! Calculate the chemical potentials of the (second) liquid phase
        call R_mix_calc(gl,Rmix)
        call dna_dni(gl,Temp, gl%rho_liq2, ChemPot_fluid4, 0)
        ChemPot_fluid4 = ChemPot_fluid4 * Rmix * Temp    !Chemical Potentials of the (second) liquid phase [J / mol]
        !----------------------------------
    End if

    if (gl%solidtype(2) .eq. 1) then
        !Hydrate phase
        !**********************************************************
        !Write the fugacities of the hydrates formers on positon 1 to nrofhydrateformers-1 of the vector fug_guest
        fug_guest(1:gl%nrofhydrateformers-1) = dexp(lnf(2:gl%nrofhydrateformers))*1.d6
        !Get the chemical potential of water in hydrate
        call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_guest,chpw)
        ChemPot_hyd = chpw              !Chemical Potential of Water in Hydrate [J / mol]
        !Get the density of the hydrate phase
        call hdrt_density(gl,temp,press*1.D6,fug_guest,occup,CiJ,rho_H)
        if (gl%solidtype(1) .gt. 0) then
            rho(3) = rho_H
        else
            rho(4) = rho_H
        end if
        !Get the composition of the hydrate phase
        call hdrt_mole_fract(gl,Temp,press*1.D6,fug_guest,occup,CiJ,x_hyd, occup_single, occup_double)
    end if

    if (gl%solidtype(1) .gt. 0) then
        !Pure solid phase
        !**********************************************************
        select case (gl%solidtype(1))
        case(1)
            ChemPot_sol = g_WaterIce(gl,Temp, press)    !Chemical Potential of Water Ice [J / mol]
            rho(4) = 1.D0 / v_WaterIce(gl,Temp,press)
        case(2)
            ChemPot_sol = g_DryIce(gl,Temp, press)     !Chemical Potential of Dry Ice [J / mol]
            rho(4) = 1.D0 / v_DryIce(gl,Temp,press)
        end select
    end if


    !----------------------------------
    ! Setting up the system of equations for the minimization of the Gibbs free energy.

    !SYSTEM OF EQUATIONS FOR BINARY MIXTURES
    !********************************************************************************************************************************
    If (gl%ncomp .eq. 2) then
        if ( (gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1) ) then !VLxLwH, LyLxLwH
            GibbsEQN(1) = ChemPot_fluid1(1) - ChemPot_fluid2(1)
            GibbsEQN(2) = ChemPot_fluid1(2) - ChemPot_fluid2(2)
            GibbsEQN(3) = ChemPot_fluid1(1) - ChemPot_fluid3(1)
            GibbsEQN(4) = ChemPot_fluid1(2) - ChemPot_fluid3(2)
            GibbsEQN(5) = ChemPot_fluid1(1) - ChemPot_hyd          !ChemPot of water is needed --> water always on position 1
        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then  !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc
            GibbsEQN(1) = ChemPot_fluid1(1) - ChemPot_fluid2(1)
            GibbsEQN(2) = ChemPot_fluid1(2) - ChemPot_fluid2(2)
            GibbsEQN(3) = ChemPot_fluid1(1) - ChemPot_hyd          !ChemPot of water is needed --> water always on position 1
            GibbsEQN(4) = ChemPot_fluid1(gl%solid_pos) - ChemPot_sol  !ChemPot_Vap(Pos_CO2) - ChemPot_sol   4
        end if
    end if
    !********************************************************************************************************************************

    !SYSTEM OF EQUATIONS FOR MIXURES WITH MORE THAN 2 COMPONENTS
    !********************************************************************************************************************************
    If (gl%ncomp .gt. 2) then
        if ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1)) then         !VLxLwH, LyLxLwH

            Do k=1, gl%ncomp
                GibbsEQN(k) = ChemPot_fluid1(k) - ChemPot_fluid2(k)
                GibbsEQN(gl%ncomp + k) = ChemPot_fluid1(k) - ChemPot_fluid3(k)
            end do
            GibbsEQN(2*gl%ncomp + 1) =  ChemPot_fluid1(1) - ChemPot_hyd    !Equality of chemical potential of water in fluid phase one and of water in hydrate phase
            Do k=1, gl%ncomp-1
                GibbsEQN(2*gl%ncomp + 1 + k) = Phasefrac(1) * x_fluid1(k) + Phasefrac(2) * x_fluid2(k) &
                    &     + Phasefrac(3) * x_fluid3(k) + Phasefrac(4) * x_hyd(k) - x_known(k)
            end do
            GibbsEQN(3*gl%ncomp + 1) = Phasefrac(1) + Phasefrac(2) + Phasefrac(3) + Phasefrac(4) - 1.D0


        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 0)) then     !VLxLwIw, LyLxLwIw, VLxLwIc, LyLxLwIc

            Do k=1, gl%ncomp
                GibbsEQN(k) = ChemPot_fluid1(k) - ChemPot_fluid2(k)
                GibbsEQN(gl%ncomp + k) = ChemPot_fluid1(k) - ChemPot_fluid3(k)
            end do
            GibbsEQN(2*gl%ncomp + 1) =  ChemPot_fluid1(gl%solid_pos) - ChemPot_Sol    !Equality of chemical potential of pure solid (water or co2) and chemical potential of the solid substance in fluid phase 1
            Do k=1, gl%ncomp-1
                GibbsEQN(2*gl%ncomp + 1 + k) = Phasefrac(1) * x_fluid1(k) + Phasefrac(2) * x_fluid2(k) &
                    &     + Phasefrac(3) * x_fluid3(k) + Phasefrac(4) * x_sol(k) - x_known(k)
            end do
            GibbsEQN(3*gl%ncomp + 1) = Phasefrac(1) + Phasefrac(2) + Phasefrac(3) + Phasefrac(4) - 1.D0


        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then     !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc

            Do k=1, gl%ncomp
                GibbsEQN(k) = ChemPot_fluid1(k) - ChemPot_fluid2(k)
            end do
            GibbsEQN(gl%ncomp + 1) = ChemPot_fluid1(1) - ChemPot_hyd    !Equality of chemical potential of water in fluid phase one and of water in hydrate phase
            GibbsEQN(gl%ncomp + 2) = ChemPot_fluid1(gl%solid_pos) - ChemPot_Sol    !Equality of chemical potential of pure solid (water or co2) and chemical potential of the solid substance in fluid phase 1

            Do k=1, gl%ncomp-1
                GibbsEQN(gl%ncomp + 2 + k) = Phasefrac(1) * x_fluid1(k) + Phasefrac(2) * x_fluid2(k) &
                    &     + Phasefrac(3) * x_hyd(k) + Phasefrac(4) * x_sol(k) - x_known(k)
            end do
            GibbsEQN(2*gl%ncomp + 2) = Phasefrac(1) + Phasefrac(2) + Phasefrac(3) + Phasefrac(4) - 1.D0


        elseif ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 0)) then     !VLxLyLz, LwLxLyLz

            Do k=1, gl%ncomp
                GibbsEQN(k) = ChemPot_fluid1(k) - ChemPot_fluid2(k)
                GibbsEQN(gl%ncomp + k) = ChemPot_fluid1(k) - ChemPot_fluid3(k)
                GibbsEQN(2*gl%ncomp + k) = ChemPot_fluid1(k) - ChemPot_fluid4(k)
            end do

            Do k=1, gl%ncomp-1
                GibbsEQN(3*gl%ncomp + k) = Phasefrac(1) * x_fluid1(k) + Phasefrac(2) * x_fluid2(k) &
                    &     + Phasefrac(3) * x_fluid3(k) + Phasefrac(4) * x_fluid4(k) - x_known(k)
            end do
            GibbsEQN(4*gl%ncomp) = Phasefrac(1) + Phasefrac(2) + Phasefrac(3) + Phasefrac(4) - 1.D0

        end if
    end if
    !********************************************************************************************************************************

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = x_known

    end subroutine SysOfEqs_solid_NC_4P
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_solid_NC_4P (gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, &
        & Phasefrac, iFlash, nr_of_fluidphases, fug_guest, mapp_index, JacMatrix, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS FOR FOUR PHASES IN EQUILIBRIUM.
    ! FOR VARIABLE DEFINITIONS AND OTHER EXPLANATIONS, SEE ROUTINE ptflash_sol_NC_4P
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !A. J輍er, February 2016







    implicit none

    type(type_gl) :: gl


    double precision:: Temp, press
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd
    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(4):: Phasefrac
    integer:: iFlash, nr_of_fluidphases
    integer, dimension(3):: mapp_index
    integer:: errval

    double precision:: rho_fluid, dChempot_dT_sol, dChempot_dp_sol
    double precision, dimension(30, 30):: dChempoti_dxj_fluid1, dChempoti_dxj_fluid2, dChempoti_dxj_fluid3, dChempoti_dxj_fluid4
    double precision, dimension(30):: dChempoti_dT_fluid1, dChempoti_dp_fluid1, dChempoti_dT_fluid2, dChempoti_dp_fluid2
    double precision, dimension(30):: dChempoti_dT_fluid3, dChempoti_dp_fluid3, dChempoti_dT_fluid4, dChempoti_dp_fluid4
    double precision, dimension(30):: d2nadnidT, dnadni
    double precision:: rhoredmix_orig, tredmix_orig, Rmix
    double precision, dimension(30):: dChempot_dxj_sol

    integer:: i, j, k

    !Variables for the mixed hydrate model
    double precision, dimension(30):: lnf, dphiidT, dphiidp, dChempot_dxj_hyd
    double precision, dimension(30,30):: dlnfidXj
    double precision, dimension(30):: fug_guest, DchpwDf
    double precision:: DchpwDT, DchpwDp, dChempot_dT_hyd, dChempot_dp_hyd
    double precision, dimension(3,30):: CiJ
    double precision, dimension(3,30,30) :: doccupidfj
    double precision, dimension(30,30):: dxidfj, dxiHdxjFluid


    JacMatrix = 0.D0
    dChempot_dT_sol = 0.D0
    dChempot_dp_sol = 0.D0
    dChempoti_dxj_fluid1 = 0.D0
    dChempoti_dT_fluid1 = 0.D0
    dChempoti_dp_fluid1 = 0.D0
    dChempoti_dxj_fluid2 = 0.D0
    dChempoti_dT_fluid2 = 0.D0
    dChempoti_dp_fluid2 = 0.D0
    dChempoti_dxj_fluid3 = 0.D0
    dChempoti_dT_fluid3 = 0.D0
    dChempoti_dp_fluid3 = 0.D0
    dChempoti_dxj_fluid4 = 0.D0
    dChempoti_dT_fluid4 = 0.D0
    dChempoti_dp_fluid4 = 0.D0
    d2nadnidT = 0.D0
    errval = 0

    !The derivative with respect to x will only be different from 0, if the solid phase is a hydrate
    dChempot_dxj_hyd = 0.D0

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    !Fluid Phase 1
    !**************************************************************************************
    ! get the properties of phase 1
    gl%molfractions = x_fluid1
    call reduced_parameters_calc(gl,Temp)
    call R_mix_calc(gl,Rmix)
    !Calculate the derivative of the chemical potential with respect to T at constant
    !X and p for fluid phase 1.
    if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
        call d2na_dnidT_P (gl,Temp, gl%rho_vap, d2nadnidT, 0)
        call dna_dni(gl,Temp, gl%rho_vap, dnadni, 0)
        dChempoti_dT_fluid1 = d2nadnidT * Rmix * Temp + dnadni * Rmix
    end if
    !Calculate the derivative of the chemical potential with respect to p at constant
    !X and T for fluid phase 1
    if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
        call d2na_dnidp_T (gl,Temp, gl%rho_vap, dChempoti_dp_fluid1, 0)
        dChempoti_dp_fluid1 = dChempoti_dp_fluid1 * Rmix * Temp
    end if
    !Calculate the derivative of the chemical potential with respect to xj at constant
    !p and T for the vapor phase
    call d2na_dnidxj_PT (gl,Temp, gl%rho_vap, dChempoti_dxj_fluid1, 0)
    dChempoti_dxj_fluid1 = dChempoti_dxj_fluid1 * Rmix * Temp
    !Get the fugacities and derivatives of the vapor phase (THIS IS NEEDED FOR THE HYDRATE MODEL)
    if (gl%solidtype(2) .eq. 1) then
        call lnf_mix(gl,Temp, gl%rho_vap, press, lnf)
        call dlnfi_dxj_TP (gl,Temp, gl%rho_vap, dlnfidXj)
        call dlnphii_dT(gl,Temp, gl%rho_vap, dphiidT)
        call dlnphii_dp(gl,Temp, gl%rho_vap, dphiidp)
    end if
    !**************************************************************************************

    !Fluid Phase 2
    !**************************************************************************************
    if (nr_of_fluidphases > 1) then
        ! get the properties of the fluid phase 2
        gl%molfractions = x_fluid2
        call reduced_parameters_calc(gl,Temp)
        call R_mix_calc(gl,Rmix)
        !Calculate the derivative of the chemical potential with respect to T at constant
        !X and p for fluid phase 2
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
            call d2na_dnidT_P (gl,Temp, gl%rho_liq1, d2nadnidT, 0)
            call dna_dni(gl,Temp, gl%rho_liq1, dnadni, 0)
            dChempoti_dT_fluid2 = d2nadnidT * Rmix * Temp + dnadni * Rmix
        end if
        !Calculate the derivative of the chemical potential with respect to p at constant
        !X and T for fluid phase 2
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
            call d2na_dnidp_T (gl,Temp, gl%rho_liq1, dChempoti_dp_fluid2, 0)
            dChempoti_dp_fluid2 = dChempoti_dp_fluid2 * Rmix * Temp
        end if
        !Calculate the derivative of the chemical potential with respect to xj at constant
        !p and T for fluid phase 2
        call d2na_dnidxj_PT (gl,Temp, gl%rho_liq1, dChempoti_dxj_fluid2, 0)
        dChempoti_dxj_fluid2 = dChempoti_dxj_fluid2 * Rmix * Temp
    End if
    !**************************************************************************************

    !Fluid Phase 3
    !**************************************************************************************
    if (nr_of_fluidphases > 2) then
        ! get the properties of the fluid phase 3
        gl%molfractions = x_fluid3
        call reduced_parameters_calc(gl,Temp)
        call R_mix_calc(gl,Rmix)
        !Calculate the derivative of the chemical potential with respect to T at constant
        !X and p for fluid phase 3
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
            call d2na_dnidT_P (gl,Temp, gl%rho_liq2, d2nadnidT, 0)
            call dna_dni(gl,Temp, gl%rho_liq2, dnadni, 0)
            dChempoti_dT_fluid3 = d2nadnidT * Rmix * Temp + dnadni * Rmix
        end if
        !Calculate the derivative of the chemical potential with respect to p at constant
        !X and T for fluid phase 3
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
            call d2na_dnidp_T (gl,Temp, gl%rho_liq2, dChempoti_dp_fluid3, 0)
            dChempoti_dp_fluid3 = dChempoti_dp_fluid3 * Rmix * Temp
        end if
        !Calculate the derivative of the chemical potential with respect to xj at constant
        !p and T for fluid phase 3
        call d2na_dnidxj_PT (gl,Temp, gl%rho_liq2, dChempoti_dxj_fluid3, 0)
        dChempoti_dxj_fluid3 = dChempoti_dxj_fluid3 * Rmix * Temp
    End if
    !**************************************************************************************

    !Fluid Phase 4
    !**************************************************************************************
    if (nr_of_fluidphases > 3) then
        ! get the properties of the fluid phase 4
        gl%molfractions = x_fluid4
        call reduced_parameters_calc(gl,Temp)
        call R_mix_calc(gl,Rmix)
        !Calculate the derivative of the chemical potential with respect to T at constant
        !X and p for fluid phase 4
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
            call d2na_dnidT_P (gl,Temp, gl%rho_liq3, d2nadnidT, 0)
            call dna_dni(gl,Temp, gl%rho_liq3, dnadni, 0)
            dChempoti_dT_fluid4 = d2nadnidT * Rmix * Temp + dnadni * Rmix
        end if
        !Calculate the derivative of the chemical potential with respect to p at constant
        !X and T for fluid phase 4
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
            call d2na_dnidp_T (gl,Temp, gl%rho_liq3, dChempoti_dp_fluid4, 0)
            dChempoti_dp_fluid4 = dChempoti_dp_fluid4 * Rmix * Temp
        end if
        !Calculate the derivative of the chemical potential with respect to xj at constant
        !p and T for fluid phase 4
        call d2na_dnidxj_PT (gl,Temp, gl%rho_liq3, dChempoti_dxj_fluid4, 0)
        dChempoti_dxj_fluid4 = dChempoti_dxj_fluid4 * Rmix * Temp
    End if
    !**************************************************************************************

    if (gl%solidtype(2) .eq. 1) then
        !Hydrate phase
        !**************************************************************************************
        ! Get the hydrate properties
        !Write the fugacities of the hydrates formers on position 1 to nrofhydrateformers-1 of the vector fug_guest
        fug_guest(1:gl%nrofhydrateformers-1) = dexp(lnf(2:gl%nrofhydrateformers))*1.d6

        !Get the derivative of the chemical potential of water in the hydrate phase with respect to the fugacity of the guest molecules
        call hdrt_Df_chem_potent_w(gl,Temp,press*1.D6,fug_guest,DchpwDf)

        !The derivative with respect to Temperature is calculated according to: dcp_dT|p,x_CO2_vap = dcp_dT|p,f_CO2 + dcp_df|T,p * dfi_dT|p, x_CO2
        ! pure hydrate: dChempot_dT_hyd = DchpwDT + DchpwDf * dphiidT(pos_CO2)*fug_CO2
        ! mixed hydrate 2015:
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
            call hdrt_DT_chem_potent_w(gl,Temp,press*1.D6,fug_guest,DchpwDT)
            dChempot_dT_hyd = DchpwDT
            Do i = 2, gl%nrofhydrateformers
                dChempot_dT_hyd = dChempot_dT_hyd + DchpwDf(i) * dphiidT(i)*dexp(lnf(i))*1.d6    ! Eq. 5.173 in Andreas' diss.
            End do
        end if

        !The derivative with respect to pressure is calculated according to: dcp_dp|T,x_CO2_vap = dcp_dp|p,f_CO2 + dcp_df|T,p * dfi_dp|T, x_CO2
        ! pure hydrate: dChempot_dp_hyd = DchpwDp + DchpwDf * fug_CO2 * (dphiidp(pos_CO2) + 1.D0 / (p*1.D6))
        ! mixed hydrate 2015:
        if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
            call hdrt_Dp_chem_potent_w(gl,Temp,press*1.D6,fug_guest,DchpwDp)
            dChempot_dp_hyd = DchpwDp
            Do i = 2, gl%nrofhydrateformers
                dChempot_dp_hyd = dChempot_dp_hyd + DchpwDf(i) * dexp(lnf(i))*1.d6 * (dphiidp(i) + 1.D0 / (press*1.D6))    ! Eq. 5.175 in diss.
            End do
        end if

        ! The derivative with respect to x_CO2 is calculated according to: dcp_dx_CO2|p,T = dcp_df|T,p * dfi_dx_CO2|T,p
        ! pure hydrate: dChempot_dxj_hyd = DchpwDf * dlnfidXj(1, pos_CO2)*fug_CO2
        ! mixed hydrate 2015:
        dChempot_dxj_hyd = 0.d0
        Do j = 1, gl%ncomp - 1 ! x(j)
            Do i = 2, gl%nrofhydrateformers
                ! mixed hydrate 2015:
                dChempot_dxj_hyd(j) = dChempot_dxj_hyd(j) + DchpwDf(i) * dlnfidxj(j, i)*dexp(lnf(i))*1.d6    ! Eq. 5.176 in Andreas' diss.
            End do
        End Do

        !Derivative of the molfraction of component i in hydrate wrt. the molfraction j in the fluid
        !If the component is not present in the hydrate phase, the derivative remains 0
        call hdrt_dxi_dfj(gl,Temp,press*1.D6,fug_guest,doccupidfj,CiJ,dxidfj)
        Do j = 1, gl%ncomp - 1
            Do i = 1, gl%nrofhydrateformers
                ! mixed hydrate 2015:
                dxiHdxjFluid(j,i) = 0.d0
                Do k = 2, gl%nrofhydrateformers    ! k-guest fugacity
                    dxiHdxjFluid(j,i) = dxiHdxjFluid(j,i) + dxidfj(k,i) * dlnfidxj(j, k)*dexp(lnf(k))*1.d6    ! Eq. 5.179 in Andreas' diss.
                End do
            End Do
        End Do

        !**************************************************************************************
    end if


    if (gl%solidtype(1) .gt. 0) then
        !Pure solid phase
        !**************************************************************************************
        select case (gl%solidtype(1))
        case(1)
            if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
                dChempot_dT_sol = dgdT_WaterIce(gl,Temp, press)   !Derivative of the Chemical Potential of Water Ice [J / mol K]
            end if
            if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
                dChempot_dp_sol = dgdp_WaterIce(gl,Temp, press)   !Derivative of the Chemical Potential of Water Ice [J / mol Pa]
            end if
        case(2)
            if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 1)) then !Derivative needed for binary mixtures and if iFlash = 1
                dChempot_dT_sol = dgdT_DryIce(gl,Temp, press)     !Derivative of the Chemical Potential of Dry Ice [J / mol K]
            end if
            if ((gl%ncomp .eq. 2) .or. (iFlash .eq. 2)) then !Derivative needed for binary mixtures and if iFlash = 2
                dChempot_dp_sol = dgdp_DryIce(gl,Temp, press)     !Derivative of the Chemical Potential of Dry Ice [J / mol Pa]
            end if
        end select
        dChempot_dxj_sol = 0.D0 !Solid is assumed to be pure --> no derivative with respect to composition
        !**************************************************************************************
    end if


    !Calculate the Jacobian matrix

    !JACOBIAN FOR BINARY MIXTURES
    !********************************************************************************************************************************
    If (gl%ncomp .eq. 2) then
        if ( (gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1) ) then !VLxLwH, LyLxLwH

            JacMatrix(1, 1) = dChempoti_dxj_fluid1(1, 1)
            JacMatrix(1, 2) = dChempoti_dxj_fluid1(1, 2)
            JacMatrix(1, 3) = dChempoti_dxj_fluid1(1, 1)
            JacMatrix(1, 4) = dChempoti_dxj_fluid1(1, 2)
            JacMatrix(1, 5) = dChempoti_dxj_fluid1(1, 1) - dChempot_dxj_hyd(1) !Derivative of chemical potential of water in the hydrate phase with respect to vapor composition

            JacMatrix(2, 1) = -dChempoti_dxj_fluid2(1, 1)
            JacMatrix(2, 2) = -dChempoti_dxj_fluid2(1, 2)
            JacMatrix(2, 3) = 0.D0
            JacMatrix(2, 4) = 0.D0
            JacMatrix(2, 5) = 0.D0

            JacMatrix(3, 1) = 0.D0
            JacMatrix(3, 2) = 0.D0
            JacMatrix(3, 3) = -dChempoti_dxj_fluid3(1, 1)
            JacMatrix(3, 4) = -dChempoti_dxj_fluid3(1, 2)
            JacMatrix(3, 5) = 0.D0

            JacMatrix(4, 1) = dChempoti_dT_fluid1(1) - dChempoti_dT_fluid2(1)
            JacMatrix(4, 2) = dChempoti_dT_fluid1(2) - dChempoti_dT_fluid2(2)
            JacMatrix(4, 3) = dChempoti_dT_fluid1(1) - dChempoti_dT_fluid3(1)
            JacMatrix(4, 4) = dChempoti_dT_fluid1(2) - dChempoti_dT_fluid3(2)
            JacMatrix(4, 5) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd         !ChemPot of water is needed --> water on position 1 at the moment

            JacMatrix(5, 1) = dChempoti_dp_fluid1(1) - dChempoti_dp_fluid2(1)
            JacMatrix(5, 2) = dChempoti_dp_fluid1(2) - dChempoti_dp_fluid2(2)
            JacMatrix(5, 3) = dChempoti_dp_fluid1(1) - dChempoti_dp_fluid3(1)
            JacMatrix(5, 4) = dChempoti_dp_fluid1(2) - dChempoti_dp_fluid3(2)
            JacMatrix(5, 5) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd         !ChemPot of water is needed --> water on position 1 at the moment

        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then  !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc

            JacMatrix(1, 1) = dChempoti_dxj_fluid1(1, 1)
            JacMatrix(1, 2) = dChempoti_dxj_fluid1(1, 2)
            JacMatrix(1, 3) = dChempoti_dxj_fluid1(1, 1) - dChempot_dxj_hyd(1)    !dChempot_dxj_hyd(2) = d琯w^H / dx_H2O^V
            JacMatrix(1, 4) = dChempoti_dxj_fluid1(1, gl%solid_pos)

            JacMatrix(2, 1) = -dChempoti_dxj_fluid2(1, 1)
            JacMatrix(2, 2) = -dChempoti_dxj_fluid2(1, 2)
            JacMatrix(2, 3) = 0.D0
            JacMatrix(2, 4) = 0.D0

            JacMatrix(3, 1) = dChempoti_dT_fluid1(1) - dChempoti_dT_fluid2(1)
            JacMatrix(3, 2) = dChempoti_dT_fluid1(2) - dChempoti_dT_fluid2(2)
            JacMatrix(3, 3) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd         !ChemPot of water is needed --> water on position 1 at the moment
            JacMatrix(3, 4) = dChempoti_dT_fluid1(gl%solid_pos) - dChempot_dT_sol   !dChempoti_dT_fluid1(solid_pos) - dChempot_dT_sol

            JacMatrix(4, 1) = dChempoti_dp_fluid1(1) - dChempoti_dp_fluid2(1)
            JacMatrix(4, 2) = dChempoti_dp_fluid1(2) - dChempoti_dp_fluid2(2)
            JacMatrix(4, 3) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd         !ChemPot of water is needed --> water on position 1 at the moment
            JacMatrix(4, 4) = dChempoti_dp_fluid1(gl%solid_pos) - dChempot_dp_sol !dChempoti_dp_fluid1(solid_pos) - dChempot_dp_sol

        end if
    end if
    !********************************************************************************************************************************

    !JACOBIAN FOR MIXTURE WITH MORE THAN 2 COMPONENTS
    !********************************************************************************************************************************
    If (gl%ncomp .gt. 2) then
        if ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 1)) then         !VLxLwH, LyLxLwH

            !i: Functions (See Sys_of_Eqs), j: independent variables (x1V, x2V, ..., xN-1V, x1Lw, x2Lw, ..., xN-1Lw, T or p or beta(4), beta(1), beta(2), beta(3)

            !In this block, the derivatives of the first ncomp equations with respect to ALL 2*ncomp + 2 independent variables are put into the Jacobian
            !-----
            Do i = 1, gl%ncomp !Loop over the Equations 猩=無x (or 無y=無x) and 猩=無w (or 無y=無w)
                Do j = 1, gl%ncomp-1 !Loop over the variables xV, xLx, and xLw (or xLy, xLx and xLw)
                    Jacmatrix(j, i) = dChempoti_dxj_fluid1(j, i)                !Derivative of first N equations with respect to xV (or xLy)
                    Jacmatrix(gl%ncomp - 1 + j, i) = -dChempoti_dxj_fluid2(j, i)   !Derivative of first N equations with respect to xLx
                    Jacmatrix(2*gl%ncomp - 2 + j, i) = 0.D0                        !Derivative of first N equations with respect to xLw

                    Jacmatrix(j, gl%ncomp + i) = dChempoti_dxj_fluid1(j, i)                             !Derivative of second N equations with respect to xV (or xLy)
                    Jacmatrix(gl%ncomp - 1 + j, gl%ncomp + i) = 0.D0                                       !Derivative of second N equations with respect to xLx
                    Jacmatrix(2*gl%ncomp - 2 + j, gl%ncomp + i) = -dChempoti_dxj_fluid3(j, i)              !Derivative of second N equations with respect to xLw
                end do
                if (iflash .eq. 1) then
                    !Derivative with respect to temperature
                    Jacmatrix(3*gl%ncomp - 3 + 1, i) = dChempoti_dT_fluid1(i) - dChempoti_dT_fluid2(i)         !Derivative of first N equations with respect to temperature
                    Jacmatrix(3*gl%ncomp - 3 + 1, gl%ncomp + i) = dChempoti_dT_fluid1(i) - dChempoti_dT_fluid3(i) !Derivative of second N equations with respect to temperature
                elseif (iflash .eq. 2) then
                    !Derivative with respect to pressure
                    Jacmatrix(3*gl%ncomp - 3 + 1, i) = dChempoti_dp_fluid1(i) - dChempoti_dp_fluid2(i)             !Derivative of first N equations with respect to pressure
                    Jacmatrix(3*gl%ncomp - 3 + 1, gl%ncomp + i) = dChempoti_dp_fluid1(i) - dChempoti_dp_fluid3(i)     !Derivative of second N equations with respect to pressure
                elseif (iflash .eq. 3) then
                    !Derivative with respect to beta
                    Jacmatrix(3*gl%ncomp - 3 + 1, i) = 0.D0            !Derivative of first N equations with respect to beta
                    Jacmatrix(3*gl%ncomp - 3 + 1, gl%ncomp + i) = 0.D0    !Derivative of second N equations with respect to beta
                end if
                !Derivatives with respect to betas for first N equations
                Jacmatrix(3*gl%ncomp - 3 + 2, i) = 0.D0
                Jacmatrix(3*gl%ncomp - 3 + 3, i) = 0.D0
                Jacmatrix(3*gl%ncomp - 3 + 4, i) = 0.D0
                !Derivatives with respect to betas for second N equations
                Jacmatrix(3*gl%ncomp - 3 + 2, gl%ncomp + i) = 0.D0
                Jacmatrix(3*gl%ncomp - 3 + 3, gl%ncomp + i) = 0.D0
                Jacmatrix(3*gl%ncomp - 3 + 4, gl%ncomp + i) = 0.D0
            end do
            !-----


            !In the following block, the derivatives of equation 2*ncomp + 1 with respect to ALL 3*ncomp + 1 independent variables are put into the Jacobian
            !-----
            Do j = 1, gl%ncomp - 1
                Jacmatrix(j, 2*gl%ncomp + 1) = dChempoti_dxj_fluid1(j, 1) - dChempot_dxj_hyd(j)    !Derivative of equation 2*ncomp+1 with respect to xV (or xLy)
                Jacmatrix(gl%ncomp - 1 + j, 2*gl%ncomp + 1) = 0.D0                                    !Derivative of equation 2*ncomp+1 with respect to xLx
                Jacmatrix(2*gl%ncomp - 2 + j, 2*gl%ncomp + 1) = 0.D0                                  !Derivative of equation 2*ncomp+1 with respect to xLw
            end do
            if (iflash .eq. 1) then
                !Derivative of equation 2*ncomp+1 with respect to temperature
                Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd
            elseif (iflash .eq. 2) then
                !Derivative of equation 2*ncomp+1 with respect to pressure
                Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd
            elseif (iflash .eq. 3) then
                !Derivative of equation 2*ncomp+1 with respect to beta
                Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1) = 0.D0
            end if
            !Derivatives with respect to betas
            Jacmatrix(3*gl%ncomp - 3 + 2, 2*gl%ncomp + 1) = 0.D0
            Jacmatrix(3*gl%ncomp - 3 + 3, 2*gl%ncomp + 1) = 0.D0
            Jacmatrix(3*gl%ncomp - 3 + 4, 2*gl%ncomp + 1) = 0.D0
            !-----


            !In the following block, the derivatives of equations 2*ncomp + 2 to 3*ncomp+1 (Material balance) with respect to ALL 2*ncomp + 2 independet variables are put into the Jacobian
            !-----
            Do i = 1, gl%ncomp-1
                Do j = 1, gl%ncomp - 1
                    Jacmatrix(j, 2*gl%ncomp + 1 + i) = 0.D0                     !Derivative of equations 2*ncomp+2 to 3*ncomp with respect to xV (or Ly)
                    Jacmatrix(gl%ncomp - 1 + j, 2*gl%ncomp + 1 + i) = 0.D0         !Derivative of equations 2*ncomp+2 to 3*ncomp with respect to xLx
                    Jacmatrix(2*gl%ncomp - 2 + j, 2*gl%ncomp + 1 + i) = 0.D0       !Derivative of equations 2*ncomp+2 to 3*ncomp with respect to xLw
                end do
                if (iflash .eq. 1) then
                    !Derivative with respect to temperature
                    Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1 + i) = 0.D0      !Derivative of equations 2*ncomp+2 to 3*ncomp with respect to temperature
                elseif (iflash .eq. 2) then
                    !Derivative with respect to pressure
                    Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1 + i) = 0.D0      !Derivative of equations 2*ncomp+2 to 3*ncomp with respect to temperature
                elseif (iflash .eq. 3) then
                    !Derivative with respect to beta
                    Jacmatrix(3*gl%ncomp - 3 + 1, 2*gl%ncomp + 1 + i) = x_fluid1(i)
                    Jacmatrix(3*gl%ncomp - 3 + 2, 2*gl%ncomp + 1 + i) = x_fluid2(i)
                    Jacmatrix(3*gl%ncomp - 3 + 3, 2*gl%ncomp + 1 + i) = x_fluid3(i)
                    Jacmatrix(3*gl%ncomp - 3 + 4, 2*gl%ncomp + 1 + i) = x_hyd(i)
                end if
                If ((iflash .eq. 1) .or. (iflash .eq. 2)) then
                    if (mapp_index(1) .eq. 1) then
                        Jacmatrix(3*gl%ncomp - 3 + 2, 2*gl%ncomp + 1 + i) = x_fluid1(i)
                    elseif (mapp_index(1) .eq. 2) then
                        Jacmatrix(3*gl%ncomp - 3 + 2, 2*gl%ncomp + 1 + i) = x_fluid2(i)
                    end if
                    if (mapp_index(2) .eq. 2) then
                        Jacmatrix(3*gl%ncomp - 3 + 3, 2*gl%ncomp + 1 + i) = x_fluid2(i)
                    elseif (mapp_index(2) .eq. 3) then
                        Jacmatrix(3*gl%ncomp - 3 + 3, 2*gl%ncomp + 1 + i) = x_fluid3(i)
                    end if
                    if (mapp_index(3) .eq. 3) then
                        Jacmatrix(3*gl%ncomp - 3 + 4, 2*gl%ncomp + 1 + i) = x_fluid3(i)
                    elseif (mapp_index(3) .eq. 4) then
                        Jacmatrix(3*gl%ncomp - 3 + 4, 2*gl%ncomp + 1 + i) = x_hyd(i)
                    end if
                end if
            end do
            !Last equation (3*ncomp+1)
            Jacmatrix(1:(3*gl%ncomp + 1) ,3*gl%ncomp + 1) = 0.D0          !Derivatives with respect to xV (xLy), xLx, and xLw
            if (iflash .eq. 1) then
                !Derivative with respect to temperature
                Jacmatrix(3*gl%ncomp - 3 + 1, 3*gl%ncomp + 1) = 0.D0
            elseif (iflash .eq. 2) then
                !Derivative with respect to pressure
                Jacmatrix(3*gl%ncomp - 3 + 1, 3*gl%ncomp + 1) = 0.D0
            elseif (iflash .eq. 3) then
                !Derivative with respect to beta
                Jacmatrix(3*gl%ncomp - 3 + 1, 3*gl%ncomp + 1) = 1.D0
            end if
            Jacmatrix(3*gl%ncomp - 3 + 2 ,3*gl%ncomp+1) = 1.D0
            Jacmatrix(3*gl%ncomp - 3 + 3 ,3*gl%ncomp+1) = 1.D0
            Jacmatrix(3*gl%ncomp - 3 + 4 ,3*gl%ncomp+1) = 1.D0
            !-----

        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 0)) then     !VLxLwIw, LyLxLwIw, VLxLwIc, LyLxLwIc
            !NOT YET IMPLEMENTED


        elseif ((gl%solidtype(1) .gt. 0) .and. (gl%solidtype(2) .eq. 1)) then     !VLwHIw, LxLwHIw, VLwHIc, LxLwHIc

            !i: Functions (See Sys_of_Eqs), j: independent variables (x1V, x2V, ..., xN-1V, x1Lw, x2Lw, ..., xN-1Lw, T or p or beta(4), beta(1), beta(2), beta(3)

            !In this block, the derivatives of the first ncomp equations with respect to ALL 2*ncomp + 2 independent variables are put into the Jacobian
            !-----
            Do i = 1, gl%ncomp !Loop over the Equations 猩=無w (or 無x=無w)
                Do j = 1, gl%ncomp-1 !Loop over the variables xV and xLw (or xLx xLw)
                    Jacmatrix(j, i) = dChempoti_dxj_fluid1(j, i)
                    Jacmatrix(gl%ncomp - 1 + j, i) = -dChempoti_dxj_fluid2(j, i)
                end do
                if (iflash .eq. 1) then
                    !Derivative with respect to temperature
                    Jacmatrix(2*gl%ncomp - 2 + 1, i) = dChempoti_dT_fluid1(i) - dChempoti_dT_fluid2(i)
                elseif (iflash .eq. 2) then
                    !Derivative with respect to pressure
                    Jacmatrix(2*gl%ncomp - 2 + 1, i) = dChempoti_dp_fluid1(i) - dChempoti_dp_fluid2(i)
                elseif (iflash .eq. 3) then
                    !Derivative with respect to beta
                    Jacmatrix(2*gl%ncomp - 2 + 1, i) = 0.D0
                end if
                !Derivatives with respect to betas
                Jacmatrix(2*gl%ncomp - 2 + 2, i) = 0.D0
                Jacmatrix(2*gl%ncomp - 2 + 3, i) = 0.D0
                Jacmatrix(2*gl%ncomp - 2 + 4, i) = 0.D0
            end do
            !-----


            !In the following block, the derivatives of equations ncomp + 1 and ncomp + 2 with respect to ALL 2*ncomp + 2 independet variables are put into the Jacobian
            !-----
            Do j = 1, gl%ncomp - 1
                Jacmatrix(j, gl%ncomp + 1) = dChempoti_dxj_fluid1(j, 1) - dChempot_dxj_hyd(j)
                Jacmatrix(gl%ncomp - 1 + j, gl%ncomp + 1) = 0.D0
                Jacmatrix(j, gl%ncomp + 2) = dChempoti_dxj_fluid1(j, gl%solid_pos)
                Jacmatrix(gl%ncomp - 1 + j, gl%ncomp + 2) = 0.D0
            end do
            if (iflash .eq. 1) then
                !Derivative with respect to temperature
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 1) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2) = dChempoti_dT_fluid1(gl%solid_pos) - dChempot_dT_sol
            elseif (iflash .eq. 2) then
                !Derivative with respect to pressure
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 1) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2) = dChempoti_dp_fluid1(gl%solid_pos) - dChempot_dp_sol
            elseif (iflash .eq. 3) then
                !Derivative with respect to beta
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 1) = 0.D0
                Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2) = 0.D0
            end if
            !Derivatives with respect to betas
            Jacmatrix(2*gl%ncomp - 2 + 2, gl%ncomp + 1) = 0.D0
            Jacmatrix(2*gl%ncomp - 2 + 2, gl%ncomp + 2) = 0.D0
            Jacmatrix(2*gl%ncomp - 2 + 3, gl%ncomp + 1) = 0.D0
            Jacmatrix(2*gl%ncomp - 2 + 3, gl%ncomp + 2) = 0.D0
            Jacmatrix(2*gl%ncomp - 2 + 4, gl%ncomp + 1) = 0.D0
            Jacmatrix(2*gl%ncomp - 2 + 4, gl%ncomp + 2) = 0.D0
            !-----


            !In the following block, the derivatives of equations ncomp + 3 to 2*ncomp+2 (Material balance) with respect to ALL 2*ncomp + 2 independet variables are put into the Jacobian
            !-----
            Do i = 1, gl%ncomp-1
                Do j = 1, gl%ncomp - 1
                    Jacmatrix(j, gl%ncomp + 2 + i) = 0.D0
                    Jacmatrix(gl%ncomp - 1 + j, gl%ncomp + 2 + i) = 0.D0
                end do
                if (iflash .eq. 1) then
                    !Derivative with respect to temperature
                    Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2 + i) = 0.D0
                elseif (iflash .eq. 2) then
                    !Derivative with respect to pressure
                    Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2 + i) = 0.D0
                elseif (iflash .eq. 3) then
                    !Derivative with respect to beta
                    Jacmatrix(2*gl%ncomp - 2 + 1, gl%ncomp + 2 + i) = x_sol(i)        !Last position beta(4), always Iw or Ic (x_sol)
                    Jacmatrix(2*gl%ncomp - 2 + 2, gl%ncomp + 2 + i) = x_fluid1(i)
                    Jacmatrix(2*gl%ncomp - 2 + 3, gl%ncomp + 2 + i) = x_fluid2(i)
                    Jacmatrix(2*gl%ncomp - 2 + 4, gl%ncomp + 2 + i) = x_hyd(i)
                end if
                If ((iflash .eq. 1) .or. (iflash .eq. 2)) then
                    if (mapp_index(1) .eq. 1) then
                        Jacmatrix(2*gl%ncomp - 2 + 2, gl%ncomp + 2 + i) = x_fluid1(i)
                    elseif (mapp_index(1) .eq. 2) then
                        Jacmatrix(2*gl%ncomp - 2 + 2, gl%ncomp + 2 + i) = x_fluid2(i)
                    end if
                    if (mapp_index(2) .eq. 2) then
                        Jacmatrix(2*gl%ncomp - 2 + 3, gl%ncomp + 2 + i) = x_fluid2(i)
                    elseif (mapp_index(2) .eq. 3) then
                        Jacmatrix(2*gl%ncomp - 2 + 3, gl%ncomp + 2 + i) = x_hyd(i)
                    end if
                    if (mapp_index(3) .eq. 3) then
                        Jacmatrix(2*gl%ncomp - 2 + 4, gl%ncomp + 2 + i) = x_hyd(i)
                    elseif (mapp_index(3) .eq. 4) then
                        Jacmatrix(2*gl%ncomp - 2 + 4, gl%ncomp + 2 + i) = x_sol(i)
                    end if
                end if
            end do
            !Last equation
            Jacmatrix(1:(2*gl%ncomp-2) ,2*gl%ncomp+2) = 0.D0
            if (iflash .eq. 1) then
                !Derivative with respect to temperature
                Jacmatrix(2*gl%ncomp - 2 + 1, 2*gl%ncomp+2) = 0.D0
            elseif (iflash .eq. 2) then
                !Derivative with respect to pressure
                Jacmatrix(2*gl%ncomp - 2 + 1, 2*gl%ncomp+2) = 0.D0
            elseif (iflash .eq. 3) then
                !Derivative with respect to beta
                Jacmatrix(2*gl%ncomp - 2 + 1, 2*gl%ncomp+2) = 1.D0
            end if
            Jacmatrix(2*gl%ncomp - 2 + 2 ,2*gl%ncomp+2) = 1.D0
            Jacmatrix(2*gl%ncomp - 2 + 3 ,2*gl%ncomp+2) = 1.D0
            Jacmatrix(2*gl%ncomp - 2 + 4 ,2*gl%ncomp+2) = 1.D0
            !-----


        elseif ((gl%solidtype(1) .eq. 0) .and. (gl%solidtype(2) .eq. 0)) then     !VLxLyLz, LwLxLyLz
            !NOT YET IMPLEMENTED

        end if
    end if
    !********************************************************************************************************************************


    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = x_known


    end subroutine Jacobi_solid_NC_4P
    !**************************************************************************


    !**************************************************************************
    module subroutine ptflash_solid_NC_3P(gl,press, Temp, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, Phasefrac, iFlash, iphase, iter, errval)
    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THREE-PHASE-EQUILIBRIA FOR MULTI COMPONENTS MIXTURES FORMING SOLID PHASES. FOLLOWING THE GIBBS' PHASE RULE
    ! THE THREE PHASE REGIONS ARE LINES IN THE T-P-DIAGRAM FOR BINARY MIXTURES AND MAY BE AREAS FOR THREE AND MORE COMPONENT MIXTURES
    !
    ! THE FOLLOWING CALCULATIONS MIGHT BE PERFORMED
    ! Case 1:   VLH / LLH   Indicators: solidtype(1) = 0, solidtype(2) = 1
    ! Case 2:   VHS / LHS   Indicators: solidtype(1) /= 0, solidtype(2) = 1
    ! Case 3:   VLS / LLS   Indicators: solidtype(1) /= 0, solidtype(2) = 0
    !
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   - THREE PHASE LINE ("Melting" Point)         : p, x Given beta_'fluid1' = 0 (VHS: beta_S = 0)              --  iFlash =  1
    !   - THREE PHASE LINE ("Melting" Point)         : p, x Given beta_'fluid1' = 0 (VHS: beta_'fluid1' = 0)       --  iFlash = -1
    !   - THREE PHASE LINE ("Freezing" Point)        : p, x Given (VLH and VHS: beta_H = 0; VLS: beta_S = 0)       --  iFlash =  2
    !   - THREE PHASE LINE ("Dew" Point)             : p, x Given beta_'fluid2' = 0 (VLH: betaL=0, L1L2H: betaL2=0)--  iFlash =  3
    !   - THREE PHASE LINE ("Melting" Point)         : T, x Given beta_'fluid1' = 0 (VHS: beta_H = 0)              --  iFlash =  4
    !   - THREE PHASE LINE ("Freezing" Point)        : T, x Given (VLH and VSH: beta_H = 0; VLS: beta_S = 0)       --  iFlash =  5
    !   - THREE PHASE LINE ("Dew" Point)             : T, x Given beta_'fluid2' = 0 (VLH: betaL=0, L1L2H: betaL2=0)--  iFlash =  6
    !   - THREE PHASE Flash: p, T and x VECTOR GIVEN                                                               --  iFlash =  7
    !   - TRIPLE POINT OF PURE SUBSTANCE             : x given, T and p unknown                                    --  iFlash =  8
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press           - Pressure [MPa]
    !   Temp            - Temperature [K]
    !   x_fluid1        - Composition of the first fluid phase. In case of VLH it is V, in case of LHS it is L, in case of LLH it is the "lighter" liquid
    !   x_fluid2        - Composition of the second fluid phase, if there is any. in all cases this is the "heavier" fluid phase
    !   x_sol           - (Pure) Solid phase composition
    !   x_hyd           - Hydrate phase composition
    !   rho_fluid1_est  - estimated density of the first fluid phase
    !   rho_fluid2_est  - estimated density of the second fluid phase
    !   iFlash          - Flash type
    !   iphase          - specify which fluid phases are assumed in case of two fluid phases:   0 - Let the density solver choose based on the gibbs energy
    !                                                                                           1 - Assume liquid / liquid equilibrium
    !                                                                                           2 - Assume vapor / liquid equilibrium
    !                     In case of only one liquid phase (e.g. VHIw)                          0 - Let the density solver choose based on the gibbs energy
    !                                                                                           1 - Assume liquid
    !                                                                                           2 - Assume vapor
    !   Phasefrac       - Phasefractions according to:  Phasefrac(1)    :   nV / n      - Molar amount of the first phase divided by overall molar amount of the mixture
    !                                                   Phasefrac(2)    :   nL1 / n     - Molar amount of the second phase divided by overall molar amount of the mixture
    !                                                   Phasefrac(3)    :   nL2 / n     - Molar amount of the third phase
    !                   The following cases of equilibria might appear: VLH / LLH   phasefrac(1) = beta_V or beta_L1, phasefrac(2) = beta_H, phasefrac(3) = beta_L2
    !                                                                   VHS / LHS   phasefrac(1) = beta_V or beta_L1, phasefrac(2) = beta_H , phasefrac(3) = beta_S
    !                                                                   VLS / LLS   phasefrac(1) = beta_V or beta_L1, phasefrac(2) = beta_S , phasefrac(3) = beta_L2

    ! OUTPUT:
    !   errval          - Error value
    !   iter            - number of iterations carried out
    !   rho             - densities of all phases found are given back according to the order of the phase eq. type: (e.g. VLH --> rho(1) = Vapor density, rho(2) = Liquid density, rho(3) = hydrate density)
    !   x_fluid1        - Composition of the first fluid phase. In case of VLH it is V, in case of LHS it is L, in case of LLH it is the "lighter" liquid
    !   x_fluid2        - Composition of the second fluid phase, if there is any. in all cases this is the "heavier" fluid phase
    !   x_sol           - (Pure) Solid phase composition
    !   x_hyd           - Hydrate phase composition
    !   Phasefrac       - Phasefractions


    !------------------0--------------------------------------------------------
    ! A. J輍er August 2012, V. Vin?- October - December 2015, S. Hielscher - April 2016







    implicit none

    type(type_gl) :: gl


    double precision, dimension(30)::  x_known
    double precision:: rhofluid1_est, rhofluid2_est
    integer:: iFlash
    integer:: iPhase

    double precision:: press, temp
    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol, x_hyd
    double precision, dimension(3):: Phasefrac

    double precision, dimension(3) :: rho
    integer:: errval, iter

    double precision, dimension(3):: Phasefrac_new
    double precision, dimension(30):: x_fluid1_new, x_fluid2_new
    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(60):: GibbsEQN, Delta_X, Var_X

    double precision:: eps_Gibbs, eps_del, max_del, stepsize, eps_Gibbs_min
    double precision:: sum_fluid1, sum_fluid2, Temp_new, press_new
    integer:: i, j, k, eqn, fluidnr

    logical:: twofluidphases !Indicates if two fluid phases are present

    character(255):: herr

    !Variables for hydrates & pure solids
    !double precision :: v_DryIce, v_WaterIce, fug_save, rho_H
    !double precision, dimension (2):: occup, Cij

    !Variables for the mixed hydrates & pure solids
    double precision :: rho_H
    double precision, dimension(30):: fug_save
    double precision, dimension(3,30):: occup, CiJ


    errval = 0
    Delta_X = 1.D0
    x_fluid1_new = 0.D0
    x_fluid2_new = 0.D0
    phasefrac_new = 0.D0
    Temp_new = 0.D0
    press_new = 0.D0
    stepsize = 1.D0
    twofluidphases = .false.
    rho = 0.D0
    rho_H = 0.D0
    fug_save = 0.D0
    JacMatrix = 0.d0
    herr = ''

    !catch NaN
    if ((press /= press) .or. (Temp /= Temp)) then
        errval = -4321
        return
    endif

    !If the maximum difference of fugacities is lower than eps_Gibbs, the calculation is finished
    eps_Gibbs = 1.d-8
    !If the relative change of the unknowns (T,p or x) is below eps_del and the maximum value of Gibbs_EQN is smaller than eps_Gibbs_min, the calculation is finished
    eps_del = 1.d-4
    eps_Gibbs_min = 1.D-4

    !Iteration counter is set back to 1
    iter = 1
    !The GibbsEQN get a high starting value
    GibbsEQN = 1.D10
    !Starting values for max_del and Var_X
    !max_del = maximum relative difference for the set of unknowns between two iterations
    max_del = 0.D0
    !Current value for all variables
    Var_X = 1.D0

    If (((gl%solidtype(1) == 0) .and. (gl%solidtype(2) == 1)) .or. ((gl%solidtype(1) > 0) .and. (gl%solidtype(2) == 0))) then
        twofluidphases = .true.
    End if

    Do i = 1, gl%ncomp
        if (dabs(x_fluid1(i)) < 1.D-14) then
            errval = -1111   !Wrong (missing) inputs to a routine
            return
        End if
        If(twofluidphases) then !In case of 2 fluid phases two startvalues for the fluid compositions are needed
            if (dabs(x_fluid2(i)) < 1.D-14) then
                errval = -1111   !Wrong (missing) inputs to a routine
                return
            End if
        End if
    End Do

    !Check for starting values
    If ((dabs(Temp) < 1.D-14) .or. (dabs(press) < 1.D-14)) then
        errval = -1111   !Wrong (missing) inputs to a routine
        return
    End if

    if ((iFlash == -1) .or. (iFlash == 1) .or. (iFlash == 4)) phasefrac(1) = 0   !see head of the routine
    if ((iFlash == 2) .or. (iFlash == 5)) phasefrac(2) = 0   !see head of the routine
    if ((iFlash == 3) .or. (iFlash == 6)) phasefrac(3) = 0   !see head of the routine

    if ((iFlash == 8) .and. (gl%ncomp > 1)) then
        !Triple point calculation only possible for pure substances
        errval = -9902
        return
    end if

    do i = 1, 60

        call SysOfEqs_solid_NC_3P(gl,press, Temp, x_known, x_fluid1, x_fluid2, x_sol,x_hyd, rhofluid1_est, rhofluid2_est, fug_save,&
            & twofluidphases, iFlash, iphase, Phasefrac, GibbsEQN, errval)

        !Break Condition
        ! This is the first exit condition for the VLE iteration!
        if (maxval(dabs(GibbsEQN)) < eps_Gibbs) then
            if (twofluidphases) then
                rho(1) = gl%rho_vap
                rho(2) = gl%rho_liq2
                if (gl%solidtype(1) == 1) then         !solid H2O
                    rho(3) = 1.D0 / v_WaterIce(gl,Temp, press)
                elseif (gl%solidtype(1) == 2) then     !solid CO2
                    rho(3) = 1.D0 / v_DryIce(gl,Temp, press)
                elseif (gl%solidtype(2) == 1) then     !Hydrate
                    call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                    rho(3) = rho_H
                end if
            else
                rho(1) = gl%rho_vap
                if (gl%solidtype(1) == 1) then         !solid H2O
                    rho(2) = 1.D0 / v_WaterIce(gl,Temp, press)
                elseif (gl%solidtype(1) == 2) then     !solid CO2
                    rho(2) = 1.D0 / v_DryIce(gl,Temp, press)
                end if
                if (gl%solidtype(2) == 1) then     !Hydrate
                    call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                    rho(3) = rho_H
                end if
            end if

            return
        end if
        if (errval /= 0) return
        call Jacobi_solid_NC_3P(gl,press, temp, x_fluid1, x_fluid2, x_sol, x_hyd, fug_save, twofluidphases, iFlash, &
            & Phasefrac, JacMatrix, errval)
        Delta_X = - GibbsEQn

        If (twofluidphases) then
            If ((gl%ncomp > 3) .or. ((gl%ncomp == 3) .and. (iFlash < 7))) then
                eqn = 2*gl%ncomp
            else
                eqn = gl%ncomp + 1
            End if
        else    !2 components or 3 component and iFlash = 7
            If ((gl%ncomp > 3) .or. ((gl%ncomp == 3) .and. (iFlash < 7))) then
                eqn = gl%ncomp + 1
            else
                eqn = 2
            End if
        End if

        call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)

        if (errval /= 0) then
            errval = -4321
            return
        endif
        ! Initialize
        sum_fluid1 = 0.D0
        sum_fluid2 = 0.D0
        iter = i

        select case (iFlash)
            ! --------------------------------------------------------------------------------
        case(-1) !Melting point
            j = 1
            press_new = press
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        temp_new = temp + Delta_X(2*gl%NCOMP-1)
                        Var_X(2*gl%NCOMP-1) = temp_new
                    else
                        temp_new = temp + Delta_X(gl%NCOMP)
                        Var_X(gl%NCOMP) = temp_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(2)
                        else
                            phasefrac_new(2) = phasefrac(2) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(2)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(2) = 0.D0
                end if
                j = j + 1
            end do

        case(1) !Melting point
            j = 1
            press_new = press
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        temp_new = temp + Delta_X(2*gl%NCOMP-1)
                        Var_X(2*gl%NCOMP-1) = temp_new
                    else
                        temp_new = temp + Delta_X(gl%NCOMP)
                        Var_X(gl%NCOMP) = temp_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(2)
                        else
                            phasefrac_new(2) = phasefrac(2) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(2)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(2) = 0.D0
                end if
                j = j + 1
            end do
            ! --------------------------------------------------------------------------------
        case(2) !Freezing point
            j = 1
            press_new = press
            phasefrac_new(2) = phasefrac(2)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        temp_new = temp + Delta_X(2*gl%NCOMP-1)
                        Var_X(2*gl%NCOMP-1) = temp_new
                    else
                        temp_new = temp + Delta_X(gl%NCOMP)
                        Var_X(gl%NCOMP) = temp_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(1) = phasefrac(1) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(1)
                        else
                            phasefrac_new(1) = phasefrac(1) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(1)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(1) = 0.D0
                end if
                j = j + 1
            end do
            !Andreas, July 2014. For two components, calculate the phase fractions separately here (three phase
            If (gl%ncomp == 2) then
                if(twofluidphases) then
                    phasefrac_new(1) = (x_known(1) - x_fluid2_new(1)) / (x_fluid1_new(1) - x_fluid2_new(1))
                else
                    phasefrac_new(1) = (x_known(1) - x_sol(1)) / (x_fluid1_new(1) - x_sol(1))
                end if
            End if
            ! --------------------------------------------------------------------------------
        case(3) !Dew point
            j = 1
            press_new = press
            phasefrac_new(3) = phasefrac(3)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        temp_new = temp + Delta_X(2*gl%NCOMP-1)
                        Var_X(2*gl%NCOMP-1) = temp_new
                    else
                        temp_new = temp + Delta_X(gl%NCOMP)
                        Var_X(gl%NCOMP) = temp_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(2)
                        else
                            phasefrac_new(2) = phasefrac(2) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(2)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (temp_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(3) = 0.D0
                end if
                j = j + 1
            end do

            ! --------------------------------------------------------------------------------
        case(4) !!Melting point
            j = 1
            Temp_new = Temp
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        press_new = press + Delta_X(2*gl%NCOMP-1) / 1.D6
                        Var_X(2*gl%NCOMP-1) = press_new
                    else
                        press_new = press + Delta_X(gl%NCOMP) / 1.D6
                        Var_X(gl%NCOMP) = press_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(2)
                        else
                            phasefrac_new(2) = phasefrac(2) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(2)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(2) = 0.D0
                end if
                j = j + 1
            end do
            ! --------------------------------------------------------------------------------
        case(5) !Freezing point
            j = 1
            Temp_new = Temp
            phasefrac_new(1) = phasefrac(1)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        press_new = press + Delta_X(2*gl%NCOMP-1) / 1.D6
                        Var_X(2*gl%NCOMP-1) = press_new
                    else
                        press_new = press + Delta_X(gl%NCOMP) / 1.D6
                        Var_X(gl%NCOMP) = press_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(1) = phasefrac(1) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(1)
                        else
                            phasefrac_new(1) = phasefrac(1) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(1)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then  !&
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(1) = 0.D0
                end if
                j = j + 1
            end do
            ! --------------------------------------------------------------------------------
        case(6) !Dew point
            j = 1
            Temp_new = Temp
            phasefrac_new(3) = phasefrac(3)
            do while (j < gl%ncomp+1)
                if (j < gl%ncomp) then
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    if(twofluidphases) then
                        x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                        press_new = press + Delta_X(2*gl%NCOMP-1) / 1.D6
                        Var_X(2*gl%NCOMP-1) = press_new
                    else
                        press_new = press + Delta_X(gl%NCOMP) / 1.D6
                        Var_X(gl%NCOMP) = press_new
                    End if
                    If (gl%ncomp > 2) then
                        if(twofluidphases) then
                            phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%NCOMP)
                            Var_X(2*gl%NCOMP) = phasefrac_new(2)
                        else
                            phasefrac_new(2) = phasefrac(2) + Delta_X(gl%NCOMP+1)
                            Var_X(gl%NCOMP+1) = phasefrac_new(2)
                        End if
                    end if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) .OR. (press_new < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(2) < 0.d0) .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(3) = 0.D0
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
                    x_fluid1_new(j) = x_fluid1(j) + Delta_X(j)
                    Var_X(j) = x_fluid1_new(j)
                    if( twofluidphases) then
                        x_fluid2_new(j) = x_fluid2(j) + Delta_X(j+gl%NCOMP-1)
                        Var_X(j+gl%NCOMP-1) = x_fluid2_new(j)
                    end if
                else
                    x_fluid1_new(gl%ncomp) = 1.D0 - sum_fluid1
                    x_fluid2_new(gl%ncomp) = 1.D0 - sum_fluid2
                    if( twofluidphases) then
                        phasefrac_new(1) = phasefrac(1) + Delta_X(2*gl%ncomp-1)
                        phasefrac_new(2) = phasefrac(2) + Delta_X(2*gl%ncomp)
                        Var_X(2*gl%ncomp-1) = phasefrac_new(1)
                        Var_X(2*gl%ncomp) = phasefrac_new(2)
                    else
                        phasefrac_new(1) = phasefrac(1) + Delta_X(gl%ncomp)
                        phasefrac_new(2) = phasefrac(2) + Delta_X(gl%ncomp+1)
                        Var_X(gl%ncomp) = phasefrac_new(1)
                        Var_X(gl%ncomp+1) = phasefrac_new(2)
                    End if
                End if
                sum_fluid1 = sum_fluid1 + x_fluid1_new(j)
                sum_fluid2 = sum_fluid2 + x_fluid2_new(j)
                ! check if the stepsize is too large.
                ! The condition: x(j)_vap AND x(j)_liq1 > 0  AND x(j)_liq2 > 0
                ! has to be fulfilled for all components
                ! Otherwise the stepsize is reduced
                if ((x_fluid1_new(j) < 0.D0) .OR. (x_fluid2_new(j) < 0.D0) &
                    & .OR. (x_fluid1_new(j) > 1.D0) .OR. (x_fluid2_new(j) > 1.D0)) then !&
                    !& .OR. (phasefrac_new(1) < 0.d0) .OR. (phasefrac_new(1) > 1.d0) .OR. (phasefrac_new(2) < 0.d0) &
                    !& .OR. (phasefrac_new(2) > 1.d0)) then
                    stepsize = 2.D0
                    Delta_X = Delta_X/stepsize
                    if (maxval(dabs(Delta_X)) < 1.D-10) then
                        errval = -3333
                        return
                    end if
                    stepsize = 1.D0
                    j = 0
                    sum_fluid1 = 0.D0
                    sum_fluid2 = 0.D0
                    phasefrac_new(1) = 0.D0
                    phasefrac_new(2) = 0.D0
                end if
                j = j + 1
            end do
            ! --------------------------------------------------------------------------------
        case(8) !Triple Point of a pure substance, Andreas July 2014
            x_fluid1_new = x_fluid1
            x_fluid2_new = x_fluid2
            temp_new = temp + Delta_X(1)
            press_new = press + Delta_X(2) * 1.D-6
            !No step reduction, usually for pure substances this should quickly converge with
            !initial estimates set properly. If it diverges most likely inputs are wrong
            if ((temp_new < 0.D0) .or. (Temp_new > 400.D0) .or. (press_new < 0.D0) .or. (press_new > 10.D0)) then
                errval = -3333
                return
            end if
            ! --------------------------------------------------------------------------------
        end select

        ! write the new values to the variables
        x_fluid1 =  x_fluid1_new
        x_fluid2 = x_fluid2_new
        Temp = Temp_new
        press = press_new
        phasefrac = phasefrac_new
        if ((iflash /= 3) .and. (iflash /= 6)) then
            phasefrac(3) = 1.D0 - phasefrac(1) - phasefrac(2)
        elseif ((iflash == 3) .or. (iflash == 6)) then
            phasefrac(1) = 1.D0 - phasefrac(3) - phasefrac(2)
        endif

        !Andreas, July 2014
        !Dont do this for pure substances!
        if (iFlash /= 8) then
            !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
            !Additionally, check whether Gibbs_EQN is smaller than eps_Gibbs_min (minumum convergence criterion) if not --> error
            max_del = 0.D0
            Do k = 1, eqn
                if(abs(delta_X(k) / Var_X(k)) > max_del) then
                    max_del = abs(delta_X(k) / Var_X(k))
                end if
            end do

            if (((max_del) < eps_del).And. (maxval(dabs(GibbsEQN)) < eps_Gibbs_min)) then

                if (twofluidphases) then
                    rho(1) = gl%rho_vap
                    rho(2) = gl%rho_liq2
                    if (gl%solidtype(1) == 1) then         !solid H2O
                        rho(3) = 1.D0 / v_WaterIce(gl,Temp, press)
                    elseif (gl%solidtype(1) == 2) then     !solid CO2
                        rho(3) = 1.D0 / v_DryIce(gl,Temp, press)
                    elseif (gl%solidtype(2) == 1) then     !Hydrate
                        call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                        rho(3) = rho_H
                        return
                    end if
                else
                    rho(1) = gl%rho_vap
                    if (gl%solidtype(1) == 1) then         !solid H2O
                        rho(2) = 1.D0 / v_WaterIce(gl,Temp, press)
                    elseif (gl%solidtype(1) == 2) then     !solid CO2
                        rho(2) = 1.D0 / v_DryIce(gl,Temp, press)
                    end if
                    if (gl%solidtype(2) == 1) then     !Hydrate
                        call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                        rho(3) = rho_H
                    end if
                end if

                return

            end if

            !Catch the case, that the phase compositions become the same
            if(twofluidphases) then
                if (maxval(dabs(x_fluid1(1:gl%ncomp)-x_fluid2(1:gl%ncomp))) < 0.0000042d0) then
                    errval = -4323
                    return
                end if
            end if

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

    end subroutine


    !**************************************************************************
    module subroutine SysOfEqs_solid_NC_3P(gl,P, T, x_known, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, fug_save, twofluidphases, iFlash, iphaseIn, Phasefrac, GibbsEQN, errval)
    !**************************************************************************

    !------------------0--------------------------------------------------------
    ! A. J輍er,  August 2012, V. Vin?- October 2015








    implicit none

    type(type_gl) :: gl


    double precision, dimension(30):: x_hyd
    double precision:: P, T
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_sol
    double precision:: rhofluid1_est, rhofluid2_est
    double precision, dimension(3):: Phasefrac
    integer:: iFlash
    logical:: twofluidphases
    integer:: iPhaseIn

    double precision, dimension(60):: GibbsEQN
    integer:: errval

    double precision, dimension(30):: Chempot_fluid1, Chempot_fluid2, lnf
    double precision:: rhoredmix_orig, tredmix_orig, d_fluid1, d_fluid2
    double precision:: Rmix
    integer:: errorflag, iPhase, i

    !!Variables for pure hydrates & pure solids
    !double precision:: fug_CO2, Chempot_hyd, Chempot_sol, fug_save
    !!double precision, dimension(30):: fug_gas
    !double precision:: chpw, g_DryIce, g_WaterIce
    !double precision, dimension (2):: occup, Cij, x_hyd2

    !Variables for the mixed hydrates & pure solids
    double precision, dimension(30):: fug_gas, fug_save, x_hyd2!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: Chempot_hyd, Chempot_sol
    double precision:: chpw, p_save
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    logical :: seasave


    errval = 0
    GibbsEQN = 0.D0
    fug_save = 0.d0
    fug_gas = 0.d0

    !Just in case, maybe not necessary
    !------------------------------
    gl%molfractions = x_known
    call reduced_parameters_calc(gl,T)
    !------------------------------
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    ! Get the density of fluid phase 1
    gl%molfractions = x_fluid1
    call reduced_parameters_calc(gl,T)
    !The iPhase is now passed by the calling routine to allow for guesses on the phases
    iPhase = iPhaseIn
    !iPhase = 0  !Let the density solver choose
    ! calculate the fluid phase 1 density
    gl%rho_vap = rhomix_calc(gl,T, P, rhofluid1_est, iPhase,0)

    if (gl%rho_vap < 1.D-12)then
        errval = -8888
        return
    end if

    !----------------------------------
    ! get the fluid phase 1 properties
    d_fluid1 = gl%rho_vap
    ! Calculate the chemical potentials for the gas phase
    call R_mix_calc(gl,Rmix)
    !call dna_dni(T, d_vap, ChemPot_vap, 0)
    !ChemPot_vap = ChemPot_vap !* Rmix * T    !Chemical Potentials of the vapor phase [J / mol]
    if((gl%seawater).and.(iphase == 1)) then
        gl%seacalc = .true.
        gl%gecarrier = .true.
        call chempot_num_reac(gl, t, d_fluid1, Chempot_fluid1, iphase)
        gl%seacalc = .false.
        gl%gecarrier = .false.
    elseif((gl%el_present).and.(iphase == 1)) then
        !p = p_calc(gl, t, dw, gl%el%solpos)
        Chempot_fluid1 = chempot_water_brine(gl, t, d_fluid1, p)
    else
        call dna_dni(gl,T, d_fluid1, ChemPot_fluid1, 0)
        ChemPot_fluid1 = Chempot_fluid1 * Rmix * T    !Chemical Potentials of the fluid phase [J / mol]
    end if
    !IF hydrates form the fugacities of the lighter fluid phase are needed as input to the hydrate model
    if (gl%solidtype(2) == 1) call lnf_mix(gl, T, d_fluid1, p, lnf)

    If(twofluidphases) then !In case of 2 fluid phases
        ! Get the density of fluid phase 2
        gl%molfractions = x_fluid2
        call reduced_parameters_calc(gl,T)
        iPhase = 1  !Liquid
        !iphase = 2     !Use this if fluid phases shall be switched --> usually phase1 = vapor or lighter liquid. Needed for example for VHIw phase boundary!
        ! calculate the fluid phase 2 density
        gl%rho_liq2 = rhomix_calc(gl,T, P, rhofluid2_est, iPhase,0)
        if (gl%rho_liq2 < 1.D-12)then
            errval = -8888
            return
        end if
        ! get the gas phase density from the module
        d_fluid2 = gl%rho_liq2
        ! Calculate the chemical potentials for the gas phase
        call R_mix_calc(gl,Rmix)
        !call dna_dni(T, d_vap, ChemPot_vap, 0)
        !ChemPot_vap = ChemPot_vap !* Rmix * T    !Chemical Potentials of the vapor phase [J / mol]
        if((gl%seawater).and.(iphase == 1)) then
            p_save = gl%sea%seap
            gl%seacalc = .true.
            gl%gecarrier = .true.
            gl%sea%seap = p_calc(gl, t, d_fluid2, 1)
            call chempot_num_reac(gl, t, d_fluid2, Chempot_fluid2, iPhase)
            gl%sea%seap = p_save
            gl%seacalc = .false.
            gl%gecarrier = .false.
        elseif((gl%el_present).and.(iphase == 1)) then
            !p = p_calc(gl, t, dw, gl%el%solpos)
            seasave = gl%gecarrier
            gl%gecarrier = .false.
            Chempot_fluid2(1) = chempot_water_brine(gl, t, d_fluid2, p)
            gl%gecarrier = seasave
        else
            call dna_dni(gl,T, d_fluid2, ChemPot_fluid2, 0)
            ChemPot_fluid2 = Chempot_fluid2 * Rmix * T    !Chemical Potentials of the fluid phase [J / mol]
        end if
    End if

    !----------------------------------
    ! get the solid phase properties
    if (gl%solidtype(1) == 2) then        !solidtype(1) = 2 --> Dry Ice
        ChemPot_sol = g_DryIce(gl,T, p)    !Chemical Potential of Dry Ice [J / mol]
    End if
    if (gl%solidtype(1) == 1) then    !solidtype(1) = 1 --> Solid water
        ChemPot_sol = g_WaterIce(gl,T, p)  !Chemical Potential of Water Ice [J / mol]
    End if
    if (gl%solidtype(2)==1) then
        !The fugacity of CO2 in the vapor phase will be used as input for the chem. pot. of water in hydrate (CO2 on position 2 in Hydrate_list)
        fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnf(2:gl%nrofhydrateformers))*1.d6
        do i = 1,gl%nrofhydrateformers-1
            if (is_infinity(fug_gas(i))) then
                errval = -15566
                return
            endif
        enddo
        !Andreas March 2014. Save fugacity of guest molecule for later use
        fug_save = fug_gas
        !Get the chemical potential of water in hydrate
        call hdrt_chem_potent_w(gl,T,p*1.d6,fug_gas,chpw)
        ChemPot_hyd = chpw              !Chemical Potential of Water in Hydrate [J / mol]

        !Get the composition of the hydrate phase
        call hdrt_mole_fract(gl,T,p*1.D6,fug_gas,occup,CiJ,x_hyd2, occup_single, occup_double)
        x_hyd = x_hyd2    !Molfractions of water in hydrate !x_hyd2(2)

    End if

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.

    !Equality of chemical potentials in all phases

    if (twofluidphases) then
        if (gl%solidtype(1) > 0) then
            GibbsEQN(1) = ChemPot_fluid1(gl%solid_pos) - ChemPot_sol
        End if
        if (gl%solidtype(2) == 1) then
            GibbsEQN(1) = ChemPot_fluid1(1) - ChemPot_hyd
        End if
        Do i = 1, gl%ncomp
            GibbsEQN(i+1) = ChemPot_fluid1(i) - ChemPot_fluid2(i)
        End Do
    else
        GibbsEQN(1) = ChemPot_fluid1(1) - ChemPot_hyd
        GibbsEQN(2) = ChemPot_fluid1(gl%solid_pos) - ChemPot_sol
    End if

    !Additional mass balance equations. These equations are necessary in the following cases:
    !1) 3 components and Flash type 1, 2, 3, 4, 5 or 6
    !2) more than 3 components in the mixture
    If ((gl%ncomp > 3) .or. ((gl%ncomp == 3) .and. (iFlash < 7))) then
        !Mass balance
        Do i = 1, gl%ncomp - 1
            if (twofluidphases) then
                if (gl%solidtype(1) > 0) then
                    !if ((iflash /= 3) .and. (iflash /= 6)) then
                    GibbsEQN(gl%ncomp+1+i) = Phasefrac(1) * x_fluid1(i) + Phasefrac(2) * x_sol(i) &
                        & + (1.D0 - Phasefrac(1) - Phasefrac(2)) * x_fluid2(i) - x_known(i)
                    !elseif ((iflash == 3) .or. (iflash == 6)) then
                    !    GibbsEQN(ncomp+1+i) = Phasefrac(3) * x_fluid1(i) + Phasefrac(2) * x_sol(i) &
                    !       & + (1.D0 - Phasefrac(3) - Phasefrac(2)) * x_fluid2(i) - x_known(i)
                    !endif
                End if
                if (gl%solidtype(2) == 1) then
                    !if ((iflash /= 3) .and. (iflash /= 6)) then
                    GibbsEQN(gl%ncomp+1+i) = Phasefrac(1) * x_fluid1(i) + Phasefrac(2) * x_hyd(i) &
                        & + (1.D0 - Phasefrac(1) - Phasefrac(2)) * x_fluid2(i) - x_known(i)
                    !elseif ((iflash == 3) .or. (iflash == 6)) then
                    !    GibbsEQN(ncomp+1+i) = Phasefrac(3) * x_fluid1(i) + Phasefrac(2) * x_hyd(i) &
                    !        & + (1.D0 - Phasefrac(3) - Phasefrac(2)) * x_fluid2(i) - x_known(i)
                    !endif
                End if
            else
                if (iflash /= -1)then
                    GibbsEQN(2+i) = Phasefrac(1) * x_sol(i) + Phasefrac(2) * x_hyd(i) &
                        & + (1.D0 - Phasefrac(1) - Phasefrac(2)) * x_fluid1(i) - x_known(i)
                elseif (iflash == -1) then
                    GibbsEQN(2+i) = (1.d0 - Phasefrac (1) - Phasefrac(2)) * x_sol(i) + Phasefrac(2) * x_hyd(i) &
                        & + Phasefrac(1) * x_fluid1(i) - x_known(i)
                endif
            End if
        End do
    End if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = x_known

    end subroutine
    !**************************************************************************


    !**************************************************************************
    module subroutine Jacobi_solid_NC_3P (gl,P, T, x_fluid1, x_fluid2, x_sol, x_hyd, fug_save, twofluidphases, iFlash, Phasefrac, JacMatrix, errval)
    !**************************************************************************

    !------------------0--------------------------------------------------------
    ! A. J輍er, August 2012, V. Vin?- October 2015







    implicit none

    type(type_gl) :: gl


    double precision:: T, p
    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol, x_hyd
    double precision, dimension(3):: Phasefrac
    integer :: iFlash
    logical :: twofluidphases
    double precision, dimension(60, 60), intent(out):: JacMatrix
    integer :: errval

    integer:: i, j, k
    double precision:: d_fluid1, d_fluid2
    double precision, dimension(30)::  z

    double precision:: rhoredmix_orig, tredmix_orig, Rmix
    double precision:: rho_fluid, dChempot_dT_sol, dChempot_dp_sol
    double precision, dimension(30, 30):: dChempoti_dxj_fluid1, dChempoti_dxj_fluid2
    double precision, dimension(30):: dChempoti_dT_fluid1, dChempoti_dP_fluid1, dChempoti_dT_fluid2,  dChempoti_dP_fluid2
    double precision, dimension(30):: dChempot_dxj_sol

    double precision, dimension(30):: d2nadnidT, dnadni

    !!Variables needed for the pure hydrate model
    !double precision:: fug_CO2, fug_save
    !double precision, dimension(30) :: doccupidfj, CiJ, dxiHdT, dxiHdp!, fug_gas
    !double precision, dimension(30,30):: dxidfj, dxiHdxjFluid
    !double precision, dimension(30):: lnf, dphiidT, dphiidp
    !double precision:: DchpwDT, DchpwDp, DchpwDf
    !double precision, dimension(30,30):: dlnfidXj

    !Variables for the mixed hydrate model
    double precision, dimension(30):: fug_gas, fug_save, DchpwDf !, fug_gas
    double precision, dimension(3,30):: CiJ
    double precision, dimension(3,30,30) :: doccupidfj
    double precision, dimension(30,30):: dxidfj, dxiHdxjFluid
    double precision, dimension(30):: lnf, dphiidT, dphiidp, dChempot_dxj_hyd
    double precision:: DchpwDT, DchpwDp, dChempot_dp_hyd, dChempot_dT_hyd
    double precision, dimension(30,30):: dlnfidXj
    double precision, dimension(30):: dxiHdT, dxiHdp, dmolfracidT, dmolfracidp   ! these quantities are not used and moreover not calculated in a correct way

    JacMatrix = 0.D0
    z = 0.D0
    dChempot_dxj_sol = 0.D0
    dChempot_dxj_hyd = 0.D0
    dChempoti_dT_fluid1 = 0.D0
    dChempoti_dT_fluid2 = 0.D0
    dChempoti_dp_fluid1 = 0.D0
    dChempoti_dp_fluid2 = 0.D0
    dxiHdxjFluid = 0.D0
    errval = 0

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    z = gl%molfractions

    ! get the vapor phase properties
    d_fluid1 = gl%rho_vap
    gl%molfractions = x_fluid1
    call reduced_parameters_calc(gl,T)
    call R_mix_calc(gl,Rmix)
    if ((iFlash == -1) .or. (iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3) .or. (iFlash == 8)) then
        !Calculate the derivative of the chemical potential with respect to T at constant
        !X and p for the first fluid phase
        call d2na_dnidT_P (gl,T, d_fluid1, d2nadnidT, 0)
        call dna_dni(gl,T, d_fluid1, dnadni, 0)
        dChempoti_dT_fluid1 = d2nadnidT * Rmix * T + dnadni * Rmix
        !The following derivatives are needed for the hydrate model
        If(gl%solidtype(2) == 1) then
            call dlnphii_dT(gl,T, d_fluid1, dphiidT)
        End if
    End if
    if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6) .or. (iFlash == 8)) then
        !Calculate the derivative of the chemical potential with respect to p at constant
        !X and T for the vapor phase
        call d2na_dnidp_T (gl,T, d_fluid1, dChempoti_dp_fluid1, 0)
        dChempoti_dp_fluid1 = dChempoti_dp_fluid1 * Rmix * T
        If(gl%solidtype(2) == 1) then
            call dlnphii_dp(gl,T, d_fluid1, dphiidp)
        End if
    End if

    if (iFlash /= 8) then
        !Calculate the derivative of the chemical potential with respect to xj at constant
        !p and T for the liquid phase
        call d2na_dnidxj_PT (gl,T, d_fluid1, dChempoti_dxj_fluid1, 0)
        dChempoti_dxj_fluid1 = dChempoti_dxj_fluid1 * Rmix * T
    end if

    !The following variables are needed for the hydrate model
    If(gl%solidtype(2) == 1) then
        !Andreas March 2013: Take the saved fugacity for the guest
        call lnf_mix(gl,T, d_fluid1, p, lnf)
        call dlnfi_dxj_TP (gl,T, d_fluid1, dlnfidXj)
    End if

    If(twofluidphases) then !In case of 2 fluid phases
        ! get the liquid phase 2 properties
        d_fluid2 = gl%rho_liq2
        gl%molfractions = x_fluid2
        call reduced_parameters_calc(gl,T)
        call R_mix_calc(gl,Rmix)
        if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3) .or. (iFlash == 8)) then
            !Calculate the derivative of the chemical potential with respect to T at constant
            !X and p for the first fluid phase
            call d2na_dnidT_P (gl,T, d_fluid2, d2nadnidT, 0)
            call dna_dni(gl,T, d_fluid2, dnadni, 0)
            dChempoti_dT_fluid2 = d2nadnidT * Rmix * T + dnadni * Rmix
        End if
        if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6) .or. (iFlash == 8)) then
            !Calculate the derivative of the chemical potential with respect to p at constant
            !X and T for the vapor phase
            call d2na_dnidp_T (gl,T, d_fluid2, dChempoti_dp_fluid2, 0)
            dChempoti_dp_fluid2 = dChempoti_dp_fluid2 * Rmix * T
        End if
        if (iFlash /= 8) then
            !Calculate the derivative of the chemical potential with respect to xj at constant
            !p and T for the liquid phase
            call d2na_dnidxj_PT (gl,T, d_fluid2, dChempoti_dxj_fluid2, 0)
            dChempoti_dxj_fluid2 = dChempoti_dxj_fluid2 * Rmix * T
        end if
    End if

    !Get the properties and derivatives of the solid phases
    if (gl%solidtype(1) > 0) then      !Pure solid substance
        if ((iFlash == -1) .or. (iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3) .or. (iFlash == 8)) then
            if (gl%solidtype(1) == 1) then
                dChempot_dT_sol = dgdT_WaterIce(gl,T, p)   !Derivative of the Chemical Potential of Dry Ice [J / mol]
            End if
            if (gl%solidtype(1) == 2) then
                dChempot_dT_sol = dgdT_DryIce(gl,T, p)     !Derivative of the Chemical Potential of Water Ice [J / mol]
            End if
        End if
        if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6) .or. (iFlash == 8)) then
            if (gl%solidtype(1) == 1) then
                dChempot_dp_sol = dgdp_WaterIce(gl,T, p)  !Derivative of the Chemical Potential of Water Ice [J / mol]
            End if
            if (gl%solidtype(1) == 2) then
                dChempot_dp_sol = dgdp_DryIce(gl,T, p)    !Derivative of the Chemical Potential of Dry Ice [J / mol]
            End if
        End if
        dChempot_dxj_sol = 0.D0     !pure solid phase has a fixed composition
    end if
    If (gl%solidtype(2) == 1) then      !Hydrate
        !The fugacity of CO2 in the vapor phase will be used as input for the chem. pot. of water in hydrate (mapping(2) = CO2, see Hydrate_list module variable)
        !fug_gas(1:nrofhydrateformers-1) = dexp(lnf(2:nrofhydrateformers))*1.d6     !Multi component hydrates
        !fug_CO2 = dexp(lnf(2)) * 1.D6
        !Andreas March 2013: Take the saved fugacity for the guest

        fug_gas = fug_save  ! 2015 - fugacitites of guests already calculated in 'SysOfEqs_solid_NC_3P'
        do i = 1,gl%nrofhydrateformers-1
            if (is_infinity(fug_gas(i))) then
                errval = -15566
                return
            endif
        enddo
        call hdrt_Df_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDf)
        call hdrt_dxi_dfj(gl,T,p*1.D6,fug_gas,doccupidfj,CiJ,dxidfj)


        if ((iFlash == -1) .or. (iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then
            !The derivative with respect to Temperature is calculated according to: dcp_dT|p,x_CO2_vap = dcp_dT|p,f_CO2 + dcp_df|T,p * dfi_dT|p, x_CO2
            call hdrt_DT_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDT)
            ! pure hydrate: dChempot_dT_hyd = DchpwDT + DchpwDf * dphiidT(2)*dexp(lnf(2))*1.d6
            ! mixed hydrate 2015:
            dChempot_dT_hyd = DchpwDT
            Do i = 2, gl%nrofhydrateformers
                dChempot_dT_hyd = dChempot_dT_hyd + DchpwDf(i) * dphiidT(i)*dexp(lnf(i))*1.d6    ! Eq. 5.173 in Andreas' diss.
            End do

            ! 2015 - dxiHdT(i) is not used, moreover its calculation is wrong - see Eq. 5.177 in Andreas' dissertation
            ! 2017 - dxiHdT(i) is calculated numerically
            !!Derivative of the molfraction of component i in hydrate wrt. T
            !!If the component is not present in the hydrate phase, the derivative remains 0
            call hdrt_Dxi_DT(gl,T,p*1.d6,fug_gas,dmolfracidT)
            Do i = 1, gl%nrofhydrateformers
                Do k = 2, gl%nrofhydrateformers
                    dxiHdT(i) = dmolfracidT(i) + dxidfj(k,i)*dphiidT(k)*dexp(lnf(k))*1.d6
                End do
            End do

        End if
        if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
            !The derivative with respect to pressure is calculated according to: dcp_dp|T,x_CO2_vap = dcp_dp|p,f_CO2 + dcp_df|T,p * dfi_dp|p, x_CO2
            call hdrt_Dp_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDp)
            ! pure hydrate: dChempot_dp_hyd = DchpwDp + DchpwDf * dexp(lnf(2))*1.d6 * (dphiidp(2) + 1.D0 / (p*1.D6))
            ! mixed hydrate 2015:
            dChempot_dp_hyd = DchpwDp
            Do i = 2, gl%nrofhydrateformers
                dChempot_dp_hyd = dChempot_dp_hyd + DchpwDf(i) * dexp(lnf(i))*1.d6 * (dphiidp(i) + 1.D0 / (p*1.D6))    ! Eq. 5.175 in diss.
            End do

            ! 2015 - dxiHdp(i) is not used, moreover its calculation is wrong - see Eq. 5.178 in Andreas' dissertation
            !!Derivative of the molfraction of component i in hydrate wrt. p
            !!If the component is not present in the hydrate phase, the derivative remains 0
            !Do i = 1, nrofhydrateformers
            !    Do k = 2, nrofhydrateformers
            !        dxiHdp(i) = dxidfj(k-1,i)*(dphiidp(k) + 1.D0 / (p*1.D6))*dexp(lnf(k))*1.d6
            !    End do
            !End do
        End if

        !The derivative with respect to x_CO2 is calculated according to: dcp_dx_CO2|p,T = dcp_df|T,p * dfi_dx_CO2|T,p
        dChempot_dxj_hyd = 0.d0
        Do j = 1, gl%ncomp - 1 ! x(j)
            Do i = 2, gl%nrofhydrateformers
                ! pure hydrate: dChempot_dxj_hyd(j) = DchpwDf * dlnfidxj(j, i)*dexp(lnf(i))*1.d6       !Only one guest
                !dChempot_dxj_sol(j) = DchpwDf(i) * dlnfidxj(2, 2)*dexp(lnf(2))*1.d6   !Multi component hydrates ... old comment - Andreas

                ! mixed hydrate 2015:
                dChempot_dxj_hyd(j) = dChempot_dxj_hyd(j) + DchpwDf(i) * dlnfidxj(j, i)*dexp(lnf(i))*1.d6    ! Eq. 5.176 in Andreas' diss.
            End do
        End Do

        !Derivative of the molfraction of component i in hydrate wrt. the molfraction j in the fluid
        !If the component is not present in the hydrate phase, the derivative remains 0
        Do j = 1, gl%ncomp - 1
            Do i = 1, gl%nrofhydrateformers
                ! pure hydrate: Do k = 2, nrofhydrateformers
                !    dxiHdxjFluid(j,i) = dxidfj(k,i)*dlnfidxj(j, k)*dexp(lnf(i))*1.d6
                !End do

                ! mixed hydrate 2015:
                dxiHdxjFluid(j,i) = 0.d0
                Do k = 2, gl%nrofhydrateformers    ! k-guest fugacity
                    dxiHdxjFluid(j,i) = dxiHdxjFluid(j,i) + dxidfj(k,i) * dlnfidxj(j, k)*dexp(lnf(k))*1.d6    ! Eq. 5.179 in Andreas' diss.
                End do
            End Do
        End Do

    End if

    !TRIPLE POINT OF PURE SUBSTANCE
    if (iflash == 8) then
        !SysofEqs:
        !GibbsEQN(1) = chempot_vap - chempot_sol
        !GibbsEQN(2) = chempot_vap - chempot_liq
        !First unknown: T, second unknown: p
        JacMatrix(1,1) = dChempoti_dT_fluid1(1) - dChempot_dT_sol
        JacMatrix(1,2) = dChempoti_dT_fluid1(1) - dChempoti_dT_fluid2(1)
        JacMatrix(2,1) = dChempoti_dp_fluid1(1) - dChempot_dp_sol
        JacMatrix(2,2) = dChempoti_dp_fluid1(1) - dChempoti_dp_fluid2(1)
        return
    end if

    !The Matrix is split into different sections (compare to notes)
    !A,B,C ..... M
    !The structure of the matrix is based on the VLL algorithm. In all cases treated here, the sections B,D,J of the VLL matrix do not exist
    !In case of 2 solid phases, the sections B,C,D,E,H,K,J do not exist
    !Furthermore the special case of bubble or dew points for three component mixtures must be considered

    !Two different cases of matrix structures might occur
    !Case 1: VLH, VLS, LLH, LLS     --> 2N equations, 2N unknowns (both fluid phase compositions, 2 out of T, p, phasefrac
    !Case 2: VHS, LHS               --> N+1, N+1 unknowns (one fluid phase composition, 2 out of T, p, phasefrac
    !Case 1:
    if (twofluidphases) then
        Do j = 1, gl%ncomp -1
            !A1     !This section contains only one equation for all case
            if (gl%solidtype(2) == 1) then
                !If hydrates form, the first equation is always equality of chem. pot of water in hydrate phase equals chem pot in the first fluid phase
                JacMatrix(j,1) = dChempoti_dxj_fluid1(j, 1) - dChempot_dxj_hyd(j)
            else
                !If no hydrates form, chem.pot for solid comp. must be equal in solid phase and first fluid phase
                JacMatrix(j,1) = dChempoti_dxj_fluid1(j, gl%solid_pos) - dChempot_dxj_sol(j)
            end if
            !C
            JacMatrix(gl%ncomp-1+j, 1) = 0.D0
            Do i = 1, gl%ncomp
                !A2
                JacMatrix(j,1+i) = dChempoti_dxj_fluid1(j, i)
                !E
                JacMatrix(gl%ncomp-1+j, 1+i) = -dChempoti_dxj_fluid2(j, i)
            End do
        End do

        !In case of "Sublimation" or "Melting" or "Dew" point calculations the derivatives of the fugacities wrt T or p are needed
        if(iFlash < 7) then
            Do i = 1, gl%ncomp
                if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then
                    !G
                    !!UNSCH烰
                    if (gl%solidtype(2) == 1) then
                        JacMatrix(2*gl%ncomp-1,1) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd
                    else
                        JacMatrix(2*gl%ncomp-1,1) = dChempoti_dT_fluid1(gl%solid_pos) - dChempot_dT_sol
                    End if
                    !H
                    JacMatrix(2*gl%ncomp-1,i+1) = dChempoti_dT_fluid1(i) - dChempoti_dT_fluid2(i)
                end if
                if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
                    !G
                    !!UNSCH烰
                    if (gl%solidtype(2) == 1) then
                        JacMatrix(2*gl%ncomp-1,1) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd
                    else
                        JacMatrix(2*gl%ncomp-1,1) = dChempoti_dp_fluid1(gl%solid_pos) - dChempot_dp_sol
                    End if
                    !H
                    JacMatrix(2*gl%ncomp-1,i+1) = dChempoti_dp_fluid1(i) - dChempoti_dp_fluid2(i)
                end if
            End do
        End if

        !If ncomp > 3 or ncomp == 3 and iflash < 7 the mass balance needs to be solved
        If ((gl%ncomp > 3) .or. ((gl%ncomp == 3) .and. (iFlash < 7))) then
            Do j = 1, gl%ncomp -1
                Do i = 1, gl%ncomp-1
                    if (iFlash == 7) then
                        !G
                        JacMatrix(2*gl%ncomp-1,1) = 0.D0
                        !H
                        JacMatrix(2*gl%ncomp-1,i+1) = 0.D0
                    End if
                    !F
                    !Sehr UNSCH烰!!!!
                    JacMatrix(2*gl%ncomp,1) = 0.D0
                    JacMatrix(2*gl%ncomp,1+i) = 0.D0
                    !I

                    if (gl%solidtype(2) == 1) then
                        if (i==j) then
                            !if ((iflash /= 3) .and. (iflash /= 6)) then
                            JacMatrix(j,gl%ncomp+1+i) = Phasefrac(1)+Phasefrac(2)*dxiHdxjFluid(j,i)   !In case of a sublimation point, this is 0
                            !elseif ((iflash == 3) .or. (iflash == 6)) then
                            !    JacMatrix(j,ncomp+1+i) = Phasefrac(3)+Phasefrac(2)*dxiHdxjFluid(j,i)
                            !endif
                        else
                            JacMatrix(j,gl%ncomp+1+i) = Phasefrac(2)*dxiHdxjFluid(j,i)
                        end if
                    else
                        if (i==j) then
                            !if ((iflash /= 3) .and. (iflash /= 6)) then
                            JacMatrix(j,gl%ncomp+1+i) = Phasefrac(1)   !In case of a sublimation point, this is 0
                            !elseif ((iflash == 3) .or. (iflash == 6)) then
                            !    JacMatrix(j,ncomp+1+i) = Phasefrac(3)
                            !endif
                        else
                            JacMatrix(j,gl%ncomp+1+i) = 0.D0
                        end if
                    End if

                    !K
                    if (i==j) then
                        !if ((iflash /= 3) .and. (iflash /= 6)) then
                        JacMatrix(j+gl%ncomp-1,gl%ncomp+1+i) = 1.D0 - Phasefrac(1) - Phasefrac(2)
                        !elseif ((iflash == 3) .and. (iflash == 6)) then
                        !    JacMatrix(j+ncomp-1,ncomp+1+i) = 1.D0 - Phasefrac(3) - Phasefrac(2)
                        !endif
                    else
                        JacMatrix(j+gl%ncomp-1,gl%ncomp+1+i) = 0.D0
                    end if
                    !L
                    if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then
                        JacMatrix(2*gl%ncomp-1,gl%ncomp+1+i) = 0.D0 !phasefrac(2)*dxiHdT
                    End if
                    if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
                        JacMatrix(2*gl%ncomp-1,gl%ncomp+1+i) = 0.D0 !phasefrac(2)*dxiHdp
                    End if
                    if (iFlash == 7) then
                        JacMatrix(2*gl%ncomp-1,gl%ncomp+1+i) = x_fluid1(i) - x_fluid2(i)
                    end if
                    !M
                    if ((iFlash == 1) .or. (iFlash == 4).or. (iFlash == 7)) then  !"sublimation point" phasefrac(1) (beta_V) = 0 (Same equation needed for iFlash = 5)
                        if (gl%solidtype(2) == 1) then
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_hyd(i) - x_fluid2(i)
                        else
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_sol(i) - x_fluid2(i)
                        end if
                    end if
                    if ((iFlash == 2) .or. (iFlash == 5)) then  !"Melting point" phasefrac(2) (gamma or ebsilon) = 0
                        if (gl%solidtype(2) == 1) then
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_fluid1(i) - x_fluid2(i)
                        else
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_fluid1(i) - x_fluid2(i)
                        end if
                    end if
                    if ((iflash == 3) .or. (iflash == 6)) then !"dew point" phasefrac(3) (beta_L2) = 0
                        if (gl%solidtype(2) == 1) then
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_hyd(i) - x_fluid1(i)
                        else
                            JacMatrix(2*gl%ncomp,gl%ncomp+1+i) = x_sol(i) - x_fluid1(i)
                        end if
                    endif
                End do
            End do
        End if

    else    !Case 2

        Do j = 1, gl%ncomp -1
            !A1
            JacMatrix(j,1) = dChempoti_dxj_fluid1(j, 1) - dChempot_dxj_hyd(j)
            !A2
            JacMatrix(j,2) = dChempoti_dxj_fluid1(j, gl%solid_pos)
        End do

        !In case of "Sublimation" or "Melting" point calculations the derivatives of the fugacities wrt T or p are needed
        if(iFlash < 7) then
            if ((iFlash == -1) .or. (iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then
                !G
                JacMatrix(gl%ncomp,1) = dChempoti_dT_fluid1(1) - dChempot_dT_hyd
                !H
                JacMatrix(gl%ncomp,2) = dChempoti_dT_fluid1(gl%solid_pos) - dChempot_dT_sol
            end if

            if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
                !G
                JacMatrix(gl%ncomp,1) = dChempoti_dp_fluid1(1) - dChempot_dp_hyd
                !H
                JacMatrix(gl%ncomp,2) = dChempoti_dp_fluid1(gl%solid_pos) - dChempot_dp_sol
            end if
        End if

        !If ncomp > 3 or ncomp == 3 and iflash < 5 the mass balance needs to be solved
        If ((gl%ncomp > 3) .or. ((gl%ncomp == 3) .and. (iFlash < 7))) then
            !F
            JacMatrix(gl%ncomp+1,1) = 0.D0
            JacMatrix(gl%ncomp+1,2) = 0.D0
            if (iFlash == 7) then
                !G
                JacMatrix(gl%ncomp-1,1) = 0.D0
                !H
                JacMatrix(gl%ncomp-1,2) = 0.D0
            End if
            Do j = 1, gl%ncomp -1
                Do i = 1, gl%ncomp-1
                    !I
                    if (iFlash /= -1 ) then
                        if (i==j) then
                            !JacMatrix(j,2+i) = Phasefrac(1)+Phasefrac(2)*dxiHdxjFluid(j,i)   !In case of a sublimation point, this is 0
                            JacMatrix(j,2+i) = (1.d0 - Phasefrac(1) - Phasefrac(2)) + Phasefrac(2)*dxiHdxjFluid(j,i)
                        else
                            JacMatrix(j,2+i) = Phasefrac(2)*dxiHdxjFluid(j,i)
                        end if
                    elseif (iFlash == -1) then
                        if (i==j) then
                            JacMatrix(j,2+i) = phasefrac(1) + Phasefrac(2) *dxiHdxjFluid(j,i)
                        else
                            JacMatrix(j,2+i) = Phasefrac(2)*dxiHdxjFluid(j,i)
                        end if
                    endif
                    !L
                    if ((iFlash == 1) .or. (iFlash == 2) .or. (iFlash == 3)) then
                        JacMatrix(gl%ncomp,2+i) = 0.D0 !phasefrac(2)*dxiHdT(j,i)
                    elseif (iFlash == -1) then
                        JacMatrix(gl%ncomp,2+i) = phasefrac(2)*dxiHdT(i)
                    End if
                    if ((iFlash == 4) .or. (iFlash == 5) .or. (iFlash == 6)) then
                        JacMatrix(gl%ncomp,2+i) = 0.D0 !phasefrac(2)*dxiHdp(j,i)
                    End if
                    if (iFlash == 7) then
                        JacMatrix(gl%ncomp,2+i) = x_hyd(i) - x_fluid1(i)
                    end if
                    !M
                    if ((iFlash == 1) .or. (iFlash == 4).or. (iFlash == 7)) then  !"Melting point" phasefrac(1) (ebsilon = 0) (Same equation needed for iFlash = 7)
                        !JacMatrix(ncomp+1,2+i) = x_sol(i) - x_fluid1(i)
                        JacMatrix(gl%ncomp+1,2+i) = x_hyd(i) - x_fluid1(i)
                    end if
                    if ((iFlash == 2) .or. (iFlash == 5)) then  !"Freezing point" phasefrac(2) (gamma) = 0
                        !JacMatrix(ncomp+1,2+i) = x_hyd(i) - x_fluid1(i)
                        JacMatrix(gl%ncomp+1,2+i) = x_sol(i) - x_fluid1(i)
                    end if
                    if (iFlash == -1) then
                        JacMatrix(gl%ncomp+1,2+i) = x_hyd(i) - x_sol(i)
                    endif
                    !if ((iflash == 3) .or. (iflash == 6)) then   !"dew point" phasefrac(3) (beta_L2) = 0
                    !    JacMatrix(ncomp+1,2+i) = x_hyd(i) - x_fluid2(i)
                    !endif
                End do
            End do
        End if
    End if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z


    end subroutine
    !**************************************************************************


    !**************************************************************************
    module subroutine ptflash_solid_NC_2P (gl,press, Temp, rho, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, &
        & iPhase, errval, iter)

    !**************************************************************************
    ! SUBROUTINE FOR THE ITERATIVE CALCULATION OF THE S-L and S-V EQUILIBRIUM
    ! THIS ROUTINE IS NOT LIMITED IN THE NUMBER OF COMPONENTS, BUT ONLY WORKS FOR
    ! 2 PHASES IN EQUILIBRIUM (OF WHICH ONE IS A SOLID PHASE)
    !
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS MAY BE PERFORMED:
    !   - SUBLIMATION / MELTING POINT: P AND x' or x'' VEXTOR GIVEN   --  iFlash = 1    (May not work for hydrates)
    !   - SUBLIMATION / MELTING POINT: T AND x' or x'' VEXTOR GIVEN   --  iFlash = 2    (May not work for hydrates)
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN --  iFlash = 3
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press           - Pressure [MPa]
    !   Temp            - Temperature [K]
    !   x_solid         - Composition of the solid phase (0 except for the position solidnr)
    !   x_hyd           - Composition of the hydrate phase
    !   x_fluid         - Composition of the fluid phase (vapor or liquid)
    !   rhofluid_est    - Estimated fluid phase density
    !   beta            - Molar "fluid" fraction (xi - xi_sol) / (xi_fluid - xi_sol)    0 --> solid phase only; 1 --> fluid phase only
    !   iFlash          - Flash mode
    !   iphase          - Specifies which fluid phase is assumed:       0 - Let the density solver choose based on the gibbs energy
    !                                                                   1 - Assume liquid
    !                                                                   2 - Assume vapor
    ! OUTPUT:
    !   rho             - densities of all phases found are given back according to the order of the phase eq. type: (e.g. VS --> rho(1) = Vapor density, rho(2) = solid density
    !                   - in principle only a two dimensional vector is needed, however, corresponding to the 3 phase routine "ptflash_solid_NC_3P" a three dimensional vector is used
    !   errval          - Error value

    !------------------0--------------------------------------------------------
    ! A. J輍er May 2012; V.Vin?November 2015







    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, beta_loc
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_fluid, x_solid, x_hyd! oder doch kein intent?!
    double precision, dimension(30):: z, x_fluid_new!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: rhofluid_est
    integer:: iFlash
    integer:: errval, iter
    integer:: iPhase
    !New output variable for phase densities, Andreas March 2014
    double precision, dimension(3) :: rho

    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(60):: GibbsEQN, Delta_X
    double precision:: sum_fluid, stepsize, Temp_new, press_new
    integer:: i, j, k, eqn
    integer, dimension(1):: maxID
    character(255):: herr

    !Variables for Solid H2O - Hydrate Equilibrium
    double precision:: ChemPot_sol

    double precision :: rho_H     ! , fug_save

    !!Variables for hydrates
    !double precision, dimension (2):: occup, Cij, x_sol
    !double precision:: fug_CO2, chpw, DchpwDf, fug_CO2_new
    !double precision, dimension(30):: lnf

    !Variables for the mixed hydrates & pure solids
    double precision :: chpw
    double precision, dimension(30):: fug_gas, fug_gas_new, x_sol, DchpwDf, lnf, fug_save
    double precision, dimension(3,30):: occup, CiJ, occup_single, occup_double


    z = x_known
    errval = 0
    Delta_X = 1.D0
    x_fluid_new = 0.D0
    Temp_new = 0.D0
    press_new = 0.D0
    stepsize = 1.D0
    rho_H = 0.D0
    fug_save = 0.D0

    iter = 1
    GibbsEQN = 1.D10
    x_sol = 0.D0

    if (iFlash < 3) then
        x_fluid = x_known
    end if

    !Solid water and hydrate in Equilibrium (ONLY TWO COMPONENTS POSSIBLE). The fugacity of the guest is the unknown!!
    if ((gl%ncomp == 2) .And. (gl%solidtype(1) == 1) .AND. (gl%solidtype(2) == 1)) then

        ChemPot_sol = g_WaterIce(gl,Temp, press)
        fug_gas = 1.D7      !Arbitrary Startvalue

        do i = 1, 30
            do k = 1,gl%nrofhydrateformers-1
                if (is_infinity(fug_gas(k))) then
                    errval = -15566
                    return
                endif
            enddo
            call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_gas,chpw)
            GibbsEQN(1) = chpw - ChemPot_sol

            !Break criterion
            if (dabs(GibbsEQN(1)) < 1.D-8) then
                rho(1) = 1.D0 / v_WaterIce(gl,Temp, press)
                !call hdrt_density(temp,press,fug_save,occup,CiJ,rho_H)
                call hdrt_density(gl,temp,press,fug_gas,occup,CiJ,rho_H)   ! mixed hydrates
                rho(2) = rho_H
                beta_loc = (z(1) - x_hyd(1))/(x_solid(1) - x_hyd(1))
                return
            end if

            call hdrt_Df_chem_potent_w(gl,Temp,press*1.D6,fug_gas,DchpwDf)
            JacMatrix(1,1) = DchpwDf(1)     ! Only one guest = pure hydrate

            Delta_X = - GibbsEQn
            eqn = 1
            call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
            if (errval /= 0) return

            ! Only one guest = pure hydrate
            do j = 1, 20
                fug_gas_new(1) = fug_gas(1) + Delta_X(1)
                if (fug_gas_new(1) < 0) then
                    Delta_X(1) = Delta_X(1) / 2.D0
                else
                    exit
                end if
                if ((j == 20) .and. (fug_gas_new(1) < 0)) errval = -2222
            end do

            fug_gas = fug_gas_new

            !Get the composition of the hydrate phase
            call hdrt_mole_fract(gl,temp,press*1.d6,fug_gas,occup,CiJ,x_sol, occup_single, occup_double)
            x_hyd(2) = x_sol(2)  !Molfractions of gas in hydrate
            x_hyd(1) = x_sol(1)  !Molfraction of water in hydrate

        end do

    elseif ((gl%ncomp > 2) .And. (gl%solidtype(1) == 1) .AND. (gl%solidtype(2) == 1)) then
        errval = -15572
    else

        do i = 1, 100
            call SysOfEqs_solid_NC_2P(gl,press, Temp, x_fluid, x_solid, x_hyd, rhofluid_est, beta_loc, fug_save, iFlash, iPhase, &
                & GibbsEQN, errval)

            if (errval /= 0) return

            !Recalculate phase fraction
            !Andreas, July 2014
            if ((gl%solidtype(1) == 0) .and. (gl%solidtype(2) == 1)) then
                beta_loc = (z(1) - x_hyd(1))/(x_fluid(1) - x_hyd(1))
            elseif ((gl%solidtype(1) > 0) .and. (gl%solidtype(2) == 0)) then
                beta_loc = (z(1) - x_solid(1))/(x_fluid(1) - x_solid(1))
            end if

            !Break Condition
            ! This is the exit condition for the iteration!
            if (maxval(dabs(GibbsEQN)) < 1.D-8) then

                rho(1) = gl%rho_vap    !as HIw equilibrium is handled above, the first phase will always be a fluid phase and the density is stored in the module variable "rho_vap"
                if (gl%solidtype(1) == 1) then         !solid H2O
                    rho(2) = 1.D0 / v_WaterIce(gl,Temp, press)
                elseif (gl%solidtype(1) == 2) then     !solid CO2
                    rho(2) = 1.D0 / v_DryIce(gl,Temp, press)
                elseif (gl%solidtype(2) == 1) then     !Hydrate
                    call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                    rho(2) = rho_H
                end if

                return
            end if

            call Jacobi_solid_NC_2P(gl,press, temp, x_fluid, x_solid, x_hyd, beta_loc, fug_save, iFlash, JacMatrix, errval)
            Delta_X = - GibbsEQn

            if (iFLash == 3) then
                eqn = (gl%ncomp - 1)
            else
                eqn = 1
            end if

            call LUdecomp(gl,eqn,JacMatrix,Delta_X,errval,herr)
            if (errval /= 0) return

            ! return the new set of variables
            sum_fluid = 0.D0

            if (maxval(dabs(Delta_X)) < 1.D-12) then

                rho(1) = gl%rho_vap    !as HIw equilibrium is handled above, the first phase will always be a fluid phase and the density is stored in the module variable "rho_vap"
                if (gl%solidtype(1) == 1) then         !solid H2O
                    rho(2) = 1.D0 / v_WaterIce(gl,Temp, press)
                elseif (gl%solidtype(1) == 2) then     !solid CO2
                    rho(2) = 1.D0 / v_DryIce(gl,Temp, press)
                elseif (gl%solidtype(2) == 1) then     !Hydrate
                    call hdrt_density(gl,temp,press,fug_save,occup,CiJ,rho_H)
                    rho(2) = rho_H
                end if

                return
            end if

            iter = i

            select case (iFlash)
                ! --------------------------------------------------------------------------------
            case(1)
                press_new = press
                x_fluid_new = x_fluid
                Temp_new = Temp + Delta_X(1)
                ! --------------------------------------------------------------------------------
            case(2)
                Temp_new = Temp
                x_fluid_new = x_fluid
                press_new = press + Delta_X(1) / 1.D6
                ! --------------------------------------------------------------------------------
            case(3)
                j = 1
                Temp_new = Temp
                press_new = press
                do while (j < gl%ncomp+1)
                    if (j < gl%ncomp) then
                        x_fluid_new(j) = x_fluid(j) + Delta_X(j)
                    else
                        x_fluid_new(gl%ncomp) = 1.D0 - sum_fluid
                    End if
                    sum_fluid = sum_fluid + x_fluid_new(j)
                    ! check if the stepsize is too large.
                    ! The 1st condition: (x(j)_vap - z(j))*(x(j)_liq - z(j)) < 0  ( z is inbetween xliq and xvap)
                    ! The 2nd condition: x(j)_vap AND x(j)_liq > 0
                    ! have to be fulfilled for all components
                    ! Otherwise the stepsize is reduced
                    if ((x_fluid_new(j) > 1) .Or. (x_fluid_new(j) < 0)) then
                        stepsize = 2.D0
                        Delta_X = Delta_X/stepsize
                        if (maxval(dabs(Delta_X)) < 1.D-10) then
                            errval = -3333
                            return
                        end if
                        stepsize = 1.D0
                        j = 0
                        sum_fluid = 0.D0
                    end if
                    j = j + 1
                end do
                !Andreas July 2014
                !if ((solidtype(1) == 0) .and. (solidtype(2) == 1)) then
                !    beta = (z(1) - x_hyd(1))/(x_fluid_new(1) - x_hyd(1))
                !elseif ((solidtype(1) > 0) .and. (solidtype(2) == 0)) then
                !    beta = (z(1) - x_solid(1))/(x_fluid_new(1) - x_solid(1))
                !end if
                ! --------------------------------------------------------------------------------
            end select

            ! write the new values to the variables
            x_fluid = x_fluid_new
            Temp = Temp_new
            press = press_new

        end do

    end if

    ! Iteration failed!
    if (i > 30) errval = -2222

    end subroutine

    !**************************************************************************

    module subroutine SysOfEqs_solid_NC_2P (gl,P, T, x_fluid, x_solid, x_hyd, rhofluid_est, beta_loc, fug_save, &
        & iFlash, iPhaseIn, GibbsEQN, errval)
    !**************************************************************************
    ! SUBROUTINE FOR SETTING UP THE SYSTEM OF EQUATIONS FOR PERFORMING PHASE
    ! EQUILIBRIUM CALCULATIONS.
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
    !   - SUBLIMATION / MELTING POINT: P AND x' or x'' VEXTOR GIVEN   --  iFlash = 1
    !   - SUBLIMATION / MELTING POINT: T AND x' or x'' VEXTOR GIVEN   --  iFlash = 2
    !   - PT-FLASH:     P, T AND x VECTOR GIVEN                       --  iFlash = 3
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   P               - Pressure
    !   T               - Temperature
    !   x_fluid         - Fluid phase composition
    !   x_solid         - Composition of the solid phase (0 except for the position solidnr)
    !   x_hyd           - Composition of the hydrate phase
    !   rhofluid_est    - Estimated vapor phase density
    !   beta            - Molar "fluid" fraction (xi - xi_sol) / (xi_fluid - xi_sol)    0 --> solid phase only; 1 --> fluid phase only
    !   iFlash          - Flash mode
    !   iphase          - Specifies which fluid phase is assumed:       0 - Let the density solver choose based on the gibbs energy
    !                                                                   1 - Assume liquid
    !                                                                   2 - Assume vapor
    ! OUTPUT:
    !   errval          - Error value
    !   GibbsEQN        - 60 x 1 matrix containing the set of equations for the Gibbs minimization algorithm
    !--------------------------------------------------------------------------
    ! A. J輍er,  May 2012






    implicit none

    type(type_gl) :: gl


    double precision:: P, T, beta_loc
    double precision, dimension(30):: x_fluid, x_solid
    double precision, dimension(30):: x_hyd
    double precision, dimension(60):: GibbsEQN
    double precision:: rhofluid_est
    integer:: iFlash
    integer:: errval
    integer:: iPhaseIn

    double precision, dimension(30):: z, ChemPot_fluid, lnf!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision::ChemPot_sol, ChemPot_hyd        !, chpw, fug_CO2, fug_save
    double precision:: rhoredmix_orig, tredmix_orig, d_fluid, dbeta_dxi
    double precision:: rhofluid, Rmix, dxk_dxj_sol,dxk_dxj_fluid
    integer:: i, j, h, errorflag, Iphase
    !!Variables for hydrates
    !double precision, dimension (2):: occup, Cij, x_sol

    !Variables for the mixed hydrates
    double precision:: chpw
    double precision, dimension(30):: x_sol, fug_gas, fug_save
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double


    errval = 0
    GibbsEQN = 0.D0
    z = 0.D0
    chpw = 0.D0
    x_sol = 0.D0
    fug_gas = 0.D0
    fug_save = 0.D0
    CiJ = 0.D0
    occup = 0.D0
    occup_single = 0.D0
    occup_double = 0.D0

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    z = gl%molfractions
    gl%molfractions = x_fluid
    call reduced_parameters_calc(gl,T)

    !The assumed phase is now passed by the calling routine
    !iPhase = 0
    iphase = iPhaseIn
    rhofluid = rhomix_calc(gl,T, P, rhofluid_est, iPhase,0)

    !Preliminary solution, save the density to BOTH global variables
    gl%rho_vap = rhofluid
    !rho_liq = rhofluid

    if (rhofluid < 1.D-12) then
        errval = -8888
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) rhofluid
        write (*,*) 'error ', errval, ' -- mixture density iteration failed'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        return
    end if

    !----------------------------------
    ! get the fluid phase properties

    if(gl%el_present) then
        Chempot_fluid = chempot_Water_Brine(gl, t, rhofluid, p)
    else
        call R_mix_calc(gl,Rmix)
        call dna_dni(gl,T, rhofluid, ChemPot_fluid, 0)
        ChemPot_fluid = Chempot_fluid * Rmix * T    !Chemical Potentials of the fluid phase [J / mol]
    end if
    !Get the fugacities of the vapor phase (THIS IS NEEDED FOR THE HYDRATE MODEL)
    if (gl%solidtype(2) == 1) call lnf_mix(gl, T, rhofluid, p, lnf)

    !----------------------------------
    ! get the solid phase properties
    if (gl%solidtype(1) == 2) then        !solidtype(1) = 2 --> Dry Ice
        ChemPot_sol = g_DryIce(gl,T, p)    !Chemical Potential of Dry Ice [J / mol]
    End if
    if (gl%solidtype(1) == 1) then    !solidtype(1) = 1 --> Solid water
        ChemPot_sol = g_WaterIce(gl,T, p)  !Chemical Potential of Water Ice [J / mol]
    End if
    if (gl%solidtype(2)==1) then
        !The fugacity of CO2 in the vapor phase will be used as input for the chem. pot. of water in hydrate (CO2 on position 2 in Hydrate_list)
        fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnf(2:gl%nrofhydrateformers))*1.d6     !Multi component hydrates
        !Andreas, March 2014. Save the guest fugacity for later use
        do i = 1,gl%nrofhydrateformers-1
            if (is_infinity(fug_gas(i))) then
                errval = -15566
                return
            endif
        enddo
        fug_save = fug_gas
        !Get the chemical potential of water in hydrate
        call hdrt_chem_potent_w(gl,T,p*1.d6,fug_gas,chpw)
        !call hdrt_chem_potent_w(T,p*1.d6,fug_gas,chpw)                             !Multi component hydrates
        ChemPot_hyd = chpw              !Chemical Potential of Water in Hydrate [J / mol]

        !Get the composition of the hydrate phase
        call hdrt_mole_fract(gl,T,p*1.D6,fug_gas,occup,CiJ,x_sol, occup_single, occup_double)
        x_hyd = x_sol !Molfractions of water in hydrate

    End if

    !----------------------------------
    ! Setting up the system of equations
    ! for the minimization of the Gibbs free energy.
    ! For two phase equilibria with hydrate and / or solid formation 5 cases with different sets of equations might occur
    !   1) number of components <= number of hydrate formers + 1 And ncomp /= 2
    !   2) number of components > number of hydrate formers + 1 And ncomp /= 2
    !   3) number of hydrate formers = 0  --> No hydrate formation but SV or SL equilibrium
    !   4) two components and water ice - hydrate equilibrium   (e.g. HIw)  !not for mixed hydrates
    !   5) two components and solid - hydrate equilibrium (e.g. HIc)        !not for mixed hydrates

    !Check if pure solid phases other than hydrate form, if not, case 1 or 2 applies
    if (gl%solidtype(1) == 0) then
        !Case 1)
        if (gl%ncomp <= gl%nrofhydrateformers+1) then
            GibbsEQN(1) = ChemPot_hyd - ChemPot_fluid(1)    !H2O always on first position on hydrate_list
            ! for the case of a p,T-flash another n-2 set of equations
            ! is needed (comes from the mass balance):
            if (iFlash == 3) then
                do j = 2, gl%ncomp - 1
                    GibbsEQN(j) = (z(j)-x_fluid(j))*(x_hyd(1)-x_fluid(1))- (z(1)-x_fluid(1))*(x_hyd(j)-x_fluid(j))
                end do
            End if

            !Case 2)
        elseif (gl%ncomp > gl%nrofhydrateformers+1) then
            GibbsEQN(1) = ChemPot_hyd - ChemPot_fluid(1)    !H2O always on first position on hydrate_list
            ! for the case of a p,T-flash another n-2 set of equations
            ! is needed (comes from the mass balance):
            if (iFlash == 3) then
                do j = 2, gl%nrofhydrateformers
                    GibbsEQN(j) = (z(j)-x_fluid(j))*(x_hyd(1)-x_fluid(1))- (z(1)-x_fluid(1))*(x_hyd(j)-x_fluid(j))
                end do
                do j = gl%nrofhydrateformers+1, gl%ncomp-1
                    GibbsEQN(j) = dlog(x_fluid(j)*z(gl%ncomp))-dlog(z(j)*x_fluid(gl%ncomp))
                end do
            End if
        End if

        !Solid phases other than hydrate form (case 3 to 5)
    else
        if ((gl%solidtype(1) > 0) .And. (gl%solidtype(2) == 0)) then !Solid phase (other than hydrate) in equilibrium with V or L phase
            GibbsEQN(1) = ChemPot_sol - ChemPot_fluid(gl%solid_pos)
            ! for the case of a p,T-flash another n-2 set of equations
            ! is needed (comes from the mass balance):
            if (iFlash == 3) then
                h = 1
                do j = 2, gl%ncomp - 1
                    if (h == gl%solid_pos) then
                        h = h + 1
                    End if
                    GibbsEQN(j) = dlog(x_fluid(h)*z(gl%ncomp))-dlog(z(h)*x_fluid(gl%ncomp))
                    h = h + 1
                end do
            End if
            !elseif ((ncomp == 2) .And. (solidtype(1) == 1) .AND. (solidtype(2) == 1)) then !Solid water and hydrate in Equilibrium (ONLY TWO COMPONENTS POSSIBLE). The fugacity of the guest is the unknown!!
            !   GibbsEQN(1) = ChemPot_hyd - ChemPot_sol
            !elseif ((ncomp == 2) .And. (solidtype(1) > 1) .AND. (solidtype(2) == 1)) then !Solid other than water and hydrate in Equilibrium (ONLY TWO COMPONENTS POSSIBLE)
            !   GibbsEQN(1) = 0.D0  !This has not to be solved, since T,p is known and the fugacity of the guest can be calculated from the solid equation
        End if

    End if

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    end subroutine
    !**************************************************************************

    !**************************************************************************
    module subroutine Jacobi_solid_NC_2P (gl,P, T, x_fluid, x_solid , x_hyd, beta_loc, fug_save, iFlash, JacMatrix, errval)
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
    !   P               - Pressure
    !   T               - Temperature
    !   solid           - indicates, which substance is solid (so far: 1 = CO2, 2 = water)
    !   solidnr         - gives the position of the fluid phase composition vector, where the substance stands, that becomes solid
    !   x_fluid         - Fluid phase composition
    !   x_solid         - Composition of the solid phase (0 except for the position solidnr)
    !   x_hyd           - Composition of the hydrate phase
    !   beta            - Molar "fluid" fraction (xi - xi_sol) / (xi_fluid - xi_sol)    0 --> solid phase only; 1 --> fluid phase only
    !   iFlash          - Flash mode
    ! OUTPUT:
    !   errval          - Error value
    !   JacMatrix       - 60 x 60 matrix containing the derivatives of all Gibbs-equations
    !                     F_i with respect to all independent variables X_i
    !--------------------------------------------------------------------------
    ! A. J輍er, March 2011







    implicit none

    type(type_gl) :: gl


    double precision:: T, p, beta_loc
    double precision, dimension(30):: x_fluid, x_solid, x_hyd
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: iFlash
    integer:: i, j, k, h
    double precision:: rho_fluid, dChempot_dT_sol, dChempot_dp_sol
    double precision, dimension(30, 30):: dChempoti_dxj
    double precision, dimension(30):: dChempoti_dT, d2nadnidT, dChempoti_dp, z, dnadni
    double precision:: rhoredmix_orig, tredmix_orig, Rmix, dummy

    !Variables needed for the hydrate model
    !double precision, dimension(30) :: lnf, dphiidT, dphiidp
    !double precision, dimension(30,30) :: dlnfidXj
    !double precision, dimension(30) :: dChempot_dxj_sol
    !double precision:: fug_CO2, DchpwDf, DchpwDT, DchpwDp, fug_save
    !double precision, dimension(30) :: doccupidfj, CiJ, fug_gas
    !double precision, dimension(30,30):: dxidfj, dxiHdxjFluid

    !Variables for the mixed hydrate model
    double precision, dimension(30):: lnf, dphiidT, dphiidp, dChempot_dxj_sol
    double precision, dimension(30,30):: dlnfidXj
    double precision, dimension(30):: fug_gas, DchpwDf, fug_save
    double precision:: DchpwDT, DchpwDp
    double precision, dimension(3,30):: CiJ
    double precision, dimension(3,30,30) :: doccupidfj
    double precision, dimension(30,30):: dxidfj, dxiHdxjFluid


    JacMatrix = 0.D0
    z = 0.D0
    dChempot_dT_sol = 0.D0
    dChempot_dp_sol = 0.D0
    dChempot_dxj_sol = 0.D0
    dChempoti_dxj = 0.D0
    dChempoti_dT = 0.D0
    dChempoti_dp = 0.D0
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    !if (T < ttp(solidnr)) then
    !if (p < ptp(solid_pos)) then
    rho_fluid = gl%rho_vap
    !else
    !    rho_fluid = rho_liq
    !end if

    z = gl%molfractions
    gl%molfractions = x_fluid
    call reduced_parameters_calc(gl,T)

    call R_mix_calc(gl,Rmix)

    select case (iFlash)
    case (1)
        !--------------------------------------------------------------------------
        ! Sublimation / Melting point calculation, p and x' or x" vector are given.
        ! T needs to be calculated
        !--------------------------------------------------------------------------
        !Calculate the derivative of the chemical potential with respect to T at constant
        !X and p for the fluid phase
        call d2na_dnidT_P (gl,T, rho_fluid, d2nadnidT, 0)
        call dna_dni(gl,T, rho_fluid, dnadni, 0)
        dChempoti_dT = d2nadnidT * Rmix * T + dnadni * Rmix
        ! get the solid phase properties
        if (gl%solidtype(1) > 0) then      !Pure solid substance
            select case (gl%solidtype(1))
            case(1)
                dChempot_dT_sol = dgdT_WaterIce(gl,T, p)  !Chemical Potential if Water Ice [J / mol K]
            case(2)
                dChempot_dT_sol = dgdT_DryIce(gl,T, p)    !Chemical Potential if Dry Ice [J / mol K]
                case default
                errval = -14444
                !write (*,*) 'error ', errval, ' -- Equation for solid does not exist'
                return
            end select
            ! Calculate the Jacobi-Matrix (which is a scalar in this case)
            JacMatrix(1, 1) = dChempot_dT_sol - dChempoti_dT(gl%solid_pos)
        else                            !Hydrate
            !The fugacity of guests in the vapor phase will be used as input for the chem. pot. of water in hydrate (mapping(2) = CO2, see Hydrate_list module variable)
            call lnf_mix(gl,T, rho_fluid, p, lnf)
            !fug_gas(1:nrofhydrateformers-1) = dexp(lnf(2:nrofhydrateformers))*1.d6     !Multi component hydrates
            fug_gas = fug_save

            !The derivative with respect to Temperature is calculated according to: dcp_dT|p,x_CO2_vap = dcp_dT|p,f_CO2 + dcp_df|T,p * dfi_dT|p, x_CO2
            call hdrt_DT_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDT)
            call hdrt_Df_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDf)
            call dlnphii_dT(gl,T, rho_fluid, dphiidT)
            ! pure hydrate: dChempot_dT_sol = DchpwDT + DchpwDf * dphiidT(2)*dexp(lnf(2))*1.d6
            ! mixed hydrate 2015:
            dChempot_dT_sol = DchpwDT
            Do i = 2, gl%nrofhydrateformers
                dChempot_dT_sol = dChempot_dT_sol + DchpwDf(i) * dphiidT(i)*dexp(lnf(i))*1.d6    ! Eq. 5.173 in Andreas' diss.
            End do

            ! Calculate the Jacobi-Matrix (which is a scalar in this case)
            JacMatrix(1, 1) = dChempot_dT_sol - dChempoti_dT(1) !Position of water, always 1
        end if

    case (2)
        !--------------------------------------------------------------------------
        ! Sublimation / Melting point calculation, T and x' or x" vector are given.
        ! p needs to be calculated
        !--------------------------------------------------------------------------
        !Calculate the derivative of the chemical potential with respect to p at constant
        !X and T for the fluid phase
        call d2na_dnidp_T (gl,T, rho_fluid, dChempoti_dp, 0)
        dChempoti_dp = dChempoti_dp * Rmix * T
        ! get the solid phase properties
        if (gl%solidtype(1) > 0) then      !Pure solid substance
            select case (gl%solidtype(1))
            case(1)
                dChempot_dp_sol = dgdp_WaterIce(gl,T, p)  !Chemical Potential of Water Ice [J / mol Pa]
            case(2)
                dChempot_dp_sol = dgdp_DryIce(gl,T, p)    !Chemical Potential of Dry Ice [J / mol Pa]
                case default
                errval = -14444
                !write (*,*) 'error ', errval, ' -- Equation for solid does not exist'
                return
            end select
            ! Calculate the Jacobi-Matrix (which is a scalar in this case)
            JacMatrix(1, 1) = dChempot_dp_sol - dChempoti_dp(gl%solid_pos)
        else
            !The fugacity of guests in the vapor phase will be used as input for the chem. pot. of water in hydrate (mapping(2) = CO2, see Hydrate_list module variable)
            call lnf_mix(gl,T, rho_fluid, p, lnf)
            !fug_gas(1:nrofhydrateformers-1) = dexp(lnf(2:nrofhydrateformers))*1.d6     !Multi component hydrates
            fug_gas = fug_save

            !The derivative with respect to pressure is calculated according to: dcp_dp|T,x_CO2_vap = dcp_dp|p,f_CO2 + dcp_df|T,p * dfi_dp|p, x_CO2
            call hdrt_Dp_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDp)
            call hdrt_Df_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDf)
            call dlnphii_dp(gl,T, rho_fluid, dphiidp)
            ! pure hydrate: dChempot_dp_sol = DchpwDp + DchpwDf * dexp(lnf(2))*1.d6 * (dphiidp(2) + 1.D0 / (p*1.D6))
            ! mixed hydrate 2015:
            dChempot_dp_sol = DchpwDp
            Do i = 2, gl%nrofhydrateformers
                dChempot_dp_sol = dChempot_dp_sol + DchpwDf(i) * dexp(lnf(i))*1.d6 * (dphiidp(i) + 1.D0 / (p*1.D6))    ! Eq. 5.175 in diss.
            End do

            ! Calculate the Jacobi-Matrix (which is a scalar in this case)
            JacMatrix(1, 1) = dChempot_dp_sol - dChempoti_dp(1) !Position of water, always 1
        end if


    case (3)
        !--------------------------------------------------------------------------
        ! p-T flash calculation, p and T and the overall x vector are given.
        ! x" or x' need to be calculated (because of the assumption, that the solid phase is pure,
        ! the composition of the solid phase does not need to be calculated)
        !--------------------------------------------------------------------------
        !Calculate the derivative of the chemical potential with respect to xj at constant
        !p and T for the fluid phase
        call d2na_dnidxj_PT (gl,T, rho_fluid, dChempoti_dxj, 0)
        dChempoti_dxj = dChempoti_dxj * Rmix * T
        ! get the solid phase properties
        if (gl%solidtype(1) > 0) then      !Pure solid substance
            dChempot_dxj_sol = 0.D0     !pure solid phase has a fixed composition
        else                            !Hydrate
            !Get the fugacities of the fluid phase (THIS IS NEEDED FOR THE HYDRATE MODEL)
            call lnf_mix(gl,T, rho_fluid, p, lnf)
            call dlnfi_dxj_TP (gl,T, rho_fluid, dlnfidXj)
            !The fugacity of guests in the vapor phase will be used as input for the chem. pot. of water in hydrate (mapping(2) = CO2, see Hydrate_list module variable)
            !fug_gas(1:nrofhydrateformers-1) = dexp(lnf(2:nrofhydrateformers))*1.d6     !Multi component hydrates
            fug_gas = fug_save
            do i = 1,gl%nrofhydrateformers-1
                if (is_infinity(fug_gas(i))) then
                    errval = -15566
                    return
                endif
            enddo
            call hdrt_Df_chem_potent_w(gl,T,p*1.D6,fug_gas,DchpwDf)
            call hdrt_dxi_dfj(gl,T,p*1.D6,fug_gas,doccupidfj,CiJ,dxidfj)

            !The derivative with respect to x_CO2 is calculated according to: dcp_dx_CO2|p,T = dcp_df|T,p * dfi_dx_CO2|T,p
            dChempot_dxj_sol = 0.d0
            Do j = 1, gl%ncomp - 1
                Do i = 2, gl%nrofhydrateformers
                    !dChempot_dxj_sol(j) = DchpwDf * dlnfidxj(j, i)*dexp(lnf(i))*1.d6       !Only one guest
                    !dChempot_dxj_sol(j) = DchpwDf(i) * dlnfidxj(2, 2)*dexp(lnf(2))*1.d6   !Multi component hydrates

                    ! mixed hydrate 2015:
                    dChempot_dxj_sol(j) = dChempot_dxj_sol(j) + DchpwDf(i) * dlnfidxj(j, i)*dexp(lnf(i))*1.d6
                End do
            End Do

            Do j = 1, gl%ncomp - 1 ! x(j) in the fluid
                Do i = 1, gl%nrofhydrateformers
                    !Do k = 2, nrofhydrateformers
                    !    dxiHdxjFluid(j,i) = dxidfj(k,i)*dlnfidxj(j, k)*dexp(lnf(i))*1.d6   !Derivative of the molfraction of component i in hydrate wrt. the molfraction j in the fluid
                    !End do

                    ! mixed hydrate 2015:
                    dxiHdxjFluid(j,i) = 0.d0
                    Do k = 2, gl%nrofhydrateformers    ! k-guest fugacity
                        dxiHdxjFluid(j,i) = dxiHdxjFluid(j,i) + dxidfj(k,i) * dlnfidxj(j, k)*dexp(lnf(k))*1.d6 !+  oder k in lnf??    ! Eq. 5.179 in Andreas' diss.
                    End do
                End Do
            End Do

        End if
        ! Calculate the Jacobi-Matrix
        ! For two phase equilibria with hydrate and / or solid formation 5 cases with different sets of Jacobi matices might occur
        !   1) number of components <= number of hydrate formers + 1 And ncomp /= 2
        !   2) number of components > number of hydrate formers + 1 And ncomp /= 2
        !   3) number of hydrate formers = 0  --> No hydrate formation but SV or SL equilibrium
        !   4) two components and water ice - hydrate equilibrium   (e.g. HIw)
        !   5) two components and solid - hydrate equilibrium (e.g. HIc)
        !open(unit=59,file="Jacobi.txt",status="replace",action="write")
        !4059 format (f16.6)
        !Check if pure solid phases other than hydrate form, if not, case 1 or 2 applies
        if (gl%solidtype(1) == 0) then
            !Case 1)
            if (gl%ncomp <= gl%nrofhydrateformers+1) then
                Do j = 1, gl%ncomp - 1
                    JacMatrix(j,1) = dChempot_dxj_sol(j) - dChempoti_dxj(j,1)
                    Do i = 2, gl%ncomp - 1
                        if (j == 1) then    !Derivative wrt. molfraction of water in fluid phase
                            !Old probably wrong
                            !JacMatrix(j,i) = (z(i)-x_fluid(i))*(dxiHdxjFluid(j,1)-1)+ &
                            !           & (x_hyd(i)-x_fluid(i))+(z(i)-x_fluid(i))*dxiHdxjFluid(1,i)
                            JacMatrix(j,i) = (z(i)-x_fluid(i))*(dxiHdxjFluid(j,1)-1)+ &
                                & (x_hyd(i)-x_fluid(i))-(z(i)-x_fluid(i))*dxiHdxjFluid(1,i)
                        else
                            if (j == i) then
                                JacMatrix(j,i) = -(x_hyd(1)-x_fluid(1))+(z(i)-x_fluid(i))*dxiHdxjFluid(j,1) &
                                    & -(z(1)-x_fluid(1))*(dxiHdxjFluid(j,i)-1)
                                !JacMatrix(j,i) = -(x_hyd(1)-x_fluid(1))+(z(i)-x_fluid(i))*dxiHdxjFluid(j,1) &
                                !           & -(z(1)-x_fluid(1))*dxiHdxjFluid(j,i)
                            End if
                            if (j /= i) then
                                JacMatrix(j,i) = (z(i)-x_fluid(i))*dxiHdxjFluid(j,1)-(z(1)-x_fluid(1))*dxiHdxjFluid(j,i)
                                !JacMatrix(j,i) = (x_hyd(i)-x_fluid(i))+(z(i)-x_fluid(i))*dxiHdxjFluid(j,1)-(z(1)-x_fluid(1))*dxiHdxjFluid(j,i)
                            End if
                        End if
                    End Do
                End Do


                !Case 2)
            elseif (gl%ncomp > gl%nrofhydrateformers+1) then
                Do j = 1, gl%ncomp - 1
                    JacMatrix(j,1) = dChempot_dxj_sol(j) - dChempoti_dxj(j,1)
                    Do i = 2, gl%nrofhydrateformers
                        if (j == 1) then    !Derivative wrt. molfraction of water in fluid phase
                            JacMatrix(j,i) = (z(i)-x_fluid(i))*(dxiHdxjFluid(j,1)-1)+ &
                                & (x_hyd(i)-x_fluid(i))+(z(i)-x_fluid(i))*dxiHdxjFluid(1,i)
                        else
                            if (j == i) then
                                JacMatrix(j,i) = -(x_hyd(1)-x_fluid(1))+(z(i)-x_fluid(i))*dxiHdxjFluid(j,1) &
                                    & -(z(1)-x_fluid(1))*(dxiHdxjFluid(j,i)-1)
                            End if
                            if ((j /= i) .or. (j > gl%nrofhydrateformers)) then
                                JacMatrix(j,i) = (z(i)-x_fluid(i))*dxiHdxjFluid(j,1)-(z(1)-x_fluid(1))*dxiHdxjFluid(j,i)
                            End if
                        End if
                    End Do
                    Do i = gl%nrofhydrateformers+1, gl%ncomp - 1
                        if (j == i) then
                            JacMatrix(j,i) = 1.d0 / x_fluid(i) + 1.D0 / x_fluid(gl%ncomp)
                        End if
                        if ((j /= i) .or. (j <= gl%nrofhydrateformers)) then
                            JacMatrix(j,i) = 1.D0 / x_fluid(gl%ncomp)
                        End if
                    End Do
                End Do
            End if

            !Solid phases other than hydrate form (case 3 to 5)
        else
            if ((gl%solidtype(1) > 0) .And. (gl%solidtype(2) == 0)) then   !Solid phase (other than hydrate) in equilibrium with V or L phase
                Do j = 1, gl%ncomp - 1
                    JacMatrix(j,1) = dChempot_dxj_sol(j) - dChempoti_dxj(j,gl%solid_pos)
                    h = 2
                    Do i = 2, gl%ncomp -1
                        if (h == gl%solid_pos) then
                            h = h + 1
                        End if
                        if (j == i) then
                            JacMatrix(j,i) = 1.d0 / x_fluid(h) + 1.D0 / x_fluid(gl%ncomp)
                        End if
                        if (j /= i) then
                            JacMatrix(j,i) = 1.D0 / x_fluid(gl%ncomp)
                        End if
                        h = h + 1
                    End do
                end do




                !elseif ((ncomp == 2) .And. (solidtype(1) == 1) .AND. (solidtype(2) == 1)) !Solid water and hydrate in Equilibrium (ONLY TWO COMPONENTS POSSIBLE). The fugacity of the guest is the unknown!!
                !    call hdrt_Df_chem_potent_w(T,p*1.D6,fug_CO2,DchpwDf)
                !    JacMatrix(1,1) = DchpwDf
                !             elseif ((ncomp == 2) .And. (solidtype(1) > 1) .AND. (solidtype(2) == 1)) !Solid other than water and hydrate in Equilibrium (ONLY TWO COMPONENTS POSSIBLE)
                !                GibbsEQN(1) = 0.D0  !This has not to be solved, since T,p is known and the fugacity of the guest can be calculated from the solid equation

            End if

        End if

        !--------------------------------------------------------------------------
        case default
        errval = -1111
        !--------------------------------------------------------------------------
    end select

    ! set the module variables back to original values
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = z

    !write(59,4059) JacMatrix
    end subroutine
    !**************************************************************************


    end submodule impl
