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


    module phase_properties
    type type_phase_properties
        double precision, dimension(30, 3) :: x_phase
        double precision, dimension(3) :: gibbs, phasefrac, rho
        double precision :: gibbs_tot
        character(5), dimension(3) :: text_phase
        integer, dimension(5) :: phasetype
        integer :: errorflag, nrofphases, hdrt_structure_flag
        integer, dimension(2) :: solidtype
    end type
    end module phase_properties

    ! module for file phasedet_sol.f90
    module phasedet_sol_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use rhomix_pt_module
    use flash_module
    use phasedet_mix_module
    use waterice_module
    use dryice_module
    use ptflash_2c_3p_module
    use hdrt_chem_pot_module
    use phasedet_mix_module
    use ancillary_equations_mix_module
    use ptflash_solids_mix_module
    use hdrt_properties_modules_module
    use reduced_parameters_calc_module


    contains

    !************************************************************************************
    subroutine PhaseDet_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified mixture
    ! This algorithms does also check if any solid phases form
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !
    ! OUTPUT:
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   x_phase(5,30)   - matrix that contains all phase compositions
    !       1   x_vap       - Vapor phase composition
    !       2   x_liq1      - (lighter) liquid phase composition
    !       3   x_liq2      - (heavier) liquid phase composition
    !       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !       5   x_hyd       - hydrate phase compsition (so far only CO2 hydrates are possible)
    !   phasefrac           - Phasefraction of the phases
    !   nrofphases  - Number of phases present
    !   errval      - Error value
    !************************************************************************************
    !Andreas Jan. 2012
    !Sebastian New Algorithm Nov. 2018

    use, intrinsic :: iso_c_binding


    use phase_properties
    use help_funcs

    implicit none

    type(type_gl) :: gl
    type(type_gl), pointer :: gl_check
    type(type_phase_properties), dimension(:), pointer :: ph_prop
    logical :: use_old_stability !temporary switch
    double precision :: press, temp
    double precision, dimension(30) :: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, dimension (5) :: phasefrac
    integer :: errval, nrofphases

    double precision :: psub, pmelt, peq
    double precision, dimension(30) :: Chempot, lnf, x_solid, x_hyd, x_fluid1, x_fluid2, x_hyd1, x_hyd2
    double precision, dimension(30) :: x_phase1, x_phase2, x_fluid
    double precision :: ChemPot_hyd, rhofluid_est, tpd_sol, rho1, rho2, vapfrac1, vapfrac2, rho1_H, rho2_H      ! , fug_CO2
    double precision :: tpd_sol_phase1, tpd_sol_phase2, tpd_phase1, tpd_phase2, rho_phase1, rho_phase2
    double precision :: tpd_sol_hyd(2) !2 dimensional for loop over s1 and s2 hydrate structure
    integer :: errval1, errval2, iFlash, i, iter, errval_fluid, nrfluidphases, k
    logical :: check_Hydsolwater

    !Indicate which phases are present
    integer, dimension(5) :: phasetype         !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !--> nrofphases = 2
    !--> phasetype(1) = 2 (light liquid)
    !--> phasetype(2) = 3 (heavy liquid)
    !E.g.: 1 phase present: liquid
    !--> nrofphases = 1
    !--> phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!

    integer :: n, j, iphase, iPhase_try, iPhase_sol
    logical :: stable, same_comp, is_VL1
    double precision, dimension(30) :: x_vap, x_liq1, x_liq2, x_VL1_V, x_VL1_L1, x_VL2_V, x_VL2_L2, x_L1L2_L1, x_L1L2_L2, x_help!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision :: rho_orig, tpd_liq1, tpd_vap, rhoVL1_V, rhoVL1_L1, rhoVL2_V, rhoVL2_L2, rhoL1L2_L1, rhoL1L2_L2
    double precision :: rho_VS_V, rho_VS_S, rho_LS_L, rho_LS_S
    double precision :: vapfrac_L1L2, vapfrac_VL1, vapfrac_VL2
    double precision:: vapfrac, rhoredmix_orig, tredmix_orig, rhovap_est, rholiq1_est, rholiq2_est
    double precision::  gVL1, gVL2, gL1L2, g_VL1L2, g_f1f2, gf1, gf2, gf3, df1, df2, df3, press_3phase, gVIw, gLIw ! G_CALC
    double precision, dimension(3):: phasefrac_VL1L2, phasefrac_3Phase_solid, rho_sol
    double precision :: ps_compi, ps_compj

    !double precision, dimension(2) :: occup, CiJ, xH
    double precision :: rho_H

    double precision, dimension (60,60) :: Mat_phasefrac_3Comp
    double precision, dimension (60) :: Vec_phasefrac_3Comp, Vec_rightside_3Comp
    integer:: eqn
    character(255):: herr

    integer :: pos_dryIce

    logical :: phys_reas1, phys_reas2   !Save the information if a test equilibrium found is physically reasonable. If yes, use this equilibrium as backup solution if no certain correct equilibrium can be found
    !This became necessary because for the mixture fluidl = 'water;co2;nitrogen' molesl = '0.5;0.49;0.01' calculated with EOS_CG for press = 5.6D0 and Temp = 268.018310851428D0
    !it is correctly predicted that solid H2O will form, but no two phase equilbrium can be found which does not lead to a negative tpd for the other fluid phase --> strange solution!
    !The different gas constants for the fluid are suspected to be the source of this error

    logical :: error_then_ice_found     !Variable to store information if ice was already checked when an error occured before
    double precision :: curvature

    integer, dimension(1) :: pos_L1, pos_L2

    !Mixed hydrates 2015
    double precision, dimension(3,30) :: occup, CiJ, occup_single, occup_double
    double precision, dimension(30) :: xH, fug_gas, fug_CO2
    integer :: hdrt_structure

    logical :: h2o_present, co2_present
    integer :: nr_checks

    errval_fluid = 0
    errval= 0

    phasefrac_VL1L2 = 0.D0
    phasefrac_3Phase_solid = 0.D0
    phasefrac = 0.D0
    tpd_phase1 = 1.D6
    tpd_phase2 = 1.D6
    phys_reas1 = .false.
    phys_reas2 = .false.

    error_then_ice_found = .false.

    call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)

    !------
    !If only one stable phase was found and no error occured, check whether solid water forms (if water is in the mixture)
    !For some mixtures with water it might happen that a phase with relatively high water amount is mispredicted to be stable
    !In these cases the following checks won't work, so this has to be checked in the beginning
    !This test catches errors of PhaseDet, however, it can also lead to problems finding phase equilibria in very water rich mixtures. Therefore, this has to be treated carefully! Andreas Feb. 2015
    !Pure Hydrates
    !If ((nrofhydrateformers > 0) .and. (nrofhydrateformers < 3) .and. (errval == 0)) then   !Water is in the mixture (on position 1)
    !Mixed Hydrates
    If ((gl%nrofhydrateformers > 0) .and. (errval == 0)) then   !Water is in the mixture (on position 1)
        if (nrofphases == 1) then
            !Check if solid water forms.
            gl%molfractions = x_phase(:,phasetype(1))
            call reduced_parameters_calc(gl,Temp)
            call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
            gl%molfractions = x
            call reduced_parameters_calc(gl,Temp)
            tpd_sol = g_WaterIce(gl,Temp,press) - Chempot(1) !consider seawater here, due to shift of chempot w
            !If solid water forms, it would still be possible that all other components are dissolved in liquid water which is in equilibrium with solid water
            !Hence, do not produce an error if the mixture consists of almost only water. The minumum water content was chosen arbitrarily!! Andreas Feb. 2015.
            !SH: Added check if phase in equlibirum with water is really liquid
            if ((tpd_sol < -1.D-14) .and. (x_phase(1,phasetype(1)) < 0.984)) then!.and. (any(phasetype(:) /= 1))  then
                errval = -101010
                !Internal error value
            end if
        end if
    end if
    !------

    if (errval == 0) then
        nrfluidphases = nrofphases
        if (nrofphases > 1) then
            Phasefrac(phasetype(1)) = vapfrac
            Phasefrac(phasetype(2)) = 1-vapfrac
        else
            Phasefrac(phasetype(1)) = 1.D0
        end if

        rhoredmix_orig = gl%rhoredmix
        tredmix_orig = gl%tredmix
    end if

    !Andreas February 2014
    !Check of the correct two phase equilibrium was found and check whether three phases coexist at equilibrium
    If ((errval == 0) .and. (nrofphases == 2)) then

        !Get overall Gibbs energy of phase equilibrium found
        !calculate light phase
        gl%molfractions = x_phase(:,phasetype(1))
        call reduced_parameters_calc(gl,Temp)
        df1 = rho(phasetype(1))
        if((gl%seawater) .and. ((phasetype(1) .eq. 2) .or. (phasetype(1) .eq. 3))) gl%seacalc = .true.
        gf1 = G_CALC(gl,Temp, df1, 0)
        gl%seacalc = .false.
        !calculate heavy phase
        gl%molfractions = x_phase(:,phasetype(2))
        call reduced_parameters_calc(gl,Temp)
        df2 = rho(phasetype(2))
        if((gl%seawater) .and. ((phasetype(2) .eq. 2) .or. (phasetype(2) .eq. 3))) gl%seacalc = .true.
        gf2 = G_CALC(gl,Temp, df2, 0)
        gl%seacalc = .false.
        g_f1f2 = gf2 + vapfrac * (gf1 - gf2)

        !Set back reducing parameters
        gl%molfractions = x
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig

        !logical to check if the tpd minimization found a new phase or if it converged to one of the already existing phases
        same_comp = .true.
        !VL1 equlibrium: is_VL1 = true, VL2 equilibrium: is_VL1 = false
        is_VL1 = .false.

        !Case: VL2 Vapor-liquid equilibrium or VL1 equilibrium
        If (phasetype(1) == 1) then
            x_vap = x_phase(:,phasetype(1))
            x_liq2 = x_phase(:,phasetype(2))
            !Check whether a second liquid phase exists
            !First identify comonents in the mixture that might be liquid (ps < p)
            iphase = 1 !Density solver prefers liquid solution
            do i = 1, gl%ncomp
                if ((temp < gl%tc(i)) .and. (press > (0.9D0 * vp_eq(gl,Temp, i))) ) then    !This component potentially forms another liquid phase
                    rho_orig = rho(phasetype(1))
                    n = 20  !Take a maximum of 20 steps of successice substitution
                    !Start the minimization of the tangent plane distant close to xi = 1 with i as the potential liquid former
                    x_liq1(i) = 0.999D0
                    do j = 1, gl%ncomp
                        if (i /= j) then
                            x_liq1(j) = (1.D0 - x_liq1(i))/(gl%ncomp-1)
                        end if
                    end do
                    !Minimize the tpd and check whether the original phases are stable
                    call Succ_Sub_tpd(gl,press, temp, x_vap, rho_orig, x_liq1, tpd_liq1, iphase, Stable, n, errval)
                    if (errval == 0) then
                        !Check whether a different composition than x_vap and x_liq1 was found
                        do j = 1, gl%ncomp
                            !if one of the mole fractions is significantly different, another composition was found
                            if ((dabs(x_vap(j) - x_liq1(j)) > 1.D-4) .and. (dabs(x_liq2(j) - x_liq1(j)) > 1.D-4)) then
                                same_comp = .false.
                                exit
                            end if
                        end do
                        !If another liquid phase was found, check whether the tpd at this composition is negative (other two phase equilibrium) or close to 0
                        if (.not. same_comp) then
                            if (tpd_liq1 < -1.D-12) then                                    !Different two phase equilibrium or three phase equilibrium

                                !Check if the original equilibrium is VL1 or VL2
                                !This is important!
                                !   In case of VL2 the other equilibrium to check for is L1L2 !!!
                                !   In case of VL1 the other equilibrium to check for is VL2 !!!
                                !Find out which component is the main component in the liquid incipient phase
                                Do j = 1, gl%ncomp
                                    if (x_liq2(j) > 0.5D0) exit
                                end do
                                !Calculate the vapor pressure for both liquid formers
                                ps_compi = vp_eq(gl,Temp, i)   !new liquid phase
                                ps_compj = vp_eq(gl,Temp, j)   !originally found liquid in equilibrium
                                !Compare the vapor pressures. The more volatile liquid is liq1, the less volatile liq2
                                if (ps_compi < ps_compj) then   !new phase is less volatile --> originally VL1 found
                                    x_help = x_liq1
                                    x_liq1 = x_liq2
                                    x_liq2 = x_help
                                    is_VL1 = .true.
                                else                            !new phase is more volatile --> originally VL2 found
                                    is_VL1 = .false.
                                end if

                                !If the tangent plane distance is "close" to 0, try to calculate a three phase equilibrium
                                !---
                                if (dabs(tpd_liq1) < 10.D0) then
                                    if (gl%ncomp > 2) then
                                        rhovap_est = 0.D0
                                        rholiq1_est = 0.D0
                                        rholiq2_est = 0.D0
                                        iFlash = 7 !pT given
                                        !Startvalue for phasefractions
                                        !If 3 components are in a three phase equilibrium, the startvalues for the phasefractions can be solved analytically
                                        if (gl%ncomp == 3) then
                                            Mat_phasefrac_3Comp(1,1) = x_vap(1)
                                            Mat_phasefrac_3Comp(2,1) = x_liq1(1)
                                            Mat_phasefrac_3Comp(3,1) = x_liq2(1)

                                            Mat_phasefrac_3Comp(1,2) = x_vap(2)
                                            Mat_phasefrac_3Comp(2,2) = x_liq1(2)
                                            Mat_phasefrac_3Comp(3,2) = x_liq2(2)

                                            Mat_phasefrac_3Comp(1,3) = 1.D0
                                            Mat_phasefrac_3Comp(2,3) = 1.D0
                                            Mat_phasefrac_3Comp(3,3) = 1.D0

                                            Vec_rightside_3Comp(1) = x(1)
                                            Vec_rightside_3Comp(2) = x(2)
                                            Vec_rightside_3Comp(3) = 1.D0

                                            Vec_phasefrac_3Comp = 0.D0

                                            !Solve the system of equations
                                            eqn = 3
                                            call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                                            if (errval /= 0) then
                                                phasefrac_VL1L2 = 1.D0 / 3.D0
                                            else
                                                phasefrac_VL1L2 = Vec_rightside_3Comp(1:3)
                                            end if
                                        else
                                            !Arbitrary start values
                                            phasefrac_VL1L2 = 1.D0 / 3.D0
                                        end if
                                        !save start values if next steps leads to same compositions of two phases
                                        x_VL1_V = x_vap
                                        x_VL1_L1 = x_liq1
                                        x_VL2_V = x_vap
                                        x_VL2_L2 = x_liq2
                                        call ptflash_NC_3P(gl,press, Temp, x, rho, x_vap, x_liq1, x_liq2, rhovap_est, &
                                            & rholiq1_est, rholiq2_est, phasefrac_VL1L2, iFlash, iter, errval)
                                        if (minval(phasefrac_VL1L2) < -1.D-16) then
                                            errval = -3333
                                        end if
                                        if (errval == 0) then
                                            phasetype(1) = 1
                                            phasetype(2) = 2
                                            phasetype(3) = 3
                                            x_Phase(:,phasetype(1)) = x_vap
                                            x_Phase(:,phasetype(2)) = x_liq1
                                            x_Phase(:,phasetype(3)) = x_liq2
                                            Phasefrac(phasetype(1)) = phasefrac_VL1L2(1)
                                            Phasefrac(phasetype(2)) = phasefrac_VL1L2(2)
                                            Phasefrac(phasetype(3)) = phasefrac_VL1L2(3)
                                            nrofphases = 3
                                            nrfluidphases = 3

                                            !Save information of the equilibrium found
                                            gl%solidtype_akt_phase = 0             !No pure solids
                                            gl%solidpos_akt_phase = 0              !No pure solids

                                            !!Get overall Gibbs energy of phase equilibrium found
                                            !!vapor phase
                                            !molfractions = x_vap
                                            !call reduced_parameters_calc(gl,Temp)
                                            !df1 = rho(phasetype(1))
                                            !gf1 = G_CALC(gl,Temp, df1, 0)
                                            !!liquid phase 1
                                            !molfractions = x_liq1
                                            !call reduced_parameters_calc(gl,Temp)
                                            !df2 = rho(phasetype(2))
                                            !gf2 = G_CALC(gl,Temp, df2, 0)
                                            !!liquid phase 2
                                            !molfractions = x_liq2
                                            !call reduced_parameters_calc(gl,Temp)
                                            !df3 = rho(phasetype(3))
                                            !gf3 = G_CALC(gl,Temp, df3, 0)
                                            !
                                            !g_VL1L2 = Phasefrac(phasetype(1))* gf1 + Phasefrac(phasetype(2)) * gf2 + Phasefrac(phasetype(3))* gf3
                                            !!Set back reducing parameters
                                            !molfractions = x
                                            !rhoredmix = rhoredmix_orig
                                            !tredmix = tredmix_orig

                                            !exit
                                            !return

                                        end if
                                        errval = 0
                                    else
                                        !rhovap_est = 0.D0
                                        !rholiq1_est = 0.D0
                                        !rholiq2_est = 0.D0
                                        !iFlash = 5 !T given, Bubble point
                                        !call ptflash_NC_3P(press_3phase, Temp, x, rho, x_vap, x_liq1, x_liq2, rhovap_est, &
                                        !   & rholiq1_est, rholiq2_est, phasefrac_VL1L2, iFlash, iter, errval)
                                    end if
                                end if
                                !---

                                if (is_VL1) then
                                    x_vap(1:gl%ncomp) = x_VL2_V(1:gl%ncomp)
                                    x_liq2(1:gl%ncomp) = x_VL2_L2(1:gl%ncomp)
                                    !Calculate VL2
                                    !Use the same variables as for VL1 (save some lines of code)
                                    iPhase_try = 5  !VLE expected
                                    rhovap_est = 0.D0
                                    rholiq1_est = 0.D0
                                    call Flash_pT_calc(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq1_est, vapfrac_VL1, iPhase_try, errval, iter)
                                    if (errval == 0) then
                                        x_VL1_V = x_vap
                                        x_VL1_L1 = x_liq2
                                        rhoVL1_V = gl%rho_vap
                                        rhoVL1_L1 = gl%rho_liq
                                        !Get overall Gibbs energy of phase equilibrium found
                                        !calculate light phase
                                        gl%molfractions = x_vap
                                        call reduced_parameters_calc(gl,Temp)
                                        df1 = rhoVL1_V
                                        gf1 = G_CALC(gl,Temp, df1, 0)
                                        !calculate heavy phase
                                        gl%molfractions = x_liq1
                                        call reduced_parameters_calc(gl,Temp)
                                        df2 = rhoVL1_L1
                                        gf2 = G_CALC(gl,Temp, df2, 0)
                                        gVL2 = gf2 + vapfrac_VL1 * (gf1 - gf2)

                                        !Set back reducing parameters
                                        gl%molfractions = x
                                        gl%rhoredmix = rhoredmix_orig
                                        gl%tredmix = tredmix_orig
                                    else
                                        !Dummy value
                                        gVL1 = 1.D10
                                        errval = 0
                                    end if

                                else
                                    x_vap(1:gl%ncomp) = x_VL1_V(1:gl%ncomp)
                                    x_liq1(1:gl%ncomp) = x_VL1_L1(1:gl%ncomp)
                                    !Calculate VL1
                                    iPhase_try = 5  !VLE expected
                                    rhovap_est = 0.D0
                                    rholiq1_est = 0.D0
                                    call Flash_pT_calc(gl,press, Temp, x, x_vap, x_liq1, rhovap_est, rholiq1_est, vapfrac_VL1, iPhase_try, errval, iter)
                                    if (errval == 0) then
                                        x_VL1_V = x_vap
                                        x_VL1_L1 = x_liq1
                                        rhoVL1_V = gl%rho_vap
                                        rhoVL1_L1 = gl%rho_liq
                                        !Get overall Gibbs energy of phase equilibrium found
                                        !calculate light phase
                                        gl%molfractions = x_vap
                                        call reduced_parameters_calc(gl,Temp)
                                        df1 = rhoVL1_V
                                        gf1 = G_CALC(gl,Temp, df1, 0)
                                        !calculate heavy phase
                                        gl%molfractions = x_liq1
                                        call reduced_parameters_calc(gl,Temp)
                                        df2 = rhoVL1_L1
                                        gf2 = G_CALC(gl,Temp, df2, 0)
                                        gVL2 = gf2 + vapfrac_VL1 * (gf1 - gf2)

                                        !Set back reducing parameters
                                        gl%molfractions = x
                                        gl%rhoredmix = rhoredmix_orig
                                        gl%tredmix = tredmix_orig
                                    else
                                        !Dummy value
                                        gVL1 = 1.D10
                                        errval = 0
                                    end if

                                end if

                                !Calculate L1L2
                                iPhase_try = 4  !LLE expected
                                rholiq1_est = 0.D0
                                rholiq2_est = 0.D0
                                call Flash_pT_calc(gl,press, Temp, x, x_liq1, x_liq2, rholiq1_est, rholiq2_est, vapfrac_L1L2, iPhase_try, errval, iter)
                                if (errval == 0) then
                                    x_L1L2_L1 = x_liq1
                                    x_L1L2_L2 = x_liq2
                                    rhoL1L2_L1 = gl%rho_vap
                                    rhoL1L2_L2 = gl%rho_liq
                                    !Get overall Gibbs energy of phase equilibrium found
                                    !calculate light phase
                                    gl%molfractions = x_liq1
                                    call reduced_parameters_calc(gl,Temp)
                                    df1 = rhoL1L2_L1
                                    gf1 = G_CALC(gl,Temp, df1, 0)
                                    !calculate heavy phase
                                    gl%molfractions = x_liq2
                                    call reduced_parameters_calc(gl,Temp)
                                    df2 = rhoL1L2_L2
                                    gf2 = G_CALC(gl,Temp, df2, 0)
                                    gL1L2 = gf2 + vapfrac_L1L2 * (gf1 - gf2)

                                    !Set back reducing parameters
                                    gl%molfractions = x
                                    gl%rhoredmix = rhoredmix_orig
                                    gl%tredmix = tredmix_orig
                                else
                                    !Dummy value
                                    gL1L2 = 1.D10
                                    errval = 0
                                end if

                                !Check which equilibrium yields the lowest gibbs energy
                                gVL2 = g_f1f2
                                if ((gL1L2 < gVL1) .and. (gL1L2 < gVL2)) then
                                    phasetype(1) = 2
                                    phasetype(2) = 3
                                    x_Phase(:,phasetype(1)) = x_liq1
                                    x_Phase(:,phasetype(2)) = x_liq2
                                    Phasefrac(phasetype(1)) = vapfrac_L1L2
                                    Phasefrac(phasetype(2)) = 1-vapfrac_L1L2
                                    rho(phasetype(1)) = rhoL1L2_L1
                                    rho(phasetype(2)) = rhoL1L2_L2
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 0             !No pure solids
                                    gl%solidpos_akt_phase = 0              !No pure solids
                                    exit
                                elseif((gVL1 < gVL2) .and. (gVL1 < gL1L2)) then
                                    phasetype(1) = 1
                                    phasetype(2) = 3
                                    x_Phase(:,phasetype(1)) = x_vap
                                    x_Phase(:,phasetype(2)) = x_liq1
                                    Phasefrac(phasetype(1)) = vapfrac_VL1
                                    Phasefrac(phasetype(2)) = 1-vapfrac_VL1
                                    rho(phasetype(1)) = rhoVL1_V
                                    rho(phasetype(2)) = rhoVL1_L1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 0             !No pure solids
                                    gl%solidpos_akt_phase = 0              !No pure solids
                                    exit
                                else
                                    !Backup if strange things happen (theoretically it is not possible that VL1 is stable, since it was found by the tpd check that the second liquid is stable)
                                    !errval = -324424
                                    !If this happens, continue the calculations with the original phases coming from PhaseDet instead of quitting with error
                                    exit
                                end if
                            else                                                            !No other stable phase found
                                !do nothing
                            end if
                        end if
                    end if
                end if

            end do


            !Liquid / liquid equilibrium found.
            !Andreas, October 2014: Check for VLL three phase equilibrium, because it can happen at elevated pressures that LL equilibria are found instead of VL equilibria, with VLLE being the correct equilibrium
            !Needs only to be checked for mixtures with more than 2 components
        else

            if (gl%ncomp > 2) then
                !Test if the "lighter" liquid could be a gas phase
                x_liq1 = x_phase(:,phasetype(1))
                x_liq2 = x_phase(:,phasetype(2))
                !Check whether a second liquid phase exists
                iphase = 1 !Density solver prefers vapor solution
                !At a first step, identify the two liquid forming components
                pos_L1 = maxloc(x_liq1)
                pos_L2 = maxloc(x_liq2)
                i = pos_L1(1)
                j = pos_L2(1)
                !Set the mole fractions of the liquid formers in the initial composition for the tpd minimization to a low value
                x_vap(pos_L1) = 1.D-5
                x_vap(pos_L2) = 1.D-5
                do k = 1, gl%ncomp
                    if ((k /= i) .and. (k /= j)) then
                        x_vap(k) = (1.D0 - x_vap(i) - x_vap(j))/(gl%ncomp-2)
                    end if
                end do
                !Minimize the tpd and check whether the original phases are stable
                rho_orig = rho(phasetype(1))    !Density of lighter liquid phase
                n = 20  !Take a maximum of 20 steps of successice substitution
                call Succ_Sub_tpd(gl,press, temp, x_liq1, rho_orig, x_vap, tpd_vap, iphase, Stable, n, errval)
                if (errval == 0) then
                    !Check whether a different composition than x_liq1 and x_liq2 was found
                    do j = 1, gl%ncomp
                        !if one of the mole fractions is significantly different, another composition was found
                        if ((dabs(x_vap(j) - x_liq1(j)) > 1.D-4) .and. (dabs(x_liq2(j) - x_vap(j)) > 1.D-4)) then
                            same_comp = .false.
                            exit
                        end if
                    end do
                    !If a vapor phase was found, check whether the tpd at this composition is close to 0, if yes, try to calculate a three phase equilibrium
                    if (.not. same_comp) then
                        if (dabs(tpd_vap) < 10.D0) then                                    !Try three phase equilibrium if tpd is small

                            rhovap_est = 0.D0
                            rholiq1_est = 0.D0
                            rholiq2_est = 0.D0
                            iFlash = 7 !pT given
                            !Startvalue for phasefractions
                            !If 3 components are in a three phase equilibrium, the startvalues for the phasefractions can be solved analytically
                            if (gl%ncomp == 3) then
                                Mat_phasefrac_3Comp(1,1) = x_vap(1)
                                Mat_phasefrac_3Comp(2,1) = x_liq1(1)
                                Mat_phasefrac_3Comp(3,1) = x_liq2(1)

                                Mat_phasefrac_3Comp(1,2) = x_vap(2)
                                Mat_phasefrac_3Comp(2,2) = x_liq1(2)
                                Mat_phasefrac_3Comp(3,2) = x_liq2(2)

                                Mat_phasefrac_3Comp(1,3) = 1.D0
                                Mat_phasefrac_3Comp(2,3) = 1.D0
                                Mat_phasefrac_3Comp(3,3) = 1.D0

                                Vec_rightside_3Comp(1) = x(1)
                                Vec_rightside_3Comp(2) = x(2)
                                Vec_rightside_3Comp(3) = 1.D0

                                Vec_phasefrac_3Comp = 0.D0

                                !Solve the system of equations
                                eqn = 3
                                call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                                if (errval /= 0) then
                                    phasefrac_VL1L2 = 1.D0 / 3.D0
                                else
                                    phasefrac_VL1L2 = Vec_rightside_3Comp(1:3)
                                end if
                            else
                                !Arbitrary start values
                                phasefrac_VL1L2 = 1.D0 / 3.D0
                            end if
                            call ptflash_NC_3P(gl,press, Temp, x, rho, x_vap, x_liq1, x_liq2, rhovap_est, &
                                & rholiq1_est, rholiq2_est, phasefrac_VL1L2, iFlash, iter, errval)
                            if (minval(phasefrac_VL1L2) < -1.D-16) then
                                errval = -3333
                            end if
                            if (errval == 0) then
                                !if (gl%N_guests > 1) then
                                !    if (.not. allocated(gl_check)) allocate(gl_check(3))
                                !    if (.not. allocated(ph_prop)) allocate(ph_prop(4))
                                !    gl_check = gl
                                !    Vec_rightside_3Comp(31:60) = 0.d0
                                !    eqn = 3
                                !    do j = 1,4!j=1:VLxH, 2: VLwH, 3: LxLwH
                                !        Vec_rightside_3Comp(1:30) = x
                                !        if (j == 1) then
                                !            ph_prop(j)%x_phase(:,1) = x_vap; ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                                !            ph_prop(j)%x_phase(:,2) = x_liq1; ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                                !            ph_prop(j)%text_phase(3) = 'hyd'; ph_prop(j)%phasetype(3) = 5
                                !        elseif (j == 2) then
                                !            ph_prop(j)%x_phase(:,1) = x_vap; ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                                !            ph_prop(j)%x_phase(:,2) = x_liq2; ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                                !            ph_prop(j)%text_phase(3) = 'hyd'; ph_prop(j)%phasetype(3) = 5
                                !        elseif (j == 3) then
                                !            ph_prop(j)%x_phase(:,1) = x_liq1; ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                                !            ph_prop(j)%x_phase(:,2) = x_liq2; ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                                !            ph_prop(j)%text_phase(3) = 'hyd'; ph_prop(j)%phasetype(3) = 5
                                !        elseif (j == 4) then
                                !            ph_prop(j)%x_phase(:,1) = x_vap; ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                                !            ph_prop(j)%x_phase(:,2) = x_liq1; ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                                !            ph_prop(j)%x_phase(:,3) = x_liq2; ph_prop(j)%text_phase(3) = 'liq2'; ph_prop(j)%phasetype(3) = 3
                                !            ph_prop(j)%rho = rho(1:3)
                                !            errval = 0
                                !        endif
                                !        if (j < 4) call ptflash_solid_NC_3P(gl_check(j),press, Temp, x, ph_prop(j)%rho, ph_prop(j)%x_phase(:,1), ph_prop(j)%x_phase(:,2), x_solid, ph_prop(j)%x_phase(:,3), rhovap_est, rholiq2_est, phasefrac, iFlash, iPhase_sol, iter, errval)
                                !        if (errval == 0) then
                                !            ph_prop(j)%errorflag = 0
                                !            Mat_phasefrac_3Comp(1,1:30) = ph_prop(j)%x_phase(:,1)
                                !            Mat_phasefrac_3Comp(2,1:30) = ph_prop(j)%x_phase(:,2)
                                !            Mat_phasefrac_3Comp(3,1:30) = ph_prop(j)%x_phase(:,3)
                                !            if (j < 4) then
                                !                call LUdecomp (gl_check(j),eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                                !            else
                                !                call LUdecomp (gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                                !            endif
                                !            ph_prop(j)%phasefrac = Vec_rightside_3Comp(1:3)
                                !            ph_prop(j)%gibbs_tot = 0.d0
                                !            do i = 1,3!i=1:phase1, 2:phase2, 3:phase3
                                !                if (j < 4) then
                                !                    gl_check(j)%molfractions = ph_prop(j)%x_phase(:,i)
                                !                    call reduced_parameters_calc(gl_check(j), 300.d0)
                                !                    if (i < 3) then
                                !                        ph_prop(j)%gibbs(i) = G_CALC(gl_check(j),Temp,ph_prop(j)%rho(i), 0)
                                !                        if (i == 1) call Chempot_CALC(gl_check(j),Temp,ph_prop(j)%rho(j), Chempot, 0)
                                !                    else
                                !                        call hdrt_Gibbs_energy(gl_check(j), ph_prop(j)%x_phase(:,3), Chempot, ph_prop(j)%gibbs(i), errval)
                                !                    endif
                                !                elseif (j == 4) then
                                !                    gl%molfractions = ph_prop(j)%x_phase(:,i)
                                !                    call reduced_parameters_calc(gl, 300.d0)
                                !                    ph_prop(j)%gibbs(i) = G_CALC(gl,Temp,rho(i), 0)
                                !                endif
                                !                ph_prop(j)%gibbs_tot = ph_prop(j)%gibbs_tot + ph_prop(j)%gibbs(i) * ph_prop(j)%phasefrac(i)
                                !            end do
                                !        else
                                !            ph_prop(j)%errorflag = errval
                                !            ph_prop(j)%gibbs_tot = 1.d60
                                !        endif
                                !    enddo
                                !endif
                                !if (errval == 0) then
                                !    if (gl%N_Guests > 1) then
                                !        j = minloc(ph_prop(:)%gibbs_tot,1)
                                !        phasetype(1:3) = ph_prop(j)%phasetype(1:3)
                                !        rho = 0.d0
                                !        rho(phasetype(1:3)) = ph_prop(j)%rho
                                !        x_Phase = 0.d0
                                !        x_Phase(:,phasetype(1)) = ph_prop(j)%x_phase(:,1)
                                !        x_Phase(:,phasetype(2)) = ph_prop(j)%x_phase(:,2)
                                !        x_Phase(:,phasetype(3)) = ph_prop(j)%x_phase(:,3)
                                !        Phasefrac = 0.d0
                                !        Phasefrac(phasetype(1:3)) = ph_prop(j)%phasefrac(1:3)
                                !        nrofphases = 3
                                !        if (j < 4) then
                                !            nrfluidphases = 2
                                !        elseif (j == 4) then
                                !            nrfluidphases = 3
                                !        endif
                                !        !Save information of the equilibrium found
                                !        gl%solidtype_akt_phase = 0             !No pure solids
                                !        gl%solidpos_akt_phase = 0              !No pure solids
                                !    else
                                !        phasetype(1) = 1
                                !        phasetype(2) = 2
                                !        phasetype(3) = 3
                                !        x_Phase(:,phasetype(1)) = x_vap
                                !        x_Phase(:,phasetype(2)) = x_liq1
                                !        x_Phase(:,phasetype(3)) = x_liq2
                                !        Phasefrac(phasetype(1)) = phasefrac_VL1L2(1)
                                !        Phasefrac(phasetype(2)) = phasefrac_VL1L2(2)
                                !        Phasefrac(phasetype(3)) = phasefrac_VL1L2(3)
                                !        nrofphases = 3
                                !        nrfluidphases = 3
                                !
                                !        !Save information of the equilibrium found
                                !        gl%solidtype_akt_phase = 0             !No pure solids
                                !        gl%solidpos_akt_phase = 0              !No pure solids
                                !    endif
                                !    return
                                !
                                !end if
                            end if
                            errval = 0

                        end if

                    end if

                end if

            end if

        end if

    end if
    use_old_stability = .false.
    if (.not.use_old_stability) then
        !########################################################################################################################################################
        !At this point only fluid phase equilibiria have been checked.
        !possible determined phase equlibiria are
        !1: V
        !2: L1
        !3: L2
        !4: VL1
        !5: VL2
        !6: L1L2
        !7: VL1L2
        !Possible phase equilibria with solid phases are:
        !   2 phase equilibria              3 phase equilibria
        ! 2 VHs1                          13  VL1Hs1
        ! 3 VHs2                          14  VL1Hs2
        ! 4 VIw                           15  VL1Iw
        ! 5 L1Hs1                         16  VL2Hs1
        ! 6 L1Hs2                         17  VL2Hs2
        ! 7 L1Iw                          18  VL2Iw
        ! 8 L2Hs1                         19  L1L2Hs1
        ! 9 L2Hs2                         20  L1L2Hs2
        !10 L2Iw                          21  L1L2Iw
        !11 Hs1Iw                         22  VHs1Iw
        !12 Hs2Iw                         23  VHs2Iw
        !                                 24  L1Hs1Iw
        !28 VIc    (if co2 is present)    25  L1Hs2Iw
        !29 L1Ic   (if co2 is present)    26  L2Hs1Iw
        !30 L2Ic   (if co2 is present)    27  L2Hs2Iw
        !31 Hs1Ic  (if co2 is present)
        !32 Hs2Ic  (if co2 is present)    33  VL1Ic    (if co2 is present)
        !                                 34  VL2Ic    (if co2 is present)
        !                                 35  L1L2Ic   (if co2 is present)
        !                                 36  VHs1Ic   (if co2 is present)
        !                                 37  VHs2Ic   (if co2 is present)
        !                                 38  L1Hs1Ic  (if co2 is present)
        !                                 39  L1Hs2Ic  (if co2 is present)
        !                                 40  L2Hs1Ic  (if co2 is present)
        !                                 41  L2Hs2Ic  (if co2 is present)
        !   27  Combinations WITHOUT CO2
        !   41  Combinations WITH CO2

        co2_present = .false.
        h2o_present = .false.
        if (count(gl%components == 'water') /= 0) h2o_present = .true.
        if (count(gl%components == 'co2') /= 0) co2_present = .true.

        if (.not. co2_present) then
            nr_checks = 26 + 1
        elseif (co2_present) then
            nr_checks = 40 + 1
        endif
        if(.not. c_associated(gl%s_p%ph_prop_handle)) then !create a new fluid and save the handle
            allocate(ph_prop(nr_checks))
            gl%s_p%ph_prop_handle = C_LOC(ph_prop)
            gl%s_p%nr_checks = nr_checks
            gl%s_p%fluid_present = .false.
        elseif (c_associated(gl%s_p%ph_prop_handle)) then !create a new fluid and save the handle
            CALL C_F_POINTER(gl%s_p%ph_prop_handle,ph_prop,[nr_checks])
        endif
        if (.not. gl%s_p%fluid_present) then

            if(.not. c_associated(gl%s_p%gl_check_handle)) then !create a new fluid and save the handle
                allocate(gl_check)
                gl%s_p%gl_check_handle = C_LOC(gl_check)
            elseif (c_associated(gl%s_p%gl_check_handle)) then !create a new fluid and save the handle
                CALL C_F_POINTER(gl%s_p%gl_check_handle,gl_check)
            endif
            !initialize variable
            !1st element
            ph_prop(1)%x_phase = 0.d0
            ph_prop(1)%gibbs = 0.d0
            ph_prop(1)%phasefrac = 0.d0
            ph_prop(1)%rho = 0.d0
            ph_prop(1)%gibbs_tot = 0.d0
            ph_prop(1)%text_phase = ''
            ph_prop(1)%phasetype = 0
            ph_prop(1)%errorflag = 0.d0
            ph_prop(1)%nrofphases = 0.d0
            ph_prop(1)%hdrt_structure_flag = 1
            ph_prop(1)%solidtype = 0
            !remaining elements
            ph_prop(:) = type_phase_properties(ph_prop(1)%x_phase, ph_prop(1)%gibbs, ph_prop(1)%phasefrac, ph_prop(1)%rho, ph_prop(1)%gibbs_tot, ph_prop(1)%text_phase, ph_prop(1)%phasetype, ph_prop(1)%errorflag, ph_prop(1)%nrofphases, ph_prop(1)%hdrt_structure_flag, ph_prop(1)%solidtype)
            !initialize start values
            !set initial estimates for the phase compositions in case the stability analysis for fluid phases failed
            if (errval /= 0) then
                errval_fluid = errval
                errval = 0
                call PTX_startvals_pT (gl,press, Temp, x, x_phase(:,1), x_phase(:,3), vapfrac, errval)
                if (h2o_present) then
                    x_phase(1,2) = 1.d-7
                    x_phase(2:30,2) = x(2:30)/(1.d0 - x_phase(1,2))
                else
                    x_phase(:,2) = x_phase(:,3)
                endif
            endif
            if ((dabs(sum(x_phase(:,1)) - 1.d0) > 1.d-14).or.(count(phasetype(:) == 1) == 0)) then
                call PTX_startvals_pT (gl,press, Temp, x, x_phase(:,1), x_liq1, vapfrac, errval)
            endif
            if ((dabs(sum(x_phase(:,2)) - 1.d0) > 1.d-14).or.(count(phasetype(:) == 2) == 0)) then
                call PTX_startvals_pT (gl,press, Temp, x, x_vap, x_phase(:,2), vapfrac, errval)
            endif
            if ((dabs(sum(x_phase(:,3)) - 1.d0) > 1.d-14).or.(count(phasetype(:) == 3) == 0)) then
                call PTX_startvals_pT (gl,press, Temp, x, x_vap, x_phase(:,3), vapfrac, errval)
            endif
            !end initialize


            if (.not. h2o_present) ph_prop(2:27)%gibbs_tot = 1.d60
            Vec_rightside_3Comp(31:60) = 0.d0
            gl%molfractions = x
            do j = 1,nr_checks!j=1:VLxH, 2: VLwH, 3: LxLwH
                gl_check = gl

                Vec_rightside_3Comp(1:30) = x
                select case (j)
                case (1) !V, L1, L2, VL1, VL2, L1L2, VL1L2
                    if (errval_fluid == 0) then
                        if (nrfluidphases == 1) then
                            x_fluid1 = 0.d0
                            ph_prop(j)%x_phase(:,1) = x_phase(:,phasetype(1));
                            if (phasetype(1) == 1) ph_prop(j)%text_phase(1) = 'vap';
                            if (phasetype(1) == 2) ph_prop(j)%text_phase(1) = 'liq1';
                            if (phasetype(1) == 3) ph_prop(j)%text_phase(1) = 'liq2';
                            ph_prop(j)%phasetype(1) = phasetype(1);
                            ph_prop(j)%phasefrac(1) = phasefrac(phasetype(1))
                            ph_prop(j)%rho(1) = rho(phasetype(1))
                        elseif (nrfluidphases == 2) then
                            x_fluid1 = 0.d0
                            x_fluid2 = 0.d0
                            ph_prop(j)%x_phase(:,1) = x_phase(:,phasetype(1));
                            if (phasetype(1) == 1) ph_prop(j)%text_phase(1) = 'vap';
                            if (phasetype(1) == 2) ph_prop(j)%text_phase(1) = 'liq1';
                            ph_prop(j)%phasetype(1) = phasetype(1);
                            ph_prop(j)%x_phase(:,2) = x_phase(:,phasetype(2));
                            if (phasetype(2) == 2) ph_prop(j)%text_phase(2) = 'liq1';
                            if (phasetype(2) == 3) ph_prop(j)%text_phase(2) = 'liq2';
                            ph_prop(j)%phasetype(2) = phasetype(2)
                            ph_prop(j)%rho(1:2) = rho(phasetype(1:2))
                            ph_prop(j)%phasefrac(1:2) = phasefrac(phasetype(1:2))
                            ph_prop(j)%rho(1:2) = rho(phasetype(1:2))
                        elseif (nrfluidphases == 3) then
                            x_fluid1 = 0.d0
                            x_fluid2 = 0.d0
                            ph_prop(j)%x_phase(:,1) = x_vap; ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                            ph_prop(j)%x_phase(:,2) = x_liq1; ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                            ph_prop(j)%x_phase(:,3) = x_liq2; ph_prop(j)%text_phase(3) = 'liq2'; ph_prop(j)%phasetype(3) = 3
                            ph_prop(j)%rho = rho(1:3)
                            ph_prop(j)%phasefrac = phasefrac(1:3)
                            ph_prop(j)%rho = rho(1:3)
                        endif
                        ph_prop(j)%nrofphases = nrfluidphases
                        errval = 0
                    else
                        ph_prop(j)%gibbs_tot = 1.d60
                        ph_prop(j)%errorflag = errval_fluid
                        cycle
                    endif
                case (2) !VHs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    ph_prop(j)%text_phase(2) = 'hyds1'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor

                case (3) !VHs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    ph_prop(j)%text_phase(2) = 'hyds2'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor

                case (4) !VIw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor

                case (5) !L1Hs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2 !x_phase(:,2) before, BS 07/2019
                    ph_prop(j)%text_phase(2) = 'hyds1'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (6) !L1Hs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    ph_prop(j)%text_phase(2) = 'hyds2'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (7) !L1Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    ph_prop(j)%nrofphases = 2
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid

                case (8) !L2Hs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    ph_prop(j)%text_phase(2) = 'hyds1'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (9) !L2Hs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    ph_prop(j)%text_phase(2) = 'hyds2'; ph_prop(j)%phasetype(2) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (10) !L2Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid

                case (11) !HIw - will not work since no fugacity of guests can be calculated
                    !if (.not. h2o_present) cycle
                    ph_prop(j)%gibbs_tot = 1.d60
                    ph_prop(j)%errorflag = -15000
                    cycle
                    !gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    !x_solid = 0.d0
                    !x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    !ph_prop(j)%text_phase(2) = 'hyds1'; ph_prop(j)%phasetype(2) = 5
                    !ph_prop(j)%hdrt_structure_flag = 1
                    !ph_prop(j)%nrofphases = 2
                    !errval = -15566
                    !gl_check%solidtype(1) = 1
                    !gl_check%solidtype(2) = 1

                case (12) !HIw - will not work since no fugacity of guests can be calculated
                    !if (.not. h2o_present) cycle
                    ph_prop(j)%gibbs_tot = 1.d60
                    ph_prop(j)%errorflag = -15000
                    cycle
                    !gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    !x_solid = 0.d0
                    !x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    !ph_prop(j)%text_phase(2) = 'hyds2'; ph_prop(j)%phasetype(2) = 5
                    !ph_prop(j)%hdrt_structure_flag = 2
                    !ph_prop(j)%nrofphases = 2
                    !errval = -15566
                    !gl_check%solidtype(1) = 1
                    !gl_check%solidtype(2) = 1

                case (13) !VL1Hs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,2); ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor/liquid

                case (14) !VL1Hs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,2); ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor/liquid

                case (15) !VL1Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    ph_prop(j)%nrofphases = 3
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,2); ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor/liquid

                case (16) !VL2Hs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor/liquid

                case (17) !VL2Hs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor/liquid

                case (18) !VL2Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor/liquid

                case (19) !L1L2Hs1
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid/liquid

                case (20) !L1L2Hs2
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 0
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid/liquid

                case (21) !L1L2Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid/liquid

                case (22) !VHs1Iw
                    if (.not. h2o_present) cycle
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    iPhase = 2 !assume vapor

                case (23) !VHs2Iw
                    if (.not. h2o_present) cycle
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    iPhase = 2 !assume vapor

                case (24) !L1Hs1Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (25) !L1Hs2Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (26) !L2Hs1Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (27) !L2Hs2Iw
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("water",gl_check%Fluidlist_hydrate)
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(1) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 1
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid


                case (28) !VIc
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor

                case (29) !L1Ic
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid

                case (30) !L2Ic
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%nrofphases = 2
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid

                case (31) !Hs1Ic - will not work since no fugacity of guests can be calculated
                    ph_prop(j)%gibbs_tot = 1.d60
                    ph_prop(j)%errorflag = -15000
                    cycle
                    !if (.not. co2_present) cycle
                    !if (gl_check%ncomp > 2) cycle
                    !gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    !x_solid = 0.d0
                    !x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    !ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    !ph_prop(j)%hdrt_structure_flag = 1
                    !ph_prop(j)%nrofphases = 2
                    !gl_check%solidtype(1) = 2
                    !gl_check%solidtype(2) = 1

                case (32) !Hs2Ic - will not work since no fugacity of guests can be calculated
                    ph_prop(j)%gibbs_tot = 1.d60
                    ph_prop(j)%errorflag = -15000
                    cycle
                    !if (.not. co2_present) cycle
                    !if (gl_check%ncomp > 2) cycle
                    !gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    !x_solid = 0.d0
                    !x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    !ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    !ph_prop(j)%hdrt_structure_flag = 2
                    !ph_prop(j)%nrofphases = 2
                    !gl_check%solidtype(1) = 2
                    !gl_check%solidtype(2) = 1

                case (33) !VL1Ic
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(2) = 'liq1'; ph_prop(j)%phasetype(2) = 2
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor/liquid

                case (34) !VL2Ic
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 2 !assume vapor/liquid

                case (35) !L1L2Ic
                    if (.not. co2_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_fluid2 = x_phase(:,3); ph_prop(j)%text_phase(2) = 'liq2'; ph_prop(j)%phasetype(2) = 3
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(3) = 'sol'; ph_prop(j)%phasetype(3) = 4
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 0
                    iPhase = 1 !assume liquid/liquid

                case (36) !VHs1Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor

                case (37) !VHs2Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,1); ph_prop(j)%text_phase(1) = 'vap'; ph_prop(j)%phasetype(1) = 1
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 2 !assume vapor

                case (38) !L1Hs1Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (39) !L1Hs2Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,2); ph_prop(j)%text_phase(1) = 'liq1'; ph_prop(j)%phasetype(1) = 2
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (40) !L2Hs1Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds1'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 1
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid

                case (41) !L2Hs1Ic
                    if (.not. co2_present) cycle
                    if (.not. h2o_present) cycle
                    gl_check%solid_pos = locate_entry("co2",gl_check%Fluidlist_hydrate)
                    x_fluid1 = x_phase(:,3); ph_prop(j)%text_phase(1) = 'liq2'; ph_prop(j)%phasetype(1) = 3
                    x_solid = 0.d0
                    x_solid(gl_check%solid_pos) = 1.d0; ph_prop(j)%text_phase(2) = 'sol'; ph_prop(j)%phasetype(2) = 4
                    ph_prop(j)%text_phase(3) = 'hyds2'; ph_prop(j)%phasetype(3) = 5
                    ph_prop(j)%hdrt_structure_flag = 2
                    ph_prop(j)%nrofphases = 3
                    gl_check%solidtype(1) = 2
                    gl_check%solidtype(2) = 1
                    iPhase = 1 !assume liquid
                    case default
                    errval = -15566
                end select
                ph_prop(j)%solidtype = gl_check%solidtype
                if (any(ph_prop(j)%phasetype(:) == 5)) then
                    call hdrt_structure_definition(gl_check,ph_prop(j)%hdrt_structure_flag)
                endif
                errval = 0
                if ((j /= 1).and.(ph_prop(j)%nrofphases == 2)) then!no phase equilibrium calc if only fluid phase is present
                    iFlash = 3
                    rhofluid_est = 0.d0
                    call ptflash_solid_NC_2P(gl_check,press, Temp, ph_prop(j)%rho, x, x_solid, x_hyd, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval, iter)
                elseif ((j /= 1).and.(ph_prop(j)%nrofphases == 3)) then
                    iFlash = 7
                    rhofluid_est = 0.d0
                    call ptflash_solid_NC_3P(gl_check,press, Temp, x, ph_prop(j)%rho, x_fluid1, x_fluid2, x_solid, x_hyd, rhovap_est, rholiq2_est, phasefrac, iFlash, iPhase, iter, errval)
                endif
                !check whether fluid phases have the given fluid state (i.e. vapor given vapor found, liquid given, liquid found)
                if ((any(ph_prop(j)%Phasetype(:) <= 3)).and.(errval == 0)) then
                    !deriv < 0 -> vapor, deriv > 0 -> liquid
                    do i = 1,ph_prop(j)%nrofphases
                        if (ph_prop(j)%phasetype(i) <= 3) then
                            gl_check%molfractions = x_phase(:,ph_prop(j)%phasetype(i))
                            call reduced_parameters_calc(gl_check, temp)
                            curvature = D2PDD2_CALC(gl_check,Temp, ph_prop(j)%rho(i), 0)*1.d6
                            gl_check%molfractions = x
                            call reduced_parameters_calc(gl_check, temp)
                            if ((ph_prop(j)%phasetype(i) == 1).and.(curvature > 0)) then!vapor
                                errval = -15566
                                exit
                            elseif (((ph_prop(j)%phasetype(i) == 2) .or. (ph_prop(j)%phasetype(i) == 3)) .and. (curvature < 0)) then!liquid
                                errval = -15566
                                exit
                            endif
                        endif
                    enddo
                endif
                if (errval == 0) then
                    if ((ph_prop(j)%nrofphases == 1).and.(j /= 1)) then
                        !no case except for 0 where this happens
                    elseif ((ph_prop(j)%nrofphases == 2).and.(j /= 1)) then
                        if ((ph_prop(j)%phasetype(1) <= 3).and.(ph_prop(j)%phasetype(2) == 4)) then!14, 24, 34
                            ph_prop(j)%x_phase(:,1) = x_fluid1
                            ph_prop(j)%x_phase(:,2) = x_solid
                        elseif ((ph_prop(j)%phasetype(1) <= 3).and.(ph_prop(j)%phasetype(2) == 5)) then!15, 25, 35
                            ph_prop(j)%x_phase(:,1) = x_fluid1
                            ph_prop(j)%x_phase(:,2) = x_hyd
                        elseif ((ph_prop(j)%phasetype(1) == 4).and.(ph_prop(j)%phasetype(2) == 5)) then!45 !not possible but for the sake of completeness
                            ph_prop(j)%x_phase(:,1) = x_solid
                            ph_prop(j)%x_phase(:,2) = x_hyd
                        endif
                        ph_prop(j)%errorflag = 0
                        ph_prop(j)%phasefrac(1) = vapfrac1
                        ph_prop(j)%phasefrac(2) = 1.d0 - vapfrac1
                        ph_prop(j)%phasefrac(3) = 0.d0

                    elseif ((ph_prop(j)%nrofphases == 3).and.(j /= 1)) then
                        if ((ph_prop(j)%phasetype(1) <= 3).and.(ph_prop(j)%phasetype(2) <= 3).and.(ph_prop(j)%phasetype(3) == 4)) then!124,134,234
                            ph_prop(j)%x_phase(:,1) = x_fluid1
                            ph_prop(j)%x_phase(:,2) = x_fluid2
                            ph_prop(j)%x_phase(:,3) = x_solid
                        elseif ((ph_prop(j)%phasetype(1) <= 3).and.(ph_prop(j)%phasetype(2) <= 3).and.(ph_prop(j)%phasetype(3) == 5)) then!125,135,235
                            ph_prop(j)%x_phase(:,1) = x_fluid1
                            ph_prop(j)%x_phase(:,2) = x_fluid2
                            ph_prop(j)%x_phase(:,3) = x_hyd
                        elseif ((ph_prop(j)%phasetype(1) <= 3).and.(ph_prop(j)%phasetype(2) == 4).and.(ph_prop(j)%phasetype(3) == 5)) then!145,245,345
                            ph_prop(j)%x_phase(:,1) = x_fluid1
                            ph_prop(j)%x_phase(:,2) = x_solid
                            ph_prop(j)%x_phase(:,3) = x_hyd
                        endif
                        ph_prop(j)%errorflag = 0
                        Mat_phasefrac_3Comp = 0.d0
                        Mat_phasefrac_3Comp(1,1:30) = ph_prop(j)%x_phase(:,1)
                        Mat_phasefrac_3Comp(2,1:30) = ph_prop(j)%x_phase(:,2)
                        Mat_phasefrac_3Comp(3,1:30) = ph_prop(j)%x_phase(:,3)
                        Vec_rightside_3Comp(1:30) = x
                        eqn = ph_prop(j)%nrofphases
                        call LUdecomp (gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                        ph_prop(j)%phasefrac = Vec_rightside_3Comp(1:ph_prop(j)%nrofphases)

                    elseif (j == 1) then
                        continue
                    else
                        ph_prop(j)%errorflag = -15000
                        cycle
                    endif

                    if ( (minval(ph_prop(j)%phasefrac(1:ph_prop(j)%nrofphases)) < 0.d0) .or. (maxval(ph_prop(j)%phasefrac(1:ph_prop(j)%nrofphases)) > 1.0d0)) then
                        ph_prop(j)%gibbs_tot = 1.d60
                        ph_prop(j)%errorflag = -15000
                    else
                        ph_prop(j)%gibbs_tot = 0.d0
                        !this is console output for debugging
                        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
                        if (ph_prop(j)%nrofphases == 2) then
                            write(*,*)"z(1)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(1,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(1,2)
                            write(*,*)"z(2)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(2,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(2,2)
                            write(*,*)"z(3)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(3,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(3,2)
                        elseif (ph_prop(j)%nrofphases == 3) then
                            write(*,*)"z(1)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(1,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(1,2) + ph_prop(j)%phasefrac(3)*ph_prop(j)%x_phase(1,3)
                            write(*,*)"z(2)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(2,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(2,2) + ph_prop(j)%phasefrac(3)*ph_prop(j)%x_phase(2,3)
                            write(*,*)"z(3)=",ph_prop(j)%phasefrac(1)*ph_prop(j)%x_phase(3,1) + ph_prop(j)%phasefrac(2)*ph_prop(j)%x_phase(3,2) + ph_prop(j)%phasefrac(3)*ph_prop(j)%x_phase(3,3)
                        endif
                        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
                        do i = 1,ph_prop(j)%nrofphases!i=1:phase1, 2:phase2, 3:phase3
                            if (j /= 1) then
                                gl_check%molfractions = ph_prop(j)%x_phase(:,i)
                                call reduced_parameters_calc(gl_check, 300.d0)
                                if (ph_prop(j)%phasetype(i) <= 3) then
                                    ph_prop(j)%gibbs(i) = G_CALC(gl_check,Temp,ph_prop(j)%rho(i), 0)
                                    if (i == 1) call Chempot_CALC(gl_check,Temp,ph_prop(j)%rho(i), Chempot, 0)
                                elseif (ph_prop(j)%phasetype(i) == 4) then
                                    if (gl_check%solidtype(1) == 2) then        !solidtype(1) = 2 --> Dry Ice
                                        ph_prop(j)%gibbs(i) = g_DryIce(gl_check,Temp, press)    !Chemical Potential of Dry Ice [J / mol]
                                    End if
                                    if (gl_check%solidtype(1) == 1) then    !solidtype(1) = 1 --> Solid water
                                        ph_prop(j)%gibbs(i) = g_WaterIce(gl_check,Temp, press)  !Chemical Potential of Water Ice [J / mol]
                                    End if
                                elseif (ph_prop(j)%phasetype(i) == 5) then
                                    call hdrt_Gibbs_energy(gl_check, ph_prop(j)%x_phase(:,i), Chempot, ph_prop(j)%gibbs(i), errval)
                                endif
                            elseif (j == 1) then
                                gl_check%molfractions = ph_prop(j)%x_phase(:,i)
                                call reduced_parameters_calc(gl_check, 300.d0)
                                ph_prop(j)%gibbs(i) = G_CALC(gl_check,Temp,ph_prop(j)%rho(i), 0)
                            endif
                        end do
                        ph_prop(j)%gibbs_tot = sum(ph_prop(j)%gibbs(1:ph_prop(j)%nrofphases) * ph_prop(j)%phasefrac(1:ph_prop(j)%nrofphases))
                        if (j == 10) then!if L1Iw and L2Iw are the same results take L2Iw
                            if (dabs(ph_prop(j)%gibbs_tot - ph_prop(7)%gibbs_tot) < 1.d-12) ph_prop(7)%gibbs_tot = 1.d60
                        elseif (j == 30) then !if L1Ic and L2Ic are the same results take L2Ic
                            if (dabs(ph_prop(j)%gibbs_tot - ph_prop(j-1)%gibbs_tot) < 1.d-12) ph_prop(j-1)%gibbs_tot = 1.d60
                        endif
                    endif
                else
                    ph_prop(j)%errorflag = errval
                    ph_prop(j)%gibbs_tot = 1.d60
                endif
            enddo
        endif
        if (any(ph_prop(:)%errorflag == 0)) then
            errval = 0
            j = minloc(ph_prop(:)%gibbs_tot,1)
            phasetype(1:ph_prop(j)%nrofphases) = ph_prop(j)%phasetype(1:ph_prop(j)%nrofphases)
            rho = 0.d0
            rho(phasetype(1:ph_prop(j)%nrofphases)) = ph_prop(j)%rho(1:ph_prop(j)%nrofphases)
            x_Phase = 0.d0
            x_Phase(:,phasetype(1:ph_prop(j)%nrofphases)) = ph_prop(j)%x_phase(:,1:ph_prop(j)%nrofphases)
            !x_Phase(:,phasetype(2)) = ph_prop(j)%x_phase(:,2)
            !x_Phase(:,phasetype(3)) = ph_prop(j)%x_phase(:,3)
            Phasefrac = 0.d0
            Phasefrac(phasetype(1:ph_prop(j)%nrofphases)) = ph_prop(j)%phasefrac(1:ph_prop(j)%nrofphases)
            nrofphases = ph_prop(j)%nrofphases

            !Save information of the equilibrium found
            if (any(phasetype(:) == 4)) then
                gl%solidpos_akt_phase = locate_entry(trim(gl%components(maxloc(x_phase(:,4),1))),gl_check%Fluidlist_hydrate)
                gl%solidtype_akt_phase = ph_prop(j)%solidtype(1)
            endif
            !Save information of the equilibrium found
            if (any(phasetype(:) == 5)) then
                gl%hdrt_structure_stable = ph_prop(j)%hdrt_structure_flag
            endif
            !if (allocated(ph_prop))deallocate(ph_prop)
            !if(allocated(gl_check))deallocate(gl_check)
        endif
        !endif
    elseif (use_old_stability) then
        !########################################################################################################################################################
        !BEGIN OF THE OLD ROUTINE:
        !If an error occurs in the phase determination for the fluid phases, try fluid phase (gas or H2O poor liquid) in equilibrium with hydrate or solid H2O
        !and solid H2O / Hydrate.
        !In case CO2 is present in the mixture, try dry ice / hydrate too
        !Hence, the following phase equilibria need to be checked:
        !Two phase:     VIw, LxIw --> LxH, VH --> IwH, IcH, (VIc, LcIc)
        !Three phase:   VLxIw --> VLxH --> VHIw,LHIw (VLxIc)
        !The following are possible, but excluded (the error is most likely due to extrapolation problems of the water equation --> solid would form. So liquid water is not checked.
        !Two phase:     VLw, LxLw, LwH, LwIc
        !Three phase:   VLwH, LxLwH
        if (errval /= 0) then
            errval_fluid = errval

            !------------------------------------------------------------------------------
            gl%solid_pos = 0
            !Check whether a fluid phase is in equilibrium with solid water
            Do i = 1, gl%ncomp
                if (gl%Fluidlist_hydrate(i) == "water") then
                    gl%solid_pos = i
                    exit
                end if
            end do
            if (gl%solid_pos == 0) return  !Quit if no water is in the mixture

            if (Temp < 273.2D0) then    !No solid water above this temperature assumed, quit with error elsewise

                !1.1) Try to calculate a vapor / solid water equilibrium first
                x_fluid(1) = 0.000001D0    !Arbitrary startvalue
                do i = 2, gl%ncomp
                    x_fluid(i) = 0.999999D0  * x(i) / (1.D0 - x(1))        !Arbitrary startvalue
                end do
                x_fluid1 = x_fluid
                x_fluid2 = x_fluid
                x_solid = 0.D0
                x_solid(gl%solid_pos) = 1.d0
                gl%solidtype(1) = 1   !solid forms
                gl%solidtype(2) = 0   !No hydrate assumed
                tpd_sol_phase1 = 1.D6
                tpd_sol_phase2 = 1.D6
                vapfrac1 = 0
                vapfrac2 = 0
                errval1 = 0
                errval2 = 0
                iFlash = 3
                rhofluid_est = 0.D0

                !Check VIw Equilibrium
                iphase = 2 !Vapor phase assumed
                call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval1, iter)
                gVIw = 1.D12
                rho_VS_V = 0.D0
                rho_VS_S = 0.D0
                !If a physically correct phase equilibrium was found, calculate the Gibbs-energy
                if (errval1 == 0) then
                    if ((vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then  !Physically plausible solution
                        rho_VS_V = rho_sol(1)
                        rho_VS_S = rho_sol(2)
                        gl%molfractions = x_fluid1
                        call reduced_parameters_calc(gl,Temp)
                        gf1 = G_CALC(gl,Temp, rho_VS_V, 0)
                        gf2 = g_WaterIce(gl,Temp, press)
                        vapfrac1 = (x(1)-x_solid(1)) / (x_fluid1(1)-x_solid(1))
                        gVIw = vapfrac1 * gf1 + (1.D0 - vapfrac1) * gf2
                        !Test if first phase is really vapor
                        curvature = d2P_drho2 (gl,Temp, rho_VS_V)
                        if (curvature > 0.D0) then  !If liquid, discard solution
                            gVIw = 1.D12
                        end if
                        !Set back reducing parameters
                        gl%molfractions = x
                        gl%rhoredmix = rhoredmix_orig
                        gl%tredmix = tredmix_orig
                    End if
                end if

                !Check LIw equilibrium
                iphase = 1 !Liquid phase assumed
                call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid2, rhofluid_est, vapfrac2, iFlash, iPhase, errval2, iter)
                gLIw = 1.D12
                !If a physically correct phase equilibrium was found, calculate the Gibbs-energy
                if (errval2 == 0) then
                    if ((vapfrac2 > 0.D0) .and. (vapfrac2 < 1.D0)) then  !Physically plausible solution
                        rho_LS_L = rho_sol(1)
                        rho_LS_S = rho_sol(2)
                        gl%molfractions = x_fluid2
                        call reduced_parameters_calc(gl,Temp)
                        gf1 = G_CALC(gl,Temp, rho_LS_L, 0)
                        gf2 = g_WaterIce(gl,Temp, press)
                        vapfrac2 = (x(1)-x_solid(1)) / (x_fluid2(1)-x_solid(1))
                        gLIw = vapfrac2 * gf1 + (1.D0 - vapfrac2) * gf2
                        !Test if first phase is really vapor
                        curvature = d2P_drho2 (gl,Temp, rho_LS_L)
                        if (curvature < 0.D0) then  !If vapor, discard solution
                            gLIw = 1.D12
                        end if
                        !Set back reducing parameters
                        gl%molfractions = x
                        gl%rhoredmix = rhoredmix_orig
                        gl%tredmix = tredmix_orig
                    End if
                end if

                !Three cases may occur:
                !1) Only one physically reasonable equilibrium found (VIw or LIw)
                !2) VIw AND LIw are physically reasonable --> find the stable one
                !3) No reasonable equilibrium found --> only solid phases left (HIw, maybe HIc)
                if ((gVIw < 0.9D12) .or. (gLIw < 0.9D12)) then       !Physically plausible VIw or/and LIw equilibrium found
                    !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                    !Save information that fluid phase in equilibrium with solid water was found
                    !Andresa July 2014
                    error_then_ice_found = .true.
                    if (gVIw < gLIw) then  !Physically plausible solution, Vapor/Ice is more stable than Liquid/Ice.
                        tpd_sol_phase1 = 1.D0
                        phasetype(1) = 1
                        rho(phasetype(1)) = rho_VS_V
                        x_phase(:,phasetype(1)) = x_fluid1
                        phasetype(2) = 4
                        rho(phasetype(2)) = rho_VS_S
                        x_phase(:,phasetype(2)) = x_solid
                        nrofphases = 2
                        errval = 0
                        Phasefrac(phasetype(1)) = vapfrac1
                        Phasefrac(phasetype(2)) = 1.D0 - Phasefrac(phasetype(1))
                        nrfluidphases = 1

                        !Save information of the equilibrium found
                        gl%solidtype_akt_phase = 1             !Solid is H2O
                        gl%solidpos_akt_phase = gl%solid_pos      !Position of H2O in the fluid vector
                    else
                        tpd_sol_phase1 = 1.D0
                        phasetype(1) = 3
                        rho(phasetype(1)) = rho_LS_L
                        x_phase(:,phasetype(1)) = x_fluid2
                        phasetype(2) = 4
                        rho(phasetype(2)) = rho_LS_S
                        x_phase(:,phasetype(2)) = x_solid
                        nrofphases = 2
                        errval = 0
                        Phasefrac(phasetype(1)) = vapfrac2
                        Phasefrac(phasetype(2)) = 1.D0 - Phasefrac(phasetype(1))
                        nrfluidphases = 1

                        !Save information of the equilibrium found
                        gl%solidtype_akt_phase = 1             !Solid is H2O
                        gl%solidpos_akt_phase = gl%solid_pos      !Position of H2O in the fluid vector
                    end if
                else
                    tpd_sol_phase1 = 1.D6
                end if

                gl%molfractions = x
                call reduced_parameters_calc(gl,Temp)

                !In principle, it would make sense to check for VLwIw or VLx equilibrium here. However, VLw is unlikely, since no two phase equilibrium was found.
                !But VLxIw must be checked.
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !THREE PHASE CHECK VLIw
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                If ((gVIw < 0.9D12) .and. (gLIw < 0.9D12) .and. (gl%ncomp > 2)) then
                    if (dabs(gVIw - gLIw) < 10.D0) then !The Gibbs energies are almost the same, check for three phase equilibrium

                        x_vap = x_fluid1
                        x_liq2 = x_fluid2

                        !Setup for the three phase flash
                        gl%solidtype(1) = 1 !Solid H2O
                        gl%solidtype(2) = 0 !No hydrate forms

                        !VLIw assumed, thus lighter phase is vapor
                        iPhase_sol = 2 !Liqhter phase is vapor

                        iFlash = 7  !T and p given

                        rhovap_est = 0.d0
                        rholiq2_est = 0.D0

                        !Startvalue for phasefractions
                        phasefrac_3Phase_solid = 1.D0 / 3.D0

                        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_vap, x_liq2, x_solid, x_hyd, rhovap_est, &
                            & rholiq2_est, phasefrac_3Phase_solid, iFlash, iPhase_sol, iter, errval)

                        if ((gl%ncomp == 3) .and. (errval == 0)) then

                            Mat_phasefrac_3Comp(1,1) = x_vap(1)
                            Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                            Mat_phasefrac_3Comp(3,1) = x_solid(1)

                            Mat_phasefrac_3Comp(1,2) = x_vap(2)
                            Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                            Mat_phasefrac_3Comp(3,2) = x_solid(2)

                            Mat_phasefrac_3Comp(1,3) = 1.D0
                            Mat_phasefrac_3Comp(2,3) = 1.D0
                            Mat_phasefrac_3Comp(3,3) = 1.D0

                            Vec_rightside_3Comp(1) = x(1)
                            Vec_rightside_3Comp(2) = x(2)
                            Vec_rightside_3Comp(3) = 1.D0

                            Vec_phasefrac_3Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                            phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)

                        end if

                        if (minval(phasefrac_3Phase_solid) < -1.D-16) then
                            errval = -3333
                        end if
                        if (errval == 0) then
                            phasetype(1) = 1
                            phasetype(2) = 3
                            phasetype(3) = 4
                            x_Phase(:,phasetype(1)) = x_vap
                            x_Phase(:,phasetype(2)) = x_liq2
                            x_Phase(:,phasetype(3)) = x_solid
                            Phasefrac(phasetype(1)) = phasefrac_3Phase_solid(1)
                            Phasefrac(phasetype(2)) = phasefrac_3Phase_solid(2)
                            Phasefrac(phasetype(3)) = phasefrac_3Phase_solid(3)
                            rho(phasetype(1)) = rho_sol(1)
                            rho(phasetype(2)) = rho_sol(2)
                            rho(phasetype(3)) = rho_sol(3)
                            nrofphases = 3

                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 1             !Solid is H2O
                            gl%solidpos_akt_phase = gl%solid_pos      !Position of H2O in the fluid vector
                            return
                        else
                            errval = 0
                        end if
                    end if
                end if
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !END OF THREE PHASE CHECK VLIw
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !1.2) Check whether hydrates may form. If yes, hydrate / solid CO2 and hydrate / solid water is possible and has to be checked
                if ((gl%nrofhydrateformers > 0) .and. (Temp < 273.2D0)) then    !Check if hydrates are possible

                    !1.2.1) If a vapor / solid H2O equilibrium was calculated, check the tpd for the correct equilbrium, try if hydrate / solid H2O leads to a lower tpd
                    !if (dabs(tpd_sol_phase1 - 1.D0) < 1.D-12) then  !This means that an equilibrium with vapor or liquid and solid H2O was found
                    if (error_then_ice_found) then  !This means that an equilibrium with vapor or liquid and solid H2O was found   !Andreas July 2014

                        !Andreas March 2015: PROCEED WITH "NORMAL ROUTINE" but first check if solid co2 forms. If yes, co2 + hydrate will be the correct equilibrium.
                        if ((gl%ncomp == 2) .and. (gl%Fluidlist_hydrate(2) == "co2")) then
                            ! Get the melting or sublimation pressure at the given temperature
                            gl%solid_pos = 2
                            if (Temp > gl%ttp(gl%solid_pos)) then
                                peq = pmelt_eq(gl,Temp, gl%solid_pos)
                            else
                                peq = psub_eq(gl,Temp, gl%solid_pos)
                            end if
                            !if (dexp(lnf(solid_pos)) > (0.9D0 *peq)) then
                            if (press > (0.9D0 *peq)) then

                                gl%molfractions = x_phase(:,phasetype(1))
                                call reduced_parameters_calc(gl,Temp)
                                call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
                                gl%molfractions = x
                                call reduced_parameters_calc(gl,Temp)

                                tpd_sol = g_DryIce(gl,Temp,press) - Chempot(gl%solid_pos)


                                if (tpd_sol < 0.D0) then
                                    !Dry ice forms, hydrate + dry ice is the correct phase equilibrium
                                    !Get the fugacity of solid CO2
                                    gl%solidpos_akt_phase = gl%solid_pos
                                    pos_dryIce = gl%solidpos_akt_phase
                                    fug_CO2(1) = fug_DryIce(gl,Temp,press, pos_dryIce) * 1.D6
                                    if (fug_CO2(1) > 1.D-12) then  !No error occured
                                        !Get hydrate phase properties
                                        !Get the chemical potential of water in hydrate
                                        call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_CO2,ChemPot_hyd)
                                        !ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                                        !Get the composition of the hydrate phase
                                        call hdrt_mole_fract(gl,Temp,press*1.D6,fug_CO2, occup, CiJ, xH,occup_single, occup_double)
                                        !x_hyd(1) = xH(1) !Molfractions of water in hydrate
                                        !x_hyd(2) = xH(2) !Molfractions of gas in hydrate
                                        x_hyd = xH

                                        phasetype(1) = 4
                                        phasetype(2) = 5
                                        nrofphases = 2
                                        !Solid CO2 properties
                                        rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                                        x_phase(:,phasetype(1)) = 0.D0
                                        x_phase(pos_dryIce,phasetype(1)) = 1.D0
                                        !Hydrate properties
                                        call hdrt_density(gl,temp,press,fug_CO2,occup,CiJ,rho_H)
                                        x_phase(:,phasetype(2)) = x_hyd
                                        rho(phasetype(2)) = rho_H

                                        phasefrac(phasetype(1)) = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                                        phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))

                                        !Save information of the equilibrium found
                                        gl%solidtype_akt_phase = 2             !Solid is co2
                                        gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector

                                        errval = 0

                                        return
                                    else

                                        errval = -7779
                                        gl%molfractions = x
                                        call reduced_parameters_calc(gl,Temp)
                                        return

                                    end if

                                end if

                            end if
                        end if
                        !!Check if hydrates form
                        !molfractions = x_phase(:,phasetype(1))
                        !call reduced_parameters_calc(gl,Temp)
                        !
                        !call lnf_mix(Temp, rho(phasetype(1)), press, lnf)
                        !!The fugacity of CO2 in the vapor phase will be used as input for the chem. pot. of water in hydrate (CO2 on position 2 in Hydrate_list)
                        !fug_CO2 = dexp(lnf(2))*1.d6
                        !!fug_gas(1:nrofhydrateformers-1) = dexp(lnf(2:nrofhydrateformers))*1.d6     !Multi component hydrates
                        !!Get the chemical potential of water in hydrate
                        !call hdrt_chem_potent_w(Temp,press*1.d6,fug_CO2,ChemPot_hyd)
                        !!call hdrt_chem_potent_w(T,p*1.d6,fug_gas,chpw)                             !Multi component hydrates
                        !!ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                        !!Get the chemical potential of water in the fluid phase
                        !call Chempot_CALC(Temp, rho(phasetype(1)), Chempot, 0)
                        !tpd_sol_phase2 = ChemPot_hyd - Chempot(1)
                        !
                        !!Get the composition of the hydrate phase
                        !call hdrt_mole_fract(Temp,press*1.D6,fug_CO2, occup, CiJ, xH, occup_ls, occup_ld, occup_sms, occup_smd)
                        !x_hyd(1) = xH(1) !Molfractions of water in hydrate
                        !x_hyd(2) = xH(2) !Molfractions of gas in hydrate
                        !
                        !!If the tangent plane distance is "close" to 0, try to calculate a three phase equilibrium (VHIW)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!THREE PHASE CHECK VHIW or LHIw
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !if ((nrofphases == 2) .and. (dabs(tpd_sol_phase2) < 10.D0) .and. (ncomp > 2)) then
                        !
                        !    x_vap = x_phase(:,phasetype(1))
                        !    x_solid = 0.D0
                        !    x_solid(1) = 1.D0
                        !
                        !    !Setup for the three phase flash
                        !    solidtype(1) = 1 !Solid water
                        !    solidtype(2) = 1 !Hydrate forms
                        !    solid_pos = 1
                        !
                        !    iPhase_sol = 2 !Liqhter phase is vapor
                        !    iFlash = 5  !T and p given
                        !
                        !    rhovap_est = 0.d0
                        !    rholiq2_est = 0.D0
                        !
                        !    !Arbitrary startvalue for phasefractions
                        !    phasefrac_3Phase_solid = 1.D0 / 3.D0
                        !
                        !    call ptflash_solid_NC_3P(press, Temp, x, rho_sol, x_vap, x_liq2, x_solid, x_hyd, rhovap_est, &
                        !            & rholiq2_est, phasefrac_3Phase_solid, iFlash, iPhase_sol, iter, errval)
                        !
                        !    if ((ncomp == 3) .and. (errval == 0)) then
                        !
                        !        Mat_phasefrac_3Comp(1,1) = x_vap(1)
                        !        Mat_phasefrac_3Comp(2,1) = x_solid(1)
                        !        Mat_phasefrac_3Comp(3,1) = x_hyd(1)
                        !
                        !        Mat_phasefrac_3Comp(1,2) = x_vap(2)
                        !        Mat_phasefrac_3Comp(2,2) = x_solid(2)
                        !        Mat_phasefrac_3Comp(3,2) = x_hyd(2)
                        !
                        !        Mat_phasefrac_3Comp(1,3) = 1.D0
                        !        Mat_phasefrac_3Comp(2,3) = 1.D0
                        !        Mat_phasefrac_3Comp(3,3) = 1.D0
                        !
                        !        Vec_rightside_3Comp(1) = x(1)
                        !        Vec_rightside_3Comp(2) = x(2)
                        !        Vec_rightside_3Comp(3) = 1.D0
                        !
                        !        Vec_phasefrac_3Comp = 0.D0
                        !
                        !        !Solve the system of equations
                        !        eqn = 3
                        !        call LUdecomp(eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                        !        phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)
                        !
                        !    end if
                        !
                        !    if (minval(phasefrac_3Phase_solid) < -1.D-16) then
                        !        errval = -3333
                        !    end if
                        !    if (errval == 0) then
                        !        phasetype(1) = 1
                        !        phasetype(2) = 4
                        !        phasetype(3) = 5
                        !        x_Phase(:,phasetype(1)) = x_vap
                        !        x_Phase(:,phasetype(2)) = x_solid
                        !        x_Phase(:,phasetype(3)) = x_hyd
                        !        Phasefrac(phasetype(1)) = phasefrac_3Phase_solid(1)
                        !        Phasefrac(phasetype(2)) = phasefrac_3Phase_solid(2)
                        !        Phasefrac(phasetype(3)) = phasefrac_3Phase_solid(3)
                        !        rho(phasetype(1)) = rho_sol(1)
                        !        rho(phasetype(2)) = rho_sol(2)
                        !        rho(phasetype(3)) = rho_sol(3)
                        !        nrofphases = 3
                        !
                        !        !Save information of the equilibrium found
                        !        solidtype_akt_phase = 1             !Solid is water
                        !        solidpos_akt_phase = solid_pos      !Water on position 1
                        !
                        !        return
                        !    else
                        !        errval = 0
                        !    end if
                        !
                        !end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!END OF THREE PHASE CHECK VHIw or LHIw
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !
                        !if (tpd_sol_phase2 < 0.D0) then !Hydrates form, VH, LH, IwH, IcH
                        !
                        !    !An equilbrium of a fluid phase (vapor) with solid water was calculated. Try to calculate a vapor / hydrate equilbrium and check if
                        !    !this equilibrium is physically reasonable. If not calculate hydrate / Solid Water
                        !    iFlash = 3
                        !    solidtype(1) = 0   !No solid forms
                        !    solidtype(2) = 1   !Hydrate assumed
                        !    x_fluid1 = x_phase(:,phasetype(1))
                        !    !VIw or LxIw equilibrium was calculated. Check whether the first phase is vapor or liquid
                        !    If (phasetype(1) == 1) then
                        !        iphase = 2  !Vapor assumed
                        !    elseif ((phasetype(1) == 2) .or. (phasetype(1) == 3)) then
                        !        iphase = 1  !Liquid assumed
                        !    else
                        !        iphase = 2  !Backup: Vapor assumed
                        !    end if
                        !    call ptflash_solid_NC_2P(press, Temp, rho_sol, x, x_solid, x_hyd1, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval1, iter)
                        !
                        !    if ((errval1 == 0) .and. (vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then
                        !        phasetype(1) = 1
                        !        rho(phasetype(1)) = rho_sol(1)!rho1
                        !        x_phase(:,phasetype(1)) = x_fluid1
                        !        !vapfrac = vapfrac1 - 10
                        !        phasetype(2) = 5
                        !        rho(phasetype(2)) = rho_sol(2)
                        !        x_phase(:,phasetype(2)) = x_hyd1
                        !        nrofphases = 2
                        !        errval = 0
                        !        Phasefrac(phasetype(1)) = (x(1)-x_hyd1(1)) / (x_fluid1(1)-x_hyd1(1))
                        !        Phasefrac(phasetype(2)) = 1.D0 - Phasefrac(phasetype(1))
                        !        nrfluidphases = 1
                        !        !Save information of the equilibrium found
                        !        solidtype_akt_phase = 0             !No pure solids
                        !        solidpos_akt_phase = 0              !No pure solids
                        !
                        !        !Vapor / hydrate equilibrium found. In case of a binary mixture check if solid CO2 forms
                        !        if ((ncomp == 2) .and. (Fluidlist_hydrate(2) == "co2")) then
                        !            !ADD SOMETHING HERE!!!
                        !            !**
                        !            !--
                        !            !**
                        !            !--
                        !        end if
                        !
                        !
                        !    else    !Try solid water + hydrate
                        !        iFlash = 3
                        !        x_solid(solid_pos) = 1.d0
                        !        solidtype(1) = 1   !solid forms
                        !        solidtype(2) = 1   !hydrate assumed
                        !        vapfrac1 = 0.D0
                        !        iphase = 0
                        !        call ptflash_solid_NC_2P(press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid, rhofluid_est, vapfrac1, iFlash, iPhase, errval1, iter)
                        !        ! Check if the result is reasonable
                        !        ! If yes, assume hydrate/ice equilibrium and check if the result is reasonable
                        !
                        !        if ((errval1 == 0) .and. (vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then
                        !            phasetype(1) = 4
                        !            rho(phasetype(1)) = rho_sol(1)!1.D0 / v_WaterIce(gl,Temp, press)
                        !            x_phase(:,phasetype(1)) = x_solid
                        !            phasetype(2) = 5
                        !            rho(phasetype(2)) = rho_sol(2)
                        !            x_phase(:,phasetype(2)) = x_hyd
                        !            nrofphases = 2
                        !            errval = 0
                        !            !vapfrac = vapfrac1 - 10
                        !            Phasefrac(phasetype(1)) = (x(1)-x_hyd(1)) / (x_solid(1)-x_hyd(1))
                        !            Phasefrac(phasetype(2)) = 1.D0 - Phasefrac(phasetype(1))
                        !            !Andreas, March 2014
                        !            !Save information of the equilibrium found
                        !            solidtype_akt_phase = 1             !Solid is water
                        !            solidpos_akt_phase = solid_pos      !Water on position 1
                        !            return
                        !        else
                        !            errval = errval_fluid
                        !            !Andreas, March 2014
                        !            return
                        !        end if
                        !    end if
                        !
                        !    molfractions = x
                        !    call reduced_parameters_calc(gl,Temp)
                        !end if

                        !If the calculation of vapor / solid water equilibrium failed, try hydrate / solid water equilibrium
                    else
                        iFlash = 3
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.d0
                        gl%solidtype(1) = 1   !solid water forms
                        gl%solidtype(2) = 1   !hydrate assumed
                        iphase = 0
                        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid, rhofluid_est, vapfrac1, iFlash, iPhase, errval1, iter)
                        ! Check if the result is reasonable
                        ! If yes, assume hydrate/ice equilibrium and check if the result is reasonable

                        !If no error occured and a physically reasonable vapor fraction was found, assume hydrate / solid water to be stable
                        if ((errval1 == 0) .and. (vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then
                            phasetype(1) = 4
                            rho(phasetype(1)) = rho_sol(1)!1.D0 / v_WaterIce(gl,Temp, press)
                            x_phase(:,phasetype(1)) = x_solid
                            phasetype(2) = 5
                            rho(phasetype(2)) = rho_sol(2)
                            x_phase(:,phasetype(2)) = x_hyd
                            nrofphases = 2
                            errval = 0
                            !vapfrac = vapfrac1 - 10
                            !Andreas March 2014
                            Phasefrac(phasetype(1)) = (x(1)-x_hyd(1)) / (x_solid(1)-x_hyd(1))
                            Phasefrac(phasetype(2)) = 1.D0 - Phasefrac(phasetype(1))
                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 1             !Solid is water
                            gl%solidpos_akt_phase = gl%solid_pos      !Water on position 1
                            return
                        else
                            !Last try.
                            !If water and co2 are present in a binary mixture, check whether dry ice forms. Andreas May 2014
                            if ((gl%ncomp == 2) .and. (gl%Fluidlist_hydrate(2) == "co2")) then

                                gl%solid_pos = 2
                                gl%solidpos_akt_phase = 2
                                ! Get the melting or sublimation pressure at the given temperature
                                if (Temp > gl%ttp(gl%solid_pos)) then
                                    peq = pmelt_eq(gl,Temp, gl%solid_pos)
                                else
                                    peq = psub_eq(gl,Temp, gl%solid_pos)
                                end if
                                !if (dexp(lnf(solid_pos)) > (0.9D0 *peq)) then
                                if (press > (0.9D0 *peq)) then

                                    !Get the fugacity of solid CO2
                                    pos_dryIce = gl%solidpos_akt_phase
                                    gl%solid_pos = gl%solidpos_akt_phase
                                    fug_CO2(1) = fug_DryIce(gl,Temp,press, pos_dryIce) * 1.D6
                                    if (fug_CO2(1) > 1.D-12) then  !No error occured
                                        !Get hydrate phase properties
                                        !Get the chemical potential of water in hydrate
                                        call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_CO2,ChemPot_hyd)
                                        !ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                                        !Get the composition of the hydrate phase
                                        call hdrt_mole_fract(gl,Temp,press*1.D6,fug_CO2, occup, CiJ, xH,occup_single, occup_double)
                                        !x_hyd(1) = xH(1) !Molfractions of water in hydrate
                                        !x_hyd(2) = xH(2) !Molfractions of gas in hydrate
                                        x_hyd = xH

                                        phasetype(1) = 4
                                        phasetype(2) = 5
                                        nrofphases = 2
                                        !Solid CO2 properties
                                        rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                                        x_phase(:,phasetype(1)) = 0.D0
                                        x_phase(pos_dryIce,phasetype(1)) = 1.D0
                                        !Hydrate properties
                                        call hdrt_density(gl,temp,press,fug_CO2,occup,CiJ,rho_H)
                                        x_phase(:,phasetype(2)) = x_hyd
                                        rho(phasetype(2)) = rho_H

                                        phasefrac(phasetype(1)) = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                                        phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))

                                        !Save information of the equilibrium found
                                        gl%solidtype_akt_phase = 2             !Solid is co2
                                        gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector

                                        errval = 0

                                        return
                                    else
                                        errval = -7779
                                        gl%molfractions = x
                                        call reduced_parameters_calc(gl,Temp)
                                        return
                                    end if

                                else
                                    errval = errval_fluid
                                    !Andreas May 2014
                                    gl%molfractions = x
                                    call reduced_parameters_calc(gl,Temp)
                                    return
                                end if

                            else
                                errval = errval_fluid
                                !Andreas March 2014
                                gl%molfractions = x
                                call reduced_parameters_calc(gl,Temp)
                                return
                            end if

                        end if

                        gl%molfractions = x
                        call reduced_parameters_calc(gl,Temp)
                        return

                    end if


                end if
            else
                errval = -15566
                return
            end if


            !------------------------------------------------------------------------------

        end if

        !SH 09/18
        !Ensure that all phasetypes have a designated value for all phases that have been found until now
        do i = 1,nrofphases
            if (phasetype(i) == 0) then
                errval = -15565
                return
            end if
        enddo

        !Andreas, March 2014
        !Up to this point the following cases may have occured
        !1) Single fluid phase in homogeneous state
        !2) two phase equilibrium of fluid phases
        !2a) VLE
        !2a) LLE
        !3) three phase equilibrium of fluid phases
        !3a) VLL
        !3b) VHIw (At the moment the algorithm returns if HIw was found)
        !4) two phase equilibrium with solid phases
        !4a) VIw / (LxIw?)
        !4b) VH  / (LxH?)
        !4c) HIw (At the moment the algorithm returns if HIw was found)

        !Check whether solid phases might occur
        !Andreas Jan 2012


        !The tangent plane distance criterion TPD(gl,x_trial) = (sum i to N)(xtrial(i)(chempot(x_trial) - chempot(x))
        !for a PURE solid phase reduces to TPD(gl,x_trial) = gsol(x_trial) - chempoti(x)
        !If TPD >= 0 no solid phase occurs, if tpd < 0, the phases found are unstable and solid phases will form

        !For hydrates the TPD criterion is also reduced to TPD(gl,x_trial) = xtrial(water)*(chempot_hyd(water) - chempot_fluid(water))
        !If TPD >= 0 hydrate forms, if tpd < 0, no hydrate forms
        !THIS CRITERION HAS TO BE TESTED!!!

        !ASSUMPTION, HAS TO BE CHECKED: If 2 fluid phases are present in a binary mixture and a solid forms, the phase with the composition
        !"closer" to the solid composition vanishes.

        !Strategy:
        !1) Check if CO2 is in the mixture and solid CO2 forms
        !1a) If yes, calculate the new equilibrium
        !1b) If no, continue with 2)
        !2) Check if H2O is in the mixture and solid H2O forms
        !2a) If yes, calculate the new equilibrium
        !2b) If no, continue with 3)
        !3) Check if Hydrates form
        !3) If yes, calculate the new equilibrium

        peq = 0.D0
        pmelt = 0.D0
        psub = 0.D0
        tpd_sol = 0.D0
        rhofluid_est = 0.D0
        gl%solid_pos = 0
        !This variable is needed, if hydrate should form according to the tpd but no valid equilibrium with hydrate was found, before hydrate solid H2O is checked.
        !in this case, Hydrate solid must be calculated without checking the tpd
        check_Hydsolwater = .false.
        !Save the original fluid phases
        x_phase1 = x_phase(:,phasetype(1))
        rho_phase1 =  rho(phasetype(1))
        if (nrofphases == 2) then
            x_phase2 = x_phase(:,phasetype(2))
            rho_phase2 =  rho(phasetype(2))
        end if

        !---------------------------------------------------------------------------------------------------------
        if (.not. error_then_ice_found) then    !Only check for solid CO2 and solid H2O if not already an equilibrium with ice was calculated
            !Step 1)
            !Check if CO2 is in the mixture and thus dry ice may form
            Do i = 1, gl%ncomp!nrofsolids
                if (gl%Fluidlist_hydrate(i) == "co2") then
                    gl%solid_pos = i
                    exit
                end if
            end do
            if (gl%solid_pos /= 0) then
                !Check if the given conditions are close to the formation conditions of dry ice
                ! This is true, if p > psub(Temp) or p > pmelt(Temp)
                !Get the fugacity of CO2 in the mixture
                !molfractions = x_phase(:,phasetype(1))
                !call reduced_parameters_calc(gl,Temp)
                !call lnf_mix(Temp, rho(phasetype(1)), press, lnf)
                !molfractions = x
                !call reduced_parameters_calc(gl,Temp)

                ! Get the melting or sublimation pressure at the given temperature
                if (Temp > gl%ttp(gl%solid_pos)) then
                    peq = pmelt_eq(gl,Temp, gl%solid_pos)
                else
                    peq = psub_eq(gl,Temp, gl%solid_pos)
                end if
                !if (dexp(lnf(solid_pos)) > (0.9D0 *peq)) then
                if (press > (0.9D0 *peq)) then
                    gl%molfractions = x_phase(:,phasetype(1))
                    call reduced_parameters_calc(gl,Temp)
                    call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
                    gl%molfractions = x
                    call reduced_parameters_calc(gl,Temp)

                    tpd_sol = g_DryIce(gl,Temp,press) - Chempot(gl%solid_pos)

                    !If the tpd is "close" to 0, try to calculate a three phase equilibrium
                    if ((nrofphases == 2) .and. (nrfluidphases == 2) .and. (dabs(tpd_sol) < 10.D0) .and. (gl%ncomp > 2)) then

                        x_vap = x_phase(:,phasetype(1))
                        x_liq2 = x_phase(:,phasetype(2))
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.D0

                        !Setup for the three phase flash
                        gl%solidtype(1) = 2 !Solid co2
                        gl%solidtype(2) = 0 !No hydrate forms
                        !solid_pos = 2

                        if (phasetype(1) == 1) then
                            iPhase_sol = 2 !Liqhter phase is vapor
                        else
                            iPhase_sol = 1 !Liqhter phase is liquid
                        end if
                        iFlash = 7  !T and p given

                        rhovap_est = 0.d0
                        rholiq2_est = 0.D0

                        !Startvalue for phasefractions
                        !If 3 components are in a three phase equilibrium, the startvalues for the phasefractions can be solved analytically
                        if (gl%ncomp == 3) then
                            !Mat_phasefrac_3Comp(1,1) = x_vap(1)
                            !Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                            !Mat_phasefrac_3Comp(3,1) = x_solid(1)
                            !
                            !Mat_phasefrac_3Comp(1,2) = x_vap(2)
                            !Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                            !Mat_phasefrac_3Comp(3,2) = x_solid(2)
                            !
                            !Mat_phasefrac_3Comp(1,3) = 1.D0
                            !Mat_phasefrac_3Comp(2,3) = 1.D0
                            !Mat_phasefrac_3Comp(3,3) = 1.D0
                            !
                            !Vec_rightside_3Comp(1) = x(1)
                            !Vec_rightside_3Comp(2) = x(2)
                            !Vec_rightside_3Comp(3) = 1.D0
                            !
                            !Vec_phasefrac_3Comp = 0.D0
                            !
                            !!Solve the system of equations
                            !eqn = 3
                            !call LUdecomp(eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                            !if (errval /= 0) then
                            !    phasefrac_3Phase_solid = 1.D0 / 3.D0
                            !else
                            !    phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)
                            !end if
                        else
                            !Arbitrary start values
                            phasefrac_3Phase_solid = 1.D0 / 3.D0
                        end if

                        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_vap, x_liq2, x_solid, x_hyd, rhovap_est, &
                            & rholiq2_est, phasefrac_3Phase_solid, iFlash, iPhase_sol, iter, errval)

                        if ((gl%ncomp == 3) .and. (errval == 0)) then

                            Mat_phasefrac_3Comp(1,1) = x_vap(1)
                            Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                            Mat_phasefrac_3Comp(3,1) = x_solid(1)

                            Mat_phasefrac_3Comp(1,2) = x_vap(2)
                            Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                            Mat_phasefrac_3Comp(3,2) = x_solid(2)

                            Mat_phasefrac_3Comp(1,3) = 1.D0
                            Mat_phasefrac_3Comp(2,3) = 1.D0
                            Mat_phasefrac_3Comp(3,3) = 1.D0

                            Vec_rightside_3Comp(1) = x(1)
                            Vec_rightside_3Comp(2) = x(2)
                            Vec_rightside_3Comp(3) = 1.D0

                            Vec_phasefrac_3Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                            phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)

                        end if

                        if (minval(phasefrac_3Phase_solid) < -1.D-16) then
                            errval = -3333
                        end if
                        if (errval == 0) then
                            !phasetype(1) = 1
                            !phasetype(2) = 3
                            phasetype(3) = 4
                            !x_Phase(:,phasetype(1)) = x_vap
                            !x_Phase(:,phasetype(2)) = x_solid
                            x_Phase(:,phasetype(3)) = x_solid
                            Phasefrac(phasetype(1)) = phasefrac_3Phase_solid(1)
                            Phasefrac(phasetype(2)) = phasefrac_3Phase_solid(2)
                            Phasefrac(phasetype(3)) = phasefrac_3Phase_solid(3)
                            rho(phasetype(1)) = rho_sol(1)
                            rho(phasetype(2)) = rho_sol(2)
                            rho(phasetype(3)) = rho_sol(3)
                            nrofphases = 3

                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 2             !Solid is co2
                            gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector
                            return
                        else
                            errval = 0
                        end if

                    end if

                    if (tpd_sol < 0) then       !Dry Ice forms
                        iFlash = 3

                        !Step 1a
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.d0
                        gl%solidtype(1) = 2   !solid forms (solid CO2)
                        gl%solidtype(2) = 0   !No hydrate assumed
                        tpd_sol_phase1 = 1.D6
                        tpd_sol_phase2 = 1.D6
                        vapfrac1 = 0
                        vapfrac2 = 0
                        errval1 = 0
                        errval2 = 0

                        ! In case of two fluid phases the correct equilibrium needs to be found
                        ! Thus both fluid phases have to be tested with the solid phase
                        ! Strategy:
                        ! 1.1) Calculate the equilibrium for both phases
                        ! 1.2) check, whether the molar vaporfraction is fine in both cases, if not, take the results with the physically correct vapor fraction
                        ! 1.3) if both vapor fractions are correct, check the tangent plane distance at the overallcomposition. The equilibrium with the lower tangent plane distance is the correct one, since the energy is lower

                        !Step 1.1)
                        x_fluid1 = x_phase(:,phasetype(1))
                        iphase = 0
                        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval, iter)
                        if(errval == 0) then
                            !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                            if ((vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then  !Physically plausible solution
                                !Check if this solution is stable or a "fake solution"
                                !Calculate the tpd for this solution at the other fluid phase as trial composion
                                !if the resulting tpd is negative, a "fake solution" was found
                                rho1 = gl%rho_vap
                                phys_reas1 = .true.
                                if (nrfluidphases == 2) then
                                    tpd_phase2 = tpd(gl,press, temp, x_fluid1, rho1, x_phase2, errval1)
                                else
                                    tpd_phase2 = 1 !dummy value
                                    errval1 = 0
                                End if
                                if ((tpd_phase2 > 0) .and. (errval1 == 0)) then
                                    !calculate the tpd at the overall composition for this solution
                                    tpd_sol_phase1 =  tpd(gl,press, temp, x_fluid1, rho1, x, errval1)
                                    if (errval1 /= 0) tpd_sol_phase1 = 1.d6
                                else
                                    tpd_sol_phase1 = 1.d6
                                    errval1 = 0
                                end if
                            else
                                phys_reas1 = .false.
                            end if
                        end if
                        gl%molfractions = x
                        call reduced_parameters_calc(gl,Temp)
                        !Andreas February 2015. "Errval == 0" statement in if case deleted. Not needed and makes problems
                        !if ((nrfluidphases == 2) .and. (errval == 0)) then
                        if (nrfluidphases == 2) then
                            x_fluid2 = x_phase(:,phasetype(2))
                            iphase = 0
                            call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid2, rhofluid_est, vapfrac2, iFlash, iPhase, errval, iter)
                            if(errval == 0) then
                                !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                                if ((vapfrac2 > 0) .and. (vapfrac2 < 1)) then  !Physically plausible solution
                                    !Check if this solution is stable or a "fake solution"
                                    !Calculate the tpd for this solution at the other fluid phase as trial composion
                                    !if the resulting tpd is negative, a "fake solution" was found
                                    rho2 = gl%rho_vap
                                    tpd_phase1 = tpd(gl,press, temp, x_fluid2, rho2, x_phase1, errval2)
                                    phys_reas2 = .true.
                                    if ((tpd_phase1 > 0) .and. (errval2 == 0)) then
                                        !calculate the tpd at the overall composition for this solution
                                        tpd_sol_phase2 =  tpd(gl,press, temp, x_fluid2, rho2, x, errval2)
                                        if (errval2 /= 0) tpd_sol_phase2 = 1.d6
                                    else
                                        tpd_sol_phase2 = 1.d6
                                        errval2 = 0
                                    end if
                                else
                                    phys_reas2 = .false.
                                end if
                            end if
                            gl%molfractions = x
                            call reduced_parameters_calc(gl,Temp)
                        end if
                        if ((errval1 == 0) .or. (errval2 == 0)) then
                            if ((tpd_sol_phase1 < 1.d6) .or. (tpd_sol_phase2 < 1.d6)) then
                                if (tpd_sol_phase1 < tpd_sol_phase2) then
                                    rho(phasetype(1)) = rho1
                                    x_phase(:,phasetype(1)) = x_fluid1
                                    vapfrac = vapfrac1
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_DryIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac1
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 2             !Solid is co2
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector
                                else if (tpd_sol_phase1 > tpd_sol_phase2) then
                                    rho(phasetype(1)) = rho2
                                    x_phase(:,phasetype(1)) = x_fluid2
                                    vapfrac = vapfrac2
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_DryIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac2
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 2             !Solid is co2
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector
                                end if
                            else
                                !No stable equilbrium found, quit with error
                                if (errval_fluid /= 0) then
                                    errval = errval_fluid
                                else
                                    errval = -2233
                                    !Andreas March 2014
                                    !The tpd indicates that dry ice forms, but no fluid phase in equilibrium can be found.
                                    !Last chance: Dry ice with hydrate
                                    if ((gl%nrofhydrateformers > 0) .and. (gl%nrofhydrateformers < 3)) then
                                        !Get the fugacity of solid CO2
                                        pos_dryIce = gl%solid_pos
                                        fug_CO2(1) = fug_DryIce(gl,Temp,press, pos_dryIce) * 1.D6
                                        if (fug_CO2(1) > 1.D-12) then  !No error occured
                                            !Get hydrate phase properties
                                            !Get the chemical potential of water in hydrate
                                            call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_CO2,ChemPot_hyd)
                                            !ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                                            !Get the composition of the hydrate phase
                                            call hdrt_mole_fract(gl,Temp,press*1.D6,fug_CO2, occup, CiJ, xH,occup_single, occup_double)
                                            !x_hyd(1) = xH(1) !Molfractions of water in hydrate
                                            !x_hyd(2) = xH(2) !Molfractions of gas in hydrate
                                            x_hyd = xH


                                            phasetype(1) = 4
                                            phasetype(2) = 5
                                            nrofphases = 2
                                            !Solid CO2 properties
                                            rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                                            x_phase(:,phasetype(1)) = 0.D0
                                            x_phase(pos_dryIce,phasetype(1)) = 1.D0
                                            !Hydrate properties
                                            call hdrt_density(gl,temp,press,fug_CO2,occup,CiJ,rho_H)
                                            x_phase(:,phasetype(2)) = x_hyd
                                            rho(phasetype(2)) = rho_H

                                            phasefrac(phasetype(1)) = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                                            phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))

                                            !Save information of the equilibrium found
                                            gl%solidtype_akt_phase = 2             !Solid is co2
                                            gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector

                                            errval = 0
                                        end if
                                    end if
                                end if

                                return
                            end if

                        end if
                    end if
                end if
            end if
            !---------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------
            !Step 2)
            !Check if H2O is in the mixture and thus water ice may form
            !If dry ice forms, this may not be checked because:
            !1) If the amount of water is higher than the solubility of the fluid phase, hydrate will form instead of solid water
            !2) If the amount of water is smaller than the solubility of the fluid phase, all water is dissolved in the fluid phase
            gl%solid_pos = 0
            Do i = 1, gl%ncomp
                if (gl%Fluidlist_hydrate(i) == "water") then
                    gl%solid_pos = i
                    exit
                end if
            end do
            if ((gl%solid_pos /= 0) .and.  (phasetype(2) < 4) .and. (gl%solidtype_akt_phase < 2)) then    !water is in the mixture and no solid co2 was found
                !Check if the given conditions are close to the formation conditions of solid H2O
                ! This is true, if p > psub(Temp) or p < pmelt(Temp)
                !Get the fugacity of H2O in the mixture
                !molfractions = x_phase(:,phasetype(1))
                !call reduced_parameters_calc(gl,Temp)
                !call lnf_mix(Temp, rho(phasetype(1)), press, lnf)
                !molfractions = x
                !call reduced_parameters_calc(gl,Temp)
                ! Get the melting AND sublimation pressure at the given temperature
                if (Temp < gl%ttp(gl%solid_pos)) then
                    pmelt = pmelt_eq(gl,Temp, gl%solid_pos)
                    if (pmelt < 1.D-14) pmelt = 1000.D0
                    psub = psub_eq(gl,Temp, gl%solid_pos)
                end if
                if ((press > (0.9D0 *psub)) .and. (press < (1.1D0 *pmelt))) then
                    gl%molfractions = x_phase(:,phasetype(1))
                    call reduced_parameters_calc(gl,Temp)
                    call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
                    gl%molfractions = x
                    call reduced_parameters_calc(gl,Temp)
                    tpd_sol = g_WaterIce(gl,Temp,press) - Chempot(gl%solid_pos)


                    !If the tpd is "close" to 0, try to calculate a three phase equilibrium
                    if ((nrofphases == 2) .and. (nrfluidphases == 2) .and. (dabs(tpd_sol) < 10.D0) .and. (gl%ncomp > 2)) then

                        x_vap = x_phase(:,phasetype(1))
                        x_liq2 = x_phase(:,phasetype(2))
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.D0

                        !Setup for the three phase flash
                        gl%solidtype(1) = 1 !Solid H2O
                        gl%solidtype(2) = 0 !No hydrate forms
                        !solid_pos = 1

                        if (phasetype(1) == 1) then
                            iPhase_sol = 2 !Liqhter phase is vapor
                        else
                            iPhase_sol = 1 !Liqhter phase is liquid
                        end if
                        iFlash = 7  !T and p given

                        rhovap_est = 0.d0
                        rholiq2_est = 0.D0

                        !Startvalue for phasefractions
                        !If 3 components are in a three phase equilibrium, the startvalues for the phasefractions can be solved analytically
                        if (gl%ncomp == 3) then

                            !Mat_phasefrac_3Comp(1,1) = x_vap(1)
                            !Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                            !Mat_phasefrac_3Comp(3,1) = x_solid(1)
                            !
                            !Mat_phasefrac_3Comp(1,2) = x_vap(2)
                            !Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                            !Mat_phasefrac_3Comp(3,2) = x_solid(2)
                            !
                            !Mat_phasefrac_3Comp(1,3) = 1.D0
                            !Mat_phasefrac_3Comp(2,3) = 1.D0
                            !Mat_phasefrac_3Comp(3,3) = 1.D0
                            !
                            !Vec_rightside_3Comp(1) = x(1)
                            !Vec_rightside_3Comp(2) = x(2)
                            !Vec_rightside_3Comp(3) = 1.D0
                            !
                            !Vec_phasefrac_3Comp = 0.D0
                            !
                            !!Solve the system of equations
                            !eqn = 3
                            !call LUdecomp(eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                            !if (errval /= 0) then
                            !    phasefrac_3Phase_solid = 1.D0 / 3.D0
                            !else
                            !    phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)
                            !end if
                        else
                            !Arbitrary start values
                            phasefrac_3Phase_solid = 1.D0 / 3.D0
                        end if

                        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_vap, x_liq2, x_solid, x_hyd, rhovap_est, &
                            & rholiq2_est, phasefrac_3Phase_solid, iFlash, iPhase_sol, iter, errval)

                        if ((gl%ncomp == 3) .and. (errval == 0)) then

                            Mat_phasefrac_3Comp(1,1) = x_vap(1)
                            Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                            Mat_phasefrac_3Comp(3,1) = x_solid(1)

                            Mat_phasefrac_3Comp(1,2) = x_vap(2)
                            Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                            Mat_phasefrac_3Comp(3,2) = x_solid(2)

                            Mat_phasefrac_3Comp(1,3) = 1.D0
                            Mat_phasefrac_3Comp(2,3) = 1.D0
                            Mat_phasefrac_3Comp(3,3) = 1.D0

                            Vec_rightside_3Comp(1) = x(1)
                            Vec_rightside_3Comp(2) = x(2)
                            Vec_rightside_3Comp(3) = 1.D0

                            Vec_phasefrac_3Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                            phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)

                        end if

                        if (minval(phasefrac_3Phase_solid) < -1.D-16) then
                            errval = -3333
                        end if

                        if (errval == 0) then
                            !phasetype(1) = 1
                            !phasetype(2) = 3
                            phasetype(3) = 4
                            !x_Phase(:,phasetype(1)) = x_vap
                            !x_Phase(:,phasetype(2)) = x_solid
                            x_Phase(:,phasetype(3)) = x_solid
                            Phasefrac(phasetype(1)) = phasefrac_3Phase_solid(1)
                            Phasefrac(phasetype(2)) = phasefrac_3Phase_solid(2)
                            Phasefrac(phasetype(3)) = phasefrac_3Phase_solid(3)
                            rho(phasetype(1)) = rho_sol(1)
                            rho(phasetype(2)) = rho_sol(2)
                            rho(phasetype(3)) = rho_sol(3)
                            nrofphases = 3

                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 1             !Solid is water
                            gl%solidpos_akt_phase = gl%solid_pos      !Position of water in the fluid vector

                            return
                        else
                            errval = 0
                        end if

                    end if

                    if (tpd_sol < 0) then
                        iFlash = 3

                        !Step 1a
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.d0
                        gl%solidtype(1) = 1   !solid forms (solid H2O)
                        gl%solidtype(2) = 0   !No hydrate assumed
                        tpd_sol_phase1 = 1.D6
                        tpd_sol_phase2 = 1.D6
                        vapfrac1 = 0
                        vapfrac2 = 0
                        errval1 = 0
                        errval2 = 0

                        ! In case of two fluid phases the correct equilibrium needs to be found
                        ! Thus both fluid phases have to be tested with the solid phase
                        ! Strategy:
                        ! 1.1) Calculate the equilibrium for both phases
                        ! 1.2) check, whether the molar vaporfraction is fine in both cases, if not, take the results with the physically correct vapor fraction
                        ! 1.3) if both vapor fractions are correct, check the tangent plane distance at the overallcomposition. The equilibrium with the lower tangent plane distance is the correct one, since the energy is lower

                        !Step 1.1)
                        x_fluid1 = x_phase(:,phasetype(1))
                        iphase = 0
                        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval, iter)
                        if(errval == 0) then
                            !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                            if ((vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then  !Physically plausible solution
                                !Check if this solution is stable or a "fake solution"
                                !Calculate the tpd for this solution at the other fluid phase as trial composion
                                !if the resulting tpd is negative, a "fake solution" was found
                                rho1 = gl%rho_vap
                                phys_reas1 = .true.
                                if (nrfluidphases == 2) then
                                    tpd_phase2 = tpd(gl,press, temp, x_fluid1, rho1, x_phase2, errval1)
                                else
                                    tpd_phase2 = 1 !dummy value
                                    errval1 = 0
                                End if
                                if ((tpd_phase2 > 0) .and. (errval1 == 0)) then
                                    !calculate the tpd at the overall composition for this solution
                                    tpd_sol_phase1 =  tpd(gl,press, temp, x_fluid1, rho1, x, errval1)
                                    if (errval1 /= 0) tpd_sol_phase1 = 1.d6
                                else
                                    tpd_sol_phase1 = 1.d6
                                    errval1 = 0
                                end if
                            else
                                phys_reas1 = .false.
                            end if
                        end if
                        gl%molfractions = x
                        call reduced_parameters_calc(gl,Temp)
                        !Andreas February 2015. "Errval == 0" statement in if case deleted. Not needed and makes problems
                        !if ((nrfluidphases == 2) .and. (errval == 0)) then
                        if (nrfluidphases == 2) then
                            x_fluid2 = x_phase(:,phasetype(2))
                            iphase = 0
                            call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid2, rhofluid_est, vapfrac2, iFlash, iPhase, errval, iter)
                            if(errval == 0) then
                                !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                                if ((vapfrac2 > 0) .and. (vapfrac2 < 1)) then  !Physically plausible solution
                                    !Check if this solution is stable or a "fake solution"
                                    !Calculate the tpd for this solution at the other fluid phase as trial composion
                                    !if the resulting tpd is negative, a "fake solution" was found
                                    rho2 = gl%rho_vap
                                    phys_reas2 = .true.
                                    tpd_phase1 = tpd(gl,press, temp, x_fluid2, rho2, x_phase1, errval2)
                                    if ((tpd_phase1 > 0) .and. (errval2 == 0)) then
                                        !calculate the tpd at the overall composition for this solution
                                        tpd_sol_phase2 =  tpd(gl,press, temp, x_fluid2, rho2, x, errval2)
                                        if (errval2 /= 0) tpd_sol_phase2 = 1.d6
                                    else
                                        tpd_sol_phase2 = 1.d6
                                        errval2 = 0
                                    end if
                                else
                                    phys_reas2 = .false.
                                end if
                            end if
                            gl%molfractions = x
                            call reduced_parameters_calc(gl,Temp)
                        end if
                        if ((errval1 == 0) .or. (errval2 == 0)) then
                            if ((tpd_sol_phase1 < 1.D6) .or. (tpd_sol_phase1 < 1.D6)) then
                                if (tpd_sol_phase1 < tpd_sol_phase2) then
                                    rho(phasetype(1)) = rho1
                                    x_phase(:,phasetype(1)) = x_fluid1
                                    !vapfrac = vapfrac1
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_WaterIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac1
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 1             !Solid is water
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of water in the fluid vector

                                else if (tpd_sol_phase1 > tpd_sol_phase2) then
                                    rho(phasetype(1)) = rho2
                                    x_phase(:,phasetype(1)) = x_fluid2
                                    vapfrac = vapfrac2
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_WaterIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac2
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 1             !Solid is water
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of water in the fluid vector

                                end if
                            else
                                !if a physically reasonable solution was found, try this solution, elsewise quit with error. Andreas May 2014
                                if (phys_reas1) then
                                    rho(phasetype(1)) = rho1
                                    x_phase(:,phasetype(1)) = x_fluid1
                                    !vapfrac = vapfrac1
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_WaterIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac1
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 1             !Solid is water
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of water in the fluid vector
                                    phys_reas1 = .false.

                                elseif (phys_reas2) then
                                    rho(phasetype(1)) = rho2
                                    x_phase(:,phasetype(1)) = x_fluid2
                                    vapfrac = vapfrac2
                                    phasetype(2) = 4
                                    rho(phasetype(2)) = 1.D0 / v_WaterIce(gl,Temp, press)
                                    x_phase(:,phasetype(2)) = x_solid
                                    nrofphases = 2
                                    phasefrac(phasetype(1)) = vapfrac2
                                    phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                    errval = 0
                                    nrfluidphases = 1
                                    !Save information of the equilibrium found
                                    gl%solidtype_akt_phase = 1             !Solid is water
                                    gl%solidpos_akt_phase = gl%solid_pos      !Position of water in the fluid vector
                                    phys_reas2 = .false.

                                else
                                    !No stable equilbrium found, quit with error
                                    if (errval_fluid /= 0) then
                                        errval = errval_fluid
                                    else
                                        errval = -2233
                                    end if
                                    return
                                end if
                            end if
                        end if
                    end if
                end if
            end if
            !---------------------------------------------------------------------------------------------------------
        end if


        !Step 3
        !If hydrates may form, check whether the temperature is in a range where this may occur
        gl%solid_pos = 0
        !Pure Hydrates
        !if ((nrofhydrateformers > 0) .and. (nrofhydrateformers < 3) .and. (Temp < 330.D0)) then
        !Mixed Hydrates
        if ((gl%nrofhydrateformers > 0) .and. (Temp < 330.D0)) then
            tpd_sol_hyd = 0.d0
            do hdrt_structure = 1,2
                check_Hydsolwater = .false. !see beginning of this routine
                call hdrt_structure_definition(gl,hdrt_structure)
                tpd_sol_phase1 = 1.D6
                tpd_sol_phase2 = 1.d6
                tpd_sol = 0.D0
                errval1 = 0
                errval2 = 0
                vapfrac1 = 0
                vapfrac2 = 0

                gl%molfractions = x_phase(:,phasetype(1))
                call reduced_parameters_calc(gl,Temp)

                call lnf_mix(gl,Temp, rho(phasetype(1)), press, lnf)
                !The fugacity of CO2 in the vapor phase will be used as input for the chem. pot. of water in hydrate (CO2 on position 2 in Hydrate_list)
                !fug_CO2 = dexp(lnf(2))*1.d6
                fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnf(2:gl%nrofhydrateformers))*1.d6     !Multi component hydrates
                !Get the chemical potential of water in hydrate
                call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_gas,ChemPot_hyd)
                !call hdrt_chem_potent_w(T,p*1.d6,fug_gas,chpw)                             !Multi component hydrates
                !ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                !Get the chemical potential of water in the fluid phase
                call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
                !tpd_sol = ChemPot_hyd - Chempot(1)
                tpd_sol_hyd(hdrt_structure) = ChemPot_hyd - Chempot(1)

                gl%molfractions = x
                call reduced_parameters_calc(gl,Temp)


                !If the tangent plane distance is "close" to 0, try to calculate a three phase equilibrium
                !---
                !if ((nrofphases == 2) .and. (dabs(tpd_sol) < 20.D0)  .and. (gl%ncomp > 2)) then
                if ((nrofphases == 2) .and. (dabs(tpd_sol_hyd(hdrt_structure)) < 20.D0)  .and. (gl%ncomp > 2)) then

                    if (phasetype(2) < 4) then          !VLwH, VLcH, LwLcH

                        x_liq2 = x_phase(:,phasetype(2))
                        !Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                        !Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                        gl%solidtype(1) = 0

                    elseif (phasetype(2) == 4) then     !VHIw, LwHIw, LcHIw, VHIc, LcHIc, LwHIc

                        x_solid = 0.D0
                        if (gl%solidtype_akt_phase == 1) then         !Solid H2O
                            x_solid(gl%solidpos_akt_phase) = 1.D0
                            gl%solidtype(1) = 1 !Solid water
                            gl%solid_pos = gl%solidpos_akt_phase
                        elseif (gl%solidtype_akt_phase == 2) then     !Solid CO2
                            x_solid(gl%solidpos_akt_phase) = 1.D0
                            gl%solidtype(1) = 2 !Solid CO2
                            gl%solid_pos = gl%solidpos_akt_phase
                        else
                            errval = -15566
                            return
                        end if
                        !Mat_phasefrac_3Comp(2,1) = x_solid(1)
                        !Mat_phasefrac_3Comp(2,2) = x_solid(2)
                    end if

                    if (phasetype(1) == 1) then
                        iPhase_sol = 2 !Liqhter phase is vapor
                    else
                        iPhase_sol = 1 !Liqhter phase is liquid
                    end if
                    iFlash = 7  !T and p given
                    gl%solidtype(2) = 1 !Hydrate forms

                    rhovap_est = 0.d0
                    rholiq2_est = 0.D0
                    x_vap = x_phase(:,phasetype(1))

                    !Arbitrary startvalue for phasefractions
                    phasefrac_3Phase_solid = 1.D0 / 3.D0

                    call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_vap, x_liq2, x_solid, x_hyd, rhovap_est, &
                        & rholiq2_est, phasefrac_3Phase_solid, iFlash, iPhase_sol, iter, errval)

                    if ((gl%ncomp == 3) .and. (errval == 0)) then

                        Mat_phasefrac_3Comp(1,1) = x_vap(1)
                        if (phasetype(2) < 4) then
                            Mat_phasefrac_3Comp(2,1) = x_liq2(1)
                        else
                            Mat_phasefrac_3Comp(2,1) = x_solid(1)
                        end if
                        Mat_phasefrac_3Comp(3,1) = x_hyd(1)

                        Mat_phasefrac_3Comp(1,2) = x_vap(2)
                        if (phasetype(2) < 4) then
                            Mat_phasefrac_3Comp(2,2) = x_liq2(2)
                        else
                            Mat_phasefrac_3Comp(2,2) = x_solid(2)
                        end if
                        Mat_phasefrac_3Comp(3,2) = x_hyd(2)

                        Mat_phasefrac_3Comp(1,3) = 1.D0
                        Mat_phasefrac_3Comp(2,3) = 1.D0
                        Mat_phasefrac_3Comp(3,3) = 1.D0

                        Vec_rightside_3Comp(1) = x(1)
                        Vec_rightside_3Comp(2) = x(2)
                        Vec_rightside_3Comp(3) = 1.D0

                        Vec_phasefrac_3Comp = 0.D0

                        !Solve the system of equations
                        eqn = 3
                        call LUdecomp(gl,eqn,Mat_phasefrac_3Comp,Vec_rightside_3Comp,errval,herr)
                        phasefrac_3Phase_solid = Vec_rightside_3Comp(1:3)

                    end if

                    if (minval(phasefrac_3Phase_solid) < -1.D-16) then
                        errval = -3333
                    end if
                    if (errval == 0) then
                        if (gl%solidtype(1) == 0) then
                            phasetype(2) = 3
                            x_Phase(:,phasetype(2)) = x_liq2
                        else
                            phasetype(2) = 4
                            x_Phase(:,phasetype(2)) = x_solid
                        end if
                        !phasetype(1) = 1
                        phasetype(3) = 5
                        x_Phase(:,phasetype(1)) = x_vap
                        x_Phase(:,phasetype(3)) = x_hyd
                        Phasefrac(phasetype(1)) = phasefrac_3Phase_solid(1)
                        Phasefrac(phasetype(2)) = phasefrac_3Phase_solid(2)
                        Phasefrac(phasetype(3)) = phasefrac_3Phase_solid(3)
                        rho(phasetype(1)) = rho_sol(1)
                        rho(phasetype(2)) = rho_sol(2)
                        rho(phasetype(3)) = rho_sol(3)
                        nrofphases = 3
                        gl%hdrt_structure_stable = hdrt_structure
                        if (hdrt_structure == 1) cycle
                        exit
                    else
                        errval = 0
                    end if

                end if

                !if (tpd_sol < 0) then
                if (tpd_sol_hyd(hdrt_structure) < 0) then
                    gl%hdrt_structure_stable = hdrt_structure
                    iFlash = 3
                    gl%solidtype(1) = 0   !No solid forms
                    gl%solidtype(2) = 1   !Hydrate assumed
                    !For hydrates the choice which phase forms together with the hydrate gets more complicated,
                    !since several cases may occur:
                    ! A fluid phase with solid CO2 in equilibrium was found
                    ! A fluid phase with solid H2O in equilibrium was found
                    ! Two fluid phases rich of the guest component(s) have been found
                    ! One guest component rich and one water rich phase found
                    ! Thus some criteria are needed, to distinguish, which phase is in equilibrium with the hydrate
                    ! Strategy:
                    !
                    ! 3.1) If the second phase is solid CO2, assume hydrate + solid CO2 in equilibrium (ONLY FOR BINARY MIXTURES!!! FOR ALL OTHER MIXTURES A THREE PHASE EQUILIBRIUM IS LIKELY)
                    ! 3.2) Calculate the equilibrium for both fluid phases (if exist)
                    ! 3.3) check, whether the molar vaporfraction is fine in both cases, if not, take the results with the physically correct vapor fraction
                    ! 3.4) if both vapor fractions are correct, check the tangent plane distance at the overallcomposition. The equilibrium with the lower tangent plane distance is the correct one, since the energy is lower
                    ! 3.5) Assume Hydrate/solid H2O if nothing else was found

                    !Step 3.1)
                    if ((gl%ncomp == 2) .and.  (phasetype(2) == 4) .and. (gl%solidtype_akt_phase == 2)) then      !Solid CO2 +  hydrate or error
                        !Get the fugacity of solid CO2
                        pos_dryIce = gl%solidpos_akt_phase
                        gl%solid_pos = gl%solidpos_akt_phase
                        fug_CO2(1) = fug_DryIce(gl,Temp,press, pos_dryIce) * 1.D6
                        if (fug_CO2(1) > 1.D-12) then  !No error occured
                            !Get hydrate phase properties
                            !Get the chemical potential of water in hydrate
                            call hdrt_chem_potent_w(gl,Temp,press*1.d6,fug_CO2,ChemPot_hyd)
                            !ChemPot_hyd              Chemical Potential of Water in Hydrate [J / mol]
                            !Get the composition of the hydrate phase
                            call hdrt_mole_fract(gl,Temp,press*1.D6,fug_CO2, occup, CiJ, xH,occup_single, occup_double)
                            !x_hyd(1) = xH(1) !Molfractions of water in hydrate
                            !x_hyd(2) = xH(2) !Molfractions of gas in hydrate
                            x_hyd = xH

                            phasetype(1) = 4
                            phasetype(2) = 5
                            nrofphases = 2
                            !Solid CO2 properties
                            rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                            x_phase(:,phasetype(1)) = 0.D0
                            x_phase(pos_dryIce,phasetype(1)) = 1.D0
                            !Hydrate properties
                            call hdrt_density(gl,temp,press,fug_CO2,occup,CiJ,rho_H)
                            x_phase(:,phasetype(2)) = x_hyd
                            rho(phasetype(2)) = rho_H

                            phasefrac(phasetype(1)) = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))

                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 2             !Solid is co2
                            gl%solidpos_akt_phase = gl%solid_pos      !Position of co2 in the fluid vector

                            errval = 0

                            return
                        end if
                    end if


                    !Step 3.2)
                    x_fluid1 = x_phase(:,phasetype(1))
                    if (phasetype(1) == 1) then
                        iphase = 2
                    elseif (phasetype(1) == 2) then
                        iphase = 1
                    else
                        iphase = 0
                    end if
                    call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd1, x_fluid1, rhofluid_est, vapfrac1, iFlash, iPhase, errval, iter)
                    if(errval == 0) then
                        !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                        if ((vapfrac1 > 0.D0) .and. (vapfrac1 < 1.D0)) then  !Physically plausible solution
                            !Check if this solution is stable or a "fake solution"
                            !Calculate the tpd for this solution at the other fluid phase as trial composion
                            !if the resulting tpd is negative, a "fake solution" was found
                            rho1 = rho_sol(1)!rho_vap
                            rho1_H = rho_sol(2)
                            if (nrfluidphases == 2) then
                                tpd_phase2 = tpd(gl,press, temp, x_fluid1, rho1, x_phase2, errval1)
                            else
                                tpd_phase2 = 1 !dummy value
                                errval1 = 0
                            End if
                            if ((tpd_phase2 > 0.D0) .and. (errval1 == 0)) then
                                !calculate the tpd at the overall composition for this solution
                                tpd_sol_phase1 =  tpd(gl,press, temp, x_fluid1, rho1, x, errval1)
                                if (errval1 /= 0) tpd_sol_phase1 = 1.d6
                            else
                                tpd_sol_phase1 = 1.d6
                                errval1 = 0
                            end if
                        end if
                    end if
                    gl%molfractions = x
                    call reduced_parameters_calc(gl,Temp)
                    if ((nrfluidphases == 2) .and. (errval == 0)) then
                        x_fluid2 = x_phase(:,phasetype(2))
                        iphase = 0
                        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd2, x_fluid2, rhofluid_est, vapfrac2, iFlash, iPhase, errval, iter)
                        if(errval == 0) then
                            !Check if the calculated equilibrium is plausible and calculate the tpd at the overall composition
                            if ((vapfrac2 > 0) .and. (vapfrac2 < 1)) then  !Physically plausible solution
                                !Check if this solution is stable or a "fake solution"
                                !Calculate the tpd for this solution at the other fluid phase as trial composion
                                !if the resulting tpd is negative, a "fake solution" was found
                                rho2 = rho_sol(1)!rho_vap
                                rho2_H = rho_sol(2)
                                tpd_phase1 = tpd(gl,press, temp, x_fluid2, rho2, x_phase1, errval2)
                                if ((tpd_phase1 > 0.D0) .and. (errval2 == 0)) then
                                    !calculate the tpd at the overall composition for this solution
                                    tpd_sol_phase2 =  tpd(gl,press, temp, x_fluid2, rho2, x, errval2)
                                    if (errval2 /= 0) tpd_sol_phase2 = 1.d6
                                else
                                    tpd_sol_phase2 = 1.d6
                                    errval2 = 0
                                end if
                            end if
                        end if
                        gl%molfractions = x
                        call reduced_parameters_calc(gl,Temp)
                    end if
                    if ((errval1 == 0) .or. (errval2 == 0)) then
                        if ((tpd_sol_phase1 < (1.d6 - 1.D-8)) .or. (tpd_sol_phase2 < (1.d6 - 1.D-8))) then
                            if (tpd_sol_phase1 < tpd_sol_phase2) then
                                phasetype(2) = 5
                                rho(phasetype(1)) = rho1
                                rho(phasetype(2)) = rho1_H
                                x_phase(:,phasetype(1)) = x_fluid1
                                x_phase(:,phasetype(2)) = x_hyd1
                                !vapfrac = vapfrac1
                                phasefrac(phasetype(1)) = vapfrac1
                                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                errval = 0
                                !Save information of the equilibrium found
                                gl%solidtype_akt_phase = 0         !No solid
                                gl%solidpos_akt_phase = 0          !No solid
                            else if (tpd_sol_phase1 > tpd_sol_phase2) then
                                phasetype(1) = phasetype(2)
                                phasetype(2) = 5
                                rho(phasetype(1)) = rho2
                                rho(phasetype(2)) = rho2_H
                                x_phase(:,phasetype(1)) = x_fluid2
                                x_phase(:,phasetype(2)) = x_hyd2
                                !vapfrac = vapfrac2
                                phasefrac(phasetype(1)) = vapfrac2
                                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                                errval = 0
                                !Save information of the equilibrium found
                                gl%solidtype_akt_phase = 0         !No solid
                                gl%solidpos_akt_phase = 0          !No solid
                            end if
                        else
                            !Obviously hydrate forms, but none of the phases checked is in equilibrium with hydrate
                            !So solid water has to form
                            check_Hydsolwater = .true.
                        end if
                        !rho(phasetype(2)) = 1000.D0
                        nrofphases = 2
                        errval1 = 0
                        errval2 = 0
                    end if

                    !Check, if hydrate and solid water are in equilibrium
                    gl%solid_pos = 1
                    if (.not. check_Hydsolwater) then
                        gl%molfractions = x_phase(:,phasetype(1))
                        call reduced_parameters_calc(gl,Temp)
                        call Chempot_CALC(gl,Temp, rho(phasetype(1)), Chempot, 0)
                        gl%molfractions = x
                        call reduced_parameters_calc(gl,Temp)
                        !tpd_sol = g_WaterIce(gl,Temp,press) - Chempot(gl%solid_pos)
                        tpd_sol_hyd(hdrt_structure) = g_WaterIce(gl,Temp,press) - Chempot(gl%solid_pos)
                    else
                        !tpd_sol = -1.D0    !Dummy value, so that the equilibrium is calculated
                        if (gl%Ncomp == 2) then ! IxH equilibrium only possible for pure hydrates since fugacity is unknown
                            tpd_sol_hyd(hdrt_structure) = -1.D0    !Dummy value, so that the equilibrium is calculated
                        else
                            nrofphases = nrfluidphases
                            errval = 0
                        endif
                    end if
                    !if ((tpd_sol < 0) .or. (vapfrac < 0) .or. (vapfrac > 1)) then
                    if (((tpd_sol_hyd(hdrt_structure) < 0) .or. (vapfrac < 0) .or. (vapfrac > 1)).and.(gl%Ncomp < 3)) then ! IxH equilibrium only possible for pure hydrates since fugacity is unknown
                        iFlash = 3
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.d0
                        gl%solidtype(1) = 1   !solid forms
                        gl%solidtype(2) = 1   !hydrate assumed
                        iphase = 0
                        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid, rhofluid_est, vapfrac, iFlash, iPhase, errval, iter)
                        if ((errval == 0) .and. (vapfrac > 0) .and. (vapfrac < 1)) then
                            phasetype(1) = 4
                            rho(phasetype(1)) = rho_sol(1)!1.D0 / v_WaterIce(gl,Temp, press)
                            x_phase(:,phasetype(1)) = x_solid
                            phasetype(2) = 5
                            rho(phasetype(2)) = rho_sol(2)
                            x_phase(:,phasetype(2)) = x_hyd
                            nrofphases = 2
                            phasefrac(phasetype(1)) = vapfrac
                            phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                            errval = 0
                            !Save information of the equilibrium found
                            gl%solidtype_akt_phase = 1             !Solid is water
                            gl%solidpos_akt_phase = gl%solid_pos      !Position of solid in fluid vector
                        else
                            !No stable equilbrium found, quit with error
                            if (errval_fluid /= 0) then
                                errval = errval_fluid
                            else
                                errval = -2233
                            end if
                            if (hdrt_structure == 1) cycle
                            return
                        end if
                    end if
                    check_Hydsolwater = .false.
                    if ((errval /= 0) .and. (dabs(tpd_sol_hyd(hdrt_structure)) < 100.d0) ) then !if close to 0 maybe no stable phase equilibrium was found
                        gl%hdrt_structure_stable = 0
                        nrofphases = nrfluidphases
                        errval = 0
                    endif
                end if
            enddo
        end if
        !END OF THE OLD ROUTINE
        !########################################################################################################################################################
    endif
    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)
    call hdrt_structure_definition(gl,gl%hdrt_structure_stable)

    end subroutine PhaseDet_sol



    ! Routines for the iterative (slow) two dimensional iteration for ph and ps flash calculations
    ! for mixtures. Solid formation of hydrate, water, and co2 is also considered
    ! Andreas, March 2014

    !************************************************************************************
    subroutine PhaseDet_ph_ps_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, z_spec, h_or_s, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for a specified mixture
    ! Formation of solid phases, i.e., hydrate, solid CO2, and solid H2O is also considered
    ! Variables:
    ! INPUT:
    !   press       - pressure [MPa]
    !   x           - Overall composition of the defined mixture
    !   z_spec      - Overall enthalpy or entropy of the mixture (z = s or h, specified in variable h_or_s)
    !   h_or_s      - Is enthalpy or entropy given? (1 = enthalpy, 2 = entropy)
    !
    ! OUTPUT:
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   phasetype   - contains the information about the phases in equilibrium
    !   x_phase(5,30)   - matrix that contains all phase compositions
    !       1   x_vap       - Vapor phase composition
    !       2   x_liq1      - (lighter) liquid phase composition
    !       3   x_liq2      - (heavier) liquid phase composition
    !       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !       5   x_hyd       - hydrate phase compsition
    !   Phasefrac           - Vector containing molar phase fractions
    !   nrofphases  - number of phases found
    !   errval      - Error value
    !   Temp        - Temperature [K]
    !************************************************************************************
    !Andreas March. 2014










    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, z_spec
    double precision, dimension(30) :: x
    double precision, dimension(5):: rho
    double precision, dimension(30,5) :: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, dimension(5):: phasefrac
    integer:: errval, nrofphases, h_or_s
    integer, dimension(5) :: phasetype

    integer :: i, j, k, mixflag, iphase
    double precision :: zmax, zmin, rhomax, rho_est
    double precision :: tredmix_orig, rhoredmix_orig
    integer, dimension(2) :: phasetype_min, phasetype_max   !These variables are only important for two component systems in which three phase equilibria exist
    !If the specified enthalpy is coincidently in the three phase area, the iterations will fail!!
    !The three phase area needs to be detected AND the corresponding temperature calculated
    integer :: nrofphases_min, nrofphases_max, iFlash, iter
    double precision, dimension(30) :: x_phase1_min, x_phase2_min, x_phase1_max, x_phase2_max, x_vap, x_liq1, x_liq2
    double precision:: rhovap_est, rholiq1_est, rholiq2_est, T_3phase, beta_V, beta_L1, beta_phase1, beta_phase2, beta_phase3
    double precision, dimension(3) :: phasefrac_2C3P
    integer :: solidtype_min, solidtype_max, solidpos_min, solidpos_max, solidtype_trial, solidpos_trial
    !double precision :: rho_phase1_min, rho_phase2_min, rho_phase1_max, rho_phase2_max

    double precision ::  d_phase, ChemPot_hyd
    double precision, dimension(3) :: z_phase

    double precision, dimension (60,60) :: Mat_phasefrac_2Comp
    double precision, dimension (60) :: Vec_phasefrac_2Comp, Vec_rightside_2Comp
    integer :: eqn, nrofsamephases
    character(255) :: herr

    logical :: bisec_3P2C
    double precision :: Ttrial, ztrial, rhofluid1_est, rhofluid2_est
    double precision, dimension(30) :: x_phase1_trial, x_phase2_trial, x_sol
    integer:: nrofphases_trial
    integer, dimension(2) :: phasetype_trial
    double precision, dimension(3) :: rho_sol
    integer:: iter_bisec

    !Variables needed for the calculation of hydrate enthalpy
    double precision, dimension(30) :: chem_pot, x_fluid, x_hyd, x_fluid1, x_fluid2     ! , fug_g
    double precision :: dchem_pot_dT_fluid, dchem_pot_dx_fluid, dfug_dT_fluid, dfug_dx_fluid, d_fluid
    double precision :: z_hyd

    integer:: pos_DryIce
    !double precision :: fug_CO2
    !double precision, dimension(2) :: occup, CiJ, xH

    !Variables neccessary for the regula falsi equilibrium pressure iteration
    !----------------------------------------------------------------------------
    integer :: Max_Iterations, Iterations
    double precision:: Tmin, Tmax, T_min_allowed, T_max_allowed, Delta_allowed
    type(type_additional_parameters)  :: parameters         ! needed to transmit T, rho_spec to Regula Falsi
    ! Mixed hydrates 2015
    double precision, dimension(3,30) :: occup, CiJ
    double precision, dimension(30) :: xH, fug_g        ! NO fug_CO2

    Max_Iterations = 50
    Delta_allowed = 1.D-6
    !parameters = 0.D0
    parameters%a_p(1) = z_spec
    parameters%a_p(2) = press
    parameters%a_p(3:32)= x
    errval = 0
    Iterations = 0
    !----------------------------------------------------------------------------

    Tmax = 0.D0
    Tmin = 0.D0
    phasetype = 0
    rho = 0.D0

    mixflag = 0

    nrofphases_min = 0
    nrofphases_max = 0
    phasetype_min = 0
    phasetype_max = 0
    x_phase1_min = 0.D0
    x_phase2_min = 0.D0
    x_phase1_max = 0.D0
    x_phase2_max = 0.D0
    !rho_phase1_min = 0.D0
    !rho_phase2_min = 0.D0
    !rho_phase1_max = 0.D0
    !rho_phase2_max = 0.D0
    bisec_3P2C = .true.

    chem_pot = 0.D0
    fug_g = 0.d0
    x_fluid = 0.d0
    dchem_pot_dT_fluid = 0.D0
    dchem_pot_dx_fluid = 0.D0
    dfug_dT_fluid = 0.d0
    dfug_dx_fluid = 0.d0
    d_fluid = 0.d0

    tredmix_orig = gl%tredmix
    rhoredmix_orig = gl%rhoredmix

    if (gl%ncomp == 1) then
        errval = -9901
        return
    end if

    if ((h_or_s < 1) .or. (h_or_s > 2)) then
        errval = -15567
    end if

    Do i = 1, gl%ncomp
        !Tmin = Tmin + tminfluid(i) * x(i)
        Tmax = Tmax + gl%tmaxfluid(i) * x(i)
    end do

    if (Tmax > 1500.D0) then
        Tmax = 1500.D0
    end if

    !Strategy:
    !1) Calculate the enthalpy at Tmax and a fairly high temperature, where no solids form (ca. 350 K) (solids = water, co2, hydrate)
    !2) Check whether the specified enthalpy is in between the calculated values
    !Case h(350K) < h_spec < h(Tmax) --> Use regula Falsi to find the correct temperature (SPECIAL CASE: 2 components, 3 phase equilibrium)
    !Case h_spec > h(Tmax) --> Quit with error
    !Case h_spec < h(350 K) --> Successively decrease the temperature in certain steps (30K ?) and try to find h(T_try) < h_spec
    !If error or T_try < 150 K, quit with error
    !Elsewise do regula falsi


    !Calculate the enthalpy at Tmax
    !********************************************************
    !Assume that only gas is present at Tmax (1 phase)
    gl%molfractions = x
    call reduced_parameters_calc(gl,Tmax)
    iphase = 2      !gas assumed
    rho_est = 0.D0
    rhomax = rhomix_calc(gl,Tmax, press, rho_est, IPHASE, mixflag)
    if (rhomax < 1.D-14) then
        errval = -8888
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig
        return
    end if
    if (h_or_s == 1) then   !Enthalpy
        zmax = h_calc(gl,Tmax, rhomax, mixflag)
    else                    !Entropy
        zmax = s_calc(gl,Tmax, rhomax, mixflag)
    end if
    if (z_spec > zmax) then
        errval = -19941
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig
        return
    end if
    nrofphases_max = 1
    phasetype_max(1) = 1
    !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
    solidtype_max = 0
    solidpos_max = 0
    !********************************************************

    !Calculate the enthalpy at Tmin
    !********************************************************
    Tmin = 340.020102010201D0  !K  Arbitrary starting temperature. No solids formation at that temperature
    do while (Tmin > 150.D0)
        call PhaseDet_sol(gl,press, Tmin, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)
        if (errval == 0) then
            if (gl%ncomp == 2) then
                nrofphases_min = nrofphases
                phasetype_min = phasetype(1:2)
                !rho_phase1_min = rho(phasetype(1))
                x_phase1_min = x_Phase(:,phasetype(1))
                if (nrofphases == 2) then
                    !rho_phase2_min = rho(phasetype(2))
                    x_phase2_min = x_phase(:,phasetype(2))
                else
                    !rho_phase2_min = 0.D0
                    x_phase2_min = 0.D0
                end if
            end if
            !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
            solidtype_min = gl%solidtype_akt_phase
            solidpos_min = gl%solidpos_akt_phase
            do i = 1, nrofphases
                if (phasetype(i) < 4) then  !Fluid phase
                    ! calculate the properties of the first phase
                    if (gl%ncomp > 1) then
                        gl%molfractions = x_phase(:,phasetype(i))
                        call reduced_parameters_calc(gl,Tmin) !Dummy temperature 300 K for the SRK
                    end if
                    d_phase = rho(phasetype(i))
                    if (h_or_s == 1) then   !Enthalpy
                        z_phase(i) = h_calc(gl,Tmin, d_phase, 0)
                    else                    !Entropy
                        z_phase(i) = s_calc(gl,Tmin, d_phase, 0)
                    end if

                elseif (phasetype(i) == 4) then     !Solid (pure)
                    if (solidtype_min == 1) then         !Water
                        if (h_or_s == 1) then   !Enthalpy
                            z_phase(i) = h_WaterIce(gl,Tmin, press)
                        else                    !Entropy
                            z_phase(i) = s_WaterIce(gl,Tmin, press)
                        end if
                    elseif (solidtype_min == 2) then     !CO2
                        if (h_or_s == 1) then   !Enthalpy
                            z_phase(i) = h_DryIce(gl,Tmin, press)
                        else                    !Entropy
                            z_phase(i) = s_DryIce(gl,Tmin, press)
                        end if
                    else                                !No solid equation available
                        errval = -9904
                    end if

                elseif (phasetype(i) == 5) then     !Hydrate

                    !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                    !calculate the entropy and enthalpy of hydrates
                    if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                        x_fluid = x_phase(:,phasetype(1))
                        d_fluid = rho(phasetype(1))
                        x_hyd = x_phase(:,phasetype(i))
                        !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                        !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                        !SH pure hdrt
                        !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                        !mix hdrt
                        call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                        if (errval == 0) then
                            if (h_or_s == 1) then   !Enthalpy
                                !SH pure hdrt
                                !call hdrt_enthalpy(Tmin, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                !mix hdrt
                                call hdrt_enthalpy(gl, Tmin, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                            else                    !Entropy
                                !SH pure hdrt
                                !call hdrt_entropy(Tmin, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                !mix hdrt
                                call hdrt_entropy(gl, Tmin, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                            end if
                        end if

                        if (errval /= 0) then
                            return
                        end if

                        z_phase(i) = z_hyd

                    else
                        if ((solidtype_min == 2) .and. (gl%ncomp == 2)) then !Dry ice in equilibrium with hydrate
                            pos_DryIce = solidpos_min
                            fug_g(1) = fug_DryIce(gl,Tmin,press, pos_dryIce) * 1.D6
                            if (fug_g(1) < 1.D-12) then
                                errval = -7779
                            else
                                !fug_CO2 = fug_g(1)
                                chem_pot(2) = g_DryIce(gl,Tmin,press)
                                call hdrt_chem_potent_w(gl,Tmin,press*1.d6,fug_g,ChemPot_hyd)
                                chem_pot(1) = ChemPot_hyd

                                x_hyd = x_phase(:,phasetype(i))

                                if (h_or_s == 1) then   !Enthalpy
                                    !SH pure hdrt
                                    !call hdrt_enthalpy(Tmin, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_enthalpy(gl, Tmin, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                else                    !Entropy
                                    !SH pure hdrt
                                    !call hdrt_entropy(Tmin, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_entropy(gl, Tmin, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                end if

                                if (errval /= 0) then
                                    return
                                end if
                                z_phase(i) = z_hyd
                            end if
                        else
                            errval = -12900
                        end if
                    end if

                end if
            end do
            if (errval == 0) then
                if (nrofphases == 3) then
                    zmin = z_phase(1) * Phasefrac(phasetype(1)) + z_phase(2) * Phasefrac(phasetype(2)) + z_phase(3) * Phasefrac(phasetype(3))
                elseif(nrofphases == 2) then
                    zmin = z_phase(1) * Phasefrac(phasetype(1)) + z_phase(2) * Phasefrac(phasetype(2))
                elseif(nrofphases == 1) then
                    zmin = z_phase(1)
                end if
                if (zmin < z_spec) then
                    exit
                else
                    Tmax = Tmin
                    nrofphases_max = nrofphases_min
                    phasetype_max = phasetype_min
                    !rho_phase1_max = rho_phase1_min
                    x_phase1_max = x_phase1_min
                    if (nrofphases == 2) then
                        !rho_phase2_max = rho_phase2_min
                        x_phase2_max = x_phase2_min
                    else
                        !rho_phase2_max = 0.D0
                        x_phase2_max = 0.D0
                    end if
                    !if (ncomp > 2) then
                    Tmin = Tmin - 5.D0
                    !else
                    !    Tmin = Tmin - 5.D0
                    !end if
                    !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                    solidtype_max = solidtype_min
                    solidpos_max = solidpos_min
                    if (Tmin < 120.D0) then
                        errval = -4405
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        return
                    end if
                end if
            else
                gl%tredmix = tredmix_orig
                gl%rhoredmix = rhoredmix_orig
                return
            end if
        else
            errval = -4405
            gl%tredmix = tredmix_orig
            gl%rhoredmix = rhoredmix_orig
            return
        end if
    end do
    !********************************************************

    if (zmin > z_spec) then
        errval = -19941
        return
    end if

    !For a binary mixture: try to find the three phase region in case hmin was found at a VLE and hmax at a LLE and hmin < h_spec < hmax
    !Several other three phase equilibrium lines are also possible:
    !For the binary system CO2 - water: VLwLc                               Only fluid phases                           phasetype: 1,2, and 3
    !                                   VLwH, LwLcH, VLcH                   2 fluid phases +  hydrate                   phasetype: two of (1,2,3) + 5
    !                                   VHIw, VHIc, LwHIw, LcHIc, LwHIc     1 fluid phase, 1 solid phase + hydrate      phasetype: one of (1,2,3) + 4 + 5
    !                                   VLwIw, VLcIc, LcLwIc                2 fluid phases + solid                      phasetype: two of (1,2,3) + 4
    !                                   (NOTE: The equilibrium IcH is not considered at the moment)
    !
    ! Strategy for finding the correct 3 phase equilibrium:
    ! 1) Check whether hmin and hmax were found at different two-phase equilibria. First check if a liquid phase at low temperature and a two phase equilibrium of different phases have been found at high temperature. If yes, do some safety bisection steps
    !   2a) If no, do not check for three phase line    --> exit
    !   2b) If yes, check whether both phases are different or only one
    !       3a) only one phase is different (e.g. 1+3 vs 2+3 --> check for VLcLw; 1+5 vs 1+4 --> check for VHIw equilbrium --> then exit
    !       3b) two phases different (e.g. 1+3 vs 2+5 --> do bisection until only only one phase different with hmin < h_spec < hmax --> go to 3a)

    if (gl%ncomp == 2) then

        !First check if different phases are present at the maximum and minimum temperature
        !If yes, do some saftey steps of bisection, Andreas, Feb 2015

        !k needs to be set to 0! Andreas, May 2018
        k = 0
        if ((nrofphases_min /= 2) .or. (nrofphases_max /= 2)) then


            !k = 0
            !Some two phase equilibria are only present in very narrow temperature ranges, e.g. VLc in the binary mixture CO2 + water. Safety first

            do while (bisec_3P2C)

                !Exit if the phases are the same or if at low and at high temperature two phases are present
                if ( ((nrofphases_min == 1) .and. (nrofphases_max == 1) .and. (phasetype_min(1) /= phasetype_max(1))) .or. &
                    &    ((nrofphases_min == 1) .and. (nrofphases_max == 2) .and. &
                    & (phasetype_min(1) /= phasetype_max(1)) .and. (phasetype_min(1) /= phasetype_max(2)))            .or. &
                    &    ((nrofphases_min == 2) .and. (nrofphases_max == 1) .and. &
                    & (phasetype_min(1) /= phasetype_max(1)) .and. (phasetype_min(2) /= phasetype_max(1)))  )          then

                else
                    exit
                end if


                !Do bisection
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                Ttrial = (Tmax + Tmin) / 2.D0
                call PhaseDet_sol(gl,press, Ttrial, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)
                if (errval == 0) then

                    nrofphases_trial = nrofphases
                    phasetype_trial = phasetype(1:2)
                    x_phase1_trial = x_Phase(:,phasetype(1))
                    if (nrofphases == 2) then
                        x_phase2_trial = x_phase(:,phasetype(2))
                    end if
                    !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                    solidtype_trial = gl%solidtype_akt_phase
                    solidpos_trial = gl%solidpos_akt_phase

                    do i = 1, nrofphases
                        ! calculate the properties of the first phase
                        if (phasetype(i) < 4) then  !Fluid phase
                            gl%molfractions = x_phase(:,phasetype(i))
                            call reduced_parameters_calc(gl,Ttrial) !Dummy temperature 300 K for the SRK
                            d_phase = rho(phasetype(i))
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(i) = h_calc(gl,Ttrial, d_phase, 0)
                            else                    !Entropy
                                z_phase(i) = s_calc(gl,Ttrial, d_phase, 0)
                            end if
                        elseif (phasetype(i) == 4) then     !Solid (pure)
                            if (solidtype_trial == 1) then         !Water
                                if (h_or_s == 1) then   !Enthalpy
                                    z_phase(i) = h_WaterIce(gl,Ttrial, press)
                                else
                                    z_phase(i) = s_WaterIce(gl,Ttrial, press)
                                end if
                            elseif (solidtype_trial == 2) then     !CO2
                                if (h_or_s == 1) then   !Enthalpy
                                    z_phase(i) = h_DryIce(gl,Ttrial, press)
                                else                    !Entropy
                                    z_phase(i) = s_DryIce(gl,Ttrial, press)
                                end if
                            else                                !No solid equation available
                                errval = -9904
                            end if
                        elseif (phasetype(i) == 5) then     !Hydrate
                            !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                            !calculate the entropy and enthalpy of hydrates
                            if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                                x_fluid = x_phase(:,phasetype(1))
                                d_fluid = rho(phasetype(1))
                                x_hyd = x_phase(:,phasetype(i))
                                !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                                !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                                !SH pure hdrt
                                !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                                !mix hdrt
                                call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                                if (errval == 0) then
                                    if (h_or_s == 1) then   !Enthalpy
                                        !SH pure hdrt
                                        !call hdrt_enthalpy(Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                        !mix hdrt
                                        call hdrt_enthalpy(gl, Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                    else                    !Entropy
                                        !SH pure hdrt
                                        !call hdrt_entropy(Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                        !mix hdrt
                                        call hdrt_entropy(gl, Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    end if
                                end if

                                if (errval /= 0) then
                                    return
                                end if

                                z_phase(i) = z_hyd

                            else
                                if ((solidtype_trial == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
                                    pos_DryIce = solidpos_trial
                                    fug_g(1) = fug_DryIce(gl,Ttrial,press, pos_dryIce) *1.D6
                                    if (fug_g(1) < 1.D-12) then
                                        errval = -7779
                                    else
                                        !fug_CO2 = fug_g(1)
                                        chem_pot(2) = g_DryIce(gl,Ttrial,press)
                                        call hdrt_chem_potent_w(gl,Ttrial,press*1.d6,fug_g,ChemPot_hyd)
                                        chem_pot(1) = ChemPot_hyd

                                        x_hyd = x_phase(:,phasetype(i))

                                        if (h_or_s == 1) then   !Enthalpy
                                            !SH pure hdrt
                                            !call hdrt_enthalpy(Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                            !mix hdrt
                                            call hdrt_enthalpy(gl, Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                        else                    !Entropy
                                            !SH pure hdrt
                                            !call hdrt_entropy(Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                            !mix hdrt
                                            call hdrt_entropy(gl, Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                        end if
                                        if (errval /= 0) then
                                            return
                                        end if

                                        z_phase(i) = z_hyd

                                    end if
                                else
                                    errval = -12900
                                end if
                            end if
                        end if
                    end do
                    if (errval == 0) then
                        if(nrofphases == 2) then
                            ztrial = z_phase(1) * Phasefrac(phasetype(1)) + z_phase(2) * Phasefrac(phasetype(2))
                        elseif(nrofphases == 1) then
                            ztrial = z_phase(1)
                        end if
                        if (ztrial < z_spec) then
                            zmin = ztrial
                            Tmin = Ttrial
                            phasetype_min = phasetype_trial
                            x_phase1_min = x_phase1_trial
                            x_phase2_min = x_phase2_trial
                            nrofphases_min = nrofphases_trial
                            !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                            solidtype_min = solidtype_trial
                            solidpos_min = solidpos_trial
                        else
                            zmax = ztrial
                            Tmax = Ttrial
                            phasetype_max = phasetype_trial
                            x_phase1_max = x_phase1_trial
                            x_phase2_max = x_phase2_trial
                            nrofphases_max = nrofphases_trial
                            !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                            solidtype_max = solidtype_trial
                            solidpos_max = solidpos_trial
                        end if
                    else
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        return
                    end if
                else
                    errval = -4405
                    gl%tredmix = tredmix_orig
                    gl%rhoredmix = rhoredmix_orig
                    return
                end if
                !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                k = k + 1
                if (k > 20) exit   !emergency exit, avoid endless iteration
            end do


        end if

        if ((nrofphases_min == 2) .and. (nrofphases_max == 2)) then


            !Check which phases are at equilibrium at hmin and hmax

            iter_bisec = 0
            bisec_3P2C = .true.
            if (gl%nrofhydrateformers < 2) then        !No hydrate forms
                iter_bisec = 0
            else                                    !Hydrate forms
                if (press < gl%p_Q_Hyd_2C(1)) then     !pressure lower than  lower quadruple point. VHIw or VLwIw may form
                    if (press > (0.95D0 * gl%p_Q_Hyd_2C(1))) then
                        iter_bisec = 8     !Do 8 steps of bisection if close to the quadruple point pressure
                    elseif (press > (0.6D0 * gl%p_Q_Hyd_2C(1))) then
                        iter_bisec = 5      !Do 5 steps of bisection if farer away from quadruple point pressure
                    elseif (press > (0.4D0 * gl%p_Q_Hyd_2C(1))) then
                        iter_bisec = 1      !Do 1 step of bisection if far away from quadruple point pressure
                    end if
                end if
                if ((gl%p_Q_Hyd_2C(2) > 1.D-12) .and. (press > gl%p_Q_Hyd_2C(2))) then !upper quadruple point exists and pressure is higher. VLwLx or LwLxH may form
                    if (press < (1.05D0 * gl%p_Q_Hyd_2C(2))) then
                        iter_bisec = 8     !Do 8 steps of bisection if close to the quadruple point pressure
                    elseif (press < (1.4D0 * gl%p_Q_Hyd_2C(2))) then
                        iter_bisec = 5      !Do 5 steps of bisection if farer away from quadruple point pressure
                    elseif (press < (1.6D0 * gl%p_Q_Hyd_2C(2))) then
                        iter_bisec = 1      !Do 1 step of bisection if far away from quadruple point pressure
                    end if
                end if
                !water - co2, Quadrupel point close to Dry ice tripel point
                if (gl%p_Q_Hyd_2C(3) > 1.D-12) then !Lower quadruple point close to dry ice triple point exists
                    if ((press > (0.8D0 * gl%p_Q_Hyd_2C(3))) .or. (press < (1.2D0 * gl%p_Q_Hyd_2C(3)))) then
                        iter_bisec = 8     !Do 8 steps of bisection if close to the quadruple point pressure
                    elseif ((press > (0.6D0 * gl%p_Q_Hyd_2C(3))) .or. (press < (1.4D0 * gl%p_Q_Hyd_2C(3)))) then
                        iter_bisec = 6      !Do 6 steps of bisection if farer away from quadruple point pressure
                    elseif ((press > (0.4D0 * gl%p_Q_Hyd_2C(3))) .or. (press < (1.6D0 * gl%p_Q_Hyd_2C(3)))) then
                        iter_bisec = 4      !Do 4 steps of bisection if far away from quadruple point pressure
                    elseif ((press > (0.2D0 * gl%p_Q_Hyd_2C(3))) .or. (press < (1.8D0 * gl%p_Q_Hyd_2C(3)))) then
                        iter_bisec = 2      !Do 2 steps of bisection if far away from quadruple point pressure
                    elseif ((press > (0.05D0 * gl%p_Q_Hyd_2C(3))) .or. (press < (2.0D0 * gl%p_Q_Hyd_2C(3)))) then
                        iter_bisec = 1      !Do 1 steps of bisection if far away from quadruple point pressure
                    end if
                end if
            end if

            !If bisection was done before, reduce the steps of bisection done here
            if (k < iter_bisec) then
                iter_bisec = iter_bisec - k
            else
                iter_bisec = 0
            end if
            k = 0
            do while (bisec_3P2C)

                nrofsamephases = 0
                do i=1,2
                    do j = 1,2
                        if (phasetype_min(i) == phasetype_max(j)) then                      !Check if the phase identifiers are the same
                            nrofsamephases = nrofsamephases + 1
                        end if
                    end do
                end do

                if (nrofsamephases == 2) then   !No special caution needs to be taken, if both phases are the same (No three phase lines inbetween)
                    exit
                end if

                ! Elsewise do at least iter_bisec steps of bisection, to find the correct solution in narrow areas, where 2 three phase lines like VHIw and VLwIw are close to each other
                ! For iter_bisec = 5 the two three phase lines are distiguishable if they have a difference of at least 30K / 2^iter_bisec = 0.94 K
                ! For iter_bisec = 10 the two three phase lines are distiguishable if they have a difference of at least 30K / 2^iter_bisec = 0.03 K
                if ((nrofsamephases > 0) .and. (k >= iter_bisec)) then
                    !If at least one phase is the same and the given number of bisection steps is reached, exit bisection method
                    bisec_3P2C = .false.
                else

                    !Do bisection until only one phase is different
                    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    Ttrial = (Tmax + Tmin) / 2.D0
                    call PhaseDet_sol(gl,press, Ttrial, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)
                    if (errval == 0) then

                        nrofphases_trial = nrofphases
                        phasetype_trial = phasetype(1:2)
                        if (nrofphases == 2) then
                            x_phase1_trial = x_Phase(:,phasetype(1))
                            x_phase2_trial = x_phase(:,phasetype(2))
                        end if
                        !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                        solidtype_trial = gl%solidtype_akt_phase
                        solidpos_trial = gl%solidpos_akt_phase

                        do i = 1, nrofphases
                            ! calculate the properties of the first phase
                            if (phasetype(i) < 4) then  !Fluid phase
                                gl%molfractions = x_phase(:,phasetype(i))
                                call reduced_parameters_calc(gl,Ttrial) !Dummy temperature 300 K for the SRK
                                d_phase = rho(phasetype(i))
                                if (h_or_s == 1) then   !Enthalpy
                                    z_phase(i) = h_calc(gl,Ttrial, d_phase, 0)
                                else                    !Entropy
                                    z_phase(i) = s_calc(gl,Ttrial, d_phase, 0)
                                end if
                            elseif (phasetype(i) == 4) then     !Solid (pure)
                                if (solidtype_trial == 1) then         !Water
                                    if (h_or_s == 1) then   !Enthalpy
                                        z_phase(i) = h_WaterIce(gl,Ttrial, press)
                                    else
                                        z_phase(i) = s_WaterIce(gl,Ttrial, press)
                                    end if
                                elseif (solidtype_trial == 2) then     !CO2
                                    if (h_or_s == 1) then   !Enthalpy
                                        z_phase(i) = h_DryIce(gl,Ttrial, press)
                                    else                    !Entropy
                                        z_phase(i) = s_DryIce(gl,Ttrial, press)
                                    end if
                                else                                !No solid equation available
                                    errval = -9904
                                end if
                            elseif (phasetype(i) == 5) then     !Hydrate
                                !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                                !calculate the entropy and enthalpy of hydrates
                                if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                                    x_fluid = x_phase(:,phasetype(1))
                                    d_fluid = rho(phasetype(1))
                                    x_hyd = x_phase(:,phasetype(i))
                                    !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                                    !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                                    !SH pure hdrt
                                    !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                                    !mix hdrt
                                    call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                                    if (errval == 0) then
                                        if (h_or_s == 1) then   !Enthalpy
                                            !SH pure hdrt
                                            !call hdrt_enthalpy(Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                            !mix hdrt
                                            call hdrt_enthalpy(gl, Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                        else                    !Entropy
                                            !SH pure hdrt
                                            !call hdrt_entropy(Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                            !mix hdrt
                                            call hdrt_entropy(gl, Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                        end if
                                    end if

                                    if (errval /= 0) then
                                        return
                                    end if

                                    z_phase(i) = z_hyd

                                else
                                    if ((solidtype_trial == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
                                        pos_DryIce = solidpos_trial
                                        fug_g(1) = fug_DryIce(gl,Ttrial,press, pos_dryIce) *1.D6
                                        if (fug_g(1) < 1.D-12) then
                                            errval = -7779
                                        else
                                            !fug_CO2 = fug_g(1)
                                            chem_pot(2) = g_DryIce(gl,Ttrial,press)
                                            call hdrt_chem_potent_w(gl,Ttrial,press*1.d6,fug_g,ChemPot_hyd)
                                            chem_pot(1) = ChemPot_hyd

                                            x_hyd = x_phase(:,phasetype(i))

                                            if (h_or_s == 1) then   !Enthalpy
                                                !SH pure hdrt
                                                !call hdrt_enthalpy(Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                                !mix hdrt
                                                call hdrt_enthalpy(gl, Ttrial, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                            else                    !Entropy
                                                !SH pure hdrt
                                                !call hdrt_entropy(Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                                !mix hdrt
                                                call hdrt_entropy(gl, Ttrial, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                            end if
                                            if (errval /= 0) then
                                                return
                                            end if

                                            z_phase(i) = z_hyd

                                        end if
                                    else
                                        errval = -12900
                                    end if
                                end if
                            end if
                        end do
                        if (errval == 0) then
                            if(nrofphases == 2) then
                                ztrial = z_phase(1) * Phasefrac(phasetype(1)) + z_phase(2) * Phasefrac(phasetype(2))
                            elseif(nrofphases == 1) then
                                ztrial = z_phase(1)
                            end if
                            if (ztrial < z_spec) then
                                zmin = ztrial
                                Tmin = Ttrial
                                phasetype_min = phasetype_trial
                                x_phase1_min = x_phase1_trial
                                x_phase2_min = x_phase2_trial
                                nrofphases_min = nrofphases_trial
                                !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                                solidtype_min = solidtype_trial
                                solidpos_min = solidpos_trial
                            else
                                zmax = ztrial
                                Tmax = Ttrial
                                phasetype_max = phasetype_trial
                                x_phase1_max = x_phase1_trial
                                x_phase2_max = x_phase2_trial
                                nrofphases_max = nrofphases_trial
                                !Save the information which solid forms (if at all) and at which position the solid former is in the fluids vector
                                solidtype_max = solidtype_trial
                                solidpos_max = solidpos_trial
                            end if
                            !If in the course of the iteration only one stable phase is found, quit this routine and continue
                            if (nrofphases == 1) then
                                exit
                            end if
                        else
                            gl%tredmix = tredmix_orig
                            gl%rhoredmix = rhoredmix_orig
                            return
                        end if
                    else
                        errval = -4405
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        return
                    end if
                    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                end if

                k = k + 1
                if (k > 100) exit   !emergency exit, avoid endless iteration
            end do

            !Check whether specified enthalpy or entropy is in the three phase region
            if ((nrofsamephases == 1) .and. (nrofphases_min == 2) .and. (nrofphases_max == 2))  then

                !Check which three phase region the enthalpy is in
                if ((phasetype_min(1) == 2) .and. (phasetype_min(2) == 3) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) then !VL and LL
                    x_vap = x_phase1_max
                    x_liq1 = x_phase1_min
                    x_liq2 = x_phase2_max
                    rhovap_est = 0.D0
                    rholiq1_est = 0.D0
                    rholiq2_est = 0.D0
                    T_3phase = (tmin + tmax) / 2    !Start value for the three phase temperature iteration
                    iFlash = 1 !p given, Bubble point
                    call ptflash_NC_3P(gl,press, T_3phase, x, rho, x_vap, x_liq1, x_liq2, rhovap_est, &
                        & rholiq1_est, rholiq2_est, phasefrac_2C3P, iFlash, iter, errval)

                    if (errval == 0) then
                        phasetype(1) = 1
                        phasetype(2) = 2
                        phasetype(3) = 3
                        x_phase(:,phasetype(1)) = x_vap
                        x_phase(:,phasetype(2)) = x_liq1
                        x_phase(:,phasetype(3)) = x_liq2
                        nrofphases = 3
                        !Compute the enthalpies of the 3 phases
                        do i = 1, nrofphases
                            ! calculate the properties of the first phase
                            gl%molfractions = x_phase(:,phasetype(i))
                            call reduced_parameters_calc(gl,T_3phase) !Dummy temperature 300 K for the SRK
                            !Fluid phase
                            d_phase = rho(phasetype(i))
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(i) = h_calc(gl,T_3phase, d_phase, 0)
                            else                    !Entropy
                                z_phase(i) = s_calc(gl,T_3phase, d_phase, 0)
                            end if
                        end do
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig

                        !Calculate the enthalpies on the phase boundaries: VL2 ->  beta_L1 = 0  -->  beta_V = (z1 - x1_L2) / (x1_V - x1_L2)
                        !Calculate the enthalpies on the phase boundaries: L1L2 -> beta_V = 0  -->  beta_L1 = (z1 - x1_L2) / (x1_L1 - x1_L2)
                        beta_V = (x(1) - x_liq2(1)) / (x_vap(1) - x_liq2(1))
                        beta_L1 = (x(1) - x_liq2(1)) / (x_liq1(1) - x_liq2(1))
                        zmin = beta_L1 * z_phase(2) + (1.D0 - beta_L1) * z_phase(3)
                        zmax = beta_V * z_phase(1) + (1.D0 - beta_V) * z_phase(3)

                        if ((zmin <= z_spec) .and. (zmax >= z_spec)) then

                            Temp = T_3phase
                            !Solve the system of equations:
                            ! x1V * beta_V + x1L1 * beta_L1 + X1L2 * beta_L2 = z1
                            ! beta_V + beta_L1 + beta_L2 = 1
                            ! beta_V * h_V + beta_L1 * h_L1 + beta_L2 * h_L2 = h_spec
                            ! for the phase fractions


                            Mat_phasefrac_2Comp(1,1) = x_vap(1)
                            Mat_phasefrac_2Comp(2,1) = x_liq1(1)
                            Mat_phasefrac_2Comp(3,1) = x_liq2(1)

                            Mat_phasefrac_2Comp(1,2) = z_phase(1)
                            Mat_phasefrac_2Comp(2,2) = z_phase(2)
                            Mat_phasefrac_2Comp(3,2) = z_phase(3)

                            Mat_phasefrac_2Comp(1,3) = 1.D0
                            Mat_phasefrac_2Comp(2,3) = 1.D0
                            Mat_phasefrac_2Comp(3,3) = 1.D0

                            Vec_rightside_2Comp(1) = x(1)
                            Vec_rightside_2Comp(2) = z_spec
                            Vec_rightside_2Comp(3) = 1.D0

                            Vec_phasefrac_2Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_2Comp,Vec_rightside_2Comp,errval,herr)
                            if (errval == 0) then
                                phasefrac(1:3) = Vec_rightside_2Comp(1:3)
                                return
                            end if
                        end if

                    end if

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !Two fluid phases + hydrate
                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                elseif  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) .or. &            !VH  / VLw      or
                    &     ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) .or. &            !LwH / VLw      or
                    &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 2)) .or. &            !VH  / VLc      or
                    &     ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) .or. &            !VH  / LcH      or
                    &     ((phasetype_min(1) == 2) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) .or. &            !VH  / LcH      or (Backup if Lc has ID 2)
                    &     ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 3)) .or. &            !LwH / LcLw     or
                    &     ((phasetype_min(1) == 2) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 3)) ) then            !LcH / LcLw

                    !in each of the 7 cases above, except for VH/LcH , phasetype_min(1) is the phase that shows up twice --> take Phasetype_max + hydrate as initial guess
                    if (phasetype_max(2) /= 5) then
                        x_fluid1 = x_phase1_max
                        x_fluid2 = x_phase2_max
                    else
                        x_fluid1 = x_phase1_max
                        x_fluid2 = x_phase1_min
                    end if
                    x_sol = 0.D0
                    x_hyd = 0.D0
                    rho_sol = 0.D0
                    rhofluid1_est = 0.D0
                    rhofluid2_est = 0.D0
                    gl%solidtype(1) = 0    !no other solid forms
                    gl%solidtype(2) = 1    !hydrate forms
                    iFlash = 1          !T given, melting point

                    if (phasetype_max(1) == 1) then
                        iphase = 2  !First phase is vapor
                    else
                        iphase = 1  !First phase is liquid
                    end if

                    T_3phase = (tmin + tmax) / 2    !Start value for the three phase temperature iteration
                    iFlash = 1 !p given, Melting point
                    call ptflash_solid_NC_3P(gl,press, T_3phase, x, rho_sol, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_2C3P, iFlash, iphase, iter, errval)

                    if (errval == 0) then
                        phasetype(1) = phasetype_max(1)
                        if (phasetype_max(2) /= 5) then
                            phasetype(2) = phasetype_max(2)
                        else
                            phasetype(2) = phasetype_min(1)
                        end if
                        phasetype(3) = 5
                        x_phase(:,phasetype(1)) = x_fluid1
                        x_phase(:,phasetype(2)) = x_fluid2
                        x_phase(:,phasetype(3)) = x_hyd
                        rho(phasetype(1)) = rho_sol(1)
                        rho(phasetype(2)) = rho_sol(2)
                        rho(phasetype(3)) = rho_sol(3)
                        nrofphases = 3
                        !Compute the enthalpies of the 3 phases
                        do i = 1, 2 !First 2 phases are fluid phases
                            ! calculate the properties of the first phase
                            gl%molfractions = x_phase(:,phasetype(i))
                            call reduced_parameters_calc(gl,T_3phase) !Dummy temperature 300 K for the SRK
                            !Fluid phase
                            d_phase = rho(phasetype(i))
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(i) = h_calc(gl,T_3phase, d_phase, 0)
                            else                    !Entropy
                                z_phase(i) = s_calc(gl,T_3phase, d_phase, 0)
                            end if
                        end do
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        !Calculate the enthalpy of the hydrate phase
                        !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                        !calculate the entropy and enthalpy of hydrates
                        if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                            x_fluid = x_phase(:,phasetype(1))
                            d_fluid = rho(phasetype(1))
                            x_hyd = x_phase(:,phasetype(3))
                            !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                            !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                            !SH pure hdrt
                            !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                            !mix hdrt
                            call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                            if (errval == 0) then
                                if (h_or_s == 1) then   !Enthalpy
                                    !SH pure hdrt
                                    !call hdrt_enthalpy(T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_enthalpy(gl, T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                else                    !Entropy
                                    !SH pure hdrt
                                    !call hdrt_entropy(T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_entropy(gl, T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                end if
                            end if

                            if (errval /= 0) then
                                return
                            end if
                            z_phase(3) = z_hyd

                        else
                            errval = -12900
                        end if

                        !Calculate the enthalpies on the phase boundaries:
                        !VH  / VLw   ->  Phase 1 + Phase 3 <-> Phase 1 + Phase 2    beta_Phase2 = 0 : beta_Phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3) <->  beta_Phase3 = 0 : beta_phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2)
                        !VH  / VLc   ->  Phase 1 + Phase 3 <-> Phase 1 + Phase 2
                        !LwH / LcLw  ->  Phase 1 + Phase 3 <-> Phase 1 + Phase 2
                        !LcH / LcLw  ->  Phase 1 + Phase 3 <-> Phase 1 + Phase 2

                        if  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) .or. &            !VH  / VLw      or
                            &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 2)) .or. &            !VH  / VLc      or (DOES THIS PHASECHANGE EXIST??)
                            &     ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 3)) .or. &            !LwH / LcLw     or
                            &     ((phasetype_min(1) == 2) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 3)) ) then            !LcH / LcLw     or

                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        !LwH / VLw   ->  Phase 2 + Phase 3 <-> Phase 1 + Phase 2     beta_Phase1 = 0 : beta_Phase2 = (z1 - x1_Phase3) / (x1_Phase2 - x1_Phase3) <->  beta_Phase3 = 0 : beta_phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2)

                        if  ( ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) ) then            !LcH / LcLw

                            beta_phase2 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(2)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase2 * z_phase(2) + (1.D0 - beta_phase2) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        !VH  / LcH   ->  Phase 1 + Phase 3 <-> Phase 2 + Phase 3    beta_Phase2 = 0 : beta_Phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3) <->  beta_Phase1 = 0 : beta_phase2 = (z1 - x1_Phase3) / (x1_Phase2 - x1_Phase3)

                        if  ( ((phasetype_min(1) == 2) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) .or. &            !VH  / LcH      or (Backup if Lc has ID 2)
                            &     ((phasetype_min(1) == 3) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) ) then            !VH  / LcH      or


                            beta_phase2 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(2)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase2 * z_phase(2) + (1.D0 - beta_phase2) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                        end if

                        if ((zmin <= z_spec) .and. (zmax >= z_spec)) then

                            Temp = T_3phase
                            !Solve the system of equations:
                            ! x1_Phase1 * beta_Phase1 + x1_Phase2 * beta_Phase2 + X1_Phase3 * beta_Phase3 = z1
                            ! beta_Phase1 + beta_Phase2 + beta_Phase3 = 1
                            ! beta_Phase1 * h_Phase1 + beta_Phase2 * h_Phase2 + beta_Phase3 * h_Phase3 = h_spec
                            ! for the phase fractions


                            Mat_phasefrac_2Comp(1,1) = x_phase(1,phasetype(1))
                            Mat_phasefrac_2Comp(2,1) = x_phase(1,phasetype(2))
                            Mat_phasefrac_2Comp(3,1) = x_phase(1,phasetype(3))

                            Mat_phasefrac_2Comp(1,2) = z_phase(1)
                            Mat_phasefrac_2Comp(2,2) = z_phase(2)
                            Mat_phasefrac_2Comp(3,2) = z_phase(3)

                            Mat_phasefrac_2Comp(1,3) = 1.D0
                            Mat_phasefrac_2Comp(2,3) = 1.D0
                            Mat_phasefrac_2Comp(3,3) = 1.D0

                            Vec_rightside_2Comp(1) = x(1)
                            Vec_rightside_2Comp(2) = z_spec
                            Vec_rightside_2Comp(3) = 1.D0

                            Vec_phasefrac_2Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_2Comp,Vec_rightside_2Comp,errval,herr)
                            if (errval == 0) then
                                !phasefrac(1:3) = Vec_rightside_2Comp(1:3)
                                phasefrac(phasetype(1)) = Vec_rightside_2Comp(1)
                                phasefrac(phasetype(2)) = Vec_rightside_2Comp(2)
                                phasefrac(phasetype(3)) = Vec_rightside_2Comp(3)
                                return
                            end if
                        end if

                    end if

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !Two fluid phases + hydrate (Special case VLcH) Andreas Feb 2015
                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                elseif  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 3) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5))   .or. &            !VLc / VH       or Andreas, Feb 2015
                    &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 2) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) ) then              !VLc / VH       or (Backup if Lc has ID 2) Andreas, Feb 2015

                    !in each of the 2 cases above, phasetype_min(1) is the phase that shows up twice --> take Phasetype_min + hydrate as initial guess
                    x_fluid1 = x_phase1_min
                    x_fluid2 = x_phase2_min

                    x_sol = 0.D0
                    x_hyd = 0.D0
                    rho_sol = 0.D0
                    rhofluid1_est = 0.D0
                    rhofluid2_est = 0.D0
                    gl%solidtype(1) = 0    !no other solid forms
                    gl%solidtype(2) = 1    !hydrate forms
                    iFlash = 1          !T given, melting point

                    if (phasetype_max(1) == 1) then
                        iphase = 2  !First phase is vapor
                    else
                        iphase = 1  !First phase is liquid
                    end if

                    T_3phase = (tmin + tmax) / 2    !Start value for the three phase temperature iteration
                    iFlash = 1 !p given, Melting point
                    call ptflash_solid_NC_3P(gl,press, T_3phase, x, rho_sol, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_2C3P, iFlash, iphase, iter, errval)

                    if (errval == 0) then
                        phasetype(1) = phasetype_min(1)
                        phasetype(2) = phasetype_min(2)
                        phasetype(3) = 5
                        x_phase(:,phasetype(1)) = x_fluid1
                        x_phase(:,phasetype(2)) = x_fluid2
                        x_phase(:,phasetype(3)) = x_hyd
                        rho(phasetype(1)) = rho_sol(1)
                        rho(phasetype(2)) = rho_sol(2)
                        rho(phasetype(3)) = rho_sol(3)
                        nrofphases = 3
                        !Compute the enthalpies of the 3 phases
                        do i = 1, 2 !First 2 phases are fluid phases
                            ! calculate the properties of the first phase
                            gl%molfractions = x_phase(:,phasetype(i))
                            call reduced_parameters_calc(gl,T_3phase) !Dummy temperature 300 K for the SRK
                            !Fluid phase
                            d_phase = rho(phasetype(i))
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(i) = h_calc(gl,T_3phase, d_phase, 0)
                            else                    !Entropy
                                z_phase(i) = s_calc(gl,T_3phase, d_phase, 0)
                            end if
                        end do
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        !Calculate the enthalpy of the hydrate phase
                        !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                        !calculate the entropy and enthalpy of hydrates
                        if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                            x_fluid = x_phase(:,phasetype(1))
                            d_fluid = rho(phasetype(1))
                            x_hyd = x_phase(:,phasetype(3))
                            !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                            !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                            !SH pure hdrt
                            !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                            !mix hdrt
                            call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                            if (errval == 0) then
                                if (h_or_s == 1) then   !Enthalpy
                                    !SH pure hdrt
                                    !call hdrt_enthalpy(T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_enthalpy(gl, T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                else                    !Entropy
                                    !SH pure hdrt
                                    !call hdrt_entropy(T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_entropy(gl, T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                end if
                            end if

                            if (errval /= 0) then
                                return
                            end if
                            z_phase(3) = z_hyd

                        else
                            errval = -12900
                        end if


                        !Calculate the enthalpies on the phase boundaries:
                        !VH  / VLc   ->  Phase 1 + Phase 3 <-> Phase 1 + Phase 2    beta_Phase2 = 0 : beta_Phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3) <->  beta_Phase3 = 0 : beta_phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2)

                        if  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 3) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5))   .or. &            !VH  / VLc      or
                            &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 2) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) ) then              !VH  / VLc      or (Backup if Lc has ID 2)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmin = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                        end if


                        if ((zmin <= z_spec) .and. (zmax >= z_spec)) then

                            Temp = T_3phase
                            !Solve the system of equations:
                            ! x1_Phase1 * beta_Phase1 + x1_Phase2 * beta_Phase2 + X1_Phase3 * beta_Phase3 = z1
                            ! beta_Phase1 + beta_Phase2 + beta_Phase3 = 1
                            ! beta_Phase1 * h_Phase1 + beta_Phase2 * h_Phase2 + beta_Phase3 * h_Phase3 = h_spec
                            ! for the phase fractions


                            Mat_phasefrac_2Comp(1,1) = x_phase(1,phasetype(1))
                            Mat_phasefrac_2Comp(2,1) = x_phase(1,phasetype(2))
                            Mat_phasefrac_2Comp(3,1) = x_phase(1,phasetype(3))

                            Mat_phasefrac_2Comp(1,2) = z_phase(1)
                            Mat_phasefrac_2Comp(2,2) = z_phase(2)
                            Mat_phasefrac_2Comp(3,2) = z_phase(3)

                            Mat_phasefrac_2Comp(1,3) = 1.D0
                            Mat_phasefrac_2Comp(2,3) = 1.D0
                            Mat_phasefrac_2Comp(3,3) = 1.D0

                            Vec_rightside_2Comp(1) = x(1)
                            Vec_rightside_2Comp(2) = z_spec
                            Vec_rightside_2Comp(3) = 1.D0

                            Vec_phasefrac_2Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_2Comp,Vec_rightside_2Comp,errval,herr)
                            if (errval == 0) then
                                !phasefrac(1:3) = Vec_rightside_2Comp(1:3)
                                phasefrac(phasetype(1)) = Vec_rightside_2Comp(1)
                                phasefrac(phasetype(2)) = Vec_rightside_2Comp(2)
                                phasefrac(phasetype(3)) = Vec_rightside_2Comp(3)
                                return
                            end if
                        end if

                    end if

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !fluid phases + solid phase + hydrate
                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                elseif  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 4)) .or. &            !VIw  / VH      or
                    &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 4)) .or. &            !VIw  / HIw     or
                    &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 4)) .or. &            !LwIw / HIw     or
                    &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 5)) .or. &            !LwH  / HIw     or
                    &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) .or. &            !VH   / HIc     or
                    &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 5)) ) then            !LcH  / HIc

                    !in each of the 4 cases above, phasetype_max(1) is the fluid phase --> take Phasetype_max(1) + solid + hydrate as initial guess

                    !No dry ice / hydrate equilibria implemented, thus VHIc and LcHIc cannot be calculated yet!!
                    !No enthalpy calculation possible for HIw equilibria at the moment!!

                    x_fluid1 = x_phase1_max

                    x_sol = 0.D0

                    if (phasetype_min(1) == 4) then
                        gl%solidtype(1) = solidtype_min
                        gl%solid_pos = solidpos_min
                    elseif(phasetype_max(2) == 4) then
                        gl%solidtype(1) = solidtype_max
                        gl%solid_pos = solidpos_max
                    else
                        errval = -15567
                        return
                    end if


                    x_sol(gl%solid_pos) = 1.D0 !solid water + hydrate
                    x_hyd = 0.D0
                    rho_sol = 0.D0
                    rhofluid1_est = 0.D0
                    rhofluid2_est = 0.D0
                    !solidtype(1) = 1    !Solid water forms
                    gl%solidtype(2) = 1    !hydrate forms
                    iFlash = 1          !T given, melting point

                    if (phasetype_max(1) == 1) then
                        iphase = 2  !First phase is vapor
                    else
                        iphase = 1  !First phase is liquid
                    end if

                    T_3phase = (tmin + tmax) / 2    !Start value for the three phase temperature iteration
                    iFlash = 1 !p given, Melting point
                    call ptflash_solid_NC_3P(gl,press, T_3phase, x, rho_sol, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_2C3P, iFlash, iphase, iter, errval)

                    if (errval == 0) then
                        phasetype(1) = phasetype_max(1)
                        phasetype(2) = 4
                        phasetype(3) = 5
                        x_phase(:,phasetype(1)) = x_fluid1
                        x_phase(:,phasetype(2)) = x_sol
                        x_phase(:,phasetype(3)) = x_hyd
                        rho(phasetype(1)) = rho_sol(1)
                        rho(phasetype(2)) = rho_sol(2)
                        rho(phasetype(3)) = rho_sol(3)
                        nrofphases = 3
                        !Compute the enthalpies of the 3 phases
                        ! calculate the properties of the first phase (fluid)
                        gl%molfractions = x_phase(:,phasetype(1))
                        call reduced_parameters_calc(gl,T_3phase) !Dummy temperature 300 K for the SRK
                        !Fluid phase
                        d_phase = rho(phasetype(1))
                        if (h_or_s == 1) then   !Enthalpy
                            z_phase(1) = h_calc(gl,T_3phase, d_phase, 0)
                        else                    !Entropy
                            z_phase(1) = s_calc(gl,T_3phase, d_phase, 0)
                        end if
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        !Calculate the enthalpy or entropy of the solid phase
                        if (gl%solidtype(1) == 1) then
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(2) = h_WaterIce(gl,T_3phase, press)
                            else                    !Entropy
                                z_phase(2) = s_WaterIce(gl,T_3phase, press)
                            end if
                        else
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(2) = h_DryIce(gl,T_3phase, press)
                            else                    !Entropy
                                z_phase(2) = s_DryIce(gl,T_3phase, press)
                            end if
                        end if
                        !Calculate the enthalpy of the hydrate phase
                        !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                        !calculate the entropy and enthalpy of hydrates
                        if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                            x_fluid = x_phase(:,phasetype(1))
                            d_fluid = rho(phasetype(1))
                            x_hyd = x_phase(:,phasetype(3))
                            !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                            !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                            !SH pure hdrt
                            !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                            !mix hdrt
                            call hdrt_ancillary_hs(gl, Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)

                            if (errval == 0) then
                                if (h_or_s == 1) then   !Enthalpy
                                    !SH pure hdrt
                                    !call hdrt_enthalpy(T_3phase, press*1.D6, x_hyd, chem_pot, fug_g,  z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_enthalpy(gl, T_3phase, press*1.D6, x_hyd, chem_pot, fug_g,  z_hyd, errval)
                                else                    !Entropy
                                    !SH pure hdrt
                                    !call hdrt_entropy(T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    !mix hdrt
                                    call hdrt_entropy(gl, T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                end if
                            end if

                            if (errval /= 0) then
                                return
                            end if

                            z_phase(3) = z_hyd

                        else
                            if ((gl%solidtype(1) == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
                                pos_DryIce = 2
                                fug_g(1) = fug_DryIce(gl,Tmin,press, pos_dryIce) *1.D6
                                if (fug_g(1) < 1.D-12) then
                                    errval = -7779
                                else
                                    !fug_CO2 = fug_g(1)
                                    chem_pot(2) = g_DryIce(gl,Tmin,press)
                                    call hdrt_chem_potent_w(gl,Tmin,press*1.d6,fug_g,ChemPot_hyd)
                                    chem_pot(1) = ChemPot_hyd

                                    x_hyd = x_phase(:,phasetype(i))

                                    if (h_or_s == 1) then   !Enthalpy
                                        !SH pure hdrt
                                        !call hdrt_enthalpy(T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                        !mix hdrt
                                        call hdrt_enthalpy(gl, T_3phase, press*1.D6, x_hyd, chem_pot, fug_g, z_hyd, errval)
                                    else                    !Entropy
                                        !SH pure hdrt
                                        !call hdrt_entropy(T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                        !mix hdrt
                                        call hdrt_entropy(gl, T_3phase, press*1.D6, x_hyd, fug_g, z_hyd, errval)
                                    end if
                                end if
                            else
                                errval = -12900
                            end if
                        end if



                        !Calculate the enthalpies on the phase boundaries:
                        !VIw  / VH   ->  Phase 1 + Phase 2 <-> Phase 1 + Phase 3    beta_Phase3 = 0 : beta_Phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2) <->  beta_Phase2 = 0 : beta_phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3)

                        if  ( ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 4)) ) then    !VIw  / VH


                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        !VIw  / HIw  ->  Phase 1 + Phase 2 <-> Phase 2 + Phase 3     beta_Phase3 = 0 : beta_Phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2) <->  beta_Phase1 = 0 : beta_phase2 = (z1 - x1_Phase3) / (x1_Phase2 - x1_Phase3)
                        !LwIw / HIw  ->  Phase 1 + Phase 2 <-> Phase 2 + Phase 3

                        if  ( ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 4)) .or. &    !VIw  / HIw     or
                            &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 4)) ) then    !LwIw / HIw

                            beta_phase2 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(2)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase2 * z_phase(2) + (1.D0 - beta_phase2) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        !LwH  / HIw  ->  Phase 1 + Phase 3 <-> Phase 2 + Phase 3    beta_Phase2 = 0 : beta_Phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3) <->  beta_Phase1 = 0 : beta_phase2 = (z1 - x1_Phase3) / (x1_Phase2 - x1_Phase3)
                        !VH   / HIc  ->  Phase 1 + Phase 3 <-> Phase 2 + Phase 3
                        !LcH  / HIc  ->  Phase 1 + Phase 3 <-> Phase 2 + Phase 3

                        if  ( ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 5)) .or. &    !LwH  / HIw     or
                            &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 5)) .or. &    !VH   / HIc     or
                            &     ((phasetype_min(1) == 4) .and. (phasetype_min(2) == 5) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 5)) ) then    !LcH  / HIc

                            beta_phase2 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(2)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase2 * z_phase(2) + (1.D0 - beta_phase2) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                        end if


                        if ((zmin <= z_spec) .and. (zmax >= z_spec)) then

                            Temp = T_3phase
                            !Solve the system of equations:
                            ! x1_Phase1 * beta_Phase1 + x1_Phase2 * beta_Phase2 + X1_Phase3 * beta_Phase3 = z1
                            ! beta_Phase1 + beta_Phase2 + beta_Phase3 = 1
                            ! beta_Phase1 * h_Phase1 + beta_Phase2 * h_Phase2 + beta_Phase3 * h_Phase3 = h_spec
                            ! for the phase fractions


                            Mat_phasefrac_2Comp(1,1) = x_phase(1,phasetype(1))
                            Mat_phasefrac_2Comp(2,1) = x_phase(1,phasetype(2))
                            Mat_phasefrac_2Comp(3,1) = x_phase(1,phasetype(3))

                            Mat_phasefrac_2Comp(1,2) = z_phase(1)
                            Mat_phasefrac_2Comp(2,2) = z_phase(2)
                            Mat_phasefrac_2Comp(3,2) = z_phase(3)

                            Mat_phasefrac_2Comp(1,3) = 1.D0
                            Mat_phasefrac_2Comp(2,3) = 1.D0
                            Mat_phasefrac_2Comp(3,3) = 1.D0

                            Vec_rightside_2Comp(1) = x(1)
                            Vec_rightside_2Comp(2) = z_spec
                            Vec_rightside_2Comp(3) = 1.D0

                            Vec_phasefrac_2Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_2Comp,Vec_rightside_2Comp,errval,herr)
                            if (errval == 0) then
                                !phasefrac(1:3) = Vec_rightside_2Comp(1:3)
                                phasefrac(phasetype(1)) = Vec_rightside_2Comp(1)
                                phasefrac(phasetype(2)) = Vec_rightside_2Comp(2)
                                phasefrac(phasetype(3)) = Vec_rightside_2Comp(3)
                                gl%solidtype_akt_phase = gl%solidtype(1)
                                gl%solidpos_akt_phase = gl%solid_pos
                                return
                            end if
                        end if

                    end if

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !two fluid phases + solid phase
                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                elseif (  ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) .or. &            !VIw / VLw      or  (VIc / VLc )
                    &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 4)) .or. &            !VIw / LwIw     or  (VIc / LcIc)
                    &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 2)) .or. &            !VIc / VLc      or  (Backup)
                    &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 4)) ) then            !VIc / LcIc         (Backup)


                    !In each of the four cases above, phasetype_min(1) is vapor
                    x_fluid1 = x_phase1_min
                    !the liquid phase is either stored in phasetype_max(1) oder phasetype_max(2)
                    if (phasetype_max(2) == 4) then
                        x_fluid2 = x_phase1_max
                    else
                        x_fluid2 = x_phase2_max
                    end if

                    x_sol = 0.D0

                    gl%solidtype(1) = solidtype_min
                    gl%solid_pos = solidpos_min

                    x_sol(gl%solid_pos) = 1.D0 !solid water + hydrate

                    x_hyd = 0.D0
                    rho_sol = 0.D0
                    rhofluid1_est = 0.D0
                    rhofluid2_est = 0.D0
                    gl%solidtype(2) = 0    !No hydrate forms
                    iFlash = 1          !T given, melting point

                    iphase = 2  !First phase is vapor

                    T_3phase = (tmin + tmax) / 2    !Start value for the three phase temperature iteration
                    iFlash = 1 !p given, Melting point
                    call ptflash_solid_NC_3P(gl,press, T_3phase, x, rho_sol, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_2C3P, iFlash, iphase, iter, errval)

                    if (errval == 0) then
                        phasetype(1) = 1
                        if (phasetype_max(2) == 4) then
                            phasetype(2) = phasetype_max(1)
                        else
                            phasetype(2) = phasetype_max(2)
                        end if
                        phasetype(3) = 4
                        x_phase(:,phasetype(1)) = x_fluid1
                        x_phase(:,phasetype(2)) = x_fluid2
                        x_phase(:,phasetype(3)) = x_sol
                        rho(phasetype(1)) = rho_sol(1)
                        rho(phasetype(2)) = rho_sol(2)
                        rho(phasetype(3)) = rho_sol(3)
                        nrofphases = 3
                        !Compute the enthalpies of the 3 phases
                        do i = 1, 2 !First 2 phases are fluid phases
                            ! calculate the properties of the first phase
                            gl%molfractions = x_phase(:,phasetype(i))
                            call reduced_parameters_calc(gl,T_3phase) !Dummy temperature 300 K for the SRK
                            !Fluid phase
                            d_phase = rho(phasetype(i))
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(i) = h_calc(gl,T_3phase, d_phase, 0)
                            else                    !Entropy
                                z_phase(i) = s_calc(gl,T_3phase, d_phase, 0)
                            end if
                        end do
                        gl%tredmix = tredmix_orig
                        gl%rhoredmix = rhoredmix_orig
                        !Calculate the enthalpy of the solid phase
                        if (gl%solidtype(1) == 1) then
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(3) = h_WaterIce(gl,T_3phase, press)
                            else
                                z_phase(3) = s_WaterIce(gl,T_3phase, press)
                            end if
                        else
                            if (h_or_s == 1) then   !Enthalpy
                                z_phase(3) = h_DryIce(gl,T_3phase, press)
                            else
                                z_phase(3) = s_DryIce(gl,T_3phase, press)
                            end if
                        end if


                        !Calculate the enthalpies on the phase boundaries:
                        !VIw / VLw  (VIc / VLc )  ->  Phase 1 + Phase 2 <-> Phase 1 + Phase 3    beta_Phase3 = 0 : beta_Phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2) <->  beta_Phase2 = 0 : beta_phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3)

                        if (  ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 3)) .or. &            !VIw / VLw      or  (VIc / VLc )
                            &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 1) .and. (phasetype_max(2) == 2)) ) then            !VIc / VLc      or  (Backup)


                            beta_phase1 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        !VIw / LwIw (VIc / LcIc)  ->  Phase 1 + Phase 3 <-> Phase 2 + Phase 3   beta_Phase3 = 0 : beta_Phase1 = (z1 - x1_Phase2) / (x1_Phase1 - x1_Phase2) <->  beta_Phase2 = 0 : beta_phase1 = (z1 - x1_Phase3) / (x1_Phase1 - x1_Phase3)

                        if (  ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 3) .and. (phasetype_max(2) == 4)) .or. &            !VIw / LwIw     or  (VIc / LcIc)
                            &     ((phasetype_min(1) == 1) .and. (phasetype_min(2) == 4) .and. (phasetype_max(1) == 2) .and. (phasetype_max(2) == 4)) ) then            !VIc / LcIc         (Backup)

                            beta_phase2 = (x(1) - x_phase(1,phasetype(3))) / (x_phase(1,phasetype(2)) - x_phase(1,phasetype(3)))
                            zmin = beta_phase2 * z_phase(2) + (1.D0 - beta_phase2) * z_phase(3)

                            beta_phase1 = (x(1) - x_phase(1,phasetype(2))) / (x_phase(1,phasetype(1)) - x_phase(1,phasetype(2)))
                            zmax = beta_phase1 * z_phase(1) + (1.D0 - beta_phase1) * z_phase(2)

                        end if


                        if ((zmin <= z_spec) .and. (zmax >= z_spec)) then

                            Temp = T_3phase
                            !Solve the system of equations:
                            ! x1_Phase1 * beta_Phase1 + x1_Phase2 * beta_Phase2 + X1_Phase3 * beta_Phase3 = z1
                            ! beta_Phase1 + beta_Phase2 + beta_Phase3 = 1
                            ! beta_Phase1 * h_Phase1 + beta_Phase2 * h_Phase2 + beta_Phase3 * h_Phase3 = h_spec
                            ! for the phase fractions


                            Mat_phasefrac_2Comp(1,1) = x_phase(1,phasetype(1))
                            Mat_phasefrac_2Comp(2,1) = x_phase(1,phasetype(2))
                            Mat_phasefrac_2Comp(3,1) = x_phase(1,phasetype(3))

                            Mat_phasefrac_2Comp(1,2) = z_phase(1)
                            Mat_phasefrac_2Comp(2,2) = z_phase(2)
                            Mat_phasefrac_2Comp(3,2) = z_phase(3)

                            Mat_phasefrac_2Comp(1,3) = 1.D0
                            Mat_phasefrac_2Comp(2,3) = 1.D0
                            Mat_phasefrac_2Comp(3,3) = 1.D0

                            Vec_rightside_2Comp(1) = x(1)
                            Vec_rightside_2Comp(2) = z_spec
                            Vec_rightside_2Comp(3) = 1.D0

                            Vec_phasefrac_2Comp = 0.D0

                            !Solve the system of equations
                            eqn = 3
                            call LUdecomp(gl,eqn,Mat_phasefrac_2Comp,Vec_rightside_2Comp,errval,herr)
                            if (errval == 0) then
                                !phasefrac(1:3) = Vec_rightside_2Comp(1:3)
                                phasefrac(phasetype(1)) = Vec_rightside_2Comp(1)
                                phasefrac(phasetype(2)) = Vec_rightside_2Comp(2)
                                phasefrac(phasetype(3)) = Vec_rightside_2Comp(3)
                                gl%solidtype_akt_phase = gl%solidtype(1)
                                gl%solidpos_akt_phase = gl%solid_pos
                                return
                            end if
                        end if

                    end if

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                end if

            end if

        end if
    end if

    T_min_allowed = Tmin
    T_max_allowed = Tmax

    if (h_or_s == 1) then   !Enthalpy
        call Regula_Falsi(gl,Enth_dif_sol, Temp, Tmin, Tmax, Delta_allowed, &
            T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errval, parameters)
    else                    !Entropy
        call Regula_Falsi(gl,Entropy_dif_sol, Temp, Tmin, Tmax, Delta_allowed, &
            T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errval, parameters)
    end if

    !Calculation successful
    if (errval == 0) then
        call PhaseDet_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)

        if ((nrofphases == 3) .and. (errval == 0)) then       !Do some more steps of Regula falsi to get better results for three phase equilibria

            !The temperature should be accurate to Delta_allowed (1.D-6). Search in the interval Temp +- 1.D-5 for a better solution
            Tmin = Temp - Delta_allowed * 1.D1
            Tmax = Temp + Delta_allowed * 1.D1
            T_min_allowed = Temp - Delta_allowed * 1.D5
            T_max_allowed = Temp + Delta_allowed * 1.D5
            Delta_allowed = 1.D-8

            if (h_or_s == 1) then   !Enthalpy
                call Regula_Falsi(gl,Enth_dif_sol, Temp, Tmin, Tmax, Delta_allowed, &
                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errval, parameters)
            else                    !Entropy
                call Regula_Falsi(gl,Entropy_dif_sol, Temp, Tmin, Tmax, Delta_allowed, &
                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errval, parameters)
            End if

            if (errval == 0) then
                call PhaseDet_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)
            else
                errval = -2222
            end if

        end if

    else
        errval = -2222
    end if

    !New module variable to save the phase information and pass it to the interface routines
    gl%phase_id = phasetype

    end



    ! This function works as zerofunction for the regula falsi routine when overall enthalpy and pressure are given
    ! Andreas March 2014
    Double Precision Function Enth_dif_sol(gl,Temp, Parameters)







    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: press             ! Inputvariables
    Double precision :: h_spec, Temp
    type(type_additional_parameters) :: parameters

    !Variables for PhaseDet
    double precision, dimension(30):: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_Phase    !x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, dimension(5) :: Phasefrac
    double precision :: d_phase
    double precision, dimension(3) :: h_phase
    integer :: errval, nrofphases, i
    integer, dimension(5) :: phasetype

    !Variables needed for the calculation of hydrate enthalpy
    double precision, dimension(30) :: chem_pot, x_fluid, x_hyd, x_fluid1, x_fluid2     ! , fug_g
    double precision :: d_fluid
    double precision :: h_hyd, ChemPot_hyd
    integer :: pos_DryIce
    !Mixed hydrate 2015:
    double precision, dimension(30) :: fug_g    ! NO fug_CO2
    !  --------------------------------------------------

    h_spec = parameters%a_p(1)
    press = parameters%a_p(2)
    x = parameters%a_p(3:32)
    call PhaseDet_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)

    if (errval == 0) then
        do i = 1, nrofphases
            ! calculate the properties of the first phase
            if (gl%ncomp > 1) then
                gl%molfractions = x_phase(:,phasetype(i))
                call reduced_parameters_calc(gl,Temp) !Dummy temperature 300 K for the SRK
            end if
            if (phasetype(i) < 4) then  !Fluid phase
                d_phase = rho(phasetype(i))
                h_phase(i) = h_calc(gl,Temp, d_phase, 0)
            elseif (phasetype(i) == 4) then     !Solid (pure)
                if (gl%solidpos_akt_phase == 1) then         !Water
                    h_phase(i) = h_WaterIce(gl,Temp, press)
                elseif (gl%solidpos_akt_phase == 2) then     !CO2
                    h_phase(i) = h_DryIce(gl,Temp, press)
                else                                !No solid equation available
                    errval = -9904
                    Enth_dif_sol = -9904
                end if
            elseif (phasetype(i) == 5) then     !Hydrate
                !Calculate the enthalpy of the hydrate phase
                !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                !calculate the entropy and enthalpy of hydrates
                if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                    x_fluid = x_phase(:,phasetype(1))
                    d_fluid = rho(phasetype(1))
                    x_hyd = x_phase(:,phasetype(i))
                    !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                    !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                    !SH pure hdrt
                    !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                    !mix hdrt
                    call hdrt_ancillary_hs(gl, Temp, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                    if (errval == 0) then
                        !SH pure hdrt
                        !call hdrt_enthalpy(Temp, press*1.D6, x_hyd, chem_pot, fug_g,  h_hyd, errval)
                        !mix hdrt
                        call hdrt_enthalpy(gl, Temp, press*1.D6, x_hyd, chem_pot, fug_g,  h_hyd, errval)
                    end if

                    if (errval /= 0) then
                        return
                    end if
                    h_phase(i) = h_hyd

                else
                    if ((gl%solidpos_akt_phase == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
                        pos_DryIce = gl%solidtype_akt_phase
                        fug_g(1) = fug_DryIce(gl,Temp, press, pos_dryIce) *1.D6
                        if (fug_g(1) < 1.D-12) then
                            errval = -7779
                        else
                            !fug_CO2 = fug_g(1)
                            chem_pot(2) = g_DryIce(gl,Temp, press)
                            call hdrt_chem_potent_w(gl,Temp, press*1.d6, fug_g, ChemPot_hyd)
                            chem_pot(1) = ChemPot_hyd

                            x_hyd = x_phase(:,phasetype(i))

                            !SH pure hdrt
                            !call hdrt_enthalpy(Temp, press*1.D6, x_hyd, chem_pot, fug_g, h_hyd, errval)
                            !mix hdrt
                            call hdrt_enthalpy(gl, Temp, press*1.D6, x_hyd, chem_pot, fug_g, h_hyd, errval)
                            if (errval /= 0) then
                                return
                            end if

                        end if
                    else
                        errval = -12900
                    end if
                end if
            end if
        end do
    else
        parameters%a_p(65) = errval
    end if


    if (errval == 0) then
        if (nrofphases == 3) then
            Enth_dif_sol = h_spec - (h_phase(1) * Phasefrac(phasetype(1)) + h_phase(2) * Phasefrac(phasetype(2)) + h_phase(3) * Phasefrac(phasetype(3)))
        elseif(nrofphases == 2) then
            Enth_dif_sol = h_spec - (h_phase(1) * Phasefrac(phasetype(1)) + h_phase(2) * Phasefrac(phasetype(2)))
        elseif(nrofphases == 1) then
            Enth_dif_sol = h_spec - h_phase(1)
        end if
    else
        parameters%a_p(65) = errval
    end if

    End Function



    ! This function works as zerofunction for the regula falsi routine when overall entropy and pressure are given
    ! Andreas September 2014
    Double Precision Function Entropy_dif_sol(gl,Temp, Parameters)







    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: press             ! Inputvariables
    Double precision :: s_spec, Temp
    type(type_additional_parameters) :: parameters

    !Variables for PhaseDet
    double precision, dimension(30):: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_Phase    !x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, dimension(5) :: Phasefrac
    double precision :: d_phase
    double precision, dimension(3) :: s_phase
    integer :: errval, nrofphases, i
    integer, dimension(5) :: phasetype

    !Variables needed for the calculation of hydrate enthalpy
    double precision, dimension(30) :: chem_pot, x_fluid, x_hyd, x_fluid1, x_fluid2     !, fug_g
    double precision :: d_fluid
    double precision :: s_hyd, ChemPot_hyd
    integer :: pos_DryIce
    !Mixed hydrate 2015:
    double precision, dimension(30) :: fug_g    ! NO fug_CO2

    !  --------------------------------------------------

    s_spec = parameters%a_p(1)
    press = parameters%a_p(2)
    x = parameters%a_p(3:32)
    call PhaseDet_sol(gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, errval)

    if (errval == 0) then
        do i = 1, nrofphases
            ! calculate the properties of the first phase
            if (gl%ncomp > 1) then
                gl%molfractions = x_phase(:,phasetype(i))
                call reduced_parameters_calc(gl,Temp) !Dummy temperature 300 K for the SRK
            end if
            if (phasetype(i) < 4) then  !Fluid phase
                d_phase = rho(phasetype(i))
                s_phase(i) = s_calc(gl,Temp, d_phase, 0)
            elseif (phasetype(i) == 4) then     !Solid (pure)
                if (gl%solidpos_akt_phase == 1) then         !Water
                    s_phase(i) = s_WaterIce(gl,Temp, press)
                elseif (gl%solidpos_akt_phase == 2) then     !CO2
                    s_phase(i) = s_DryIce(gl,Temp, press)
                else                                !No solid equation available
                    errval = -9904
                    Entropy_dif_sol = -9904
                end if
            elseif (phasetype(i) == 5) then     !Hydrate
                !Calculate the enthalpy of the hydrate phase
                !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                !calculate the entropy and enthalpy of hydrates
                if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
                    x_fluid = x_phase(:,phasetype(1))
                    d_fluid = rho(phasetype(1))
                    x_hyd = x_phase(:,phasetype(i))
                    !Ancillary equation has to be called before enthalpy, entropy, or Gibbs energy of hydrates can be calculated because
                    !the chemical potential of all components in the hydrate phase as well as the fugacity of the guest is needed
                    !SH pure hdrt
                    !call hdrt_ancillary_hs(Tmin, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                    !mix hdrt
                    call hdrt_ancillary_hs(gl, Temp, d_fluid, press *1.D6, x_fluid, chem_pot, fug_g, errval)
                    if (errval == 0) then
                        !SH pure hdrt
                        !call hdrt_entropy(Temp, press*1.D6, x_hyd, fug_g, s_hyd, errval)
                        !mix hdrt
                        call hdrt_entropy(gl, Temp, press*1.D6, x_hyd, fug_g, s_hyd, errval)
                    end if

                    if (errval /= 0) then
                        return
                    end if
                    s_phase(i) = s_hyd

                else
                    if ((gl%solidpos_akt_phase == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
                        pos_DryIce = gl%solidtype_akt_phase
                        fug_g(1) = fug_DryIce(gl,Temp, press, pos_dryIce) *1.D6
                        if (fug_g(1) < 1.D-12) then
                            errval = -7779
                        else
                            !fug_CO2 = fug_g(1)
                            chem_pot(2) = g_DryIce(gl,Temp, press)
                            call hdrt_chem_potent_w(gl,Temp, press*1.d6, fug_g, ChemPot_hyd)
                            chem_pot(1) = ChemPot_hyd

                            x_hyd = x_phase(:,phasetype(i))

                            !SH pure hdrt
                            !call hdrt_entropy(Temp, press*1.D6, x_hyd, fug_g, s_hyd, errval)
                            !mix hdrt
                            call hdrt_entropy(gl, Temp, press*1.D6, x_hyd, fug_g, s_hyd, errval)

                            if (errval /= 0) then
                                return
                            end if
                        end if
                    else
                        errval = -12900
                    end if
                end if
            end if
        end do
    else
        parameters%a_p(65) = errval
    end if


    if (errval == 0) then
        if (nrofphases == 3) then
            Entropy_dif_sol = s_spec - (s_phase(1) * Phasefrac(phasetype(1)) + s_phase(2) * Phasefrac(phasetype(2)) + s_phase(3) * Phasefrac(phasetype(3)))
        elseif(nrofphases == 2) then
            Entropy_dif_sol = s_spec - (s_phase(1) * Phasefrac(phasetype(1)) + s_phase(2) * Phasefrac(phasetype(2)))
        elseif(nrofphases == 1) then
            Entropy_dif_sol = s_spec - s_phase(1)
        end if
    else
        parameters%a_p(65) = errval
    end if

    End Function


    !************************************************************************************
    subroutine Flash_PhaseBoundary_sol (gl,press, Temp, x, rho, x_Phase, phasetype, Phasefrac, nrofphases, iFlash, errval)
    !DEC$ ATTRIBUTES DLLEXPORT :: Flash_PhaseBoundary_sol
    !************************************************************************************
    !Subroutine for calculation of the temperature or pressure at which the first solid
    !phase occurs
    !The routine is basically a bisection method to find a temperature interval where the
    !first solid phase in the system forms, followed by a precise calculation of the
    !respective two or three phase equilibrium
    !Andreas, July 2014
    !--------------------------------------------------------------------------
    ! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS MAY BE PERFORMED:
    !   - SUBLIMATION / MELTING POINT: P GIVEN    --  iFlash = 1
    !   - SUBLIMATION / MELTING POINT: T GIVEN    --  iFlash = 2
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT: (depending on chosen iFlash)
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !   iFlash      - see above
    !
    ! OUTPUT: (depending on chosen iFlash)
    !   Temp        - Temperature [K]
    !   press       - Pressure [MPa]
    !
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   x_phase(5,30)   - matrix that contains all phase compositions
    !       1   x_vap       - Vapor phase composition
    !       2   x_liq1      - (lighter) liquid phase composition
    !       3   x_liq2      - (heavier) liquid phase composition
    !       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !       5   x_hyd       - hydrate phase compsition (so far only CO2 hydrates are possible)
    !   phasefrac           - Phasefraction of the phases
    !   nrofphases  - Number of phases present
    !   errval      - Error value
    !************************************************************************************






    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision, dimension(30):: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5):: x_Phase !x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    integer:: errval, iFlash, nrofphases

    integer:: i

    logical :: hdrt_check               !Variable indicates if hydrates may form (water + hydrate former in the mixture)
    logical :: dry_ice_check            !Variable indicates if dry may form (CO2 in the mixture)
    logical :: water_ice_check          !Variable indicates if solid water may form (water in the mixture)

    integer :: sol_pos_dryIce
    integer :: sol_pos_waterIce

    !Following variables are needed to store the temperature, phasetype, and nrofphases at the "low" and "high" test point as well as the "trial" test point
    !Inbetween these test points the point of solid formation should occur
    double precision :: T_high, T_low, T_trial, step_temp
    integer, dimension(5) :: phasetype_high, phasetype_low, phasetype_trial
    integer :: nrofphases_high, nrofphases_low, nrofphases_trial
    double precision, dimension(30) :: x_phase1_low, x_phase2_low, x_phase3_low, x_phase1_high, x_phase2_high, x_phase3_high
    double precision, dimension(5) :: phasefrac_low, phasefrac_high, phasefrac_calc

    !Variables needed for ptflash_solid_NC_2P
    double precision, dimension(30) :: x_hyd, x_fluid, x_solid
    double precision :: rhofluid_est, beta_loc
    integer :: iphase, iter
    double precision, dimension(3) :: rho_sol

    !Variables needed for ptflash_solid_NC_3P
    double precision, dimension(30) :: x_fluid1, x_fluid2
    double precision :: rhofluid1_est, rhofluid2_est, rhovap_est


    hdrt_check = .false.
    dry_ice_check = .false.
    water_ice_check = .false.

    sol_pos_dryIce = 0
    sol_pos_waterIce = 0

    T_high = 0.D0
    T_low = 0.D0
    x_phase1_low = 0.D0
    x_phase2_low = 0.D0
    x_phase3_low = 0.D0
    x_phase1_high = 0.D0
    x_phase2_high = 0.D0
    x_phase3_high = 0.D0
    phasefrac_low = 0.D0
    phasefrac_high = 0.d0

    rhofluid_est = 0.D0
    beta_loc = 0.d0
    iphase = 0
    rho_sol = 0.D0

    x_fluid1 = 0.d0
    x_fluid2 = 0.D0
    x_solid = 0.D0
    x_hyd = 0.D0
    rhovap_est = 0.D0
    rhofluid1_est = 0.D0
    rhofluid2_est = 0.D0

    if (iFlash == 1) then
        !Search for the temperature where the first solid (at the moment: dry ice, solid H2O, or hydrate) forms
        !Strategy:
        ! 1) Check, whether the solids can even be predicted for the mixture! At the moment: If no water or CO2 is present, quit with error.
        ! 2) Chose arbitrary elevated start temperature at which no solids should form under reasonable pressures and calculate pT-Flash with solids
        ! 3) Reduce temperature in 5 K steps until a solid phase is predicted, if no solid is predicted down to 120 K quit with error
        ! 4) Do a bisection until the temperature interval is small enough (arbitrarily chosen here)
        ! 5) CALCULATE THE PHASE BOUNDARY!
        !   Different scenarios can occur now (recall, in a mixture solid cannot be the only phase --> at low temp only more than one phase present):
        !   5a) high: 1 phase, low: 2 phase --> two phase boundary
        !   5b) high: 2 phase, low: 2 phase --> three phase boundary
        !   Only more than 2 components:
        !   5c) high: 1 phase, low: 3 phase --> three phase boundary    (Very unlikely case)
        !   5d) high: 2 phase, low: 3 phase --> three phase boundary
        !   5e) high: 3 phase, low: 3 phase --> four phase boundary    (Not implemented yet)
        !   --> IF ((nrofphases_high == 1) .and. (nrofphases_low == 2)) THEN two phase ELSE three phase

        !1)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        !At the moment only one hydrate former can be considered in the mixture. So only check hydrates of nrofhydrateformers == 2 (nrofhydrateformers == 1 will never happen, water is counted too)
        if ((gl%nrofhydrateformers > 0) .and. (gl%nrofhydrateformers < 3)) then
            hdrt_check = .true.
            water_ice_check = .true.
        end if

        if (.not. water_ice_check) then
            !Check if solid water has to be checked
            Do i = 1, gl%ncomp
                if (gl%Fluidlist_hydrate(i) == "water") then
                    sol_pos_waterIce = i
                    water_ice_check = .true.
                    exit
                end if
            end do
        end if

        !Check if solid co2 has to be checked
        Do i = 1, gl%ncomp
            if (gl%Fluidlist_hydrate(i) == "co2") then
                sol_pos_dryIce = i
                dry_ice_check = .true.
                exit
            end if
        end do

        !No model for solid formers available, quit with error
        If ((.not. hdrt_check) .and. (.not. water_ice_check) .and. (.not. dry_ice_check)) then
            errval = -12902
            return
        End if
        If ((.not. hdrt_check) .and. (water_ice_check) .and. (.not. dry_ice_check)) then    !If only solid water can form --> quit if pressure is above maximum pressure
            if (press > gl%pmax_Waterice) then
                errval = -12902
                return
            end if
        End if
        !1) END+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        !2)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        !Select start temperature and calculate pt-Flash with solids, store information at that point
        if (press < 10.D0) then
            T_high = 300.1D0  !K  Arbitrary starting temperature. No solid formation at that temperature
        else
            T_high = 340.1D0
        end if
        call PhaseDet_sol(gl,press, T_high, x, rho, x_Phase, phasetype_high, phasefrac_high, nrofphases_high, errval)
        If (errval /= 0) then
            return
        end if
        !Check if solid is predicted at this temperature, if yes, the temperature is too low (not very likely, error is more likely) --> quit with error
        if (maxval(phasetype_high) > 3) then
            errval = -15570
            return
        end if
        x_phase1_high = x_Phase(:,phasetype_high(1))
        if (nrofphases_high > 1) then
            x_phase2_high = x_Phase(:,phasetype_high(2))
        end if
        if (nrofphases_high > 2) then
            x_phase3_high = x_Phase(:,phasetype_high(3))
        end if
        !2) END+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        !3)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        !Search for the temperature where the first solid phase emerges. If no solid is predicted down to 120 K, quit with error
        T_low = T_high - 5.D0
        do while (T_low > 150.D0)
            call PhaseDet_sol(gl,press, T_low, x, rho, x_Phase, phasetype_low, phasefrac_low, nrofphases_low, errval)
            If (errval /= 0) then
                return
            end if
            !Check if solid is predicted at this temperature, if yes, leave this loop, if not, update upper boundary
            if (maxval(phasetype_low) > 3) then
                x_phase1_low = x_Phase(:,phasetype_low(1))
                if (nrofphases_low > 1) then
                    x_phase2_low = x_Phase(:,phasetype_low(2))
                end if
                if (nrofphases_low > 2) then
                    x_phase3_low = x_Phase(:,phasetype_low(3))
                end if
                exit
            else
                phasetype_high = phasetype_low
                x_phase1_high = x_Phase(:,phasetype_high(1))
                nrofphases_high = nrofphases_low
                if (nrofphases_high > 1) then
                    x_phase2_high = x_Phase(:,phasetype_high(2))
                end if
                if (nrofphases_high > 2) then
                    x_phase3_high = x_Phase(:,phasetype_high(3))
                end if
                T_high = T_low
                phasefrac_high = phasefrac_low
            end if
            !Stepsizes adjustable
            !Andreas Feb 2015
            !No hydrates: do always 5 K steps
            step_temp = 5.D0
            if (hdrt_check) then                                                                                ! SH: if very low water content decrease temperature faster
                if ((gl%T_Q_Hyd_2C(2) > 1.D-12) .and. ((T_low > 295.D0) .or. (T_low > (gl%T_Q_Hyd_2C(2) + 5.D0)) ) .or. (x(1) < 0.01d0 ) ) then
                    step_temp = 5.D0

                    !Check whether the given pressure is close to a quadruple point pressure and the temperature as well. If yes, do small steps. Andreas Feb 2015
                    !This could be improved in this future, priliminary version
                elseif ((press < gl%p_Q_Hyd_2C(2)) .and. (press > (gl%p_Q_Hyd_2C(2) - 0.5D0)) .and. (T_low < (gl%T_Q_Hyd_2C(2) + 1.5D0)) .and. (T_low > (gl%T_Q_Hyd_2C(2) - 1.5D0))) then !Upper quadruple point
                    step_temp = 0.05D0

                elseif ((press < gl%p_Q_Hyd_2C(1)) .and. (press > (gl%p_Q_Hyd_2C(1) - 0.5D0)) .and. (T_low < (gl%T_Q_Hyd_2C(1) + 1.5D0)) .and. (T_low > (gl%T_Q_Hyd_2C(1) - 1.5D0))) then !Upper quadruple point
                    step_temp = 0.05D0

                else
                    step_temp = 1.D0
                end if
            end if
            T_low = T_low - step_temp
            if (T_low < 150.D0) then
                errval = -15570
                return
            end if
            !If no dry ice in mixture, quit earlier (Arbitrarily set boundary)
            if (.not.dry_ice_check) then
                if (T_low < 190.D0) then
                    errval = -15570
                    return
                End if
            end if
        end do
        !3) END+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        !4)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        !Do bisection until the temperature interval is small enough (arbitrarily chosen: 4 steps at the moment, worst case: temperature is known within 2.5/2^4 K= 0,15 K)
        Do i = 1, 4
            T_trial = (T_high + T_low) / 2.D0
            call PhaseDet_sol(gl,press, T_trial, x, rho, x_Phase, phasetype_trial, Phasefrac, nrofphases_trial, errval)
            If (errval /= 0) then
                if (i > 1) then
                    !Andreas, September 2014
                    !If at least one bisection worked, proceed without improved start values
                    exit
                else
                    !Andreas, September 2014
                    !If the first bisection fails, quit with error
                    return
                end if
            end if
            !Check if solid is predicted at this temperature, if yes, replace low temperature, if not, replace high temperature
            if (maxval(phasetype_trial) > 3) then
                T_low = T_trial
                phasetype_low = phasetype_trial
                nrofphases_low = nrofphases_trial
                phasefrac_low = phasefrac
                x_phase1_low = x_Phase(:,phasetype_trial(1))
                if (nrofphases_trial > 1) then
                    x_phase2_low = x_Phase(:,phasetype_trial(2))
                end if
                if (nrofphases_trial > 2) then
                    x_phase3_low = x_Phase(:,phasetype_trial(3))
                end if
            else
                T_high = T_trial
                phasetype_high = phasetype_trial
                nrofphases_high = nrofphases_trial
                phasefrac_high = phasefrac
                x_phase1_high = x_Phase(:,phasetype_trial(1))
                if (nrofphases_trial > 1) then
                    x_phase2_high = x_Phase(:,phasetype_trial(2))
                end if
                if (nrofphases_trial > 2) then
                    x_phase3_high = x_Phase(:,phasetype_trial(3))
                end if
            end if
        End do
        !4) END+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        !5)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        !Calculate the phase boundary
        If ((nrofphases_high == 1) .and. (nrofphases_low == 2)) then
            !If one phase shows up twice --> Two phase boundary, if not --> three phase boundary

            Temp = (T_low + T_high) / 2.D0

            !TWO PHASE BOUNDARY (VH, VIw, VIc, LH, LIw, LIc possible)
            if (phasetype_low(1) == phasetype_high(1)) then !second phase at low temperature should always be solid
                !Check if fluid phase is vapor or liquid:
                If (phasetype_low(1) == 1) then
                    iphase = 2  !Vapor solution for density solver
                elseif((phasetype_low(1) == 2) .or. (phasetype_low(1) == 3)) then
                    iphase = 1  !Liquid solution for density solver
                else
                    !Quit with error
                    errval = -15570
                    return
                end if

                !Check which solid forms
                if (phasetype_low(2) == 4) then
                    !solid H2O or solid CO2
                    gl%solidtype(1) = gl%solidtype_akt_phase
                    gl%solidtype(2) = 0   !No hydrate forms
                    gl%solid_pos = gl%solidpos_akt_phase
                    x_solid = 0.D0
                    x_solid(gl%solid_pos) = 1.D0
                elseif (phasetype_low(2) == 5) then
                    gl%solidtype(1) = 0
                    gl%solidtype(2) = 1   !hydrate forms
                    gl%solid_pos = 1
                end if

                call ptflash_solid_NC_2P (gl,press, Temp, rho_sol, x, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, &
                    & iPhase, errval, iter)
                if (errval /= 0) then
                    return
                end if

                nrofphases = 2
                !First phase always fluid phase
                phasetype(1) = phasetype_low(1)
                rho(phasetype(1)) = rho_sol(1)
                !Second phase hydrate or solid
                phasetype(2) = phasetype_low(2)
                rho(phasetype(2)) = rho_sol(2)
                phasefrac(phasetype(1)) = 1.D0
                Phasefrac(phasetype(2)) = 0.D0
                x_Phase(:,phasetype(1)) = x_fluid
                if (phasetype(2) == 4) then
                    x_Phase(:,phasetype(2)) = x_solid
                elseif (phasetype(2) == 5) then
                    x_Phase(:,phasetype(2)) = x_hyd
                else
                    !Quit with error
                    errval = -15570
                    return
                end if

            else

                !THREE PHASE BOUNDARY. (Most likely LcHIc)

                !Initial estimates for phases

                !Phase 1, phase at high temperature will definitely be a fluid phase (vapor or light liquid)
                x_fluid1 = x_phase1_high
                phasetype(1) = phasetype_high(1)
                if (phasetype_high(1) == 1) then     !Phase(1) = Vapor
                    iphase = 2  !First phase is vapor
                elseif((phasetype_high(1) == 2) .or. (phasetype_high(1) == 3)) then    !Phase(1) = liquid
                    iphase = 1  !First phase is liquid
                else
                    !Quit with error
                    errval = -15570
                    return
                end if

                !Phase 2
                if (phasetype_low(1) < 4) then   !TWO FLUID PHASES
                    x_fluid2 = x_phase1_low
                    phasetype(2) = 3
                else
                    if (phasetype_low(1) == 4) then
                        !Initial estimates for solid phase
                        gl%solid_pos = gl%solidpos_akt_phase
                        gl%solidtype(1) = gl%solidtype_akt_phase
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.D0
                        phasetype(2) = 4
                    elseif (phasetype_low(1) == 5) then
                        gl%solidtype(2) = 1   !hydrate forms
                        phasetype(2) = 5
                    else
                        !Quit with error
                        errval = -15570
                        return
                    end if
                end if

                !Phase 3
                if (phasetype_low(2) == 4) then
                    !Initial estimates for solid phase
                    gl%solid_pos = gl%solidpos_akt_phase
                    gl%solidtype(1) = gl%solidtype_akt_phase
                    x_solid = 0.D0
                    x_solid(gl%solid_pos) = 1.D0
                    phasetype(3) = 4
                    !Initial estimates for hydrate phase
                elseif(phasetype_low(2) == 5) then
                    gl%solidtype(2) = 1   !hydrate forms
                    phasetype(3) = 5
                else
                    !Quit with error
                    errval = -15570
                    return
                end if

                !Initial estimate for phase fraction.
                phasefrac(1) = 0.5D0
                phasefrac(3) = 0.5D0

                !p given, melting point
                iFlash = 2

                call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_fluid1, x_fluid2, x_solid, x_hyd, rhofluid1_est, &
                    & rhofluid2_est, Phasefrac, iFlash, iphase, iter, errval)

                if (errval /= 0) then
                    return
                end if

                nrofphases = 3
                !First phase
                x_Phase(:,phasetype(1)) = x_fluid1
                rho(phasetype(1)) = rho_sol(1)

                !Second phase liquid
                if (phasetype(2) == 3) then
                    x_Phase(:,phasetype(2)) = x_fluid2
                elseif (phasetype(2) == 4) then
                    x_Phase(:,phasetype(2)) = x_solid
                elseif (phasetype(2) == 5) then
                    x_Phase(:,phasetype(2)) = x_hyd
                else
                    !Quit with error
                    errval = -15570
                    return
                end if
                rho(phasetype(2)) = rho_sol(2)

                !Third phase solid
                if (phasetype(3) == 4) then
                    x_Phase(:,phasetype(3)) = x_solid
                elseif (phasetype(3) == 5) then
                    x_Phase(:,phasetype(3)) = x_hyd
                else
                    !Quit with error
                    errval = -15570
                    return
                end if
                rho(phasetype(3)) = rho_sol(3)

            end if

        else
            !THREE PHASE BOUNDARY
            ! Possible: VLH, LLH, VHIw, VLIw, LLIw?, VLIc, LLIc?, LwHIc, LwHIw
            ! Rather other solids before these three phase lines appear: VHIc, LcHIc, VLcH

            !Initial temperature guess
            Temp = (T_high + T_low) / 2.D0
            iFlash = 2      !THREE PHASE LINE ("Freezing" Point)        : p, x Given (VLH and VSH: beta_H = 0; VLS: beta_I = 0)

            rhovap_est = 0.D0
            rhofluid1_est = 0.D0
            rhofluid2_est = 0.D0

            !If at the low temperature 3 phases are present --> already take these as initial guess and calculate flash with beta_sol = 0
            if (nrofphases_low == 3) then

                !Initial estimates for fluid phase (First phase can be treated here, since first phase is always a fluid phase)
                x_fluid1 = x_phase1_low
                if (phasetype_low(1) == 1) then     !Phase(1) = Vapor
                    iphase = 2  !First phase is vapor
                elseif((phasetype_low(1) == 2) .or. (phasetype_low(1) == 3)) then    !Phase(1) = liquid
                    iphase = 1  !First phase is liquid
                else
                    !Quit with error
                    errval = -15570
                    return
                end if

                !Initial estimate for phase fraction.
                !NOTE: On the freezing line, the phase fraction of solid will be 0. To simplify calculations, the SECOND phasefraction in ptflash_solid_NC_3P will always be the solid
                phasefrac(1)= phasefrac_low(1)
                phasefrac(3) = 1.D0 - phasefrac(1)

                if (phasetype_low(2) < 4) then   !TWO FLUID PHASES
                    x_fluid2 = x_phase2_low
                    if (phasetype_low(3) == 4) then     !Phase(3) = Ice
                        !Initial estimates for solid phase
                        gl%solid_pos = gl%solidpos_akt_phase
                        gl%solidtype(1) = gl%solidtype_akt_phase
                        x_solid = 0.D0
                        x_solid(gl%solid_pos) = 1.D0
                        !Initial estimates for hydrate phase
                        gl%solidtype(2) = 0   !No hydrate forms
                    elseif(phasetype_low(3) == 5) then    !Phase(3) = Hydrate
                        !Initial estimates for hydrate phase
                        gl%solidtype(2) = 1   !hydrate forms
                        !No solid forms
                        gl%solidtype(1) = 0
                    else
                        !Quit with error
                        errval = -15570
                        return
                    end if

                else                             !Solid + hydrate + fluid
                    !Initial estimates for solid phase
                    gl%solid_pos = gl%solidpos_akt_phase
                    gl%solidtype(1) = gl%solidtype_akt_phase
                    x_solid = 0.D0
                    x_solid(gl%solid_pos) = 1.D0
                    !Initial estimates for hydrate phase
                    gl%solidtype(2) = 1   !hydrate forms
                end if

                !p given, melting point
                iFlash = 2

                call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_fluid1, x_fluid2, x_solid, x_hyd, rhofluid1_est, &
                    & rhofluid2_est, Phasefrac, iFlash, iphase, iter, errval)

                if (errval /= 0) then
                    return
                end if

                nrofphases = 3
                !First phase
                phasetype(1) = phasetype_low(1)
                x_Phase(:,phasetype(1)) = x_fluid1
                rho(phasetype(1)) = rho_sol(1)

                !Second phase liquid
                phasetype(2) = 3
                x_Phase(:,phasetype(2)) = x_fluid2
                rho(phasetype(2)) = rho_sol(2)

                !Third phase solid
                if (gl%solidtype(2) == 1) then         !hydrate
                    phasetype(3) = 5
                    x_Phase(:,phasetype(3)) = x_hyd
                else                                !Ice
                    phasetype(3) = 4
                    x_Phase(:,phasetype(3)) = x_solid
                end if
                rho(phasetype(3)) = rho_sol(3)

            elseif (nrofphases_low == 2) then
                !If at the low temperature 2 phases are present, distuingish betweeen
                !   1) Two fluid phases + hydrate
                !   2) Two fluid phases + solid


                !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !two fluid phases + solid phase
                !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                if     (  ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 4) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 3)) .or. &            !VIw / VLw      or  (VIc / VLc )
                    &     ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 4) .and. (phasetype_high(1) == 3) .and. (phasetype_high(2) == 4)) .or. &            !VIw / LwIw     or  (VIc / LcIc)
                    &     ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 4) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 2)) .or. &            !VIc / VLc      or  (Backup)
                    &     ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 4) .and. (phasetype_high(1) == 2) .and. (phasetype_high(2) == 4)) ) then            !VIc / LcIc         (Backup)


                    !In each of the four cases above, phasetype_min(1) is vapor
                    x_fluid1 = x_phase1_low
                    !the liquid phase is either stored in phasetype_max(1) oder phasetype_max(2)
                    if (phasetype_high(2) == 4) then
                        x_fluid2 = x_phase1_high
                    else
                        x_fluid2 = x_phase2_high
                    end if

                    gl%solidtype(1) = gl%solidtype_akt_phase
                    gl%solidtype(2) = 0    !No hydrate forms
                    gl%solid_pos = gl%solidpos_akt_phase

                    x_solid = 0.D0
                    x_solid(gl%solid_pos) = 1.D0 !solid water + hydrate

                    iphase = 2  !First phase is vapor

                    !Initial estimate for phase fraction.
                    !NOTE: On the freezing line, the phase fraction of solid will be 0. To simplify calculations, the SECOND phasefraction in ptflash_solid_NC_3P will always be the solid
                    phasefrac_calc(1)= phasefrac_low(1)
                    phasefrac_calc(3) = 1.D0 - phasefrac_low(1)


                    !p given, melting point
                    iFlash = 2

                    call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_fluid1, x_fluid2, x_solid, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_calc, iFlash, iphase, iter, errval)

                    if (errval /= 0) then
                        return
                    end if

                    nrofphases = 3
                    !First phase always vapor
                    phasetype(1) = 1
                    !Second phase liquid
                    phasetype(2) = 3
                    !Third phase solid
                    phasetype(3) = 4
                    phasefrac(phasetype(1)) = phasefrac_calc(1)
                    Phasefrac(phasetype(2)) = phasefrac_calc(3)
                    Phasefrac(phasetype(3)) = 0.D0
                    x_Phase(:,phasetype(1)) = x_fluid1
                    x_Phase(:,phasetype(2)) = x_fluid2
                    x_Phase(:,phasetype(3)) = x_solid
                    rho(phasetype(1)) = rho_sol(1)
                    rho(phasetype(2)) = rho_sol(2)
                    rho(phasetype(3)) = rho_sol(3)

                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !Two fluid phases + hydrate
                    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                elseif  ( ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 3)) .or. &            !VH  / VLw      or
                    &     ((phasetype_low(1) == 3) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 3)) .or. &            !LwH / VLw      or
                    &     ((phasetype_low(1) == 1) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 2)) .or. &            !VH  / VLc      or (DOES THIS PHASECHANGE EXIST??)
                    &     ((phasetype_low(1) == 2) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 5)) .or. &            !VH  / LcH      or (Backup if Lc has ID 2)
                    &     ((phasetype_low(1) == 3) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 1) .and. (phasetype_high(2) == 5)) .or. &            !VH  / LcH      or
                    &     ((phasetype_low(1) == 3) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 2) .and. (phasetype_high(2) == 3)) .or. &            !LwH / LcLw     or
                    &     ((phasetype_low(1) == 2) .and. (phasetype_low(2) == 5) .and. (phasetype_high(1) == 2) .and. (phasetype_high(2) == 3)) ) then            !LcH / LcLw


                    !in each of the 7 cases above, except for VH/LcH , phasetype_low(1) is the phase that shows up twice --> take Phasetype_high + hydrate as initial guess
                    if (phasetype_high(2) /= 5) then
                        x_fluid1 = x_phase1_high
                        x_fluid2 = x_phase2_high
                    else
                        x_fluid1 = x_phase1_high
                        x_fluid2 = x_phase1_low
                    end if

                    gl%solidtype(1) = 0    !no other solid forms
                    gl%solidtype(2) = 1    !hydrate forms

                    if (phasetype_high(1) == 1) then
                        iphase = 2  !First phase is vapor
                    else
                        iphase = 1  !First phase is liquid
                    end if

                    !Initial estimate for phase fraction.
                    !NOTE: On the freezing line, the phase fraction of solid will be 0. To simplify calculations, the SECOND phasefraction in ptflash_solid_NC_3P will always be the solid
                    phasefrac_calc(1)= phasefrac_low(1)
                    phasefrac_calc(3) = 1.D0 - phasefrac_low(1)

                    !p given, melting point
                    iFlash = 2

                    call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_fluid1, x_fluid2, x_solid, x_hyd, rhofluid1_est, &
                        & rhofluid2_est, phasefrac_calc, iFlash, iphase, iter, errval)

                    nrofphases = 3
                    !First phase is liquid or vapor
                    if (iphase == 2) then
                        phasetype(1) = 1    !Vapor
                    else
                        phasetype(1) = 2    !Liquid
                    end if
                    rho(phasetype(1)) = rho_sol(1)
                    !Second phase liquid
                    phasetype(2) = 3
                    rho(phasetype(2)) = rho_sol(2)
                    !Third phase hydrate
                    phasetype(3) = 5
                    rho(phasetype(3)) = rho_sol(3)
                    phasefrac(phasetype(1)) = phasefrac_calc(1)
                    Phasefrac(phasetype(2)) = phasefrac_calc(3)
                    Phasefrac(phasetype(3)) = 0.D0
                    x_Phase(:,phasetype(1)) = x_fluid1
                    x_Phase(:,phasetype(2)) = x_fluid2
                    x_Phase(:,phasetype(3)) = x_hyd

                    if (errval /= 0) then
                        return
                    end if

                else
                    !Quit with error
                    errval = -15570
                    return
                end if

            else

                !Quit with error
                errval = -15570
                return

            end if

        end if
        !5) END+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    elseif (iFlash == 2) then
        !Given temperature and unknown pressure is not implemented yet
        errval = -12901
        return
    else
        !iFlash does not exist
        errval = -12901
        return
    end if

    end subroutine
    !************************************************************************************


    end module phasedet_sol_module
