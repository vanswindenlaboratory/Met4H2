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

    ! module for file phasedet_mix.f90
    submodule (phasedet_mix_module) impl
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use rhomix_pt_module
    use flash_module
    use vle_derivs_module
    use phasenv_vbased_module
    use reduced_parameters_calc_module

    contains





    !************************************************************************************
    module subroutine PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified mixture
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !           (page 136 - 138)
    !--------------------------------------------------------------------------
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
    !   vapfrac        - Molar vapor fraction
    !   nrofphases  - Number of phases present
    !   errval      - Error value
    !************************************************************************************
    !Andreas Aug. 2011

    !--------------------Reference-----------------
    ! J. G. Gernert, A. Jaeger, R. Span, Fluid Phase Equiliria, 375:209-218 (2014)









    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x
    double precision, dimension(5),intent(out):: rho
    double precision, dimension(30,5),intent(out):: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision, intent(out):: vapfrac
    integer, intent(out):: errval, nrofphases

    double precision:: tpd_vap, tpd_liq2, del_g, curvature
    double precision:: rhovap_est, rholiq_est, sum, rho_orig
    double precision, dimension(30):: x_vap_start, x_liq2_start, K_val_start
    double precision, dimension(30):: x_vap,  x_liq2
    integer:: n, iFlash, iter, i, j, IFound
    logical:: converged
    ! for phase envelope calculation:
    double precision:: T_spec, p_spec, T_return(6), rho_return(6,2), x_return(30,6), help, help_x(30), step
    double precision:: x_vapold(30), x_liqold(30), help1, help2

    integer:: points_found
    !Variables needed for the stability analysis
    double precision, dimension(30):: x_found1, x_found2, x_trial
    logical:: StablePhase1, StablePhase2
    double precision:: rho_trial, rho_phase, rhomix_dummy

    !safty check
    double precision:: gibbsvap, gibbsliq, gibbs2phase, gibbshom, rhotest !g_calc

    !Indicate which phases are present
    integer, dimension(5), intent(out):: phasetype         !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    !E.g.: 1 phase present: liquid
    !-- >  nrofphases = 1
    !-- >  phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!
    integer :: iphase   !Indicates which phase is expected in the succ_sub_tpd routine, Andreas February 2014


    iFlash = 5
    vapfrac = 0.5D0
    sum = 0.D0
    x_trial = 0.D0
    x_found1 = 0.D0
    x_found2 = 0.D0
    rho = 0.d0
    rho_orig = 0.d0
    curvature = 0.d0
    x_Phase = 0.D0
    phasetype = 0
    iphase = 0
    nrofphases = 0
    errval = 0
    K_val_start= 0.d0
    !---------------------------------------------------------------------------
    !1) Assume, that 2 phases are present, Gernert et al. (2014), Section 3.1.1
    !---------------------------------------------------------------------------
    !Generate starting values for the phase split
    call PTX_startvals_pT (gl,press, Temp, x, x_vap, x_liq2, vapfrac, errval)
    x_vap_start = x_vap
    x_liq2_start = x_liq2

    !Gernert et al. (2014), Eq. (11)
    Do i = 1, gl%ncomp
        K_val_start(i) =  x_vap_start(i) / x_liq2_start(i)
    End do
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !2a) If the startvalues are fine, do n steps of successive substitution
    !    and check whether the Gibbs-Energy of the split system is lower than the
    !    Gibbs-Energy of the initial system
    !---------------------------------------------------------------------------
    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)
    ! calculate the density of the original phase, needed for tangent plane distance
    rho_orig = rhomix_calc(gl,Temp, Press, 0.D0, 0,0)
    if (rho_orig < 1.d-14) then
        errval = -8888
        !write (*,*)'error ', errval, ' -- overall mixture density iteration failed in subroutine "PhaseDet"'
        return
    end if
    ! check if the original phase is vapor like or liquid like
    curvature = d2P_drho2 (gl,Temp, rho_orig)
    ! negative curvature means: vapor like phase, positive curvature: liquid like phase
    ! write the original density on the correct position in the return vector
    ! THIS WILL BE OVERWRITTEN IF THE SYSTEM IS FOUND UNSTABLE !!!
    if (curvature > 0.d0) then
        rho(3) = rho_orig
        phasetype(1) = 3
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(1) = 3
        x_phase(:,3) = x
    else
        rho(1) = rho_orig
        phasetype(1) = 1
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(1) = 1
        x_phase(:,1) = x
    end if

    n = 3
    converged = .false.
    rhovap_est = 0.D0
    rholiq_est = 0.D0
    call Succ_Sub(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash, errval, n, converged)
    !Error handling was missing here, A. Jäger, March 2016
    !Error handling now leads to errors for SRK + PR S. Hielscher, May 2016
    !Error for SRK and PR resolved, May 2016, A. Jäger
    if (errval /= 0) then
        if (maxval(gl%molfractions) > 0.9d0) then
            if ((Temp > gl%tc(maxloc(gl%molfractions,1))) .and. (press * gl%molfractions(minloc(gl%molfractions(1:gl%ncomp),1)) < gl%ptp(minloc(gl%molfractions(1:gl%ncomp),1)) )) then
                errval = 0
                nrofphases = 1
                vapfrac = 1.d0
                return
            endif
        else
            return
        endif
    end if

    !Andreas, May 2014
    !if ((vapfrac <= 0) .or. (vapfrac >= 1)) then
    if ((vapfrac < 1.D-14) .or. ((vapfrac - 1.D0) > 1.D-14)) then
        StablePhase1 = .false.
        StablePhase2 = .false.
        !Generate trial phases for the tangent plane distance

        !Heavy trial phase
        !Gernert et al. (2014), Eq. (23); from dew point (x" == x_spec) to heavy phase
        Do i = 1, gl%ncomp
            x_trial(i) = x(i) / K_val_start(i)
            sum = sum + x_trial(i)
        End do
        x_trial = x_trial / sum
        sum = 0.D0
        x_found1 = x_trial
        n = 20
        call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found1, tpd_liq2, iphase, StablePhase1, n, errval)
        !Error handling was missing here, A. Jäger, March 2016
        if (errval /= 0) then
            if (maxval(gl%molfractions) > 0.9d0) then
                if ((Temp > gl%tc(maxloc(gl%molfractions,1))) .and. (press * gl%molfractions(minloc(gl%molfractions(1:gl%ncomp),1)) < gl%ptp(minloc(gl%molfractions(1:gl%ncomp),1)) )) then
                    errval = 0
                    nrofphases = 1
                    vapfrac = 1.d0
                    return
                endif
            else
                return
            endif
        end if

        !Light trial phase
        Do i = 1, gl%ncomp
            x_trial(i) = x(i) * K_val_start(i)
            sum = sum + x_trial(i)
        End do
        x_trial = x_trial / sum
        sum = 0.D0
        x_found2 = x_trial
        n = 20
        call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found2, tpd_vap, iphase, StablePhase2, n, errval)
        !Error handling was missing here, A. Jäger, March 2016
        if (errval /= 0) then
            if (maxval(gl%molfractions) > 0.9d0) then
                if ((Temp > gl%tc(maxloc(gl%molfractions,1))) .and. (press * gl%molfractions(minloc(gl%molfractions(1:gl%ncomp),1)) < gl%ptp(minloc(gl%molfractions(1:gl%ncomp),1)) )) then
                    errval = 0
                    nrofphases = 1
                    vapfrac = 1.d0
                    return
                endif
            else
                return
            endif
        end if

        !Check whether an instability was found. If not exit, if yes choose the trial phase with the lowest tpd
        If(StablePhase1 .And. StablePhase2) then
            nrofphases = 1
            if (curvature > 0.d0) then
                x_liq2 = x
                x_vap = 0.d0
                !Andreas March 2013
                vapfrac = 0.D0
            else
                x_vap = x
                x_liq2 = 0.d0
                !Andreas March 2013
                vapfrac = 1.D0
            end if
        else
            nrofphases = 2
            if(tpd_liq2 < tpd_vap) then
                x_vap = x_found2!x
                x_liq2 = x_found1
            else
                x_liq2 = x_found1!x
                x_vap = x_found2
            end if
        end if
    else
        tpd_vap = tpd(gl,press, temp, x, rho_orig, x_vap, errval)        !Calculate the tangent plane distance with the vapor phase as trial phase
        tpd_liq2 = tpd(gl,press, temp, x, rho_orig, x_liq2, errval)      !Calculate the tangent plane distance with the liquid phase as trial phase
        !Gernert et al. (2014), Eq. (20)
        del_g = (1-vapfrac) * tpd_liq2 + vapfrac * tpd_vap        !Calculate the change of the Gibbs-Energy

        !If the Gibbs Energy decreased, then the system splits in two phases
        !Calculate the phase equilibrium
        If ((del_g < -1.D-12) .or. (tpd_vap < -1.D-12) .or. (tpd_liq2 < -1.D-12)) then
            nrofphases = 2
            !call ptflash(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash, 0, errval, iter)
        else
            StablePhase1 = .false.
            StablePhase2 = .false.
            !Generate trial phases for the tangent plane distance

            !Heavy trial phase
            Do i = 1, gl%ncomp
                x_trial(i) = x(i) / K_val_start(i)
                sum = sum + x_trial(i)
            End do
            x_trial = x_trial / sum
            sum = 0.D0
            x_found1 = x_trial
            n = 20
            call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found1, tpd_liq2, iphase, StablePhase1, n, errval)
            !Error handling was missing here, A. Jäger, March 2016
            if (errval /= 0) then
                return
            end if

            ! check if the phase found by succ_sub_tpd is equivalent to the original phase
            if (.NOT. Stablephase1 .and. maxval(abs(x - x_found1)) < 1.d-4) then
                gl%molfractions = x_found1
                call reduced_parameters_calc(gl,temp)
                rho_trial = rhomix_calc(gl,press, temp, 0.d0, 0, 0)
                if (abs(rho_trial - rho_orig) < 1.d0) Stablephase1 = .true.
                gl%molfractions = x
                call reduced_parameters_calc(gl,temp)
            end if

            !Light trial phase
            Do i = 1, gl%ncomp
                x_trial(i) = x(i) * K_val_start(i)
                sum = sum + x_trial(i)
            End do
            x_trial = x_trial / sum
            sum = 0.D0
            x_found2 = x_trial
            n = 20
            call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found2, tpd_vap, iphase, StablePhase2, n, errval)
            !Error handling was missing here, A. Jäger, March 2016
            if (errval /= 0) then
                return
            end if

            ! check if the phase found by succ_sub_tpd is equivalent to the original phase
            if (.NOT. Stablephase2 .and. maxval(abs(x - x_found1)) < 1.d-4) then
                gl%molfractions = x_found1
                call reduced_parameters_calc(gl,temp)
                rho_trial = rhomix_calc(gl,press, temp, 0.d0, 0, 0)
                if (abs(rho_trial - rho_orig) < 1.d0) Stablephase1 = .true.
                gl%molfractions = x
                call reduced_parameters_calc(gl,temp)
            end if

            !Check whether an instability was found. If not exit, if yes choose the trial phase with the lowest tpd
            If(StablePhase1 .And. StablePhase2) then
                nrofphases = 1
                if (curvature > 0.d0) then
                    x_liq2 = x
                    x_vap = 0.d0
                    !Andreas March 2013
                    vapfrac = 0.D0
                else
                    x_vap = x
                    x_liq2 = 0.d0
                    !Andreas March 2013
                    vapfrac = 1.D0
                end if
            else
                nrofphases = 2
                if(tpd_liq2 < tpd_vap) then
                    x_vap = x_found2!x
                    x_liq2 = x_found1
                else
                    x_liq2 = x_found1!x
                    x_vap = x_found2
                end if
            end if
        end if
    end if

    errval = 0

    if (nrofphases == 2) then
        ! continue with Newton-Raphson iteration of the phase compositions, because the stability algorithm found the original (homogeneous) phase to be unstable.
        call Flash_pT(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, errval, iter)
        !Andreas Jäger, December 2018
        !Back-up solution - if the first attempt to calculate the phase equilibrium fails, try again with the start values coming from the PTX_startvals_pT routine.
        if (errval /= 0) then
            errval = 0
            x_vap = x_vap_start
            x_liq2 = x_liq2_start
            !Check if one of the start values for the phases is equal to the overall composition (simplified by only checking the first component).
            !If yes, slightly change the phase composition by:
            !assuming a phase fraction very close to 0 (if the liquid phase is equal to the overall composition) or 1 (if the vapor phase is equal to the overall composition)
            if (dabs(x_vap(1) - x(1)) < 1.D-14) then
                !Vapor composition is very close to the overall composition. Assume beta_V = 0.999.
                !It is: beta_V = (z(i) - x(i)_L) / (x(i)_V - x(i)_L) --> x(i)_V = (z(i) - x(i)_L) / beta_V + x(i)_L
                !Change the vapor composition accordingly
                sum = 0.D0
                do i = 1, gl%ncomp
                    x_vap(i) = (x(i) - x_liq2(i)) / 0.999D0 + x_liq2(i)
                    sum = sum + x_vap(i)
                end do
                !Normalize x_vap
                x_vap = x_vap / sum
            elseif (dabs(x_liq2(1) - x(1)) < 1.D-14) then
                !Liquid composition is very close to the overall composition. Assume beta_L = 0.999.
                !It is: beta_L = (z(i) - x(i)_V) / (x(i)_L - x(i)_V) --> x(i)_L = (z(i) - x(i)_V) / beta_L + x(i)_V
                !Change the liquid composition accordingly
                sum = 0.D0
                do i = 1, gl%ncomp
                    x_liq2(i) = (x(i) - x_vap(i)) / 0.999D0 + x_vap(i)
                    sum = sum + x_liq2(i)
                end do
                !Normalize x_vap
                x_liq2 = x_liq2 / sum
            end if
            if (errval == 0) then
                call Flash_pT(gl, press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, errval, iter)
            end if

            if (errval /= 0) then
                !Second backup method before trying to calculate the whole phase envelope (which is very slow and should be avoided, if possible)
                !Directly try to calculate the bubble point temperature and the dew point temperature at the given composition. Two cases may occur:
                !
                !1) Either the bubble or the dew point calculation failed: Proceed with phase envelope calculation
                !
                !
                !2) Both calculations succeeded: The phase envelope calculation is no longer needed. Three option can occur now:
                !
                !   a) If the specified temperature lies between the bubble point and the dew point temperature, the mixture is in the two phase region.
                !      Try once again pT-flash with the compositions on the phase boundary. If this fails, quit with error
                !   b) The temperature is below the bubble point temperature: The phase is homogeneous liquid
                !   c) The temperature is above the dew point temperature: The phase is homogeneous vapor

                !Try to calculate the dew point temperature
                errval = 0
                x_vap = 0.D0
                x_liq2 = 0.D0
                T_return = 0.D0
                iFlash = 2      !Dew point at given p
                x_vap = x
                call Flash_PhaseBoundary(gl,press, T_return(1), x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, iFlash, 0, errval, iter)
                if (errval == 0) then
                    !Calculation of dew point temperature successful.
                    !Save composition of the liquid phase
                    x_found2 = x_liq2
                    !Try to calculate the bubble point temperature
                    errval = 0
                    x_vap = 0.D0
                    x_liq2 = 0.D0
                    iFlash = 1      !Bubble point at given p
                    x_liq2 = x
                    call Flash_PhaseBoundary(gl,press, T_return(2), x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, iFlash, 0, errval, iter)
                    if (errval == 0) then
                        !Calculation of bubble point temperature successful.
                        !Save composition of the vapor phase
                        x_found1 = x_vap
                        !Check whether the given temperature is below the bubble point temperature, between the dew and bubble point temperature, or above the dew point temperature
                        if (Temp < T_return(2)) then
                            !Homogeneous liquid
                            nrofphases = 1
                        elseif (Temp > T_return(1)) then
                            !Homogeneous vapor
                            nrofphases = 1
                        else
                            !Temperature is in between the dew point temperature and the bubble point temperature. Try a flash calculation with the phase compositions at the dew and bubble point as initial values
                            errval = 0
                            x_vap = x_found1
                            x_liq2 = x_found2
                            call Flash_pT(gl, press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, errval, iter)
                        end if
                    end if
                end if
                iflash = 5  !Set back iFlash
                T_return = 0.D0 !Set back T_return
            end if

        end if
        if ((errval == 0) .and. (nrofphases == 2)) then
            gl%molfractions = x_vap
            call reduced_parameters_calc(gl,Temp)
            rhotest = gl%rho_vap
            gibbsvap = G_CALC(gl,Temp,rhotest, 0)
            gl%molfractions = x_liq2
            call reduced_parameters_calc(gl,Temp)
            rhotest = gl%rho_liq
            if((gl%seawater) .or. gl%el_present ) gl%gecarrier = .true.
            gibbsliq = G_CALC(gl,Temp,rhotest, 0)
            gl%gecarrier = .false.
            !wrong        gibbs2phase = gibbsvap + vapfrac * (gibbsliq - gibbsvap)
            gibbs2phase = gibbsliq + vapfrac * (gibbsvap - gibbsliq)
            gl%molfractions = x
            call reduced_parameters_calc(gl,Temp)
            rhotest = rho_orig
            gibbshom = G_CALC(gl,Temp,rhotest, 0)
            if((gibbshom - gibbs2phase) < 0.d0) nrofphases = 1
        end if
    end if

    ! if the Michelsen-Method fails, try the (very slow!) phase envelope calculation
    if (errval /= 0) then
        p_spec = press
        T_spec = 0.d0
        ! find all points on the phase envelope at the given pressure
        ! modified 08.2012 by J.Gernert: now uses the new ptdiag routine that tries to construct the whole phase envelope
        call ptdiag(gl,x, p_spec, T_spec, t_return, rho_return, x_return, points_found, errval)
        if (points_found < 1) then
            if (errval == 0) errval = -2222
            return !phasenv failed
        end if
        ! sort the points with respect to temperature
        do i = 1, points_found-1
            if ((T_return(i) > T_return(i+1)) .and. (T_return(i+1) > 0.d0)) then
                help = T_return(i)
                T_return(i) = T_return(i+1)
                T_return(i+1) = help
                help1 = rho_return(i,1)
                rho_return(i,1) = rho_return(i+1,1)
                rho_return(i+1,1) = help1
                help2 = rho_return(i,2)
                rho_return(i,2) = rho_return(i+1,2)
                rho_return(i+1,2) = help2
                help_x = x_return(:,i)
                x_return(:,i) = x_return(:,i+1)
                x_return(:,i+1) = help_x
            end if
        end do
        if ((T_return(1) >  T_return(2)) .and. (T_return(2) > 0.d0)) then
            help = T_return(1)
            T_return(1) = T_return(2)
            T_return(2) = help
            help1 = rho_return(1,1)
            rho_return(1,1) = rho_return(2,1)
            rho_return(2,1) = help1
            help2 = rho_return(i,2)
            rho_return(1,1) = rho_return(2,2)
            rho_return(2,2) = help2
            help_x = x_return(:,1)
            x_return(:,1) = x_return(:,2)
            x_return(:,2) = help_x
        end if

        ! check whether the system is in the two-phase region
        select case (points_found)
            !--------------------------------------------------------------------------------
        case(3) !unusual (open) phase envelope. the two phase regions are BELOW point 1 and BETWEEN point 2 and 3
            !--------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------
            if ((Temp < T_return(1)) .OR. ((Temp > T_return(2)) .AND. (Temp < T_return(3)))) then ! two-phase region
                !--------------------------------------------------------------------------------
                nrofphases = 2
                ! try calculating the pt-flash using the values on the phase boundary as initials
                if (Temp < T_return(1)) then
                    help = T_return(1)
                    help_x = x_return(:,1)
                else
                    help = T_return(3)
                    help_x = x_return(:,3)
                end if
                !--------------------------------------------------------------------------------
            else ! single phase
                !--------------------------------------------------------------------------------
                nrofphases = 1
                return
            end if
            !--------------------------------------------------------------------------------
        case (2) ! closed phase envelope (typical for natural gases)
            !--------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------
            if ((T_return(1) < Temp) .AND. (Temp < T_return(2))) then  ! two-phase region
                !--------------------------------------------------------------------------------
                help = T_return(2)
                help_x = x_return(:,2)
                !--------------------------------------------------------------------------------
            else !single phase
                !--------------------------------------------------------------------------------
                nrofphases = 1
                return
            end if
            !--------------------------------------------------------------------------------
        case (1) ! open phase envelope, one point at the designated pressure
            !--------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------
            if (Temp < T_return(1)) then  ! two-phase region
                !--------------------------------------------------------------------------------
                help = T_return(1)
                help_x = x_return(:,1)
                !--------------------------------------------------------------------------------
            else !single phase
                !--------------------------------------------------------------------------------
                nrofphases = 1
                return
            end if
        end select

        ! calculate the pt-flash, using the values on the phase boundary as initials
        x_liq2 = help_x
        x_vap = x
        errval = 0
        ! make 3 steps with successive substitution first
        call Succ_Sub(gl,press, temp, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, 5, errval, 3, converged)
        ! continue with Newton, using the results from succ. subst. as initials
        if (errval == 0) then
            call Flash_pT(gl,press, temp, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, errval, iter)
        end if
        if (errval /= 0) then
            ! if the flash calculation fails, try to go step by step from the phase boundary into the 2-phase region
            ! the compositions from the previous step are used as initials for the next step
            step = (help - Temp)/5.d0
            x_liq2 = help_x
            x_vap = x
            do while (help > (Temp + step))
                help = help - step
                x_liqold = x_liq2
                x_vapold = x_vap
                ! make 3 steps with successive substitution first
                call Succ_Sub(gl,press, help, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, 5, errval, 3, converged)
                ! continue with Newton, using the results from succ. subst. as initials
                call Flash_pT(gl,press, help, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, errval, iter)
                if (errval /= 0) then
                    help = help + step
                    step = step / 2.d0
                    x_vap = x_vapold
                    x_liq2 = x_liqold
                    if (step < 1.d-2) return    ! return if the step size becomes too small - >  phaseDet failed!
                end if
            end do
            ! now we have initial values close to the point we are looking for
            ! make 3 steps with successive substitution first
            call Succ_Sub(gl,press, temp, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, 5, errval, 3, converged)
            ! continue with Newton, using the results from succ. subst. as initials
            call Flash_pT(gl,press, temp, x, x_vap, x_liq2, 0.d0, 0.d0, vapfrac, errval, iter)
            !Andreas September 2012: error handle is needed here!!!
            if (errval /= 0) return
        end if
    end if

    if (nrofphases > 1) then
        !Determine which phases are present
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = gl%rho_vap
        If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .or. (gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
            rhomix_dummy = 0.d0
            do i = 1,gl%ncomp
                rhomix_dummy = gl%molfractions(i) * 1.d0 / gl%rhored (i) + rhomix_dummy
            enddo
            rhomix_dummy = 1.d0 / rhomix_dummy
        else
            rhomix_dummy = gl%rhoredmix
        endif
        if (rho_phase <= gl%factor_decide4vap * rhomix_dummy) then
            curvature = -1.d0
        else
            curvature = d2P_drho2 (gl,Temp, rho_phase)
        endif
        !do i = 1,30
        !    do j = 1,bad_nrmixcomp
        !        if (trim(components(i)) == trim(bad_mixcomp(j))) then
        !            bad_mixmodel = .true.
        !            call rho_test(Temp, press, rho_phase, 2, IFound,1)
        !            if (IFound == 1) then
        !                curvature = -1.d0
        !                exit
        !            endif
        !        endif
        !    enddo
        !enddo
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        if (curvature > 0.d0) then      !Light liquid phase
            rho(2) = rho_phase
            phasetype(1) = 2
            !New module variable to save the phase information and pass it to the interface routines
            !Andreas April 2013
            gl%phase_id(1) = 2
            x_phase(:,2) = x_vap
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            !New module variable to save the phase information and pass it to the interface routines
            !Andreas April 2013
            gl%phase_id(1) = 1
            x_phase(:,1) = x_vap
        end if
        gl%molfractions = x_liq2
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = gl%rho_liq
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        !Some cases might occur, where a (correct) liquid phase solution has a negative curvature and therefore
        !would be indicated as vapor. If this happens, the correct vapor solution would be overwritten. To prevent this the solution
        !will ALWAYS be stored at the third position as heavy liquid.
        rho(3) = rho_phase
        phasetype(2) = 3
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(2) = 3
        x_phase(:,3) = x_liq2
        ! make sure that the liquid with the higher density is always in position 3
        if (rho(3) < rho(2)) then
            !help = rho(3)
            !help_x = x_phase(:,3)
            rho(2:3) = rho(3:2:-1)
            x_phase(:,2:3) = x_phase(:,3:2:-1)
            !rho(2) = help
            !x_phase(:,2) = help_x
            vapfrac = 1.d0 - vapfrac  !Monika, Jan. 2017
        end if
    End if

    gl%molfractions = x
    gl%nphases = nrofphases
    call reduced_parameters_calc(gl,Temp)

    end subroutine PhaseDet

    !************************************************************************************
    module subroutine PhaseDet2(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    !************************************************************************************
    ! THIS ROUTINE IS STILL NEEDED BY THE SUBROUTINE phasenv IN ORDER TO GENERATE START VALUES!!!
    ! Subroutine for determining how many phases are present for the specified mixture
    ! THE ALGORITHM IS BASED ON THE FOLLOWING PUBLICATION:
    !--------------------------------------------------------------------------
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !           (page 136 - 138)
    !--------------------------------------------------------------------------
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
    !   vapfrac        - Molar vapor fraction
    !   errval      - Error value
    !************************************************************************************
    !Andreas Aug. 2011






    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x
    double precision, dimension(5):: rho
    double precision, dimension(30,5):: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision:: vapfrac
    integer:: errval, nrofphases

    double precision:: tpd_vap, tpd_liq2, del_g, curvature
    double precision:: rhovap_est, rholiq_est, sum, rho_orig
    double precision, dimension(30):: x_vap_start, x_liq2_start, K_val_start
    double precision, dimension(30):: x_vap,  x_liq2, x_sol
    integer:: n, iFlash, iter, i
    logical:: converged

    !Variables needed for the stability analysis
    double precision, dimension(30):: x_found1, x_found2, x_trial
    logical:: StablePhase1, StablePhase2
    double precision:: rho_phase

    !Indicate which phases are present
    integer, dimension(5) :: phasetype         !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    !E.g.: 1 phase present: liquid
    !-- >  nrofphases = 1
    !-- >  phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!
    integer :: iphase   !Indicates which phase is expected in the succ_sub_tpd routine, Andreas February 2014


    iFlash = 5
    vapfrac = 0.5D0
    sum = 0.D0
    x_trial = 0.D0
    x_found1 = 0.D0
    x_found2 = 0.D0
    rho = 0.d0
    rho_orig = 0.d0
    curvature = 0.d0
    x_Phase = 0.D0
    x_sol = 0.D0
    phasetype = 0
    nrofphases = 1
    errval = 0
    iphase = 0
    !---------------------------------------------------------------------------
    !1) Assume, that 2 phases are present
    !---------------------------------------------------------------------------
    !Generate starting values for the phase split
    call PTX_startvals_pT (gl,press, Temp, x, x_vap, x_liq2, vapfrac, errval)
    x_vap_start = x_vap
    x_liq2_start = x_liq2
    Do i = 1, gl%ncomp
        K_val_start(i) =  x_vap_start(i) / x_liq2_start(i)
    End do
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !2a) If the startvalues are fine, do n steps of successive substitution
    !    and check whether the Gibbs-Energy of the split system is lower than the
    !    Gibbs-Energy of the initial system
    !---------------------------------------------------------------------------
    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)
    ! calculate the density of the original phase
    rho_orig = rhomix_calc(gl,Temp, Press, 0.D0, 0,0)
    if (rho_orig == 0.d0) then
        errval = -8888
        !write (*,*)'error ', errval, ' -- overall mixture density iteration failed in subroutine "PhaseDet"'
        return
    end if
    ! check if the original phase is vapor like or liquid like
    curvature = d2P_drho2 (gl,Temp, rho_orig)
    ! negative curvature means: vapor like phase, positive curvature: liquid like phase
    ! write the original density on the correct position in the return vector
    ! THIS WILL BE OVERWRITTEN IF THE SYSTEM IS FOUND UNSTABLE !!!
    if (curvature > 0.d0) then
        rho(3) = rho_orig
        phasetype(1) = 3
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(1) = 3
        x_phase(:,3) = x
    else
        rho(1) = rho_orig
        phasetype(1) = 1
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(1) = 1
        x_phase(:,1) = x
    end if

    n = 3
    converged = .false.
    rhovap_est = 0.D0
    rholiq_est = 0.D0
    call Succ_Sub(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash, errval, n, converged)

    !Andreas, May 2014
    !if ((vapfrac <= 0) .or. (vapfrac >= 1)) then
    if ((vapfrac < 1.D-14) .or. ((vapfrac - 1.D0) > 1.D-14)) then
        StablePhase1 = .false.
        StablePhase2 = .false.
        !Generate trial phases for the tangent plane distance

        !Heavy trial phase
        Do i = 1, gl%ncomp
            x_trial(i) = x(i) / K_val_start(i)
            sum = sum + x_trial(i)
        End do
        x_trial = x_trial / sum
        sum = 0.D0
        x_found1 = x_trial
        n = 20
        call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found1, tpd_liq2, iphase, StablePhase1, n, errval)

        !Light trial phase
        Do i = 1, gl%ncomp
            x_trial(i) = x(i) * K_val_start(i)
            sum = sum + x_trial(i)
        End do
        x_trial = x_trial / sum
        sum = 0.D0
        x_found2 = x_trial
        n = 20
        call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found2, tpd_vap, iphase, StablePhase2, n, errval)


        !Check whether an instability was found. If not exit, if yes choose the trial phase with the lowest tpd
        If(StablePhase1 .And. StablePhase2) then
            nrofphases = 1
            if (curvature > 0.d0) then
                x_liq2 = x
                x_vap = 0.d0
            else
                x_vap = x
                x_liq2 = 0.d0
            end if
        else
            nrofphases = 2
            if(tpd_liq2 < tpd_vap) then
                x_vap = x_found2!x
                x_liq2 = x_found1
            else
                x_liq2 = x_found1!x
                x_vap = x_found2
            end if
        end if
    else
        tpd_vap = tpd(gl,press, temp, x, rho_orig, x_vap, errval)        !Calculate the tangent plane distance with the vapor phase as trial phase
        tpd_liq2 = tpd(gl,press, temp, x, rho_orig, x_liq2, errval)      !Calculate the tangent plane distance with the liquid phase as trial phase
        del_g = (1-vapfrac) * tpd_liq2 + vapfrac * tpd_vap        !Calculate the change of the Gibbs-Energy

        !If the Gibbs Energy decreased, then the system splits in two phases
        !Calculate the phase equilibrium
        If ((del_g < -1.D-12) .or. (tpd_vap < -1.D-12) .or. (tpd_liq2 < -1.D-12)) then
            nrofphases = 2
            !call ptflash(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash, 0, errval, iter)
        else
            StablePhase1 = .false.
            StablePhase2 = .false.
            !Generate trial phases for the tangent plane distance

            !Heavy trial phase
            Do i = 1, gl%ncomp
                x_trial(i) = x(i) / K_val_start(i)
                sum = sum + x_trial(i)
            End do
            x_trial = x_trial / sum
            sum = 0.D0
            x_found1 = x_trial
            n = 20
            call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found1, tpd_liq2, iphase, StablePhase1, n, errval)

            !Light trial phase
            Do i = 1, gl%ncomp
                x_trial(i) = x(i) * K_val_start(i)
                sum = sum + x_trial(i)
            End do
            x_trial = x_trial / sum
            sum = 0.D0
            x_found2 = x_trial
            n = 20
            call Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_found2, tpd_vap, iphase, StablePhase2, n, errval)

            !Check whether an instability was found. If not exit, if yes choose the trial phase with the lowest tpd
            If(StablePhase1 .And. StablePhase2) then
                nrofphases = 1
                if (curvature > 0.d0) then
                    x_liq2 = x
                    x_vap = 0.d0
                else
                    x_vap = x
                    x_liq2 = 0.d0
                end if
            else
                nrofphases = 2
                if(tpd_liq2 < tpd_vap) then
                    x_vap = x_found2!x
                    x_liq2 = x_found1
                else
                    x_liq2 = x_found1!x
                    x_vap = x_found2
                end if
            end if
        end if
    end if

    errval = 0
    if (nrofphases == 2) then
        ! do 3 steps of successive substitution
        !call Succ_Sub(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 3, converged)
        ! continue with Newton
        call Flash_pT(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, errval, iter)
        !Andreas September 2012: error handle is needed here!!!
        if (errval /= 0) return
    end if

    if (nrofphases > 1) then
        !Determine which phases are present
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = gl%rho_vap
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        if (curvature > 0.d0) then      !Light liquid phase
            rho(2) = rho_phase
            phasetype(1) = 2
            !New module variable to save the phase information and pass it to the interface routines
            !Andreas April 2013
            gl%phase_id(1) = 2
            x_phase(:,2) = x_vap
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            !New module variable to save the phase information and pass it to the interface routines
            !Andreas April 2013
            gl%phase_id(1) = 1
            x_phase(:,1) = x_vap
        end if
        gl%molfractions = x_liq2
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = gl%rho_liq
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        !Some cases might occur, where a (correct) liquid phase solution has a negative curvature and therefore
        !would be indicated as vapor. If this happens, the correct vapor solution would be overwritten. To prevent this the solution
        !will ALWAYS be stored at the third position as heavy liquid.
        rho(3) = rho_phase
        phasetype(2) = 3
        !New module variable to save the phase information and pass it to the interface routines
        !Andreas April 2013
        gl%phase_id(2) = 3
        x_phase(:,3) = x_liq2
    End if

    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)


    end subroutine PhaseDet2


    !************************************************************************************
    double precision module function tpd(gl,press, temp, x, rho_orig, x_trial, errval)
    !************************************************************************************
    ! Function for the calculation of the tangent plane distance of a trial phase x_trial
    !--------------------------------------------------------------------------
    !           Kunz, O. et al.
    !           The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !           GERG TM15, 2007
    !           (page 136 - 138)
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !   rho_orig    - density of the orignal phase
    !   x_trial     - Trial phase composition
    !************************************************************************************
    !Andreas Aug. 2011





    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, rho_orig
    double precision, dimension (30):: x, x_trial
    integer :: errval

    double precision, dimension (30) :: lnf_trial, lnf
    double precision:: rho_trial, rho

    integer:: i, iPhase

    tpd = 0.D0
    rho_trial = 0.D0
    rho = 0.D0
    iPhase = 0
    errval = 0

    ! Get the density of the trial phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    gl%molfractions = x_trial
    call reduced_parameters_calc(gl,Temp)
    ! calculate the liquid phase 2 density
    rho_trial = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)

    if (rho_trial == 0.D0)then
        errval = -8888
        !write (*,*) 'error ', errval, ' -- trial phase mixture density iteration failed in function "tpd"'
        return
    end if

    call lnf_mix(gl,Temp, rho_trial, press, lnf_trial)

    ! Get the density of the overall phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)
    ! calculate the liquid phase 2 density
    if (rho_orig > 0.d0) then
        rho = rho_orig
    else
        rho = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)
    end if

    if (rho == 0.D0)then
        errval = -8888
        !write (*,*) 'error ', errval, ' -- overall mixture density iteration failed in function "tpd"'
        return
    end if

    call lnf_mix(gl,Temp, rho, press, lnf)

    !Gernert et al. (2014), Eq. (21 or 22)
    Do i = 1, gl%Ncomp
        tpd = tpd + x_trial(i) * (lnf_trial(i) - lnf(i))
    end do

    End function tpd


    !************************************************************************************
    module subroutine d_tpd_dxj(gl,press, temp, x, x_trial, dtpddxj, errval)
    !************************************************************************************
    ! Subroutine for the calculation of the first derivative of the tangent plane distance of a trial phase x_trial
    ! with respect to each molfraction xj. The "last" molfraction xN is replaced by xN = 1 - sum(xi)
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !   x_trial     - Trial phase composition
    !
    ! OUTPUT:
    !   dtpddxj     - Vector containing the derivatives
    !************************************************************************************
    !Andreas Aug. 2011





    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x, x_trial
    !double precision, dimension(30), intent(inout):: x_trial
    integer:: errval
    double precision, dimension(30):: dtpddxj

    double precision, dimension(30):: lnf_trial, lnf
    double precision, dimension(30,30):: dlnfi_dxj_TP_trial
    double precision:: rho_trial, rho

    integer:: i, j, iPhase

    rho_trial = 0.D0
    rho = 0.D0
    dtpddxj = 0.D0
    iPhase = 0


    ! Get the density of the trial phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    gl%molfractions = x_trial
    call reduced_parameters_calc(gl,Temp)
    ! calculate the liquid phase 2 density
    rho_trial = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)

    if (rho_trial == 0.D0)then
        errval = -8888
        !write (*,*) 'error ', errval, ' -- trial phase mixture density iteration failed in function "tpd"'
        return
    end if

    call lnf_mix(gl,Temp, rho_trial, press, lnf_trial)
    call dlnfi_dxj_TP (gl,Temp, rho_trial, dlnfi_dxj_TP_trial)

    ! Get the density of the overall phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)
    ! calculate the liquid phase 2 density
    rho = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)

    if (rho == 0.D0)then
        errval = -8888
        !write (*,*) 'error ', errval, ' -- overall mixture density iteration failed in function "tpd"'
        return
    end if

    call lnf_mix(gl,Temp, rho, press, lnf)

    Do j = 1, gl%Ncomp - 1
        dtpddxj = lnf_trial(j) - lnf_trial(gl%Ncomp)+ lnf(gl%ncomp) - lnf(j)
        Do i = 1, gl%Ncomp
            dtpddxj(j) = dtpddxj(j) + x_trial(i) * dlnfi_dxj_TP_trial(j,i)
        end do
    End Do

    End subroutine d_tpd_dxj


    !************************************************************************************
    module subroutine Succ_Sub_tpd(gl,press, temp, x, rho_orig, x_trial, tpd_min, iphase, Stable, NrofIter, errval)
    !************************************************************************************
    !Successive Substitution routine for finding a minimum of the tangent plane distance
    !function
    !Added iphase (estimated phase for density solver). Andreas February 2014
    !--------------------------------------------------------------------------
    !This algorithm is according to
    !Michelson, M.L. ; Mollerup, J.M.
    !"Thermodynamic Models: Fundamentals & Computational Aspects"
    !Tie-Line Publications, Denmark 2004
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !   rho_orig    - density of the original phase
    !   x_trial     - Trial phase composition
    !   i
    !
    ! OUTPUT:
    !
    !   tpd_min     - Vector containing the derivatives
    !   x_trial     - composition, that minimizes the tpd
    !************************************************************************************
    !Andreas Aug. 2011







    implicit none

    type(type_gl) :: gl


    double precision:: press, temp
    double precision, dimension(30):: x
    double precision:: rho_orig
    double precision, dimension(30):: x_trial
    integer :: NrofIter
    integer :: errval
    double precision :: tpd_min
    logical :: Stable

    double precision :: rho, rho_trial, sum, rhoredmix_orig, Tredmix_orig, tpd_trial
    double precision, dimension(30):: x_trial_old, dev_x, dev_x_orig, x_trial_min
    double precision, dimension(30):: fugcoeff_trial, fugcoeff
    integer:: j, i, iPhase

    x_trial_old = 0.D0
    !iPhase = 0
    sum = 0.D0
    Stable = .true.
    i = 0
    errval = 0
    tpd_min = 0.D0
    x_trial_min = 0.d0

    ! Get the density of the overall phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    gl%molfractions = x
    call reduced_parameters_calc(gl,temp)
    rhoredmix_orig = gl%rhoredmix
    Tredmix_orig = gl%tredmix
    ! calculate the liquid phase 2 density
    if (rho_orig > 0.d0) then
        rho = rho_orig
    else
        rho = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)
    end if

    !Andreas Jäger, November 2018. Changed double comparison from rho == 0.D0 to rho < 1.D-14
    !if (rho == 0.D0)then
    if (rho < 1.D-14) then
        errval = -8888
        !write (*,*) 'error ', errval, ' -- overall mixture density iteration failed in function "tpd"'
        gl%molfractions = x
        gl%rhoredmix = rhoredmix_orig
        gl%tredmix = tredmix_orig
        return
    end if

    call FUGCO_CALC_MIX(gl,Temp,rho, fugcoeff)

    ! Get the density of the trial phase
    ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
    ! Check the tangent plane distance at the initial point
    ! Check for stability of the given mixture
    tpd_trial = tpd(gl,press, temp, x, rho, x_trial, errval)
    !if (tpd_trial < -1.D-12) then
    tpd_min = tpd_trial
    x_trial_min = x_trial
    !end if

    Do while (i <=  NrofIter)

        !if (tpd_min < -1.D-12 .and. stable) then !Achtung
        if (tpd_min < -1.D-4 .and. stable) then
            stable = .false.
            !unstable -- >  carry out 3 more iteration steps
            i = Nrofiter - 3
        End if

        ! Get the density of the trial phase
        ! Take the most promising solution of the density (solution with lowest Gibbs-Energy)
        gl%molfractions = x_trial
        call reduced_parameters_calc(gl,Temp)
        ! calculate the liquid phase 2 density
        rho_trial = rhomix_calc(gl,Temp, Press, 0.D0, iPhase,0)
        !Andreas Jäger, November 2018. Changed double comparison from rho_trial == 0.D0 to dabs(rho_trial) < 1.D-14
        !    if (rho_trial == 0.D0)then
        if (rho_trial < 1.D-14) then
            errval = -8888
            !write (*,*) 'error ', errval, ' -- trial phase mixture density iteration failed in function "tpd"'
            gl%molfractions = x
            gl%rhoredmix = rhoredmix_orig
            gl%tredmix = tredmix_orig
            return
        end if
        call FUGCO_CALC_MIX(gl,Temp,rho_trial, fugcoeff_trial)

        !Gernert et al. (2014), Eq. (24)
        Do j = 1, gl%Ncomp
            x_trial(j) = fugcoeff(j) * x(j) / fugcoeff_trial(j)
            sum = sum + x_trial(j)
        End do

        x_trial = x_trial / sum

        ! Check the tangent plane distance at the initial point
        ! Check for stability of the given mixture
        tpd_trial = tpd(gl,press, temp, x, rho, x_trial, errval)
        if (tpd_min > tpd_trial) then
            tpd_min = tpd_trial
            x_trial_min = x_trial
        end if

        !    !Check for stability of the given mixture
        !    tpd_min = tpd(gl,press, temp, x, rho, x_trial, errval)
        !    if (tpd_min < -1.D-12 .and. stable) then
        !        stable = .false.
        !        !unstable -- >  carry out 3 more iteration steps
        !        i = Nrofiter - 3
        !    End if

        Do j = 1, gl%Ncomp
            !In dev_x the relative change of the molefractions compared to the last step is saved
            dev_x(j) = (x_trial(j) - x_trial_old(j)) / x_trial(j)
            !In dev_x_orig the "distance" from the original molfraction is stored
            dev_x_orig(j) = (x_trial(j) - x(j)) / x(j)
        End do

        !Exit criterion: if the molfractions don't change any more
        !OR if the composition gets too close to the overall composition
        if ((maxval(dabs(dev_x)) < 1.D-10) .or. (maxval(dabs(dev_x_orig)) < 1.D-4)) then
            gl%molfractions = x
            gl%rhoredmix = rhoredmix_orig
            gl%tredmix = tredmix_orig
            x_trial = x_trial_min
            return
        end if
        x_trial_old = x_trial
        sum = 0.D0
        i = i + 1
    end do

    !if (i > NrofIter) errval = -2222

    gl%molfractions = x
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    x_trial = x_trial_min

    End subroutine Succ_Sub_tpd


    !************************************************************************************
    module subroutine PhaseDet_td(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, rho_spec, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified mixture
    ! For this calculation the temperature and the overall density are given
    ! Variables:
    ! INPUT:
    !   Temp        - Temperature [K]
    !   x           - Overall composition of the defined mixture
    !   rho_spec    - Overall density
    !
    ! OUTPUT:
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   x_phase(5,30)   - matrix that contains all phase compositions
    !       1   x_vap       - Vapor phase composition
    !       2   x_liq1      - (lighter) liquid phase composition
    !       3   x_liq2      - (heavier) liquid phase composition
    !       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !       5   x_hyd       - hydrate phase compsition (so far only CO2 hydrates are possible)
    !   vapfrac        - Molar vapor fraction
    !   nrofphases  - number of phases found
    !   errval      - Error value
    !   press       - Pressure [MPa]
    !************************************************************************************
    !Andreas Dec. 2012
    !THIS IS A PRELIMINARY VERSION. ON THE LONG RUN, A FASTER ROUTINE WITH T AND RHO AS INPUT WILL BE DEVELOPED







    implicit none

    type(type_gl) :: gl


    double precision:: press, temp, rho_spec
    double precision, dimension(30) :: x
    double precision, dimension(5):: rho
    double precision, dimension(30,5) :: x_Phase!x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision:: vapfrac
    double precision:: dens_residuum
    integer:: errval, nrofphases
    integer, dimension(5) :: phasetype

    integer :: i
    double precision :: rho_cal

    !Variables neccessary for the regula falsi equilibrium pressure iteration
    !----------------------------------------------------------------------------
    integer :: Max_Iterations, Iterations
    double precision:: pmin, pmax, p_min_allowed, p_max_allowed, Delta_allowed
    type(type_additional_parameters) :: parameters         ! needed to transmit T, rho_spec to Regula Falsi
    Max_Iterations = 50
    Delta_allowed = 1.D-8
    !parameters = 0.D0
    parameters%a_p(1) = rho_spec
    parameters%a_p(2) = Temp
    parameters%a_p(3:32)= x
    errval = 0
    rho = 0.d0
    Iterations = 0
    !----------------------------------------------------------------------------
    if ((gl%savebounds_p .eqv. .true.) .and. ((dabs(gl%pmin_old) > 1.d-12) .and. (dabs(gl%pmax_old)  > 1.d-12))) then
        pmin = gl%pmin_old
        pmax = gl%pmax_old
    elseif ((gl%savebounds_p .eqv. .false.) .or. ((dabs(gl%pmin_old) <= 1.d-12) .and. (dabs(gl%pmax_old) <= 1.d-12))) then
        pmax = 0.D0
        pmin = 0.D0

        Do i = 1, gl%ncomp
            pmin = 1D-6 !MPa This is an arbitrary value
            pmax = pmax + gl%pmaxfluid(i) * x(i)
        end do

        phasetype = 0

        !Check if the calculation of phase equilibrium at the boundaries works (pmin and pmax)
        Do i = 1, 20
            !Try to calculate the overall density at the pressure minimum. If the calculation fails,
            !increase the pressure up to 20 times
            call PhaseDet(gl,pmin, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval == 0) then
                exit
            else
                pmin = pmin * 2.0D0
            end if
        End do
        !If after 5 steps PhaseDet still did not converge quit with error
        if (errval /= 0) return
        !Calculate the density at the minimum pressure
        if(nrofphases == 2) then
            rho_cal = 1.D0/(1.D0/rho(phasetype(1))*vapfrac + 1.D0/rho(phasetype(2)) * (1.D0-vapfrac))
        else
            rho_cal = rho(phasetype(1))
        end if
        !Check whether the density is smaller than the density specified. If not quit with error
        if (rho_cal > rho_spec) then
            errval = -9943
            return
        end if

        Do i = 1, 20
            !Try to calculate the overall density at the pressure maximum. If the calculation fails,
            !decrease the pressure up to 20 times
            call PhaseDet(gl,pmax, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval == 0) then
                exit
            else
                pmax = pmax * 0.90D0
            end if
        End do
        !If after 5 steps PhaseDet still did not converge quit with error
        if (errval /= 0) return
        !Calculate the density at the maximum pressure
        if(nrofphases == 2) then
            rho_cal = 1.D0/(1.D0/rho(phasetype(1))*vapfrac + 1.D0/rho(phasetype(2)) * (1.D0-vapfrac))
        else
            rho_cal = rho(phasetype(1))
        end if
        !Check whether the density is bigger than the density specified. If not quit with error
        if (rho_cal < rho_spec) then
            errval = -9943
            return
        end if
    end if
    !SH 01/2020
    !try first if trival solution is the correct density
    press = p_calc(gl, temp, rho_spec, 0)
    if (press < 0.d0) then
        p_min_allowed = pmin
        p_max_allowed = pmax
        call Regula_Falsi(gl,Density_dif, press, pmin, pmax, Delta_allowed, &
            p_min_allowed,p_max_allowed, Max_Iterations, Iterations, errval, parameters)
    else
        errval = 0
    endif

    !Check if quick calculation was succesful, i.e. specified density == calculated density
    if (errval == 0) then
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval == 0) then !SH 08/16: errval was not catched before
            !Check if the solution is correct. If not, quit with error
            !Andreas, Feb 2015
            if(nrofphases == 2) then
                dens_residuum = rho_spec- 1.D0/(1.D0/rho(phasetype(1))*vapfrac + 1.D0/rho(phasetype(2)) * (1.D0-vapfrac))
            else
                dens_residuum = rho_spec - rho(phasetype(1))
            end if
            !if quick solution is not the correct density try to determine pressure iteratively with regula falsi
            if (dabs(dens_residuum/rho_spec) > 1.D-3) then  !0.1 percent deviation in the density allowed. arbitrarily set
                errval = 0
                p_min_allowed = pmin
                p_max_allowed = pmax
                call Regula_Falsi(gl,Density_dif, press, pmin, pmax, Delta_allowed, &
                    p_min_allowed,p_max_allowed, Max_Iterations, Iterations, errval, parameters)

            end if
        else
            errval = -4407
        endif
    else
        errval = -4407
    end if




    if ((errval == 0).and.(gl%savebounds_p .eqv. .true.)) then
        gl%pmin_old = press * gl%min_factor
        gl%pmax_old = press
    end if

    !Calculation successful
    if (errval == 0) then
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval == 0) then !SH 08/16: errval was not catched before
            !Check if the solution is correct. If not, quit with error
            !Andreas, Feb 2015
            if(nrofphases == 2) then
                dens_residuum = rho_spec- 1.D0/(1.D0/rho(phasetype(1))*vapfrac + 1.D0/rho(phasetype(2)) * (1.D0-vapfrac))
            else
                dens_residuum = rho_spec - rho(phasetype(1))
            end if
            if (dabs(dens_residuum/rho_spec) > 1.D-3) then  !0.1 percent deviation in the density allowed. arbitrarily set
                errval = -4407
            end if
        else
            errval = -4407
        endif
    else
        errval = -4407
    end if

    !New module variable to save the phase information and pass it to the interface routines
    !Andreas April 2013
    gl%phase_id = phasetype

    end subroutine PhaseDet_td



    !This function works as zerofunction for the regula falsi routine when overall Density and Temperature are given
    !in the single phase region     Andreas Nov 2011
    Double Precision module Function Density_dif(gl,press, Parameters)

    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: press             ! Inputvariables
    Double precision :: rho_spec, Temp
    type(type_additional_parameters) :: parameters


    !Variables for PhaseDet
    double precision, dimension(30):: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_Phase    !x_vap, x_liq1, x_liq2, x_sol, x_hyd
    double precision :: vapfrac
    integer :: errval, nrofphases
    integer, dimension(5) :: phasetype
    !  --------------------------------------------------
    Density_dif = 0.d0
    rho_spec = parameters%a_p(1)
    Temp = parameters%a_p(2)
    x = parameters%a_p(3:32)
    call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    if (errval == 0) then
        if(nrofphases == 2) then
            Density_dif = rho_spec- 1.D0/(1.D0/rho(phasetype(1))*vapfrac + 1.D0/rho(phasetype(2)) * (1.D0-vapfrac))
        else
            Density_dif = rho_spec - rho(phasetype(1))
        end if
    else
        parameters%a_p(65) = errval
    end if

    End Function



    !************************************************************************************
    module subroutine PhaseDet_ps(gl,press, Temp, x, rho, x_phase, phasetype, vapfrac, s_spec, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified mixture
    ! For this calculation the pressure and the enthalpy are given
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   s_spec      - Entropy [J / mol K]
    !   x           - Overall composition of the defined mixture
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   x_vap       - Vapor phase composition
    !   x_liq1      - (lighter) liquid phase composition
    !   x_liq2      - (heavier) liquid phase composition
    !   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !   x_hyd       - hydrate phase compsition (so far only CO2 hydrates are possible)
    !   vapfrac        - Molar vapor fraction
    !   errval      - Error value
    !************************************************************************************
    !Andreas Nov. 2011








    implicit none

    type(type_gl) :: gl


    double precision :: press, s_spec, vapfrac, temp
    double precision, dimension(30) :: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_phase
    double precision, dimension(30):: x_vap, x_liq2
    integer :: errval, nrofphases

    double precision:: rhovap_est, rholiq_est  !S_CALC
    double precision:: rho_phase, curvature, T_spec
    integer:: iter, nr_pts_found, i

    !Variables necessary for the determination of the entropy limits
    double precision:: dens_vap, dens_liq, s_vap, s_liq, s_max, s_min, Tmax, Tmin
    double precision, dimension(30,5):: x_phase_pt
    integer, dimension(5) :: phasetype_pt
    double precision, dimension(5):: rho_pt

    !Variable for checking the calculated entropy
    double precision:: s_check

    !Indicate which phases are present
    integer, dimension(5) :: phasetype         !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    !E.g.: 1 phase present: liquid
    !-- >  nrofphases = 1
    !-- >  phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!


    !Variables neccessary for the regula falsi equilibrium temperature iteration
    !----------------------------------------------------------------------------
    integer :: Max_Iterations, Iterations, errorflag    ! transmission variable
    double precision:: T_min, T_max, T_min_allowed, T_max_allowed, Delta_allowed
    !warnings (Moni)
    !double precision:: Subp_dif, Meltp_dif
    type(type_additional_parameters) :: parameters         ! needed to transmit T, p to Regula Falsi
    Max_Iterations = 200
    Delta_allowed = 1.D-6
    !parameters = 0.D0
    parameters%a_p(1) = s_spec
    parameters%a_p(2) = press
    errorflag = 0
    Iterations = 0
    rho = 0.d0
    !----------------------------------------------------------------------------

    rhovap_est = 0.D0
    rholiq_est = 0.D0
    nr_pts_found = 0
    T_spec = 0.D0
    Tmin = 0.D0
    Tmax = 0.D0
    phasetype = 0

    !NEW ALGORITHM.
    !Andreas Feb. 2012
    !1) Check whether the desired entropy is valid within the range of Tmin to Tmax:
    !Calculate "overall" minimum and maximum temperature for the mixture
    Do i = 1, gl%ncomp
        Tmin = Tmin + gl%tminfluid(i) * x(i)
        Tmax = Tmax + gl%tmaxfluid(i) * x(i)
    end do

    !---------- this is a save but not very flexible boundary! J.G. 09.2012
    !Tmin = maxval(tminfluid)    !Uncommented in December 2013, Andreas
    !Tmin = minval(tminfluid)

    !Determine the minimum entropy value possible at the given pressure
    !Therefore let phasedet determine how many phases are present
    call PhaseDet(gl,press, Tmin, x, rho_pt, x_Phase_pt, phasetype_pt, vapfrac, nrofphases, errval)
    if (errval /= 0) then
        !If phasedet fails at the minimum temperature, try to raise the temperature and recalculate the equilibrium
        Do i = 1, 10
            Tmin = Tmin + 3.D0*i
            call PhaseDet(gl,press, Tmin, x, rho_pt, x_Phase_pt, phasetype_pt, vapfrac, nrofphases, errval)
            if (errval == 0) exit
        End Do
        if (errval /= 0) return !Minimum temperature iteration failed
    end if
    if (nrofphases == 1) then
        !For simplification, call the existing phase vapor
        gl%molfractions = x_Phase_pt(:,phasetype_pt(1))
        call reduced_parameters_calc(gl,Tmin)
        dens_vap = rho_pt(phasetype_pt(1))
        s_min = S_CALC(gl,Tmin,dens_vap,0)
    else
        !Vapour phase
        gl%molfractions = x_Phase_pt(:,phasetype_pt(1))
        call reduced_parameters_calc(gl,Tmin)
        dens_vap = rho_pt(phasetype_pt(1))
        s_vap = S_CALC(gl,Tmin,dens_vap,0)

        !Liquid phase
        gl%molfractions = x_Phase_pt(:,phasetype_pt(2))
        call reduced_parameters_calc(gl,Tmin)
        dens_liq = rho_pt(phasetype_pt(2))
        s_liq = S_CALC(gl,Tmin,dens_liq,0)

        !Calculate overall entropy
        s_min =  vapfrac * s_vap + (1.D0 - vapfrac) * s_liq
    end if

    !Determine the maximum entropy value possible at the given pressure
    !As the temperature is very high, a single phase (vapor) is assumed, so no phaseequilibrium
    !calculations are necessary
    gl%molfractions = x
    call reduced_parameters_calc(gl,Tmax)
    rho_phase = rhomix_calc(gl,Tmax, press, 0.D0 ,0,0)
    if (rho_phase == 0.D0) then
        errval = -8888
        return
    end if
    s_max = S_CALC(gl,Tmax, rho_phase, 0)

    !2) Check whether the specified entropy is between the maximum and minimum value. If not -- >  Calculation not possible
    if ((s_spec < s_min) .or. (s_spec > s_max)) then
        errval = -9942
        return
    end if

    !3) Calculate the specified entropy using the regula falsi routine
    !   The calculation might fail, if the enthalpy specified is in the two phase region
    T_min = Tmin
    T_min_allowed = Tmin
    T_max = Tmax
    T_max_allowed = Tmax
    call Regula_Falsi(gl,Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    rho(1) = parameters%a_p(3)
    !Andreas December 2013
    !The regula falsi might find a wrong solution, but accept it as correct. So the correctness of the solution has to be checked
    rho_phase = rho(1)
    if (temp > 0.D0) then
        s_check = S_CALC(gl,Temp, rho_phase, 0)
        if (dabs(s_spec-s_check) > 1.D-1) then
            errorflag = -9941
        end if
    else
        errorflag = -9941
    end if

    !4) Check, whether the solution found is stable or not
    ! Andreas, December 2013
    ! It is assumed that the Regula falsi method will find the correct solution, if the mixture is stable
    ! So an error (Iteration failed or Enthalpy not correct) will ONLY occur in the two phase region
    ! But the Regula Falsi can find a metastable, mathematically correct solution in the two phase region
    ! So the mixture has to be checked for stability, even if no error occurs
    if(errorflag == 0) then
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval /= 0) then
            return
        end if
    else
        nrofphases = 2
        !Start values for the following ps flash need to be generated
        !Assumption: 2 phases in equilibrium
        !Andreas, Dez. 2013
        !Try startvalue from the Regula Falsi first
        if (temp > 0.D0) then
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval /= 0) then
                return
            end if
        else
            nrofphases = 1
        end if
        !2 phases expected, if only one phase is found, try linear interpolation of the temperature
        if (nrofphases == 1) then

            !Linear interpolation of Temp
            Temp = (s_spec - s_min) / (s_max - s_min) * (Tmax - Tmin) + Tmin
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval /= 0) then
                return
            end if
            !If an error in the regula falsi occured, two phases are expected. So if the mixture is predicted to be stable return with error here
            !Andreas December 2013
            if (nrofphases == 1) then
                errval = -4406
                return
            end if

        end if
    end if

    if(nrofphases == 1) then
        !Stable phase, solution found!! No changes necessary
        x_vap = x_Phase(:,phasetype(1))
        if (rho(1) == 0.D0) rho(1) = rho(phasetype(1))
    else
        !this means the solution is in the two phase region
        x_vap = x_phase(:,phasetype(1))
        x_liq2 = x_phase(:,phasetype(2))
        rhovap_est = 0.D0
        rholiq_est = 0.D0
        call Flash_ps(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
        if(errval == 0) then
            nrofphases = 2
            rho(1) = gl%rho_vap
            rho(3) = gl%rho_liq
            !        else   ! commented out the block below, made things worse!  J.G. 10.2012
            !            !If the first try failed, try to use linear starting value for the temperature
            !            Temp = (s_spec - s_min) / (s_max - s_min) * (Tmax - Tmin) + Tmin
            !            call PhaseDet(press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            !            x_vap = x_phase(:,phasetype(1))
            !            x_liq2 = x_phase(:,phasetype(2))
            !            rhovap_est = 0.D0
            !            rholiq_est = 0.D0
            !            call Flash_ps(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
            !            if (errval == 0) then
            !                nrofphases = 2
            !                rho(1) = rho_vap
            !                rho(3) = rho_liq
            !            else
            !                return
            !            end if
        end if
    End if


    !!1) Try to calculate the bubble point
    !iFlash = 1
    !Temp = 0.D0
    !x_vap = 0.D0
    !x_liq2 = 0.D0
    !call Flash_PhaseBoundary_calc(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash,&
    !                        & 0, 0, errval, iter)
    !if(errval == 0) then
    !    !Calculate the enthalpy for the bubble point
    !    rho_return(1) = rho_liq
    !    T_return(1) = Temp
    !    s_return(1) = S_CALC(gl,T_return(1), rho_return(1), 0)
    !    !Try to calculate the dew point
    !    iFlash = 2
    !    Temp = 0.D0
    !    x_vap = 0.D0
    !    x_liq2 = 0.D0
    !    nr_pts_found = 1
    !    call Flash_PhaseBoundary_calc(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash,&
    !                        & 0, 0, errval, iter)
    !    if(errval == 0) then
    !        !Calculate the enthalpy for the dew point
    !        rho_return(2) = rho_vap
    !        T_return(2) = Temp
    !        s_return(2) = S_CALC(gl,T_return(2), rho_return(2), 0)
    !        nr_pts_found = 2
    !    End if
    !End if
    !Temp = 0.D0
    !x_vap = 0.D0
    !x_liq2 = 0.D0
    !!2) Check, whether either the bubblepoint OR the dewpoint calculation failed
    !if(nr_pts_found < 2) then !If true, at least one of the calculations failed
    !    !Set back the points found
    !    rho_return = 0.D0
    !    T_return = 0.D0
    !    s_return = 0.D0
    !    nr_pts_found = 0
    !    !Start calculating the phase envelope
    !    !call phasenv(x, press, T_return, rho_return, nr_pts_found, errval)
    !    call phasenv(x, press, T_spec, T_return, rho_return, x_calc, nr_pts_found, errval)
    !    if (errval == 0) then
    !        do i = 1, nr_pts_found-1
    !            if ((T_return(i) > T_return(i+1)) .and. (T_return(i+1) > 0.d0)) then
    !                help = T_return(i)
    !                T_return(i) = T_return(i+1)
    !                T_return(i+1) = help
    !                help = rho_return(i)
    !                rho_return(i) = rho_return(i+1)
    !                rho_return(i+1) = help
    !            end if
    !        end do
    !        if ((T_return(1) >  T_return(2)) .and. (T_return(2) > 0.d0)) then
    !            help = T_return(1)
    !            T_return(1) = T_return(2)
    !            T_return(2) = help
    !            help = rho_return(1)
    !            rho_return(1) = rho_return(2)
    !            rho_return(2) = help
    !        end if
    !        !Calculate the enthalpy for all points found
    !        Do i = 1, nr_pts_found
    !                !Calculate the enthalpy for each point found
    !                s_return(i) = S_CALC(gl,T_return(i), rho_return(i), 0)
    !        End Do
    !    else
    !        write(*,*) "The envelope calculations failed."
    !        return
    !    End if
    !End if
    !
    !!3) Determine whether the given enthalpy is in a two phase region or in the single phase and calculate
    !select case (nr_pts_found)
    !    case(3)
    !    !This means the two phase regions are BELOW point 1 and BETWEEN point 2 and 3
    !        if (s_spec < s_return(1)) then
    !            !Two phase
    !            !Startvalue for the temperature:
    !            Temp = (T_return(2) - T_return(1))/(s_return(2) - s_return(1))*(s_spec - s_return(1)) + T_return(1)
    !            call Flash_ps(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
    !            nrofphases = 2
    !            rho(1) = rho_vap
    !            rho(3) = rho_liq
    !        elseif ((s_return(3) > s_spec) .And. (s_return(2) < s_spec)) then
    !            !Two phase
    !            !Startvalue for the temperature:
    !            Temp = (s_spec - s_return(2)) / (s_return(3) - s_return(2)) * (T_return(3) - T_return(2)) + T_return(2)
    !            call Flash_ps(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
    !            nrofphases = 2
    !            rho(1) = rho_vap
    !            rho(3) = rho_liq
    !        else
    !            !single phase
    !            if(s_spec < s_return(2)) then !single phase region between the points 1 and 2
    !                T_min = T_return(1)
    !                T_min_allowed = T_return(1)
    !                T_max = T_return(2)
    !                T_max_allowed = T_return(2)
    !                call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                rho(3) = parameters(3)
    !                !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                rho_vap = rho(3)
    !            else                        !single phase region is beyond point 3
    !                T_min = T_return(3)
    !                T_min_allowed = T_return(3)
    !                T_max = minval(tmaxfluid(1:ncomp))
    !                T_max_allowed = minval(tmaxfluid(1:ncomp))
    !                call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                rho(1) = parameters(3)
    !                !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                rho_vap = rho(1)
    !            End if
    !            nrofphases = 1
    !        End if
    !    case(2)
    !    !This means the two phase reagion is BETWEEN the two points found
    !        if ((s_return(2) > s_spec) .And. (s_return(1) < s_spec)) then
    !            !Two phase
    !            !Startvalue for the temperature:
    !            Temp = (s_spec - s_return(1)) / (s_return(2) - s_return(1)) * (T_return(2) - T_return(1)) + T_return(1)
    !            call Flash_ps(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
    !            nrofphases = 2
    !            rho(1) = rho_vap
    !            rho(3) = rho_liq
    !        else
    !            !Single phase
    !            if(s_spec < s_return(1)) then !single phase region below point 1
    !                T_min = maxval(tminfluid)
    !                T_min_allowed = maxval(tminfluid)
    !                T_max = T_return(1)
    !                T_max_allowed = T_return(1)
    !                call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                rho(3) = parameters(3)
    !                !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                rho_vap = rho(3)
    !            else                        !single phase region is above point 2
    !                T_min = T_return(2)
    !                T_min_allowed = T_return(2)
    !                T_max = minval(tmaxfluid(1:ncomp))
    !                T_max_allowed = minval(tmaxfluid(1:ncomp))
    !                call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                rho(1) = parameters(3)
    !                !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                rho_vap = rho(1)
    !            End if
    !            nrofphases = 1
    !        End if
    !    case(1)
    !    !This means the two phase region is BELOW the point found
    !        if (s_spec < s_return(1)) then
    !            !Two phase
    !            ! calculate a density and enthalpy in the two-phase region close to the point on the
    !            ! phase envelope in order to extrapolate a start value for the temperature
    !            T_help = T_return(1)-10.d0
    !            call PhaseDet(press, T_help, x, rho, x_vap, x_liq1, x_liq2, x_sol, x_hyd, vapfrac, nrofphases, errval)
    !            if (errval /= 0) return
    !            molfractions = x_vap
    !            call reduced_parameters_calc
    !            d_vap = rhomix_calc(gl,T_help, press, 0.d0, 0,0)
    !            s_1 = s_CALC(gl,T_help, d_vap, 0)
    !            molfractions = x_liq2
    !            call reduced_parameters_calc
    !            d_liq = rhomix_calc(gl,T_help, press, 0.d0, 0,0)
    !            s_2 = s_CALC(gl,T_help, d_liq, 0)
    !            s_help = s_1*vapfrac + s_2*(1.d0 - vapfrac)
    !            molfractions = x
    !            call reduced_parameters_calc
    !            x_vap = 0.d0; x_liq2 = 0.d0
    !            Temp = (T_help - T_return(1))/(s_help - s_return(1))*(s_spec - s_return(1)) + T_return(1)
    !            call Flash_ps(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, s_spec, 0, errval, iter)
    !            nrofphases = 2
    !            rho(1) = rho_vap
    !            rho(3) = rho_liq
    !        else
    !            !Single phase
    !            T_min = T_return(1)
    !            T_min_allowed = T_return(1)
    !            T_max = minval(tmaxfluid(1:ncomp))
    !            T_max_allowed = minval(tmaxfluid(1:ncomp))
    !            call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !            nrofphases = 1
    !            rho(1) = parameters(3)
    !            !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !            !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !            rho_vap = rho(1)
    !        End if
    !    case(0)
    !    !This means the point is supercritical
    !    !Single phase
    !    T_min = maxval(tminfluid(1:ncomp))
    !    T_min_allowed = maxval(tminfluid(1:ncomp))
    !    T_max = minval(tmaxfluid(1:ncomp))
    !    T_max_allowed = minval(tmaxfluid(1:ncomp))
    !    call Regula_Falsi(Entropy_dif, Temp, T_min, T_max, Delta_allowed, &
    !        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !    nrofphases = 1
    !    rho(1) = parameters(3)
    !    !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !    !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !    rho_vap = rho(3)
    !End select

    if (nrofphases == 1) then

        gl%molfractions = x
        call reduced_parameters_calc(gl,Temp)
        rho_phase = rho(1)
        ! check if the original phase is vapor like or liquid like
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        ! write the original density on the correct position in the return vector
        if (curvature > 0.d0) then
            rho(3) = rho_phase
            phasetype(1) = 3
            x_phase(:,3) = x
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            x_phase(:,1) = x
        end if

    elseif(nrofphases == 2) then

        !Determine which phases are present
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = rho(1)
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        if (curvature > 0.d0) then      !Light liquid phase
            rho(2) = rho_phase
            phasetype(1) = 2
            x_phase(:,2) = x_vap
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            x_phase(:,1) = x_vap
        end if

        gl%molfractions = x_liq2
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = rho(3)
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        !Some cases might occur, where a (correct) liquid phase solution has a negative curvature and therefore
        !would be indicated as vapor. If this happens, the correct vapor solution would be overwritten. To prevent this the solution
        !will ALWAYS be stored at the third position as heavy liquid.

        rho(3) = rho_phase
        phasetype(2) = 3
        x_phase(:,3) = x_liq2

        !    if (curvature > 0.d0) then
        !        rho(3) = rho_phase
        !        phasetype(2) = 3
        !        x_phase(:,3) = x_liq2
        !    else
        !        rho(1) = rho_phase
        !        phasetype(2) = 1
        !        x_phase(:,1) = x_liq2
        !    end if

    End if


    !New module variable to save the phase information and pass it to the interface routines
    !Andreas April 2013
    gl%phase_id = phasetype

    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)

    end subroutine PhaseDet_ps


    !This function works as zerofunction for the regula falsi routine when Entropy and pressure are given and the respective temperature is needed
    !in the single phase region     Andreas Nov 2011
    Double Precision module Function Entropy_dif(gl,Temp, Parameters)


    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: Temp, press, s_spec  ! Inputvariables
    Double Precision ::  dens             ! Functions for calculation of the Enthalpy and Density      !S_CALC
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------
    s_spec = parameters%a_p(1)
    press = parameters%a_p(2)
    dens = rhomix_calc(gl,Temp, Press, 0.D0, 0,0)
    Entropy_dif = s_spec-S_CALC(gl,Temp, dens, 0)
    parameters%a_p(3) = dens
    !if (dens == 0.d0) parameters(65) = -999

    End Function Entropy_dif


    !************************************************************************************
    module subroutine PhaseDet_ph(gl,press, Temp, x, rho, x_phase, phasetype, vapfrac, h_spec, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified mixture
    ! For this calculation the pressure and the enthalpy are given
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   h_spec      - Enthalpy [J / mol]
    !   x           - Overall composition of the defined mixture
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - densities of all phases found, vector, length 5 (vapor, (light) liquid, (heavy) liquid, solid, hydrate)
    !   x_phase(5,30)   - matrix that contains all phase compositions
    !       1   x_vap       - Vapor phase composition
    !       2   x_liq1      - (lighter) liquid phase composition
    !       3   x_liq2      - (heavier) liquid phase composition
    !       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !       5   x_hyd       - hydrate phase compsition (so far only CO2 hydrates are possible)
    !   vapfrac        - Molar vapor fraction
    !   errval      - Error value
    !************************************************************************************
    !Andreas Nov. 2011








    implicit none

    type(type_gl) :: gl


    double precision :: press, h_spec, vapfrac, temp
    double precision, dimension(30) :: x
    double precision, dimension(5) :: rho
    double precision, dimension(30,5) :: x_phase
    double precision, dimension(30):: x_vap, x_liq2
    integer :: errval, nrofphases

    double precision:: rhovap_est, rholiq_est  !H_CALC
    !warnings (Moni)
    !double precision:: rho_help
    double precision:: rho_phase, curvature, T_spec
    integer:: iter, nr_pts_found, i

    !Variables necessary for the determination of the enthalpy limits
    double precision:: dens_vap, dens_liq, h_vap, h_liq, h_max, h_min, Tmax, Tmin
    double precision, dimension(30,5):: x_phase_pt
    integer, dimension(5) :: phasetype_pt
    double precision, dimension(5):: rho_pt

    !Variable for checking the calculated enthalpy
    double precision:: h_check


    !Indicate which phases are present
    integer, dimension(5) :: phasetype         !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    !E.g.: 1 phase present: liquid
    !-- >  nrofphases = 1
    !-- >  phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!

    !Variables neccessary for the regula falsi equilibrium temperature iteration
    !----------------------------------------------------------------------------
    integer :: Max_Iterations, Iterations, errorflag    ! transmission variable
    double precision:: T_min, T_max, T_min_allowed, T_max_allowed, Delta_allowed
    !warnings (Moni)
    !double precision:: psub_eq, Meltp_dif, Subp_dif, pmelt_eq
    type(type_additional_parameters) :: parameters         ! needed to transmit T, p to Regula Falsi
    Max_Iterations = 200
    Delta_allowed = 1.D-6
    !parameters = 0.D0
    parameters%a_p(1) = h_spec
    parameters%a_p(2) = press
    errorflag = 0
    rho = 0.d0
    Iterations = 0
    !----------------------------------------------------------------------------

    rhovap_est = 0.D0
    rholiq_est = 0.D0
    nr_pts_found = 0
    T_spec = 0.D0
    Tmin = 0.D0
    Tmax = 0.D0
    phasetype = 0
    h_check = 0.D0

    !NEW ALGORITHM.
    !Andreas Feb. 2012
    !1) Check whether the desired enthalpy is valid within the range of Tmin to Tmax:
    !Calculate "overall" minimum and maximum temperature for the mixture
    Do i = 1, gl%ncomp
        Tmin = Tmin + gl%tminfluid(i) * x(i)
        Tmax = Tmax + gl%tmaxfluid(i) * x(i)
    end do

    
    ! -------------------------------------------------
    ! Sven and Robin : Advanced temperature scaling 
    ! -------------------------------------------------
    if(Tmin/2d0 > minval(gl%tminfluid(1:gl%ncomp))) then
        Tmin =  Tmin / 2d0
    else
        Tmin = minval(gl%tminfluid(1:gl%ncomp))
    end if

    if(Tmax/2d0 > minval(gl%tminfluid(1:gl%ncomp))) then
        Tmax =  Tmax / 2d0
    else
        Tmax = minval(gl%tminfluid(1:gl%ncomp))
    end if
    ! -------------------------------------------------
    
    !---------- this is a save but not very flexible boundary! J.G. 09.2012
    !Tmin = maxval(tminfluid)   !Uncommented and changed in December 2013, Andreas
    !Tmin = minval(tminfluid)

    !Determine the minimum enthalpy value possible at the given pressure
    !Therefore let phasedet determine how many phases are present
    call PhaseDet2(gl,press, Tmin, x, rho_pt, x_Phase_pt, phasetype_pt, vapfrac, nrofphases, errval)
    if (errval /= 0) then
        !If phasedet fails at the minimum temperature, try to raise the temperature and recalculate the equilibrium
        Do i = 1, 10
            Tmin = Tmin + 3.D0*i
            call PhaseDet2(gl,press, Tmin, x, rho_pt, x_Phase_pt, phasetype_pt, vapfrac, nrofphases, errval)
            if (errval == 0) exit
        End Do
        if (errval /= 0) return !Minimum temperature iteration failed
    end if
    if (nrofphases == 1) then
        !For simplification, call the existing phase vapor
        gl%molfractions = x_Phase_pt(:,phasetype_pt(1))
        call reduced_parameters_calc(gl,Tmin)
        dens_vap = rho_pt(phasetype_pt(1))
        h_min = H_CALC(gl,Tmin,dens_vap,0)
    else
        !Vapour phase
        gl%molfractions = x_Phase_pt(:,phasetype_pt(1))
        call reduced_parameters_calc(gl,Tmin)
        dens_vap = rho_pt(phasetype_pt(1))
        h_vap = H_CALC(gl,Tmin,dens_vap,0)

        !Liquid phase
        gl%molfractions = x_Phase_pt(:,phasetype_pt(2))
        call reduced_parameters_calc(gl,Tmin)
        dens_liq = rho_pt(phasetype_pt(2))
        h_liq = H_CALC(gl,Tmin,dens_liq,0)

        !Calculate overall enthalpy
        h_min =  vapfrac * h_vap + (1.D0 - vapfrac) * h_liq
    end if

    !Determine the maximum enthalpy value possible at the given pressure
    !As the temperature is very high, a single phase (vapor) is assumed, so no phaseequilibrium
    !calculations are necessary
    gl%molfractions = x
    do i = 1,12
        call reduced_parameters_calc(gl,Tmax)
        rho_phase = rhomix_calc(gl,Tmax, press, 0.D0 ,0,0)
        if (rho_phase > 0) exit
        Tmax = Tmax - 100.d0
    end do
    if (rho_phase == 0.D0) then
        errval = -8888
        return
    end if
    h_max = H_CALC(gl,Tmax, rho_phase, 0)

    !2) Check whether the specified enthalpy is between the maximum and minimum value. If not -- >  Calculation not possible
    if ((h_spec < h_min) .or. (h_spec > h_max)) then
        errval = -9941
        return
    end if

    !3) Calculate the specified enthalpy using the regula falsi routine
    !   The calculation might fail, if the enthalpy specified is in the two phase region
    T_min = Tmin
    T_min_allowed = Tmin
    T_max = Tmax
    T_max_allowed = Tmax
    call Regula_Falsi(gl,Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !SH: Delete after Validation is completed (2016-10-24)
    !if (errorflag == 3) errorflag = 0
    rho(1) = parameters%a_p(3)
    !Andreas December 2013
    !The regula falsi might find a wrong solution, but accept it as correct. So the correctness of the solution has to be checked
    rho_phase = rho(1)
    if (temp > 0.D0) then
        h_check = H_CALC(gl,Temp, rho_phase, 0)
        if (dabs(h_spec-h_check) > 1.D-1) then
            errorflag = -4405
        end if
    else
        errorflag = -4405
    end if


    !4) Check, whether the solution found is stable or not
    ! Andreas, December 2013
    ! It is assumed that the Regula falsi method will find the correct solution, if the mixture is stable
    ! So an error (Iteration failed or Enthalpy not correct) will ONLY occur in the two phase region
    ! But the Regula Falsi can find a metastable, mathematically correct solution in the two phase region
    ! So the mixture has to be checked for stability, even if no error occurs
    if(errorflag == 0) then
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval /= 0) then
            return
        end if
    else
        nrofphases = 2
        !Start values for the following ph flash need to be generated
        !Assumption: 2 phases in equilibrium
        !Andreas, Dez. 2013
        !Try startvalue from the Regula Falsi first
        if (temp > 0.D0) then
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval /= 0) then
                return
            end if
        else
            nrofphases = 1
        end if

        !2 phases expected, if only one phase is found, try linear interpolation of the temperature
        if (nrofphases == 1) then

            !If a stable phase was found, try linear interpolation of Temp
            Temp = (h_spec - h_min) / (h_max - h_min) * (Tmax - Tmin) + Tmin
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            if (errval /= 0) then
                return
            end if
            !If an error in the regula falsi occured, two phases are expected. So if the mixture is predicted to be stable return with error here
            !Andreas December 2013
            if (nrofphases == 1) then
                errval = -4405
                return
            end if

        end if
    end if

    if (nrofphases == 1) then
        !Stable phase, solution found!! No changes necessary
        x_vap = x_Phase(:,phasetype(1))
        if (rho(1) == 0.D0) rho(1) = rho(phasetype(1))
    else
        !this means the solution is in the two phase region
        x_vap = x_phase(:,phasetype(1))
        x_liq2 = x_phase(:,phasetype(2))
        rhovap_est = 0.D0
        rholiq_est = 0.D0
        call Flash_ph(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
        if(errval == 0) then
            nrofphases = 2
            rho(1) = gl%rho_vap
            rho(3) = gl%rho_liq
        else
            !If the first try failed, try to use linear starting value for the temperature
            Temp = (h_spec - h_min) / (h_max - h_min) * (Tmax - Tmin) + Tmin
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
            x_vap = x_phase(:,phasetype(1))
            x_liq2 = x_phase(:,phasetype(2))
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            call Flash_ph(gl,press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
            if (errval == 0) then
                nrofphases = 2
                rho(1) = gl%rho_vap
                rho(3) = gl%rho_liq
            else
                return
            end if
        end if
    End if

    !    !1) Try to calculate the bubble point
    !    iFlash = 1
    !    Temp = 0.D0
    !    x_vap = 0.D0
    !    x_liq2 = 0.D0
    !    call Flash_PhaseBoundary_calc(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash,&
    !                            & 0, 0, errval, iter)
    !    if(errval == 0) then
    !        !Calculate the enthalpy for the bubble point
    !        rho_return(1) = rho_liq
    !        t_return(1) = Temp
    !        h_return(1) = H_CALC(gl,t_return(1), rho_return(1), 0)
    !        !Try to calculate the dew point
    !        iFlash = 2
    !        Temp = 0.D0
    !        x_vap = 0.D0
    !        x_liq2 = 0.D0
    !        nr_pts_found = 1
    !        call Flash_PhaseBoundary_calc(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, iFlash,&
    !                            & 0, 0, errval, iter)
    !        if(errval == 0) then
    !            !Calculate the enthalpy for the dew point
    !            rho_return(2) = rho_vap
    !            t_return(2) = Temp
    !            h_return(2) = H_CALC(gl,t_return(2), rho_return(2), 0)
    !            nr_pts_found = 2
    !        End if
    !    End if
    !    Temp = 0.D0
    !    x_vap = 0.D0
    !    x_liq2 = 0.D0
    !    !2) Check, whether either the bubblepoint OR the dewpoint calculation failed
    !    if(nr_pts_found < 2) then !If true, at least one of the calculations failed
    !        !Set back the points found
    !        rho_return = 0.D0
    !        t_return = 0.D0
    !        h_return = 0.D0
    !        nr_pts_found = 0
    !        !Start calculating the phase envelope
    !        !call phasenv(x, press, t_return, rho_return, nr_pts_found, errval)
    !        call phasenv(x, press, T_spec, t_return, rho_return, x_calc, nr_pts_found, errval)
    !        if (errval == 0) then
    !            do i = 1, nr_pts_found-1
    !                if ((t_return(i) > t_return(i+1)) .and. (t_return(i+1) > 0.d0)) then
    !                    help = t_return(i)
    !                    t_return(i) = t_return(i+1)
    !                    t_return(i+1) = help
    !                    help = rho_return(i)
    !                    rho_return(i) = rho_return(i+1)
    !                    rho_return(i+1) = help
    !                end if
    !            end do
    !            if ((t_return(1) >  t_return(2)) .and. (t_return(2) > 0.d0)) then
    !                help = t_return(1)
    !                t_return(1) = t_return(2)
    !                t_return(2) = help
    !                help = rho_return(1)
    !                rho_return(1) = rho_return(2)
    !                rho_return(2) = help
    !            end if
    !            !Calculate the enthalpy for all points found
    !            Do i = 1, nr_pts_found
    !                    !Calculate the enthalpy for each point found
    !                    h_return(i) = H_CALC(gl,t_return(i), rho_return(i), 0)
    !            End Do
    !        else
    !            write(*,*) "The envelope calculations failed."
    !            return
    !        End if
    !    End if
    !
    !    !3) Determine whether the given enthalpy is in a two phase region or in the single phase and calculate
    !    select case (nr_pts_found)
    !        case(3)
    !        !This means the two phase regions are BELOW point 1 and BETWEEN point 2 and 3
    !            if (h_spec < h_return(1)) then
    !                !Two phase
    !                !Startvalue for the temperature:
    !                Temp = (t_return(2) - t_return(1))/(h_return(2) - h_return(1))*(h_spec - h_return(1)) + t_return(1)
    !                call Flash_ph(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
    !                nrofphases = 2
    !                rho(1) = rho_vap
    !                rho(3) = rho_liq
    !            elseif ((h_return(3) > h_spec) .And. (h_return(2) < h_spec)) then
    !                !Two phase
    !                !Startvalue for the temperature:
    !                Temp = (h_spec - h_return(2)) / (h_return(3) - h_return(2)) * (t_return(3) - t_return(2)) + t_return(2)
    !                call Flash_ph(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
    !                nrofphases = 2
    !                rho(1) = rho_vap
    !                rho(3) = rho_liq
    !            else
    !                !single phase
    !                if(h_spec < h_return(2)) then !single phase region between the points 1 and 2
    !                    T_min = t_return(1)
    !                    T_min_allowed = t_return(1)
    !                    T_max = t_return(2)
    !                    T_max_allowed = t_return(2)
    !                    call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                    rho(3) = parameters(3)
    !                    !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                    !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                    rho_vap = rho(3)
    !                else                        !single phase region is beyond point 3
    !                    T_min = t_return(3)
    !                    T_min_allowed = t_return(3)
    !                    T_max = minval(tmaxfluid(1:ncomp))
    !                    T_max_allowed = minval(tmaxfluid(1:ncomp))
    !                    call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                    rho(1) = parameters(3)
    !                    !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                    !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                    rho_vap = rho(1)
    !                End if
    !                nrofphases = 1
    !            End if
    !        case(2)
    !        !This means the two phase region is BETWEEN the two points found
    !            if ((h_return(2) > h_spec) .And. (h_return(1) < h_spec)) then
    !                !Two phase
    !                !Startvalue for the temperature:
    !                Temp = (h_spec - h_return(1)) / (h_return(2) - h_return(1)) * (t_return(2) - t_return(1)) + t_return(1)
    !                call Flash_ph(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
    !                nrofphases = 2
    !                rho(1) = rho_vap
    !                rho(3) = rho_liq
    !            else
    !                !Single phase
    !                if(h_spec < h_return(1)) then !single phase region below point 1
    !                    T_min = maxval(tminfluid)
    !                    T_min_allowed = maxval(tminfluid)
    !                    T_max = t_return(1)
    !                    T_max_allowed = t_return(1)
    !                    call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                    rho(3) = parameters(3)
    !                    !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                    !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                    rho_vap = rho(3)
    !                else                        !single phase region is above point 2
    !                    T_min = t_return(2)
    !                    T_min_allowed = t_return(2)
    !                    T_max = minval(tmaxfluid(1:ncomp))
    !                    T_max_allowed = minval(tmaxfluid(1:ncomp))
    !                    call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                        T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                    rho(1) = parameters(3)
    !                    !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                    !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                    rho_vap = rho(1)
    !                End if
    !                nrofphases = 1
    !            End if
    !        case(1)
    !        !This means the two phase region is BELOW the point found
    !            if (h_spec < h_return(1)) then
    !                !Two phase
    !                ! calculate a density and enthalpy in the two-phase region close to the point on the
    !                ! phase envelope in order to extrapolate a start value for the temperature
    !                ! If the generated starting values is bad, reduce the temperature up to 20 times by 10 K
    !                Do i = 1, 40
    !                    T_help = t_return(1)-10.d0
    !                    !call PhaseDet(press, T_help, x, rho, x_vap, x_liq1, x_liq2, x_sol, x_hyd, vapfrac, nrofphases, errval)
    !                    call PhaseDet(press, T_help, x, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
    !                    if (errval /= 0) return
    !                    molfractions = x_Phase(:,phasetype(1))
    !                    call reduced_parameters_calc
    !                    d_vap = rho(phasetype(1))
    !                    h_1 = h_CALC(gl,T_help, d_vap, 0)
    !                    molfractions = x_Phase(:,phasetype(2))
    !                    call reduced_parameters_calc
    !                    d_liq = rho(phasetype(2))
    !                    h_2 = h_CALC(gl,T_help, d_liq, 0)
    !                    h_help = h_1*vapfrac + h_2*(1.d0 - vapfrac)
    !                    molfractions = x
    !                    call reduced_parameters_calc
    !                    x_vap = x_Phase(:,phasetype(1))
    !                    x_liq2 = x_Phase(:,phasetype(2))
    !                    Temp = (T_help - t_return(1))/(h_help - h_return(1))*(h_spec - h_return(1)) + t_return(1)
    !                    call Flash_ph(press, Temp, x, x_vap, x_liq2, rhovap_est, rholiq_est, vapfrac, h_spec, 0, errval, iter)
    !                    if (errval == 0) then
    !                       exit
    !                    else
    !                        t_return(1) = T_help
    !                        h_return(1) = h_help
    !                        x_vap = 0.D0
    !                        x_liq2 =0.D0
    !                    end if
    !                end do
    !                nrofphases = 2
    !                rho(1) = rho_vap
    !                rho(3) = rho_liq
    !            else
    !                !Single phase
    !                T_min = t_return(1)
    !                T_min_allowed = t_return(1)
    !                T_max = minval(tmaxfluid(1:ncomp))
    !                T_max_allowed = minval(tmaxfluid(1:ncomp))
    !                call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !                    T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !                nrofphases = 1
    !                rho(1) = parameters(3)
    !                !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !                !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !                rho_vap = rho(1)
    !            End if
    !        case(0)
    !        !This means the point is supercritical
    !        !Single phase
    !        T_min = maxval(tminfluid(1:ncomp))
    !        T_min_allowed = maxval(tminfluid(1:ncomp))
    !        T_max = minval(tmaxfluid(1:ncomp))
    !        T_max_allowed = minval(tmaxfluid(1:ncomp))
    !        call Regula_Falsi(Enthalpy_dif, Temp, T_min, T_max, Delta_allowed, &
    !            T_min_allowed,T_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !        nrofphases = 1
    !        rho(1) = parameters(3)
    !        !The density is stored in the module variable rho_vap in EVERY SINGLE PHASE CASE!!!
    !        !If it is really a vapor phase or rather a liquid phase will be determined at the end!!
    !        rho_vap = rho(1)
    !    End select
    !

    if (nrofphases == 1) then

        gl%molfractions = x
        call reduced_parameters_calc(gl,Temp)
        rho_phase = rho(1)
        ! check if the original phase is vapor like or liquid like
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        ! write the original density on the correct position in the return vector
        if (curvature > 0.d0) then
            rho(3) = rho_phase
            phasetype(1) = 3
            x_phase(:,3) = x
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            x_phase(:,1) = x
        end if

    elseif(nrofphases == 2) then

        !Determine which phases are present
        gl%molfractions = x_vap
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = rho(1)
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        if (curvature > 0.d0) then      !Light liquid phase
            rho(2) = rho_phase
            phasetype(1) = 2
            x_phase(:,2) = x_vap
        else
            rho(1) = rho_phase
            phasetype(1) = 1
            x_phase(:,1) = x_vap
        end if

        gl%molfractions = x_liq2
        call reduced_parameters_calc(gl,Temp)
        ! check if the phase is vapor (like) or liquid (like)
        rho_phase = rho(3)
        curvature = d2P_drho2 (gl,Temp, rho_phase)
        ! negative curvature means: vapor like phase, positive curvature: liquid like phase
        !Some cases might occur, where a (correct) liquid phase solution has a negative curvature and therefore
        !would be indicated as vapor. If this happens, the correct vapor solution would be overwritten. To prevent this the solution
        !will ALWAYS be stored at the third position as heavy liquid.

        rho(3) = rho_phase
        phasetype(2) = 3
        x_phase(:,3) = x_liq2

        !    if (curvature > 0.d0) then
        !        rho(3) = rho_phase
        !        phasetype(2) = 3
        !        x_phase(:,3) = x_liq2
        !    else
        !        rho(1) = rho_phase
        !        phasetype(2) = 1
        !        x_phase(:,1) = x_liq2
        !    end if

    End if

    !New module variable to save the phase information and pass it to the interface routines
    !Andreas April 2013
    gl%phase_id = phasetype

    gl%molfractions = x
    call reduced_parameters_calc(gl,Temp)

    end subroutine PhaseDet_ph


    !This function works as zerofunction for the regula falsi routine when Enthalpy and pressure are given and the respective temperature is needed
    !in the single phase region     Andreas Nov 2011
    Double Precision module Function Enthalpy_dif(gl,Temp, Parameters)



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: Temp, press, h_spec  ! Inputvariables
    Double Precision :: dens             ! Functions for calculation of the Enthalpy and Density !H_CALC
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    press = parameters%a_p(2)
    h_spec = parameters%a_p(1)
    dens = rhomix_calc(gl,Temp, Press, 0.D0, 0,0)
    parameters%a_p(3) = dens
    Enthalpy_dif = h_spec-H_CALC(gl,Temp, dens, 0)

    End Function Enthalpy_dif


    end submodule impl
