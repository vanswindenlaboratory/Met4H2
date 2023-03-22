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

    ! module for file rhomix_pt.f90
    submodule (rhomix_pt_module) impl
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use ancillary_equations_mix_module
    use cubic_eos_module
    use setup_module
    use fnrderivs_module
    contains




    !**************************************************************************
    !           --------------------------------------------------
    !           Function for the iterative calculation of the density
    !           of a mixture as a function of Temperature and pressure
    !
    !           J. Gernert, Boulder, Aug. 2010
    !           A. Jäger March 2011 (modified)
    !           --------------------------------------------------
    !**************************************************************************
    !
    !**************************************************************************
    double precision module function rhomix_calc(gl,TEMPERATURE, PRESSURE, RHO_EST_GIVEN, IPHASE, nrsubst)
    !**************************************************************************
    !DEC$ ATTRIBUTES DLLEXPORT :: rhomix_calc
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! PRESSURE    - P   MPa
    ! RHO_EST_GIVEN - Initial value for the density iteration, can be zero if no inital is known
    ! IPHASE      - Indicates, for which phase the density has to be calculated
    !               1: liquid phase
    !               2: vapor phase
    !               0: no information on the phase
    !--------------------------------------------------------------------------


    implicit none

    type(type_gl) :: gl


    double precision :: TEMPERATURE, PRESSURE, RHO_EST_GIVEN
    integer:: IPHASE,nrsubst

    integer:: IFound, errval, j, k, o
    double precision:: rho_est, rho_old, rho_new, Zerofunction, Deriv
    double precision:: p_min, p_max, T, p, rholiq, rhovap, rho_save, gibbs_vap, gibbs_liq, delta_rho! , saddle_point
    double precision:: gibbs_vap_id, gibbs_liq_id, gibbs_diff_id, Rmix
    double precision:: slope, curvature, enh_fac
    logical:: min_found, max_found, extrem, quick_liq
    integer::  i
    !Variables neccessary for the regula falsi iteration
    !----------------------------------------------------------------------------
    integer :: Max_Iterations, Iterations, errorflag    ! transmission variable
    double precision:: rhomin,rhomax, rhomin_allowed, rhomax_allowed, Delta_allowed, rho_root, rho_atmin, rho_atmax, factorpress_inv
    type(type_additional_parameters) :: parameters         ! needed to transmit T, p to Regula Falsi

    !SRK
    Double precision:: rholiqtest, rhovaptest
    !PR, Stefan Feb 2014
    Double precision:: rhol_est, rhov_est
    !Check if density is vapor phase
    Double precision:: rhomix_dummy,ptest,rhoc_comp

    factorpress_inv = 1.d0/gl%factorpress

    Max_Iterations = 30
    Delta_allowed = 1.D-6
    !parameters = 0.D0

    errorflag = 0
    Iterations = 0
    rho_root = 0.d0

    T = TEMPERATURE
    P = pressure
    !Catch NaN
    if ((t /= t).or.(p /= p)) then
        rhomix_calc = -8888.d0
        return
    endif
    parameters%a_p(1) = T
    parameters%a_p(2) = nrsubst
    extrem = .false.
    min_found = .false.
    max_found = .false.
    quick_liq = .false.
    rhomix_calc = 0.d0

    rholiq = 0.d0
    rho_atmin = 0.D0
    rho_atmax = 0.D0

    p_min = 0.D0
    p_max = 0.D0

    !----------------------------------------------------------------------------------
    if ((IPHASE == 1) .or. (IPHASE == 0)) then
        !----------------------------------------------------------------------------------
        if ((gl%mix_Type == 1).OR.(gl%mix_Type == 4) .OR. (gl%mix_Type == 6) .OR. (gl%mix_Type == 7) .OR. (gl%mix_type == 11) .OR. (gl%mix_type == 12) .OR. (gl%mix_type == 13) &
            &  .OR. (gl%mix_type == 110) .OR. (gl%mix_type == 111) .OR. (gl%mix_type == 120) .OR. (gl%mix_type == 121)) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012

            !-----------------------------------------------------
            ! Step 1: Try calculating the density using initial values
            !-----------------------------------------------------
            if (rho_est_given > 0.d0) then
                rho_est = rho_est_given
            else
                rho_est = PSRK(gl,T, P, 1, nrsubst)
                !rho_est = rho_est*1.09d0 ! TE DANGER

                ! in case a pure fluid is treated, unreasonable results of the SRK can be corrected, if PHASETYPE is known  RS, 09/2014

                if (nrsubst /= 0) then ! pure fluid        MT 02/2016: DO we need it here? Already tested before...

                    if ((iphase == 1).and.(T < 0.98*gl%TC(nrsubst)) .and. (gl%mix_Type == 1)) Then ! Liquid clearly below Tc; ancillaries only available for Helmholtz EOS
                        rhol_est = dl_eq(gl,T, nrsubst)
                        if (rho_est < 0.98*rhol_est) rho_est = rhol_est/0.95D0  ! Lower limit for Regula Falsi corresponds to saturated liquid density of anc. eq.
                    end if
                    if ((iphase == 1) .and. rho_est < gl%rhoc(nrsubst)) then
                        rho_est = gl%tc(nrsubst)/TEMPERATURE*gl%rhoc(nrsubst)
                    end  if
                end if

            end if

            ! set the boundaries for the Regula Falsi
            rhomax = 1.10D0*rho_est
            rhomin = 0.90D0*rho_est

            ! more narrow interval for pure liquid    RS, 09/2014
            if ((nrsubst == 1).and.(iphase == 1).and.(gl%mix_type/=6)) then
                rhomax = 1.05D0*rho_est
                rhomin = 0.95D0*rho_est
            end if

            rhomax_allowed = 2.D0*rho_est  !MT, for Lennard-Jones (original value was 1.5D0)
            rhomin_allowed = 0.7D0*rho_est

            if(gl%seacalc) then
                rhomax = 1.07d0 * rho_est
                rhomin = 0.84d0 * rho_est
                rhomax_allowed = 1.15d0 * rho_est
                rhomin_allowed = 0.6d0 * rho_est
            end if

            ptest = P_CALC(gl, temperature, rhomax, nrsubst)
            if (ptest < 1.d-12) then
                rhomax = (rhomin_allowed+rhomax_allowed)/2.d0
                rhomin = rhomax * 0.8d0
            else if (isnan(ptest)) then
                rhomax = rhomax*0.7d0
                rhomin = rhomax * 0.5d0
                rhomin_allowed = rhomin * 0.1d0
            end if

            ! call the iteration routine for rhomix(T,P)
            rholiq = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)

            !Validation failed, SH, 01/2015
            !! For pure fluids no additional stability tests are required = >  Quick exit!  RS, 09/2014
            !
            !If (nrsubst == 1) Then
            !    rhomix_calc = rholiq
            !    return
            !end if

            !if(gl%seacalc) then
            !    continue
            !end if

            ! check if the result s physically reasonable
            if (rholiq > 0.d0) then
                if(gl%seacalc) gl%sea%seap_save = gl%sea%seap
                call rho_test (gl,T, p, rholiq, 1, IFound,nrsubst)
                if(gl%seacalc) gl%sea%seap = gl%sea%seap_save
                if (IFound < 1) then
                    ! if an initial estimate for rho is given, check if the value found is in the vapor phase
                    if (rho_est_given > 0.d0) then
                        call rho_test (gl,T, p, rholiq, 2, IFound,nrsubst)
                        if (iFound > 0) then
                            rhomix_calc = rholiq    ! quick exit
                            return
                        end if
                    end if
                    rholiq = 0.d0
                    quick_liq = .false.
                else
                    quick_liq = .true.
                    if (IPHASE == 1) then
                        rhomix_calc = rholiq    ! quick exit
                        return
                    end if
                end if
            end if

            !if(gl%seacalc) then
            !    continue
            !end if

            !-----------------------------------------------------
            ! Step 2: if Step 1 failed, find the minimum with the highest density, using Newton
            !-----------------------------------------------------
            if (dabs(rholiq) < 1.d-14) then
                ! choose a very high density, where the fluid is definity a liquid
                if (gl%mix_type == 6) then !SH: This is a heuristic added by Theresa for SAFT which leads to errors for hydrates
                    if (gl%rhoredmix /= 1.d0) then
                        rhomax = gl%rhoredmix*4.2d0
                    else if (nrsubst == 0) then
                        rhoc_comp = 0.d0
                        do j = 1, gl%ncomp
                            rhoc_comp = rhoc_comp + gl%molfractions(j)/gl%rhoc(j)
                        end do
                        rhoc_comp = 1.d0/rhoc_comp
                        rhomax = 3.2d0*rhoc_comp
                    end if
                else
                    rhomax = gl%rhoredmix*4.2d0
                endif
                if (nrsubst /= 0) rhomax = gl%rhoc(nrsubst)*4.2d0
                if (gl%mix_type == 19) then
                    !The one-fluid model is non-corresponding states based. Instead of rhoredmix, the covolume of the SRK is used here to define the maximum density
                    rhomax = 0.D0
                    do i=1,gl%ncomp
                        rhomax = rhomax + gl%molfractions(i) * 0.08664D0*8.3144589D0* gl%tc(i) / (gl%pc(i)*1.D6/gl%factorpress)
                    end do
                    rhomax = 1.D0/rhomax * 1.2d0
                end if
                !SAFT
                ! first, check if the EOS still behaves reasonably at such high densities
                ! This means: slope and curvature must be positive
                slope = DPDD_CALC(gl,T, rhomax, nrsubst)*1d6*factorpress_inv
                curvature = D2PDD2_CALC (gl,T, rhomax,nrsubst)*1d6*factorpress_inv
                ! if this is not the case, reduce the density rhomax until the three criteria are fulfilled
                ! this is necessary for the newton algorithm to work! - >  J. Gernert, April 2012
                i = 0
                do while ((slope < 0.d0) .OR. (curvature < 0.d0))
                    rhomax = rhomax*0.95
                    slope = DPDD_CALC(gl,T, rhomax, nrsubst)*1d6*factorpress_inv
                    curvature = D2PDD2_CALC (gl,T, rhomax,nrsubst)*1d6*factorpress_inv
                    i = i + 1
                    if (i > 10) exit
                end do
                rho_old = rhomax
                !rho_old = maxval(rhomaxfluid)
                rho_new = rho_old
                errval = 0
                rhomin = 0.d0
                do k = 1, 50
                    ZeroFunction = DPDD_CALC(gl,T, rho_old, nrsubst)*1.d6*factorpress_inv
                    ! break criterion for Newton
                    if (abs(ZeroFunction) < 1.d-8) exit !changed from 1.d-6 to 1.d-8 - J.G. 2.2012
                    Deriv = D2PDD2_CALC (gl,T, rho_old,nrsubst)*1.d6*factorpress_inv
                    if (gl%mix_type == 6) then !SH: This is a heuristic added by Theresa for SAFT which leads to errors for hydrates
                        delta_rho = Zerofunction/Deriv/2.d0
                    else
                        delta_rho = Zerofunction/Deriv
                    endif
                    rho_new = rho_old - delta_rho
                    j = 1
                    do while ((rho_new < 0.d0) .OR. ((rho_new > rhomax) .and. (gl%mix_type /= 19))) !Added an exception for one-fluid mixing rules here. Andreas, July 2016
                        j = j + 1
                        rho_new = rho_old - delta_rho/j
                        if (j > 10) then
                            errval = -8888
                            exit
                        end if
                    end do
                    rho_old = rho_new
                    if (errval < 0) exit
                end do
                rhomin = rho_new
                rho_atmin = rhomin

                !if(gl%seacalc) then
                !    continue
                !end if


                !-----------------------------------------------------
                ! Step 3: check if an extremum was found
                !-----------------------------------------------------
                if ((k < 50) .AND. (errval >= 0)) then
                    extrem = .true.
                    Deriv = D2PDD2_CALC(gl,T, rhomin, nrsubst)*1.d6*factorpress_inv
                    ! check if it is really a minimum or a saddle point
                    ! Since it is not possible to hit a saddle point exactly several criterions are checked to determine if
                    ! the point found is close to a saddle point
                    ! 1) First derivative = 0
                    ! 2) Second derivative = 0
                    ! 3) The slope becomes positive left to the point

                    if (Deriv > 0.d0) then  !Deriv > 0 -- >  possible Minimum
                        min_found = .true.
                        if (abs(Deriv) < 1.D-6) then
                            ! check for possible saddle point by looking at the slope left from the point
                            do i = 1, 10
                                rho_new = rho_new*(1.d0-i*0.01d0)
                                slope = DPDD_CALC(gl,T, rho_new, nrsubst)*1.d6*factorpress_inv
                                if (slope > 1.d-8) then
                                    min_found = .false. ! not a minimum but a saddle point
                                    exit
                                end if !
                            end do
                        end if
                    end if

                end if

                ! If a maximum was found instead of a minimum, try searching for the minimum at higher densities again
                if ((extrem) .and. (.not. min_found)) then
                    rhomin = rhomin*1.001
                    rhomax_allowed = rhomax
                    rhomin_allowed = rhomin
                    call Regula_Falsi(gl,dp_drho_zero, rho_root, rhomin, rhomax, Delta_allowed, rhomin_allowed, rhomax_allowed, &
                        & Max_iterations, Iterations, Errorflag, Parameters)
                    ! check if the root is the correct maximum
                    if (rho_root > 0.d0) then
                        rhomin = rho_root
                        Deriv = D2PDD2_CALC(gl,T, rhomin, nrsubst)*1.d6*factorpress_inv
                        ! check if it is really a minimum or a saddle point
                        ! Since it is not possible to hit a saddle point exactly several criterions are checked to determine if
                        ! the point found is close to a saddle point
                        ! 1) First derivative = 0
                        ! 2) Second derivative = 0
                        ! 3) The slope becomes positive left to the point
                        if (Deriv > 0.d0) then  !Deriv > 0 -- >  possible Minimum
                            min_found = .true.
                            if (abs(Deriv) < 1.D-6) then
                                ! check for possible saddle point by looking at the slope left from the point
                                do i = 1, 10
                                    rho_new = rho_new*(1.d0-i*0.01d0)
                                    slope = DPDD_CALC(gl,T, rho_new, nrsubst)*1.d6*factorpress_inv
                                    if (slope > 1.d-8) then
                                        min_found = .false. ! not a minimum but a saddle point
                                        exit
                                    end if !
                                end do
                            end if
                        end if
                    end if
                end if

                !if(gl%seacalc) then
                !    continue
                !end if


                !-----------------------------------------------------
                ! Step 4: check if the density is between the start density and the density at the maximum
                !-----------------------------------------------------
                if (min_found) then
                    p_min = p_CALC(gl,T, rhomin,nrsubst)
                    if (p_min < p) then
                        ! now an interval is found that encloses the density, start regula falsi
                        rhomax_allowed = rhomax
                        rhomin_allowed = rhomin
                        rholiq = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                        !Andreas Dec 2011: If Iphase is 1 then no vapor phase is tried to be found -- >  Set maxfound true, such that the program does not try to find a vapor phase later
                        if (Iphase == 1) max_found = .true.
                    end if
                end if
            end if

        elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then            !SRK is used
            !THIS IF-STATEMENT SEEMS TO NOT MAKE SENSE. IF NRSUBST IS FOR EXAMPLE 1 (PURE COMPONENT) AND THE COMPONENT IS OIL, NOTHING IS DONE.
            !ANYHOW IN BOTH CASES THE SAME CALL TO THE SAME ROUTINE IS MADE. ANDREAS, MAY 2016
            if (nrsubst /= 0) then
                if (gl%components(nrsubst) /= "oil") rholiq = rho_SRK(gl,T, P, 1, nrsubst)
            else
                rholiq = rho_SRK(gl,T, P, 1, nrsubst)
            end if
            !If an estimate for the density is given, check if the density found is near that value!
            !Andreas, Nov. 2012
            if (rho_est_given > 0.d0) then
                rholiqtest = rho_SRK(gl,T, P, 2, nrsubst)
                if (abs(rholiqtest - rholiq) > 1.D-5) then
                    if (abs(rholiqtest - rho_est_given) < abs(rholiq - rho_est_given)) then
                        rholiq = rholiqtest
                    end if
                end if
                rhomix_calc = rholiq
                return
            end if
            !rhomix_calc = rholiq
            !Changes made for using the whole density solver for the SRK as well
            !That means the density solver decides for mixtures which density is the more reasonable one
            !Andreas June 2013
            !return !For testing only

            !Check if the density is a liquid density (d2pdrho2 > 0)
            Deriv = D2PDD2_CALC(gl,T, rholiq, nrsubst)*1.d6*factorpress_inv
            if (Deriv < 0.D0) then
                rholiq = 0.D0
                min_found = .false.
            end if

            ! Stefan Feb 2014
        elseif ((gl%mix_type == 3) .or. (gl%mix_type == 31)) then            !PR is used
            rholiq = rho_PR(gl,T, P, 1, nrsubst)
            !If an estimate for the density is given, check if the density found is near that value!
            if (rho_est_given > 0.d0) then
                rholiqtest = rho_PR(gl,T, P, 2, nrsubst)
                if (abs(rholiqtest - rholiq) > 1.D-5) then
                    if (abs(rholiqtest - rho_est_given) < abs(rholiq - rho_est_given)) then
                        rholiq = rholiqtest
                    end if
                end if
                rhomix_calc = rholiq
                return
            end if
            !rhomix_calc = rholiq
            !Changes made for using the whole density solver for the PR as well
            !That means the density solver decides for mixtures which density is the more reasonable one
            !return !For testing only

            !Check if the density is a liquid density (d2pdrho2 > 0)
            Deriv = D2PDD2_CALC(gl,T, rholiq, nrsubst)*1.d6*factorpress_inv
            if (Deriv < 0.D0) then
                rholiq = 0.D0
                min_found = .false.
            end if

        End if
    End if

    rhovap = 0.d0
    !----------------------------------------------------------------------------------
    if ((IPHASE == 2) .OR. (IPHASE == 0)) then
        !----------------------------------------------------------------------------------
        if ((gl%mix_Type == 1).OR.(gl%mix_Type == 4) .OR.(gl%mix_Type == 6) .OR. (gl%mix_Type == 7) .OR. (gl%mix_type == 11) .OR. (gl%mix_type == 12) .OR. (gl%mix_type == 13) &
            &  .OR. (gl%mix_type == 110) .OR. (gl%mix_type == 111) .OR. (gl%mix_type == 120) .OR. (gl%mix_type == 121)) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012

            !-----------------------------------------------------
            ! Step 1: Try calculating the density using initial values
            !-----------------------------------------------------
            if (rho_est_given > 0.d0) then
                rho_est = rho_est_given
            else
                rho_est = PSRK(gl,T, P, 2, nrsubst)

                ! in case a pure fluid is treated, unreasonable results of the SRK can be corrected, if PHASETYPE is known  RS, 09/2014

                if (nrsubst == 1) then ! pure fluid
                    if ((iphase == 2).and. (T < 0.98*gl%TC(nrsubst)) .and. (gl%mix_Type == 1)) Then  ! Vapour clearly below Tc
                        rhov_est = dv_eq(gl,T, nrsubst)
                        if ((rho_est > 1.02*rhov_est) .and. (rhov_est > 1.d-13)) rho_est = rhov_est/1.1D0  ! Upper limit for Regula Falsi corresponds to saturated vapour density of anc. eq.
                    end if

                end if

            end if


            ! set the boundaries for the Regula Falsi
            rhomax = 1.10D0*rho_est
            rhomin = 0.90D0*rho_est
            rhomax_allowed = 1.5D0*rho_est
            rhomin_allowed = 0.7D0*rho_est


            if (isnan(rho_est)) then
                rho_est = 1.d-4
                rhomax = 1.10D0*rho_est
                rhomin = 0.90D0*rho_est
                rhomax_allowed = 1.D6*rho_est
                rhomin_allowed = 1.D-6 *rho_est
            else if (nrsubst /= 0) then
                if (rho_est > gl%rhoc(nrsubst)) then
                    rho_est = gl%rhoc(nrsubst)*temperature/gl%tc(nrsubst)
                    rhomax = 1.10D0*rho_est
                    rhomin = 0.90D0*rho_est
                    rhomax_allowed = 1.D6*rho_est
                    rhomin_allowed = 1.D-6 *rho_est
                end if
            end if


            ! call the iteration routine for rhomix(T,P)
            rhovap = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
            ! check if the result s physically reasonable
            if (rhovap > 0.d0) then
                call rho_test (gl,T, p, rhovap, 2, IFound,nrsubst)
                if (IFound < 1) then
                    rhovap = 0.d0
                else
                    if (quick_liq) then
                        min_found = .true.
                        max_found = .true.
                    end if
                    if (min_found) max_found = .true.
                end if
            end if

            !-----------------------------------------------------
            ! Step 2: if Step 1 failed, find the maximum with the lowest density, using Newton
            !-----------------------------------------------------
            if (rhovap <= 1.d-14) then
                rhomin = 1.d-9!*factorpress_inv
                rho_old = rhomin
                rho_new = rho_old
                errval = 0
                rhomax = 0.d0
                extrem = .false.

                !estimation of critical density
                if (nrsubst == 0) then
                    rhoc_comp = 0.d0
                    do j = 1, gl%ncomp
                        rhoc_comp = rhoc_comp + gl%molfractions(j)/gl%rhoc(j)
                    end do
                    rhoc_comp = 1.d0/rhoc_comp * 1.5d0
                else
                    rhoc_comp = gl%rhoc(nrsubst)
                end if

                do k = 1, 30
                    ZeroFunction = DPDD_CALC(gl,T, rho_old,nrsubst)*1.d6*factorpress_inv
                    ! Abbruchkrit.
                    if (abs(ZeroFunction) < 1.d-8) exit
                    Deriv = D2PDD2_CALC (gl,T, rho_old,nrsubst)*1.d6*factorpress_inv
                    !Hier sollte vielleicht abgefangen werden, dass die Krümmung mit der neuen Dichte positiv sein kann und man dann in der Flüssigphase ist (also schon zu weit rechts)
                    rho_new = rho_old - Zerofunction/Deriv
                    ! Sicherheitsabfrage
                    j = 1
                    do while ((rho_new < 0.d0) .OR. ((rho_new > rhoc_comp) .and. (gl%mix_type /= 19)))   !Added an exception here for one-fluid mixture model, Andreas July 2016
                        j = j +1
                        rho_new = rho_old - Zerofunction/Deriv/j
                        if (j > 10) then
                            errval = -8888
                            exit
                        end if
                    end do
                    rho_old = rho_new
                    if (errval < 0) exit
                end do
                rhomax = rho_new
                rho_atmax = rhomax

                !if(gl%seacalc) then
                !    continue
                !end if



                !-----------------------------------------------------
                ! Step 3: check if an extremum was found
                !-----------------------------------------------------
                if ((k < 30) .AND. (errval >= 0)) then
                    extrem = .true.
                    Deriv = D2PDD2_CALC(gl,T, rhomax, nrsubst)*1.d6*factorpress_inv
                    ! check if it is really a maximum or a saddle point
                    ! Since it is not possible to hit a saddle point exactly several criterions are checked to determine if
                    ! the point found is close to a saddle point
                    ! 1) First derivative = 0
                    ! 2) Second derivative = 0
                    ! 3) The slope becomes positive to the right of the point
                    if (Deriv < 0.d0) then  !Deriv > 0 -- >  possible Maximum
                        max_found = .true.
                        if (abs(Deriv) < 1.D-6) then
                            ! check for possible saddle point by looking at the slope right from the point
                            do i = 1, 10
                                rho_new = rho_new*(1.d0 + i*0.01d0)
                                slope = DPDD_CALC(gl,T, rho_new, nrsubst)*1.d6*factorpress_inv
                                if (slope > 1.d-8) then
                                    max_found = .false. ! not a maximum but a saddle point
                                    exit
                                end if !
                            end do
                        end if
                    end if
                end if

                ! If a minimum was found instead of a maximum, try searching for the maximum at lower densities again
                if ((extrem) .and. (.not. max_found)) then
                    rhomax = rhomax/1.001d0
                    rhomax_allowed = rhomax
                    rhomin_allowed = rhomin
                    call Regula_Falsi(gl,dp_drho_zero, rho_root, rhomin, rhomax, Delta_allowed, rhomin_allowed, rhomax_allowed, &
                        & Max_iterations, Iterations, Errorflag, Parameters)
                    ! check if the root is the correct maximum
                    if (rho_root > 0.d0) then
                        rhomax = rho_root
                        Deriv = D2PDD2_CALC(gl,T, rhomax, nrsubst)*1.d6*factorpress_inv
                        ! check if it is really a maximum or a saddle point
                        ! Since it is not possible to hit a saddle point exactly several criterions are checked to determine if
                        ! the point found is close to a saddle point
                        ! 1) First derivative = 0
                        ! 2) Second derivative = 0
                        ! 3) The slope becomes positive to the right of the point
                        if (Deriv < 0.d0) then  !Deriv > 0 -- >  possible Maximum
                            max_found = .true.
                            if (abs(Deriv) < 1.D-6) then
                                ! check for possible saddle point by looking at the slope right from the point
                                do i = 1, 10
                                    rho_new = rho_new*(1.d0 + i*0.01d0)
                                    slope = DPDD_CALC(gl,T, rho_new, nrsubst)*1.d6*factorpress_inv
                                    if (slope > 1.d-8) then
                                        max_found = .false. ! not a maximum but a saddle point
                                        exit
                                    end if !
                                end do
                            end if
                        end if
                    end if
                end if

                !if(gl%seacalc) then
                !    continue
                !end if

                !-----------------------------------------------------
                ! Step 4: check if the density is between the start density and the density at the maximum
                !-----------------------------------------------------
                if (max_found) then
                    p_max = p_CALC(gl,T, rhomax, nrsubst)
                    if (quick_liq) min_found = .true.
                    if (p_max > p) then
                        ! now an interval is found that encloses the density, start regula falsi
                        rhomax_allowed = rhomax
                        rhomin_allowed = rhomin
                        rhovap = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                        if (Iphase == 2) min_found = .true.
                    end if
                end if
            end if

        elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then          !SRK is used

            rhovap = rho_SRK(gl,T, P, 2, nrsubst)
            !If an estimate for the density is given, check if the density found is near that value!
            !Andreas, Nov. 2012
            if (rho_est_given > 0.d0) then
                rhovaptest = rho_SRK(gl,T, P, 1, nrsubst)
                if (abs(rhovaptest - rhovap) > 1.D-5) then
                    if (abs(rhovaptest - rho_est_given) < abs(rhovap - rho_est_given)) then
                        rhovap = rhovaptest
                    end if
                end if
                rhomix_calc = rhovap
                return
            end if
            !Changes made for using the whole density solver for the SRK as well
            !That means the density solver decides for mixtures which density is the more reasonable one
            !Andreas June 2013
            !rhomix_calc = rhovap
            !return !For testing only


            rhomix_dummy = 0.d0
            do i = 1,gl%ncomp
                rhomix_dummy = gl%molfractions(i) * 1.d0 / gl%rhored (i) + rhomix_dummy
            enddo
            rhomix_dummy = 1.d0 / rhomix_dummy

            !Check if the density is a vapor density (d2pdrho2 > 0)
            !if density of current Phase less than factor_decide4vap  * rhoredmix -> gas (factor_decide4vap in modules) else curvature = d2P_drho2 (Temp, rho_phase)
            if (rhovap <= gl%factor_decide4vap * rhomix_dummy) then
                Deriv = -1.d0
            else
                Deriv = D2PDD2_CALC(gl,T, rhovap, nrsubst)*1.d6*factorpress_inv
            endif

            if (Deriv > 0.D0) then
                rhovap = 0.D0
                min_found = .false.
            end if

            if((rhovap > 0.D0) .and. (rholiq >  0.D0)) then
                min_found = .true.
                max_found = .true.
            end if

            !Stefan Feb 2014
        elseif ((gl%mix_type == 3) .or. (gl%mix_type == 31)) then          !PR is used

            rhovap = rho_PR(gl,T, P, 2, nrsubst)
            !If an estimate for the density is given, check if the density found is near that value!
            if (rho_est_given > 0.d0) then
                rhovaptest = rho_PR(gl,T, P, 1, nrsubst)
                if (abs(rhovaptest - rhovap) > 1.D-5) then
                    if (abs(rhovaptest - rho_est_given) < abs(rhovap - rho_est_given)) then
                        rhovap = rhovaptest
                    end if
                end if
                rhomix_calc = rhovap
                return
            end if
            !Changes made for using the whole density solver for the SRK as well
            !That means the density solver decides for mixtures which density is the more reasonable one
            !rhomix_calc = rhovap
            !return !For testing only

            !Check if the density is a vapor density (d2pdrho2 > 0)
            rhomix_dummy = 0.d0
            do i = 1,gl%ncomp
                rhomix_dummy = gl%molfractions(i) * 1.d0 / gl%rhored (i) + rhomix_dummy
            enddo
            rhomix_dummy = 1.d0 / rhomix_dummy

            !Check if the density is a vapor density (d2pdrho2 > 0)
            !if density of current Phase less than factor_decide4vap  * rhoredmix -> gas (factor_decide4vap in modules) else curvature = d2P_drho2 (Temp, rho_phase)
            if (rhovap <= gl%factor_decide4vap * rhomix_dummy) then
                Deriv = -1.d0
            else
                Deriv = D2PDD2_CALC(gl,T, rhovap, nrsubst)*1.d6*factorpress_inv
            endif
            if (Deriv > 0.D0) then
                rhovap = 0.D0
                min_found = .false.
            end if

            if((rhovap > 0.D0) .and. (rholiq >  0.D0)) then
                min_found = .true.
                max_found = .true.
            end if

        End if

    End if

    !-----------------------------------------------------
    ! Step 5: Check, if either the density iteration for liquid or vapor failed and check the extrema found
    !-----------------------------------------------------
    !Special case for cubic equations of state. It may happen for cubic equations of state that if the user wants to calculate a liquid density only a vapor density is found
    !and vice versa. In this case, the algorithm went the same way as for Helmholtz EOS and LKP before, meaning the density might be iterated with regula falsi which produced an error
    !However, this whole procedure is unnecessary because for cubics the vapor and liquid solution can be readily obtained by solving the cubic EOS. To handle this here in a more simple way the
    !handling of cubic EOS is separated from the other more complex EOS-
    !Andreas, May 2016
    if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then    !SRK with linear or quadratic mixing rules for b
        If (iphase == 0) then         !User wants to calculate both density solutions.
            !In an extreme case, it might happen that the vapor and the liquid solution are rejected because of wrong curvature
            !Thus, the densities are recalculated here and the curvature check is omitted
            if ((rhovap <= 1.D-15) .and. (rholiq <= 1.D-15)) then
                rhovap = rho_SRK(gl,T, P, 2, nrsubst)
                rholiq = rho_SRK(gl,T, P, 1, nrsubst)
            end if
        elseif (iphase == 1) then     !User wants liquid density
            !It might happen that the liquid solution is rejected because of wrong curvature
            !Thus, the liquid density is recalculated here and the curvature check is omitted
            if (rholiq <= 1.D-15) then
                rholiq = rho_SRK(gl,T, P, 1, nrsubst)
                !Theoretically it should not happen that no density is returned by the routine rho_SRK, because if a liquid density is specified and not found, the found vapor density is returned instead
                !However, for safety it is checked if a density was found and if not, the other density is tried to calculate
                if (rholiq < 1.D-15) then
                    rholiq = rho_SRK(gl,T, P, 2, nrsubst)
                end if
                rhovap = 0.D0
            end if
        elseif (iphase == 2) then     !User wants gas density
            !It might happen that the gas solution is rejected because of wrong curvature
            !Thus, the gas density is recalculated here and the curvature check is omitted
            if (rhovap <= 1.D-15) then
                rhovap = rho_SRK(gl,T, P, 2, nrsubst)
                !Theoretically it should not happen that no density is returned by the routine rho_SRK, because if a gas density is specified and not found, the found liquid density is returned instead
                !However, for safety it is checked if a density was found and if not, the other density is tried to calculate
                if (rhovap < 1.D-15) then
                    rhovap = rho_SRK(gl,T, P, 1, nrsubst)
                end if
                rholiq = 0.D0
            end if
        end if

    else if((gl%mix_type == 3) .or. (gl%mix_type == 31)) then !PR with linear or quadratic mixing rules for b

        If (iphase == 0) then         !User wants to calculate both density solutions.
            !In an extreme case, it might happen that the vapor and the liquid solution are rejected because of wrong curvature
            !Thus, the densities are recalculated here and the curvature check is omitted
            if ((rhovap <= 1.D-15) .and. (rholiq <= 1.D-15)) then
                rhovap = rho_PR(gl,T, P, 2, nrsubst)
                rholiq = rho_PR(gl,T, P, 1, nrsubst)
            end if
        elseif (iphase == 1) then     !User wants liquid density
            !It might happen that the liquid solution is rejected because of wrong curvature
            !Thus, the liquid density is recalculated here and the curvature check is omitted
            if (rholiq <= 1.D-15) then
                rholiq = rho_PR(gl,T, P, 1, nrsubst)
                !Theoretically it should not happen that no density is returned by the routine rho_PR, because if a liquid density is specified and not found, the found vapor density is returned instead
                !However, for safety it is checked if a density was found and if not, the other density is tried to calculate
                if (rholiq < 1.D-15) then
                    rholiq = rho_PR(gl,T, P, 2, nrsubst)
                end if
                rhovap = 0.D0
            end if
        elseif (iphase == 2) then     !User wants gas density
            !It might happen that the gas solution is rejected because of wrong curvature
            !Thus, the gas density is recalculated here and the curvature check is omitted
            if (rhovap <= 1.D-15) then
                rhovap = rho_PR(gl,T, P, 2, nrsubst)
                !Theoretically it should not happen that no density is returned by the routine rho_PR, because if a gas density is specified and not found, the found liquid density is returned instead
                !However, for safety it is checked if a density was found and if not, the other density is tried to calculate
                if (rhovap < 1.D-15) then
                    rhovap = rho_PR(gl,T, P, 1, nrsubst)
                end if
                rholiq = 0.D0
            end if
        end if

    else if ((gl%mix_Type == 1).OR.(gl%mix_Type == 4) .OR.(gl%mix_Type == 6) .OR. (gl%mix_Type == 7) .OR. (gl%mix_type == 11) .OR. (gl%mix_type == 12) .OR. (gl%mix_type == 13)) then

        if ((rhovap <= 0.D0) .OR. (rholiq <= 0.d0)) then
            ! Check, if no extremum was found (supercritical isotherms)
            ! in that case, try regula falsi one more time with the whole density range
            if ((.not. min_found) .and. (.not. max_found) .and. (rhovap <= 0.D0) .and. (rholiq <= 0.d0)) then
                !Andreas April 2014, lower limit for rhomin adjusted
                !rhomin = 1.d-3*factorpress_inv !Moni
                rhomin = 1.d-8*factorpress_inv
                rhomax = maxval(gl%rhomaxfluid) !rhoredmix*4.d0
                if (nrsubst == 0 .and. (dabs(gl%rhoredmix - 1.d0) < 1.d-12)) then
                    rhoc_comp = 0.d0
                    do j = 1, gl%ncomp
                        rhoc_comp = rhoc_comp + gl%molfractions(j)/gl%rhoc(j)
                    end do
                    enh_fac = 4.d0
                    rhomax = 1.d0/rhoc_comp*enh_fac
                end if
                rhomax_allowed = rhomax
                rhomin_allowed = rhomin
                !check if better starting point can be calculated, Monika March 2018
                !Monika/Sebastian: October 2018: only for supercritical states
                if (T > gl%tredmix) then
                    rho_est = PSRK(gl, T, P, 2, nrsubst)
                    if (rho_est > 1.d-12) then
                        rhomin = rho_est * 0.8d0
                        rhomax = rho_est * 1.2d0
                    end if
                endif
                if (nrsubst /= 0) then
                    if ((gl%components(nrsubst) == "oil") .and. (gl%mix_type == 2)) then
                        rhomin = 400.d0
                        rhomin_allowed = 400.d0
                    else if ((gl%components(nrsubst) == "pec7") .and. (gl%mix_type == 1) .and. iphase == 2) then
                        rhomin = 1.d-12
                        rhomin_allowed = 1.d-12
                        rhomax = 200.d0
                        rhomax_allowed = 2000.d0
                    end if
                end if
                rhovap = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                !error catching for PC-SAFT, Monika 2017-11-14
                if (rhovap .lt. 1.d-12) then
                    o = 1
                    enh_fac = 1.d0 ! TE 180126
                    do while (rhovap .lt. 1.d-12)
                        o = o + 1
                        enh_fac = enh_fac * 0.9d0
                        rhomax = rhomax*enh_fac
                        rhomax_allowed = rhomax
                        rhovap = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                        if (o == 10) exit
                    end do
                end if
                ! check again, if still no density is found, the solver fails!
                if (rhovap <= 0.d0) return
            else if (max_found .and. (.not. min_found)) then
                ! in this case, the search for the maximum and rholiq failed. Try to find
                ! the minimum using fixed density steps
                rho_old = gl%rhoredmix*4.d0
                if (gl%mix_type == 19) then
                    !The one-fluid model is non-corresponding states based. Instead of rhoredmix, the covolume of the SRK is used here to define the maximum density
                    !Andreas Jäger, July 2016
                    rho_old = 0.D0
                    do i=1,gl%ncomp
                        rho_old = rho_old + gl%molfractions(i) * 0.08664D0*8.3144589D0* gl%tc(i) / (gl%pc(i)*1.D6/gl%factorpress)
                    end do
                    rho_old = 1.D0/rho_old * 1.2
                end if
                rho_new = rho_old
                j = 1
                do while (rho_new >= rhomax)
                    !Reinstoffe und krit. Pkt. Monika
                    if(((gl%tredmix-t)/gl%Tredmix < 0.02d0) .and. ((gl%tredmix-t)/gl%Tredmix > 0.d0) .and. (gl%ncomp == 1)) then
                        rho_new = rho_old*(1.D0 - j*0.001D0)
                    else
                        rho_new = rho_old*(1.D0 - j*0.02D0)
                    end if
                    Zerofunction = dpdd_CALC(gl,T, rho_new ,nrsubst)*1.d6*factorpress_inv
                    if (Zerofunction < 0.d0) then
                        min_found = .true.
                        exit
                    end if
                    j = j + 1
                end do
                ! if an interval is found where the derivative dp_drho becomes zero, use Regula Falsi to find the root
                rho_root = 0.D0
                if (min_found) then
                    rhomax = rho_new + rho_old*0.02D0
                    rhomin = rho_new
                    rhomax_allowed = rhomax
                    rhomin_allowed = rhomin
                    call Regula_Falsi(gl,dp_drho_zero, rho_root, rhomin, rhomax, Delta_allowed, rhomin_allowed, rhomax_allowed, &
                        & Max_iterations, Iterations, Errorflag, Parameters)
                else
                    return ! no minimum was found, solver fails!
                end if
                if (rho_root > 0.d0) then! check if Regula succeeded
                    p_min = p_CALC(gl,T, rho_root ,nrsubst)
                    if (p_min < p) then
                        ! now an interval is found that encloses the density, start regula falsi (AGAIN!!!)
                        rhomin = rho_root
                        rhomax = rho_old
                        rhomax_allowed = rhomax
                        rhomin_allowed = rhomin
                        rholiq = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                        ! check again, if still no density is found, the solver fails!
                        if (rholiq <= 0.d0) return
                    end if
                end if
            else if ((.not. max_found) .and. (min_found)) then
                ! in this case, the search for the maximum and rhovap failed. Try to find
                ! the maximum using fixed density steps
                !rho_old = rhoredmix*4.d0
                rho_old = maxval(gl%rhomaxfluid)
                if (gl%mix_type == 19) then
                    !The one-fluid model is non-corresponding states based. Instead of rhoredmix, the covolume of the SRK is used here to define the maximum density
                    !Andreas Jäger, July 2016
                    rho_old = 0.D0
                    do i=1,gl%ncomp
                        rho_old = rho_old + gl%molfractions(i) * 0.08664D0*8.3144589D0* gl%tc(i) / (gl%pc(i)*1.D6/gl%factorpress)
                    end do
                    rho_old = 1.D0/rho_old * 1.2
                end if
                rho_new = 0.D0
                j = 1
                do while (rho_new <= rho_old)
                    rho_new = rho_old*(j*0.02D0)
                    Zerofunction = dpdd_CALC(gl,T, rho_new,nrsubst)*1.d6*factorpress_inv
                    if (Zerofunction < 0.d0) then
                        max_found = .true.
                        exit
                    end if
                    j = j + 1
                end do
                ! if an interval is found where the derivative dp_drho becomes zero, use Regula Falsi to find the root
                rho_root = 0.D0
                if (max_found) then
                    rhomax = rho_new
                    rhomin = rho_new - rho_old*0.02D0
                    if (rhomin <= 0.d0) rhomin = 0.1d0*factorpress_inv  ! changed from if (rhomin < 0.d0) to avoid evaluation at rho == 0 - J.G. 11.2012
                    rhomax_allowed = rhomax
                    rhomin_allowed = rhomin
                    call Regula_Falsi(gl,dp_drho_zero, rho_root, rhomin, rhomax, Delta_allowed, rhomin_allowed, rhomax_allowed, &
                        & Max_iterations, Iterations, Errorflag, Parameters)
                else
                    return ! no minimum was found, solver fails!
                end if
                if (rho_root > 0.d0) then
                    p_max = p_CALC(gl,T, rho_root,nrsubst)
                    if (p_max > p) then
                        ! now an interval is found that encloses the density, start regula falsi (AGAIN!!!)
                        rhomin = 1.D-1*factorpress_inv
                        rhomax = rho_root
                        rhomax_allowed = rhomax
                        rhomin_allowed = rhomin
                        rholiq = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                        ! check again, if still no density is found, the solver fails!
                        if (rholiq <= 0.d0) return
                    end if
                end if

                !   IMPORTANT, if quick liquid density or quick vapor density was found, the pressure at the minimum and/or maximum will not be calculated
                !              thus for these cases the test if pmax < p < pmin is missleading and MUST NOT be used
            else if ((min_found) .and. (max_found) .and. (p_max < p) .and. (p_min > p) .and. &
                & (p_max > 0.D0) .and. (quick_liq .eqv. .false.)) then
                !Andreas, November 2013.
                !In case a minimum and a maximum was found, with p_minimum > p > p_maximum (which means that no valid density solution exists!) accept the wrong solution
                !This composition has to be in the two phase region then

                !Search the root between rhomin and rhomax. THIS IS IN PRINCIPLE NOT A VALID ROOT. BUT IF NO OTHER SOLUTION EXISTS, THIS ONE IS TAKEN
                rhomax = rho_atmin
                rhomin = rho_atmax
                rhomax_allowed = rhomax
                rhomin_allowed = rhomin
                rholiq = rhomix_iter(gl,T, p, rhomin, rhomax, rhomin_allowed, rhomax_allowed,nrsubst)
                !If no root was found -- >  exit
                if (rholiq <= 0.d0) return
                !Check if at least the criterion dpdd > 0 is satisfied. If not, exit with error (a root with dpdd < 0 is not accepted)
                Deriv = dpdd_CALC(gl,T, rholiq ,nrsubst)*factorpress_inv
                if (Deriv < 0.D0) return
                !Check if the density is vapor or liquid like
                Deriv = D2PDD2_CALC(gl,T, rholiq, nrsubst)*1.d6*factorpress_inv
                if (Deriv < 0) then
                    rhovap = rholiq
                    rholiq = 0.D0
                else
                    rhovap = 0.D0
                end if
            end if
        end if

    end if
    !-----------------------------------------------------
    ! Step 6: check if one of the density searches failed, then the other one is taken
    !-----------------------------------------------------
    if ((rhovap <= 0.d0) .and. (rholiq > 0.d0)) then
        rhomix_calc = rholiq
        return
    else if ((rholiq <= 0.d0) .and. (rhovap > 0.d0)) then
        rhomix_calc = rhovap
        return
    else if (rholiq + rhovap <= 0.d0) then
        !Allow negative density solutions if the pressure is negative as well. Needed for evaluation of critical points at negative pressure. Andreas, August 2016
        If (p > 0.D0) then
            rhomix_calc = 0.d0
        else
            if (rholiq < 0.D0) then
                rhomix_calc = rholiq
            else
                rhomix_calc = rhovap
            end if
        end if
        return
    end if

    !-----------------------------------------------------
    ! Step 7: In case of T,p-flash and two density solutions for one phase composition the more likely
    !         solution is the one with the lower Gibbs energy value
    !-----------------------------------------------------
    if ((rhovap == 0.D0) .OR. (rholiq == 0.D0)) return
    if (IPHASE == 0) then
        If ((abs(rholiq - rhovap)/rhovap) > 1.d-5) then
            gibbs_liq = 1.D12
            gibbs_vap = 1.D12
            gibbs_liq_id = 1.D12
            gibbs_vap_id = 1.D12
            rho_save = 0.d0

            !Old version, Andreas June 2013
            !    if (rholiq > 0.D0) gibbs_liq = G_CALC(gl,T,rholiq,nrsubst) ! ganz heikle angelegenheit!!! momentan weiss ich aber nix besseres, klappt meistens....
            !    if (rhovap > 0.D0) gibbs_vap = G_CALC(gl,T,rhovap,nrsubst)
            !    gibbs_diff = gibbs_liq-gibbs_vap
            !
            !    if (gibbs_liq < gibbs_vap) then !the liquid density is the more likely solution
            !        rho_save = rholiq
            !    else if ((gibbs_vap <= gibbs_liq).And.(gibbs_vap < 1.D12)) then !the vapor density is the more likely solution
            !        rho_save = rhovap
            !!    else
            !!        rhomix_calc = 0.D0 ! both density iterations failed
            !    end if


            !The ideal parts of the Gibbs energy at constant composition and temperature cancel out, except for the +ln(del) part
            !Thus not the whole ideal part is needed for the gibbs energy comparison
            !Andreas June 2013
            if (rholiq > 0.D0) gibbs_liq_id = GR_CALC(gl,T,rholiq,nrsubst)
            if (rhovap > 0.D0) gibbs_vap_id = GR_CALC(gl,T,rhovap,nrsubst)
            call R_mix_calc(gl,Rmix)
            do i = 1, gl%ncomp
                if (rholiq > 0.D0) gibbs_liq_id = gibbs_liq_id + gl%molfractions(i) * dlog(rholiq / gl%rhored(i)) * Rmix*T
                if (rhovap > 0.D0) gibbs_vap_id = gibbs_vap_id + gl%molfractions(i) * dlog(rhovap / gl%rhored(i)) * Rmix*T
            end do
            gibbs_diff_id = gibbs_liq_id - gibbs_vap_id

            if (gibbs_diff_id < 0.D0) then !the liquid density is the more likely solution
                rho_save = rholiq
            else if ((gibbs_diff_id >= 0.D0).And.(gibbs_vap_id < 1.D12)) then !the vapor density is the more likely solution
                rho_save = rhovap
                !    else
                !        rhomix_calc = 0.D0 ! both density iterations failed
            end if


            !check the densities for physically wrong solutions
            !This does not make sense for cubic equations and gives misleading results, Andreas June 2013
            if ((gl%mix_type == 1) .or. (gl%mix_type == 11) .or. (gl%mix_type == 12) .or. (gl%mix_type == 13)) then !Helmholtz Eos
                if (rho_save > 0.d0) then
                    if (rho_save < rholiq) then
                        call rho_test (gl,T, p, rho_save, 2, IFound,nrsubst) ! check if rho_save is really a vapor phase density
                        if (IFound == 0) then
                            rhomix_calc = rholiq
                        else
                            rhomix_calc = rho_save
                        end if
                        return
                    else
                        if(gl%seacalc) gl%sea%seap_save = gl%sea%seap
                        call rho_test(gl,T, p, rho_save, 1, IFound,nrsubst) ! check if rho_save is really a liquid phase density
                        if(gl%seacalc) gl%sea%seap = gl%sea%seap_save
                        if (IFound == 0) then
                            rhomix_calc = rhovap
                        else
                            rhomix_calc = rho_save
                        end if
                        return
                    end if
                end if
            end if
            rhomix_calc = rho_save
        else    !If the same density was found twice, it does not matter which density is picked
            rhomix_calc = rhovap
        End if
    end if

    end function rhomix_calc
    !**************************************************************************


    !**************************************************************************
    Double Precision module Function dp_drho_zero(gl,rho, parameters)
    !**************************************************************************



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: Temp  ! Inputvariables
    Double Precision :: rho
    type(type_additional_parameters) :: parameters
    integer::nrsubst
    !  --------------------------------------------------

    Temp = parameters%a_p(1)
    nrsubst = nint(parameters%a_p(2))
    dp_drho_zero = dpdd_CALC(gl,Temp, rho, nrsubst)

    End Function dp_drho_zero
    !**************************************************************************


    !**************************************************************************
    module subroutine rho_test (gl,T, p, rho, IPhase, IFound, nrsubst)
    !**************************************************************************
    ! This subroutine checks, if a density found by the density solver is a physically
    ! correct solution or a solution in a metastable loop.
    !-------------------------------------------------------------------------
    ! Variables:
    !   T - Temperature [K]    - INPUT
    !   p - pressure [MPa]     - INPUT
    !  rho - density [mol/m^3] - INPUT
    ! IPhase - indicates the phase, the density should have
    !       1: liquid phase
    !       2: vapor phase    - INPUT
    !  IFound - indicates the search results:
    !       0: wrong density;
    !       1: correct density - OUTPUT
    !
    ! J. Gernert, Aug. 2011
    !-------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: T, p, rho
    integer:: IPhase, nrsubst
    integer:: IFound

    double precision:: p_old, p_new, rho_new, dpdd1, ddpddd, step, rhoc_comp
    integer:: i, j
    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives


    IFound = 0

    ! has already been checked befroe enrtering this routine!
    !! calculate the derivative d(p)/d(rho)_T
    !dpdd = dpdd_CALC(gl,T,rho, 0)
    !
    !! 1. criterion, valid for all valid phases: d(p)/d(rho)_T must be positive!
    !if (dpdd < 0.d0) return

    !-------------------------------------------------------------------------
    ! check for gas phase density
    if (IPhase == 2) then
        step = 0.1d0/gl%factorpress
        p_old = p
        do i = 0, 4
            rho_new = rho*(1.d0 - i*step)  ! decrease the density

            ! calculate ALL needed derivatives at once
            GETDERR = (/0,1,1,0,0,0,0,1,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
            if (nrsubst == 0) then
                CALL MIXDERIVSFNR(gl,T,rho_new,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
            else
                CALL FNRDERIVS(gl,T,rho_new,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
            end if

            p_new = p_calc_from_der(gl,T, rho_new, nrsubst, FNRDER)!-1.d0/dpdd * rho_new
            ! 2. criterion: the pressure must always be positive
            if (p_new < 0.d0) return
            ! 3. criterion: the gas phase pressure must always decrease with density
            if ((p_new - p_old > 1.d-12) .and. (i > 0)) return
            ! 4. criterion: the slope of an isotherm in the gas phase of a p-rho-plot must always be positive
            dpdd1 = dpdd_CALC_from_der(gl,T,rho_new, nrsubst, FNRDER)
            if (dpdd1 < 0.d0) return
            !Sebastian: For mixture models that have bad mixing rules not all physical cirteria have to be met (curvature)
            ! 5. criterion: the curvature of an isotherm in the gas phase of a p-rho-plot must always be negative
            if (gl%mix_type == 6) gl%bad_mixmodel = .true.    !SAFT
            if (.not. gl%bad_mixmodel) then
                ddpddd = D2PDD2_CALC_from_der(gl, t, rho_new, nrsubst, FNRDER)
                if (ddpddd > 0.d0) return
            endif
            p_old = p_new
        end do
        ! if all criteria are fulfilled in the search area the densityis (likely) a physically correct solution
        IFound = 1

        !-------------------------------------------------------------------------
        ! check for liquid phase density
    else if (IPhase == 1) then
        step = 0.1d0/gl%factorpress
        p_old = p

        !if (gl%seacalc) then
        !    step = 0.05d0 /gl%factorpress
        !end if

        do i = 0, 4

            rho_new = rho*(1.d0 + i*step)  ! increase the density
            ! If a cubic equation is used in the multi-fluid mixture model, it has to be checked here that the density does not exceed the maximum density allowed for cubic equations of state (the inverse of the covolume b)
            ! Andreas Jäger, May 2018
            !--
            do j = 1, gl%ncomp
                if (gl%Eq_type(j) == 2) then   !SRK
                    if (rho_new > (1.D0 / gl%bi_SRK(j))) then
                        !Adjust the stepsize such that the density never gets larger than the inverse of the covolume
                        rho_new = (1.D0 / gl%bi_SRK(j)) * 0.95D0   !Set rho to 95 % of the inverse covolume
                        !Preliminary solution: Exit if this happens once and take solution as correct
                        IFound = 1
                        return
                    end if
                elseif (gl%Eq_type(j) == 3) then !Peng-Robinson
                    if (rho_new > (1.D0 / gl%bi_PR(j))) then
                        !Adjust the stepsize such that the density never gets larger than the inverse of the covolume
                        rho_new = (1.D0 / gl%bi_PR(j)) * 0.95D0   !Set rho to 95 % of the inverse covolume
                        !Preliminary solution: Exit if this happens once and take solution as correct
                        IFound = 1
                        return
                    end if
                end if
            end do
            !--

            ! calculate ALL needed derivatives at once
            GETDERR = (/0,1,1,0,0,0,0,1,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
            if (nrsubst == 0) then
                CALL MIXDERIVSFNR(gl,T,rho_new,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
            else
                CALL FNRDERIVS(gl,T,rho_new,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
            end if

            ! 2. criterion: the pressure must always be positive
            p_new = p_calc_from_der(gl,T, rho_new, nrsubst,FNRDER)!-1.d0/dpdd * rho_new
            if (p_new < 0.d0) return
            !    gl%sea%seap = p_new
            !end if
            if ( (p_new < 0.d0) .and. (gl%seacalc)) then
                continue
            elseif((p_new < 0.d0)) then
                return
            end if
            ! 3. criterion: the liquid phase pressure must always increase with density
            if ((p_new - p_old < 0.d0) .and. (i > 0)) return
            ! 4. criterion: the slope of an isotherm in the liquid phase of a p-rho-plot must always be positive
            dpdd1 = dpdd_CALC_from_der(gl,T,rho_new, nrsubst, FNRDER)
            if (dpdd1 < 0.d0) return
            !Sebastian: Not used anymore
            !THERESA, in case of a lower quality EOS, this criterion does not hold!
            !if (.not.twophasecalc) then ! module parameter, density solution doesn't have to strictly fullfil all physical requirements (curvature)
            !Sebastian: For mixture models that have bad mixing rules not all physical cirteria have to be met (curvature)
            ! 5. criterion: the curvature of an isotherm in the liquid phase of a p-rho-plot must always be positive
            if (gl%mix_type == 6) gl%bad_mixmodel = .true.    !SAFT
            if (.not. gl%bad_mixmodel) then
                ddpddd = D2PDD2_CALC_from_der(gl, t, rho_new, nrsubst, FNRDER)
                if (ddpddd < 0.d0) return
            else if (nrsubst == 0) then
                rhoc_comp = 0.d0
                do j = 1, gl%ncomp
                    rhoc_comp = rhoc_comp + gl%molfractions(j)/gl%rhoc(j)
                end do
                rhoc_comp = 1.d0/rhoc_comp
                if (rho < rhoc_comp) return
            end if
            p_old = p_new
        end do
        if (nrsubst /= 0) then
            if (rho < gl%rhoc(nrsubst)) return
        end if
        ! if all criteria are fulfilled in the search area the density is (likely) a physically correct solution
        IFound = 1
    end if

    endsubroutine rho_test
    !**************************************************************************

    !**************************************************************************
    Double Precision module Function rhomix_iter(gl,T, p, rho_min, rho_max, rho_min_allowed, rho_max_allowed,nrsubst)
    !**************************************************************************
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !---------------------------------------------------------
    Double Precision :: T, p, rho_min, rho_max, rho_min_allowed, rho_max_allowed! Inputvariablen
    Double Precision :: rho_root               ! transmission variable
    Double Precision :: Delta_allowed          ! transmission variable
    integer :: Max_Iterations, Iterations, nrsubst, errorflag    ! transmission variable

    !Integer :: i    !Andreas March 2011

    type(type_additional_parameters) :: parameters    ! needed to transmit T, p to Regula Falsi
    !parameters = 0.d0
    Parameters%a_p(1) = T
    Parameters%a_p(2) = p
    Parameters%a_p(3) = nrsubst
    !Catch NaN
    if ((T /= T) .or. (p /= p)) then
        rhomix_iter = -8888
        return
    endif
    
    Delta_allowed = 1.d-8/gl%factorpress    ! Iteration to 1.d-8 in density was reasonable for specific densities.
    !Delta_allowed = 1.d-11/factorpress    ! for oils
    !Delta_allowed = 1.d-6/factorpress      ! With densities in mol/m3 a more relaxed iteration criterion seems reasonable
    !Return to smaller Delta_allowed value because Validation failed by 1.d-9 % deviation               SH,01/2015
    Max_Iterations = 50

    call Regula_Falsi(gl,pressure_dif, rho_root, rho_min, rho_max, Delta_allowed, &
        rho_min_allowed, rho_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
    !if ((errorflag == 1) .or. (errorflag == 2)) then
    !    rhomix_iter = -8888
    !else
    rhomix_iter = rho_root
    !endif

    if (errorflag == 1 .and. nrsubst /= 0) then
        if (gl%components(nrsubst) == "pec5".or. gl%components(nrsubst) == "pec7") then
            Delta_allowed = 1.d-13
            call Regula_Falsi(gl,pressure_dif, rho_root, rho_min, rho_max, Delta_allowed, &
                rho_min_allowed, rho_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
        end if

    else if (errorflag == 2 .and. nrsubst /= 0) then
        if (gl%components(nrsubst) == "pec5".or. gl%components(nrsubst) == "pec7") then
            rho_max = rho_max*1.5d0
            rho_max_allowed = rho_max_allowed*1.5d0
            call Regula_Falsi(gl,pressure_dif, rho_root, rho_min, rho_max, Delta_allowed, &
                rho_min_allowed, rho_max_allowed, Max_Iterations, Iterations, errorflag, parameters)
        end if
    end if



    rhomix_iter = rho_root

    End Function rhomix_iter
    !**************************************************************************


    !**************************************************************************
    Double Precision module Function pressure_dif(gl,rho, parameters)
    !**************************************************************************



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: rho, T, p! Inputvariablen
    !Double Precision :: P_CALC                  ! function allocating calculated pressure
    type(type_additional_parameters) :: parameters
    integer::nrsubst
    !  --------------------------------------------------

    T = parameters%a_p(1)
    p = parameters%a_p(2)
    nrsubst = int(parameters%a_p(3))
    pressure_dif = p-P_CALC(gl,T, rho, nrsubst)

    End Function pressure_dif
    !**************************************************************************


    end submodule impl
