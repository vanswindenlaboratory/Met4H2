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



    ! module for file ptx-diag.f90
    module ptx_diag_module
    !global use inclusion
    use module_all_types
    use module_isoch_therm
    use flash_pure_module
    use flash_module
    use pc_saft_ARX_derivs_module
    use fnrderivs_module
    use phasedet_mix_module
    use spline_module
    use setup_module
    use fniderivs_module
    use reduced_parameters_calc_module

    contains



    !

    !--------------------------------------------------------------------------------
    !subroutine pxdiag_old(Temp, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, errval)
    !!--------------------------------------------------------------------------------
    !! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    !! TEMPERATURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    !!--------------------------------------------------------------------------------
    !! J. Gernert, March 2012
    !! S.Hielscher und S.Herrig, June 2016, revised Nov 2016



    !
    !implicit none
    !!--------------------------------------------------------------------------------
    !integer, parameter:: maxi = 300 ! maximum number of calculated points, lengt of the return vectors!
    !double precision, intent(in):: Temp     ! specified temperature for the p-x diagram
    !double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    !double precision, intent(out):: x_points(maxi,4)   ! this array is filled according to:
    !                                                    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_vap(2)
    !                                                    ! x_points(i,3) = x_liq(1), x_points(i,4) = x_liq(2)
    !character(255), intent(in):: fileout    ! optional path for output file
    !integer:: points, errval                ! number of calculated points, error code
    !!--------------------------------------------------------------------------------
    !
    !! for calling pure fluid routines
    !double precision:: psat_1, psat_2, psat_high
    !integer::mixflag, errPsat
    !! for calling PhaseDet
    !double precision:: rho(5), x_Phase(30,5), x_known(30)
    !integer::nrofphases, phasetype(5)
    !! for calling Flash_pt and Flash_phase_boundary:
    !double precision::x_vap(30), x_liq(30), beta
    !integer::iter
    !! for calling the spline routines:
    !integer, parameter:: n = 10 ! number of points used for the spline interpolation of the critical point
    !double precision, dimension (n) :: xi(n), yp(n), bp(n), cp(n), dp(n), yx(n), bx(n), cx(n), dx(n)   ! parameters of the spline interpolation
    !double precision:: ispline  ! spline function
    !! internal variables:
    !double precision::x_old, x_crit(2),pstart, deltap, deltap_max, dpdx, press, d_vap, d_liq, rho1, rho2, rhocrit, delta, T, x_backup, press_backup, deltap_backup, delta_backup
    !double precision, dimension(maxi):: T_pts_right, p_pts_right, rhovap_pts_right, rholiq_pts_right
    !double precision:: x_pts_right(maxi,4)
    !logical::TwoPhaseFound, azeo
    !integer:: iFlash, i, j, k, m, small, large, error, count, lastpoint
    !character(255):: errMSG
    !
    !! output file format
    !1000 format(f10.5, f10.4, f10.6, f10.6, f10.6, f10.6, f10.3, f10.3)
    !open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
    !
    !p_points = 0.d0; T_points = 0.d0; x_points = 0.d0; rhovap_points = 0.d0; rholiq_points = 0.d0;
    !p_pts_right = 0.d0; T_pts_right = 0.d0; x_pts_right = 0.d0; rhovap_pts_right = 0.d0; rholiq_pts_right = 0.d0;
    !x_vap = 0.d0; x_liq = 0.d0; points = 1; errval = 0; errMSG = ''
    !d_vap = 0.d0; d_liq = 0.d0
    !azeo = .false.
    !T = Temp
    !deltap_max = 5.d-1  !maximum pressure step
    !psat_high = 0.d0
    !!--------------------------------------------------------------------------------
    !! STEP 0: CHECK FOR BINARY MIXTURE
    !!--------------------------------------------------------------------------------
    !if (ncomp /= 2) then
    !    errval = -9956
    !    return
    !end if
    !
    !!--------------------------------------------------------------------------------
    !! STEP 1: START VALUES
    !!--------------------------------------------------------------------------------
    !! get the vapor pressures of the pure fluids / if one fluid is supercritical an estimation is used
    !if (T < tc(1)) then
    !    psat_1 = vp_eq(gl,T, 1)
    !else
    !    psat_1 = pc(1)*10.d0**(-2.333333d0*(1.d0+accen(1))*(tc(1)/T-1.d0))
    !end if
    !if (T < tc(2)) then
    !    psat_2 = vp_eq(gl,T, 2)
    !else
    !    psat_2 = pc(2)*10.d0**(-2.333333d0*(1.d0+accen(2))*(tc(2)/T-1.d0))
    !end if
    !molfractions(1) = 0.999d0
    !! use the lower vapor pressure as start value and set the composition accordingly
    !if (psat_1 < psat_2) then
    !    pstart = psat_1
    !    small = 2
    !    large = 1
    !else
    !    pstart = psat_2
    !    molfractions(1) = 1.d0 - molfractions(1)
    !    small = 1
    !    large = 2
    !end if
    !
    !mixflag = large
    !! calculate the pure fluid as the first point
    !x_known(large) = 1.d0
    !x_known(small) = 1.d0 - x_known(large)
    !press = pstart
    !call Flash_Pure_PhaseBoundary(press, t, d_vap, d_liq, 1, errPsat, iter, mixflag)
    !if (d_liq*d_vap*press > 0.d0) then
    !    ! write the results to module variables
    !    T_points(points) = T
    !    p_points(points) = press
    !    x_points(points, 1) = x_known(1)
    !    x_points(points, 2) = x_known(2)
    !    x_points(points, 3) = x_known(1)
    !    x_points(points, 4) = x_known(2)
    !    rhovap_points(points) = d_vap
    !    rholiq_points(points) = d_liq
    !end if
    !
    !
    !press = pstart*1.01d0
    !molfractions(2) = 1.d0 - molfractions(1)
    !x_known = molfractions
    !!call reduced_parameters_calc(gl,300.d0)
    !
    !!--------------------------------------------------------------------------------
    !! STEP 2: FIRST MIXTURE VLE POINT
    !!--------------------------------------------------------------------------------
    !points = points + 1
    !TwoPhaseFound = .false.
    !call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !
    !count = 0
    !if (nrofphases < 2) then
    !    !press = pstart*1.01d0
    !    do while(nrofphases < 2)
    !        press = (press + pstart)/2.d0
    !        count = count + 1
    !        if (count >  5) then
    !            errval = -2222
    !            exit
    !        end if
    !        call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !    end do
    !end if
    !
    !if ((nrofphases < 2) .OR. (errval /= 0)) then
    !    if (phasetype(1) == 3) then
    !        iFlash = 3
    !    else
    !        iFlash = 4
    !    end if
    !    call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, 0.d0, 0.d0, beta, iFlash,&
    !                        &  0, errval, iter)
    !else
    !    x_liq = x_Phase(:,phasetype(2))
    !    x_vap = x_Phase(:,phasetype(1))
    !    rho_liq = rho(phasetype(2))
    !    rho_vap = rho(phasetype(1))
    !end if
    !
    !
    !
    !if (errval == 0) then
    !    ! write the results to module variables
    !    T_points(points) = T
    !    p_points(points) = press
    !    x_points(points, 1) = x_vap(1)
    !    x_points(points, 2) = x_vap(2)
    !    x_points(points, 3) = x_liq(1)
    !    x_points(points, 4) = x_liq(2)
    !    rhovap_points(points) = rho_vap
    !    rholiq_points(points) = rho_liq
    !    TwoPhaseFound = .true.
    !else
    !    p_points = errval
    !    goto 111
    !end if
    !
    !!--------------------------------------------------------------------------------
    !! STEP 3: PRESSURE CONTROLLED STEPS !!e.g. methane + ethane
    !!--------------------------------------------------------------------------------
    !x_old = (x_vap(small) + x_liq(small))/2.d0
    !dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
    !deltap = press*0.1d0     !abs(psat_1 - psat_2)/50.d0
    !if (deltap > deltap_max) deltap = deltap_max    !maximum pressure step
    !
    !do while (points < maxi-n)
    !    !create backup of xliq(left component)
    !    x_backup = x_points(points,small+2)
    !    deltap_backup = deltap
    !    press_backup = press
    !    press = press + deltap
    !
    !    x_known(small) = (press - pstart)/dpdx !look for new composition in the 2-phase-region
    !
    !    do while (x_known(small) >= 0.999d0) !in vicinity of the pure component
    !        press = press - deltap
    !        deltap = deltap/2.d0
    !        press = press + deltap
    !        x_known(small) = (press - pstart)/dpdx
    !        if (deltap < 1.d-5) then
    !            errval = -3333
    !            p_points = errval
    !            goto 111
    !        end if
    !    end do
    !    if (x_known(small) <= 0.d0) x_known(small) = 0.001d0 !This could happen for very small dpdx (e.g. for azeotropic mixtures)
    !    ! -> small change in pressure then results in a large composition step that might exceed x = 1!
    !    x_known(large) = 1.d0 - x_known(small)
    !    ! calculate the pT-flash
    !    call Flash_pT(press, T, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
    !    nrofphases = 2
    !
    !    !this will prevent large steps in composition -> reduces pressure step
    !    do while ((errval == 0).and.(abs(x_backup-x_liq(small)) > 0.1d0))
    !        press = press_backup
    !        deltap = deltap / 2.d0
    !        press = press + deltap
    !        x_known(small) = (press - pstart)/dpdx
    !        x_known(large) = 1.d0 - x_known(small)
    !        call Flash_pT(press, T, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
    !        nrofphases = 2
    !    end do
    !    !deltap = deltap_backup
    !
    !
    !    if (errval /= 0) then
    !        call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !        if ((nrofphases > 1) .and. (errval == 0)) then
    !            x_liq = x_Phase(:,phasetype(2))
    !            x_vap = x_Phase(:,phasetype(1))
    !            rho_liq = rho(phasetype(2))
    !            rho_vap = rho(phasetype(1))
    !        end if
    !    end if
    !    if ((errval == 0) .AND. (nrofphases > 1)) then
    !        deltap = deltap_backup
    !        TwoPhaseFound = .true.
    !        points = points + 1
    !        ! write the results to module variables
    !        T_points(points) = T
    !        p_points(points) = press
    !        if (rho_vap > rho_liq) then
    !            x_points(points, 1) = x_liq(1)
    !            x_points(points, 2) = x_liq(2)
    !            x_points(points, 3) = x_vap(1)
    !            x_points(points, 4) = x_vap(2)
    !            rhovap_points(points) = rho_liq
    !            rholiq_points(points) = rho_vap
    !        else
    !            x_points(points, 1) = x_vap(1)
    !            x_points(points, 2) = x_vap(2)
    !            x_points(points, 3) = x_liq(1)
    !            x_points(points, 4) = x_liq(2)
    !            rhovap_points(points) = rho_vap
    !            rholiq_points(points) = rho_liq
    !        end if
    !        ! use the mean value of x_vap and x_liq as overall composition of the next step
    !        x_old = (x_vap(small) + x_liq(small))/2.d0
    !        dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
    !        deltap = press*0.1d0
    !        if (deltap > deltap_max) deltap = deltap_max    !maximum pressure step
    !    else if (TwoPhaseFound) then
    !        ! try reducing the pressure step
    !        press = press - deltap
    !        deltap = deltap/2.d0
    !        if (deltap < 1.d-3) exit    ! jump to STEP 5
    !    end if
    !    !check if close to the critical point
    !    rhocrit = (rhovap_points(points) + rholiq_points(points))/2.d0
    !    if ((rhocrit > 0.d0) .AND. ((rhocrit - rhovap_points(points))/rhocrit < 0.1)) exit  ! jump to STEP 5
    !
    !    !--------------------------------------------------------------------------------
    !    ! STEP 4: CHECK IF THE COMPOSITION RANGE WAS COMPLETELY CROSSED
    !    !--------------------------------------------------------------------------------
    !    if (x_known(small) > 0.998) then  ! Das Zweiphasengebiet wurde durchlaufen
    !        x_known(small) = 1.d0
    !        mixflag = small
    !        ! calculate the pure fluid as the last point
    !        x_known(large) = 1.d0 - x_known(small)
    !        call Flash_Pure_PhaseBoundary(press, t, d_vap, d_liq, 1, errPsat, iter, mixflag)
    !        !call VLEpure(T,press,d_vap,d_liq,errPsat, mixflag)
    !        if (d_liq*d_vap*press > 0.d0) then
    !            ! write the results to module variables
    !            T_points(points) = T
    !            p_points(points) = press
    !            x_points(points, 1) = x_known(1)
    !            x_points(points, 2) = x_known(2)
    !            x_points(points, 3) = x_known(1)
    !            x_points(points, 4) = x_known(2)
    !            rhovap_points(points) = d_vap
    !            rholiq_points(points) = d_liq
    !        end if
    !        goto 111
    !    end if
    !    !write(13,1000) p_points(points), T_points(points), x_points(points, 1), x_points(points, 2), x_points(points, 3), &
    !    !            & x_points(points, 4), rholiq_points(points)/1000.d0, rhovap_points(points)/1000.d0
    !end do
    !
    !!----------------------------------------------------------------------------------------------
    !! STEP 5: COMPOSITION CONTROLLED STEPS
    !!----------------------------------------------------------------------------------------------
    !delta = 1.d0
    !do while (points < maxi-n)
    !    ! estimation of the crit. point:
    !    x_crit(small) = (x_vap(small) + x_liq(small))/2.d0
    !    x_crit(large) = 1.d0 - x_crit(small)
    !    if (x_vap(small) < x_liq(small)) azeo = .true. !azeotropic point was already passed
    !    !new step: 1/delta of the distance to the crit. point
    !    if (azeo .eqv. .true.) then
    !                    !composition on bubble line of comp 1
    !        if ((x_points(points,small+2) - x_crit(small)) < 0) then
    !            x_known(small) = x_points(points,small+2) + 0.025d0
    !        else
    !            x_known(small) = x_points(points,small+2) + (x_points(points,small+2) - x_crit(small))/delta
    !        endif
    !        delta_backup = delta
    !        do while ((x_known(small) > 1.d0) .or. (x_known(small) < 0.d0))
    !            delta_backup = 0.5d0 * delta_backup
    !            if (delta_backup < 1) then
    !                x_known(small) = x_points(points,small+2) + (x_points(points,small+2) - x_crit(small))*delta_backup
    !            else
    !                x_known(small) = x_points(points,small+2) + (x_points(points,small+2) - x_crit(small))/delta_backup
    !            endif
    !        end do
    !    else            !composition on dew line of comp 1
    !        x_known(small) = x_points(points,small) + (x_crit(small) - x_points(points,small))/delta
    !    end if
    !    x_known(large) = 1.d0 - x_known(small)
    !    ! calculate a T,x'-Flash
    !    iFlash = 3
    !    x_liq = 0.d0
    !    x_vap = 0.d0
    !    call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_points(points), rholiq_points(points), beta, iFlash,&
    !                        &  0, errval, iter)
    !    if ((errval /= 0).and.(azeo .eqv. .false.)) then
    !        do while (errval /= 0)
    !            delta = delta * 2.d0
    !            x_known(small) = x_points(points,small+2) + (x_crit(small) - x_points(points,small+2))/delta
    !            x_known(large) = 1.d0 - x_known(small)
    !            x_liq = 0.d0
    !            x_vap = 0.d0
    !            call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_points(points), rholiq_points(points), beta, iFlash,&
    !                            &  0, errval, iter)
    !            if (delta > 33.d0) then !try decreasing composition step 5 times. if not successful: Jump to STEP 7
    !                !Get the higher pure vapor pressure
    !                Select case (small)
    !                Case (1)
    !                    if (T < tc(1)) then
    !                        psat_1 = vp_eq(gl,T, 1)
    !                    else
    !                        psat_1 = 0.d0
    !                    endif
    !                    psat_high = psat_1
    !                Case (2)
    !                    if (T < tc(2)) then
    !                        psat_2 = vp_eq(gl,T, 2)
    !                    else
    !                        psat_2 = 0.d0
    !                    endif
    !                    psat_high = psat_2
    !                end select
    !                exit
    !            endif
    !        enddo
    !        delta = 1.d0
    !        if ((errval /= 0).and.(azeo .eqv. .false.)) then
    !            if (x_points(points, 1) - x_points(points, 3) < 5.d-3) errval = 0
    !            exit   ! Jump to STEP 7
    !        endif
    !    endif
    !    points = points + 1
    !    ! write the results to module variables
    !    T_points(points) = T
    !    p_points(points) = press
    !    !if (p_points(points-1)-p_points(points) > 1.e-14) azeo = .true. !azeotropic point was already passed
    !    if (rho_vap > rho_liq) then
    !        x_points(points, 1) = x_liq(1)
    !        x_points(points, 2) = x_liq(2)
    !        x_points(points, 3) = x_vap(1)
    !        x_points(points, 4) = x_vap(2)
    !        rhovap_points(points) = rho_liq
    !        rholiq_points(points) = rho_vap
    !    else
    !        x_points(points, 1) = x_vap(1)
    !        x_points(points, 2) = x_vap(2)
    !        x_points(points, 3) = x_liq(1)
    !        x_points(points, 4) = x_liq(2)
    !        rhovap_points(points) = rho_vap
    !        rholiq_points(points) = rho_liq
    !    end if
    !    !write(*,*)abs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points)
    !    !BREAK CRITERION: Very close to the critical point or no change in the densities
    !    if ((dabs(rho_vap - rho_liq)/rho_liq < 5.d-2) .OR. &
    !        & ((dabs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 1.d-3).and.(azeo .eqv. .false.))) then
    !        !& ((dabs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 1.d-3).and.(x_vap(1) - x_liq(1)  < 5.d-3).and.(azeo == .false.))) then
    !        !& ((x_vap(1)-x_liq(1))/x_liq(1) < 1.d-2) .and. (azeo == .false.)) then
    !        !No change in density could also occur near an azeotrope, e.g. so2 + chlorine with positive azeotrope at 273 K
    !        !This is not checked again, if an azeotrope was already detected. We assume that multiple azeotropes do not excist!
    !        if (azeo .eqv. .false.) then
    !
    !            !Get the higher pure vapor pressure
    !            Select case (small)
    !                Case (1)
    !                    if (T < tc(1)) then
    !                        psat_1 = vp_eq(gl,T, 1)
    !                    else
    !                        psat_1 = 0.d0
    !                    endif
    !                    psat_high = psat_1
    !                Case (2)
    !                    if (T < tc(2)) then
    !                        psat_2 = vp_eq(gl,T, 2)
    !                    else
    !                        psat_2 = 0.d0
    !                    endif
    !                    psat_high = psat_2
    !                end select
    !
    !            !If the current pressure is higher than the higher vapor pressure of the pure components, a positive azeotrope exists
    !            !Here we assume that an azeotrope would not occur, if one component is supercritical... which is hopefully true... In fact, it is not. But it's found to be too exotic
    !            if ((press > psat_high) .AND. (psat_1*psat_2 > 0.d0)) then
    !
    !                azeo = .true. !Azeotrope was detected
    !                !Composition is pushed over the azeotropic one
    !                x_known(small) = x_known(small) + (1.d0 - x_known(small))/(10.d0) !Bigger composition step!
    !                x_known(large) = 1.d0 - x_known(small)
    !                ! calculate a T,x'-Flash
    !                iFlash = 3
    !                call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_points(points), rholiq_points(points), beta, iFlash,&
    !                                    &  0, errval, iter)
    !                if (errval /= 0) exit   ! Jump to STEP 7
    !                points = points + 1
    !                ! write the results to module variables
    !                T_points(points) = T
    !                p_points(points) = press
    !                if (rho_vap > rho_liq) then
    !                    x_points(points, 1) = x_liq(1)
    !                    x_points(points, 2) = x_liq(2)
    !                    x_points(points, 3) = x_vap(1)
    !                    x_points(points, 4) = x_vap(2)
    !                    rhovap_points(points) = rho_liq
    !                    rholiq_points(points) = rho_vap
    !                else
    !                    x_points(points, 1) = x_vap(1)
    !                    x_points(points, 2) = x_vap(2)
    !                    x_points(points, 3) = x_liq(1)
    !                    x_points(points, 4) = x_liq(2)
    !                    rhovap_points(points) = rho_vap
    !                    rholiq_points(points) = rho_liq
    !                end if
    !
    !            !If two real pure vapor pressures do exist, but the current pressure is not higher than the higher pure vapor pressure, a former azetropic mixture has splitted in a left and right VLE plot
    !            !e.g. ethane + co2 at 293 K
    !            elseif ((press < psat_high) .AND. (psat_1*psat_2 > 0.d0)) then
    !                exit
    !            !if relative change in rhovap from previous to current step is small a critical point is approached -> close loop if composition difference too large
    !            !e.g. methane + ethane @ 250 K
    !            elseif ((dabs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 1.d-3) .and. (x_vap(1) - x_liq(1)  > 5.d-3)) then
    !                T_points(points) = T_points(points - 1)
    !                p_points(points) = p_points(points - 1)
    !                x_points(points, 1) = 0.5d0 * (x_points(points - 1, 1) + x_points(points - 1, 3))
    !                x_points(points, 2) = 0.5d0 * (x_points(points - 1, 2) + x_points(points - 1, 4))
    !                x_points(points, 3) = 0.5d0 * (x_points(points - 1, 1) + x_points(points - 1, 3))
    !                x_points(points, 4) = 0.5d0 * (x_points(points - 1, 2) + x_points(points - 1, 4))
    !                rhovap_points(points) = rhovap_points(points - 1)
    !                rholiq_points(points) = rholiq_points(points - 1)
    !                exit
    !            end if
    !
    !        !end if
    !        else
    !            exit !Jump to STEP 7
    !        end if
    !
    !    elseif ((azeo .eqv. .true.).and.(x_liq(small) > 0.996d0)) then
    !        exit
    !    end if
    !
    !end do
    !
    !!--------------------------------------------------------------------------------
    !! STEP 6: CHECK IF AZEOTROPE WAS FOUND AND THE COMPOSITION RANGE WAS COMPLETELY CROSSED
    !!--------------------------------------------------------------------------------
    !if (azeo .eqv. .true.) then
    !    if (x_known(small) > 0.995d0) then  ! Das Zweiphasengebiet wurde durchlaufen
    !        x_known(small) = 1.d0
    !        mixflag = small
    !        ! calculate the pure fluid as the last point
    !        x_known(large) = 1.d0 - x_known(small)
    !        call Flash_Pure_PhaseBoundary(press, t, d_vap, d_liq, 1, errPsat, iter, mixflag)
    !        !call VLEpure(T,press,d_vap,d_liq,errPsat, mixflag)
    !        if ((d_liq*d_vap*press > 0.d0).and.(dabs((press - p_points(points))/p_points(points)) < 1.d-1)) then
    !            ! write the results to module variables
    !            T_points(points) = T
    !            p_points(points) = press
    !            x_points(points, 1) = x_known(1)
    !            x_points(points, 2) = x_known(2)
    !            x_points(points, 3) = x_known(1)
    !            x_points(points, 4) = x_known(2)
    !            rhovap_points(points) = d_vap
    !            rholiq_points(points) = d_liq
    !        end if
    !        goto 111
    !    end if
    !end if
    !
    !
    !!--------------------------------------------------------------------------------
    !! STEP 7: FINAL ESTIMATE OF THE CRITICAL POINT (SPLINE INTERPOLATION)
    !!--------------------------------------------------------------------------------
    !! use the last n/2 points for the generation of the spline parameters
    !! the pressure is used as a function of rholiq and rhovap in order to interpolate
    !! the critical density, so actually n points are used
    !
    !!--------------------------------------------------------------------------------
    !! first, check if the steps above lead close to the critical point. If not,
    !! the two-phase region may not have a critical point and an estimate would be
    !! meaningless.
    !if (dabs(rho_vap - rho_liq)/rho_liq > 1.d-1) then
    !    errMSG = '*** No critical point found at this temperature ***'
    !    !goto 111
    !    goto 110
    !end if
    !!--------------------------------------------------------------------------------
    !
    !k = n/2
    !j = points-k
    !xi(1:k) = rhovap_points(j+1:points)
    !yp(1:k) = p_points(j+1:points)
    !yx(1:k) = x_points(j+1:points,1)
    !do m = 1, k
    !    xi(m+k) = rholiq_points(points+1-m)
    !    yp(m+k) = p_points(points+1-m)
    !    yx(m+k) = x_points(points+1-m,3)
    !end do
    !
    !! generate the spline coefficients
    !call spline (xi, yp, bp, cp, dp, n)
    !! generate the spline coefficients
    !call spline (xi, yx, bx, cx, dx, n)
    !
    !d_vap = xi(k)
    !d_liq = xi(k+1)
    !delta = (d_liq - d_vap)/20.d0
    !! do 10 more density controlled steps toward the critical point
    !do j = 1, 10
    !    rho1 = d_vap + real(j)*delta
    !    rho2 = d_liq - real(j)*delta
    !    press = ispline(gl,rho1, xi, yp, bp, cp, dp, n)
    !    x_vap(1) = ispline(gl,rho1, xi, yx, bx, cx, dx, n)
    !    x_vap(2) = 1.d0 - x_vap(1)
    !    x_liq(1) = ispline(gl,rho2, xi, yx, bx, cx, dx, n)
    !    x_liq(2) = 1.d0 - x_liq(1)
    !    ! write the results to module variables
    !    points = points + 1
    !    T_points(points) = T
    !    p_points(points) = press
    !    x_points(points, 1) = x_vap(1)
    !    x_points(points, 2) = x_vap(2)
    !    x_points(points, 3) = x_liq(1)
    !    x_points(points, 4) = x_liq(2)
    !    rhovap_points(points) = d_vap + real(j)*delta
    !    rholiq_points(points) = d_liq - real(j)*delta
    !
    !end do
    !
    !!--------------------------------------------------------------------------------
    !! STEP 8: CHECK IF AN AZEOTROPIC MIXTURE HAS SPLITTED INTO A LEFT AND RIGHT VLE-PLOT
    !!--------------------------------------------------------------------------------
    !110 continue
    !
    !!does the check with the change in vapor density make sense? hope the current critera is better
    !!if ( (dabs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 5.d-3) .and. (azeo == .false.) .and. (psat_1*psat_2 > 0.d0) ) then
    !!if (((x_vap(1) - x_liq(1))/x_liq(1) <  1.d-2) .and. (azeo == .false.) .and. (psat_1*psat_2 > 0.d0) .and. (psat_high > 1.d-12)) then
    !if (((x_points(points, 1) - x_points(points, 3))/x_points(points, 3) <  2.d-2) .and. (azeo .eqv. .false.) .and. (psat_1*psat_2 > 0.d0) .and. (psat_high > 1.d-12)) then
    !    lastpoint = points
    !    points = 1
    !    ! calculate the pure fluid as the first point of the right VLE plot
    !    x_known(small) = 1.d0
    !    x_known(large) = 1.d0 - x_known(small)
    !    mixflag = small
    !    pstart = psat_high
    !    d_vap = 0.d0
    !    d_liq = 0.d0
    !    call Flash_Pure_PhaseBoundary(pstart, t, d_vap, d_liq, 1, errPsat, iter, mixflag)
    !    ! write the results to module variables
    !    T_pts_right(points) = T
    !    p_pts_right(points) = pstart
    !    x_pts_right(points, 1) = x_known(1)
    !    x_pts_right(points, 2) = x_known(2)
    !    x_pts_right(points, 3) = x_known(1)
    !    x_pts_right(points, 4) = x_known(2)
    !    rhovap_pts_right(points) = d_vap
    !    rholiq_pts_right(points) = d_liq
    !
    !    continue
    !    !--------------------------------------------------------------------------------
    !    !Step 8a: FIRST MIXTURE VLE POINT FROM RIGHT
    !    !--------------------------------------------------------------------------------
    !    press = pstart*1.01d0
    !    x_known(small) = 0.999d0
    !    x_known(large) = 1.d0 - x_known(small)
    !
    !    points = points + 1
    !
    !    !T-x' flash -> x_known = x'
    !    iFlash = 3
    !    call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_pts_right(points), rholiq_pts_right(points), beta, iFlash,&
    !                        &  0, errval, iter)
    !    T_pts_right(points) = T
    !    p_pts_right(points) = press
    !    x_pts_right(points, 1) = x_liq(1)
    !    x_pts_right(points, 2) = x_liq(2)
    !    x_pts_right(points, 3) = x_vap(1)
    !    x_pts_right(points, 4) = x_vap(2)
    !
    !    !Start of new section!
    !
    !    !--------------------------------------------------------------------------------
    !    ! STEP 8b: PRESSURE CONTROLLED STEPS - FROM THE RIGHT TO THE LEFT
    !    !--------------------------------------------------------------------------------
    !    x_old = (x_vap(large) + x_liq(large))/2.d0 ! now the composition of the component on the right starts with small values ... very confusing
    !    dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
    !    deltap = press*0.025d0     !abs(psat_1 - psat_2)/50.d0
    !    if (deltap > deltap_max) deltap = deltap_max    !maximum pressure step
    !
    !    do while (points < maxi-n-lastpoint)
    !        !x_backup = x_points(points,large+2)
    !        !deltap_backup = deltap
    !        !press_backup = press
    !        press = press + deltap
    !
    !        x_known(large) = (press - pstart)/dpdx !look for new composition in the 2-phase-region
    !
    !        !do while (x_known(large) >= 0.999d0) !in vicinity of the pure component
    !        !    press = press - deltap
    !        !    deltap = deltap/2.d0
    !        !    press = press + deltap
    !        !    x_known(large) = (press - pstart)/dpdx
    !        !    if (deltap < 1.d-5) then
    !        !        errval = -3333
    !        !        p_points = errval
    !        !        goto 111
    !        !    end if
    !        !end do
    !        !if (x_known(large) <= 0.d0) x_known(large) = 0.001d0 !This could happen for very small dpdx (e.g. for azeotropic mixtures)
    !        ! -> small change in pressure then results in a large composition step that might exceed x = 1!
    !        x_known(small) = 1.d0 - x_known(large)
    !        ! calculate the pT-flash
    !        call Flash_pT(press, T, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
    !        nrofphases = 2
    !
    !        !do while ((errval == 0).and.(abs(x_backup-x_liq(large)) > 0.1d0))
    !        !    press = press_backup
    !        !    deltap = deltap / 2.d0
    !        !    press = press + deltap
    !        !    x_known(large) = (press - pstart)/dpdx
    !        !    x_known(small) = 1.d0 - x_known(large)
    !        !    call Flash_pT(press, T, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
    !        !    nrofphases = 2
    !        !end do
    !        !deltap = deltap_backup
    !
    !
    !        if (errval /= 0) then
    !            call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !            if ((nrofphases > 1) .and. (errval == 0)) then
    !                x_liq = x_Phase(:,phasetype(2))
    !                x_vap = x_Phase(:,phasetype(1))
    !                rho_liq = rho(phasetype(2))
    !                rho_vap = rho(phasetype(1))
    !            end if
    !        end if
    !        if ((errval == 0) .AND. (nrofphases > 1)) then
    !            !deltap = deltap_backup
    !            TwoPhaseFound = .true.
    !            points = points + 1
    !            ! write the results to module variables
    !            T_pts_right(points) = T
    !            p_pts_right(points) = press
    !            if (rho_vap > rho_liq) then
    !                x_pts_right(points, 1) = x_liq(1)
    !                x_pts_right(points, 2) = x_liq(2)
    !                x_pts_right(points, 3) = x_vap(1)
    !                x_pts_right(points, 4) = x_vap(2)
    !                rhovap_pts_right(points) = rho_liq
    !                rholiq_pts_right(points) = rho_vap
    !            else
    !                x_pts_right(points, 1) = x_vap(1)
    !                x_pts_right(points, 2) = x_vap(2)
    !                x_pts_right(points, 3) = x_liq(1)
    !                x_pts_right(points, 4) = x_liq(2)
    !                rhovap_pts_right(points) = rho_vap
    !                rholiq_pts_right(points) = rho_liq
    !            end if
    !            ! use the mean value of x_vap and x_liq as overall composition of the next step
    !            x_old = (x_vap(large) + x_liq(large))/2.d0
    !            dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
    !            deltap = press*0.025d0
    !            if (deltap > deltap_max) deltap = deltap_max    !maximum pressure step
    !        else if (TwoPhaseFound) then
    !            ! try reducing the pressure step
    !            press = press - deltap
    !            deltap = deltap/2.d0
    !            if (deltap < 5.d-3) exit    ! jump to STEP 5
    !        end if
    !        !check if close to the critical point
    !        rhocrit = (rhovap_pts_right(points) + rholiq_pts_right(points))/2.d0
    !        if ((rhocrit > 0.d0) .AND. ((rhocrit - rhovap_pts_right(points))/rhocrit < 0.1)) exit  ! jump to STEP 5
    !    enddo
    !
    !    !----------------------------------------------------------------------------------------------
    !    ! STEP 8c: COMPOSITION CONTROLLED STEPS
    !    !----------------------------------------------------------------------------------------------
    !    delta = 1.d0
    !    do while (points < maxi-n-lastpoint)
    !        ! estimation of the crit. point:
    !        x_crit(large) = (x_vap(large) + x_liq(large))/2.d0
    !        x_crit(small) = 1.d0 - x_crit(large)
    !        !if (x_vap(large) < x_liq(large)) azeo = .true. !azeotropic point was already passed
    !        !new step: 1/delta of the distance to the crit. point
    !        !if (azeo == .true.) then
    !        !                !composition on bubble line of comp 1
    !        !    if ((x_points(points,large+2) - x_crit(large)) < 0) then
    !        !        x_known(large) = x_points(points,large+2) + 0.025d0
    !        !    else
    !                !why did we choose the dewline for the next step?
    !                !x_known(small) = x_pts_right(points,small+2) - (x_pts_right(points,small+2) - x_crit(small))/delta
    !                !possible fix?
    !                x_known(small) = x_pts_right(points,small) - (x_pts_right(points,small+2) - x_crit(small))/delta
    !        !    endif
    !        !    delta_backup = delta
    !        !    do while ((x_known(large) > 1.d0) .or. (x_known(large) < 0.d0))
    !        !        delta_backup = 0.5d0 * delta_backup
    !        !        if (delta_backup < 1) then
    !        !            x_known(large) = x_points(points,large+2) + (x_points(points,large+2) - x_crit(large))*delta_backup
    !        !        else
    !        !            x_known(large) = x_points(points,large+2) + (x_points(points,large+2) - x_crit(large))/delta_backup
    !        !        endif
    !        !    end do
    !        !else            !composition on dew line of comp 1
    !        !    x_known(large) = x_points(points,large) + (x_crit(large) - x_points(points,large))/delta
    !        !end if
    !        x_known(large) = 1.d0 - x_known(small)
    !        ! calculate a T,x'-Flash
    !        iFlash = 3
    !        x_liq = 0.d0
    !        x_vap = 0.d0
    !        call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_pts_right(points), rholiq_pts_right(points), beta, iFlash,&
    !                            &  0, errval, iter)
    !        do while (errval /= 0)
    !            delta = delta * 2.d0
    !            x_known(small) = x_pts_right(points,small+2) - (x_pts_right(points,small+2) - x_crit(small))/delta
    !            x_known(large) = 1.d0 - x_known(small)
    !            call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_pts_right(points), rholiq_pts_right(points), beta, iFlash,&
    !                            &  0, errval, iter)
    !            if (press - p_pts_right(points) < 0) errval = 1
    !        enddo
    !        points = points + 1
    !        ! write the results to module variables
    !        T_pts_right(points) = T
    !        p_pts_right(points) = press
    !        !if (p_points(points-1)-p_points(points) > 1.e-14) azeo = .true. !azeotropic point was already passed
    !        if (rho_vap > rho_liq) then
    !            x_pts_right(points, 1) = x_liq(1)
    !            x_pts_right(points, 2) = x_liq(2)
    !            x_pts_right(points, 3) = x_vap(1)
    !            x_pts_right(points, 4) = x_vap(2)
    !            rhovap_pts_right(points) = rho_liq
    !            rholiq_pts_right(points) = rho_vap
    !        else
    !            x_pts_right(points, 1) = x_vap(1)
    !            x_pts_right(points, 2) = x_vap(2)
    !            x_pts_right(points, 3) = x_liq(1)
    !            x_pts_right(points, 4) = x_liq(2)
    !            rhovap_pts_right(points) = rho_vap
    !            rholiq_pts_right(points) = rho_liq
    !        end if
    !        !write(*,*)abs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points)
    !        !BREAK CRITERION: Very close to the critical point or no change in the densities
    !        !if ((abs(rho_vap - rho_liq)/rho_liq < 5.d-2) .OR. &
    !        !    & ((abs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 5.d-3))) then
    !        if ((x_liq(1) - x_vap(1) < 1.d-5) .or. (delta > 31))  then
    !            !!No change in density could also occur near an azeotrope, e.g. so2 + chlorine with positive azeotrope at 273 K
    !            !!This is not checked again, if an azeotrope was already detected. We assume that multiple azeotropes do not excist!
    !            !if (azeo == .false.) then
    !            !
    !            !    !Get the higher pure vapor pressure
    !            !    Select case (large)
    !            !        Case (1)
    !            !            psat_1 = vp_eq(gl,T, 1)
    !            !            psat_high = psat_1
    !            !        Case (2)
    !            !            psat_2 = vp_eq(gl,T, 2)
    !            !            psat_high = psat_2
    !            !        end select
    !            !
    !            !    !If the current pressure is higher than the higher vapor pressure of the pure components, a positive azeotrope exists
    !            !    !Here we assume that an azeotrope would not occur, if one component is supercritical... which is hopefully true... In fact, it is not. But it's found to be too exotic
    !            !    if ((press > psat_high) .AND. (psat_1*psat_2 > 0.d0)) then
    !            !
    !            !        azeo = .true. !Azeotrope was detected
    !            !        !Composition is pushed over the azeotropic one
    !            !        x_known(large) = x_known(large) + (1.d0 - x_known(large))/(10.d0) !Bigger composition step!
    !            !        x_known(small) = 1.d0 - x_known(large)
    !            !        ! calculate a T,x'-Flash
    !            !        iFlash = 3
    !            !        call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_points(points), rholiq_points(points), beta, iFlash,&
    !            !                            &  0, errval, iter)
    !            !        if (errval /= 0) exit   ! Jump to STEP 7
    !            !        points = points + 1
    !            !        ! write the results to module variables
    !            !        T_points(points) = T
    !            !        p_points(points) = press
    !            !        if (rho_vap > rho_liq) then
    !            !            x_points(points, 1) = x_liq(1)
    !            !            x_points(points, 2) = x_liq(2)
    !            !            x_points(points, 3) = x_vap(1)
    !            !            x_points(points, 4) = x_vap(2)
    !            !            rhovap_points(points) = rho_liq
    !            !            rholiq_points(points) = rho_vap
    !            !        else
    !            !            x_points(points, 1) = x_vap(1)
    !            !            x_points(points, 2) = x_vap(2)
    !            !            x_points(points, 3) = x_liq(1)
    !            !            x_points(points, 4) = x_liq(2)
    !            !            rhovap_points(points) = rho_vap
    !            !            rholiq_points(points) = rho_liq
    !            !        end if
    !            !
    !            !    !If two real pure vapor pressures do exist, but the current pressure is not higher than the higher pure vapor pressure, a former azetropic mixture has splitted in a left and right VLE plot
    !            !    !e.g. ethane + co2 at 293 K
    !            !    elseif ((press < psat_high) .AND. (psat_1*psat_2 > 0.d0)) then
    !            !        exit
    !            !    end if
    !            !
    !            !!end if
    !            !else
    !                T_points(lastpoint + 1) = -1.d0
    !                p_points(lastpoint + 1) = -1.d0
    !                x_points(lastpoint + 1,:) = -1.d0
    !                rhovap_points(lastpoint + 1) = -1.d0
    !                rholiq_points(lastpoint + 1) = -1.d0
    !                j = points
    !                do i = lastpoint+2,lastpoint + points + 1
    !                    T_points(i) = T_pts_right(j)
    !                    p_points(i) = p_pts_right(j)
    !                    x_points(i,:) = x_pts_right(j,:)
    !                    rhovap_points(i) = rhovap_pts_right(j)
    !                    rholiq_points(i) = rholiq_pts_right(j)
    !                    j = j - 1
    !                enddo
    !                exit !Jump to STEP 7
    !            !end if
    !
    !        !elseif ((azeo == .true.).and.(x_liq(large) > 0.996d0)) then
    !        !    exit
    !        end if
    !
    !    end do
    !
    !!end of new section!
    !
    !
    !points = points + lastpoint + 1
    !
    !end if
    !
    !!--------------------------------------------------------------------------------
    !! CHECK IF WRITING THE DATA TO A FILE IS WISHED
    !!--------------------------------------------------------------------------------
    !111 continue
    !if (fileout /= '') then
    !    open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
    !    if (error /= 0) return
    !    write(13,*)'!---------------------------------------------------------------------------------'
    !    write(13,*)'! P-X DIAGRAM DATA FILE FOR THE BINARY SYSTEM'
    !    write(13,*)'! ', trim(gl%components(1)),'-',trim(gl%components(2))
    !    write(13,*)'!'
    !    write(13,'(A20, F8.3)')'! TEMPERATURE [K]: ', T
    !    write(13,*)'!---------------------------------------------------------------------------------'
    !    write(13,*)'!'
    !    write(13,*)'! p/MPa    ','T/K       ','x',trim(gl%components(1)),'    x',trim(gl%components(2)),'    y',trim(gl%components(1)), &
    !            &'    y',trim(gl%components(2)),'   dl/mol/l ', 'dv/mol/l'
    !    do i = 1, points
    !        write(13,1000) p_points(i), T_points(i), x_points(i, 3), x_points(i, 4), x_points(i, 1), x_points(i, 2), &
    !                & rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
    !    end do
    !    if (errMSG /= '') write(13,*)'! ', errMSG
    !    close(13)
    !end if
    !
    !end subroutine pxdiag_old



    !--------------------------------------------------------------------------------
    !subroutine txdiag_old(press, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, switchfluids, errval)
    !!--------------------------------------------------------------------------------
    !! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    !! PRESSURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    !!--------------------------------------------------------------------------------
    !! C. Witzke March 2015, based on: J. Gernert, March 2012
    !! S.Hielscher und S.Herrig, June 2016



    !
    !implicit none
    !!--------------------------------------------------------------------------------
    !integer, parameter:: maxi = 300 ! maximum number of calculated points, lengt of the return vectors!
    !double precision, intent(in):: press     ! specified temperature for the p-x diagram
    !double precision, intent(out):: p_points(maxi), T_points(maxi), rhovap_points(maxi), rholiq_points(maxi)    ! return vectors of the calculated points
    !double precision, intent(out):: x_points(maxi,4)   ! this array is filled according to:
    !                                                    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_vap(2)
    !                                                    ! x_points(i,3) = x_liq(1), x_points(i,4) = x_liq(2)
    !character(255), intent(in):: fileout    ! optional path for output file
    !integer:: points, switchfluids, errval                ! number of calculated points, error code
    !!--------------------------------------------------------------------------------
    !
    !! for calling pure fluid routines
    !double precision:: psat_1, psat_2
    !integer::mixflag, errPsat
    !! for calling PhaseDet:
    !double precision:: rho(5), x_Phase(30,5), x_known(30)
    !integer::nrofphases, phasetype(5)
    !! for calling Flash_pt and Flash_phase_boundary:
    !double precision::x_vap(30), x_liq(30), beta
    !integer::iter
    !! for calling the spline routines:
    !integer, parameter:: n = 10 ! number of points used for the spline interpolation of the critical point
    !double precision, dimension (n) :: xi(n), yp(n), bp(n), cp(n), dp(n), yx(n), bx(n), cx(n), dx(n)   ! parameters of the spline interpolation
    !double precision:: ispline  ! spline function
    !! internal variables:
    !double precision::x_old, x_crit(2),pstart, deltap, deltap_max, dpdx, d_vap, d_liq, rho1, rho2, rhocrit, delta, T
    !logical::TwoPhaseFound
    !integer:: iFlash, i, j, k, m, small, large, error, count
    !character(255):: errMSG
    !
    !!Witzke
    !double precision :: deltaT, deltaT_max, dTdx, p, prop2, T_EOS, Tboil_1, Tboil_2, Tstart, Temp,Tsave
    !character(255) :: fluids, moles, path, EOS_indicator, input
    !
    !! output file format
    !1000 format(f10.5, f10.4, f10.6, f10.6, f10.6, f10.6, f10.3, f10.3)
    !open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
    !
    !p_points = 0.d0; T_points = 0.d0; x_points = 0.d0; rhovap_points = 0.d0; rholiq_points = 0.d0;
    !x_vap = 0.d0; x_liq = 0.d0; points = 0; errval = 0; errMSG = ''
    !p = press
    !deltaT_max = 1.d0  !maximum temperature step
    !
    !!--------------------------------------------------------------------------------
    !! STEP 0: CHECK FOR BINARY MIXTURE
    !!--------------------------------------------------------------------------------
    !if (ncomp /= 2) then
    !    errval = -9956
    !    return
    !end if
    !
    !!--------------------------------------------------------------------------------
    !! STEP 1: START VALUES
    !!--------------------------------------------------------------------------------
    !! get the boiling temperatures of the pure fluids
    !
    !    input = 'pliq'
    !    prop2 = 0.d0
    !
    !    fluids = 'ethane'
    !    moles = '1.d0'
    !    EOS_indicator = '1;1'
    !    path = 'Z:\Aktuelle_Fluid_Dateien\GERG-2008\'
    !    Tboil_1 = T_EOS (gl,input, p, prop2, fluids, moles, EOS_indicator, path)
    !
    !
    !    fluids = 'propane'
    !    Tboil_2 = T_EOS (gl,input, p, prop2, fluids, moles, EOS_indicator, path)
    !
    !    fluids = 'ethane;propane'
    !    EOS_indicator = '1;1;1'
    !    moles = '0.5d0;0.5d0'
    !    Temp = T_EOS (gl,input, p, prop2, fluids, moles, EOS_indicator, path)
    !    Temp = 0.d0
    !
    !!if (T < tc(1)) then
    !!    psat_1 = vp_eq(gl,T, 1)
    !!else
    !!    psat_1 = pc(1)*10.d0**(-2.333333d0*(1.d0+accen(1))*(tc(1)/T-1.d0))
    !!end if
    !!if (T < tc(2)) then
    !!    psat_2 = vp_eq(gl,T, 2)
    !!else
    !!    psat_2 = pc(2)*10.d0**(-2.333333d0*(1.d0+accen(2))*(tc(2)/T-1.d0))
    !!end if
    !
    !molfractions(1) = 0.999d0
    !
    !! use the higher boiling temperature as start value and set the composition accordingly
    !if (Tboil_1 < Tboil_2) then
    !    Tstart = Tboil_2
    !    small = 1
    !    large = 2
    !    molfractions(1) = 1.d0 - molfractions(1)
    !    Tsave =  Tboil_2-Tboil_1
    !else
    !    Tstart = Tboil_1
    !    small = 2
    !    large = 1
    !    Tsave =  Tboil_1-Tboil_2
    !end if
    !
    !
    !
    !Temp = Tstart*0.99d0
    !molfractions(2) = 1.d0 - molfractions(1)
    !x_known = molfractions
    !call reduced_parameters_calc(gl,300.d0)
    !
    !!--------------------------------------------------------------------------------
    !! STEP 2: FIRST VLE POINT
    !!--------------------------------------------------------------------------------
    !points = 1
    !TwoPhaseFound = .false.
    !call PhaseDet(p, Temp, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !
    !count = 0
    !if (nrofphases < 2) then
    !    do while(nrofphases < 2)
    !        Temp = (Temp + Tstart)/2.d0
    !        count = count + 1
    !        call PhaseDet(p, Temp, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !        !if (count >  5) then
    !        !    errval = -2222
    !        !    exit
    !        !end if
    !    end do
    !end if
    !
    !
    !!if ((nrofphases < 2) .OR. (errval /= 0)) then
    !!    if (phasetype(1) == 3) then
    !!        iFlash = 3
    !!    else
    !!        iFlash = 4
    !!    end if
    !!    call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, 0.d0, 0.d0, beta, iFlash,&
    !!                        &  0, errval, iter)
    !!else
    !    x_liq = x_Phase(:,phasetype(2))
    !    x_vap = x_Phase(:,phasetype(1))
    !    rho_liq = rho(phasetype(2))
    !    rho_vap = rho(phasetype(1))
    !!end if
    !
    !
    !if (errval == 0) then
    !    ! write the results to module variables
    !    T_points(points) = Temp
    !    p_points(points) = p
    !    x_points(points, 1) = x_vap(1)
    !    x_points(points, 2) = x_vap(2)
    !    x_points(points, 3) = x_liq(1)
    !    x_points(points, 4) = x_liq(2)
    !    rhovap_points(points) = rho_vap
    !    rholiq_points(points) = rho_liq
    !    TwoPhaseFound = .true.
    !else
    !    write(errMSG,*) errval
    !    errMSG = trim(errMSG)//'  *** Calculation of the first point failed ***'
    !    goto 111
    !end if
    !write(13,1000) p_points(points), T_points(points), x_points(points, 1), x_points(points, 2), x_points(points, 3), &
    !                & x_points(points, 4), rholiq_points(points)/1000.d0, rhovap_points(points)/1000.d0
    !
    !!--------------------------------------------------------------------------------
    !! STEP 3: TEMPERATURE CONTROLLED STEPS
    !!--------------------------------------------------------------------------------
    !x_old = (x_vap(small) + x_liq(small))/2.d0
    !dTdx = (Temp - Tstart)/x_old    !calculate the slope of the temperature with the overall x
    !deltaT = Temp*0.1d0
    !if (deltaT > deltaT_max) deltaT = deltaT_max    !maximum temperature step
    !
    !do while (points < maxi-n)
    !    Temp = Temp - deltaT
    !    x_known(small) = (Temp - Tstart)/dTdx
    !    !do while (x_known(small) >= 0.999d0)
    !    !    press = press - deltap
    !    !    deltap = deltap/2.d0
    !    !    press = press + deltap
    !    !    x_known(small) = (press - pstart)/dpdx
    !    !    if (deltap < 1.d-5) then
    !    !        errval = -3333
    !    !        write(errMSG,*) errval
    !    !        errMSG = trim(errMSG)//'  *** Calculation of the second point failed ***'
    !    !
    !    !        goto 111
    !    !    end if
    !    !end do
    !    !if (x_known(small) <= 0.d0) x_known(small) = 0.001d0
    !    x_known(large) = 1.d0 - x_known(small)
    !    ! calculate the pT-flash
    !    call Flash_pT(p, Temp, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
    !    nrofphases = 2
    !    !if (errval /= 0) then
    !    !    call PhaseDet(p, Temp, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
    !    !    if ((nrofphases > 1) .and. (errval == 0)) then
    !    !        x_liq = x_Phase(:,phasetype(2))
    !    !        x_vap = x_Phase(:,phasetype(1))
    !    !        rho_liq = rho(phasetype(2))
    !    !        rho_vap = rho(phasetype(1))
    !    !    end if
    !    !end if
    !    if ((errval == 0) .AND. (nrofphases > 1)) then
    !        TwoPhaseFound = .true.
    !        points = points + 1
    !        ! write the results to module variables
    !        T_points(points) = Temp
    !        p_points(points) = p
    !        if (rho_vap > rho_liq) then
    !            x_points(points, 1) = x_liq(1)
    !            x_points(points, 2) = x_liq(2)
    !            x_points(points, 3) = x_vap(1)
    !            x_points(points, 4) = x_vap(2)
    !            rhovap_points(points) = rho_liq
    !            rholiq_points(points) = rho_vap
    !        else
    !            x_points(points, 1) = x_vap(1)
    !            x_points(points, 2) = x_vap(2)
    !            x_points(points, 3) = x_liq(1)
    !            x_points(points, 4) = x_liq(2)
    !            rhovap_points(points) = rho_vap
    !            rholiq_points(points) = rho_liq
    !        end if
    !        !use the mean value of x_vap and x_liq as overall composition of the next step
    !        x_old = (x_vap(small) + x_liq(small))/2.d0
    !        dTdx = (Temp - Tstart)/x_old    !calculate the slope of the temperature with the overall x
    !        deltaT = Temp*0.1d0
    !        if (deltaT > deltaT_max) deltaT = deltaT_max    !maximum pressure step
    !        else if (TwoPhaseFound) then
    !        ! try reducing the temperature step
    !        Temp = Temp - deltaT
    !        deltaT = deltaT/2.d0
    !        !if (deltap < 1.d-5) exit    ! jump to STEP 5
    !    end if
    !    !check if close to the critical point
    !    !rhocrit = (rhovap_points(points) + rholiq_points(points))/2.d0
    !    !if ((rhocrit > 0.d0) .AND. ((rhocrit - rhovap_points(points))/rhocrit < 0.1)) exit  ! jump to STEP 5
    !
    !    !--------------------------------------------------------------------------------
    !    ! STEP 4: CHECK IF THE COMPOSITION RANGE WAS COMPLETELY CROSSED
    !    !--------------------------------------------------------------------------------
    !    if (x_known(small) > 0.998) then  ! Das Zweiphasengebiet wurde durchlaufen
    !        x_known(small) = 1.d0
    !        mixflag = small
    !        ! calculate the pure fluid as the last point
    !        x_known(large) = 1.d0 - x_known(small)
    !        Temp = Tstart-Tsave+0.1
    !        d_vap = 0.d0
    !        d_liq = 0.d0
    !        call Flash_Pure_PhaseBoundary(press, Temp, d_vap, d_liq, 2, errPsat, iter, 1)
    !        !call VLEpure(T,press,d_vap,d_liq,errPsat, mixflag)
    !        if (d_liq*d_vap*Temp > 0.d0) then
    !            ! write the results to module variables
    !            T_points(points) = Temp
    !            p_points(points) = p
    !            x_points(points, 1) = x_known(1)
    !            x_points(points, 2) = x_known(2)
    !            x_points(points, 3) = x_known(1)
    !            x_points(points, 4) = x_known(2)
    !            rhovap_points(points) = d_vap
    !            rholiq_points(points) = d_liq
    !        end if
    !        goto 111
    !    end if
    !    write(13,1000) p_points(points), T_points(points), x_points(points, 1), x_points(points, 2), x_points(points, 3), &
    !                & x_points(points, 4), rholiq_points(points)/1000.d0, rhovap_points(points)/1000.d0
    !end do
    !
    !!!--------------------------------------------------------------------------------
    !!! STEP 5: COMPOSITION CONTROLLED STEPS
    !!!--------------------------------------------------------------------------------
    !!delta = 10.d0
    !!do while (points < maxi-n)
    !!    ! estimation of the crit. point:
    !!    x_crit(1) = (x_vap(1) + x_liq(1))/2.d0
    !!    x_crit(2) = 1.d0 - x_crit(1)
    !!    !new step: 1/delta of the distance to the crit. point
    !!    x_known(1) = x_points(points,1) + (x_crit(1) - x_points(points,1))/delta
    !!    x_known(2) = 1.d0 - x_known(1)
    !!    ! calculate a T,x'-Flash
    !!    iFlash = 3
    !!    call Flash_PhaseBoundary(press, T, x_known, x_vap, x_liq, rhovap_points(points), rholiq_points(points), beta, iFlash,&
    !!                        &  0, errval, iter)
    !!    if (errval /= 0) exit   ! Jump to STEP 7
    !!    points = points + 1
    !!    ! write the results to module variables
    !!    T_points(points) = T
    !!    p_points(points) = press
    !!    if (rho_vap > rho_liq) then
    !!        x_points(points, 1) = x_liq(1)
    !!        x_points(points, 2) = x_liq(2)
    !!        x_points(points, 3) = x_vap(1)
    !!        x_points(points, 4) = x_vap(2)
    !!        rhovap_points(points) = rho_liq
    !!        rholiq_points(points) = rho_vap
    !!    else
    !!        x_points(points, 1) = x_vap(1)
    !!        x_points(points, 2) = x_vap(2)
    !!        x_points(points, 3) = x_liq(1)
    !!        x_points(points, 4) = x_liq(2)
    !!        rhovap_points(points) = rho_vap
    !!        rholiq_points(points) = rho_liq
    !!    end if
    !!    !--------------------------------------------------------------------------------
    !!    !BREAK CRITERION: Very close to the critical point or no change in the densities
    !!    if ((abs(rho_vap - rho_liq)/rho_liq < 5.d-2) .OR. &
    !!      & (abs(rhovap_points(points) - rhovap_points(points-1))/rhovap_points(points) < 5.d-3)) then
    !!        exit    ! Jump to STEP 7
    !!    end if
    !!    write(13,1000) p_points(points), T_points(points), x_points(points, 1), x_points(points, 2), x_points(points, 3), &
    !!                & x_points(points, 4), rholiq_points(points)/1000.d0, rhovap_points(points)/1000.d0
    !!end do
    !
    !!!--------------------------------------------------------------------------------
    !!! STEP 7: FINAL ESTIMATE OF THE CRITICAL POINT (SPLINE INTERPOLATION)
    !!!--------------------------------------------------------------------------------
    !!! use the last n/2 points for the generation of the spline parameters
    !!! the pressure is used as a function of rholiq and rhovap in order to interpolate
    !!! the critical density, so actually n points are used
    !!
    !!!--------------------------------------------------------------------------------
    !!! first, check if the steps above lead close to the critical point. If not,
    !!! the two-phase region may not have a critical point and an estimate would be
    !!! meaningless.
    !!if (abs(rho_vap - rho_liq)/rho_liq > 1.d-1) then
    !!    errMSG = '*** No critical point found at this temperature ***'
    !!    goto 111
    !!end if
    !!!--------------------------------------------------------------------------------
    !!
    !!k = n/2
    !!j = points-k
    !!xi(1:k) = rhovap_points(j+1:points)
    !!yp(1:k) = p_points(j+1:points)
    !!yx(1:k) = x_points(j+1:points,1)
    !!do m = 1, k
    !!    xi(m+k) = rholiq_points(points+1-m)
    !!    yp(m+k) = p_points(points+1-m)
    !!    yx(m+k) = x_points(points+1-m,3)
    !!end do
    !!
    !!! generate the spline coefficients
    !!call spline (xi, yp, bp, cp, dp, n)
    !!! generate the spline coefficients
    !!call spline (xi, yx, bx, cx, dx, n)
    !!
    !!d_vap = xi(k)
    !!d_liq = xi(k+1)
    !!delta = (d_liq - d_vap)/20.d0
    !!! do 10 more density controlled steps toward the critical point
    !!do j = 1, 10
    !!    rho1 = d_vap + real(j)*delta
    !!    rho2 = d_liq - real(j)*delta
    !!    press = ispline(gl,rho1, xi, yp, bp, cp, dp, n)
    !!    x_vap(1) = ispline(gl,rho1, xi, yx, bx, cx, dx, n)
    !!    x_vap(2) = 1.d0 - x_vap(1)
    !!    x_liq(1) = ispline(gl,rho2, xi, yx, bx, cx, dx, n)
    !!    x_liq(2) = 1.d0 - x_liq(1)
    !!    ! write the results to module variables
    !!    points = points + 1
    !!    T_points(points) = T
    !!    p_points(points) = press
    !!    x_points(points, 1) = x_vap(1)
    !!    x_points(points, 2) = x_vap(2)
    !!    x_points(points, 3) = x_liq(1)
    !!    x_points(points, 4) = x_liq(2)
    !!    rhovap_points(points) = d_vap + real(j)*delta
    !!    rholiq_points(points) = d_liq - real(j)*delta
    !!
    !!end do
    !
    !!--------------------------------------------------------------------------------
    !! CHECK IF WRITING THE DATA TO A FILE IS WISHED
    !!--------------------------------------------------------------------------------
    !111 continue
    !if (fileout /= '') then
    !    open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
    !    if (error /= 0) return
    !    write(13,*)'!---------------------------------------------------------------------------------'
    !    write(13,*)'! T-X DIAGRAM DATA FILE FOR THE BINARY SYSTEM'
    !    write(13,*)'! ', trim(gl%components(1)),'-',trim(gl%components(2))
    !    write(13,*)'!'
    !    write(13,'(A20, F8.3)')'! PRESSURE [MPa]: ', p
    !    write(13,*)'!---------------------------------------------------------------------------------'
    !    write(13,*)'!'
    !    write(13,*)'! p/MPa    ','T/K       ','y',trim(gl%components(1)),'    y',trim(gl%components(2)),'    x',trim(gl%components(1)), &
    !            &'    x',trim(gl%components(2)),'   dl/mol/l ', 'dv/mol/l'
    !    do i = 1, points
    !        write(13,1000) p_points(i), T_points(i), x_points(i, 1), x_points(i, 2), x_points(i, 3), x_points(i, 4), &
    !                & rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
    !    end do
    !    if (errMSG /= '') write(13,*)'! ', errMSG
    !    close(13)
    !end if
    !
    !end subroutine txdiag_old

    !--------------------------------------------------------------------------------
    subroutine pxdiag(gl,Temp, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, switchfluids, errval)
    !--------------------------------------------------------------------------------
    ! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    ! TEMPERATURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    !--------------------------------------------------------------------------------
    ! J. Gernert, March 2012
    ! S.Hielscher und S.Herrig, June 2016, revised Nov 2016





    implicit none

    type(type_gl) :: gl

    !--------------------------------------------------------------------------------
    integer, parameter:: maxi = 300 ! maximum number of calculated points, lengt of the return vectors!
    double precision, intent(in):: Temp     ! specified temperature for the p-x diagram
    double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision, intent(out):: x_points(maxi,4)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_vap(2)
    ! x_points(i,3) = x_liq(1), x_points(i,4) = x_liq(2)
    character(255), intent(in):: fileout    ! optional path for output file
    !switchfluids = 1: fluid order was switched -> this has to be undone by swapping point order and molefractions x1 -> x2, x2 -> x1
    !switchfluids = 2: former azeotropic mixture, splited in two vle regions (co2-ethane T > 290 K)
    !swtichfluids = 3: Mostly relevant for open VLE regions -> start with component with smaller vapor pressure, because this component will most likely be attached to y axis while the other component is detached from y axis
    !switchsupercrit: indicates if the point of minimal x(1) on the dewline for a binary mixtures where one component is supercritical was passed
    integer:: points, switchfluids, switchsupercrit, errval                ! number of calculated points, error code
    !--------------------------------------------------------------------------------

    ! for calling pure fluid routines
    double precision:: psat_1, psat_2, psat_high
    integer::mixflag, errPsat
    ! for calling PhaseDet
    double precision:: rho(5), x_Phase(30,5), x_known(30)
    integer::nrofphases, phasetype(5)
    ! for calling Flash_pt and Flash_phase_boundary:
    double precision::x_vap(30), x_liq(30), beta_lc
    integer::iter, save_calc, force_spline !0: error, don't save 1: save, 2: finish with spline
    ! for calling the spline routines:
    integer, parameter:: n = 10 ! number of points used for the spline interpolation of the critical point
    double precision, dimension (n) :: xi(n), yp(n), bp(n), cp(n), dp(n), yx(n), bx(n), cx(n), dx(n)   ! parameters of the spline interpolation
    ! internal variables:
    double precision::x_old, x_crit(2),pstart, deltap, deltap_max, dpdx(2), dpdx_old(2), press, d_vap, d_liq, rho1, rho2, rhocrit, delta, T, x_backup, press_backup, deltap_backup, delta_backup, x_step
    double precision, dimension(:), allocatable :: T_pts_right, p_pts_right, rhovap_pts_right, rholiq_pts_right
    double precision, dimension(:,:), allocatable :: x_pts_right
    logical::TwoPhaseFound, azeo
    !                                                                   counter for open two phase region
    integer:: iFlash, i, j, k, m, small, large, error, count, lastpoint, open_count, Flash_type, iphase_try, points_copy
    character(255):: errMSG

    ! output file format
1000 format(f10.5, f10.4, f10.6, f10.6, f10.6, f10.6, f10.3, f10.3)


    if (.not. allocated(T_pts_right)) allocate(T_pts_right(maxi))
    if (.not. allocated(p_pts_right)) allocate(p_pts_right(maxi))
    if (.not. allocated(rhovap_pts_right)) allocate(rhovap_pts_right(maxi))
    if (.not. allocated(rholiq_pts_right)) allocate(rholiq_pts_right(maxi))
    if (.not. allocated(x_pts_right)) allocate(x_pts_right(maxi,30))

    p_points = 0.d0; T_points = 0.d0; x_points = 0.d0; rhovap_points = 0.d0; rholiq_points = 0.d0;
    p_pts_right = 0.d0; T_pts_right = 0.d0; x_pts_right = 0.d0; rhovap_pts_right = 0.d0; rholiq_pts_right = 0.d0;
    x_vap = 0.d0; x_liq = 0.d0; errval = 0; errMSG = '';points_copy=0
    ! points = 1 moved to interface because of fromer azeotrop
    d_vap = 0.d0; d_liq = 0.d0
    azeo = .false.
    T = Temp
    deltap_max = 1.d0  !maximum pressure step changed from 0.5d0 to 1.d0 10/19 TN
    psat_high = 0.d0
    switchsupercrit = 0
    force_spline = 0
    save_calc = 0
    dpdx = 0.d0
    dpdx_old = 0.d0
    open_count = 0
    Flash_type = 1 !(1 pT Flash, 2 Tx' or Tx'' Flash)
    !switchfluids = 0 TN 04.19 wieder auskommentiert, da schon in interface_routines intialisiert.

    !in case the calculation starts also from the otherside to make sure to stay below maxi
    if (points>1)then
        points_copy=points
        points = 1
    end if

    !--------------------------------------------------------------------------------
    ! STEP 1: START VALUES
    !--------------------------------------------------------------------------------
    ! get the vapor pressures of the pure fluids
    psat_1 = 0.d0
    psat_2 = 0.d0

    psat_1 = vp_eq(gl,T, 1)
    if (T < gl%tc(2)) then
        psat_2 = vp_eq(gl,T, 2)
    end if
    if ((psat_1 > psat_2) .and. (psat_2 > 1.d-10) .and. (switchfluids == 0) ) then !.and. (psat_1/psat_2) > 1.d2 !Abfrage rausgenommen TN. 19.3.19
        switchfluids = 3
        return
    endif
    ! calculate and save the pure fluid as the first point
    call Flash_Pure_PhaseBoundary(gl,psat_1, t, d_vap, d_liq, 1, errPsat, iter, 1)
    T_points(points) = T
    p_points(points) = psat_1
    pstart = psat_1
    x_points(points, 1) = 1.d0 !xvap1
    x_points(points, 2) = 0.d0 !xvap2
    x_points(points, 3) = 1.d0 !xliq1
    x_points(points, 4) = 0.d0 !xliq2
    rhovap_points(points) = d_vap
    rholiq_points(points) = d_liq
    points = points + 1
    !--------------------------------------------------------------------------------
    ! STEP 2: FIRST MIXTURE VLE POINT
    !--------------------------------------------------------------------------------
    !if (T/gl%tc(1) >= 0.95d0) then
    !    x_step = 5.d-3
    !else
    !    x_step = 1.d-2
    !endif


    if (T/gl%tc(1) >= 0.95d0) then
        x_step = -1.00606d-1 * (T/gl%tc(1)) + 1.00576d-1 !Function to reduce the stepsize depending on the distance to tc at 0.95 x_step~5.d-3
    else
        x_step = 1.d-2
    endif


    x_known(1) = 1.d0
    x_known(2) = 1.d0 - x_known(1)
    iFlash = 3 ! start with T x' Flash
    ! calculate the pure fluid as the first point
    press = pstart
    do while (points < maxi-points_copy-n)
        if (Flash_type ==1)then    !pt based Flash first  !TN
            if (open_count <= 1) then
                save_calc = 0
                x_known(1) = x_known(1) - x_step
                !quit if the whole molesrange was passed
                if ((x_known(1) < 1.d-10) .and. (points>1)) then
                    if ((psat_2 > 1.d-10) .and. (press/psat_2 >0.95)) then

                        ! calculate and save the pure fluid as the first point
                        d_vap = rhovap_points(points-1)
                        d_liq = rholiq_points(points-1)
                        call Flash_Pure_PhaseBoundary(gl,psat_2, t, d_vap, d_liq, 1, errPsat, iter, 2)
                        if (errPsat==0)then
                            T_points(points) = T
                            p_points(points) = psat_2
                            !pstart = psat_1
                            x_points(points, 1) = 0.d0 !xvap1
                            x_points(points, 2) = 1.d0 !xvap2
                            x_points(points, 3) = 0.d0 !xliq1
                            x_points(points, 4) = 1.d0 !xliq2
                            rhovap_points(points) = d_vap
                            rholiq_points(points) = d_liq
                        else
                            points = points - 1
                        endif
                        exit
                    else
                        points = points - 1
                        save_calc = 2
                    endif
                    !exit
                endif
                x_known(2) = 1.d0 - x_known(1)
                x_liq = 0.d0
                x_vap = 0.d0
                !Why not use the compositions of the vapor phase and liquid phase as initial estimates here?
                !Added, if it does not work for some reasons, delete lateron, Andreas Jäger, September 2017
                !----
                if(points > 2) then
                    if (iflash == 3) then
                        x_vap(1) = x_points(points-1, 1)   !xvap1
                        x_vap(2) = x_points(points-1, 2)   !xvap2
                        x_liq(1) = x_known(1)              !xliq1
                        x_liq(2) = x_known(2)              !xliq2
                    elseif (iflash == 4) then
                        x_vap(1) = x_known(1)              !xvap1
                        x_vap(2) = x_known(2)              !xvap2
                        x_liq(1) = x_points(points-1, 3)   !xliq1
                        x_liq(2) = x_points(points-1, 4)   !xliq2
                    endif
                end if
                !----
                call Flash_PhaseBoundary(gl,press, T, x_known, x_vap, x_liq, rhovap_points(points-1), rholiq_points(points-1), beta_lc, iFlash,&
                    &  0, errval, iter)
                if (errval /= 0) then
                    call Flash_PhaseBoundary(gl,press, T, x_known, x_vap, x_liq, 0.d0, 0.d0, beta_lc, iFlash,&
                        &  0, errval, iter)
                endif
                !Catch NaN
                if (press /= press) then
                    switchfluids = 0
                    if (errval == 0) errval = -898964
                    return
                endif
                !everything fine, save variables and continue calculations
                if (errval == 0) save_calc = 1

                !reset a false set of switchsupercrit = 1
                !Example: When the composition of bubble line is the input composition, the composition on the dewline is calculated. If one component is overcritical it could happen that the amount of
                !component two increases on the dewline at lower pressures and decreases at higher pressures. If the flag switch supercrit is one, but the amount of component two still increases the flag is set back to 0
                if ((save_calc == 1) .and. (switchsupercrit == 1)) then
                    !if (iflash == 3) then !T x' Flash
                    if (x_vap(1) - x_points(points - 1, 1)  < 0.d0) then
                        switchsupercrit = 0
                    endif
                    !elseif (iflash == 4) then !T x'' Flash
                    !    !if (x_points(points - 1, 1) - x_points(points, 1) < 0.d0) then
                    !    !    switchsupercrit = 0
                    !    !endif
                    !    continue
                    !endif
                elseif ((save_calc == 1) .and. (iflash == 3) .and. (x_vap(1) - x_points(points-1,1) > 0.d0)) then
                    switchsupercrit = 1
                endif

                !the dewline (or bubble line in rare cases) might change the direction of molar composition from positive to negative direction
                !when this happens switchsupercrit is set to 1 (from 0)
                !if the composition on the dewline was the input composition, the composition on the bubbleline is now the input composition
                !when the pressure now is decreasing, calculations are finished with spline interpolations
                if ((switchsupercrit == 1) .and. (press - p_points(points-1) < 0.d0) .and. (errval == 0)) then
                    save_calc = 2 !quit calculations and go to spline calculations
                    !errval = -10 !arbitrary error value
                endif

                !if the given temperature is >0.95*Tc it (hopefully) does not happen, that the slope of bubble and dewline do not change
                !if this happens, the previous step was wrong and saved values are discarded and calculations are finished with spline interpolations
                if (((T/gl%tc(1) >= 0.95d0).and.(gl%tc(2) < 1.d-12)) .and. (dpdx(1)*dpdx_old(1) < 0.d0).and.(dpdx(2)*dpdx_old(2) < 0.d0)) then
                    save_calc = 2
                    points = points - 1
                endif

                !if right component is supercritical it might happen that the current composition step is out of the two phase region
                !-> try changing form dew to bubble line or vice versa
                if ((errval /= 0) .and. (psat_2 < 1.d-10) .and. (switchsupercrit == 0)) then
                    if (iflash == 3) then ! T x' Flash -> T x'' Flash
                        if ((rholiq_points(points-1)-rhovap_points(points-1))/(rhovap_points(points-1)) < 1.d-1) then
                            save_calc = 2
                        else
                            iflash = 4
                            x_known(1) = x_points(points-1,1) !xvap1
                            x_known(2) = x_points(points-1,2) !xvap2
                        endif
                    elseif (iflash == 4) then! T x'' Flash -> T x' flash
                        iflash = 3
                        x_known(1) = x_points(points-1,3) !xliq1
                        x_known(2) = x_points(points-1,4) !xliq2
                    endif
                    switchsupercrit = switchsupercrit + 1
                    cycle
                elseif ((errval /= 0) .and. (switchsupercrit == 1) .and. (dabs(rhovap_points(points-1) - rholiq_points(points-1))/rholiq_points(points-1) <= 1.d-1)) then
                    save_calc = 2
                elseif ((errval /= 0) .and. (switchsupercrit == 1) .and. (dabs(rhovap_points(points-1) - rholiq_points(points-1))/rholiq_points(points-1) > 1.d-1)) then
                    press = p_points(points-1)
                    x_vap(1) = x_points(points-1, 1)   !xvap1
                    x_vap(2) = x_points(points-1, 2)   !xvap2
                    x_liq(1) = x_points(points-1, 3)   !xliq1
                    x_liq(2) = x_points(points-1, 4)   !xliq2
                    if (points == 2) then
                        open_count = 2
                    else
                        open_count = 3
                    endif
                    cycle


                    !elseif ((dabs(rhovap_points(points-1) - rholiq_points(points-1))/rhovap_points(points-1) < 1.d-1) .and. (temp/minval(gl%tc(1:gl%ncomp)) > 0.85) .and. (errval /= 0))

                elseif ((dabs(gl%rho_vap - gl%rho_liq)/gl%rho_vap < 1.d-1) .and. (temp/minval(gl%tc(1:gl%ncomp)) > 0.85) .and. (errval /= 0)) then   !RB/TN Test ob Temp, nahe kritischer Temperatur bezieht sich jetzt nur auf die Tc der vorhandenen Komponenten und missachtet genullte Werte
                    if (points >= 2) then
                        switchfluids = 2
                        errval = 0
                        save_calc = 2
                    endif
                endif
            elseif (open_count >= 2) then
                save_calc = 0
                !first pressure based vle point
                if (open_count == 2) then
                    TwoPhaseFound = .false.
                    x_known(1) = 0.999d0
                    x_known(2) = 1.d0 - x_known(1)
                    press = 1.01d0 * psat_1
                    call PhaseDet(gl,press, T, x_known, rho, x_Phase, phasetype, beta_lc, nrofphases, errval)

                    count = 0
                    if (nrofphases < 2) then
                        !press = pstart*1.01d0
                        do while(nrofphases < 2)
                            press = (press + pstart)/2.d0
                            count = count + 1
                            if (count >  10) then
                                errval = -2222
                                exit
                            end if
                            call PhaseDet(gl,press, T, x_known, rho, x_Phase, phasetype, beta_lc, nrofphases, errval)
                        end do
                    end if
                    if (errval == 0) then  !2019/03/14 Moni, error handling added
                        x_liq = x_Phase(:,phasetype(2))
                        x_vap = x_Phase(:,phasetype(1))
                        gl%rho_liq = rho(phasetype(2))
                        gl%rho_vap = rho(phasetype(1))
                    end if

                elseif (open_count > 2) then

                    !if (Flash_type == 1)then !follwing pressure based points
                    x_old = (x_vap(2) + x_liq(2))/2.d0
                    dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
                    deltap = press*0.1d0     !abs(psat_1 - psat_2)/50.d0
                    if (deltap > deltap_max) deltap = deltap_max    !maximum pressure step
                    press = press + deltap
                    x_known(2) = (press - pstart)/dpdx(1) !look for new composition in the 2-phase-region
                    !Very likely not happen, delete?
                    do while (x_known(2) >= 0.999d0) !in vicinity of the pure component
                        press = press - deltap
                        deltap = deltap/2.d0
                        press = press + deltap
                        x_known(2) = (press - pstart)/dpdx(1)
                        if (deltap < 1.d-5) then
                            errval = -3333
                            !Andreas Jäger, September 2017. Case occurs for at least for new mixture model (12). Changed that all pressures are set to the error value in order to at least save the previous calculations
                            !p_points = errval
                            exit
                        end if
                    end do
                    !end delete?
                    x_known(1) = 1.d0 - x_known(2)
                    ! calculate the pT-flash
                    call Flash_pT(gl,press, T, x_known, x_vap, x_liq, gl%rho_vap, gl%rho_liq, beta_lc, errval, iter)
                    nrofphases = 2

                    !in case of unlucky x_known values try again
                    if (errval/=0)then
                        do i=1, 10
                            x_known(2) = x_known(2) * 1.001
                            x_known(1) = 1.d0 - x_known(2)
                            call Flash_pT(gl,press, T, x_known, x_vap, x_liq, gl%rho_vap, gl%rho_liq, beta_lc, errval, iter)
                            nrofphases = 2
                            if (errval==0) exit

                        end do
                    end if


                endif
            end if
        elseif (Flash_Type == 2)then    !Phaseboundary based Flash  !TN

            if ((x_known(1) < 1.d-10) .and. (points>1)) then
                if ((psat_2 > 1.d-10) .and. (press/psat_2 >0.95)) then
                    errval = 0
                    ! calculate and save the pure fluid as the first point
                    d_vap = rhovap_points(points-1)
                    d_liq = rholiq_points(points-1)
                    call Flash_Pure_PhaseBoundary(gl,psat_2, t, d_vap, d_liq, 1, errPsat, iter, 2)
                    if (errPsat==0)then
                        T_points(points) = T
                        p_points(points) = psat_2
                        !pstart = psat_1
                        x_points(points, 1) = 0.d0 !xvap1
                        x_points(points, 2) = 1.d0 !xvap2
                        x_points(points, 3) = 0.d0 !xliq1
                        x_points(points, 4) = 1.d0 !xliq2
                        rhovap_points(points) = d_vap
                        rholiq_points(points) = d_liq
                    else
                        points = points - 1
                    endif
                    exit
                else
                    points = points - 1
                    save_calc = 2
                endif
                !exit
            endif

            !If (T/gl%tc(1) >= 0.95d0) then
            !    x_step = 2.d-3
            !else
            !    x_step = 1.d-2
            !endif

            dpdx(1) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,1)-x_points(points-2,1)) !dewline
            dpdx(2) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,3)-x_points(points-2,3)) !bubbleline


            if ((abs((x_vap(1)-x_liq(1))) < 5.d-3))then
                If (T/gl%tc(1) >= 0.95d0) then
                    x_step = (-1.00606d-1 * (T/gl%tc(1)) + 1.00576d-1) / 2.d0
                else
                    x_step = 2.d-3
                endif
            end if

            x_known = 0.d0

            !check which line is less steep and calculate the next composition step there   !TN
            if (abs(dpdx(1)) < abs(dpdx(2))) then !on dewline
                iFlash = 4
                !iPhase_try = 5
                x_known(2) = x_points(points-1,2) - x_step
                x_known(1) = 1.d0 - x_known(2)
                x_vap = x_known
                x_liq(1) = x_points(points-1, 3)
                x_liq(2) = x_points(points-1, 4)
            else    !on bubbleline
                iFlash = 3
                !iPhase_try = 1
                x_known(2) = x_points(points-1,4) + x_step
                x_known(1) = 1.d0 - x_known(2)
                x_liq = x_known
                x_vap(1) = x_points(points-1, 1)
                x_vap(2) = x_points(points-1, 2)
            end if

            press = p_points(points-1)

            iPhase_try = 0  !not sure, but works so far !TN


            call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, gl%rho_vap, gl%rho_liq, beta_lc, iFlash,&
                & iPhase_try, nrofphases, errval, iter)

        end if



        if ((errval == 0) .and. (save_calc /= 2)) then

            if ((points > 2) .and. ((((((press-p_points(points-1)) < 1.d-10)) .or. ((x_liq(1) - x_vap(1)) < 1.d-10))) .and. ((T > gl%tc(1)) .or. T > gl%tc(2))) .or. (abs((press-p_points(points-1))/press) < 1.d-6)) then !stop criteria if ((press is going down again) .or. (the phases are switch) .and. one component is supercritical (for azeotrops)) .or. if (pressure change is too small) !TN

                if (open_count < 3) then
                    open_count = open_count + 1
                    cycle
                else
                    save_calc = 2
                    nrofphases = 0
                    force_spline = 0

                    call Phasedet(gl,press, T, x_known, rho, x_Phase, phasetype, beta_lc, nrofphases, errval)
                    if (nrofphases == 1) then
                        !pressure is above VLE Region -> if sign of slopes of dew and bubble lines are different, closed phase envelope -> spline calculation
                        dpdx(1) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,1)-x_points(points-2,1)) !dewline
                        dpdx(2) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,3)-x_points(points-2,3)) !bubble line
                        if (dpdx(1)/dpdx(2) < 0 ) then
                            !errval = -1
                            save_calc = 2 !quit calculations and go to spline calculations  !TN
                            force_spline = 1
                        endif
                    endif
                endif
            else

                save_calc = 1
                !open_count = open_count + 1

            endif

        elseif (errval /= 0) then

            if (Flash_type == 1)then
                !if the pt based flash fails try phaseboundary based flash !TN
                save_calc = 2

            elseif (Flash_type == 2)then

                save_calc = 2
                nrofphases = 0
                force_spline = 0

                call Phasedet(gl,press, T, x_known, rho, x_Phase, phasetype, beta_lc, nrofphases, errval)
                if (nrofphases == 1) then
                    !pressure is above VLE Region -> if sign of slopes of dew and bubble lines are different, closed phase envelope -> spline calculation
                    dpdx(1) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,1)-x_points(points-2,1)) !dewline
                    dpdx(2) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,3)-x_points(points-2,3)) !bubble line
                    if (dpdx(1)/dpdx(2) < 0 ) then
                        !errval = -1
                        save_calc = 2
                        force_spline = 1
                    endif
                endif
            end if

        endif


        !if calculations failed twice (after switching from dew to bubbleline or vice versa) cancel calculations
        if ((save_calc == 0) .and. (switchsupercrit == 0)) then
            return
            !calculation successful
        elseif (save_calc == 1) then!(errval == 0) then!((x_points(points-1,1) - x_vap(1) > 0.d0) .or. (x_points(points-1,3) - x_liq(1) > 0.d0))) then
            if (open_count <= 1) then
                !if step was to large on the dewline, change flash from T-x' to T-x''
                !-------
                !Nothing changed, but for some mixtures this section needs to be commented. Andreas Jäger, September 2017
                if ((x_points(points-1,1) - x_vap(1)) >= 0.07d0) then   !changed from 0.05 to 0.07 10/19 TN
                    iflash = 4
                    x_known(1) = x_points(points-1,1) !xvap1
                    x_known(2) = x_points(points-1,2) !xvap2
                    !if open vle region, the step will be also very large -> should happen when composition on bubble line or dewline is given either way
                    open_count = open_count + 1
                    cycle
                    !if step was to large on the bubbleline, change flash from T-x'' to T-x'
                elseif ((x_points(points-1,3) - x_liq(1)) > 0.07d0) then    !changed from 0.05 to 0.07 10/19 TN
                    iflash = 3
                    x_known(1) = x_points(points-1,3) !xliq1
                    x_known(2) = x_points(points-1,4) !xliq2
                    !if open vle region, the step will be also very large -> should happen when composition on bubble line or dewline is given either way
                    open_count = open_count + 1
                    cycle
                endif
                if ((open_count == 1) .and. ((dabs(x_points(points-1,3) - x_liq(1)) < 1.d-4))) then
                    open_count = open_count + 1
                    cycle
                endif
                !-------
                open_count = 0
            endif

            ! write the results to module variables
            T_points(points) = T
            p_points(points) = press
            x_points(points, 1) = x_vap(1)  !xvap1
            x_points(points, 2) = x_vap(2)  !xvap2
            x_points(points, 3) = x_liq(1)  !xliq1199
            x_points(points, 4) = x_liq(2)  !xliq2
            if (gl%rho_vap < gl%rho_liq) then
                rhovap_points(points) = gl%rho_vap
                rholiq_points(points) = gl%rho_liq
            else
                rhovap_points(points) = gl%rho_liq
                rholiq_points(points) = gl%rho_vap
            endif

            points = points + 1

            !save slope of curves
            dpdx_old = dpdx
            dpdx(1) = -(p_points(points-1) - p_points(points-2))/dabs(x_points(points-1,1) - x_points(points-2,1)) !slope of dewline
            dpdx(2) = -(p_points(points-1) - p_points(points-2))/dabs(x_points(points-1,3) - x_points(points-2,3)) !slope of bubble line

            !error occurred
        elseif (save_calc == 2) then

            if((Flash_type ==1) .and. (switchfluids /= 2))then

                if (points > 2)then
                    Flash_type = 2
                    save_calc = 0
                else
                    exit
                endif
            else
                !TN: Does this make sense here?
                !if (psat_2 > 1.d-10) then
                !    ! calculate and save the pure fluid as the first point
                !    call Flash_Pure_PhaseBoundary(gl,psat_2, t, d_vap, d_liq, 1, errPsat, iter, 2)
                !    T_points(points) = T
                !    p_points(points) = psat_2
                !    !pstart = psat_1
                !    x_points(points, 1) = 0.d0 !xvap1
                !    x_points(points, 2) = 1.d0 !xvap2
                !    x_points(points, 3) = 0.d0 !xliq1
                !    x_points(points, 4) = 1.d0 !xliq2
                !    rhovap_points(points) = d_vap
                !    rholiq_points(points) = d_liq
                !
                !    exit
                !else


                if (points > (n/2 + 1)) then    !Elsewise, the code crashes, Andreas Jäger, September 2017
                    points = points - 1

                    !for less points spline does not makes sense
                    if (points<n/2)then
                        exit
                    end if

                    !--------------------------------------------------------------------------------
                    ! FINAL ESTIMATE OF THE CRITICAL POINT (SPLINE INTERPOLATION)
                    !--------------------------------------------------------------------------------
                    ! use the last n/2 points for the generation of the spline parameters
                    ! the pressure is used as a function of rholiq and rhovap in order to interpolate
                    ! the critical density, so actually n points are used

                    !--------------------------------------------------------------------------------
                    ! first, check if the steps above are close to the critical point. If not,
                    ! the two-phase region may not have a critical point and an estimate would be
                    ! meaningless.
                    if ((dabs(rhovap_points(points) - rholiq_points(points))/rholiq_points(points) > 3.d-1) .or. (dabs(x_points(points,1) - x_points(points,3)) > 0.05d0 ) .or. (force_spline == 0)) then
                        errMSG = '*** No critical point found at this temperature ***'
                        !goto 111
                        !goto 110
                        exit
                    end if
                    k = n/2
                    j = points-k
                    xi(1:k) = rhovap_points(j+1:points)
                    yp(1:k) = p_points(j+1:points)
                    yx(1:k) = x_points(j+1:points,1)
                    xi(k+1:2*k) = rholiq_points(points:j+1:-1)
                    if (maxloc(xi,1) /= n) then !density must be in strictly increasing order
                        errval = 0
                        exit
                    endif
                    yp(k+1:2*k) = p_points(points:j+1:-1)
                    yx(k+1:2*k) = x_points(points:j+1:-1,3)
                    ! generate the spline coefficients
                    call spline (gl,xi, yp, bp, cp, dp, n)
                    ! generate the spline coefficients
                    call spline (gl,xi, yx, bx, cx, dx, n)

                    !Catch NaN
                    if (any(isnan(xi)) .or. any(isnan(yp)) .or. any(isnan(bp)) .or. any(isnan(cp)) .or. any(isnan(dp)) .or. any(isnan(yx)) .or. any(isnan(bx)) .or. any(isnan(cx)) .or. any(isnan(dx)))then
                        !errval = -5532 !noch hinzufügen
                        exit
                    end if

                    d_vap = xi(k)
                    d_liq = xi(k+1)
                    delta = (d_liq - d_vap)/20.d0
                    ! do 10 more density controlled steps toward the critical point
                    do j = 1, 10
                        rho1 = d_vap + real(j)*delta
                        rho2 = d_liq - real(j)*delta
                        press = ispline(gl,rho1, xi, yp, bp, cp, dp, n)
                        x_vap(1) = ispline(gl,rho1, xi, yx, bx, cx, dx, n)
                        x_vap(2) = 1.d0 - x_vap(1)
                        x_liq(1) = ispline(gl,rho2, xi, yx, bx, cx, dx, n)
                        x_liq(2) = 1.d0 - x_liq(1)
                        !catch case when dew or bubble line is very steep -> spline interpolation is not working
                        !if ((x_vap(1) < 0.d0) .or. (x_vap(1) > 1.d0) .or.(x_liq(1) < 0.d0) .or. (x_liq(1) > 1.d0)) then
                        !    errval = 0
                        !    exit !cancel interpolation
                        !endif



                        ! write the results to module variables
                        points = points + 1
                        T_points(points) = T
                        p_points(points) = press
                        x_points(points, 1) = x_vap(1)
                        x_points(points, 2) = x_vap(2)
                        x_points(points, 3) = x_liq(1)
                        x_points(points, 4) = x_liq(2)
                        rhovap_points(points) = d_vap + real(j)*delta
                        rholiq_points(points) = d_liq - real(j)*delta

                    end do
                    errval = 0
                    exit
                end if

                !end if
                exit
            end if

        end if

    end do

    if ((points == maxi-points_copy-n) .and. (save_calc == 1)) points = points - 1
    !press = pstart*1.01d0
    !molfractions(2) = 1.d0 - molfractions(1)
    !x_known = molfractions
    !call reduced_parameters_calc(gl,300.d0)

#IF DEFINED(_DEBUG)     
    !--------------------------------------------------------------------------------
    ! CHECK IF WRITING THE DATA TO A FILE IS WISHED
    !--------------------------------------------------------------------------------
    !111 continue
    if (fileout /= '') then
        open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
        if (error /= 0) return
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'! P-X DIAGRAM DATA FILE FOR THE BINARY SYSTEM'
        write(13,*)'! ', trim(gl%components(1)),'-',trim(gl%components(2))
        write(13,*)'!'
        write(13,'(A20, F8.3)')'! TEMPERATURE [K]: ', T
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'!'
        write(13,*)'! p/MPa    ','T/K       ','x',trim(gl%components(1)),'    x',trim(gl%components(2)),'    y',trim(gl%components(1)), &
            &'    y',trim(gl%components(2)),'   dl/mol/l ', 'dv/mol/l'
        do i = 1, points
            write(13,1000) p_points(i), T_points(i), x_points(i, 3), x_points(i, 4), x_points(i, 1), x_points(i, 2), &
                & rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
        end do
        if (errMSG /= '') write(13,*)'! ', errMSG
        close(13)
    end if
#ENDIF

    end subroutine pxdiag


    subroutine txdiag(gl,press, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, switchfluids, errval)
    !--------------------------------------------------------------------------------
    ! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    ! PRESSURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    !--------------------------------------------------------------------------------
    ! J. Gernert, March 2012
    ! S.Hielscher und T.Neumann, July 2017





    implicit none

    type(type_gl) :: gl

    !--------------------------------------------------------------------------------
    integer, parameter:: maxi = 300 ! maximum number of calculated points, lengt of the return vectors!(300)
    double precision, intent(in):: press     ! specified pressure for the T-x diagram
    double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision, intent(out):: x_points(maxi,4)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_vap(2)
    ! x_points(i,3) = x_liq(1), x_points(i,4) = x_liq(2)
    character(255), intent(in):: fileout    ! optional path for output file
    !switchfluids = 1: fluid order was switched -> this has to be undone by swapping point order and molefractions x1 -> x2, x2 -> x1
    !switchfluids = 2: former azeotropic mixture, splited in two vle regions (co2-ethane T > 290 K)
    !swtichfluids = 3: Mostly relevant for open VLE regions -> start with component with higher saturation temperature, because this component will most likely be attached to y axis while the other component is detached from y axis
    !switchsupercrit: indicates if the point of minimal x(1) on the dewline for a binary mixtures where one component is supercritical was passed
    integer:: points, switchfluids, switchsupercrit, errval                ! number of calculated points, error code
    !--------------------------------------------------------------------------------

    ! for calling pure fluid routines
    double precision:: tsat_1, tsat_2, Tsat_high
    integer::mixflag, errPsat,stepsize_counter
    ! for calling PhaseDet
    double precision:: rho(5), x_Phase(30,5), x_known(30)
    integer::nrofphases, phasetype(5)
    ! for calling Flash_pt and Flash_phase_boundary:
    double precision::x_vap(30), x_liq(30), beta_lc
    integer::iter, save_calc, force_spline !0: error, don't save 1: save, 2: finish with spline, 3: finish without spline
    ! for calling the spline routines:
    integer, parameter:: n = 10 ! number of points used for the spline interpolation of the critical point
    double precision, dimension (n) :: xi(n), yT(n), bT(n), cT(n), dT(n), yx(n), bx(n), cx(n), dx(n)   ! parameters of the spline interpolation
    ! internal variables:
    double precision::x_old, x_crit(2),Tstart, deltap, deltaT_max, dTdx(2),dTdx_point(2), dTdx_old(2), d_vap, d_liq, rho1, rho2, rhocrit, delta, Temp, x_backup, press_backup, deltap_backup, delta_backup, x_step
    double precision, dimension(maxi):: T_pts_right, p_pts_right, rhovap_pts_right, rholiq_pts_right
    double precision:: x_pts_right(maxi,30)
    logical::TwoPhaseFound, azeo
    !                                                                   counter for open two phase region
    integer:: iFlash,iFlash2, i, j, k, m, small, large, error, dens_count, lastpoint, open_count
    character(255):: errMSG

    ! output file format
1000 format(f10.5, f10.4, f10.6, f10.6, f10.6, f10.6, f10.3, f10.3)


    p_points = 0.d0; T_points = 0.d0; x_points = 0.d0; rhovap_points = 0.d0; rholiq_points = 0.d0;
    p_pts_right = 0.d0; T_pts_right = 0.d0; x_pts_right = 0.d0; rhovap_pts_right = 0.d0; rholiq_pts_right = 0.d0;
    x_vap = 0.d0; x_liq = 0.d0; points = 1; errval = 0; errMSG = ''
    d_vap = 0.d0; d_liq = 0.d0
    azeo = .false.
    deltaT_max = 0.5d0  !maximum pressure step
    tsat_high = 0.d0
    switchsupercrit = 0
    force_spline = 0
    save_calc = 0
    dtdx = 0.d0
    dtdx_old = 0.d0
    open_count = 0
    dens_count = 0
    stepsize_counter = 0

    !--------------------------------------------------------------------------------
    ! STEP 1: START VALUES
    !--------------------------------------------------------------------------------
    ! get the saturation temperature of the pure fluids
    tsat_1 = 0.d0
    tsat_2 = 0.d0

    tsat_1 = Estimate_Tsat(gl,press, 1)
    if (press < gl%pc(2)) then
        tsat_2 = Estimate_Tsat(gl,press, 2)
    end if
    !relevant for open VLE regions, check head
    if ((tsat_1 < tsat_2) .and. (tsat_2 > 1.d-10) .and. (switchfluids == 0)) then ! .and. (tsat_2/tsat_1) > 1.d2) then
        switchfluids = 3
        return
    endif
    !Check if both components are below their critical pressures
    if ((press < gl%pc(1)) .and. (press < gl%pc(2))) then
        ! calculate and save the pure fluid as the first point
        call Flash_Pure_PhaseBoundary(gl,press, tsat_1, d_vap, d_liq, 2, errPsat, iter, 1)
        T_points(points) = tsat_1
        Tstart = Tsat_1
        p_points(points) = press
        x_points(points, 1) = 1.d0 !xvap1
        x_points(points, 2) = 0.d0 !xvap2
        x_points(points, 3) = 1.d0 !xliq1
        x_points(points, 4) = 0.d0 !xliq2

        x_known(1) = 1.d0
        x_known(2) = 1.d0 - x_known(1)
    else
        !Iterate startpoint:
        !if pressure is higher than pc of both components, start with x(1)=0.9 and check for two phase region by performing p x'' Flash
        !if homogeneous lower amount of (1) by 0.1 and so on.
        outer: do i = 90,10,-10
            x_vap = 0.d0
            x_liq = 0.d0
            x_known(1) = i * 1.d-2
            x_known(2) = 1.d0 - x_known(1)
            iflash = 2
            call Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, 0.d0, 0.d0, beta_lc, iFlash, 0, errval, iter)
            if (errval == 0) then
                !after two phase region was found, raise amount of (1) by 0.01 until homogeneous again. The last heterogeneous point is taken for startvalues
                do j = 1,10
                    x_known(1) = x_known(1) + 1.d-2
                    x_known(2) = 1.d0 - x_known(1)
                    x_vap = 0.d0
                    x_liq = 0.d0
                    call Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, 0.d0, 0.d0, beta_lc, iFlash, 0, errval, iter)
                    if (errval /= 0) then
                        x_known(1) = x_known(1) - 1.d-2
                        x_known(2) = 1.d0 - x_known(1)
                        x_vap = 0.d0
                        x_liq = 0.d0
                        errval = 0
                        call Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, 0.d0, 0.d0, beta_lc, iFlash, 0, errval, iter)
                        T_points(points) = temp
                        p_points(points) = press
                        x_points(points, 1) = x_vap(1) !xvap1
                        x_points(points, 2) = x_vap(2) !xvap2
                        x_points(points, 3) = x_liq(1) !xliq1
                        x_points(points, 4) = x_liq(2) !xliq2
                        x_known = x_liq
                        exit outer
                    endif
                end do
            endif
        end do outer
    endif


    rhovap_points(points) = d_vap
    rholiq_points(points) = d_liq
    points = points + 1
    !--------------------------------------------------------------------------------
    ! STEP 2: FIRST MIXTURE VLE POINT
    !--------------------------------------------------------------------------------
    If (press/gl%pc(1) >= 0.95d0) then
        x_step = 5.d-3
    else
        x_step = 1.d-2
    endif



    iFlash = 1 ! start with p x' Flash

    Temp = T_points(points)
    do while (points < maxi-n)
        !if (open_count <= 1) then
        save_calc = 0
        !when approaching the critical point abs(drho/dx) gets very large -> reduce x_step (MAYBE DELETED IN THE FUTURE)
        !if (points > 2) then
        !    if ((dens_count == 0) .and. (dabs((rhovap_points(points-1) - rhovap_points(points-2))/x_step) > 1.d5)) then
        !        x_step = x_step / 10.d0
        !        dens_count = dens_count + 1
        !    endif
        !endif
        !if (points>100)then
        !x_step = 2.d-5
        !endif
        !
        !        if (points>500)then
        !x_step = 1.d-5
        !endif



        do while((stepsize_counter < 8) .and. ((x_known(1) - x_step)<0.d0))     !make sure that x_step is small enough so that x_known stays within its limits TN 5.2.19
            x_step=x_step/2.d0
            stepsize_counter =stepsize_counter+1
        end do


        x_known(1) = x_known(1) - x_step


        !quit if the whole molesrange was passed
        if (x_known(1) < 1.d-10) then
            if ((Tsat_2 > 1.d-10) .and. (Tsat_2 < T_points(points-1)) .and. (T_points(points-1)-Tsat_2 <10)) then !added second test for not closed t,x diagrams TN 5.2.19
                ! calculate and save the pure fluid as the first point

                d_vap = rhovap_points(points-1)
                d_liq = rholiq_points(points-1)
                call Flash_Pure_PhaseBoundary(gl,press, Tsat_2, d_vap, d_liq, 2, errPsat, iter, 2)
                if (errPsat==0)then
                    T_points(points) = Tsat_2
                    p_points(points) = press
                    !pstart = psat_1
                    x_points(points, 1) = 0.d0 !xvap1
                    x_points(points, 2) = 1.d0 !xvap2
                    x_points(points, 3) = 0.d0 !xliq1
                    x_points(points, 4) = 1.d0 !xliq2
                    rhovap_points(points) = d_vap
                    rholiq_points(points) = d_liq
                else
                    points = points - 1
                endif
            else
                points = points - 1
            endif
            exit
        endif
        x_known(2) = 1.d0 - x_known(1)
        x_liq = 0.d0
        x_vap = 0.d0
        call Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, rhovap_points(points-1), rholiq_points(points-1), beta_lc, iFlash,&
            &  0, errval, iter)
        if ((points > 2).and.(errval == 0)) then
            if ((iflash == 1).and.((x_vap(1) - x_points(points-1,1))/(x_points(points-1,1)-x_points(points-2,1)) < 0.d0) ) then
                errval = -10 !abitray errorvalue to calculate flash again without startvalues
            elseif ((iflash == 2).and.((x_liq(1) - x_points(points-1,3))/(x_points(points-1,3)-x_points(points-2,3)) < 0.d0) ) then
                errval = -10 !abitray errorvalue to calculate flash again without startvalues
            endif
        endif
        if (errval /= 0) then
            call Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, 0.d0, 0.d0, beta_lc, iFlash,&
                &  0, errval, iter)
        endif
        !Catch NaN
        if (Temp /= Temp) then
            switchfluids = 0
            if (errval == 0) errval = -898964
            return
        endif



        if (iflash==1)then
            iflash2=2
        else
            iflash2=1
        end if
        dTdx_point(2) = -(Temp - T_points(points-1))/dabs(x_vap(1) - x_points(points-1,1)) !slope of dewline
        dTdx_point(1) = -(Temp - T_points(points-1))/dabs(x_liq(1) - x_points(points-1,3)) !slope of bubble line

        if ((((abs(Temp-T_points(points-1))>10.d0).and.(abs(dTdx_point(iflash))>1000.d0)) .or. (errval /= 0)) .and. (stepsize_counter < 8) .and. (open_count > 0))then !reduce stepsize by half up to to 8 times if temperature step and dtdx to big TN 5.2.19

            if (dabs(dTdx_point(iflash))<=dabs(dTdx_point(iflash2)))then

                if( iflash == 1)then
                    x_known(1) = x_points(points-1,3) !xliq1
                    x_known(2) = x_points(points-1,4) !xliq2
                elseif (iFlash==2)then
                    x_known(1) = x_points(points-1,1) !xvap1
                    x_known(2) = x_points(points-1,2) !xvap2
                end if

                x_step=x_step/2.d0
                stepsize_counter =stepsize_counter+1


                cycle
            else
                if (stepsize_counter<4)then     !even if the flash is calculated on the saturation line with the higher dtdx first try to reduce stepsize 4 times and then swith flash

                    if( iflash == 1)then
                        x_known(1) = x_points(points-1,3) !xliq1
                        x_known(2) = x_points(points-1,4) !xliq2
                    elseif (iFlash==2)then
                        x_known(1) = x_points(points-1,1) !xvap1
                        x_known(2) = x_points(points-1,2) !xvap2
                    end if

                    x_step=x_step/2.d0
                    stepsize_counter =stepsize_counter+1


                    cycle

                else
                    if (open_count<1)then

                        if (iflash==1)then
                            !change flash from p-x' to p-x''
                            iflash = 2

                            x_known(1) = x_points(points-1,1) !xvap1
                            x_known(2) = x_points(points-1,2) !xvap2

                            open_count = open_count + 1
                            cycle

                        elseif (iflash==2) then
                            !change flash from p-x'' to p-x'
                            iflash = 1

                            x_known(1) = x_points(points-1,3) !xliq1
                            x_known(2) = x_points(points-1,4) !xliq2

                            open_count = open_count + 1
                            cycle

                        end if
                        open_count = 0
                    end if



                end if

            end if

        end if


        !everything fine, save variables and continue calculations
        if (errval == 0) save_calc = 1

        !reset a false set of switchsupercrit = 1
        !Example: When the composition of bubble line is the input composition, the composition on the dewline is calculated. If one component is overcritical it could happen that the amount of
        !component two increases on the dewline at higher temperatures and decreases at lower temperatures. If the flag switch supercrit is one, but the amount of component two still increases the flag is set back to 0
        if ((save_calc == 1) .and. (switchsupercrit == 1)) then
            if (x_vap(1) - x_points(points - 1, 1)  < 0.d0) then
                switchsupercrit = 0
            endif
        elseif ((save_calc == 1) .and. (iflash == 1) .and. (x_vap(1) - x_points(points-1,1) > 0.d0)) then
            switchsupercrit = 1
        endif

        !the dewline (or bubble line in rare cases) might change the direction of molar composition from positive to negative direction
        !when this happens switchsupercrit is set to 1 (from 0)
        !if the composition on the dewline was the input composition, the composition on the bubbleline is now the input composition
        !when the pressure now is decreasing, calculations are finished with spline interpolations
        if ((switchsupercrit == 1) .and. (press - p_points(points-1) < 0.d0) .and. (errval == 0)) then
            save_calc = 3 !quit calculations and go to spline calculations
            !errval = -10 !arbitrary error value
        endif

        !if the given temperature is >0.95*Tc it (hopefully) does not happen, that the slope of bubble and dewline do not change
        !if this happens, the previous step was wrong and saved values are discarded and calculations are finished with spline interpolations
        !if (((T/tc(1) >= 0.95d0).and.(tc(2) < 1.d-12)) .and. (dpdx(1)*dpdx_old(1) < 0.d0).and.(dpdx(2)*dpdx_old(2) < 0.d0)) then
        !    save_calc = 2
        !    points = points - 1
        !endif

        !if right component is supercritical it might happen that the current composition step is out of the two phase region
        !-> try changing form dew to bubble line or vice versa
        if (((errval /= 0) .and. (Tsat_2 < 1.d-10) .and. (switchsupercrit == 0)))then
            if (iflash == 1) then ! p x' Flash -> p x'' Flash
                if ((rholiq_points(points-1)-rhovap_points(points-1))/(rhovap_points(points-1)) < 1.d-1) then
                    save_calc = 2
                else
                    iflash = 2
                    x_known(1) = x_points(points-1,1) !xvap1
                    x_known(2) = x_points(points-1,2) !xvap2
                endif
            elseif (iflash == 2) then! p x'' Flash -> p x' flash
                iflash = 1
                x_known(1) = x_points(points-1,3) !xliq1
                x_known(2) = x_points(points-1,4) !xliq2
            endif
            switchsupercrit = switchsupercrit + 1
            cycle
        elseif ((errval /= 0) .and. (switchsupercrit == 1) .and. (dabs(rhovap_points(points-1) - rholiq_points(points-1))/rholiq_points(points-1) <= 1.d-1)) then
            save_calc = 2
        elseif ((errval /= 0) .and. (switchsupercrit == 1) .and. (dabs(rhovap_points(points-1) - rholiq_points(points-1))/rholiq_points(points-1) > 1.d-1)) then
            if (dabs(x_vap(1) - x_liq(1)) > 5.d-2) then
                Temp = T_points(points-1)
                x_vap(1) = x_points(points-1, 1)   !xvap1
                x_vap(2) = x_points(points-1, 2)   !xvap2
                x_liq(1) = x_points(points-1, 3)   !xliq1
                x_liq(2) = x_points(points-1, 4)   !xliq2
                if (points == 2) then
                    open_count = 2
                else
                    open_count = 3
                endif
                cycle
            else
                points = points - 1
                errval = 0
                exit
            endif
            !elseif ((dabs(gl%rho_vap - gl%rho_liq)/gl%rho_vap < 1.d-1) .and. (press/minval(gl%pc) > 0.85) .and. (errval /= 0)) then        !deleted because it caused problems, maybe one has to rethink here TN 5.2.19
            !    switchfluids = 2
            !    errval = 0
            !    save_calc = 2
        elseif ((errval /= 0)) then        ! TN 5.2.19

            errval = 0
            save_calc = 2
        endif



        if ((abs(Temp-T_points(points-1))<5.d0) .and. (stepsize_counter>0))then     !if temperature step is to small start increasing stepsize again

            x_step=x_step*2.d0
            stepsize_counter=stepsize_counter-1

        end if


        !elseif (open_count >= 2) then
        !    save_calc = 0
        !    !first pressure based vle point
        !    if (open_count == 2) then
        !        TwoPhaseFound = .false.
        !        x_known(1) = 0.999d0
        !        x_known(2) = 1.d0 - x_known(1)
        !        press = 1.01d0 * psat_1
        !        call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
        !
        !        count = 0
        !        if (nrofphases < 2) then
        !            !press = pstart*1.01d0
        !            do while(nrofphases < 2)
        !                press = (press + pstart)/2.d0
        !                count = count + 1
        !                if (count >  10) then
        !                    errval = -2222
        !                    exit
        !                end if
        !                call PhaseDet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
        !            end do
        !        end if
        !        x_liq = x_Phase(:,phasetype(2))
        !        x_vap = x_Phase(:,phasetype(1))
        !        rho_liq = rho(phasetype(2))
        !        rho_vap = rho(phasetype(1))
        !    !follwing pressure based points
        !    elseif (open_count > 2) then
        !        x_old = (x_vap(2) + x_liq(2))/2.d0
        !        dpdx = (press - pstart)/x_old    !calculate the slope of the pressure with the overall x
        !        deltap = press*0.1d0     !abs(psat_1 - psat_2)/50.d0
        !        if (deltap > deltaT_max) deltap = deltaT_max    !maximum pressure step
        !        press = press + deltap
        !        x_known(2) = (press - pstart)/dpdx(1) !look for new composition in the 2-phase-region
        !        !Very likely not happen, delete?
        !        do while (x_known(2) >= 0.999d0) !in vicinity of the pure component
        !            press = press - deltap
        !            deltap = deltap/2.d0
        !            press = press + deltap
        !            x_known(2) = (press - pstart)/dpdx(1)
        !            if (deltap < 1.d-5) then
        !                errval = -3333
        !                p_points = errval
        !                exit
        !            end if
        !        end do
        !        !end delete?
        !        x_known(1) = 1.d0 - x_known(2)
        !        ! calculate the pT-flash
        !        call Flash_pT(press, T, x_known, x_vap, x_liq, rho_vap, rho_liq, beta, errval, iter)
        !        nrofphases = 2
        !        if (errval /= 0) then
        !            nrofphases = 0
        !            call Phasedet(press, T, x_known, rho, x_Phase, phasetype, beta, nrofphases, errval)
        !            if (nrofphases == 1) then
        !                !pressure is above VLE Region -> if sign of slopes of dew and bubble lines are different, closed phase envelope -> spline calculation
        !                dpdx(1) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,1)-x_points(points-2,1)) !dewline
        !                dpdx(2) = (p_points(points-1)-p_points(points-2))/(x_points(points-1,3)-x_points(points-2,3)) !bubble line
        !                if (dpdx(1)/dpdx(2) < 0 ) then
        !                    errval = -1
        !                    save_calc = 2
        !                    force_spline = 1
        !                endif
        !            endif
        !        endif
        !
        !
        !    endif
        !    if (errval == 0) then
        !        save_calc = 1
        !        open_count = open_count + 1
        !    endif
        !endif

        !if calculations failed twice (after switching from dew to bubbleline or vice versa) cancel calculations
        if ((save_calc == 0) .and. (switchsupercrit == 0)) then
            continue
            return

            !calculation successful
        elseif (save_calc == 1) then!(errval == 0) then!((x_points(points-1,1) - x_vap(1) > 0.d0) .or. (x_points(points-1,3) - x_liq(1) > 0.d0))) then
            if (open_count < 1) then

                if(abs(dTdx_point(iflash)) >= abs(dTdx_point(iflash2)))then !TN 5.2.19: this is new, changed the criteria for changing the flash calculation type, always calculate on the saturation line with the smaller dtdx (old code and criteria below)
                    if (iflash==1)then
                        !change flash from p-x' to p-x''
                        iflash = 2

                        x_known(1) = x_points(points-1,1) !xvap1
                        x_known(2) = x_points(points-1,2) !xvap2

                        open_count = open_count + 1
                        cycle

                    elseif (iflash==2) then
                        !change flash from p-x'' to p-x'
                        iflash = 1

                        x_known(1) = x_points(points-1,3) !xliq1
                        x_known(2) = x_points(points-1,4) !xliq2

                        open_count = open_count + 1
                        cycle

                    endif

                    open_count = 0

                end if

                open_count = 0

            else if(open_count == 1)then
                if (((iFlash == 2) .and. ((dabs(x_points(points-1,3) - x_liq(1)) < 1.d-4))) .or. ((iFlash == 1) .and. ((dabs(x_points(points-1,1) - x_vap(1)) < 1.d-4))))  then
                    open_count = open_count + 1
                    cycle
                endif
                open_count = 0

                !!if step was to large on the dewline, change flash from p-x' to p-x''
                !if ((x_points(points-1,1) - x_vap(1)) >= 0.05d0) then
                !    iflash = 2
                !    x_known(1) = x_points(points-1,1) !xvap1
                !    x_known(2) = x_points(points-1,2) !xvap2
                !    !if open vle region, the step will be also very large -> should happen when composition on bubble line or dewline is given either way
                !    open_count = open_count + 1
                !    cycle
                !    !if step was to large on the bubbleline, change flash from p-x'' to p-x'
                !elseif ((x_points(points-1,3) - x_liq(1)) > 0.05d0) then
                !    iflash = 1
                !    x_known(1) = x_points(points-1,3) !xliq1
                !    x_known(2) = x_points(points-1,4) !xliq2
                !    !if open vle region, the step will be also very large -> should happen when composition on bubble line or dewline is given either way
                !    open_count = open_count + 1
                !    cycle
                !endif
                !if ((open_count == 1) .and. ((dabs(x_points(points-1,3) - x_liq(1)) < 1.d-4))) then
                !    open_count = open_count + 1
                !    cycle
                !endif
                !open_count = 0

            endif


            ! write the results to module variables
            T_points(points) = Temp
            p_points(points) = press
            x_points(points, 1) = x_vap(1)  !xvap1
            x_points(points, 2) = x_vap(2)  !xvap2
            x_points(points, 3) = x_liq(1)  !xliq1
            x_points(points, 4) = x_liq(2)  !xliq2
            if (gl%rho_vap < gl%rho_liq) then
                rhovap_points(points) = gl%rho_vap
                rholiq_points(points) = gl%rho_liq
            else
                rhovap_points(points) = gl%rho_liq
                rholiq_points(points) = gl%rho_vap
            endif
            points = points + 1

            !save slope of curves
            dTdx_old = dTdx
            dTdx(1) = -(T_points(points-1) - T_points(points-2))/dabs(x_points(points-1,1) - x_points(points-2,1)) !slope of dewline
            dTdx(2) = -(T_points(points-1) - T_points(points-2))/dabs(x_points(points-1,3) - x_points(points-2,3)) !slope of bubble line

            !error occurred
        elseif (save_calc == 2) then
            points = points - 1

            !for less points spline does not makes sense
            if (points<n/2)then
                exit
            end if

            !--------------------------------------------------------------------------------
            ! FINAL ESTIMATE OF THE CRITICAL POINT (SPLINE INTERPOLATION)
            !--------------------------------------------------------------------------------
            ! use the last n/2 points for the generation of the spline parameters
            ! the temperature is used as a function of rholiq and rhovap in order to interpolate
            ! the critical density, so actually n points are used

            !--------------------------------------------------------------------------------
            ! first, check if the steps above are close to the critical point. If not,
            ! the two-phase region may not have a critical point and an estimate would be
            ! meaningless.
            if ((dabs(rhovap_points(points) - rholiq_points(points))/rholiq_points(points) > 3.d-1) .or. (dabs(x_points(points,1) - x_points(points,3)) > 0.05d0 ) .or. (force_spline == 0))  then
                errMSG = '*** No critical point found at this temperature ***'
                !goto 111
                !goto 110
                exit
            end if
            k = n/2
            j = points-k
            xi(1:k) = rhovap_points(j+1:points)
            yT(1:k) = T_points(j+1:points)
            yx(1:k) = x_points(j+1:points,1)
            xi(k+1:2*k) = rholiq_points(points:j+1:-1)
            if (maxloc(xi,1) /= n) then !density must be in strictly increasing order
                errval = 0
                exit
            endif
            yT(k+1:2*k) = T_points(points:j+1:-1)
            yx(k+1:2*k) = x_points(points:j+1:-1,3)
            ! generate the spline coefficients
            call spline (gl,xi, yT, bT, cT, dT, n)
            ! generate the spline coefficients
            call spline (gl,xi, yx, bx, cx, dx, n)

            !Catch NaN
            if (any(isnan(xi)) .or. any(isnan(yT)) .or. any(isnan(bT)) .or. any(isnan(cT)) .or. any(isnan(dT)) .or. any(isnan(yx)) .or. any(isnan(bx)) .or. any(isnan(cx)) .or. any(isnan(dx)))then
                !errval = -5532 !noch hinzufügen
                exit
            end if

            d_vap = xi(k)
            d_liq = xi(k+1)
            delta = (d_liq - d_vap)/20.d0
            ! do 10 more density controlled steps toward the critical point
            do j = 1, 10
                rho1 = d_vap + real(j)*delta
                rho2 = d_liq - real(j)*delta
                Temp = ispline(gl,rho1, xi, yT, bT, cT, dT, n)
                x_vap(1) = ispline(gl,rho1, xi, yx, bx, cx, dx, n)
                x_vap(2) = 1.d0 - x_vap(1)
                x_liq(1) = ispline(gl,rho2, xi, yx, bx, cx, dx, n)
                x_liq(2) = 1.d0 - x_liq(1)
                points = points + 1
                T_points(points) = Temp
                p_points(points) = press
                x_points(points, 1) = x_vap(1)
                x_points(points, 2) = x_vap(2)
                x_points(points, 3) = x_liq(1)
                x_points(points, 4) = x_liq(2)
                rhovap_points(points) = d_vap + real(j)*delta
                rholiq_points(points) = d_liq - real(j)*delta

            end do
            errval = 0
            exit
        elseif (save_calc == 3) then
            points = points - 1
        end if

    end do

    if ((points == maxi-n) .and. (save_calc == 1)) points = points - 1
    !press = pstart*1.01d0
    !molfractions(2) = 1.d0 - molfractions(1)
    !x_known = molfractions
    !call reduced_parameters_calc(gl,300.d0)

#IF DEFINED(_DEBUG) 
    !--------------------------------------------------------------------------------
    ! CHECK IF WRITING THE DATA TO A FILE IS WISHED
    !--------------------------------------------------------------------------------
    !111 continue
    if (fileout /= '') then
        open(unit = 13,file=fileout, status='unknown', action='write', iostat=error)
        if (error /= 0) return
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'! P-X DIAGRAM DATA FILE FOR THE BINARY SYSTEM'
        write(13,*)'! ', trim(gl%components(1)),'-',trim(gl%components(2))
        write(13,*)'!'
        write(13,'(A20, F8.3)')'! pressure [MPa]: ', press
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'!'
        write(13,*)'! p/MPa    ','T/K       ','x',trim(gl%components(1)),'    x',trim(gl%components(2)),'    y',trim(gl%components(1)), &
            &'    y',trim(gl%components(2)),'   dl/mol/l ', 'dv/mol/l'
        do i = 1, points
            write(13,1000) p_points(i), T_points(i), x_points(i, 3), x_points(i, 4), x_points(i, 1), x_points(i, 2), &
                & rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
        end do
        if (errMSG /= '') write(13,*)'! ', errMSG
        close(13)
    end if
#ENDIF

    end subroutine txdiag


    subroutine pxdiag_db(gl,Temp, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, errval)
    !--------------------------------------------------------------------------------
    ! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    ! TEMPERATURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    ! Simplified routine that only tries to calculate bubble and dew points at given
    ! overall compositions and temperatures.
    !--------------------------------------------------------------------------------
    ! A. Jäger, August 2017





    implicit none

    type(type_gl) :: gl

    !--------------------------------------------------------------------------------
    integer, parameter:: maxi = 410 ! maximum number of calculated points, lengt of the return vectors!
    double precision, intent(in):: Temp     ! specified temperature for the p-x diagram
    double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision, intent(out):: x_points(maxi,2)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_liq(1)
    character(255), intent(in):: fileout    ! optional path for output file
    integer:: points,errval                ! number of calculated points, error code
    !--------------------------------------------------------------------------------

    double precision, dimension(2):: psat_fluid   !Saturation pressures of fluid 1 and 2
    logical, dimension(2) :: supercritical
    double precision, dimension(maxi) :: T_points_bubble, T_points_dew
    double precision, dimension(maxi) :: p_points_bubble, p_points_dew
    double precision, dimension(maxi) :: rho_points_bubble, rho_points_dew
    double precision, dimension(maxi) :: x_points_bubble, x_points_dew
    integer:: points_dew, points_bubble

    !Variables for pure vapor pressure iteration
    double precision :: d_vap, d_liq   ![mol/m³]
    integer :: iFlash, nrsubst
    integer :: iter

    integer:: start_fluid, second_fluid   !Indicates which fluid has the lower vapor pressure. Start with this fluid and save other fluid in second_fluid

    double precision, dimension(30):: z_overall     !overall composition
    double precision:: step

    !Variables for mixture vapor pressure iteration
    double precision:: press, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: Nr_x_given, iPhase_try

    integer:: i, j, k
    logical:: err_in_dew, err_in_bub

    logical:: Bubble_pts
    !--------------------------------------------------------------------------------

    !Variable specifies if bubble or dew points should be calculated along the phase boundary
    Bubble_pts = .true.

    !If one component is below its triple point temperature, quit with error
    if (Temp < gl%ttp(1)) then
        errval = -9981
        return
    end if
    if (Temp < gl%ttp(2)) then
        errval = -9981
        return
    end if

    !Initialize variables
    supercritical = .false.
    points = 0
    points_dew = 1
    points_bubble = 1
    z_overall = 0.D0
    step = 1.D0 / 400.d0

    !--------------------------------------------------------------------------------
    ! STEP 1: PURE FLUID START VALUES
    !--------------------------------------------------------------------------------
    ! get the vapor pressures of the pure fluids
    psat_fluid = 0.d0

    !Check if fluid 1 is supercritical and if not, calculate its vapor pressure at given temperature
    if (temp < gl%tc(1)) then
        psat_fluid(1) = vp_eq(gl,Temp, 1)
    else
        supercritical(1) = .true.
    end if

    !Check if fluid 2 is supercritical and if not, calculate its vapor pressure at given temperature
    if (temp < gl%tc(2)) then
        psat_fluid(2) = vp_eq(gl,Temp, 2)
    else
        supercritical(2) = .true.
    end if

    if (supercritical(1) .and. supercritical(2)) then   !Both fluids supercritical, either no px-diagram calculation possible or if possible not implemented yet
        errval = -9981
        return
    elseif(supercritical(1) == .true.) then
        start_fluid = 2
        second_fluid = 1
    elseif(supercritical(2) == .true.) then
        start_fluid = 1
        second_fluid = 2
    else
        if (psat_fluid(1) < psat_fluid(2)) then
            start_fluid = 1
            second_fluid = 2
        else
            start_fluid = 2
            second_fluid = 1
        end if
    end if

    !start_fluid = 2
    !second_fluid = 1

    !start_fluid = 1
    !second_fluid = 2

    if (supercritical(1) == .false.) then
        !Get the vapor pressure of fluid 1
        !----------------------------------
        d_vap = 0.D0
        d_liq = 0.D0
        iFlash = 1
        nrsubst = 1
        call Flash_Pure_PhaseBoundary(gl,psat_fluid(1), Temp, d_vap, d_liq, iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 1) then
            points = points + 1
            T_points(points) = Temp
            p_points(points) = psat_fluid(1)
            x_points(points, 1) = 1.d0 !xvap1
            x_points(points, 2) = 1.d0 !xliq1
            rhovap_points(points) = d_vap
            rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if

    if (supercritical(2) == .false.) then
        !Get the vapor pressure of fluid 2
        !----------------------------------
        d_vap = 0.D0
        d_liq = 0.D0
        iFlash = 1
        nrsubst = 2
        call Flash_Pure_PhaseBoundary(gl,psat_fluid(2), Temp, d_vap, d_liq, iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 2) then
            points = points + 1
            T_points(points) = Temp
            p_points(points) = psat_fluid(2)
            x_points(points, 1) = 1.d0 !xvap1
            x_points(points, 2) = 1.d0 !xliq1
            rhovap_points(points) = d_vap
            rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if


    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! STEP 2: FIRST BUBBLE AND DEW POINT OF THE MIXTURE
    !--------------------------------------------------------------------------------

    if (Bubble_pts) then
        !Generate first composition close to start fluid
        z_overall(start_fluid) = 1.D0 - step
        z_overall(second_fluid) = 0.D0 + step

        do i = 1, 20    !Try up to 20 times calculating the first bubble point of the mixture
            !Calculate the first mixture bubble point at given composition
            !-------------------------------------------------------------
            gl%molfractions = z_overall
            call reduced_parameters_calc(gl,Temp)
            press = 0.D0
            x_vap = 0.D0
            x_liq = 0.D0
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            vapfrac = 0
            iFlash = 3      !Bubble point calculation at given temperature
            iPhase_try = 5  !Try VLE
            Nr_x_given = 0
            call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                & iPhase_try, Nr_x_given, errval, iter)
            if (errval /= 0) then
                !if (i == 10) then
                !    return
                !end if
                !z_overall(start_fluid) = 1.D0 - step/(i*2.D0)
                !z_overall(second_fluid) = 0.D0 + step/(i*2.D0)
            else
                points_bubble = points_bubble + 1
                T_points_bubble(points_bubble) = Temp
                p_points_bubble(points_bubble) = press
                rho_points_bubble(points_bubble) = gl%rho_liq
                x_points_bubble(points_bubble) = x_liq(1)
                T_points_dew (points_bubble) = Temp
                p_points_dew(points_bubble) = press
                rho_points_dew(points_bubble) = gl%rho_vap
                x_points_dew(points_bubble) = x_vap(1)
                !Save for output
                T_points(points_bubble) = Temp
                p_points(points_bubble) = press
                x_points(points_bubble, 1) = x_vap(1) !xvap1
                x_points(points_bubble, 2) = x_liq(1) !xliq1
                rhovap_points(points_bubble) = gl%rho_vap
                rholiq_points(points_bubble) = gl%rho_liq
                points = points_bubble
                exit
            end if
            !-------------------------------------------------------------
        end do

    else
        !Reset overall composition
        z_overall(start_fluid) = 1.D0 - step
        z_overall(second_fluid) = 0.D0 + step

        do i = 1, 20    !Try up to 20 times calculating the first dew point of the mixture
            !Calculate the first mixture dew point at given composition
            !-------------------------------------------------------------
            gl%molfractions = z_overall
            call reduced_parameters_calc(gl,Temp)
            press = 0.D0
            x_vap = 0.D0
            x_liq = 0.D0
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            vapfrac = 0
            iFlash = 4      !Dew point calculation at given temperature
            iPhase_try = 5  !Try VLE
            Nr_x_given = 0
            call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                & iPhase_try, Nr_x_given, errval, iter)

            if (errval /= 0) then
                !if (i == 10) then
                return
                !end if
                !z_overall(start_fluid) = 1.D0 - step/(i*2.D0)
                !z_overall(second_fluid) = 0.D0 + step/(i*2.D0)
            else
                points_dew = points_dew + 1
                T_points_dew(points_dew) = Temp
                p_points_dew(points_dew) = press
                rho_points_dew(points_dew) = gl%rho_vap
                x_points_dew(points_dew) = x_vap(1)
                T_points_bubble (points_dew) = Temp
                p_points_bubble(points_dew) = press
                rho_points_bubble(points_dew) = gl%rho_liq
                x_points_bubble(points_dew) = x_liq(1)
                !Save for output
                T_points(points_dew) = Temp
                p_points(points_dew) = press
                x_points(points_dew, 1) = x_vap(1) !xvap1
                x_points(points_dew, 2) = x_liq(1) !xliq1
                rhovap_points(points_dew) = gl%rho_vap
                rholiq_points(points_dew) = gl%rho_liq
                points = points_dew
                exit
            end if
        end do
        !-------------------------------------------------------------
    end if


    !--------------------------------------------------------------------------------
    ! STEP 3: MARCH ALONG THE PHASE BOUNDARIES BY CHANGING THE OVERALL COMPOSITION
    !--------------------------------------------------------------------------------
    k = 2
    do j = 2, 400

        if (bubble_pts) then
            !Generate first composition close to start fluid
            z_overall(start_fluid) = 1.D0 - k*step
            z_overall(second_fluid) = 0.D0 + k*step

            if (z_overall(start_fluid) < step) then
                exit
            end if
            err_in_dew = .false.
            err_in_bub = .false.

            do i = 1, 20    !Try up to 20 times calculating the next bubble point of the mixture
                !Calculate the first mixture bubble point at given composition
                !-------------------------------------------------------------
                gl%molfractions = z_overall
                call reduced_parameters_calc(gl,Temp)
                press = p_points_bubble(points_bubble)
                x_vap = 0.D0
                x_liq = z_overall
                x_vap(1) = x_points_dew(points_bubble)
                x_vap(2) = 1.D0 - x_points_dew(points_bubble)
                rhovap_est = rho_points_dew(points_bubble)
                rholiq_est = rho_points_bubble(points_bubble)
                vapfrac = 0
                iFlash = 3      !Bubble point calculation at given temperature
                iPhase_try = 5  !Try VLE
                Nr_x_given = 0
                call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                    & iPhase_try, Nr_x_given, errval, iter)
                if (errval /= 0) then
                    !err_in_bub = .true.
                    !if ((i == 20) .or. (dabs(z_overall(1) - x_points_bubble(points_bubble)) < 1.D-10) ) then
                    !    return
                    !end if
                    !z_overall(start_fluid) = z_overall(start_fluid) + step/20.D0
                    !z_overall(second_fluid) = z_overall(second_fluid) - step/20.D0
                    return
                else
                    points_bubble = points_bubble + 1
                    T_points_bubble(points_bubble) = Temp
                    p_points_bubble(points_bubble) = press
                    rho_points_bubble(points_bubble) = gl%rho_liq
                    x_points_bubble(points_bubble) = x_liq(1)
                    T_points_dew (points_bubble) = Temp
                    p_points_dew(points_bubble) = press
                    rho_points_dew(points_bubble) = gl%rho_vap
                    x_points_dew(points_bubble) = x_vap(1)
                    !Save for output
                    T_points(points_bubble) = Temp
                    p_points(points_bubble) = press
                    x_points(points_bubble, 1) = x_vap(1) !xvap1
                    x_points(points_bubble, 2) = x_liq(1) !xliq1
                    rhovap_points(points_bubble) = gl%rho_vap
                    rholiq_points(points_bubble) = gl%rho_liq
                    points = points_bubble
                    exit
                end if
                !-------------------------------------------------------------
            end do

        else

            !Reset overall composition
            z_overall(start_fluid) = 1.D0 - k*step
            z_overall(second_fluid) = 0.D0 + k*step

            err_in_dew = .false.
            err_in_bub = .false.

            do i = 1, 20    !Try up to 20 times calculating the next dew point of the mixture
                !Calculate the first mixture dew point at given composition
                !-------------------------------------------------------------
                gl%molfractions = z_overall
                call reduced_parameters_calc(gl,Temp)
                press = p_points_dew(points_dew)
                x_vap = z_overall
                x_liq = 0.D0
                x_liq(1) = x_points_bubble(points_dew)
                x_liq(2) = 1.D0 - x_points_bubble(points_dew)
                rhovap_est = rho_points_dew(points_dew)
                rholiq_est = rho_points_bubble(points_dew)
                vapfrac = 0
                iFlash = 4      !Dew point calculation at given temperature
                iPhase_try = 5  !Try VLE
                Nr_x_given = 0
                call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                    & iPhase_try, Nr_x_given, errval, iter)

                if (errval /= 0) then
                    !err_in_dew = .true.
                    !if ((i == 20) .or. (dabs(z_overall(1) - x_points_dew(points_dew)) < 1.D-10)) then
                    !    return
                    !end if
                    !z_overall(start_fluid) = z_overall(start_fluid) + step/20.D0
                    !z_overall(second_fluid) = z_overall(second_fluid) - step/20.D0
                    return
                else
                    points_dew = points_dew + 1
                    T_points_dew(points_dew) = Temp
                    p_points_dew(points_dew) = press
                    rho_points_dew(points_dew) = gl%rho_vap
                    x_points_dew(points_dew) = x_vap(1)
                    T_points_bubble (points_dew) = Temp
                    p_points_bubble(points_dew) = press
                    rho_points_bubble(points_dew) = gl%rho_liq
                    x_points_bubble(points_dew) = x_liq(1)
                    !Save for output
                    T_points(points_dew) = Temp
                    p_points(points_dew) = press
                    x_points(points_dew, 1) = x_vap(1) !xvap1
                    x_points(points_dew, 2) = x_liq(1) !xliq1
                    rhovap_points(points_dew) = gl%rho_vap
                    rholiq_points(points_dew) = gl%rho_liq
                    points = points_dew
                    exit
                end if
            end do

        end if

        !Check whether an error occured in the bubble or dew point calculation. If yes, repeat
        if (err_in_dew .or. err_in_bub) then
            k = k   !Redo step
        else
            k = k+1
        end if

    end do
    !--------------------------------------------------------------------------------

    !pause

    end subroutine


    subroutine Txdiag_db(gl,press, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, errval)
    !--------------------------------------------------------------------------------
    ! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    ! PRESSURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    ! Simplified routine that only tries to calculate bubble and dew points at given
    ! overall compositions and pressures.
    !--------------------------------------------------------------------------------
    ! A. Jäger, January 2018





    implicit none

    type(type_gl) :: gl

    !--------------------------------------------------------------------------------
    integer, parameter:: maxi = 300 ! maximum number of calculated points, lengt of the return vectors!
    double precision, intent(in):: press     ! specified temperature for the p-x diagram
    double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision, intent(out):: x_points(maxi,2)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_liq(1)
    character(255), intent(in):: fileout    ! optional path for output file
    integer:: points,errval                ! number of calculated points, error code
    !--------------------------------------------------------------------------------

    double precision, dimension(2):: Tsat_fluid   !Saturation temperatures of fluid 1 and 2
    logical, dimension(2) :: supercritical
    double precision, dimension(maxi) :: T_points_bubble, T_points_dew
    double precision, dimension(maxi) :: p_points_bubble, p_points_dew
    double precision, dimension(maxi) :: rho_points_bubble, rho_points_dew
    double precision, dimension(maxi) :: x_points_bubble, x_points_dew
    integer:: points_dew, points_bubble

    !Variables for pure vapor pressure iteration
    double precision :: d_vap, d_liq   ![mol/m³]
    integer :: iFlash, nrsubst
    integer :: iter

    integer:: start_fluid, second_fluid   !Indicates which fluid has the lower vapor pressure. Start with this fluid and save other fluid in second_fluid

    double precision, dimension(30):: z_overall     !overall composition
    double precision:: step

    !Variables for mixture vapor pressure iteration
    double precision:: Temp, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: Nr_x_given, iPhase_try

    integer:: i, j, k
    logical:: err_in_dew, err_in_bub
    logical:: Bubble_pts
    !--------------------------------------------------------------------------------

    !Variable specifies if bubble or dew points should be calculated along the phase boundary
    Bubble_pts = .false.

    !If one component is below its triple point pressure, quit with error
    if (press < gl%ptp(1)) then
        errval = -9981
        return
    end if
    if (press < gl%ptp(2)) then
        errval = -9981
        return
    end if

    !Initialize variables
    supercritical = .false.
    points = 0
    points_dew = 1
    points_bubble = 1
    z_overall = 0.D0
    step = 1.D0 / 400.d0

    !--------------------------------------------------------------------------------
    ! STEP 1: PURE FLUID START VALUES
    !--------------------------------------------------------------------------------
    ! get the vapor pressures of the pure fluids
    Tsat_fluid = 0.d0

    !Check if fluid 1 is supercritical and if not, calculate its boiling temperature at given pressure
    if (press < gl%pc(1)) then
        Tsat_fluid(1) = Estimate_Tsat(gl,press, 1)
    else
        supercritical(1) = .true.
    end if

    !Check if fluid 2 is supercritical and if not, calculate its boiling temperature at given pressure
    if (press < gl%pc(2)) then
        Tsat_fluid(2) = Estimate_Tsat(gl,press, 2)
    else
        supercritical(2) = .true.
    end if

    if (supercritical(1) .and. supercritical(2)) then   !Both fluids supercritical, either no px-diagram calculation possible or if possible not implemented yet
        errval = -9981
        return
    elseif(supercritical(1) == .true.) then
        start_fluid = 2
        second_fluid = 1
    elseif(supercritical(2) == .true.) then
        start_fluid = 1
        second_fluid = 2
    else
        if (Tsat_fluid(1) < Tsat_fluid(2)) then
            start_fluid = 2
            second_fluid = 1
        else
            start_fluid = 1
            second_fluid = 2
        end if
    end if

    !start_fluid = 2
    !second_fluid = 1

    !start_fluid = 1
    !second_fluid = 2

    if (supercritical(1) == .false.) then
        !Get the boiling temperature of fluid 1
        !----------------------------------
        d_vap = 0.D0
        d_liq = 0.D0
        iFlash = 2
        nrsubst = 1
        call Flash_Pure_PhaseBoundary(gl,press, Tsat_fluid(1), d_vap, d_liq, iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 1) then
            points = points + 1
            T_points(points) = Tsat_fluid(1)
            p_points(points) = press
            x_points(points, 1) = 1.d0 !xvap1
            x_points(points, 2) = 1.d0 !xliq1
            rhovap_points(points) = d_vap
            rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if

    if (supercritical(2) == .false.) then
        !Get the vapor pressure of fluid 2
        !----------------------------------
        d_vap = 0.D0
        d_liq = 0.D0
        iFlash = 2
        nrsubst = 2
        call Flash_Pure_PhaseBoundary(gl,press, Tsat_fluid(2), d_vap, d_liq, iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 2) then
            points = points + 1
            T_points(points) = Tsat_fluid(2)
            p_points(points) = press
            x_points(points, 1) = 1.d0 !xvap1
            x_points(points, 2) = 1.d0 !xliq1
            rhovap_points(points) = d_vap
            rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if


    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! STEP 2: FIRST BUBBLE AND DEW POINT OF THE MIXTURE
    !--------------------------------------------------------------------------------

    if (Bubble_pts) then
        !Generate first composition close to start fluid
        z_overall(start_fluid) = 1.D0 - step
        z_overall(second_fluid) = 0.D0 + step

        do i = 1, 20    !Try up to 20 times calculating the first bubble point of the mixture
            !Calculate the first mixture bubble point at given composition
            !-------------------------------------------------------------
            gl%molfractions = z_overall
            call reduced_parameters_calc(gl,Temp)
            Temp = 0.D0
            x_vap = 0.D0
            x_liq = 0.D0
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            vapfrac = 0
            iFlash = 1      !Bubble point calculation at given temperature
            iPhase_try = 5  !Try VLE
            Nr_x_given = 0
            call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                & iPhase_try, Nr_x_given, errval, iter)
            if (errval /= 0) then
                if (i == 10) then
                    return
                end if
                z_overall(start_fluid) = 1.D0 - step/(i*2.D0)
                z_overall(second_fluid) = 0.D0 + step/(i*2.D0)
            else
                points_bubble = points_bubble + 1
                T_points_bubble(points_bubble) = Temp
                p_points_bubble(points_bubble) = press
                rho_points_bubble(points_bubble) = gl%rho_liq
                x_points_bubble(points_bubble) = x_liq(1)
                T_points_dew (400+points_bubble) = Temp
                p_points_dew(400+points_bubble) = press
                rho_points_dew(400+points_bubble) = gl%rho_vap
                x_points_dew(400+points_bubble) = x_vap(1)
                !Save for output
                T_points(points_bubble) = Temp
                p_points(points_bubble) = press
                x_points(points_bubble, 1) = x_vap(1) !xvap1
                x_points(points_bubble, 2) = x_liq(1) !xliq1
                rhovap_points(points_bubble) = gl%rho_vap
                rholiq_points(points_bubble) = gl%rho_liq
                points = points_bubble
                exit
            end if
            !-------------------------------------------------------------
        end do

    else

        !Reset overall composition
        z_overall(start_fluid) = 1.D0 - step
        z_overall(second_fluid) = 0.D0 + step

        do i = 1, 20    !Try up to 20 times calculating the first dew point of the mixture
            !Calculate the first mixture dew point at given composition
            !-------------------------------------------------------------
            gl%molfractions = z_overall
            call reduced_parameters_calc(gl,Temp)
            temp = 0.D0
            x_vap = 0.D0
            x_liq = 0.D0
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            vapfrac = 0
            iFlash = 2      !Dew point calculation at given temperature
            iPhase_try = 5  !Try VLE
            Nr_x_given = 0
            call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                & iPhase_try, Nr_x_given, errval, iter)

            if (errval /= 0) then
                if (i == 10) then
                    return
                end if
                z_overall(start_fluid) = 1.D0 - step/(i*2.D0)
                z_overall(second_fluid) = 0.D0 + step/(i*2.D0)
            else
                points_dew = points_dew + 1
                T_points_dew(points_dew) = Temp
                p_points_dew(points_dew) = press
                rho_points_dew(points_dew) = gl%rho_vap
                x_points_dew(points_dew) = x_vap(1)
                T_points_bubble (400+points_dew) = Temp
                p_points_bubble(400+points_dew) = press
                rho_points_bubble(400+points_dew) = gl%rho_liq
                x_points_bubble(400+points_dew) = x_liq(1)
                exit
            end if
        end do
        !-------------------------------------------------------------
    end if


    !--------------------------------------------------------------------------------
    ! STEP 3: MARCH ALONG THE PHASE BOUNDARIES BY CHANGING THE OVERALL COMPOSITION
    !--------------------------------------------------------------------------------
    k = 2
    do j = 2, 500

        if (bubble_pts) then
            !Generate first composition close to start fluid
            z_overall(start_fluid) = 1.D0 - k*step
            z_overall(second_fluid) = 0.D0 + k*step

            if (z_overall(start_fluid) < step) then
                exit
            end if
            err_in_dew = .false.
            err_in_bub = .false.

            do i = 1, 20    !Try up to 20 times calculating the next bubble point of the mixture
                !Calculate the first mixture bubble point at given composition
                !-------------------------------------------------------------
                gl%molfractions = z_overall
                call reduced_parameters_calc(gl,Temp)
                Temp = T_points_bubble(points_bubble)
                x_vap = 0.D0
                x_liq = z_overall
                x_vap(1) = x_points_dew(400+points_bubble)
                x_vap(2) = 1.D0 - x_points_dew(400+points_bubble)
                rhovap_est = rho_points_dew(400+points_bubble)
                rholiq_est = rho_points_bubble(points_bubble)
                vapfrac = 0
                iFlash = 1      !Bubble point calculation at given temperature
                iPhase_try = 5  !Try VLE
                Nr_x_given = 0
                call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                    & iPhase_try, Nr_x_given, errval, iter)
                if (errval /= 0) then
                    !err_in_bub = .true.
                    !if ((i == 20) .or. (dabs(z_overall(1) - x_points_bubble(points_bubble)) < 1.D-10) ) then
                    !    return
                    !end if
                    !z_overall(start_fluid) = z_overall(start_fluid) + step/20.D0
                    !z_overall(second_fluid) = z_overall(second_fluid) - step/20.D0
                    return
                else
                    points_bubble = points_bubble + 1
                    T_points_bubble(points_bubble) = Temp
                    p_points_bubble(points_bubble) = press
                    rho_points_bubble(points_bubble) = gl%rho_liq
                    x_points_bubble(points_bubble) = x_liq(1)
                    T_points_dew (400+points_bubble) = Temp
                    p_points_dew(400+points_bubble) = press
                    rho_points_dew(400+points_bubble) = gl%rho_vap
                    x_points_dew(400+points_bubble) = x_vap(1)
                    !Save for output
                    T_points(points_bubble) = Temp
                    p_points(points_bubble) = press
                    x_points(points_bubble, 1) = x_vap(1) !xvap1
                    x_points(points_bubble, 2) = x_liq(1) !xliq1
                    rhovap_points(points_bubble) = gl%rho_vap
                    rholiq_points(points_bubble) = gl%rho_liq
                    points = points_bubble
                    exit
                end if
                !-------------------------------------------------------------
            end do

        else
            !Reset overall composition
            z_overall(start_fluid) = 1.D0 - j*step
            z_overall(second_fluid) = 0.D0 + j*step

            do i = 1, 20    !Try up to 20 times calculating the next dew point of the mixture
                !Calculate the first mixture dew point at given composition
                !-------------------------------------------------------------
                gl%molfractions = z_overall
                call reduced_parameters_calc(gl,Temp)
                Temp = T_points_dew(points_dew)
                x_vap = z_overall
                x_liq = 0.D0
                x_liq(1) = x_points_bubble(400+points_dew)
                x_liq(2) = 1.D0 - x_points_bubble(400+points_dew)
                rhovap_est = rho_points_dew(points_dew)
                rholiq_est = rho_points_bubble(400+points_dew)
                vapfrac = 0
                iFlash = 2      !Dew point calculation at given temperature
                iPhase_try = 5  !Try VLE
                Nr_x_given = 0
                call Flash_PhaseBoundary_calc(gl, press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                    & iPhase_try, Nr_x_given, errval, iter)

                if (errval /= 0) then
                    !err_in_dew = .true.
                    !if ((i == 20) .or. (dabs(z_overall(1) - x_points_dew(points_dew)) < 1.D-10)) then
                    !    return
                    !end if
                    !z_overall(start_fluid) = z_overall(start_fluid) + step/20.D0
                    !z_overall(second_fluid) = z_overall(second_fluid) - step/20.D0
                    return
                else
                    points_dew = points_dew + 1
                    T_points_dew(points_dew) = Temp
                    p_points_dew(points_dew) = press
                    rho_points_dew(points_dew) = gl%rho_vap
                    x_points_dew(points_dew) = x_vap(1)
                    T_points_bubble (400+points_dew) = Temp
                    p_points_bubble(400+points_dew) = press
                    rho_points_bubble(400+points_dew) = gl%rho_liq
                    x_points_bubble(400+points_dew) = x_liq(1)
                    !Save for output
                    T_points(points_dew) = Temp
                    p_points(points_dew) = press
                    x_points(points_dew, 1) = x_vap(1) !xvap1
                    x_points(points_dew, 2) = x_liq(1) !xliq1
                    rhovap_points(points_dew) = gl%rho_vap
                    rholiq_points(points_dew) = gl%rho_liq
                    points = points_dew
                    exit
                end if
            end do

        end if

        !Check whether an error occured in the bubble or dew point calculation. If yes, repeat
        if (err_in_dew .or. err_in_bub) then
            k = k   !Redo step
        else
            k = k+1
        end if

    end do
    !--------------------------------------------------------------------------------

    !pause

    end subroutine



    !Andreas Jäger, October 2018
    !Routine for the calculation of all needed derivatives for isochoric tracing of isotherms (px-diagrams)
    subroutine Derivs_isoch_therm(gl, isoch_therm, Temp, Dens, errval)
    !--------------------------------------------------------------------------------
    ! This subroutine calculates properties of mixture models formulated
    ! in the dimensionless Helmholtz energy alpha, which are needed for
    ! tracing of isotherms in isochoric thermodynamics
    ! The algorithm and method is described in:
    ! Bell, I.H.; Deiters U.K.: "On the Construction of Binary Mixture p-x and T-x Diagrams from Isochoric Thermodynamics"
    !
    ! Input:    Temperature     [K]
    !           Density         [mol/m³]
    !           Mole fractions  (contained in type gl variable)
    ! Output:   All properties contained in type type_isoch_therm
    !           errval: specifies if an error occured during the calculations
    !--------------------------------------------------------------------------------
    ! A. Jäger, October 2018




    implicit none

    !Input and output variables of the routine
    type(type_gl) :: gl
    type(type_isoch_therm) :: isoch_therm
    double precision:: Temp
    double precision:: Dens
    integer:: errval

    !Partial densities
    double precision, dimension(:), allocatable :: rho_i

    !Variables for FNRDERIVS, MIXDERIVSFNR, and MIXDERIVSFNI
    integer, dimension(nderivs):: GETDER                            ! array specifier to indicate, which derivative is needed for the residual part
    double precision, dimension(:,:), allocatable :: FNRDER         ! array with the computed values for the derivatives of the residual part (pure fluids under mixture conditions)
    double precision, dimension(nderivs)::MIXDERFNR                 ! array with the computed values for the derivatives of the residual part (mixture)
    double precision, dimension(nderivs)::DEPFUNCDER                ! array with the computed values for the derivatives of the departure function
    double precision, dimension(:,:,:), allocatable :: DEPFUNCBIN     ! array with the computed values for the derivatives of the binary specific departure function
    integer, dimension(nderivsi)::GETDERI                                 ! array specifier to indicate which derivative is needed for the ideal part
    double precision, dimension(nderivsi)::MIXDERFNI                      ! array with computed values for the derivatives of the ideal part

    !Variables for the derivatives of the reducing functions with respect to mole fractions
    !First derivative with respect to xi
    !double precision, dimension(30)::dTred_dxi, drhored_dxi
    double precision, dimension(:), allocatable :: dvred_dxi
    integer::i, j, m, n
    double precision:: xi, xj, rhoci, rhocj, tci, tcj
    ! parameters for the 1st loop
    double precision:: betatji, betarhoji, gammatji, gammarhoji
    double precision:: ct_ji, crho_ji
    double precision:: dftji_dxi, dfrhoji_dxi
    ! parameters for the 2nd loop
    double precision:: betatij, betarhoij, gammatij, gammarhoij
    double precision:: ct_ij, crho_ij
    double precision:: dftij_dxi, dfrhoij_dxi
    !Second derivative with respect to xi and xj
    !double precision, dimension(30,30):: d2Tred_dxidxj, d2rhored_dxidxj
    double precision, dimension(:,:), allocatable :: d2vred_dxidxj
    double precision::d2fvij_dxidxj, d2fTij_dxidxj, cv_ij, help
    !double precision, dimension(30)::d2Tred_dxi2, d2rhored_dxi2   ! Return vectors, maximum 30 components
    !double precision, dimension(30)::d2vred_dxi2                    ! internal vector for vred
    ! parameters for the 1st loop
    double precision:: d2ftji_dxi2, d2frhoji_dxi2
    ! parameters for the 2nd loop
    double precision:: d2ftij_dxi2, d2frhoij_dxi2

    !Variables for the derivatives of the reducing functions with respect to rho_i and rho_j
    double precision, dimension(:,:), allocatable :: kron_delta
    double precision, dimension(:,:), allocatable :: dxm_drhoi            !Help variable
    double precision, dimension(:,:,:), allocatable :: d2xm_drhoi_drhoj   !Help variable

    !Variables for (preliminary) subroutine calls in order to get the derivatives of the reduced residual Helmholtz energy with respect to compositions
    double precision, dimension(:), allocatable :: DERIVFXDEL, DERIVFXTAU
    double precision, dimension(:), allocatable :: DERIVFX, DERIVFX2
    double precision, dimension(:,:), allocatable :: DERIVFXiXj
    double precision, dimension (nderivs,gl%ncomp) :: DERIVFX_ALL
    double precision, dimension (nderivs,gl%ncomp,gl%ncomp) :: DERIVFX2_ALL


    !Allocate variables
    allocate(rho_i(gl%ncomp))
    allocate(dvred_dxi(gl%ncomp))
    allocate(d2vred_dxidxj(gl%ncomp,gl%ncomp))
    allocate(kron_delta(gl%ncomp,gl%ncomp))
    allocate(dxm_drhoi(gl%ncomp,gl%ncomp))
    allocate(d2xm_drhoi_drhoj(gl%ncomp,gl%ncomp,gl%ncomp))
    allocate(FNRDER(nderivs,gl%ncomp))
    allocate(DEPFUNCBIN(nderivs,gl%ncomp,gl%ncomp))

    !Allocate preliminarily used variables for subroutine calls
    allocate(DERIVFXDEL(30))
    allocate(DERIVFXTAU(30))
    allocate(DERIVFX(30))
    allocate(DERIVFX2(30))
    allocate(DERIVFXiXj(30,30))

    !Initialize variables
    MIXDERFNR = 0.D0
    isoch_therm%alpha_it = 0.D0
    isoch_therm%alpha0_it = 0.D0
    isoch_therm%alphar_it = 0.D0
    isoch_therm%dalphar_dxi_it = 0.D0
    isoch_therm%d2alphar_dxiddel_it = 0.D0
    isoch_therm%d2alphar_dxidtau_it = 0.D0
    isoch_therm%d2alphar_dxidxj_it = 0.D0

    !Check if the specified mixture model is implemented and can be calculated with this routine
    !At the moment, only the standard multi-fluid mixture model (mix_type == 1) and PCP-SAFT (mix_type == 6) work here, since for the other models the required variables are not available
    if ((gl%mix_type == 11) .or. (gl%mix_type == 12) .or. (gl%mix_type == 13) .or. (gl%mix_type == 2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22) .or. (gl%mix_type == 3) .or. (gl%mix_type == 31) .or. (gl%mix_type == 4)) then
        errval = -20000 !Internal error in isochoric tracing at the moment
        return
    end if

    !Fill the Kronecker delta
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            if (j == i) then
                kron_delta(i,j) = 1.D0
            else
                kron_delta(i,j) = 0.D0
            end if
        end do
    end do

    !Calculate energy densities for all components in the phase
    do i = 1, gl%ncomp
        rho_i(i) = Dens * gl%molfractions(i)
    end do

    !Calculate the inverse reduced temperature tau
    isoch_therm%tau_it = gl%tredmix / Temp
    !Calculate the reduced density delta
    isoch_therm%delta_it = Dens / gl%rhoredmix

    !Get the residual dimensionless Helmholtz energies and derivatives with respect to delta and tau
    !------------------------------------------------------------------------------------------------
    GETDER = 0
    GETDER(1) = 1   !alphar
    GETDER(2) = 1   !dalphar_ddelta         * delta
    GETDER(3) = 1   !d2alphar_ddelta2       * delta^2
    GETDER(4) = 1   !dalphar_dtau           * tau
    GETDER(5) = 1   !d2alphar_dtau2         * tau^2
    GETDER(6) = 1   !d2alphar_dtauddelta    * tau * delta
    !For multi-fluid mixture models, calculate the pure fluid parts first and then use them to get the mixture properties.
    if ((gl%Mix_type == 1) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then !Multi-fluid mixture model and models based on this model (MIXTURE MODELS 11 and 19 NOT YET IMPLEMENTED)

        !sum of the residual Helmholtz energies and its derivatives of the pure fluids
        FNRDER = 0.D0
        DEPFUNCDER = 0.d0
        DEPFUNCBIN = 0.d0
        do i = 1, gl%NCOMP

            call FNRDERIVS(gl,Temp, Dens, GETDER, FNRDER(:,i), i) ! call the routine for the derivatives of fluid i
            !Calculate the derivatives of the residual part for the mixture (residual Helmholtz energy and all first and second derivatives)
            do j = 1, 6
                MIXDERFNR(j) = MIXDERFNR(j) + gl%MOLFRACTIONS(i)*FNRDER(j,i) ! sum of derivatives over all pure fluids, multiplied with the respective mole fraction
            end do

        end do

        !Calculate the departure function
        if (gl%Mix_type == 1) then  !Usual multi-fluid mixture model

            !sum of the binary departure functions
            do i = 1, gl%ncomp - 1
                do j = i + 1, gl%ncomp
                    if (dabs(gl%Fij(i, j)) > 1.d-14) then    ! check if departure function exists for this binary combination
                        call DEPFUNCFNR (gl, Temp, Dens, GETDER, DEPFUNCBIN(:,i,j), i, j)   ! calls the routine for the calculation of the binary specific departure function and its derivatives
                        DEPFUNCBIN(:,j,i) = DEPFUNCBIN(:,i,j)
                        do m = 1, 6                      ! loop through all derivatives returned from DEPFUNCFNR
                            DEPFUNCDER(m) = DEPFUNCDER(m) + gl%molfractions(i) * gl%molfractions(j) * gl%Fij(i, j) * DEPFUNCBIN(m,i,j)     ! sum over all binary departure functions, multiplied with the mole fractions and the weighing factor Fij
                        end do
                    end if
                end do
            end do

        elseif ((gl%Mix_type == 12) .or. (gl%Mix_type == 13)) then !Multi-fluid mixture model with UNIFAC or COSMO-SAC in the departure function

            !excess based departure function and derivatives
            call DEPFUNC_GE_BASED (gl,Temp, Dens, GETDER, DEPFUNCDER)

        end if

        ! add the two parts of the residual Helmholtz energy for mixtures
        do m = 1, 6
            MIXDERFNR(m) = MIXDERFNR(m) + DEPFUNCDER(m)
        end do

    else    !All other models (Cubics, LKP, PCP-SAFT, etc). Directly get the mixture quantity because the pure fluid quantities are not needed for the composition derivatives

        call MIXDERIVSFNR (gl,Temp, Dens, GETDER, MIXDERFNR)

    end if
    isoch_therm%alphar_it = MIXDERFNR
    !Test
    !call MIXDERIVSFNR (gl,Temp, Dens, GETDER, MIXDERFNR)
    !------------------------------------------------------------------------------------------------

    !Get the ideal dimensionless Helmholtz energies and derivatives with respect to delta and tau
    !------------------------------------------------------------------------------------------------
    GETDERI = 0
    GETDERI(1) = 1   !alpha0
    GETDERI(2) = 1   !dalpha0_ddelta         * delta
    GETDERI(3) = 1   !d2alpha0_ddelta2       * delta^2
    GETDERI(4) = 1   !d2alpha0_dtauddelta    * tau * delta
    GETDERI(5) = 1   !dalpha0_dtau           * tau
    GETDERI(6) = 1   !d2alpha0_dtau2         * tau^2
    call MIXDERIVSFNI(gl,Temp, Dens, GETDERI, MIXDERFNI)
    isoch_therm%alpha0_it(1:6) = MIXDERFNI(1:6)
    !------------------------------------------------------------------------------------------------

    !Calculate the dimensionless Helmholtz energies and derivatives with respect to delta and tau
    !------------------------------------------------------------------------------------------------
    isoch_therm%alpha_it(1) = isoch_therm%alpha0_it(1) + isoch_therm%alphar_it(1)
    isoch_therm%alpha_it(2) = isoch_therm%alpha0_it(2) + isoch_therm%alphar_it(2)
    isoch_therm%alpha_it(3) = isoch_therm%alpha0_it(3) + isoch_therm%alphar_it(3)
    isoch_therm%alpha_it(4) = isoch_therm%alpha0_it(5) + isoch_therm%alphar_it(4)
    isoch_therm%alpha_it(5) = isoch_therm%alpha0_it(6) + isoch_therm%alphar_it(5)
    isoch_therm%alpha_it(6) = isoch_therm%alpha0_it(4) + isoch_therm%alphar_it(6)
    !------------------------------------------------------------------------------------------------




    !Calculate the derivative of the reducing function Tr and rhor with respect to xi and xj
    !------------------------------------------------------------------------------------------------
    !dTred_dxi = 0.D0
    !drhored_dxi = 0.D0
    dvred_dxi = 0.D0
    !d2Tred_dxidxj = 0.D0
    d2vred_dxidxj = 0.D0
    !d2rhored_dxidxj = 0.D0

    If ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then     !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016),  Mixtype = 12 and Mixtype = 13 are for excess based departure function mixing rules

        !First derivatives with respect to xi
        !-------------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, gl%NCOMP
            rhoci = gl%rhoc(i)
            tci = gl%tc(i)
            xi = gl%MOLFRACTIONS(i)
            isoch_therm%dTr_dxi_it(i) = 2.D0*xi*tci
            dvred_dxi(i) = 2.D0*xi/rhoci
            ! 1st loop for j < i
            do j = 1, i-1
                xj = gl%MOLFRACTIONS(j)
                rhocj = gl%rhoc(j)
                tcj = gl%tc(j)
                betatji = gl%RFBETAT(j, i)
                betarhoji = gl%RFBETARHO(j, i)
                gammatji = gl%RFGAMMAT(j, i)
                gammarhoji = gl%RFGAMMARHO(j, i)
                ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
                crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
                dftji_dxi = xj*(xj + xi)/(betatji**2*xj + xi) + &
                    & xj*xi/(betatji**2*xj + xi)*(1.D0 - (xj + xi)/(betatji**2*xj + xi))
                !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi:
                dfrhoji_dxi = xj*(xj + xi)/(betarhoji**2*xj + xi) + &
                    & xj*xi/(betarhoji**2*xj + xi)*(1.D0 - (xj + xi)/(betarhoji**2*xj + xi))
                ! 1st derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
                isoch_therm%dTr_dxi_it(i) = isoch_therm%dTr_dxi_it(i) +  ct_ji * dftji_dxi
                ! 1st derivative of rhored w.r.t xi - 1st summation for j = 1 to i-1
                dvred_dxi(i) = dvred_dxi(i) +  crho_ji * dfrhoji_dxi
            end do
            !2nd loop for j > i
            do j = i + 1, gl%NCOMP
                xj = gl%MOLFRACTIONS(j)
                rhocj = gl%rhoc(j)
                tcj = gl%tc(j)
                betatij = gl%RFBETAT(i, j)
                betarhoij = gl%RFBETARHO(i, j)
                gammatij = gl%RFGAMMAT(i, j)
                gammarhoij = gl%RFGAMMARHO(i, j)
                ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
                crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
                dftij_dxi = xj*(xi + xj)/(betatij**2*xi + xj) + &
                    & xi*xj/(betatij**2*xi + xj)*(1 - betatij**2*(xi + xj)/(betatij**2*xi + xj))
                !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi:
                dfrhoij_dxi = xj*(xi + xj)/(betarhoij**2*xi + xj) + &
                    & xi*xj/(betarhoij**2*xi + xj)*(1 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj))
                ! 1st derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
                isoch_therm%dTr_dxi_it(i) = isoch_therm%dTr_dxi_it(i) +  ct_ij * dftij_dxi
                ! 1st derivative of rhored w.r.t xi - 2nd summation for j = i+1 to ncomp
                dvred_dxi(i) = dvred_dxi(i) +  crho_ij * dfrhoij_dxi
            end do
            isoch_therm%drhor_dxi_it(i) = - gl%rhoredmix**2*dvred_dxi(i) ! drhored/dxi = -rhored^2*(d(1/rhored)dxi) - the routine above actually calculates dvred/dxi ...
        end do
        !-------------------------------------------------------------------------------------------------------------------------------------------

        !Second derivatives with respect to xi and xj
        !-------------------------------------------------------------------------------------------------------------------------------------------
        !Mixed derivatives
        do i = 1, (gl%ncomp - 1)
            do j = (i + 1), gl%ncomp
                betatij = gl%RFBETAT(i,j)
                betarhoij = gl%RFBETARHO(i,j)
                gammatij = gl%RFGAMMAT(i,j)
                gammarhoij = gl%RFGAMMARHO(i,j)
                xi = gl%MOLFRACTIONS(i)
                xj = gl%MOLFRACTIONS(j)
                tci = gl%tc(i)
                tcj = gl%tc(j)
                rhoci = gl%rhoc(i)
                rhocj = gl%rhoc(j)
                !critical term, temperature
                cT_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0
                !critical term, volume
                cv_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci**(1.D0/3.D0)) + 1.D0/(rhocj**(1.D0/3.D0)))**3
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
                d2fTij_dxidxj = (xi + xj)/(betatij**2*xi + xj) + xj/(betatij**2*xi + xj)*(1.D0 - (xi + xj)/(betatij**2*xi + xj)) &
                    & + xi/(betatij**2*xi + xj)*(1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) &
                    & - xi*xj/((betatij**2*xi + xj)**2)*(1.D0 + betatij**2 - 2.D0*betatij**2*(xi + xj)/(betatij**2*xi + xj))
                !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi and xj:
                help = 1.d0/(betarhoij**2*xi + xj)
                d2fvij_dxidxj = (xi + xj)*help &
                    & + xj*help*(1.D0 - (xi + xj)*help) &
                    & + xi*help*(1.D0 - betarhoij**2*(xi + xj)*help) &
                    & - xi*xj/((betarhoij**2*xi + xj)**2)*(1.D0 + betarhoij**2 &
                    & - 2.D0*betarhoij**2*(xi + xj)*help)
                !2nd derivative of Tred w.r.t xi and xj
                isoch_therm%d2Tr_dxidxj_it(i,j) = cT_ij*d2fTij_dxidxj
                isoch_therm%d2Tr_dxidxj_it(j,i) = isoch_therm%d2Tr_dxidxj_it(i,j)
                !2nd derivative of vred w.r.t xi and xj
                d2vred_dxidxj(i,j) = cv_ij*d2fvij_dxidxj
                d2vred_dxidxj(j,i) = d2vred_dxidxj(i,j)
                ! calculate d^2(rhored)/d(xi)d(xj) from d^2(1/rhored)/d(xi)d(xj)
                isoch_therm%d2rhor_dxidxj_it(i,j) = 2.D0/gl%rhoredmix*isoch_therm%drhor_dxi_it(i)*isoch_therm%drhor_dxi_it(j) - gl%rhoredmix**2*d2vred_dxidxj(i,j)
                isoch_therm%d2rhor_dxidxj_it(j,i) = isoch_therm%d2rhor_dxidxj_it(i,j)
            end do
        end do
        !Derivatives with respect to xi^2
        do i = 1, gl%NCOMP
            rhoci = gl%rhoc(i)
            tci = gl%tc(i)
            xi = gl%MOLFRACTIONS(i)
            isoch_therm%d2Tr_dxidxj_it(i,i) = 2.D0*tci
            d2vred_dxidxj(i,i) = 2.D0/rhoci
            ! 1st loop for j < i
            do j = 1, i-1
                xj = gl%MOLFRACTIONS(j)
                rhocj = gl%rhoc(j)
                tcj = gl%tc(j)
                betatji = gl%RFBETAT(j, i)
                betarhoji = gl%RFBETARHO(j, i)
                gammatji = gl%RFGAMMAT(j, i)
                gammarhoji = gl%RFGAMMARHO(j, i)
                ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
                crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
                d2ftji_dxi2 = (1.D0 - (xj + xi)/(betatji**2*xj + xi)) / (betatji**2*xj + xi) * &
                    & (2.D0*xj - xj*xi*2.D0/(betatji**2*xj + xi))
                !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:
                d2frhoji_dxi2 = (1.D0 - (xj + xi)/(betarhoji**2*xj + xi)) / (betarhoji**2*xj + xi) * &
                    & (2.D0*xj - xj*xi*2.D0/(betarhoji**2*xj + xi))
                ! 2nd derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
                isoch_therm%d2Tr_dxidxj_it(i,i) = isoch_therm%d2Tr_dxidxj_it(i,i) +  ct_ji * d2ftji_dxi2
                ! 2nd derivative of vred w.r.t xi - 1st summation for j = 1 to i-1
                d2vred_dxidxj(i,i) = d2vred_dxidxj(i,i) +  crho_ji * d2frhoji_dxi2
            end do
            !2nd loop for j > i
            do j = i + 1, gl%NCOMP
                xj = gl%MOLFRACTIONS(j)
                rhocj = gl%rhoc(j)
                tcj = gl%tc(j)
                betatij = gl%RFBETAT(i, j)
                betarhoij = gl%RFBETARHO(i, j)
                gammatij = gl%RFGAMMAT(i, j)
                gammarhoij = gl%RFGAMMARHO(i, j)
                ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
                crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
                d2ftij_dxi2 = (1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) / (betatij**2*xi + xj) * &
                    & (2.D0*xj - xi*xj*2.D0*betatij**2/(betatij**2*xi + xj))
                !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:
                d2frhoij_dxi2 = (1.D0 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj)) / (betarhoij**2*xi + xj) * &
                    & (2.D0*xj - xi*xj*2.D0*betarhoij**2/(betarhoij**2*xi + xj))
                ! 2nd derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
                isoch_therm%d2Tr_dxidxj_it(i,i) = isoch_therm%d2Tr_dxidxj_it(i,i) +  ct_ij * d2ftij_dxi2
                ! 2nd derivative of vred w.r.t xi - 2nd summation for j = i+1 to ncomp
                d2vred_dxidxj(i,i) = d2vred_dxidxj(i,i) +  crho_ij * d2frhoij_dxi2
            end do
            ! the routine above actually calculates d^2vred/dxi^2, therefore must be transformed to:
            ! d^2rhored/dxi^2 = 2/rhored*(drhored/dxi)^2 - rhored^2*(d^2rhored/dxi^2)
            isoch_therm%d2rhor_dxidxj_it(i,i) = 2.D0/gl%rhoredmix*(isoch_therm%drhor_dxi_it(i))**2 - gl%rhoredmix**2*d2vred_dxidxj(i,i)
        end do
        !-------------------------------------------------------------------------------------------------------------------------------------------

    else if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .OR. (gl%Mix_type == 3) .OR. (gl%Mix_type == 31) .or. (gl%mix_type == 6) .or. (gl%mix_type == 19)) then

        isoch_therm%dTr_dxi_it = 0.D0
        isoch_therm%drhor_dxi_it = 0.D0
        isoch_therm%d2Tr_dxidxj_it = 0.D0
        isoch_therm%d2rhor_dxidxj_it = 0.D0

    else if (gl%Mix_type == 4) then

        !LKP, NOT YET IMPLEMENTED!!
        !For implementing the LKP here, the derivatives of the reducing functions of the LKP are needed with assuming all x being independent variables.

    end if
    !------------------------------------------------------------------------------------------------



    !Calculate the derivative of the reducing function Tr and rhor with respect to rho_i and rho_j
    !------------------------------------------------------------------------------------------------
    isoch_therm%dTr_drhoi_it = 0.D0
    isoch_therm%drhor_drhoi_it = 0.D0
    isoch_therm%d2Tr_drhoidrhoj_it = 0.D0
    isoch_therm%d2rhor_drhoidrhoj_it = 0.D0
    If ((gl%Mix_type == 1) .or. (gl%Mix_type == 4) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then
        !Calculate some help variables
        !Calculate dxm_drhoi_constrhoj
        do m = 1, gl%ncomp
            do i = 1, gl%ncomp
                dxm_drhoi(m,i) = (kron_delta(m,i) - gl%molfractions(m)) / Dens
                do j = i, gl%ncomp
                    d2xm_drhoi_drhoj(m,i,j) = (2.D0 * gl%molfractions(m) - kron_delta(m,j) - kron_delta(m,i)) / Dens**2
                    d2xm_drhoi_drhoj(m,j,i) = d2xm_drhoi_drhoj(m,i,j)
                end do
            end do
        end do
        !-------------------------------------
        !First derivatives
        do i = 1, gl%ncomp
            do m = 1, gl%ncomp
                isoch_therm%dTr_drhoi_it(i) = isoch_therm%dTr_drhoi_it(i) + isoch_therm%dTr_dxi_it(m) * dxm_drhoi(m,i)
                isoch_therm%drhor_drhoi_it(i) = isoch_therm%drhor_drhoi_it(i) + isoch_therm%drhor_dxi_it(m) * dxm_drhoi(m,i)
            end do
        end do
        !-------------------------------------
        !Second derivatives
        do i = 1, gl%ncomp
            do j = i, gl%ncomp
                do m = 1, gl%ncomp
                    isoch_therm%d2Tr_drhoidrhoj_it(i,j) = isoch_therm%d2Tr_drhoidrhoj_it(i,j) + isoch_therm%dTr_dxi_it(m) * d2xm_drhoi_drhoj(m,i,j)
                    isoch_therm%d2rhor_drhoidrhoj_it(i,j) = isoch_therm%d2rhor_drhoidrhoj_it(i,j) + isoch_therm%drhor_dxi_it(m) * d2xm_drhoi_drhoj(m,i,j)
                    do n = 1, gl%ncomp
                        isoch_therm%d2Tr_drhoidrhoj_it(i,j) = isoch_therm%d2Tr_drhoidrhoj_it(i,j) + isoch_therm%d2Tr_dxidxj_it(n,m) *dxm_drhoi(n,j) * dxm_drhoi(m,i)
                        isoch_therm%d2rhor_drhoidrhoj_it(i,j) = isoch_therm%d2rhor_drhoidrhoj_it(i,j) + isoch_therm%d2rhor_dxidxj_it(n,m) *dxm_drhoi(n,j) * dxm_drhoi(m,i)
                    end do
                end do
                isoch_therm%d2Tr_drhoidrhoj_it(j,i) = isoch_therm%d2Tr_drhoidrhoj_it(i,j)
                isoch_therm%d2rhor_drhoidrhoj_it(j,i) = isoch_therm%d2rhor_drhoidrhoj_it(i,j)
            end do
        end do
    end if
    !------------------------------------------------------------------------------------------------




    !Calculate the x-derivatives of the residual dimensionless Helmholtz energies
    !------------------------------------------------------------------------------------------------
    if ((gl%Mix_type == 1) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then !Multi-fluid mixture model and models based on this model (MIXTURE MODELS 11 and 19 NOT YET IMPLEMENTED)

        !Calculate first x-derivative
        Do i = 1, gl%ncomp
            !First derivative with respect to xi
            isoch_therm%dalphar_dxi_it(i) = FNRDER(1,i)
            !Second derivative with respect to xi and delta MULTIPLIED WITH DELTA!
            isoch_therm%d2alphar_dxiddel_it(i) = FNRDER(2,i)
            !Second derivative with respect to xi and tau MULTIPLIED WITH TAU!
            isoch_therm%d2alphar_dxidtau_it(i) = FNRDER(4,i)

            !Add contribution of the departure function for mix_type == 1
            if (gl%Mix_type == 1) then
                do j = 1, gl%ncomp
                    !Calculate first x-derivative
                    isoch_therm%dalphar_dxi_it(i) = isoch_therm%dalphar_dxi_it(i) + gl%molfractions(j) * gl%Fij(i, j) * DEPFUNCBIN(1,i,j)
                    isoch_therm%d2alphar_dxiddel_it(i) = isoch_therm%d2alphar_dxiddel_it(i) + gl%molfractions(j) * gl%Fij(i, j) * DEPFUNCBIN(2,i,j)
                    isoch_therm%d2alphar_dxidtau_it(i) = isoch_therm%d2alphar_dxidtau_it(i) + gl%molfractions(j) * gl%Fij(i, j) * DEPFUNCBIN(4,i,j)
                    !Calculate second x-derivative
                    isoch_therm%d2alphar_dxidxj_it(i,j) = gl%Fij(i, j) * DEPFUNCBIN(1,i,j)
                end do
            end if

        end do

        if ((gl%Mix_type == 12) .or. (gl%Mix_type == 13)) then

            !THIS SIMPLIFIED VERSION DOES NOT WORK AT THE MOMENT, BECAUSE THERE IS NO WAY (?) OF TRANSFORMING X-DERIVATIVES WITH LAST XN REPLACED TO X-DERIVATIVES WITH ALL X BEING INDEPENDENT
            !!----------------------------------------------------------------
            !! "slow" solution for the moment
            !! Call the functions that calculate the derivatives
            !! Andreas Jäger, November 2018
            !
            !! get the derivatives d(ar)/d(xi) for all xi
            !call dar_dxi (gl, Temp, dens, DERIVFX)
            !
            !! get the derivatives del*d^2(ar)/d(xi)d(del) for all xi
            !call d2ar_dxiddel (gl, Temp, Dens, DERIVFXDEL)
            !
            !! get the derivatives tau*d^2(ar)/d(xi)d(tau) for all xi
            !call d2ar_dxidtau (gl, Temp, dens, DERIVFXTAU)
            !
            !! get the derivatives d^2(ar)/d(xi)d(xj) for all xi, xj
            !call d2ar_dxidxj (gl, Temp, dens, DERIVFXiXj)
            !call d2ar_dxi2 (gl, Temp, dens, DERIVFX2)
            !
            !!Do a variable transformation from dependend mole fractions to independent mole fractions and write the results in the type variables
            !do i = 1, gl%ncomp
            !    !First x-derivatives
            !    isoch_therm%dalphar_dxi_it(i) = DERIVFX(i)
            !    isoch_therm%d2alphar_dxiddel_it(i) = DERIVFXDEL(i)
            !    isoch_therm%d2alphar_dxidtau_it(i) = DERIVFXTAU(i)
            !    do j = 1, gl%ncomp
            !        !Second x-derivatives
            !        if (i == j) then
            !            isoch_therm%d2alphar_dxidxj_it(i,j) = DERIVFX2(i)
            !        else
            !            isoch_therm%d2alphar_dxidxj_it(i,j) = DERIVFXiXj(i,j)
            !        end if
            !    end do
            !end do
            !!----------------------------------------------------------------

        end if


    else    !All other models

        !PCP-SAFT
        if (gl%Mix_type == 6) then

            !Get the first derivatives with respect to xi
            GETDER = 0
            GETDER(1) = 1   !alphar                                 with respect to xi
            GETDER(2) = 1   !dalphar_ddelta         * delta         with respect to xi
            GETDER(4) = 1   !dalphar_dtau           * tau           with respect to xi
            call ARX1DERIVS(gl, Temp, Dens, GETDER, DERIVFX_ALL)

            !Get the second derivatives with respect to xi
            GETDER = 0
            GETDER(1) = 1   !alphar                                 with respect to xi and xj
            call ARX2DERIVS(gl, Temp, Dens, GETDER, DERIVFX2_ALL)

            do i = 1, gl%ncomp
                !First x-derivatives
                isoch_therm%dalphar_dxi_it(i) = DERIVFX_ALL(1,i)
                isoch_therm%d2alphar_dxiddel_it(i) = DERIVFX_ALL(2,i)
                isoch_therm%d2alphar_dxidtau_it(i) = DERIVFX_ALL(4,i)
                do j = 1, gl%ncomp
                    !Second x-derivatives
                    isoch_therm%d2alphar_dxidxj_it(i,j) = DERIVFX2_ALL(1,i,j)
                end do
            end do

        end if


        !THIS SIMPLIFIED VERSION DOES NOT WORK AT THE MOMENT, BECAUSE THERE IS NO WAY (?) OF TRANSFORMING X-DERIVATIVES WITH LAST XN REPLACED TO X-DERIVATIVES WITH ALL X BEING INDEPENDENT
        !!----------------------------------------------------------------
        !! "slow" solution for the moment
        !! Call the functions that calculate the derivatives
        !! Andreas Jäger, November 2018
        !
        !! get the derivatives d(ar)/d(xi) for all xi
        !call dar_dxi (gl, Temp, dens, DERIVFX)
        !
        !! get the derivatives del*d^2(ar)/d(xi)d(del) for all xi
        !call d2ar_dxiddel (gl, Temp, Dens, DERIVFXDEL)
        !
        !! get the derivatives tau*d^2(ar)/d(xi)d(tau) for all xi
        !call d2ar_dxidtau (gl, Temp, dens, DERIVFXTAU)
        !
        !! get the derivatives d^2(ar)/d(xi)d(xj) for all xi, xj
        !call d2ar_dxidxj (gl, Temp, dens, DERIVFXiXj)
        !call d2ar_dxi2 (gl, Temp, dens, DERIVFX2)
        !
        !!Do a variable transformation from dependend mole fractions to independent mole fractions and write the results in the type variables
        !do i = 1, gl%ncomp
        !    !First x-derivatives
        !    isoch_therm%dalphar_dxi_it(i) = DERIVFX(i)
        !    isoch_therm%d2alphar_dxiddel_it(i) = DERIVFXDEL(i)
        !    isoch_therm%d2alphar_dxidtau_it(i) = DERIVFXTAU(i)
        !    do j = 1, gl%ncomp
        !        !Second x-derivatives
        !        if (i == j) then
        !            isoch_therm%d2alphar_dxidxj_it(i,j) = DERIVFX2(i)
        !        else
        !            isoch_therm%d2alphar_dxidxj_it(i,j) = DERIVFXiXj(i,j)
        !        end if
        !    end do
        !end do
        !!----------------------------------------------------------------


    end if
    !------------------------------------------------------------------------------------------------


    end subroutine Derivs_isoch_therm





    !Subroutine for tracing isotherms with isochoric thermodynamics
    !Andreas Jäger, October 2018
    subroutine px_diag_isochoric(gl, Temp, spec_deriv, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, errval)
    !--------------------------------------------------------------------------------
    ! This subroutine calculates an isotherm in a px-diagram with isochoric thermodynamics
    ! The basic algorithm and method is described in:
    ! Bell, I.H.; Deiters U.K.: "On the Construction of Binary Mixture p-x and T-x Diagrams from Isochoric Thermodynamics", AIChe Journal 64(7), 2745-2757, 2018.
    !
    !
    ! Input:    Temperature     [K]
    !           spec_deriv       -   Indicate which molar concentration is used for tracing, in case both fluids are supercritical concerning the temperature  (1 = drhoi / drho_1_vap, 2 = drhoi / drho_2_vap, 3 = drhoi / drho_1_liq, 4 = drhoi / drho_2_liq)
    !
    ! Output:   p_points        -   Array (10001) containing the pressure of the points
    !           T_points        -   Array (10001) containing the temperature of the points
    !           x_points        -   Matrix (10001,2) containing the molfraction of component 1 in the vapor (:,1) and the molfraction of component 1 in the liquid (:,2)
    !           rho_vap_points  -   Array (10001) containing the densities of the vapor phase at the calculated points
    !           rho_liq_points  -   Array (10001) containing the densities of the liquid phase at the calculated points
    !           points          -   Number of points calculated
    !           fileout         -   If a path to a file is supplied here, the results are written to this file
    !           errval: specifies if an error occured during the calculations
    !--------------------------------------------------------------------------------
    ! A. Jäger, October 2018





    implicit none

    !Input and output variables of the routine
    type(type_gl) :: gl
    double precision:: Temp
    integer :: spec_deriv                                   !specifier which derivative is needed (1 = drhoi / drho_1_vap, 2 = drhoi / drho_2_vap, 3 = drhoi / drho_1_liq, 4 = drhoi / drho_2_liq)
    integer, parameter:: maxi = 300 ! maximum number of calculated points, length of the return vectors!
    double precision, dimension(maxi) :: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision :: x_points(maxi,2)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_liq(2)
    character(255) :: fileout    ! optional path for output file
    integer:: errval
    integer :: points

    !!Create type in which all needed properties of one phase for the isochoric tracing are stored
    !type(type_isoch_therm), dimension(2) :: props_phase
    !
    !!Derivatives of the residual Helmholtz energy density (Psi^r) with respect to tau, delta, x and help variables
    !double precision :: dPsir_ddelta
    !double precision :: dPsir_dtau
    !double precision, dimension(:), allocatable :: dPsir_dxm
    !double precision :: d2Psir_ddelta2
    !double precision :: d2Psir_dtauddelta
    !double precision :: d2Psir_dtau2
    !double precision, dimension(:), allocatable :: d2Psir_dxmddel
    !double precision, dimension(:), allocatable :: d2Psir_dxmdtau
    !double precision, dimension(:,:), allocatable :: d2Psir_dxmdxn
    !!Help variables
    !double precision, dimension(:,:), allocatable :: kron_delta
    !double precision, dimension(:), allocatable :: drhorTr_dxm
    !double precision, dimension(:,:), allocatable :: d2rhorTr_dxmdxn
    !double precision, dimension(:), allocatable :: ddelta_drhoi
    !double precision, dimension(:), allocatable :: dtau_drhoi
    !double precision, dimension(:,:), allocatable :: dxm_drhoi
    !double precision, dimension(:), allocatable :: d2tau_drhoi_dT
    !double precision :: dtau_dT
    !double precision, dimension(:,:), allocatable :: d2tau_drhoidrhoj
    !double precision, dimension(:,:), allocatable :: d2delta_drhoidrhoj
    !double precision, dimension(:,:,:), allocatable :: d2xm_drhoidrhoj
    !
    !!Derivatives of the residual Helmholtz energy density (Psi^r) with respect to independent variables
    !!double precision, dimension(:), allocatable :: dPsi0_drhoi              !First derivative of the ideal part of Psi with respect to rho_i
    !double precision, dimension(:), allocatable :: dPsir_drhoi              !First derivative of the residual part of Psi with respect to rho_i
    !!double precision, dimension(:), allocatable :: d2Psi0_drhoi_dT          !First derivative of the ideal part of Psi with respect to rho_i
    !double precision, dimension(:), allocatable :: d2Psir_drhoidT          !First derivative of the residual part of Psi with respect to rho_i
    !!double precision, dimension(:,:), allocatable :: d2Psi0_drhoi_drhoj     !First derivative of the ideal part of Psi with respect to rho_i
    !double precision, dimension(:,:), allocatable :: d2Psir_drhoidrhoj     !First derivative of the residual part of Psi with respect to rho_i
    !
    !!Mixed derivatives of the residual Helmholtz energy density with respect to rhoi, delta, tau, and xi
    !double precision, dimension(:), allocatable :: ddPsirddelta_drhoj
    !double precision, dimension(:), allocatable :: ddPsirdtau_drhoj
    !double precision, dimension(:,:), allocatable :: ddPsirdxm_drhoj

    !Other variables
    double precision :: Dens
    double precision :: R_mix
    integer :: a, i, j, k, m, n
    double precision, dimension(2):: psat_fluid   !Saturation pressures of fluid 1 and 2
    logical, dimension(2) :: supercritical
    double precision, dimension(30) :: z_overall
    integer:: start_fluid, second_fluid   !Indicates which fluid has the lower vapor pressure. Start with this fluid and save other fluid in second_fluid

    !!Variables for tracing
    !double precision, dimension(:,:,:), allocatable:: Hessian       !Last entry: Number of the phase
    double precision, dimension(2) :: dens_phase                    !densities of the phases: dens_phase(1): vapor phase, dens_phase(2): liquid phase
    double precision, dimension(:,:), allocatable :: dens_i_phase   !Molar concentratons of the phases: dens_i_phase(1:ncomp,1): vapor phase, dens_i_phase(1:ncomp,2): liquid phase
    double precision, dimension(:,:), allocatable :: x_phase        !compositions of the phases: x_phase(1:ncomp,1): vapor phase,x_phase(1:ncomp,2): liquid phase
    !double precision, dimension(:,:), allocatable :: drhoi_dp       !Derivative of molar concentrations of the different phases along the phase boundary. drhoi_dp(1:ncomp,1): vapor phase,drhoi_dp(1:ncomp,2): liquid phase

    !Variables for pure vapor pressure iteration
    double precision, dimension(2) :: d_vap_fluid, d_liq_fluid
    integer :: iFlash, nrsubst
    integer :: iter

    !Variables for mixture vapor pressure iteration
    double precision:: press, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: Nr_x_given, iPhase_try

    !!Variables for Gauss-algorithm. ADJUST THE SIZE DYNAMICALLY!!!!
    !double precision, dimension(60,60):: MatrixA
    !integer:: rankA
    !double precision, dimension(60):: vectorb
    !double precision, dimension(60):: vectorx
    !double precision:: det_A

    !Variables for calculating drhoi_drhoj
    integer:: nr_of_eqs
    !double precision:: rhoi_march
    double precision, dimension(4) :: y_vector              !Four equations need to be solved in the ODE
    double precision, dimension(10) :: ind_var              !Other independent variables that might be needed to evaluate the derivatives of y with respect to t (here: ind_var(1) = temperature)
    !double precision, dimension(4) :: dy_dt

    !Variables needed for the Runge-Kutta-Fehlberg method
    double precision :: t_march
    double precision :: t_min
    double precision :: t_max
    double precision :: step_min
    double precision :: eps_max
    double precision :: fac_rel
    double precision, dimension (1000) :: t_res
    double precision, dimension (4, 1000) :: y_res
    double precision, dimension (10, 1000):: res_var
    integer :: nr_of_points
    logical :: post_proc
    logical :: pressure_march                               !If pressure_march == .true., march in pressure (used if one component is supercritical)
    logical :: molar_conc_march                             !If molar_conc_march == .true., march in the specified (defined in spec_deriv) molar concentration
    logical :: dew_failed


    !!Test variables
    !!double precision, dimension(30) :: ndTred_dni_1, ndrhored_dni_1
    !double precision, dimension(30) :: Chempot_r
    !integer:: OIR



    !do a = 1, 2
    !    !Allocate generated type props_phase and define size of variables for both phases
    !    if(.not.allocated(props_phase(a)%alpha_it))       allocate(props_phase(a)%alpha_it(nderivs))
    !    if(.not.allocated(props_phase(a)%alpha0_it))      allocate(props_phase(a)%alpha0_it(nderivs))
    !    if(.not.allocated(props_phase(a)%alphar_it))      allocate(props_phase(a)%alphar_it(nderivs))
    !
    !    if(.not.allocated(props_phase(a)%dTr_drhoi_it))           allocate(props_phase(a)%dTr_drhoi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%dTr_dxi_it))             allocate(props_phase(a)%dTr_dxi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2Tr_drhoidrhoj_it))     allocate(props_phase(a)%d2Tr_drhoidrhoj_it(gl%ncomp,gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2Tr_dxidxj_it))         allocate(props_phase(a)%d2Tr_dxidxj_it(gl%ncomp,gl%ncomp))
    !
    !    if(.not.allocated(props_phase(a)%drhor_drhoi_it))         allocate(props_phase(a)%drhor_drhoi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%drhor_dxi_it))           allocate(props_phase(a)%drhor_dxi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2rhor_drhoidrhoj_it))   allocate(props_phase(a)%d2rhor_drhoidrhoj_it(gl%ncomp,gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2rhor_dxidxj_it))       allocate(props_phase(a)%d2rhor_dxidxj_it(gl%ncomp,gl%ncomp))
    !
    !    if(.not.allocated(props_phase(a)%dalpha_dxi_it))          allocate(props_phase(a)%dalpha_dxi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alpha_dxiddel_it))     allocate(props_phase(a)%d2alpha_dxiddel_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alpha_dxidtau_it))     allocate(props_phase(a)%d2alpha_dxidtau_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alpha_dxidxj_it))      allocate(props_phase(a)%d2alpha_dxidxj_it(gl%ncomp,gl%ncomp))
    !
    !    if(.not.allocated(props_phase(a)%dalphar_dxi_it))          allocate(props_phase(a)%dalphar_dxi_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alphar_dxiddel_it))     allocate(props_phase(a)%d2alphar_dxiddel_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alphar_dxidtau_it))     allocate(props_phase(a)%d2alphar_dxidtau_it(gl%ncomp))
    !    if(.not.allocated(props_phase(a)%d2alphar_dxidxj_it))      allocate(props_phase(a)%d2alphar_dxidxj_it(gl%ncomp,gl%ncomp))
    !end do
    !
    !allocate(dPsir_dxm(gl%ncomp))
    !allocate(d2Psir_dxmddel(gl%ncomp))
    !allocate(d2Psir_dxmdtau(gl%ncomp))
    !allocate(d2Psir_dxmdxn(gl%ncomp,gl%ncomp))
    !allocate(kron_delta(gl%ncomp,gl%ncomp))
    !allocate(drhorTr_dxm(gl%ncomp))
    !allocate(d2rhorTr_dxmdxn(gl%ncomp,gl%ncomp))
    !allocate(ddelta_drhoi(gl%ncomp))
    !allocate(dtau_drhoi(gl%ncomp))
    !allocate(dxm_drhoi(gl%ncomp,gl%ncomp))
    !allocate(d2tau_drhoi_dT(gl%ncomp))
    !allocate(d2tau_drhoidrhoj(gl%ncomp,gl%ncomp))
    !allocate(d2delta_drhoidrhoj(gl%ncomp,gl%ncomp))
    !allocate(d2xm_drhoidrhoj(gl%ncomp,gl%ncomp,gl%ncomp))
    !
    !!allocate(dPsi0_drhoi(gl%ncomp))
    !allocate(dPsir_drhoi(gl%ncomp))
    !!allocate(d2Psi0_drhoidrhoj(gl%ncomp,gl%ncomp))
    !allocate(d2Psir_drhoidrhoj(gl%ncomp,gl%ncomp))
    !allocate(d2Psir_drhoidT(gl%ncomp))
    !
    !allocate(ddPsirddelta_drhoj(gl%ncomp))
    !allocate(ddPsirdtau_drhoj(gl%ncomp))
    !allocate(ddPsirdxm_drhoj(gl%ncomp,gl%ncomp))
    !
    !allocate(Hessian(gl%ncomp,gl%ncomp,2))
    allocate(dens_i_phase(gl%ncomp,2))
    allocate(x_phase(gl%ncomp,2))
    !allocate(drhoi_dp(gl%ncomp,2))

1000 format(f12.7, f10.4, f10.6, f10.6, f10.6, f10.6, f10.3, f10.3)

    !Initialize variables
    pressure_march = .false.
    molar_conc_march = .false.
    dew_failed = .false.
    errval = 0


    !ndTred_dni = 0.D0
    !ndrhored_dni = 0.D0

    !!Fill the Kronecker delta
    !do i = 1, gl%ncomp
    !    do j = 1, gl%ncomp
    !        if (j == i) then
    !            kron_delta(i,j) = 1.D0
    !        else
    !            kron_delta(i,j) = 0.D0
    !        end if
    !    end do
    !end do


    !!Get the mixture gas constant
    !call R_mix_calc(gl, R_mix)


    !If one component is below its triple point temperature, quit with error
    if (Temp < gl%ttp(1)) then
        errval = -9981
        return
    end if
    if (Temp < gl%ttp(2)) then
        errval = -9981
        return
    end if

    !Initialize variables
    supercritical = .false.
    !--------------------------------------------------------------------------------
    ! STEP 1: PURE FLUID START VALUES
    !--------------------------------------------------------------------------------
    ! get the vapor pressures of the pure fluids
    psat_fluid = 0.d0

    !Check if fluid 1 is supercritical and if not, calculate its vapor pressure at given temperature
    if (temp < gl%tc(1)) then
        psat_fluid(1) = vp_eq(gl,Temp, 1)
    else
        supercritical(1) = .true.
    end if

    !Check if fluid 2 is supercritical and if not, calculate its vapor pressure at given temperature
    if (temp < gl%tc(2)) then
        psat_fluid(2) = vp_eq(gl,Temp, 2)
    else
        supercritical(2) = .true.
    end if

    if (supercritical(1) .and. supercritical(2)) then   !Both fluids supercritical, either no px-diagram calculation possible or if possible not implemented yet
        errval = -9981
        return
    elseif(supercritical(1) == .true.) then
        !If the first component is supercritical, start at pure second fluid
        start_fluid = 2
        second_fluid = 1
        !Set marching variable:
        pressure_march = .true.
        molar_conc_march = .false.
    elseif(supercritical(2) == .true.) then
        !If the second component is supercritical, start at pure first fluid
        start_fluid = 1
        second_fluid = 2
        pressure_march = .true.
        molar_conc_march = .false.
        !TODO
    else
        !if (psat_fluid(1) < psat_fluid(2)) then     !Start marching at the pure fluid with the lower vapor pressure (arbitrarily chosen)
        !    start_fluid = 1
        !    second_fluid = 2
        !else
        !    start_fluid = 2
        !    second_fluid = 1
        !end if
        pressure_march = .false.
        molar_conc_march = .true.
        !Check if spec_deriv is between 1 and 4, if not, quit with error
        if ((spec_deriv < 1) .or. (spec_deriv > 4)) then
            errval = -9955
            return
        end if
        !spec_deriv = 3      !Specify which variable is the marching variable (1 = drho_1_vap, 2 = drho_2_vap, 3 = drho_1_liq, 4 = drho_2_liq) --> This needs to be done by the user now
        if ((spec_deriv == 2) .or. (spec_deriv == 4)) then
            !Component 1 is start fluid (This means that marching is done in the molar density of the second fluid starting with rho_i_secondfluid close to 0!!)
            start_fluid = 1
            second_fluid = 2
        else
            !Component 2 is start fluid (This means that marching is done in the molar density of the first fluid starting with rho_i_firstfluid close to 0!!)
            start_fluid = 2
            second_fluid = 1
        end if
    end if

    if (supercritical(1) == .false.) then
        !Get the vapor pressure of fluid 1
        !----------------------------------
        d_vap_fluid(1) = 0.D0
        d_liq_fluid(1) = 0.D0
        iFlash = 1
        nrsubst = 1
        call Flash_Pure_PhaseBoundary(gl,psat_fluid(1), Temp, d_vap_fluid(1), d_liq_fluid(1), iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 1) then
            !points = points + 1
            !T_points(points) = Temp
            !p_points(points) = psat_fluid(1)
            !x_points(points, 1) = 1.d0 !xvap1
            !x_points(points, 2) = 1.d0 !xliq1
            !rhovap_points(points) = d_vap
            !rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if

    if (supercritical(2) == .false.) then
        !Get the vapor pressure of fluid 2
        !----------------------------------
        d_vap_fluid(2) = 0.D0
        d_liq_fluid(2) = 0.D0
        iFlash = 1
        nrsubst = 2
        call Flash_Pure_PhaseBoundary(gl,psat_fluid(2), Temp, d_vap_fluid(2), d_liq_fluid(2), iFlash, errval, iter, nrsubst)
        if (errval /= 0) then
            return
        end if
        if (start_fluid == 2) then
            !points = points + 1
            !T_points(points) = Temp
            !p_points(points) = psat_fluid(2)
            !x_points(points, 1) = 1.d0 !xvap1
            !x_points(points, 2) = 1.d0 !xliq1
            !rhovap_points(points) = d_vap
            !rholiq_points(points) = d_liq
        end if
        !----------------------------------
    end if
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! STEP 2: FIRST DEW POINT OF THE MIXTURE
    !--------------------------------------------------------------------------------

    !Set overall composition
    z_overall(start_fluid) = 0.99999D0
    z_overall(second_fluid) = 0.00001D0

    do i = 1, 20    !Try up to 20 times calculating the first dew point of the mixture
        !Calculate the first mixture dew point at given composition
        !-------------------------------------------------------------
        gl%molfractions = z_overall
        call reduced_parameters_calc(gl, Temp)
        press = 0.D0
        x_vap = 0.D0
        x_liq = 0.D0
        rhovap_est = 0.D0
        rholiq_est = 0.D0
        vapfrac = 0
        iFlash = 4      !Dew point calculation at given temperature
        iPhase_try = 5  !Try VLE
        Nr_x_given = 0
        call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
            & iPhase_try, Nr_x_given, errval, iter)

        if (errval /= 0) then
            if (i == 20) then
                !If the calculation of a dew point failed, try to calculate a bubble point instead
                dew_failed = .true.
            end if
            z_overall(start_fluid) = z_overall(start_fluid) - 0.001D0
            z_overall(second_fluid) = 1.D0 - z_overall(start_fluid)
        else
            p_points(1) = press
            T_points(1) = Temp
            x_points(i, 1) = x_vap(1)
            x_points(i, 2) = x_liq(1)
            exit
        end if
    end do
    !-------------------------------------------------------------

    if (dew_failed == .true.) then

        !Reset overall composition
        z_overall(start_fluid) = 0.99D0
        z_overall(second_fluid) = 0.01D0

        do i = 1, 20    !Try up to 20 times calculating the first bubble point of the mixture
            !Calculate the first mixture bubble point at given composition
            !-------------------------------------------------------------
            gl%molfractions = z_overall
            call reduced_parameters_calc(gl, Temp)
            press = 0.D0
            x_vap = 0.D0
            x_liq = 0.D0
            rhovap_est = 0.D0
            rholiq_est = 0.D0
            vapfrac = 0
            iFlash = 3      !Bubble point calculation at given temperature
            iPhase_try = 5  !Try VLE
            Nr_x_given = 0
            call Flash_PhaseBoundary_calc(gl,press, Temp, z_overall, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                & iPhase_try, Nr_x_given, errval, iter)

            if (errval /= 0) then
                if (i == 20) then
                    !if the dew and bubble point calculation failed, quit with error
                    return
                end if
                z_overall(start_fluid) = z_overall(start_fluid) - 0.001D0
                z_overall(second_fluid) = 1.D0 - z_overall(start_fluid)
            else
                p_points(1) = press
                T_points(1) = Temp
                x_points(i, 1) = x_vap(1)
                x_points(i, 2) = x_liq(1)
                exit
            end if
        end do

    end if
    !-------------------------------------------------------------

    !--------------------------------------------------------------------------------



    !Save densities of the phase
    dens_phase(1) = gl%rho_vap
    dens_phase(2) = gl%rho_liq

    !Save compositions of the phases and calculate molar concentrations
    !Vapor phase
    x_phase(start_fluid, 1) = x_vap(start_fluid)                        !Phase composition
    x_phase(second_fluid, 1) = x_vap(second_fluid)                      !Phase composition
    dens_i_phase(start_fluid, 1) = x_vap(start_fluid) * dens_phase(1)   !Molar concentration
    dens_i_phase(second_fluid, 1) = x_vap(second_fluid) * dens_phase(1) !Molar concentration
    !Liquid phase
    x_phase(start_fluid, 2) = x_liq(start_fluid)                        !Phase composition
    x_phase(second_fluid, 2) = x_liq(second_fluid)                      !Phase composition
    dens_i_phase(start_fluid, 2) = x_liq(start_fluid) * dens_phase(2)   !Molar concentration
    dens_i_phase(second_fluid, 2) = x_liq(second_fluid) * dens_phase(2) !Molar concentration


    !--------------------------------------------------------------------------------
    ! STEP 3: Call the Runge-Kutta-Fehlberg method to march along the phase boundary
    !--------------------------------------------------------------------------------
    If (molar_conc_march .eqv. .true.) then
        nr_of_eqs = 4       !Four equations to be solved in the ODE (see Eq. 40 in the article of Bell and Deiters)
        !The marching will be carried out with the molar concentration being the independent variable from the pure start_fluid to the second_fluid
        !This means, that marching starts at a molar concentration (of vapor or liquid) of the second fluid close to 0 and ends close to the density (of vapor or liquid) of the second fluid
        if (spec_deriv < 3) then
            !Marching in the vapor phase molar concentration of second_fluid
            t_min = dens_i_phase(second_fluid, 1)
            t_max = d_vap_fluid(second_fluid) - 1.D-6
        else
            !Marching in the liquid phase molar concentration of second_fluid
            t_min = dens_i_phase(second_fluid, 2)
            t_max = d_liq_fluid(second_fluid) - 1.D-6
        end if
        y_vector(1) = dens_i_phase(1, 1)
        y_vector(2) = dens_i_phase(2, 1)
        y_vector(3) = dens_i_phase(1, 2)
        y_vector(4) = dens_i_phase(2, 2)
        step_min = 1.D-4                !Minimum step allowed (in mol/m³)
        eps_max = 1.D-8
        fac_rel = 1.D0!0.95D0
        errval = 0
        ind_var = 0.D0
        ind_var(1) = Temp               !An additional independent variable of drhoi_drhoj is the temperature
        ind_var(2) = spec_deriv         !Specify the marching variable
        post_proc = .true.              !Use post-processing routine. Here: calculate improved values for rho_i by solving the phase equilibrium condition with Newton-Raphson method.
        call RKF45_ODE_INT(gl, drhoi_drhoj, ind_var, post_proc, nr_of_eqs, y_vector, t_min, t_max, step_min, eps_max, fac_rel, nr_of_points, t_res, y_res, res_var, errval)
        if ((errval /= 0) .and. (nr_of_points < 2)) return
    end if

    if (pressure_march .eqv. .true.) then
        nr_of_eqs = 4       !Four equations to be solved in the ODE (see Eq. 39 in the article of Bell and Deiters)
        !The marching will be carried out in pressure from the pure start_fluid to an arbitrarily defined high pressure
        t_min = p_points(1) * 1.D6
        t_max = 1000 * 1.D6 !Upper pressure boundary: 1 GPA (Arbitrarily picked)
        y_vector(1) = dens_i_phase(1, 1)
        y_vector(2) = dens_i_phase(2, 1)
        y_vector(3) = dens_i_phase(1, 2)
        y_vector(4) = dens_i_phase(2, 2)
        step_min = 1                   !Minimum step allowed (in Pa)
        eps_max = 1.D-7!1.D-8
        fac_rel = 1.D0!0.95D0
        errval = 0
        ind_var = 0.D0
        ind_var(1) = Temp               !An additional independent variable of drhoi_drhoj is the temperature
        ind_var(2) = 5                  !Specify the marching variable. 5 = March in pressure
        post_proc = .true.              !Use post-processing routine. Here: calculate improved values for rho_i by solving the phase equilibrium condition with Newton-Raphson method.
        call RKF45_ODE_INT(gl, drhoi_drhoj, ind_var, post_proc, nr_of_eqs, y_vector, t_min, t_max, step_min, eps_max, fac_rel, nr_of_points, t_res, y_res, res_var, errval)
        if ((errval /= 0) .and. (nr_of_points < 2)) return
    end if
    !--------------------------------------------------------------------------------



    !--------------------------------------------------------------------------------
    ! STEP 5: Write the results in the output variables of the subroutine
    !--------------------------------------------------------------------------------
    !Write the points from fluid 1 to fluid 2 (p plotted over x_component2)
    if (molar_conc_march .eqv. .true.) then
        points = nr_of_points + 1
        p_points(1) = psat_fluid(1)
        T_points(1) = Temp
        rhovap_points(1) = d_vap_fluid(1)
        rholiq_points(1) = d_liq_fluid(1)
        x_points(1, 1) = 1.D0
        x_points(1, 2) = 1.D0
        if (start_fluid == 1) then
            !If tracing starts at fluid 1, go from 2 to nr_of_points
            do i = 2, nr_of_points
                p_points(i) = res_var(1,i) / 1.D6
                T_points(i) = Temp
                rhovap_points(i) = y_res(1,i) + y_res(2,i)
                rholiq_points(i) = y_res(3,i) + y_res(4,i)
                x_points(i, 1) = y_res(1,i) / (y_res(1,i) + y_res(2,i))
                x_points(i, 2) = y_res(3,i) / (y_res(3,i) + y_res(4,i))
            end do
        else
            !If tracing starts at fluid 2, go from nr_of_points to 2
            do i = 2, nr_of_points
                p_points(nr_of_points+2-i) = res_var(1,i) / 1.D6
                T_points(nr_of_points+2-i) = Temp
                rhovap_points(nr_of_points+2-i) = y_res(1,i) + y_res(2,i)
                rholiq_points(nr_of_points+2-i) = y_res(3,i) + y_res(4,i)
                x_points(nr_of_points+2-i, 1) = y_res(1,i) / (y_res(1,i) + y_res(2,i))
                x_points(nr_of_points+2-i, 2) = y_res(3,i) / (y_res(3,i) + y_res(4,i))
            end do
        end if
        p_points(nr_of_points + 1) = psat_fluid(2)
        T_points(nr_of_points + 1) = Temp
        rhovap_points(nr_of_points + 1) = d_vap_fluid(2)
        rholiq_points(nr_of_points + 1) = d_liq_fluid(2)
        x_points(nr_of_points + 1, 1) = 0.D0
        x_points(nr_of_points + 1, 2) = 0.D0
    end if
    if (pressure_march .eqv. .true.) then
        points = nr_of_points
        p_points(1) = psat_fluid(start_fluid)
        T_points(1) = Temp
        rhovap_points(1) = d_vap_fluid(start_fluid)
        rholiq_points(1) = d_liq_fluid(start_fluid)
        if (start_fluid == 2) then
            x_points(1, 1) = 0.D0
            x_points(1, 2) = 0.D0
        else
            x_points(1, 1) = 1.D0
            x_points(1, 2) = 1.D0
        end if
        !if (start_fluid == 2) then
        !If tracing starts at fluid 2, go from 2 to nr_of_points
        do i = 2, nr_of_points
            p_points(i) = res_var(1,i) / 1.D6
            T_points(i) = Temp
            rhovap_points(i) = y_res(1,i) + y_res(2,i)
            rholiq_points(i) = y_res(3,i) + y_res(4,i)
            x_points(i, 1) = y_res(1,i) / (y_res(1,i) + y_res(2,i))
            x_points(i, 2) = y_res(3,i) / (y_res(3,i) + y_res(4,i))
        end do
        !else
        !    !If tracing starts at fluid 1, go from nr_of_points to 2
        !    do i = 2, nr_of_points
        !        p_points(nr_of_points+2-i) = res_var(1,i) / 1.D6
        !        T_points(nr_of_points+2-i) = Temp
        !        rhovap_points(nr_of_points+2-i) = y_res(1,i) + y_res(2,i)
        !        rholiq_points(nr_of_points+2-i) = y_res(3,i) + y_res(4,i)
        !        x_points(nr_of_points+2-i, 1) = y_res(1,i) / (y_res(1,i) + y_res(2,i))
        !        x_points(nr_of_points+2-i, 2) = y_res(3,i) / (y_res(3,i) + y_res(4,i))
        !    end do
        !end if
    end if
    !--------------------------------------------------------------------------------


#IF DEFINED(_DEBUG) 
    !--------------------------------------------------------------------------------
    ! CHECK IF WRITING THE DATA TO A FILE IS WISHED
    !--------------------------------------------------------------------------------
    !111 continue
    if (fileout /= '') then
        open(unit = 13,file=fileout, status='unknown', action='write', iostat=errval)
        if (errval /= 0) return
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'! P-X DIAGRAM DATA FILE FOR THE BINARY SYSTEM'
        write(13,*)'! ', trim(gl%components(1)),'-',trim(gl%components(2))
        write(13,*)'!'
        write(13,'(A20, F8.3)')'! TEMPERATURE [K]: ', Temp
        write(13,*)'!---------------------------------------------------------------------------------'
        write(13,*)'!'
        write(13,*)'! p/MPa    ','T/K       ','x',trim(gl%components(1)),'    y',trim(gl%components(1)),'   dl/mol/l ', 'dv/mol/l'
        do i = 1, points
            write(13,1000) p_points(i), T_points(i),  x_points(i, 1), x_points(i, 2), rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
        end do
        !if (errMSG /= '') write(13,*)'! ', errMSG
        close(13)
    end if
    !--------------------------------------------------------------------------------
#ENDIF

    !!Compute the needed derivatives and the Hessians of the phases in equilibrium
    !do a = 1, 2
    !
    !    !Set density and composition of the phase and get all needed derivatives
    !    Dens = dens_phase(a)
    !    gl%molfractions = 0.D0
    !    gl%molfractions(1:2) = x_phase(:,a)
    !    call reduced_parameters_calc(gl,Temp)
    !
    !    call Derivs_isoch_therm(gl, props_phase(a), Temp, Dens, errval)
    !
    !
    !    !Calculate the first derivative of Psi with respect to the molar concentration (i.e. the chemical potentials)
    !    !-------------------------------------------------------------------------------------------------------------
    !    !Calculate help variables
    !    ! ddelta/drhoi, dtau/drhoi, dxm/drhoi, d(rhor*Tr)/dxm
    !    do i = 1, gl%ncomp
    !        ddelta_drhoi(i) = (gl%rhoredmix - Dens * props_phase(a)%drhor_drhoi_it(i)) / gl%rhoredmix**2
    !        dtau_drhoi(i) = props_phase(a)%dTr_drhoi_it(i) / Temp
    !        drhorTr_dxm(i) = gl%tredmix * props_phase(a)%drhor_dxi_it(i) + gl%rhoredmix * props_phase(a)%dTr_dxi_it(i)
    !        do m = 1, gl%ncomp
    !            dxm_drhoi(m,i) = (Dens * kron_delta(m,i) - Dens * gl%molfractions(m)) / Dens**2
    !        end do
    !    end do
    !
    !    !Calculate the derivative of the residual Helmholtz energy density with respect to delta
    !    !dPsir/ddelta
    !    dPsir_ddelta = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1))
    !
    !    !Calculate the derivative of the residual Helmholtz energy density with respect to tau
    !    !dPsir/dtau
    !    dPsir_dtau = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1))
    !
    !    !Calculate the derivative of the residual Helmholtz energy density with respect to xm
    !    !dPsir/dxm
    !    do m = 1, gl%ncomp
    !        dPsir_dxm(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * (props_phase(a)%alphar_it(1) * drhorTr_dxm(m) + gl%rhoredmix * gl%tredmix * props_phase(a)%dalphar_dxi_it(m))
    !    end do
    !
    !    !Calculate the first derivatives of Psi_r with respect to the molar concentrations (rho_i)
    !    !It is (dPsi_r / drho_i = mue_r_i)
    !    Do i = 1, gl%ncomp
    !        dPsir_drhoi(i) = dPsir_ddelta * ddelta_drhoi(i) + dPsir_dtau * dtau_drhoi(i)
    !        do m = 1, gl%ncomp
    !            dPsir_drhoi(i) = dPsir_drhoi(i) + dPsir_dxm(m) * dxm_drhoi(m,i)
    !        end do
    !    End do
    !    !Test
    !    OIR = 2
    !    call dna_dni (gl, Temp, Dens, Chempot_r, OIR)
    !    Chempot_r = Chempot_r * R_mix * Temp
    !    !-------------------------------------------------------------------------------------------------------------
    !
    !
    !
    !
    !    !Calculate the second derivatives of Psi with respect to the molar concentration rho_i and rho_j
    !    !-------------------------------------------------------------------------------------------------------------
    !    !Calculate help variables
    !    ! d2delta/drhoidrhoj, d2tau/drhoidrhoj, d2xm/drhoidrhoj, d(rhor*Tr)/dxmdxn
    !    do i = 1, gl%ncomp
    !        do j = i, gl%ncomp
    !            d2delta_drhoidrhoj(i,j) = (gl%rhoredmix**2 * (-props_phase(a)%drhor_drhoi_it(j) - props_phase(a)%drhor_drhoi_it(i) - Dens * props_phase(a)%d2rhor_drhoidrhoj_it(i,j)) + Dens * props_phase(a)%drhor_drhoi_it(i) * 2.D0 * gl%rhoredmix * props_phase(a)%drhor_drhoi_it(j)) / gl%rhoredmix**4
    !            d2tau_drhoidrhoj(i,j) = props_phase(a)%d2Tr_drhoidrhoj_it(i,j) / Temp
    !            d2rhorTr_dxmdxn(i,j) = gl%tredmix * props_phase(a)%d2rhor_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(i) * props_phase(a)%drhor_dxi_it(j) + gl%rhoredmix * props_phase(a)%d2Tr_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(j) * props_phase(a)%drhor_dxi_it(i)
    !            do m = 1, gl%ncomp
    !                d2xm_drhoidrhoj(m,i,j) = (Dens**2 * (kron_delta(i,m) - kron_delta(j,m)) - (Dens * kron_delta(i,m) - gl%molfractions(m) * Dens) * 2.D0 * Dens) / Dens**4
    !                if (i /= j) then
    !                    d2xm_drhoidrhoj(m,j,i) = d2xm_drhoidrhoj(m,i,j)
    !                end if
    !            end do
    !            if (i /= j) then
    !                d2delta_drhoidrhoj(j,i) = d2delta_drhoidrhoj(i,j)
    !                d2tau_drhoidrhoj(j,i) = d2tau_drhoidrhoj(i,j)
    !                d2rhorTr_dxmdxn(j,i) = d2rhorTr_dxmdxn(i,j)
    !            end if
    !        end do
    !    end do
    !
    !    !Calculate the second derivative of the residual Helmholtz energy density with respect to delta and delta
    !    !d2Psir/ddelta^2
    !    d2Psir_ddelta2 = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(3) + 2.D0 * props_phase(a)%alphar_it(2)) / props_phase(a)%delta_it
    !
    !    !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and delta
    !    !d2Psir/dtauddelta
    !    d2Psir_dtauddelta = gl%rhoredmix * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1) - props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(6))
    !
    !    !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and tau
    !    !d2Psir/dtau^2
    !    d2Psir_dtau2 = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**3 * (props_phase(a)%alphar_it(5) - 2.D0 * props_phase(a)%alphar_it(4) + 2.D0 * props_phase(a)%alphar_it(1))
    !
    !    do m = 1, gl%ncomp
    !        !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and delta
    !        !d2Psir/dxmddelta
    !        d2Psir_dxmddel(m) = R_mix / props_phase(a)%tau_it * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxiddel_it(m) + props_phase(a)%dalphar_dxi_it(m)))
    !
    !        !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and tau
    !        !d2Psir/dxmdtau
    !        d2Psir_dxmdtau(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it**2 * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxidtau_it(m) - props_phase(a)%dalphar_dxi_it(m)))
    !
    !        do n = m, gl%ncomp
    !            !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and xn
    !            !d2Psir/dxmdxn
    !            d2Psir_dxmdxn(m,n) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * &
    !                               & (props_phase(a)%alphar_it(1) * d2rhorTr_dxmdxn(m,n) + props_phase(a)%dalphar_dxi_it(n) * drhorTr_dxm(m) + &
    !                               &  gl%rhoredmix * gl%tredmix * props_phase(a)%d2alphar_dxidxj_it(m,n) + props_phase(a)%dalphar_dxi_it(m) * drhorTr_dxm(n))
    !            if (n /= m) then
    !                d2Psir_dxmdxn(n,m) = d2Psir_dxmdxn(m,n)
    !            end if
    !        end do
    !    end do
    !
    !
    !    !Calculate the second derivatives of the residual Helmholtz energy density with respect to delta, tau, xm and rho_j
    !    !d(dPsir/ddelta)/drho_j, d(dPsir/dtau)/drho_j, d(dPsir/dxm)/drho_j
    !    do j = 1, gl%ncomp
    !        ddPsirddelta_drhoj(j) = d2Psir_ddelta2 * ddelta_drhoi(j) + d2Psir_dtauddelta * dtau_drhoi(j)
    !        ddPsirdtau_drhoj(j) = d2Psir_dtauddelta * ddelta_drhoi(j) + d2Psir_dtau2 * dtau_drhoi(j)
    !        do m = 1, gl%ncomp
    !            ddPsirddelta_drhoj(j) = ddPsirddelta_drhoj(j) + d2Psir_dxmddel(m) * dxm_drhoi(m,j)
    !            ddPsirdtau_drhoj(j) = ddPsirdtau_drhoj(j) + d2Psir_dxmdtau(m) * dxm_drhoi(m,j)
    !            ddPsirdxm_drhoj(m,j) = d2Psir_dxmddel(m) * ddelta_drhoi(j) + d2Psir_dxmdtau(m) * dtau_drhoi(j)
    !            do n = 1, gl%ncomp
    !                ddPsirdxm_drhoj(m,j) = ddPsirdxm_drhoj(m,j) + d2Psir_dxmdxn(m,n) * dxm_drhoi(n,j)
    !            end do
    !        end do
    !    end do
    !
    !
    !    !Calculate the second derivatives of Psi_r with respect to the molar concentrations rho_i and rho_j
    !    !(d2Psi_r / (drho_i drho_j))
    !    Do i = 1, gl%ncomp
    !        do j = i, gl%ncomp
    !            d2Psir_drhoidrhoj(i,j) = dPsir_ddelta * d2delta_drhoidrhoj(i,j) + ddPsirddelta_drhoj(j) * ddelta_drhoi(i) + dPsir_dtau * d2tau_drhoidrhoj(i,j) + ddPsirdtau_drhoj(j) * dtau_drhoi(i)
    !            do m = 1, gl%ncomp
    !                d2Psir_drhoidrhoj(i,j) = d2Psir_drhoidrhoj(i,j) + dPsir_dxm(m) * d2xm_drhoidrhoj(m,i,j) + ddPsirdxm_drhoj(m,j) * dxm_drhoi(m,i)
    !            end do
    !            if (i /= j) then
    !                d2Psir_drhoidrhoj(j,i) = d2Psir_drhoidrhoj(i,j)
    !            end if
    !        end do
    !    End do
    !    !-------------------------------------------------------------------------------------------------------------
    !
    !    !Calculate the Hessian matrix
    !    do i = 1, gl%ncomp
    !        do j = i, gl%ncomp
    !            if (i == j) then
    !                Hessian(i,j,a) = d2Psir_drhoidrhoj(i,j) + R_mix * Temp / gl%molfractions(i) / Dens
    !            else
    !                Hessian(i,j,a) = d2Psir_drhoidrhoj(i,j)
    !                Hessian(j,i,a) = Hessian(i,j,a)
    !            end if
    !        end do
    !    end do
    !
    !
    !    !Calculate the second derivatives of Psi with respect to the molar concentration rho_i and T
    !    !-------------------------------------------------------------------------------------------------------------
    !    ! Calculate help variables
    !    ! dtau/dT, d2tau/drhoidT (Variables not explicitly given in paper by Bell and Deiters?)
    !    dtau_dT = - gl%tredmix / Temp**2
    !    do i = 1, gl%ncomp
    !        d2tau_drhoi_dT(i) = -props_phase(a)%dTr_drhoi_it(i) / Temp**2
    !    end do
    !
    !    !Calculate the second derivatives of Psi_r with respect to T and molar concentration rhoi
    !    !(d2Psi_r / (drho_i dT))
    !    do i = 1, gl%ncomp
    !        d2Psir_drhoidT(i) = dPsir_dtau * d2tau_drhoi_dT(i) + d2Psir_dtauddelta * ddelta_drhoi(i) * dtau_dT + &
    !                          & d2Psir_dtau2 * dtau_drhoi(i) * dtau_dT
    !        do m = 1, gl%ncomp
    !            d2Psir_drhoidT(i) = d2Psir_drhoidT(i) + d2Psir_dxmdtau(m) * dxm_drhoi(m,i) * dtau_dT
    !        end do
    !    end do
    !    !-------------------------------------------------------------------------------------------------------------
    !
    !end do
    !
    !
    !!Calculate the derivative of the molar concentrations with respect to pressure along the bubble line
    !!Eq. (19) in the article of Bell and Deiters
    !!-------------------------------------------------------------------------------------------------------------
    !!Set up right side of the system of equations
    !do i = 1, gl%ncomp
    !    vectorb(i) = 1.D0
    !end do
    !!Set up Matrix A (FOR 2 COMPONENTS ONLY)
    !MatrixA(1,1) = Hessian(1,1,2) * dens_i_phase(1,1) + Hessian(1,2,2) * dens_i_phase(2,1)
    !MatrixA(1,2) = Hessian(2,1,2) * dens_i_phase(1,1) + Hessian(2,2,2) * dens_i_phase(2,1)
    !MatrixA(2,1) = Hessian(1,1,2) * dens_i_phase(1,2) + Hessian(1,2,2) * dens_i_phase(2,2)
    !MatrixA(2,2) = Hessian(2,1,2) * dens_i_phase(1,2) + Hessian(2,2,2) * dens_i_phase(2,2)
    !rankA = gl%ncomp
    !det_A = 0.D0
    !call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errval)
    !drhoi_dp(1,2) = vectorx(1)
    !drhoi_dp(2,2) = vectorx(2)
    !!-------------------------------------------------------------------------------------------------------------
    !
    !
    !!Calculate the derivative of the molar concentrations with respect to pressure along the dew line
    !!Eq. (20) in the article of Bell and Deiters
    !!-------------------------------------------------------------------------------------------------------------
    !!Set up right side of the system of equations (FOR 2 COMPONENTS ONLY)
    !vectorb(1) = Hessian(1,1,2) * drhoi_dp(1,2) + Hessian(1,2,2) * drhoi_dp(2,2)
    !vectorb(2) = Hessian(2,1,2) * drhoi_dp(1,2) + Hessian(2,2,2) * drhoi_dp(2,2)
    !!Set up Matrix A (Matrix A is simply the Hessian of the vapor phase in this case)
    !MatrixA(1,1) = Hessian(1,1,1)
    !MatrixA(1,2) = Hessian(1,2,1)
    !MatrixA(2,1) = Hessian(2,1,1)
    !MatrixA(2,2) = Hessian(2,2,1)
    !rankA = gl%ncomp
    !det_A = 0.D0
    !call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errval)
    !drhoi_dp(1,1) = vectorx(1)
    !drhoi_dp(2,1) = vectorx(2)
    !!-------------------------------------------------------------------------------------------------------------
    !
    !ind_var = 0.D0
    !ind_var(1) = Temp       !First independent variable needed is the temperature
    !nr_of_eqs = 4
    !rhoi_march = dens_i_phase(1,2)
    !y_vector(1) = dens_i_phase(1,1)
    !y_vector(2) = dens_i_phase(2,1)
    !y_vector(3) = dens_i_phase(1,2)
    !y_vector(4) = dens_i_phase(2,2)
    !dy_dt(:) = drhoi_drhoj(gl, nr_of_eqs, rhoi_march, y_vector, ind_var, errval)



    end subroutine px_diag_isochoric




    !RUNGE-KUTTA-FEHLBERG METHOD FOR SOLVING A FIRST ORDER ORDINARY DIFFERENTIAL EQUATION
    !GIVEN A STARTING VALUE y0(t0)
    !Andreas Jäger, November 2018
    !!--------------------------------------------------------------------------------
    subroutine RKF45_ODE_INT(gl, y_deriv, ind_var, post_proc, nr_of_eqs, y_min, t_min, t_max, step_min, eps_max, fac_rel, nr_of_points, t_res, y_res, res_var, errorflag)
    !This subroutines approximates an ordinary (system of) partial differental equation(s) (ODE)
    !by means of an adaptive Runge-Kutta-Fehlberg 4th/5th order method
    !
    !INPUT:     y_deriv     -   An external function which is passed to the routine. The function must return a vector of "nr_of_eqs" derivatives y' at given inputs y and t (y' = f(t,y))
    !           ind_var     -   Other independent variables that might be needed to evaluate the derivatives of y with respect to t (for the external function y_deriv)
    !           post_proc   -   Boolean that indicates if a post processing function shall be used in order to improve the obtained values for y(t):
    !                           post_proc = true: USE post-processing routine, post_proc = false: DO NOT USE post-processing routine
    !                           As it is implemented, the function y_deriv also serves as possible post-processing routine. Any RKF-function might be supplemented with an option to improve values for y instead of calculating y' = f(t,y), if post_proc = true
    !           y_min       -   Vector which contains the start values y_0(t_0) (here: y_min(t_min))
    !           nr_of_eqs   -   The number of equations the ODE consists of
    !           t_min       -   Starting value of the variable t
    !           t_max       -   End value of the variable t
    !           step_min    -   This is the minimum step size allowed. The algorithm starts with a multiple of the allowed minimum step size. It terminates, when the adaptive step size correction reduces the step below step_min
    !           eps_max     -   maximum error allowed (between 4th and 5th order Runge-Kutta) (around 10^-10)
    !           fac_rel     -   Relaxation factor for the step (something close to one. Smaller than one will lead to more careful steps)
    !
    !OUTPUT:    nr_of_points    -   nr_of_points calculated during the iteration process
    !           t_res           -   vector with results for the variable t at each iteration point
    !           y_res           -   Matrix containing the values of the "nr_of_eqs" function(s) y(t)
    !           errorflag       -   In case of an error, this variable is set to a value not equal to 0
    !!--------------------------------------------------------------------------------



    implicit none

    interface
    function y_deriv(gl, nr_of_eqs, t_march, y_vector, ind_var, out_var, post_proc, errorflag)  !The ODE (y'(T) = f(t,y)) which must be passed to the subroutine as an external function
    use module_all_types
    implicit none
    type(type_gl) :: gl
    integer :: nr_of_eqs                                    !Number of equations the ODE consists of
    double precision, dimension(nr_of_eqs) :: y_vector      !Values of y(t_march)
    double precision :: t_march                             !Marching variable t
    double precision, dimension(10) :: ind_var              !Other independent variables that might be needed to evaluate the derivatives of y with respect to t
    double precision, dimension(10) :: out_var              !Other (optional) output variables
    logical :: post_proc                                    !True: Run function as post-processing function, false: run function in normal mode
    integer :: errorflag                                    !Error value
    double precision, dimension(nr_of_eqs) :: y_deriv       !Output vector of derivatives
    end function
    end interface



    type(type_gl) :: gl

    !Input variables
    integer :: nr_of_eqs
    double precision, dimension(10) :: ind_var          !Other independent variables that might be needed to evaluate the derivatives of y with respect to t (for the external function y_deriv)
    logical :: post_proc
    double precision, dimension(nr_of_eqs) :: y_min
    double precision :: t_min
    double precision :: t_max
    double precision :: step_min                        !Minimum step size allowed
    double precision :: eps_max
    double precision :: fac_rel

    !Output variables
    integer :: nr_of_points
    double precision, dimension (1000) :: t_res           !List of variables t at which y has been evaluated
    double precision, dimension (4,1000) :: y_res         !Results for y evaluated at t
    double precision, dimension (10,1000) :: res_var       !Additional result variables
    integer :: errorflag

    !Other variables
    double precision :: t_march                         !actual variable of marching variable t
    double precision, dimension(nr_of_eqs) :: y_march   !Actual value of the variables y according to RK method of fourth order
    double precision, dimension(nr_of_eqs) :: z_march   !Actual value of the variables z according to RK method of fifth order
    double precision, dimension(nr_of_eqs) :: dummy_march   !dummy variable for post-processing
    integer :: i, j                                     !Iteration number
    double precision :: step_h                          !Step size for changing t
    double precision :: step_max                        !Maximum step size
    double precision, dimension(nr_of_eqs) :: var_k1    !Variable "k1" of the Runge-Kutta-Fehlberg method
    double precision, dimension(nr_of_eqs) :: var_k2    !Variable "k2" of the Runge-Kutta-Fehlberg method
    double precision, dimension(nr_of_eqs) :: var_k3    !Variable "k3" of the Runge-Kutta-Fehlberg method
    double precision, dimension(nr_of_eqs) :: var_k4    !Variable "k4" of the Runge-Kutta-Fehlberg method
    double precision, dimension(nr_of_eqs) :: var_k5    !Variable "k5" of the Runge-Kutta-Fehlberg method
    double precision, dimension(nr_of_eqs) :: var_k6    !Variable "k6" of the Runge-Kutta-Fehlberg method
    double precision :: t_k2                            !Changed variable t (t = t + 0.25*h) for calculating k2
    double precision :: t_k3                            !Changed variable t (t = t + 0.375*h) for calculating k3
    double precision :: t_k4                            !Changed variable t (t = t + 12/13*h) for calculating k4
    double precision :: t_k5                            !Changed variable t (t = t + h) for calculating k5
    double precision :: t_k6                            !Changed variable t (t = t + 0.5*h) for calculating k6
    double precision, dimension(nr_of_eqs) :: y_k2      !Changed variable y (y = y + 0.25 * k1) for calculating k2
    double precision, dimension(nr_of_eqs) :: y_k3      !Changed variable y (y = y + 3/32*k1 + 9/32*k2) for calculating k3
    double precision, dimension(nr_of_eqs) :: y_k4      !Changed variable y (y = y + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3) for calculating k4
    double precision, dimension(nr_of_eqs) :: y_k5      !Changed variable y (y = y + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4) for calculating k5
    double precision, dimension(nr_of_eqs) :: y_k6      !Changed variable y (y = y - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5) for calculating k6
    double precision:: s_stepmod                        !Stepsize control s that depends on the difference between fourth (y_march) and fifth (z_march) order Runge-Kutta
    double precision:: error_zy                         !Difference between 4th and 5th order Runge-Kutta
    logical :: step_accepted                            !The step of the Runge-Kutta method is accepted if the error is smaller than eps_max
    integer :: count_step_min                           !Count how many times the step was smaller than step_min. If this happens a certain number of times, stop the algorithm
    double precision, dimension(10) :: out_var          !Other (optional) output variables
    logical :: post_proc_RKF                            !For the RKF steps itself, post-processing needs to be disabled -> this variable is always .false.

    !Allocate variables
    !allocate(t_res(10000))
    !allocate(y_res(nr_of_eqs,10000))
    !allocate(res_var(10, 10000))

    !Initialize variables
    out_var = 0.D0
    count_step_min = 0

    !Set the marching variable to the minimum value
    t_march = t_min
    !Set the values of y and z to the start values
    y_march(:) = y_min(:)
    z_march(:) = y_min(:)
    !Begin with a step size 100 times the minimum step size
    step_h = step_min * 100.D0
    !Define the maximum step size
    step_max = (t_max - t_min) / 200.D0     !Set maximum step size such that a minimum of 200 steps (arbitrarily picked) are carried out
    !Write the solution to the marching and result variables
    t_res(1) = t_march          !Save the marching variable at step i
    y_res(:,1) = z_march(:)     !Use fifth order solution z
    y_march(:) = z_march(:)     !Use fifth order solution z
    !y_res(:,1) = y_march(:)    !Use fourth order solution y
    !y_march(:) = y_march(:)    !Use fourth order solution y


    !!Check whether the specified step size allows for coming from t_min to t_max in 10000 steps, if not, quit with error
    !if ((t_min + 10000 * step_min) < t_max) then
    !    errorflag = -20000
    !    return
    !end if

    !No post-processing in RKF routine
    post_proc_RKF = .false.

    i = 1
    do while ((t_march + step_h) < t_max)

        i = i + 1           !Increase iteration count
        !Break criterion of i is larger than the maximum number of iteration steps
        if (i > 1000) then
            return
        end if

        nr_of_points = i

        !Set step accepted by default to false
        step_accepted = .false.
        do while (step_accepted == .false.)

            !Calculate the derivatives at t and y and store in k1 variable
            var_k1(:) = step_h * y_deriv(gl, nr_of_eqs, t_march, y_march, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Calcuate the derivatives at t_k2 and y_k2 and store in k2 variable
            t_k2 = t_march + 0.25D0 * step_h
            y_k2(:) = y_march(:) + 0.25D0 * var_k1(:)
            var_k2(:) = step_h * y_deriv(gl, nr_of_eqs, t_k2, y_k2, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Calcuate the derivatives at t_k3 and y_k3 and store in k3 variable
            t_k3 = t_march + 0.375D0 * step_h
            y_k3(:) = y_march(:) + 3.D0/32.D0 * var_k1(:) + 9.D0/32.D0 * var_k2(:)
            var_k3(:) = step_h * y_deriv(gl, nr_of_eqs, t_k3, y_k3, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Calcuate the derivatives at t_k4 and y_k4 and store in k4 variable
            t_k4 = t_march + 12.D0/13.D0 * step_h
            y_k4(:) = y_march(:) + 1932.D0/2197.D0 * var_k1(:) - 7200.D0/2197.D0 * var_k2(:) + 7296.D0/2197.D0 * var_k3(:)
            var_k4(:) = step_h * y_deriv(gl, nr_of_eqs, t_k4, y_k4, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Calcuate the derivatives at t_k5 and y_k5 and store in k5 variable
            t_k5 = t_march + step_h
            y_k5(:) = y_march(:) + 439.D0/216.D0 * var_k1(:) - 8.D0 * var_k2(:) + 3680.D0/513.D0 * var_k3(:) - 845.D0/4104.D0 * var_k4(:)
            var_k5(:) = step_h * y_deriv(gl, nr_of_eqs, t_k5, y_k5, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if
            !Calcuate the derivatives at t_k5 and y_k5 and store in k5 variable
            t_k6 = t_march + 0.5D0 * step_h
            y_k6(:) = y_march(:) - 8.D0/27.D0 * var_k1(:) + 2.D0 * var_k2(:) - 3544.D0/2565.D0 * var_k3(:) + 1859.D0/4104.D0 * var_k4(:) - 11.D0/40.D0 * var_k5(:)
            var_k6(:) = step_h * y_deriv(gl, nr_of_eqs, t_k6, y_k6, ind_var, out_var, post_proc_RKF, errorflag)
            if (errorflag /= 0) then
                return
            end if

            !Fourth order Runge-Kutta method
            y_march(:) = y_march(:) + 25.D0/216.D0 * var_k1(:) + 1408.D0/2565.D0 * var_k3(:) + 2197.D0/4104.D0 * var_k4(:) - 1.D0/5.D0 * var_k5(:)

            !Fifth order Runge-Kutta method
            z_march(:) = z_march(:) + 16.D0/135.D0 * var_k1(:) + 6656.D0/12825.D0 * var_k3(:) + 28561.D0/56430.D0 * var_k4(:) - 9.D0/50.D0 * var_k5(:) + 2.D0/55.D0 * var_k6(:)

            !Calculate the relative error (difference between z_march and y_march divided by z_march. In case of more than one equation, compute the lengths of the relative difference vector)
            error_zy = 0.D0
            if (count(dabs(z_march(:)) > 0.D0) == nr_of_eqs) then
                !Calculate relative differences if all z_march are not equal to 0
                do j = 1, nr_of_eqs
                    error_zy = error_zy + ((z_march(j) - y_march(j))/z_march(j))**2
                end do
            else
                !Calculate absolute differences if any z_march(j) is equal to 0
                do j = 1, nr_of_eqs
                    error_zy = error_zy + (z_march(j) - y_march(j))**2
                end do
            end if
            error_zy = error_zy**0.5

            !Adjust step size
            if (error_zy > eps_max) then
                s_stepmod = fac_rel * (eps_max/error_zy)**0.3
                step_h = step_h * s_stepmod
                step_accepted = .false.
                !Set back marching variables
                y_march(:) = y_res(:,i-1)
                z_march(:) = y_res(:,i-1)
                !t_march = t_res(i-1) + step_h
            else
                step_accepted = .true.
            end if

        end do

        !Write the solution to the marching and result variables
        t_res(i) = t_march + step_h         !Save the marching variable at step i
        y_res(:,i) = z_march(:)             !Use fifth order solution z
        y_march(:) = z_march(:)             !Use fifth order solution z
        !y_res(:,1) = y_march(:)            !Use fourth order solution y
        !y_march(:) = y_march(:)            !Use fourth order solution y

        !Increase marching variable
        t_march = t_march + step_h

        !Post processing function. If applicable, get better values for y after the RKF-step here (i.e. solve the phase equilibrium conditions)
        !This is for example needed for drawing px (or Tx) diagrams
        if (post_proc == .true.) then
            dummy_march(:) = y_deriv(gl, nr_of_eqs, t_march, y_march, ind_var, out_var, post_proc, errorflag)
            if (errorflag /= 0) then
                return
            end if
            res_var(:,i) = out_var(:)
            y_res(:,i) = y_march(:)
        end if

        !Increase the step size, if the last step was accepted
        if ((error_zy < eps_max) .and. (error_zy > 0.D0)) then
            s_stepmod = fac_rel * (eps_max/error_zy)**0.2
            !If the step size modification factor gets too large, decrease it
            if (s_stepmod > 5.D0) then
                s_stepmod = 5.D0
            end if
            step_h = step_h * s_stepmod
        end if

        !Test if the step is bigger than the maximum step size allowed. If yes, set the step size to the maximum allowed step size
        if (step_h > step_max) then
            step_h = step_max
        end if

        !If the step is smaller than the minimum step specified, return
        if (step_h < step_min) then
            count_step_min = count_step_min + 1
            if (count_step_min > 3) then
                return
            end if
        end if

    end do


    end subroutine RKF45_ODE_INT



    !Function for calculating the derivatives of the molar concentration rho_i with respect to a certain molar concentration
    !rho_j (marching variable) or pressure (marching variable) along the phase boundary for a binary mixture and for two phases
    !Andreas Jäger, November 2018
    function drhoi_drhoj(gl, nr_of_eqs, var_march, y_vector, ind_var, out_var, post_proc, errval)
    !--------------------------------------------------------------------------------
    ! This function calculates the following derivatives along the phase boundary:
    !
    ! If ind_var(2) = 1,2,3,4: --> Marching in a molar concentration
    ! Phase 1: drho1/drhoi, drho2/drhoi
    ! Phase 2: drho1/drhoi, drho2/drhoi
    ! where rho_i = x_i * rho_phase
    ! If ind_var(2) = 5: --> Marching in pressure
    ! Phase 1: drho1/p, drho2/p
    ! Phase 2: drho1/p, drho2/p
    !
    ! The basic algorithm and method is described in:
    ! Bell, I.H.; Deiters U.K.: "On the Construction of Binary Mixture p-x and T-x Diagrams from Isochoric Thermodynamics", AIChe Journal 64(7), 2745-2757, 2018.
    ! The derivatives calculated in this function correspond to Eq. (40) for ind_var(2) = 1,2,3,4 and Eq. (39) for ind_var = 5 of the article
    !
    ! This function can be used as an input to the Runge-Kutta-Fehlberg method RKF45_ODE_INT()
    !
    ! Input:    nr_of_eqs       -   number of equations (here: for a binary mixture 4 equations are needed)
    !           var_march       -   specified value for marching variable (rho_i or pressure) (see, ind_var(2) = 1 -> drho_1_vap, ind_var(2) = 2 -> drho_2_vap, ind_var(2) = 3 -> drho_1_liq, ind_var(2) = 4 -> drho_2_liq, ind_var(2) = 5 -> pressure)
    !           y_vector        -   the values for rho_1_vap, rho_2_vap, rho_1_liq, and rho_2_liq
    !           ind_var         -   Other independent variables that might be needed to evaluate the derivatives of y with respect to t
    !                               here:       INPUT:
    !                                           ind_var(1) = temperature
    !                                           ind_var(2) = specifier which derivative is needed (1 = drhoi / drho_1_vap, 2 = drhoi / drho_2_vap, 3 = drhoi / drho_1_liq, 4 = drhoi / drho_2_liq, 5 = drho_i / dp)
    !
    !           out_var         -   Additional output variables. Here: out_var(1) = pressure
    !           post_proc       -   True: Run function as post-processing function, false: run function in normal mode. (here: Normal mode solves y' = f(t,y), post-processing solves the the phase equilibrium condition for the molar concentration rho1_vap, rho2_vap, rho1_liq, and rho2_liq)
    !
    ! Output:   drhoi_drhoj  -  The function itself, returning the derivatives as follows:
    !                               ind_var(2) = 1,2,3,4:
    !                               drhoi_drhoj(1) = drho1_vap/drho_specified (see ind_var(2))
    !                               drhoi_drhoj(2) = drho2_vap/drho_specified (see ind_var(2))
    !                               drhoi_drhoj(3) = drho1_liq/drho_specified (see ind_var(2))
    !                               drhoi_drhoj(4) = drho2_liq/drho_specified (see ind_var(2))
    !                               ind_var(2) = 5:
    !                               drhoi_drhoj(1) = drho1_vap/dp
    !                               drhoi_drhoj(2) = drho2_vap/dp
    !                               drhoi_drhoj(3) = drho1_liq/dp
    !                               drhoi_drhoj(4) = drho2_liq/dp
    !
    !           errval          -   specifies if an error occured during the calculations
    !--------------------------------------------------------------------------------
    ! A. Jäger, November 2018




    implicit none

    !Input and output variables of the routine
    type(type_gl) :: gl
    integer:: nr_of_eqs
    double precision:: var_march
    double precision, dimension(nr_of_eqs) :: y_vector
    double precision, dimension(10) :: ind_var              !Other independent variables that might be needed to evaluate the derivatives of y with respect to t (here: ind_var(1) = temperature)
    double precision, dimension(10) :: out_var              !Other (optional) output variables
    logical :: post_proc
    integer:: errval
    double precision, dimension(nr_of_eqs) :: drhoi_drhoj   !(the function itself)

    !Create type in which all needed properties of one phase for the isochoric tracing are stored
    type(type_isoch_therm), dimension(2) :: props_phase

    !Derivatives of the residual Helmholtz energy density (Psi^r) with respect to tau, delta, x and help variables
    double precision :: Psir
    double precision :: dPsir_ddelta
    double precision :: dPsir_dtau
    double precision, dimension(:), allocatable :: dPsir_dxm
    double precision :: d2Psir_ddelta2
    double precision :: d2Psir_dtauddelta
    double precision :: d2Psir_dtau2
    double precision, dimension(:), allocatable :: d2Psir_dxmddel
    double precision, dimension(:), allocatable :: d2Psir_dxmdtau
    double precision, dimension(:,:), allocatable :: d2Psir_dxmdxn
    !Help variables
    double precision, dimension(:,:), allocatable :: kron_delta
    double precision, dimension(:), allocatable :: drhorTr_dxm
    double precision, dimension(:,:), allocatable :: d2rhorTr_dxmdxn
    double precision, dimension(:), allocatable :: ddelta_drhoi
    double precision, dimension(:), allocatable :: dtau_drhoi
    double precision, dimension(:,:), allocatable :: dxm_drhoi
    double precision, dimension(:), allocatable :: d2tau_drhoi_dT
    double precision :: dtau_dT
    double precision, dimension(:,:), allocatable :: d2tau_drhoidrhoj
    double precision, dimension(:,:), allocatable :: d2delta_drhoidrhoj
    double precision, dimension(:,:,:), allocatable :: d2xm_drhoidrhoj

    !Derivatives of the residual Helmholtz energy density (Psi^r) with respect to independent variables
    !double precision, dimension(:), allocatable :: dPsi0_drhoi              !First derivative of the ideal part of Psi with respect to rho_i
    double precision, dimension(:), allocatable :: dPsir_drhoi              !First derivative of the residual part of Psi with respect to rho_i
    !double precision, dimension(:), allocatable :: d2Psi0_drhoi_dT          !Second derivative of the ideal part of Psi with respect to rho_i and T
    double precision, dimension(:), allocatable :: d2Psir_drhoidT          !Second derivative of the residual part of Psi with respect to rho_i and T
    !double precision, dimension(:,:), allocatable :: d2Psi0_drhoi_drhoj     !Second derivative of the ideal part of Psi with respect to rho_i and rho_j
    double precision, dimension(:,:), allocatable :: d2Psir_drhoidrhoj     !Second derivative of the residual part of Psi with respect to rho_i and rho_j

    !Mixed derivatives of the residual Helmholtz energy density with respect to rhoi, delta, tau, and xi
    double precision, dimension(:), allocatable :: ddPsirddelta_drhoj
    double precision, dimension(:), allocatable :: ddPsirdtau_drhoj
    double precision, dimension(:,:), allocatable :: ddPsirdxm_drhoj

    !Other variables
    double precision :: Dens, Temp
    double precision :: R_mix
    double precision, dimension(2) :: press_phase
    integer :: a, i, j, k, m, n

    !Variables for tracing
    double precision, dimension(:,:,:), allocatable:: Hessian       !Last entry: Number of the phase
    double precision, dimension(2) :: dens_phase                    !densities of the phases: dens_phase(1): vapor phase, dens_phase(2): liquid phase
    double precision, dimension(:,:), allocatable :: dens_i_phase   !Molar concentratons of the phases: dens_i_phase(1:ncomp,1): vapor phase, dens_i_phase(1:ncomp,2): liquid phase
    double precision, dimension(:,:), allocatable :: x_phase        !compositions of the phases: x_phase(1:ncomp,1): vapor phase,x_phase(1:ncomp,2): liquid phase
    double precision, dimension(:,:), allocatable :: drhoi_dp       !Derivative of molar concentrations of the different phases along the phase boundary. drhoi_dp(1:ncomp,1): vapor phase,drhoi_dp(1:ncomp,2): liquid phase
    integer :: spec_deriv                                           !specifier which derivative is needed (1 = drhoi / drho_1_vap, 2 = drhoi / drho_2_vap, 3 = drhoi / drho_1_liq, 4 = drhoi / drho_2_liq)

    !Variables for Gauss-algorithm. ADJUST THE SIZE DYNAMICALLY!!!!
    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision, dimension(60):: vectorb
    double precision, dimension(60):: vectorx
    double precision:: det_A

    !Variables for the "post-processing" step
    double precision, dimension(4) :: Sys_of_Eqs    !Residual values of the system of equations
    double precision, dimension(4,4) :: Jac_Mat     !Jacobian Matrix
    double precision, dimension(4) :: Delta_X       !Newton-Raphson step
    double precision, dimension(:,:), allocatable :: dp_drhoj                       !Derivative of p with respect to the molar concentration rho_j at constant T and rho_K. dp_drhoj(1:ncomp,1) : vapor phase, dp_drhoj(1:ncomp,2) : liquid phase
    double precision, dimension(:,:), allocatable :: dPsir_drhoi_phase              !Variable for saving the first derivative of the residual part of Psi with respect to rho_i
    double precision, dimension(:,:,:), allocatable :: d2Psir_drhoidrhoj_phase      !Variable for saving the second derivative of the residual part of Psi with respect to rho_i and rho_j
    double precision, dimension(2) :: R_mix_phase     !Variable that is necessary because it might happen that different equations of state use different values for the universal gas constant. In these cases, the gas constant is mixed with the mole fractions
    double precision :: eps_res                     !Maximum residuum allowed
    double precision :: eps_deltax                  !Maximum change of iteration variables allowed



    do a = 1, 2
        !Allocate generated type props_phase and define size of variables for both phases
        if(.not.allocated(props_phase(a)%alpha_it))       allocate(props_phase(a)%alpha_it(nderivs))
        if(.not.allocated(props_phase(a)%alpha0_it))      allocate(props_phase(a)%alpha0_it(nderivs))
        if(.not.allocated(props_phase(a)%alphar_it))      allocate(props_phase(a)%alphar_it(nderivs))

        if(.not.allocated(props_phase(a)%dTr_drhoi_it))           allocate(props_phase(a)%dTr_drhoi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%dTr_dxi_it))             allocate(props_phase(a)%dTr_dxi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2Tr_drhoidrhoj_it))     allocate(props_phase(a)%d2Tr_drhoidrhoj_it(gl%ncomp,gl%ncomp))
        if(.not.allocated(props_phase(a)%d2Tr_dxidxj_it))         allocate(props_phase(a)%d2Tr_dxidxj_it(gl%ncomp,gl%ncomp))

        if(.not.allocated(props_phase(a)%drhor_drhoi_it))         allocate(props_phase(a)%drhor_drhoi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%drhor_dxi_it))           allocate(props_phase(a)%drhor_dxi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2rhor_drhoidrhoj_it))   allocate(props_phase(a)%d2rhor_drhoidrhoj_it(gl%ncomp,gl%ncomp))
        if(.not.allocated(props_phase(a)%d2rhor_dxidxj_it))       allocate(props_phase(a)%d2rhor_dxidxj_it(gl%ncomp,gl%ncomp))

        if(.not.allocated(props_phase(a)%dalpha_dxi_it))          allocate(props_phase(a)%dalpha_dxi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alpha_dxiddel_it))     allocate(props_phase(a)%d2alpha_dxiddel_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alpha_dxidtau_it))     allocate(props_phase(a)%d2alpha_dxidtau_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alpha_dxidxj_it))      allocate(props_phase(a)%d2alpha_dxidxj_it(gl%ncomp,gl%ncomp))

        if(.not.allocated(props_phase(a)%dalphar_dxi_it))          allocate(props_phase(a)%dalphar_dxi_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alphar_dxiddel_it))     allocate(props_phase(a)%d2alphar_dxiddel_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alphar_dxidtau_it))     allocate(props_phase(a)%d2alphar_dxidtau_it(gl%ncomp))
        if(.not.allocated(props_phase(a)%d2alphar_dxidxj_it))      allocate(props_phase(a)%d2alphar_dxidxj_it(gl%ncomp,gl%ncomp))
    end do

    allocate(dPsir_dxm(gl%ncomp))
    allocate(d2Psir_dxmddel(gl%ncomp))
    allocate(d2Psir_dxmdtau(gl%ncomp))
    allocate(d2Psir_dxmdxn(gl%ncomp,gl%ncomp))
    allocate(kron_delta(gl%ncomp,gl%ncomp))
    allocate(drhorTr_dxm(gl%ncomp))
    allocate(d2rhorTr_dxmdxn(gl%ncomp,gl%ncomp))
    allocate(ddelta_drhoi(gl%ncomp))
    allocate(dtau_drhoi(gl%ncomp))
    allocate(dxm_drhoi(gl%ncomp,gl%ncomp))
    allocate(d2tau_drhoi_dT(gl%ncomp))
    allocate(d2tau_drhoidrhoj(gl%ncomp,gl%ncomp))
    allocate(d2delta_drhoidrhoj(gl%ncomp,gl%ncomp))
    allocate(d2xm_drhoidrhoj(gl%ncomp,gl%ncomp,gl%ncomp))

    !allocate(dPsi0_drhoi(gl%ncomp))
    allocate(dPsir_drhoi(gl%ncomp))
    !allocate(d2Psi0_drhoidrhoj(gl%ncomp,gl%ncomp))
    allocate(d2Psir_drhoidrhoj(gl%ncomp,gl%ncomp))
    allocate(d2Psir_drhoidT(gl%ncomp))

    allocate(ddPsirddelta_drhoj(gl%ncomp))
    allocate(ddPsirdtau_drhoj(gl%ncomp))
    allocate(ddPsirdxm_drhoj(gl%ncomp,gl%ncomp))

    allocate(Hessian(gl%ncomp,gl%ncomp,2))
    allocate(dens_i_phase(gl%ncomp,2))
    allocate(x_phase(gl%ncomp,2))
    allocate(drhoi_dp(gl%ncomp,2))

    allocate(dp_drhoj(gl%ncomp,2))
    allocate(dPsir_drhoi_phase(gl%ncomp,2))
    allocate(d2Psir_drhoidrhoj_phase(gl%ncomp,gl%ncomp,2))

    !Initialize variables
    errval = 0
    !ndTred_dni = 0.D0
    !ndrhored_dni = 0.D0

    !Fill the Kronecker delta
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            if (j == i) then
                kron_delta(i,j) = 1.D0
            else
                kron_delta(i,j) = 0.D0
            end if
        end do
    end do

    !Check whether the specified derivative is an valid input and set spec_deriv accordingly
    if (ind_var(2) < 0.D0 .or. ind_var(2) > 5.D0) then
        errval = -20000
        return
    else
        spec_deriv = ind_var(2)
    end if

    !Write the temperature in the temp variable
    Temp = ind_var(1)


    If (post_proc .eqv. .false.) then   !"Normal routine (y'(t) = f(t,y(t))) is executed

        !Write the y-function in the corresponding variables
        dens_i_phase(1,1) = y_vector(1)
        dens_i_phase(2,1) = y_vector(2)
        dens_i_phase(1,2) = y_vector(3)
        dens_i_phase(2,2) = y_vector(4)

        !Density of the vapor phase (phase 1)
        dens_phase(1) = y_vector(1) + y_vector(2)
        !Calculate the composition of the vapor phase (phase 1)
        x_phase(1,1) = y_vector(1) / dens_phase(1)
        x_phase(2,1) = y_vector(2) / dens_phase(1)

        !Density of the liquid phase (phase 2)
        dens_phase(2) = y_vector(3) + y_vector(4)
        !Calculate the composition of the vapor phase (phase 1)
        x_phase(1,2) = y_vector(3) / dens_phase(2)
        x_phase(2,2) = y_vector(4) / dens_phase(2)


        !Compute the needed derivatives and the Hessians of the phases in equilibrium
        !1: Vapor phase, 2: Liquid phase
        do a = 1, 2

            !Set density and composition of the phase and get all needed derivatives
            Dens = dens_phase(a)
            gl%molfractions = 0.D0
            gl%molfractions(1:2) = x_phase(:,a)
            call reduced_parameters_calc(gl,Temp)

            !Get the mixture gas constant
            call R_mix_calc(gl, R_mix)

            call Derivs_isoch_therm(gl, props_phase(a), Temp, Dens, errval)

            !Calculate the residual Helmholtz energy density Psi^r = a^r * rho
            Psir = Dens * R_mix * Temp * props_phase(a)%alphar_it(1)

            !Calculate the first derivative of Psi with respect to the molar concentration (i.e. the chemical potentials)
            !-------------------------------------------------------------------------------------------------------------
            !Calculate help variables
            ! ddelta/drhoi, dtau/drhoi, dxm/drhoi, d(rhor*Tr)/dxm
            do i = 1, gl%ncomp
                ddelta_drhoi(i) = (gl%rhoredmix - Dens * props_phase(a)%drhor_drhoi_it(i)) / gl%rhoredmix**2
                dtau_drhoi(i) = props_phase(a)%dTr_drhoi_it(i) / Temp
                drhorTr_dxm(i) = gl%tredmix * props_phase(a)%drhor_dxi_it(i) + gl%rhoredmix * props_phase(a)%dTr_dxi_it(i)
                do m = 1, gl%ncomp
                    dxm_drhoi(m,i) = (Dens * kron_delta(m,i) - Dens * gl%molfractions(m)) / Dens**2
                end do
            end do

            !Calculate the derivative of the residual Helmholtz energy density with respect to delta
            !dPsir/ddelta
            dPsir_ddelta = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1))

            !Calculate the derivative of the residual Helmholtz energy density with respect to tau
            !dPsir/dtau
            dPsir_dtau = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1))

            !Calculate the derivative of the residual Helmholtz energy density with respect to xm
            !dPsir/dxm
            do m = 1, gl%ncomp
                dPsir_dxm(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * (props_phase(a)%alphar_it(1) * drhorTr_dxm(m) + gl%rhoredmix * gl%tredmix * props_phase(a)%dalphar_dxi_it(m))
            end do

            !Calculate the first derivatives of Psi_r with respect to the molar concentrations (rho_i)
            !It is (dPsi_r / drho_i = mue_r_i)
            Do i = 1, gl%ncomp
                dPsir_drhoi(i) = dPsir_ddelta * ddelta_drhoi(i) + dPsir_dtau * dtau_drhoi(i)
                do m = 1, gl%ncomp
                    dPsir_drhoi(i) = dPsir_drhoi(i) + dPsir_dxm(m) * dxm_drhoi(m,i)
                end do
            End do
            !-------------------------------------------------------------------------------------------------------------




            !Calculate the second derivatives of Psi with respect to the molar concentration rho_i and rho_j
            !-------------------------------------------------------------------------------------------------------------
            !Calculate help variables
            ! d2delta/drhoidrhoj, d2tau/drhoidrhoj, d2xm/drhoidrhoj, d(rhor*Tr)/dxmdxn
            do i = 1, gl%ncomp
                do j = i, gl%ncomp
                    d2delta_drhoidrhoj(i,j) = (gl%rhoredmix**2 * (-props_phase(a)%drhor_drhoi_it(j) - props_phase(a)%drhor_drhoi_it(i) - Dens * props_phase(a)%d2rhor_drhoidrhoj_it(i,j)) + Dens * props_phase(a)%drhor_drhoi_it(i) * 2.D0 * gl%rhoredmix * props_phase(a)%drhor_drhoi_it(j)) / gl%rhoredmix**4
                    d2tau_drhoidrhoj(i,j) = props_phase(a)%d2Tr_drhoidrhoj_it(i,j) / Temp
                    d2rhorTr_dxmdxn(i,j) = gl%tredmix * props_phase(a)%d2rhor_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(i) * props_phase(a)%drhor_dxi_it(j) + gl%rhoredmix * props_phase(a)%d2Tr_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(j) * props_phase(a)%drhor_dxi_it(i)
                    do m = 1, gl%ncomp
                        d2xm_drhoidrhoj(m,i,j) = (Dens**2 * (kron_delta(i,m) - kron_delta(j,m)) - (Dens * kron_delta(i,m) - gl%molfractions(m) * Dens) * 2.D0 * Dens) / Dens**4
                        if (i /= j) then
                            d2xm_drhoidrhoj(m,j,i) = d2xm_drhoidrhoj(m,i,j)
                        end if
                    end do
                    if (i /= j) then
                        d2delta_drhoidrhoj(j,i) = d2delta_drhoidrhoj(i,j)
                        d2tau_drhoidrhoj(j,i) = d2tau_drhoidrhoj(i,j)
                        d2rhorTr_dxmdxn(j,i) = d2rhorTr_dxmdxn(i,j)
                    end if
                end do
            end do

            !Calculate the second derivative of the residual Helmholtz energy density with respect to delta and delta
            !d2Psir/ddelta^2
            d2Psir_ddelta2 = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(3) + 2.D0 * props_phase(a)%alphar_it(2)) / props_phase(a)%delta_it

            !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and delta
            !d2Psir/dtauddelta
            d2Psir_dtauddelta = gl%rhoredmix * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1) - props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(6))

            !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and tau
            !d2Psir/dtau^2
            d2Psir_dtau2 = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**3 * (props_phase(a)%alphar_it(5) - 2.D0 * props_phase(a)%alphar_it(4) + 2.D0 * props_phase(a)%alphar_it(1))

            do m = 1, gl%ncomp
                !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and delta
                !d2Psir/dxmddelta
                d2Psir_dxmddel(m) = R_mix / props_phase(a)%tau_it * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxiddel_it(m) + props_phase(a)%dalphar_dxi_it(m)))

                !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and tau
                !d2Psir/dxmdtau
                d2Psir_dxmdtau(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it**2 * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxidtau_it(m) - props_phase(a)%dalphar_dxi_it(m)))

                do n = m, gl%ncomp
                    !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and xn
                    !d2Psir/dxmdxn
                    d2Psir_dxmdxn(m,n) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * &
                        & (props_phase(a)%alphar_it(1) * d2rhorTr_dxmdxn(m,n) + props_phase(a)%dalphar_dxi_it(n) * drhorTr_dxm(m) + &
                        &  gl%rhoredmix * gl%tredmix * props_phase(a)%d2alphar_dxidxj_it(m,n) + props_phase(a)%dalphar_dxi_it(m) * drhorTr_dxm(n))
                    if (n /= m) then
                        d2Psir_dxmdxn(n,m) = d2Psir_dxmdxn(m,n)
                    end if
                end do
            end do


            !Calculate the second derivatives of the residual Helmholtz energy density with respect to delta, tau, xm and rho_j
            !d(dPsir/ddelta)/drho_j, d(dPsir/dtau)/drho_j, d(dPsir/dxm)/drho_j
            do j = 1, gl%ncomp
                ddPsirddelta_drhoj(j) = d2Psir_ddelta2 * ddelta_drhoi(j) + d2Psir_dtauddelta * dtau_drhoi(j)
                ddPsirdtau_drhoj(j) = d2Psir_dtauddelta * ddelta_drhoi(j) + d2Psir_dtau2 * dtau_drhoi(j)
                do m = 1, gl%ncomp
                    ddPsirddelta_drhoj(j) = ddPsirddelta_drhoj(j) + d2Psir_dxmddel(m) * dxm_drhoi(m,j)
                    ddPsirdtau_drhoj(j) = ddPsirdtau_drhoj(j) + d2Psir_dxmdtau(m) * dxm_drhoi(m,j)
                    ddPsirdxm_drhoj(m,j) = d2Psir_dxmddel(m) * ddelta_drhoi(j) + d2Psir_dxmdtau(m) * dtau_drhoi(j)
                    do n = 1, gl%ncomp
                        ddPsirdxm_drhoj(m,j) = ddPsirdxm_drhoj(m,j) + d2Psir_dxmdxn(m,n) * dxm_drhoi(n,j)
                    end do
                end do
            end do


            !Calculate the second derivatives of Psi_r with respect to the molar concentrations rho_i and rho_j
            !(d2Psi_r / (drho_i drho_j))
            Do i = 1, gl%ncomp
                do j = i, gl%ncomp
                    d2Psir_drhoidrhoj(i,j) = dPsir_ddelta * d2delta_drhoidrhoj(i,j) + ddPsirddelta_drhoj(j) * ddelta_drhoi(i) + dPsir_dtau * d2tau_drhoidrhoj(i,j) + ddPsirdtau_drhoj(j) * dtau_drhoi(i)
                    do m = 1, gl%ncomp
                        d2Psir_drhoidrhoj(i,j) = d2Psir_drhoidrhoj(i,j) + dPsir_dxm(m) * d2xm_drhoidrhoj(m,i,j) + ddPsirdxm_drhoj(m,j) * dxm_drhoi(m,i)
                    end do
                    if (i /= j) then
                        d2Psir_drhoidrhoj(j,i) = d2Psir_drhoidrhoj(i,j)
                    end if
                end do
            End do
            !-------------------------------------------------------------------------------------------------------------

            !Calculate the Hessian matrix
            do i = 1, gl%ncomp
                do j = i, gl%ncomp
                    if (i == j) then
                        Hessian(i,j,a) = d2Psir_drhoidrhoj(i,j) + R_mix * Temp / gl%molfractions(i) / Dens
                    else
                        Hessian(i,j,a) = d2Psir_drhoidrhoj(i,j)
                        Hessian(j,i,a) = Hessian(i,j,a)
                    end if
                end do
            end do


            !This derivative is (at the moment) not needed in this routine
            !!Calculate the second derivatives of Psi with respect to the molar concentration rho_i and T
            !!-------------------------------------------------------------------------------------------------------------
            !! Calculate help variables
            !! dtau/dT, d2tau/drhoidT (Variables not explicitly given in paper by Bell and Deiters?)
            !dtau_dT = - gl%tredmix / Temp**2
            !do i = 1, gl%ncomp
            !    d2tau_drhoi_dT(i) = -props_phase(a)%dTr_drhoi_it(i) / Temp**2
            !end do
            !
            !!Calculate the second derivatives of Psi_r with respect to T and molar concentration rhoi
            !!(d2Psi_r / (drho_i dT))
            !do i = 1, gl%ncomp
            !    d2Psir_drhoidT(i) = dPsir_dtau * d2tau_drhoi_dT(i) + d2Psir_dtauddelta * ddelta_drhoi(i) * dtau_dT + &
            !                      & d2Psir_dtau2 * dtau_drhoi(i) * dtau_dT
            !    do m = 1, gl%ncomp
            !        d2Psir_drhoidT(i) = d2Psir_drhoidT(i) + d2Psir_dxmdtau(m) * dxm_drhoi(m,i) * dtau_dT
            !    end do
            !end do
            !!-------------------------------------------------------------------------------------------------------------

            !Calculate the pressure of the phase

            press_phase(a) = dens_phase(a) * R_mix * Temp - Psir
            do i = 1, gl%ncomp
                press_phase(a) = press_phase(a) + dens_i_phase(i,a) * dPsir_drhoi(i)
            end do
            !Test
            !press_phase(a) = dens_phase(a) * R_mix * Temp * (1.D0 + props_phase(a)%alphar_it(2))

        end do


        !Calculate the derivative of the molar concentrations with respect to pressure along the bubble line
        !Eq. (19) in the article of Bell and Deiters
        !-------------------------------------------------------------------------------------------------------------
        !Set up right side of the system of equations
        do i = 1, gl%ncomp
            vectorb(i) = 1.D0
        end do
        !Set up Matrix A (FOR 2 COMPONENTS ONLY)
        MatrixA(1,1) = Hessian(1,1,2) * dens_i_phase(1,1) + Hessian(1,2,2) * dens_i_phase(2,1)
        MatrixA(1,2) = Hessian(2,1,2) * dens_i_phase(1,1) + Hessian(2,2,2) * dens_i_phase(2,1)
        MatrixA(2,1) = Hessian(1,1,2) * dens_i_phase(1,2) + Hessian(1,2,2) * dens_i_phase(2,2)
        MatrixA(2,2) = Hessian(2,1,2) * dens_i_phase(1,2) + Hessian(2,2,2) * dens_i_phase(2,2)
        rankA = gl%ncomp
        det_A = 0.D0
        call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errval)
        drhoi_dp(1,2) = vectorx(1)
        drhoi_dp(2,2) = vectorx(2)
        !-------------------------------------------------------------------------------------------------------------


        !Calculate the derivative of the molar concentrations with respect to pressure along the dew line
        !Eq. (20) in the article of Bell and Deiters
        !-------------------------------------------------------------------------------------------------------------
        !Set up right side of the system of equations (FOR 2 COMPONENTS ONLY)
        vectorb(1) = Hessian(1,1,2) * drhoi_dp(1,2) + Hessian(1,2,2) * drhoi_dp(2,2)
        vectorb(2) = Hessian(2,1,2) * drhoi_dp(1,2) + Hessian(2,2,2) * drhoi_dp(2,2)
        !Set up Matrix A (Matrix A is simply the Hessian of the vapor phase in this case)
        MatrixA(1,1) = Hessian(1,1,1)
        MatrixA(1,2) = Hessian(1,2,1)
        MatrixA(2,1) = Hessian(2,1,1)
        MatrixA(2,2) = Hessian(2,2,1)
        rankA = gl%ncomp
        det_A = 0.D0
        call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errval)
        drhoi_dp(1,1) = vectorx(1)
        drhoi_dp(2,1) = vectorx(2)
        !-------------------------------------------------------------------------------------------------------------


        !Write the needed derivatives in the function output
        if ((spec_deriv == 1) .or. (spec_deriv == 2)) then
            !Derivatives with respect to a molar concentration of the gas phase
            drhoi_drhoj(1) = drhoi_dp(1,1) / drhoi_dp(spec_deriv,1)      !drho1_vap/drhoi
            drhoi_drhoj(2) = drhoi_dp(2,1) / drhoi_dp(spec_deriv,1)      !drho2_vap/drhoi
            drhoi_drhoj(3) = drhoi_dp(1,2) / drhoi_dp(spec_deriv,1)      !drho1_liq/drhoi
            drhoi_drhoj(4) = drhoi_dp(2,2) / drhoi_dp(spec_deriv,1)      !drho2_liq/drhoi
        else if ((spec_deriv == 3) .or. (spec_deriv == 4)) then
            !Derivatives with respect to a molar concentration of the liquid phase
            drhoi_drhoj(1) = drhoi_dp(1,1) / drhoi_dp(spec_deriv-2,2)      !drho1_vap/drhoi
            drhoi_drhoj(2) = drhoi_dp(2,1) / drhoi_dp(spec_deriv-2,2)      !drho2_vap/drhoi
            drhoi_drhoj(3) = drhoi_dp(1,2) / drhoi_dp(spec_deriv-2,2)      !drho1_liq/drhoi
            drhoi_drhoj(4) = drhoi_dp(2,2) / drhoi_dp(spec_deriv-2,2)      !drho2_liq/drhoi
        else if (spec_deriv == 5) then
            !Derivatives with respect topressire
            drhoi_drhoj(1) = drhoi_dp(1,1)      !drho1_vap/dp
            drhoi_drhoj(2) = drhoi_dp(2,1)      !drho2_vap/dp
            drhoi_drhoj(3) = drhoi_dp(1,2)      !drho1_liq/dp
            drhoi_drhoj(4) = drhoi_dp(2,2)      !drho2_liq/dp
        end if

        !Write pressure to output file (The pressures are most likely not exactly equal, arbitrarily take pressure of vapor phase)
        out_var(1) = press_phase(1)

    else    !Do a "polishing" step, as explained in the supplementary material, sections 2.1 and section 2.2 of the article of Bell and Deiters (2018)

        !Newton-Raphson iteration method for solving the phase equilibrium conditions formulated in the Helmholtz energy density Psi
        eps_res = 1.D-8     !Convergence criterion, arbitrary value (this is maybe too small for differences in pressure in Pa. Thus, mostly the second convergence criterion will be used)
        eps_deltax = 1.D-10 !Second convergence criterion, arbitrary value

        !-----------------------------------------------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------------------------------------
        !Do a specified number of iterations
        Do k = 1, 30

            !Write the y-function in the corresponding variables
            dens_i_phase(1,1) = y_vector(1)
            dens_i_phase(2,1) = y_vector(2)
            dens_i_phase(1,2) = y_vector(3)
            dens_i_phase(2,2) = y_vector(4)

            !Density of the vapor phase (phase 1)
            dens_phase(1) = y_vector(1) + y_vector(2)
            !Calculate the composition of the vapor phase (phase 1)
            x_phase(1,1) = y_vector(1) / dens_phase(1)
            x_phase(2,1) = y_vector(2) / dens_phase(1)

            !Density of the liquid phase (phase 2)
            dens_phase(2) = y_vector(3) + y_vector(4)
            !Calculate the composition of the vapor phase (phase 1)
            x_phase(1,2) = y_vector(3) / dens_phase(2)
            x_phase(2,2) = y_vector(4) / dens_phase(2)

            !Calculate the needed derivatives for setting up the system of equations and the Jacobian matrix
            do a = 1, 2

                !Set density and composition of the phase and get all needed derivatives
                Dens = dens_phase(a)
                gl%molfractions = 0.D0
                gl%molfractions(1:2) = x_phase(:,a)
                call reduced_parameters_calc(gl,Temp)

                !Get the mixture gas constant
                call R_mix_calc(gl, R_mix)
                R_mix_phase(a) = R_mix      !Save R for the phase in case different R are used with the different EOS

                call Derivs_isoch_therm(gl, props_phase(a), Temp, Dens, errval)
                if (errval /= 0) then
                    return
                end if

                !Calculate the residual Helmholtz energy density Psi^r = a^r * rho
                Psir = Dens * R_mix * Temp * props_phase(a)%alphar_it(1)

                !Calculate the first derivative of Psi with respect to the molar concentration (i.e. the chemical potentials)
                !-------------------------------------------------------------------------------------------------------------
                !Calculate help variables
                ! ddelta/drhoi, dtau/drhoi, dxm/drhoi, d(rhor*Tr)/dxm
                do i = 1, gl%ncomp
                    ddelta_drhoi(i) = (gl%rhoredmix - Dens * props_phase(a)%drhor_drhoi_it(i)) / gl%rhoredmix**2
                    dtau_drhoi(i) = props_phase(a)%dTr_drhoi_it(i) / Temp
                    drhorTr_dxm(i) = gl%tredmix * props_phase(a)%drhor_dxi_it(i) + gl%rhoredmix * props_phase(a)%dTr_dxi_it(i)
                    do m = 1, gl%ncomp
                        dxm_drhoi(m,i) = (Dens * kron_delta(m,i) - Dens * gl%molfractions(m)) / Dens**2
                    end do
                end do

                !Calculate the derivative of the residual Helmholtz energy density with respect to delta
                !dPsir/ddelta
                dPsir_ddelta = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1))

                !Calculate the derivative of the residual Helmholtz energy density with respect to tau
                !dPsir/dtau
                dPsir_dtau = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1))

                !Calculate the derivative of the residual Helmholtz energy density with respect to xm
                !dPsir/dxm
                do m = 1, gl%ncomp
                    dPsir_dxm(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * (props_phase(a)%alphar_it(1) * drhorTr_dxm(m) + gl%rhoredmix * gl%tredmix * props_phase(a)%dalphar_dxi_it(m))
                end do

                !Calculate the first derivatives of Psi_r with respect to the molar concentrations (rho_i)
                !It is (dPsi_r / drho_i = mue_r_i)
                Do i = 1, gl%ncomp
                    dPsir_drhoi(i) = dPsir_ddelta * ddelta_drhoi(i) + dPsir_dtau * dtau_drhoi(i)
                    do m = 1, gl%ncomp
                        dPsir_drhoi(i) = dPsir_drhoi(i) + dPsir_dxm(m) * dxm_drhoi(m,i)
                    end do
                End do
                !-------------------------------------------------------------------------------------------------------------


                !Calculate the second derivatives of Psi with respect to the molar concentration rho_i and rho_j
                !-------------------------------------------------------------------------------------------------------------
                !Calculate help variables
                ! d2delta/drhoidrhoj, d2tau/drhoidrhoj, d2xm/drhoidrhoj, d(rhor*Tr)/dxmdxn
                do i = 1, gl%ncomp
                    do j = i, gl%ncomp
                        d2delta_drhoidrhoj(i,j) = (gl%rhoredmix**2 * (-props_phase(a)%drhor_drhoi_it(j) - props_phase(a)%drhor_drhoi_it(i) - Dens * props_phase(a)%d2rhor_drhoidrhoj_it(i,j)) + Dens * props_phase(a)%drhor_drhoi_it(i) * 2.D0 * gl%rhoredmix * props_phase(a)%drhor_drhoi_it(j)) / gl%rhoredmix**4
                        d2tau_drhoidrhoj(i,j) = props_phase(a)%d2Tr_drhoidrhoj_it(i,j) / Temp
                        d2rhorTr_dxmdxn(i,j) = gl%tredmix * props_phase(a)%d2rhor_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(i) * props_phase(a)%drhor_dxi_it(j) + gl%rhoredmix * props_phase(a)%d2Tr_dxidxj_it(i,j) + props_phase(a)%dTr_dxi_it(j) * props_phase(a)%drhor_dxi_it(i)
                        do m = 1, gl%ncomp
                            d2xm_drhoidrhoj(m,i,j) = (Dens**2 * (kron_delta(i,m) - kron_delta(j,m)) - (Dens * kron_delta(i,m) - gl%molfractions(m) * Dens) * 2.D0 * Dens) / Dens**4
                            if (i /= j) then
                                d2xm_drhoidrhoj(m,j,i) = d2xm_drhoidrhoj(m,i,j)
                            end if
                        end do
                        if (i /= j) then
                            d2delta_drhoidrhoj(j,i) = d2delta_drhoidrhoj(i,j)
                            d2tau_drhoidrhoj(j,i) = d2tau_drhoidrhoj(i,j)
                            d2rhorTr_dxmdxn(j,i) = d2rhorTr_dxmdxn(i,j)
                        end if
                    end do
                end do

                !Calculate the second derivative of the residual Helmholtz energy density with respect to delta and delta
                !d2Psir/ddelta^2
                d2Psir_ddelta2 = gl%rhoredmix * R_mix * Temp * (props_phase(a)%alphar_it(3) + 2.D0 * props_phase(a)%alphar_it(2)) / props_phase(a)%delta_it

                !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and delta
                !d2Psir/dtauddelta
                d2Psir_dtauddelta = gl%rhoredmix * R_mix * gl%tredmix / props_phase(a)%tau_it**2 * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1) - props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(6))

                !Calculate the second derivative of the residual Helmholtz energy density with respect to tau and tau
                !d2Psir/dtau^2
                d2Psir_dtau2 = gl%rhoredmix * props_phase(a)%delta_it * R_mix * gl%tredmix / props_phase(a)%tau_it**3 * (props_phase(a)%alphar_it(5) - 2.D0 * props_phase(a)%alphar_it(4) + 2.D0 * props_phase(a)%alphar_it(1))

                do m = 1, gl%ncomp
                    !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and delta
                    !d2Psir/dxmddelta
                    d2Psir_dxmddel(m) = R_mix / props_phase(a)%tau_it * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(2) + props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxiddel_it(m) + props_phase(a)%dalphar_dxi_it(m)))

                    !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and tau
                    !d2Psir/dxmdtau
                    d2Psir_dxmdtau(m) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it**2 * (drhorTr_dxm(m) * (props_phase(a)%alphar_it(4) - props_phase(a)%alphar_it(1)) + gl%rhoredmix * gl%tredmix * (props_phase(a)%d2alphar_dxidtau_it(m) - props_phase(a)%dalphar_dxi_it(m)))

                    do n = m, gl%ncomp
                        !Calculate the second derivative of the residual Helmholtz energy density with respect to xm and xn
                        !d2Psir/dxmdxn
                        d2Psir_dxmdxn(m,n) = props_phase(a)%delta_it * R_mix / props_phase(a)%tau_it * &
                            & (props_phase(a)%alphar_it(1) * d2rhorTr_dxmdxn(m,n) + props_phase(a)%dalphar_dxi_it(n) * drhorTr_dxm(m) + &
                            &  gl%rhoredmix * gl%tredmix * props_phase(a)%d2alphar_dxidxj_it(m,n) + props_phase(a)%dalphar_dxi_it(m) * drhorTr_dxm(n))
                        if (n /= m) then
                            d2Psir_dxmdxn(n,m) = d2Psir_dxmdxn(m,n)
                        end if
                    end do
                end do


                !Calculate the second derivatives of the residual Helmholtz energy density with respect to delta, tau, xm and rho_j
                !d(dPsir/ddelta)/drho_j, d(dPsir/dtau)/drho_j, d(dPsir/dxm)/drho_j
                do j = 1, gl%ncomp
                    ddPsirddelta_drhoj(j) = d2Psir_ddelta2 * ddelta_drhoi(j) + d2Psir_dtauddelta * dtau_drhoi(j)
                    ddPsirdtau_drhoj(j) = d2Psir_dtauddelta * ddelta_drhoi(j) + d2Psir_dtau2 * dtau_drhoi(j)
                    do m = 1, gl%ncomp
                        ddPsirddelta_drhoj(j) = ddPsirddelta_drhoj(j) + d2Psir_dxmddel(m) * dxm_drhoi(m,j)
                        ddPsirdtau_drhoj(j) = ddPsirdtau_drhoj(j) + d2Psir_dxmdtau(m) * dxm_drhoi(m,j)
                        ddPsirdxm_drhoj(m,j) = d2Psir_dxmddel(m) * ddelta_drhoi(j) + d2Psir_dxmdtau(m) * dtau_drhoi(j)
                        do n = 1, gl%ncomp
                            ddPsirdxm_drhoj(m,j) = ddPsirdxm_drhoj(m,j) + d2Psir_dxmdxn(m,n) * dxm_drhoi(n,j)
                        end do
                    end do
                end do


                !Calculate the second derivatives of Psi_r with respect to the molar concentrations rho_i and rho_j
                !(d2Psi_r / (drho_i drho_j))
                Do i = 1, gl%ncomp
                    do j = i, gl%ncomp
                        d2Psir_drhoidrhoj(i,j) = dPsir_ddelta * d2delta_drhoidrhoj(i,j) + ddPsirddelta_drhoj(j) * ddelta_drhoi(i) + dPsir_dtau * d2tau_drhoidrhoj(i,j) + ddPsirdtau_drhoj(j) * dtau_drhoi(i)
                        do m = 1, gl%ncomp
                            d2Psir_drhoidrhoj(i,j) = d2Psir_drhoidrhoj(i,j) + dPsir_dxm(m) * d2xm_drhoidrhoj(m,i,j) + ddPsirdxm_drhoj(m,j) * dxm_drhoi(m,i)
                        end do
                        if (i /= j) then
                            d2Psir_drhoidrhoj(j,i) = d2Psir_drhoidrhoj(i,j)
                        end if
                    end do
                End do
                !-------------------------------------------------------------------------------------------------------------

                !Save variables for the phases for setting up the system of equations
                dPsir_drhoi_phase(:,a) = dPsir_drhoi(:)
                d2Psir_drhoidrhoj_phase(:,:,a) = d2Psir_drhoidrhoj(:,:)

                !Calculate the pressure of the phase
                press_phase(a) = dens_phase(a) * R_mix * Temp - Psir
                do i = 1, gl%ncomp
                    press_phase(a) = press_phase(a) + dens_i_phase(i,a) * dPsir_drhoi(i)
                end do
                !Test
                !press_phase(a) = dens_phase(a) * R_mix * Temp * (1.D0 + props_phase(a)%alphar_it(2))

                !Calculate the derivative of the pressure of the phase with respect to the molar concentracions rho_i
                do j = 1, gl%ncomp
                    dp_drhoj(j,a) = R_mix * Temp
                    do i = 1, gl%ncomp
                        dp_drhoj(j,a) = dp_drhoj(j,a) + dens_i_phase(i,a) * d2Psir_drhoidrhoj(i,j)
                    end do
                end do

            end do

            !Write pressure to output file (arbitrarily take pressure of vapor phase)
            out_var(1) = press_phase(1)

            !Set up the system of equations according to the specified marching variable:
            if (spec_deriv < 5) then
                !Sys_of_Eqs(1) = ChemPot_1_liq - Chempot1_vap = 0
                !Sys_of_Eqs(2) = ChemPot_2_liq - Chempot2_vap = 0
                !Sys_of_Eqs(3) = p_liq - p_vap = 0
                !Sys_of_Eqs(4) = rho(specified component, specified phase) - rho_1(specified)
                !******************************
                Sys_of_Eqs(1) = dPsir_drhoi_phase(1,2) + R_mix_phase(2) * Temp * dlog(dens_i_phase(1,2)) - dPsir_drhoi_phase(1,1) - R_mix_phase(1) * Temp * dlog(dens_i_phase(1,1))
                Sys_of_Eqs(2) = dPsir_drhoi_phase(2,2) + R_mix_phase(2) * Temp * dlog(dens_i_phase(2,2)) - dPsir_drhoi_phase(2,1) - R_mix_phase(1) * Temp * dlog(dens_i_phase(2,1))
                Sys_of_Eqs(3) = press_phase(2) - press_phase(1)
                Select case (spec_deriv)
                    !NOTE: THIS EQUATION IS IN PRINCIPLE NOT NECESSARY, SHOULD ALWAYS BE 0. However, this simplifies the implementation of the algorithm
                Case (1)
                    !Molar concentration of first component in vapor is the specified variable
                    Sys_of_Eqs(4) = dens_i_phase(1,1) - var_march
                Case (2)
                    !Molar concentration of second component in vapor is the specified variable
                    Sys_of_Eqs(4) = dens_i_phase(2,1) - var_march
                Case (3)
                    !Molar concentration of first component in liquid is the specified variable
                    Sys_of_Eqs(4) = dens_i_phase(1,2) - var_march
                Case (4)
                    !Molar concentration of second component in liquid is the specified variable
                    Sys_of_Eqs(4) = dens_i_phase(2,2) - var_march
                end select
                !******************************
            else if (spec_deriv == 5) then
                !Sys_of_Eqs(1) = ChemPot_1_liq - Chempot1_vap = 0
                !Sys_of_Eqs(2) = ChemPot_2_liq - Chempot2_vap = 0
                !Sys_of_Eqs(3) = p_liq - p_spec = 0
                !Sys_of_Eqs(4) = p_vap - p_spec = 0
                !******************************
                Sys_of_Eqs(1) = dPsir_drhoi_phase(1,2) + R_mix_phase(2) * Temp * dlog(dens_i_phase(1,2)) - dPsir_drhoi_phase(1,1) - R_mix_phase(1) * Temp * dlog(dens_i_phase(1,1))
                Sys_of_Eqs(2) = dPsir_drhoi_phase(2,2) + R_mix_phase(2) * Temp * dlog(dens_i_phase(2,2)) - dPsir_drhoi_phase(2,1) - R_mix_phase(1) * Temp * dlog(dens_i_phase(2,1))
                Sys_of_Eqs(3) = press_phase(2) - var_march
                Sys_of_Eqs(4) = press_phase(1) - var_march
                !******************************
            end if

            !Check for convergence
            if (maxval(dabs(Sys_of_Eqs)) < eps_res) then
                drhoi_drhoj = 0.D0      !Derivatives are not calculated if a post-processing step is taken
                exit
            end if

            !Set up the Jacobian matrix
            !(Note: order of variables is different here compared to the article of Bell and Deiters (2018): (HERE: rho1_vap, rho2_vap, rho1_liq, rho2_liq)
            if (spec_deriv < 5) then
                !******************************
                Jac_Mat(1,1) = -d2Psir_drhoidrhoj_phase(1,1,1) - R_mix_phase(1) * Temp / dens_i_phase(1,1)
                Jac_Mat(1,2) = -d2Psir_drhoidrhoj_phase(1,2,1)
                Jac_Mat(1,3) = d2Psir_drhoidrhoj_phase(1,1,2) + R_mix_phase(2) * Temp / dens_i_phase(1,2)
                Jac_Mat(1,4) = d2Psir_drhoidrhoj_phase(1,2,2)

                Jac_Mat(2,1) = -d2Psir_drhoidrhoj_phase(2,1,1)
                Jac_Mat(2,2) = -d2Psir_drhoidrhoj_phase(2,2,1) - R_mix_phase(1) * Temp / dens_i_phase(2,1)
                Jac_Mat(2,3) = d2Psir_drhoidrhoj_phase(2,1,2)
                Jac_Mat(2,4) = d2Psir_drhoidrhoj_phase(2,2,2) + R_mix_phase(2) * Temp / dens_i_phase(2,2)

                Jac_Mat(3,1) = -dp_drhoj(1,1)
                Jac_Mat(3,2) = -dp_drhoj(2,1)
                Jac_Mat(3,3) = dp_drhoj(1,2)
                Jac_Mat(3,4) = dp_drhoj(2,2)

                Jac_Mat(4,1) = 0.D0
                Jac_Mat(4,2) = 0.D0
                Jac_Mat(4,3) = 0.D0
                Jac_Mat(4,4) = 0.D0
                Jac_Mat(4,spec_deriv) = 1.D0
                !******************************
            else if (spec_deriv == 5) then
                !******************************
                Jac_Mat(1,1) = -d2Psir_drhoidrhoj_phase(1,1,1) - R_mix_phase(1) * Temp / dens_i_phase(1,1)
                Jac_Mat(1,2) = -d2Psir_drhoidrhoj_phase(1,2,1)
                Jac_Mat(1,3) = d2Psir_drhoidrhoj_phase(1,1,2) + R_mix_phase(2) * Temp / dens_i_phase(1,2)
                Jac_Mat(1,4) = d2Psir_drhoidrhoj_phase(1,2,2)

                Jac_Mat(2,1) = -d2Psir_drhoidrhoj_phase(2,1,1)
                Jac_Mat(2,2) = -d2Psir_drhoidrhoj_phase(2,2,1) - R_mix_phase(1) * Temp / dens_i_phase(2,1)
                Jac_Mat(2,3) = d2Psir_drhoidrhoj_phase(2,1,2)
                Jac_Mat(2,4) = d2Psir_drhoidrhoj_phase(2,2,2) + R_mix_phase(2) * Temp / dens_i_phase(2,2)

                Jac_Mat(3,1) = 0.D0
                Jac_Mat(3,2) = 0.D0
                Jac_Mat(3,3) = dp_drhoj(1,2)
                Jac_Mat(3,4) = dp_drhoj(2,2)

                Jac_Mat(4,1) = dp_drhoj(1,1)
                Jac_Mat(4,2) = dp_drhoj(2,1)
                Jac_Mat(4,3) = 0.D0
                Jac_Mat(4,4) = 0.D0
                !******************************
            end if

            !Solve the system of equations
            vectorb(1:4) = -Sys_of_Eqs(:)
            !Set up Matrix A (Matrix A is simply the Hessian of the vapor phase in this case)
            MatrixA(1:4,1:4) = Jac_Mat(:,:)
            rankA = 4
            det_A = 0.D0
            call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errval)
            Delta_X(:) = vectorx(1:4)

            !Check for convergence
            if (maxval(dabs(Delta_X)) < eps_res) then
                drhoi_drhoj = 0.D0      !Derivatives are not calculated if a post-processing step is taken
                exit
            end if

            !Execute the step
            y_vector = y_vector + Delta_X

        end do

        if (k > 30) errval = -2222
        drhoi_drhoj = 0.D0      !Derivatives are not calculated if a post-processing step is taken

        !-----------------------------------------------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------------------------------------

    end if

    end function drhoi_drhoj




    end module ptx_diag_module



    !This module contains test functions for the Runge-Kutta-Fehlberg Method
    !The module with the "contains"-statement is necessary in order to create
    !functions that are able to return a vector of values
    !Andreas Jäger, November 2018
    module module_RKF_functions
    use module_all_types
    use module_isoch_therm
    use pc_saft_ARX_derivs_module
    use ptx_diag_module

    implicit none
    contains

    !!--------------------------------------------------------------------------------
    !Test function for the Runge-Kutta-Fehlberg method.
    function RKF_test_function1(gl, nr_of_eqs, t_march, y_vector, ind_var, out_var, post_proc, errorflag)
    !!--------------------------------------------------------------------------------
    !! THIS SUBROUTINE IS A TEST FUNCTION FOR THE RUNGE-KUTTA-FEHLBERG METHOD
    !! FOR SOLVING AN ORDINARY DIFFERENTIAL EQUATION OF FIRST ORDER OF THE FORM
    !! dy/dt = f(t,y)
    !! HENCE, THIS FUNCTION TAKES THE ACTUAL VALUES Y AS INPUT AND CALCULATES THE
    !! DERIVATIVES OF Y WITH RESPECT TO t AT THE SPECIFIED POSITION t
    !! THE TEST FUNCTION IS:
    !! dy/dt = y
    !! SO OBVIOUSLY THE SOLUTION IS
    !! y = e^t + C
    !!
    !! INPUT:   nr_of_eqs           - The number of equations the ODE consists of
    !!          t_march             - value of the "marching" variable t
    !!          y_vector            - vector of size 1 containing the values of y(t)
    !!          ind_var             - Other independent variables that might be needed to evaluate the derivatives of y with respect to t (here: none)
    !!          out_var             - Other (additional) output variables (here: none)
    !!          post_proc           - True: Run function as post-processing function, false: run function in normal mode. Normal mode should solve y' = f(t,y) (here: no post-processing necessary and thus not implemented)
    !! OUTPUT:  RKF_test_function1  - Derivative of y with respect to t (here: y' = y)
    !!          errorflag           - indicates if an error occured in the calculations
    !!--------------------------------------------------------------------------------
    !Andreas Jäger, November 2018



    implicit none

    type(type_gl) :: gl

    !Input variables
    integer :: nr_of_eqs
    double precision, dimension(nr_of_eqs) :: y_vector
    double precision :: t_march
    double precision, dimension(10) :: ind_var              !Other independent variables that might be needed to evaluate the derivatives of y with respect to t
    double precision, dimension(10) :: out_var              !Other (optional) output variables
    logical :: post_proc
    !Output variables
    integer :: errorflag
    double precision, dimension(nr_of_eqs) :: RKF_test_function1   !(the function itself)

    errorflag = 0

    !y'(t) = y(t)
    RKF_test_function1(:) = y_vector(:)

    end function RKF_test_function1






    end module


