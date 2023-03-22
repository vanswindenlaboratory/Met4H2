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

    ! module for file phasedet_sol_pure.f90
    module phasedet_sol_pure_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use rhomix_pt_module
    use flash_pure_module
    use waterice_module
    use dryice_module
    use seawater_module
    use phasedet_pure_module
    use ancillary_equations_mix_module
    use ptflash_solids_mix_module

    contains






    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !Andreas Februar 2014
    !PhaseDet pure routines with solid formation
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !************************************************************************************
    subroutine PhaseDet_pure_tp_sol (gl,press, Temp, d, rho, phasetype, phasefrac, nrofphases, errval)
    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !
    ! OUTPUT:
    !   rho         - density [mol/m³] of all phases
    !   d           - density [mol/m³] overall

    !   The following two variables indicate the phase type. For a pure substance,
    !   a maximum of 3 phases can be in equilibrium
    !   phasetype   - phase identifier
    !   phasefrac   - Molar phase fraction

    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************








    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    integer :: nrofphases
    integer :: errval
    integer :: fluidnr

    double precision:: pmelt, psub, psat, temptrip, presstrip
    double precision:: d2pdd2, gice, gsea, gvap

    double precision, dimension(30):: x_known, x_solid, x_hyd, x_fluid, x_fluid1, x_fluid2, x_sol
    double precision :: rhofluid_est, beta_loc, rhovap_est, rholiq_est
    integer :: iflash, iPhase, iter

    double precision :: ptrip, ttrip, dw, dsat_w, rho_l_melt_w, tfreeze   !seawater

    double precision, dimension(3) :: rho_sol

    !double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_sol, x_hyd

    phasefrac = 0.D0
    phasetype = 0
    rho = 0.D0

    fluidnr = 1

    x_known = 0.D0
    x_solid = 0.D0
    x_hyd = 0.D0
    x_fluid = 0.D0
    x_known(1) = 1.D0
    rhofluid_est = 0.D0
    tfreeze = 256.164d0

    !In case two melting pressures at given temperature exist for solid water, the higher one is saved in the variable pmelt_high
    gl%pmelt_high = 0.D0

    pmelt = 0.D0
    psub = 0.D0


    if (gl%solidtype(1) == 1) then      !WATER

        gl%solidtype(1) = 1 !Solid water has indicator "1"
        gl%solid_pos = 1 !Water on position 1 in the fluid vector
        gl%solidtype_akt_phase = 1
        gl%solidpos_akt_phase = 1

        if((gl%seawater) .or. (gl%el_present)) then
            !gl%seacalc = .true.
            !call sea_trip(gl, temp, press, d, ttrip, ptrip)

            iFlash = 8 !for triple point
            x_known(1) = 1.D0
            x_fluid1 = 0.d0
            x_fluid2 = 0.d0
            x_sol = 0.d0
            x_fluid1(1) = 1.D0
            x_fluid2(1) = 1.D0
            x_hyd = 0.D0
            x_sol(1) = 1.D0

            rhovap_est = 0.D0
            rholiq_est = 0.D0

            iFlash = 8
            iphase = 2
            temptrip = 250.d0
            presstrip = 500.d-6


            if(gl%seawater) gl%seacalc = .true.
            if(gl%el_present) gl%gecarrier = .true.

            call ptflash_solid_NC_3P(gl,presstrip, Temptrip, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhovap_est, &
                & rholiq_est, Phasefrac, iFlash, iphase, iter, errval)
            gl%seacalc  =.false.
            gl%gecarrier = .false.

            gl%ttp(fluidnr) = temptrip
            gl%ptp(fluidnr) = presstrip
            ! pfreeze =
            ! psub_sea =
            ! tsub_sea =
            ! tfreeze =
        end if

        ! Since the melting curve has a minimum temperature for water, check up front if the
        ! given temperature and pressure combination is in the solid region

        if ((temp < 256.164d0) .and. (temp < tfreeze)) then  !AREA BELOW THE TEMPERATURE MINIMUM, WHERE SOLIDS MAY FORM (256.164 K, lower temperature limit for ICE I)

            !Check the limits of the solid EOS
            if (press > gl%pmax_WaterIce) then
                errval = -19932
                return
            end if

            !Calculate the sublimation pressure at specified temperature
            psub = psub_eq(gl,Temp, fluidnr)
            !if(gl%seawater) psub = psub_sea(gl, t, p, d)
            !Check, whether the specified pressure is close to the sublimation pressure
            !If it is sufficiently close, recalculate the sublimation pressure with the solid water equation
            if ((dabs(press - psub) / psub) < 0.05D0) then
                iflash = 2  !T given, compute p
                iPhase = 2  !fluid phase is vapor
                call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
                if (errval /= 0) then
                    return
                end if
            end if

            if (press > psub) then  !Solid phase
                phasetype(1) = 4
                rho(phasetype(1)) = 1.d0 / v_WaterIce(gl,Temp,press)
            else                                    !vapor phase
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
            end if

        elseif (Temp < (gl%ttp(fluidnr) + 1.D-12)) then       !AREA FROM 256.164 K TO TRIPLE TEMP. ICE I, ICE III or ICE IV may form

            !Calculate the sublimation pressure at specified temperature
            psub = psub_eq(gl,Temp, fluidnr)
            !if(gl%seawater) psub = psub_sea(gl, t, p, d)
            !Check, whether the specified pressure is close to the sublimation pressure
            !If it is sufficiently close, recalculate the sublimation pressure with the solid water equation
            if ((dabs(press - psub) / psub) < 0.05D0) then
                iflash = 2  !T given, compute p
                iPhase = 2  !fluid phase is vapor
                call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
                if (errval /= 0) then
                    return
                end if
            end if

            pmelt = pmelt_eq(gl,Temp, fluidnr)
            !if(gl%seawater) then pmelt = pfreeze !really necessary,
            !Check, whether the specified pressure is close to the melting pressure
            !If it is sufficiently close, recalculate the melting pressure with the solid water equation
            if ((dabs(press - pmelt) / pmelt) < 0.05D0) then

                iflash = 2  !T given, compute p
                iPhase = 1  !fluid phase is liquid
                call ptflash_solid_NC_2P(gl,pmelt, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
                if (errval /= 0) then
                    return
                end if
                if(gl%seawater) then
                    pmelt = freeze_press(gl, temp, pmelt, 56000.d0)
                    rho_sol = 1.d0 / v_waterice(gl, temp, pmelt)
                end if

            end if

            if ((press > psub) .and. (press < pmelt)) then  !Solid phase
                phasetype(1) = 4
                rho(phasetype(1)) = 1.d0 / v_WaterIce(gl,Temp,press)
            elseif ((press > gl%pmelt_high) .and. (gl%pmelt_high > 0.D0)) then
                errval = -19915
                return
            else !Fluid phase
                d = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
                if(gl%seawater) dw = d_sea_calc(gl, temp, press, d)
                !d2pdd2 = D2PDD2_CALC (gl,Temp, rho(1), fluidnr)
                if (d < gl%rhoc(1)) then
                    !vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    d = 0.D0
                else
                    !liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    !if(gl%seawater) rho(phasetype(1)) = dw
                    d = 0.D0
                end if
            end if


        else    !Fluid phase or solid phase possible. ICE V, ICE VI or ICE VII may form

            if (temp < 715.D0) then !715 K upper temperature limit for ice VII equation
                pmelt = pmelt_eq(gl,Temp, fluidnr)
                !if(gl%seawater) pmelt = freeze_press(gl, temp, pmelt, rhomix_calc(gl, temp, press, 0.d0, 1, 1))
                !seawater pmelt here not necessary?
                if ((press > pmelt) .and. (pmelt > 0.D0)) then
                    errval = -19915
                    return
                end if
            end if

            if (press < gl%pmaxfluid(fluidnr)) then
                d = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
                !if(gl%seawater) dw = d_sea_calc(gl, temp, press, d)    !not necessary as p_max_sea_ < p_crit_water
                !d2pdd2 = D2PDD2_CALC (gl,Temp, rho(1), fluidnr)
                if (d < gl%rhoc(1)) then
                    !vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    d = 0.D0
                else
                    !liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    !if(gl%seawater) rho(phasetype(1)) = dw
                    d = 0.D0
                end if
            else
                errval = -9932
                return
            end if
        end if


    elseif (gl%solidtype(1) == 2) then         !CO2

        gl%solidtype(1) = 2 !Solid CO2 has indicator "2"
        gl%solid_pos = 1 !CO2 on position 1 in the fluid vector
        gl%solidtype_akt_phase = 2
        gl%solidpos_akt_phase = 1

        if (temp < gl%tmin_DryIce) then
            errval = -19912
            return
        end if

        if (Temp < (gl%ttp(fluidnr) + 1.D-12)) then           !ARE BELOW THE TRIPLE POINT. VAPOR OR ICE POSSIBLE
            psub = psub_eq(gl,Temp, fluidnr)
            !Check, whether the specified pressure is close to the sublimation pressure
            !If it is sufficiently close, recalculate the sublimation pressure with the solid co2 equation
            if ((dabs(press - psub) / psub) < 0.05D0) then
                iflash = 2  !T given, compute p
                iPhase = 2  !fluid phase is vapor
                call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
                if (errval /= 0) then
                    return
                end if
            end if

            if (press > psub) then  !Solid phase
                if (press < gl%pmax_DryIce) then
                    phasetype(1) = 4
                    rho(phasetype(1)) = 1.d0 / v_DryIce(gl,Temp,press)
                else
                    errval = -19932
                    return
                end if
            else                                    !Fluid phase
                !vapor phase
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
            end if
        else                !ARE ABOVE THE TRIPLE POINT. VAPOR, LIQUID, OR ICE POSSIBLE

            if (temp < 300.D0) then !300 K upper temperature limit for ice equation
                pmelt = pmelt_eq(gl,Temp, fluidnr)

                !Check, whether the specified pressure is close to the melting pressure
                !If it is sufficiently close, recalculate the melting pressure with the solid water equation
                if ((dabs(press - pmelt) / pmelt) < 0.05D0) then

                    !Check if the solid water equation is valid for the specified temperature
                    !Elsewise use the pressure calculated with the auxiliary equation
                    if (temp < 300.D0) then

                        iflash = 2  !T given, compute p
                        iPhase = 1  !fluid phase is liquid
                        call ptflash_solid_NC_2P(gl,pmelt, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
                        if (errval /= 0) then
                            return
                        end if


                    end if

                end if

                !if ((press > pmelt) .and. (pmelt > 0.D0)) then
                !    errval = -19915
                !    return
                !end if

            end if

            if ((press > pmelt) .and. (pmelt > 0.D0)) then  !Solid phase
                if (press < gl%pmax_DryIce) then
                    phasetype(1) = 4
                    rho(phasetype(1)) = 1.d0 / v_DryIce(gl,Temp,press)
                else
                    errval = -19932
                    return
                end if
            else                                    !Fluid phase
                d = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
                !d2pdd2 = D2PDD2_CALC (gl,Temp, rho(1), fluidnr)
                if (d < gl%rhoc(1)) then
                    !vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    d = 0.D0
                else
                    !liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    d = 0.D0
                end if
            end if

        end if

    else !No equation for the solid phase of this substance available

        gl%solidtype_akt_phase = 0
        gl%solidpos_akt_phase = 1

        if (Temp < gl%ttp(fluidnr)) then
            psub = psub_eq(gl,Temp, fluidnr)
            if ((psub > 0.D0) .and. (press > psub)) then  !Solid phase (only if auxiliary equation for psub exists, check if solid forms)
                errval = -19915
                phasetype(1) = 4
                return
            else                                    !Fluid phase
                d = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
                !d2pdd2 = D2PDD2_CALC (gl,Temp, rho(1), fluidnr)
                if (d < gl%rhoc(1)) then
                    !vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    d = 0.D0
                else
                    !liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    d = 0.D0
                end if
            end if
        else
            if (Temp < gl%pmeltmaxtemp(fluidnr)) then
                pmelt = pmelt_eq(gl,Temp, fluidnr)
            end if
            if ((pmelt > 0.D0) .and. (press > pmelt)) then  !Solid phase (only if auxiliary equation for pmelt exists, check if solid forms)
                errval = -19915
                phasetype(1) = 4
                return
            else                                    !Fluid phase
                d = rhomix_calc(gl,Temp, press, 0.d0, 0, fluidnr)
                !d2pdd2 = D2PDD2_CALC (gl,Temp, rho(1), fluidnr)
                if (d < gl%rhoc(1)) then
                    !vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    d = 0.D0
                else
                    !liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    d = 0.D0
                end if
            end if
        end if


    end if



    !!MB#####################################
    !!Checking the phasetype of the seawater,which can be only "phasetype(1) = 3", with respect to the determined phasetype of water
    !!modified bs
    !if (gl%seawater == .true.) then   !Das ist sowieso nur f reinstoffe, daher muss man nicht auf die Position "fluidnr" achten, oder?
    !
    !    call phasecheck_sea(gl, press, temp, d, rho, phasetype, phasefrac, nrofphases, errval)
    !    !phasefrac(phasetype(1)) = 1.D0
    !    !nrofphases = 1
    !    !if (phasetype(1) == 4) then
    !    !    gsea = chempot_sea_calc(gl,temp, press, rhomix_calc(gl, temp, press, 0.d0, 1, fluidnr))
    !    !    gice = g_WaterIce(gl,temp, press) / gl%wm(1)
    !    !    if ( gsea .le. gice ) then
    !    !        phasetype(1) = 3
    !    !        rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 1, fluidnr) !muss das nicht die Dichte von Seewasser sein?
    !    !        d = 0.D0
    !    !    else
    !    !        phasetype(1) = 4
    !    !        !errval = -19915
    !    !        !return
    !    !    endif
    !    !
    !    !elseif (phasetype(1) == 1) then
    !    !
    !    !    gsea = chempot_sea_calc(gl,temp, press, rhomix_calc(gl, temp, press, 0.d0, 1, fluidnr))
    !    !    gvap = g_calc(gl,temp, rho(phasetype(1)), 1) / gl%wm(1)
    !    !
    !    !    if ( gsea .le. gvap ) then !hardcoded nrsubst=1, because of water which is due to the hydrate always put on the first place
    !    !        phasetype(1) = 3
    !    !        rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 1, fluidnr)
    !    !        d = 0.D0
    !    !    else
    !    !        phasetype(1) = 1
    !    !        !errval = -19915 !MB muss noch der entsprechende fehlercode f vapor anstatt solid rein
    !    !        return
    !    !    endif
    !    !endif
    !
    !    !if(gl%seacalc)then
    !    !rho(phasetype(1)) = d_sea_calc(gl, temp, press, rho(phasetype(1)))
    !    !end if
    !endif



    phasefrac(phasetype(1)) = 1.D0
    nrofphases = 1
    d = rho(phasetype(1))

    end subroutine
    !************************************************************************************


    !************************************************************************************
    subroutine PhaseDet_pure_td_sol (gl,press, Temp, d, rho, phasetype, phasefrac, nrofphases, errval)
    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   d           - density [mol/m³]
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - density [mol/m³] of all phases
    !
    !   The following two variables indicate the phase type. For a pure substance,
    !   a maximum of 3 phases can be in equilibrium
    !   phasetype   - phase identifier
    !   phasefrac   - Molar phase fraction

    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************








    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    integer :: nrofphases
    integer :: errval
    integer :: fluidnr

    double precision:: pmelt, psub, psat

    double precision :: rhomax_dryIce, rhomax_waterice

    double precision, dimension(30):: x_known, x_solid, x_hyd, x_fluid
    double precision :: rhofluid_est, beta_loc
    double precision :: rho_s_sub, rho_v_sub, rho_s_melt, rho_l_melt, rho_v_sat, rho_l_sat, rho_l_high
    integer :: iflash, iPhase, iter
    double precision ::d_w_new, d_water, tfreeze, pfreeze, p_sea, temptrip, presstrip!seawater variables
    double precision :: ptrip, ttrip, dw, dsat_w, rho_l_melt_w, icefrac, rho_l_high_water, rho_l_melt_water, vapfrac  !seawater
    double precision, dimension(3) :: rho_sol


    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol
    double precision ::  rhovap_est, rholiq_est



    phasefrac = 0.D0
    phasetype = 0
    rho = 0.D0

    fluidnr = 1

    x_known = 0.D0
    x_solid = 0.D0
    x_hyd = 0.D0
    x_fluid = 0.D0
    x_known(1) = 1.D0
    rhofluid_est = 0.D0

    !In case two melting pressures at given temperature exist for solid water, the higher one is saved in the variable pmelt_high
    gl%pmelt_high = 0.D0

    pmelt = 0.D0
    psub = 0.D0

    rho_s_sub = 0.D0
    rho_v_sub = 0.D0
    rho_s_melt = 0.D0
    rho_l_melt = 0.D0
    rho_v_sat = 0.D0
    rho_l_sat = 0.D0
    rho_l_high = 0.D0
    tfreeze = 256.164d0 !for freezeing temp of water, if decision more simple with this initialization


    if (gl%solidtype(1) == 1) then      !WATER

        gl%solidtype(1) = 1 !Solid water has indicator "1"
        gl%solid_pos = 1 !Water on position 1 in the fluid vector
        gl%solidtype_akt_phase = 1
        gl%solidpos_akt_phase = 1
        ttrip = gl%ttp(fluidnr)
        ptrip = gl%ptp(fluidnr)

        ! Since the melting curve has a minimum temperature for water, check up front if the
        ! given temperature and pressure combination is in the solid region

        if(gl%seawater .or. gl%el_present) then


            errval = -12800
            return

            iFlash = 8 !for triple point
            x_known(1) = 1.D0
            x_fluid1 = 0.d0
            x_fluid2 = 0.d0
            x_sol = 0.d0
            x_fluid1(1) = 1.D0
            x_fluid2(1) = 1.D0
            x_hyd = 0.D0
            x_sol(1) = 1.D0

            rhovap_est = 0.D0
            rholiq_est = 0.D0

            iFlash = 8
            iphase = 2
            temptrip = 265.d0
            presstrip = 500.d-6

            call ptflash_solid_NC_3P(gl,presstrip, Temptrip, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhovap_est, &
                & rholiq_est, Phasefrac, iFlash, iphase, iter, errval)

            !gl%seacalc = .false.
            gl%ttp(fluidnr) = temptrip
            gl%ptp(fluidnr) = presstrip

            !call sea_trip(gl, 273.2d0, 0.6d-3, 56000.d0, ttrip, ptrip)
            !gl%ttp(fluidnr) = ttrip
            !gl%ptp(fluidnr) = ptrip
        end if

        if ((temp < 256.164d0) .and. (temp < tfreeze)) then  !AREA BELOW THE TEMPERATURE MINIMUM, WHERE SOLIDS MAY FORM (256.164 K, lower temperature limit for ICE I)

            !Check the limits of the solid EOS
            rhomax_waterice = 1.D0 / v_waterIce(gl,Temp, gl%pmax_WaterIce)
            if (d > rhomax_waterice) then
                errval = -19914
                return
            end if

            !Get an inital estimate for the sublimation pressure
            psub = psub_eq(gl,Temp, fluidnr)

            !Calculate p_sub for seawater?! really necessary?

            !Calculate the sublimation pressure with the solid water equation
            iflash = 2  !T given, compute p
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Density of solid water at the sublimation line
            rho_s_sub = rho_sol(2)!1.D0 / v_WaterIce(gl,Temp,psub)
            !Density of gaseous water at the sublimation line
            rho_v_sub = rho_sol(1)!rhomix_calc(gl,Temp, psub, 0.d0, 2, fluidnr)

            if (d < rho_v_sub) then             !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif (d > rho_s_sub) then         !solid
                press = p_WaterIce(gl,Temp, d, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            else                                !sublimation region
                press = psub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_s_sub))/((1.D0/rho_v_sub)-(1.D0/rho_s_sub))
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                rho(phasetype(2)) = rho_s_sub
                nrofphases = 2
            end if


        elseif (Temp < (ttrip + 1.D-12))then       !AREA FROM 256.164 K TO TRIPLE TEMP. ICE I, ICE III or ICE V may form   !ice and vapor

            !Get an inital estimate for the sublimation pressure
            psub = psub_eq(gl,Temp, fluidnr)
            !Calculate the sublimation pressure with the solid water equation
            iflash = 2  !T given, compute p
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            rho_v_sub = rho_sol(1)!rhomix_calc(gl,Temp, psub, 0.d0, 2, fluidnr)
            rho_s_sub = rho_sol(2)!1.D0 / v_WaterIce(gl,Temp, psub)

            !Get an inital estimate for the melting pressure
            pmelt = pmelt_eq(gl,Temp, fluidnr)
            !if(gl%seawater) pmelt = freeze_press(gl, temp, press, 56000.d0)
            iflash = 2  !T given, compute p
            iPhase = 1  !fluid phase is liquid              !here comes the seawater!!??
            call ptflash_solid_NC_2P(gl,pmelt, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            rho_s_melt = rho_sol(2)!1.D0 / v_WaterIce(gl,Temp, pmelt)
            rho_l_melt = rho_sol(1)!rhomix_calc(gl,Temp, pmelt, 0.d0, 1, fluidnr)  seawater density???!!
            if(gl%seawater) then
                pmelt = freeze_press(gl, temp, pmelt, rho_l_melt)
                rho_l_melt_water =  rhomix_calc(gl, temp, pmelt, 0.d0,1,1)
                rho_l_melt = d_sea_calc(gl, temp, pmelt, rhomix_calc(gl, temp, pmelt, 0.d0,1,1))
                rho_s_melt = 1.D0 / v_WaterIce(gl,Temp, pmelt)
            end if
            rho_l_high = rhomix_calc(gl,Temp, gl%pmelt_high, 0.d0, 0, fluidnr)

            if (d < rho_v_sub) then                                     !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > rho_v_sub) .and. (d < rho_s_sub)) then         !sublimation region
                press = psub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_s_sub))/((1.D0/rho_v_sub)-(1.D0/rho_s_sub))
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                rho(phasetype(2)) = rho_s_sub
                nrofphases = 2
            elseif ((d > rho_s_sub) .and. (d < rho_s_melt)) then        !solid Ih
                press = p_WaterIce(gl,Temp, d, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > rho_s_melt) .and. (d < rho_l_melt)) then         !melting region SLE
                press = pmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_s_melt))/((1.D0/rho_l_melt)-(1.D0/rho_s_melt))    !seawater mix here
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                rho(phasetype(2)) = rho_s_melt
                nrofphases = 2
                if(gl%seawater) then
                    icefrac = icefrac_sea(gl, temp, press, rho_l_melt_water)
                    phasefrac(phasetype(2)) = icefrac
                    phasefrac(phasetype(1)) = 1.d0 -icefrac
                end if
            elseif ((d > rho_l_melt) .and. (d < rho_l_high)) then         !liquid region
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 3
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d           !seawaterdensity here
            else
                errval = -19914
                return
            end if


        elseif  (Temp < gl%tc(fluidnr)) then       !AREA FROM TRIPLE TEMP TO CRITICAL TEMP. ICE V, or ICE IV may form

            !Compute the vapor pressure and the dew and bubble densities
            iflash = 1  !T given
            call Flash_Pure_PhaseBoundary(gl,psat, Temp, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if(gl%seawater) then
                psat = boil_press(gl, temp, 0.1d0, 56000.d0)
                rho_v_sat = rhomix_calc(gl, temp, psat, 0.d0, 2, 1)
                dsat_w = rhomix_calc(gl, temp, psat, 0.d0, 1, 1)
                rho_l_sat = d_sea_calc(gl, temp, psat, dsat_w)
            end if
            if (errval /= 0) then
                return
            end if

            !Compute maximum pressure and the corresponding density
            pmelt = pmelt_eq(gl,Temp, fluidnr)
            rho_l_high = rhomix_calc(gl,Temp, pmelt, 0.d0, 0, fluidnr)  !seawater
            if(gl%seawater) then
                pmelt = freeze_press(gl, temp, pmelt, rho_l_sat)
                rho_l_high_water = rhomix_calc(gl,temp, pmelt, 0.d0, 1, 1)
                rho_l_high = d_sea_calc(gl, temp, pmelt, rho_l_high_water)
            end if

            !In case density calculation failes, take dummy value, Andreas September 2014
            if (rho_l_high < 1.D-14) then
                rho_l_high = gl%rhomaxfluid(fluidnr)
            end if

            if (d < rho_v_sat) then                                     !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > rho_v_sat) .and. (d < rho_l_sat)) then         !two phase (VLE) region
                press = psat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_l_sat))/((1.D0/rho_v_sat)-(1.D0/rho_l_sat))       !seawatermix here
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                rho(phasetype(2)) = rho_l_sat
                nrofphases = 2
                if(gl%seawater) then
                    vapfrac = vapfrac_sea(gl, temp, psat, dsat_w)
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.d0 - vapfrac
                end if
            elseif ((d > rho_l_sat) .and. (d < rho_l_high)) then        !liquid
                press = P_CALC(gl,Temp, d, fluidnr)
                if(gl%seawater) press = p_w_new(gl, temp, 0.1d0, d)
                phasetype(1) = 3
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
                if(gl%seawater) temp = temp !seawater here
            else
                errval = -19914
                return
            end if


        else     !Fluid phase or solid phase possible. ICE VII may form

            if (temp < 715.D0) then

                pmelt = pmelt_eq(gl,Temp, fluidnr)
                rho_l_high = rhomix_calc(gl,Temp, pmelt, 0.d0, 0, fluidnr)      !seawater

                if (rho_l_high > 1.D-12) then   !Reasonable liquid density found

                    if (d > rho_l_high) then
                        errval = -19914
                        return
                    end if
                end if

            end if

            !Fluid phase
            press = P_CALC(gl,Temp, d, fluidnr)             !seawater

            !d2pdd2 = D2PDD2_CALC (gl,Temp, d, fluidnr)
            if (d < gl%rhoc(1)) then
                !vapor phase
                phasetype(1) = 1
                rho(phasetype(1)) = d
            else
                !liquid phase
                phasetype(1) = 3
                rho(phasetype(1)) = d               !seawater
                if(gl%seawater) press = p_w_new(gl, temp, press, d)
            end if
            phasefrac(phasetype(1)) = 1.D0
            nrofphases = 1

        end if


    elseif (gl%solidtype(1) == 2) then     !CO2

        gl%solidtype(1) = 2 !Solid CO2 has indicator "2"
        gl%solid_pos = 1 !CO2 on position 1 in the fluid vector
        gl%solidtype_akt_phase = 2
        gl%solidpos_akt_phase = 1

        if (temp < gl%tmin_DryIce) then
            errval = -19912
            return
        end if

        if ((Temp < (gl%ttp(fluidnr) + 1.D-12)) .or. (Temp < (ttrip + 1.d-12))) then           !BELOW THE TRIPLE POINT. VAPOR OR ICE POSSIBLE

            !Check the limits of the solid EOS
            rhomax_dryice = 1.D0 / v_DryIce(gl,Temp, gl%pmax_DryIce)
            if (d > rhomax_dryice) then
                errval = -19914
                return
            end if

            !Get inital estimate for sublimation pressure
            psub = psub_eq(gl,Temp, fluidnr)
            !Calculate the sublimation pressure with the dry ice equation
            iflash = 2  !T given, compute p
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Density of solid water at the sublimation line
            rho_s_sub = rho_sol(2)!1.D0 / v_DryIce(gl,Temp,psub)
            !Density of gaseous water at the sublimation line
            rho_v_sub = rho_sol(1)!rhomix_calc(gl,Temp, psub, 0.d0, 2, fluidnr)

            if (d < rho_v_sub) then             !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif (d > rho_s_sub) then         !solid
                press = p_DryIce(gl,Temp, d, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            else                                !sublimation region SVE
                press = psub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_s_sub))/((1.D0/rho_v_sub)-(1.D0/rho_s_sub))
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                rho(phasetype(2)) = rho_s_sub
                nrofphases = 2
            end if

        elseif (Temp < gl%tc(fluidnr)) then               !TEMPERATURE BETWEEN TRIPLE POINT AND CRITICAL TEMPERATURE. VAPOR, LIQUID, OR SOLID POSSIBLE

            !Check the limits of the solid EOS
            rhomax_dryice = 1.D0 / v_DryIce(gl,Temp, gl%pmax_DryIce)
            if (d > rhomax_dryice) then
                errval = -19914
                return
            end if

            !Compute the vapor pressure and the dew and bubble densities
            iflash = 1  !T given
            call Flash_Pure_PhaseBoundary(gl,psat, Temp, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)

            if (errval /= 0) then
                return
            end if

            !Initial Estimate for the melting pressure
            pmelt = pmelt_eq(gl,Temp, fluidnr)
            !Compute melting pressure with the dry ice equation (ASSUME THE SOLID EOS VALID UP TO CRITICAL TEMPERATURE. THIS IS NOT EXACTLY TRUE, BUT WILL GIVE REASONABLE ESTIMATES FOR THE MELTING PRESSURE)
            iflash = 2  !T given, compute p
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if
            rho_s_melt = rho_sol(2)!1.D0 / v_DryIce(gl,Temp, pmelt)
            rho_l_melt = rho_sol(1)!rhomix_calc(gl,Temp, pmelt, 0.d0, 1, fluidnr)


            if (d < rho_v_sat) then                                     !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > (rho_v_sat-1.D-12)) .and. (d < rho_l_sat)) then         !two phase (VLE) region
                press = psat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_l_sat))/((1.D0/rho_v_sat)-(1.D0/rho_l_sat))
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                rho(phasetype(2)) = rho_l_sat
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
            elseif ((d > (rho_l_sat-1.D-12)) .and. (d < rho_l_melt)) then        !liquid
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 3
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > (rho_l_melt-1.D-12)) .and. (d < rho_s_melt)) then        !two phase SLE region
                press = pmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_s_melt))/((1.D0/rho_l_melt)-(1.D0/rho_s_melt))
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                rho(phasetype(2)) = rho_s_melt
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
            else                                                                !Solid CO2
                press = p_DryIce(gl,Temp, d, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            end if

        else                                    !TEMPERATURE ABOVE THE CRITICAL TEMPERATURE. ONLY FLUID SOLUTIONS POSSIBLE

            !Fluid phase
            press = P_CALC(gl,Temp, d, fluidnr)
            !d2pdd2 = D2PDD2_CALC (gl,Temp, d, fluidnr)
            if (d < gl%rhoc(1)) then
                !vapor phase
                phasetype(1) = 1
                rho(phasetype(1)) = d
            else
                !liquid phase
                phasetype(1) = 3
                rho(phasetype(1)) = d
            end if
            phasefrac(phasetype(1)) = 1.D0
            nrofphases = 1


        end if

    else !No equation for the solid phase of this substance available

        gl%solidtype_akt_phase = 0
        gl%solidpos_akt_phase = 1

        if (Temp < gl%ttp(fluidnr)) then           !BELOW THE TRIPLE POINT. VAPOR OR ICE POSSIBLE

            !Get sublimation pressure
            psub = psub_eq(gl,Temp, fluidnr)

            if (psub > 0.D0) then
                !Density of gaseous water at the sublimation line
                rho_v_sub = rhomix_calc(gl,Temp, psub, 0.d0, 2, fluidnr)
            end if

            if (d < rho_v_sub) then             !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            else                                !Solid forms
                errval = -19915
                phasetype(1) = 4
                return
            end if

        elseif (Temp < gl%tc(fluidnr)) then               !TEMPERATURE BETWEEN TRIPLE POINT AND CRITICAL TEMPERATURE. VAPOR, LIQUID, OR SOLID POSSIBLE

            !Compute the vapor pressure and the dew and bubble densities
            iflash = 1  !T given
            call Flash_Pure_PhaseBoundary(gl,psat, Temp, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)

            if (errval /= 0) then
                return
            end if

            !Compute the melting pressure
            pmelt = pmelt_eq(gl,Temp, fluidnr)
            rho_l_melt = rhomix_calc(gl,Temp, pmelt, 0.d0, 1, fluidnr)


            if (d < rho_v_sat) then                                     !vapor
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 1
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            elseif ((d > (rho_v_sat-1.D-12)) .and. (d < rho_l_sat)) then         !two phase (VLE) region
                press = psat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = ((1.D0/d)-(1.D0/rho_l_sat))/((1.D0/rho_v_sat)-(1.D0/rho_l_sat))
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                rho(phasetype(2)) = rho_l_sat
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
            elseif ((d > (rho_l_sat-1.D-12)) .and. (d < rho_l_melt)) then        !liquid
                press = P_CALC(gl,Temp, d, fluidnr)
                phasetype(1) = 3
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
                rho(phasetype(1)) = d
            else                                                                !Solid
                errval = -19915
                phasetype(1) = 4
                return
            end if

        else                                    !TEMPERATURE ABOVE THE CRITICAL TEMPERATURE. ONLY FLUID SOLUTIONS POSSIBLE

            !If melting pressure equation exists at the given temperature, check whether solid forms
            if (Temp < gl%pmeltmaxtemp(fluidnr)) then
                !Compute the melting pressure
                pmelt = pmelt_eq(gl,Temp, fluidnr)
                rho_l_melt = rhomix_calc(gl,Temp, pmelt, 0.d0, 1, fluidnr)
                if (d > rho_l_melt) then
                    errval = -19915
                    phasetype(1) = 4
                    return
                end if
            end if

            !Fluid phase
            press = P_CALC(gl,Temp, d, fluidnr)
            !d2pdd2 = D2PDD2_CALC (gl,Temp, d, fluidnr)
            if (d < gl%rhoc(1)) then
                !vapor phase
                phasetype(1) = 1
                rho(phasetype(1)) = d
            else
                !liquid phase
                phasetype(1) = 3
                rho(phasetype(1)) = d
            end if
            phasefrac(phasetype(1)) = 1.D0
            nrofphases = 1


        end if


    end if

    end subroutine
    !************************************************************************************



    !************************************************************************************
    subroutine PhaseDet_pure_ph_sol (gl,press, Temp, d, rho, phasetype, phasefrac, h_spec, nrofphases, errval)
    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   h_spec      - overall specified enthalpy [J/mol]
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - density [mol/m³] of all phases
    !   d           - overall density [mol/m³]
    !
    !   The following two variables indicate the phase type. For a pure substance,
    !   a maximum of 3 phases can be in equilibrium
    !   phasetype   - phase identifier
    !   phasefrac   - Molar phase fraction

    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************








    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    double precision :: h_spec

    integer :: nrofphases
    integer :: errval
    integer :: fluidnr

    double precision:: pmelt, psub, psat
    double precision:: Tmelt, Tsub, Tsat

    double precision :: rhomax_dryIce, rhomax_waterice, p_tmeltmin, p_tmeltmax_ICE_VII, rhomax, rhomin
    double precision :: hmin, hmax

    double precision, dimension(30):: x_known, x_solid, x_hyd, x_fluid, x_fluid1, x_fluid2, x_sol
    double precision :: rhofluid_est, beta_loc
    double precision :: h_s_sub, h_v_sub, h_s_melt, h_l_melt, h_v_sat, h_l_sat, h_l_high
    double precision :: rho_s_sub, rho_v_sub, rho_s_melt, rho_l_melt, rho_v_sat, rho_l_sat, rho_l_high
    integer :: iflash, iPhase, iter

    double precision, dimension(3) :: rho_sol

    ! needed to transmit T, p to Regula Falsi:
    double precision :: t_min
    double precision :: t_max
    double precision :: t_min_allowed
    double precision :: t_max_allowed
    double Precision :: DeltaH_allowed, press_rho_l_melt_water, press_rho_l_sat
    type(type_additional_parameters) :: parameters
    integer :: max_iterations
    integer :: iterations
    double precision ::  rhovap_est, rholiq_est
    double precision :: ptrip, ttrip, dw, dsat_w, rho_l_melt_w, icefrac, &
        & rho_l_melt_water, mix_sal, vapfrac,temptrip, presstrip   !seawater

    psat= 0.d0
    Tsat= 0.d0
    !parameters = 0.d0
    h_s_sub = 0.D0
    h_v_sub = 0.D0
    h_s_melt = 0.D0
    h_l_melt = 0.D0
    h_v_sat = 0.D0
    h_l_sat = 0.D0
    h_l_high = 0.D0

    rho_sol = 0.D0
    x_known = 0.D0
    x_solid = 0.D0
    x_hyd = 0.D0
    x_fluid = 0.D0
    x_solid(1) = 1.D0
    x_fluid(1) = 1.D0
    x_known(1) = 1.D0
    rhofluid_est = 0.D0

    rho_v_sat = 0.d0
    rho_l_sat = 0.d0

    hmin = 0.D0
    hmax = 0.D0

    fluidnr = 1

    phasetype = 0

    !Variables for the enthalpy iteration of the fluid phase at given pressure
    t_min = 0.D0
    t_max = 0.D0
    t_min_allowed = 0.D0
    t_max_allowed = 0.D0
    DeltaH_allowed = 1.d-8
    max_iterations = 200

    if(gl%seawater) mix_sal = gl%sea%salinity

    if (gl%solidtype(1) == 1) then     !WATER

        gl%solidtype(1) = 1 !Solid water has indicator "1"
        gl%solid_pos = 1 !Water on position 1 in the fluid vector
        gl%solidtype_akt_phase = 1
        gl%solidpos_akt_phase = 1

        p_tmeltmin = 170.184416372448D0                 !Pressure at temperature minimum of Ice Ih. Calculated with the dry ice equation
        p_tmeltmax_ICE_VII = 20617.7780087102D0         !Pressure maximum of ICE VII equation (out of bounds of water equation!)

        if(gl%seawater .or. gl%el_present) then

            errval = -12800
            return


            iFlash = 8 !for triple point
            x_known(1) = 1.D0
            x_fluid1 = 0.d0
            x_fluid2 = 0.d0
            x_sol = 0.d0
            x_fluid1(1) = 1.D0
            x_fluid2(1) = 1.D0
            x_hyd = 0.D0
            x_sol(1) = 1.D0

            rhovap_est = 0.D0
            rholiq_est = 0.D0

            iFlash = 8
            iphase = 2
            temptrip = 265.d0
            presstrip = 500.d-6

            call ptflash_solid_NC_3P(gl,presstrip, Temptrip, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhovap_est, &
                & rholiq_est, Phasefrac, iFlash, iphase, iter, errval)

            !gl%seacalc = .false.
            gl%ttp(fluidnr) = temptrip
            gl%ptp(fluidnr) = presstrip

            !call sea_trip(gl, temp, press, d, ttrip, ptrip)
            !gl%ttp(fluidnr) = ttrip
            !gl%ptp(fluidnr) = ptrip
            !! pfreeze =
            ! tfreeze =
        end if

        !Check the boundaries of the solid equation
        hmin = h_WaterIce(gl,gl%Tmin_Waterice, press)
        if (h_spec < hmin) then
            errval = -19941
            return
        end if

        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        hmax = H_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
        if (h_spec > hmax) then
            errval = -19941
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature from the ancillary equation
            psub = press
            Tsub = tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then
                return
            end if
            !Calculate the sublimation temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Tsub, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in sublimation equilibrium
            rho_s_sub = rho_sol(2) !1.D0 / v_WaterIce(gl,Tsub, psub)
            h_s_sub = h_WaterIce(gl,Tsub, psub)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_v_sub = rho_sol(1) !rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
            h_v_sub = H_CALC(gl,Tsub, rho_v_sub, fluidnr)

            if (h_spec <= h_s_sub) then              !solid water Ih region
                Temp = T_WaterIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (h_spec > h_v_sub) then        !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tsub
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Sublimation region, SVE
                Temp = Tsub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (h_spec-h_s_sub)/(h_v_sub-h_s_sub)
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_sub
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            !if(gl%seawater) Tmelt = freeze_temp(gl, tmelt, pmelt, 56000.d0)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of liquid in melting equilibrium
            rho_l_melt = rho_sol(1) !rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
            h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)


            rho_s_melt = rho_sol(2)
            !
            !if(gl%seawater) then
            !Tmelt = freeze_temp(gl, tmelt, pmelt, 56000.d0)
            !rho_l_melt = rhomix_calc(gl, tmelt, pmelt, 55500.d0)
            !h_l_helt = h_sea_calc(gl, tmelt, pmelt, rho_l_melt)
            !rho_l_melt = d_sea_calc(gl, tmelt, pmelt, rho_l_melt)
            !rho_s_melt = 1.d0 / v_waterice(gl, tmelt, pmelt)
            !end if
            !
            h_s_melt = h_WaterIce(gl,Tmelt, pmelt)

            !Get the density and enthalpy of solid water in melting equilibrium
            !rho_s_melt = rho_sol(2) !1.D0 / v_WaterIce(gl,Tmelt, pmelt)


            if(gl%seawater) then
                Tmelt = freeze_temp(gl, tmelt, pmelt, 56000.d0)
                rho_s_melt = 1.D0 / v_WaterIce(gl,Tmelt, pmelt)
                h_s_melt = h_WaterIce(gl,Tmelt, pmelt)
                rho_l_melt = rhomix_calc(gl, tmelt, pmelt, 0.d0, 1,1 )
                h_l_melt = h_sea_calc(gl, tmelt, pmelt, rho_l_melt)
                !rho_l_melt = d_sea_calc(gl, tmelt, pmelt, rho_l_melt)
            end if

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            !if(gl%seawater) then
            !    tsat = boil_temp(gl, temp, psat, rho_l_sat)
            !end if
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            h_l_sat = H_CALC(gl,Tsat, rho_l_sat, fluidnr)
            if(gl%seawater) then
                tsat = boil_temp(gl, tsat, psat, rho_l_sat)
                rho_l_sat = rhomix_calc(gl, tsat, psat, 0.d0, 1, 1)
                rho_v_sat = rhomix_calc(gl, tsat, psat, 0.d0, 2, 1)
                h_l_sat = h_sea_calc(gl, tsat,psat, rho_l_sat)
                !rho_l_sat = d_sea_calc(gl, tsat, psat, rho_l_sat)
            end if
            !Get the enthalpy of saturated vapor
            h_v_sat = H_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if (h_spec < h_s_melt) then              !solid water Ih region
                Temp = T_WaterIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((h_spec >= h_s_melt) .and. (h_spec <= h_l_melt)) then       !SLE region
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (h_spec-h_s_melt)/(h_l_melt-h_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))

                if(gl%seawater) then
                    icefrac = icefrac_sea(gl, temp, press, rho_l_melt_water)
                    phasefrac(phasetype(2)) = icefrac
                    phasefrac(phasetype(1)) = 1.d0 - icefrac
                    mix_sal = freeze_sal(gl, temp, press, rho_l_melt_water)
                    !definie new props here (density, ...)
                end if

            elseif ((h_spec > h_l_melt) .and. (h_spec < h_l_sat)) then        !Liquid region
                !determine temperature from p and h
                !parameters = 0.d0
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tmelt
                t_max = Tsat
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((h_spec >= h_l_sat) .and. (h_spec <= h_v_sat)) then       !SLE region
                Temp = Tsat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (h_spec-h_l_sat)/(h_v_sat-h_l_sat)
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))

                if(gl%seawater) then
                    vapfrac = vapfrac_sea(gl, temp, press, rho_l_sat)
                    phasefrac(phasetype(2)) = vapfrac
                    phasefrac(phasetype(1)) = 1.d0 - vapfrac
                    mix_sal = boil_sal(gl, temp, press,rho_l_sat)
                end if

            elseif (h_spec > h_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < p_tmeltmin) then        !Region from the maximum pressure where the Ice Ih equation is valid down to the critical pressure

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in melting equilibrium
            rho_s_melt = rho_sol(2) !1.D0 / v_WaterIce(gl,Tmelt, pmelt)
            h_s_melt = h_WaterIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_l_melt = rho_sol(1) !rhomix_calc(gl,Tmelt, pmelt, 0.d0, 2, fluidnr)
            h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            if(gl%seawater) then
                Tmelt = freeze_temp(gl, tmelt, pmelt, rho_l_melt)
                rho_l_melt = rhomix_calc(gl, tmelt, pmelt, 0.d0, 1, 1)
            end if
            !seawater also here!!

            if (h_spec < h_s_melt) then              !solid water Ih region
                Temp = T_WaterIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (h_spec > h_l_melt) then        !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Melting region, SLE
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (h_spec-h_s_melt)/(h_l_melt-h_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
                if(gl%seawater) then
                    icefrac = icefrac_sea(gl, temp, press, rho_l_melt_water)
                    phasefrac(phasetype(2)) = icefrac
                    phasefrac(phasetype(1)) = 1.d0 - icefrac
                    mix_sal = freeze_sal(gl, temp, press,rho_l_melt_water)
                    !definie new props here (density, ...)
                end if
            end if

        elseif (press < p_tmeltmax_ICE_VII) then        !pressure is smaller than the maximum pressure for ICE VII (Fluid EOS anyways only valid up to 1000 MPa, this is already catched in the beginning)

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if


            !Get the density and enthalpy of liquid in melting equilibrium
            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
            h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            if (h_spec < h_l_melt) then             !solid water region or melting region. No equation available
                errval = -19915
                phasetype(1) = 4
                return
            else                                    !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        end if



    elseif (gl%solidtype(1) == 2) then         !CO2

        gl%solidtype(1) = 2 !Solid CO2 has indicator "2"
        gl%solid_pos = 1 !CO2 on position 1 in the fluid vector
        gl%solidtype_akt_phase = 2
        gl%solidpos_akt_phase = 1

        !Check the boundaries of the solid equation
        hmin = h_DryIce(gl,gl%Tmin_DryIce, press)
        if (h_spec < hmin) then
            errval = -19941
            return
        end if

        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        hmax = H_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
        if (h_spec > hmax) then
            errval = -19941
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature
            psub = press
            Tsub =  tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then
                return
            end if
            !Calculate the sublimation temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Tsub, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid CO2 in sublimation equilibrium
            rho_s_sub = rho_sol(2) !1.D0 / v_DryIce(gl,Tsub, psub)
            h_s_sub = h_DryIce(gl,Tsub, psub)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_v_sub = rho_sol(1) !rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
            h_v_sub = H_CALC(gl,Tsub, rho_v_sub, fluidnr)

            if (h_spec < h_s_sub) then              !solid region
                Temp = T_DryIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (h_spec > h_v_sub) then        !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tsub
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Sublimation region, SVE
                Temp = Tsub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (h_spec-h_s_sub)/(h_v_sub-h_s_sub)
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_sub
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid co2 in melting equilibrium
            rho_s_melt = rho_sol(2) !1.D0 / v_DryIce(gl,Tmelt, pmelt)
            h_s_melt = h_DryIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of liquid in melting equilibrium
            rho_l_melt = rho_sol(1) !rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
            h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            h_l_sat = H_CALC(gl,Tsat, rho_l_sat, fluidnr)
            !Get the enthalpy of saturated vapor
            h_v_sat = H_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if (h_spec < h_s_melt) then              !solid region
                Temp = T_DryIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((h_spec >= h_s_melt) .and. (h_spec <= h_l_melt)) then       !SLE region
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (h_spec-h_s_melt)/(h_l_melt-h_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif ((h_spec > h_l_melt) .and. (h_spec < h_l_sat)) then        !Liquid region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tmelt
                t_max = Tsat
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((h_spec >= h_l_sat) .and. (h_spec <= h_v_sat)) then       !SLE region
                Temp = Tsat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (h_spec-h_l_sat)/(h_v_sat-h_l_sat)
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif (h_spec > h_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < gl%pmax_Dryice) then        !pressure is smaller than the maximum pressure for the dry ice equation

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in melting equilibrium
            rho_s_melt = rho_sol(2) !1.D0 / v_DryIce(gl,Tmelt, pmelt)
            h_s_melt = h_DryIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_l_melt = rho_sol(1) !rhomix_calc(gl,Tmelt, pmelt, 0.d0, 2, fluidnr)
            h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            if (h_spec < h_s_melt) then              !solid region
                Temp = T_DryIce_ph(gl,press, h_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (h_spec > h_l_melt) then        !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Melting region, SLE
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (h_spec-h_s_melt)/(h_l_melt-h_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        else        !Dry Ice equation not valid for that region (500 MPa, approximately Tcrit)

            !Get the enthalpy of liquid CO2 at the limit of validity of the dry ice equation

            !Check the boundaries of the fluid equation
            rhomin = rhomix_calc(gl,gl%tc(fluidnr), press, 0.d0, 0, fluidnr)
            if (rhomax < 1.D-12) then
                errval = -8888
                return
            end if
            hmin = H_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
            if (h_spec < hmax) then
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            end if

            !determine temperature from p and h
            parameters%a_p(1) = press
            parameters%a_p(2) = h_spec
            t_min = gl%tc(fluidnr)
            t_max = gl%tmaxfluid(fluidnr)
            t_min_allowed = t_min
            t_max_allowed = t_max
            call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                &                       Max_iterations, Iterations, errval, Parameters)

            if (errval /= 0) then
                return
            end if
            !Determine if the phase is liquid or vapor like
            d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
            if (d < gl%rhoc(fluidnr)) then
                phasetype(1) = 1
            else
                phasetype(1) = 3
            end if
            rho(phasetype(1)) = d
            if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                errval = -1235
                return
            end if
            !d = rho(phasetype(1))
            phasefrac(phasetype(1)) = 1.D0
            nrofphases = 1


        end if


    else !No equation for the solid phase of this substance available

        gl%solidtype_akt_phase = 0
        gl%solidpos_akt_phase = 1
        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        hmax = H_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
        if (h_spec > hmax) then
            errval = -19941
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature
            psub = press
            Tsub =  tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then       !Do not exit on error. If no sublimation equation exists, assume fluid
                Tsub = 0.D0
                h_v_sub = 0.D0          !Dummy value
            else
                !Get the density and enthalpy of vapor in sublimation equilibrium
                rho_v_sub = rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
                h_v_sub = H_CALC(gl,Tsub, rho_v_sub, fluidnr)
            end if



            if ((h_spec < h_v_sub) .and. (Tsub > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            else                                    !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                if (Tsub > 0.D0) then
                    t_min = Tsub
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then    !Do not exit on error. If no sublimation equation exists, assume fluid
                Tmelt = 0.D0
                h_l_melt = h_spec-1.D0         !Dummy value
            else
                !Get the density and enthalpy of liquid in melting equilibrium
                rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
                h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)
            end if

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            h_l_sat = H_CALC(gl,Tsat, rho_l_sat, fluidnr)
            !Get the enthalpy of saturated vapor
            h_v_sat = H_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if ((h_spec < h_l_melt) .and. (Tmelt > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            elseif ((h_spec > h_l_melt) .and. (h_spec < h_l_sat)) then        !Liquid region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                if (Tmelt > 0.D0) then
                    t_min = Tmelt
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = Tsat
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((h_spec >= h_l_sat) .and. (h_spec <= h_v_sat)) then       !SLE region
                Temp = Tsat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (h_spec-h_l_sat)/(h_v_sat-h_l_sat)
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif (h_spec > h_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        else        !Region above the critical pressure

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then    !Do not exit on error. If no sublimation equation exists, assume fluid
                Tmelt = 0.D0
                h_l_melt = h_spec-1.D0         !Dummy value
            else
                !Get the density and enthalpy of liquid in melting equilibrium
                rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
                h_l_melt = H_CALC(gl,Tmelt, rho_l_melt, fluidnr)
            end if


            if ((h_spec < h_l_melt) .and. (Tmelt > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            else     !liquid or supercritical region

                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = h_spec
                if (Tmelt > 0.D0) then
                    t_min = Tmelt
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(h_spec - H_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1

            end if


        end if


    end if

    if(gl%seawater) gl%sea%salinity = mix_sal

    end subroutine
    !************************************************************************************




    !************************************************************************************
    subroutine PhaseDet_pure_ps_sol (gl,press, Temp, d, rho, phasetype, phasefrac, s_spec, nrofphases, errval)
    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   s_spec      - overall specified entropy [J/molK]
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - density [mol/m³] of all phases
    !   d           - overall density [mol/m³]
    !
    !   The following two variables indicate the phase type. For a pure substance,
    !   a maximum of 3 phases can be in equilibrium
    !   phasetype   - phase identifier
    !   phasefrac   - Molar phase fraction

    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************








    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    double precision :: s_spec

    integer :: nrofphases
    integer :: errval
    integer :: fluidnr

    double precision:: pmelt, psub, psat
    double precision:: Tmelt, Tsub, Tsat

    double precision :: rhomax_dryIce, rhomax_waterice, p_tmeltmin, p_tmeltmax_ICE_VII, rhomax, rhomin
    double precision :: smin, smax

    double precision, dimension(30):: x_known, x_solid, x_hyd, x_fluid
    double precision :: rhofluid_est, beta_loc
    double precision :: s_s_sub, s_v_sub, s_s_melt, s_l_melt, s_v_sat, s_l_sat, s_l_high
    double precision :: rho_s_sub, rho_v_sub, rho_s_melt, rho_l_melt, rho_v_sat, rho_l_sat, rho_l_high
    integer :: iflash, iPhase, iter

    double precision, dimension(3) :: rho_sol

    ! needed to transmit T, p to Regula Falsi:
    double precision :: t_min
    double precision :: t_max
    double precision :: t_min_allowed
    double precision :: t_max_allowed
    double Precision :: DeltaS_allowed
    type(type_additional_parameters) :: parameters
    integer :: max_iterations
    integer :: iterations

    double precision :: ptrip, ttrip, dw, dsat_w, rho_l_melt_w  !seawater

    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol
    double precision ::  rhovap_est, rholiq_est, temptrip, presstrip

    s_s_sub = 0.D0
    s_v_sub = 0.D0
    s_s_melt = 0.D0
    s_l_melt = 0.D0
    s_v_sat = 0.D0
    s_l_sat = 0.D0
    s_l_high = 0.D0

    rho_sol = 0.D0
    x_known = 0.D0
    x_solid = 0.D0
    x_hyd = 0.D0
    x_fluid = 0.D0
    x_known(1) = 1.D0
    x_solid(1) = 1.D0
    x_fluid(1) = 1.D0
    rhofluid_est = 0.D0

    rho_v_sat = 0.d0
    rho_l_sat = 0.d0

    smin = 0.D0
    smax = 0.D0
    tsat = 0.d0
    psat = 0.d0

    fluidnr = 1

    phasetype = 0

    !Variables for the enthalpy iteration of the fluid phase at given pressure
    t_min = 0.D0
    t_max = 0.D0
    t_min_allowed = 0.D0
    t_max_allowed = 0.D0
    DeltaS_allowed = 1.d-8
    max_iterations = 200



    if (gl%solidtype(1) == 1) then         !WATER

        gl%solidtype(1) = 1 !Solid water has indicator "1"
        gl%solid_pos = 1 !Water on position 1 in the fluid vector
        gl%solidtype_akt_phase = 1
        gl%solidpos_akt_phase = 1

        p_tmeltmin = 170.184416372448D0                 !Pressure at temperature minimum of Ice Ih. Calculated with the water ice equation
        p_tmeltmax_ICE_VII = 20617.7780087102D0         !Pressure maximum of ICE VII equation (out of bounds of water equation!)

        if(gl%seawater .or. gl%el_present) then

            errval = -12800
            return

            iFlash = 8 !for triple point
            x_known(1) = 1.D0
            x_fluid1 = 0.d0
            x_fluid2 = 0.d0
            x_sol = 0.d0
            x_fluid1(1) = 1.D0
            x_fluid2(1) = 1.D0
            x_hyd = 0.D0
            x_sol(1) = 1.D0

            rhovap_est = 0.D0
            rholiq_est = 0.D0

            iFlash = 8
            iphase = 2
            temptrip = 265.d0
            presstrip = 500.d-6

            call ptflash_solid_NC_3P(gl,presstrip, Temptrip, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhovap_est, &
                & rholiq_est, Phasefrac, iFlash, iphase, iter, errval)

            !gl%seacalc = .false.
            gl%ttp(fluidnr) = temptrip
            gl%ptp(fluidnr) = presstrip


            !call sea_trip(gl, temp, press, d, ttrip, ptrip)
            !gl%ttp(fluidnr) = ttrip
            !gl%ptp(fluidnr) = ptrip
            ! pfreeze =
            ! tfreeze =
        end if

        !Check the boundaries of the solid equation

        smin = s_WaterIce(gl,gl%Tmin_Waterice, press)
        if (s_spec < smin) then
            errval = -9942
            return
        end if

        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        smax = S_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr) !water vapor?!
        if (s_spec > smax) then
            errval = -9942
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature
            psub = press
            Tsub =  tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then
                return
            end if
            !Calculate the sublimation temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Tsub, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in sublimation equilibrium
            rho_s_sub = 1.D0 / v_WaterIce(gl,Tsub, psub)
            s_s_sub = s_WaterIce(gl,Tsub, psub)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_v_sub = rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
            s_v_sub = S_CALC(gl,Tsub, rho_v_sub, fluidnr)

            if (s_spec < s_s_sub) then              !solid water Ih region
                Temp = T_WaterIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (s_spec > s_v_sub) then        !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tsub
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Sublimation region, SVE
                Temp = Tsub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (s_spec-s_s_sub)/(s_v_sub-s_s_sub)
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_sub
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if(gl%seawater) Tmelt = freeze_temp(gl, 250.d0, press, 56000.d0)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if(gl%seawater) Tmelt = freeze_temp(gl, 250.d0, press, 56000.d0)
            if (errval /= 0) return
            if (Tmelt < 1.D-12) errval = -2224

            !Get the density and enthalpy of solid water in melting equilibrium
            rho_s_melt = 1.D0 / v_WaterIce(gl,Tmelt, pmelt)
            s_s_melt = s_WaterIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of liquid in melting equilibrium
            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr) !seawater here
            s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)        !seawater here
            if(gl%seawater) then
                rho_l_melt_w = rho_l_melt
                rho_l_melt = d_sea_calc(gl, Tmelt, pmelt, rho_l_melt)
                s_l_melt = s_sea_calc(gl, Tmelt, pmelt, rho_l_melt_w)
            end if

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            s_l_sat = S_CALC(gl,Tsat, rho_l_sat, fluidnr)   !seawater here
            if(gl%seawater) then
                tsat = boil_temp(gl, temp, press, rho_l_sat)
                rho_l_sat = rhomix_calc(gl, tsat, press, 0.d0, 1, 1)
                rho_v_sat = rhomix_calc(gl, tsat, press, 0.d0, 2, 1)
                s_l_sat = s_sea_calc(gl, tsat, press, rho_l_sat)
                rho_l_sat = d_sea_calc(gl, tsat, press, rho_l_sat)
            end if
            !Get the enthalpy of saturated vapor
            s_v_sat = S_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if (s_spec < s_s_melt) then              !solid water Ih region
                Temp = T_WaterIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((s_spec >= s_s_melt) .and. (s_spec <= s_l_melt)) then       !SLE region
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (s_spec-s_s_melt)/(s_l_melt-s_s_melt)     !seawater here
                rho(phasetype(1)) = rho_l_melt                                       !seawater here
                phasetype(2) = 4                                                     !seawater here
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))             !seawater here
                nrofphases = 2                                                       !seawater here
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif ((s_spec > s_l_melt) .and. (s_spec < s_l_sat)) then        !Liquid region
                !determine temperature from p and h
                !parameters = 0.d0
                parameters%a_p(1) = press                                           !seawater here
                parameters%a_p(2) = s_spec                                          !seawater here
                t_min = Tmelt                                                   !seawater here
                t_max = Tsat                                                    !seawater here
                t_min_allowed = t_min                                           !seawater here
                t_max_allowed = t_max                                           !seawater here
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((s_spec >= s_l_sat) .and. (s_spec <= s_v_sat)) then       !SLE region
                Temp = Tsat                                                         !seawater here
                phasetype(1) = 1                                                    !seawater here
                phasefrac(phasetype(1)) = (s_spec-s_l_sat)/(s_v_sat-s_l_sat)        !seawater here
                rho(phasetype(1)) = rho_v_sat                                       !seawater here
                phasetype(2) = 3                                                    !seawater here
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))            !seawater here
                nrofphases = 2                                                      !seawater here
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif (s_spec > s_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < p_tmeltmin) then        !Region from the maximum pressure where the Ice Ih equation is valid down to the critical pressure

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if(gl%seawater) Tmelt = freeze_temp(gl, tmelt, pmelt, rho_l_melt)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in melting equilibrium
            rho_s_melt = 1.D0 / v_WaterIce(gl,Tmelt, pmelt)
            s_s_melt = s_WaterIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 2, fluidnr)
            s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)
            if(gl%seawater) then
                s_l_melt = s_sea_calc(gl, tmelt, pmelt, rho_l_melt)
                rho_l_melt = d_sea_calc(gl, tmelt, pmelt, rho_l_melt)
            end if

            if (s_spec < s_s_melt) then              !solid water Ih region
                Temp = T_WaterIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_WaterIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (s_spec > s_l_melt) then        !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Melting region, SLE
                Temp = Tmelt                                                !seawater here
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (s_spec-s_s_melt)/(s_l_melt-s_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        elseif (press < p_tmeltmax_ICE_VII) then        !pressure is smaller than the maximum pressure for ICE VII (Fluid EOS anyways only valid up to 1000 MPa, this is already catched in the beginning)

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the melting temperature with the solid water equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of liquid in melting equilibrium
            if(gl%seawater) tmelt = freeze_temp(gl, tmelt, pmelt, rho_l_melt)

            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)   !seawater here
            s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)              !seawater here

            if(gl%seawater) then
                s_l_melt = s_sea_calc(gl, tmelt, pmelt, rho_l_melt)
                rho_l_melt = d_sea_calc(gl, tmelt, pmelt, rho_l_melt)
            end if
            if (s_spec < s_l_melt) then             !solid water region or melting region. No equation available
                errval = -19915
                phasetype(1) = 4
                return
            else                                    !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then                          !Assuming seawater pcrit as water pcrit, as eos not defined at that point!
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        end if



    elseif (gl%solidtype(1) == 2) then             !CO2

        gl%solidtype(1) = 2 !Solid CO2 has indicator "2"
        gl%solid_pos = 1 !CO2 on position 1 in the fluid vector
        gl%solidtype_akt_phase = 2
        gl%solidpos_akt_phase = 1

        !Check the boundaries of the solid equation
        smin = s_DryIce(gl,gl%Tmin_DryIce, press)
        if (s_spec < smin) then
            errval = -19941
            return
        end if

        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        smax = S_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
        if (s_spec > smax) then
            errval = -19941
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature
            psub = press
            Tsub =  tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then
                return
            end if
            !Calculate the sublimation temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 2  !fluid phase is vapor
            call ptflash_solid_NC_2P(gl,psub, Tsub, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid CO2 in sublimation equilibrium
            rho_s_sub = 1.D0 / v_DryIce(gl,Tsub, psub)
            s_s_sub = s_DryIce(gl,Tsub, psub)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_v_sub = rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
            s_v_sub = S_CALC(gl,Tsub, rho_v_sub, fluidnr)

            if (s_spec < s_s_sub) then              !solid region
                Temp = T_DryIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (s_spec > s_v_sub) then        !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tsub
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Sublimation region, SVE
                Temp = Tsub
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (s_spec-s_s_sub)/(s_v_sub-s_s_sub)
                rho(phasetype(1)) = rho_v_sub
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_sub
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the sublimation temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid co2 in melting equilibrium
            rho_s_melt = 1.D0 / v_DryIce(gl,Tmelt, pmelt)
            s_s_melt = s_DryIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of liquid in melting equilibrium
            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
            s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            s_l_sat = S_CALC(gl,Tsat, rho_l_sat, fluidnr)
            !Get the enthalpy of saturated vapor
            s_v_sat = S_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if (s_spec < s_s_melt) then              !solid region
                Temp = T_DryIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((s_spec >= s_s_melt) .and. (s_spec <= s_l_melt)) then       !SLE region
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (s_spec-s_s_melt)/(s_l_melt-s_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif ((s_spec > s_l_melt) .and. (s_spec < s_l_sat)) then        !Liquid region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tmelt
                t_max = Tsat
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((s_spec >= s_l_sat) .and. (s_spec <= s_v_sat)) then       !SLE region
                Temp = Tsat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (s_spec-s_l_sat)/(s_v_sat-s_l_sat)
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif (s_spec > s_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < gl%pmax_Dryice) then        !pressure is smaller than the maximum pressure for the dry ice equation

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            !Calculate the sublimation temperature with the solid CO2 equation
            iflash = 1  !p given, compute T
            iPhase = 1  !fluid phase is liquid
            call ptflash_solid_NC_2P(gl,pmelt, Tmelt, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
            if (errval /= 0) then
                return
            end if

            !Get the density and enthalpy of solid water in melting equilibrium
            rho_s_melt = 1.D0 / v_DryIce(gl,Tmelt, pmelt)
            s_s_melt = s_DryIce(gl,Tmelt, pmelt)
            !Get the density and enthalpy of vapor in sublimation equilibrium
            rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 2, fluidnr)
            s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)

            if (s_spec < s_s_melt) then              !solid region
                Temp = T_DryIce_ps(gl,press, s_spec, errval)
                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 4
                rho(phasetype(1)) = 1.D0 / v_DryIce(gl,Temp, press)
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif (s_spec > s_l_melt) then        !liquid / supercritical region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tmelt
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            else                                    !Melting region, SLE
                Temp = Tmelt
                phasetype(1) = 3
                phasefrac(phasetype(1)) = (s_spec-s_s_melt)/(s_l_melt-s_s_melt)
                rho(phasetype(1)) = rho_l_melt
                phasetype(2) = 4
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_s_melt
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            end if

        else        !Dry Ice equation not valid for that region (500 MPa, approximately Tcrit)

            !Get the enthalpy of liquid CO2 at the limit of validity of the dry ice equation

            !Check the boundaries of the fluid equation
            rhomin = rhomix_calc(gl,gl%tc(fluidnr), press, 0.d0, 0, fluidnr)
            if (rhomax < 1.D-12) then
                errval = -8888
                return
            end if
            smin = S_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
            if (s_spec < smax) then
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            end if

            !determine temperature from p and h
            parameters%a_p(1) = press
            parameters%a_p(2) = s_spec
            t_min = gl%tc(fluidnr)
            t_max = gl%tmaxfluid(fluidnr)
            t_min_allowed = t_min
            t_max_allowed = t_max
            call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                &                       Max_iterations, Iterations, errval, Parameters)

            if (errval /= 0) then
                return
            end if
            !Determine if the phase is liquid or vapor like
            d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
            if (d < gl%rhoc(fluidnr)) then
                phasetype(1) = 1
            else
                phasetype(1) = 3
            end if
            rho(phasetype(1)) = d
            if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                errval = -1235
                return
            end if
            !d = rho(phasetype(1))
            phasefrac(phasetype(1)) = 1.D0
            nrofphases = 1


        end if


    else !No equation for the solid phase of this substance available

        gl%solidtype_akt_phase = 0
        gl%solidpos_akt_phase = 1
        !Check the boundaries of the fluid equation
        rhomax = rhomix_calc(gl,gl%tmaxfluid(fluidnr), press, 0.d0, 2, fluidnr)
        if (rhomax < 1.D-12) then
            errval = -8888
            return
        end if
        smax = S_CALC(gl,gl%tmaxfluid(fluidnr), rhomax, fluidnr)
        if (s_spec > smax) then
            errval = -19941
            return
        end if
        if (press > gl%pmaxfluid(fluidnr)) then
            errval = -9932
            return
        end if

        If (press < gl%ptp(fluidnr)) then          !pressure below triple point pressure

            !Calculate the sublimation temperature
            psub = press
            Tsub =  tsub_eq(gl,psub, fluidnr, errval)
            if (errval /= 0) then       !Do not exit on error. If no sublimation equation exists, assume fluid
                Tsub = 0.D0
                s_v_sub = 0.D0          !Dummy value
            else
                !Get the density and enthalpy of vapor in sublimation equilibrium
                rho_v_sub = rhomix_calc(gl,Tsub, psub, 0.d0, 2, fluidnr)
                s_v_sub = S_CALC(gl,Tsub, rho_v_sub, fluidnr)
            end if



            if ((s_spec < s_v_sub) .and. (Tsub > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            else                                    !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                if (Tsub > 0.D0) then
                    t_min = Tsub
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        elseif (press < gl%pc(fluidnr)) then

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then    !Do not exit on error. If no sublimation equation exists, assume fluid
                Tmelt = 0.D0
                s_l_melt = s_spec-1.D0         !Dummy value
            else
                !Get the density and enthalpy of liquid in melting equilibrium
                rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
                s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)
            end if

            !Calculate the vapor temperature
            psat = press
            iflash = 2  !p given
            call Flash_Pure_PhaseBoundary(gl,psat, Tsat, rho_v_sat, rho_l_sat, iFlash, errval, iter, fluidnr)
            if (errval /= 0) then
                return
            end if

            !Get the enthalpy of saturated liquid
            s_l_sat = S_CALC(gl,Tsat, rho_l_sat, fluidnr)
            !Get the enthalpy of saturated vapor
            s_v_sat = S_CALC(gl,Tsat, rho_v_sat, fluidnr)

            if ((s_spec < s_l_melt) .and. (Tmelt > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            elseif ((s_spec > s_l_melt) .and. (s_spec < s_l_sat)) then        !Liquid region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                if (Tmelt > 0.D0) then
                    t_min = Tmelt
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = Tsat
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 3
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 1, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            elseif ((s_spec >= s_l_sat) .and. (s_spec <= s_v_sat)) then       !SLE region
                Temp = Tsat
                phasetype(1) = 1
                phasefrac(phasetype(1)) = (s_spec-s_l_sat)/(s_v_sat-s_l_sat)
                rho(phasetype(1)) = rho_v_sat
                phasetype(2) = 3
                phasefrac(phasetype(2)) = 1.D0 - phasefrac(phasetype(1))
                nrofphases = 2
                rho(phasetype(2)) = rho_l_sat
                d = 1.D0/(phasefrac(phasetype(1))/rho(phasetype(1))+phasefrac(phasetype(2))/rho(phasetype(2)))
            elseif (s_spec > s_v_sat) then  !Gas region
                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                t_min = Tsat
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                phasetype(1) = 1
                rho(phasetype(1)) = rhomix_calc(gl,temp, press, 0.d0, 2, fluidnr)
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1
            end if

        else        !Region above the critical pressure

            !Calculate the melting temperature
            pmelt = press
            Tmelt = tmelt_eq(gl,pmelt, fluidnr)
            if (Tmelt < 1.D-12) then    !Do not exit on error. If no sublimation equation exists, assume fluid
                Tmelt = 0.D0
                s_l_melt = s_spec-1.D0         !Dummy value
            else
                !Get the density and enthalpy of liquid in melting equilibrium
                rho_l_melt = rhomix_calc(gl,Tmelt, pmelt, 0.d0, 1, fluidnr)
                s_l_melt = S_CALC(gl,Tmelt, rho_l_melt, fluidnr)
            end if


            if ((s_spec < s_l_melt) .and. (Tmelt > 0.D0)) then              !solid region
                errval = -19915
                phasetype(1) = 4
                nrofphases = 1
                return
            else     !liquid or supercritical region

                !determine temperature from p and h
                parameters%a_p(1) = press
                parameters%a_p(2) = s_spec
                if (Tmelt > 0.D0) then
                    t_min = Tmelt
                else
                    t_min = gl%tminfluid(fluidnr)
                end if
                t_max = gl%tmaxfluid(fluidnr)
                t_min_allowed = t_min
                t_max_allowed = t_max
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, errval, Parameters)

                if (errval /= 0) then
                    return
                end if
                !Determine if the phase is liquid or vapor like
                d = rhomix_calc(gl,temp, press, 0.d0, 0, fluidnr)      !Calculate density
                if (d < gl%rhoc(fluidnr)) then
                    phasetype(1) = 1
                else
                    phasetype(1) = 3
                end if
                rho(phasetype(1)) = d
                if (dabs(s_spec - S_CALC(gl,temp, rho(phasetype(1)), fluidnr)) > 1.D-4) then
                    errval = -1235
                    return
                end if
                !d = rho(phasetype(1))
                phasefrac(phasetype(1)) = 1.D0
                nrofphases = 1

            end if


        end if


    end if

    end subroutine
    !************************************************************************************


    !************************************************************************************
    subroutine Flash_Pure_PhaseBoundary_sol (gl,press, Temp, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errval)
    !DEC$ ATTRIBUTES DLLEXPORT :: Flash_Pure_PhaseBoundary_sol
    !************************************************************************************
    !Subroutine for calculation states on the sublimation or melting line for pure substances
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
    !   iFlash      - see above
    !   iPhase      - Fluid phase specifier: 1: Liquid, 2: Vapor
    !
    ! OUTPUT: (depending on chosen iFlash)
    !   Temp        - Temperature [K]
    !   press       - Pressure [MPa]
    !
    !   rho         - density [mol/m³] of all phases
    !   d           - overall density [mol/m³]
    !
    !   The following two variables indicate the phase type. For a pure substance,
    !   a maximum of 3 phases can be in equilibrium
    !   phasetype   - phase identifier
    !   phasefrac   - Molar phase fraction

    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************






    implicit none

    type(type_gl) :: gl


    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    integer:: errval, iFlash, nrofphases, iPhase

    double precision :: psub, Tsub, pmelt, Tmelt
    integer :: fluidnr

    !Variables needed for ptflash_solid_NC_2P
    double precision, dimension(30) :: x_known, x_hyd, x_solid, x_fluid
    double precision:: beta_loc, rhofluid_est
    integer:: iter

    double precision :: Tmin, Tmax, pmin, pmax

    double precision, dimension(3) :: rho_sol

    fluidnr = 1
    errval = 0
    d = 0.d0
    phasetype = 0

    !!WATER
    !!MLT line
    !Tmin = 256.164d0                !K
    !Tmax = ttp(fluidnr)             !K
    !pmax = 170.184416372448D0       !MPa
    !pmin = ptp(1)                   !MPa
    !!SBL line
    !Tmin = 50.D0                    !K
    !Tmax = ttp(fluidnr)             !K
    !pmin = 0.D0                     !MPa
    !pmax = ptp(fluidnr)             !MPa

    !!CO2
    !!MLT line
    !Tmin = ttp(fluidnr)             !K
    !Tmax = tc(fluidnr)              !K
    !pmax = 500.D0                   !MPa
    !pmin = ptp(fluidnr)             !MPa
    !!SBL line
    !Tmin = 80.D0                    !K
    !Tmax = ttp(fluidnr)             !K
    !pmin = 5.65D-012                !MPa
    !pmax = ptp(fluidnr)             !Mpa

    !!OTHER
    !!MLT line
    !Tmin = ttp(fluidnr)
    !Tmax = pmeltmaxtemp(fluidnr)
    !pmax = pmelt_eq(gl,Tmax, fluidnr)
    !pmin = ptp(fluidnr)
    !!SBL line
    !Tmin = 50.D0
    !Tmax = ttp(fluidnr)
    !pmin = pmelt_eq(gl,Tmin, fluidnr)
    !pmax = ptp(fluidnr)


    !check the temperature or pressure inputs vor validity
    if (gl%solidtype(1) == 1) then         !Water
        if (iFlash == 1) then !p given
            if (iphase == 1) then !melting line
                pmin = gl%ptp(fluidnr)
                pmax = gl%pmax_Waterice !170.184416372448D0
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2230
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                pmin = 0.D0
                pmax = gl%ptp(fluidnr)
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2232
                    return
                end if
            else
                errval = -1111
                return
            end if
        elseif (iFlash == 2) then !T given
            if (iphase == 1) then !melting line
                Tmin = 256.164d0
                Tmax = gl%ttp(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2231
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                Tmin = 50.d0
                Tmax = gl%ttp(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2233
                    return
                end if
            else
                errval = -1111
                return
            end if
        else
            errval = -1111
            return
        end if
    elseif (gl%solidtype(1) == 2) then     !CO2
        if (iFlash == 1) then !p given
            if (iphase == 1) then !melting line
                pmin = gl%ptp(fluidnr)
                pmax = 500.D0
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2230
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                pmin = 5.65D-012
                pmax = gl%ptp(fluidnr)
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2232
                    return
                end if
            else
                errval = -1111
                return
            end if
        elseif (iFlash == 2) then !T given
            if (iphase == 1) then !melting line
                Tmin = gl%ttp(fluidnr)
                Tmax = gl%tc(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2231
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                Tmin = 80.d0
                Tmax = gl%ttp(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2233
                    return
                end if
            else
                errval = -1111
                return
            end if
        else
            errval = -1111
            return
        end if
    else                                !other substance
        if (iFlash == 1) then !p given
            if (iphase == 1) then !melting line
                pmin = gl%ptp(fluidnr)
                pmax = pmelt_eq(gl,gl%pmeltmaxtemp(fluidnr), fluidnr)
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2230
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                pmin = 0.D0
                pmax = gl%ptp(fluidnr)
                if ((press < pmin) .or. (press > pmax)) then
                    errval = -2232
                    return
                end if
            else
                errval = -1111
                return
            end if
        elseif (iFlash == 2) then !T given
            if (iphase == 1) then !melting line
                Tmin = gl%ttp(fluidnr)
                Tmax = gl%pmeltmaxtemp(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2231
                    return
                end if
            elseif (iphase == 2) then !sublimation line
                Tmin = 50.d0
                Tmax = gl%ttp(fluidnr)
                if ((Temp < Tmin) .or. (Temp > Tmax)) then
                    errval = -2233
                    return
                end if
            else
                errval = -1111
                return
            end if
        else
            errval = -1111
            return
        end if
    end if



    if (iFlash == 1) then
        !Get an inital estimate for the sublimation or melting temperature
        if (iPhase == 2) then
            Tsub = Tsub_eq(gl,press, fluidnr, errval)
            if (errval /= 0) return
            Temp = Tsub
            phasetype(1) = 1    !First phase is vapor
        elseif (iPhase == 1) then
            Tmelt = Tmelt_eq(gl,press, fluidnr)
            if (Tmelt < 1.D-12) then
                errval = -2224
                return
            end if
            Temp= Tmelt
            phasetype(1) = 3    !First phase is liquid
        else
            errval = -1111 !wrong input
        end if
    elseif (Iflash == 2) then
        !Get an inital estimate for the sublimation or melting pressure
        if (iPhase == 2) then
            psub = psub_eq(gl,Temp, fluidnr)
            if (psub < 1.D-12) then
                errval = -2225
                return
            end if
            press = psub
            phasetype(1) = 1    !First phase is vapor
        elseif (iPhase == 1) then
            pmelt = pmelt_eq(gl,Temp, fluidnr)
            if (pmelt < 1.D-12) then
                errval = -2226
                return
            end if
            press = pmelt
            phasetype(1) = 3    !First phase is liquid
        else
            errval = -1111 !wrong input
        end if
    else
        errval = -1111 !wrong input
    end if

    !Second phase is solid
    phasetype(2) = 4

    if ((gl%solidtype(1) == 1) .or. (gl%solidtype(1) == 2)) then
        !Calculate the equilibrium pressure
        x_known = 0.D0
        x_solid = 0.d0
        x_hyd = 0.D0
        x_fluid = 0.D0
        rhofluid_est = 0.D0
        beta_loc = 0.d0
        x_known(1) = 1.D0
        x_solid(1) = 1.D0
        x_fluid(1) = 1.D0

        call ptflash_solid_NC_2P(gl,press, Temp, rho_sol, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, iPhase, errval, iter)
        if (errval /= 0) then
            return
        end if
    end if

    nrofphases = 2

    !Density of solid phase equilibrium conditions
    if (gl%solidtype(1) == 1) then         !Water
        rho(phasetype(2)) = rho_sol(2)!1.D0 / v_WaterIce(gl,Temp,press)
    elseif (gl%solidtype(1) == 2) then     !CO2
        rho(phasetype(2)) = rho_sol(2)!1.D0 / v_DryIce(gl,Temp,press)
    else                                !No solid equation available
        rho(phasetype(2)) = 0.D0
        !errval = -9904
    end if
    !Density of gaseous water at the sublimation line
    rho(phasetype(1)) = rho_sol(1)!rhomix_calc(gl,Temp, press, 0.d0, iPhase, fluidnr)



    end subroutine
    !************************************************************************************


    !************************************************************************************
    subroutine phasecheck_sea(gl,press, Temp, d, rho, phasetype, phasefrac, nrofphases, errval)



    implicit none
    type(type_gl) :: gl

    double precision :: press
    double precision :: temp
    double precision :: d
    double precision, dimension(5) :: rho
    double precision, dimension(5) :: phasefrac     !phasefrac(1): Phase fraction of phase 1
    !phasefrac(2): Phase fraction of phase 2
    !phasefrac(3): Phase fraction of phase 3
    integer, dimension(5) :: phasetype              !phasetype = 1: vapor
    !phasetype = 2: liquid phase 1
    !phasetype = 3: liquid phase 2
    !phasetype = 4: solid phase
    !phasetype = 5: hydrate phase
    integer :: nrofphases
    integer :: errval
    integer :: fluidnr


    double precision :: gice, gsea, gvap

    fluidnr = 1


    phasefrac(phasetype(1)) = 1.D0
    nrofphases = 1
    if (phasetype(1) == 4) then
        gsea = chempot_sea_calc(gl,temp, press, rhomix_calc(gl, temp, press, 0.d0, 1, fluidnr))
        gl%seacalc = .false.
        gice = g_WaterIce(gl,temp, press)
        gl%seacalc = .true.
        if ( gsea .le. gice ) then
            phasetype(1) = 3
            rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 1, fluidnr) !muss das nicht die Dichte von Seewasser sein?
            d = 0.D0
        elseif( abs(gsea - gice) .le. 1.d-14) then
            phasetype(1) = 3
            phasetype(2) = 4
            !errval = -19915
            !return
        else
            phasetype(1) = 4
        endif

    elseif (phasetype(1) == 1) then

        gsea = chempot_sea_calc(gl,temp, press, rhomix_calc(gl, temp, press, 0.d0, 1, fluidnr))
        gl%seacalc = .false.
        gvap = g_calc(gl,temp, rho(phasetype(1)), 1)
        gl%seacalc = .true.

        if ( gsea .lt. gvap ) then !hardcoded nrsubst=1, because of water which is due to the hydrate always put on the first place
            phasetype(1) = 3
            rho(phasetype(1)) = rhomix_calc(gl,Temp, press, 0.d0, 1, fluidnr)
            d = 0.D0
        elseif(abs(gsea -gvap) .le. 1.d-14) then
            phasetype(1) = 1
            phasetype(2) = 3
            !errval = -19915 !MB muss noch der entsprechende fehlercode f vapor anstatt solid rein
            !return
        else
            phasetype(1) = 1
        endif
    endif

    end subroutine


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    end module phasedet_sol_pure_module
