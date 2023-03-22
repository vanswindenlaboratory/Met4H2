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




    !DEC$ IF DEFINED(HDRT)
    submodule (hdrt_phaselines_calc) impl
    use module_all_types
    use module_regula_falsi
    use rhomix_pt_module
    use flash_module
    use phasedet_mix_module
    use hdrt_chem_pot_module
    use ptflash_solids_mix_module
    use hdrt_properties_modules_module
    use reduced_parameters_calc_module
    use hdrt_main_module
    use ptflash_2c_3p_module
    use waterice_module
    use dryice_module

    use module_hdrt_phaselines
    use omp_lib

    contains
    !******************************************************************************
    module subroutine quadruplepoints (gl, path, Q_calc, Q_point, Qexist, Temp_Q, press_Q, &
        & x, x_ph1, x_ph2, x_ph3, x_ph4)
    !******************************************************************************
    ! Calculation of all existing quadruple points for pure hydrate



    implicit none
    type(type_gl) :: gl
    character(len=255) :: path

    double precision, dimension(30) :: x, lnfi !<-lnfi is currently used in the commented part ONLY
    double precision:: Temp, Dens
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est  !, fug_gas
    integer :: error, i, solidnr, solid !<- solidnr and solid are used in the commented part ONLY
    integer :: iPhase, iFlash, iter, iPhase_sol

    !Pressures and Temperatures at each Quadrupelpoint
    character(len=6), dimension(4)  :: Q_point       !Quadruple point
    double precision, dimension (6) :: press_Q, Temp_Q
    double precision, dimension(4,30) :: x_ph1, x_ph2, x_ph3, x_ph4, x_phX
    logical, dimension (4):: Qexist
    integer, dimension(5) :: error_Q, Q_calc

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double !<- CiJ and occup are used in the commented part ONLY
    double precision, dimension(30):: xH, fug_gas!, occup_ls, occup_ld, occup_sms, occup_smd !<- xH and fug_gas are used in the commented part ONLY

    !New quadruple point routine
    double precision, dimension(4):: rho_qpts, Phasefrac
    double precision, dimension(30):: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    integer:: Phasefrac_0

    !//////////////////////////////////////////////////////////////////

    ! Q_calc defines which Q-point shall be calculated
    ! 1-VLwHIw, 2-VLwLcH, 3-VLxHIx, 4-LwLxHIx, 5-all existing Q-points and stored in txt file

    error = 0
    if ((x(1) < 1.d-12) .and. (x(2) < 1.d-12)) then
        x = 0.D0
        Phasefrac = 0.D0
    else
        phasefrac(1) = 0.4d0
        phasefrac(2) = 0.4d0
        phasefrac(3) = 0.2d0
        phasefrac(4) = 0.0d0
    endif
    x_ph1 = 0.d0
    x_ph2 = 0.d0
    x_ph3 = 0.D0
    x_ph4 = 0.D0
    x_phX = 0.d0        ! only help variable = 0
    beta = 0.D0
    rhoph1_est = 0.D0
    rhoph2_est = 0.D0
    rhoph3_est = 0.D0
    rhophX_est = 0.D0   ! only help variable = 0
    !rhofluid_est = 0.D0
    iFlash = 5
    Qexist = .false.
    Q_point = ''
    press_Q = 0.d0
    Temp_Q = 0.d0
    error_Q = 0

    rhofluid1_est = 0.D0
    rhofluid2_est = 0.D0
    rhofluid3_est = 0.D0
    rhofluid4_est = 0.D0

    !The iPhase_sol can be used to assume fluid phases in two and three phase equilibrium calculation
    ! 0 : Let the algorithms choose the phases
    ! 1 : Liquid phase of the hydrate guest is assumed (e.g. LwH or LwLcH)
    ! 2 : Vapor phase of the hydrate guest is assumed (e.g. VH or VLwH
    !NOTE: IN case of two fluid phases, the heavier phase is always assumed as liquid phase (no vapor / vapor equilibrium)
    iPhase_sol = 0

    ! In the following section quadruple points need to be calculated. Depending on the guest molecule, different quadruple
    ! points will need to be found
    ! In general, all Hydrates should have the Quadruple point: VLwHIw
    ! Depending on the guest, the quadruple points VLwLxH, VLcHIx and LwLcHIx can be calculated.
    ! In this case x stands for the guest molecule and Lx for a liquid phase and Ix for a solid phase
    ! of this guest
    ! The algorithms implemented so far only work for binary mixtures!!!
    ! FURTHERMORE WATER HAS ALWAYS TO BE COMPONENT 1!!!
    !----------------------------------------------------------------------------------------------

    !---------------------------------------------
    if ((Q_calc(1) == 1).or.(Q_calc(5) == 1)) then
        ! VLwHIw - ph1 = vap, ph2 = liq2, ph3 = hyd, ph4 = sol
        !Calculate the quadruple point VLwHIw (Q1)
        Q_point(1) = 'VLwHIw'
        !Choose startvalues depending on the hydrate former
        !Set arbitrary values first
        temp_Q(1) = 272.D0
        press_Q(1) = 5.0D0
        if (gl%Fluidlist_hydrate(2) == "co2") then
            temp_Q(1) = 271.246D0
            press_Q(1) = 1.017105D0
            !if ((any(gl%components == 'ethanol')).or.(any(gl%components == 'methanol'))) press_Q(1) = 2.D0
        elseif (gl%Fluidlist_hydrate(2) == "nitrogen") then
            temp_Q(1) = 272.D0
            press_Q(1) = 15.D0
        elseif (gl%Fluidlist_hydrate(2) == "oxygen") then
            temp_Q(1) = 272.D0
            press_Q(1) = 11.D0
        elseif (gl%Fluidlist_hydrate(2) == "argon") then
            temp_Q(1) = 273.D0
            press_Q(1) = 11.D0
        elseif (gl%Fluidlist_hydrate(2) == "CO") then
            temp_Q(1) = 272.D0
            press_Q(1) = 13.6D0
        elseif (gl%Fluidlist_hydrate(2) == "methane") then
            temp_Q(1) = 272.D0
            press_Q(1) = 2.5D0
        elseif (gl%Fluidlist_hydrate(2) == "ethane") then
            temp_Q(1) = 272.D0
            press_Q(1) = 0.45D0
        elseif (gl%Fluidlist_hydrate(2) == "propane") then
            temp_Q(1) = 272.D0
            press_Q(1) = 0.16D0
        end if
        !
        !!Old quadruple point routine
        !solid = 2
        !solidnr = 1
        !x_ph4(1,1) = 1.d0                 ! x_sol(1) = 1.D0
        !x_ph4(2,1) = 0.d0                 ! x_sol(2) = 0.D0
        !x_ph1(1,1) = 0.00064468D0       ! x_vap_Q(1,1) = 0.00064468D0
        !x_ph1(1,2) = 0.99935532D0       ! x_vap_Q(1,2) = 0.99935532D0
        !x_ph2(1,1) = 0.98468113D0       ! x_liq2_Q(1,1) = 0.98468113D0
        !x_ph2(1,2) = 0.01531887D0       ! x_liq2_Q(1,2) = 0.01531887D0
        !call ptflash_sol_2C_4P(press_Q(1), Temp_Q(1), solid, solidnr, 2, x, x_ph4(1,:), x_ph1(1,:), x_phX(1,:), &
        !                     & x_ph2(1,:), rhoph1_est, rhophX_est, rhoph2_est, error_Q(1), iter)




        iFlash = 1          !p given, Temperature iterated
        iphase = 2          !First fluid phase is vapor
        Phasefrac_0 = 4     !beta_Iw = 0
        gl%solidtype(1) = 1    !Solid water forms
        gl%solidtype(2) = 1    !Hydrate forms
        gl%solid_pos = 1       !Water on position 1
        x_sol(1) = 1.d0                 ! x_sol(1) = 1.D0
        x_sol(2) = 0.d0                 ! x_sol(2) = 0.D0
        x_fluid1(1) = 0.00064468D0       ! x_vap_Q(1,1) = 0.00064468D0
        x_fluid1(2) = 0.99935532D0       ! x_vap_Q(1,2) = 0.99935532D0
        x_fluid2(1) = 0.98468113D0       ! x_liq2_Q(1,1) = 0.98468113D0
        x_fluid2(2) = 0.01531887D0       ! x_liq2_Q(1,2) = 0.01531887D0
        !if ((any(gl%components == 'ethanol')).or.(any(gl%components == 'methanol'))) then
        !    x_fluid1(1) = 0.000551d0!0.9d0 * x_fluid1(1)
        !    x_fluid1(2) = 0.998750d0!0.9d0 * x_fluid1(2)
        !    x_fluid1(3) = 1.d0 - sum(x_fluid1(1:2))
        !    x_fluid2(1) = 0.826458d0!0.9d0 * x_fluid2(1)
        !    x_fluid2(2) = 0.008635633d0!0.9d0 * x_fluid2(2)
        !    x_fluid2(3) = 1.d0 - sum(x_fluid2(1:2))
        !endif
        x_hyd = 0.D0

        !Test of quadruple point calculation for three components
        !x(1) = 0.1D0
        !x(2) = 0.89D0
        !x(3) = 0.01D0
        !x_sol(1) = 1.d0                 ! x_sol(1) = 1.D0
        !x_sol(2) = 0.d0                 ! x_sol(2) = 0.D0
        !x_fluid1(1) = 0.00064468D0       ! x_vap_Q(1,1) = 0.00064468D0
        !x_fluid1(2) = 0.98935532D0       ! x_vap_Q(1,2) = 0.99935532D0
        !x_fluid1(3) = 0.01
        !x_fluid2(1) = 0.98468113D0       ! x_liq2_Q(1,1) = 0.98468113D0
        !x_fluid2(2) = 0.01531787D0       ! x_liq2_Q(1,2) = 0.01531887D0
        !x_fluid2(3) = 0.000001
        !Phasefrac(1) = 0.45D0
        !Phasefrac(2) = 0.45D0
        !Phasefrac(3) = 0.1D0
        !Phasefrac(4) = 0.D0
        !press_Q(1) = 0.99D0
        !Temp_Q(1) = 271.2D0
        call ptflash_solid_NC_4P(gl, press_Q(1), Temp_Q(1), x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
            & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error_Q(1), iter)
        x_ph1(1,:) = x_fluid1
        x_ph2(1,:) = x_fluid2
        x_ph3(1,:) = x_hyd
        x_ph4(1,:) = x_sol

        !Calculate the Chem. Pot and Fugacities of the phases
        if (error_Q(1) == 0) then
            !Needed for old quadruple point routine
            !IPhase = 0
            !molfractions = x_ph1(1,:)
            !call reduced_parameters_calc(gl,Temp_Q(1))
            !Dens = rhomix_calc(gl,Temp_Q(1), press_Q(1), 0.D0, IPhase,0)
            !call lnf_mix(gl,Temp_Q(1), Dens, press_Q(1), lnfi)

            !
            !!Get the composition of the hydrate phase
            !! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
            !! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
            !fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
            !call hdrt_mole_fract(gl, Temp_Q(1),press_Q(1)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            !x_ph3(1,1) = xH(1)
            !x_ph3(1,2) = xH(2)
            !quadruple point exists and was found
            Qexist(1) = .true.
        else
            write(12,*) 'Failed to calculate quadruple point VLwHIw'
            return
        End if

    end if      ! end of VLwHIw

    !---------------------------------------------
    if ((Q_calc(2) == 1).or.(Q_calc(5) == 1)) then
        ! VLwLxH - ph1 = vap, ph2 = liq2, ph3 = liq1, ph4 = hyd
        !Quadruple point for systems where liquid liquid equilibria form
        if ((gl%Fluidlist_hydrate(2) == "co2") .or. (gl%Fluidlist_hydrate(2) == "ethane") .or. (gl%Fluidlist_hydrate(2) == "propane")) then
            !Calculate Q2 (VLwLxH)
            !Vapor Liquid water hydrate and another liquid in equilibrium
            if (gl%Fluidlist_hydrate(2) == "co2") then
                Q_point(2) = 'VLwLcH'
                press_Q(2) = 4.5D0
                temp_Q(2) = 283.0D0
            elseif (gl%Fluidlist_hydrate(2) == "ethane") then
                Q_point(2) = 'VLwLeH'
                press_Q(2) = 3.3D0
                temp_Q(2) = 288.0D0
            elseif (gl%Fluidlist_hydrate(2) == "propane") then
                Q_point(2) = 'VLwLpH'
                press_Q(2) = 0.56D0
                temp_Q(2) = 278.5D0
            end if
            !Old routine to calculate quadruple point
            !solid = 0
            !solidnr = 1
            !x_phX(2,1) = 0.d0     !x_sol(1) = 0.D0
            !x_phX(2,2) = 1.d0     !x_sol(2) = 1.D0
            !x_ph1(2,1) = 0.001D0    !x_vap
            !x_ph1(2,2) = 0.999D0
            !x_ph3(2,1) = 0.002D0    !x_liq1
            !x_ph3(2,2) = 0.998D0
            !x_ph2(2,1) = 0.98D0     !x_liq2
            !x_ph2(2,2) = 0.02D0
            !call ptflash_sol_2C_4P(press_Q(2), Temp_Q(2), solid, solidnr, 2, x, x_phX(2,:), x_ph1(2,:), x_ph3(2,:), &
            !                     & x_ph2(2,:), rhoph1_est, rhoph3_est, rhoph2_est, error_Q(2), iter)


            iFlash = 1          !p given, Temperature iterated
            iphase = 2          !First fluid phase is vapor
            Phasefrac_0 = 3     !beta_Lc = 0
            gl%solidtype(1) = 0    !No pure solid forms
            gl%solidtype(2) = 1    !Hydrate forms
            gl%solid_pos = 1       !Dummy value, no solid forms

            x_fluid1(1) = 0.001D0       ! x_vap_Q(1,1) = 0.00064468D0
            x_fluid1(2) = 0.999D0       ! x_vap_Q(1,2) = 0.99935532D0

            !From Andy (possibly wrong?):
            x_fluid2(1) = 0.002D0       ! x_liq2_Q(1,1) = 0.98468113D0
            x_fluid2(2) = 0.998D0       ! x_liq2_Q(1,2) = 0.01531887D0
            x_fluid3(1) = 0.98D0        ! x_liq2_Q(1,1) = 0.98468113D0
            x_fluid3(2) = 0.02D0        ! x_liq2_Q(1,2) = 0.01531887D0
            !SH
            !x_fluid2(1) = 0.98D0      ! x_liq2_Q(1,1) = 0.98468113D0
            !x_fluid2(2) = 0.02D0      ! x_liq2_Q(1,2) = 0.01531887D0
            !x_fluid3(1) = 0.002D0       ! x_liq2_Q(1,1) = 0.98468113D0
            !x_fluid3(2) = 0.998D0       ! x_liq2_Q(1,2) = 0.01531887D0



            x_hyd = 0.D0
            x_sol = 0.D0
            call ptflash_solid_NC_4P(gl, press_Q(2), Temp_Q(2), x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, &
                & x_hyd, rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, &
                & iphase, Phasefrac_0, error_Q(2), iter)
            x_ph1(2,:) = x_fluid1
            x_ph2(2,:) = x_fluid3
            x_ph3(2,:) = x_fluid2
            !x_ph2(2,:) = x_fluid2
            !x_ph3(2,:) = x_fluid3
            x_ph4(2,:) = x_hyd

            if (error_Q(2) == 0) then
                !Needed for old quadruple point routine
                !Calculate the Chem. Pot and Fugacities of the phases
                !IPhase = 0
                !molfractions = x_ph1(2,:)
                !call reduced_parameters_calc(gl,Temp_Q(2))
                !Dens = rhomix_calc(gl,Temp_Q(2), press_Q(2), 0.D0, IPhase,0)
                !call lnf_mix(gl,Temp_Q(2), Dens, press_Q(2), lnfi)

                !
                !!Get the composition of the hydrate phase
                !! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
                !! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
                !fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
                !call hdrt_mole_fract(gl, Temp_Q(2),press_Q(2)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
                !x_ph4(2,1) = xH(1)
                !x_ph4(2,2) = xH(2)
                !quadruple point exists and was found
                Qexist(2) = .true.
            else
                write(12,*) 'Failed to calculate quadruple point VLwLxH'
                return
            End if
        End if
    end if      ! end of VLwLxH

    !---------------------------------------------
    if ((Q_calc(3) == 1).or.(Q_calc(4) == 1).or.(Q_calc(5) == 1)) then
        ! VLxHIx - ph1 = vap, ph2 = liq1, ph3 = hyd, ph4 = sol
        !Quadruple point with hydrate and solid guest
        if (gl%Fluidlist_hydrate(2) == "co2") then
            !Calculate Q3 (VLxHIx)
            !Vapor, liquid guest, Hydrate and solid guest in equilibrium
            if (gl%Fluidlist_hydrate(2) == "co2") then
                Q_point(3) = 'VLcHIc'
                press_Q(3) = 0.51791D0
                temp_Q(3) = 216.59D0
            End if
            !Old quadruple point routine
            !solid = 1
            !solidnr = 2
            !x_ph4(3,1) = 0.d0           !x_sol(1) = 0.D0
            !x_ph4(3,2) = 1.d0           !x_sol(2) = 1.D0
            !x_ph1(3,1) = 0.00000358D0   !x_vap_Q(3,1) = 0.00000358D0
            !x_ph1(3,2) = 0.99999642D0   !x_vap_Q(3,2) = 0.99999642D0
            !x_ph2(3,1) = 0.000052D0     !x_liq1_Q(3,1) = 0.000052D0
            !x_ph2(3,2) = 0.999948D0     !x_liq1_Q(3,2) = 0.999948D0
            !x_phX(3,1) = 0.d0           !x_liq2(1) = 0.D0
            !x_phX(3,2) = 0.d0           !x_liq2(2) = 0.D0
            !call ptflash_sol_2C_4P(press_Q(3), Temp_Q(3), solid, solidnr, 2, x, x_ph4(3,:), x_ph1(3,:), x_ph2(3,:), &
            !                     & x_phX(3,:), rhoph1_est, rhoph2_est, rhophX_est, error_Q(3), iter)



            iFlash = 1          !p given, Temperature iterated
            iphase = 2          !First fluid phase is vapor
            Phasefrac_0 = 4     !beta_Ic = 0
            gl%solidtype(1) = 2    !Dry Ice forms
            gl%solidtype(2) = 1    !Hydrate forms
            gl%solid_pos = 2       !Position of CO2 in fluid vector
            x_fluid1(1) = 0.00000358D0       ! x_vap_Q(1,1) = 0.00064468D0
            x_fluid1(2) = 0.99999642D0       ! x_vap_Q(1,2) = 0.99935532D0
            x_fluid2(1) = 0.000052D0       ! x_liq2_Q(1,1) = 0.98468113D0
            x_fluid2(2) = 0.999948D0       ! x_liq2_Q(1,2) = 0.01531887D0
            x_hyd = 0.D0
            x_sol(1) = 0.D0
            x_sol(2) = 1.D0
            call ptflash_solid_NC_4P(gl,press_Q(3), Temp_Q(3), x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, &
                & x_hyd, rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, &
                & iphase, Phasefrac_0, error_Q(3), iter)
            x_ph1(3,:) = x_fluid1
            x_ph2(3,:) = x_fluid2
            x_ph3(3,:) = x_hyd
            x_ph4(3,:) = x_sol

            If (error_Q(3) == 0) then
                !Needed for old quadruple point routine
                !Calculate the Chem. Pot and Fugacities of the phases
                !IPhase = 0
                !molfractions = x_ph1(3,:)
                !call reduced_parameters_calc(gl,Temp_Q(3))
                !Dens = rhomix_calc(gl,Temp_Q(3), press_Q(3), 0.D0, IPhase,0)
                !call lnf_mix(gl,Temp_Q(3), Dens, press_Q(3), lnfi)

                !
                !!Get the composition of the hydrate phase
                !! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
                !! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
                !fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
                !call hdrt_mole_fract(gl, Temp_Q(3),press_Q(3)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
                !x_ph3(3,1) = xH(1)
                !x_ph3(3,2) = xH(2)
                !quadruple point exists and was found
                Qexist(3) = .true.
            else
                write(12,*) 'Failed to calculate quadruple point VLxHIx'
                return
            End if
        end if
    end if      ! end of VLxHIx

    !---------------------------------------------
    if ((Q_calc(4) == 1).or.(Q_calc(5) == 1)) then
        ! LwLxHIx - ph1 = liq2, ph2 = liq1, ph3 = hyd, ph4 = sol
        if (gl%Fluidlist_hydrate(2) == "co2") then
            !Calculate Q4 (LwLxHIx)
            !Liquid water, liquid guest, Hydrate and solid guest in equilibrium
            if (gl%Fluidlist_hydrate(2) == "co2") then
                Q_point(4) = 'LwLcHIc'
                press_Q(4) = 450.D0
                temp_Q(4) = 280.D0
            end if

            !Old quadruple point calculation routine
            !solid = 1
            !solidnr = 2
            !x_ph4(4,1) = 0.D0       !x_sol(1) = 0.D0
            !x_ph4(4,2) = 1.D0       !x_sol(2) = 1.D0
            !x_ph2(4,1) = 0.005D0    !x_liq1_Q(4,1) = 0.005D0
            !x_ph2(4,2) = 0.995D0    !x_liq1_Q(4,2) = 0.995D0
            !x_ph1(4,1) = 0.950D0    !x_liq2_Q(4,1) = 0.950D0
            !x_ph1(4,2) = 0.050D0    !x_liq2_Q(4,2) = 0.050D0

            !x_phX(4,1) = 0.d0       !x_vap(1) = 0.D0
            !x_phX(4,2) = 0.d0       !x_vap(2) = 0.D0
            !call ptflash_sol_2C_4P(press_Q(4), Temp_Q(4), solid, solidnr, 2, x, x_ph4(4,:), x_ph1(4,:), x_ph2(4,:), &
            !                     & x_phX(4,:), rhophX_est, rhoph2_est, rhoph1_est, error_Q(4), iter)


            iFlash = 1          !p given, Temperature iterated
            iphase = 1          !First fluid phase is liquid
            Phasefrac_0 = 4     !beta_Ic = 0
            gl%solidtype(1) = 2    !Dry Ice forms
            gl%solidtype(2) = 1    !Hydrate forms
            gl%solid_pos = 2       !Position of CO2 in fluid vector
            x_fluid1(1) = 0.005D0       ! x_vap_Q(1,1) = 0.00064468D0
            x_fluid1(2) = 0.995D0       ! x_vap_Q(1,2) = 0.99935532D0
            x_fluid2(1) = 0.950D0       ! x_liq2_Q(1,1) = 0.98468113D0
            x_fluid2(2) = 0.050D0       ! x_liq2_Q(1,2) = 0.01531887D0
            x_hyd = 0.D0
            x_sol(1) = 0.D0
            x_sol(2) = 1.D0
            call ptflash_solid_NC_4P(gl,press_Q(4), Temp_Q(4), x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, &
                & x_hyd, rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, &
                & iphase, Phasefrac_0, error_Q(4), iter)
            x_ph1(4,:) = x_fluid2
            x_ph2(4,:) = x_fluid1
            x_ph3(4,:) = x_hyd
            x_ph4(4,:) = x_sol

            If (error_Q(4) == 0) then
                !Needed for old quadruple point calculation routine
                !!Calculate the Chem. Pot and Fugacities of the phases
                !IPhase = 0
                !molfractions = x_ph2(4,:)
                !call reduced_parameters_calc(gl,Temp_Q(4))
                !Dens = rhomix_calc(gl,Temp_Q(4), press_Q(4), 0.D0, IPhase,0)
                !call lnf_mix(gl,Temp_Q(4), Dens, press_Q(4), lnfi)

                !
                !!Get the composition of the hydrate phase
                !! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
                !! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
                !fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
                !call hdrt_mole_fract(gl, Temp_Q(4),press_Q(4)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
                !x_ph3(4,1) = xH(1)
                !x_ph3(4,2) = xH(2)
                !quadruple point exists and was found
                Qexist(4) = .true.
            else
                write(12,*) 'Failed to calculate quadruple point LwLxHIx'
                return
            End if
        End if
    end if      ! end of LwLxHIx

    !----------------------------------------------------------------------------------------------
    ! If all existing quadruple points are calculated,
    ! the subroutine writes them in the Q-points.txt file
    if (Q_calc(5) == 1) then

3024    format (a6,' ',f8.4,' ',f10.6,' ',f16.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
            &,' ',f12.10,' ',f12.10)

        open(unit=12, file=trim(path) // 'Q-points.txt', status='unknown', action='write', iostat=error)

        do i = 1,4
            if (Qexist(i)) then
                write(12,3024) Q_point(i), Temp_Q(i), press_Q(i), x_ph1(i,1), x_ph1(i,2), x_ph2(i,1), x_ph2(i,2), &
                    & x_ph3(i,1), x_ph3(i,2), x_ph4(i,1), x_ph4(i,2)
            end if
        end do

        close(12)
    end if
    !-------------------------------------------------------------------------------------------

    end subroutine quadruplepoints
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    module subroutine threephase_equil(gl, pel, Eqtype, cur_3ph,iFlash,Temp, press, x , x_ph1, x_ph2, x_ph3, &
        &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
    !&   Q_point, Qexist,Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    !******************************************************************************
    ! Calculation of a single three-phase equilibrium point by calling 'ptflash_solid_NC_3P'
    ! iTp = 1 ... temperature given, iTp = 2 ... pressure given




    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    character(6), intent(in) :: Eqtype        !Type of Equilibrium, which should be calculated
    integer, intent(in) :: iflash

    !Variables at the three-phase equilibrium
    double precision :: Temp, press
    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3, x_phX
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: Dens
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est
    double precision, dimension(5):: rho
    double precision, dimension(3):: rho_sol
    integer :: error
    integer :: iPhase, iter, iPhase_sol!, iFlash
    integer :: cur_3ph

    !Variables at the Quadruple point
    !character(len=6), dimension(4)  :: Q_point
    !double precision, dimension (4) :: press_Q, Temp_Q
    !double precision, dimension(4,30) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    !logical, dimension (4):: Qexist

    double precision, dimension(3) :: Phasefrac

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(30):: xH, fug_gas

    !//////////////////////////////////////////////////////////////////


    error = 0
    !x = 0.D0
    beta = 0.D0
    rhoph1_est = rho_ph1
    rhoph2_est = rho_ph2
    rhoph3_est = rho_ph3
    rhophX_est = 0.D0   ! only help variable = 0
    x_phX = 0.d0
    Phasefrac = 0.D0

    CiJ = 0.d0
    occup = 0.d0
    xH = 0.d0
    fug_gas = 0.d0

    !pel%phl4%Q_point_mix(:,pel%phl4%Q_map(EqTypeStart, 1))
    !pel%phl4%Qexist_mix(:,pel%phl4%Q_map(EqTypeStart, 1))
    !!pel%phl4%Temp_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1))
    !!pel%phl4%press_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1))
    !pel%phl4%x_Q_ph1_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1))
    !pel%phl4%x_Q_ph2_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1))
    !pel%phl4%x_Q_ph3_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1))
    !pel%phl4%x_Q_ph4_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)))

    select case (Eqtype)
    case ('VLwH')   ! ph1 = vap, ph2 = liq2, ph3 = hyd (phX = sol)

        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(1)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(1)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph2 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(1,:)    !vap  (Q = VLwHIw)
            !x_ph2 = x_Q_ph2(1,:)    !liq2 (Q = VLwHIw)
            x_phX = 0.D0
        end if

        gl%solidtype(1) = 0 !No pure solid
        gl%solidtype(2) = 1 !Hydrate forms

        iPhase_sol = 2 !Lighter phase is vapor
        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_ph2, x_phX, x_ph3, rhoph1_est, &
            & rhoph2_est, Phasefrac, iFlash, iPhase_sol, iter, error)
        if (error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph3 = xH  ! Mixed hydrate

            rho_ph1 = gl%rho_vap
            rho_ph2 = gl%rho_liq2
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph3)
        end if


    case ('VHIw')   ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq2)
        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(1)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(1)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(1,:)    !vap (Q = VLwHIw)
        end if
        x_ph3 = 0.D0        ! water ice
        x_ph3(1) = 1.D0     ! water ice

        gl%solidtype(1) = 1 !Solid water forms
        gl%solidtype(2) = 1 !Hydrate forms
        gl%solid_pos = 1    !Position of the pure solid former in the fluid vector (Water = 1)

        iPhase_sol = 2 !Liqhter phase is vapor
        !iflash = 1
        Phasefrac = 0.5d0
        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_phX, x_ph3, x_ph2, rhoph1_est, &
            & rhophX_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph2 = xH  ! Mixed hydrate

            rho_ph1 = gl%rho_vap
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph2)
            rho_ph3 = 1.D0 / dgdp_WaterIce(gl,Temp,press)
        end if


    case ('LwHIw')  ! ph1 = liq2, ph2 = hyd, ph3 = sol (phX = liq1)

        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(1)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(1)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph2(1,:)    !liq2 (Q = VLwHIw)
        end if
        x_ph3 = 0.D0
        x_ph3(1) = 1.D0

        gl%solidtype(1) = 1 !Solid water forms
        gl%solidtype(2) = 1 !Hydrate forms
        gl%solid_pos = 1    !Position of the pure solid former in the fluid vector (Water = 1)

        iPhase_sol = 1 !Liqhter (only fluid) phase is liquid

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_phX, x_ph3, x_ph2, rhoph1_est, &
            & rhophX_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 1
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph2 = xH  ! Mixed hydrate

            rho_ph1 = Dens
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph2)
            rho_ph3 = 1.D0 / dgdp_WaterIce(gl,Temp,press)
        end if


    case ('LwLxH') ! ph1 = liq1, ph2 = liq2, ph3 = hyd (phX = sol) ... it is LxLwH in fact
        if (pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1) < 0) return
        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(2)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(2)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph2 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph3_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph3(2,:)    ! liq1 (Q = VLwLxH)
            !x_ph2 = x_Q_ph2(2,:)    ! liq2 (Q = VLwLxH)
        End if

        x_phX = 0.D0

        gl%solidtype(1) = 0 !No pure solid
        gl%solidtype(2) = 1 !Hydrate forms

        iPhase_sol = 1 !Liqhter phase is liquid

        !iflash = 3 !betaLw = 0 is calculated (phasefraction of first fluid phase = 0)

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_ph2, x_phX, x_ph3, rhoph1_est, &
            & rhoph2_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 1
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase, 0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph3 = xH  ! Mixed hydrate

            rho_ph1 = gl%rho_vap
            rho_ph2 = gl%rho_liq2
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph3)
        end if


    case ('VLxH')   ! ph1 = vap, ph2 = liq1, ph3 = hyd (phX = sol or liq2)

        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(2)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(2)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph2 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph3_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(2,:)      ! vap  (Q = VLwLxH)
            !x_ph2 = x_Q_ph3(2,:)      ! liq1 (Q = VLwLxH)
        End if

        x_phX = 0.D0

        gl%solidtype(1) = 0 !No pure solid
        gl%solidtype(2) = 1 !Hydrate forms

        !iFlash = 3 !betaLx = 0 is calculated (phasefraction of second fluid phase = 0)

        iPhase_sol = 2 !Liqhter phase is vapor

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_ph2, x_phX, x_ph3, rhoph1_est, &
            & rhophX_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if (error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase, 0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph3 = xH

            rho_ph1 = gl%rho_vap
            rho_ph2 = gl%rho_liq2
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph3)
        end if


    case ('VHIx')   ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq1)

        ! Mixed hydrates 2015
        if (gl%Fluidlist_hydrate(2) /= 'co2') then
            !write(*,*) 'Error - equilibrium with DRY ICE => CO2 must be the 2nd component!'
            call messages(gl, pel,4, trim(''), trim(''),0,0,0)
        end if

        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(3)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(3)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(3,:)    ! vap  (Q = VLxHIx)
        End if
        x_ph3 = 0.D0
        x_ph3(2) = 1.D0  ! CO2 as 2nd component!

        gl%solidtype(1) = 2 !Pure solid forms (at the moment CO2)
        gl%solidtype(2) = 1 !Hydrate forms
        gl%solid_pos = 2    !Position of the pure solid former in the fluid vector (CO2 = 2)

        iPhase_sol = 2 !Liqhter phase is vapor

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_phX, x_ph3, x_ph2, rhoph1_est, &
            & rhophX_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase, 0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph2 = xH

            rho_ph1 = gl%rho_vap
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph2)
            if (gl%Fluidlist_hydrate(2) == 'co2') then
                rho_ph3 = 1.d0/dgdp_DryIce(gl,Temp,press) ! only for CO2
            else
                rho_ph3 = -9999.d0
            end if

        end if


    case ('LwHIx') ! ph1 = liq2, ph2 = hyd, ph3 = sol (phX = vap)

        ! Mixed hydrates 2015
        if (gl%Fluidlist_hydrate(2) /= 'co2') then
            !write(*,*) 'Error - equilibrium with DRY ICE => CO2 must be the 2nd component!'
            call messages(gl, pel,4, trim(''), trim(''),0,0,0)
        end if

        !Starting point is Quadruple Point 4
        if ((pel%phl4(4)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(4)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(4,:)  ! liq2  (Q = LwLxHIx)
        end if
        x_ph3 = 0.D0
        x_ph3(2) = 1.D0

        gl%solidtype(1) = 2 !Pure solid forms (at the moment CO2)
        gl%solidtype(2) = 1 !Hydrate forms
        gl%solid_pos = 2    !Position of the pure solid former in the fluid vector (CO2 = 2)

        iPhase_sol = 2 !Liqhter (only fluid) phase is liquid

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_phX, x_ph3, x_ph2, rhoph1_est, &
            & rhoph1_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 1
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase, 0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph2(1) = xH(1)
            x_ph2(2) = xH(2)

            rho_ph1 = Dens
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph2)
            if (gl%Fluidlist_hydrate(2) == 'co2') then
                rho_ph3 = 1.d0/dgdp_DryIce(gl,Temp,press) ! only for CO2
            else
                rho_ph3 = -9999.d0
            end if
        end if


    case ('LxHIx')  ! ph1 = liq1, ph2 = hyd, ph3 = sol (phX = vap)

        ! Mixed hydrates 2015
        if (gl%Fluidlist_hydrate(2) /= 'co2') then
            !write(*,*) 'Error - equilibrium with DRY ICE => CO2 must be the 2nd component!'
            call messages(gl, pel,4, trim(''), trim(''),0,0,0)
        end if

        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(3)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(3)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph2(3,:)  ! liq1  (Q = VLxHIx)
        End if
        x_ph3 = 0.D0
        x_ph3(2) = 1.D0

        gl%solidtype(1) = 2 !Pure solid forms (at the moment CO2)
        gl%solidtype(2) = 1 !Hydrate forms
        gl%solid_pos = 2    !Position of the pure solid former in the fluid vector (CO2 = 2)

        iPhase_sol = 1 ! 2: density is from vapor phase?

        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_phX, x_ph3, x_ph2, rhoph1_est, &
            & rhophX_est, Phasefrac, iFlash, iPhase_sol, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 1
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            !Get the composition of the hydrate phase
            fug_gas(1:gl%nrofhydrateformers-1) = dexp(lnfi(2:gl%nrofhydrateformers))*1.d6
            call hdrt_mole_fract(gl, temp,press*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
            x_ph2(1) = xH(1)
            x_ph2(2) = xH(2)

            rho_ph1 = Dens
            call hdrt_density(gl, temp,press*1.d6,fug_gas,occup,CiJ,rho_ph2)
            if (gl%Fluidlist_hydrate(2) == 'co2') then
                rho_ph3 = 1.d0/dgdp_DryIce(gl,Temp,press) ! only for CO2
            else
                rho_ph3 = -9999.d0
            end if
        end if


    case ('VLxLw')  ! ph1 = vap, ph2 = liq1, ph3 = liq2 (phX = vap)


        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(2)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(2)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph2 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph3_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph3 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(2,:)        ! vap (Q = VLwLxH)
            !x_ph2 = x_Q_ph3(2,:)        ! liq1
            !x_ph3 = x_Q_ph2(2,:)        ! liq2
        End if

        !UNCOMMENT THE FOLLOWING SECTION TO CALCULATE VLwLx INDEPENDENTLY OF QUADRUPLE POINT VLxLwH
        !-------------------------------------------------------
        !press = 3.D0 !CO2=3.5;C2H6=3;C3H8=0.4;
        !x_liq2(1) = 0.99999D0
        !x_liq2(2) = 0.00001D0
        !x_liq1(1) = 0.01D0
        !x_liq1(2) = 0.99D0
        !x_vap(1) = 0.001D0
        !x_vap(2) = 0.999D0
        !Temp = 278.D0     !CO2=283.D0;C2H6=;C3H8=278.D90
        !--------------------------------------------------------

        call ptflash_NC_3P(gl,press, Temp, x, rho, x_ph1, x_ph2, x_ph3, rhoph1_est, &
            & rhoph2_est, rhoph3_est, Phasefrac, iFlash, iter, error)

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rho(1)
            call Chempot_Calc(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            rho_ph1 = rho(1)
            rho_ph2 = rho(2)
            rho_ph3 = rho(3)
            occup = 5.555
        end if


    case ('VLwIw')  ! ph1 = vap, ph2 = liq2, ph3 = sol (phX = liq1)
        !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
        if ((pel%phl4(1)%Qexist_mix(pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1)) ).and.(dabs(x_ph1(1))<1.d-12)) then
            !if ((Qexist(1)).and.(dabs(x_ph1(1))<1.d-12)) then
            x_ph1 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            x_ph2 = pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(pel%phl3(cur_3ph)%EqTypeStart)%Q_map(1))
            !x_ph1 = x_Q_ph1(1,:)        ! vap (Q = VLwHIw)
            !x_ph2 = x_Q_ph2(1,:)        ! liq2
        End if
        x_ph3 = 0.D0
        x_ph3(1) = 1.D0

        !UNCOMMENT THE FOLLOWING SECTION TO CALCULATE VLwIw INDEPENDENTLY OF QUADRUPLE POINT VLwHIw
        !-------------------------------------------------------
        !min_press = 1.D0 !CO2=1.5;N2=;CO=20;Ar=12;O2=20;CH4=3;C2H6=1;C3H8=0.25;
        !x_liq2(1) = 0.9999976197D0
        !x_liq2(2) = 0.0000023803D0
        !x_vap(1) = 0.05D0
        !x_vap(2) = 0.95D0
        !Temp = 273.1597D0
        !press = 0.007D0
        !--------------------------------------------------------

        gl%solidtype(1) = 1 !Pure solid forms (at the moment H2O)
        gl%solidtype(2) = 0 !Hydrate forms
        gl%solid_pos = 1    !Position of the pure solid former in the fluid vector (H2O = 1)

        iPhase_sol = 2 !Liqhter (only fluid) phase is liquid
        call ptflash_solid_NC_3P(gl,press, Temp, x, rho_sol, x_ph1, x_ph2, x_ph3, x_phX, rhoph1_est, &
            & rhoph2_est, Phasefrac, iFlash, iPhase_sol, iter, error)
        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            IPhase = 2
            gl%molfractions = x_ph1
            call reduced_parameters_calc(gl,Temp)
            Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_Calc(gl,Temp,Dens, Chempot, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi)

            rho_ph1 = gl%rho_vap
            rho_ph2 = gl%rho_liq2
            rho_ph3 = 1.D0 / dgdp_WaterIce(gl,Temp,press)
            occup = 5.555
        end if

    end select
    gl%molfractions = x
    end subroutine threephase_equil
    !******************************************************************************
    !******************************************************************************

    module subroutine phaselines_init(phl3)

    !use module_hdrt_phaselines

    implicit none

    type(type_hdrt_3ph_lines) :: phl3


    !pel%phl3(cur_3ph)%np_3ph = pel%phl3(cur_3ph)%loops + 1
    phl3%np_3ph = phl3%loops + 1



    !in three component system with water-guest1-guest2 exist threephase lines with different phasefraction = 0
    !for a fixed overall composition there are two threephase lines the calculation is possible
    !also startvalues for press, Temp and compositions are set
    !allocate arrays for values of threephase lines
    !if (.not. allocated(pel%phl3(cur_3ph)%press_tr_mix) ) allocate(pel%phl3(cur_3ph)%press_tr_mix(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%temp_tr_mix) ) allocate(pel%phl3(cur_3ph)%temp_tr_mix(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%chempot) ) allocate(pel%phl3(cur_3ph)%chempot(pel%phl3(cur_3ph)%np_3ph,30))
    !if (.not. allocated(pel%phl3(cur_3ph)%lnfi) ) allocate(pel%phl3(cur_3ph)%lnfi(pel%phl3(cur_3ph)%np_3ph,30))
    !if (.not. allocated(pel%phl3(cur_3ph)%x_ph) ) allocate(pel%phl3(cur_3ph)%x_ph(pel%phl3(cur_3ph)%np_3ph,3,30))
    !if (.not. allocated(pel%phl3(cur_3ph)%rho_tr_mix) ) allocate(pel%phl3(cur_3ph)%rho_tr_mix(pel%phl3(cur_3ph)%np_3ph,3))
    !if (.not. allocated(pel%phl3(cur_3ph)%occup) ) allocate(pel%phl3(cur_3ph)%occup(pel%phl3(cur_3ph)%np_3ph,3,30))
    !if (.not. allocated(pel%phl3(cur_3ph)%hyd_nr) ) allocate(pel%phl3(cur_3ph)%hyd_nr(pel%phl3(cur_3ph)%np_3ph,30))
    !if (.not. allocated(pel%phl3(cur_3ph)%s_hyd) ) allocate(pel%phl3(cur_3ph)%s_hyd(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%h_hyd) ) allocate(pel%phl3(cur_3ph)%h_hyd(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%h_melt) ) allocate(pel%phl3(cur_3ph)%h_melt(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%tpd_sol) ) allocate(pel%phl3(cur_3ph)%tpd_sol(pel%phl3(cur_3ph)%np_3ph))
    !if (.not. allocated(pel%phl3(cur_3ph)%a0_other_struct) ) allocate(pel%phl3(cur_3ph)%a0_other_struct(pel%phl3(cur_3ph)%np_3ph))
    if (.not. allocated(phl3%press_tr_mix) ) allocate(phl3%press_tr_mix(phl3%np_3ph))
    if (.not. allocated(phl3%temp_tr_mix) ) allocate(phl3%temp_tr_mix(phl3%np_3ph))
    if (.not. allocated(phl3%chempot) ) allocate(phl3%chempot(phl3%np_3ph,30))
    if (.not. allocated(phl3%lnfi) ) allocate(phl3%lnfi(phl3%np_3ph,30))
    if (.not. allocated(phl3%x_ph) ) allocate(phl3%x_ph(phl3%np_3ph,3,30))
    if (.not. allocated(phl3%rho_tr_mix) ) allocate(phl3%rho_tr_mix(phl3%np_3ph,3))
    if (.not. allocated(phl3%occup) ) allocate(phl3%occup(phl3%np_3ph,3,30))
    if (.not. allocated(phl3%hyd_nr) ) allocate(phl3%hyd_nr(phl3%np_3ph,30))
    if (.not. allocated(phl3%s_hyd) ) allocate(phl3%s_hyd(phl3%np_3ph))
    if (.not. allocated(phl3%h_hyd) ) allocate(phl3%h_hyd(phl3%np_3ph))
    if (.not. allocated(phl3%h_melt) ) allocate(phl3%h_melt(phl3%np_3ph))
    if (.not. allocated(phl3%tpd_sol) ) allocate(phl3%tpd_sol(phl3%np_3ph))
    if (.not. allocated(phl3%a0_other_struct) ) allocate(phl3%a0_other_struct(phl3%np_3ph))
    phl3%press_tr_mix = 0.d0
    phl3%temp_tr_mix = 0.d0
    phl3%rho_tr_mix = 0.d0
    phl3%pathout = ""
    phl3%x_ph = 0.d0

    end subroutine phaselines_init

    module subroutine phaselines_out(gl, pel, phaseline, pathIn, errorflag)

    !use module_all_types
    !use module_hdrt_phaselines

    implicit none

    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel

    integer :: threephlineunit, fourphlineunit, eqtype_ind, phaseline, cur_3ph
    integer :: i, k, errorflag, row
    character(255) :: path, pathIn

    if (phaseline == 3) then
        path = trim(pathIn) // trim(pel%CompsHdrtsAll) // '\three_phase_lines\'
        do cur_3ph = 1,pel%n_3ph_lines * pel%phl3(1)%nflashs
            if (count(pel%phl3(cur_3ph)%press_tr_mix(:) > 1.d-14) == 0) cycle
            pel%phl3(cur_3ph)%pathout = trim(path) // trim(trim(pel%phl3(cur_3ph)%trline_mix))
            pel%phl3(cur_3ph)%pathout = trim(pel%phl3(cur_3ph)%pathout) // '_z='
            do i = 1,gl%ncomp-1
                write(pel%phl3(cur_3ph)%pathout,'(A,F5.3,A)')trim(pel%phl3(cur_3ph)%pathout),gl%molfractions(i),trim('-')
            enddo
            if (gl%n_guests > 1 .or. any(gl%components(:) == 'ethanol')) then
                write(pel%phl3(cur_3ph)%pathout,'(A,F5.3,A,A)')trim(pel%phl3(cur_3ph)%pathout),gl%molfractions(gl%ncomp),trim(pel%phl3(cur_3ph)%strFlash),trim('.txt')
            else
                write(pel%phl3(cur_3ph)%pathout,'(A,F5.3,A,A)')trim(pel%phl3(cur_3ph)%pathout),gl%molfractions(gl%ncomp),trim('.txt')
            endif
            if (pel%phl3(cur_3ph)%pathout /= "") then
                open(newunit=threephlineunit, file=trim(pel%phl3(cur_3ph)%pathout), status='unknown', action='write', iostat=errorflag)
                if (errorflag /= 0) then
                    call messages(gl, pel,7, trim(''), trim(''),0,0,0)
                    cycle
                endif

                if (pel%print_head .eqv. .true.) then
                    if (gl%N_guests == 1) then
                        if ((pel%phl3(cur_3ph)%Eqtype == 'VLwH').or.(pel%phl3(cur_3ph)%Eqtype == 'VHIw')) then
                            write(threephlineunit,3011)'Eqtype','Temp','press','Chempot1','Chempot2','lnfi1','lnfi2','xph11','xph12','xph21','xph22','xph31','xph32','rhoph1','rhoph2','rhoph3','occupS','occupL','hydnr','hmelt','s_hyd','h_hyd','tpd_other','latt_par'
                        else
                            write(threephlineunit,3012)'Eqtype','Temp','press','Chempot1','Chempot2','lnfi1','lnfi2','xph11','xph12','xph21','xph22','xph31','xph32','rhoph1','rhoph2','rhoph3','occupS','occupL','hydnr','s_hyd','h_hyd','tpd_other','latt_par'
                        endif
                    elseif (gl%N_guests == 2) then
                        write(threephlineunit,3013)'Eqtype','Temp','press','Chempot1','Chempot2','Chempot3','lnfi1','lnfi2','lnfi3','xph11','xph12','xph13','xph21','xph22','xph23','xph31','xph32','xph33', &
                            & 'rhoph1','rhoph2','rhoph3','occupS1','occupS2','occupL1','occupL2','hydnr1','hydnr2','s_hyd','h_hyd'
                    endif
                endif
                !close VLxLw phase envelope
                if (gl%n_guests == 2) then
                    if ((any(gl%components == 'co2')) .and. (any(gl%components == 'methane'))) then
                        if ((pel%phl3(cur_3ph)%trline_mix(1:5) == 'VLcLw').and.(2-mod(cur_3ph,2) == 2)) then
                            row = maxloc(pel%phl3(cur_3ph-1)%press_tr_mix(:),1)
                            if (row < 100) then
                                pel%phl3(cur_3ph)%temp_tr_mix(count(pel%phl3(cur_3ph)%temp_tr_mix(:) /= 0.d0)+1) = pel%phl3(cur_3ph-1)%temp_tr_mix( row )
                                pel%phl3(cur_3ph)%press_tr_mix(count(pel%phl3(cur_3ph)%press_tr_mix(:) /= 0.d0)+1) = pel%phl3(cur_3ph-1)%press_tr_mix( row )
                            endif
                        endif
                    endif
                endif

                do i = 1, count(pel%phl3(cur_3ph)%press_tr_mix(:) /= 0.d0)!pel%phl3(cur_3ph)%np_3ph
                    !updated for mixed hydrate for up to - 2 - guests
                    if (gl%N_guests == 1) then
                        if ((pel%phl3(cur_3ph)%Eqtype == 'VLwH').or.(pel%phl3(cur_3ph)%Eqtype == 'VHIw')) then
                            ! Enthalpy of hydrate formation
                            write(threephlineunit,3021) pel%phl3(cur_3ph)%Eqtype, pel%phl3(cur_3ph)%temp_tr_mix(i), pel%phl3(cur_3ph)%press_tr_mix(i), pel%phl3(cur_3ph)%chempot(i,1), pel%phl3(cur_3ph)%chempot(i,2), &
                                & pel%phl3(cur_3ph)%lnfi(i,1), pel%phl3(cur_3ph)%lnfi(i,2), pel%phl3(cur_3ph)%x_ph(i,1,1), pel%phl3(cur_3ph)%x_ph(i,1,2), &
                                & pel%phl3(cur_3ph)%x_ph(i,2,1), pel%phl3(cur_3ph)%x_ph(i,2,2), pel%phl3(cur_3ph)%x_ph(i,3,1), pel%phl3(cur_3ph)%x_ph(i,3,2), &
                                & pel%phl3(cur_3ph)%rho_tr_mix(i,1), pel%phl3(cur_3ph)%rho_tr_mix(i,2), pel%phl3(cur_3ph)%rho_tr_mix(i,3), &
                                & pel%phl3(cur_3ph)%occup(i,1,1),pel%phl3(cur_3ph)%occup(i,2,1), pel%phl3(cur_3ph)%hyd_nr(i,1), pel%phl3(cur_3ph)%h_melt(i), pel%phl3(cur_3ph)%s_hyd(i), &
                                & pel%phl3(cur_3ph)%h_hyd(i), pel%phl3(cur_3ph)%tpd_sol(i), pel%phl3(cur_3ph)%a0_other_struct(i)      !Added the tangent plane distance of the hydrate former in the "other" phase here as well as the lattice parameter at ref. conditions
                        else
                            write(threephlineunit,3022) pel%phl3(cur_3ph)%Eqtype, pel%phl3(cur_3ph)%temp_tr_mix(i), pel%phl3(cur_3ph)%press_tr_mix(i), pel%phl3(cur_3ph)%chempot(i,1), pel%phl3(cur_3ph)%chempot(i,2), &
                                & pel%phl3(cur_3ph)%lnfi(i,1), pel%phl3(cur_3ph)%lnfi(i,2), pel%phl3(cur_3ph)%x_ph(i,1,1), pel%phl3(cur_3ph)%x_ph(i,1,2), &
                                & pel%phl3(cur_3ph)%x_ph(i,2,1), pel%phl3(cur_3ph)%x_ph(i,2,2), pel%phl3(cur_3ph)%x_ph(i,3,1), pel%phl3(cur_3ph)%x_ph(i,3,2), &
                                & pel%phl3(cur_3ph)%rho_tr_mix(i,1), pel%phl3(cur_3ph)%rho_tr_mix(i,2), pel%phl3(cur_3ph)%rho_tr_mix(i,3), &
                                & pel%phl3(cur_3ph)%occup(i,1,1), pel%phl3(cur_3ph)%occup(i,2,1), pel%phl3(cur_3ph)%hyd_nr(i,1), pel%phl3(cur_3ph)%s_hyd(i), pel%phl3(cur_3ph)%h_hyd(i), &
                                & pel%phl3(cur_3ph)%tpd_sol(i), pel%phl3(cur_3ph)%a0_other_struct(i)     !Added the tangent plane distance of the hydrate former in the "other" phase here as well as the lattice parameter at ref. conditions
                        end if
                    elseif (gl%N_guests == 2) then
                        write(threephlineunit,3023) pel%phl3(cur_3ph)%Eqtype, pel%phl3(cur_3ph)%temp_tr_mix(i), pel%phl3(cur_3ph)%press_tr_mix(i), pel%phl3(cur_3ph)%chempot(i,1), pel%phl3(cur_3ph)%chempot(i,2), pel%phl3(cur_3ph)%chempot(i,3), &
                            & pel%phl3(cur_3ph)%lnfi(i,1), pel%phl3(cur_3ph)%lnfi(i,2), pel%phl3(cur_3ph)%lnfi(i,3), &
                            & pel%phl3(cur_3ph)%x_ph(i,1,1), pel%phl3(cur_3ph)%x_ph(i,1,2), pel%phl3(cur_3ph)%x_ph(i,1,3), &
                            & pel%phl3(cur_3ph)%x_ph(i,2,1), pel%phl3(cur_3ph)%x_ph(i,2,2), pel%phl3(cur_3ph)%x_ph(i,2,3), &
                            & pel%phl3(cur_3ph)%x_ph(i,3,1), pel%phl3(cur_3ph)%x_ph(i,3,2), pel%phl3(cur_3ph)%x_ph(i,3,3), &
                            & pel%phl3(cur_3ph)%rho_tr_mix(i,1), pel%phl3(cur_3ph)%rho_tr_mix(i,2), pel%phl3(cur_3ph)%rho_tr_mix(i,3), &
                            & pel%phl3(cur_3ph)%occup(i,1,1),pel%phl3(cur_3ph)%occup(i,1,2), pel%phl3(cur_3ph)%occup(i,2,1), pel%phl3(cur_3ph)%occup(i,2,2), &
                            & pel%phl3(cur_3ph)%hyd_nr(i,1), pel%phl3(cur_3ph)%hyd_nr(i,2), pel%phl3(cur_3ph)%s_hyd(i), pel%phl3(cur_3ph)%h_hyd(i)
                    end if
                enddo
                close(threephlineunit)
            endif
        enddo
    elseif (phaseline == 4) then
        do eqtype_ind = 1,4
            if (any(pel%phl4(eqtype_ind)%Qexist_mix(:)) .eqv. .true.) then
                path = trim(pathIn) // trim(pel%CompsHdrtsAll)
                if (pel%phl4(eqtype_ind)%Qexist_mix(1)) then
                    path = trim(path) // '\four_phase_lines\' // trim(pel%phl4(eqtype_ind)%Q_point_mix(1))  !outputpath
                elseif (pel%phl4(eqtype_ind)%Qexist_mix(200)) then
                    path = trim(path) // '\four_phase_lines\' // trim(pel%phl4(eqtype_ind)%Q_point_mix(200))  !outputpath
                endif
                path = trim(path) // '_z='
                do i = 1,gl%ncomp-1
                    write(path,'(A,F4.2,A)')trim(path),gl%molfractions(i),trim('-')
                enddo
                write(path,'(A,F4.2,A)')trim(path),gl%molfractions(gl%ncomp),trim('.txt')
                !only open and write to file, when at least one point was saved
                open(newunit=fourphlineunit, file=trim(path), status='unknown', action='write', iostat=errorflag)
                if (errorflag /= 0) call messages(gl, pel,410,trim(''),trim(''),0,0,0)!write(*,*)'Error in fourphaselines: path does not exist. No File will be written for 4-phase line, but calculation for 3-phase line will continue.'

                if (pel%print_head .eqv. .true.) then
                    write(fourphlineunit,3013)'Qpointtype','TempQ','pressQ','z1','z2','z3','xQph11','xQph12','xQph13','xQph21','xQph22','xQph23','xQph31','xQph32','xQph33','xQph41','xQph42','xQph43','phasefrac1','phasefrac2','phasefrac3','phasefrac4','rho1','rho2','rho3','rho4'
                endif
                i = count(pel%phl4(eqtype_ind)%press_q_mix(:) /= 0.d0)
                do k = 1, i
                    if ((eqtype_ind == 1) .or. (eqtype_ind == 3)) then
                        if (k == pel%phl4(eqtype_ind)%Q_map(1))  write(fourphlineunit,'(A8)')'physical'
                    elseif ((eqtype_ind == 2) .or. (eqtype_ind == 4)) then
                        if (k == pel%phl4(eqtype_ind)%Q_map(1)) write(fourphlineunit,'(A8)')'physical'
                    endif
                    if ((k == i) .and. (dabs(pel%phl4(eqtype_ind)%Temp_Q_mix(k)) < 1.d-14 )) then
                        write(fourphlineunit,3024) pel%phl4(eqtype_ind)%Q_point_mix(200), pel%phl4(eqtype_ind)%Temp_Q_mix(200), pel%phl4(eqtype_ind)%press_Q_mix(200), gl%molfractions(1), gl%molfractions(2), gl%molfractions(3), pel%phl4(eqtype_ind)%x_Q_ph1_mix(1,200), pel%phl4(eqtype_ind)%x_Q_ph1_mix(2,200),pel%phl4(eqtype_ind)%x_Q_ph1_mix(3,200), &
                            & pel%phl4(eqtype_ind)%x_Q_ph2_mix(1,200), pel%phl4(eqtype_ind)%x_Q_ph2_mix(2,200), pel%phl4(eqtype_ind)%x_Q_ph2_mix(3,200), pel%phl4(eqtype_ind)%x_Q_ph3_mix(1,200), pel%phl4(eqtype_ind)%x_Q_ph3_mix(2,200), pel%phl4(eqtype_ind)%x_Q_ph3_mix(3,200), &
                            &  pel%phl4(eqtype_ind)%x_Q_ph4_mix(1,200), pel%phl4(eqtype_ind)%x_Q_ph4_mix(2,200), pel%phl4(eqtype_ind)%x_Q_ph4_mix(3,200), pel%phl4(eqtype_ind)%phasefrac_mix(1,200), pel%phl4(eqtype_ind)%phasefrac_mix(2,200), &
                            & pel%phl4(eqtype_ind)%phasefrac_mix(3,k), pel%phl4(eqtype_ind)%phasefrac_mix(4,k), pel%phl4(eqtype_ind)%rho_Q_mix(1,200),pel%phl4(eqtype_ind)%rho_Q_mix(2,200),pel%phl4(eqtype_ind)%rho_Q_mix(3,200),pel%phl4(eqtype_ind)%rho_Q_mix(4,200)
                    else
                        write(fourphlineunit,3024) pel%phl4(eqtype_ind)%Q_point_mix(k), pel%phl4(eqtype_ind)%Temp_Q_mix(k), pel%phl4(eqtype_ind)%press_Q_mix(k), gl%molfractions(1), gl%molfractions(2), gl%molfractions(3), pel%phl4(eqtype_ind)%x_Q_ph1_mix(1,k), pel%phl4(eqtype_ind)%x_Q_ph1_mix(2,k),pel%phl4(eqtype_ind)%x_Q_ph1_mix(3,k), &
                            & pel%phl4(eqtype_ind)%x_Q_ph2_mix(1,k), pel%phl4(eqtype_ind)%x_Q_ph2_mix(2,k), pel%phl4(eqtype_ind)%x_Q_ph2_mix(3,k), pel%phl4(eqtype_ind)%x_Q_ph3_mix(1,k), pel%phl4(eqtype_ind)%x_Q_ph3_mix(2,k), pel%phl4(eqtype_ind)%x_Q_ph3_mix(3,k), &
                            &  pel%phl4(eqtype_ind)%x_Q_ph4_mix(1,k), pel%phl4(eqtype_ind)%x_Q_ph4_mix(2,k), pel%phl4(eqtype_ind)%x_Q_ph4_mix(3,k), pel%phl4(eqtype_ind)%phasefrac_mix(1,k), pel%phl4(eqtype_ind)%phasefrac_mix(2,k), &
                            & pel%phl4(eqtype_ind)%phasefrac_mix(3,k), pel%phl4(eqtype_ind)%phasefrac_mix(4,k), pel%phl4(eqtype_ind)%rho_Q_mix(1,k),pel%phl4(eqtype_ind)%rho_Q_mix(2,k),pel%phl4(eqtype_ind)%rho_Q_mix(3,k),pel%phl4(eqtype_ind)%rho_Q_mix(4,k)
                    endif
                    if ((eqtype_ind == 1) .or. (eqtype_ind == 3)) then
                        if (k == pel%phl4(eqtype_ind)%Q_map(3)) write(fourphlineunit,'(A8)')'physical'
                    elseif ((eqtype_ind == 2) .or. (eqtype_ind == 4)) then
                        if (k == pel%phl4(eqtype_ind)%Q_map(4)) write(fourphlineunit,'(A8)')'physical'
                    endif
                enddo
                close(fourphlineunit)
            endif
        enddo
    endif

    !Changed format statements for tpd and lattice parameter output in order to examine stability of pure hydrate structure. Andreas Jäger November 2017
3011 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ', a8, ' ', a12, ' ', a12, ' ', a12, ' ', a12, ' ', a12)
    !Changed format statements for tpd and lattice parameter output in order to examine stability of pure hydrate structure. Andreas Jäger November 2017
3012 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ', a8, ' ', a12, ' ', a8 ' ', a12, ' ', a12, ' ', a12, ' ', a12)
3021 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ', f8.4, ' ', f12.4, ' ', f12.6, ' ', f12.4, ' ', f12.6, ' ', f8.2)

3022 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ', f8.4, ' ', f12.6, ' ', f12.4, ' ', f12.6, ' ', f8.2)

3013 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ',a12,' ',a12,' ', a8,' ', a8, ' ', a12, ' ', a12)

3023 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ', f8.4,' ', f8.4, ' ', f12.6, ' ', f12.4)

3024 format (a6,' ',f8.4,' ',f10.6,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f12.10,' ',f15.10,' ',f15.10,' ',f15.10,' ',f15.10' ',f12.6' ',f12.6' ',f12.6' ',f12.6)

    end subroutine phaselines_out


    module subroutine phaselines_expert(gl, pel, pathIn, additional_args, errorflag)
    !use module_all_types
    !use module_hdrt_phaselines

    implicit none

    type(type_gl), dimension(:) :: gl
    type(type_hdrt_pheq) :: pel

    double precision, dimension(30) :: xIn
    double precision, dimension(3,30) :: x_ph
    double precision, dimension(3) :: rho
    double precision :: p_start, p_end, press, press_step
    double precision :: T_start, T_end, Temp, Temp_step

    ! iTp = 1 ... temperature given, iTp = 2 ... pressure given
    integer :: iFlash, iTp
    integer :: errorflag

    character(255) :: pathIn
    character(256) :: additional_args
    character(30), dimension(30):: args_split


    integer :: cur_3ph, k, j, i

    if (gl(1)%N_Guests == 1) then
        pel%phl3%nflashs = 1
    elseif (gl(1)%N_Guests == 2) then
        pel%phl3%nflashs = 2
        do i = pel%n_3ph_lines,1,-1
            !pel%phl3%trline_mix(i*2-1:i*2) = pel%phl3%trline_mix(i)
            pel%phl3(i*2-1)%trline_mix = pel%phl3(i)%trline_mix
            pel%phl3(i*2)%trline_mix = pel%phl3(i)%trline_mix
        enddo
    endif

    k = 1

    do cur_3ph = 1,pel%n_3ph_lines * pel%phl3(1)%nflashs
        call phaselines_init(pel%phl3(cur_3ph))
    enddo
    call split_char(additional_args, args_split, k, ";")
    do k = 1,size(args_split,1)

        if (index(args_split(k),'pstart') > 0) then
            iTp = 2
            read(args_split(k)(index(args_split(k),'=')+1:),*)p_start
            pel%phl3%press_tr_mix(1) = p_start
        elseif (index(args_split(k),'pend') > 0) then
            read(args_split(k)(index(args_split(k),'=')+1:),*)p_end
            !pel%phl3%press_tr_mix(pel%phl3(1)%np_3ph) = p_end
        elseif (index(args_split(k),'Temp') > 0) then
            read(args_split(k)(index(args_split(k),'=')+1:),*)temp
            pel%PHL3%temp_tr_mix(1) = temp
        elseif (index(args_split(k),'Tstart') > 0) then
            iTp = 1
            read(args_split(k)(index(args_split(k),'=')+1:),*)T_start
            pel%phl3%temp_tr_mix(1) = T_start
        elseif (index(args_split(k),'Tend') > 0) then
            read(args_split(k)(index(args_split(k),'=')+1:),*)T_end
            !pel%phl3%press_tr_mix(pel%phl3(1)%np_3ph) = p_end
        elseif (index(args_split(k),'press') > 0) then
            read(args_split(k)(index(args_split(k),'=')+1:),*)press
            pel%PHL3%press_tr_mix(1) = press
        elseif (index(args_split(k),'xph1') > 0) then
            j = 1
            call split_dbl(args_split(k)(index(args_split(k),'=')+1:),pel%PHL3(1)%X_PH(1,1,:),j,":")
            if (dabs(sum(pel%PHL3(1)%X_PH(1,1,:))-1.d0) > 1.d-14) pel%PHL3(1)%X_PH(1,1,count(pel%PHL3(1)%X_PH(1,1,:) > 0.d0)+1) = 1.d0 - sum(pel%PHL3(1)%X_PH(1,1,1:count(pel%PHL3(1)%X_PH(1,1,:) > 0.d0)))
        elseif (index(args_split(k),'xph2') > 0) then
            j = 1
            call split_dbl(args_split(k)(index(args_split(k),'=')+1:),pel%PHL3(1)%X_PH(1,2,:),j,":")
            if (dabs(sum(pel%PHL3(1)%X_PH(1,2,:))-1.d0) > 1.d-14) pel%PHL3(1)%X_PH(1,2,count(pel%PHL3(1)%X_PH(1,2,:) > 0.d0)+1) = 1.d0 - sum(pel%PHL3(1)%X_PH(1,2,1:count(pel%PHL3(1)%X_PH(1,2,:) > 0.d0)))
        elseif (index(args_split(k),'xph3') > 0) then
            j = 1
            call split_dbl(args_split(k)(index(args_split(k),'=')+1:),pel%PHL3(1)%X_PH(1,3,:),j,":")
            if (dabs(sum(pel%PHL3(1)%X_PH(1,3,:))-1.d0) > 1.d-14) pel%PHL3(1)%X_PH(1,3,count(pel%PHL3(1)%X_PH(1,3,:) > 0.d0)+1) = 1.d0 - sum(pel%PHL3(1)%X_PH(1,2,1:count(pel%PHL3(1)%X_PH(1,3,:) > 0.d0)))
        else
        endif
    enddo
    do i = 1,4
        pel%phl4(i)%Q_map = 1
    enddo
    do cur_3ph = 1,pel%n_3ph_lines * pel%phl3(1)%nflashs
        xIn = gl(1)%molfractions
        x_ph = pel%PHL3(1)%X_PH(1,:,:)
        rho = 0.d0
        temp = pel%PHL3(cur_3ph)%temp_tr_mix(1)
        pel%phl3(cur_3ph)%EQTYPESTART = 1
        pel%phl3(cur_3ph)%Eqtype = pel%PHL3(cur_3ph)%trline_mix
        press = pel%PHL3(cur_3ph)%press_tr_mix(1)
        if (iTp == 1) then
            Temp_step = (t_end - pel%PHL3(cur_3ph)%temp_tr_mix(1)) / pel%PHL3(cur_3ph)%loops
        elseif (iTp == 2) then
            press_step = (p_end - pel%PHL3(cur_3ph)%press_tr_mix(1)) / pel%PHL3(cur_3ph)%loops
        endif
        do i = 1, pel%PHL3(cur_3ph)%np_3ph
            select case (trim(pel%PHL3(cur_3ph)%trline_mix))
            case ('VLwH')   ! ph1 = vap, ph2 = liq1, ph3 = hyd (phX = sol or liq2)
                if (cur_3ph == 1) then
                    if (iTp == 1) then
                        iflash = 5
                    elseif (iTp == 2) then
                        iflash = 2
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaH=0")
                elseif (cur_3ph == 2) then
                    if (iTp == 1) then
                        iflash = 1
                    elseif (iTp == 2) then
                        iflash = 2
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
                endif
            case ('VLxH')   ! ph1 = vap, ph2 = liq1, ph3 = hyd (phX = sol or liq2)
                if (cur_3ph == 1) then
                    if (iTp == 1) then
                        iflash = 4
                    elseif (iTp == 2) then
                        iflash = 3
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaLx=0")
                elseif (cur_3ph == 2) then
                    if (iTp == 1) then
                        iflash = 1
                    elseif (iTp == 2) then
                        iflash = 2
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
                endif
            case ('VHIx')   ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq1)
                if (cur_3ph == 1) then
                    if (iTp == 1) then
                        iflash = 4
                    elseif (iTp == 2) then
                        iflash = 1
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaIc=0")
                elseif (cur_3ph == 2) then
                    if (iTp == 1) then
                        write(*,*) "Unsupported input"
                        pause
                        stop
                    elseif (iTp == 2) then
                        iflash = -1
                    endif
                    pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
                endif
            end select
            call threephase_equil(gl(1), pel, pel%PHL3(cur_3ph)%trline_mix, cur_3ph, iFlash, temp, press, xIn , x_ph(1,:), x_ph(2,:), x_ph(3,:), rho(1), &
                rho(2), rho(3), pel%PHL3(cur_3ph)%Chempot(i,:), pel%PHL3(cur_3ph)%lnfi(i,:), pel%PHL3(cur_3ph)%occup(i,:,:), errorflag)!
                write(*,*)"w_C2H6O^L = ", x_ph(2,3) * gl(1)%wm(3) / (x_ph(2,1) * gl(1)%wm(1) + x_ph(2,2) * gl(1)%wm(2) + x_ph(2,3) * gl(1)%wm(3))
            if (errorflag == 0) then
                pel%PHL3(cur_3ph)%X_PH(i,1:3,:) = x_ph(1:3,:)
                pel%PHL3(cur_3ph)%RHO_TR_MIX(i,:) = rho
                pel%PHL3(cur_3ph)%temp_tr_mix(i) = temp
                pel%PHL3(cur_3ph)%press_tr_mix(i) = press
            endif
            if (iTp == 1) then
                temp = temp + temp_step
            elseif (iTp == 2) then
                press = press + press_step
            endif
        enddo
        pel%phl3(cur_3ph)%trline_mix = 'EP' // trim(pel%phl3(cur_3ph)%trline_mix)
    enddo
    call phaselines_out(gl(1), pel, 3, pathIn, errorflag)

    end subroutine phaselines_expert



    !******************************************************************************
    module subroutine phaselines (gl, pel, pathIn, x_vap_given, errorflag)
    !******************************************************************************
    ! Subroutine calculates three- and four-phase equilibrium lines and calls 'threephase_equil'

    implicit none
    type(type_gl), dimension(:) :: gl
    type(type_gl), dimension(:), allocatable :: gl_parallel
    type(type_hdrt_pheq) :: pel
    type(type_hdrt_pheq), dimension(:), allocatable :: pel_parallel

    !> current three phase line local variable
    integer :: cur_3ph
    !> unit for outputfiles
    integer :: threephlineunit, fourphlineunit

    character(255) :: pathIn, path, filename, dummy
    character(6) :: EqtypeIn, Eqtype        !Type of Equilibrium, which should be calculated
    logical :: Eqtype_found, VLwH_done


    !Variables at/for the three-phase equilibrium
    double precision :: Temp, press, presstemp
    double precision, dimension(30) :: xIn, x, x_ph1, x_ph2, x_ph3, x_phX, x_vap_given
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi
    double precision:: Dens, press_step, max_press, min_press
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est

    integer :: trloop, i, k, iTp, iflash

    integer :: row
    integer :: iPhase, iPhase_sol
    logical :: write_3ph
    !>Indicate which Qpoint is used for start values
    !!EqType is for quadruple line, iQ is for startpoint with desired phasefractions (plural!) = 0
    integer :: EqTypeStart, EqTypeEnd, iQ

    ! Mixed hydrates 2015
    double precision, dimension(30) :: xH
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(3) :: Phasefrac

    ! Properties for enthalpy of formation and hydration number
    integer :: ii, J, refin_VHIw
    double precision :: h_vap, h_hyd, s_hyd, h_liq, h_melt, h_phase, h_sol
    double precision, dimension(30) :: fug_g, hyd_nr!, occup_ls, occup_ld, occup_sms, occup_smd

    integer :: errorflag, eqtype_ind, avail_cores

    !New variables for tpd calculation
    double precision:: tpd_sol, ChemPot_hyd, a0_other_struct
    integer:: hdrt_structure_flag

    !//////////////////////////////////////////////////////////////////
    !Initialization of local variables
    write_3ph = .false.
    dummy = ''
    x = 0.d0
    x_ph1 = 0.d0
    x_ph2 = 0.d0
    x_ph3 = 0.d0
    x_phX = 0.d0
    !//////////////////////////////////////////////////////////////////
    !allocate(pel_parallel(pel%n_3ph_lines))
    avail_cores =  OMP_GET_NUM_PROCS()
    !avail_cores = 1
    !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
    write(*,*) "Performing calculations on up to",avail_cores," cores."
    !DEC$ END IF ! WO_WRITE_TO_CONSOLE
    call OMP_SET_NUM_THREADS(avail_cores)
    call OMP_SET_DYNAMIC(.false.)
    call OMP_SET_NESTED(.false.)

    !allocate (spltstr)

    if (errorflag /= 0) return


    x = gl(1)%molfractions

    !indicates on which position the pressure for specified ratio of guest2/guest3 in the vapor phase is: Start Value: -1
    do i = 1,size(pel%phl4)
        pel%phl4(i)%Q_map = -1
    enddo

    if (gl(1)%N_Guests == 1) then
        pel%phl3%nflashs = 1
        do i = 1,size(pel%phl4)
            if (pel%phl4(i)%Qexist_mix(1) .eqv. .true.) pel%phl4(i)%Q_map(:) = 1
        enddo
        !for ternary mixtures (2 guests) fourphaseline needs to be calculated (connects quadruple points of binary mixtures with water + guest1 and water + guest2
    elseif (gl(1)%N_Guests == 2) then
        do i = pel%n_3ph_lines,1,-1
            !pel%phl3%trline_mix(i*2-1:i*2) = pel%phl3%trline_mix(i)
            pel%phl3(i*2-1)%trline_mix = pel%phl3(i)%trline_mix
            pel%phl3(i*2)%trline_mix = pel%phl3(i)%trline_mix
        enddo
        pel%phl3%nflashs = 2
    endif

    if (pel%n_3ph_lines * pel%phl3(1)%nflashs < 4) then
        allocate(gl_parallel(4))
        do i = 1,4
            gl_parallel(i) = gl(1)
        enddo
    else
        allocate(gl_parallel(pel%n_3ph_lines * pel%phl3(1)%nflashs))
        do i = 1,pel%n_3ph_lines * pel%phl3(1)%nflashs
            gl_parallel(i) = gl(1)
            !pel_parallel(i) = pel
        enddo
    endif
    if (gl(1)%N_Guests == 2) then
        !$omp parallel do private(eqtype_ind) firstprivate(x, x_vap_given) lastprivate(errorflag)  schedule(static,1)
        do eqtype_ind = 1,4
            call phaselines_4(gl_parallel(eqtype_ind),pel%phl4(eqtype_ind),x, x_vap_given, eqtype_ind, errorflag)
        end do
        !$omp end parallel do

        if (errorflag /= 0) write(*,*)'phaselines_4: errorflag = ',errorflag

        call phaselines_out(gl(1), pel, 4, pathIn, errorflag)

    endif

    !$omp parallel do private(cur_3ph) lastprivate(errorflag)  schedule(static,1)
    do cur_3ph = 1,pel%n_3ph_lines * pel%phl3(1)%nflashs
        !threephaselineloop: do trloop = 1,pel%n_3ph_lines
        !flashloop: do k = 1,pel%phl3(trloop*2)%nflashs
        !cur_3ph = cur_3ph + 1
        call phaselines_3(gl_parallel(cur_3ph), pel, x, x_vap_given, cur_3ph, errorflag)
        !end do flashloop
        !end do threephaselineloop
    enddo
    !$omp end parallel do
    if (errorflag == 0) then
        call phaselines_out(gl(1), pel, 3, pathIn, errorflag)
    endif

    !Changed format statements for tpd and lattice parameter output in order to examine stability of pure hydrate structure. Andreas Jäger November 2017
3011 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ', a8, ' ', a12, ' ', a12, ' ', a12, ' ', a12, ' ', a12)
    !Changed format statements for tpd and lattice parameter output in order to examine stability of pure hydrate structure. Andreas Jäger November 2017
3012 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ', a8, ' ', a12, ' ', a8 ' ', a12, ' ', a12, ' ', a12, ' ', a12)
3021 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ', f8.4, ' ', f12.4, ' ', f12.6, ' ', f12.4, ' ', f12.6, ' ', f8.2)

3022 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ', f8.4, ' ', f12.6, ' ', f12.4, ' ', f12.6, ' ', f8.2)

3013 format (a6,' ',a8,' ',a12,' ',a16,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12,' ',a12 &
        &,' ',a12,' ',a12,' ',a12,' ',a8,' ',a8,' ',a8,' ',a12,' ',a12,' ',a12,' ',a12,' ', a8,' ', a8, ' ', a12, ' ', a12)

3023 format (a6,' ',f8.4,' ',f12.6,' ',f16.8,' ',f12.4,' ',f12.4,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',f8.2,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ', f8.4,' ', f8.4, ' ', f12.6, ' ', f12.4)

3024 format (a6,' ',f8.4,' ',f10.6,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f12.10,' ',f15.10,' ',f15.10,' ',f15.10,' ',f15.10)
    !    xQph41     xQph42      xQph43  phasefrac1  phasefrac2  phasefrac3  phasefrac4



    end subroutine phaselines

    module subroutine phaselines_3 (gl, pel, xIn, x_vap_given, cur_3ph, errorflag)
    !******************************************************************************
    ! Subroutine calculates three-phase equilibrium lines and calls 'threephase_equil'




    implicit none
    type(type_gl) :: gl
    type(type_gl), dimension(:), allocatable :: gl_parallel
    type(type_hdrt_pheq) :: pel
    type(type_hdrt_pheq), dimension(:), allocatable :: pel_parallel

    !> current three phase line local variable
    integer :: cur_3ph

    character(255) :: path, filename, dummy
    character(6) :: EqtypeIn, Eqtype        !Type of Equilibrium, which should be calculated
    logical :: Eqtype_found, VLwH_done


    !Variables at/for the three-phase equilibrium
    double precision :: Temp, press, presstemp
    double precision, dimension(30) :: xIn, x, x_ph1, x_ph2, x_ph3, x_phX, x_vap_given
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi
    double precision:: Dens, press_step, max_press, min_press
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est

    integer :: trloop, i, k, iTp, iflash

    integer, dimension(1) :: row
    integer :: iPhase, iPhase_sol
    logical :: write_3ph
    !>Indicate which Qpoint is used for start values
    !!EqType is for quadruple line, iQ is for startpoint with desired phasefractions (plural!) = 0
    integer :: EqTypeStart, EqTypeEnd, iQ

    ! Mixed hydrates 2015
    double precision, dimension(30) :: xH
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(3) :: Phasefrac

    ! Properties for enthalpy of formation and hydration number
    integer :: ii, J, refin_VHIw
    double precision :: h_vap, h_hyd, s_hyd, h_liq, h_melt, h_phase, h_sol
    double precision, dimension(30) :: fug_g, hyd_nr!, occup_ls, occup_ld, occup_sms, occup_smd

    integer :: errorflag, eqtype_ind

    !New variables for tpd calculation
    double precision:: tpd_sol, ChemPot_hyd, a0_other_struct
    integer:: hdrt_structure_flag
    !//////////////////////////////////////////////////////////////////
    !Initialization of local variables
    filename = ''
    dummy = ''
    x_ph1 = 0.d0
    x_ph2 = 0.d0
    x_ph3 = 0.d0
    x_phX = 0.d0
    beta = 0.D0
    Phasefrac = 0.D0

    !The iPhase_sol can be used to assume fluid phases in two and three phase equilibrium calculation
    ! 0 : Let the algorithms choose the phases
    ! 1 : Liquid phase of the hydrate guest is assumed (e.g. LwH or LwLcH)
    ! 2 : Vapor phase of the hydrate guest is assumed (e.g. VH or VLwH
    !NOTE: IN case of two fluid phases, the heavier phase is always assumed as liquid phase (no vapor / vapor equilibrium)
    iPhase_sol = 0
    !!count three phase lines that need to be calculated
    !pel%n_3ph_lines = 1
    !do i = 1,len(trim(pel%threephaselinestr))
    !    if ((pel%threephaselinestr(i:i) == ';') .and. (i < len(trim(pel%threephaselinestr))))then
    !        pel%threephaselinestr(i:i) = ','
    !        pel%n_3ph_lines = pel%n_3ph_lines + 1
    !    elseif (pel%threephaselinestr(i:i) == ' ') then
    !        pel%threephaselinestr(i:) = pel%threephaselinestr(i+1:)
    !    endif
    !enddo

    !n_3ph_lines = 6 ! 6 different 3 phase lines with 2 different betas -> 2*6 = 12
    !allocate(spltstr.vector(pel%n_3ph_lines*pel%phl3%nflashs))
    !spltstr.splitstring = pel%threephaselinestr
    !call splitstr(spltstr)
    !pel%phl3%trline_mix = spltstr.vector
    !if (gl%N_guests >= 2) then
    !    ii = 1
    !    do i = 1,pel%n_3ph_lines
    !        pel%phl3%trline_mix(ii) = spltstr.vector(i)
    !        pel%phl3%trline_mix(ii + 1) = spltstr.vector(i)
    !        ii = ii + 2
    !    enddo
    !endif
    !deallocate(spltstr.vector)
    VLwH_done = .false.
    !cur_3ph = 1
    !call setup(gl,inptorig, temp, press, fluidl_mod , molesl_mod, moles_mod, path_old, eos_indicator_mod, errorflag)
    Eqtype = ''
    iFlash = 0
    temp = 0.d0
    press = 0.d0
    if (pel%CompositionType /= 2) x_ph1 = 0.d0
    x_ph2 = 0.d0
    x_ph3 = 0.d0
    rho_ph1 = 0.d0
    rho_ph2 = 0.d0
    rho_ph3 = 0.d0
    Chempot = 0.d0
    lnfi = 0.d0
    occup = 0.d0
    errorflag = 0
    write_3ph = .false.
    Eqtype_found = .false.
    !//////////////////////////////////////////////////////////////////
    x = xIn

    EqtypeIn = trim(pel%phl3(cur_3ph)%trline_mix)

    Eqtype = EqtypeIn
    !Write the equation type in a more general form and
    !set EqTypeStart for the correct fourphaseline to start from
    !-----------------------------------------------------------
    If ((EqtypeIn == 'LwLcH') .or. (EqtypeIn == 'LwLpH') .or. (EqtypeIn == 'LwLeH')) then
        Eqtype = 'LwLxH'
        EqTypeStart = 2
        EqtypeEnd = 4
        Eqtype_found = .true.

    ElseIf ((EqtypeIn == 'VLcH') .or. (EqtypeIn == 'VLpH') .or. (EqtypeIn == 'VLeH')) then
        Eqtype = 'VLxH'
        !special case for co2: try to start from VLcHIc fourphaseline to calculate VLcH threephase line
        if ((any(gl%components(:) == 'co2')).and.(any(gl%components(:) == 'methane'))) then
            EqTypeStart = 3
            EqtypeEnd = 2
            Eqtype_found = .true.
        endif
        !for all other guest components
        if (Eqtype_found .neqv. .true.) then
            EqTypeStart = 2
            EqtypeEnd = 3
            Eqtype_found = .true.
        endif

    Elseif (EqtypeIn == 'VLwH') then
        EqTypeStart = 1
        EqtypeEnd = 2
        Eqtype_found = .true.

    ElseIf (EqtypeIn == 'VHIc') then
        Eqtype = 'VHIx'
        EqTypeStart = 3
        EqtypeEnd = 0
        Eqtype_found = .true.

    ElseIf (EqtypeIn == 'VHIw') then
        EqTypeStart = 1
        EqtypeEnd = 0
        Eqtype_found = .true.

    ElseIf (EqtypeIn == 'LcHIc') then
        Eqtype = 'LxHIx'
        EqTypeStart = 3
        EqtypeEnd = 4
        Eqtype_found = .true.

    ElseIf (EqtypeIn == 'LwHIc') then
        Eqtype = 'LwHIx'
        EqTypeStart = 4
        EqtypeEnd = 0
        Eqtype_found = .true.

    ElseIf ((EqtypeIn == 'VLcLw').or.(EqtypeIn == 'VLpLw').or.(EqtypeIn == 'VLeLw')) then
        Eqtype = 'VLxLw'
        EqTypeStart = 2
        EqtypeEnd = 0
        Eqtype_found = .true.

    ElseIf (EqtypeIn == 'VLwIw') then
        EqTypeStart = 1
        EqtypeEnd = 0
        Eqtype_found = .true.
    Elseif (EqtypeIn == 'LwHIw') then
        !pure hydrates only
        EqTypeStart = 1
        EqtypeEnd = 0
        Eqtype_found = .true.
    endif
    pel%phl3(cur_3ph)%EqTypeStart = EqTypeStart
    !-----------------------------------------------------------
    !threephase equilibrium does not exist in this system OR is not implemented (yet)
    if (EqType_found .eqv. .false.) then
        call messages(gl, pel,5, trim(EqtypeIn), trim(''),0,0,0)
        !pause
        return
        !fourphase equilibrium found (=startvalues) for requested three phase line
        !startvalues
    elseif (EqType_found .eqv. .true.) then
        call phaselines_init(pel%phl3(cur_3ph))

        rhoph1_est = 0.D0
        rhoph2_est = 0.D0
        rhoph3_est = 0.D0
        rhophX_est = 0.D0
        rho_ph1 = 0.d0
        rho_ph2 = 0.d0
        rho_ph3 = 0.d0
        !default values for flash calculation
        iQ = 1
        pel%phl3(cur_3ph)%strFlash = ''
        iFlash = 2
        max_press = 0.d0
        min_press = 0.d0
        press_step = 0.d0
        refin_VHIw = 0
        if (gl%N_guests == 1) then
            k = 1
        elseif (gl%N_guests == 2) then
            k = 2 - mod(cur_3ph,2)
        endif

        !set the startvalue depending on eqtype and phaseboundary
        select case (trim(Eqtype))
            !
        case ('VLwH')   ! ph1 = vap, ph2 = liq2, ph3 = hyd (phX = sol)
            if (k == 1) then
                iflash = 2
                !flexible determination of startvalues:
                !for VLwH_betaH=0 the hydrate and iced water phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaH=0")
            elseif (k == 2) then
                iflash = 3
                !flexible determination of startvalues:
                !for VLwH_betaLw=0 the liquid water and iced water phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaLw=0")
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(1)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    Temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    !vap  (Q = VLwHIw)
                    x_ph2 = pel%phl4(EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    !liq2 (Q = VLwHIw)
                    x_phX = 0.D0
                else
                    return
                end if
                if (k == 1) then !betaH = 0
                    call phaseemerge_detection(pel,EqTypeEnd, (/3,4/), iQ)
                    if (iQ == 6) then
                        call phaseemerge_detection(pel,EqTypeEnd,(/1,4/),iQ)
                    endif
                    if (iQ == 6) then
                        EqTypeEnd = 4
                        call phaseemerge_detection(pel,EqTypeEnd,(/3,4/),iQ)
                    endif
                elseif (k == 2) then!betaLw=0
                    call phaseemerge_detection(pel,EqTypeEnd, (/2,3/), iQ)
                    if (iQ == 6) then
                        call phaseemerge_detection(pel,EqTypeEnd,(/1,4/),iQ)
                    endif
                    if (iQ == 6) then
                        EqTypeEnd = 4
                        call phaseemerge_detection(pel,EqTypeEnd,(/1,4/),iQ)
                    endif
                endif
                if (iQ < 6) then
                    max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                else
                    if ((gl%N_Guests == 1).and.(pel%phl4(2)%Qexist_mix(1))) then
                        EqTypeEnd = 2
                        max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                    else
                        max_press = 1000.d0
                    endif
                endif
                press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
            endif
        case ('VHIw')   ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq2)
            if (k == 1) then
                iflash = 2
                !flexible determination of startvalues:
                !for VHIw_betaH=0 the hydrate and liquid water phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,3/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaH=0")
            elseif (k == 2) then
                iflash = 1
                iQ = 3
                !flexible determination of startvalues:
                !for VHIw_betaIw=0 the iced water and liquid water phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaIw=0")
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(1)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(1)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(1)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(1)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    !vap (Q = VLwHIw)
                else
                    return
                end if

                x_ph3 = 0.D0
                x_ph3(1) = 1.D0
                min_press = 0.009D0
                press_step = (min_press - press) / pel%phl3(cur_3ph)%loops

                refin_VHIw = 1
            endif

        case ('LwHIw')   ! ph1 = liq2, ph2 = hyd, ph3 = sol (phX = liq1)
            if (gl%n_guests /= 1) then
                errorflag = -10051
            endif
            if (errorflag == 0) then
                !Only for pure hydrates ATM
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(EqTypeStart)%Q_map(1) == -1) then
                    errorflag = -10050
                    write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
                else
                    if (pel%phl4(1)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(1))) then
                        temp = pel%phl4(1)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(1))
                        press = pel%phl4(1)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(1))
                        x_ph1 = pel%phl4(1)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(1))    !liq2 (Q = VLwHIw)
                    else
                        return
                    end if
                    x_ph3 = 0.D0
                    x_ph3(1) = 1.D0
                    min_press = 280.D0
                    press_step = (min_press - press) / pel%phl3(cur_3ph)%loops
                endif
            endif
        case ('LwLxH')  ! ph1 = liq1, ph2 = liq2, ph3 = hyd (phX = sol) ... it is LxLwH in fact
            if (k == 1) then
                iflash = 2
                !flexible determination of startvalues:
                !for LwLxH_betaH=0 the hydrate and vapor phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) then!if no root on VLLH line for LwLcH betaH = 0 found, take root for betaLw = 0
                    call phaseemerge_detection(pel,EqTypeStart, (/1,2/), iQ)
                    if (iQ == 6) then
                        LLH_parachute: do i = 1,size(pel%phl3,1)
                            if ((trim(pel%phl3(i)%trline_mix(1:4)) == 'VLcH') .or. (trim(pel%phl3(i)%trline_mix(1:4)) == 'VLeH') .or. (trim(pel%phl3(i)%trline_mix(1:4)) == 'VLpH')) then
                                if (dabs(maxval(pel%phl3(i)%press_tr_mix(:)) - maxval(pel%phl4(2)%press_Q_mix(:)))/maxval(pel%phl4(2)%press_Q_mix(:)) < 5.d-2) then
                                    iQ = 5
                                    pel%phl4(EqTypeStart)%Q_map(iQ) = maxloc(pel%phl4(2)%press_Q_mix(:),1)
                                    exit LLH_parachute
                                endif
                            endif
                        enddo LLH_parachute
                    endif
                    if (iQ == 6) then !if still no startvalues try highest pressure of VLxLwH line as startpoint
                        iQ = 5
                        pel%phl4(EqTypeStart)%Q_map(iQ) = maxloc(pel%phl4(2)%press_Q_mix(:),1)
                    endif
                endif
                pel%phl3(cur_3ph)%strFlash = trim("_betaH=0")
            elseif (k == 2) then
                iflash = 3
                !flexible determination of startvalues:
                !for LwLxH_betaLw=0 the liquid water and vapor phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,2/), iQ)
                if (iQ == 6) then!if no root on VLLH line for LwLcH betaLw = 0 found, take root for betaH = 0
                    call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                endif
                pel%phl3(cur_3ph)%strFlash = trim("_betaLw=0")
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(2)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    Temp = pel%phl4(2)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(2)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(2)%x_Q_ph3_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    ! liq1 (Q = VLwLxH)
                    x_ph2 = pel%phl4(2)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    ! liq2 (Q = VLwLxH)
                End if
                x_phX = 0.D0
                if (k == 1) then!betaH = 0
                    call phaseemerge_detection(pel,EqTypeEnd,(/3,4/),iQ)
                elseif (k == 2) then!betaLw = 0
                    call phaseemerge_detection(pel,EqTypeEnd,(/1,4/),iQ)
                endif
                if (pel%phl4(EqTypeEnd)%Q_map(iQ) > 0) then
                    max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                else
                    max_press = 1000.d0
                endif
                press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
            endif
        case ('VLxH')   ! ph1 = vap, ph2 = liq1, ph3 = hyd (phX = sol or liq2)
            if (EqtypeStart == 2) then
                if (k == 1) then
                    iflash = 3
                    !flexible determination of startvalues:
                    !for VLxH_betaLx=0 the liquid guestX rich and liquid water phases emerge
                    call phaseemerge_detection(pel,EqTypeStart,(/2,3/),iQ)
                    if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart,(/3,4/),iQ)
                    pel%phl3(cur_3ph)%strFlash = trim("_betaLx=0")
                elseif (k == 2) then
                    iflash = 1
                    !flexible determination of startvalues:
                    !for VLxH_betaV=0 the vapor and liquid water phases emerge
                    call phaseemerge_detection(pel,EqTypeStart,(/1,2/),iQ)
                    if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart,(/3,4/),iQ)
                    pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
                endif
            elseif (EqtypeStart == 3) then
                if (k == 1) then
                    iflash = 3
                    !flexible determination of startvalues:
                    !for VLxH_betaLx=0 the liquid guestX rich and liquid water phases emerge
                    call phaseemerge_detection(pel,EqTypeStart,(/2,4/),iQ)
                    if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart,(/1,4/),iQ)
                    pel%phl3(cur_3ph)%strFlash = trim("_betaLx=0")
                elseif (k == 2) then
                    iflash = 1
                    !flexible determination of startvalues:
                    !for VLxH_betaV=0 the vapor and liquid water phases emerge
                    call phaseemerge_detection(pel,EqTypeStart,(/1,4/),iQ)
                    pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
                endif
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(EqTypeStart)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    Temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    if (EqtypeStart == 2) then
                        x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))      ! vap  (Q = VLwLxH)
                        x_ph2 = pel%phl4(EqTypeStart)%x_Q_ph3_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))      ! liq1 (Q = VLwLxH)
                    elseif (EqtypeStart == 3) then
                        x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))      ! vap  (Q = VLxHIc)
                        x_ph2 = pel%phl4(EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))      ! liq1 (Q = VLxHIc)
                    endif
                End if

                x_phX = 0.D0
                if (EqTypeStart == 2) then
                    if (k == 1) then !betaLx = 0
                        call phaseemerge_detection(pel,EqTypeEnd, (/2,4/), iQ)
                    elseif (k == 2) then !betaV = 0
                        call phaseemerge_detection(pel,EqTypeEnd, (/1,4/), iQ)
                    endif
                    if (iQ < 6) then
                        min_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                    else
                        !min_press = 0.5d0
                        min_press = 0.01d0
                    endif
                    press_step = (min_press - press) / pel%phl3(cur_3ph)%loops
                elseif (EqTypeStart == 3) then
                    if (k == 1) then !betaLx = 0
                        call phaseemerge_detection(pel,EqTypeEnd, (/2,3/), iQ)
                    elseif (k == 2) then !betaV = 0
                        call phaseemerge_detection(pel,EqTypeEnd, (/1,2/), iQ)
                    endif
                    if (iQ < 6) then
                        max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                    else
                        !max_press = 10.d0
                        if (gl%N_guests == 1) then
                            max_press = pel%phl4(2)%press_Q_mix(1)
                        else
                            max_press = maxval(pel%phl4(EqTypeStart)%press_Q_mix(:)) * 2.d0
                        endif
                    endif
                    press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
                endif
            endif

        case ('VHIx')   ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq1)
            if (k == 1) then
                iflash = 1
                iQ = 1
                !flexible determination of startvalues:
                !for VHIx_betaIc=0 the solid co2 and liquid guestX rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaIc=0")
            elseif (k == 2) then
                iflash = -1
                iQ = 3
                !flexible determination of startvalues:
                !for VHIx_betaIc=0 the vapor and liquid guestX rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,2/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
            endif
            !Only for pure hydrates ATM
            !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                if (pel%phl4(EqTypeStart)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    if ((gl%molfractions(3) > 0.065d0 ).and.(k == 2)) press = press * 0.95d0
                    x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))    ! vap  (Q = VLxHIx)
                End if
                x_ph3 = 0.D0
                x_ph3(2) = 1.D0

                min_press = 0.01D0
                press_step = (min_press - press) / pel%phl3(cur_3ph)%loops
            endif

        case ('LwHIx')
            if (k == 1) then
                iflash = -1
                iQ = 3
                !flexible determination of startvalues:
                !for LwHIx_betaLw=0 the liquid water and liquid guestX rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,2/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaLw=0")
            elseif (k == 2) then
                iflash = 2
                !flexible determination of startvalues:
                !for LwHIx_betaH=0 the hydrate and liquid guestX rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,3/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaH=0")
            endif
            ! ph1 = liq2, ph2 = hyd, ph3 = sol (phX = vap)
            !Only for pure hydrates ATM
            !Starting point is Quadruple Point 4
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,'(I)')errorflag
            else
                if (pel%phl4(4)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(4)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(4)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(4)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))  ! liq2  (Q = LwLxHIx)
                else
                    return
                end if
                x_ph3 = 0.D0
                x_ph3(2) = 1.D0
                max_press = 1500.D0
                press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
            endif

        case ('LxHIx')  ! ph1 = liq1, ph2 = hyd, ph3 = sol (phX = vap)
            if (k == 1) then
                iflash = 1
                iQ = 2
                !flexible determination of startvalues:
                !for LxHIx_betaIx=0 the solid co2 and liquid water rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaIx=0")
            elseif (k == 2) then
                iflash = -1
                iQ = 3
                !flexible determination of startvalues:
                !for LxHIx_betaLx=0 the liquid guestX rich and liquid water rich phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,2/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaLx=0")
            endif
            !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
            if ((pel%phl4(EqTypeStart)%Q_map(iQ) == -1) .and. (k == 1)) then !beta Ix = 0
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            elseif (((pel%phl4(EqTypeStart)%Q_map(iQ) == -1) .and. (k == 2))) then!betaLx = 0
                call phaseemerge_detection(pel,EqTypeEnd, (/1,2/), iQ)
                if (pel%phl4(EqTypeEnd)%Q_map(iQ) /= -1) then
                    EqTypeStart = 4
                    EqTypeEnd = 3
                endif
                if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                    errorflag = -10050
                    write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
                elseif (pel%phl4(EqTypeStart)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))  ! liq1  (Q = LxLwHIx)
                else
                    return
                End if

                x_ph3 = 0.D0
                x_ph3(2) = 1.D0
                min_press = 5.5d0
                press_step = (min_press - press) / pel%phl3(cur_3ph)%loops
            else
                if (pel%phl4(EqTypeStart)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))  ! liq1  (Q = VLxHIx)
                else
                    return
                End if
                x_ph3 = 0.D0
                x_ph3(2) = 1.D0
                if (k == 1) then!beta Ix = 0
                    call phaseemerge_detection(pel,EqTypeEnd, (/1,4/), iQ)
                elseif (k == 2) then!beta Lx = 0
                    call phaseemerge_detection(pel,EqTypeEnd, (/1,2/), iQ)
                endif
                if (iQ < 6) then
                    max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(iQ))
                else
                    if (gl%N_guests == 1) then
                        EqTypeEnd = 4
                        max_press = pel%phl4(EqTypeEnd)%press_Q_mix(pel%phl4(EqTypeEnd)%Q_map(1))
                    else
                        max_press = 1000.d0
                    endif
                endif
                press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
            endif

        case ('VLxLw')  ! ph1 = vap, ph2 = liq1, ph3 = liq2 (phX = vap)
            if (k == 1) then
                iflash = 2
                !flexible determination of startvalues:
                !for VLxLw_betaLx=0 the liquid guestX rich and hydrate phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaLx=0")
            elseif (k == 2) then
                iflash = 1
                !flexible determination of startvalues:
                !for VLxLw_betaV=0 the vapor and hydrate phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/1,4/), iQ)
                if (iQ == 6) call phaseemerge_detection(pel,EqTypeStart,(/3,4/),iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaV=0")
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !Only for pure hydrates ATM
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                !First estimate: always take V
                if (pel%phl4(EqTypeStart)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(EqTypeStart)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(EqTypeStart)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(EqTypeStart)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))        ! vap (Q = VLwLxH)
                    x_ph2 = pel%phl4(EqTypeStart)%x_Q_ph3_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))        ! liq1
                    x_ph3 = pel%phl4(EqTypeStart)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))        ! liq2
                    !Temp = 283.078715104786d0
                    !press = 4.48014869363662d0
                    !x_ph1(1) = 6.710531532570730D-004
                    !x_ph1(2) = 0.999328946846743d0
                    !x_ph2(1) = 2.347845147541576D-003
                    !x_ph2(2) = 0.997652154852458d0
                    !x_ph3(1) = 0.971602961042455d0
                    !x_ph3(2) = 2.839703895754497D-002
                else
                    return
                End if
                if (gl%n_guests == 1) then
                    max_press = 0.9*gl%pc(2)
                else
                    !max_press = 1.15*pc(2)!6.D0
                    if ((any(gl%components == 'propane')) .and. (any(gl%components == 'nitrogen')) .and. (k == 2)) then
                        max_press = 3.d0
                    elseif (pel%phl4(EqTypeStart)%fourphase_physical .eqv. .true.) then
                        max_press = 1.25*maxval(pel%phl4(EqTypeStart)%press_Q_mix(:))
                    elseif (pel%phl4(EqTypeStart)%fourphase_physical .eqv. .false.) then
                        max_press = 1.1*maxval(pel%phl4(EqTypeStart)%press_Q_mix(:))
                    endif
                endif
                !UNCOMMENT THE FOLLOWING SECTION TO CALCULATE VLwLx INDEPENDENTLY OF QUADRUPLE POINT VLxLwH
                !-------------------------------------------------------
                !press = 3.D0 !CO2=3.5;C2H6=3;C3H8=0.4;
                !x_liq2(1) = 0.99999D0
                !x_liq2(2) = 0.00001D0
                !x_liq1(1) = 0.01D0
                !x_liq1(2) = 0.99D0
                !x_vap(1) = 0.001D0
                !x_vap(2) = 0.999D0
                !Temp = 278.D0     !CO2=283.D0;C2H6=;C3H8=278.D90
                !--------------------------------------------------------

                press_step = (max_press - press) / pel%phl3(cur_3ph)%loops
            endif
        case ('VLwIw')     ! ph1 = vap, ph2 = liq2, ph3 = sol (phX = liq1)
            if (k == 1) then
                iflash = 2
                !flexible determination of startvalues:
                !for VLxLw_betaIw=0 the solid water and hydrate phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/3,4/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaIw=0")
            elseif (k == 2) then
                iflash = 3
                !flexible determination of startvalues:
                !for VLxLw_betaLw=0 the liquid water and hydrate phases emerge
                call phaseemerge_detection(pel,EqTypeStart, (/2,3/), iQ)
                pel%phl3(cur_3ph)%strFlash = trim("_betaLw=0")
            endif
            if (pel%phl4(EqTypeStart)%Q_map(iQ) == -1) then
                errorflag = -10050
                write(pel%phl3(cur_3ph)%trline_mix,*)errorflag
            else
                !Only for pure hydrates ATM
                !The first value will be calculated at a given startvalue for temperature, pressure and phase composition
                if (pel%phl4(1)%Qexist_mix(pel%phl4(EqTypeStart)%Q_map(iQ))) then
                    temp = pel%phl4(1)%Temp_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    press = pel%phl4(1)%press_Q_mix(pel%phl4(EqTypeStart)%Q_map(iQ))
                    x_ph1 = pel%phl4(1)%x_Q_ph1_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))        ! vap (Q = VLwHIw)
                    x_ph2 = pel%phl4(1)%x_Q_ph2_mix(:,pel%phl4(EqTypeStart)%Q_map(iQ))        ! liq2
                End if
                x_ph3 = 0.D0
                x_ph3(1) = 1.D0

                !min_press = 0.0007D0
                min_press = 0.005D0
                !UNCOMMENT THE FOLLOWING SECTION TO CALCULATE VLwIw INDEPENDENTLY OF QUADRUPLE POINT VLwHIw
                !-------------------------------------------------------
                !min_press = 1.D0 !CO2=1.5;N2=;CO=20;Ar=12;O2=20;CH4=3;C2H6=1;C3H8=0.25;
                !x_liq2(1) = 0.9999976197D0
                !x_liq2(2) = 0.0000023803D0
                !x_vap(1) = 0.05D0
                !x_vap(2) = 0.95D0
                !Temp = 273.1597D0
                !press = 0.007D0
                !--------------------------------------------------------
                press_step = (min_press - press) / pel%phl3(cur_3ph)%loops
            endif

        end select  ! Eqtype


        !no fourphase equilibrium found or calculated to get startvalues for requested three phase line / area -> jump to next threephaseline
        if ((errorflag /= 0) .and. (EqType_found .eqv. .true.)) then
            call messages(gl, pel,6, trim(EqtypeIn) // dummy, trim(pel%phl3(cur_3ph)%strFlash) // dummy,0,0,0)
            !cur_3ph = cur_3ph + 1
            errorflag = 0
            !cycle
            return
        endif

        if (gl%n_guests == 1) pel%phl3(cur_3ph)%strFlash = ""
        ! -----------------------------------
        ! Calculation of the three-phase line
        ! -----------------------------------
98      i = 1
        pressureloop: Do i = 1, pel%phl3(cur_3ph)%loops + 1
            if (pel%show_progress .eqv. .true.) then
                !write(*,'(''+'',A,A,I5,A4,I5,A)')trim(Eqtype),': calculating point ',i, ' of ',loops+1,' ' // trim(pel%phl3(cur_3ph)%strFlash(2:)) // '          '
                call messages(gl, pel,100, trim(Eqtype) // dummy, trim(pel%phl3(cur_3ph)%strFlash) // dummy, i, pel%phl3(cur_3ph)%loops,0)
            elseif ((pel%show_progress .eqv. .false.) .and. (i == 1)) then
                !write(*,*)'Calculation of ' // trim(Eqtype) // ' ' // trim(pel%phl3(cur_3ph)%strFlash(2:)) // '         '
                call messages(gl, pel,8, trim(Eqtype) // dummy,  trim(pel%phl3(cur_3ph)%strFlash) // dummy,0,0,0)
            endif
            !#############################################
            !BEGIN parameter study water content variation
            !if (trim(eqtype) == 'VLwH') then
            !    x_vap_given = 0.d0
            !    x_vap_given(2) = x(2)
            !    x_vap_given(3) = x(3)
            !    open(unit=999,file='D:\test.txt')
            !    do ii = 20,70
            !        x(1) = ii/100.d0
            !        x(2) = (1.d0 - x(1))*x_vap_given(2)/(x_vap_given(2)+x_vap_given(3))
            !        x(3) = (1.d0 - x(1))*x_vap_given(3)/(x_vap_given(2)+x_vap_given(3))
            !        call threephase_equil(gl,pel,Eqtype,cur_3ph, iFlash,Temp, press, x, x_ph1, x_ph2, x_ph3, &
            !                        &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, errorflag, &
            !                        &   Q_point_mix(:,Q_map(EqTypeStart, 1)), Qexist_mix(:,Q_map(EqTypeStart, 1)),Temp_Q_mix(:,Q_map(EqTypeStart, 1)), press_Q_mix(:,Q_map(EqTypeStart, 1)), &
            !                        &   x_Q_ph1_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph2_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph3_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph4_mix(:,:,Q_map(EqTypeStart, 1)))
            !        call threephase_equil_iter(gl,Eqtype,cur_3ph, iFlash,Temp, press, x, x_ph1, x_ph2, x_ph3, &
            !                &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, errorflag, &
            !                &   Q_point_mix(:,Q_map(EqTypeStart, 1)), Qexist_mix(:,Q_map(EqTypeStart, 1)),Temp_Q_mix(:,Q_map(EqTypeStart, 1)), press_Q_mix(:,Q_map(EqTypeStart, 1)), &
            !                &   x_Q_ph1_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph2_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph3_mix(:,:,Q_map(EqTypeStart, 1)), x_Q_ph4_mix(:,:,Q_map(EqTypeStart, 1)),x_vap_given)
            !        write(999,'(F10.8,2X,F10.8,2X,F10.8,2X,F10.8,2X,F10.8,2X,F10.8,2X,F12.6,2X,F10.6,2X,I6)')x(1),x(2),x(3),x_ph1(1),x_ph1(2),x_ph1(3),Temp,press,errorflag
            !    enddo
            !    close(999)
            !endif
            !END   parameter study water content variation
            !#############################################
            call threephase_equil(gl, pel, Eqtype, cur_3ph, iFlash,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, errorflag)
            !in case the ratio of the guests in the gasphase is given, the overall compositions needs to be adjusted until the ratio of the guests is reached
            if ((pel%CompositionType == 2) .and. (Eqtype(1:1) == 'V') .and. (errorflag == 0)) then
                call threephase_equil_iter(gl, pel, Eqtype, cur_3ph, iFlash, Temp, press, x, x_ph1, x_ph2, x_ph3, &
                    &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, errorflag, &
                    !&   pel%phl4%Q_point_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%Qexist_mix(:,pel%phl4%Q_map(EqTypeStart, 1)),pel%phl4%Temp_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%press_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), &
                    !&   pel%phl4%x_Q_ph1_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph2_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph3_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph4_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1))&
                    & x_vap_given, 0)
            endif
            if (gl%n_guests >= 2) then
                !if no error -> save results
                if ((errorflag == 0) .and. ((dabs((rho_ph1-rho_ph2)/rho_ph1) >= 3.d-1).or.(trim(Eqtype) /= 'VLxLw'))) then !.and. (dabs((rho_ph1-rho_ph3)/rho_ph1) >= 2.d-1) .and. (dabs((rho_ph3 - rho_ph2)/rho_ph2) >= 2.d-2)) ) then
                    if (i > 1) then
                        if (dabs(pel%phl3(cur_3ph)%temp_tr_mix(i-1) - temp) / pel%phl3(cur_3ph)%temp_tr_mix(i-1) > 0.15d0 ) then
                            !write(*,*)'Large jump in calculated temperature (> 15%). Values will not be saved.'
                            call messages(gl, pel,9,trim(''),trim(''),0,0,0)
                        else
                            write_3ph = .true.
                        endif
                    else
                        write_3ph = .true.
                    endif
                    !if error maybe closed phase envelope -> go one pressure step back an increase pressure with small steps
                elseif (((errorflag /= -1111).and.(errorflag /= 0)) .or. ((dabs((rho_ph1 - rho_ph2)/rho_ph1) < 3.d-1) .and. ((rho_ph1 > 0.d0) .and. (rho_ph2 > 0.d0)) .and. ((trim(Eqtype) /= 'VLwH').or.(trim(Eqtype) /= 'LwHIw').or.(trim(Eqtype) /= 'LwHIx')))) then !.or. (dabs((rho_ph1 - rho_ph3)/rho_ph1) < 2.d-1) .or. (dabs((rho_ph3 - rho_ph2)/rho_ph2) < 2.d-2))) then
                    presstemp = press
                    if (errorflag /= 0) press = press - press_step
                    if ((dabs((rho_ph1 - rho_ph2)/rho_ph1) < 3.d-1) .or. (errorflag /= 0)) then
                        if (errorflag /= 0) errorflag = 0
                        do while ((dabs((rho_ph1 - rho_ph2)/rho_ph1) > 2.d-2) .and. (errorflag == 0))
                            press = press + 5.d-3 * press_step
                            call threephase_equil(gl, pel, Eqtype, cur_3ph, iFlash,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                                &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, errorflag)!, &
                            !&   pel%phl4%Q_point_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%Qexist_mix(:,pel%phl4%Q_map(EqTypeStart, 1)),pel%phl4%Temp_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%press_Q_mix(:,pel%phl4%Q_map(EqTypeStart, 1)), &
                            !&   pel%phl4%x_Q_ph1_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph2_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph3_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)), pel%phl4%x_Q_ph4_mix(:,:,pel%phl4%Q_map(EqTypeStart, 1)))

                            if ((dabs((rho_ph1 - rho_ph2)/rho_ph1) > 1.d0).and.(errorflag == 0).and.(trim(Eqtype) /= 'VLxH')) then
                                errorflag = -10200
                                exit
                            endif
                            if ((press_step > 0.d0).and.(presstemp <= press)) then
                                errorflag = -10200
                                exit
                            elseif ((press_step < 0.d0).and.(presstemp >= press)) then
                                errorflag = -10200
                                exit
                            endif
                        end do

                        if (errorflag == 0) then
                            if (i > 1) then
                                if (dabs(pel%phl3(cur_3ph)%temp_tr_mix(i-1) - temp) / pel%phl3(cur_3ph)%temp_tr_mix(i-1) > 0.15d0 ) then
                                    !write(*,*)'Large jump in calculated temperature (> 15%). Values will not be saved.'
                                    call messages(gl, pel,9,trim(''),trim(''),0,0,0)
                                else
                                    write_3ph = .true.
                                endif
                            else
                                write_3ph = .true.
                            endif
                        elseif (errorflag /= 0) then
                            write_3ph = .false.
                        endif
                    endif
                else
                    pel%phl3(cur_3ph)%temp_tr_mix(i) = errorflag
                    errorflag = 0
                end if
            elseif (gl%n_guests == 1) then
                if (errorflag == 0) write_3ph = .true.
            endif


            if (write_3ph .eqv. .true.) then
                ! Calculation of hydration number
                ! -------------------------------
                do J = 1,gl%N_guests
                    hyd_nr(J) = 0.d0
                    do ii = 1,gl%N_cavi
                        hyd_nr(J) = hyd_nr(J) + gl%v_cavi(ii) * occup(ii,J)
                    end do
                    hyd_nr(J) = 1.d0/hyd_nr(J)
                end do


                ! Calculation of enthalpy and entropy
                ! -----------------------------------
                ! Calculate the Chem. Pot, Fugacities and mole fract. of the phases
                IPhase = 2
                gl%molfractions = x_ph1
                call reduced_parameters_calc(gl,Temp)
                Dens = rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
                call Chempot_CALC(gl,Temp,Dens, Chempot, 0)
                call lnf_mix(gl,Temp, Dens, press, lnfi)
                fug_g(1:gl%N_guests) = dexp(lnfi(2:gl%N_hdrts)) * 1.D6
                call hdrt_mole_fract(gl, temp,press*1.D6,fug_g,occup,CiJ,xH,occup_single, occup_double)
                ! Enthalpy of hydrates
                call hdrt_enthalpy(gl, Temp, press*1.D6, xH, Chempot, fug_g, h_hyd, errorflag)
                ! Entropy of hydrates
                call hdrt_entropy(gl, temp, press*1.D6, xH, fug_g, s_hyd, errorflag)

                !Calculate the tpd of the hydrate in the other structure. Do this for "pure" hydrates only in order to check their stability
                !Andreas Jäger, November 2017
                if (gl%N_guests == 1) then
                    !Change hydrate structure to "other" structure
                    if (gl%a0_hdrt(1) < 15.D0) then !Check which structure is the natural structure for this hydrate
                        hdrt_structure_flag = 2 !Set to structure 2
                    else
                        hdrt_structure_flag = 1 !Set to structure 1
                    end if
                    call hdrt_structure_definition(gl, hdrt_structure_flag)
                    !Calculate the tpd with the changed hydrate structure along the three-phase line
                    call hdrt_chem_potent_w(gl, Temp,press*1.d6,fug_g,ChemPot_hyd)
                    tpd_sol = ChemPot_hyd - Chempot(1) !Note: for the real tpd, this would need to be multiplied by the mole fraction of water in the hydrate phase. However, as only the sign is important for deciding if the originally calculated hydrate is stable or not, the mole fraction is omitted here
                    !Save lattice parameter in other structure
                    a0_other_struct = gl%a0_hdrt(1)
                    !Set back hydrate structure
                    if (hdrt_structure_flag == 1) then
                        hdrt_structure_flag = 2
                    else
                        hdrt_structure_flag = 1
                    end if
                    call hdrt_structure_definition(gl, hdrt_structure_flag)
                end if

                ! Calculation of enthalpy of melting for pure hydrate @ VLwH and VHIw
                ! -------------------------------------------------------------------
                if ((Eqtype == 'VLwH').and.(gl%N_guests < 2)) then ! ph1 = vap, ph2 = liq2, ph3 = hyd (phX = sol)

                    !Get the enthalpy of the vapor phase
                    ! Note: there is a problem with reduced parameters
                    !       call molfractions & recuded parameters before H_CALC
                    gl%molfractions = x_ph1
                    call reduced_parameters_calc(gl,Temp)
                    h_vap = H_CALC(gl,Temp,Dens, 0)

                    !Get the enthalpy of the liquid phase
                    iPhase = 1
                    gl%molfractions = x_ph2
                    call reduced_parameters_calc(gl,Temp)
                    Dens = rhomix_calc(gl,Temp, press, 0.D0, iPhase,0)
                    h_liq = H_CALC(gl,Temp,Dens, 0)

                    x = xH
                    Phasefrac(2) = (x(1) - x_ph2(1)) / (x_ph1(1) - x_ph2(1))
                    Phasefrac(1) = 1.D0 - Phasefrac(2)

                    !Calculate the enthalpy of the fluid phases (VLw)
                    h_phase = Phasefrac(1) * h_liq + Phasefrac(2) * h_vap

                    !Calculate the enthalpy of melting for hydrates, based on the total number of moles in the
                    !fluid phase h_melt = h_melt * (1 + hyd_nr)
                    h_melt = (h_phase - h_hyd) * (1.D0 + hyd_nr(1))


                elseif ((Eqtype == 'VHIw').and.(gl%N_guests < 2)) then ! ph1 = vap, ph2 = hyd, ph3 = sol (phX = liq2)

                    !Get the enthalpy of the vapor phase
                    ! Note: there is a problem with reduced parameters
                    !       call molfractions & recuded parameters before H_CALC
                    gl%molfractions = x_ph1
                    call reduced_parameters_calc(gl,Temp)
                    h_vap = H_CALC(gl,Temp,Dens, 0)

                    !Get the enthalpy of the solid water phase
                    h_sol = h_WaterIce(gl,Temp,press)

                    !Calculate the amount of vapor and liquid in the fluid phase
                    x = xH
                    Phasefrac(2) = (x(1) - x_ph3(1)) / (x_ph1(1) - x_ph3(1))
                    Phasefrac(1) = 1.D0 - Phasefrac(2)

                    !Calculate the enthalpy of the fluid phases (VLw)
                    h_phase = Phasefrac(1) * h_sol + Phasefrac(2) * h_vap

                    !Calculate the enthalpy of melting for hydrates, based on the total number of moles in the
                    !fluid phase h_melt = h_melt * (1 + hyd_nr)
                    h_melt = (h_phase - h_hyd) * (1.D0 + hyd_nr(1))

                end if

                !saving values for python goes here
                pel%phl3(cur_3ph)%Eqtype = EqTypeIn
                pel%phl3(cur_3ph)%temp_tr_mix(i) = temp
                pel%phl3(cur_3ph)%press_tr_mix(i) = press
                pel%phl3(cur_3ph)%chempot(i,:) = chempot
                pel%phl3(cur_3ph)%lnfi(i,:) = lnfi
                pel%phl3(cur_3ph)%x_ph(i,1,:) = x_ph1
                pel%phl3(cur_3ph)%x_ph(i,2,:) = x_ph2
                pel%phl3(cur_3ph)%x_ph(i,3,:) = x_ph3
                pel%phl3(cur_3ph)%rho_tr_mix(i,1) = rho_ph1
                pel%phl3(cur_3ph)%rho_tr_mix(i,2) = rho_ph2
                pel%phl3(cur_3ph)%rho_tr_mix(i,3) = rho_ph3
                pel%phl3(cur_3ph)%occup(i,:,:) = occup
                pel%phl3(cur_3ph)%hyd_nr(i,:) = hyd_nr
                pel%phl3(cur_3ph)%s_hyd(i) = s_hyd
                pel%phl3(cur_3ph)%h_hyd(i) = h_hyd
                pel%phl3(cur_3ph)%h_melt(i) = h_melt
                pel%phl3(cur_3ph)%tpd_sol(i) = tpd_sol
                pel%phl3(cur_3ph)%a0_other_struct(i) = a0_other_struct

                write_3ph = .false.
                press = press + press_step
                if (dabs((rho_ph1 - rho_ph2)/rho_ph1) <= 5.d-2) exit
            elseif (write_3ph .eqv. .false.) then
                exit
            endif

        End do pressureloop

        if (pel%show_progress .eqv. .true.) then
            call messages(gl, pel,99,trim(''),trim(''),0,0,0)
        endif


        errorflag = 0
        ! -----------------------------------

        !k =  1
    end if



    !-------------------------------------------------------------------------------------------

    end subroutine phaselines_3
    !******************************************************************************

    !******************************************************************************
    module subroutine phaseemerge_detection(pel, EqType, phases, iQ)
    !******************************************************************************
    !
    !


    implicit none
    type(type_hdrt_pheq) :: pel
    integer, intent(in) :: phases(2), EqType
    integer, intent(out) :: iQ
    do iQ = 1,6
        if ( ((pel%phl4(Eqtype)%phaseemerge_mix(iQ,1) == phases(1)) .and. (pel%phl4(Eqtype)%phaseemerge_mix(iQ,2) == phases(2))) .or. ((pel%phl4(Eqtype)%phaseemerge_mix(iQ,1) == phases(2)) .and. (pel%phl4(Eqtype)%phaseemerge_mix(iQ,2) == phases(1))) ) exit
    enddo
    if (iQ == 7) iQ = 6

    end subroutine
    !******************************************************************************


    !******************************************************************************
    module subroutine threephase_equil_iter(gl, pel, Eqtype,cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, rho_ph1, rho_ph2, rho_ph3,  &
        & Chempot, lnfi, occup, error, &!Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4, &
        x_vap_given, bisect_mode)
    !******************************************************************************
    ! This subroutine iterates the overall composition to match the given ratio of
    ! molefractions for (atm) two components of the gasphase



    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    character(len=6), intent(in) :: Eqtype        !Type of Equilibrium, which should be calculated
    integer, intent(in) :: iTp

    !Variables at the three-phase equilibrium
    double precision :: Temp, press
    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3, x_vap_given, x_lower, x_upper
    double precision :: rho_ph1, rho_ph2, rho_ph3, x_diff_exit, x_change, x_diff_prev
    double precision, dimension(30) :: Chempot, lnfi


    integer :: cur_3ph
    integer :: error, itermax, iterations
    integer :: iPhase, iFlash, bisect_mode

    !Variables at the Quadruple point
    !character(len=6), dimension(4), intent(in)  :: Q_point
    !double precision, dimension (4), intent(in) :: press_Q, Temp_Q
    !double precision, dimension(4,30), intent(in) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: occup


    !//////////////////////////////////////////////////////////////////

    x_diff_exit = 1.3d-4                !convergence crterion in percent/100
    !x_diff_exit = 1.0d-6               !convergence crterion in percent/100
    itermax = 100                       !exit criterion
    iterations = 1                      !startvalue
    if (minval(x(1:gl%n_hdrts)) <= 5.d-4) then
        x_change = 0.99*minval(x(1:gl%n_hdrts))
    elseif (minval(x(1:gl%n_hdrts)) <= 1.d-3) then
        x_change = 0.99*minval(x(1:gl%n_hdrts))
    else
        x_change = 1.d-3                    !value by how much the molfraction of comp 1 and 2 are changed to obtain the upper limit
    endif

    call threephase_equil(gl,pel, Eqtype,cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, &
        &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
    !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)

    if ((dabs((x_ph1(2)/x_ph1(3) - x_vap_given(2)/x_vap_given(3))/(x_vap_given(2)/x_vap_given(3))) < x_diff_exit)) then!
        continue
        return
    else
        if (error == 0) then
            if (gl%n_Guests == 2) then

                if (x_ph1(2)/x_ph1(3) < x_vap_given(2)/x_vap_given(3)) then
                    x_lower = x
                    do while (x_ph1(2)/x_ph1(3) < x_vap_given(2)/x_vap_given(3))
                        x(2) = x(2) + x_change
                        x(3) = x(3) - x_change
                        if (x(3) < 1.d-12) then
                            error = -10150
                            return
                        endif
                        gl%moles_hdrt = x
                        call hdrt_ref_a_gw0_hw0_COC(gl)
                        call threephase_equil(gl,pel,Eqtype, cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                            &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
                        !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
                    end do
                    x_upper = x
                elseif (x_ph1(2)/x_ph1(3) > x_vap_given(2)/x_vap_given(3)) then
                    x_upper = x
                    do while (x_ph1(2)/x_ph1(3) > x_vap_given(2)/x_vap_given(3))
                        x(2) = x(2) - x_change
                        x(3) = x(3) + x_change
                        if (x(2) < 1.d-12) then
                            error = -10150
                            return
                        endif
                        gl%moles_hdrt = x
                        call hdrt_ref_a_gw0_hw0_COC(gl)
                        call threephase_equil(gl,pel,Eqtype,cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                            &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
                        !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
                    end do
                    x_lower = x
                endif


                bisect: do
                    if (bisect_mode == 0) then                !MB The bisect_mode determines the mixrule of x_lower and x_upper (The need of this comes from falsecodes of specific "numbers" which Jacobi_solid_NC_2P "don't like" - the coefficents are getting at some point twice as big as should be)
                        x = (x_lower + x_upper) / 2.d0
                    elseif (iterations /= 1) then
                        x = (x_lower + x_upper) / 2.d0
                    elseif (bisect_mode == 1) then
                        x = 0.25d0 * x_lower + 0.75d0 * x_upper
                    elseif (bisect_mode == 2) then
                        x = 0.75d0 * x_lower + 0.25d0 * x_upper
                    endif
                    call hdrt_ref_a_gw0_hw0_COC(gl)
                    call threephase_equil(gl,pel,Eqtype,cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                        &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
                    !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
                    !convergence criterion:
                    if ((dabs((x_ph1(2)/x_ph1(3) - x_vap_given(2)/x_vap_given(3))/(x_vap_given(2)/x_vap_given(3))) < x_diff_exit)) then
                        continue
                        exit
                        !exit criterion:
                    elseif (error /= 0) then
                        call messages(gl, pel,10,trim(''),trim(''),error,0,0)
                        exit
                    elseif (iterations > itermax) then
                        call messages(gl, pel,11,trim(''),trim(''),0,0,0)
                        pause
                        exit
                    endif
                    x_diff_prev = dabs((x_ph1(2)/x_ph1(3) - x_vap_given(2)/x_vap_given(3))/(x_vap_given(2)/x_vap_given(3)))

                    if (x_ph1(2)/x_ph1(3) < x_vap_given(2)/x_vap_given(3)) then
                        x_lower = x
                    elseif (x_ph1(2)/x_ph1(3) > x_vap_given(2)/x_vap_given(3)) then
                        x_upper = x
                    endif
                    iterations = iterations + 1
                enddo bisect
            endif
        endif
    endif

    end subroutine threephase_equil_iter
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    module subroutine phaselines_4(gl, phl4, x, x_vap_given, eqtype_ind, error)
    ! Subroutine calculates four-phase equilibrium lines and calls 'threephase_equil'




    implicit none
    type(type_gl) :: gl
    type(type_hdrt_4ph_lines) :: phl4
    type(type_hdrt_pheq), allocatable :: pel
    !type(type_splitstr), allocatable :: spltstr
    character(255) :: dummy

    !Variables at/for the four-phase equilibrium
    double precision :: Temp, press, press_start
    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3, x_vap_given
    double precision :: first_molfrac_ratio
    double precision :: phasefrac_prev
    integer :: root_dir, root_sysofeqs
    !>0 by default, shows number of found roots in the physical region
    integer :: inphysical
    logical :: abovephysical
    logical :: last_run


    double precision:: Dens, press_intermediate, press_after, press_before, press_diff_exit, mat_bal_exit, press_factor
    !1: pressure difference from one quadruple point of binary system to the other devided by 100
    !2: estimate pressure where "real" phasefraction should occur -> devide the difference from start pressur to estimated pressure of physical region in 100 parts
    !double precision, dimension(2):: press_step
    double precision :: press_step

    !>Trend error variable, loop index i, loop index k, loop index iterations, upper limit for iterations
    integer :: error, i, k, l, iterations1, iterations2, itermax, i_step
    !>Indicator for current Equilibrium type. 4-ph equilibrium depends on involved mixtrues. For CO2-CH4:
    !!1: VLwHIw, 2: VLcLwH, 3: VLcHIc 4: LwLcHI
    !!Equilibrium type is stored in Q_point_mix
    integer :: eqtype_ind!, eqtype_start = 1, eqtype_end = 1
    !>Indicates which quaduple point is used to generate start values for 4-ph equilibrium
    integer, dimension(2) :: pure_Q_flag
    !>Flash routine variables
    integer :: iPhase, iFlash, iter
    !>Write channel of the outputfile
    integer :: fourphlineunit


    !>Variables for solving mass balance
    double precision, dimension(60,60) :: mat_matbal
    double precision, dimension(60,30) :: vec_compphafra, vec_compphafra_prev
    double precision :: phafra_2ndsmall, press_prev
    integer :: root_between, roots
    !double precision ::
    character(255):: herr
    integer :: ierr, neq

    !New quadruple point routine
    double precision, dimension(4):: rho_qpts, Phasefrac
    double precision, dimension(30):: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    integer:: Phasefrac_0
    !//////////////////////////////////////////////////////////////////
    if (.not. any(phl4%Qexist_mix)) then
        error = 0
    elseif (any(phl4%Qexist_mix)) then
        dummy = ''
        press_factor = 0.1d0
        root_between = 0
        roots = 0
        inphysical = 0
        root_dir = 0
        root_sysofeqs = 0
        vec_compphafra = 0.d0
        vec_compphafra_prev = -9.d10
        last_run =  .false.
        abovephysical = .false.
        if (.not. allocated(pel)) allocate(pel)
        !allocate(spltstr)
        error = 0

        phl4%phaseemerge_mix = 0
        pel%phl4 = phl4
        if (phl4%show_progress .eqv. .true.) then
            !write(*,'(A)')' Init 4phaseline calc'
            call messages(gl, pel,400,trim(''),trim(''),0,0,0)
        endif

        if ((phl4%Qexist_mix(1) .eqv. .True.) .or. (phl4%Qexist_mix(np_4ph) .eqv. .True.)) then
            if ((any(gl%components == 'propane')) .and. (any(gl%components == 'nitrogen')))  then
                pure_Q_flag(1) = np_4ph
                pure_Q_flag(2) = 1
                i_step = 1
            elseif (phl4%Qexist_mix(1) .eqv. .True.) then
                pure_Q_flag(1) = 1
                pure_Q_flag(2) = np_4ph
                i_step = 1
            elseif (phl4%Qexist_mix(np_4ph) .eqv. .True.) then
                pure_Q_flag(1) = np_4ph
                pure_Q_flag(2) = 1
                i_step = 1
            endif




            !Set startvalues
            !default values for loop to iterate ratio of guest components
            press_diff_exit = 1.d-6             !convergence crterion
            mat_bal_exit = 1.d-2             !convergence crterion
            itermax = 500                       !break criterion
            iterations1 = 1                      !startvalue

            !Startvalues
            Temp = phl4%Temp_Q_mix(pure_Q_flag(1))
            press = phl4%press_Q_mix(pure_Q_flag(1))
            !dummy value for previous pressure
            press_prev = -9.d0
            !VLwHIw
            if (eqtype_ind == 1) then
                if ((any(gl%components == 'propane')) .and. (any(gl%components == 'nitrogen'))) then
                    press_step = 0.01d0
                else
                    press_step = 0.05d0*(phl4%press_Q_mix(pure_Q_flag(2)) - phl4%press_Q_mix(pure_Q_flag(1)))
                endif
                press_factor = 0.1d0
                x_fluid1 = phl4%x_Q_ph1_mix(:,pure_Q_flag(1))       !composition fluid phase 1
                x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))       !composition fluid phase 2
                if (pure_Q_flag(1) == 1) then
                    x_fluid1(2) = x_fluid1(2) - 0.001d0                    !change fluidvector to have a small amout of
                    x_fluid1(3) = 0.001d0                                  ! guest2 in molesfrac
                    x_fluid2(2) = x_fluid2(2) - 1.d-5                      !change fluidvector to have a small amout of
                    x_fluid2(3) = 1.d-5                                    ! guest2 in molesfrac
                elseif (pure_Q_flag(1) == np_4ph) then
                    x_fluid1(3) = x_fluid1(3) - 0.001d0                    !change fluidvector to have a small amout of
                    x_fluid1(2) = 0.001d0                                  ! guest2 in molesfrac
                    x_fluid2(3) = x_fluid2(3) - 5.d-8!1.d-5                      !change fluidvector to have a small amout of
                    x_fluid2(2) = 5.d-8!1.d-5                              ! guest2 in molesfrac
                endif
                x_fluid3 = 0.d0                                        !composition fluid phase 3
                x_fluid4 = 0.d0                                        !composition fluid phase 4
                x_sol = 0.d0                                           !initialize solid phase
                x_sol(1) = 1.d0                                        !solidphase consits only of water
                Phasefrac_0 = 2                                        !phasefrac(4) = 0 -> "VLwHIw" Lw is emerging phase
                iFlash = 1                                             !p given
                iphase = 2                                             !first fluid phase is a gasphase
                gl%solidtype(1) = 1                                       !solid phase is Iw
                gl%solidtype(2) = 1                                       !hydrate phase exists
                gl%solid_pos = 1                                          !Iw -> solid phase -> water on position 1 in components vector
                rhofluid1_est = 0.d0                                   !dummy value for density estimations
                rhofluid2_est = 0.d0                                   !dummy value for density estimations
                rhofluid3_est = 0.d0                                   !dummy value for density estimations
                rhofluid4_est = 0.d0                                   !dummy value for density estimations
                Phasefrac(1) = 0.4d0                                   !dummy value for phasefrac
                Phasefrac(2) = 0.4d0                                   !dummy value for phasefrac
                Phasefrac(3) = 0.2d0                                   !dummy value for phasefrac
                Phasefrac(4) = 0.d0                                    !Iw is emerging phase

                !VLwLxH
            elseif (eqtype_ind == 2) then
                !Pressure step
                !press_step = 1.d-3
                if ((phl4%Qexist_mix(1) .eqv. .True.) .and. (phl4%Qexist_mix(np_4ph)) .eqv. .True.) then
                    press_step = (phl4%press_Q_mix(np_4ph) - phl4%press_Q_mix(1))/49
                else
                    if (trim(gl%components(2)) == 'co') then
                        press_step = 1.d-1
                    elseif ((any(gl%components == 'propane')) .and. (any(gl%components == 'oxygen'))) then
                        press_step = 2.d-1
                    elseif ((any(gl%components == 'propane')) .and. (any(gl%components == 'nitrogen'))) then
                        press_step = 1.5d-1
                    else
                        press_step = 5.d-2
                    endif
                endif

                if ((any(gl%components == 'propane')) .and. (any(gl%components == 'nitrogen'))) then
                    press_factor = 1.d0
                else
                    press_factor = 0.4d0
                endif
                x_fluid1 = phl4%x_Q_ph1_mix(:,pure_Q_flag(1))       !composition fluid phase 1
                if (pure_Q_flag(1) == 1) then
                    !vapor phase
                    x_fluid1(2) = x_fluid1(2) - 0.001d0                !change fluidvector to have a small amout of
                    x_fluid1(3) = 0.001d0                              ! guest2 in molesfrac
                    !liquid water like phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(2) = x_fluid2(2) - 1.d-6   !change fluidvector to have a small amout of
                    x_fluid2(3) = 1.d-6                 ! guest2 in molesfrac
                    !liquid x phase
                    x_fluid3 = phl4%x_Q_ph3_mix(:,pure_Q_flag(1))       !composition fluid phase 2
                    x_fluid3(2) = x_fluid3(2) - 1.d-2   !change fluidvector to have a small amout of
                    x_fluid3(3) = 1.d-2                 ! guest2 in molesfrac
                    x_fluid4 = 0.d0                     !composition fluid phase 4
                elseif (pure_Q_flag(1) == np_4ph) then
                    !vapor phase
                    x_fluid1(3) = x_fluid1(3) - 0.001d0                !change fluidvector to have a small amout of
                    x_fluid1(2) = 0.001d0                              ! guest2 in molesfrac
                    !liquid water like phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(3) = x_fluid2(3) - 1.d-6   !change fluidvector to have a small amout of
                    x_fluid2(2) = 1.d-6                 ! guest2 in molesfrac
                    !liquid x phase
                    x_fluid3 = phl4%x_Q_ph3_mix(:,pure_Q_flag(1))       !composition fluid phase 2
                    x_fluid3(3) = x_fluid3(3) - 1.d-2   !change fluidvector to have a small amout of
                    x_fluid3(2) = 1.d-2                 ! guest2 in molesfrac
                    x_fluid4 = 0.d0                     !composition fluid phase 4
                endif
                x_sol = 0.d0                        !initialize solid phase
                !x_sol(1) = 1.d0                     !solidphase consits only of water
                !Phasefrac_0 = 1                     !phasefrac(1) = 0 -> "VLwLcH" V is emerging phase
                !Phasefrac_0 = 2                     !phasefrac(2) = 0 -> "VLwLcH" Lw is emerging phase
                !Phasefrac_0 = 3                     !phasefrac(3) = 0 -> "VLwLcH" Lc is emerging phase
                Phasefrac_0 = 4                     !phasefrac(4) = 0 -> "VLwLcH" H is emerging phase
                if ((any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                    iFlash = 2                          !T given
                else
                    iFlash = 1                          !p given
                endif
                iphase = 2                          !first fluid phase is a gasphase
                gl%solidtype(1) = 0                    !solid phase does not exist
                gl%solidtype(2) = 1                    !hydrate phase exists
                gl%solid_pos = 0                       !Iw does not exist
                rhofluid1_est = 0.d0                !dummy value for density estimations
                rhofluid2_est = 0.d0                !dummy value for density estimations
                rhofluid3_est = 0.d0                !dummy value for density estimations
                rhofluid4_est = 0.d0                !dummy value for density estimations
                Phasefrac(1) = 0.8d0                !dummy value for phasefrac
                Phasefrac(2) = 0.1d0                !dummy value for phasefrac
                Phasefrac(3) = 0.1d0                !Lc is emerging phase
                Phasefrac(4) = 0.0d0                !dummy value for phasefrac

                !VLcHIc
            elseif (eqtype_ind == 3) then
                press_diff_exit = 1.d-5 !slighty lower convergence criterion, higher criterion lead to errors
                press_factor = 0.8d0  !this fourphase line is rather large
                if ((trim(gl%components(2)) == 'co2') .and. (trim(gl%components(3)) == 'nitrogen')) then
                    press_step = 2.d-1
                else
                    !Pressure step
                    press_step = 5.d-2
                endif

                x_fluid1 = phl4%x_Q_ph1_mix(:,pure_Q_flag(1))       !composition fluid phase 1
                if (pure_Q_flag(1) == 1) then
                    !vapor phase
                    x_fluid1(2) = x_fluid1(2) - 0.1d0                !change fluidvector to have a small amout of
                    x_fluid1(3) = 0.1d0                              ! guest2 in molesfrac
                    !liquid guest1 rich phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(2) = x_fluid2(2) - 1.d-2   !change fluidvector to have a small amout of
                    x_fluid2(3) = 1.d-2                 ! guest2 in molesfrac

                elseif (pure_Q_flag(1) == np_4ph) then
                    !vapor phase
                    x_fluid1(3) = x_fluid1(3) - 0.1d0                !change fluidvector to have a small amout of
                    x_fluid1(2) = 0.1d0                              ! guest2 in molesfrac
                    !liquid guest1 rich phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(3) = x_fluid2(3) - 1.d-2   !change fluidvector to have a small amout of
                    x_fluid2(2) = 1.d-2                 ! guest2 in molesfrac
                endif
                x_fluid3 = 0.d0                                        !composition fluid phase 3
                x_fluid4 = 0.d0                                        !composition fluid phase 4

                x_sol = 0.d0                        !initialize solid phase
                x_sol(2) = 1.d0                     !solidphase consits only of co2
                !Phasefrac_0 = 1                     !phasefrac(1) = 0 -> "VLcHIc" V is emerging phase
                Phasefrac_0 = 2                     !phasefrac(2) = 0 -> "VLcHIc" Lx is emerging phase
                !Phasefrac_0 = 3                     !phasefrac(3) = 0 -> "VLcHIc" H is emerging phase
                !Phasefrac_0 = 4                     !phasefrac(4) = 0 -> "VLcHIc" Ic is emerging phase
                if ((any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                    iFlash = 2                          !T given
                else
                    iFlash = 1                          !p given
                endif
                iphase = 2                          !first fluid phase is a gasphase
                gl%solidtype(1) = 2                    !solid phase with co2
                gl%solidtype(2) = 1                    !hydrate phase exists
                gl%solid_pos = 2                       !Ic -> solid phase -> co2 on position 2 in components vector
                rhofluid1_est = 0.d0                !dummy value for density estimations
                rhofluid2_est = 0.d0                !dummy value for density estimations
                rhofluid3_est = 0.d0                !dummy value for density estimations
                rhofluid4_est = 0.d0                !dummy value for density estimations
                Phasefrac(1) = 0.8d0                !dummy value for phasefrac
                Phasefrac(2) = 0.1d0                !dummy value for phasefrac
                Phasefrac(3) = 0.1d0                !Lc is emerging phase
                Phasefrac(4) = 0.0d0                !dummy value for phasefrac

                !LxLwHIy
            elseif (eqtype_ind == 4) then
                press_diff_exit = 1.d-3 !slighty lower convergence criterion, higher criterion lead to errors
                !Pressure step
                press_step = 1.d1
                press_factor = 1.d0

                x_fluid1 = phl4%x_Q_ph1_mix(:,pure_Q_flag(1))       !composition fluid phase 1
                if (pure_Q_flag(1) == 1) then
                    !liquid water rich phase
                    x_fluid1(2) = x_fluid1(2) - 1.d-2                !change fluidvector to have a small amout of
                    x_fluid1(3) = 1.d-2                              ! guest2 in molesfrac
                    !liquid guest rich phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(2) = x_fluid2(2) - 1.d-2   !change fluidvector to have a small amout of
                    x_fluid2(3) = 1.d-2                 ! guest2 in molesfrac

                elseif (pure_Q_flag(1) == np_4ph) then
                    !liquid water rich phase
                    x_fluid1(2) = 0.01*x_fluid1(3)!1.d-2                              ! guest2 in molesfrac
                    x_fluid1(3) = 0.99*x_fluid1(3) !- 1.d-2                !change fluidvector to have a small amout of
                    x_fluid1 = x_fluid1/sum(x_fluid1)
                    !liquid guest rich phase
                    x_fluid2 = phl4%x_Q_ph2_mix(:,pure_Q_flag(1))   !composition fluid phase 2
                    x_fluid2(2) = 0.01*x_fluid2(3)!1.d-2                 ! guest2 in molesfrac
                    x_fluid2(3) = 0.99*x_fluid2(3) !- 1.d-2   !change fluidvector to have a small amout of
                    x_fluid2 = x_fluid2/sum(x_fluid2)
                endif
                x_fluid3 = 0.d0                                        !composition fluid phase 3
                x_fluid4 = 0.d0                                        !composition fluid phase 4

                x_sol = 0.d0                        !initialize solid phase
                x_sol(2) = 1.d0                     !solidphase consits only of co2
                !Phasefrac_0 = 1                     !phasefrac(1) = 0 -> "LxLwHIy" Lw is emerging phase
                !Phasefrac_0 = 2                     !phasefrac(2) = 0 -> "LxLwHIy" Lx is emerging phase
                Phasefrac_0 = 3                     !phasefrac(3) = 0 -> "LxLwHIy" H is emerging phase
                !Phasefrac_0 = 4                     !phasefrac(4) = 0 -> "LxLwHIy" Iy is emerging phase
                iFlash = 1                          !p given
                iphase = 2                          !first fluid phase is a gasphase
                gl%solidtype(1) = 2                    !solid phase with co2
                gl%solidtype(2) = 1                    !hydrate phase exists
                gl%solid_pos = 2                       !Ic -> solid phase -> co2 on position 2 in components vector
                rhofluid1_est = 0.d0                !dummy value for density estimations
                rhofluid2_est = 0.d0                !dummy value for density estimations
                rhofluid3_est = 0.d0                !dummy value for density estimations
                rhofluid4_est = 0.d0                !dummy value for density estimations
                Phasefrac(1) = 0.8d0                !dummy value for phasefrac
                Phasefrac(2) = 0.1d0                !dummy value for phasefrac
                Phasefrac(3) = 0.1d0                !Lc is emerging phase
                Phasefrac(4) = 0.0d0                !dummy value for phasefrac

            endif

            !if overall composition is given the "given" ratio of guests in the vapor phase has to be generated with the given values
            !to generate better start values for the three phase lines
            if (pel%CompositionType == 1) then
                x_vap_given(2) = x(2)/x(1)
                x_vap_given(3) = x(3)/x(1)
            endif
            !calculate first massbalance with half press_step
            press_start = phl4%press_Q_mix(pure_Q_flag(1))

            !startvalues for finding points with two phasefractions = 0
            vec_compphafra = -9.d10; vec_compphafra_prev = -9.d10
            root_dir = 0 ! positive if root is at higher pressure and vice versa
            root_between = 0 !root is between last and current pressure step
            inphysical = 0 !is zero below first root of physical reasonable phasefractions and the value is number of found roots
            abovephysical = .false. !if region of physical reasonable phasefractions already passed
            press_intermediate = 0.d0

            !calculate one point with three components to get start values for the compositions and phasefractions
            !press = press_start + press_step * press_factor!0.5d0
            if ((any(gl%components == 'propane')) .and. (any(gl%components == 'oxygen'))) then
                press = press_start + press_step
            elseif ((any(gl%components == 'ethane')) .and. (any(gl%components == 'co2'))) then
                if (iflash == 1) then
                    press = press_start + press_step * 0.05d0
                elseif (iflash == 2) then
                    if (eqtype_ind == 2) then
                        temp = temp + 0.1d0
                    elseif (eqtype_ind == 3) then
                        temp = temp - 0.1d0
                    endif
                endif
            elseif ((any(gl%components == 'co')) .and. (any(gl%components == 'co2'))) then
                if (iflash == 1) then
                    press = press_start + press_step * 0.05d0
                elseif (iflash == 2) then
                    if (eqtype_ind == 2) then
                        temp = temp + 0.1d0
                    elseif (eqtype_ind == 3) then
                        temp = temp - 0.5d0
                    endif
                endif
            else
                press = press_start + press_step * 0.1d0
            endif
            call hdrt_ref_a_gw0_hw0_COC(gl)
            call ptflash_solid_NC_4P(gl,press, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)
            i = 1
            do while ((error /= 0).and.(iflash == 1))
                !press = press_start + press_step * press_factor!0.5d0
                press = press_start + press_step * 0.1d0
                call ptflash_solid_NC_4P(gl,press, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                    & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)
                i = i + 1
                if (i > itermax) then
                    if (eqtype_ind == 4) then
                        phl4%Qexist_mix(pure_q_flag(1)) = .false.
                        error = 0
                        return
                    else
                        return
                    endif
                endif
            enddo
            call massbalance_4ph(gl,press, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysical + 1, root_dir, root_between, root_sysofeqs)
            if (iflash == 1) then
                press_prev = press
                press = press_start
            elseif (iflash == 2) then
                if ((eqtype_ind == 2).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                    temp = temp + 0.1d0
                endif
            endif
            i = 1
            fourphaselineloop: do
                if ((phl4%fourphase_physical .eqv. .false.).or.(inphysical /= 0)) i = i + i_step
                if ((i > np_4ph).or.(i == 0)) exit
                if ((phl4%fourphase_physical .eqv. .true.) .and. (abovephysical .eqv. .true.)) exit
                root_between = 0
                root_dir = 0
                if (iflash == 1) then
                    !before physical region: inphysical = 0
                    !in physical region: 0 <= phasefrac <= 1 (do smaller steps)
                    !after physical region: inphysical /= 0 minval(phasefrac) < 0
                    if (inphysical == 0) then
                        if ((eqtype_ind == 1).and.(any(gl%components == 'propane')) .and. (any(gl%components == 'oxygen')) .and. press < 1.5d0) then
                            press = press + press_step * 0.1d0
                        elseif ((eqtype_ind == 3).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                            press = press + press_step * 0.1d0
                        elseif ((eqtype_ind == 1).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'co'))) then
                            press = press + press_step * press_factor
                        else
                            press = press + press_step
                        endif
                    elseif ((inphysical /= 0) .and. (abovephysical .eqv. .false.)) then
                        !press = press + press_step(2) * press_factor
                        press = press + press_step * press_factor
                    elseif ((inphysical /= 0) .and. (abovephysical .eqv. .true.)) then
                        !press = press + press_step(1)
                        press = press + press_step
                    endif
                elseif (iflash == 2) then
                    if ((eqtype_ind == 3).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'methane'))) then
                        temp = temp - 1.d0
                    elseif ((eqtype_ind == 2).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                        temp = temp + 0.05d0
                        if (dabs(phl4%Temp_Q_mix(i-1) - phl4%Temp_Q_mix(200)) < 0.5d0) then
                            iflash = 1
                            press_step = -press_step
                            i = i - i_step
                            cycle
                        endif
                    elseif ((eqtype_ind == 3).and.(any(gl%components == 'co2')) .and. (any(gl%components == 'ethane'))) then
                        temp = temp - 0.1d0
                    endif
                    if (i > 2) press_step =  press - phl4%press_Q_mix(i-2)
                endif

                !EXIT CRITERIONS - start
                !check if the Qpoint for the binary mixtures water + guest2 exists and if the current pressure is above the pressure for the binary system water + guest2
                !-> this is the final state point and is written to the Q matrices
                if (phl4%Qexist_mix(pure_Q_flag(2))) then
                    !check if the current pressure is above the pressure of binary mixture water + guest2 if press_step is positive
                    if ( (iflash == 1).and.( ((press_step > 0.d0) .and. (press > phl4%press_Q_mix(pure_Q_flag(2)))) .or. ((press_step < 0.d0) .and. (press < phl4%press_Q_mix( pure_Q_flag(2)))) ).and.(i < np_4ph)) then
                        call phaselines_4_move_q_point(phl4, i, pure_Q_flag(2))
                        !write pure quadruple point directly after mixed hydrate data
                        if (abovephysical .eqv. .false.) then
                            if (inphysical == 0) then
                                continue
                            elseif (inphysical == 1) then
                                if ((phl4%phaseemerge_mix(inphysical,1) == 2) .and. (phl4%phaseemerge_mix(inphysical,2) == 3)) then
                                    phl4%phaseemerge_mix(2,1) = 3
                                    phl4%phaseemerge_mix(2,2) = 4
                                    phl4%phaseemerge_mix(3,1) = 2
                                    phl4%phaseemerge_mix(3,2) = 4
                                elseif ((phl4%phaseemerge_mix(inphysical,1) == 2) .and. (phl4%phaseemerge_mix(inphysical,2) == 4)) then
                                    phl4%phaseemerge_mix(2,1) = 3
                                    phl4%phaseemerge_mix(2,2) = 4
                                    phl4%phaseemerge_mix(3,1) = 2
                                    phl4%phaseemerge_mix(3,2) = 3
                                endif
                                phl4%Q_map(2) = phl4%Q_map(inphysical)
                                phl4%Q_map(3) = phl4%Q_map(inphysical)
                            elseif (inphysical == 2) then
                                if ((phl4%phaseemerge_mix(inphysical-1,1) == 2) .and. (phl4%phaseemerge_mix(inphysical-1,2) == 3)) then
                                    phl4%phaseemerge_mix(3,1) = 2
                                    phl4%phaseemerge_mix(3,2) = 4
                                elseif ((phl4%phaseemerge_mix(inphysical-1,1) == 2) .and. (phl4%phaseemerge_mix(inphysical-1,2) == 4)) then
                                    phl4%phaseemerge_mix(3,1) = 2
                                    phl4%phaseemerge_mix(3,2) = 3
                                endif
                                phl4%Q_map(3) = phl4%Q_map(inphysical)
                            endif
                        endif
                        last_run = .true.
                        i = i + i_step
                        exit
                    endif
                    !if VLLH fourphaseline look for density difference, if small -> exit
                elseif ((eqtype_ind == 2).and.(phl4%fourphase_physical .eqv. .false.)) then
                    if (dabs((rho_qpts(1) - rho_qpts(3))/rho_qpts(3)) < 0.2d0) then
                        !if (dabs((rho_qpts(1) - rho_qpts(3))/rho_qpts(3)) < 0.05d0) then
                        !write(*,*)'Densities of Vapor and Guest rich phase very close. Calculation ends.'
                        call messages(gl, pel,401,trim(''),trim(''),0,0,0)
                        !when the quadruple point of guest#2 was used for startvalues of 4-phaseline it has to be moved to 1st place in module variables
                        if (pure_Q_flag(1) == np_4ph) then
                            call phaselines_4_move_q_point(phl4, 1, pure_q_flag(1))
                        endif
                        exit
                    endif
                elseif ((eqtype_ind == 3).and.(iflash == 1)) then
                    if (x_fluid1(1) < 3.d-7 ) then
                        iflash = 2
                        if (i /= 0) i = i - i_step
                        cycle
                    elseif ((eqtype_ind == 3).and.(dabs(rho_qpts(1) - rho_qpts(2))/rho_qpts(2) < 0.2d0)) then
                        call messages(gl, pel,401,trim(''),trim(''),0,0,0)
                        if (pure_Q_flag(1) == np_4ph) then
                            call phaselines_4_move_q_point(phl4, i, pure_q_flag(1))
                        endif
                        exit
                    endif
                elseif ((eqtype_ind == 3).and.(iflash == 2)) then
                    if (x_fluid1(1) < 1.d-17 ) then
                        if (pure_Q_flag(1) == np_4ph) then
                            call phaselines_4_move_q_point(phl4, i, pure_Q_flag(1))
                        endif
                        exit
                    endif
                elseif (eqtype_ind == 4) then
                    !if amount of guest#2 in water rich phase is very low
                    if (x_fluid1(3) < 1.d-10) then
                        i = i - i_step
                        exit
                        !if the densities of guest rich and dry ice phase come close to each other
                    elseif (dabs(rho_qpts(2)-rho_qpts(4))/rho_qpts(2) < 2.d-2) then
                        i = i - i_step
                        exit
                    endif
                    !if only physical reasonable phasefractions should be saved and current pressure is above this region
                elseif ((abovephysical .eqv. .true.) .and. (phl4%fourphase_physical .eqv. .true.)) then
                    exit
                endif
                !EXIT CRITERIONS - end
                if (pel%show_progress .eqv. .true.) then
                    !write(*,'(''+'',A,I5)')Q_point_mix(Eqtype_ind,pure_Q_flag(1)) // ': calculating point ',i!, ' of ',100
                    call messages(gl, pel,402,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),i,0,0)
                elseif (pel%show_progress .eqv. .false.) then
                    if (i == 2) call messages(gl, pel,403,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),0,0,0)!write(*,*)'Calculation of ' // trim(Q_point_mix(Eqtype_ind,pure_Q_flag(1)))
                endif

                call hdrt_ref_a_gw0_hw0_COC(gl)

                !if fourphase line with physical reasonable values for phasefraction needs to be calculated boundaries have to be determined
                call ptflash_solid_NC_4P(gl,press, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                    & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)!               inphysical + 1, change in phasefraction sign needs to be checked at the conditions for nex root
                if (error /= 0) exit
                if (abovephysical .eqv. .false.) call massbalance_4ph(gl,press, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysical + 1, root_dir, root_between, root_sysofeqs)
                press_prev = press
                !ITERATIONS FOR PHYSICAL REASONABLE REGION OF FOURPHASELINES - START
                !if already in the physical region first pressure step above this region will result in unphysical betas ->
                ! beta < 0 or beta > 1 calculation is stopped when only the fourphase line with physical correct betas needs to be calculated (fourphase_physical .eqv. .true.)
                !iterate pressure at root:
                if ((root_between > 0).and.(abovephysical .eqv. .false.)) then
                    if (iflash == 1) then
                        roots = root_between
                        inphysical = inphysical + 1
                        if (pel%show_progress .eqv. .true.) call messages(gl, pel,404,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),i,inphysical,0)! write(*,'(''+'',A,A,I5,A,I2)')Q_point_mix(Eqtype_ind,pure_Q_flag(1)),': calculating point ', i, ' Finding physical solutions',inphysical

                        !set boundaries
                        press_after = press
                        !set the  lower limit slightly above the previous root
                        press_before = press - press_step
                        !if two roots are >>very<< close together and the second root is currently between current and last step,
                        !it could happen that press_step(2) > press-press_intermediate (press_intermediate has the value with
                        !the last iteratet pressure at root) in that case the first root is between press_before and press_after
                        if (press_step > 0) then
                            if (press_before < press_intermediate) press_before = press_intermediate
                        elseif ((press_step < 0) .and. (press_intermediate > 1.d-14) )then
                            if (press_before > press_intermediate) press_before = press_intermediate
                        endif
                        !if in the region of physical reasonable values for phasefraction (0 <= beta <= 1) iterate the roots between previous and current pressure step if present
                        !if the first root (inphysical = 1) is determined after finding a root pressure is decreased to check whether the first root was really found.
                        iterations1 = 0
                        do
                            iterations1 = iterations1 + 1
                            if (pel%show_progress .eqv. .true.) call messages(gl, pel,405,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),i,inphysical,iterations1)!write(*,'(''+'',A,A,I5,A,I2,A,I3)')Q_point_mix(Eqtype_ind,pure_Q_flag(1)),': calculating point ', i, ' Finding physical solutions',inphysical,' Iteration',iterations1
                            if (iterations1 > itermax) then
                                error = -10000; exit
                            endif
                            if (root_between <= 1) then
                                press_intermediate = 0.5d0 * (press_after + press_before)
                            elseif (root_between >= 2) then
                                press_intermediate = 0.99d0 * press_before + 0.01d0 * press_after
                            else
                                !write(*,*)'More than two roots found, code needs to be updated!'
                                call messages(gl, pel,406,trim(''),trim(''),0,0,0)
                            endif
                            root_between = 0
                            call ptflash_solid_NC_4P(gl,press_intermediate, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                                & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)
                            if (abovephysical .eqv. .false.) call massbalance_4ph(gl,press_intermediate, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysical, root_dir, root_between, root_sysofeqs)
                            press_prev = press_intermediate
                            if (root_dir == -1) then
                                press_after = press_intermediate
                            elseif (root_dir == 1) then
                                press_before = press_intermediate
                            endif
                            !exit criterion: pressure difference of upper and lower boundary small and 0 <= beta <= 1
                            if ((dabs(press_after - press_before)/press_intermediate < press_diff_exit).and.((minval(vec_compphafra(:,1)) >= -1.d-14).or.(minval(vec_compphafra(:,2)) >= -1.d-14).or.(minval(vec_compphafra(:,3)) >= -1.d-14).or.(minval(vec_compphafra(:,4)) >= -1.d-14))) then
                                !if there is a root with a lower pressure press_intermediate should be smaller then press upper -> search again for root with lower boundaries
                                if ((dabs(press_intermediate - press_after))/press_intermediate > press_diff_exit) then
                                    press_after = press_intermediate
                                    press_before = press_intermediate - press_step * 0.5d0
                                    !otherwise: finally the root was found
                                else
                                    if (inphysical > 1) then
                                        !check if the same root was found in the previous step
                                        if (dabs((phl4%press_Q_mix(phl4%Q_map( inphysical - 1 )) - press_intermediate)/press_intermediate) < 1.d-6) then !same root was found
                                            inphysical = inphysical - 1
                                            !if (fourphase_physical .eqv. .true.) i = i - 1
                                            if (pel%show_progress .eqv. .true.)  call messages(gl, pel,407,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),i,inphysical,0)!write(*,'(''+'',A,A,I5,A,I2,A)')Q_point_mix(Eqtype_ind,pure_Q_flag(1)),': calculating point ', i, ' Found',inphysical ,' physical solutions    '
                                            call ptflash_solid_NC_4P(gl,press, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                                                & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)
                                            !abovephysical = .true.
                                            root_between = 0
                                            exit
                                        endif
                                    endif
                                    !check if after physical reasonable region  (press_intermediate+press_step(2)*press_factor)
                                    if ( (eqtype_ind == 1) .or. (eqtype_ind == 3) ) then
                                        if (inphysical >= 3) abovephysical = .true.
                                    elseif ( (eqtype_ind == 2) .or. (eqtype_ind == 4) ) then
                                        if (inphysical >= 4) abovephysical = .true.
                                    endif

                                    !determine the two emerging phases on the four phase line
                                    !vec_compphafra contains the phasefractions depending on which phase is set as an emerging phase
                                    !root_sysofeqs indicates the system of equations (see massbalance_4ph) where the change in the sign of the respective phasefraction occured
                                    !              this system contains the two emerging phasefractions and will be used later to generate startvalues for threephaselines
                                    !1st emerging phase is the smallest value in the array, where the sign change occured
                                    !2nd emerging phase is the second smallest value
                                    phl4%phaseemerge_mix(inphysical,1) = minloc(vec_compphafra(:,root_sysofeqs),1)
                                    phafra_2ndsmall = 1.d9
                                    do k = 1,4
                                        if ((phafra_2ndsmall > vec_compphafra(k,root_sysofeqs)).and.(k /= phl4%phaseemerge_mix(inphysical,1))) then
                                            phafra_2ndsmall = vec_compphafra(k,root_sysofeqs)
                                            phl4%phaseemerge_mix(inphysical,2) = k
                                        endif
                                    enddo
                                    !in case the root is in positve pressure direction push the pressure so far, that the sign of the relevant phasefraction changes
                                    !this has to be done, because otherwise in the next step the same root will be found
                                    if (root_dir == 1) then
                                        root_between = 0
                                        press = press_intermediate
                                        do
                                            !press = press + press_step(2) * 1.d-3
                                            if (eqtype_ind /= 4) then
                                                press = press + press_step * 1.d-6
                                            elseif (eqtype_ind == 4) then
                                                press = press + press_step
                                            endif
                                            call ptflash_solid_NC_4P(gl,press, Temp, x, rho_qpts, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
                                                & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, error, iter)
                                            call massbalance_4ph(gl,press, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysical, root_dir, root_between, root_sysofeqs)
                                            if (root_dir == -1) exit
                                        enddo
                                    endif
                                    !backspace(200)
                                    !write(200,'(I3,A1,F16.8,A1,F16.8,A1,F16.8,A1,F16.8,A2,F16.8,A1,F16.8,A1,F16.8,A1,F16.8,A2,F16.8,A1,F16.8,A1,F16.8,A1,F16.8,A2,F16.8,A1,F16.8,A1,F16.8,A1,F16.8,A2)')i,",",vec_compphafra(1,1),",",vec_compphafra(2,1),",",vec_compphafra(3,1),",",vec_compphafra(4,1),",,",vec_compphafra(1,2),",",vec_compphafra(2,2),",",vec_compphafra(3,2),",",vec_compphafra(4,2),",,",vec_compphafra(1,3),",",vec_compphafra(2,3),",",vec_compphafra(3,3),",",vec_compphafra(4,3),",,",vec_compphafra(1,4),",",vec_compphafra(2,4),",",vec_compphafra(3,4),",",vec_compphafra(4,4),",*"
                                    phl4%Q_map(inphysical) = i
                                    press = press_intermediate! + press_step(2) * press_factor
                                    exit
                                endif
                            endif
                        enddo
                    elseif (iflash == 2) then
                        inphysical = inphysical + 1
                        if (pel%show_progress .eqv. .true.) call messages(gl, pel,404,trim(phl4%Q_point_mix(pure_Q_flag(1))) // dummy,trim(''),i,inphysical,0)! write(*,'(''+'',A,A,I5,A,I2)')Q_point_mix(Eqtype_ind,pure_Q_flag(1)),': calculating point ', i, ' Finding physical solutions',inphysical
                        phl4%Q_map(inphysical) = i
                        root_between = 0
                        phl4%phaseemerge_mix(inphysical,1) = minloc(vec_compphafra(:,root_sysofeqs),1)
                        phafra_2ndsmall = 1.d9
                        do k = 1,4
                            if ((phafra_2ndsmall > vec_compphafra(k,root_sysofeqs)).and.(k /= phl4%phaseemerge_mix(inphysical,1))) then
                                phafra_2ndsmall = vec_compphafra(k,root_sysofeqs)
                                phl4%phaseemerge_mix(inphysical,2) = k
                            endif
                        enddo
                        if (eqtype_ind < 4) then
                            if ( ((minval(vec_compphafra(:,1)) < 0.d0) .or. (maxval(vec_compphafra(:,1)) > 1.d0)) .and. ((minval(vec_compphafra(:,2)) < 0.d0) .or. (maxval(vec_compphafra(:,2)) > 1.d0))) abovephysical = .true.
                        elseif (eqtype_ind == 4) then
                            if ( ((minval(vec_compphafra(:,1)) < 0.d0) .or. (maxval(vec_compphafra(:,1)) > 1.d0)) .and. ((minval(vec_compphafra(:,2)) < 0.d0) .or. (maxval(vec_compphafra(:,2)) > 1.d0)) &
                                & .and. ((minval(vec_compphafra(:,3)) < 0.d0) .or. (maxval(vec_compphafra(:,3)) > 1.d0)) ) abovephysical = .true.
                        endif
                    endif
                endif
                !ITERATIONS FOR PHYSICAL REASONABLE REGION OF FOURPHASELINES - END

                !((if only physical reasonable phasefractions have to be saved) and (pressure is in this range)) or (complete fourphaseline is calculated)
                !-> save data
                if ((error == 0).and.(i > 0)) then
                    phl4%Q_point_mix(i) = phl4%Q_point_mix(pure_Q_flag(1))
                    phl4%Qexist_mix(i) = phl4%Qexist_mix(pure_Q_flag(1))
                    phl4%Temp_Q_mix(i) = Temp
                    phl4%rho_q_mix(:,i) = rho_qpts
                    if (root_between > 0) then
                        phl4%press_Q_mix(i) = press_intermediate
                        phl4%phasefrac_mix(:,i) = vec_compphafra(1:4,root_sysofeqs)
                    else
                        phl4%press_Q_mix(i) = press
                        phl4%phasefrac_mix(:,i) = phasefrac
                    endif

                    if ((eqtype_ind == 1).or.(eqtype_ind == 3).or.(eqtype_ind == 4)) then
                        phl4%x_Q_ph1_mix(:,i) = x_fluid1
                        phl4%x_Q_ph2_mix(:,i) = x_fluid2
                        phl4%x_Q_ph3_mix(:,i) = x_hyd
                        phl4%x_Q_ph4_mix(:,i) = x_sol
                    elseif (eqtype_ind == 2) then
                        phl4%x_Q_ph1_mix(:,i) = x_fluid1
                        phl4%x_Q_ph2_mix(:,i) = x_fluid2
                        phl4%x_Q_ph3_mix(:,i) = x_fluid3
                        phl4%x_Q_ph4_mix(:,i) = x_hyd
                    endif

                    if (i == 2) first_molfrac_ratio = x_fluid1(2)/x_fluid1(3) !save first ratio of guests: if > 1 more of guest 1 at the beginning

                elseif (error /= 0) then
                    !write(*,*)
                    !write(*,*)'Convergence Problem in fourphaselines in ptflash_sol_NC_4P.'
                    !write(*,*)'It is very likely that there exists no Quadruple point for the binary system of'
                    !write(*,*)'water and one guest. Calculation of threephase lines starting from ' // Q_point_mix(Eqtype_ind,pure_Q_flag(1)) // 'line'
                    !write(*,*)'maybe is not possible. Molfractions of Guests from last successfully calculated'
                    !write(*,*)'point on 4-ph-line:'
                    if (pel%CompositionType == 1) then
                        !write(*,*)x_Q_ph1_mix(eqtype_ind,2:3,i-1)*(1.d0-x(1))
                        call messages(gl, pel,408,trim(''),trim(''),eqtype_ind,i,pure_Q_flag(1))
                    elseif( pel%CompositionType == 2) then
                        !write(*,*)x_Q_ph1_mix(eqtype_ind,2:3,i-1)/(1.d0-x_Q_ph1_mix(eqtype_ind,2:3,i-1))
                        call messages(gl, pel,409,trim(''),trim(''),eqtype_ind,i,pure_Q_flag(1))
                    endif
                    !write(*,*)
                    if (abovephysical .eqv. .false.) phl4%Q_map(1) = i - i_step
                    if (eqtype_ind == 1) then
                        error = -10100 ! for VLHI there exists a quadruple point for co2-methane -> force calculation for other fourphaselines until error, but continue threephaseline calculation
                    else
                        error = 0
                    endif

                    exit
                endif
                press_prev = phl4%press_Q_mix(i)
            enddo fourphaselineloop
            if (pel%show_progress .eqv. .true.) then
                !write(*,*)
                call messages(gl, pel,99,trim(''),trim(''),0,0,0)
            endif

            if (error /= 0) then
                i = i - i_step
            endif

        else
            if (eqtype_ind == 1) then
                call messages(gl, pel,411,trim(''),trim(''),0,0,0)
                !write(*,*)'Quadruple point VLwHIw not present for both hydrate former.'
            elseif (eqtype_ind == 2) then
                call messages(gl, pel,412,trim(''),trim(''),0,0,0)
                !write(*,*)'Quadruple point VLwLxH not present for both hydrate former.'
            elseif (eqtype_ind == 3) then
                call messages(gl, pel,413,trim(''),trim(''),0,0,0)
                !write(*,*)'Quadruple point VLxHIc not present for both hydrate former.'
            elseif (eqtype_ind == 4) then
                call messages(gl, pel,414,trim(''),trim(''),0,0,0)
                !write(*,*)'Quadruple point LxLyHIc not present for both hydrate former.'
            endif
        endif
        if (pel%show_progress .eqv. .true.) then
            !write(*,*)
            call messages(gl, pel,99,trim(''),trim(''),0,0,0)
        endif
        !enddo Eqtypeloop
        continue
        !close(200)
    endif
    end subroutine phaselines_4

    module subroutine phaselines_4_move_q_point(phl4, i_in, pure_Q_flag)

    implicit none
    type(type_hdrt_4ph_lines) :: phl4
    integer :: pure_Q_flag,i,i_in

    i = i_in
    if (i > 2) then
        if ( (phl4%press_Q_mix(i-1) - phl4%press_Q_mix(i-2))/(phl4%press_Q_mix(pure_Q_flag) - phl4%press_Q_mix(i-1)) > 0.d0) then
            continue
        else
            i = 1
        endif
    else
        continue
    endif
    phl4%Q_point_mix(i) = phl4%Q_point_mix(pure_Q_flag)
    phl4%Qexist_mix(i) = phl4%Qexist_mix(pure_Q_flag)
    phl4%Temp_Q_mix(i) = phl4%Temp_Q_mix(pure_Q_flag)
    phl4%press_Q_mix(i) = phl4%press_Q_mix(pure_Q_flag)
    phl4%x_Q_ph1_mix(:,i) = phl4%x_Q_ph1_mix(:,pure_Q_flag)
    phl4%x_Q_ph2_mix(:,i) = phl4%x_Q_ph2_mix(:,pure_Q_flag)
    phl4%x_Q_ph3_mix(:,i) = phl4%x_Q_ph3_mix(:,pure_Q_flag)
    phl4%x_Q_ph4_mix(:,i) = phl4%x_Q_ph4_mix(:,pure_Q_flag)
    phl4%Q_point_mix(pure_Q_flag)     = ''
    phl4%Qexist_mix(pure_Q_flag)      = .false.
    phl4%Temp_Q_mix(pure_Q_flag)      = 0.d0
    phl4%press_Q_mix(pure_Q_flag)     = 0.d0
    phl4%x_Q_ph1_mix(:,pure_Q_flag)   = 0.d0
    phl4%x_Q_ph2_mix(:,pure_Q_flag)   = 0.d0
    phl4%x_Q_ph3_mix(:,pure_Q_flag)   = 0.d0
    phl4%x_Q_ph4_mix(:,pure_Q_flag)   = 0.d0
    end subroutine phaselines_4_move_q_point
    !******************************************************************************
    !******************************************************************************

    !This routine solves the massbalance for a fourphase equilibrium with different phasefractions = 0
    !inputs:
    !press                 -  pressure corresponding to current values of phasefractions
    !press_prev            -  pressure corresponding to previous values of phasefractions
    ! x                    -  overall composition
    ! x_fluidj             -  composition of the fluid phase j
    ! x_sol                -  composition of the solid phase
    ! x_hyd                -  composition of the hydrate phase
    ! vec_compphafra       -  vector with composition at the beginning and phasefractions after solving the sysofeqs
    ! vec_compphafra_prev  -  vector with composition at the beginning and phasefractions after solving the sysofeqs of the previous run of this sub
    ! eqtype_ind           -  depending on the equilibrium type different system of equations are solved -> see below in code
    ! root_dir             -  direction of the root: 1 if root is in positive direction of pressure step
    !                                               -1 if root is in negative direction of pressure step
    ! root_between         -  if the root is between the previous and current step the root is between this pressure range
    module subroutine massbalance_4ph(gl, press, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysicalx, root_dirx, root_between, root_sysofeqs)

    implicit none
    type(type_gl) :: gl
    double precision :: press, press_prev, press_sol(4,3,2)
    double precision, dimension(60,60) :: mat_matbal
    double precision, dimension(60,30) :: vec_compphafra, vec_compphafra_prev
    double precision, dimension(30):: x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd
    integer :: eqtype_ind, inphysicalx, inphysical, nsysofeqs, betas_good
    character(255):: herr
    !integer variables for loops etc
    integer :: ierr, neq, l, ll, i, ii, j
    !sign_changed is a vector that indicates if if the sign of beta_ph_i of sysofeqs l has changed from previous to current step
    integer :: root_between, root_dirx, root_dir, root_dir_changed, sign_changed(4), root_sysofeqs

    !create local copies
    root_dir = root_dirx
    inphysical = inphysicalx
    root_dir_changed = 0
    neq = 3

    if ((eqtype_ind == 1).or.(eqtype_ind == 2).or.(eqtype_ind == 3)) then  !VLxHIy, VLxLyH
        nsysofeqs = 2
    elseif (eqtype_ind == 4) then !LxLyHIz
        nsysofeqs = 3
    endif
    !create vector that contains the overall composition
    !each column contains the overall composition for a system of equations
    do l = 1,nsysofeqs
        vec_compphafra(:,l) = 0.d0
        vec_compphafra(1:30,l) = x(1:30)
    enddo

    if (eqtype_ind == 1) then !VLHI

        !possible check with this sysofeqs:
        !betaLw = 0 + betaV = 0
        !betaLw = 0 + betaH = 0
        !betaLw = 0 + betaIw = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_hyd
        mat_matbal(3,1:30) = x_sol
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,1),ierr,herr)
        vec_compphafra(4,1) = vec_compphafra(3,1) ! betaI
        vec_compphafra(3,1) = vec_compphafra(2,1) ! betaH
        vec_compphafra(2,1) = 0.d0                ! betaLw
        !vec_compphafra(:,2) = vec_compphafra(:,1)
        !vec_compphafra(:,3) = vec_compphafra(:,1)
        !new
        !betaH = 0 + betaV = 0
        !betaH = 0 + betaLw = 0
        !betaH = 0 + betaIw = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_fluid2
        mat_matbal(3,1:30) = x_sol
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,2),ierr,herr)
        vec_compphafra(4,2) = vec_compphafra(3,2)  !betaI
        vec_compphafra(3,2) = 0.d0 !betaH

    elseif (eqtype_ind == 2) then !VLwLcH
        !possible check with this sysofeqs:
        !betaH = 0 + betaV = 0
        !betaH = 0 + betaLw = 0
        !betaH = 0 + betaLx = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_fluid2
        mat_matbal(3,1:30) = x_fluid3
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,1),ierr,herr)

        !possible check with this sysofeqs:
        !betaLw = 0 + betaV = 0
        !betaLw = 0 + betaLx = 0
        !betaLw = 0 + betaH = 0
        mat_matbal(1,1:30) = x_fluid1!V
        mat_matbal(2,1:30) = x_fluid3!Lx
        mat_matbal(3,1:30) = x_hyd
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,2),ierr,herr)
        vec_compphafra(4,2) = vec_compphafra(3,2)  !betaH
        vec_compphafra(3,2) = vec_compphafra(2,2)  !betaLx
        vec_compphafra(2,2) = 0.d0                 !betaLw

    elseif (eqtype_ind == 3) then !VLxHIc

        !possible check with this sysofeqs:
        !betaLx = 0 + betaV = 0
        !betaLx = 0 + betaH = 0
        !betaLx = 0 + betaIx = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_hyd
        mat_matbal(3,1:30) = x_sol
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,1),ierr,herr)
        vec_compphafra(4,1) = vec_compphafra(3,1)  !Ic
        vec_compphafra(3,1) = vec_compphafra(2,1)  !H
        vec_compphafra(2,1) = 0.d0                 !Lx

        !betaIx = 0 + betaV = 0
        !betaIx = 0 + betaLx = 0
        !betaIx = 0 + betaH = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_fluid2
        mat_matbal(3,1:30) = x_hyd
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,2),ierr,herr)

    elseif (eqtype_ind == 4) then !LxLwHIy

        !possible check with this sysofeqs:
        !betaIy = 0 + betaLx = 0
        !betaIy = 0 + betaLw = 0
        !betaIy = 0 + betaH = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_fluid2
        mat_matbal(3,1:30) = x_hyd
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,1),ierr,herr)

        !possible check with this sysofeqs:
        !betaLw = 0 + betaLx = 0
        !betaLw = 0 + betaH = 0
        !betaLw = 0 + betaIy = 0
        mat_matbal(1,1:30) = x_fluid2
        mat_matbal(2,1:30) = x_hyd
        mat_matbal(3,1:30) = x_sol
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,2),ierr,herr)
        vec_compphafra(4,2) = vec_compphafra(3,2) !Iy
        vec_compphafra(3,2) = vec_compphafra(2,2) !H
        vec_compphafra(2,2) = vec_compphafra(1,2) !Lx
        vec_compphafra(1,2) = 0.d0                !Lw
        !possible check with this sysofeqs:
        !betaH = 0 + betaLw = 0
        !betaH = 0 + betaLx = 0
        !betaH = 0 + betaIy = 0
        mat_matbal(1,1:30) = x_fluid1
        mat_matbal(2,1:30) = x_fluid2
        mat_matbal(3,1:30) = x_sol
        call LUdecomp (gl,neq,mat_matbal,vec_compphafra(:,3),ierr,herr)
        vec_compphafra(4,3) = vec_compphafra(3,3) !Iy
        vec_compphafra(3,3) = 0.d0                !H

    endif


    !first step is very close to quadruple point of one guest
    !the values for phasefractions are saved and used for comparision in the next steps
    if (maxval(vec_compphafra_prev) < -8.5d10) then
        vec_compphafra_prev = vec_compphafra
    else

        do l = 1,4
            !sign_changed is a vector that indicates if if the sign of beta_ph_i of sysofeqs l has changed from previous to current step
            sign_changed = 0
            do i = 1,4
                !if sign changed the ratio is negative
                if ((vec_compphafra(i,l)/vec_compphafra_prev(i,l)) < 0) then
                    sign_changed(i) = 1
                endif
            enddo
            !if sign has changed for at least one phasefraction there migh be a root
            if ((sum(sign_changed) > 0).and.(sum(sign_changed) < 3)) then
                do i = 1,4
                    !if the sign did not change and the other phasefractions are within 0 < beta < 1 a root was found
                    if ((sign_changed(i) == 0).and.(((vec_compphafra(i,l) > 1.d-14).and. (vec_compphafra(i,l) < 1.d0)).or.((minval(vec_compphafra_prev(:,l)) > -1.d-14).and.(maxval(vec_compphafra_prev(:,l)) < 1.d0)))) then
                        root_between = root_between + 1
                    endif
                enddo
                !if no root was found perfom a linear interpolation between current and previous pressure and the corresponding phasefractions
                if (root_between == 0) then
                    if (press_prev > 0.d0) then
                        press_sol = 0.d0
                        sysofeqs: do ll = 1,4
                            j = 0
                            if (minval(vec_compphafra(:,ll)) > -8.9d10) then
                                phasefracloop: do ii = 1,4
                                    if (dabs(vec_compphafra(ii,ll)) > 1.d-14) then
                                        j = j + 1
                                        !in press_sol (#1,#2,#3) the pressures where beta_ph_ii = 0 or beta_ph_ii = 1 for the phase ii and the sysofeqs ll are stored:
                                        !#1 : corresponds to the system of equations, e.g. for VLwHIw only press_sol(1,:,:) and press_sol(2,:,:) is used
                                        !#2 : corresponds to the phases, e.g. for VLxLyH in press_sol(1,1,1) the pressure where beta_ph_V = 0 for sysofeqs 1 is stored
                                        !#3 : beta_ph_ii for sysofeqs = ll = 0 or 1      in press_sol(1,1,2) the pressure where beta_ph_V = 1 for sysofeqs 1 is stored
                                        !                                                in press_sol(1,2,1) the pressure where beta_ph_Lx = 0 for sysofeqs 1 is stored
                                        !                                                in press_sol(1,3,1) the pressure where beta_ph_Ly = 0 for sysofeqs 1 is stored
                                        !                                                in press_sol(1,4,1) the pressure where beta_ph_H = 0 for sysofeqs 1 is stored
                                        !                                                in press_sol(2,4,1) the pressure where beta_ph_H = 0 for sysofeqs 2 is stored
                                        press_sol(ll,j,1) = (- vec_compphafra(ii,ll) + (vec_compphafra(ii,ll) - vec_compphafra_prev(ii,ll)) / (press-press_prev) * press ) /((vec_compphafra(ii,ll) - vec_compphafra_prev(ii,ll)) / (press-press_prev))
                                        press_sol(ll,j,2) = (1 - vec_compphafra(ii,ll) + (vec_compphafra(ii,ll) - vec_compphafra_prev(ii,ll)) / (press-press_prev) * press ) /((vec_compphafra(ii,ll) - vec_compphafra_prev(ii,ll)) / (press-press_prev))
                                        !don't look at solutions before previous pressure value because it was already found in the previous calculations that there is no root
                                        !if ((press-press_prev) < 0) then !press_step < 0
                                        !    if (press_sol(ll,j,1) > press_prev) press_sol(ll,j,1) = press_prev
                                        !    if (press_sol(ll,j,2) > press_prev) press_sol(ll,j,2) = press_prev
                                        !elseif ((press-press_prev) > 0) then !press_step > 0
                                        !    if (press_sol(ll,j,1) < press_prev) press_sol(ll,j,1) = press_prev
                                        !    if (press_sol(ll,j,2) < press_prev) press_sol(ll,j,2) = press_prev
                                        !endif
                                        if ((press-press_prev) > 0) then !press_step > 0
                                            if (((maxval(press_sol(ll,j,:)) > press_prev).and.(maxval(press_sol(ll,j,:)) < press)) .or. ((minval(press_sol(ll,j,:)) > press_prev).and.(minval(press_sol(ll,j,:)) < press))) then
                                                root_between = root_between + 1
                                            else
                                                root_between = 0
                                                exit phasefracloop
                                            endif
                                        elseif ((press-press_prev) < 0) then !press_step < 0
                                            if (((maxval(press_sol(ll,j,:)) < press_prev).and.(maxval(press_sol(ll,j,:)) > press)) .or. ((minval(press_sol(ll,j,:)) < press_prev).and.(minval(press_sol(ll,j,:)) > press))) then
                                                root_between = root_between + 1
                                            else
                                                root_between = 0
                                                exit phasefracloop
                                            endif
                                        endif
                                    endif
                                enddo phasefracloop
                                !a root is found if all pressure values are within previous and current pressure steps
                                !if ((press-press_prev) < 0) then !press_step < 0
                                !    if (minval(press_sol(ll,:,:)) >= press) root_between = root_between + 2
                                !elseif ((press-press_prev) > 0) then !press_step > 0
                                !    if (maxval(press_sol(ll,:,:)) <= press) root_between = root_between + 2
                                !endif
                            endif
                        enddo sysofeqs
                    endif
                endif
                !if the sign has changed for all phasefractions there must be a root, no further checks required
            elseif (sum(sign_changed) >= 3) then
                root_between = root_between + 1
            endif
            !root was found but root_directions has to be changed
            if ((root_between > 0) .and.(root_dir_changed == 0)) then
                if (root_dir == 0) then
                    root_dir = -1 !root is in negative direction
                else
                    root_dir = -root_dir
                endif
                root_dir_changed = 1
                root_sysofeqs = l
            endif
        enddo
        root_dirx = root_dir
        vec_compphafra_prev = vec_compphafra
    endif
    end subroutine massbalance_4ph


    !********************************************************************************************
    module subroutine GenerateFitFile(gl, pel, Filepath,Datenfilepath,typeofcalc,dT_ab_Q,dp_ab_Q, &
        &          Q_point, Qexist,Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    !********************************************************************************************




    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    !Path for the Daten and CPO files + UNC file with results of the work_form fitting
    character(len=255) :: Datenfilepath, CPOfilepath, UNCfilepath, UNCfilepath_Dp, CompsHdrtsAll
    character(len=6) :: dummy                           !used to jump over the not needed lines of the inputfile
    character(len=255) :: author
    character(len=6) :: EqtypeIn                        !Transform the eqtype to a general form: Le,Lp,Lc -> Lg (g = guest)

    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3
    double precision, dimension(30) :: x_exp_ph1, x_exp_ph2, x_exp_ph3      ! experimental compositions
    double precision, dimension(30) :: Chempot_ph1, lnfi_ph1, Chempot_ph2, lnfi_ph2
    double precision, dimension(30) :: Chempot_ph3, lnfi_ph3

    double precision:: Temp, press
    double precision:: rho_ph1, rho_ph2, rho_ph3

    integer :: iTp, i, j, error, i3p, in_scr
    double precision :: weight

    !Variables needed for PhaseDet
    double precision, dimension(5) :: rho_out

    !Variables for the uncertainty estimates
    integer :: typeofcalc
    character(len=255) :: Filepath
    double precision, dimension(2) :: dT_ab, dp_ab, dT_ab_Q, dp_ab_Q
    character(len=255) :: Describ_path = ''   !Path for Hydrate model description file
    character(12), dimension(30) :: COMPONENTSX

    ! Variables for the regula falsi - optimization of dT or dp
    double precision, dimension(2) :: Work_Q, Work_Q_2 , Work_Q_3 ! [J] work of formation at the Q-point for T and p shifted
    double precision, External :: WorkForm_diff         ! external function for difference from fitted work of formation
    double Precision, dimension(65) :: parameters       ! needed to transmit T, p, etc. to Regula Falsi
    double Precision, dimension(2) :: Tp, Work_Tp, Work_Tp_diff, Tp_min, Tp_max, Delta_allowed, &
        Tp_min_allowed, Tp_max_allowed, Chempot_Tp, p_Lapl_Tp     !, lnphi_Tp
    integer :: Iterations, Max_Iterations
    integer :: Errorflag
    double precision ::  u_cpo, Chempot_exp
    !integer :: ph_eq_type
    integer, dimension(2) :: error_RF, iter_RF, ph_eq_type

    ! Work of formation at the shifted T or p
    double precision:: work_form, p_Lapl

    !Variables at the Quadruple point
    character(len=6), dimension(4), intent(in)  :: Q_point
    double precision, dimension (6), intent(in) :: press_Q, Temp_Q
    double precision, dimension(4,30), intent(in) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist

    !Variables at the three-phase equilibrium
    !double precision, dimension(30) :: x_ph1, x_ph2, x_ph3
    !double precision :: rho_ph1, rho_ph2, rho_ph3
    !double precision, dimension(30) :: Chempot, lnfi
    double precision, dimension(2) :: Temp_3eq, press_3eq, Chempot_3eq   ! , lnfi_3eq, occup
    integer, dimension(2) :: error_3p

    ! Variables for work check
    integer :: error_check, dT_VHIw_dp_VLwH, Dp_write
    double Precision, dimension(3) :: Tp3, work3

    ! Variables for mixed hydrates
    double precision, dimension(2,30) :: lnphi_Tp
    double precision, dimension(30) :: fug_gas
    double precision, dimension(2,30) :: lnfi_3eq
    double precision, dimension(3,30) :: occup

    !//////////////////////////////////////////////////////////////

    COMPONENTSX = gl%Fluidlist_hydrate     ! defined in module solids

    error = 0
    x = 0.D0    ! overall composition
    x_ph1 = 0.D0
    x_ph2 = 0.D0
    x_ph3 = 0.D0

7007 format(a8,e15.8,a9,e15.8)

    pel%CompsHdrtsAll = trim(COMPONENTSX(2))
    if (gl%nrofhydrateformers > 2) then
        do i = 3,gl%nrofhydrateformers
            pel%CompsHdrtsAll = trim(pel%CompsHdrtsAll) // '_' // trim(COMPONENTSX(i))
        end do
    end if

    Select Case (typeofcalc)
    Case(11)
        !The z_ChemPot_out file type is .cpo
        CPOfilepath = trim(Filepath) // trim(pel%CompsHdrtsAll) // '.cpo'
        dT_ab = (/0.d0, 0.d0/)
        dp_ab = (/0.d0, 0.d0/)
    Case(12)
        ! Calculation at decreased T
        CPOfilepath = trim(Filepath) // trim(pel%CompsHdrtsAll) // '_T.cpo'
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define slope (dT_a) for T-variation. [0]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dT_ab(1)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define constant term (dT_b) for T-variation in K. [0.5]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dT_ab(2)
        dp_ab = (/0.d0, 0.d0/)
    Case(13)
        ! Calculation at increased p
        CPOfilepath = trim(Filepath) // trim(pel%CompsHdrtsAll) // '_p.cpo'
        dT_ab = (/0.d0, 0.d0/)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define slope (dp_a) for p-variation. [0]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dp_ab(1)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define constant term (dp_b) for p-variation in MPa. [0.5]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dp_ab(2)
    Case(14)
        ! Optimization of dT_ab and dp_ab from work of formation
        CPOfilepath = trim(Filepath) // trim(pel%CompsHdrtsAll) // '_workF.cpo'
        UNCfilepath = trim(Filepath) // trim(pel%CompsHdrtsAll) // '_dT_dp.unc'
        ! write Laplace pressure in 'gas_dT_dp_Dp.unc'
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Write Laplace pressure in file gas_dT_dp_Dp.unc? [1-NO]:'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) Dp_write
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        UNCfilepath_Dp = trim(Filepath) // trim(pel%CompsHdrtsAll) // '_dT_dp_Dp.unc'
        dT_ab = (/0.d0, 0.d0/)
        dp_ab = (/0.d0, 0.d0/)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

        write(*,*) 'Definition of the work of formation:'
        write(*,*) ' 1 - work form evaluated at the Q-point for given dT and dp,'
        write(*,*) ' 2 - work from at given three-phase equilibrium point,'
        write(*,*) ' 3 - work from defined in the code.'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) in_scr
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE

        if (in_scr == 1) then
            ! Work of formation at the lower Q-point VLwHIw
            Temp = gl%T_Q_Hyd_2C(1)
            press = gl%p_Q_Hyd_2C(1)
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

            write(*,*) 'Definition of 2-phase equilibrium used at Q-point for work_form:'
            write(*,*) ' 1 - dT and dp steps calculated from stable VIw or VLw,'
            write(*,*) ' 2 - dT-step and dp-step calculated ONLY from VLw,'
            write(*,*) ' 3 - dT-step and dp-step calculated ONLY from VIw.'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) dT_VHIw_dp_VLwH
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) ' '
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE

        elseif (in_scr == 2) then
            ! Work of formation at the three-equilibrium point
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) 'Please define equilibrium type [VHIw, VLwH, LwLxH]:'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) EqtypeIn
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) ' - Please define phase equilibrium temperature [K]:'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) Temp
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) ' - Please provide ESTIMATE for pressure [MPa]:'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) press

            iTp = 1
            call threephase_equil(gl,pel,EqtypeIn,1, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, &
                &   rho_ph1, rho_ph2, rho_ph3, Chempot_ph1, lnfi_ph1, occup, error_3p(iTp))!, &
            !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

            write(*,*) ' - Calculated equilibrium pressure [MPa] p = ', press
            write(*,*) ' '
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE

            dT_VHIw_dp_VLwH = 1     ! stable 2-phase equil. VLw or VIw will be found

        elseif (in_scr == 3) then
            ! Work of formation defined HERE
            !*******************************************************
            Work_Q(1) = 4.5d-17     ! work form for temperature step
            Work_Q(1) = 4.5d-17     ! work form for pressure step
            Work_Q(2) = Work_Q(1)
            !*******************************************************

8007        format(a12,e15.8,a14,e15.8)
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,8007) 'work_f(dT)= ', Work_Q(1), ', work_f(dp)= ', Work_Q(2)
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        end if


        if (in_scr /= 3) then
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

            write(*,*) 'Please define TEMPERATURE scatter in K for work_form. [0 or 10]'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) dT_ab(2)
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) 'Please define PRESSURE scatter in MPa for work_form. [0.02]'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) dp_ab(2)

            ! Work of formation at the lower Q-point VLwHIw
            !call WorkForm_Qpoint(dT_VHIw_dp_VLwH,Temp,press,dT_ab,dp_ab,Work_Q)

            if (in_scr == 2) then
                Work_Q(2) = Work_Q(1)   ! same work_form for T-step and p-step

            elseif ((in_scr == 1).and.(dT_VHIw_dp_VLwH == 1)) then

                ! work_form as geometrical average for VLw and VIw
                !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
                write(*,*) ' '
                write(*,*) ' Shall work_form be taken as geom. average of w(VLw) and w(VIw)? [1-Yes]'
                !DEC$ END IF ! WO_WRITE_TO_CONSOLE
                read(*,*) i

                if (i == 1) then
                    ! Work of formation at the lower Q-point VLwHIw
                    !call WorkForm_Qpoint(2,Temp,press,dT_ab,dp_ab,Work_Q_2)   ! dT-step and dp-step calculated ONLY from VLw
                    !call WorkForm_Qpoint(3,Temp,press,dT_ab,dp_ab,Work_Q_3)   ! dT-step and dp-step calculated ONLY from VIw

                    Work_Q(1) = dsqrt(Work_Q_2(1)*Work_Q_3(1))
                    Work_Q(2) = dsqrt(Work_Q_2(2)*Work_Q_3(2))
                end if


            end if
        end if

8008    format(a20,f10.3,a11,f10.4)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        write(*,8008) 'Given point: T[K] = ',Temp,', p[MPa] = ',press
        write(*,*) ' '
        write(*,*) ' Work of form. at dT [J] = ', Work_Q(1)
        write(*,*) ' Work of form. at dp [J] = ', Work_Q(2)
        write(*,*) ' '
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        !Work_Q(1) = Work_Q(1)/2.d0 + Work_Q(2)/2.d0
        !Work_Q(2) = Work_Q(1)
        !write(*,*) ' Average work form.  [J] = ', Work_Q(2)
        !write(*,*) ' '

    End select

    dT_ab_Q = dT_ab
    dp_ab_Q = dp_ab

    ! Quadruple points stored in parameters for Regula Falsi
    parameters(20) = Temp_Q(1)
    parameters(21) = Temp_Q(2)
    parameters(22) = Temp_Q(3)
    parameters(23) = Temp_Q(4)

    parameters(24) = press_Q(1)
    parameters(25) = press_Q(2)
    parameters(26) = press_Q(3)
    parameters(27) = press_Q(4)

3011 format (a6)
3012 format (a255)
3013 format (a)
    !CPO output:Equil, weight, Temp, press, chempot_H2O, Chempot_gas, lnf_H2O, lnf_gas, ...
    !           x_H2O(phase1), x_gas(phase1), x_H2O(phase2), x_gas(phase2), x_H2O(phase3), x_gas(phase3), dens(phase1), dens(phase2)
3020 format (a6,' ',f8.4,' ',f8.4,' ',f9.4,' ',f14.6,' ',f14.6,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10 ,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2)
    !CPO output:Equil, weight, Temp, press, error, error, ...
3019 format (a6,' ',f8.4,' ',f8.4,' ',f9.4,' ',i10,' ',i10,' ',i10,' ',i10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10 ,' ',f12.10,' ',f12.10,' ',i10,' ',i10)
    !CPO at shifted T or p - last column = work of formation for given dT or dp    ... work_form at last position
3021 format (a6,' ',f8.4,' ',f8.4,' ',f9.4,' ',f14.6,' ',f14.6,' ',f12.8,' ',f12.8,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10 ,' ',f12.10,' ',f12.10,' ',f8.2,' ',f8.2,' ',e16.8)
    !CPO from fitting of the work of formation
    !            Equil, weight, Temp, press, T_3eq, chempot_3eq(dT), lnfi_3eq(dT), T(dT), Work_T(dT), chempot_H2O(dT), lnf_gas(dT),...
    !            p_3eq, chempot_3eq(dp), lnfi_3eq(dp),  p(dp), Work_p, chempot_H2O(dp), lnf_gas(dp)
3022 format (a6,' ',f8.4,' ',f8.4,' ',f9.4,' ',f9.4,' ',f14.6,' ',f12.8 &
        &,' ',f10.5,' ',f14.6,' ',f12.8)
    ! UNC Output:Equil, weight, Temp, press, chempot_exp, T(dT), Work_T, chempot_H2O(dT), iter, error_T, p(dp), Work_p, chempot_H2O(dT), iter, error_p, u_cpo
3023 format (a5,' ',f8.4,' ',f8.4,' ',f10.5,' ',f14.6,' ',f9.4,' ',e16.8,' ',f14.6,' ',i4,' ',i6,' ',f10.5,' ',e16.8' ',f14.6,' ',i4,' ',i6,' ',f14.6)
    ! UNC with Laplace pressure Dp
    ! UNC Output:Equil, weight, Temp, press, chempot_exp, T(dT), Work_T, chempot_H2O(dT), iter, error_T, Dp(T), p(dp), Work_p, chempot_H2O(dT), iter, error_p, Dp(p), u_cpo
3024 format (a5,' ',f8.4,' ',f8.4,' ',f10.5,' ',f14.6,' ',f9.4,' ',e16.8,' ',f14.6,' ',i4,' ',i6,' ',f12.7,' ',f10.5,' ',e16.8' ',f14.6,' ',i4,' ',i6,' ',f12.7,' ',f14.6)


    !Open daten file and output file
    !----------------------------------------------------------------------------------
    !!The path, were the Hydrate data file is stored
    open(unit=10, file=trim(Datenfilepath), status='old', action='read', iostat=error)
    if (error /= 0) then
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Data file not found, calculation aborted'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        return
    end if
    !The path, were the calculated chemical potentials will be written to
    open(unit=11, file=trim(CPOfilepath), status='unknown', action='write', iostat=error)
    if (error /= 0) then
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'An error occured when trying to generate the CPO file'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        return
    end if
    !The path, were the results of the work of formation fitting will be stored
    if (typeofcalc == 4) then
        open(unit=12, file=trim(UNCfilepath), status='unknown', action='write', iostat=error)
        if (Dp_write /= 1) then
            open(unit=22, file=trim(UNCfilepath_Dp), status='unknown', action='write', iostat=error)
        end if
        if (error /= 0) then
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) 'An error occured when trying to generate the UNC file'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            return
        else
7008        format (a8, f10.5, a10, f8.5, a11, e15.8)
            write(12,7008)  '! Temp= ', Temp,  ' dT_ab(2)= ', dT_ab(2), '  workQ(dT)= ',Work_Q(1)
            write(12,7008)  '! pres= ', press, ' dp_ab(2)= ', dp_ab(2), '  workQ(dp)= ',Work_Q(2)
7068        format (a,a)
            write(12,7068) '! Eq., weight, T_exp, p_exp, cpo_exp, T(wfQ), wf(dT), cpo(dT), iter, error_dT, ', &
                & 'p(wfQ), wf(dp), cpo(dp), iter, error_p, u(cpo)'
            write(11,7068) '! Eq., weight, T_exp, p_exp, T_3eq, cpo(Teq), lnfi(Teq)', &
                &  ', p_3eq, cpo(peq), lnfi(peq)'
            if (Dp_write /= 1) then
                write(22,7008)  '! Temp= ', Temp,  ' dT_ab(2)= ', dT_ab(2), '  workQ(dT)= ',Work_Q(1)
                write(22,7008)  '! pres= ', press, ' dp_ab(2)= ', dp_ab(2), '  workQ(dp)= ',Work_Q(2)
                write(22,7068) '! Eq., weight, T_exp, p_exp, cpo_exp, T(wfQ), wf(dT), cpo(dT), iter, error_dT, Dp(T), ', &
                    & 'p(wfQ), wf(dp), cpo(dp), iter, error_p, Dp(p), u(cpo)'
            end if
        end if
    end if
    !----------------------------------------------------------------------------------



    if (error == 0) then
        Do while (.true.)
            !Go through file, until first row, which has no "!" as first sign
            !Read in Authors name
            !Read in first 6 signs to determine which phase equilibrium data is given
            !Go back and read full line
            !Try to calculate the two phase equilibria without Hydrates at every three phase data point
            !Write results in output textfile
            Read(10,3011) dummy
            dummy = trim(dummy)
            !End of file reached
            if (trim(dummy) == '!End') exit
            !Begin reading the data
            if (dummy(1:1) /= '!') then
                backspace(10)
                Read(10,3012) author
                !Read all data corresponding to that author
                write(11,3013) '! ' // trim(author)
                if (typeofcalc == 4) then
                    write(12,3013) '! ' // trim(author)
                    if (Dp_write /= 1) then
                        write(22,3013) '! ' // trim(author)
                    end if
                end if
                Do while (.true.)
                    Read(10,*) dummy
                    dummy = trim(dummy)
                    !Transform the eqtype
                    EqtypeIn = dummy
                    If ((EqtypeIn == 'LwLcH') .or. (EqtypeIn == 'LwLpH').or. (EqtypeIn == 'LwLeH')) then
                        EqtypeIn = 'LwLxH'
                    end if
                    If ((EqtypeIn == 'VLcH') .or. (EqtypeIn == 'VLpH') .or. (EqtypeIn == 'VLeH')) then
                        EqtypeIn = 'VLxH'
                    end if
                    backspace(10)

                    !----------------------------------------------------------
                    ! Calculation of the 2-phase equilibria along 3-phase lines
                    if ((EqtypeIn == 'VLwH').or.(EqtypeIn == 'VHIw').or.(EqtypeIn == 'VLxH').or.(EqtypeIn == 'LwLxH')) then

                        ! Type of 3-phase equilibria as integer 'ph_eq_type'
                        if (EqtypeIn == 'VLwH') then
                            ph_eq_type(1) = 1
                        else if (EqtypeIn == 'VHIw') then
                            ph_eq_type(1) = 2
                        else if (EqtypeIn == 'VLxH') then
                            ph_eq_type(1) = 3
                        else if (EqtypeIn == 'LwLxH') then
                            ph_eq_type(1) = 4
                        end if


                        Read(10,*) dummy, weight, Temp, press, x(1), x(2), x_exp_ph1(1), x_exp_ph1(2), &
                            &   x_exp_ph2(1), x_exp_ph2(2), x_exp_ph3(1), x_exp_ph3(2)
                        if (x_exp_ph1(1) /= 0.d0) then
                            x_ph1 = x_exp_ph1
                            x_ph2 = x_exp_ph2
                            x_ph3 = x_exp_ph3
                        end if

                        ! Added for CPO uncertainty evaluation
                        if (typeofcalc < 14 ) then
                            ! Exp. point shifted for CPO
                            Temp = Temp - (dabs(gl%T_Q_Hyd_2C(1)-Temp)*dT_ab(1) + dT_ab(2)) ! exp. T decreased
                            press = press + (dabs(gl%p_Q_Hyd_2C(1)-press)*dp_ab(1) + dp_ab(2)) ! exp. p increased
                            x_ph1 = x_exp_ph1
                            x_ph2 = x_exp_ph2
                            x_ph3 = x_exp_ph3

                            call cpo_2phase_at_3phase (ph_eq_type(1),Temp,press,x,Chempot_ph1, lnfi_ph1, &
                                &           x_ph1, x_ph2, x_ph3, rho_out, error)


                            ! Data upload in a cpo file
                            if(error == 0) then
                                ! Calculation of the work of formation
                                if ((typeofcalc == 12).or.(typeofcalc==13)) then
                                    ! Mixed hydrates 2015
                                    fug_gas(1:gl%nrofhydrateformers-1) = 1.d6*dexp(lnfi_ph1(2:gl%nrofhydrateformers))
                                    !call WorkForm(Temp,press*1.d6,Chempot_ph1(1),fug_gas,work_form,p_Lapl)
                                    write(11,3021) dummy, weight, Temp, press, Chempot_ph1(1), Chempot_ph1(2), &
                                        & lnfi_ph1(1), lnfi_ph1(2), x_ph1(1), x_ph1(2), &
                                        & x_ph2(1), x_ph2(2), x_ph3(1), x_ph3(2), rho_out(1), rho_out(2), work_form
                                else
                                    write(11,3020) dummy, weight, Temp, press, Chempot_ph1(1), Chempot_ph1(2), &
                                        & lnfi_ph1(1), lnfi_ph1(2), x_ph1(1), x_ph1(2), &
                                        & x_ph2(1), x_ph2(2), x_ph3(1), x_ph3(2), rho_out(1), rho_out(2)
                                end if
                            else
                                !write(11,*) dummy,' ',error
                                write(11,3019) dummy, weight, Temp, press, error, error, &
                                    & error, error, 0.d0, 0.d0, &
                                    & 0.d0, 0.d0, 0.d0, 0.d0, error, error
                                error = 0
                            end if

                            !----------------------------------------------------------
                            ! Solution of dT and dp from the work of formation approach
                        elseif (typeofcalc == 14 ) then

                            if (EqtypeIn == 'VLxH') then
                                weight = 0.d0
                            end if

                            if (weight /= 0.d0) then

                                ! Three-phase equilibrium at press (1) or temp (2)
                                do iTp = 1,2  ! temperature or pressure fitting
                                    Temp_3eq(iTp) = Temp
                                    press_3eq(iTp) = press

                                    ! Input to 'threephase_equil' - type of flash calculation
                                    if (iTp == 1) then      ! T varied, p constant
                                        i3p = 2
                                    elseif (iTp == 2) then  ! p varied, T constant
                                        i3p = 1
                                    end if

                                    x = 0.d0; x_ph1 = 0.d0; x_ph2 = 0.d0; x_ph3 = 0.d0    ! new equil type - initial composition will be taken from Q-points in 'threephase_equil'

                                    call threephase_equil(gl,pel,EqtypeIn,1,i3p,Temp_3eq(iTp), press_3eq(iTp),x , x_ph1, x_ph2, x_ph3, &
                                        &   rho_ph1, rho_ph2, rho_ph3, Chempot_ph1, lnfi_ph1, occup, error_3p(iTp))!, &
                                    !&   Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)

                                    Chempot_3eq(iTp) = Chempot_ph1(1)
                                    !lnfi_3eq(iTp) = lnfi_ph1(2)
                                    ! Mixed hydrates 2015
                                    lnfi_3eq(iTp,1:gl%nrofhydrateformers-1) = lnfi_ph1(2:gl%nrofhydrateformers)
                                end do

                                !*****************************
                                ! Regula_Falsi input variables
                                Tp = (/Temp, press/)
                                if (dT_ab_Q(2) > 0.d0) then
                                    Tp_min = (/Temp_3eq(1)-8.d0*dT_ab_Q(2), press_3eq(2)/)
                                    Tp_min_allowed = (/Temp_3eq(1)-10.d0*dT_ab_Q(2), press_3eq(2)/)
                                else
                                    Tp_min = (/Temp_3eq(1)-20.d0, press_3eq(2)/)
                                    Tp_min_allowed = (/Temp_3eq(1)-30.d0, press_3eq(2)/)
                                end if

                                !if (press_3eq(2) > 50.d0) then
                                !    Tp_max = (/Temp_3eq(1), press_3eq(2)+50.d0 /)
                                !    Tp_max_allowed = (/Temp_3eq(1), press_3eq(2)+130.d0/)
                                !else
                                !    Tp_max = (/Temp_3eq(1), press_3eq(2)*1.3d0/)
                                !    Tp_max_allowed = (/Temp_3eq(1), press_3eq(2)*3.d0/)
                                !end if

                                if (press_3eq(2) > 50.d0) then
                                    Tp_max = (/Temp_3eq(1), press_3eq(2)*5.d0 /)
                                    Tp_max_allowed = (/Temp_3eq(1), press_3eq(2)+150.d0/)
                                else
                                    Tp_max = (/Temp_3eq(1), press_3eq(2)*5.0d0/)
                                    Tp_max_allowed = (/Temp_3eq(1), press_3eq(2)*10.d0/)
                                end if

                                Delta_allowed = (/0.0001d0, press/100.d0*0.001d0/)
                                Delta_allowed = (/1.d-4, 1.d-7/)
                                Max_Iterations = 100

                                parameters(3) = ph_eq_type(1)   ! phase equilibrium type 1-VLwH, 2-VHIw, 3-VLxH, 4-LwLxH
                                parameters(4) = Temp            ! [K]
                                parameters(5) = press           ! [MPa]

                                do iTp = 1,2  ! temperature or pressure fitting

                                    parameters(1) = iTp             ! 1-Temp [K], 2-press [MPa] calculation
                                    parameters(2) = Work_Q(iTp)     ! [J] work of formation in the Q-point at T-lowered (1), p-increased (2)

                                    ! Initial estimate for root Tp(iTp) must lie within the limits
                                    Tp(iTp) = Tp_min(iTp)/2.d0 + Tp_max(iTp)/2.d0

                                    if (error_3p(iTp) == 0) then
                                        !call Regula_Falsi(WorkForm_diff,Tp(iTp),Tp_min(iTp),Tp_max(gl,iTp),Delta_allowed(iTp), &
                                        !     &  Tp_min_allowed(iTp),Tp_max_allowed(iTp),Max_Iterations,Iterations,errorflag,parameters)

                                        !Work_Tp_diff(iTp) = WorkForm_diff    ! [J] difference of the work of formation from hydrate and fluid models
                                        Work_Tp(iTp) = parameters(6)          ! [J] work of formation
                                        Chempot_Tp(iTp) = parameters(7)       ! [J/mol] chemical potential of water
                                        p_Lapl_Tp(iTp)  = parameters(9)                ! [Pa] Laplace pressure

                                        !lnphi_Tp(iTp) = dlog(parameters(8)*1.d-6)      ! [MPa] log(fugacity of gas)
                                        ! Mixed hydrates
                                        lnphi_Tp(iTp,1:gl%nrofhydrateformers-1) = dlog(parameters(10:9+gl%nrofhydrateformers-1)*1.d-6)  ! [MPa] log(fugacity of gas)

                                        error_RF(iTp) = errorflag
                                        iter_RF(iTp) = Iterations

                                        parameters(63+iTp) = parameters(64)   ! error from the 'cpo_2phase_at_3phase' subroutine
                                        ! parameters(64) ... error from the 'cpo_2phase_at_3phase' subroutine - T optimization
                                        ! parameters(65) ... error from the 'cpo_2phase_at_3phase' subroutine - p optimization

                                    else
                                        error_RF(iTp) = error_3p(iTp)  ! bad 3-phase equil, i.e. bad input conditions for the Regula_Falsi
                                        Tp(iTp) = 0.d0
                                        Work_Tp(iTp) = 0.d0
                                        iter_RF(iTp) = 0.d0
                                        Chempot_Tp(iTp) = 0.d0
                                        lnphi_Tp(iTp,:) = 0.d0
                                    end if

                                end do

                                !dT_ab(2) = Temp_3eq(1) - Tp(1)         ! [K] temperature step for given work of formation
                                !dp_ab(2) = press_3eq(2) - Tp(2)        ! [MPa] pressure step for given work of formation
                                !
                                !dT_ab(2) = Temp - Tp(1)         ! [K] temperature step for given work of formation
                                !dp_ab(2) = press - Tp(2)        ! [MPa] pressure step for given work of formation

                                ! CPO combined uncertainty
                                if ((error_RF(1) == 0).and.(error_RF(2) == 0)) then

                                    ! Calculation of CPO at the experimental point Temp, press for converged points
                                    call cpo_2phase_at_3phase (ph_eq_type(1),Temp,press,x,Chempot_ph1, lnfi_ph1, &
                                        &           x_ph1, x_ph2, x_ph3, rho_out, error)
                                    Chempot_exp = Chempot_ph1(1)
                                    !u_cpo = dsqrt( (Chempot_Tp(1)-Chempot_exp)**2 + (Chempot_Tp(2)-Chempot_exp)**2  )


                                    ! Chemical potential at three-phase equilibrium calculated from the model
                                    Chempot_ph1(1) = Chempot_3eq(1)  ! chemical potential at equilibrium point press = p_experimental, T - varied
                                    Chempot_ph1(2) = Chempot_3eq(2)  ! chemical potential at equilibrium point Temp = T_experimental, p - varied

                                    u_cpo = dsqrt( (Chempot_Tp(1)-Chempot_ph1(1))**2 + (Chempot_Tp(2)-Chempot_ph1(2))**2  )

                                else    ! if error occurs - data point is not taken into account
                                    u_cpo = 0.d0
                                    weight = 0.d0
                                end if

                                ! Data upload in a CPO_workF file
                                write(11,3022) dummy, weight, Temp, press, &
                                    & Temp_3eq(1),  Chempot_3eq(1), lnfi_3eq(1,1), &
                                    & press_3eq(2), Chempot_3eq(2), lnfi_3eq(2,1)
                                ! Data upload in a UNC file
                                write(12,3023) dummy, weight, Temp, press, Chempot_exp, Tp(1), Work_Tp(1), Chempot_Tp(1), iter_RF(1), error_RF(1), &
                                    &  Tp(2), Work_Tp(2), Chempot_Tp(2), iter_RF(2), error_RF(2), u_cpo
                                ! Data upload in a UNC file with Laplace pressure
                                if (Dp_write /= 1) then
                                    write(22,3024) dummy, weight, Temp, press, Chempot_exp, Tp(1), Work_Tp(1), Chempot_Tp(1), iter_RF(1), error_RF(1), p_Lapl_Tp(1)/1.d6, &
                                        &  Tp(2), Work_Tp(2), Chempot_Tp(2), iter_RF(2), error_RF(2), p_Lapl_Tp(2)/1.d6, u_cpo
                                end if

                            else    ! weight = 0

                                ! Data upload in a CPO file
                                write(11,3022) dummy, weight, Temp, press, &
                                    & 0.d0,  0.d0, 0.d0, &
                                    & 0.d0, 0.d0, 0.d0
                                ! Data upload in a UNC file
                                write(12,3023) dummy, weight, Temp, press, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                    &  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
                                ! Data upload in a UNC file with Laplace pressure
                                write(22,3024) dummy, weight, Temp, press, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                    &  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0

                            end if  ! end weight /= 0

                            ph_eq_type(2) = ph_eq_type(1)

                        end if ! end of EqtypeIn
                        !----------------------------------------------------------

                    elseif (EqtypeIn == '!Next') then
                        exit
                    else
                        Read(10,3011) dummy
                    end if

                end do
            end if

        end do
    end if
    close(10)   ! HDT Daten file
    write(11,3013) 'End'    ! CPO file
    close(11)
    if (typeofcalc == 4) then
        write(12,3013) 'End'    ! UNC file
        close(12)
        if (Dp_write /= 1) then
            write(22,3013) 'End'    ! UNC file
            close(22)
        end if
    end if

    end subroutine GenerateFitFile
    !******************************************************************************
    !******************************************************************************






    !******************************************************************************
    module subroutine cpo_2phase_at_3phase (ph_eq_type, Temp, press, x, &
        &           Chempot_ph1, &
        &           lnfi_ph1, x_ph1, x_ph2, &
        &           x_ph3, rho_out, error)
    !******************************************************************************
    ! Calculation of a single two-phase equilibrium point at the three-phase equilibrium line




    implicit none
    type(type_gl), allocatable :: gl
    integer :: ph_eq_type

    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3 ! composition 1-vap, 2-liq, 3-sol
    double precision, dimension(30) :: Chempot_ph1, Chempot_ph2, Chempot_ph3
    double precision, dimension(30) :: lnfi_ph1, lnfi_ph2, lnfi_ph3

    double precision:: Temp, press, Dens
    double precision:: beta, rho_ph1, rho_ph2, rho_ph3
    integer :: iPhase, iFlash, iter, i, j, error, iPhase_sol
    double precision:: weight

    !Variables needed for PhaseDet
    double precision, dimension(5):: rho, rho_out
    double precision, dimension(3):: rho2P
    integer:: nrofphases
    double precision, dimension(30,5):: x_Phase
    integer, dimension(5) :: phasetype                  !phasetype contains the phase indicator number

    !The iPhase_sol can be used to assume fluid phases in two and three phase equilibrium calculation
    ! 0 : Let the algorithms choose the phases
    ! 1 : Liquid phase of the hydrate guest is assumed (e.g. LwH or LwLcH)
    ! 2 : Vapor phase of the hydrate guest is assumed (e.g. VH or VLwH
    !NOTE: IN case of two fluid phases, the heavier phase is always assumed as liquid phase (no vapor / vapor equilibrium)
    iPhase_sol = 0

    error = 0
    beta = 0.D0
    rho_ph1 = 0.D0
    rho_ph2 = 0.D0
    rho_ph3 = 0.D0

    x_ph1 = 0.D0
    x_ph2 = 0.D0
    x_ph3 = 0.D0

    Select Case (ph_eq_type)

    case (1)    !('VLwH')
        !ph1 = vap, ph2 = liq2, ph3 = hyd
        !Three phase equilibrium Vapor (V), Liquid water with dissolved component(s) (Lw) and Hydrate (H)

        !       A guest component rich vapor phase and a water rich liquid phase are expected
        !       If no overall composition is given: Set the overall composition to a value that is most likely inbetween the two phases in equilibrium:
        if ((x(1) == 0.D0) .OR. (x(2) == 0.D0)) then
            x = 0.D0
            ! Mixed hydrates 2015
            do i = 1,gl%ncomp
                x(i) = 1.D0/gl%ncomp
            end do
        End if
        gl%molfractions = x
        call reduced_parameters_calc(gl,Temp)
        !                       Call PhaseDet (let the program determine the stable phase equilibrium)
        !                       Try to calculate the VLw equilibrium
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, beta, nrofphases, error)

        if (nrofphases == 1) error = -6666

        if(error == 0) then
            !Calculate the Chem. Pot and Fugacities of the phases
            x_ph1 = x_Phase(:,phasetype(1))     ! vap
            x_ph2 = x_Phase(:,phasetype(2))     ! liq2
            Dens = rho(phasetype(1))
            gl%molfractions = x_ph1                ! vap
            call reduced_parameters_calc(gl,Temp)
            call Chempot_CALC(gl,Temp,Dens, Chempot_ph1, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi_ph1)

            rho_out(1) = rho(phasetype(1))
            rho_out(2) = rho(phasetype(2))
            !write(11,3020) dummy, weight, Temp, press, Chempot_vap(1), Chempot_vap(2), &
            !            & lnfi_vap(1), lnfi_vap(2), x_vap(1), x_vap(2), &
            !            & x_liq2(1), x_liq2(2), x_hyd(1), x_hyd(2), rho(phasetype(1)), rho(phasetype(2))
        end if

    case (2)    !('VHIw')
        !ph1 = vap, ph2 = hyd, ph3 = sol
        !Three phase equilibrium Vapor (V), Hydrate (H) and Solid H2O (Iw)
        if (x(1) == 0.D0) then
            x = 0.d0
            ! Mixed hydrates 2015
            do i = 1,gl%ncomp
                x(i) = 1.d0/gl%ncomp
            end do
        end if
        gl%molfractions = x
        call reduced_parameters_calc(gl,Temp)
        !H2O = 1, CO2 = 2 (IN OUR CASE, IN GENERAL ALL SOLID FORMERS ARE NUMBERED CONSECUTIVELY)
        gl%solidtype(1) = 1
        gl%solidtype(2) = 0    !No Hydrate for the twophase equilibrium
        gl%solid_pos = 1        !Water on position 1
        !IFlash = 3 means, that T and p are given a a two phase equilibrium is calculated
        iFlash = 3
        !Initial values for the vapor phase
        !x_ph1(1) = 0.01D0   ! vap
        !x_ph1(2) = 0.99D0   ! vap
        !x_ph3(1) = 1.D0     ! sol
        !x_ph3(2) = 0.D0     ! sol
        ! Mixed hydrates 2015
        x_ph1(1) = 0.01D0   ! vap (water - low amount in vapor phase)
        do i = 2,gl%ncomp
            x_ph1(i) = (1.D0-x_ph1(1))/(gl%ncomp-1)  ! vap (all components except water)
        end do
        x_ph3 = 0.D0        ! sol (all components except water)
        x_ph3(1) = 1.D0     ! sol (water ice)

        iPhase_sol = 2 !Liqhter (only fluid) phase is vapor

        !Try to calculate the VIw equilibrium
        call ptflash_solid_NC_2P(gl,press, Temp, rho2P, x, x_ph3, x_ph2, x_ph1, &
            & rho_ph1, beta, iFlash, iPhase_sol, error, iter)

        if(error == 0) then
            IPhase = 2
            gl%molfractions = x_ph1    ! vap
            call reduced_parameters_calc(gl,Temp)
            Dens = gl%rho_vap  ! vap!rhomix_calc(gl,Temp, press, 0.D0, IPhase,0)
            call Chempot_CALC(gl,Temp,Dens, Chempot_ph1, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi_ph1)
            rho_ph3 = 1.D0 / dgdp_WaterIce(gl,Temp,press)

            rho_out(1) = Dens
            rho_out(2) = rho_ph3    ! solid = ice
            !write(11,3020) dummy, weight, Temp, press, Chempot_vap(1), Chempot_vap(2), &
            !            & lnfi_vap(1), lnfi_vap(2), x_vap(1), x_vap(2), &
            !            & x_hyd(1), x_hyd(2), x_sol(1), x_sol(2), Dens, rho_sol
        end if

    case (3)    !('VLxH')
        !ph1 = vap, ph2 = liq1, ph3 = hyd
        !Three phase equilibrium Vapor (V), Hydrate (H) and liquid component other than water (e.g. CO2: Lc)
        !        ------------------------------------------------------------------------------------
        !         A carbon dioxide rich vapor phase and a water rich liquid are expected
        !         If no overall composition is given: Set the overall composition to a value that is most likely inbetween the two phases in equilibrium:
        if ((x(1) == 0.D0) .OR. (x(2) == 0.D0)) then
            x = 0.D0
            x(1) = 0.0001D0     ! low amount of water
            ! Mixed hydrates 2015
            do i = 2,gl%ncomp
                x(i) = (1.d0-x(1))/gl%ncomp
            end do
        End if
        !Try to find a two phase equilibrium with the given overall composition
        !If the mixture does not split into two phases, change the overall composition and try again
        Do j = 1, 20
            gl%molfractions = x
            call reduced_parameters_calc(gl,Temp)
            !                       Call PhaseDet (let the program determine the stable phase equilibrium)
            call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, beta, nrofphases, error)
            if (nrofphases == 2) exit
            !x(2) = x(2) - 0.0002D0
            !x(1) = 1.D0 - x(2)
            ! Mixed hydrates 2015
            x(1) = x(1) + 0.0002D0
            do i = 2,gl%ncomp
                x(i) = (1.d0-x(1))/gl%ncomp
            end do
        End Do
        if (error == 0) then
            if (nrofphases == 1) then
                error = -6666
            else
                x_ph1 = x_Phase(:,phasetype(1))     ! vap
                x_ph2 = x_Phase(:,phasetype(2))     ! liq1
            end if
        end if

        !if ((x_ph2(2) < 0.5D0) .or. (x_ph1(2) < 0.5D0)) error = -6667   ! pure hydrates (water + 1 gas)
        if ((x_ph2(1) > 0.5D0) .or. (x_ph1(1) > 0.5D0)) error = -6667   ! mixed hydrates (error = high content of water)

        if(error == 0) then
            gl%molfractions = x_ph1    ! vap
            call reduced_parameters_calc(gl,Temp)
            Dens = rho(phasetype(1))
            call Chempot_CALC(gl,Temp,Dens, Chempot_ph1, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi_ph1)

            rho_out(1) = rho(phasetype(1))
            rho_out(2) = rho(phasetype(2))
            !write(11,3020) dummy, weight, Temp, press, Chempot_vap(1), Chempot_vap(2), &
            !            & lnfi_vap(1), lnfi_vap(2), x_vap(1), x_vap(2), &
            !            & x_liq1(1), x_liq1(2), x_hyd(1), x_hyd(2), rho(phasetype(1)), rho(phasetype(2))
        end if

    case (4)    !('LwLxH')
        !ph1 = liq2, ph2 = liq1, ph3 = hyd
        !Three phase equilibrium Liquid Water (Lw), Hydrate (H) and liquid component other than water (e.g. CO2: Lc)
        !       Two liquids that are rich of water and the guest component are expected
        !       If no overall composition is given: Set the overall composition to a value that is most likely inbetween the two phases in equilibrium:
        if ((x(1) == 0.D0) .OR. (x(2) == 0.D0)) then
            x = 0.D0
            ! Pure hydrates
            !x(1) = 0.5d0    !0.7D0
            !x(2) = 0.5d0    !0.3D0
            ! Mixed hydrates 2015
            do i = 1,gl%ncomp
                x(i) = 1.D0/gl%ncomp
            end do
        End if
        gl%molfractions = x
        call reduced_parameters_calc(gl,Temp)
        !                       Call PhaseDet (let the program determine the stable phase equilibrium)
        call PhaseDet(gl,press, Temp, x, rho, x_Phase, phasetype, beta, nrofphases, error)

        if (nrofphases == 1) error = -6666

        if(error == 0) then
            x_ph2 = x_Phase(:,phasetype(1))    ! liq1
            x_ph1 = x_Phase(:,phasetype(2))    ! liq2
            gl%molfractions = x_ph2    ! liq1
            call reduced_parameters_calc(gl,Temp)
            Dens = rho(phasetype(1))
            call Chempot_CALC(gl,Temp,Dens, Chempot_ph1, 0)
            call lnf_mix(gl,Temp, Dens, press, lnfi_ph1)

            rho_out(1) = rho(phasetype(1))
            rho_out(2) = rho(phasetype(2))
            !write(11,3020) dummy, weight, Temp, press, Chempot_liq1(1), Chempot_liq1(2), &
            !            & lnfi_liq1(1), lnfi_liq1(2), x_liq2(1), x_liq2(2), &
            !            & x_liq1(1), x_liq1(2), x_hyd(1), x_hyd(2), rho(phasetype(1)), rho(phasetype(2))
        end if

        case default
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'unknown phase equilibria type'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        error = 300
    end select

    end subroutine cpo_2phase_at_3phase
    !******************************************************************************
    !******************************************************************************




    !!******************************************************************************
    !subroutine QudruplePoints (path, Q_calc, Q_point, Qexist, Temp_Q, press_Q, &
    !                           & x_ph1, x_ph2, x_ph3, x_ph4)
    !!******************************************************************************
    !! Calculation of all existing quadruple points for pure hydrate
    !




    !
    !implicit none
    !
    !character(len=255) :: path, filename
    !character(len=6) :: EqtypeIn, Eqtype        !Type of Equilibrium, which should be calculated
    !
    !double precision, dimension(30) :: x, Chempot, lnfi
    !double precision:: Temp, press, Dens, press_step, max_press, min_press, rhomix_calc
    !double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est  !, fug_gas
    !double precision, dimension(5):: rho
    !double precision, dimension(3):: rho_sol
    !integer :: error, solidnr, solid, i
    !integer :: iPhase, iFlash, iter, loops, iPhase_sol
    !
    !!Pressures and Temperatures at each Quadrupelpoint
    !character(len=6), dimension(4)  :: Q_point       !Quadruple point
    !double precision, dimension (6) :: press_Q, Temp_Q
    !double precision, dimension(4,30) :: x_ph1, x_ph2, x_ph3, x_ph4, x_phX
    !logical, dimension (4):: Qexist
    !integer, dimension(5) :: error_Q, Q_calc
    !
    !double precision, dimension(3) :: Phasefrac
    !
    !! Mixed hydrates 2015
    !double precision, dimension(3,30):: CiJ, occup
    !double precision, dimension(30):: xH, fug_gas
    !
    !!//////////////////////////////////////////////////////////////////
    !
    !! Q_calc defines which Q-point shall be calculated
    !! 1-VLwHIw, 2-VLwLcH, 3-VLxHIx, 4-LwLxHIx, 5-all existing Q-points and stored in txt file
    !
    !error = 0
    !x = 0.D0
    !x_ph1 = 0.d0
    !x_ph2 = 0.d0
    !x_ph3 = 0.D0
    !x_ph4 = 0.D0
    !x_phX = 0.d0        ! only help variable = 0
    !beta = 0.D0
    !rhoph1_est = 0.D0
    !rhoph2_est = 0.D0
    !rhoph3_est = 0.D0
    !rhophX_est = 0.D0   ! only help variable = 0
    !!rhofluid_est = 0.D0
    !iFlash = 5
    !Phasefrac = 0.D0
    !Qexist = .false.
    !press_Q = 0.d0
    !Temp_Q = 0.d0
    !Q_point = ''
    !
    !!The iPhase_sol can be used to assume fluid phases in two and three phase equilibrium calculation
    !! 0 : Let the algorithms choose the phases
    !! 1 : Liquid phase of the hydrate guest is assumed (e.g. LwH or LwLcH)
    !! 2 : Vapor phase of the hydrate guest is assumed (e.g. VH or VLwH
    !!NOTE: IN case of two fluid phases, the heavier phase is always assumed as liquid phase (no vapor / vapor equilibrium)
    !iPhase_sol = 0
    !
    !! In the following section quadruple points need to be calculated. Depending on the guest molecule, different quadruple
    !! points will need to be found
    !! In general, all Hydrates should have the Quadruple point: VLwHIw
    !! Depending on the guest, the quadruple points VLwLxH, VLcHIx and LwLcHIx can be calculated.
    !! In this case x stands for the guest molecule and Lx for a liquid phase and Ix for a solid phase
    !! of this guest
    !! The algorithms implemented so far only work for binary mixtures!!!
    !! FURTHERMORE WATER HAS ALWAYS TO BE COMPONENT 1!!!
    !!----------------------------------------------------------------------------------------------
    !
    !!---------------------------------------------
    !if ((Q_calc(1) == 1).or.(Q_calc(5) == 1)) then
    !! VLwHIw - ph1 = vap, ph2 = liq2, ph3 = hyd, ph4 = sol
    !!Calculate the quadruple point VLwHIw (Q1)
    !Q_point(1) = 'VLwHIw'
    !!Choose startvalues depending on the hydrate former
    !!Set arbitrary values first
    !temp_Q(1) = 272.D0
    !press_Q(1) = 5.0D0
    !if (Fluidlist_hydrate(2) == "co2") then
    !temp_Q(1) = 271.246D0
    !press_Q(1) = 1.017105D0
    !elseif (Fluidlist_hydrate(2) == "nitrogen") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 15.D0
    !elseif (Fluidlist_hydrate(2) == "oxygen") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 11.D0
    !elseif (Fluidlist_hydrate(2) == "argon") then
    !temp_Q(1) = 273.D0
    !press_Q(1) = 11.D0
    !elseif (Fluidlist_hydrate(2) == "CO") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 13.6D0
    !elseif (Fluidlist_hydrate(2) == "methane") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 2.5D0
    !elseif (Fluidlist_hydrate(2) == "ethane") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 0.45D0
    !elseif (Fluidlist_hydrate(2) == "propane") then
    !temp_Q(1) = 272.D0
    !press_Q(1) = 0.16D0
    !end if
    !
    !solid = 2
    !solidnr = 1
    !x_ph4(1,1) = 1.d0                 ! x_sol(1) = 1.D0
    !x_ph4(2,1) = 0.d0                 ! x_sol(2) = 0.D0
    !x_ph1(1,1) = 0.00064468D0       ! x_vap_Q(1,1) = 0.00064468D0
    !x_ph1(1,2) = 0.99935532D0       ! x_vap_Q(1,2) = 0.99935532D0
    !x_ph2(1,1) = 0.98468113D0       ! x_liq2_Q(1,1) = 0.98468113D0
    !x_ph2(1,2) = 0.01531887D0       ! x_liq2_Q(1,2) = 0.01531887D0
    !call ptflash_sol_2C_4P(press_Q(1), Temp_Q(1), solid, solidnr, 2, x, x_ph4(1,:), x_ph1(1,:), x_phX(1,:), &
    !                     & x_ph2(1,:), rhoph1_est, rhophX_est, rhoph2_est, error_Q(1), iter)
    !
    !    !Calculate the Chem. Pot and Fugacities of the phases
    !    if (error_Q(1) == 0) then
    !        IPhase = 0
    !        molfractions = x_ph1(1,:)
    !        call reduced_parameters_calc(gl,Temp_Q(1))
    !        Dens = rhomix_calc(gl,Temp_Q(1), press_Q(1), 0.D0, IPhase,0)
    !        call lnf_mix(gl,Temp_Q(1), Dens, press_Q(1), lnfi)
    !
    !        !Get the composition of the hydrate phase
    !        ! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
    !        ! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
    !        fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
    !        call hdrt_mole_fract(gl, Temp_Q(1),press_Q(1)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
    !        x_ph3(1,1) = xH(1)
    !        x_ph3(1,2) = xH(2)
    !        !quadruple point exists and was found
    !        Qexist(1) = .true.
    !    else
    !        write(12,*) 'Failed to calculate quadruple point VLwHIw'
    !        return
    !    End if
    !
    !end if      ! end of VLwHIw
    !
    !!---------------------------------------------
    !if ((Q_calc(2) == 1).or.(Q_calc(5) == 1)) then
    !! VLwLxH - ph1 = vap, ph2 = liq2, ph3 = liq1, ph4 = hyd
    !!Quadruple point for systems where liquid liquid equilibria form
    !if ((Fluidlist_hydrate(2) == "co2") .or. (Fluidlist_hydrate(2) == "ethane") .or. (Fluidlist_hydrate(2) == "propane")) then
    !    !Calculate Q2 (VLwLxH)
    !    !Vapor Liquid water hydrate and another liquid in equilibrium
    !    if (Fluidlist_hydrate(2) == "co2") then
    !        Q_point(2) = 'VLwLcH'
    !        press_Q(2) = 4.5D0
    !        temp_Q(2) = 283.0D0
    !    elseif (Fluidlist_hydrate(2) == "ethane") then
    !        Q_point(2) = 'VLwLeH'
    !        press_Q(2) = 3.3D0
    !        temp_Q(2) = 288.0D0
    !    elseif (Fluidlist_hydrate(2) == "propane") then
    !        Q_point(2) = 'VLwLpH'
    !        press_Q(2) = 0.56D0
    !        temp_Q(2) = 278.5D0
    !    end if
    !    solid = 0
    !    solidnr = 1
    !    x_phX(2,1) = 0.d0     !x_sol(1) = 0.D0
    !    x_phX(2,2) = 1.d0     !x_sol(2) = 1.D0
    !    x_ph1(2,1) = 0.001D0    !x_vap
    !    x_ph1(2,2) = 0.999D0
    !    x_ph3(2,1) = 0.002D0    !x_liq1
    !    x_ph3(2,2) = 0.998D0
    !    x_ph2(2,1) = 0.98D0     !x_liq2
    !    x_ph2(2,2) = 0.02D0
    !    call ptflash_sol_2C_4P(press_Q(2), Temp_Q(2), solid, solidnr, 2, x, x_phX(2,:), x_ph1(2,:), x_ph3(2,:), &
    !                         & x_ph2(2,:), rhoph1_est, rhoph3_est, rhoph2_est, error_Q(2), iter)
    !
    !        if (error_Q(2) == 0) then
    !            !Calculate the Chem. Pot and Fugacities of the phases
    !            IPhase = 0
    !            molfractions = x_ph1(2,:)
    !            call reduced_parameters_calc(gl,Temp_Q(2))
    !            Dens = rhomix_calc(gl,Temp_Q(2), press_Q(2), 0.D0, IPhase,0)
    !            call lnf_mix(gl,Temp_Q(2), Dens, press_Q(2), lnfi)
    !
    !            !Get the composition of the hydrate phase
    !            ! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
    !            ! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
    !            fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
    !            call hdrt_mole_fract(gl, Temp_Q(2),press_Q(2)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
    !            x_ph4(2,1) = xH(1)
    !            x_ph4(2,2) = xH(2)
    !            !quadruple point exists and was found
    !            Qexist(2) = .true.
    !        else
    !            write(12,*) 'Failed to calculate quadruple point VLwLxH'
    !            return
    !        End if
    !End if
    !end if      ! end of VLwLxH
    !
    !!---------------------------------------------
    !if ((Q_calc(3) == 1).or.(Q_calc(4) == 1).or.(Q_calc(5) == 1)) then
    !! VLxHIx - ph1 = vap, ph2 = liq1, ph3 = hyd, ph4 = sol
    !!Quadruple point with hydrate and solid guest
    !if (Fluidlist_hydrate(2) == "co2") then
    !    !Calculate Q3 (VLxHIx)
    !    !Vapor, liquid guest, Hydrate and solid guest in equilibrium
    !    if (Fluidlist_hydrate(2) == "co2") then
    !        Q_point(3) = 'VLcHIc'
    !        press_Q(3) = 0.51791D0
    !        temp_Q(3) = 216.59D0
    !    End if
    !    solid = 1
    !    solidnr = 2
    !    x_ph4(3,1) = 0.d0           !x_sol(1) = 0.D0
    !    x_ph4(3,2) = 1.d0           !x_sol(2) = 1.D0
    !    x_ph1(3,1) = 0.00000358D0   !x_vap_Q(3,1) = 0.00000358D0
    !    x_ph1(3,2) = 0.99999642D0   !x_vap_Q(3,2) = 0.99999642D0
    !    x_ph2(3,1) = 0.000052D0     !x_liq1_Q(3,1) = 0.000052D0
    !    x_ph2(3,2) = 0.999948D0     !x_liq1_Q(3,2) = 0.999948D0
    !    x_phX(3,1) = 0.d0           !x_liq2(1) = 0.D0
    !    x_phX(3,2) = 0.d0           !x_liq2(2) = 0.D0
    !    call ptflash_sol_2C_4P(press_Q(3), Temp_Q(3), solid, solidnr, 2, x, x_ph4(3,:), x_ph1(3,:), x_ph2(3,:), &
    !                         & x_phX(3,:), rhoph1_est, rhoph2_est, rhophX_est, error_Q(3), iter)
    !
    !        If (error_Q(3) == 0) then
    !            !Calculate the Chem. Pot and Fugacities of the phases
    !            IPhase = 0
    !            molfractions = x_ph1(3,:)
    !            call reduced_parameters_calc(gl,Temp_Q(3))
    !            Dens = rhomix_calc(gl,Temp_Q(3), press_Q(3), 0.D0, IPhase,0)
    !            call lnf_mix(gl,Temp_Q(3), Dens, press_Q(3), lnfi)
    !
    !            !Get the composition of the hydrate phase
    !            ! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
    !            ! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
    !            fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
    !            call hdrt_mole_fract(gl, Temp_Q(3),press_Q(3)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
    !            x_ph3(3,1) = xH(1)
    !            x_ph3(3,2) = xH(2)
    !            !quadruple point exists and was found
    !            Qexist(3) = .true.
    !        else
    !            write(12,*) 'Failed to calculate quadruple point VLxHIx'
    !            return
    !        End if
    !end if
    !end if      ! end of VLxHIx
    !
    !!---------------------------------------------
    !if ((Q_calc(4) == 1).or.(Q_calc(5) == 1)) then
    !! LwLxHIx - ph1 = liq2, ph2 = liq1, ph3 = hyd, ph4 = sol
    !if (Fluidlist_hydrate(2) == "co2") then
    !!Calculate Q4 (LwLxHIx)
    !!Liquid water, liquid guest, Hydrate and solid guest in equilibrium
    !    if (Fluidlist_hydrate(2) == "co2") then
    !        Q_point(4) = 'LwLcHIc'
    !        press_Q(4) = 450.D0
    !        temp_Q(4) = 280.D0
    !    end if
    !    solid = 1
    !    solidnr = 2
    !    x_ph4(4,1) = 0.D0       !x_sol(1) = 0.D0
    !    x_ph4(4,2) = 1.D0       !x_sol(2) = 1.D0
    !    x_ph2(4,1) = 0.005D0    !x_liq1_Q(4,1) = 0.005D0
    !    x_ph2(4,2) = 0.995D0    !x_liq1_Q(4,2) = 0.995D0
    !    x_ph1(4,1) = 0.950D0    !x_liq2_Q(4,1) = 0.950D0
    !    x_ph1(4,2) = 0.050D0    !x_liq2_Q(4,2) = 0.050D0
    !    x_phX(4,1) = 0.d0       !x_vap(1) = 0.D0
    !    x_phX(4,2) = 0.d0       !x_vap(2) = 0.D0
    !    call ptflash_sol_2C_4P(press_Q(4), Temp_Q(4), solid, solidnr, 2, x, x_ph4(4,:), x_ph1(4,:), x_ph2(4,:), &
    !                         & x_phX(4,:), rhophX_est, rhoph2_est, rhoph1_est, error_Q(4), iter)
    !
    !        If (error_Q(4) == 0) then
    !            !Calculate the Chem. Pot and Fugacities of the phases
    !            IPhase = 0
    !            molfractions = x_ph2(4,:)
    !            call reduced_parameters_calc(gl,Temp_Q(4))
    !            Dens = rhomix_calc(gl,Temp_Q(4), press_Q(4), 0.D0, IPhase,0)
    !            call lnf_mix(gl,Temp_Q(4), Dens, press_Q(4), lnfi)
    !
    !            !Get the composition of the hydrate phase
    !            ! Pure hydrate: fug_gas = dexp(lnfi(2)) * 1.D6
    !            ! Mixed hydrate: ... however this does not really apply in quadruple point subr. (ncomp = 2)
    !            fug_gas(1:nrofhydrateformers-1) = dexp(lnfi(2:nrofhydrateformers))*1.d6
    !            call hdrt_mole_fract(gl, Temp_Q(4),press_Q(4)*1.D6,fug_gas,occup,CiJ,xH,occup_single, occup_double)
    !            x_ph3(4,1) = xH(1)
    !            x_ph3(4,2) = xH(2)
    !            !quadruple point exists and was found
    !            Qexist(4) = .true.
    !        else
    !            write(12,*) 'Failed to calculate quadruple point LwLxHIx'
    !            return
    !        End if
    !End if
    !end if      ! end of LwLxHIx
    !
    !!----------------------------------------------------------------------------------------------
    !! If all existing quadruple points are calculated,
    !! the subroutine writes them in the Q-points.txt file
    !if (Q_calc(5) == 1) then
    !
    !3024    format (a6,' ',f8.4,' ',f10.6,' ',f16.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
    !              &,' ',f12.10,' ',f12.10)
    !
    !open(unit=12, file=trim(path) // 'Q-points.txt', status='unknown', action='write', iostat=error)
    !
    !    do i = 1,4
    !        if (Qexist(i)) then
    !            write(12,3024) Q_point(i), Temp_Q(i), press_Q(i), x_ph1(i,1), x_ph1(i,2), x_ph2(i,1), x_ph2(i,2), &
    !                            & x_ph3(i,1), x_ph3(i,2), x_ph4(i,1), x_ph4(i,2)
    !        end if
    !    end do
    !
    !close(12)
    !end if
    !!-------------------------------------------------------------------------------------------
    !
    !end subroutine QudruplePoints
    !!******************************************************************************
    !!******************************************************************************
    end submodule impl
    !DEC$ END IF