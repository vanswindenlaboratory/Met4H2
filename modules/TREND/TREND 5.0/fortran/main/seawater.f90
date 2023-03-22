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

    ! module for file seawater.f90
    module seawater_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use rhomix_pt_module
    use waterice_module
    use setup_module
    use ancillary_equations_mix_module
    use reduced_parameters_calc_module
    use electrolytes

    contains




    !sewater equation


    Double Precision Function Density_water_diff(gl, p_w, parameters)



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: d_w, p_w        ! Inputvariablen
    Double Precision :: D_sea, T, P_new, M_salt, d_w_new, rho_water, wm_sea
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------


    T = parameters%a_p(1)
    D_sea = parameters%a_p(2)
    wm_sea = parameters%a_p(3)


    d_w_new = rhomix_calc(gl, T, p_w, 0.d0, 1, 1)


    Density_water_diff = (d_w_new * gl%wm(1)) - (1 / ((1 / (D_sea * wm_sea)) - v_saline(gl, T, p_w)))!von Max

    end function Density_water_diff
    !------------------------------------------------------------------------------------------


    !************************************************************************************
    !Function for iteration temperature from ps input

    Double Precision Function entropy_diff_sea(gl, titer, Parameters)



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: d_w, p_w        ! Inputvariablen
    Double Precision :: D_sea, T, P_new, M_salt, d_w_new, rho_water, wm_sea, d_water, p, s_in, s_new, titer
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------


    p = parameters%a_p(1)
    s_in = parameters%a_p(2)
    wm_sea = parameters%a_p(3)

    d_water = rhomix_calc(gl, titer, p, 0.d0, 1, 1)

    s_new = s_sea_calc(gl, titer, p, d_water)

    entropy_diff_sea = s_new - s_in
    !Density_water_diff = (d_w_new) - (1 / ((1 / (D_sea)) - (v_saline(gl, T, p_w) * M_salt)))


    end function entropy_diff_sea
    !------------------------------------------------------------------------------------------------


    !_----------------------------------------------------------------------------------------------
    Double Precision Function ps_input_sea(gl, T, P, s_in)



    implicit none

    type(type_gl) :: gl


    double precision :: s_min
    double precision :: s_max
    double precision :: s_min_allowed
    double precision :: s_max_allowed
    double Precision :: Delta_allowed
    type(type_additional_parameters) :: param_sea
    double Precision :: T, P, D_sea, d_w, p_w, s_in, titer
    integer :: max_iterations
    integer :: iterations
    integer :: errorflag, iter



    s_min = -10.d0
    s_max = 12000.d0
    s_min_allowed = 1.3d0 * s_min
    s_max_allowed = 1.3d0 * s_max

    !parameters_sea = 0.d0
    errorflag=0

    !Parameters for iteration:
    Delta_allowed = 1.d-11
    Max_Iterations = 100
    param_sea%a_p(1) = p
    param_sea%a_p(2) = s_in
    param_sea%a_p(3) = molar_sea(gl)
    titer = t

    CALL Regula_Falsi(gl, entropy_diff_sea, titer, s_min, s_max, Delta_allowed,&
        & s_min_allowed, s_max_allowed, Max_iterations, Iterations, errorflag, param_sea)


    if(errorflag .eq. 0) then
        ps_input_sea = titer
    else
        ps_input_sea = errorflag
    end if


    End Function ps_input_sea
    !----------------------------------------------------------------------------------------------------------------------------------




    !-------------------------------------------------------------------------------------------
    DOUBLE PRECISION FUNCTION G_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, g_water, press

    !!gl%seacalc = .false.
    d = d!rhomix_calc(gl, t, p, 0.d0, 1, 1)
    g_water = g_calc(gl, T, D, 1)
    press = p!p_calc(gl, t, d, 1)
    !gl%seacalc = .true.
    g_water = g_water / gl%wm(1)


    G_SEA_CALC = g_water + g_saline(gl, T, P)

    !call convert_single_prop(gl, 2, g_sea_calc, 9999)

    END FUNCTION G_SEA_CALC
    !********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION D_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, d_water

    !d = rhomix_calc(gl, t, p, 0.d0, 1, 1)  ! wie war das mit reduced parameters?!

    d_water = d * gl%wm(1)

    d_water = rhomix_calc(gl, t, p, 0.d0, 1, 1) * gl%wm(1)

    D_SEA_CALC = 1/((1 / d_water) + v_saline(gl, T, P))

    !call convert_single_prop(gl, 1, d_sea_calc, 9999)

    END FUNCTION D_SEA_CALC
    !********************************************************************************


    !!********************************************************************************
    !DOUBLE PRECISION FUNCTION V_SEA_CALC(gl, T, P, D)
    !!achtung, nur spezifisches Volumen des Salzteils !!!!


    !implicit none
    !
    !type(type_gl) :: gl
    !
    !DOUBLE PRECISION :: T, P, D, SAL, mol_w
    !
    !!V_SEA_CALC = v_saline(gl, T, P)
    !
    !END FUNCTION V_SEA_CALC
    !!********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION U_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, u_water, u_saline

    !gl%seacalc = .false.

    u_water = u_calc (gl, T, D, 1)

    !gl%seacalc = .true.
    u_water = u_water / gl%wm(1)

    u_saline =  g_saline(gl, T, P) - (T * DGDT_saline(gl, T, P)) - ((P * 1000000.d0) * v_saline(gl, T, P))

    U_SEA_CALC = u_water + u_saline

    !call convert_single_prop(gl, 2, u_sea_calc, 9999)

    END FUNCTION U_SEA_CALC
    !********************************************************************************


    !********************************************************************************
    double precision function u_saline(gl,t,p,d)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, u_water

    u_saline =  g_saline(gl, T, P) - (T * DGDT_saline(gl, T, P)) - ((P * 1000000.d0) * v_saline(gl, T, P))

    end function
    !********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION S_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, s_water

    ! gl%seacalc = .false.

    s_water = s_calc (gl, T, D, 1)

    !gl%seacalc = .true.

    s_water = s_water / gl%wm(1)

    S_SEA_CALC = s_water - DGDT_saline(gl, T, P)

    !call convert_single_prop(gl, 2, s_sea_calc, 9999)

    END FUNCTION S_SEA_CALC
    !********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION H_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, h_water
    logical :: seawater

    !gl%seacalc = .false.

    h_water = h_calc (gl, T, D, 1)
    !
    ! gl%seacalc = .true.

    h_water = h_water / gl%wm(1)
    !
    H_SEA_CALC = h_water + h_saline(gl, T, P)

    !call convert_single_prop(gl, 2, h_sea_calc, 9999)

    END FUNCTION H_SEA_CALC
    !********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION h_saline(gl, T, P)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P

    h_saline =  g_saline(gl, T, P) - (T * DGDT_saline(gl, T, P))

    END FUNCTION h_saline
    !********************************************************************************


    !********************************************************************************
    DOUBLE PRECISION FUNCTION A_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, a_water, a_saline,g, v, press
    logical :: seawater

    gl%seacalc = .false.

    a_water = a_calc (gl, T, D, 1)

    !gl%seacalc = .true.

    a_water = a_water / gl%wm(1)

    g = g_saline(gl, T, P)

    v = v_saline(gl, t, p)

    press = p * 1.d6


    a_saline = g - press*v

    A_SEA_CALC = a_water + a_saline

    !call convert_single_prop(gl, 2, a_sea_calc, 9999)


    END FUNCTION A_SEA_CALC
    !**************************************************************************************


    !**************************************************************************************
    double precision function a_saline(gl, t, p)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, v, press, gp, g

    gp = v_saline(gl, t, p)
    g = g_saline(gl,t , p)

    a_saline = g - (p*1.d6*gp)

    end function a_Saline
    !**************************************************************************************


    !***************************************************************************************
    DOUBLE PRECISION FUNCTION CP_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, cp_water
    logical :: seawater

    !gl%seacalc = .false.

    cp_water = cp_calc (gl, T, D, 1)

    !!gl%seacalc = .true.

    cp_water = cp_water / gl%wm(1)

    CP_SEA_CALC = cp_water + (- T) * D2GDT2_saline(gl, T, P)

    !call convert_single_prop(gl, 2, cp_sea_calc, 9999)

    END FUNCTION CP_SEA_CALC

    DOUBLE PRECISION FUNCTION cp_Saline(gl, T, P)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D
    logical :: seawater

    cp_Saline = (- T) * D2GDT2_saline(gl, T, P)

    END FUNCTION cp_saline


    DOUBLE PRECISION FUNCTION cv_Saline(gl, T,P)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D
    logical :: seawater

    cv_Saline = (T) * (d2gdtdp_saline(gl, t, p)**2.d0 - d2gdt2_saline(gl, t, p) * d2gdp2_saline(gl, t, p)) / d2gdp2_saline(gl, t, p)

    END FUNCTION cv_saline
    !***************************************************************

    !***************************************************************
    DOUBLE PRECISION FUNCTION WS_SEA_CALC(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, dgdp, d2gdt2, d2gdp2, d2gdtdp, ws_saline, ws_water
    logical :: seawater, seacalc
    !D2GDP2_w = D2GDP2_water(gl, T, D, 1) / gl%wm(1)
    !D2GDP2_s = D2GDP2_saline(gl, T, P)


    !DGDP_sea = v_saline(gl, T, P) + (1 / (d * gl%wm(1)))
    !D2GDT2_sea = (D2GDT2_water(gl, T, D, 1) / gl%wm(1)) + D2GDT2_saline(gl, T, P)
    !D2GDTDP_sea = (D2GDTDP_water(gl, T, D, 1) / gl%wm(1)) + D2GDTDP_saline(gl, T, P)
    !D2GDP2_sea = (D2GDP2_water(gl, T, D, 1) / gl%wm(1)) + D2GDP2_saline(gl, T, P)

    dgdp = DGDP_sea(gl, T,p, D)!/gl%sea%wm_sea
    d2gdt2 = D2GDT2_sea(gl, T,p, D)!/gl%sea%wm_sea
    d2gdtdp = D2GDTDP_sea(gl, T, p,D)!/gl%sea%wm_sea
    d2gdp2 = D2GDP2_sea(gl, T, P, D)!/gl%sea%wm_sea

    WS_SEA_CALC = DGDP * sqrt(D2GDT2 / ((D2GDTDP ** 2.d0) - (D2GDT2 * D2GDP2)))    !das wäre die spezifische variante

    ws_sea_calc = DGDP  * sqrt((D2GDT2 / ((D2GDTDP ** 2.d0) - (D2GDT2 * D2GDP2)))/gl%wm(1))   !das hier funktioniert in molar.....

    !gl%seacalc = .false.

    ws_water = ws_calc(gl, t, d,1)

    !gl%seacalc = .true.

    dgdp = v_saline(gl, t, p)
    d2gdt2 = d2gdt2_saline(gl, t, p)
    d2gdtdp = d2gdtdp_saline(gl, t, p)
    d2gdp2 = d2gdp2_saline(gl, t, p)

    ws_saline = dgdp * sqrt( (d2gdt2 / ((d2gdtdp**2.d0) - (d2gdt2 * d2gdp2))))

    ws_saline = ws_saline * ws_saline

    dgdp = DGDP_sea(gl, T,p, D) - v_saline(gl, t, p)
    d2gdt2 = D2GDT2_sea(gl, T,p, D) - d2gdt2_saline(gl, t, p)
    d2gdtdp = D2GDTDP_sea(gl, T, p,D) - d2gdtdp_saline(gl, t, p)
    d2gdp2 = D2GDP2_sea(gl, T, P, D) - d2gdp2_saline(gl, t, p)

    ws_water = dgdp * sqrt( (d2gdt2 / ((d2gdtdp**2.d0) - (d2gdt2 * d2gdp2))))


    !Ws_sea_calc =( ( dgdp**2.d0 ) * d2gdt2 ) / (d2gdtdp - (d2gdt2 * d2gdp2))

    !    ws_sea_calc = ((d2gdtdp/(d2gdtdp-(d2gdt2*d2gdp2)))*(dgdp**2))**0.5

    END FUNCTION WS_SEA_CALC
    !************************************************************************************************


    !************************************************************************************************
    DOUBLE PRECISION FUNCTION WS_SEA_Mix(gl, T, P, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, P, D, dgdp&
        &, d2gdt2, d2gdp2, d2gdtdp, dgdpl, d2gdt2l, d2gdtdpl, d2gdp2l, wmmix

    logical :: seawater, seacalc

    call wm_mix_calc(gl, wmmix)

    dgdpl = DGDP_mix(gl, T,p, D) / wmmix
    d2gdt2l = D2GDT2_mix(gl, T,p, D)/ wmmix!/gl%sea%wm_sea
    d2gdtdpl = D2GDTDP_mix(gl, T, p,D)!* wmmix!/gl%sea%wm_sea
    d2gdp2l = D2GDP2_mix(gl, T, P, D)!* wmmix!/gl%sea%wm_sea



    dgdp = dgdpl + v_saline(gl, t, p)
    d2gdt2 = d2gdt2l + d2gdt2_saline(gl, t, p)
    d2gdtdp = d2gdtdpl + d2gdtdp_Saline(gl, t, p)
    d2gdp2 = d2gdp2l + d2gdp2_saline(gl, t, p)

    WS_SEA_mix = DGDP * sqrt(D2GDT2 / ((D2GDTDP ** 2.d0) - (D2GDT2 * D2GDP2)))    !das wäre die spezifische variante

    ws_sea_mix = DGDP  * sqrt((D2GDT2 / ((D2GDTDP ** 2.d0) - (D2GDT2 * D2GDP2)))/wmmix)   !das hier funktioniert in molar.....



    END FUNCTION WS_SEA_mix
    !*******************************************************************************************


    !****************************************************************************************
    DOUBLE PRECISION FUNCTION chempot_SEA_CALC(gl, T, P, D)
    !calculates chemical potential of water in seawater


    implicit none

    type(type_gl) :: gl
    logical :: seawater
    DOUBLE PRECISION :: T, P, D, d_w, dgds, chemw
    double precision, dimension(30) :: chempotv
    !double precision, dimension(30) :: chempot

    !Chempot_CALC(gl,T,D, Chempot, 0)  !hardcoded oir=0 for water
    dgds = DGDS_saline(gl, t, p)
    ! gl%seacalc = .false.
    d_w = d!rhomix_Calc(gl, t, p, 0.d0, 1, 1)
    chemw = (g_sea_calc(gl, t,p, d_w)) !/ 0.018015268d0)!gl%wm(1))
    !gl%seacalc = .true.
    chempot_SEA_CALC = chemw - (gl%sea%salinity * DGDS) !+ chempot_water


    END FUNCTION chempot_SEA_CALC
    !*****************************************************************************

    !***********************************************************************
    double precision function d3gds2dp_saline(gl, t, p)



    implicit none

    type(type_gl) :: gl
    double precision :: numfac, salsave, salp, sall, d2gdsdph, d2gdsdpl, t, p

    salsave = gl%sea%salinity

    numfac = 1.d-7

    salp = gl%sea%salinity * (1.d0 + numfac)

    sall = gl%sea%salinity * (1.d0 - numfac)

    gl%sea%salinity = salp

    d2gdsdph = d2gdsdp_saline(gl, t, p)

    gl%sea%salinity = sall

    d2gdsdpl = d2gdsdp_saline(gl, t, p)

    d3gds2dp_saline = (d2gdsdph - d2gdsdpl) / (salp - sall)

    gl%sea%salinity = salsave

    end function d3gds2dp_saline
    !**********************************************************************

    !***********************************************************************
    double precision function d2gds2_saline_num(gl, t, p)



    implicit none

    type(type_gl) :: gl
    double precision :: numfac, salsave, salp, sall, d2gds2h, d2gds2l, t, p

    salsave = gl%sea%salinity

    numfac = 1.d-4

    salp = gl%sea%salinity * (1.d0 + numfac)

    sall = gl%sea%salinity * (1.d0 - numfac)

    gl%sea%salinity = salp

    d2gds2h = dgds_saline(gl, t, p)

    gl%sea%salinity = sall

    d2gds2l = dgds_saline(gl, t, p)

    d2gds2_saline_num = (d2gds2h - d2gds2l) / (salp - sall)

    gl%sea%salinity = salsave

    end function d2gds2_saline_num
    !**********************************************************************

    !****************************************************************************
    double precision function d2gdsdp_saline(gl, t, p)


    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_IJK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K

    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_IJK = 0.d0

    do i = 2,7
        do j = 0,6
            do k = 1,5
                GIJK_sum_IJK = GIJK_sum_IJK + ( i * k * gl%sea%GIJK_list( i + 1, j + 1, k + 1) * x**(i-2) * y**j * z**(k-1) )
            end do !k
        end do !j
    end do !i

    d2gdsdp_saline = GIJK_sum_IJK / (2.d0 * sal_r * p_r)

    end function d2gdsdp_saline
    !*************************************************************************

    !*******************************************************************
    !Original
    Double Precision Function d_w_new(gl, T, P, D_sea)



    implicit none

    type(type_gl) :: gl


    double precision :: p_min
    double precision :: p_max
    double precision :: p_min_allowed
    double precision :: p_max_allowed
    double Precision :: Delta_allowed
    double Precision :: T, P, D_sea, d_w, p_w, p_w_new
    type(type_additional_parameters) :: param_sea
    integer :: max_iterations
    integer :: iterations
    integer :: errorflag, iter


    !MB das d-7 ist nicht eine Begrenzung der Genauigkeit nur eine Erweiterung des Gültigkeitsbereches der einer mit der
    p_min = gl%sea%pmin_seawater !gl%sea%pmin_seawater
    p_max = 1.d0 * gl%sea%pmax_seawater !+ 1.d-7                !MB Weil die ungenauigkeit wegen der Umrechnung mit der Molmasse bei d-8 liegt, wird hier d-7 dazu gerechnet
    p_min_allowed = gl%sea%pmin_seawater ! * gl%sea%pmax_seawater) !gl%sea%pmin_seawater                !MB d-8 hat auch noch funktioniert, aber vorsichtshalber d-7
    p_max_allowed = 1.d0 * gl%sea%pmax_seawater !+ 1.d-7        !MB Beim TD input wird mit der Molmasse (zwischen mlar und spez.) umgerechnet wegen kontrolwerten

    !parameters_sea = 0.d0
    Delta_allowed = 0.d0
    errorflag=0

    !Parameters for iteration:
    Delta_allowed = 1.d-11
    Max_Iterations = 100
    param_sea%a_p(1) = T
    param_sea%a_p(2) = D_Sea
    param_sea%a_p(3) = molar_sea(gl)


    CALL Regula_Falsi(gl, Density_water_diff, p_w, p_min, p_max, Delta_allowed,&
        & p_min_allowed, p_max_allowed, Max_iterations, Iterations, errorflag, param_sea)


    !MB hier einen geeigneten Fehlerabfang implementieren
    if (errorflag == 0) then
        p_w_new = p_w
        gl%sea%seap = p_w
        d_w_new = rhomix_calc(gl, T, p_w, 0.d0, 1, 1)
    else
        if (errorflag == 1) then
            d_w_new = -1.d0
            gl%sea%Regula_temp = d_w_new
        elseif (errorflag == 2) then
            d_w_new = -2.d0
            gl%sea%Regula_temp = d_w_new
        elseif (errorflag == 3) then
            d_w_new = -3.d0
            gl%sea%Regula_temp = d_w_new
        else
            d_w_new = -123
            gl%sea%Regula_temp = d_w_new
        endif
        !gl%Regula_temp = T
        !gl%regula_press = p
        gl%sea%regula_rho = rhomix_calc(gl, T, P, 0.d0, 1, 1)
        !pause
    endif

    End Function d_w_new
    !_----------------------------------------------------------------------------


    !-----------------------------------------------------------------------------
    Double Precision Function p_w_new(gl, T, P, D_sea)




    implicit none

    type(type_gl) :: gl


    ! needed to transmit T, p to Regula Falsi:
    double precision :: p_min
    double precision :: p_max
    double precision :: p_min_allowed
    double precision :: p_max_allowed
    double Precision :: Delta_allowed
    double Precision :: T, P, D_sea, d_w, p_w
    type(type_additional_parameters) :: param_sea
    integer :: max_iterations
    integer :: iterations
    integer :: errorflag, iter

    !MB das d-7 ist nicht eine Begrenzung der Genauigkeit nur eine Erweiterung des Gültigkeitsbereches der einer mit der
    p_min = gl%sea%pmin_seawater
    p_max = gl%sea%pmax_seawater !+ 1.d-7                !MB Weil die ungenauigkeit wegen der Umrechnung mit der Molmasse bei d-8 liegt, wird hier d-7 dazu gerechnet
    p_min_allowed = gl%sea%pmin_seawater *2.d0               !MB d-8 hat auch noch funktioniert, aber vorsichtshalber d-7
    p_max_allowed = gl%sea%pmax_seawater * 2.d0!+ 1.d-7        !MB Beim TD input wird mit der Molmasse (zwischen mlar und spez.) umgerechnet wegen kontrolwerten

    !parameters_sea = 0.d0
    Delta_allowed = 0.d0
    errorflag=0

    !Parameters for iteration:
    Delta_allowed = 1.d-11
    Max_Iterations = 100
    param_sea%a_p(1) = T
    param_sea%a_p(2) = D_Sea
    param_sea%a_p(3) = molar_sea(gl)

    if(p_w == 0.d0) p_w = 0.1d0

    CALL Regula_Falsi(gl, Density_water_diff, p_w, p_min, p_max, Delta_allowed,&
        & p_min_allowed, p_max_allowed, Max_iterations, Iterations, errorflag, param_sea)

    if (errorflag == 0) then
        p_w_new = p_w
        !d_w_new = rhomix_calc(gl, T, p_w, 0.d0, 1, 1)
    else
        !write(*,*) "fehler in regulappp"
        !pause
    endif

    End Function p_w_new


    subroutine set_coef_Seawater (gl)
    !----------------------------------------------------------------------------
    !This Subroutine sets the coefficients of the Gibbs-Seawater Equation as a global type variable (list)
    !-----------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl
    gl%sea%GIJK_list = 0.d0

    !Die 41 Gibbs-Koefizienten für den reinen Wasseranteil
    gl%sea%GIJK_list(1, 3, 1) = -12357.785933039d0               !Anfang-Koeffizienten für Gibbs-wasseranteil, which ARE NOT USED in this TREND implementation but
    gl%sea%GIJK_list(1, 4, 1) = 736.741204151612d0               !for the sake of completeness listed
    gl%sea%GIJK_list(1, 5, 1) = -148.185936433658d0              !Koefizienten
    gl%sea%GIJK_list(1, 6, 1) = 58.0259125842571d0               !Koefizienten
    gl%sea%GIJK_list(1, 7, 1) = -18.9843846514172d0              !Koefizienten
    gl%sea%GIJK_list(1, 8, 1) = 3.05081646487967d0               !Koefizienten
    gl%sea%GIJK_list(1, 1, 2) = 100015.695367145d0               !Koefizienten
    gl%sea%GIJK_list(1, 2, 2) = -270.983805184062d0              !Koefizienten
    gl%sea%GIJK_list(1, 3, 2) = 1455.0364540468d0                !Koefizienten
    gl%sea%GIJK_list(1, 4, 2) = -672.50778314507d0               !Koefizienten
    gl%sea%GIJK_list(1, 5, 2) = 397.968445406972d0               !Koefizienten
    gl%sea%GIJK_list(1, 6, 2) = -194.618310617595d0              !Koefizienten
    gl%sea%GIJK_list(1, 7, 2) = 63.5113936641785d0               !Koefizienten
    gl%sea%GIJK_list(1, 8, 2) = -9.63108119393062d0              !Koefizienten
    gl%sea%GIJK_list(1, 1, 3) = -2544.5765420363d0               !Koefizienten
    gl%sea%GIJK_list(1, 2, 3) = 776.153611613101d0               !Koefizienten
    gl%sea%GIJK_list(1, 3, 3) = -756.558385769359d0              !Koefizienten
    gl%sea%GIJK_list(1, 4, 3) = 499.360390819152d0               !Koefizienten
    gl%sea%GIJK_list(1, 5, 3) = -301.815380621876d0              !Koefizienten
    gl%sea%GIJK_list(1, 6, 3) = 120.520654902025d0               !Koefizienten
    gl%sea%GIJK_list(1, 7, 3) = -22.2897317140459d0              !Koefizienten
    gl%sea%GIJK_list(1, 1, 4) = 284.517778446287d0               !Koefizienten
    gl%sea%GIJK_list(1, 2, 4) = -196.51255088122d0               !Koefizienten
    gl%sea%GIJK_list(1, 3, 4) = 273.479662323528d0               !Koefizienten
    gl%sea%GIJK_list(1, 4, 4) = -239.545330654412d0              !Koefizienten
    gl%sea%GIJK_list(1, 5, 4) = 152.196371733841d0               !Koefizienten
    gl%sea%GIJK_list(1, 6, 4) = -55.2723052340152d0              !Koefizienten
    gl%sea%GIJK_list(1, 7, 4) = 8.17060541818112d0               !Koefizienten
    gl%sea%GIJK_list(1, 1, 5) = -33.3146754253611d0              !Koefizienten
    gl%sea%GIJK_list(1, 2, 5) = 28.9796526294175d0               !Koefizienten
    gl%sea%GIJK_list(1, 3, 5) = -55.5604063817218d0              !Koefizienten
    gl%sea%GIJK_list(1, 4, 5) = 48.8012518593872d0               !Koefizienten
    gl%sea%GIJK_list(1, 5, 5) = -26.3748377232802d0              !Koefizienten
    gl%sea%GIJK_list(1, 6, 5) = 6.48190668077221d0               !Koefizienten
    gl%sea%GIJK_list(1, 1, 6) = 4.20263108803084d0               !Koefizienten
    gl%sea%GIJK_list(1, 2, 6) = -2.13290083518327d0              !Koefizienten
    gl%sea%GIJK_list(1, 3, 6) = 4.34420671917197d0               !Koefizienten
    gl%sea%GIJK_list(1, 4, 6) = -1.66307106208905d0              !Koefizienten
    gl%sea%GIJK_list(1, 1, 7) = -0.546428511471039d0             !Ende-Koefizienten
    !IAPWS-95 reference state condition
    !energy = 0 and entropy = 0 at the triple point:
    !gl%sea%GIJK_list(1, 1, 1) = 101.342743139672d0
    !gl%sea%GIJK_list(1, 2, 1) = 5.90578348518236d0

    !Oder doch lieber alternativ diese Werte?
    !!quadruple precision values (D.G.Wright 21 July 2008)
    !!gl%sea%GIJK_list(0, 0, 0) = 1.013427431396741480431228220832E2
    !!gl%sea%GIJK_list(0, 1, 0) = 5.905783479094018366702121889468E0
    gl%sea%GIJK_list(1, 1, 1) = 101.342743139674d0
    gl%sea%GIJK_list(1, 2, 1) = 5.90578347909402d0

    !Die 64 Gibbs-Koefizienten für den Salzgehaltanteil
    !Coefficients from the IAPWS Release 2008
    !-----------------------------------------------
    !Logarithmic terms from Reference Composition 2008, computed 01 Aug 2007
    !Deep-Sea Research I 55(2008)50-72.
    gl%sea%GIJK_list(2, 1, 1) = 5812.81456626732d0
    gl%sea%GIJK_list(2, 2, 1) = 851.226734946706d0

    !Seawater Reference State defined by WG127, computed 03 March 2008
    gl%sea%GIJK_list(3, 1, 1) = 1416.27648484197d0  !computed from a quadruple-precision implementation
    gl%sea%GIJK_list(3, 2, 1) = 168.072408311545d0

    !Thermal and colligative properties at 101325 Pa, computed 01 Aug 2007
    gl%sea%GIJK_list(4, 1, 1) = -2432.14662381794d0
    gl%sea%GIJK_list(5, 1, 1) = 2025.80115603697d0
    gl%sea%GIJK_list(6, 1, 1) = -1091.66841042967d0
    gl%sea%GIJK_list(7, 1, 1) = 374.60123787784d0
    gl%sea%GIJK_list(8, 1, 1) = -48.5891069025409d0
    gl%sea%GIJK_list(4, 2, 1) = -493.407510141682d0
    gl%sea%GIJK_list(5, 2, 1) = 543.835333000098d0
    gl%sea%GIJK_list(6, 2, 1) = -196.028306689776d0
    gl%sea%GIJK_list(7, 2, 1) = 36.7571622995805d0
    gl%sea%GIJK_list(3, 3, 1) = 880.031352997204d0
    gl%sea%GIJK_list(4, 3, 1) = -43.0664675978042d0
    gl%sea%GIJK_list(5, 3, 1) = -68.5572509204491d0
    gl%sea%GIJK_list(3, 4, 1) = -225.267649263401d0
    gl%sea%GIJK_list(4, 4, 1) = -10.0227370861875d0
    gl%sea%GIJK_list(5, 4, 1) = 49.3667694856254d0
    gl%sea%GIJK_list(3, 5, 1) = 91.4260447751259d0
    gl%sea%GIJK_list(4, 5, 1) = 0.875600661808945d0
    gl%sea%GIJK_list(5, 5, 1) = -17.1397577419788d0
    gl%sea%GIJK_list(3, 6, 1) = -21.6603240875311d0
    gl%sea%GIJK_list(5, 6, 1) = 2.49697009569508d0
    gl%sea%GIJK_list(3, 7, 1) = 2.13016970847183d0

    !coefficients of the pressure part of the 2003 Gibbs function
    !Progress in Oceanography 58(2003)43-114
    gl%sea%GIJK_list(3, 1, 2) = -3310.49154044839d0
    gl%sea%GIJK_list(4, 1, 2) = 199.459603073901d0
    gl%sea%GIJK_list(5, 1, 2) = -54.7919133532887d0
    gl%sea%GIJK_list(6, 1, 2) = 36.0284195611086d0
    gl%sea%GIJK_list(3, 2, 2) = 729.116529735046d0
    gl%sea%GIJK_list(4, 2, 2) = -175.292041186547d0
    gl%sea%GIJK_list(5, 2, 2) = -22.6683558512829d0
    gl%sea%GIJK_list(3, 3, 2) = -860.764303783977d0
    gl%sea%GIJK_list(4, 3, 2) = 383.058066002476d0
    gl%sea%GIJK_list(3, 4, 2) = 694.244814133268d0
    gl%sea%GIJK_list(4, 4, 2) = -460.319931801257d0
    gl%sea%GIJK_list(3, 5, 2) = -297.728741987187d0
    gl%sea%GIJK_list(4, 5, 2) = 234.565187611355d0

    gl%sea%GIJK_list(3, 1, 3) = 384.794152978599d0
    gl%sea%GIJK_list(4, 1, 3) = -52.2940909281335d0
    gl%sea%GIJK_list(5, 1, 3) = -4.08193978912261d0
    gl%sea%GIJK_list(3, 2, 3) = -343.956902961561d0
    gl%sea%GIJK_list(4, 2, 3) = 83.1923927801819d0
    gl%sea%GIJK_list(3, 3, 3) = 337.409530269367d0
    gl%sea%GIJK_list(4, 3, 3) = -54.1917262517112d0
    gl%sea%GIJK_list(3, 4, 3) = -204.889641964903d0
    gl%sea%GIJK_list(3, 5, 3) = 74.726141138756d0

    gl%sea%GIJK_list(3, 1, 4) = -96.5324320107458d0
    gl%sea%GIJK_list(4, 1, 4) = 68.0444942726459d0
    gl%sea%GIJK_list(5, 1, 4) = -30.1755111971161d0
    gl%sea%GIJK_list(3, 2, 4) = 124.687671116248d0
    gl%sea%GIJK_list(4, 2, 4) = -29.483064349429d0
    gl%sea%GIJK_list(3, 3, 4) = -178.314556207638d0
    gl%sea%GIJK_list(4, 3, 4) = 25.6398487389914d0
    gl%sea%GIJK_list(3, 4, 4) = 113.561697840594d0
    gl%sea%GIJK_list(3, 5, 4) = -36.4872919001588d0

    gl%sea%GIJK_list(3, 1, 5) = 15.8408172766824d0
    gl%sea%GIJK_list(4, 1, 5) = -3.41251932441282d0
    gl%sea%GIJK_list(3, 2, 5) = -31.656964386073d0
    gl%sea%GIJK_list(3, 3, 5) = 44.2040358308d0
    gl%sea%GIJK_list(3, 4, 5) = -11.1282734326413d0

    gl%sea%GIJK_list(3, 1, 6) = -2.62480156590992d0
    gl%sea%GIJK_list(3, 2, 6) = 7.04658803315449d0
    gl%sea%GIJK_list(3, 3, 6) = -7.92001547211682d0

    end subroutine set_coef_Seawater

    subroutine seawater_limits(gl)



    implicit none

    type(type_gl) :: gl

    gl%sea%tmin_seawater = 261.0d0
    gl%sea%tmax_seawater = 353.d0
    gl%sea%seapmin_seawater = 0.d0
    gl%sea%pmax_seawater = 100.0d0 * 1.1d0 !MB Inreased for the uncertainty of an "TD" input compared to an "TP" input by 100MPa (.d-8 would be also possible)
    gl%sea%salmin_seawater = 0.d0
    gl%sea%salmax_seawater = 0.12d0 !kg/kg


    end subroutine seawater_limits


    !*******************************************************************************************

    subroutine check_seawater_limits(gl, temp, press, errorflag)



    implicit none

    type(type_gl) :: gl

    double precision:: temp, press
    integer:: errorflag

    errorflag = 0


    if (gl%sea%salinity .gt. 0.12d0) then
        errorflag = -9936.d0

    elseif((gl%sea%salinity .gt. 0.05d0)) then

        if((temp .lt. 261.d0) .or. (temp .gt. 353.d0)) then
            errorflag = -9936.d0
        elseif((press .lt. 0.09d0) .or. (press .gt. 0.102d0)) then
            errorflag = -9936.d0
        end if

    elseif(gl%sea%salinity .gt. 0.042d0) then

        if(temp .gt. 313.d0) then
            errorflag = -9936.d0
        elseif(press .gt. 0.102d0) then
            errorflag = -9936.d0
        end if

    elseif(temp .gt. 313.d0) then
        if(gl%sea%salinity .gt. 0.042d0) then
            errorflag = -9936.d0
        elseif((press .le. 0.1d0) .and. (100.d0 .le. press)) then
            errorflag = -9936.d0
        end if
    end if

    end subroutine check_seawater_limits

    !*********************************************************************************************

    subroutine setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)
    !----------------------------------------------------------------------------
    !This Subroutine sets the reduced variables/parameters of seawater for the calculation routines of the saline part
    !-----------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    !double precision, dimension(8,8,7) :: GIJK_list
    double precision :: sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea, sal_n_molar, sal_r_molar, x_molar, sal_molar

    sal_n = 0.03516504d0
    sal_r = sal_n * (40.0d0 / 35.0d0)


    sal_n_molar = 0.020479686d0
    sal_r_molar = 0.023456358d0 ! mol/mol

    t_r = 40.d0 ! Kelvin
    p_r = 100000000.d0 !Pa
    sal = gl%sea%salinity                           !Salinity in kg/kg
    x = sqrt(sal/sal_r)                         !reduced salinity
    !sal_molar = 1.d0
    !x_molar = sqrt(sal_molar / sal_r_molar)
    y = (T - 273.15d0) / t_r                    !reduced temperature
    z = (P - 0.101325d0) / (p_r / 1000000.d0)   !reduced Pressure

    !The Molar Mass of Seawater
    wm_sea = molar_sea(gl)
    ! gl%sea%sal_m  = sal_molar
    !gl%molcoef = wm_sea / gl%wm(1)
    end subroutine setvariables
    !----------------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------------
    double precision Function molar_sea(gl)
    !----------------------------------------------------------------------------
    !This Subroutine calculates the molar mass of seawater with respect to the salinity
    !-----------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    !double precision, dimension(8,8,7) :: GIJK_list
    double precision :: sal, M_salt, mass_w, mass_s, x_s, x_w

    sal = gl%sea%salinity

    !The Molar Mass of Seawater

    M_salt = 0.0314038218d0 !Molar Mass of Salt
    molar_sea = ((1.d0 / ((1.d0 - sal) / sal) + 1.d0) / ((sal / ((1.d0 - sal) * M_salt)) + (1.d0 / gl%wm(1))))
    gl%sea%wm_sea = molar_sea

    ! ende max

    mass_w = 1.d0 - sal
    mass_s = sal

    gl%sea%x_s = (mass_s / M_salt) / ( (mass_s / M_salt ) + (mass_w /gl%wm(1) ) )

    gl%sea%x_w = 1.d0 - gl%sea%x_s

    molar_sea = gl%sea%x_w * gl%wm(1) + gl%sea%x_s * M_salt

    end function molar_sea
    !***************************************************************************************


    !***************************************************************************************
    double precision Function molality_sea(gl)
    !----------------------------------------------------------------------------
    !This Subroutine calculates the molar mass of seawater with respect to the salinity
    !-----------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    !double precision, dimension(8,8,7) :: GIJK_list
    double precision :: sal, M_salt, mass_w, mass_s, x_s, x_w

    M_salt = 0.0314038218d0 !Molar Mass of Salt

    molality_sea = gl%sea%salinity / ( ( 1.d0 - gl%sea%salinity) *M_salt )


    end function molality_sea





    DOUBLE PRECISION FUNCTION g_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the Gibbs Free Energy of the Saline Part of Seawater
    !   in J / kg
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision :: GIJK_SUMS, GIJK_SUMB
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea, log_x
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUMS = 0.d0
    GIJK_SUMB = 0.d0
    do K = 0, 5
        do J = 0, 6
            GIJK_SUMS = 0.D0
            do I = 2, 7

                GIJK_SUMS = GIJK_SUMS + gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** I)

            end do
            GIJK_SUMB = GIJK_SUMB + ((GIJK_SUMS) * (y ** J) * (z** K))
        end do
    end do
    log_x = 0.d0
    log_x = log(x)
    g_saline = 0.d0
    g_saline = (gl%sea%GIJK_list(2, 1, 1) + gl%sea%GIJK_list(2, 2, 1) * y) * (x **2.d0) * (log_x) + GIJK_SUMB

    !MB The "check_zerosalinity" function is relevant in case the Salinity (so also x) are equal "0" (to be precize "< 1.d-14").
    !MB If it is 0 we become, thus of the log(x), a "NaN" value from FORTRAN. In this case this function sets the saline function equal to "0.d0"
    !MB Otherwise it gives the write copy of the saline function back (check the "check_zerosalinity" function out for better understanding).
    g_saline = check_zerosalinity (gl, g_saline)

    !g_saline = g_saline / molality_sea(gl)

    !call convert_single_prop(gl, 2, g_saline, 9999)

    End function g_saline
    !_---------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------
    DOUBLE PRECISION FUNCTION v_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (dG/dp)S,T
    !Function for calculating the specific Volume of the Saline Part of Seawater
    !   in m^3 / kg
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_IJK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    double precision, dimension(8,8,7) :: GIJK_list

    call set_coef_Seawater (gl)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_IJK = 0.d0
    do K = 1, 6, 1
        do J = 0, 7, 1
            do I = 2, 7, 1

                GIJK_SUM_IJK = GIJK_SUM_IJK + gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** (I)) * (y ** J) * K * (z ** (K - 1))

            end do
        end do
    end do


    v_saline = (1 / p_r) * GIJK_SUM_IJK

    ! v_saline = v_saline* 3.14038218d-2
    !call convert_single_prop(gl, 2, v_saline, 9999)

    End function v_saline

    DOUBLE PRECISION FUNCTION DGDT_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (dG/dt)S,p
    !Function for calculating the 1st derivative of Gibbs-Energy wrt Temperature of the Saline Part of Seawater
    !   in J / kg
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision :: GIJK_SUM_I, GIJK_SUM_JK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_I = 0.d0
    GIJK_SUM_JK = 0.d0
    do K = 0, 6, 1
        do J = 1, 7, 1
            GIJK_SUM_I = 0.D0
            do I = 2, 7, 1

                GIJK_SUM_I = GIJK_SUM_I + gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** I)

            end do
            GIJK_SUM_JK = GIJK_SUM_JK + ((GIJK_SUM_I + (gl%sea%GIJK_list(2, J + 1, K + 1) * (x ** 2) * (log(x)))) * (y ** (J - 1)) * (z** K) * J)
        end do
    end do

    DGDT_saline = (1 / t_r) * GIJK_SUM_JK

    !MB The "check_zerosalinity" function is relevant in case the Salinity (so also x) are equal "0" (to be precize "< 1.d-14").
    !MB If it is 0 we become, thus of the log(x), a "NaN" value from FORTRAN. In this case this function sets the saline function equal to "0.d0"
    !MB Otherwise it gives the write copy of the saline function back (check the "check_zerosalinity" function out for better understanding).
    DGDT_saline = check_zerosalinity (gl, DGDT_saline)

    !call convert_single_prop(gl, 2, Dgdt_saline, 9999)

    End function DGDT_saline



    DOUBLE PRECISION FUNCTION D2GDT2_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (d2G/dt2)S,p
    !Function for calculating the 2st derivative of Gibbs-Energy wrt Temperature2 of the Saline Part of Seawater
    !   in J / (kg K^2)
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_IJK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_IJK = 0.d0
    do K = 0, 6, 1
        do J = 2, 7, 1
            do I = 2, 7, 1

                GIJK_SUM_IJK = GIJK_SUM_IJK + (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** I) * (y ** (J - 2)) * (z ** K) * J * (J - 1))

            end do
        end do
    end do


    D2GDT2_saline = (1 / (t_r ** 2)) * GIJK_SUM_IJK

    !call convert_single_prop(gl, 2, D2gdt2_saline, 9999)

    End function D2GDT2_saline



    DOUBLE PRECISION FUNCTION D2GDP2_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (d2G/dp2)S,p
    !Function for calculating the 2st derivative of Gibbs-Energy wrt Pressure2 of the Saline Part of Seawater
    !   in m^3 / (kg Pa)
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_IJK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_IJK = 0.d0
    do K = 2, 6, 1
        do J = 0, 7, 1
            do I = 2, 7, 1

                GIJK_SUM_IJK = GIJK_SUM_IJK + (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** I) * (y ** J) * (z ** (K - 2)) * K * (K - 1))

            end do
        end do
    end do


    D2GDP2_saline = (1 / (p_r ** 2)) * GIJK_SUM_IJK

    !call convert_single_prop(gl, 2, d2gdp2_saline, 9999)

    End function D2GDP2_saline
    !--------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION FUNCTION D2GDTDP_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (d2G/dtdp)S,p
    !Function for calculating the 2st derivative of Gibbs-Energy wrt Temperature and Pressure of the Saline Part of Seawater
    !   in m^3 / (kg K)
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_IJK
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_IJK = 0.d0

    do K = 1, 6
        do J = 1, 7
            do I = 2, 7

                GIJK_SUM_IJK = GIJK_SUM_IJK + (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** I) * (y ** (J - 1)) * (z ** (K - 1)) * J * K)

            end do
        end do
    end do


    D2GDTDP_saline = (1 / (t_r * p_r)) * GIJK_SUM_IJK

    !call convert_single_prop(gl, 2, d2gdtdp_saline, 9999)

    End function D2GDTDP_saline
    !_------------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION FUNCTION DGDS_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (dG/ds)S,p
    !Function for calculating the 1st derivative of Gibbs-Energy wrt Salinity of the Saline Part of Seawater
    !   in m^3 / (kg K)
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_JK, GIJK_SUM_I
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_JK = 0.d0
    GIJK_SUM_I = 0.d0
    do K = 0, 5
        do J = 0, 6
            GIJK_SUM_I = 0.D0
            do I = 2, 7

                GIJK_SUM_I = GIJK_SUM_I + I * (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** (I - 2.d0)))

            end do
            GIJK_SUM_JK = GIJK_SUM_JK + ((GIJK_SUM_I + (gl%sea%GIJK_list(2, J + 1, K + 1) * (2 * log(x) + 1))) * (y ** J) * (z** K))
        end do
    end do


    DGDS_saline = (1 / (2 * sal_r)) * GIJK_SUM_JK

    !MB The "check_zerosalinity" function is relevant in case the Salinity (so also x) are equal "0" (to be precize "< 1.d-14").
    !MB If it is 0 we become, thus of the log(x), a "NaN" value from FORTRAN. In this case this function sets the saline function equal to "0.d0"
    !MB Otherwise it gives the write copy of the saline function back (check the "check_zerosalinity" function out for better understanding).
    DGDS_saline = check_zerosalinity (gl, DGDS_saline)

    !call convert_single_prop(gl, 2, dgds_saline, 9999)

    End function DGDS_saline
    !***************************************************************************


    !_------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION FUNCTION D2GDS2_saline(gl, t, p)
    !----------------------------------------------------------------------------------
    !   (d2G/ds2)S,p
    !Function for calculating the 2st derivative of Gibbs-Energy wrt Salinity of the Saline Part of Seawater
    !   in m^3 / (kg K)
    !   Input Parameters:
    !
    !   SAL   -   Salinity in kg / kg
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !----------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    double precision :: GIJK_SUM_JK, GIJK_SUM_I
    double precision ::  sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea
    integer :: I, J, K
    !double precision, dimension(8,8,7) :: GIJK_list

    !call setcoeff_G (gl, GIJK_list)
    call setvariables (gl, sal, sal_n, sal_r, t, t_r, p, p_r, x, y, z, wm_sea)

    GIJK_SUM_JK = 0.d0
    GIJK_SUM_I = 0.d0
    do K = 0, 5
        do J = 0, 6
            GIJK_SUM_I = 0.D0
            do I = 2, 7

                ! GIJK_SUM_I = GIJK_SUM_I + I * (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** (I - 2.d0))) 1rst deriv

                GIJK_SUM_I = GIJK_SUM_I + I * ((0.5d0 - 0.5*(i))) * (1.d0 / sal_r)  * (gl%sea%GIJK_list(I + 1, J + 1, K + 1) * (x ** ( -1.5d0 + 0.5d0 * i)))

            end do
            GIJK_SUM_JK = GIJK_SUM_JK + ((GIJK_SUM_I + (gl%sea%GIJK_list(2, J + 1, K + 1) * (2.d0 * (1.d0/(2.d0 * (sal + (2.d0)**0.5d0 * (sal)**0.5d0))) ) ))) * (y ** J) * (z** K)
        end do
    end do


    D2GDS2_saline = (1 / (2 * sal_r)) * GIJK_SUM_JK

    !MB The "check_zerosalinity" function is relevant in case the Salinity (so also x) are equal "0" (to be precize "< 1.d-14").
    !MB If it is 0 we become, thus of the log(x), a "NaN" value from FORTRAN. In this case this function sets the saline function equal to "0.d0"
    !MB Otherwise it gives the write copy of the saline function back (check the "check_zerosalinity" function out for better understanding).
    D2GDS2_saline = check_zerosalinity (gl, D2GDS2_saline)

    !call convert_single_prop(gl, 2, dgds_saline, 9999)

    End function D2GDS2_saline
    !***************************************************************************


    !***************************************************************************
    DOUBLE PRECISION FUNCTION D2GDP2_water(gl, T, D, nrsubst)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D
    integer:: nrsubst

    D2GDP2_water = (- ((compt_calc (gl, T, D, 1)) / (d*gl%wm(1)))) / 1.d6

    END FUNCTION D2GDP2_water
    !************************************************************************

    !***********************************************************************
    DOUBLE PRECISION FUNCTION D2GDP2_sea(gl, T, P, D)


    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, p
    integer:: nrsubst


    D2GDP2_sea = D2GDP2_water(gl, T, D, 1)

    D2GDP2_sea = D2GDP2_sea + d2gdp2_saline(gl, t, p)

    !call convert_single_prop(gl, 2, D2GDP2_sea, 9999)

    D2GDP2_sea = D2GDP2_sea !* 1.d6

    END FUNCTION D2GDP2_sea
    !*************************************************************************

    !***************************************************************************
    DOUBLE PRECISION FUNCTION D2GDTDP_water(gl, T, D, nrsubst)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water, D2FDTDD_water
    integer:: nrsubst
    logical :: seawater

    !gl%seacalc = .false.

    D2FDTDD_water = (VOLEXP_CALC(gl,T,D, 1)) / (((d*gl%wm(1)) ** 2.d0) * compt_calc(gl,T,D, 1)*1.d-6)

    D2GDTDP_water = (d*gl%wm(1)) * (compt_calc(gl,T,D, 1)*1.d-6) * (D2FDTDD_water)

    !!gl%seacalc = .true.

    END FUNCTION D2GDTDP_water
    !********************************************************************************************

    !********************************************************************************************
    DOUBLE PRECISION FUNCTION D2GDTDP_sea(gl, T, p,D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, p
    integer:: nrsubst
    logical :: seawater

    D2GDTDP_sea = D2GDTDP_water(gl, T, d, 1) + d2gdtdp_saline(gl, t, p)

    !call convert_single_prop(gl, 2, D2GDTDP_sea, 9999)

    END FUNCTION D2GDTDP_sea
    !***********************************************************************

    !********************************************************************
    DOUBLE PRECISION FUNCTION D2GDT2_water(gl, T, D, nrsubst)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water
    integer:: nrsubst
    logical :: seawater

    !gl%seacalc = .false.

    D2GDT2_water = - (cp_CALC(gl,T,D, 1) / gl%wm(1)/ T)

    ! gl%seacalc = .true.

    END FUNCTION D2GDT2_water
    !************************************************************

    !**************************************************************
    DOUBLE PRECISION FUNCTION D2GDT2_sea(gl, T,p, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.

    D2GDT2_sea = D2GDT2_water(gl, T, D, 1)

    ! gl%seacalc = .true.

    D2GDT2_sea = D2GDT2_sea + D2GDT2_saline(gl, t, p)

    !call convert_single_prop(gl, 2, D2GDT2_sea, 9999)

    END FUNCTION D2GDT2_sea

    DOUBLE PRECISION FUNCTION DGDT_sea(gl, T,p, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.

    DGDT_sea = (- s_calc(gl, T, D, 1) / gl%wm(1) )

    ! gl%seacalc = .true.

    DGDT_sea = DGDT_sea + dgdt_saline(gl, t, p)

    !call convert_single_prop(gl, 2, DGDT_sea, 9999)

    END FUNCTION DGDT_sea
    !***********************************************************************

    !**********************************************************************
    DOUBLE PRECISION FUNCTION DGDP_sea(gl, T,p, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water, dw, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.

    dw = d * gl%wm(1)!rhomix_calc(gl, t, p, 0.d0, 0, 1)*gl%wm(1)

    dw = rhomix_calc(gl, t, p, 0.d0, 1, 1)*gl%wm(1)

    !gl%seacalc = .true.

    DGDP_sea = 1.d0/(dw) + v_saline(gl, t, p)

    !call convert_single_prop(gl, 2, DGDP_sea, 9999)

    END FUNCTION DGDP_sea
    !***********************************************************************

    !*********************************************************************
    double Precision Function check_zerosalinity (gl, saline_sum)



    implicit none

    type(type_gl) :: gl
    double precision :: saline_sum

    if (dabs(gl%sea%salinity) <= 1.d-14) then
        check_zerosalinity = 0.d0
    else
        check_zerosalinity = saline_sum
    endif

    end function check_zerosalinity


    !Functions for saturation props of seawater (s, t, p)

    double precision function diff_boilsal(gl,salb,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_boil, d, g_vap, t, p, salb, dvap
    type(type_additional_parameters) :: parameters

    !sald = parameters(3)

    gl%sea%salinity = salb

    t = parameters%a_p(1)
    p = parameters%a_p(2)
    dvap = rhomix_calc(gl, t, p, 0.d0, 2, 1)
    d = parameters%a_p(3)
    parameters%a_p(4) = dvap

    g_boil =chempot_sea_Calc(gl, t, p, d)

    !call convert_single_prop(gl, 2, g_boil, 9999)

    gl%seacalc = .false.
    g_vap = g_calc(gl, t, dvap, 1)
    !gl%seacalc = .true.

    diff_boilsal = g_boil - g_vap

    !if((dabs(diff_boilsal)) .le. (1.d-14)) then
    !return
    !end if

    end function
    !_------------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    double precision function boil_sal(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_boil, salmin, salmax, delta_allowed, t, p, salsave, d, salb
    type(type_additional_parameters) :: parameters
    integer :: maxiter, errorflag, iterations

    salsave = gl%sea%salinity
    salb = gl%sea%salinity
    !parameters = 0.d0
    Delta_allowed = 1.d-13
    maxiter = 100
    parameters%a_p (1) = T
    parameters%a_p (2) = p
    parameters%a_p (3) = d

    call Regula_Falsi(gl,diff_boilsal, salb, 1.d-10, 0.2d0, Delta_allowed, 0.d0, 0.2d0, &
        &                       maxiter, Iterations, Errorflag, parameters)

    if(errorflag .eq. 0)then
        boil_sal = salb
    else
        boil_sal = errorflag
    end if

    gl%sea%salinity = salsave

    end function
    !_------------------------------------------------------------------------------------------------------------


    !****************************
    double precision function diff_freezesal(gl,sald,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, d, g_ice, t, p, sald, gt
    type(type_additional_parameters) :: parameters


    gl%sea%salinity = sald

    t = parameters%a_p(1)
    p = parameters%a_p(2)
    d = parameters%a_p(3)

    g_freeze = chempot_sea_calc(gl, t, p, d)
    gt = g_sea_calc(gl, t, p, d)

    !call convert_single_prop(gl, 2, g_freeze, 9999)

    gl%seacalc = .false.
    g_ice = g_waterIce(gl, t, p)
    !gl%seacalc = .true.

    diff_freezesal = g_freeze - g_ice
    gt = gt - g_ice

    end function
    !_------------------------------------------------------------------------------------------------------------

    !_------------------------------------------------------------------------------------------------------------
    double precision function freeze_sal(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, salf
    type(type_additional_parameters) :: parameters
    integer :: iterations, maxiter, errorflag

    salsave = gl%sea%salinity

    !parameters = 0.d0
    Delta_allowed = 1.d-13
    maxiter = 100
    parameters%a_p(1) = T
    parameters%a_p(2) = p
    parameters%a_p(3) = d

    call Regula_Falsi(gl,diff_freezesal, salf, 0.d0, 0.2d0, Delta_allowed, 0.d0, 0.2d0, &
        &                       maxiter, Iterations, Errorflag,parameters)

    if(errorflag .eq. 0)then
        freeze_sal = salf
    else
        freeze_sal = errorflag
    end if

    gl%sea%salinity = salsave

    end function
    !_------------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    double precision function diff_freezetemp(gl,tf,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, d, g_ice, t, p, tf, dliq, test
    type(type_additional_parameters) :: parameters


    !gl%sea%salinity = sald

    !t = parameters(1)
    p = parameters%a_p(2)
    d = parameters%a_p(3)
    dliq = rhomix_calc(gl, tf, p, 0.d0, 1, 1)!parameters(4)

    g_freeze = chempot_sea_calc(gl, tf, p, dliq)

    !call convert_single_prop(gl, 2, g_freeze, 9999)

    gl%seacalc = .false.
    g_ice = g_waterIce(gl, tf, p)
    !gl%seacalc = .true.

    diff_freezetemp = g_freeze - g_ice

    test = freeze_press(gl, tf, gl%sea%seap, d)


    end function
    !_------------------------------------------------------------------------------------------------------------



    !_------------------------------------------------------------------------------------------------------------
    double precision function freeze_temp(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tf
    type(type_additional_parameters) :: parameters
    integer :: iterations, maxiter, errorflag

    !salsave = gl%sea%salinity

    !parameters = 0.d0
    Delta_allowed = 1.d-13
    maxiter = 100
    ! parameters(1) = T
    parameters%a_p(2) = p
    parameters%a_p(3) = d
    parameters%a_p(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    tf = t

    call Regula_Falsi(gl,diff_freezetemp, tf, 250.d0, 300.d0, Delta_allowed, 200.d0, 500.d0, &
        &                       maxiter, Iterations, Errorflag,parameters)

    if(errorflag .eq. 0)then
        freeze_temp = tf
    else
        freeze_temp = errorflag
    end if

    ! gl%sea%salinity = salsave

    end function
    !_------------------------------------------------------------------------------------------------------------

    !_------------------------------------------------------------------------------------------------------------
    double precision function diff_boiltemp(gl,tb,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, d, g_ice, t, p, tb, dliq, g_boil, dvap, g_vap
    type(type_additional_parameters) :: parameters

    !gl%sea%salinity = sald

    !t = parameters(1)
    p = parameters%a_p(2)
    !d = parameters(3)
    dliq = rhomix_calc(gl, tb, p, 0.d0, 1, 1)!parameters(4)
    dvap = rhomix_calc(gl, tb, p, 0.d0, 2, 1)!vapor density

    g_boil = chempot_sea_calc(gl, tb, p, dliq)

    !call convert_single_prop(gl, 2, g_freeze, 9999)

    gl%seacalc = .false.
    g_vap = g_calc(gl, tb, dvap, 1)
    !gl%seacalc = .true.

    diff_boiltemp = g_boil - g_vap

    end function
    !_------------------------------------------------------------------------------------------------------------



    !_------------------------------------------------------------------------------------------------------------
    double precision function boil_temp(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tb
    type(type_additional_parameters) :: parameters

    integer :: iterations, maxiter, errorflag

    !salsave = gl%sea%salinity

    !parameters = 0.d0
    Delta_allowed = 1.d-12
    maxiter = 1000
    ! parameters(1) = T
    parameters%a_p(2) = p
    !parameters(3) = d
    parameters%a_p(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    tb = t

    call Regula_Falsi(gl,diff_boiltemp, tb, 250.d0, 500.d0, Delta_allowed, 200.d0, 800.d0, &
        &                       maxiter, Iterations, Errorflag,parameters)

    if(errorflag .eq. 0)then
        boil_temp = tb
    else
        boil_temp = tb
    end if

    ! gl%sea%salinity = salsave

    end function
    !_------------------------------------------------------------------------------------------------------------
    !
    !
    !
    !_------------------------------------------------------------------------------------------------------------
    double precision function diff_freezepress(gl,pf,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, d, g_ice, t, pf, dliq
    type(type_additional_parameters) :: parameters

    !gl%sea%salinity = sald

    t = parameters%a_p(1)
    !p = parameters(2)
    d = parameters%a_p(3)
    dliq = rhomix_calc(gl, t, pf, 0.d0, 1, 1)!parameters(4)

    g_freeze = chempot_sea_calc(gl, t, pf, dliq)

    !call convert_single_prop(gl, 2, g_freeze, 9999)

    gl%seacalc = .false.
    g_ice = g_waterIce(gl, t, pf)
    !gl%seacalc = .true.

    diff_freezepress = g_freeze - g_ice

    end function
    !***************************************************************************

    !**************************************************************************
    double precision function freeze_press(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tf, pf, pmin, pmax
    type(type_additional_parameters) :: parameters
    integer :: iterations, maxiter, errorflag
    logical :: seacalcorig

    !salsave = gl%sea%salinity

    seacalcorig = gl%seacalc

    !parameters = 0.d0
    Delta_allowed = 1.d-11
    maxiter = 100
    parameters%a_p(1) = T
    !parameters(2) = p
    parameters%a_p(3) = d
    parameters%a_p(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    pf = p
    pmin = 0.1d0* pmelt_eq(gl,T, 1)
    pmax = 10.d0*pmelt_eq(gl, t, 1)

    call Regula_Falsi(gl,diff_freezepress, pf, pmin, pmax, Delta_allowed,0.1d0*pmin, 10.d0*pmax, &
        &                       maxiter, Iterations, Errorflag,parameters)

    if((errorflag .eq. 0) .or. (errorflag .eq. 3))then
        freeze_press = pf
        !gl%sea%seap = pf
    else
        freeze_press = errorflag
    end if

    ! gl%sea%salinity = salsave

    gl%seacalc = seacalcorig

    end function
    !************************************************************************

    !************************************************************************
    double precision function diff_boilpress(gl,pb,parameters)



    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, d, g_vap, t, p, tb, dliq, pb, g_boil, dvap
    type(type_additional_parameters) :: parameters

    !gl%sea%salinity = sald

    t = parameters%a_p(1)
    !p = parameters(2)
    d = parameters%a_p(3)
    dliq = rhomix_calc(gl, t, pb, 0.d0, 1, 1)!parameters(4)
    dvap = rhomix_calc(gl, t, pb, 0.5d0, 2, 1)!vapor density
    !if(dvap .ge. 17000.d0) then
    !    continue
    !end if!dvap = 17000.d0
    g_boil = chempot_sea_calc(gl, t, pb, dliq)

    !call convert_single_prop(gl, 2, g_freeze, 9999)

    gl%seacalc = .false.
    g_vap = g_calc(gl, t, dvap, 1)
    !gl%seacalc = .true.

    diff_boilpress = g_boil - g_vap

    end function
    !_------------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    double precision function boil_press(gl, t, p, d)


    implicit none

    type(type_gl) :: gl
    double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tb, pb, p_max
    type(type_additional_parameters) :: parameters
    integer :: iterations, maxiter, errorflag

    !salsave = gl%sea%salinity

    !parameters = 0.d0
    Delta_allowed = 1.d-11
    maxiter = 100
    parameters%a_p(1) = T
    !parameters(2) = p
    parameters%a_p(3) = d
    parameters%a_p(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    pb = p
    p_max = vp_eq(gl,T,1)

    ! call Regula_Falsi(gl,diff_boilpress, pb, 1.d-10, 0.1d-2, Delta_allowed, 00.d0, 50.d0, &
    !     &                       maxiter, Iterations, Errorflag,Parameters)

    call Regula_Falsi(gl,diff_boilpress, pb, 2.d-4, 1.1d0*p_max, Delta_allowed, 0.1d-5, 3.d0*p_max, &
        &                       maxiter, Iterations, Errorflag,parameters)

    if((errorflag .eq. 0) .or. (errorflag .eq.3))then
        boil_press = pb
    else
        boil_press = errorflag
    end if

    ! gl%sea%salinity = salsave

    end function
    !_------------------------------------------------------------------------------------------------------------



    !!********************************************************************
    !!tripe point pressure
    !now realised with solid_3p
    !
    !subroutine sea_trip(gl, t, p, d, ttrip, ptrip)


    !implicit none
    !
    !type(type_gl) :: gl
    !double precision :: diff, delta_allowed, t, p,  d, tb, ttrip, ptrip
    !double precision :: parameters(65)
    !integer :: iterations, maxiter, errorflag
    !
    !
    !parameters = 0.d0
    !Delta_allowed = 1.d-11
    !maxiter = 100
    !parameters(1) = T
    !parameters(2) = p
    !parameters(3) = d
    !parameters(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    !ttrip = t
    !
    !call Regula_Falsi(gl,diff_triptemp, ttrip, 262.d0, 273.17d0, Delta_allowed, 260.d0, 273.2d0, &
    !    &                       maxiter, Iterations, Errorflag,Parameters)
    !
    !ptrip = parameters(2)
    !ttrip = parameters(1)
    !gl%seacalc = .false.
    !
    !end subroutine sea_trip
    !_------------------------------------------------------------------------------------------------------------


    !!_------------------------------------------------------------------------------------------------------------
    !double precision function diff_trippress(gl,pb,parameters)
    !


    !implicit none
    !
    !type(type_gl) :: gl
    !double precision :: diff, g_freeze, d, g_vap,diff1, diff2, g_ice, g_waterice, t, p,  chempot_sea_calc, tb,  rhomix_calc, dliq, pb, g_boil, dvap, trip_temp
    !double precision :: parameters(65)
    !
    !!gl%sea%salinity = sald
    !
    !t = parameters(1)
    !
    !dliq = rhomix_calc(gl, t, pb, 0.d0, 1, 1)!parameters(4)
    !dvap = rhomix_calc(gl, t, pb, 0.d0, 2, 1)!vapor density
    !
    !
    !gl%seacalc = .false.
    !g_vap = g_calc(gl, t, dvap, 1)
    !g_ice = g_waterice(gl, t, pb)
    !gl%seacalc = .false.
    !
    !
    !diff_trippress = g_vap - g_ice
    !
    !if (abs(diff_trippress) .le. 1.d-10) then
    !    continue
    !end if
    !
    !end function diff_trippress
    !**********************************************************

    !!**********************
    !!triple point temp
    !double precision function trip_temp(gl, t, p, d)


    !implicit none
    !
    !type(type_gl) :: gl
    !double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tb
    !double precision :: parameters(65)
    !integer :: iterations, maxiter, errorflag
    !
    !!salsave = gl%sea%salinity
    !
    !parameters = 0.d0
    !Delta_allowed = 1.d-11
    !maxiter = 100
    !! parameters(1) = T
    !parameters(2) = p
    !parameters(3) = d
    !parameters(4) = rhomix_calc(gl, t, p, 0.d0, 1, 1)
    !tb = t
    !
    !call Regula_Falsi(gl,diff_triptemp, tb, 255.d0, 275.d0, Delta_allowed, 230.d0, 280.d0, &
    !    &                       maxiter, Iterations, Errorflag,Parameters)
    !
    !if(errorflag .eq. 0)then
    !    trip_temp = tb
    !else
    !    trip_temp = t * 0.98d0
    !end if
    !
    !! gl%sea%salinity = salsave
    !
    !end function
    !_------------------------------------------------------------------------------------------------------------


    !_------------------------------------------------------------------------------------------------------------
    !double precision function diff_triptemp(gl,tb,parameters)
    !


    !implicit none
    !
    !type(type_gl) :: gl
    !double precision :: diff, g_freeze, d, g_vap,diff1, diff2,  g_waterice, t, p,  chempot_sea_calc, tb,  rhomix_calc, dliq, pb, g_boil, dvap, chempot_sea, g_ice, ptrip
    !double precision :: parameters(65)
    !
    !!gl%sea%salinity = sald
    !
    !!t = parameters(1)
    !p = parameters(2)
    !d = parameters(3)
    !dliq = rhomix_calc(gl, tb, p, 0.d0, 1, 1)!parameters(4)
    !dvap = rhomix_calc(gl, tb, p, 0.d0, 2, 1)!vapor density
    !
    !chempot_sea = chempot_sea_calc(gl, tb, p, dliq)
    !!g_boil = g_calc(gl, t, dliq, 1)
    !!call convert_single_prop(gl, 2, g_freeze, 9999)
    !
    !gl%seacalc = .false.
    !call g_sea_trip_diff(gl, tb, parameters, g_ice, ptrip)
    !
    !!gl%seacalc = .true.
    !dliq = rhomix_calc(gl, tb, ptrip, 0.d0, 1, 1)
    !chempot_sea = chempot_sea_calc(gl, tb, ptrip, dliq)
    !
    !diff_triptemp = chempot_sea - g_ice
    !!diff2 = g_boil - g_ice
    !
    !!    diff_triptemp = diff1! - diff2
    !
    !
    !end function diff_triptemp
    !_------------------------------------------------------------------------------------------------------------



    !!_------------------------------------------------------------------------------------------------------------
    !subroutine g_sea_trip_diff(gl, tb, parameters, g_ice, ptrip)


    !implicit none
    !
    !type(type_gl) :: gl
    !double precision :: diff, g_freeze, salmin, salmax, delta_allowed, t, p, salsave, d, tb, g_waterice, g_ice, ptrip
    !double precision :: parameters(65)
    !integer :: iterations, maxiter, errorflag
    !
    !!salsave = gl%sea%salinity
    !
    !
    !Delta_allowed = 1.d-11
    !maxiter = 100
    !parameters(1) = tb
    !ptrip = parameters(2)
    !!if((ptrip .le.0.3d-3) .or. (ptrip .ge. 0.
    !
    !call Regula_Falsi(gl,diff_trippress, ptrip, 0.3d-3, 0.62d-3, Delta_allowed, 0.25d-4, 0.65d-2, maxiter, Iterations, Errorflag,Parameters)
    !
    !if((errorflag .eq. 0) .or. (errorflag .eq. 3))then
    !    g_ice = g_waterice(gl, tb, ptrip)
    !    parameters(2) = ptrip
    !else
    !    tb = t * 0.98d0
    !end if
    !
    !! gl%sea%salinity = salsave
    !
    !end subroutine g_sea_trip_diff
    !******************************************************************************
    !******************************************************************************


    !****************************************************************************
    double precision function vapfrac_sea(gl, t, p, d)



    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, sb, sum

    sb = boil_sal(gl, t, p, d)

    vapfrac_sea = 1.d0 - (gl%sea%salinity / sb)             !mass fraction of water vapor

    sum = 0.d0

    sum = vapfrac_sea / gl%wm(1) +  ( (1.d0 - vapfrac_sea) / molar_sea(gl))

    vapfrac_sea = (vapfrac_sea / gl%wm(1)) / sum        !molar fraction of water vapor

    end function vapfrac_sea
    !*******************************************************************************


    !****************************************************************************
    double precision function icefrac_sea(gl, t, p, d)



    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, sf, sum

    sf = freeze_sal(gl, t, p, d)

    icefrac_sea = 1.d0 - (gl%sea%salinity / sf)             !mass fraction of water vapor

    sum = 0.d0

    sum = icefrac_sea / gl%wm(1) +  ( (1.d0 - icefrac_sea) / molar_sea(gl))

    icefrac_sea = (icefrac_sea / gl%wm(1)) / sum        !molar fraction of water vapor

    end function icefrac_sea
    !*******************************************************************************
    !********************************************************************************


    !Special routines for seawater mixtures in liquid phase as ws needs to be calculated over gibbs energy, simple addition not possible!! (BS 09/2018)
    !
    double precision function dgdp_mix(gl, t, p, d)


    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water, dw, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.

    dw = d!rhomix_calc(gl, t, p, 0.d0, 1, 0)

    !gl%seacalc = .true.

    DGDP_mix = 1.d0/(dw)

    !call convert_single_prop(gl, 2, DGDP_mix, 2)

    END FUNCTION DGDP_mix
    !**********************************************************************

    !**********************************************************************
    DOUBLE PRECISION FUNCTION D2GDT2_mix(gl, T, p,D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.

    D2GDT2_mix = - (cp_CALC(gl,T,D, 0)) / T !/ gl%wm(1)/ T)

    gl%seacalc = .true.

    END FUNCTION D2GDT2_mix
    !************************************************************

    !*************************************************************
    DOUBLE PRECISION FUNCTION D2GDTDP_mix(gl, T,p, D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, D2GDP2_water, D2FDTDD_mix, wmmix, p
    integer:: nrsubst
    logical :: seawater

    gl%seacalc = .false.
    call wm_mix_calc(gl, wmmix)
    !d = rhomix_calc(gl, t, p, 0.d0, 1, 0)
    D2FDTDD_mix = (VOLEXP_CALC(gl,T,D, 0)) / (((d*wmmix) ** 2.d0) * compt_calc(gl,T,D, 0)*1.d-6)

    D2GDTDP_mix = (d*wmmix) * (compt_calc(gl,T,D, 0)*1.d-6) * (D2FDTDD_mix)

    !gl%seacalc = .true.

    END FUNCTION D2GDTDP_mix
    !************************************************************************

    !**********************************************************************
    DOUBLE PRECISION FUNCTION D2GDP2_mix(gl, T, p,D)



    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION :: T, D, wmmix,p
    integer:: nrsubst

    call wm_mix_calc(gl, wmmix)
    D2GDP2_mix = (- ((compt_calc (gl, T, D, 0)) / (d*wmmix))) / 1.d6

    END FUNCTION D2GDP2_mix
    !*****************************************************************************


    !*****************************************************************************
    !function for molar density of seawater and other brines
    !Benedikt 2018
    !double precision function brine_dens(gl, t, p, d)
    !
    !
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !double precision, intent(in) :: t, p, d
    !double precision ::  v_helm, v_salt, v_brine, x_salt, x_helm, wm_mix_helm, wm_water, wm_brine, d_helm
    !!
    !
    !d_helm = d
    !
    !v_salt = v_saline(gl, T, p)! * gl%sea%wm_sea!gl%sea%wm_sea
    !
    !call wm_mix_calc(gl, wm_mix_helm)
    !
    !d_helm = d_helm * wm_mix_helm
    !
    !v_helm = 1.d0 / d_helm
    !
    !x_helm = 1.d0
    !
    !x_salt = gl%molfractions(1)
    !
    !v_brine = v_helm * x_helm + v_salt * x_salt
    !
    !brine_dens = 1.d0 / v_brine
    !
    !wm_water = gl%wm(1)
    !
    !gl%wm(1) = gl%sea%wm_sea
    !
    !call wm_mix_calc(gl, wm_brine)
    !
    !brine_dens = brine_dens/wm_brine
    !
    !gl%wm(1) = wm_Water
    !
    !
    !end function brine_dens
    !!********************************************************************************
    !!********************************************************************************

    subroutine seawater_remap(gl, fluids_in, moles_in, x_phase_in, chempot_in, phasetype)


    implicit none

    type(type_gl) :: gl

    character(30), dimension(30) :: fluids_in, fluids_out
    double precision, dimension(30) :: moles_in, moles_out
    double precision, dimension(30,5) :: x_phase_in, x_phase_out, prop_in, prop_out, chempot_in, chempot_out
    integer :: i, seaphases
    integer, dimension(5) :: phasetype

    fluids_out = ' '
    moles_out = 0.d0
    x_phase_out = 0.d0
    ! prop_out = 0.d0

    do i = 1,gl%salpos-1

        fluids_out(i) = fluids_in(i)
        moles_out(i) = moles_in(i)
        x_phase_out(i,:) = x_phase_in(i,:)
        chempot_out(i,:) = chempot_in(i,:)
        ! prop_out(i,:) = prop_in(i,:)

    end do

    fluids_out(gl%salpos-1) = 'seawater'
    fluids_out(gl%salpos) = 'salinity'
    moles_out(gl%salpos) = gl%sea%salinity
    if(phasetype(1)==2 .or. phasetype(2)==2) then
    x_phase_out(gl%salpos,2) = gl%sea%salinity
    end if
    if(phasetype(1)==3 .or. phasetype(2)==3 .or. phasetype(3)==3) then
    x_phase_out(gl%salpos,3) = gl%sea%salinity
    end if
    chempot_out(gl%salpos,:) = 0.d0
    ! prop_out(gl%salpos,:) = 0.d0


    do i = gl%salpos+1,gl%ncomp+1

        fluids_out(i) = fluids_in(i-1)
        moles_out(i) = moles_in(i-1)
        x_phase_out(i,:) = x_phase_in(i-1,:)
        chempot_out(i,:) = chempot_in(i-1,:)
        !        prop_out(i,:) = prop_in(i-1,:)

    end do

    fluids_in = fluids_out
    moles_in = moles_out
    x_phase_in = x_phase_out
    chempot_in = chempot_out
    !  prop_in = prop_out


    end subroutine !seawater_remap
    !********************************************************************
    !********************************************************************


   



    end module seawater_module
