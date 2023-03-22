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
    ! module for file gibbsderivs.f90
    module gibbsderivs_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use rhomix_pt_module
    use seawater_module
    use reduced_parameters_calc_module
    use fnrderivs_module
    use fniderivs_module
    use electrolytes
    use phasedet_sol_module
    use phasedet_sol_pure_module
    !use interface_support_module

    contains

    !**************************************************************************
    !Describption of modelflags
    ! 1 : seawater
    ! 2 : Pitzer-electrolyte model




    !****************************************************************************************
    double precision function A10_sea(gl, t, d, p)
    !_-------------------------------------------------------------------------------------
    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018

    implicit none

    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdt, dgdp, d_w, d, t_water, p_water, d_water, d_mixsave
    logical :: seaflagorig, seasave

    d_mixsave = d
    seasave = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1)
    !gl%gecarrier = seasave

    !seaflagorig = gl%gecarrier

    !gl%gecarrier = .false.

    d_w = d

    R = gl%req(1) / gl%sea%wm_sea
    g = g_Sea_calc(gl,t,p, d_w)
    dgdt =  DGDT_sea(gl, t, p, d_w)
    dgdp = DGDP_sea(gl, t, p, d_w)

    p = gl%sea%seap

    A10_sea = a10_gibbs_trans(gl,t, p,  R, g, dgdt, dgdp) !* gl%molfractions(1)

    gl%gecarrier = seasave

    d = d_mixsave

    end function
    !***********************************************************************************************************


    !****************************************************************************************
    double precision function A10_el(gl, t, d, p)
    !_-------------------------------------------------------------------------------------
    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018

    implicit none

    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdt, dgdp, d_w, d, t_water, p_water, d_water, d_mixsave
    logical :: seaflagorig, seasave
    integer :: pos

    pos = gl%el%solpos
    p = gl%el%press

    d_mixsave = d
    seasave = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos)
    gl%gecarrier = seasave

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    d_w = d

    R = gl%req(pos)
    g = g_brine(gl,t,d_w, p)
    dgdt = - s_brine(gl, t, d_w,p)
    dgdp = dbrine_dp(gl, t, d_w,p)

    A10_el = a10_gibbs_trans(gl,t, p,  R, g, dgdt, dgdp) !* gl%molfractions(1)

    gl%gecarrier = seaflagorig

    d = d_mixsave

    end function
    !***********************************************************************************************************



    !***********************************************************************************************************
    double precision function A01_sea(gl, t, d, p)
    !--------------------------------------------------------------------------------------------
    !function for calculating  A01 from seawater gibbs equations for saline part
    !Benedikt July 2018

    implicit none

    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdp, d_w, p_test, dens_sea, d, salsave, t_water, p_water, d_water, d_ard, d_mixsave
    logical :: seaflagorig, seasave
    integer :: ncomp

    d_mixsave = d
    seasave = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1)
    gl%gecarrier = seasave

    salsave = gl%sea%salinity

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d_w = d

    gl%sea%wm_sea = molar_sea(gl)

    R = gl%req(1) / gl%sea%wm_sea
    g = g_sea_calc(gl, t, p, d_w)
    dgdp = dgdp_sea(gl, t, p, d_w)

    A01_sea = a01_gibbs_trans(gl, t, p, dgdp, R,g) !*gl%molfractions(1)

    gl%gecarrier = seaflagorig

    gl%sea%salinity = salsave

    d = d_mixsave
    end function
    !*******************************************************************************************


    !***********************************************************************************************************
    double precision function A01_el(gl, t, d, p)
    !--------------------------------------------------------------------------------------------
    !function for calculating  A01 from seawater gibbs equations for saline part
    !Benedikt July 2018

    implicit none
    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdp, d_w, p_test, dens_sea, d, salsave, t_water, p_water, d_water, d_ard, d_mixsave
    logical :: seaflagorig, seasave
    integer :: ncomp, pos

    pos = gl%el%solpos

    p = gl%el%press

    d_mixsave = d
    seasave = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos)
    gl%gecarrier = seasave

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d_w = d

    R = gl%req(pos)
    g = g_brine(gl, t,d_w, p)
    dgdp = dbrine_dp(gl, t, d_w, p)

    A01_el = a01_gibbs_trans(gl, t, p, dgdp, R,g) !*gl%molfractions(1)

    gl%gecarrier = seaflagorig
    d = d_mixsave
    end function
    !*******************************************************************************************



    !****************************************************************************************
    double precision function A11_sea(gl, t, d, p)
    !-----------------------------------------------------------------------------------------
    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, R, d2gdt2, d_w, d, a01, a20, a02, d_mixsave
    logical :: seaflagorig, seasave

    d_mixsave = d

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier  

    gl%gecarrier = .false.

    R = gl%req(1)/ gl%sea%wm_sea

    d_w = d!rhomix_calc(gl , t, p, 0.d0, 1, 1)

    d2gdt2 =  D2GDT2_sea(gl, t, p, d_w)

    a01 = A01_sea(gl, t,d, p)

    a20 = A20_sea(gl, t, d, p)

    a02 = A02_sea(gl, t, d, p)

    A11_sea = a11_gibbs_trans(gl, t, R, d2gdt2, a01, a20, a02) 

    gl%gecarrier = seaflagorig

    d = d_mixsave
    end function A11_sea
    !**************************************************************************************************************

    !#################################################################################################
    !#################################################################################################


    !****************************************************************************************
    double precision function A11_el(gl, t, d, p)
    !-----------------------------------------------------------------------------------------
    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018

    implicit none

    type(type_gl) :: gl

    double precision :: t, p, R, d2gdt2, d_w, d, a01, a20, a02, d_mixsave
    logical :: seaflagorig, seasave
    integer :: pos

    pos = gl%el%solpos
    p = gl%el%press
    d_mixsave = d

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    R = gl%req(pos)

    d_w = d
    d2gdt2 =  d2brine_dt2(gl, t, d_w, p)

    a01 = A01_el(gl, t,d, p)

    a20 = A20_el(gl, t, d, p)

    a02 = A02_el(gl, t, d, p)

    A11_el = a11_gibbs_trans(gl, t, R, d2gdt2, a01, a20, a02) 

    gl%gecarrier = seaflagorig

    d = d_mixsave
    end function A11_el
    !**************************************************************************************************************


    !*************************************************************************************************************
    double precision function A20_sea(gl, t, d, p)

    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl

    double precision :: t, p, R, g, d, d_w, d2gdt2, d2gdp2, d2gdtdp, d_mixsave
    logical :: seaflagorig

    d_mixsave = d
    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    R = gl%req(1) / gl%sea%wm_sea

    d_w  = d

    g = g_Sea_calc(gl,t,p, d_w)
    d2gdt2 =  D2GDT2_sea(gl, t, p, d_w)
    d2gdp2 = D2GDP2_sea(gl, t, p, d_w)
    d2gdtdp = D2GDTDP_sea(gl, t, p, d_w)

    A20_sea = a20_gibbs_trans(gl,t, p, R, d2gdt2, d2gdp2, d2gdtdp) 

    gl%gecarrier = seaflagorig

    d = d_mixsave

    end function
    !********************************************************************************************



    !*************************************************************************************************************
    double precision function A20_el(gl, t, d, p)

    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl

    double precision :: t, p, R, g, d, d_w, d2gdt2, d2gdp2, d2gdtdp, d_mixsave
    logical :: seaflagorig
    integer :: pos

    pos = gl%el%solpos
    p = gl%el%press
    d_mixsave = d
    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    R = gl%req(pos)

    d_w  = d

    g = g_brine(gl,t, d_w, p)
    d2gdt2 =  d2brine_dt2(gl, t, d_w, p)
    d2gdp2 = d2_brine_dp2(gl, t, d_w, p)
    d2gdtdp = d2_brine_dtdp(gl, t, d_w, p)

    A20_el = a20_gibbs_trans(gl,t, p, R, d2gdt2, d2gdp2, d2gdtdp) 

    gl%gecarrier = seaflagorig

    d = d_mixsave

    end function
    !********************************************************************************************



    !********************************************************************************************
    double precision function A02_sea(gl, t, d, p)
    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdp, d2gdp2, d2gdt2, d2gdtdp, d_w, d, t_water, p_water, d_water, d_mixsave
    logical :: seaflagorig

    d_mixsave = d

    seaflagorig = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    d_w  = d!rhomix_Calc(gl, t, p, 0.d0, 1, 1)

    !call partial_reac_press(gl, t, p, d, 1, t_water, p_water, d_water)

    R = gl%req(1) / gl%sea%wm_sea!gl%wm(1)
    g = g_sea_calc(gl,t,p, d_w)
    d2gdt2 =  D2GDT2_sea(gl, t, p, d_w)
    d2gdp2 = D2GDP2_sea(gl, t, p, d_w)
    dgdp = dgdp_sea(gl, T, p, d_w)
    d2gdtdp = D2GDTDP_sea(gl, t, p, d_w)

    A02_sea = a02_gibbs_trans(gl, t, p, dgdp, R, d2gdp2, d2gdt2, d2gdtdp, g) 

    gl%gecarrier = seaflagorig

    d = d_mixsave

    end function
    !*************************************************************************************************************


    !********************************************************************************************
    double precision function A02_el(gl, t, d, p)

    !function for calculating  A10 from seawater gibbs equations for saline part
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdp, d2gdp2, d2gdt2, d2gdtdp, d_w, d, t_water, p_water, d_water, d_mixsave
    logical :: seaflagorig
    integer :: pos

    pos = gl%el%solpos
    p = gl%el%press
    d_mixsave = d

    seaflagorig = gl%gecarrier
    gl%gecarrier = .false.
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos)
    gl%gecarrier = seaflagorig

    seaflagorig = gl%gecarrier

    gl%gecarrier = .false.

    d_w  = d!rhomix_Calc(gl, t, p, 0.d0, 1, 1)

    !call partial_reac_press(gl, t, p, d, 1, t_water, p_water, d_water)

    R = gl%req(gl%el%solpos)! / gl%sea%wm_sea!gl%wm(1)
    g = g_brine(gl,t,d_w, p)
    d2gdt2 =  d2brine_dt2(gl, t, d_w, p)!cp und T d2brine_dt2
    d2gdp2 = d2_brine_dp2(gl, t, d_w, p)! aus v brine
    dgdp = dbrine_dp(gl, T, d_w, p)
    d2gdtdp = d2_brine_dtdp(gl, t, d_w,p)

    A02_el = a02_gibbs_trans(gl, t, p, dgdp, R, d2gdp2, d2gdt2, d2gdtdp, g) 
    gl%gecarrier = seaflagorig

    d = d_mixsave

    end function
    !*************************************************************************************************************


    !*******************************************************************************************************
    !----------------------------------------------------------------------------------------
    !Gerenral functions for converting gibbs to helmholtz
    !Props have to be given to this functions via the gibbs models with futher routines
    !----------------------------------------------------------------------------------------
    !!Benedikt July 2018

    double precision function a_gibbs_trans(gl,t, p, R, g, dgdp)
    !function for calculating  the reduced helmholtz energy from gibbs
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, g, dgdp, R

    a_gibbs_trans = ( g- p*1.d6*dgdp ) / R / T

    end function a_gibbs_trans
    !****************************************************************************************

    !*****************************************************************************************
    double precision function a10_gibbs_trans(gl,t, p, R, g, dgdt, dgdp)
    !function for calculating  tau alpha tau from gibbs equations, differs from A10 !!!
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, g, dgdp, R, dgdt

    a10_gibbs_trans =  (( g - p*1.d6 * dgdp ) / (R * t) )- ( t * (dgdt / ( R * t) ) )

    end function a10_gibbs_trans
    !*******************************************************************************************

    !*******************************************************************************************
    double precision function a01_gibbs_trans(gl, t, p, dgdp, R,g)
    !founction for calculating A01 from gibbs eqations
    !Benedikt July 2018

    implicit none
    type(type_gl) :: gl
    double precision :: t, p, g, dgdp, R

    p = gl%gepress

    a01_gibbs_trans =(p * 1.d6*dgdp) / (R *t) - 1.d0

    end function
    !*************************************************************************************************


    !*************************************************************************************************
    double precision function a20_gibbs_trans(gl,t, p, R, d2gdt2, d2gdp2, d2gdtdp)
    !founction for calculating tau^2 alpha tautau from gibbs eqations c_v / R  , differs from A20!!
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p, g, dgdp, R, d2gdp2, d2gdt2, d2gdtdp, dgdt

    a20_gibbs_trans = - ( T * ( (d2gdtdp**2.d0 - (d2gdt2 * d2gdp2) )* d2gdp2**(-1.d0) ) ) / R

    end function
    !***************************************************************************************************


    !***************************************************************************************************
    double precision function a02_gibbs_trans(gl, t, p, dgdp, R, d2gdp2, d2gdt2, d2gdtdp, g)
    !founction for calculating A02 from gibbs eqations
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p,  g, dgdp, R, d2gdp2, d2gdt2, d2gdtdp, a01

    a01  = a01_gibbs_trans(gl, t, p, dgdp, R,g)
    a02_gibbs_trans =  ( ((dgdp**2.d0) / (- d2gdp2 * R * T ))  - 2.d0 * a01  ) -1.d0    

    end function
    !**********************************************************************************************


    !**********************************************************************************************
    double precision function a11_gibbs_trans(gl, t, R, d2gdt2, a01, a20, a02)
    !function for calculating A11 from gibbs eqations
    !Benedikt July 2018


    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, g, dgdp, R, d2gdp2, d2gdtdp, dgdt, d2gdt2, a01, a02, a20

    a11_gibbs_trans = sqrt( (  ((-t * d2gdt2)/R)+ a20 ) * (1.d0+2.d0*a01+a02) ) + 1.d0 + a01  

    end function
    !****************************************************************************************************
    !****************************************************************************************************


    !*************************************************************************
    !sub for getting pressure_reac and reduced parameters for seawater
    subroutine partial_reac_press(gl, t, p, d, ncomp, t_w, p_w, d_w)
    !-------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d, delta, tau, d_w, t_w, p_w, rhored_orig, tred_orig, tc, rhoc, torig, dorig
    integer :: ncomp
    logical :: calcflag

    rhored_orig = gl%rhoredmix
    tred_orig = gl%tredmix
    calcflag = gl%gecarrier
    torig = t
    dorig = d

    delta = d / gl%rhoredmix
    tau = gl%tredmix / t

    d_w = delta * gl%rhored(1)
    t_w = gl%tc(1) / tau

    gl%rhoredmix = gl%rhored(1)
    gl%tredmix = gl%tc(1)

    !call reduced_parameters_calc(gl, 300.d0)

    gl%gecarrier = .false.

    p_w = p_calc(gl, t_w, d_w, 1)

    gl%gecarrier = calcflag

    gl%rhoredmix = rhored_orig
    gl%tredmix = tred_orig
    t = torig
    d = dorig

    end subroutine !partial_reac_press(gl, t, p, d, ncomp)
    !*********************************************************************************************


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!NEW FUNCTIONS FOR THE GE CARRIER  changed to previous functions for seawater
    !




    !**************************************************************
    !
    !Functions for seawater derivatives

    !-----------------------------------------------------------------------------
    double precision function A_sea(gl, t,d, p)
    !----------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl

    double precision :: t, p, R, g, dgdp, test, d_w, d, t_water, p_water, d_water, ai, ar, d_mixsave, a
    double precision, dimension(nderivs)::FNRDER, FNIDER, FNRDER_w, FNIDER_w
    integer, dimension(nderivs) :: Getderr
    integer, dimension(10) :: getderi
    logical :: seaflagorig

    d_mixsave = d                                   !from old props_Derivs
    seaflagorig = gl%gecarrier                            !from old props_Derivs
    !from old props_Derivs
    gl%gecarrier = .false.                            !from old props_Derivs
    d = rhomix_calc(gl, t, gl%sea%seap, 0.d0, 1, 1) !from old props_Derivs
    gl%gecarrier = seaflagorig                            !from old props_Derivs

    seaflagorig = gl%gecarrier  !raus?

    gl%gecarrier = .false.
    gl%sea%wm_sea = molar_sea(gl)
    R = gl%req(1) / gl%sea%wm_sea
    d_w = d

    !call partial_reac_press(gl, t, p, d, 1, t_water, p_water, d_water)
    g = G_SEA_CALC(gl, T, p, d_w)
    dgdp = dgdp_sea(gl, t, p, d_w)

    A_sea = a_gibbs_trans(gl, t, p, R, g, dgdp) !* gl%molfractions(1) !now here *molfrac

    d = d_mixsave

    gl%gecarrier = seaflagorig

    end function A_sea
    !****************************************************************************************


    !-----------------------------------------------------------------------------
    double precision function A_el(gl, t,d, p)
    !----------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, R, g, dgdp, test, d_w, d, t_water, p_water, d_water, ai, ar, d_mixsave, deybe, a
    double precision, dimension(nderivs)::FNRDER, FNIDER, FNRDER_w, FNIDER_w
    integer, dimension(nderivs) :: Getderr
    integer, dimension(10) :: getderi
    logical :: seaflagorig
    integer :: pos

    pos = gl%el%solpos
    p = gl%el%press
    d_mixsave = d                                   !from old props_Derivs
    seaflagorig = gl%gecarrier                            !from old props_Derivs
    !from old props_Derivs
    gl%gecarrier = .false.                            !from old props_Derivs
    d = rhomix_calc(gl, t, p, 0.d0, 1, pos) !from old props_Derivs
    gl%gecarrier = seaflagorig                            !from old props_Derivs

    seaflagorig = gl%gecarrier  !raus?

    gl%gecarrier = .false.

    R = gl%req(1)
    d_w = d



    g = g_brine(gl, T, d_w, p)
    dgdp = dbrine_dp(gl, t, d_w, p)

    A_el = a_gibbs_trans(gl, t, p, R, g, dgdp) !* gl%molfractions(1)

    d = d_mixsave

    gl%gecarrier = seaflagorig
    end function A_el
    !****************************************************************************************


    !****************************************************************************************
    !new for carrier
    double precision function a_gibbs(gl,t, d, p)
    !function for calculating  the reduced helmholtz energy from gibbs
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d, g, dgdp

    select case(gl%modelflag)
    case(1)
        a_gibbs = A_sea(gl, t, d, p)*gl%molfractions(1)
    case(2)
        a_gibbs = A_el(gl, t, d, p)*gl%molfractions(1)
    end select

    end function a_gibbs
    !****************************************************************************************



    !#################################################################################################
    !#################################################################################################


    !*****************************************************************************************
    double precision function a10_gibbs(gl, t, d, p)
    !function for calculating  tau alpha tau from gibbs equations, differs from A10 !!!
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d

    select case(gl%modelflag)
    case(1)
        a10_gibbs = a10_sea(gl,t, d, p)*gl%molfractions(1)
    case(2)
        a10_gibbs = a10_el(gl, t, d, p)*gl%molfractions(1)
    end select

    end function a10_gibbs
    !*******************************************************************************************

    !*******************************************************************************************
    double precision function a01_gibbs(gl, t, d, p)
    !founction for calculating A01 from gibbs eqations
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d

    select case(gl%modelflag)
    case(1)
        a01_gibbs = a01_sea(gl,t, d, p)*gl%molfractions(1)
    case(2)
        a01_gibbs = a01_el(gl, t, d, p)*gl%molfractions(1)
    end select

    end function
    !*************************************************************************************************


    !*************************************************************************************************
    double precision function a20_gibbs(gl, t, d, p)
    !founction for calculating tau^2 alpha tautau from gibbs eqations c_v / R  , differs from A20!!
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d

    select case(gl%modelflag)
    case(1)
        a20_gibbs = a20_sea(gl, t, d, p)*gl%molfractions(1)
    case(2)
        a20_gibbs = a20_el(gl, t, d, p)*gl%molfractions(1)
    end select

    end function
    !***************************************************************************************************


    !***************************************************************************************************
    double precision function a02_gibbs(gl, t, d, p)
    !founction for calculating A02 from gibbs eqations
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d

    select case(gl%modelflag)
    case(1)
        a02_gibbs =  a02_sea(gl, t, d, p)*gl%molfractions(1)
    case(2)
        a02_gibbs = a02_el(gl,t, d, p)*gl%molfractions(1)
    end select

    end function
    !**********************************************************************************************


    !**********************************************************************************************
    double precision function a11_gibbs(gl, t, d, p)
    !function for calculating A11 from gibbs eqations
    !Benedikt July 2018
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d

    select case(gl%modelflag)
    case(1)
        a11_gibbs = a11_sea(gl,t, d, p)*gl%molfractions(1)
    case(2)
        a11_gibbs = a11_el(gl, t, d, p)*gl%molfractions(1)
    end select

    end function
    !****************************************************************************************************
    !****************************************************************************************************


    !*****************************************************************************************************
    !Routine for the cancellation of double calculated parts for in Helmholtz-Gibbs mixutres
    !Benedikt 2020
    subroutine solvent_cor(gl, t, p, GETDERI, GETDERR, corr_der_i, corr_der_r)
    implicit none
    type(type_gl) :: gl
    double precision :: t, d, dw, p
    double precision, dimension(nderivs):: FNRDER_w,  corr_der_r
    double precision, dimension(nderivsi) :: corr_der_i, FNIDER_w
    integer, dimension(nderivs) :: Getderr
    integer, dimension(10) :: getderi

    dw = rhomix_calc(gl, t, p, 0.d0, 1, 1)

    CALL FNRDERIVS(gl,T,Dw,GETDERR,FNRDER_w,1)   !subroutine calculates derivatives of residual part
    CALL FNIDERIVS(gl,T,Dw,GETDERI,FNIDER_w,1)   !subroutine calculates derivatives of ideal part

    corr_der_i = fnider_w * gl%molfractions(1)
    corr_der_r = fnrder_w * gl%molfractions(1)


    end subroutine
    !******************************************************************************************************
    !******************************************************************************************************


    !*****************************************************************************
    !function for molar density of seawater and other brines
    !Benedikt 2018
    double precision function brine_dens(gl, t, p, d)
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d
    double precision ::  v_helm, v_salt, v_brine, x_brine, x_helm, wm_mix_helm, wm_water, wm_brine_mix, d_helm, prop_out, vw, vb, dv, g, a
    !

    d_helm = d
 
    a = a_calc(gl, t, d_helm, 0)
    g = g_calc(gl, t, d_helm, 0)
    prop_out = (a-g)/(-p*1.d6)
    prop_out = 1.d0/prop_out

    !end if
    brine_dens = prop_out
    end function brine_dens
    !********************************************************************************


    !*******************************************************************
    !subroutine for chemical potential of brine mixtures
    !Benedikt, 01/2019
    !*******************************************************************
    subroutine chempot_num_reac(gl, t, d, chempot_num, n_phase)
    implicit none

    type(type_gl) :: gl

    double precision :: t, d, p, sal, sal_p, sal_m, numdiff, sal_save, salinity_shifted, Rmix, sum_moles, numdiff_abs, d_save, gm, gp, mm, mp, gbrine, chempot_factor
    double precision, dimension (30) :: chempot_num, moles_p, moles_m, mol_save, a_p, a_m, moles_norm, moles_var_p, moles_var_m, chempot, diff, moles_var_2p, moles_var_2m, moles_2m,&
        & moles_2p, a_2m, a_2p, numdif_vec, h_m, s_m, g_diff, gibbs_point, substance

    integer :: i, n_phase
    logical :: calcflag

    numdiff_abs = 1.d-4
    mol_save = gl%molfractions
    !sal_save = gl%sea%salinity
    moles_m = gl%molfractions
    moles_p = gl%molfractions
    chempot_num = 0.d0
    chempot = 0.d0
    a_p = 0.d0
    a_m = 0.d0
    a_2m = 0.d0
    a_2p = 0.d0
    numdif_vec = 0.d0
    moles_var_p = 0.d0
    moles_var_m = 0.d0
    !calcflag = gl%seacalc
    d_save = d
    chempot_factor = 1.d0

    if(gl%seawater) then
        calcflag = gl%seacalc
        sal_save = gl%sea%salinity
        p = gl%sea%seap
        !elseif(gl%el_present) then
        !end if

        if((gl%ncomp == 1) .and. (n_phase == 1)) then

            if(gl%modelflag == 1) then

                chempot_num(1) = chempot_SEA_CALC(gl, T, gl%sea%seap, D) * gl%wm(1)

            else if(gl%modelflag == 2 ) then
                calcflag = gl%gecarrier
                gl%gecarrier = .false.
                chempot_num(1) = chempot_water_brine(gl, t, d, p)
                gl%gecarrier = calcflag
            end if

        elseif(n_phase == 1) then
            numdif_vec = 1.d-3!1.d-3

            !p = gl%sea%seap
            call reduced_parameters_calc(gl, t)
            !p = p_calc(gl, t,d,0)
            !1rst component
            do i=1,gl%ncomp !positive loop
                chempot_factor = 1.d0
                moles_p = gl%molfractions
                numdiff = numdiff_abs
                !numdif_vec(i) = dmod( gl%molfractions(i),10.d0 )* numdiff_abs
                moles_p(i) = gl%molfractions(i) * (1.d0 + numdif_vec(i))
                if((gl%seawater) .and. (n_phase == 1)) then
                    gl%seacalc = .true.
                    gl%gecarrier = gl%seacalc
                elseif((gl%el_present) .and. (n_phase == 1)) then
                    gl%gecarrier = .true.
                    !salinity = gl%el%molality(1)
                    !gl%el%x_salt = gl%molfractions(1) *  gl%wm(1) * gl%el%molality(1)
                end if
                call molnorm(gl, moles_p, moles_norm, salinity_shifted, sum_moles)
                if((gl%seacalc) .and. (i==1)) then
                    gl%sea%salinity = salinity_shifted
                elseif((gl%el_present)) then
                    if(i == gl%el%solpos) then
                        gl%el%molality(1) = salinity_shifted
                        chempot_Factor = gl%wm(gl%el%solpos)
                    end if
                end if
                !moles_var_p(i) = sum_moles - 1.d0!moles_norm(i)*sum_moles
                moles_var_p(i) = moles_norm(i)*sum_moles
                gl%molfractions = moles_norm
                call R_mix_calc(gl, Rmix)
                call reduced_parameters_calc(gl, t)

                d = rhomix_calc(gl, t, p, 0.d0, n_phase, 0)
                if (i==1) then
                    a_p(i) = g_calc(gl, t, d, 0)*sum_moles*chempot_factor
                else
                    a_p(i) = g_calc(gl, t, d, 0)*sum_moles
                end if
                gl%molfractions = mol_save
                if(gl%seawater) then
                    gl%sea%salinity = sal_save
                    gl%seacalc = calcflag
                    gl%gecarrier = gl%seacalc
                elseif(gl%el_present) then
                    gl%el%molality(1) = sal_save
                    gl%gecarrier = calcflag
                end if
            end do

            do i=1,gl%ncomp !negative loop
                chempot_factor = 1.d0
                moles_m = gl%molfractions
                !gibbs_point(i) = g_calc(gl, t, d, 0)
                !numdiff = gl%molfractions(i) *numdiff_abs
                numdiff = numdiff_abs
                moles_m(i) = gl%molfractions(i) * (1.d0 - numdif_vec(i))
                if((gl%seawater) .and. (n_phase == 1)) then
                    gl%seacalc = .true.
                    gl%gecarrier = gl%seacalc
                elseif((gl%el_present) .and. (n_phase == 1)) then
                    gl%gecarrier = .true.
                end if
                call molnorm(gl, moles_m, moles_norm, salinity_shifted, sum_moles)
                substance(i) = sum_moles
                if((gl%seacalc) .and. (i==1)) then
                    gl%sea%salinity = salinity_shifted
                elseif((gl%el_present)) then
                    if(i == gl%el%solpos) then
                        gl%el%molality(1) = salinity_shifted
                        chempot_Factor = gl%wm(gl%el%solpos)
                    end if
                end if
                !moles_var_m(i) = sum_moles - 1.d0!moles_norm(i)*sum_moles
                moles_var_m(i)= moles_norm(i)*sum_moles
                gl%molfractions = moles_norm
                call R_mix_calc(gl, Rmix)
                call reduced_parameters_calc(gl, t)
                d = rhomix_calc(gl, t, p, 0.d0, n_phase, 0)
                if(i==1) then
                    a_m(i) = g_calc(gl, t, d, 0)*sum_moles*chempot_factor
                else
                    a_m(i) = g_calc(gl, t, d, 0)*sum_moles
                end if
                gl%molfractions = mol_save
                if(gl%seawater) then
                    gl%sea%salinity = sal_save
                    gl%seacalc = calcflag
                    gl%gecarrier = gl%seacalc
                elseif(gl%el_present) then
                    gl%el%molality(1) = sal_save
                    gl%gecarrier = calcflag
                end if
            end do

            chempot_num(1:gl%ncomp) = ( a_p(1:gl%ncomp) - a_m(1:gl%ncomp) ) / ( moles_var_p(1:gl%ncomp) - moles_var_m(1:gl%ncomp) )!( moles_var_p(1:gl%ncomp) + moles_var_m(1:gl%ncomp))!

            gl%molfractions = mol_save

            call reduced_parameters_calc(gl, t)
            d = d_save
        end if
   ! elseif(gl%el_present)then
    !    p = gl%gepress
        !call partial_molar_gibbs(gl, t, d, p, chempot)
        !call reduced_parameters_calc(gl, t)
        !call R_mix_calc(gl, Rmix)
        !p = p_calc
        !d = rhomix_calc(gl, t, p, 0.d0, n_phase,0)
        !call chempot_calc(gl, t, d,chempot,0)
        !
        !chempot_num = chempot * Rmix * t
    end if

    !call partial_molar_gibbs(gl, t, d, p, chempot_num)!, errorflag)


    !end if

    end subroutine !chempot_num_reac
    !*****************************************************************************
    !****************************************************************************


    !******************************************************************************
    !subroutine for norming the composition and ajusting the reactive components
    !Benedikt, 01/2019
    subroutine molnorm(gl, moles_in, moles_norm, salinity_shifted, sum_moles)

    implicit none

    type(type_gl) :: gl

    double precision :: salinity, salinity_shifted, sum_moles, wm_sea, n_water, n_salt
    double precision, dimension(30) :: moles_in, moles_norm

    moles_norm = 0.d0

    sum_moles = sum(moles_in(:))

    if(gl%seacalc) then

        wm_sea = molar_sea(gl)

        n_water = gl%molfractions(1) * sum_moles

        n_salt = n_water * gl%sea%x_s

        n_water = n_water * gl%sea%x_w

        salinity_shifted = (n_salt * gl%sea%wm_salt) / ( (n_salt * gl%sea%wm_salt) + ( n_water * gl%wm(1) ))
        write(*,*) salinity_shifted
    elseif(gl%el_present) then

        !moles_norm = moles_in / sum_moles

        n_water = gl%molfractions(1) * sum_moles

        !n_salt = gl%el%x_salt
        n_salt = gl%el%n2
        !gl%el%mol_in = gl%el%molality(1)

        salinity_shifted = n_salt / (n_Water * gl%wm(1))

    end if

    moles_norm = moles_in / sum_moles

    end subroutine !molnorm
    !*************************************************************************
    
    
    !*********************************************************************************
    subroutine chempot_brine(gl, t, d, chempot)
    !chempot for brines
    implicit none
    type(type_gl)::gl
    double precision ::t, d, p, add_part_brine, R, mw, m, xw, nw
    double precision, dimension(30):: chempot
    integer::i
    mw = gl%wm(1)
    nw = 1.d0/mw
    r = gl%req(1)
    m = gl%el%molality(1)
    xw = gl%molfractions(1)
    xw = 1.d0!xw*(1.d0 - (m/(nw+m)))
    
    call chempot_calc(gl, t, d, chempot, 0)
    p = gl%gepress
    add_part_brine = (- osmotic_coeff(gl, t, d, p)*dlog(1.d0 + Mw*m)+ 2.d0*m*Mw)*R*T
    !add_part_brine = R*T*2.d0*m*M
    
    chempot(1) = chempot(1) + xw* (add_part_brine)
    
    end subroutine
    !*********************************************************************************
    !*********************************************************************************
    


    !!**************************************************************************
    !subroutine partial_molar_gibbs(gl, t, d, p, chempot)
    !!function for the calcuation of the partial molar gibbs energy of mixture; Benedikt 10/2020
    !implicit none
    !type(type_gl)::gl
    !double precision :: t, d, p, g, xs, xsp, xsm, numdiff, sum_x_dg_dx, minfrac, saltsave, x_salt, x_wb, x_1_orig
    !double precision, dimension(2) :: salt_adjusted, molality_adjusted !(1) +; (2) -
    !double precision, dimension(30) :: chempot, x_orig, x_p, x_m, dg_dx, gp, gm, x_m_save, x_p_save
    !integer :: i, j
    !logical :: absdiff, brine_flag
    !
    !
    !absdiff = .false.
    !x_orig = gl%molfractions
    !chempot = 0.d0
    !numdiff = 1.d-4
    !x_p_save = 0.d0
    !x_m_save = 0.d0
    !gp = 0.d0
    !gm = 0.d0
    !dg_dx = 0.d0
    !brine_flag = gl%gecarrier
    !if(gl%el_present) gl%gecarrier = .true.
    !
    !!avoid problems when one comp has a high concentration (x_com*(1+numdiff) would exceed 1)
    !minfrac = minval(x_orig(1:gl%ncomp))
    !if(minfrac .le. (numdiff/10.d0)) then
    !    numdiff = minfrac / 100.d0
    !    absdiff = .true.
    !end if
    !
    !g = g_calc(gl, t, d, 0)
    !
    !if(gl%el_present) then
    !    saltsave = gl%el%molality(1)
    !    x_1_orig = x_orig(1)
    !    x_salt = gl%el%molality(1) / (1.d0/gl%wm(1) + gl%el%molality(1) )
    !    x_wb = 1.d0 - x_salt
    !    call adjust_salt(gl, x_orig, saltsave, salt_adjusted, molality_adjusted)
    !    gl%el%molality(1) = molality_adjusted(1)
    !    gp(1) = g_calc(gl, t, d, 0)          !2, as water is saved in 1st position
    !    gl%el%molality(1) = molality_adjusted(2)
    !    gm(1) = g_calc(gl, t, d, 0)
    !    gl%el%molality(1) = saltsave
    !    dg_dx(1) = (gp(1) - gm(1)) / (salt_adjusted(1) - salt_adjusted(2) )
    !    x_p_save(1) = salt_adjusted(1)
    !    x_m_save(1) = salt_adjusted(2)
    !    !x_orig(1) = x_Salt
    !    gl%el%molality(1) = saltsave
    !
    !    !loops up to ncomp are needed, as a pseudopure substance is involved
    !    do i=2,gl%ncomp !loop for upper gibbs energy
    !        x_p = x_orig
    !        if(absdiff == .true.) then
    !            x_p(i) = x_orig(i) + numdiff
    !        else
    !            x_p(i) = x_orig(i) * (1.d0 + numdiff)
    !        end if
    !        x_p(1) = 1.d0 - sum(x_p(2:gl%ncomp))
    !        gl%molfractions = x_p
    !        x_p_save(i) = x_p(i)
    !        call reduced_parameters_calc(gl, t)
    !        d = rhomix_calc(gl, t, p,0.d0, 0, 0)
    !        gp(i) = g_calc(gl, t, d,0)
    !    end do !i loop, upper
    !
    !    do i=2,gl%ncomp   !minus loop
    !        x_m = x_orig
    !        if (absdiff) then
    !            x_m(i) = x_orig(i) - numdiff
    !        else
    !            x_m(i) = x_orig(i) * (1.d0 - numdiff)
    !        end if
    !        x_m(1) = 1.d0 - sum(x_m(2:gl%ncomp))
    !        gl%molfractions = x_m
    !        x_m_save(i) = x_m(i)
    !        call reduced_parameters_calc(gl, t)
    !        d = rhomix_calc(gl, t, p,0.d0, 0, 0)
    !        gm(i) = g_calc(gl, t, d,0)
    !    end do !i  lower loop
    !
    !    !---------------  Variaton loops done!
    !    do i=1,gl%ncomp
    !        dg_dx(i) = ( gp(i) - gm(i) ) / (x_p_save(i) - x_m_save(i))
    !    end do
    !
    !    sum_x_dg_dx = 0.d0
    !    sum_x_dg_dx = x_salt * dg_dx(1)
    !    do i=2,gl%ncomp  !loop for x_i * dg/dxi
    !        sum_x_dg_dx = sum_x_dg_dx + x_orig(i) * dg_dx(i)
    !    end do !i
    !
    !    chempot(1) = g - sum_x_dg_dx
    !    do i=2,gl%ncomp  !loop for calculation of the chemical potentail with partial molar gibbs up to ncomp-1
    !        chempot(i) = g + dg_dx(i) - sum_x_dg_dx
    !    end do !i
    !
    !
    !    !end if
    !
    !
    !
    !
    !    !---------------------------------------------------------------------------------
    !    !normal way, similar to chempot in mixtures (no pseudo-pure comps)
    !elseif(.not. gl%el_present) then
    !    !    if(gl%ncomp == 1) then
    !    !    if(gl%el_present) then
    !    !        gl%gecarrier = .true.
    !    !        g = g_calc(gl, t, d, 0)
    !    !        call adjust_salt(gl, x_orig, saltsave, salt_adjusted, molality_adjusted)
    !    !        gl%el%molality(1) = molality_adjusted(1)
    !    !        gp(2) = g_calc(gl, t, d,1)          !2, as water is saved in 1st position
    !    !        gl%el%molality(1) = molality_adjusted(2)
    !    !        gm(2) = g_calc(gl, t, d,1)
    !    !        gl%el%molality(1) = saltsave
    !    !        dg_dx(2) = (gp(2) - gm(2)) / (salt_adjusted(1) - salt_adjusted(2) )
    !    !        sum_x_dg_dx = x_salt * dg_dx(2)
    !    !        gl%gecarrier = brine_flag
    !    !    end if
    !    !else
    !    !normal procedure for common mixture, no pseudo-pure substances, etc.
    !    do i=1,gl%ncomp-1   !plus loop
    !        x_p = x_orig
    !        if(gl%el_present) gl%el%molality(1) = saltsave
    !        if(absdiff == .true.) then
    !            x_p(i) = x_orig(i) + numdiff
    !        else
    !            x_p(i) = x_orig(i) * (1.d0 + numdiff)
    !        end if
    !        x_p(gl%ncomp) = 1.d0 - sum(x_p(1:gl%ncomp-1))
    !        gl%molfractions = x_p
    !        x_p_save(i) = x_p(i)
    !        call reduced_parameters_calc(gl, t)
    !        d = rhomix_calc(gl, t, p,0.d0, 0, 0)
    !        !if((gl%el_present) .and. (i==1)) then
    !        !    gl%el%molality(1) = molality_adjusted(1)
    !        !end if
    !        gp(i) = g_calc(gl, t, d,0)
    !    end do !i
    !
    !    do i=1,gl%ncomp-1   !minus loop
    !        x_m = x_orig
    !        if(gl%el_present) gl%el%molality(1) = saltsave
    !        if (absdiff) then
    !            x_m(i) = x_orig(i) - numdiff
    !        else
    !            x_m(i) = x_orig(i) * (1.d0 - numdiff)
    !        end if
    !        x_m(gl%ncomp) = 1.d0 - sum(x_m(1:gl%ncomp-1))
    !        gl%molfractions = x_m
    !        x_m_save(i) = x_m(i)
    !        call reduced_parameters_calc(gl, t)
    !        d = rhomix_calc(gl, t, p,0.d0, 0, 0)
    !        !if((gl%el_present) .and. (i==1)) then
    !        !    gl%el%molality(1) = molality_adjusted(2)
    !        !end if
    !        gm(i) = g_calc(gl, t, d,0)
    !    end do !i
    !
    !    gl%molfractions = x_orig
    !
    !    do i=1,gl%ncomp -1
    !        dg_dx(i) = ( gp(i) - gm(i) ) / (x_p_save(i) - x_m_save(i))
    !    end do
    !
    !    sum_x_dg_dx = 0.d0
    !    do i=1,gl%ncomp-1  !loop for x_i * dg/dxi
    !        sum_x_dg_dx = sum_x_dg_dx + x_orig(i) * dg_dx(i)
    !    end do !i
    !
    !    do i=1,gl%ncomp-1  !loop for calculation of the chemical potentail with partial molar gibbs up to ncomp-1
    !        chempot(i) = g + dg_dx(i) - sum_x_dg_dx
    !    end do !i
    !    chempot(gl%ncomp) = g - sum_x_dg_dx
    !end if
    !!last (or second) component with different formulation (analog to x_n = 1-sum(x_i...N-1)
    !
    !
    !gl%gecarrier = brine_flag
    !end subroutine !partial molar gibbs
    !!****************************************************************************
    !!****************************************************************************

    !****************************************************************************
    subroutine adjust_salt(gl, x_in, mol_salt_orig, salt_adjusted, molality_adjusted)
    implicit none
    type(type_gl) :: gl
    double precision :: mol_salt_orig, numdiff
    double precision, dimension(2) :: salt_adjusted, molality_adjusted
    double precision, dimension(30) :: x_in

    numdiff = 1.d-4

    if(gl%seawater) then
        molality_adjusted(1) = mol_salt_orig * (1.d0 + numdiff)
        molality_adjusted(2) = mol_salt_orig * (1.d0 - numdiff)
    elseif(gl%el_present) then
        molality_adjusted(1) = mol_salt_orig * (1.d0 + numdiff) !/ (1.d0/gl%wm(1) + x_salt_orig * (1.d0 + numdiff))
        molality_adjusted(2) = mol_salt_orig * (1.d0 - numdiff) !/ (1.d0/gl%wm(1) + x_salt_orig * (1.d0 - numdiff))
        salt_adjusted(1) = molality_adjusted(1) / ( 1.d0/gl%wm(1) + molality_adjusted(1))
        salt_adjusted(2) = molality_adjusted(2) / ( 1.d0/gl%wm(1) + molality_adjusted(2))
    end if

    end subroutine !adjust_salt


    end module gibbsderivs_module
