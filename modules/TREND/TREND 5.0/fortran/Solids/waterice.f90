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

    ! module for file waterice.f90
    module waterice_module
    !global use inclusion
    use module_all_types
    use module_regula_falsi
    use module_regula_falsi_support

    contains






    Subroutine Set_coef_H20(gl,REOS_or_GERG)




    implicit none

    type(type_gl) :: gl


    integer:: REOS_or_GERG    !1 = REOS --> Reference equations are used
    !2 = GERG --> GERG eq. of state are used

    !Molar weight in kg/mol
    gl%M_H20 = 18.01528D-3

    gl%Ttr_water = 273.16D0
    gl%ptr_water = 611.657D0!4771007894D0
    gl%p0_water = 101325.D0

    if (REOS_or_GERG == 1) then
        !Reference point is linked to the reference point of the IAPWS-95
        gl%g00_waterice = -0.632020233335886D6
        !Reference point is linked to the reference point of the IAPWS-95
        gl%s0 = -0.332733756492168D4 !3333.18160308627D0
        !    elseif (REOS_or_GERG == 2) then
        !!        !REFERENCE POINT SET ACCORDING TO THE FLUID FILE!!
        !!        !Reference point is linked to the reference point of the GERG 2004
        !!!!!!!!!!g00_waterice = -644416.505137970D0 ! (old value)
        !        g00_waterice = -637966.096687019D0
        !!        !Reference point is linked to the reference point of the GERG 2004
        !!!!!!!!!!s0 = -3374.56552904248D0 ! (old value)
        !        s0 = -3350.95149898262D0
    elseif (REOS_or_GERG == 2) then
        !ORIGINAL REFERENCE POINT OF THE GERG 2004!!
        !Reference point is linked to the reference point of the GERG 2004
        gl%g00_waterice = -3185921.68104710D0
        !Reference point is linked to the reference point of the GERG 2004
        gl%s0 = -10312.6873872297D0
    elseif (REOS_or_GERG == 3) then
        !REFERENCE POINT ACCORDING TO THE FLUID FILE OF WATER, but for the SRK EOS
        gl%g00_waterice = 1780890.87494544D0
        !REFERENCE POINT ACCORDING TO THE FLUID FILE OF WATER, but for the SRK EOS
        gl%s0 = 5567.63044407941D0
    end if

    gl%g01_waterice = 0.655022213658955D0
    gl%g02_waterice = -0.189369929326131D-7
    gl%g03_waterice = 0.339746123271053D-14
    gl%g04_waterice = -0.556464869058991D-21

    gl%t1 = (0.368017112855051D-1, 0.510878114959572D-1)
    gl%r1 = (0.447050716285388D2, 0.656876847463481D2)
    gl%t2 = (0.337315741065416D0, 0.335449415919309D0)
    gl%r20 = (-0.725974574329220D2, -0.781008427112870D2)
    gl%r21 = (-0.557107698030123D-4, 0.464578634580806D-4)
    gl%r22 = (0.234801409215913D-10, -0.285651142904972D-10)

    End subroutine

    DOUBLE PRECISION FUNCTION g_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the Gibbs Free Energy (chemical potential) of pure water
    !   ice in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0
    double precision ::Real_part

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    gl%g0_waterice = gl%g00_waterice + gl%g01_waterice * (pi_loc - pi_loc_0) + gl%g02_waterice * (pi_loc - pi_loc_0)**2 + gl%g03_waterice * (pi_loc - pi_loc_0)**3 + gl%g04_waterice * (pi_loc - pi_loc_0)**4
    gl%r2 = gl%r20 + gl%r21 * (pi_loc - pi_loc_0) + gl%r22 * (pi_loc - pi_loc_0)**2

    Real_part =  real((gl%r1 * ((gl%t1 - tau) * cdlog(gl%t1 - tau) + (gl%t1 + tau) * cdlog(gl%t1 + tau) - 2.D0 * gl%t1 * cdlog(gl%t1) - tau**2 / gl%t1) + &
        &  gl%r2 * ((gl%t2 - tau) * cdlog(gl%t2 - tau) + (gl%t2 + tau) * cdlog(gl%t2 + tau) - 2.D0 * gl%t2 * cdlog(gl%t2) - tau**2 / gl%t2)))

    g_WaterIce = (gl%g0_waterice - gl%s0 * T + gl%Ttr_water * Real_part) * gl%M_H20


    End function

    DOUBLE PRECISION FUNCTION dgdp_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of the Gibbs Free Energy with respect to
    !   pressure p at constant temperature T in [J / mol] / [Pa]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0, g0p
    double complex :: r2p
    double precision ::Real_part

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    g0p = (gl%g01_waterice  + gl%g02_waterice * 2.D0 * (pi_loc - pi_loc_0) + gl%g03_waterice * 3.D0 * (pi_loc - pi_loc_0)**2 + gl%g04_waterice * 4.D0 * (pi_loc - pi_loc_0)**3) / gl%ptr_water
    r2p = (gl%r21  + gl%r22 * 2.D0 * (pi_loc - pi_loc_0)) / gl%ptr_water

    Real_part =  real(r2p * ((gl%t2 - tau) * cdlog(gl%t2 - tau) + (gl%t2 + tau) * cdlog(gl%t2 + tau) - 2.D0 * gl%t2 * cdlog(gl%t2) - tau**2 / gl%t2))

    dgdp_WaterIce = (g0p + gl%Ttr_water * Real_part) * gl%M_H20


    End function

    DOUBLE PRECISION FUNCTION dgdT_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of the Gibbs Free Energy with respect to
    !   Temperature T at constant pressure p in [J / mol] / [K]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0
    double precision ::Real_part

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    gl%r2 = gl%r20 + gl%r21 * (pi_loc - pi_loc_0) + gl%r22 * (pi_loc - pi_loc_0)**2

    Real_part =  real(gl%r1 * (-cdlog(gl%t1 - tau) + cdlog(gl%t1 + tau) - 2.D0 * tau / gl%t1) + &
        &  gl%r2 * (-cdlog(gl%t2 - tau) + cdlog(gl%t2 + tau) - 2.D0 * tau / gl%t2))

    dgdT_WaterIce = (-gl%s0 + Real_part) * gl%M_H20


    End function

    DOUBLE PRECISION FUNCTION d2gdT2_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of the Gibbs Free Energy with respect to
    !   Temperature T at constant pressure p in [J / mol] / [K]^2
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0
    double precision ::Real_part

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    gl%r2 = gl%r20 + gl%r21 * (pi_loc - pi_loc_0) + gl%r22 * (pi_loc - pi_loc_0)**2

    Real_part =  real(gl%r1 * (1.D0 / (gl%t1 - tau) + 1.D0 / (gl%t1 + tau) - 2.D0  / gl%t1) + &
        &  gl%r2 * (1.D0 / (gl%t2 - tau) + 1.D0 / (gl%t2 + tau) - 2.D0  / gl%t2))

    d2gdT2_WaterIce = (Real_part / gl%Ttr_water) * gl%M_H20


    End function



    DOUBLE PRECISION FUNCTION d2gdTdp_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of the Gibbs Free Energy with respect to
    !   Temperature T and pressure p  in [J / mol] / [K] [Pa]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0
    double precision :: Real_part
    double complex :: r2p

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    r2p =  (gl%r21  + gl%r22 * 2.D0 * (pi_loc - pi_loc_0)) / gl%ptr_water

    Real_part =  real(r2p * (-cdlog(gl%t2 - tau) + cdlog(gl%t2 + tau) - 2.D0 * tau / gl%t2))

    d2gdTdp_WaterIce =  Real_part * gl%M_H20


    End function


    DOUBLE PRECISION FUNCTION d2gdp2_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of the Gibbs Free Energy with respect to
    !   pressure p at constant temperature T in [J / mol] / [Pa]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision :: T, p
    double precision :: tau, pi_loc, pi_loc_0, g0pp
    double complex :: r2pp
    double precision ::Real_part

    !call Set_coef_H20()

    tau = T / gl%Ttr_water
    pi_loc = p *1.D6 / gl%ptr_water
    pi_loc_0 = gl%p0_water / gl%ptr_water

    g0pp = (gl%g02_waterice * 2.D0 + gl%g03_waterice * 6.D0 * (pi_loc - pi_loc_0) + gl%g04_waterice * 12.D0 * (pi_loc - pi_loc_0)**2) / gl%ptr_water**2
    r2pp = (gl%r22 * 2.D0) / gl%ptr_water**2

    Real_part =  real(r2pp * ((gl%t2 - tau) * cdlog(gl%t2 - tau) + (gl%t2 + tau) * cdlog(gl%t2 + tau) - 2.D0 * gl%t2 * cdlog(gl%t2) - tau**2 / gl%t2))

    d2gdp2_WaterIce = (g0pp + gl%Ttr_water * Real_part) * gl%M_H20


    End function



    DOUBLE PRECISION FUNCTION v_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar volume of pure solid water
    !   in  m³ / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    v_WaterIce = dgdp_WaterIce(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION s_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar entropy of pure pure solid water
    !   in J / mol K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    s_WaterIce = -dgdT_WaterIce(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION cp_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar heat capacity of pure solid water
    !   in J / mol K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    cp_WaterIce = -T*d2gdT2_WaterIce(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION h_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar enthalpy of pure solid water
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    h_WaterIce = g_waterIce(gl,T,p)-T*dgdT_WaterIce(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION u_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar energy of pure solid water
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    u_WaterIce = g_WaterIce(gl,T,p)-T*dgdT_WaterIce(gl,T,p)-p*1.D6*dgdp_waterice(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION f_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar Helmholtz energy of pure solid water
    !   in J / mol K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    f_WaterIce = g_WaterIce(gl,T,p)-p*1.D6*dgdp_waterice(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION kappa_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the isothermal compressibility of pure water ice
    !   in 1/Pa
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    kappa_WaterIce = -d2gdp2_WaterIce(gl,T,p)/ dgdp_waterice(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION alpha_WaterIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the expansion coefficient of pure Dry ice
    !   in 1 / K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision :: T, p

    alpha_WaterIce = d2gdTdp_WaterIce(gl,T,p)/ dgdp_waterice(gl,T,p)

    END function


    !Andreas February 2014
    Double precision function p_WaterIce(gl,T,rho,errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the pressure in MPa of solid water at given T [K] and
    !   rho [mol/m³]
    !
    !   Input Parameters:
    !
    !   T       -   Temperature in K
    !   rho     -   density in mol/m³
    !
    !----------------------------------------------------------------------------------


    implicit none

    type(type_gl) :: gl

    double precision:: T , rho

    double precision :: p_root, p_min, p_max, p_min_allowed, p_max_allowed
    double precision :: delta_allowed
    type(type_additional_parameters) :: parameters
    integer:: max_iterations, iterations, errval

    !parameters = 0.d0
    Parameters%a_p(1) = T
    Parameters%a_p(2) = rho

    Delta_allowed = 1.d-8
    Max_Iterations = 30

    !pressure interval
    p_min = 0.D0
    p_max = 211.D0
    p_min_allowed = 0.D0
    p_max_allowed = 211.D0

    call Regula_Falsi(gl,rho_diff_WaterIce, p_root, p_min, p_max, Delta_allowed, &
        p_min_allowed, p_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(p_root) < 1.d-10) then
        errval = -18867
    end if

    if ( dabs( (1.D0 / v_WaterIce(gl,T,p_root)) - rho) > 1.d-6 ) then
        errval = -18867
    end if

    p_WaterIce = p_root

    return

    End function


    double precision function rho_diff_WaterIce(gl,p, parameters)
    !----------------------------------------------------------------------------------
    !   Function for the iterative solution of the pressure at given T and rho for dry
    !   ice
    !   Input Parameters:
    !
    !   T       -   Temperature in K
    !   rho     -   density in mol/m³
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: p    ! Inputvariablen
    Double Precision :: T, rho    !
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    T = parameters%a_p(1)
    rho = parameters%a_p(2)
    rho_diff_WaterIce = rho-(1.D0/v_WaterIce(gl,T,p))

    return
    End function



    Double precision function T_WaterIce_ph(gl,p, h, errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the temperature in K of water ice at given p [MPa] and
    !   h [J/mol]
    !
    !   Input Parameters:
    !
    !   p       -   pressure in MPa
    !   h       -   enthalpy in J/mol
    !
    !----------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: p , h

    double precision :: T_root, T_min, T_max, T_min_allowed, T_max_allowed
    double precision :: delta_allowed
    type(type_additional_parameters) :: parameters
    integer:: max_iterations, iterations, errval

    !parameters = 0.d0
    Parameters%a_p(1) = p
    Parameters%a_p(2) = h

    Delta_allowed = 1.d-4
    Max_Iterations = 30

    !Temperature interval
    T_min = gl%Tmin_Waterice
    T_min_allowed = T_min
    T_max = 273.2D0          !Upper temperature limit of water ice equation (about triple point)
    T_max_allowed = T_max

    call Regula_Falsi(gl,h_diff_WaterIce, T_root, T_min, T_max, Delta_allowed, &
        T_min_allowed, T_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(T_root) < 1.d-10) then
        errval = -18867
    end if

    if (dabs(h_WaterIce(gl,T_root, p) - h) > 1.d-2) then
        errval = -18867
    end if

    T_WaterIce_ph = T_root

    return

    End function


    double precision function h_diff_WaterIce(gl,T, parameters)
    !----------------------------------------------------------------------------------
    !   Function for the iterative solution of the temperature at given h and p for dry
    !   ice
    !   Input Parameters:
    !
    !   p       -   pressure in MPa
    !   h       -   enthalpy in J/mol
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: T  ! Inputvariablen
    Double Precision :: h, p      !
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    h = parameters%a_p(2)
    h_diff_WaterIce = h - h_WaterIce(gl,T,p)

    return

    End function


    Double precision function T_WaterIce_ps(gl,p, s, errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the temperature in K of water ice at given p [MPa] and
    !   s [J/molK]
    !
    !   Input Parameters:
    !
    !   p       -   pressure in MPa
    !   s       -   entropy in J/molK
    !
    !----------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: p, s

    double precision :: T_root, T_min, T_max, T_min_allowed, T_max_allowed
    double precision :: delta_allowed
    type(type_additional_parameters) :: parameters
    integer:: max_iterations, iterations, errval

    !parameters = 0.d0
    Parameters%a_p(1) = p
    Parameters%a_p(2) = s

    Delta_allowed = 1.d-4
    Max_Iterations = 30

    !Temperature interval
    T_min = gl%Tmin_Waterice
    T_min_allowed = T_min
    T_max = 273.2D0          !Upper temperature limit of water ice equation (about triple point)
    T_max_allowed = T_max

    call Regula_Falsi(gl,s_diff_WaterIce, T_root, T_min, T_max, Delta_allowed, &
        T_min_allowed, T_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(T_root) < 1.d-10) then
        errval = -18867
    end if

    if (dabs(s_WaterIce(gl,T_root, p) - s) > 1.d-3) then
        errval = -18867
    end if

    T_WaterIce_ps = T_root

    return

    End function


    double precision function s_diff_WaterIce(gl,T, parameters)
    !----------------------------------------------------------------------------------
    !   Function for the iterative solution of the temperature at given s and p for dry
    !   ice
    !   Input Parameters:
    !
    !   p       -   pressure in MPa
    !   s       -   entropy in J/molK
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: T ! Inputvariablen
    Double Precision :: s, p      !
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    s = parameters%a_p(2)
    s_diff_WaterIce = s - s_WaterIce(gl,T,p)

    return

    End function


    end module waterice_module
