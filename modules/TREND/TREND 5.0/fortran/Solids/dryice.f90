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

    ! module for file dryice.f90
    module dryice_module
    !global use inclusion
    use module_all_types
    use module_regula_falsi
    use rhomix_pt_module
    use rhomix_pt_module
    use ancillary_equations_mix_module
    use fnrderivs_module

    contains





    !   Andreas Jäger Dec 2011



    Subroutine Set_coef_CO2(gl,REOS_or_GERG)




    implicit none

    type(type_gl) :: gl


    integer, intent(in):: REOS_or_GERG       !1 = REOS --> Reference equations are used
    !2 = GERG --> GERG eq. of state are used

    !Reference point has been set arbitrarily
    gl%T0_dryice = 150.D0
    gl%p0_dryice = 101325.D0

    if (REOS_or_GERG == 1) then
        !Reference point fitted to the reference eq. for CO2 Span and Wagner (1996)
        gl%g0 = -2.6385478D0!-2.63854779500635D0
        !Reference point fitted to the reference eq. for CO2 Span and Wagner (1996)
        gl%g1 = 4.5088732D0!4.50887318966592D0
        !elseif (REOS_or_GERG == 2) then
        !    !REFERENCE POINT SET ACCORDING TO THE FLUID FILE!!
        !    !Reference point fitted to the GERG 2004
        !    g0 = -2.63776975693671D0
        !    !Reference point fitted to the GERG 2004
        !    g1 = 4.50822211352089D0
    elseif (REOS_or_GERG == 2) then
        !ORIGINAL REFERENCE POINT OF THE GERG 2004!!
        !Reference point is linked to the reference point of the GERG 2004
        gl%g0 = -6.02360919073367D0
        !Reference point is linked to the reference point of the GERG 2004
        gl%g1 = 1.90068705897955D1
    elseif (REOS_or_GERG == 3) then
        !REFERENCE POINT ACCORDING TO THE FLUID FILE OF CO2, but for the SRK EOS
        gl%g0 = -2.98568108789825D0
        !REFERENCE POINT ACCORDING TO THE FLUID FILE OF CO2, but for the SRK EOS
        gl%g1 = 3.54365418008556D0
    end if

    gl%g2 = -2.0109135D0!-2.01091349805817D0
    gl%g3 = -2.7976237D0!-2.79762374688068D0
    gl%g4 = 2.6427834D-1!2.64278343584356D-1
    gl%g5 = 3.8259935D0!3.82599347359082D0
    gl%g6 = 3.1711996D-1!3.1711995750868D-1
    gl%g7 = 2.2087195D-3!2.20871949316193D-3
    gl%g8 = -1.1289668D0!-1.12896682263017D0
    gl%g9 = 9.2923982D-3!9.29239821442875D-3
    gl%g10 = 3.3914617D3!3.39146173274561D3

    gl%exp_n = 7.D0

    gl%g0_a = 3.9993365D-2
    gl%g1_a = 2.3945101D-3
    gl%g2_a = 3.2839467D-1
    gl%g3_a = 5.7918471D-2
    gl%g4_a = 2.3945101D-3
    gl%g5_a = -2.6531689D-3
    gl%g6_a = 1.6419734D-1
    gl%g7_a = 1.7594802D-1
    gl%g8_a = 2.6531689D-3

    gl%g0_k = 2.2690751D-1
    gl%g1_k = -7.5019750D-2
    gl%g2_k = 2.6442913D-1

    gl%const_pi = 4.D0 * atan(1.D0)

    gl%R_CO2 = 8.314472D0

    End subroutine



    DOUBLE PRECISION FUNCTION g_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the Gibbs Free Energy (chemical potential) of pure Dry ice
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    double precision:: f_alpha, K_p
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    f_alpha = f_alpha_calc(gl,theta)
    K_p = K_p_calc(gl,theta)

    g_DryIce= (gl%g0 + gl%g1 * del_theta + gl%g2 * del_theta**2 + &
        & gl%g3 * (dlog((theta**2+gl%g4**2)/(1.D0+gl%g4**2)) - 2.D0*theta/gl%g4*(atan(theta/gl%g4)-atan(1.D0/gl%g4))) + &
        & gl%g5 * (dlog((theta**2+gl%g6**2)/(1.D0+gl%g6**2)) - 2.D0*theta/gl%g6*(atan(theta/gl%g6)-atan(1.D0/gl%g6))) + &
        & gl%g7 * del_pi * (dexp(f_alpha) + K_p*gl%g8) + gl%g9 * K_p * ((pi_lc + gl%g10)**((gl%exp_n-1.D0)/gl%exp_n) - &
        & (1.D0 + gl%g10)**((gl%exp_n-1.D0)/gl%exp_n))) * gl%R_CO2 * gl%T0_dryice

    End function



    DOUBLE PRECISION FUNCTION f_alpha_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the integral of the thermal expansion coefficient alpha
    !   over T at constant p
    !   Input Parameters:
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    f_alpha_calc = gl%g0_a * (theta**2 - 1.D0) &
        & + gl%g1_a*dlog((theta**2-gl%g2_a*theta+gl%g3_a) / (1.D0 - gl%g2_a + gl%g3_a))&
        & + gl%g4_a*dlog((theta**2+gl%g2_a*theta+gl%g3_a) / (1.D0 + gl%g2_a + gl%g3_a)) &
        & + gl%g5_a*(atan((theta-gl%g6_a)/gl%g7_a)-atan((1.D0-gl%g6_a)/gl%g7_a)) &
        & + gl%g8_a*(atan((theta+gl%g6_a)/gl%g7_a)-atan((1.D0+gl%g6_a)/gl%g7_a))

    End function

    DOUBLE PRECISION FUNCTION df_alpha_dtheta_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of f_alpha with respect to reduced Temperture Theta
    !   at constant p
    !   Input Parameters:
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    df_alpha_dtheta_calc = 2.D0 * gl%g0_a * theta + &
        & gl%g1_a * (2.D0*theta-gl%g2_a) / (theta**2-gl%g2_a*theta+gl%g3_a) + &
        & gl%g4_a * (2.D0*theta+gl%g2_a) / (theta**2+gl%g2_a*theta+gl%g3_a) + &
        & gl%g5_a / gl%g7_a / (1.D0 + ((theta - gl%g6_a) / gl%g7_a)**2) + &
        & gl%g8_a / gl%g7_a / (1.D0 + ((theta + gl%g6_a) / gl%g7_a)**2)

    End function

    DOUBLE PRECISION FUNCTION d2f_alpha_dtheta2_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of f_alpha with respect to the reduced Temperture Theta
    !   at constant p
    !   Input Parameters:
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    d2f_alpha_dtheta2_calc = 2.D0 * gl%g0_a + &
        & gl%g1_a * (2.D0*(theta**2-gl%g2_a*theta+gl%g3_a) - (2.D0*theta-gl%g2_a)**2)/(theta**2-gl%g2_a*theta+gl%g3_a)**2 + &
        & gl%g4_a * (2.D0*(theta**2+gl%g2_a*theta+gl%g3_a) - (2.D0*theta+gl%g2_a)**2)/(theta**2+gl%g2_a*theta+gl%g3_a)**2 - &
        & gl%g5_a/gl%g7_a**2/(1.D0 + ((theta - gl%g6_a) / gl%g7_a)**2)**2*2.D0*(theta-gl%g6_a)/gl%g7_a - &
        & gl%g8_a/gl%g7_a**2/(1.D0 + ((theta + gl%g6_a) / gl%g7_a)**2)**2*2.D0*(theta+gl%g6_a)/gl%g7_a
    End function

    DOUBLE PRECISION FUNCTION K_p_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the integral of the temperature dependence of the compressibility
    !   Input Parameters:
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    K_p_calc = gl%g0_k * theta**2 + gl%g1_k * theta + gl%g2_k

    End function


    DOUBLE PRECISION FUNCTION dK_p_dtheta_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of Kp with respect to the reduced temperature
    !   theta
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    dK_p_dtheta_calc = 2.D0 * gl%g0_k * theta + gl%g1_k

    End function

    DOUBLE PRECISION FUNCTION d2K_p_dtheta2_calc(gl,theta)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of Kp with respect to the reduced temperature
    !   theta
    !
    !   theta   -   reduced Temperature (Temp / T0)
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: theta

    d2K_p_dtheta2_calc = 2.D0 * gl%g0_k

    End function

    !
    !
    DOUBLE PRECISION FUNCTION dgdp_DryIce(gl,T,p)
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


    Double precision:: T, p
    double precision:: f_alpha
    double precision:: K_p
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    f_alpha = f_alpha_calc(gl,theta)
    K_p = K_p_calc(gl,theta)


    dgdp_DryIce = gl%R_CO2 * gl%T0_dryice / gl%p0_dryice * (gl%g7*(dexp(f_alpha)+gl%g8*K_p)+gl%g9*K_p*(gl%exp_n-1.D0)/gl%exp_n*(pi_lc+gl%g10)**(-1.D0/gl%exp_n))

    END function

    DOUBLE PRECISION FUNCTION d2gdpdT_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of the Gibbs Free Energy with respect to
    !   pressure p and temperature T in [J / mol] / [Pa K]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    double precision:: df_alpha_dtheta, f_alpha
    double precision:: dK_p_dtheta
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    df_alpha_dtheta = df_alpha_dtheta_calc(gl,theta)
    dK_p_dtheta = dK_p_dtheta_calc(gl,theta)
    f_alpha = f_alpha_calc(gl,theta)


    d2gdpdT_DryIce = gl%R_CO2 / gl%p0_dryice * (gl%g7*(dexp(f_alpha)*df_alpha_dtheta+gl%g8*dK_p_dtheta)+ &
        & gl%g9*dK_p_dtheta*(gl%exp_n-1.D0)/gl%exp_n*(pi_lc+gl%g10)**(-1.D0/gl%exp_n))

    END function

    DOUBLE PRECISION FUNCTION d2gdp2_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of the Gibbs Free Energy with respect to
    !   pressure p at constant temperature T in [J / mol] / [Pa^2]
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    double precision:: K_p
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    K_p = K_p_calc(gl,theta)

    d2gdp2_DryIce = gl%R_CO2 * gl%T0_dryice / gl%p0_dryice**2 * (gl%g9*K_p*(1.D0-gl%exp_n)/gl%exp_n**2*(pi_lc+gl%g10)**(-(gl%exp_n + 1.D0)/gl%exp_n))

    END function

    !
    DOUBLE PRECISION FUNCTION dgdT_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the derivative of the Gibbs Free Energy (chemical potential) of pure Dry ice
    !   with respect to temperature T at constant pressure p in J / mol / K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    double precision:: f_alpha, df_alpha_dtheta
    double precision:: K_p, dK_p_dtheta
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    f_alpha = f_alpha_calc(gl,theta)
    df_alpha_dtheta = df_alpha_dtheta_calc(gl,theta)
    K_p = K_p_calc(gl,theta)
    dK_p_dtheta = dK_p_dtheta_calc(gl,theta)


    dgdT_DryIce = gl%R_CO2 * (gl%g1 + gl%g2 * 2.D0 * del_theta - &
        & 2.D0 * gl%g3 / gl%g4 * (atan(theta/gl%g4) - atan(1.D0/gl%g4)) - &
        & 2.D0 * gl%g5 / gl%g6 * (atan(theta/gl%g6) - atan(1.D0/gl%g6)) + &
        & gl%g7 * del_pi* (dexp(f_alpha)*df_alpha_dtheta + gl%g8*dK_p_dtheta) + &
        & gl%g9 * dK_p_dtheta * ((pi_lc+gl%g10)**((gl%exp_n - 1.D0)/gl%exp_n)-(1.D0+gl%g10)**((gl%exp_n - 1.D0)/gl%exp_n)))


    END function


    DOUBLE PRECISION FUNCTION d2gdT2_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the second derivative of the Gibbs Free Energy (chemical potential) of pure Dry ice
    !   with respect to temperature T at constant pressure p in J / mol / K^2
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    double precision:: f_alpha, df_alpha_dtheta
    double precision:: d2f_alpha_dtheta2
    double precision:: d2K_p_dtheta2
    double precision:: theta, pi_lc, del_theta, del_pi

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    f_alpha = f_alpha_calc(gl,theta)
    df_alpha_dtheta = df_alpha_dtheta_calc(gl,theta)
    d2f_alpha_dtheta2 = d2f_alpha_dtheta2_calc(gl,theta)
    d2K_p_dtheta2 = d2K_p_dtheta2_calc(gl,theta)

    d2gdT2_DryIce = gl%R_CO2 / gl%T0_dryice * (2.D0 * gl%g2 - 2.D0 * gl%g3 / (theta**2 + gl%g4**2) - 2.D0 * gl%g5 / (theta**2 + gl%g6**2) + &
        & gl%g7 * del_pi * (dexp(f_alpha)*d2f_alpha_dtheta2 + dexp(f_alpha)*df_alpha_dtheta**2 + gl%g8* d2K_p_dtheta2) + &
        & gl%g9 * d2K_p_dtheta2 * ((pi_lc+gl%g10)**((gl%exp_n-1.D0)/gl%exp_n)-(1.D0+gl%g10)**((gl%exp_n-1.D0)/gl%exp_n)))

    END function


    !
    DOUBLE PRECISION FUNCTION v_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the molar volume of pure Dry ice
    !   in m³ / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    v_DryIce = dgdp_DryIce(gl,T,p)

    END function


    DOUBLE PRECISION FUNCTION s_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the Entropy of pure Dry ice
    !   in J / mol K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    s_DryIce = -dgdT_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION h_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the Enthalpy of pure Dry ice
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    h_DryIce = g_DryIce(gl,T,p)- T * dgdT_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION u_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the inner energy of pure Dry ice
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    u_DryIce = g_DryIce(gl,T,p)- T * dgdT_DryIce(gl,T,p) - p*dgdp_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION f_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the helmholtz energy of pure Dry ice
    !   in J / mol
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    f_DryIce = g_DryIce(gl,T,p) - p*dgdp_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION cp_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the isobaric heat capacity of pure Dry ice
    !   in J / mol K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    cp_DryIce = -T*d2gdT2_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION alpha_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the thermal expansion coefficient of pure Dry ice
    !   in 1 / K
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    alpha_DryIce = d2gdpdT_DryIce(gl,T,p) / dgdp_DryIce(gl,T,p)

    END function

    DOUBLE PRECISION FUNCTION kappa_DryIce(gl,T,p)
    !----------------------------------------------------------------------------------
    !   Function for calculating the isothermal compressibility of pure Dry ice
    !   in 1 / Pa
    !   Input Parameters:
    !
    !   T   -   Temperature in K
    !   p   -   pressure in MPa
    !
    !----------------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    Double precision:: T, p

    kappa_DryIce = -d2gdp2_DryIce(gl,T,p) / dgdp_DryIce(gl,T,p)

    END function


    !Andreas February 2014
    Double precision function p_DryIce(gl,T,rho,errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the pressure in MPa of dry ice at given T [K] and
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
    parameters%a_p(1) = T
    parameters%a_p(2) = rho

    Delta_allowed = 1.d-8
    Max_Iterations = 30

    !pressure interval
    p_min = 0.D0
    p_max = 500.D0
    p_min_allowed = 0.D0
    p_max_allowed = 500.D0

    call Regula_Falsi(gl,rho_diff_DryIce, p_root, p_min, p_max, Delta_allowed, &
        p_min_allowed, p_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(p_root) < 1.d-10) then
        errval = -18867
    end if

    if ( dabs( (1.D0 / v_DryIce(gl,T,p_root)) - rho) > 1.d-6 ) then
        errval = -18867
    end if

    p_DryIce = p_root

    return

    End function


    double precision function rho_diff_DryIce(gl,p, parameters)
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
    Double Precision :: p   ! Inputvariablen
    Double Precision :: T, rho    !
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    T = parameters%a_p(1)
    rho = parameters%a_p(2)
    rho_diff_DryIce = rho-(1.D0/v_DryIce(gl,T,p))

    return

    End function


    Double precision function T_DryIce_ph(gl,p, h, errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the temperature in K of dry ice at given p [MPa] and
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
    integer:: max_iterations, iterations, errval
        type(type_additional_parameters) :: parameters


    !parameters = 0.d0
    parameters%a_p(1) = p
    parameters%a_p(2) = h

    Delta_allowed = 1.d-4
    Max_Iterations = 30

    !Temperature interval
    T_min = gl%Tmin_Dryice
    T_min_allowed = T_min

    if (p < 0.6D0) then         !Divide the temperature interval appoximately at the triple point
        T_max = 220.D0          !Upper temperature limit of dry ice equation
    else
        T_max = 300.D0
    end if

    T_max_allowed = T_max

    call Regula_Falsi(gl,h_diff_DryIce, T_root, T_min, T_max, Delta_allowed, &
        T_min_allowed, T_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(T_root) < 1.d-10) then
        errval = -18867
    end if

    if (dabs(h_DryIce(gl,T_root, p) - h) > 1.d-2) then
        errval = -18867
    end if

    T_DryIce_ph = T_root

    return

    End function


    double precision function h_diff_DryIce(gl,T, parameters)
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
    Double Precision :: T   ! Inputvariablen
    Double Precision :: h, p      !
        type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    h = parameters%a_p(2)
    h_diff_DryIce = h - h_DryIce(gl,T,p)

    return

    End function


    Double precision function T_DryIce_ps(gl,p, s, errval)
    !----------------------------------------------------------------------------------
    !   Function for calculating the temperature in K of dry ice at given p [MPa] and
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
    integer:: max_iterations, iterations, errval
        type(type_additional_parameters) :: parameters


    !parameters = 0.d0
    parameters%a_p(1) = p
    parameters%a_p(2) = s

    Delta_allowed = 1.d-4
    Max_Iterations = 30

    !Temperature interval
    T_min = gl%Tmin_Dryice
    T_min_allowed = T_min

    if (p < 0.6D0) then         !Divide the temperature interval appoximately at the triple point
        T_max = 220.D0          !Upper temperature limit of dry ice equation
    else
        T_max = 300.D0
    end if

    T_max_allowed = T_max

    call Regula_Falsi(gl,s_diff_DryIce, T_root, T_min, T_max, Delta_allowed, &
        T_min_allowed, T_max_allowed, Max_Iterations, Iterations, errval, parameters)

    if (dabs(T_root) < 1.d-10) then
        errval = -18867
    end if

    if (dabs(s_DryIce(gl,T_root, p) - s) > 1.d-2) then
        errval = -18867
    end if

    T_DryIce_ps = T_root

    return

    End function


    double precision function s_diff_DryIce(gl,T, parameters)
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
    Double Precision :: T   ! Inputvariablen
    Double Precision :: s, p      !
        type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    s = parameters%a_p(2)
    s_diff_DryIce = s - s_DryIce(gl,T,p)

    return

    End function





    DOUBLE PRECISION FUNCTION fug_DryIce(gl,T,p, pos_dryIce)
    !Andreas March 2014
    !----------------------------------------------------------------------------------
    !   Function for calculating the fugacity of pure Dry ice
    !   in MPa
    !   Input Parameters:
    !
    !   T           -   Temperature in K
    !   p           -   pressure in MPa
    !   pos_dryIce  -   position of CO2 in fluid vector
    !
    !----------------------------------------------------------------------------------









    implicit none

    type(type_gl) :: gl


    Double precision:: T, p
    integer:: pos_dryIce

    double precision:: f_alpha
    double precision:: K_p
    double precision:: theta, pi_lc, del_theta, del_pi, pi_sat

    double precision:: ln_phi, ln_phi_psat, phi_sat, psat
    double precision:: d_fluid

    double precision :: AR,D_ARD, rhoredmixorg, tredmixorg
    integer, dimension(nderivs) :: GETDERR
    !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    double precision, dimension(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives

    integer:: errval, iphase

    GETDERR = 0
    GETDERR(1:2) = 1

    ln_phi = 0.D0
    ln_phi_psat = 0.D0

    theta = T / gl%T0_dryice
    pi_lc = p *1.D6 / gl%p0_dryice
    del_theta = theta - 1.D0
    del_pi = pi_lc - 1.D0

    f_alpha = f_alpha_calc(gl,theta)
    K_p = K_p_calc(gl,theta)

    if (T < gl%Tmin_Dryice) then
        errval = -2233
        fug_DryIce = errval
        return
    end if
    if (T > 300.D0) then
        errval = -2231
        fug_DryIce = errval
        return
    end if

    ! Get the melting or sublimation pressure at the given temperature
    ! At the moment the pressure is calculated from the ancillary equation. This could be improved by computing the pressure from the dry ice equation
    if (T > gl%ttp(pos_dryIce)) then
        psat = pmelt_eq(gl,T, pos_dryIce)
        iphase = 1  !For later calculation of the fluid phase density --> liquid
    else
        psat = psub_eq(gl,T, pos_dryIce)
        iphase = 2  !For later calculation of the fluid phase density --> vapor
    end if

    if (psat < 1.D-14) then
        errval = -2225
        fug_DryIce = errval
        return
    end if

    pi_sat = psat *1.D6 / gl%p0_dryice

    !Calculate the fugacity coefficient of PURE CO2
    !----

    phi_sat = 0.D0         !Initialize
    errval = 0

    rhoredmixorg = gl%rhoredmix
    tredmixorg = gl%tredmix

    !This is necessary, since elsewise tau and del get reduced with the wrong parameters!!
    gl%rhoredmix = gl%rhored(pos_dryIce)
    gl%tredmix = gl%tred(pos_dryIce)

    d_fluid = rhomix_calc(gl,T, psat, 0.D0, IPHASE, pos_dryIce)
    if (d_fluid < 1.D-14) then
        errval = -8888
        fug_DryIce = errval
        gl%tredmix = tredmixorg
        gl%rhoredmix = rhoredmixorg
        return
    end if
    CALL FNRDERIVS(gl,T,d_fluid,GETDERR,FNRDER,pos_dryIce)   !subroutine calculates derivatives of residual part

    AR = FNRDER(1)
    D_ARD = FNRDER(2)

    ! safety caution, exp(710) would crash
    if ((AR + D_ARD) > 700.d0) then
        errval = -1444
        fug_DryIce = errval
        gl%tredmix = tredmixorg
        gl%rhoredmix = rhoredmixorg
        return
    else if (isnan(AR)) then   !((AR + D_ARD) == NaN) then
        errval = -9876
        fug_DryIce = errval
        gl%tredmix = tredmixorg
        gl%rhoredmix = rhoredmixorg
        return
    else
        phi_sat = dEXP(AR + D_ARD) / (1.D0 + D_ARD)
    end if

    gl%tredmix = tredmixorg
    gl%rhoredmix = rhoredmixorg
    !----

    !Old version
    !ln_phi_psat = dlog(phi_sat)
    !ln_phi = (   ln_phi_psat + R_CO2 * T0 / p0 * ( g7*(dexp(f_alpha)+g8*K_p)*(p-psat)+g9*K_p*p0*(pi_lc+g10)**((exp_n-1.D0)/exp_n) )  ) / R_CO2 / T
    !fug_DryIce = dexp(ln_phi) * p

    !New version, Andreas October 2014
    fug_DryIce = phi_sat*psat*dexp(gl%T0_dryice / gl%p0_dryice / T * ( gl%g7*(dexp(f_alpha)+gl%g8*K_p)*(p-psat)*1.D6+gl%g9*K_p*gl%p0*( (pi_lc+gl%g10)**((gl%exp_n-1.D0)/gl%exp_n)-(pi_sat+gl%g10)**((gl%exp_n-1.D0)/gl%exp_n) ) )  )

    END function




    !******************************************************************************
    subroutine dfug_DryIce_dT(gl,T, p, pos_dryIce, dfugdT, errval)
    !******************************************************************************
    !
    !  Temperature-derivative of the dry ice fugacity [Pa / T]
    !  Andreas, March 2014




    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision :: t, p
    integer:: pos_dryIce
    ! Output arguments
    double precision:: dfugdT

    double precision:: T1, T2, fug_1, fug_2
    integer:: errval
    !//////////////////////////////////////////////////////////////////

    T1 = T*(1.d0+1.d-4)
    T2 = T*(1.d0-1.d-4)

    fug_2 = fug_DryIce(gl,T2, p, pos_DryIce) *1.D6 !Fugacity in Pa
    if (fug_2 < -1.D12) then                         !If an error occurs: return
        errval = -7779
        return
    end if
    fug_1 = fug_DryIce(gl,T1, p, pos_DryIce) *1.D6 !Fugacity in Pa
    if (fug_1 < -1.D12) then                         !If an error occurs: return
        errval = -7779
        return
    end if

    dfugdT = (fug_2 - fug_1) / (T2 - T1)        ! Pa/K

    end subroutine
    !******************************************************************************
    !******************************************************************************


    end module dryice_module
