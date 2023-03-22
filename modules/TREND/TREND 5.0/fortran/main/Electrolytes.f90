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

    !file contianing the calculations for Electrolytes with fundamental EOS
    module electrolytes

    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use rhomix_pt_module
    use setup_module
    use transport_module
    !use seawater_module
    !use interface_support_module


    contains


    !----------------------------------------------------------------------
    !Subroutine for the calculation of gE for Electrolytes, based on Pitzer's Equations. E.g. Rogers and Pitzer 1982, Rowland 2013,
    !subroutine gE_Pitzer(gl, temp, dw, p, gE)
    double precision function gE_pitzer(gl, temp, dw, p)
    !-------------------------------------------------------------------------------
    implicit none

    type(type_gl) :: gl
    !type (electrolytes_vars) :: el

    integer :: i, salt

    double precision :: dens_sol, Ms, de, Ax, Ion, temp, p, mo, ge, Cmx, pred, pred2, tred, tred2, beta0, beta1, beta2, R, dw, press, h, m, B, C, t, v, lny, osmo!, f_param

    !call electro_init(gl%el)
    salt = 1

    m = gl%el%molality(1)

    press = p
    t = temp
    v = gl%el%vx + gl%el%vm

    Ax =Deybe_hueckel(gl, t, dw, p)
    Ion = ionic_Strength(gl)

    h =  dlog(1.d0 + gl%el%b(gl%el%salt)*Ion**0.5d0)

    B = B_mx(gl, temp, press, Ion)
    C = C_mx(gl, temp, press, gl%el%zm, gl%el%zx)

    gE_pitzer = -Ax * (4.d0*Ion / gl%el%b(gl%el%salt)) * h + 2.d0 * gl%el%vm*gl%el%vx * ( m**2.d0 * B + m**3.d0 *gl%el%vm*gl%el%zm*C )
    lny = ln_el_activity(gl, t, dw, p)
    osmo = osmotic_coeff(gl, t, dw, p)
    gE_pitzer =  v*m*(lny+1.d0-osmo)
   
    end function gE_pitzer! for Electrolytes
    !--------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    !Function for the parameter B in Pitzers eqs (eq6, PPB84)
    double precision function dB_mx_dn(gl, t, p, I)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,I,a, beta0, beta1

    a = gl%el%alpha1
    beta0 = f_param(gl,1,222, t, p)
    beta1 = f_param(gl,2,222, t, p)

    dB_mx_dn = beta0 +  beta1 * dexp(- a * I**0.5d0)

    end function dB_mx_dn
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------



    !------------------------------------------------------------------------------
    !Function for the parameter B in Pitzers eqs (eq6, PPB84)
    double precision function B_mx(gl, t, p, I)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,I,a, beta0, beta1

    a = gl%el%alpha1
    beta0 = f_param(gl,1,222, t, p)
    beta1 = f_param(gl,2,222, t, p)

    B_mx = beta0 + 2.d0* beta1 *(1.d0 - (1.d0 + a*I**0.5d0)*dexp(-a*I**0.5d0))/(a**2.d0*I)

    end function B_mx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !Function for the parameter C in Pitzers eqs (eq7, PPB84)
    double precision function C_mx(gl, t, p, zm, zx)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,zm, zx

    C_mx = f_param(gl,4,222, t, p)/(2.d0*(dabs(zm*zx))**0.5d0)

    end function C_mx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    !Function for the parameter BV in Pitzers eqs
    double precision function BVmx(gl, t, p, I)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,I,a, beta0

    a = gl%el%alpha1
    beta0 = d_f_param_dp(gl,t,p,1)
    !beta1 = f_param(gl,2,222, t, p)    !no pressure dependency

    BVmx = beta0 !+ 2.d0* beta1 *(1.d0 - (1.d0 + a*I**0.5d0)dexp(-a*I**0.5d0))/(a**2.d0*I)

    end function BVmx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !Function for the parameter CV in Pitzers eqs
    double precision function CVmx(gl, t, p, zm, zx)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,zm, zx

    CVmx = d_f_param_dp(gl,t,p,4)/(2.d0*(dabs(zm*zx))**0.5d0)

    end function CVmx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !Function for the parameter BL in Pitzers eqs (eq. 16, PPB84)
    double precision function BLmx(gl, t, p, I)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,I,a, beta0, beta1

    a = gl%el%alpha1
    beta0 = df_param_dT(gl,1, t,p)
    beta1 = df_param_dT(gl,2, t,p)     !no pressure dependency

    BLmx = beta0 + 2.d0* beta1 *(1.d0 - (1.d0 + a*I**0.5d0)*dexp(-a*I**0.5d0))/(a**2.d0*I)

    end function BLmx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !Function for the parameter CL in Pitzers eqs (eq. 16, PPB84)
    double precision function CLmx(gl, t, p, zm, zx)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,zm, zx

    CLmx = df_param_dt(gl,4,t,p)/(2.d0*(dabs(zm*zx))**0.5d0)

    end function CLmx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    !Function for the parameter BJ in Pitzers eqs (eq. 22, PPB84)
    double precision function BJ_mx(gl, t, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,I,a, beta0, beta1,g, gr, Ir

    I = Ionic_strength(gl)
    a = gl%el%alpha1
    beta0 = df2_param_dt2(gl,1,t,p)
    beta1 = df2_param_dt2(gl,2,t,p)     !no pressure dependency

    !Ir = gl%el%m_ref

    g = (1.d0 - (1.d0 + a*I**0.5d0)*dexp(-a*I**0.5d0))/(a**2.d0*I)

    BJ_mx = beta0 + beta1 *2.d0* g

    end function BJ_mx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !Function for the parameter CJ in Pitzers eqs (eq. 22, PPB84)
    double precision function CJ_mx(gl, t, p, zm, zx)
    implicit none
    type(type_gl) :: gl
    double precision :: t,p,zm, zx

    CJ_mx = df2_param_dt2(gl,4,t,p)/(2.d0*(dabs(zm*zx))**0.5d0)

    end function CJ_mx
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    double precision function beta0_func(gl, t, p)
    implicit none
    type (type_gl) :: gl
    double precision :: a, temp, press, p, t, beta0, beta1, beta2

    temp = t
    press = p

    select case(gl%el%betaflag)

    case(1)
        !beta0 = ( f_param(gl,1,1,temp, press) +  1.d-2 * f_param(gl,1,2,temp, press)*gl%el%pred + 1.d-4 * f_param(gl,1,2,temp, press)*gl%el%pred2 ) / mo
        beta0_func = f_param(gl,1,222, temp, press)
    case(2) !Book Pitzer
        beta0_func = f_param(gl,1,222, temp, press)

    end select



    end function beta0_func
    !------------------------------------------------------------------------------



    !---------------------------------------------------------------------------
    double precision function beta1_func(gl, t, p)
    implicit none
    type (type_gl) :: gl
    double precision :: a, temp, press, beta0, beta1, beta2, t, p, mo

    temp = t
    press = p

    select case(gl%el%betaflag)

    case(1)
        beta1_func = f_param(gl,2,222, temp, press)
        !beta1_func = ( f_param(gl,2,1,temp, press) +  1.d-2 * f_param(gl,2,2,temp, press)*gl%el%pred + 1.d-4 * f_param(gl,2,3,temp, press)*gl%el%pred2 ) / mo
    case(2) !Book Pitzer
        beta1_func = f_param(gl,2,222, temp, press)
    end select


    end function beta1_func
    !------------------------------------------------------------------------------
    !******************************************************************************



    !---------------------------------------------------------------------------------
    !function for d gE dp´, first pressure derivative
    double precision function dgE_dp(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, pp, pm, dep, dem, ax, axp, axm, Ms, dax_dp, d_ln_de_dp, compt, dpdd, numdiff, dwm, dwp, deybe_part, virial_part, b_test, b_part, c_part
    double precision :: b0p, b0m, b1m, b1p, cmxp, cmxm, ionic_strenght, ion_strength, db0_dp, db1_dp, cmx, beta2, dcmx_dp, de, d, beta0_ref, dap, da_dp, h!, d_f_param_dp
    integer :: i,j

    dax_dp = av(gl, t, dw, p)  !dA/dp /(RT)


    ion_strength = ionic_strength(gl)

    db0_dp = BVmx(gl, t, p, ion_strength)

    dcmx_dp = CVmx(gl, t, p, gl%el%zm, gl%el%zx)


    h = (1.d0 /( 2.d0 * gl%el%b(1))) *dlog( 1.d0 + gl%el%b(1)*(Ion_strength**(0.5d0)) )


    deybe_part = ( gl%el%vm + gl%el%vx ) * dabs( gl%el%zm * gl%el%zx ) * (dAx_dp) * h

    b_part = 2.d0 * gl%el%vm * gl%el%vx * gl%req(gl%el%solpos) * T*  gl%el%molality(1) *(db0_dp )
    c_part = 2.d0 * gl%el%vm * gl%el%vx * gl%req(gl%el%solpos) * T* gl%el%molality(1)**2.d0 * gl%el%vm * gl%el%zm * dCmx_dp

    dge_dp =  deybe_part + b_part + c_part

    dgE_dp = dgE_dp

    end function dgE_dp
    !-------------------------------------------------------------------------------
    !*****************************************************************************************


    !-------------------------------------------------------------------------------
    !Derivative of Pitzer Excess part wrt tempemerature
    double precision function dgE_dt(gl, t, dw, p)
    !----------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t,p, dw, cp, Ax, axp, axm,zmzx, dwp, dwm, dd_dt, de, dep,dem, da_dt, de_dt, v, h, tp, tm, ge, L_H, ion_strength, get
    double precision :: d_b1_dt, d_B0_dt, d_cmx_dt, db_dt, beta0, beta1, Cmx, part1, part2, gem, gep
    integer :: i

    L_H = h_Ex(gl, t, dw, p)*gl%el%molality(1)!multiplied by moles of salt, as reffered to moles of salt

    ge = gE_pitzer(gl, t, dw,p)!reffered to moles of water

    dge_dt = ( L_H - ge)/gl%el%molality(gl%el%salt)
    !gives S/R reffered to mole of salt (molar prop)
    end function dgE_dt
    !-----------------------------------------------------------------------------------
    !*********************************************************************************************************

    !-------------------------------------------------------------------------------
    !Derivative of Pitzer Excess part wrt tempemerature
    double precision function dse_dn(gl, t, dw, p)
    !----------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t,p, dw, dh, dg
    integer :: i

    dh = dh_ex_dn(gl, t,dw, p)!/(gl%req(1)*T)
    dg = dge_dn(gl, t, dw, p)

    dse_dn = (dh -dg) / T

    end function dse_dn
    !-----------------------------------------------------------------------------------
    !*********************************************************************************************************


    !-------------------------------------------------------------------------------
    !Derivative of Pitzer Excess part wrt tempemerature
    double precision function dgE_dt_num(gl, t, dw, p)
    !----------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t,p, dw, cp, Ax, axp, axm,zmzx, dwp, dwm, dd_dt, de, dep,dem, da_dt, de_dt, v, h, tp, tm, ge, L_H, ion_strength, get
    double precision :: d_b1_dt, d_B0_dt, d_cmx_dt, db_dt, beta0, beta1, Cmx, part1, part2, gem, gep
    integer :: i

    tp = t*(1.d0+1.d-6)
    dwp = rhomix_calc(gl, tp, p, 0.d0, 0, 1)

    gep = gE_pitzer(gl, tp, dwp,p) * Tp !* gl%req(1)

    tm = t*(1.d0-1.d-6)
    dwm = rhomix_calc(gl, tm, p, 0.d0, 0, 1)

    gem = gE_pitzer(gl, tm, dwm,p)* Tm!* gl%req(1)

    dge_dt_num = -(gep-gem)/(tp-tm)

    end function dgE_dt_num
    !-----------------------------------------------------------------------------------
    !*********************************************************************************************************


    !-------------------------------------------------------------------------------
    !function for the enthalpy of the brine, excess part only
    double precision function dh_ex_dn(gl, t, dw, p)
    !-------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, da_dt, v, d_b0_dt, d_b1_dt, d_cmx_dt, beta0, h, I, db_dt, zmzx, R

    R = gl%req(1)

    da_dt =  ah(gl, t, dw,p)*0.5d0

    v = gl%el%vx + gl%el%vm

    I= Ionic_strength(gl)

    d_b0_dt = df_param_dt(gl,1,t, p)
    d_b1_dt = df_param_dt(gl,2, t, p)
    d_Cmx_Dt = df_param_dt(gl, 4, t, p) / 2.d0

    h = I**1.5d0/(1.d0+gl%el%b(1)*I**0.5d0)

    dB_dt = (d_b0_dt +  d_b1_dt * dexp(- gl%el%alpha1 * I**0.5d0))

    zmzx = dabs(gl%el%zm * gl%el%zx)

    dh_ex_dn =  zmzx * da_dt * h - 2.d0 * gl%el%vm * gl%el%vx* T* ( gl%el%molality(gl%el%salt)**2.d0 * dB_dT + 2.d0 * gl%el%molality(gl%el%salt)**3.d0 * (gl%el%vm * gl%el%zm) * d_Cmx_Dt )
    dh_ex_dn = dh_ex_dn * R*T

    end function dh_ex_dn
    !------------------------------------------------------------------------------
    !******************************************************************************************




    !-------------------------------------------------------------------------------
    !function for the enthalpy of the brine, excess part only
    double precision function h_ex(gl, t, dw, p)
    !-------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, da_dt, v, d_b0_dt, d_b1_dt, d_cmx_dt, beta0, h, ion_strength, db_dt, zmzx

    da_dt = ah(gl, t, dw,p)

    v = gl%el%vx + gl%el%vm

    d_b0_dt = df_param_dt(gl,1,t, p)
    d_b1_dt = df_param_dt(gl,2, t, p)
    d_Cmx_Dt = df_param_dt(gl, 4, t, p) / 2.d0

    beta0 = f_param(gl, 1, 222, t, p)

    h = dlog(1.d0 + gl%el%b(gl%el%salt) * Ionic_strength(gl)**0.5d0 ) / ( 2.d0 * gl%el%b(gl%el%salt) )

    ion_strength = Ionic_strength(gl)

    dB_dt = (d_b0_dt +  d_b1_dt * g(gl, gl%el%alpha1, Ion_strength**0.5d0))

    zmzx = dabs(gl%el%zm * gl%el%zx)

    h_ex = v * zmzx * da_dt * h - 2.d0 * gl%el%vm * gl%el%vx  * t* ( gl%el%molality(gl%el%salt) * dB_dT + gl%el%molality(gl%el%salt)**2.d0 * (gl%el%vm * gl%el%zm) * d_Cmx_Dt )


    end function h_ex
    !------------------------------------------------------------------------------
    !******************************************************************************************




    !-------------------------------------------------------------------------------------------------
    !Function for the second derivative of the Gibbs excess energy for calculation of Cpex (cp = -T*d2g/dT2)
    !Result is given back as reduced prop (dimensionless)
    double precision function d2ge_dt2(gl, t, dw, p)
    !--------------------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cex, v, Bjmx, cjmx, h, db_dt, zmzx, ajf, cp02,hr, dm, dm2!, cp0_salt
    integer :: i

    ajf = Ac_num(gl, t, dw, p)

    h = dlog(1.d0 + gl%el%b(gl%el%salt)*Ionic_strength(gl)**0.5d0)/(2.d0*gl%el%b(gl%el%salt))

    v = gl%el%vm + gl%el%vx

    bjmx = BJ_MX(gl, t, p)
    cjmx = cj_mx(gl, t, p, gl%el%zm, gl%el%zx)

    dm = gl%el%molality(1) !- gl%el%m_ref
    dm2 = gl%el%molality(1)**2.d0 !- gl%el%m_ref**2.d0
    zmzx = dabs(gl%el%zm*gl%el%zx)

    cEx =  v*zmzx * Ajf * h*gl%req(1) - 2.d0*gl%el%vm*gl%el%vx * gl%req(gl%el%solpos) * t**2.d0 *(gl%el%molality(gl%el%salt)*bjmx + &
        & gl%el%molality(gl%el%salt)**2.d0*gl%el%vm*gl%el%zm*cjmx)

    d2ge_dt2 = cex
    !Not devided by R but reffered to moles of salt
    end function d2ge_dt2
    !.------------------------------------------------------------------------------------------------------
    !********************************************************************************************************


    !-------------------------------------------------------------------------------------------------
    !Function for the second pressure derivative of Gibbs excess energy (d2G/dp2) for calculation of isothermal compressibility
    double precision function d2ge_dp2(gl, t, dw, p)
    !--------------------------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p,ak, m, v, d_b0_dpp, d_b1_dpp, d_cmx_dpp, h, ion_strength, db_dpp, zmzx
    integer :: i

    ak = ak_num(gl, t, dw, p)!*1.d-6

    m = gl%el%molality(1)
    v = gl%el%vx + gl%el%vm

    d_b0_dpp = d2f_param_dp2(gl, t, p, 1)!/T
    d_b1_dpp = d2f_param_dp2(gl, t, p, 2)!/T
    d_Cmx_dpp = d2f_param_dp2(gl, t, p, 4) / (2.d0)

    !beta0 = f_param(gl, 1, 222, t, p)

    h = dlog(1.d0 + gl%el%b(gl%el%salt) * Ionic_strength(gl)**0.5d0 ) / ( 2.d0 * gl%el%b(gl%el%salt) )

    ion_strength = Ionic_strength(gl)

    dB_dpp = (d_b0_dpp +  d_b1_dpp * g(gl, gl%el%alpha1, Ion_strength**0.5d0))

    zmzx = dabs(gl%el%zm * gl%el%zx)

    d2ge_dp2 = v*zmzx*ak * h + 2.d0*  gl%el%vm * gl%el%vx  * t* gl%req(1) * ( gl%el%molality(1)**2.d0 * db_dpp + gl%el%molality(1)**3.d0 * d_Cmx_dpp )


    end function d2ge_dp2
    !.------------------------------------------------------------------------------------------------------
    !*********************************************************************************************

    !-------------------------------------------------------------------------------------------------
    !Function for the second (mixed) derivative of the Gibbs excess energy used for the isobaric expansion
    double precision function d2ge_dtdp(gl, t, dw, p)

    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, v, d2a_dtdp, d_b0_dtdp, d_b1_dtdp, d_Cmx_Dtdp, h, ion_strength, dB_dtdp, zmzx, dhe_dp, v_ex, b0_tp, b0_p
    integer :: i


    !v_ex = dge_dp(gl, t, dw, p)

    d2a_dtdp = deybe_dtdp(gl, t, dw,p)!cm^3*kg^0.5/(mol^1.5*K)

    v = gl%el%vx + gl%el%vm
    b0_tp = d2f_param_dtdp(gl,1,t, p)
    b0_p = d_f_param_dp(gl, t, p, 1)/T
    d_b0_dtdp = d2f_param_dtdp(gl,1,t, p) + d_f_param_dp(gl, t, p, 1)/T  !kg/(mol*MPa*K)

    d_Cmx_Dtdp = d2f_param_dtdp(gl, 4, t, p) / 2.d0 + d_f_param_dp(gl, t, p, 4) / (2.d0*T)!kg/(mol*MPa*K)

    !beta0 = f_param(gl, 1, 222, t, p)

    h = dlog(1.d0 + gl%el%b(gl%el%salt) * Ionic_strength(gl)**0.5d0 ) / ( 2.d0 * gl%el%b(gl%el%salt) )

    ion_strength = Ionic_strength(gl)

    dB_dtdp = (d_b0_dtdp )!+  d_b1_dtdp * g(gl, gl%el%alpha1, Ion_strength**0.5d0))

    zmzx = dabs(gl%el%zm * gl%el%zx)

    d2ge_dtdp = v*zmzx*d2a_dtdp * h + 2.d0*  gl%el%vm * gl%el%vx  * t* gl%req(1) * ( gl%el%molality(1)**2.d0 * db_dtdp + gl%el%molality(1)**3.d0 * d_Cmx_Dtdp )
    !all together is cm^3/(mol*K)
    end function d2ge_dtdp
    !.------------------------------------------------------------------------------------------------------





    !--------------------------------------------------------------------------------------------------
    !calculates factor dBdt (also named BL) for Piters Eq, e.g. 16 PPB 84
    double precision function dB_dt(gl, t, p)
    !------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: d_b0_dt, t, p, d_b1_dt, d_cmx_dt, ion_strength, h, beta0

    d_b0_dt = df_param_dt(gl,1,t, p)
    d_b1_dt = df_param_dt(gl,2, t, p)
    !d_Cmx_Dt = df_param_dt(gl, 4, t, p) / 2.d0

    beta0 = f_param(gl, 1, 222, t, p)

    h = dlog(1.d0 + gl%el%b(gl%el%salt) * Ionic_strength(gl)**0.5d0 ) / ( 2.d0 * gl%el%b(gl%el%salt) )

    ion_strength = Ionic_strength(gl)

    dB_dt = (d_b0_dt +  d_b1_dt * g(gl, gl%el%alpha1, Ion_strength**0.5d0))

    end function dB_dt
    !-------------------------------------------------------------------------------------------------
    !************************************************************************************************


    !-------------------------------------------------------------------------------
    !Function for the Gibbs energy of the derived by pressure "pure salt", Pitzer 1984, eq.41
    !Volume of the "pure" salt at reference molality; corrected with excess function to real concentration
    double precision function V0_2(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, s0_water, gex_ref, dw, addpart_dp, H0_2v, eos_ref_adjust, m_save, vex_ref, v0_Water, M_salt, M_w,  n_water, Vm_ref, m_ref, n_ref, n_salt
    integer :: i

    M_w = gl%wm(1)
    M_salt = 0.05844277d0
    v0_water = 1/dw                     !molar volume of water(mol/m^3)
    n_water = 1/M_w                     !moles per kg of water
    n_salt = gl%el%molality(1)/M_salt
    Vm_ref = V_ref(gl, t, dw, p) * 1.d-6
    !m_ref = 5.550825d0
    m_ref = gl%el%m_ref
    n_ref = m_ref /M_salt

    m_save = gl%el%molality(1)
    gl%el%molality(1) = m_ref
    vex_ref  = dge_dp(gl, t, dw, p)  !/m_ref
    gl%el%molality(1) = m_save

    V0_2 = Vm_ref / m_ref - (n_water / m_ref) * v0_water - vex_ref

    end function V0_2
    !------------------------------------------------------------------
    !****************************************************************************************************************


    !-------------------------------------------------------------------------------
    !Function for the Gibbs energy of the derived by pressure "pure salt", RP82, eq.23
    !Volume of the "pure" salt also known as V_02, at reference molality, correction of conentration in other routines
    double precision function V_ref(gl, t, dw, p)
    !-------------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, s0_water, gex_ref, dw, addpart_dp, H0_2v, eos_ref_adjust, m_save, vex_ref, v0_Water, m_ref, p0
    integer :: i
    gl%seacalc = .false.

    m_ref = gl%el%m_ref
    p0 = 0.1d0

    v_ref = gl%el%g0_param(1,1) + gl%el%g0_param(1,2) * T + gl%el%g0_param(1,3)*T**2.d0  + gl%el%g0_param(1,4)*t**3.d0 + gl%el%g0_param(1,5)*(p-p0) &
        & + gl%el%g0_param(1,6)*(p-p0)*T + gl%el%g0_param(1,7)*(p-p0)*t**2.d0 + gl%el%g0_param(1,8)*(p**2.d0-p0**2.d0) + gl%el%g0_param(1,9)*(p**2.d0-p0**2.d0)*T

    end function V_ref
    !------------------------------------------------------------------
    !********************************************************************************


    !-------------------------------------------------------------------------------
    !Function for the Gibbs energy of the derived by pressure "pure salt", RP82, eq.23
    !Volume of the "pure" salt also known as V_02, at reference molality, correction of conentration in other routines
    double precision function dV_ref_dt(gl, t, dw, p)
    !-------------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, dw, p0
    integer :: i
    ! gl%seacalc = .false.

    p0 = 0.1d0

    dV_ref_dt = gl%el%g0_param(1,2) + gl%el%g0_param(1,3)*T*2.d0  + 3.d0*gl%el%g0_param(1,4)*t**2.d0 +  &
        &  gl%el%g0_param(1,6)*(p-p0) + gl%el%g0_param(1,7)*(p-p0)*t*2.d0 + gl%el%g0_param(1,9)*(p**2.d0-p0**2.d0)
    
    !dV_ref_dt = gl%el%g0_param(1,2) + gl%el%g0_param(1,3)*T*2.d0  + 2.d0*gl%el%g0_param(1,4)*t**2.d0 +  &
    !    &  gl%el%g0_param(1,6)*(p) + gl%el%g0_param(1,7)*(p)*t*2.d0 + gl%el%g0_param(1,9)*(p**2.d0)
    !
    !    v_ref =  gl%el%g0_param(1,2)  + gl%el%g0_param(1,3)*T*2.d0  + gl%el%g0_param(1,4)*3.d0*t**2.d0  &
    !    & + gl%el%g0_param(1,6)*(p-p0) + gl%el%g0_param(1,7)*(p-p0)*t*2.d0  + gl%el%g0_param(1,9)*(p**2.d0-p0**2.d0)

    end function dV_ref_dt
    !------------------------------------------------------------------
    !********************************************************************************


    !-------------------------------------------------------------------------------
    !Function for the Gibbs energy of the derived by pressure "pure salt", Pitzer 1984, eq.41
    !Volume of the "pure" salt at reference molality; corrected with excess function to real concentration
    double precision function dV0_2_dt(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, s0_water, gex_ref, dw, dVm_ref_dt, volexp, m_save, d2gw_dtdp, v0_Water, M_salt, M_w,  n_water, dvex_ref_dt, m_ref, n_ref, n_salt
    integer :: i

    M_w = gl%wm(1)
    M_salt = 0.05844277d0
    v0_water = 1/(dw*M_w)                    !molar volume of water(mol/m^3)
    n_water = 1/M_w                     !moles per kg of water
    dVm_ref_dt = dV_ref_dt(gl, t, dw, p) * 1.d-6
    !m_ref = 5.550825d0
    m_ref = gl%el%m_ref
    !n_ref = m_ref /M_salt

    volexp = VOLEXP_CALC(gl,T,Dw, 1)
    d2gw_dtdp = volexp * v0_water

    m_save = gl%el%molality(1)
    gl%el%molality(1) = m_ref
    dvex_ref_dt  = d2ge_dtdp(gl, t, dw, p) * 1.d-6 /m_ref
    gl%el%molality(1) = m_save

    dV0_2_dt = dVm_ref_dt - dvex_ref_dt - (d2gw_dtdp / m_ref)

    end function dV0_2_dt
    !------------------------------------------------------------------
    !****************************************************************************************************************




    !-----------------------------------------------------------------------------------------------
    !Function for the apparent volume of brines, Rogers & Pitzer 1982 eqs.
    double precision function v_app(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, vw, vex_ref, v_ex, M_s, M_w , m_ref, m, v02

    M_s = 0.05844277d0
    M_w = gl%wm(1)
    vw = 1.d0 / (dw*m_w)
    m_ref = gl%el%m_ref
    m = gl%el%molality(1)
    gl%el%molality(1) = m_ref

    vex_ref = dge_dp(gl, t, dw, p) *1.d-6!/ m_ref
    gl%el%molality(1) = m
    v_ex = dge_dp(gl, t, dw, p) *1.d-6!/ m

    V02 =  V_ref(gl, t, dw, p)*1.d-6

    v_app = (- (vw / m_ref)  + v_ex - vex_ref + v02)
    ! apparent molar volume of the brine
    end function v_app
    !-------------------------------------------------------------------------------------------------
    !*************************************************************************************************

    !-------------------------------------------------------------------------------
    !Function for the Gibbs energy of the derived by pressure "pure salt", RP82, eq.23
    !Volume of the "pure" salt also known as V_02, at reference molality, correction of conentration in other routines
    double precision function dV_ref_dp(gl, t, dw, p)
    !-------------------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    double precision :: t, p, dw, vex_ref, p0
    integer :: i
    !gl%seacalc = .false.

    !m_ref = gl%el%m_ref
    p0 = 0.1d0

    vex_ref = gl%el%g0_param(1,5)*(1.d0-p0) + gl%el%g0_param(1,6)*(1.d0-p0)*T + gl%el%g0_param(1,7)*(1.d0-p0)*t**2.d0 + gl%el%g0_param(1,8)*(2.d0*p-p0**2.d0) + gl%el%g0_param(1,9)*(2.d0*p - p0**2.d0)*T
    !vex_ref = gl%el%g0_param(1,5)*(1.d0) + gl%el%g0_param(1,6)*(1.d0)*T + gl%el%g0_param(1,7)*(1.d0)*t**2.d0 + gl%el%g0_param(1,8)*(2.d0*p) + gl%el%g0_param(1,9)*(2.d0*p)*T
    dv_ref_dp = vex_ref !/m_ref

    end function dV_ref_dp
    !------------------------------------------------------------------
    !********************************************************************************




    !-----------------------------------------------------------------------------------------------
    !Function for the apparent volume of brines, Rogers & Pitzer 1982 eqs.
    double precision function dv_app_dp(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, vw, vex_ref, v_ex, M_s, M_w , m_ref, m, v02, d2gw_dp2, comp

    M_s = 0.05844277d0
    M_w = gl%wm(1)
    vw = 1.d0 / (dw*m_w)
    comp = COMPT_CALC(gl,T,Dw, 1)
    d2gw_dp2 = - comp * vw*1.d-6

    m_ref = gl%el%m_ref
    m = gl%el%molality(1)
    gl%el%molality(1) = m_ref
    vex_ref = d2ge_dp2(gl, t, dw, p) *1.d-12/ m_ref
    gl%el%molality(1) = m
    v_ex = d2ge_dp2(gl, t, dw, p) *1.d-12 / m
    V02 =  dV_ref_dp(gl, t, dw, p)*1.d-12
    !V02 =- (d2gw_dp2 / m_ref) - vex_ref + v02
    dv_app_dp = (- (d2gw_dp2 / m_ref)  + v_ex -vex_ref + v02)


    end function dv_app_dp
    !-------------------------------------------------------------------------------------------------
    !*************************************************************************************************

    !-----------------------------------------------------------------------------------------------
    !Function for the apparent volume of brines, Rogers & Pitzer 1982 eqs.
    double precision function dv_app_dp_num(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, vw, dwp, dwm, pp, pm, vapp_p, vapp_m

    pp = p*(1.d0+1.d-5)
    pm = p*(1.d0-1.d-5)

    dwp = rhomix_Calc(gl, t, pp, 0.d0, 1, gl%el%solpos)
    dwm = rhomix_Calc(gl, t, pm, 0.d0, 1, gl%el%solpos)

    vapp_p =v_app(gl, t, dwp, pp)
    vapp_m =v_app(gl, t, dwm, pm)

    dv_app_dp_num = (vapp_p - vapp_m)/(pp-pm) *1.d6

    end function dv_app_dp_num
    !-------------------------------------------------------------------------------------------------
    !*************************************************************************************************



    !-----------------------------------------------------------------------------------------------
    !Function for the apparent volume of brines, Rogers & Pitzer 1982 eqs.
    double precision function d2_brine_dp2(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, vw, vex_ref, v_ex, M_s, M_w , m, d2g_dp2_app, d2gw_dp2, comp

    M_s = 0.05844277d0
    M_w = gl%wm(1)
    m = gl%el%molality(1)

    vw = 1.d0 / (dw*m_w)
    comp = COMPT_CALC(gl,T,Dw, 1)
    d2gw_dp2 = - comp * vw*1.d-6 !MPa to Pa!!

    d2g_dp2_app =  dv_app_dp(gl, t, dw, p)!*1.d-6
    d2_brine_dp2 = ( d2gw_dp2 + m * d2g_dp2_app) / (1.d0 + m*M_s)
    d2_brine_dp2 = (d2gw_dp2 + m * d2g_dp2_app) / (1.d0/M_w + m)

    end function d2_brine_dp2
    !-------------------------------------------------------------------------------------------------
    !*************************************************************************************************



    !-----------------------------------------------------------------------------------------------
    !Function for the apparent volume of brines, Rogers & Pitzer 1982 eqs.
    double precision function d2_brine_dtdp(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, vw, dv02dt, M_s, M_w, m, d2gw_dtdp, d2g_dtdp_app, dge_dtdp, vol

    M_s = 0.05844277d0
    M_w = gl%wm(1)
    m = gl%el%molality(1)

    vw = 1.d0 / (dw*m_w)
    vol = volexp_calc(gl,T,Dw, 1)
    d2gw_dtdp = vol * vw

    dge_dtdp = d2ge_dtdp(gl, t, dw, p)*1.d-6/m

    dv02dt = dV0_2_dt(gl, t, dw, p)!*1.d-6

    d2g_dtdp_app =  m *( dv02dt + dge_dtdp)
    d2_brine_dtdp = ( d2gw_dtdp + d2g_dtdp_app ) / (1.d0 + m*M_s)
    d2_brine_dtdp = ( d2gw_dtdp + d2g_dtdp_app ) / (1.d0/M_w + m)

    end function d2_brine_dtdp
    !-------------------------------------------------------------------------------------------------
    !*************************************************************************************************

    !*************************************************************************************************
    double precision function comp_brine(gl, t, dw, p)
    implicit none
    type(type_gl)::gl
    double precision:: t, dw, p, g_pp, g_p

    g_pp = d2_brine_dp2(gl, t, dw, p)

    g_p = dbrine_dp(gl, t, dw, p)

    comp_brine = - g_pp/g_p

    end function !comp_brine
    !-----------------------------------------------------------------------------------------------
    !***********************************************************************************************

    !*************************************************************************************************
    double precision function expan_brine(gl, t, dw, p)
    implicit none
    type(type_gl)::gl
    double precision:: t, dw, p, g_tp, g_p

    g_tp = d2_brine_dtdp(gl, t, dw, p)

    g_p = dbrine_dp(gl, t, dw, p)

    expan_brine = g_tp/g_p

    end function !expan_brine
    !-----------------------------------------------------------------------------------------------
    !***********************************************************************************************

    

    !--------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------------------------------------------------------------------------------
    !     FUNCTIONS FOR THE PROPERITES OF THE BRINE
    !--------------------------------------------------------------------------------------
    !-----------------------------------------------------------------
    !Gibbs energy of the brine
    double precision function g_brine(gl, t, dw, p, nw_in, n_salt_in)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, G_salt, n_salt, g_ex, v, nw, g_water, R, mix_part, n2, g02_salt, g02_w, M_2, m, xw, x_salt, g_brinen, g_brinex!, g_exold, mix_old
    double precision, optional :: nw_in, n_Salt_in

    M_2 =  0.05844277d0
    nw = 1.d0 / gl%wm(gl%el%solpos)
    R = gl%req(1)
    v = gl%el%vm + gl%el%vx
    n_salt = gl%el%molality(1)
    m = gl%el%molality(1)

    if(present(nw_in))then
        xw = nw_in
    end if
    if(present(n_salt_in)) then
        x_Salt = n_salt_in
    end if

    g_water = g_calc(gl, t, dw, 1)
    G_salt = G0_2(gl, t, dw, p)
    !dw = rhomix_calc(gl, t, p, 0.d0, 1, gl%el%solpos)
    g_ex = ge_pitzer(gl,t, dw, p)*gl%req(1) * T
    mix_part = R * T * n_salt * v*(1.d0 - dlog(gl%el%molality(1)))
    g02_w =G0_2_split_wex(gl, t, dw, p)
    g02_Salt = G0_2_split_salt(gl, t, dw, p)

    g_brine = (nw*g_water + n_salt * g02_Salt + n_salt*g02_w +  g_ex  -  mix_part)/(nw + n_salt)

    ! molar Gibbs energy of the brine, refferend to moles of brine!!
    end function g_brine
    !------------------------------------------------------------------
    !*****************************************************************


    !-----------------------------------------------------------------
    !Gibbs energy of the salt
    double precision function G0_2(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, G_salt, n_salt, g_ex, v, nw, g_water, R, mol_save, t0, p0, G_salt_ref, S_salt_ref, nr, dw_ref, g_w_ref, s_w_ref, g_w, ge_tpm_ref, se_tpm_ref, ge_m_ref
    double precision :: cp_salt_ref_int, v_salt_ref_int, Mw

    G_salt_ref = gl%el%g0_param(1,10)
    S_salt_ref = gl%el%g0_param(1,11)
    t0 = 298.15d0
    p0 = 0.1d0
    nr = gl%el%m_ref

    Mw = gl%wm(1)
    nw = 1.d0/gl%wm(1)
    !n_salt = gl%el%molality(1)
    !n_salt = gl%el%x_salt

    !call reduced_parameters_calc(gl,T)
    dw_ref =  rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)!55344.5576839165d0
    !dw_ref = 55344.5576839165d0
    g_w_ref = g_calc(gl, t0, dw_ref, 1)
    s_w_ref = s_calc(gl, t0, dw_ref, 1)

    g_w  = g_calc(gl, t, dw, 1)
    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    ge_tpm_ref = ge_pitzer(gl, t0, dw_ref, p0)*gl%req(1)*t0
    se_tpm_ref = dgE_dt(gl, t0, dw_ref, p0)*gl%req(1)!*nr
    ge_m_ref = ge_pitzer(gl, t, dw, p)*gl%req(1)*t

    gl%el%molality(1) = mol_save

    cp_salt_ref_int = t * ( (t**2.d0 - t0**2.d0)*gl%el%cp0_param(3)/6.d0 + (t-t0)*gl%el%cp0_param(2)/2.d0 +gl%el%cp0_param(1)*dlog(t/t0) &
        & - (2.d0 *t0**3.d0*gl%el%cp0_param(3) + 3.d0 * T0**2.d0 *gl%el%cp0_param(2) + 6.d0 * gl%el%cp0_param(1))/(6.d0 * t0) +  &
        &   (2.d0 *t0**3.d0*gl%el%cp0_param(3)+ 3.d0 * T0**2.d0 *gl%el%cp0_param(2) + 6.d0 * gl%el%cp0_param(1))/(6.d0 * t) )


    v_salt_ref_int = (p**3.d0 - p0**3.d0)*(t**2.d0*gl%el%g0_param(1,9)/3.d0 + gl%el%g0_param(1,8)/3.d0 ) + &
        &(p**2.d0 - p0**2.d0)*(t**2.d0*gl%el%g0_param(1,7)/2.d0 + T*gl%el%g0_param(1,6)/2.d0 + gl%el%g0_param(1,5)/2.d0 ) + &
        &(p-p0)* ( t**3.d0*gl%el%g0_param(1,4) - T**2.d0*p0**2.d0*gl%el%g0_param(1,9) - t**2.d0*p0*gl%el%g0_param(1,7) + T**2.d0*gl%el%g0_param(1,3) + &
        &  T*gl%el%g0_param(1,2) - p0**2.d0*gl%el%g0_param(1,8) - p0*gl%el%g0_param(1,5) + gl%el%g0_param(1,1) )

    G0_2 = G_salt_ref + nw / nr *(g_w_ref - g_w) - (ge_tpm_ref - ge_m_ref)/nr - (t-t0)*(S_salt_ref + nw/nr * s_w_ref + se_tpm_ref/nr) &
        & - cp_salt_ref_int + v_salt_ref_int

    end function G0_2
    !------------------------------------------------------------------
    !*****************************************************************


    !-----------------------------------------------------------------
    !Gibbs energy of the salt
    double precision function G0_2_split_salt(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, G_salt, n_salt, g_ex, v, nw, g_water, R, mol_save, t0, p0, G_salt_ref, S_salt_ref, nr, dw_ref, g_w_ref, s_w_ref, g_w, ge_tpm_ref, se_tpm_ref, ge_m_ref
    double precision :: cp_salt_ref_int, v_salt_ref_int, Mw

    G_salt_ref = gl%el%g0_param(1,10)
    S_salt_ref = gl%el%g0_param(1,11)
    t0 = 298.15d0
    p0 = 0.1d0
    nr = gl%el%m_ref

    Mw = gl%wm(1)
    nw = 1.d0/gl%wm(1)
    !n_salt = gl%el%molality(1)
    !n_salt = gl%el%x_salt

    !call reduced_parameters_calc(gl,T)
    dw_ref =  rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)
    dw_ref = 55344.5576839165d0
    g_w_ref = g_calc(gl, t0, dw_ref, 1) 
    s_w_ref = s_calc(gl, t0, dw_ref, 1)

    g_w  = g_calc(gl, t, dw, 1)
    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    ge_tpm_ref = ge_pitzer(gl, t0, dw_ref, p0)*gl%req(1)*t0 
    se_tpm_ref = dgE_dt(gl, t0, dw_ref, p0)*gl%req(1)
    ge_m_ref = ge_pitzer(gl, t, dw, p)*gl%req(1)*t

    gl%el%molality(1) = mol_save

    cp_salt_ref_int = t * ( (t**2.d0 - t0**2.d0)*gl%el%cp0_param(3)/6.d0 + (t-t0)*gl%el%cp0_param(2)/2.d0 +gl%el%cp0_param(1)*dlog(t/t0) &
        & - (2.d0 *t0**3.d0*gl%el%cp0_param(3) + 3.d0 * T0**2.d0 *gl%el%cp0_param(2) + 6.d0 * gl%el%cp0_param(1))/(6.d0 * t0) +  &
        &   (2.d0 *t0**3.d0*gl%el%cp0_param(3)+ 3.d0 * T0**2.d0 *gl%el%cp0_param(2) + 6.d0 * gl%el%cp0_param(1))/(6.d0 * t) )


    v_salt_ref_int = (p**3.d0 - p0**3.d0)*(t**2.d0*gl%el%g0_param(1,9)/3.d0 + gl%el%g0_param(1,8)/3.d0 ) + &
        &(p**2.d0 - p0**2.d0)*(t**2.d0*gl%el%g0_param(1,7)/2.d0 + T*gl%el%g0_param(1,6)/2.d0 + gl%el%g0_param(1,5)/2.d0 ) + &
        &(p-p0)* ( t**3.d0*gl%el%g0_param(1,4) - T**2.d0*p0**2.d0*gl%el%g0_param(1,9) - t**2.d0*p0*gl%el%g0_param(1,7) + T**2.d0*gl%el%g0_param(1,3) + &
        &  T*gl%el%g0_param(1,2) - p0**2.d0*gl%el%g0_param(1,8) - p0*gl%el%g0_param(1,5) + gl%el%g0_param(1,1) )

    G0_2_split_salt = G_salt_ref  - (t-t0)*(S_salt_ref ) - cp_salt_ref_int + v_salt_ref_int

    end function G0_2_split_salt
    !------------------------------------------------------------------
    !*****************************************************************




    !-----------------------------------------------------------------
    !Gibbs energy of the salt
    double precision function G0_2_split_wex(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, G_salt, n_salt, g_ex, v, nw, g_water, R, mol_save, t0, p0, G_salt_ref, S_salt_ref, nr, dw_ref, g_w_ref, s_w_ref, g_w, ge_tpm_ref, se_tpm_ref, ge_m_ref
    double precision :: cp_salt_ref_int, v_salt_ref_int, Mw

    G_salt_ref = gl%el%g0_param(1,10) !
    S_salt_ref = gl%el%g0_param(1,11)
    t0 = 298.15d0
    p0 = 0.1d0
    nr = gl%el%m_ref

    Mw = gl%wm(1)
    nw = 1.d0/gl%wm(1)


    !call reduced_parameters_calc(gl,T)
    dw_ref =  rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)!55344.5576839165d0
    !dw_ref = 55344.5576839165d0
    g_w_ref = g_calc(gl, t0, dw_ref, 1)
    s_w_ref = s_calc(gl, t0, dw_ref, 1)

    g_w  = g_calc(gl, t, dw, 1)
    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    ge_tpm_ref = ge_pitzer(gl, t0, dw_ref, p0)*gl%req(1)*t0 !ge_pitzer(gl, t0, dw_ref, p0)
    se_tpm_ref = dgE_dt(gl, t0, dw_ref, p0)*gl%req(1)!*nr
    ge_m_ref = ge_pitzer(gl, t, dw, p)*gl%req(1)*t

    gl%el%molality(1) = mol_save



    G0_2_split_wex =  nw / nr *(g_w_ref - g_w) - (ge_tpm_ref - ge_m_ref)/nr - (t-t0)*( nw/nr * s_w_ref + se_tpm_ref/nr)


    end function G0_2_split_wex
    !------------------------------------------------------------------
    !*****************************************************************




    !-----------------------------------------------------------------
    !Gibbs energy of the brine
    double precision function s_brine(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, s_salt, n_salt, s_ex, v, nw, s_water, R

    s_salt = S0_2(gl, t, dw, p)
    n_salt = gl%el%molality(1)
    s_ex = dgE_dt(gl,t, dw, p)*gl%req(1)
    v = gl%el%vm + gl%el%vx

    s_water = s_calc(gl, t, dw, 1)

    nw = 1.d0 / gl%wm(gl%el%solpos)
    R = gl%req(1)

    s_brine = (nw * s_water + n_salt * s_salt + s_ex - R * n_salt * v*(1.d0 - dlog(gl%el%molality(1))) )/ (nw + n_salt)

    end function s_brine
    !------------------------------------------------------------------
    !*****************************************************************


    !-----------------------------------------------------------------
    !Gibbs energy of the brine
    double precision function h_brine(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, g, s

    g = g_brine(gl, t, dw, p)
    s = s_brine(gl, t, dw, p)

    h_brine = g + T*s

    end function h_brine
    !------------------------------------------------------------------
    !*****************************************************************



    !-----------------------------------------------------------------------------------------
    !*****************************************************************************************
    !-----------------------------------------------------------------------------------------



    !--------------!--------------!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------!-------------------------------------------------------------------------------
    !--------------!Functions for Deybe-Hueckel and its derivatives for the needed properties
    !--------------!-------------------------------------------------------------------------------
    !--------------!--------------!************************************************************************************

    !-------------------------------------------------------------------------------
    !!function for Debye-Hückel term
    double precision function Deybe_hueckel(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: dens_sol, t, dw, p, Ms, de, ax, d_in

    !dens_sol = rhomix_calc(gl, t, p, 0.d0, 1, gl%el%solpos)!not good for numerical derivatives :) !be aware of redfunc in mixtures
    d_in  = dw
    !dw = 55344.5576839165d0

    dens_sol = dw


    Ms = gl%wm(1)!gl%wm(gl%el%solpos)

    if(.not. allocated(gl%de)) allocate(gl%de)
    if(.not. allocated(gl%el)) allocate(gl%el)

    !call initialize_Electrolytes(gl) !only for fitting preprocessing
    call initialize_dew(gl)

    call dielectric_swap(gl, 1)
    call dielectric_swap(gl, 2)

    de = DE_CALC(gl,T, Dens_sol,1)! gl%el%solpos)

    call dielectric_swap(gl,3)

    !Debye_Hückel Parameter

    Ax = 1.d0/3.d0 * (2.d0 * gl%el%pi * gl%el%Na * dens_sol * Ms )**0.5d0 * (gl%el%e_const**2.d0 / (gl%el%pi * 4.d0 * gl%el%el_permit * de* gl%el%kb * T))**(3.d0/2.d0)

    dw = d_in

    Deybe_hueckel = ax

    end function deybe_hueckel
    !----------------------------------------------------------------------------
    !****************************************************************************


    !----------------------------------------------------------------------------
    !first derivative of deybe-hückel wrt p   !cm^3 kg^0.5 / mol^-1.5
    double precision function Av(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, numdiff, MS, ion_strength, de, ax, dddp, compt, pp, pm, d_ln_de_dp, dwp, dwm, dpdd, dep, dem, dw_in
    integer :: i

    dw_in = dw

    !dw = 55344.5576839165d0

    numdiff = 1.d-7

    !Ms = gl%wm(gl%el%solpos)

    !   Ion_strength = ionic_strength(gl)

    !de = DE_CALC(gl,t, dw, gl%el%solpos)

    Ax = Deybe_hueckel(gl, t, dw, p)

    DPDD=DPDD_CALC(gl,T,Dw,gl%el%solpos)

    compt = 1.d0/Dw/DPDD_CALC(gl,T,Dw, gl%el%solpos)
    compt = compt / 10.d0 !MPa - bar soll so

    dwp = dw*(1.d0 + numdiff)
    pp = p_calc(gl, t, dwp, gl%el%solpos)

    if(.not. allocated(gl%de)) allocate(gl%de)
    call initialize_dew(gl)
    call dielectric_swap(gl, 1)
    call dielectric_swap(gl, 2)

    de = de_Calc(gl,t, dw, gl%el%solpos)

    dep = de_calc(gl, t, dwp, gl%el%solpos)

    dwm = dw * (1.d0-numdiff)

    pm = p_calc(gl, t, dwm, gl%el%solpos)

    dem = de_calc(gl, t, dwm, gl%el%solpos)

    call dielectric_swap(gl,3)

    !dep = dep-dem
    d_ln_de_dp = (dep-dem) / (dwp-dwm) * (1.d0 / DpDD) !ok, checked against paper
    d_ln_de_dp = d_ln_De_dp / 10.d0

    av = 2.d0 * gl%REq(1) * t * ( 3.d0/(de) * d_ln_de_dp - compt) * ax * 10.d0
    dw = dw_in
    end function Av
    !---------------------------------------------------------------------------
    !*************************************************************************************


    !-------------------------------------------------------------------------------------
    !first derivative of Deybe-hückel wrt t/RT
    double precision function Ah(gl, t, dw, p)
    implicit none
    type(Type_gl) :: gl
    double precision :: t, dw, p, tp, tm, dwp, dwm, dep, dem, de_Dt, de, ax, dd_dt, da_dt, ap, am

    if(.not. allocated(gl%de)) allocate(gl%de)
    call initialize_dew(gl)
    call dielectric_swap(gl, 1)
    call dielectric_swap(gl, 2)

    tp = t* (1.d0+1.d-6)
    tm = t*(1.d0-1.d-6)

    dwp = rhomix_Calc(gl, tp, p, 0.d0, 1, gl%el%solpos)
    dwm = rhomix_Calc(gl, tm, p, 0.d0, 1, gl%el%solpos)

    !ap = deybe_hueckel(gl, tp, dwp, p)
    !am = deybe_hueckel(gl, tm, dwm, p)

    dep = de_calc(gl, tp, dwp, gl%el%solpos)

    dem = de_calc(gl, tm, dwm, gl%el%solpos)

    de_dt = (dep - dem ) / (tp - tm)

    de = de_calc(gl, t, dw, gl%el%solpos)

    call dielectric_swap(gl,3)

    ax = deybe_hueckel(gl, t, dw, p)

    dd_dT = DDDT_Calc(gl,t, dw, gl%el%solpos)

    da_dt = - 6.d0 * ax * ( 1.d0 + (( (T * de_dt) / de) - ( (T * dd_dt) / (3.d0 * dw ) ) ))!* gl%req(gl%el%solpos) * t

    !da_dt = (ap -am) / (tp-tm)

    ah=da_dt

    end function ah
    !-----------------------------------------------------------------------------------------
    !******************************************************************************************


    !-------------------------------------------------------------------------------------
    !first derivative of Deybe-hückel wrt t, numerical deriv
    double precision function Ah_num(gl, t, dw, p)
    implicit none
    type(Type_gl) :: gl
    double precision :: t, dw, p, tp, tm, dwp, dwm, dep, dem, de_Dt, de, ax, dd_dt, da_dt, ap, am

    if(.not. allocated(gl%de)) allocate(gl%de)
    call initialize_dew(gl)
    call dielectric_swap(gl, 1)
    call dielectric_swap(gl, 2)

    tp = t* (1.d0+1.d-5)
    tm = t*(1.d0-1.d-5)

    dwp = rhomix_Calc(gl, tp, p, 0.d0, 1, gl%el%solpos)
    dwm = rhomix_Calc(gl, tm, p, 0.d0, 1, gl%el%solpos)

    ap = deybe_hueckel(gl, tp, dwp, p)
    am = deybe_hueckel(gl, tm, dwm, p)

    da_dt = (ap -am) / (tp-tm)

    ah_num=da_dt

    ah_num = da_dt * 4.d0 * t  !acutally Ah/RT

    end function ah_num
    !-----------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------
    !Derivative of Deybe-Hueckel for calculation of heat capacities corresponds to second temp deriv
    double precision function Ac_num(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, dwp, dwm, tp, tm, d2a_dt2, ah, ap, a2p, am, a2m, t2p, t2m, dw2m, dw2p, a, tdiff, da_dt, R

    tdiff = t*1.d-4

    tp = t+tdiff
    tm = t-tdiff
    R = gl%req(1)

    t2p = t+2.d0*tdiff
    t2m  = t-2.d0 * tdiff

    dw = rhomix_calc(gl, t, p, 0.0d0, 1, 0)

    dwp = rhomix_calc(gl, tp, p, 0.0d0, 1, 0)
    dwm = rhomix_calc(gl, tm, p, 0.0d0, 1, 0)
    dw2p = rhomix_calc(gl, t2p, p, 0.0d0, 1, 0)
    dw2m = rhomix_calc(gl, t2m, p, 0.0d0, 1, 0)

    a = deybe_hueckel(gl, t, dw, p)

    ap = deybe_hueckel(gl, tp, dwp, p)
    am = deybe_hueckel(gl, tm, dwm, p)

    a2p = deybe_hueckel(gl, t2p, dw2p, p)
    a2m = deybe_hueckel(gl, t2m, dw2m, p)

    d2a_dt2 = (-a2p + 16.d0 * ap - 30.d0*a + 16.d0 * am - a2m) / (12.d0*tdiff**2.d0 )

    da_dt = (ap-am)/(tp-tm)
    d2a_dt2 = (ap - 2.d0*a + am) / (tdiff**2.d0)
    Ac_num = 8.d0 * gl%req(1) * t * da_dt + 4.d0 * gl%req(1) * t**2.d0 * d2a_dt2


    ac_num = ac_num /gl%req(1)

    end function Ac_num
    !---------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------
    !Derivative of Deybe-Hueckel for calculation of functions with second pressure deriv corresponds to second pressure deriv
    double precision function Ak_num(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, R,dwp, dwm, tp, tm, d2a_dt2, ah, ap, a2p, am, a2m, t2p, t2m, dw2m, dw2p, a, tdiff, da_dt, pp, pm, p2p, p2m, da2_dp2, ac_num, pdiff, d2a_dp2, dadp, rt, d2a_dp2s, numdiff, pdiffp, pdiffm
    double precision :: de, a_v, dep, dem, d_pp, dpdd, de_dp, d2e_dp2, compt, diff
    R=gl%req(1)
    diff = 1.d-2
    !
    pdiff = diff

    pp = p + diff
    pm = p - diff
    p2p = p+2.d0*diff
    p2m = p-2.d0*diff
    dwp = rhomix_calc(gl, t, pp, 0.0d0, 1, 1)
    dwm = rhomix_calc(gl, t, pm, 0.0d0, 1, 1)
    dw2p = rhomix_calc(gl, t, p2p, 0.0d0, 1, 1)
    dw2m = rhomix_calc(gl, t, p2m, 0.0d0, 1, 1)
    a = deybe_hueckel(gl, t, dw, p)!* 4.d0 !*R*T!*4.d0
    ap = deybe_hueckel(gl, t, dwp, pp)!* 4.d0 !*R*T!*4.d0
    am = deybe_hueckel(gl, t, dwm, pm)!* 4.d0 !*R*T!*4.d0
    d2a_dp2s = (ap-am)/(pp-pm)!/4.d0
    d2a_dp2s = (ap - 2.d0*a  + am)/((diff**2.d0))
    a2p = deybe_hueckel(gl, t, dw2p, p2p)!* 4.d0
    a2m = deybe_hueckel(gl, t, dw2m, p2m)!* 4.d0
    rt = gl%req(1) * T*4.d0
    d2a_dp2s = (-a2p + 16.d0 * ap - 30.d0*a + 16.d0 * am - a2m) / (12.d0*diff**2.d0 ) *(-rt)
    ak_num = d2a_dp2s
    
    end function Ak_num
    !---------------------------------------------------------------------------------------
    !***********************************************************************************************

    !--------------------------------------------------------------------------------------
    !function for the second pressure derivative of the Deybe hueckel slope
    double precision function a_vv(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, pp, pm, p2p, p2m, ap, am, a2p, a2m, numdiff, r, add, dwp, dwm, dw2p, dw2m, a, avd

    r = gl%req(gl%el%solpos)

    numdiff = 1.d-6

    add = p*numdiff

    pp = p*(1.d0+numdiff)
    pm = p*(1.d0-numdiff)

    p2p = p*(1.d0+2.d0*numdiff)
    p2m = p*(1.d0+2.d0*numdiff)

    dwp = rhomix_calc(gl, t, pp, 0.d0, 1, gl%el%solpos)
    dwm = rhomix_calc(gl, t, pm, 0.d0, 1, gl%el%solpos)

    dw2p = rhomix_calc(gl, t, p2p, 0.d0, 1, gl%el%solpos)
    dw2m = rhomix_calc(gl, t, p2m, 0.d0, 1, gl%el%solpos)

    a = deybe_hueckel(gl, t, dw, p)*4.d0*R*T!*4.d0*r

    ap = deybe_hueckel(gl, t, dwp, pp)*4.d0*R*T!*4.d0*r
    am = deybe_hueckel(gl, t, dwm, pm)*4.d0*R*T!*4.d0*r

    a2p = deybe_hueckel(gl, t, dw2p, p2p)*4.d0*R*T!*4.d0*r
    a2m = deybe_hueckel(gl, t, dw2m, p2m)*4.d0*R*T

    a_vv = (ap-2.d0*a + am)/(add)**2.d0 !*r*T



    end function
    !------------------------------------------------------------------------------------
    !************************************************************************************


    !--------------------------------------------------------------------------------------
    !function for the second pressure derivative of the Deybe hueckel slope
    double precision function deybe_dtdp(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, tp, tm, dw, p, pp, pm, p2p, p2m, ap, am, a2p, a2m, numdiff, r, add, dwp, dwm, dw2p, dw2m, a, avd
    double precision ::dwp_p, dwm_p, dwp_T, dwm_T, ap_p, am_p, ap_T, am_T

    numdiff = 1.d-4

    Tp = T*1.d0+numdiff
    Tm = T*1.d0-numdiff

    dwp_T = rhomix_calc(gl, tp, p, 0.d0, 1, gl%el%solpos)
    dwm_T = rhomix_calc(gl, tm, p, 0.d0, 1, gl%el%solpos)

    ap_T = av(gl, tp, dwp_T, p)
    am_T = av(gl, tm, dwm_T, p)

    deybe_dtdp = (ap_T - am_T )/(tp-tm)

    end function deybe_dtdp
    !------------------------------------------------------------------------------------
    !************************************************************************************




    !----------!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !------------------------------------------------------------------------------------
    !-----Functions for the calculation and handling of the dielectric constant of pure water, also handels initisalisation and swapping within mixtures
    !----------------------------------------------------------------------------------------------------
    !*****************************************************************************************************************


    !---------------------------------------------------------------------------
    !function for the dielectric constant of water, hardcoded parameters
    double precision function de_water(gl,t, dw)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p

    if(.not. allocated(gl%de)) allocate(gl%de)

    call dielectric_swap(gl, 1)
    call dielectric_swap(gl, 2)

    de_water = DE_CALC(gl,T, dw, gl%el%solpos)

    call dielectric_swap(gl,3)

    end function
    !---------------------------------------------------------------------------
    !***************************************************************************


    !---------------------------------------------------------------------------------
    !function for param derivative
    double precision function d_f_param_dp(gl, temp, press, i)
    implicit none
    type(type_gl) :: gl
    double precision :: temp, press, p, p_red
    integer :: i!, functionflag

    p = press

    !p_red = gl%el%pred

    select case(gl%el%functionflag)

    case(1)
        p=p*10.d0
        d_f_param_dp = gl%el%lit_param(i,3) + 2.d0 * gl%el%lit_param(i,4)*p + 3.d0 * gl%el%lit_param(i,5)*p**2.d0 +  &
            & (  gl%el%lit_param(i,8) + 2.d0 * gl%el%lit_param(i,9)*p + 3.d0 * gl%el%lit_param(i,10)*p**2.d0 )*temp + &
            & (  gl%el%lit_param(i,12) + 2.d0*gl%el%lit_param(i,13)*p)*temp**2.d0 + &
            & ( gl%el%lit_param(i,15) + 2.d0 * gl%el%lit_param(i,16)*p + 3.d0 * gl%el%lit_param(i,17)*p**2.d0 )/(temp - 227.d0) &
            & + ( gl%el%lit_param(i,19)+ 2.d0 * gl%el%lit_param(i,20)*p + 3.d0*gl%el%lit_param(i,21)*p**2.d0 )/(680.d0 -temp)
        p = p/10.d0
    case(2) 


        d_f_param_dp = gl%el%pitzer_param(i,3) + 2.d0 * gl%el%pitzer_param(i,4)*p + 3.d0 * gl%el%pitzer_param(i,5)*p**2.d0 +  &
            & (  gl%el%pitzer_param(i,8) + 2.d0 * gl%el%pitzer_param(i,9)*p + 3.d0 * gl%el%pitzer_param(i,10)*p**2.d0 )*temp + &
            & (  gl%el%pitzer_param(i,12) + 2.d0*gl%el%pitzer_param(i,13)*p)*temp**2.d0 + &
            & ( gl%el%pitzer_param(i,15) + 2.d0 * gl%el%pitzer_param(i,16)*p + 3.d0 * gl%el%pitzer_param(i,17)*p**2.d0 )/(temp - 227.d0) &
            & + ( gl%el%pitzer_param(i,19)+ 2.d0 * gl%el%pitzer_param(i,20)*p + 3.d0*gl%el%pitzer_param(i,21)*p**2.d0 )/(680.d0 -temp)


    end select

    end function d_f_param_dp

    !------------------------------------------------------------------------------------------
    !******************************************************************************************

    !---------------------------------------------------------------------------------
    !function for param derivative
    double precision function d2f_param_dp2(gl, temp, press, i)
    implicit none
    type(type_gl) :: gl
    double precision :: temp, press, p, p_red
    integer :: i!, functionflag

    p = press

    !p_red = gl%el%pred

    select case(gl%el%functionflag)

    case(1)
        p=p*10.d0
        d2f_param_dp2 =  2.d0 * gl%el%lit_param(i,4) + 6.d0 * gl%el%lit_param(i,5)*p +  &
            & ( 2.d0 * gl%el%lit_param(i,9) + 6.d0 * gl%el%lit_param(i,10)*p)*temp + &
            & ( 2.d0*gl%el%lit_param(i,13))*temp**2.d0 + &
            & ( 2.d0 * gl%el%lit_param(i,16) + 6.d0 * gl%el%lit_param(i,17)*p)/(temp - 227.d0) &
            & + ( 2.d0 * gl%el%lit_param(i,20) + 6.d0*gl%el%lit_param(i,21)*p)/(680.d0 -temp)
        p=p/10.d0
    case(2) 


        d2f_param_dp2 =  2.d0 * gl%el%pitzer_param(i,4) + 6.d0 * gl%el%pitzer_param(i,5)*p +  &
            & ( 2.d0 * gl%el%pitzer_param(i,9) + 6.d0 * gl%el%pitzer_param(i,10)*p)*temp + &
            & ( 2.d0*gl%el%pitzer_param(i,13))*temp**2.d0 + &
            & ( 2.d0 * gl%el%pitzer_param(i,16) + 6.d0 * gl%el%pitzer_param(i,17)*p)/(temp - 227.d0) &
            & + ( 2.d0 * gl%el%pitzer_param(i,20) + 6.d0*gl%el%pitzer_param(i,21)*p)/(680.d0 -temp)


    end select

    end function d2f_param_dp2

    !------------------------------------------------------------------------------------------
    !******************************************************************************************



    !--------------------------------------------------------------------------------
    ! Parameterfunction for temperature and pressure dependent parameters (eeq. 3.6 in Rowland 2013)
    double precision function f_param(gl, i, j,temp, p)
    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el
    integer :: i, j
    double precision :: temp, p, p_red

    p_red = gl%el%pred

    select case(gl%el%functionflag)

    case(1) ! lit_params_pitzer
        p = p*10.d0
        f_param = gl%el%lit_param(i,1) / Temp + gl%el%lit_param(i,2) + gl%el%lit_param(i,3)*p + gl%el%lit_param(i,4)*p**2.d0 + gl%el%lit_param(i,5)*p**3.d0 + gl%el%lit_param(i,6)*dlog(temp) &
            & +( gl%el%lit_param(i,7) + gl%el%lit_param(i,8)*p + gl%el%lit_param(i,9)*p**2.d0 + gl%el%lit_param(i,10)*p**3.d0 )*temp + ( gl%el%lit_param(i,11) + gl%el%lit_param(i,12)*p + &
            & gl%el%lit_param(i,13)*p**2.d0)*temp**2.d0 + (gl%el%lit_param(i,14) + gl%el%lit_param(i,15)*p +gl%el%lit_param(i,16)*p**2.d0 + gl%el%lit_param(i,17)*p**3.d0 )/(temp - 227.d0) &
            & + (gl%el%lit_param(i,18) + gl%el%lit_param(i,19)*p + gl%el%lit_param(i,20)*p**2.d0 + gl%el%lit_param(i,21)*p**3.d0)/(680.d0 -temp)
        p = p/10.d0
    case(2) !fittted
        !P = p*10.d0
        !
        !p = p / p_red

        f_param = gl%el%pitzer_param(i,1) / Temp + gl%el%pitzer_param(i,2) + gl%el%pitzer_param(i,3)*p + gl%el%pitzer_param(i,4)*p**2.d0 + gl%el%pitzer_param(i,5)*p**3.d0 + gl%el%pitzer_param(i,6)*dlog(temp) &
            & +( gl%el%pitzer_param(i,7) + gl%el%pitzer_param(i,8)*p + gl%el%pitzer_param(i,9)*p**2.d0 + gl%el%pitzer_param(i,10)*p**3.d0 )*temp + ( gl%el%pitzer_param(i,11) + gl%el%pitzer_param(i,12)*p + &
            & gl%el%pitzer_param(i,13)*p**2.d0)*temp**2.d0 + (gl%el%pitzer_param(i,14) + gl%el%pitzer_param(i,15)*p +gl%el%pitzer_param(i,16)*p**2.d0 + gl%el%pitzer_param(i,17)*p**3.d0 )/(temp - 227.d0) &
            & + (gl%el%pitzer_param(i,18) + gl%el%pitzer_param(i,19)*p + gl%el%pitzer_param(i,20)*p**2.d0 + gl%el%pitzer_param(i,21)*p**3.d0)/(680.d0 -temp)

        !p = p * p_red
        !
        !p = p/10.d0
    end select

    end function f_param
    !-----------------------------------------------------------------------------
    !********************************************************************************


    !--------------------------------------------------------------------------------
    !Derivative of Parameterfunction for temperature and pressure dependent parameters  wrt temperature )
    double precision function df_param_dt(gl, i, temp, p)
    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el
    integer :: i, j
    double precision :: temp, p, p_red

    p_red = gl%el%pred

    select case(gl%el%functionflag)

    case(1) ! lit pitzer
        P = p*10.d0
        df_param_dt = (- gl%el%lit_param(i,1)) / Temp**2.d0 + gl%el%lit_param(i,6) / temp &
            & +( gl%el%lit_param(i,7) + gl%el%lit_param(i,8)*p + gl%el%lit_param(i,9)*p**2.d0 + gl%el%lit_param(i,10)*p**3.d0 ) + &
            & 2.d0 *( gl%el%lit_param(i,11) + gl%el%lit_param(i,12)*p + gl%el%lit_param(i,13)*p**2.d0)*temp + &
            &( -1.d0*(gl%el%lit_param(i,14) + gl%el%lit_param(i,15)*p +gl%el%lit_param(i,16)*p**2.d0 + gl%el%lit_param(i,17)*p**3.d0 )/(temp - 227.d0)**2.d0 )&
            & + (gl%el%lit_param(i,18) + gl%el%lit_param(i,19)*p + gl%el%lit_param(i,20)*p**2.d0 + gl%el%lit_param(i,21)*p**3.d0)/(680.d0 -temp)**2.d0
        p = p/10.d0
    case(2) !fitted
        !P = p*10.d0
        !
        !p = p / p_red

        df_param_dt = (- gl%el%pitzer_param(i,1)) / Temp**2.d0 + gl%el%pitzer_param(i,6) / temp &
            & +( gl%el%pitzer_param(i,7) + gl%el%pitzer_param(i,8)*p + gl%el%pitzer_param(i,9)*p**2.d0 + gl%el%pitzer_param(i,10)*p**3.d0 ) + &
            & 2.d0 *( gl%el%pitzer_param(i,11) + gl%el%pitzer_param(i,12)*p + gl%el%pitzer_param(i,13)*p**2.d0)*temp + &
            &( -1.d0*(gl%el%pitzer_param(i,14) + gl%el%pitzer_param(i,15)*p +gl%el%pitzer_param(i,16)*p**2.d0 + gl%el%pitzer_param(i,17)*p**3.d0 )/(temp - 227.d0)**2.d0 )&
            & + (gl%el%pitzer_param(i,18) + gl%el%pitzer_param(i,19)*p + gl%el%pitzer_param(i,20)*p**2.d0 + gl%el%pitzer_param(i,21)*p**3.d0)/(680.d0 -temp)**2.d0
        !
        !p = p * p_red
        !
        !p = p/10.d0
    end select

    end function df_param_dt
    !-----------------------------------------------------------------------------
    !*****************************************************************************


    !--------------------------------------------------------------------------------
    !Derivative of Parameterfunction for temperature and pressure dependent parameters  wrt temperature )
    double precision function d2f_param_dtdp(gl, i, temp, p)
    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el
    integer :: i, j
    double precision :: temp, p, p_red

    p_red = gl%el%pred

    select case(gl%el%functionflag)

    case(1) ! lit pitzer
        P = p*10.d0
        d2f_param_dtdp = ( gl%el%lit_param(i,8) + gl%el%lit_param(i,9)*p*2.d0 + 3.d0*gl%el%lit_param(i,10)*p**2.d0 ) + &
            & 2.d0 *( gl%el%lit_param(i,12) + gl%el%lit_param(i,13)*p*2.d0)*temp + &
            &( -1.d0*( gl%el%lit_param(i,15) +gl%el%lit_param(i,16)*p*2.d0 + 2.d0*gl%el%lit_param(i,17)*p**2.d0 )/(temp - 227.d0)**2.d0 )&
            & + ( gl%el%lit_param(i,19) + gl%el%lit_param(i,20)*p*2.d0 + 3.d0*gl%el%lit_param(i,21)*p**2.d0)/(680.d0 -temp)**2.d0
        p = p/10.d0
    case(2) !fitted


        d2f_param_dtdp = ( gl%el%pitzer_param(i,8) + gl%el%pitzer_param(i,9)*p*2.d0 + 3.d0*gl%el%pitzer_param(i,10)*p**2.d0 ) + &
            & 2.d0 *( gl%el%pitzer_param(i,12) + gl%el%pitzer_param(i,13)*p*2.d0)*temp + &
            &( -1.d0*( gl%el%pitzer_param(i,15) +gl%el%pitzer_param(i,16)*p*2.d0 + 2.d0*gl%el%pitzer_param(i,17)*p**2.d0 )/(temp - 227.d0)**2.d0 )&
            & + ( gl%el%pitzer_param(i,19) + gl%el%pitzer_param(i,20)*p*2.d0 + 3.d0*gl%el%pitzer_param(i,21)*p**2.d0)/(680.d0 -temp)**2.d0

    end select

    end function d2f_param_dtdp
    !-----------------------------------------------------------------------------
    !*****************************************************************************


    !--------------------------------------------------------------------------------
    !Derivative of Parameterfunction for temperature and pressure dependent parameters  wrt temperature**2 (second deriv) )
    double precision function df2_param_dt2(gl, i, temp, p)
    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el
    integer :: i, j
    double precision :: temp, p, dfdt1



    select case(gl%el%functionflag)

    case(1) ! pitzer
        df2_param_dt2 = (2.d0 *gl%el%lit_param(i,1)) / Temp**3.d0 + (-1.d0) * gl%el%lit_param(i,6) / temp**2.d0 &
            & + 2.d0 *( gl%el%lit_param(i,11) + gl%el%lit_param(i,12)*p + gl%el%lit_param(i,13)*p**2.d0) + &
            &( 2.d0*(gl%el%lit_param(i,14) + gl%el%lit_param(i,15)*p +gl%el%lit_param(i,16)*p**2.d0 + gl%el%lit_param(i,17)*p**3.d0 )/(temp - 227.d0)**3.d0 )&
            & + (2.d0* (gl%el%lit_param(i,18) + gl%el%lit_param(i,19)*p + gl%el%lit_param(i,20)*p**2.d0 + gl%el%lit_param(i,21)*p**3.d0))/(680.d0 -temp)**3.d0
        dfdt1 = 2.d0*df_param_dt(gl, i, temp,p)/temp
        df2_param_dt2 = df2_param_dt2 + dfdt1

    case(2) !fitted

        df2_param_dt2 = (2.d0 *gl%el%pitzer_param(i,1)) / Temp**3.d0 + (-1.d0) * gl%el%pitzer_param(i,6) / temp**2.d0 &
            & + 2.d0 *( gl%el%pitzer_param(i,11) + gl%el%pitzer_param(i,12)*p + gl%el%pitzer_param(i,13)*p**2.d0) + &
            &( 2.d0*(gl%el%pitzer_param(i,14) + gl%el%pitzer_param(i,15)*p +gl%el%pitzer_param(i,16)*p**2.d0 + gl%el%pitzer_param(i,17)*p**3.d0 )/(temp - 227.d0)**3.d0 )&
            & + (2.d0* (gl%el%pitzer_param(i,18) + gl%el%pitzer_param(i,19)*p + gl%el%pitzer_param(i,20)*p**2.d0 + gl%el%pitzer_param(i,21)*p**3.d0))/(680.d0 -temp)**3.d0

        dfdt1 = 2.d0*df_param_dt(gl, i, temp,p)/temp
        df2_param_dt2 = df2_param_dt2 + dfdt1
    end select

    end function df2_param_dt2
    !-----------------------------------------------------------------------------
    !*******************************************************************************

    !----------------------------------------------------------------------------
    !Fuction for coefficient of B_MX (see Rowland 2013, chap. 3.2), or Pitzer 84
    ! has nothing to do with Gibbs energy!!!!!!!
    double precision function g(gl, alpha, sqrt_Ion_strength)
    implicit none
    type (type_gl) :: gl
    double precision :: alpha, sqrt_Ion_strength, dummy
    !g does not have to do anything with the Gibbs energy. Its just a function for easier programming
    g = 2.d0*(1.d0 - (1.d0 + alpha * sqrt_Ion_strength)*dexp(- alpha * sqrt_Ion_strength))/(alpha * sqrt_Ion_strength)**2.d0

    end function g
    !-----------------------------------------------------------------------------
    !*******************************************************************************


    !---------------------------------------------------------------------------------------------
    !function for the calcutation of the brine density
    double precision function brine_dens_el(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw ,p, mol_wm, v_w, v_02, v_ex, v_t, rho_tot, Y, m, m_salt, v_ref_mr, v_ex_mr, M_2, deybe, mass, v_phi

    dw = rhomix_calc(gl, t, p, 0.d0, 1, gl%el%solpos) !* gl%wm(1)
    v_w = 1.d0/dW
    m = gl%el%molality(1)
    M_2 =  0.05844277d0

    v_phi = v_app(gl, t, dw, p)

    v_t = (m * v_phi + v_w/gl%wm(1) ) / (1.d0 + m*M_2)

    brine_dens_el = 1.d0 / v_t

    end function brine_dens_el
    !------------------------------------------------------------------------------------------
    !******************************************************************************************


    !---------------------------------------------------------------------------------------------
    !function for the calcutation of the molar volume of the brine (dg_brine/dp)_T
    double precision function dbrine_dp(gl, t, dw, p)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw ,p, mol_wm, v_w, v_02, v_ex, v_t, rho_tot, Y, m, m_salt, v_ref_mr, v_ex_mr, M_2, deybe, mass, v_phi, nw, n2

    dw = rhomix_calc(gl, t, p, 0.d0, 1, gl%el%solpos) !* gl%wm(1)
    v_w = 1.d0/dW
    m = gl%el%molality(1)
    n2 = gl%el%mol_in
    M_2 =  0.05844277d0
    nw = 1.d0/gl%wm(1)

    v_phi = v_app(gl, t, dw, p)

    dbrine_dp = (m* v_phi + v_w/gl%wm(1) ) / (nw + m)
    dbrine_dp = (n2* v_phi + v_w/gl%wm(1) ) / (nw + n2)
    if(gl%el%mol_in .le. 1.d-6) then
        dbrine_dp = v_w
    end if
    !molar volume of the brine
    end function dbrine_dp
    !------------------------------------------------------------------------------------------
    !******************************************************************************************

    !------------------------------------------------------------------------------------------
    double precision function cp_brine(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_w, cp_salt, deybe, app_cp, m

    m = gl%el%molality(1)

    app_cp = apparent_cp(gl, t, dw, p)

    cp_w = cp_calc(gl, t, dw, 1)*gl%wm(1)

    cp_brine = (cp_w + app_cp) / (1.d0/gl%wm(1) + m)


    end function cp_brine
    !------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    double precision function dg2_brine_dt2(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_w, cp_salt, deybe, app_cp, m, cp_b

    cp_b = brine_cp(gl, t, dw, p)

    dg2_brine_dt2 = - cp_b / T

    end function dg2_brine_dt2
    !------------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------
    double precision function apparent_cp(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_ex, cp_salt, deybe,cp_ref, m_save

    cp_ex = d2ge_dt2(gl, t, dw, p)!*gl%el%molality(1)!*gl%req(1)


    cp_salt = cp0_salt(gl, t, dw, p)

    apparent_cp = (cp_ex + cp_salt)

    end function apparent_cp
    !------------------------------------------------------------------------------------------------

    !********************************************************************************************
    !********************************************************************************************
    double precision function cp0_salt(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_ex_ref, m_ref, M_salt, M_water, n_water, n_salt, cp_ex, m, n_ref, cp_w, cp0_mr, deybe
    logical :: gemod

    !m_ref = 5.5508d0
    m_ref = gl%el%m_ref
    m = gl%el%molality(1)
    M_salt = 0.05844277d0

    gl%el%molality(1) = m_ref

    cp_ex_ref = d2ge_dt2(gl, t, dw, p)!*m_ref!*gl%req(1)

    gl%el%molality(1) = m


    gemod = gl%gecarrier
    gl%gecarrier = .false.
    cp_w = cp_calc(gl, t, dw, gl%el%solpos) / gl%wm(1)
    gl%gecarrier = gemod

    cp0_mr = cp_ref(gl, t, dw, p)!/n_ref

    cp0_salt = cp0_mr -  (cp_w / m_ref) -  cp_ex_ref


    end function cp0_salt
    !------------------------------------------------------------------------------------------------


    !********************************************************************************************
    !********************************************************************************************
    double precision function cp_ref(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_t_ref, cp0_ref, v_Ref, p0, cp_p

    p0 = 0.1d0

    cp_ref = gl%el%cp0_param(1) + gl%el%cp0_param(2) * T + gl%el%cp0_param(3) *t**2.d0

    !presspart from V0_ref:
    cp_p = 2.d0/3.d0 * T * gl%el%g0_param(1,9)*(p**3.d0 - p0**3.d0) + T * gl%el%g0_param(1,7)*(p**2.d0 - p0**2.d0) &
        & + (p - p0)* (6.d0*T**2.d0*gl%el%g0_param(1,4) - 2.d0*T*p0**2.d0*gl%el%g0_param(1,9) - 2.d0*T*p0*gl%el%g0_param(1,7) + 2.d0*T*gl%el%g0_param(1,3))

    cp_ref = cp_ref + cp_p

    end function cp_ref
    !------------------------------------------------------------------------------------------------
    !*********************************************************************************************


    !********************************************************************************************
    !********************************************************************************************
    double precision function hm_ref(gl, t, dw, p)
    !
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_t_ref, cp0_ref, v_Ref, p0, cp_p, t0

    p0 = 0.1d0
    t0 = 298.15d0

    hm_ref = (t**3.d0 - t0**3.d0) * (gl%el%cp0_param(3)/3.d0 + 2.d0 * (p - p0) * gl%el%g0_param(1,4) ) + (T-T0)*gl%el%cp0_param(1) &
        & + (T**2.d0 - T0**2.d0) * (gl%el%cp0_param(2)/2.d0 + gl%el%g0_param(1,9)/3.d0 * p**3.d0 + gl%el%g0_param(1,7)/2.d0 * p**2.d0 - p*p0**2.d0*gl%el%g0_param(1,9) &
        &  -p*p0*gl%el%g0_param(1,7) + p*gl%el%g0_param(1,3) + 2.d0/3.d0 * p0**3.d0*gl%el%g0_param(1,9) + p0**2.d0/2.d0 * gl%el%g0_param(1,7) - p0*gl%el%g0_param(1,3) )

    !g0 = G_02(gl, t, dw, p)
    !s0 = s_02(gl, t, dw, p)
    end function hm_ref
    !------------------------------------------------------------------------------------------------
    !*********************************************************************************************



    !********************************************************************************************
    !********************************************************************************************
    double precision function h0_2(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, cp_t_ref, cp0_ref, v_Ref, p0, cp_p, mols, nr, hw, h_ex_mr, href_m

    hw = h_calc(gl, t, dw, 1)/gl%wm(1)

    nr = gl%el%molality(1)

    mols = gl%el%molality(1)
    gl%el%molality(1) = gl%el%m_ref

    h_ex_mr = h_ex(gl, t, dw, p)
    gl%el%molality(1) = mols
    href_m = hm_ref(gl, t, dw, p)

    h0_2 = href_m - hw/nr - h_ex_mr

    end function h0_2
    !------------------------------------------------------------------------------------------------
    !*********************************************************************************************



    !********************************************************************************************
    !********************************************************************************************
    double precision function s0_2(gl, t, dw, p)
    !function for the apparant molal heat capacity (isobaric) cp
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, h, g

    h = h0_2(gl, t, dw, p)
    g = g0_2(gl, t, dw, p)

    s0_2 = (h-g)/T/gl%req(1)

    end function s0_2
    !----------------------------------------------------------------------------------------------
    !*********************************************************************************************


    !***********************************************************************************************
    !**********************************************************************************************
    double precision function brine_cp(gl, t, dw, p)
    !function for the heat capacity of 1 kg brine (not 1kg w + m moles salt)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, app_cp, cp_w, moles_water, M_salt
    logical :: gemod

    app_cp = apparent_cp(gl, t, dw, p)

    gemod = gl%gecarrier
    gl%gecarrier = .false.
    cp_w = cp_calc(gl, t, dw, 1)/gl%wm(1)
    gl%gecarrier = gemod

    moles_water = 1.d0/gl%wm(1)

    M_salt = 0.05844277d0

    brine_cp = ( cp_w + gl%el%molality (1) * app_cp) / (moles_water + gl%el%molality(1))
    !brine_cp = brine_cp_spec(gl, t, dw, p)
    !brine_cp = brine_cp*wm_brine(gl)
    end function brine_cp
    !---------------------------------------------------------------------------------------------------
    !**********************************************************************************************


    !***********************************************************************************************
    !**********************************************************************************************
    double precision function brine_cp_spec(gl, t, dw, p)
    !function for the heat capacity of 1 kg brine (not 1kg w + m moles salt)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, app_cp, cp_w, moles_water, M_salt
    logical :: gemod

    app_cp = apparent_cp(gl, t, dw, p)

    gemod = gl%gecarrier
    gl%gecarrier = .false.
    cp_w = cp_calc(gl, t, dw, 1)/gl%wm(1)
    gl%gecarrier = gemod

    moles_water = 1/gl%wm(1)

    M_salt = 0.05844277d0
    brine_cp_spec = ( cp_w + gl%el%molality (1) * app_cp) / (1.d0 + gl%el%molality(1)*M_salt)


    end function brine_cp_spec
    !---------------------------------------------------------------------------------------------------
    !**********************************************************************************************



    !***********************************************************************************************
    !**********************************************************************************************
    double precision function d2brine_dt2(gl, t, dw, p)
    !function for the heat capacity of 1 kg brine (not 1kg w + m moles salt)
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p

    d2brine_dt2 = - brine_cp(gl, t, dw, p) / T

    end function d2brine_dt2
    !---------------------------------------------------------------------------------------------------
    !**********************************************************************************************


    !----------------------------------------------------------------------
    !Subroutine for the calculation of gE for Electrolytes, based on Pitzer's Equations. E.g. Rogers and Pitzer 1982, Rowland 2013,
    !subroutine gE_Pitzer(gl, temp, dw, p, gE)
    double precision function dgE_dn(gl, temp, dw, p)
    !-------------------------------------------------------------------------------
    implicit none

    type(type_gl) :: gl
    !type (electrolytes_vars) :: el

    integer :: salt

    double precision :: dens_sol, Ms, de, Ax, I, temp, p, mo, ge, Cmx, pred, pred2, tred, tred2, beta0, beta1, beta2, R, dw, press, h, m, B, C, t!, f_param

    !call electro_init(gl%el)
    salt = 1
    R = gl%req(1)
    m = gl%el%molality(1)

    press = p
    t = temp

    Ax =Deybe_hueckel(gl, t, dw, p)
    I= ionic_Strength(gl)

    h =  (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )

    B = dB_mx_dn(gl, temp, press, I)
    C = C_mx(gl, temp, press, gl%el%zm, gl%el%zx)

    dgE_dn = -2.d0* Ax  * h + 2.d0 * gl%el%vm*gl%el%vx * ( m**2.d0 * B + 2.d0*m**3.d0 *gl%el%vm*gl%el%zm*C )
    dge_dn = dge_dn*R*t
    !Gibbs excess energy of the brine, reffered to moles of water in the brine and divided by RT!!
    end function dgE_dn! for Electrolytes
    !--------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------





    !------------------------------------------------------------------------------------------
    !function for the chemical potential of water in brines according tp Pitzer; refferd to mol
    double precision function chempot_water_brine(gl, t, dw, p)
    !see eq.24 PPB84 for more information
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, Ax, beta0, beta1,m, Cmx, v, I, g_water, t0, p0, nr, mw, nw, n_salt, dw_ref, g_w_ref, s_w_ref, g_w, mol_save, ge_tpm_ref, se_tpm_ref, ge_m_ref, G02_part, S_salt_ref, water_Ref, ex_part, s_salt
    double precision :: mm, mp, gm, gp, gbrine, sal_save, n2_in, nw_p, nw_m, ln_part, ge_part, R, n2, mix_num_p, mix_num_m, mix_num, gep, gem, ge, xs, xsp, xsm
    double precision, dimension(30) :: chempot
    logical :: gesave

    v = 2.d0
    r = gl%req(1)
    gesave = gl%gecarrier
    gl%gecarrier = .false.
    Ax = deybe_hueckel(gl, t, dw, p)
    !beta0 = beta0_func(gl, t, p)
    beta0 = f_param(gl,1,222, t, p)
    !beta1 = beta1_func(gl, t, p)
    beta1 = f_param(gl,2,222, t, p)
    Cmx = C_mx(gl, t, p, gl%el%zm,gl%el%zx)!/(2.d0*dabs(gl%el%zx*gl%el%zm)**0.5d0)
    I = Ionic_strength(gl)
    m = gl%el%molality(1)

    ln_part = gl%wm(1) *(- v * m)
    ge_part = gl%wm(1) * ( 2.d0 * Ax * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * m**2.d0 * ( beta0 + beta1 * dexp(- gl%el%alpha1 * I**0.5d0) +2.d0 * gl%el%vm * gl%el%zm * m* Cmx ) )

    ln_part = ln_part * gl%req(1) * T
    ge_part = ge_part * gl%req(1) * T

    chempot_water_brine = gl%wm(1) * ( - v * m + 2.d0 * Ax * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * m**2.d0 * ( beta0 + beta1 * dexp(- gl%el%alpha1 * I**0.5d0) +2.d0 * gl%el%vm * gl%el%zm * m* Cmx ) )

    chempot_water_brine = chempot_water_brine * gl%req(1) * T

    g_water = g_calc(gl, t, dw, 1)

    chempot_Water_brine = chempot_water_brine + g_water

    t0 = 298.15d0
    p0 = 0.1d0
    nr = gl%el%m_ref
    S_salt_ref = gl%el%g0_param(1,11)

    Mw = gl%wm(1)
    nw = 1.d0/gl%wm(1)
    n_salt = gl%el%molality(1)
    !n_salt = gl%el%x_salt

    !call reduced_parameters_calc(gl,T)
    dw_ref = rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)
    g_w_ref = g_calc(gl, t0, dw_ref, 1)
    s_w_ref = s_calc(gl, t0, dw_ref, 1)

    g_w  = g_calc(gl, t, dw, 1)
    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    ge_tpm_ref = dge_dn(gl, t0, dw_ref, p0)!*mw
    se_tpm_ref = dsE_dn(gl, t0, dw_ref, p0)
    ge_m_ref = dge_dn(gl, t, dw, p)!*mw

    gl%el%molality(1) = mol_save

    !s_salt = -(t-t0) * mw*n_salt*s_salt_ref
    water_ref = n_salt/nr*mw * (g_w_ref - g_w ) - (t-t0)*Mw*n_salt * s_w_ref/nr
    ex_part =  ( Mw*(ge_tpm_ref - ge_m_ref) /nr - mw*(t-t0)*se_tpm_ref/nr )*mw*n_salt!**2.d0
    G02_part = (water_ref + ex_part)

    chempot_Water_brine = chempot_Water_brine + G02_part
   
    end function chempot_water_brine
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------










    !
    !------------------------------------------------------------------------------------------
    !function for the chemical potential of water in brines according tp Pitzer
    double precision function dchempot_water_brine_dT(gl, t, dw, p)
    !see eq.24 PPB84 for more information
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, Ax, beta0, beta1, Cmx, v, I, g_water,tp, tm, chempotp, chempotm, dp, dm, axp, axm, beta0p, beta1p, cmxp, beta0m, beta1m, cmxm, gwp, gwm, numdiff
    double precision :: t0, p0, nr, s_salt_ref, Mw, nw, n_salt, dw_ref, g_w_ref, s_w_ref, mol_save,numdiffref, t0p, t0m, d0m, d0p, ge_tpm_ref, se_tpm_ref, ge_m_refp, ge_m_refm, g02_partp, g02_partm, ge_m_ref, water_ref, ex_part

    numdiff = 1.d-6
    numdiffref = 1.d-6

    v = 2.d0
    I = Ionic_strength(gl)

    t0 = 298.15d0
    p0 = 0.1d0
    nr = gl%el%m_ref
    S_salt_ref = gl%el%g0_param(1,11)

    Mw = gl%wm(1)
    nw = 1.d0/gl%wm(1)
    n_salt = gl%el%molality(1)
    !n_salt = gl%el%x_salt

    !call reduced_parameters_calc(gl,T)
    dw_ref = rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)!55344.5576839165d0 !rhomix_calc(gl, t0, p0, 0.d0, 1, gl%el%solpos)
    g_w_ref = g_calc(gl, t0, dw_ref, 1) !Auf Variablen speichern!
    s_w_ref = s_calc(gl, t0, dw_ref, 1)


    !upper T
    tp = T*(1.d0 + numdiff)
    t0p = t0*(1.d0+numdiffref)
    !dp = rhomix_calc(gl, tp, p, 0.d0, 2, gl%el%solpos)
    dp = rhomix_calc(gl, tp, p, 0.d0, 1, gl%el%solpos)
    d0p = rhomix_calc(gl, t0p, p, 0.d0, 1, gl%el%solpos)
    Axp = deybe_hueckel(gl, tp, dp, p)
    beta0p = f_param(gl,1,222, tp, p)
    beta1p = f_param(gl,2,222, tp, p)
    Cmxp = C_mx(gl, tp, p, gl%el%zm,gl%el%zx)!/(2.d0*dabs(gl%el%zx*gl%el%zm)**0.5d0)

    chempotp = gl%wm(1) * ( - v * gl%el%molality(gl%el%salt) + 2.d0 * Axp * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * gl%el%molality(gl%el%salt)**2.d0 * ( beta0p + beta1p * dexp(- gl%el%alpha1 * I**0.5d0) &
        & +2.d0 * gl%el%vm * gl%el%zm * gl%el%molality(gl%el%salt)* Cmxp ) )

    chempotp = chempotp * gl%req(1) * Tp

    gwp = g_calc(gl, tp, dp, 1)

    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    !g_wp  = g_calc(gl, t, dw, 1)
    !mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    ge_tpm_ref = dge_dn(gl, t0p, d0p, p0)!*mw
    se_tpm_ref = dsE_dn(gl, t0p, d0p, p0)!*mw
    ge_m_ref = dge_dn(gl, tp, dp, p)!*mw

    gl%el%molality(1) = mol_save

    !s_salt = -(t-t0) * mw*n_salt*s_salt_ref
    water_ref = n_salt/nr*mw * (g_w_ref - gwp ) - (t-t0)*Mw*n_salt * s_w_ref/nr
    ex_part =  ( Mw*(ge_tpm_ref - ge_m_ref) /nr - mw*(t-t0)*se_tpm_ref/nr )*mw*n_salt
    G02_partp = -(water_ref + ex_part)

    chempotp = chempotp + gwp + g02_partp


    !------------
    !lower T
    tm = T*(1.d0 - numdiff)
    t0m = t0*(1.d0-numdiffref)
    dm = rhomix_calc(gl, tm, p, 0.d0, 1, gl%el%solpos)
    d0m = rhomix_calc(gl, t0m, p, 0.d0, 1, gl%el%solpos)
    Axm = deybe_hueckel(gl, tm, dm, p)
    beta0m = f_param(gl,1,222, tm, p)
    beta1m = f_param(gl,2,222, tm, p)
    Cmxm = C_mx(gl, tm, p, gl%el%zm,gl%el%zx)!/(2.d0*dabs(gl%el%zx*gl%el%zm)**0.5d0)

    chempotm = gl%wm(1) * ( - v * gl%el%molality(gl%el%salt) + 2.d0 * Axm * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * gl%el%molality(gl%el%salt)**2.d0 * ( beta0m + beta1m * dexp(- gl%el%alpha1 * I**0.5d0) &
        & +2.d0 * gl%el%vm * gl%el%zm * gl%el%molality(gl%el%salt)* Cmxm) )

    chempotm = chempotm * gl%req(1) * Tm

    gwm = g_calc(gl, tm, dm, 1)

    mol_save = gl%el%molality(1)
    gl%el%molality(1) = nr
    !g_wm  = g_calc(gl, tm, dw, 1)
    !mol_save = gl%el%molality(1)
    !gl%el%molality(1) = nr
    ge_tpm_ref = dge_dn(gl, t0m, d0m, p0)!*mw
    se_tpm_ref = dsE_dn(gl, t0m, d0m, p0)!*mw
    ge_m_ref = dge_dn(gl, tm, dm, p)!*mw

    gl%el%molality(1) = mol_save

    !s_salt = -(t-t0) * mw*n_salt*s_salt_ref
    water_ref = n_salt/nr*mw * (g_w_ref - gwm ) - (t-t0)*Mw*n_salt * s_w_ref/nr
    ex_part =  ( Mw*(ge_tpm_ref - ge_m_ref) /nr - mw*(t-t0)*se_tpm_ref/nr )*mw*n_salt
    G02_partm = -(water_ref + ex_part)

    ex_part = (g02_partp-g02_partm)/(tp-tm)

    chempotm = chempotm + gwm + g02_partm
    !------------------
    !num diff:

    dchempot_water_brine_dT = (chempotp-chempotm)/(tp-tm)

    end function dchempot_water_brine_dT
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------------------------------
    !function for the chemical potential of water in brines according tp Pitzer
    double precision function dchempot_water_brine_dp(gl, t, dw, p)
    !see eq.24 PPB84 for more information
    implicit none
    type(type_gl) :: gl
    double precision :: t, dw, p, Ax, beta0, beta1, Cmx, v, I, g_water,pp, pm, chempotp, chempotm, dp, dm, axp, axm, beta0p, beta1p, cmxp, beta0m, beta1m, cmxm, gwp, gwm, numdiff

    numdiff = 1.d-6

    v = 2.d0
    I = Ionic_strength(gl)

    !upper T
    pp = p*(1.d0 + numdiff)
    dp = rhomix_calc(gl, t, pp, 0.d0, 2, gl%el%solpos)
    Axp = deybe_hueckel(gl, t, dp, pp)
    beta0p = f_param(gl,1,222, t, pp)
    beta1p = f_param(gl,2,222, t, pp)
    Cmxp = C_mx(gl, t, pp, gl%el%zm,gl%el%zx)/(2.d0*dabs(gl%el%zx*gl%el%zm)**0.5d0)

    chempotp = gl%wm(1) * ( - v * gl%el%molality(gl%el%salt) + 2.d0 * Axp * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * gl%el%molality(gl%el%salt)**2.d0 * ( beta0p + beta1p * dexp(- gl%el%alpha1 * I**0.5d0) &
        & +2.d0 * gl%el%vm * gl%el%zm * gl%el%molality(gl%el%salt)* Cmxp ) )

    chempotp = chempotp * gl%req(1) * T

    gwp = g_calc(gl, t, dp, 0)

    chempotp = chempotp + gwp
    !------------
    !lower T
    pm = p*(1.d0 - numdiff)
    dm = rhomix_calc(gl, t, pm, 0.d0, 2, gl%el%solpos)
    Axm = deybe_hueckel(gl, t, dm, pm)
    beta0m = f_param(gl,1,222, t, pm)
    beta1m = f_param(gl,2,222, t, pm)
    Cmxm = C_mx(gl, t, pm, gl%el%zm,gl%el%zx)/(2.d0*dabs(gl%el%zx*gl%el%zm)**0.5d0)

    chempotm = gl%wm(1) * ( - v * gl%el%molality(gl%el%salt) + 2.d0 * Axm * ( (I**1.5d0) / ( 1.d0 + gl%el%b(gl%el%salt) * I**0.5d0 )) &
        & -2.d0 * gl%el%vm * gl%el%vx * gl%el%molality(gl%el%salt)**2.d0 * ( beta0m + beta1m * dexp(- gl%el%alpha1 * I**0.5d0) &
        & +2.d0 * gl%el%vm * gl%el%zm * gl%el%molality(gl%el%salt)* Cmxm ) )

    chempotm = chempotm * gl%req(1) * T

    gwm = g_calc(gl, t, dm, 0)

    chempotm = chempotm + gwm
    !------------------
    !num diff:

    dchempot_water_brine_dp = (chempotp-chempotm)/(pp-pm)

    end function dchempot_water_brine_dp
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------



    !------------------------------------------------------------------
    !Function for the calculation of the ionic strength, depeding on salt and molality
    double precision function ionic_strength(gl)
    !------------------------------------------
    implicit none
    type (type_gl) :: gl
    double precision :: m, z
    integer :: i

    ionic_strength = 0.d0

    do i = 0,gl%el%n_salts !has to be zeor, change for higher order salts
        ionic_strength = ionic_strength  + gl%el%molality(gl%el%n_salts) * gl%el%z_salt(gl%el%n_salts)
    end do

    ionic_strength = ionic_strength * 0.5d0

    end function ionic_strength
    !-------------------------------------------------------------------------------
    !***********************************************************************************



    !**********************************************************************************
    double precision function osmotic_coeff(gl, t, d, p)
    !Function for the calculaton of the osmotic coefficient (Pitzer etal 1984), eq(3)
    implicit none
    type(type_gl) :: gl
    double precision :: t, d, p, phi, v, R, I, A, b0, b1, Cphi, m, z

    z = abs(gl%el%zm * gl%el%zx)
    A = deybe_hueckel(gl, t, d, p)
    I = Ionic_strength(gl)
    b0 = beta0_func(gl, t, p)
    b1 = beta1_func(gl, t, p)*dexp(- gl%el%alpha1*(I**0.5d0))
    Cphi = C_mx(gl, t, p, gl%el%zm, gl%el%zx)*(2.d0*z**0.5d0)
    m = gl%el%molality(1)
    v=2.d0

    phi = 1.d0 - z * A * (I**0.5d0 / (1.d0 + gl%el%b(1)*I**0.5d0) ) + m * ( (2.d0 * gl%el%vm*gl%el%vx)/ v ) * ( b0 + b1 ) &
        & + m**2.d0 * ( (2.d0 * (gl%el%vm*gl%el%vx)**1.5d0)/ v ) * Cphi

    osmotic_coeff = phi

    end function osmotic_coeff
    !*****************************************************************************************
    !****************************************************************************************


    !**********************************************************************************
    double precision function ln_el_activity(gl, t, d, p)
    !Function for the calculaton of the osmotic coefficient (Pitzer etal 1984), eq(3)
    implicit none
    type(type_gl) :: gl
    double precision :: t, d, p, y, v, R, I, A, b0, b1, Cphi, m, z, b, vm, vx, p1, p2, p3, f, fb1

    z = abs(gl%el%zm * gl%el%zx)
    A = deybe_hueckel(gl, t, d, p)
    I = Ionic_strength(gl)
    b0 = beta0_func(gl, t, p)
    b1 = beta1_func(gl, t, p)
    Cphi = C_mx(gl, t, p, gl%el%zm, gl%el%zx)*(2.d0*z**0.5d0)
    m = gl%el%molality(1)
    b=gl%el%b(1)
    vm = gl%el%vm
    vx = gl%el%vx
    v=2.d0
    f = (I**0.5d0 / (1.d0 + b*I**0.5d0)) + (2.d0/b)*dlog(1.d0 + b*I**0.5d0)
    !fb0
    fb1 = 2.d0/((gl%el%alpha1**2.d0) *I) * (1.d0 - (1.d0+gl%el%alpha1*(I**0.5d0) - (0.5d0*(gl%el%alpha1**2.d0)*I) )*dexp(-gl%el%alpha1*(I**0.5d0)) )

    p1 = - z * A * f
    p2 = + m * ( (2.d0 * gl%el%vm*gl%el%vx)/ v ) * ( 2.d0*b0 +  (2.d0*b1 / (gl%el%alpha1**2.d0 * I)) * (1.d0 - (1.d0 + gl%el%alpha1*I**0.5d0 - 0.5d0*(gl%el%alpha1**2.d0 * I) )*dexp(-gl%el%alpha1* I**0.5d0)) )
    p3 = + 0.5d0*3.d0*(m**2.d0) *( ((2.d0 * (vm*vx)**1.5d0)/ v)*Cphi )

    !y = - z * A * ( (I**0.5d0 / (1.d0 + b*I**0.5d0)) + (2.d0/b)*dlog(1.d0 + b*I**0.5d0) ) &
    !    & + m * ( (2.d0 * gl%el%vm*gl%el%vx)/ v ) * ( 2.d0*b0 +  (2.d0*b1 / (gl%el%alpha1**2.d0 * I)) * (1.d0 - (1.d0 + gl%el%alpha1*I**0.5d0 - 0.5d0*(gl%el%alpha1**2.d0 * I) )*dexp(-gl%el%alpha1* I**0.5d0)) ) &
    !    & + 0.5d0*3.d0*m**2.d0 *( ((2.d0 * (vm*vx)**1.5d0)/ v)*Cphi )
    y = p1+p2+p3

    ln_el_activity = y
    end function ln_el_activity
    !*****************************************************************************************
    !****************************************************************************************



    !**********************************************************************************
    double precision function dge_w_dnw(gl, t, d, p)
    !Function for the calculaton of the osmotic coefficient (Pitzer etal 1984), eq(3)
    implicit none
    type(type_gl) :: gl
    double precision :: t, d, p, v, R, osmo

    v = gl%el%vm + gl%el%vx
    osmo = osmotic_coeff(gl, t, d, p) - 1.d0
    R = gl%req(1)
    !osmo = 0.935d0 - 1.d0
    dge_w_dnw = osmo * (-v*R*T)

    end function dge_w_dnw
    !*****************************************************************************************
    !****************************************************************************************

    !**********************************************************************************
    double precision function dge_w_dnsalt(gl, t, d, p)
    !Function for the calculaton of the osmotic coefficient (Pitzer etal 1984), eq(3)
    implicit none
    type(type_gl) :: gl
    double precision :: t, d, p, v, R, activ

    v = gl%el%vm + gl%el%vx
    activ = ln_el_activity(gl, t, d, p)
    R = gl%req(1)

    dge_w_dnsalt = activ * (v*R*T)

    end function dge_w_dnsalt
    !*****************************************************************************************
    !****************************************************************************************

    !-----------------------------------------------------------------------------------------
    double precision Function molar_el(gl)
    !----------------------------------------------------------------------------
    !This Subroutine calculates the molar mass of seawater with respect to the salinity
    !-----------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl

    !double precision, dimension(8,8,7) :: GIJK_list
    double precision :: sal, M_salt, mass_w, mass_s, x_s, x_w

    sal = gl%el%molality(1)

    M_salt = 0.05844277d0 !Molar Mass of Salt

    mass_w = 1.d0
    mass_s = sal*M_salt

    x_s = (mass_s / M_salt) / ( (mass_s / M_salt ) + (mass_w /gl%wm(1) ) )

    x_w = 1.d0 - x_s

    molar_el = x_w * gl%wm(1) + x_s * M_salt

    end function !molar_el
    !***************************************************************************************





    !------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !---------------------------------------------------------------------------------
    !--------------Needed Routines for the setup of electrolytes within the setup rotuine in TREND
    !------------------------------------------------------------------------------------------
    !*********************************************************************************************


    !-----------------------------------------------------------------------
    !function for electrolytes eqtypes
    !looking for eq type of salts and saves the positions, cross out these eq types and lines, moves other up
    subroutine electrolytes_seteq(gl, eqtype, ne, errorflag)
    !---------------------------------------------------------------
    implicit none
    type (type_gl) :: gl
    integer, dimension(30) :: eqtype
    integer :: ne, eleven, i, checkcount, errorflag

    gl%el%eq_types_in = eqtype
    gl%el%n_salts = count(eqtype(:) .eq. 15)
    checkcount = 0
    do i=1,ne
        if(eqtype(i) == 15) then
            if(checkcount == 0) then
                gl%el%saltpos = i
                gl%el_present = .true.
                eqtype(i) = 1
            else
                errorflag = -99999.d0 !up to now only one salt is available
            end if
        end if
    end do
    if(gl%el_present) then
        do i=gl%el%saltpos, ne-1
            eqtype(i) = eqtype(i+1)
        end do
        eqtype(ne) = 0
    end if

    ne = ne - gl%el%n_salts
    if((gl%el%n_salts) .gt. 1) then
        ne = 31 !with the intention th get a setup error!
    end if
    gl%el%Eq_type_el = eqtype
    end subroutine electrolytes_seteq
    !-------------------------------------------------------------
    !******************************************************************************


    !-------------------------------------------------------------
    !creats the fluid array without the salts for the normal setup routnies
    subroutine electrolytes_setfluids(gl, fluids, nf)
    !-----------------------------------------------
    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el
    character(30), dimension(30) :: fluids
    integer :: nf, i, j, k

    gl%el%fluids_input = fluids

    if(gl%el%n_salts .gt. 1) then
        k = 1
        do j=1,nf+1
            if((j .lt. 30) .and. (j .lt. gl%el%saltpos)) then
                fluids(j) = fluids(j)
            elseif(j == gl%el%saltpos) then
                gl%el%salts(k) = fluids(j)
                k = k + 1
            elseif((j .lt. 30) .and. (j .ge. gl%el%saltpos) ) then
                fluids(j) = fluids(j+1)
            elseif((j .eq. 30) .and. (j .ge. gl%salpos)) then
                fluids(j) = ' '
                !else
                !   fluids(j) = fluids(j)
            end if
        end do
    elseif(gl%el%n_salts .eq. 1) then
        k = 1
        do j=1,nf+1
            if((j .lt. 30) .and. (j .lt. gl%el%saltpos)) then
                fluids(j) = fluids(j)
            elseif(j == gl%el%saltpos) then
                gl%el%salts(k) = fluids(j)
                fluids(j) = fluids(j+1)
                gl%el%solpos = j-1
                k = k+1
            elseif((j .lt. 30) .and. (j .ge. gl%el%saltpos) ) then
                fluids(j) = fluids(j+1)
            elseif((j .eq. 30) .and. (j .ge. gl%salpos)) then
                fluids(j) = ' '
                !else
                !   fluids(j) = fluids(j)
            end if
        end do
    end if

    gl%el%fluids_wo_Salts = fluids

    gl%el%n_fluids = nf - gl%el%n_salts

    end subroutine electrolytes_setfluids
    !--------------------------------------------------------------
    !*****************************************************************


    !------------------------------------------------------------------------
    !Routine for adjusting the moles array, sum of moles woud exceed 1
    !Saves the salt molalities and creats a vector without salts
    subroutine electrolytes_adjustmoles(gl, moles, nm)
    !-------------------------------------------------------------------------
    implicit none
    type(type_gl) :: gl

    double precision, dimension(30) :: moles
    integer :: nm, i, j, k
    gl%el%molality = 0.d0
    gl%el%moles_input = moles
    k = 1
    do j=1,30
        if((j .lt. 30) .and. (j .lt. gl%el%saltpos)) then
            moles(j) = moles(j)
        elseif(j == gl%el%saltpos) then
            gl%el%molality(k) = moles(j)
            k = k + 1
            moles(j) = moles(j+1)
        elseif((j .lt. 30) .and. (j .ge. gl%el%saltpos) ) then
            moles(j) = moles(j+1)
        elseif((j .eq. 30) .and. (j .ge. gl%salpos)) then
            moles(j) = 0.d0
            !else
            !   fluids(j) = fluids(j)
        end if
    end do

    !nm = nm - gl%el%n_salts

    end subroutine electrolytes_adjustmoles
    !---------------------------------------------------------------
    !***************************************************************


    !---------------------------------------------------------------
    !Adjusts fluids and moles arrays, as salts are crossed out
    subroutine setelectrolytes(gl, fluids, eqtype, molfraction, nf)
    !----------------------------------------------------------------
    implicit none
    type(type_gl) :: gl
    !type(electrolytes_vars) :: el

    character(30), dimension(30) :: fluids, saltlist
    integer :: water, j, t, nf
    integer, dimension(30) :: eqtype
    double precision, dimension(30) :: molfraction, moles
    character(30) :: tf

    saltlist(:) = ' '
    saltlist(1) = 'nacl'
    saltlist(2) = 'kcl'
    saltlist(3) = 'mgcl'

    moles = molfraction

    do j = 1, count(saltlist(:) /= '')
        tf = saltlist(j)
        if(index(tf," ") .eq. 1) exit
        if(gl%el%n_salts .eq. nf) exit
        water = count(fluids(:) == 'water')
        if(water /= 1) exit
        do t = 1,nf
            if (tf .eq. fluids(t)) then
                gl%el%n_salts = gl%el%n_salts  + 1
                gl%el%salts(gl%el%n_salts) = fluids(t)
                gl%el%molality(gl%el%n_salts) = molfraction(t)
                gl%el%mapping(gl%el%salt) = t
                if(t .lt. nf) then
                    fluids(t) = fluids(t+1)
                    moles(t) = molfraction(t+1)
                else
                    fluids(t) = ' '
                    moles(t) = 0.d0
                end if
            end if
        end do
    end do

    molfraction = moles

    do j = 1, nf
        if (gl%el%mapping(j) .ne. j) then
            gl%el%resorted = .TRUE.
            !call initialize_ideal_gas_coefficients(gl)
        end if
    end do


    end subroutine
    !-------------------------------------------------------------------------
    !**************************************************************************

    double precision function wm_brine(gl)
    implicit none
    type(type_gl) :: gl
    double precision :: n_w, n_s, n_tot, xw, xs

    n_w = 1.d0/gl%wm(1)
    n_s = gl%el%molality(1)
    n_tot = n_w + n_s

    xw = n_w / n_tot
    xs = 1.d0 - xw

    wm_brine = xw*gl%wm(1) + xs*gl%el%wm_salt

    end function

    !---------------------------------------------------------------
    subroutine normalize_fluids(gl)

    implicit none
    type (type_gl) :: gl
    !type (electrolytes_vars) :: el

    integer :: n_fluids, i
    double precision :: sum_fluids

    n_fluids = gl%ncomp - gl%el%n_salts

    sum_fluids = 0.d0

    do i=1,n_fluids

        sum_fluids = sum_fluids + gl%molfractions(i)

    end do

    !do i=1,n_fluids

    gl%molfractions = gl%molfractions / sum_fluids

    end subroutine normalize_fluids
    !-------------------------------------------------------------



    !---------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine el_remap(gl, fluids_in, moles_in,eq_types_in, x_phase_in, chempot_in, phasetype)


    implicit none

    type(type_gl) :: gl

    character(30), dimension(30) :: fluids_in, fluids_out
    double precision, dimension(30) :: moles_in, moles_out
    integer, dimension(5) :: phasetype
    integer, dimension(30) :: eq_types_out, eq_types_in
    double precision, dimension(30,5) :: x_phase_in, x_phase_out, prop_in, prop_out, chempot_in, chempot_out
    integer :: i, seaphases

    fluids_out = ' '
    moles_out = 0.d0
    x_phase_out = 0.d0
    chempot_out = 0.d0
    ! prop_out = 0.d0

    do i = 1,gl%el%solpos

        fluids_out(i) = fluids_in(i)
        moles_out(i) = moles_in(i)
        x_phase_out(i,:) = x_phase_in(i,:)
        chempot_out(i,:) = chempot_in(i,:)
        ! prop_out(i,:) = prop_in(i,:)

    end do
    fluids_out = gl%el%fluids_input
    moles_out = gl%el%moles_input
    eq_types_out = gl%el%eq_types_in
    !fluids_out(gl%el%solpos) = 'water'
    !fluids_out(gl%el%solpos+1) = 'NaCl'
    !moles_out(gl%el%solpos+1) = gl%el%molality(1)
    if(phasetype(1)==2 .or. phasetype(2)==2 .or. phasetype(3)==2) then
        x_phase_out(gl%el%solpos+1,2) = gl%el%molality(1)
    end if
    if(phasetype(1)==3 .or. phasetype(2)==3 .or. phasetype(3)==3) then
        x_phase_out(gl%el%solpos+1,3) = gl%el%molality(1)
    end if
    chempot_out(gl%el%solpos+1,:) = 0.d0
    ! prop_out(gl%salpos,:) = 0.d0


    do i = gl%el%solpos+2,gl%ncomp+1

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
    eq_types_in = eq_types_out
    !  prop_in = prop_out


    end subroutine !el_remap
    !**************************************************************************
    !**************************************************************************

    subroutine check_el_limits(gl, t, p, input, errorflag)
    implicit none
    type(type_gl)::gl
    double precision :: t, p
    character(12)::input
    integer :: errorflag
    logical :: valid_input

    valid_input = .false.
    if(input .eq. 'tp') then
        valid_input = .true.
    elseif(input .eq. 'tp+') then
        valid_input = .true.
    elseif(input .eq. 'pmltl+') then
        valid_input = .true.
    end if


    if(t>374.d0)then
        errorflag = -9936
        return
    elseif(p>50.d0)then
        errorflag = -9936
        return
    elseif(gl%el%molality(1)> 6.d0)then
        errorflag = -9936
        return
    elseif(.not. valid_input) then
        errorflag = -12910
        return
    elseif(gl%ncomp .gt. 2) then
        errorflag = -5223
        return
    end if

    end subroutine



    !----------------------------------------------------------
    subroutine dielectric_swap(gl, flag)
    implicit none
    type(type_gl) :: gl
    integer :: flag

    if(.not. allocated(gl%de)) allocate(gl%de)
    if(.not. allocated(gl%dew)) allocate(gl%dew)
    if(.not. allocated(gl%des)) allocate(gl%des)
    call initialize_dew(gl)


    if (flag == 1) then

        gl%des%de_read = gl%de%de_read
        gl%des%demodel = gl%de%demodel
        gl%des%tmin_de = gl%de%tmin_de
        gl%des%tmax_de = gl%de%tmax_de
        gl%des%pmax_de = gl%de%pmax_de
        gl%des%rhomax_de = gl%de%rhomax_de
        gl%des%ref_temp_de = gl%de%ref_temp_de
        gl%des%ref_dens_de = gl%de%ref_dens_de
        gl%des%ref_press_de = gl%de%ref_press_de
        gl%des%term_num1_de = gl%de%term_num1_de
        gl%des%term_num2_de = gl%de%term_num2_de
        gl%des%term_num3_de = gl%de%term_num3_de
        gl%des%term_num4_de = gl%de%term_num4_de
        gl%des%term_num5_de = gl%de%term_num5_de
        gl%des%term_num6_de = gl%de%term_num6_de
        gl%des%coeffde = gl%de%coeffde
        gl%des%texpde = gl%de%texpde
        gl%des%dexpde = gl%de%dexpde
        gl%des%pexpde = gl%de%pexpde

    elseif(flag == 2) then
        gl%de%de_read = gl%dew%de_read
        gl%de%demodel = gl%dew%demodel
        gl%de%tmin_de = gl%dew%tmin_de
        gl%de%tmax_de = gl%dew%tmax_de
        gl%de%pmax_de = gl%dew%pmax_de
        gl%de%rhomax_de = gl%dew%rhomax_de
        gl%de%ref_temp_de = gl%dew%ref_temp_de
        gl%de%ref_dens_de = gl%dew%ref_dens_de
        gl%de%ref_press_de = gl%dew%ref_press_de
        gl%de%term_num1_de = gl%dew%term_num1_de
        gl%de%term_num2_de = gl%dew%term_num2_de
        gl%de%term_num3_de = gl%dew%term_num3_de
        gl%de%term_num4_de = gl%dew%term_num4_de
        gl%de%term_num5_de = gl%dew%term_num5_de
        gl%de%term_num6_de = gl%dew%term_num6_de
        gl%de%coeffde = gl%dew%coeffde
        gl%de%texpde = gl%dew%texpde
        gl%de%dexpde = gl%dew%dexpde
        gl%de%pexpde = gl%dew%pexpde
    elseif(flag==3) then
        gl%de%de_read = gl%des%de_read
        gl%de%demodel = gl%des%demodel
        gl%de%tmin_de = gl%des%tmin_de
        gl%de%tmax_de = gl%des%tmax_de
        gl%de%pmax_de = gl%des%pmax_de
        gl%de%rhomax_de = gl%des%rhomax_de
        gl%de%ref_temp_de = gl%des%ref_temp_de
        gl%de%ref_dens_de = gl%des%ref_dens_de
        gl%de%ref_press_de = gl%des%ref_press_de
        gl%de%term_num1_de = gl%des%term_num1_de
        gl%de%term_num2_de = gl%des%term_num2_de
        gl%de%term_num3_de = gl%des%term_num3_de
        gl%de%term_num4_de = gl%des%term_num4_de
        gl%de%term_num5_de = gl%des%term_num5_de
        gl%de%term_num6_de = gl%des%term_num6_de
        gl%de%coeffde = gl%des%coeffde
        gl%de%texpde = gl%des%texpde
        gl%de%dexpde = gl%des%dexpde
        gl%de%pexpde = gl%des%pexpde
    end if
    end subroutine dielectric_swap
    !--------------------------------------------------------------------------
    !**********************************************************************************


    !---------------------------------------------------------------------------------
    !Initialises the parameters for the electrolyte model
    subroutine initialize_Electrolytes(gl) !für allcoate_subtypees, dahin kopieren
    !--------------------------------------------------
    implicit none
    !type(electrolytes_vars) :: el
    type(type_gl) :: gl
    character(255) :: path_el
    character(16), dimension(100) :: params_read_str
    integer :: i, index
    double precision, dimension(100) :: params_read

    !call initialize(gl)
    if(.not. allocated(gl%dew)) allocate(gl%dew)
    if (.not. allocated(gl%des)) allocate(gl%des)
    !constants
    gl%el%pi = 3.14159265d0
    gl%el%e_const = 1.602176634d-19
    gl%el%el_permit = 8.8541878128d-12 !Vacuum permitivity
    gl%el%Na = 6.02214076d23
    gl%el%kb = 1.380649d-23!Boltzmann
    gl%el%mixedelectrolytes = .false.
    gl%el%brine = .true.

    gl%el%solpos = 0
    gl%el%n_salts = 0
    gl%el%salt = 0
    gl%el%salts = ' '
    gl%el%z_salt(:) = 0.d0
    gl%el%tred = 0.d0
    gl%el%tred2 = 0.d0
    gl%el%tr(:) = 0.d0
    gl%el%pr(:) = 0.d0
    gl%el%pred = 1.d0 !10.d0
    gl%el%pred2 = 0.d0
    gl%el%mapping(:) = 0
    gl%el%resorted = .false.
    gl%el%b = 1.2d0
    gl%el%alpha1 = 2.d0
    gl%el%el_permit = 8.8541878128d-12
    gl%el%ref_molality = 0.d0
    gl%el%g0_param = 0.d0

    !hardcoeded for first tests with nacl
    gl%el%n_salts = 1
    gl%el%z_salt(1) = 1
    gl%el%zm = -1.d0
    gl%el%salt = 1
    gl%el%vm = 1.d0
    gl%el%vx = 1.d0
    gl%el%zm = 1.d0
    gl%el%zx = 1.d0
    gl%el%wm_salt = 0.05844277d0

    gl%el%m_ref = 6.d0


    gl%el%bijk(1:4,1:3,1:7) = 0.d0


    gl%seacalc = .false.
    gl%modelflag = 2
    gl%el%modelflag = 2

    ! the following ones should be set to 2
    gl%el%betaflag = 2
    gl%el%functionflag = 2!1 lit params Pitzer for testing; 2: fitted params


    gl%el%eq_type_el = 0

    !gl%el%g0(1) = 42.d0 !Gibbs Energy of Salt in Standard State
    !gl%el%g0_water = -52.2583216989442d0  !(G(water, 293.15K, 0.101325 MPa))

    !call initialize_Saltparams(gl%el)

    !  if (.not. gl%el%fitting) then
    !path_el = trim(gl%path)//'brines'//'nacl.fld'
    !open (2, file = path_el,  status = 'old', action='read', iostat = index)
    !
    !do i = 1,100
    !    if (index /=0) exit
    !    read(2,*) params_read_str(i)
    !    if (trim(params_read_str(i)) =='@end') exit
    !    read(2,*) params_read(i)
    !end do


    !      !NaCL;
    if(gl%el%functionflag ==2) then
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gl%el%pitzer_param(:,:) = 0.d0     !%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gl%el%g0_param(1,1) =	0.0000000000d+00
        gl%el%g0_param(1,2) =	1.7292128557d+00
        gl%el%g0_param(1,3) =	-5.5243652443d-03
        gl%el%g0_param(1,4) =	6.1892782279d-06
        gl%el%g0_param(1,5) =	0.0000000000d+00
        gl%el%g0_param(1,6) =	-3.2538942207d-04
        gl%el%g0_param(1,7) =	4.7502666930d-07
        gl%el%g0_param(1,8) =	7.4820007582d-04
        gl%el%g0_param(1,9) =	-2.0876532494d-06
        gl%el%cp0_param(1) =	5.0687618589d+02
        gl%el%cp0_param(2) =	2.1771592349d+00
        gl%el%cp0_param(3) =	-4.5843046020d-03
        gl%el%pitzer_param(1,1) =	0.0000000000d+00
        gl%el%pitzer_param(1,2) = 	0.0000000000d+00
        gl%el%pitzer_param(1,3) = 	-9.8837557195d-04
        gl%el%pitzer_param(1,4) = 	-2.4883715573d-05
        gl%el%pitzer_param(1,5) = 	5.4655833666d-07
        gl%el%pitzer_param(1,6) = 	-3.8988471012d-01
        gl%el%pitzer_param(1,7) = 	1.2004575621d-02
        gl%el%pitzer_param(1,8) = 	4.1436390055d-06
        gl%el%pitzer_param(1,9) = 	1.0497822524d-07
        gl%el%pitzer_param(1,10) =	-6.3136082504d-09
        gl%el%pitzer_param(1,11) =	-1.4929938605d-05
        gl%el%pitzer_param(1,12) =	-5.4777559831d-09
        gl%el%pitzer_param(1,13) =	-1.1032477478d-10
        gl%el%pitzer_param(1,14) =	2.7452521144d+00
        gl%el%pitzer_param(1,15) =	2.5977287106d-02
        gl%el%pitzer_param(1,16) =	2.1987584555d-04
        gl%el%pitzer_param(1,17) =	-1.4381675814d-05
        gl%el%pitzer_param(1,18) =	0.0000000000d+00
        gl%el%pitzer_param(1,19) =	0.0000000000d+00
        gl%el%pitzer_param(1,20) =	0.0000000000d+00
        gl%el%pitzer_param(1,21) =	5.8645385261d-04
        gl%el%pitzer_param(2,1) = 	0.0000000000d+00
        gl%el%pitzer_param(2,2) = 	4.9792052035d-01
        gl%el%pitzer_param(2,7) = 	1.2763010540d-05
        gl%el%pitzer_param(2,14) =	-4.4647592597d+00
        gl%el%pitzer_param(4,1) = 	0.0000000000d+00
        gl%el%pitzer_param(4,2) = 	-4.1794074615d+00
        gl%el%pitzer_param(4,3) = 	9.2859065193d-04
        gl%el%pitzer_param(4,6) = 	1.0367844997d+00
        gl%el%pitzer_param(4,7) = 	-7.9694527602d-03
        gl%el%pitzer_param(4,8) = 	-5.0177652463d-06
        gl%el%pitzer_param(4,11) =	6.6796152160d-06
        gl%el%pitzer_param(4,12) =	7.2968037320d-09
        gl%el%pitzer_param(4,14) =	-5.4509584252d-01
        gl%el%pitzer_param(4,15) =	-6.7064023011d-03
        gl%el%pitzer_param(4,18) =	2.5114732220d+01
        gl%el%pitzer_param(4,19) =	0.0000000000d+00



        !*******  REF STATE PARAMS ****************
        gl%el%g0_param(1,10) = -393.1d3 !DGf0 (J/mol)

        gl%el%g0_param(1,11) = 115.05d0  !S0 J/mol K

        !*******************************************
    elseif(gl%el%functionflag ==1) then

        gl%el%lit_param = 0.d0
      
    end if





    !initialize dewater
    gl%dew%de_read = .false.
    gl%dew%demodel = ' '
    gl%dew%tmin_de = 0.d0
    gl%dew%tmax_de = 0.d0
    gl%dew%pmax_de = 0.d0
    gl%dew%rhomax_de = 0.d0
    gl%dew%ref_temp_de = 0.d0
    gl%dew%ref_dens_de = 0.d0
    gl%dew%ref_press_de = 0.d0
    gl%dew%term_num1_de = 0
    gl%dew%term_num2_de = 0
    gl%dew%term_num3_de = 0
    gl%dew%term_num4_de = 0
    gl%dew%term_num5_de = 0
    gl%dew%term_num6_de = 0
    gl%dew%coeffde = 0.d0
    gl%dew%texpde = 0.d0
    gl%dew%dexpde = 0.d0
    gl%dew%pexpde = 0.d0

    !water params hardcoded here!

    gl%dew%de_read = .true.
    gl%dew%demodel(1) = 'DE2'
    gl%dew%tmin_de(1) = 273.16d0
    gl%dew%tmax_de = 13500.d0
    gl%dew%pmax_de = 0.d0
    gl%dew%tmin_de(1) = 273.16d0
    gl%dew%tmax_de(1) = 13500.d0
    gl%dew%pmax_de = 0.d0
    gl%dew%rhomax_de = 0.d0
    gl%dew%ref_temp_de(1) = 647.096d0
    gl%dew%ref_dens_de(1) = 17.8737279956d0
    gl%dew%ref_press_de(1) = 1.d0
    gl%dew%term_num1_de(1) = 11
    gl%dew%term_num2_de(1) = 1
    gl%dew%term_num3_de = 0
    gl%dew%term_num4_de = 0
    gl%dew%term_num5_de = 0
    gl%dew%term_num6_de = 0
    gl%dew%coeffde(1,1) = 0.978224486826d0
    gl%dew%coeffde(1,2) = -0.957771379375d0
    gl%dew%coeffde(1,3) = 0.237511794148d0
    gl%dew%coeffde(1,4) = 0.714692244396d0
    gl%dew%coeffde(1,5) =-0.298217036956d0
    gl%dew%coeffde(1,6) =-0.108863472196d0
    gl%dew%coeffde(1,7) = 0.0949327488264d0
    gl%dew%coeffde(1,8) =-0.00980469816509d0
    gl%dew%coeffde(1,9) = 0.16516763497d-4
    gl%dew%coeffde(1,10) = 0.937359795772d-4
    gl%dew%coeffde(1,11) =-0.12317921872d-9
    gl%dew%coeffde(1,12) = 0.00196096504426d0

    gl%dew%texpde(1,1)= 0.25d0
    gl%dew%texpde(1,2)= 1.0d0
    gl%dew%texpde(1,3)= 2.5d0
    gl%dew%texpde(1,4)= 1.5d0
    gl%dew%texpde(1,5)= 1.5d0
    gl%dew%texpde(1,6)= 2.5d0
    gl%dew%texpde(1,7)= 2.0d0
    gl%dew%texpde(1,8)= 2.0d0
    gl%dew%texpde(1,9)= 5.0d0
    gl%dew%texpde(1,10)= 0.5d0
    gl%dew%texpde(1,11)= 10.0d0
    gl%dew%texpde(1,12)= 228.0d0

    gl%dew%dexpde(1,1) = 1.0d0
    gl%dew%dexpde(1,2) = 1.0d0
    gl%dew%dexpde(1,3) = 1.0d0
    gl%dew%dexpde(1,4) = 2.0d0
    gl%dew%dexpde(1,5) = 3.0d0
    gl%dew%dexpde(1,6) = 3.0d0
    gl%dew%dexpde(1,7) = 4.0d0
    gl%dew%dexpde(1,8) = 5.0d0
    gl%dew%dexpde(1,9) = 6.0d0
    gl%dew%dexpde(1,10) = 7.0d0
    gl%dew%dexpde(1,11) =10.0d0
    gl%dew%dexpde(1,12) = 1.0d0

    gl%dew%pexpde(1,12) = 1.2d0
    !
    !
    !initialize desave
    gl%des%de_read = .false.
    gl%des%demodel = ' '
    gl%des%tmin_de = 0.d0
    gl%des%tmax_de = 0.d0
    gl%des%pmax_de = 0.d0
    gl%des%rhomax_de = 0.d0
    gl%des%ref_temp_de = 0.d0
    gl%des%ref_dens_de = 0.d0
    gl%des%ref_press_de = 0.d0
    gl%des%term_num1_de = 0
    gl%des%term_num2_de = 0
    gl%des%term_num3_de = 0
    gl%des%term_num4_de = 0
    gl%des%term_num5_de = 0
    gl%des%term_num6_de = 0
    gl%des%coeffde = 0.d0
    gl%des%texpde = 0.d0
    gl%des%dexpde = 0.d0
    gl%des%pexpde = 0.d0


    end subroutine
    !----------------------------------------------------------
    !***********************************************************************************************************



    !-------------------------------------------------------------------
    !subroutine read_el_params(gl)
    !!----------------------
    !implicit none
    !type(type_gl) :: gl
    !character(255) :: path_el
    !character(16), dimension(100) :: params_read_str
    !integer :: i, index
    !double precision, dimension(100) :: params_read
    !
    !path_el = trim(gl%path)//'brines\'//'nacl.fld'
    !open (2, file = path_el,  status = 'old', action='read', iostat = index)
    !params_read = 0.d0
    !do i = 1,100
    !    if (index /=0) exit
    !    read(2,*) params_read_str(i)
    !    if (trim(params_read_str(i)) =='@end') exit
    !    read(2,*) params_read(i)
    !end do
    !close(2)
    !
    !!write params to variables here!!
    !gl%el%g0_param(1,1:9) = params_read(1:9)
    !gl%el%cp0_param(1:3) = params_read(10:12)
    !gl%el%pitzer_param(1,1:21) = params_read(13:33)
    !gl%el%pitzer_param(2,1:2) = params_read(34:35)
    !gl%el%pitzer_param(2,7) = params_read(36)
    !gl%el%pitzer_param(2,14) = params_read(37)
    !gl%el%pitzer_param(4,1:3) = params_read(38:40)
    !gl%el%pitzer_param(4,6:8) = params_read(41:43)
    !gl%el%pitzer_param(4,14:15) = params_read(44:45)
    !gl%el%pitzer_param(4,18:19) = params_read(46:47)
    !
    !gl%el%g0_param(1,10:11) = params_read(48:49)
    !end subroutine read_el_params
    !*******************************************************************************
    !******************************************************************************


    !-------------------------------------------------------------------------------------
    !
    subroutine initialize_dew(gl)
    implicit none
    type(type_gl) :: gl
    gl%dew%de_read = .true.
    gl%dew%demodel(1) = 'DE2'
    gl%dew%tmin_de(1) = 273.16d0
    gl%dew%tmax_de = 13500.d0
    gl%dew%pmax_de = 0.d0
    gl%dew%tmin_de(1) = 273.16d0
    gl%dew%tmax_de(1) = 13500.d0
    gl%dew%pmax_de = 0.d0
    gl%dew%rhomax_de = 0.d0
    gl%dew%ref_temp_de(1) = 647.096d0
    gl%dew%ref_dens_de(1) = 17.8737279956d0
    gl%dew%ref_press_de(1) = 1.d0
    gl%dew%term_num1_de(1) = 11
    gl%dew%term_num2_de(1) = 1
    gl%dew%term_num3_de = 0
    gl%dew%term_num4_de = 0
    gl%dew%term_num5_de = 0
    gl%dew%term_num6_de = 0
    gl%dew%coeffde(1,1) = 0.978224486826d0
    gl%dew%coeffde(1,2) = -0.957771379375d0
    gl%dew%coeffde(1,3) = 0.237511794148d0
    gl%dew%coeffde(1,4) = 0.714692244396d0
    gl%dew%coeffde(1,5) =-0.298217036956d0
    gl%dew%coeffde(1,6) =-0.108863472196d0
    gl%dew%coeffde(1,7) = 0.0949327488264d0
    gl%dew%coeffde(1,8) =-0.00980469816509d0
    gl%dew%coeffde(1,9) = 0.16516763497d-4
    gl%dew%coeffde(1,10) = 0.937359795772d-4
    gl%dew%coeffde(1,11) =-0.12317921872d-9
    gl%dew%coeffde(1,12) = 0.00196096504426d0

    gl%dew%texpde(1,1)= 0.25d0
    gl%dew%texpde(1,2)= 1.0d0
    gl%dew%texpde(1,3)= 2.5d0
    gl%dew%texpde(1,4)= 1.5d0
    gl%dew%texpde(1,5)= 1.5d0
    gl%dew%texpde(1,6)= 2.5d0
    gl%dew%texpde(1,7)= 2.0d0
    gl%dew%texpde(1,8)= 2.0d0
    gl%dew%texpde(1,9)= 5.0d0
    gl%dew%texpde(1,10)= 0.5d0
    gl%dew%texpde(1,11)= 10.0d0
    gl%dew%texpde(1,12)= 228.0d0

    gl%dew%dexpde(1,1) = 1.0d0
    gl%dew%dexpde(1,2) = 1.0d0
    gl%dew%dexpde(1,3) = 1.0d0
    gl%dew%dexpde(1,4) = 2.0d0
    gl%dew%dexpde(1,5) = 3.0d0
    gl%dew%dexpde(1,6) = 3.0d0
    gl%dew%dexpde(1,7) = 4.0d0
    gl%dew%dexpde(1,8) = 5.0d0
    gl%dew%dexpde(1,9) = 6.0d0
    gl%dew%dexpde(1,10) = 7.0d0
    gl%dew%dexpde(1,11) =10.0d0
    gl%dew%dexpde(1,12) = 1.0d0

    gl%dew%pexpde(1,12) = 1.2d0
    end subroutine
    !----------------------------------------------------------------------------------------------
    !***********************************************************************************************

    end module electrolytes