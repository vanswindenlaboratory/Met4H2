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
    !ATAN
    !You should have received a copy of the GNU General Public License
    !along with this program.  If not, see < http://www.gnu.org/licenses/ >.
    !
    !*************************************************************************************



    SUBMODULE (calc_functions) impl

    use module_all_types
    use ancillary_equations_mix_module
    use flash_pure_module
    use setup_module
    use seawater_module
    use fniderivs_module
    use gibbsderivs_module
    use vle_derivs_module
    use reduced_parameters_calc_module
    use module_regula_falsi
    use module_regula_falsi_support


    contains
    !*************************************************************************************
    !				TREND Version 3.0
    !		   Thermodynamic Reference & Engineering Data
    !
    !- software for the calculation of thermodynamic and other properties -
    !
    !Copyright (C) 2016,  Prof. Dr.-Ing. R.Span
    !                     Lehrstuhl fuer Thermodynamik
    !                     Ruhr-Universitaet Bochum
    !                     Universitaetsstr. 150
    !                     D-44892 Bochum
    !
    !Cite as: Span, R.; Eckermann, T.; Herrig, S.; Hielscher, S.; Jäger, A.; Thol, M. (2016):
    !	     TREND. Thermodynamic Reference and Engineering Data 3.0. Lehrstuhl fuer
    !	     Thermodynamik, Ruhr-Universitaet Bochum.
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




    !-----------------------------------------
    !Function for the acentric factor
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION OMEGA_CALC(gl,nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: PSAT,TSAT,RhoV,RhoL
    INTEGER::errval,nrsubst,iter,iFlash

    iFlash = 1
    TSAT = 0.D0

    TSAT = 0.7D0*gl%TRED(nrsubst)

    CALL Flash_Pure_PhaseBoundary (gl,PSAT, Tsat, RhoV, RhoL, iFlash, errval, iter, nrsubst)

    OMEGA_CALC = 0.D0

    if( errval == 0 ) then

        OMEGA_CALC = -DLOG10(PSAT/gl%PC(nrsubst))-1.D0

    else

        OMEGA_CALC = gl%accen(nrsubst)

    end if

    END FUNCTION OMEGA_CALC




    !----------------------------------------------
    ! PURE FLUID SATURATED VAPOR ENTHALPY [J/mol]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION HSAT_CALC(gl,TSAT, IPHASE, nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: TSAT, PSAT, RHOV, RHOL, ERROR1, ERROR2
    integer:: ERR_PSAT, IPHASE, nrsubst, iter

    ERROR1 = -444.D0    ! ERROR: TEMPERATURE OUTSIDE THE TWO PHASE REGION
    ERROR2 = -999.D0    ! ERROR: TRIED TO CALL PURE FLUID ROUTINE WITH MIXTURE PARAMETER

    PSAT = 0.D0
    RHOV = 0.D0
    RHOL = 0.D0
    ERR_PSAT = 0
    iter = 0


    call Flash_Pure_PhaseBoundary(gl,psat, tsat, RHOV, RHOL, 1, err_Psat, iter, nrsubst)
    !CALL VLEPURE(TSAT, PSAT, RHOV, RHOL, ERR_PSAT, nrsubst)

    IF (ERR_PSAT > 0) THEN  ! VLE CALCULATION NOT SUCCEEDED
        HSAT_CALC = ERROR1
        RETURN
    END IF

    IF (nrsubst > 0) THEN
        IF (IPHASE == -1) THEN       ! IPHASE -1: SATURATED LIQUID
            HSAT_CALC = H_CALC(gl,TSAT, RHOL, nrsubst)
        ELSEIF (IPHASE == -2) THEN   ! IPHASE -2: SATURATED VAPOR
            HSAT_CALC = H_CALC(gl,TSAT, RHOV, nrsubst)
        END IF
    ELSE
        HSAT_CALC = ERROR2
    END IF

    END FUNCTION HSAT_CALC





    !-------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DSPIN_CALC (gl,T,spinprop,nrsubst) !mol/m³
    !-------------------------------------------------------------------------
    !Function calculates the spinodales with respect to the temperature or pressure
    !Input parameters:
    !- prop: if sprinprop = 1: "liq" = with respect to the temperature
    !        if sprinprop = 2: "vap" = with respect to the temperature
    !- nrsubst: number of fluid
    !M. Kluge, M. Thol, T. Eckermann 2015/01
    !-------------------------------------------------------------------------

    implicit none

    type(type_gl) :: gl


    double precision:: T, rhomax, slope, curvature, factorpress_inv, rhoold, rhomin, p, dvap, dliq, dstart, dend
    integer:: nrsubst, i, max_iterations, errDspin, iterations, errtemp, errval, iter, nsteps, spinprop
    double precision:: dstart_min, dstart_max, delta_allowed, Dmin_allowed, Dmax_allowed, Dspinfound
    double precision:: tstart_min, tstart_max, tmin_allowed,tmax_allowed, Tempfound, tliq, step, rhoright, rholeft
    type(type_additional_parameters) :: ResDspin_param

    factorpress_inv = 1.d0/gl%factorpress
    errval=0
    DSPIN_CALC = 0.d0
    max_iterations=50
    errDspin=0
    iterations=0
    !ResDspin_param = 0.d0
    rhomax = 0.d0
    rhomin = 0.d0

    !with respect to the temperature
    if(nrsubst==0) then
        DSPIN_CALC = -5668d0
        return
    else
        if (spinprop == 1) then     ! liquid phase, temperature as input parameter

            !start with staring values of ancillaries - if available
            rhomax = dl_eq(gl,T, nrsubst)

            ! choose a very high density, where the fluid is definitly a liquid if no starting value is given
            if(rhomax .lt. 1.d14) then
                rhomax = gl%rhoredmix*4.d0
                if (gl%rhoredmix == 1.d0) rhomax = 4.d0 * maxval(gl%rhoc)
            end if

            !reduce density rhomax until the the slope is negative or the curvature is positive (minimal turning point)
            do i = 1, 50 !prevent endless running!--------------------------
                rhoold=rhomax
                rhomax = rhomax*0.95d0
                slope = DPDD_CALC(gl,T, rhomax, nrsubst)*1d6*factorpress_inv

                if (slope < 0.d0) exit
            end do

            if (i > 50) return

            !search density for current temperature where slope is zero
            dstart_min=rhomax
            dstart_max=rhoold
            delta_allowed=1.d-8
            Dmin_allowed=rhomax*0.9d0
            Dmax_allowed=rhoold*1.1d0
            !ResDspin_param=0.d0
            ResDspin_param%a_p(2)=T     !start value
            ResDspin_param%a_p(3)=nrsubst

            CALL Regula_Falsi(gl,Res_DSPIN,Dspinfound,dstart_min,dstart_max,Delta_allowed,&
                Dmin_allowed,Dmax_allowed, Max_iterations,Iterations, errDspin,ResDspin_param)

            if (.not. gl%startvaluespin) then
                ! calculate the saturated liquid and vapor density/temperature
                p = 0.d0; dvap = 0.d0; dliq = 0.d0
                call Flash_Pure_PhaseBoundary(gl,p, T, dvap, dliq, 1, errval, iter, nrsubst)

                dstart = dliq
                if (Dspinfound > gl%rhoc(nrsubst)) then
                    dend = Dspinfound
                else
                    dend = gl%rhoc(nrsubst)
                end if

                !check if there's a saddle point
                do

                    curvature = D2PDD2_CALC (gl,T, dstart,nrsubst)*1d6*factorpress_inv
                    if (curvature < 0)then
                        DSPIN_CALC=-5667.d0  !saddle point found
                        return
                    end if
                    dstart = dstart*0.995d0
                    if (dstart <= dend) exit
                end do
                continue
            else
                gl%startvaluespin = .false.
            end if

        elseif (spinprop == 2) then     !vapor phase, temperature as input parameter

            !start with staring values of ancillaries - if available
            rhomin = dv_eq(gl,T, nrsubst)

            ! choose a very low density, where the fluid is definitly a vapor if no starting value is given
            if (rhomin .lt. 1.d-14) then
                !rhomin = 1.d-14
                rhomin = 1.d-6  !Robin and Moni, 2018/04 the upper value is way too small
            end if

            ! increase density rhomin until the the slope is negative or the curvature is positive (maximal turning point)
            do i = 1, 500
                rhoold=rhomin
                rhomin = rhomin*1.2d0       !Robin and Moni, 2018/04 changed from 10% to 20%
                slope = DPDD_CALC(gl,T, rhomin, nrsubst)*1d6*factorpress_inv
                curvature = D2PDD2_CALC (gl,T, rhomin,nrsubst)*1d6*factorpress_inv
                if ((slope < 0.d0)) exit ! .OR. (curvature > 0.d0)) exit
            end do

            if (i > 500) return

            !search density for current temperature where slope is zero
            dstart_min=rhoold
            dstart_max=rhomin
            delta_allowed=1.d-8
            Dmin_allowed=rhoold*0.9d0
            Dmax_allowed=rhomin*1.1d0
            !ResDspin_param=0.d0 !slope has to be zero!
            ResDspin_param%a_p(2)=T !start value
            ResDspin_param%a_p(3)=nrsubst

            CALL Regula_Falsi(gl,Res_DSPIN,Dspinfound,dstart_min,dstart_max,Delta_allowed,&
                Dmin_allowed,Dmax_allowed, Max_iterations,Iterations, errDspin,ResDspin_param)

        end if
    end if
    DSPIN_CALC=Dspinfound
    if (abs(DSPIN_CALC) .lt. 1.d-14) DSPIN_CALC = -5669.d0
    return

    end function























    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !Functions calculate properties from derivatives of EQS
    !Vector for the ideal part
    !       getderi(/ai,aid*d,aidd*d**2,aidt*d*t,ait*t,aitt*t**2/)
    !       1 = is calculated
    !       0 = is not calculated

    !Vector for the residual part
    !       getderr(/ar,ard*d,ardd*d**2,art*t,artt*t,ardt*d*t,ardtt*d*t**2,arddd*d**3/)
    !       1 = is calculated
    !       0 = is not calculated


    !----------------------------------------------------------------------------------
    !Function for calculating pressure [MPa]
    !----------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION P_CALC(gl,T,D, nrsubst)
    !!DEC$ ATTRIBUTES DLLEXPORT :: P_CALC

    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "P_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, fnrder_w              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD, Rmix, D_ARD_sea, p_helm, d_helm, p_helm_save, d_save, d_mixsave, ar_helm, ar_sea, ar_diff !,wmmix
    logical :: seacalcorig

    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    ! d_ard_sea = 0.d0

    !d_mixsave = d

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)


    P_CALC = 0.D0

    P_CALC = (1.d0 + D_ARD ) * D * Rmix * T * 1.D-6 *gl%factorpress


    END FUNCTION P_CALC
    !---------------------------------------------------------------------

    !----------------------------------------------------------------------------------
    !Function for calculating pressure [MPa]
    !----------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION P_CALC_FROM_DER(gl,T,D, nrsubst, FNRDER)
    !!DEC$ ATTRIBUTES DLLEXPORT :: P_CALC

    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD, Rmix, D_ARD_sea, p_helm, d_helm, p_helm_save, d_save, d_mixsave, ar_helm, ar_sea, ar_diff !,wmmix

    !    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    ! d_ard_sea = 0.d0

    !d_mixsave = d

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
    else
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    P_CALC_FROM_DER = (1.d0 + D_ARD ) * D * Rmix * T * 1.D-6 *gl%factorpress
    END FUNCTION P_CALC_FROM_DER
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    !for reactive mixtures to avoid recursive calls
    !Benedikt 12/2018
    !------------------------------------------------------
    !DOUBLE PRECISION module FUNCTION P_CALC_reac(gl,T,D, nrsubst)
    !!!!!!DEC$ ATTRIBUTES DLLEXPORT :: P_CALC


    !

    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !
    !!ASSEMBLY_INTERFACE(NAME = "P_CALC")
    !
    !integer,DIMENSION(nderivs) :: GETDERR
    !DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, fnrder_w              !Vector gives the values of the calculated derivatives
    !integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    !DOUBLE PRECISION :: T,D,D_ARD, Rmix, D_ARD_sea, p_helm !,wmmix
    !logical :: seacalcorig
    !
    !GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    !d_ard_sea = 0.d0
    !
    !
    !if (nrsubst == 0) then
    !    !call wm_mix_calc(wmmix)                 !subroutine calculates molmass of fluidmixture
    !    call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
    !    !Rmix = 8.314472D0
    !    CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    !else
    !    CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    !    !wmmix = wm(nrsubst)
    !    Rmix = gl%REQ(nrsubst)
    !    !Rmix = 8.314472D0
    !end if
    !
    !D_ARD = FNRDER(2)
    !
    !P_CALC_reac = 0.D0
    !
    !P_CALC_reac = (1.d0 + D_ARD) * D * Rmix * T * 1.D-6 *gl%factorpress
    !
    !
    !END FUNCTION P_CALC_reac
    !---------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------
    !Function for calculating compressibilty
    ! Calculates the real gas factor z as a function of temperature and density
    !----------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Z_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl

    integer,DIMENSION(nderivs) :: GETDERR
    integer, dimension(nderivsi) :: getderi

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr              !Vector gives the values of the calculated derivatives
    double precision, dimension(nderivsi) :: cor_deri                   !needed for ge carrier
    integer:: nrsubst, i                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD, ardel, d_mixsave, ar_del_ge, a01_ge!, a01_gibbs !,Z_CALC
    logical :: seasave

    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    Z_CALC = 0.D0

    if(gl%seawater) gl%gepress = gl%sea%seap

    d_mixsave = d
    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    D_ARD = FNRDER(2)

    if( (gl%gecarrier) .and. (nrsubst .eq. 1)) then

        ardel = a01_gibbs(gl,t, d, gl%gepress)
        z_calc = (1.d0 + D_ARD + ardel)

    elseif((gl%gecarrier) .and. (nrsubst .eq. 0)) then

        gl%gemix = .true.
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)
        ardel = a01_gibbs(gl,t, d, gl%gepress) - cor_derr(2)

        Z_CALC = (1.d0 + D_ARD + ardel)
        gl%gemix = .false.
    else
        Z_CALC = (1.D0 + D_ARD)
    end if
    d = d_mixsave
    END FUNCTION Z_CALC


    !----------------------------------------------------
    !Function for calculating derivative of pressure WR D in [MPa/mol*m^3]
    !----------------------------------------------------
    DOUBLE PRECISION module FUNCTION DPDD_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, fnrder_w              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION T,D,D_ARD,DD_ARDD,wmmix, Rmix, dd_ardd_sea, D_ARD_sea, d_mixsave
    logical :: seasave

    GETDERR = (/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    DPDD_CALC = 0.D0
    d_ard_Sea = 0.d0
    dd_ardd_sea = 0.d0
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if
    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)


    DPDD_CALC = (1.D0 + 2.D0 * (D_ARD ) + DD_ARDD ) * Rmix * T * 1.D-6 * gl%factorpress

    END FUNCTION DPDD_CALC
    !**************************************************************************************


    !----------------------------------------------------
    !Function for calculating derivative of pressure WR D in [MPa/mol*m^3]
    !----------------------------------------------------
    DOUBLE PRECISION module FUNCTION DPDD_CALC_FROM_DER(gl,T,D, nrsubst, FNRDER)

    implicit none

    type(type_gl) :: gl

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION T,D,D_ARD,DD_ARDD,wmmix, Rmix


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
    else
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)


    DPDD_CALC_FROM_DER = (1.D0 + 2.D0 * (D_ARD ) + DD_ARDD ) * Rmix * T * 1.D-6 * gl%factorpress

    END FUNCTION DPDD_CALC_FROM_DER
    !**************************************************************************************


    !**************************************************************************************
    DOUBLE PRECISION module FUNCTION D2PDD2_CALC (gl,T, D, nrsubst)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE PRESSURE WITH RESPECT TO DENSITY
    ! AT CONSTANT TEMPERATURE AND MOLE NUMBERS / CONSTANT TEMPERATURE AND MOLAR COMPOSITION AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ddPddD      - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: T, D

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: dar_dd, d2ar_dd2, d3ar_dd3, Rmix
    integer:: nrsubst


    getder_res = (/0,1,1,0,0,0,0,1,0,0,0,0,0,0,0/)
    der_res = 0.D0


    D2PDD2_CALC = 0.D0                             !Initialize return variable
    dar_dd = 0.D0                       !Initialize
    d2ar_dd2 = 0.D0               !Initialize
    d3ar_dd3= 0.D0        !Initialize


    IF (nrsubst == 0) then
        Call MIXDERIVSFNR (gl,T, D, getder_res, der_res)
        Call R_mix_calc(gl,Rmix)

    else
        Call FNRDERIVS (gl,T, D, getder_res, der_res,nrsubst)
        Rmix = gl%REQ(nrsubst)
    END IF

    dar_dd = der_res(2)                 !Derivative of the residual Helmholtz with respect to delta
    d2ar_dd2 = der_res(3)         !Second derivative of the residual Helmholtz with respect to delta
    d3ar_dd3 = der_res(8) !Third derivative of the residual Helmholtz with respect to delta

    D2PDD2_CALC = T *  Rmix / d * (2.d0*dar_dd + 4.d0 * d2ar_dd2 + d3ar_dd3)* 1.D-6 * gl%factorpress
    end function D2PDD2_CALC

    !**************************************************************************************
    DOUBLE PRECISION module FUNCTION D2PDD2_CALC_FROM_DER (gl,T, D, nrsubst, der_res)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE PRESSURE WITH RESPECT TO DENSITY
    ! AT CONSTANT TEMPERATURE AND MOLE NUMBERS / CONSTANT TEMPERATURE AND MOLAR COMPOSITION AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ddPddD      - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: T, D

    double precision, dimension(nderivs):: der_res
    double precision:: dar_dd, d2ar_dd2, d3ar_dd3, Rmix
    integer:: nrsubst


    IF (nrsubst == 0) then
        Call R_mix_calc(gl,Rmix)

    else
        Rmix = gl%REQ(nrsubst)
    END IF

    dar_dd = der_res(2)                 !Derivative of the residual Helmholtz with respect to delta
    d2ar_dd2 = der_res(3)         !Second derivative of the residual Helmholtz with respect to delta
    d3ar_dd3 = der_res(8) !Third derivative of the residual Helmholtz with respect to delta

    D2PDD2_CALC_FROM_DER = T *  Rmix / d * (2.d0*dar_dd + 4.d0 * d2ar_dd2 + d3ar_dd3)* 1.D-6 * gl%factorpress
    end function D2PDD2_CALC_FROM_DER

    !--------------------------------------------------------------------
    !Function for calculating derivative of pressure WR T at constant RHO
    !--------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DPDT_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD,DT_ARDT,wmmix, Rmix

    GETDERR = (/0,1,0,0,0,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DT_ARDT = FNRDER(6)

    DPDT_CALC = 0.D0
    DPDT_CALC = (1.D0 + D_ARD - DT_ARDT) * Rmix * D * 1.D-6 * gl%factorpress

    END FUNCTION DPDT_CALC


    !------------------------------------------------------------------
    !Function for calculating derivative of density WR T in [mol*m^3/T]
    !------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DDDT_CALC(gl,T,D, nrsubst)



    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: T, D, DPDT, DPDD, rmix, wmmix
    INTEGER:: nrsubst

    DDDT_CALC = 0.D0

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
    else
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    DPDT=DPDT_CALC(gl,T,D, nrsubst)
    DPDD=DPDD_CALC(gl,T,D, nrsubst)

    DDDT_CALC = -DPDT/DPDD

    END FUNCTION DDDT_CALC


    !--------------------------------------------------------------------
    !Function for calculating 2nd derivative of pressure WR T^2 [MPa/K^2]
    !--------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION D2PDTT_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "D2PDTT_CALC_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TTD_ARTTD,wmmix, Rmix

    GETDERR=(/0,0,0,0,0,0,1,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),residual

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    TTD_ARTTD = FNRDER(7)

    D2PDTT_CALC = 0.D0
    D2PDTT_CALC = TTD_ARTTD * Rmix * D * 1.d-6 * gl%factorpress / T

    END FUNCTION D2PDTT_CALC


    !--------------------------------------------------------------------------
    !Function for calculating mixed derivative of pressure WR T and D [J/MOL-K]
    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION D2PDDT_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "D2PDDT_CALC_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                          !Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                   !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD,DD_ARDD,DT_ARDT,DDT_ARDDT,wmmix, Rmix

    GETDERR = (/0,1,1,0,0,1,0,0,0,1,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),residual

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    DT_ARDT = FNRDER(6)
    DDT_ARDDT = FNRDER(10)

    D2PDDT_CALC = 0.D0
    D2PDDT_CALC = Rmix * (1.d0 + 2.d0 * D_ARD + DD_ARDD - 2.d0 * DT_ARDT - DDT_ARDDT)


    END FUNCTION D2PDDT_CALC



    !---------------------------------------------
    !Function for calculating fugacity coefficient
    !---------------------------------------------
    module Subroutine FUGCO_CALC_MIX(gl,T,D, fugcoef_mix)

    !DEC$ ATTRIBUTES DLLEXPORT :: FUGCO_CALC_MIX







    implicit none

    type(type_gl) :: gl



    double precision :: T, D
    double precision, dimension(30), intent(out):: fugcoef_mix

    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                     !Vector gives the values of the calculated derivatives
    integer:: i, oir
    DOUBLE PRECISION :: D_ARD
    double precision, dimension(30):: chempot_res

    GETDERR = (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    fugcoef_mix = 0.d0

    CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    D_ARD = FNRDER(2)

    oir = 2                    !The residual chemical potential will be calculated (0 = overall, 1 = ideal)
    call dna_dni(gl,T, D, chempot_res, OIR)

    Do i = 1, gl%NCOMP
        if ((1.D0 + D_ARD) < 0.d0) then
            fugcoef_mix(i) = 1.d-10
        else
            if ((chempot_res(i) - dLOG(1.D0 + D_ARD)) < 700.d0) then
                fugcoef_mix(i) = dEXP(chempot_res(i) - dLOG(1.D0 + D_ARD))
            else ! unreasonable high values, dexp evaluation would crash
                fugcoef_mix(i) = 1.d+50
            end if
        end if
        if (fugcoef_mix(i) == 0.d0) fugcoef_mix(i) = 1.d-10
    End Do


    END Subroutine FUGCO_CALC_MIX

    !---------------------------------------------
    !Function for calculating fugacity coefficient
    !---------------------------------------------
    DOUBLE PRECISION module FUNCTION FUGCOPURE_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl



    double precision :: T, D

    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst
    DOUBLE PRECISION :: AR,D_ARD


    GETDERR = (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part

    AR = FNRDER(1)
    D_ARD = FNRDER(2)

    ! safety caution, exp(710) would crash
    if ((AR + D_ARD) > 700.d0) then
        !Andreas Oktober 2012
        !Changed this, since for mixtures in some calculations only the fugacity coefficient of a pure compound in the mixture is nedded
        !If this is the case, the other substances may have unreasonable fugacity coefficients, but the routine may not return here!!!
        FUGCOPURE_CALC = 1.d+100
        !errval = -1444
        !return
    else if (isnan(AR)) then   !((AR + D_ARD) == NaN) then
        FUGCOPURE_CALC = -9876
        return
    else if (1.D0 + D_ARD == 0.d0) then   ! / 0
        FUGCOPURE_CALC = -9876
        return
    else
        FUGCOPURE_CALC = dEXP(AR + D_ARD) / (1.D0 + D_ARD)
    end if

    END function FUGCOPURE_CALC

    !---------------------------------------------
    !Function for the residual chemical potential
    !---------------------------------------------
    DOUBLE PRECISION module FUNCTION CPOTR_CALC(gl,T,D,nrsubst)






    implicit none

    type(type_gl) :: gl



    double precision :: T, D

    integer:: oir, nrsubst, i
    DOUBLE PRECISION :: CPOTR
    double precision, dimension(30):: chempot_res

    CPOTR_CALC = 0.d0
    oir = 2                    !The residual chemical potential will be calculated (0 = overall, 1 = ideal, 2 = residual)
    chempot_res=0.d0
    call dna_dni(gl,T, D, chempot_res, OIR)

    CPOTR_CALC = 0.d0

    Do i = 1, gl%NCOMP
        CPOTR_CALC = chempot_res(i)
    end do

    END function CPOTR_CALC

    !---------------------------------------------
    !Function for the ideal chemical potential
    !---------------------------------------------
    DOUBLE PRECISION module FUNCTION CPOTI_CALC(gl,T,D,nrsubst)






    implicit none

    type(type_gl) :: gl


    double precision :: T, D

    integer:: oir, nrsubst
    DOUBLE PRECISION :: CPOTI
    double precision, dimension(30):: chempot_res

    CPOTI_CALC = 0.d0
    oir = 1                    !The residual chemical potential will be calculated (0 = overall, 1 = ideal, 2 = residual)
    chempot_res=0.d0
    call dna_dni(gl,T, D, chempot_res, OIR)
    CPOTI = chempot_res(nrsubst)
    CPOTI_CALC = CPOTI

    END function CPOTI_CALC



    DOUBLE PRECISION module FUNCTION DlnFugcoef_pureDD (gl,T, D, i)






    implicit none

    type(type_gl) :: gl


    double precision :: T, D, del
    integer, dimension (nderivs) :: GETDERR
    double precision, dimension (nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    double precision :: D_ARD, DD_ARDD
    integer :: i

    GETDERR = (/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    call FNRDERIVS(gl,T,D,GETDERR,FNRDER,i)       !subroutine calculates derivatives of residual part


    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    del = D/gl%rhored(i)

    DlnFugcoef_pureDD = 1.d0 / gl%rhored(i) / del * (2.d0*D_ARD + DD_ARDD - (D_ARD+DD_ARDD)/(1.d0+D_ARD))


    end function DlnFugcoef_pureDD


    ! DOUBLE PRECISION module FUNCTION DFugcoef_pureDD (T, D, i)
    !


    !
    ! implicit none
    !
    ! double precision :: T, D, del
    ! integer, dimension (nderivs) :: GETDERR = (/1,1,1,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    ! double precision, dimension (nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    ! double precision :: D_ARD, DD_ARDD, AR
    ! integer :: i
    !
    ! call FNRDERIVS(T,D,GETDERR,FNRDER,i)       !subroutine calculates derivatives of residual part
    !
    ! AR = FNRDER(1)
    ! D_ARD = FNRDER(2)
    ! DD_ARDD = FNRDER(3)
    ! del = D/rhored(i)
    !
    ! DFugcoef_pureDD = (1.d0 / rhored(i) / del * (2.d0*D_ARD + DD_ARDD - (D_ARD+DD_ARDD)/(1.d0+D_ARD))) * &
    !                    dEXP(AR + D_ARD)/(1.D0 + D_ARD)
    !
    !
    ! end function

    DOUBLE PRECISION module FUNCTION DFugcoefpureDD_CALC (gl,T, D, nrsubst)






    implicit none

    type(type_gl) :: gl


    double precision :: T, D, del
    integer, dimension (nderivs) :: GETDERR
    double precision, dimension (nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    double precision :: D_ARD, DD_ARDD, AR
    integer :: nrsubst

    GETDERR = (/1,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    call FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)       !subroutine calculates derivatives of residual part

    AR = FNRDER(1)
    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    !del = D/rhored(nrsubst)

    ! safety caution, exp(710) would crash
    if ((AR + D_ARD) > 700.d0) then
        !Andreas Oktober 2012
        !Changed this, since for mixtures in some calculations only the fugacity coefficient of a pure compound in the mixture is nedded
        !If this is the case, the other substances may have unreasonable fugacity coefficients, but the routine may not return here!!!
        DFugcoefpureDD_CALC = 1.d+100
        !errval = -1445
        !return
    else
        !DFugcoefpureDD_CALC = (1.d0 / rhored(nrsubst) / del * (2.d0*D_ARD + DD_ARDD - (D_ARD+DD_ARDD)/(1.d0+D_ARD))) * &
        !                    dEXP(AR + D_ARD)/(1.D0 + D_ARD)
        DFugcoefpureDD_CALC = (1.d0 / D * (2.d0*D_ARD + DD_ARDD - (D_ARD+DD_ARDD)/(1.d0+D_ARD))) * &
            dEXP(AR + D_ARD)/(1.D0 + D_ARD)
    end if

    end function DFugcoefpureDD_CALC


    DOUBLE PRECISION module FUNCTION DFugcoef_pureDT (gl,T, D, i)






    implicit none

    type(type_gl) :: gl


    double precision :: T, D, tau
    integer, dimension (nderivs) :: GETDERR
    double precision, dimension (nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    double precision :: D_ARD, T_ART, DT_ARDT, AR
    integer :: i

    GETDERR = (/1,1,0,1,0,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    call FNRDERIVS(gl,T,D,GETDERR,FNRDER,i)       !subroutine calculates derivatives of residual part

    AR = FNRDER(1)
    D_ARD = FNRDER(2)
    T_ART = FNRDER(4)
    DT_ARDT = FNRDER(6)
    tau = gl%Tred(i)/T

    ! safety caution, exp(710) would crash
    if ((AR + D_ARD) > 700.d0) then
        !Andreas Oktober 2012
        !Only changed this else statement, since this routine is already ok for mixtures
        DFugcoef_pureDT = 1.d+100
        !errval = -1446
        !        return
    else
        DFugcoef_pureDT = (-1.d0 / T * (T_ART + DT_ARDT - (DT_ARDT)/(1.d0+D_ARD))) * &
            dEXP(AR + D_ARD)/(1.D0 + D_ARD)
    end if

    end function DFugcoef_pureDT


    !------------------------------------------------------------------------
    !Function for calculating fugacity coefficient as a function of T and P
    !------------------------------------------------------------------------
    module Subroutine FUGCO_CALC_PURE_TP(gl,T,P, fugcoef_pure)







    implicit none

    type(type_gl) :: gl



    double precision :: T, P
    double precision, dimension(30), intent(out)::fugcoef_pure

    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer :: i
    DOUBLE PRECISION :: AR,D_ARD, rhoredmixorg, tredmixorg, D

    GETDERR = (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    fugcoef_pure = 0.D0         !Initialize


    rhoredmixorg = gl%rhoredmix
    tredmixorg = gl%tredmix

    Do i = 1, gl%NCOMP

        !Andreas Dec. 2010. This is necessary, since elsewise tau and del get reduced with the wrong parameters!!
        gl%rhoredmix = gl%rhored(i)
        gl%tredmix = gl%tred(i)

        if (gl%rho_TP_pure(i) > 0) then              !Module Variable rho_TP_Pure. The density of each pure component at T and p is stored in
            D = gl%rho_TP_pure(i)
        else                                   !If the module variable was not set, the densities have to be calculated in this routine
            D = rhomix_calc(gl,T, P, 0.d0, 0, i)
        end if

        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,i)   !subroutine calculates derivatives of residual part

        AR = FNRDER(1)
        D_ARD = FNRDER(2)

        !fugcoef_pure(i) = dEXP(AR + D_ARD - dLOG(1.D0 + D_ARD))
        fugcoef_pure(i) = dEXP(AR + D_ARD)/(1.D0 + D_ARD)

    End Do

    gl%tredmix = tredmixorg
    gl%rhoredmix = rhoredmixorg

    END Subroutine FUGCO_CALC_PURE_TP


    !-----------------------------------------
    !Function for calculating entropy [J/(mol-K)]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION S_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "S_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,DIMENSION(nderivsi) :: getderi
    integer:: nrsubst, i                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr              !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,DIMENSION(nderivsi) :: FNIDER, cor_deri              !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION :: T,D,AI,T_AIT,AR,T_ART,wmmix, Rmix,  a_ge, atau, ms, d_mixsave, a_helm, a_gi, a_tau_ge, a, a01_gibbs!, a_gibbs
    logical :: seasave

    GETDERR=(/1,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    getderi=(/1,0,0,0,1,0,0,0,0,0/)     !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AI = FNIDER(1)
    T_AIT = FNIDER(5)

    AR = FNRDER(1)
    T_ART = FNRDER(4)

    S_CALC = 0.D0
    if(gl%seawater) gl%gepress = gl%sea%seap

    if( (gl%gecarrier) .and. (nrsubst .eq. 1)) then

        a_tau_ge = a10_gibbs(gl, t, d, gl%gepress)
        a = a_gibbs(gl, t, d, gl%gepress)

        S_CALC = ( -a + a_tau_ge) * Rmix
        d = d_mixsave

    elseif((gl%gecarrier) .and. (nrsubst .eq. 0) ) then

        gl%gemix = .true.

        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)

        a_tau_ge = a10_gibbs(gl, t, d, gl%gepress) - cor_deri(5) - cor_derr(4)
        a_ge = a_gibbs(gl, t, d, gl%gepress) - cor_deri(1) - cor_derr(1)
        
        S_CALC = (T_AIT + T_ART - AI - AR + a_tau_ge - a_ge) * Rmix
        d = d_mixsave
        gl%gemix = .false.
    else
        S_CALC = (T_AIT + T_ART - AI - AR) * Rmix
    end if

    end FUNCTION S_CALC



    !-----------------------------------------
    !Function for internal energy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION U_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "U_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst, i                                ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr              !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri              !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,T_AIT,T_ART,wmmix, Rmix, atau, diff, us, u_sea_calc, d_mixsave, a_tau,a01_helm
    !logical :: seasave


    GETDERR=(/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,1,0,0,0,0,0/)     !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    T_AIT = FNIDER(5)
    T_ART = FNRDER(4)

    U_CALC = 0.D0

    if(gl%seawater) gl%gepress = gl%sea%seap

    if( (gl%gecarrier) .and. (nrsubst .eq. 1) ) then

        atau = a10_gibbs(gl,t, d, gl%gepress)
        U_CALC = (atau ) * Rmix * T


    elseif((gl%gecarrier) .and. (nrsubst .eq. 0)) then
        !gl%gemix = .true.
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)

        atau = a10_gibbs(gl,t, d, gl%gepress) - cor_deri(5) - cor_derr(4)

        U_CALC = (T_AIT + T_ART + atau) * Rmix * T
        !gl%gemix = .false.
    else
        U_CALC = (T_AIT + T_ART) * Rmix * T
    end if
    d=d_mixsave

    END FUNCTION U_CALC



    !-------------------------------------------
    !Function for isochoric heat capacity [J/(mol-K)]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION CV_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "CV_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst, i                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr             !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri              !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TT_AITT,TT_ARTT,wmmix, Rmix, att, d_mixsave
    logical :: seasave

    GETDERR=(/0,0,0,0,1,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)     !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    TT_AITT = FNIDER(6)
    TT_ARTT = FNRDER(5)

    CV_CALC = 0.D0
    if(gl%seawater) gl%gepress = gl%sea%seap
    !if((gl%seacalc) .and. (nrsubst .eq. 1)) then
    if( (gl%gecarrier) .and. (nrsubst .eq. 1) ) then

        att = a20_gibbs(gl, t, d, gl%gepress)
        cv_calc = -(TT_AITT + TT_ARTT + att) * Rmix

    elseif( (gl%gecarrier) .and. (nrsubst .eq. 0)) then

        gl%gemix = .true.
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)
        att = a20_gibbs(gl, t, d, gl%gepress) - cor_deri(6) - cor_derr(5)
        CV_CALC = - (att + TT_AITT + TT_ARTT ) * Rmix

        gl%gemix = .false.
    else
        CV_CALC = - (TT_AITT + TT_ARTT ) * Rmix
    end if
    d=d_mixsave

    END FUNCTION CV_CALC



    !-----------------------------------------
    !Function for enthalpy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION H_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "H_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst, i                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,T_AIT,T_ART,D_ARD,wmmix, Rmix
    double precision :: atau, ardel, wm_sea, hs, rho, d_mixsave!, p_calc, p
    logical :: seasave


    GETDERR=(/0,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,1,0,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if
    T_AIT = FNIDER(5)

    D_ARD = FNRDER(2)
    T_ART = FNRDER(4)

    H_CALC = 0.D0
    if(gl%seawater) gl%gepress = gl%sea%seap

    if((gl%gecarrier) .and. (nrsubst .eq. 1) ) then

        atau = a10_gibbs(gl, t, d, gl%gepress)
        ardel = a01_gibbs(gl,t,d, gl%gepress)

        h_calc = (1.d0 + atau + ardel) * Rmix * T

    elseif( (gl%gecarrier) .and. (nrsubst .eq. 0) ) then
        gl%gemix = .true.
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)
        atau = a10_gibbs(gl, t, d, gl%gepress)-cor_deri(5) - cor_derr(4)
        ardel = a01_gibbs(gl,t,d, gl%gepress) - cor_derr(2)

        h_calc = (1.d0 + atau + T_AIT + T_ART + ardel + D_ARD ) * Rmix * T

        gl%gemix = .false.
    else
        h_calc = (1.d0 + T_AIT + T_ART + D_ARD ) * Rmix * T
    end if
    d=d_mixsave

    END FUNCTION H_CALC




    !-----------------------------------------
    !Function for isobaric heat capacity [J/(mol-K)]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION CP_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "CP_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst    , i                             ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TT_AITT,TT_ARTT,D_ARD,DT_ARDT,DD_ARDD,wmmix, Rmix, att, ardel, ardd, ardtdr, cp_sea, cp_saline, d_mixsave, p!, a01_gibbs
    logical :: seasave

    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),ideal
    GETDERR=(/0,1,1,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),residual
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    TT_AITT = FNIDER(6)

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    CP_CALC = 0.D0
    if(gl%seawater) gl%gepress = gl%sea%seap
    ! if((gl%seacalc) .and. (nrsubst .eq. 1))then
    if((gl%gecarrier) .and. (nrsubst .eq. 1) ) then

        att = a20_gibbs(gl, t, d, gl%gepress)
        ardel = a01_gibbs(gl, t, d, gl%gepress)
        ardd = a02_gibbs(gl, t, d, gl%gepress)
        ardtdr = a11_gibbs(gl, t, d, gl%gepress)

        cp_calc =  (-(att) + (1.D0 + ardel - ardtdr) **2.d0 / (1.D0 + 2.D0 * ardel + ardd)) * Rmix

    elseif( (gl%gecarrier) .and. (nrsubst .eq. 0) ) then
        !d_mixsave = d
        !seasave = gl%seacalc

        gl%gemix = .true.
        !gibbsderivs
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)

        att = a20_gibbs(gl, t, d, gl%gepress) - cor_deri(6) - cor_derr(5)
        ardel = a01_gibbs(gl, t, d, gl%gepress) - cor_derr(2)
        ardd = a02_gibbs(gl, t, d, gl%gepress) - cor_derr(3)
        ardtdr = a11_gibbs(gl, t, d, gl%gepress) - cor_derr(6)
        !
        CP_CALC = (-( att + TT_AITT + TT_ARTT ) + (1.D0  + ardel - ardtdr + D_ARD - DT_ARDT) ** 2/(1.D0 + 2.d0*ardel + ardd + 2.D0 * D_ARD + DD_ARDD)) * Rmix

        gl%gemix =.false.
    else
        CP_CALC = (-(TT_AITT + TT_ARTT ) + (1.D0 + D_ARD - DT_ARDT ) ** 2/(1.D0 + 2.D0 * D_ARD + DD_ARDD )) * Rmix
    end if
    d=d_mixsave

    END FUNCTION CP_CALC


    !-------------------------------------------
    !Function for phase identifier parameter [-]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION PIP_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl

    !ASSEMBLY_INTERFACE(NAME = "PIP_CALC")

    DOUBLE PRECISION :: T, D, D2PDDT, DPDT, D2PDD2, DPDD
    ! DOUBLE PRECISION :: D2PDDT_CALC, DPDT_CALC, D2PDD2_CALC, DPDD_CALC
    INTEGER :: nrsubst

    D2PDDT = D2PDDT_CALC(gl,T,D, nrsubst)
    DPDT = DPDT_CALC(gl,T,D, nrsubst) * gl%factor
    D2PDD2 = D2PDD2_CALC(gl,T,D, nrsubst) * gl%factor * gl%factortrans
    DPDD = DPDD_CALC(gl,T,D, nrsubst) * gl%factortrans

    PIP_CALC = 0.D0
    PIP_CALC = 2.d0 - (D2PDDT / DPDT - D2PDD2 / DPDD) * D / gl%factor

    END FUNCTION PIP_CALC

    !-------------------------------------------
    !Function for phase identifier parameter with respect to density [-]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION D_PIP_DD_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst,i
    double precision, dimension(2) :: stenc,dense
    double precision :: num_del

    num_del = 1.D-6
    dense(1) = D + num_del
    dense(2) = D - num_del

    do i=1,2
        stenc(i) =  PIP_CALC(gl,T,dense(i), nrsubst)
    end do

    D_PIP_DD_CALC = (stenc(1)-stenc(2))/(2d0*num_del)

    END FUNCTION D_PIP_DD_CALC



    !-------------------------------------------
    !Function for phase identifier parameter with respect to density^2 [-]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION D2_PIP_D2D_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst,i
    double precision, dimension(2) :: stenc,dense
    double precision :: num_del

    num_del = 1.D-6
    dense(1) = D + num_del
    dense(2) = D - num_del

    do i=1,2
        stenc(i) =  D_PIP_DD_CALC(gl,T,dense(i), nrsubst)
    end do

    D2_PIP_D2D_CALC = (stenc(1)-stenc(2))/(2d0*num_del)

    END FUNCTION D2_PIP_D2D_CALC

    !-------------------------------------------
    !Function for phase identifier parameter with respect to temperature [-]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION D_PIP_DT_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst,i
    double precision, dimension(2) :: stenc,temp
    double precision :: num_del

    num_del = 1.D-6
    temp(1) = T + num_del
    temp(2) = T - num_del

    do i=1,2
        stenc(i) =  PIP_CALC(gl,temp(i),D, nrsubst)
    end do

    D_PIP_DT_CALC = (stenc(1)-stenc(2))/(2d0*num_del)

    END FUNCTION D_PIP_DT_CALC


    !-------------------------------------------
    !Function for Grueneisen coefficient [-]
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION GRUEN_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl

    !ASSEMBLY_INTERFACE(NAME = "GRUEN_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: TT_AITT,TT_ARTT,D_ARD,DT_ARDT, Rmix, test

    DOUBLE PRECISION :: T, D, CV, DPDT
    !DOUBLE PRECISION :: CV_CALC, DPDT_CALC
    INTEGER :: nrsubst

    GETDERR=(/0,1,0,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),residual
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),ideal

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        Rmix = gl%REQ(nrsubst)
    end if

    TT_AITT = FNIDER(6)

    D_ARD = FNRDER(2)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    GRUEN_CALC=(1.d0 + D_ARD - DT_ARDT) / (-1.d0 *(TT_AITT + TT_ARTT))

    END FUNCTION GRUEN_CALC


    !-------------------------------------------
    !Function for Grueneisen coefficient [-]  w/0  cv0
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION GRUEN_WO_CV0_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl

    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION :: TT_ARTT,D_ARD,DT_ARDT, Rmix, test
    DOUBLE PRECISION :: T, D, CV, DPDT
    INTEGER :: nrsubst

    GETDERR=(/0,1,0,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),residual

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if


    D_ARD = FNRDER(2)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    GRUEN_WO_CV0_CALC=(1.d0 + D_ARD - DT_ARDT) / (-1.d0 *(TT_ARTT))

    END FUNCTION GRUEN_WO_CV0_CALC

    !-------------------------------------------
    !Function for Grueneisen coefficient [-]  w/0  cv0
    !-------------------------------------------
    DOUBLE PRECISION module FUNCTION GRUEN_WO_CV0_CALC_DT(gl,T,D, nrsubst)
    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    double precision :: del_num  = 1d-6

    GRUEN_WO_CV0_CALC_DT = (GRUEN_WO_CV0_CALC(gl,T+del_num,D, nrsubst) - GRUEN_WO_CV0_CALC(gl,T-del_num,D, nrsubst))/(2d0*del_num)

    end function



    !-----------------------------------------
    !Function for ideal-Gas isobaric heat capacity [J/(mol-K)]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION CP0_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl


    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TT_AITT,wmmix, Rmix, CV0_CALC

    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1),ideal


    D = 1.D-12

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)!gl%cp0red(nrsubst)
    end if

    TT_AITT = FNIDER(6)

    CP0_CALC = 0.D0
    CV0_CALC = -TT_AITT * Rmix
    CP0_CALC = -(TT_AITT - 1.D0) * Rmix

    END FUNCTION CP0_CALC


    !-----------------------------------------
    !Function for gibbs energy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION G_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "G_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst  , i                               ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr               !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri               !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AI,AR,D_ARD,wmmix, Rmix, ardel, a, salb, freeze_temp, a_spec, molality_sea, molality, d_mixsave, a_test, a_log, d_mixcalc, ar_diff, ar_helm, ar_gi


    logical :: seasave

    GETDERR=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/1,0,0,0,0,0,0,0,0,0/)     !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AI = FNIDER(1)

    AR = FNRDER(1)
    D_ARD = FNRDER(2)

    ar_helm = d_ard

    G_CALC = 0.D0
    if(gl%seawater) gl%gepress = gl%sea%seap

    if( (gl%gecarrier) .and. (nrsubst .eq. 1)) then

        ardel = a01_gibbs(gl, t, d, gl%gepress)
        a = a_gibbs(gl, t, d, gl%gepress)
        g_calc = (1.d0 + a + ardel) * Rmix * T

    elseif( (gl%gecarrier) .and. (nrsubst .eq. 0)) then

        gl%gemix = .true.
        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)
        ardel = a01_gibbs(gl, t, d, gl%gepress) - cor_derr(2)
        a = a_gibbs(gl, t, d, gl%gepress) - cor_deri(1) - cor_derr(1)

        g_calc = (1.d0 + a + ardel + AI + AR + D_ARD ) * Rmix * T

        gl%gemix = .false.
        d = d_mixsave
    else
        G_CALC = (1.D0 + AI + AR + D_ARD) * Rmix * T
    end if
    d=d_mixsave

    END FUNCTION G_CALC


    !-----------------------------------------
    !Function for residual gibbs energy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION GR_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "GR_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER               !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AR,D_ARD,wmmix, Rmix

    GETDERR=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AR = FNRDER(1)
    D_ARD = FNRDER(2)

    GR_CALC = 0.D0
    GR_CALC = (AR + D_ARD) * Rmix * T

    END FUNCTION GR_CALC


    !!-----------------------------------------------------------------------------
    !!Alternative Function for calculating the reduced Gibbs-Energy of mixing [-]. Andreas Jäger, Jan. 2011
    !!The Gibbs-Energy of Mixing is calculated by subtracting the Gibbs-Energy of the pure components at
    !!The SAME Pressure and Temperature from the total Gibbs-Energy
    !!-----------------------------------------------------------------------------
    ! DOUBLE PRECISION module FUNCTION G_MIX_CALC(T,D)
    !



    !
    ! IMPLICIT NONE
    !
    !   DOUBLE PRECISION :: T,D, wmmix, P, P_CALC
    !   double precision, dimension(30):: fugcoeff_mix, fugcoeff_pure
    !   integer:: i
    !
    !   Call FUGCO_CALC_MIX(T,D, fugcoeff_mix)
    !
    !   ! check if pure fluid fugacity coefficients are already available as mudule variables
    !   ! the module variable "fugco_pure" from the module "module_VLE" is given a value in
    !   ! the routine "newvars" in the flash algorithm package
    !   if (fugco_pure(1) == 0.D0) then
    !        P = P_CALC(gl,T, D, 0)
    !        Call FUGCO_CALC_PURE_TP(T,P, fugcoeff_pure)
    !   else
    !        fugcoeff_pure = fugco_pure
    !   end if
    !
    !   call wm_mix_calc(wmmix)      ! subroutine calculates molmass of fluidmixture
    !
    !   G_MIX_CALC = 0.d0 !Initialize Gmix
    !
    !   Do i =1, NCOMP - 1
    !
    !    G_MIX_CALC = G_MIX_CALC + molfractions(i)* dlog(molfractions(i) * fugcoeff_mix(i) * &
    !            & fugcoeff_pure(NCOMP) / (molfractions(NCOMP) * fugcoeff_mix(NCOMP) * fugcoeff_pure(i)))
    !
    !   End do
    !
    !   G_MIX_CALC = G_MIX_CALC + dlog(molfractions(NCOMP)* fugcoeff_mix(NCOMP) /  fugcoeff_pure(NCOMP))
    !
    ! RETURN
    ! END
    !
    !!-----------------------------------------------------------------------------
    !!Alternative Function for calculating the reduced Excess Gibbs-Energy  [-]. Andreas Jäger, Jan. 2011
    !!-----------------------------------------------------------------------------
    ! DOUBLE PRECISION module FUNCTION G_Excess_CALC(gl,T,D)
    !


    !
    ! IMPLICIT NONE
    !
    !   DOUBLE PRECISION :: T,D, wmmix, P, P_CALC
    !   double precision, dimension(30):: fugcoeff_mix, fugcoeff_pure
    !   integer:: i
    !
    !
    !   Call FUGCO_CALC_MIX(T,D, fugcoeff_mix)
    !
    !   P = P_CALC(gl,T, D, 0)
    !
    !   Call FUGCO_CALC_PURE_TP(T,P, fugcoeff_pure)
    !
    !   call wm_mix_calc(wmmix)      ! subroutine calculates molmass of fluidmixture
    !
    !   Do i =1, NCOMP - 1
    !
    !    G_Excess_CALC = G_Excess_CALC + molfractions(i)* dlog(fugcoeff_mix(i) * &
    !            & fugcoeff_pure(NCOMP) / (fugcoeff_mix(NCOMP) * fugcoeff_pure(i)))
    !
    !   End do
    !
    !   G_Excess_CALC = G_Excess_CALC + dlog(fugcoeff_mix(NCOMP) /  fugcoeff_pure(NCOMP))
    !
    ! RETURN
    ! END FUNCTION G_Excess_CALC

    !-----------------------------------------
    !Function for the chemical potential [J/mol]
    !-----------------------------------------
    module Subroutine Chempot_CALC(gl,T,D, Chempot, oir)
    !!DEC$ ATTRIBUTES DLLEXPORT :: Chempot_CALC


    implicit none

    type(type_gl) :: gl


    double precision:: T, D
    integer :: oir
    double precision, dimension(30), intent(out):: chempot

    double precision::  Rmix, wmmix, p, a, v, dgds, d_w

    integer:: i

    !oir = 0                    !The overall chemical potential will be calculated (1 = ideal, 2 = residual)
    chempot = 0.d0                !Initialize chempot
    call R_mix_calc(gl,Rmix)
    call wm_mix_calc(gl,wmmix)

    call dna_dni(gl,T, D, CHEMPOT, OIR)

    if(gl%seawater) gl%gepress = gl%sea%seap
    p = gl%gepress
    Do i = 1, gl%NCOMP

        if(gl%gecarrier .and. gl%seawater) then
            call chempot_num_reac(gl, t, d, Chempot, 1)
        elseif(gl%gecarrier .and. gl%el_present) then
            !call partial_molar_gibbs(gl, t, d, p, chempot)
            call chempot_num_reac(gl, t, d, Chempot, 1)
        else
            Chempot(i)= Chempot(i) * Rmix * T
        end if
    end do

    end subroutine Chempot_CALC



    !-----------------------------------------
    !Function for speed of sound [m/s]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION WS_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "WS_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst       , i                          ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD,DD_ARDD,TT_ARTT,DT_ARDT,TT_AITT,wmmix, Rmix,RADICAL, ardel, a2del, ataudel,  a2tau, wm_water, rad1, &
        & rad2, ws_sea_calc, ws_sea_mix, d_mixsave, wm_w
    logical :: seasave


    GETDERR=(/0,1,1,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    TT_AITT = FNIDER(6)

    WS_CALC = 0.D0

    if(gl%seawater) gl%gepress = gl%sea%seap

    if( (gl%gecarrier) .and. (nrsubst .eq. 1) ) then

        ardel = a01_gibbs(gl, t, d, gl%gepress)
        a2del = a02_gibbs(gl, t, d, gl%gepress)
        ataudel = a11_gibbs(gl,t, d, gl%gepress)
        a2tau = a20_gibbs(gl, t, d, gl%gepress)

        if(gl%seawater) then
            wmmix = gl%sea%wm_sea
        elseif(gl%el_present) then
            wmmix = wm_brine(gl)
        end if

        RADICAL = (1.D0  +2.d0*ardel + a2del  - (1.D0 +  ardel - ataudel )**2.d0/(a2tau)) * Rmix/wmmix * T


    elseif((gl%gecarrier) .and. (nrsubst .eq. 0)) then

        gl%gemix =.true.
        wm_w = gl%wm(1)
        if(gl%seawater) then
            gl%wm(1) = gl%sea%wm_sea
        elseif(gl%el_present) then
            gl%wm(1) = wm_brine(gl)
        end if
        call wm_mix_calc(gl, wmmix)
        gl%wm(1) = wm_w

        call solvent_cor(gl, t, gl%gepress, getderi, getderr, cor_deri, cor_derr)

        ardel = a01_gibbs(gl, t, d, gl%gepress) - cor_derr(2)
        a2del = a02_gibbs(gl, t, d, gl%gepress) - cor_derr(3)
        ataudel = a11_gibbs(gl,t, d, gl%gepress) - cor_derr(6)
        a2tau = a20_gibbs(gl, t, d, gl%gepress) - cor_derr(5) - cor_deri(6)

        RADICAL = (1.D0 + 2.D0 *( D_ARD + ardel) + DD_ARDD +a2del - (1.D0 + D_ARD + ardel - DT_ARDT - ataudel )**2.d0/(TT_AITT + TT_ARTT + a2tau )) * Rmix/wmmix * T
        gl%gemix = .false.

    else
        RADICAL = (1.D0 + 2.D0 * D_ARD + DD_ARDD - (1.D0 + D_ARD - DT_ARDT)**2/(TT_AITT + TT_ARTT)) * Rmix/wmmix * T
    end if

    IF (RADICAL < 0.D0) THEN
        WS_CALC = -9971
    ELSE
        WS_CALC = RADICAL**0.5d0

    END IF

    d=d_mixsave

    END FUNCTION WS_CALC




    !-----------------------------------------
    !Function for joule - thomson coefficient
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION JTCO_CALC(gl,T,D,nrsubst)



    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TT_AITT,D_ARD,DD_ARDD,TT_ARTT,DT_ARDT,wmmix, Rmix

    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERR=(/0,1,1,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    TT_AITT = FNIDER(6)

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)


    JTCO_CALC = 0.D0
    JTCO_CALC = - (D_ARD + DD_ARDD + DT_ARDT)/((1.D0 + D_ARD - DT_ARDT)**2 - (TT_AITT + TT_ARTT) * &
        (1.D0 + 2.D0 * D_ARD + DD_ARDD)) /( Rmix* D) * gl%factortrans !/wmmix (müsste hinter Rmix)

    END FUNCTION JTCO_CALC



    !----------------------------------------------
    !Function for Second thermal virial coefficient [m^3/mol] - geändert von [cc/mol] - J.G. 11.2011
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION B_CALC(gl,T, nrsubst)


    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD, rhoredmix_orig, tredmix_orig

    GETDERR=(/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    If ((gl%mix_Type == 1) .or. (gl%mix_Type == 12) .or. (gl%mix_Type == 13)) then
        D = 0.D0!1.D-12     rho = 0 works for Helmholtz EOS only, at the moment
    else
        D = 1.D-12
    end if

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of ideal part
    else
        gl%rhoredmix = gl%rhored(nrsubst)
        gl%tredmix = gl%tred(nrsubst)
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    D_ARD = FNRDER(2)

    B_CALC = 0.D0
    If ((gl%mix_Type == 1) .or. (gl%mix_Type == 12) .or. (gl%mix_Type == 13)) then
        B_CALC = D_ARD/ gl%rhoredmix !*1.d+6    !Andreas, September 2015
    else
        B_CALC = D_ARD/ d !*1.d+6
    end if

    gl%rhoredmix =rhoredmix_orig
    gl%tredmix = tredmix_orig

    END FUNCTION B_CALC

    !----------------------------------------------
    !Function for the first derivative of the second thermal virial coefficient WRT T [m^3/mol/T]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION DBDT_CALC(gl,T, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,DT_ARDT, rhoredmix_orig, tredmix_orig

    GETDERR=(/0,0,0,0,0,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    If (gl%mix_Type == 1) then
        D = 0.D0!1.D-12     rho = 0 works for Helmholtz EOS only, at the moment
    else
        D = 1.D-12
    end if

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part of the mixture
    else
        gl%rhoredmix = gl%rhored(nrsubst)
        gl%tredmix = gl%tred(nrsubst)
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    DT_ARDT = FNRDER(6)

    DBDT_CALC = 0.D0

    If (gl%mix_Type == 1) then
        DBDT_CALC = (-1.d0)* DT_ARDT/ gl%rhoredmix * gl%tredmix / T**2
    else
        DBDT_CALC = (-1.d0)* DT_ARDT/ d / T
    end if

    gl%rhoredmix =rhoredmix_orig
    gl%tredmix = tredmix_orig

    END FUNCTION DBDT_CALC

    !---------------------------------------------------------
    !Function for Third thermal virial coefficient [m^3/mol]^2
    !---------------------------------------------------------
    DOUBLE PRECISION module FUNCTION C_CALC(gl,T, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,DD_ARDD, rhoredmix_orig, tredmix_orig

    GETDERR=(/0,0,1,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    If (gl%mix_Type == 1) then
        D = 0.D0!1.D-10     rho = 0 works for Helmholtz EOS only, at the moment
    else
        D = 1.D-10
    end if

    !Andreas, September 2015
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix

    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        gl%rhoredmix = gl%rhored(nrsubst)
        gl%tredmix = gl%tred(nrsubst)
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    DD_ARDD = FNRDER(3)

    C_CALC = 0.D0

    If (gl%mix_Type == 1) then
        C_CALC = DD_ARDD/gl%rhoredmix**2!*(1.D+6)**2.D0
    else
        C_CALC = DD_ARDD/D**2!*(1.D+6)**2.D0
    end if

    gl%rhoredmix =rhoredmix_orig
    gl%tredmix = tredmix_orig

    END FUNCTION C_CALC


    !----------------------------------------------
    !Function for the first derivative of the third thermal virial coefficient WRT T [m^6/mol^2/T]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION DCDT_CALC(gl,T, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,DDT_ARDT, rhoredmix_orig, tredmix_orig

    GETDERR=(/0,0,0,0,0,0,0,0,0,1,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (gl%mix_Type == 1) then
        D = 0.D0!1.D-10     rho = 0 works for Helmholtz EOS only, at the moment
    else
        D = 1.D-10
    end if

    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix


    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of ideal part
    else
        gl%rhoredmix = gl%rhored(nrsubst)
        gl%tredmix = gl%tred(nrsubst)
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    DDT_ARDT = FNRDER(10)

    DCDT_CALC = 0.D0

    if (gl%Mix_type == 1) then
        DCDT_CALC = (-1.d0)* DDT_ARDT/ gl%rhoredmix**2 * gl%tredmix / T**2
    else
        DCDT_CALC = (-1.d0)* DDT_ARDT/d**2/T
    end if

    gl%rhoredmix =rhoredmix_orig
    gl%tredmix = tredmix_orig

    END FUNCTION DCDT_CALC



    !---------------------------------------------------------
    !Function for fouth thermal virial coefficient [m^3/mol]^3
    !---------------------------------------------------------
    DOUBLE PRECISION module FUNCTION D_CALC(gl,T, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,DDD_ARDDD

    GETDERR=(/0,0,0,0,0,0,0,1,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    If (gl%mix_Type == 1) then
        D = 0.D0!1.D-6     rho = 0 works for Helmholtz EOS only, at the moment
    else
        D = 1.D-6
    end if

    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    DDD_ARDDD = FNRDER(8)

    D_CALC = 0.D0

    if (gl%mix_type == 1) then
        D_CALC = DDD_ARDDD/gl%rhoredmix**3/2.d0
    else
        D_CALC = DDD_ARDDD/D**3/2.d0
    end if

    END FUNCTION D_CALC


    !---------------------------------------------------------
    !Function for the volume expansivity alpha [1/K]
    !---------------------------------------------------------
    DOUBLE PRECISION module FUNCTION VOLEXP_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: T, D, DPDD, DPDT
    INTEGER:: nrsubst


    DPDD=DPDD_CALC(gl,T,D,nrsubst)
    DPDT=DPDT_CALC(gl,T,D,nrsubst)

    VOLEXP_CALC = 1.d0/D/DPDD*DPDT

    END FUNCTION VOLEXP_CALC


    !-----------------------------------------------------------
    !Function for the isothermal compressibility kappa T [1/MPa]
    !-----------------------------------------------------------
    DOUBLE PRECISION module FUNCTION COMPT_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: T, D, DPDD
    INTEGER:: nrsubst


    DPDD=DPDD_CALC(gl,T,D,nrsubst)

    COMPT_CALC = 1.d0/D/DPDD

    END FUNCTION COMPT_CALC

    !-----------------------------------------------------------
    !Function for the isentropic compressibility kappa s [1/MPa]
    !-----------------------------------------------------------
    DOUBLE PRECISION module FUNCTION COMPS_CALC(gl,T,D, nrsubst)



    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD,DD_ARDD,TT_ARTT,DT_ARDT,TT_AITT,wmmix, Rmix,RADICAL

    GETDERR=(/0,1,1,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)



    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    TT_AITT = FNIDER(6)

    RADICAL = (1.D0 + 2.D0 * D_ARD + DD_ARDD - (1.D0 + D_ARD - DT_ARDT)**2/(TT_AITT + TT_ARTT)) * Rmix * T

    COMPS_CALC = 1.d0/D/RADICAL * 1.d6 / gl%factorpress

    END FUNCTION COMPS_CALC


    !-----------------------------------------------------------------------
    !Function for the isothermal throtteling coefficient delta T [m³/mol]
    !-----------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION THROT_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD,DD_ARDD,DT_ARDT,wmmix, Rmix

    GETDERR=(/0,1,1,0,0,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    DT_ARDT = FNRDER(6)

    THROT_CALC = (1.d0 - (1.d0 + D_ARD - DT_ARDT)/(1.d0 + 2.d0 * D_ARD + DD_ARDD)) / D

    END FUNCTION THROT_CALC


    !-----------------------------------------------------------
    !Function for the isentropic expansion coefficient ks [-]
    !-----------------------------------------------------------
    DOUBLE PRECISION module FUNCTION EXPANS_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,D_ARD,DD_ARDD,TT_ARTT,DT_ARDT,TT_AITT,wmmix, Rmix,RADICAL,P

    GETDERR=(/0,1,1,0,1,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DD_ARDD = FNRDER(3)
    TT_ARTT = FNRDER(5)
    DT_ARDT = FNRDER(6)

    TT_AITT = FNIDER(6)

    RADICAL = (1.D0 + 2.D0 * D_ARD + DD_ARDD - (1.D0 + D_ARD - DT_ARDT)**2/(TT_AITT + TT_ARTT)) * Rmix * T
    P = (1.D0 + D_ARD) * D * Rmix * T

    EXPANS_CALC = D*RADICAL/P

    END FUNCTION EXPANS_CALC


    !-----------------------------------------------------------
    !Function for the isothermal expansion coefficient kT [-]
    !-----------------------------------------------------------
    DOUBLE PRECISION module FUNCTION EXPANT_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: T, D, P,DPDD
    INTEGER:: nrsubst

    P = P_CALC(gl,T,D,nrsubst)
    DPDD = DPDD_CALC(gl,T,D,nrsubst)

    EXPANT_CALC = D/P*DPDD

    END FUNCTION EXPANT_CALC


    !----------------------------------------------------------------------------------------------------------------------------------------
    !  Calculate the derivative of the internal energy with respect to the volume at constant temperature [(dU/dV) @ T=const] DUDV [J/mol/m³]
    !----------------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DUDV_CALC(gl,T,D, nrsubst)


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: T, D, P,DPDT
    INTEGER:: nrsubst


    P = P_CALC(gl,T,D,nrsubst)
    DPDT = DPDT_CALC(gl,T,D,nrsubst)

    DUDV_CALC = -P + T * DPDT

    END FUNCTION DUDV_CALC


    !----------------------------------------------
    !Function for Helmholtz energy [J/mol]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION A_CALC(gl,T, D, nrsubst)





    implicit none

    type(type_gl) :: gl


    !ASSEMBLY_INTERFACE(NAME = "A_CALC")

    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst  , i                               ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER, cor_derr               !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER, cor_deri               !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AI,AR,wmmix, Rmix, a_ge, Ms, Msea, chem, d_mixsave, a_mix, a_gi, a_diff
    logical :: seasave

    GETDERR=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/1,0,0,0,0,0,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    d_mixsave = d

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AI = FNIDER(1)
    AR = FNRDER(1)

    A_CALC = 0.D0

    a_mix = AI + AR

    if(gl%seawater) gl%gepress = gl%sea%seap


    if((gl%gecarrier ) .and. (nrsubst .eq. 1))then


        a_calc =  a_gibbs(gl, t, d, gl%gepress) * Rmix * T

    elseif( (gl%gecarrier) .and. (nrsubst .eq. 0) ) then
        gl%gemix = .true.

        call solvent_cor(gl, t,  gl%gepress, getderi, getderr, cor_deri, cor_derr)

        a_ge = a_gibbs(gl, t, d, gl%gepress) - cor_deri(1) - cor_derr(1)

        a_calc = ( a_ge + AI + AR) * Rmix * T

        gl%gemix =.false.
    else

        a_calc = (ai + ar) * Rmix * T
    end if
    d=d_mixsave

    END FUNCTION A_CALC

    !----------------------------------------------
    !Function for the ideal Helmholtz energy [J/mol]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION AI_CALC(gl,T, D, nrsubst)




    implicit none

    type(type_gl) :: gl


    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AI,Rmix,wmmix

    GETDERI=(/1,0,0,0,0,0,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AI = FNIDER(1)

    AI_CALC = 0.D0
    AI_CALC = AI * Rmix*T

    END FUNCTION AI_CALC


    !----------------------------------------------
    !Function for Second acoustic virial coefficient
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION BETA_CALC(gl,T, nrsubst)

    !   USE module_general_eos_parameters




    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,KSI,TN,DN, CP0,D_ARD,DT_ARDT,DTT_ARDTT,wmmix, Rmix

    GETDERR=(/0,1,0,0,0,1,1,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,0,1,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    D = 1.D-12

    TN = gl%tredmix/T
    DN = D/gl%rhoredmix

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERI,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    DT_ARDT = FNRDER(6)
    DTT_ARDTT = FNRDER(7)

    CP0 = (1.D0 - FNIDER(6))*Rmix
    KSI = CP0/(CP0 - Rmix)

    BETA_CALC = 0.D0
    BETA_CALC = (2.D0 * D_ARD / DN - 2.d0 * (KSI - 1.d0) / KSI / DN * DT_ARDT + (KSI - 1.D0)**2 / KSI  / DN * DTT_ARDTT) / gl%rhoredmix

    ! BETA_CALC = (2.D0 * D_ARD - 2.d0 * (KSI - 1.d0) / KSI / DN * DT_ARDT + (KSI - 1.D0)**2 / KSI * TN**2 * DTT_ARDTT) / rhoredmix

    END FUNCTION BETA_CALC


    !-----------------------------------------
    !Function for reference enthalpy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION H_REF(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: GETDERI
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,T_AIT,T_ART,D_ARD, rhored_orig, tred_orig, Rmix!,wmmix

    GETDERR=(/0,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    GETDERI=(/0,0,0,0,1,0,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    ! the reference state of each component needs to be calculated at the reduced pure fluid parameters
    rhored_orig = gl%rhoredmix
    gl%rhoredmix = gl%rhored(nrsubst)
    tred_orig = gl%tredmix
    gl%tredmix = gl%tred(nrsubst)

    CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part

    !wmmix = wm(nrsubst)
    Rmix = gl%REQ(nrsubst)

    T_AIT = FNIDER(5)

    D_ARD = FNRDER(2)
    T_ART = FNRDER(4)

    H_REF = 0.D0
    H_REF = (1 + T_AIT + T_ART + D_ARD) * Rmix * T!/wmmix

    gl%rhoredmix = rhored_orig
    gl%tredmix = tred_orig

    END FUNCTION H_REF


    !-----------------------------------------
    !Function for calculating the reference entropy in [J/(mol-K)]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION S_REF(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer,dimension(nderivsi) :: getderi
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER              !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AI,T_AIT,AR,T_ART, rhored_orig, tred_orig, Rmix!,wmmix

    GETDERR=(/1,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)
    getderi=(/1,0,0,0,1,0,0,0,0,0/)     !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    ! the reference state of each component needs to be calculated at the reduced pure fluid parameters
    rhored_orig = gl%rhoredmix
    gl%rhoredmix = gl%rhored(nrsubst)
    tred_orig = gl%tredmix
    gl%tredmix = gl%tred(nrsubst)

    CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    CALL FNIDERIVS(gl,T,D,GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
    !wmmix = wm(nrsubst)
    Rmix = gl%REQ(nrsubst)

    AI = FNIDER(1)
    T_AIT = FNIDER(5)

    AR = FNRDER(1)
    T_ART = FNRDER(4)

    S_REF = 0.D0
    S_REF = (T_AIT + T_ART - AI - AR) * Rmix!/wmmix

    gl%rhoredmix = rhored_orig
    gl%tredmix = tred_orig

    end FUNCTION S_REF

    !-----------------------------------------
    !Function for residual internal energy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION UR_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,T_ART,wmmix, Rmix

    GETDERR=(/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    T_ART = FNRDER(4)

    UR_CALC = 0.D0
    UR_CALC = T_ART * Rmix * T

    END FUNCTION UR_CALC


    !-----------------------------------------
    !Function for calculating residual entropy [J/(mol-K)]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION SR_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    DOUBLE PRECISION :: T,D,AR,T_ART,wmmix, Rmix

    GETDERR=(/1,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AR = FNRDER(1)
    T_ART = FNRDER(4)

    SR_CALC = 0.D0
    SR_CALC = (T_ART - AR) * Rmix

    end FUNCTION SR_CALC


    !--------------------------------------------------------
    !Function for residual isochoric heat capacity [J/(mol-K)]
    !---------------------------------------------------------
    DOUBLE PRECISION module FUNCTION CVR_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,TT_ARTT,wmmix, Rmix

    GETDERR=(/1,0,0,0,1,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    TT_ARTT = FNRDER(5)

    CVR_CALC = 0.D0
    CVR_CALC = - TT_ARTT * Rmix

    END FUNCTION CVR_CALC


    !-----------------------------------------
    !Function for residual enthalpy [J/mol]
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION HR_CALC(gl,T,D, nrsubst)






    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,T_ART,D_ARD,wmmix, Rmix

    GETDERR=(/0,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)
    T_ART = FNRDER(4)

    HR_CALC = 0.D0
    HR_CALC = (T_ART + D_ARD) * Rmix * T

    END FUNCTION HR_CALC


    !----------------------------------------------
    !Function for residual Helmholtz energy [J/mol]
    !----------------------------------------------
    DOUBLE PRECISION module FUNCTION AR_CALC(gl,T, D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                !Vector gives the values of the calculated derivatives

    DOUBLE PRECISION :: T,D,AR,wmmix, Rmix

    GETDERR=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        wmmix = gl%wm(nrsubst)
        Rmix = gl%REQ(nrsubst)
    end if

    AR = FNRDER(1)

    AR_CALC = 0.D0
    AR_CALC = AR * Rmix*T

    END FUNCTION AR_CALC

    !----------------------------------------------------------------------------------
    !Function for calculating residual pressure [MPa]
    !----------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION PR_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD, Rmix!,wmmix

    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)

    PR_CALC = 0.D0
    PR_CALC = D_ARD * D * Rmix * T * 1.D-6 * gl%factorpress

    END FUNCTION PR_CALC


    !----------------------------------------------------------------------------------
    !Function for calculating the residual alpha : A00 [-]
    !-----------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A00_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,AR, Rmix

    GETDERR = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    AR = FNRDER(1)

    A00_CALC = 0.D0
    A00_CALC = AR

    END FUNCTION A00_CALC

    !-----------------------------------------------------------------------------------------------------
    !Function for calculating the (1st derivative of the residual alpha with respect to tau) * tau : A10 [-]
    !-------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A10_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    !integer,DIMENSION(15) :: GETDERR = (/1,1,1,1,1,1,1,1,1,1/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated


    DOUBLE PRECISION :: T,D,T_ART, Rmix

    GETDERR = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        !CALL FNRDERIVS_NUM(T,D,GETDERR,FNRDERNUM,nrsubst)   !subroutine calculates derivatives of residual part

        !do i=1,10
        !diff(i) = FNRDER(i) - FNRDERNUM(i)
        !!-0.268321837611150 Moni
        !end do


        Rmix = gl%REQ(nrsubst)
    end if

    T_ART = FNRDER(4)

    A10_CALC = 0.D0
    A10_CALC = T_ART

    END FUNCTION A10_CALC

    !---------------------------------------------------------------------------------------------------------
    !Function for calculating the (1st derivative of the residual alpha with respect to delta) * delta : A01 [-]
    !-----------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A01_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER          !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,D_ARD, Rmix

    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    D_ARD = FNRDER(2)

    A01_CALC = 0.D0
    A01_CALC = D_ARD

    END FUNCTION A01_CALC

    !-----------------------------------------------------------------------------------------------------------------------------
    !Function for calculating the (1st mixed derivative of the residual alpha with respect to tau and delta) * tau * delta : A11 [-]
    !-------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A11_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,DT_ARDT, Rmix

    GETDERR = (/0,0,0,0,0,1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    DT_ARDT = FNRDER(6)

    A11_CALC = 0.D0
    A11_CALC = DT_ARDT

    END FUNCTION A11_CALC

    !---------------------------------------------------------------------------------------------------------
    !Function for calculating the (2nd derivative of the residual alpha with respect to tau^2) * tau^2 : A20 [-]
    !-----------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A20_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,TT_ARTT, Rmix

    GETDERR = (/0,0,0,0,1,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    TT_ARTT = FNRDER(5)

    A20_CALC = 0.D0
    A20_CALC = TT_ARTT

    END FUNCTION A20_CALC

    !-------------------------------------------------------------------------------------------------------------
    !Function for calculating the (2nd derivative of the residual alpha with respect to delta^2) * delta^2 : A02 [-]
    !---------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A02_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,DD_ARDD, Rmix

    GETDERR = (/0,0,1,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    DD_ARDD = FNRDER(3)

    A02_CALC = 0.D0
    A02_CALC = DD_ARDD

    END FUNCTION A02_CALC

    !--------------------------------------------------------------------------------------------------------------------------------
    !Function for calculating the (2nd mixed derivative of the residual alpha with respect to tau and delta^2) * tau * delta^2: A12 [-]
    !----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A12_CALC(gl,T,D, nrsubst)





    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,DDT_ARDDT, Rmix

    GETDERR = (/0,0,0,0,0,0,0,0,0,1,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    DDT_ARDDT = FNRDER(10)

    A12_CALC = 0.D0
    A12_CALC = DDT_ARDDT

    END FUNCTION A12_CALC


    !-----------------------------------------------------------------------------------------------------------
    !Function for calculating the (3rd derivative of the residual alpha with respect to delta) * delta^3 : A03 [-]
    !-------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A03_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,DDD_ARDDD, Rmix

    GETDERR = (/0,0,0,0,0,0,0,1,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    DDD_ARDDD = FNRDER(8)

    A03_CALC = 0.D0
    A03_CALC = DDD_ARDDD

    END FUNCTION A03_CALC


    !--------------------------------------------------------------------------------------------------------
    !Function for calculating the (3rd derivative of the residual alpha with respect to tau) * tau^3: A30 [-]
    !--------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A30_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR
    integer:: nrsubst                                                            ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,TTT_ARTTT, Rmix

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER                                     !Vector gives the values of the calculated derivatives

    GETDERR = (/0,0,0,0,0,0,0,0,1,0,0,0,0,0,0/)                   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)     ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    TTT_ARTTT = FNRDER(9)

    A30_CALC = 0.D0
    A30_CALC = TTT_ARTTT

    END FUNCTION A30_CALC

    DOUBLE PRECISION module FUNCTION A40_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    double precision :: num_del = 1D-6


    A40_CALC =  (   A30_CALC(gl,T+num_del,D, nrsubst)   -  A30_CALC(gl,T-num_del,D, nrsubst) ) /(2d0*num_del)



    end FUNCTION A40_CALC

    !----------------------------------------------------------------------------------------------------------------------------------
    !Function for calculating the (2nd mixed derivative of the residual alpha with respect to tau^2 and delta) * tau^2 * delta: A21 [-]
    !----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION A21_CALC(gl,T,D, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer,DIMENSION(nderivs) :: GETDERR

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: nrsubst,I                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T,D,DTT_ARDTT, Rmix

    GETDERR = (/0,0,0,0,0,0,1,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)


    if (nrsubst == 0) then
        call R_mix_calc(gl,Rmix)                   ! calculates the "mixture molar gas constant"
        CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if

    DTT_ARDTT = FNRDER(7)

    A21_CALC = 0.D0
    A21_CALC = DTT_ARDTT

    END FUNCTION A21_CALC



    !**************************************************************************************
    ! functions needed internally in case a property is to be calculated from ph or ps
    ! Theresa Wiens, Dez. 2011
    !**************************************************************************************

    DOUBLE PRECISION module FUNCTION S_CALC_TP(gl,T,P)



    implicit none

    type(type_gl) :: gl


    double precision :: T,P
    double precision :: density
    !integer :: errormsg
    integer :: nrsubst

    !errormsg = 0
    nrsubst = 1

    ! iterate density from given T and p
    density = rhomix_calc(gl,T, P, 0.d0, 0, nrsubst)

    if (density > 0) then ! successfull
        S_CALC_TP = S_CALC(gl,T,density,nrsubst)

        gl%rho_it = density    ! save iterated density for later usage
        gl%p_it = p            !control parameter
    else                                        ! density iteration failed
        S_CALC_TP = -44444.d0
    end if


    END FUNCTION S_CALC_TP

    DOUBLE PRECISION module FUNCTION H_CALC_TP(gl,T,P)



    implicit none

    type(type_gl) :: gl


    double precision :: T,P
    double precision :: density
    !integer :: errormsg
    integer :: nrsubst

    !errormsg = 0
    nrsubst = 1


    ! iterate density from given T and p
    density = rhomix_calc(gl,T, P, 0.d0, 0, nrsubst)

    if (density > 0) then ! successfull
        H_CALC_TP = H_CALC(gl,T,density,nrsubst)

        gl%rho_it = density    ! save iterated density for later usage
        gl%p_it = p            !control parameter

    else                                        ! density iteration failed
        H_CALC_TP = -44444.d0
    end if


    END FUNCTION H_CALC_TP



    DOUBLE PRECISION module FUNCTION RIEM_CALC (gl,t,d,nrsubst)
    !F_STDCALL DOUBLE PRECISION module FUNCTION RIEM_CALC (gl,input, prop1, prop2, fluids, moles, path)
    !!DEC$ ATTRIBUTES DLLEXPORT :: RIEM_CALC








    implicit none

    type(type_gl) :: gl


    !define input variables
    character (12) :: input
    double precision :: t, d
    character (255) :: fluids
    character (255) :: moles
    character (255) :: path
    character (255) :: EOS_indicator

    double precision :: A01, A02, A03, A20, A30, A11, A12, A21
    double precision :: xna,sum, term
    double precision :: tau, del
    integer, dimension(nderivs)::GETDERIVS
    integer, dimension(nderivsi)::GETDERIVSID
    DOUBLE PRECISION, dimension(nderivs)::FNRDER
    DOUBLE PRECISION, dimension(nderivsi)::FNIDER
    integer::nrsubst
    double precision::A20_2,A11_2,A21_2,A12_2, tau2, del2, tau3, del3!exponential terms

    integer:: i

    GETDERIVSID=(/1,1,1,1,1,1,1,1,1,1/)
    GETDERIVS=(/1,1,1,1,1,1,1,1,1,1,0,0,0,0,0/)

    FNRDER=0.d0
    FNIDER=0.d0

    if (nrsubst == 0) then
        CALL MIXDERIVSFNR(gl,T,D,GETDERIVS,FNRDER)   !subroutine calculates derivatives of residual part
        CALL MIXDERIVSFNI(gl,T,D,GETDERIVSID,FNIDER)   !subroutine calculates derivatives of ideal part
    else
        CALL FNRDERIVS(gl,T,D,GETDERIVS,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        CALL FNIDERIVS(gl,T,D,GETDERIVSID,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
    end if


    call reduced_parameters_calc(gl,300.d0)
    tau=gl%tredmix/t
    del=d/gl%rhoredmix

    tau2=tau*tau
    tau3=tau2*tau
    del2=del*del
    del3=del2*del

    !ideal and residual helmholtz energies derivatives
    A01=(FNIDER(2)+ FNRDER(2))/del
    A02=(FNIDER(3)+ FNRDER(3))/del2
    A03=(FNIDER(8)+FNRDER(8))/del3
    A11=(FNIDER(4)+ FNRDER(6))/tau/del     !FNIDER(4):=1ST MIXED DERIVATIVE WITH RESPECT TO D AND T
    A12=(FNIDER(10)+ FNRDER(10))/tau/del2
    A21=(FNIDER(7)+FNRDER(7))/tau2/del
    A20=(FNIDER(6)+ FNRDER(5))/tau2
    A30=(FNIDER(9)+ FNRDER(9))/tau3

    !terms in riemc
    A20_2=A20*A20
    A11_2=A11*A11
    A21_2=A21*A21
    A12_2=A12*A12
    term=(2.d0*A01+del*A02)*(2.d0*A01+del*A02)

    sum=-2.d0*       A20_2*A01&
        +4.d0*del*   A30*A01*A11 +4.d0*del   *A20*A11_2&
        -4.d0*del*   A20*A01*A21 -2.d0*del2*A01*A21_2&
        -4.d0*del*   A20_2*A02  +2.d0*del2*A30*A11*A02&
        -5.d0*del2*A20*A21*A02 -     del3*A21_2*A02&
        +2.d0*del2*A30*A01*A12 +4.d0*del2*A20*A11*A12&
        +    del3*A30*A02*A12 +     del3*A20*A12_2&
        -   del2*A20_2*A03  -     del3*A20*A21*A03

    xna=6.0221367d+23!Avogadro-Constant

    riem_calc = 0.d0

    riem_calc=1.d0/(2.d0*xna*del2*gl%rhoredmix*A20_2*term)*sum
    !unit: m³

    return
    end function



    !-----------------------------------------
    !Function for the Boyle temperature
    !criterion: B = 0
    !-----------------------------------------
    DOUBLE PRECISION module FUNCTION TBoyle_calc(gl,fluidnr)








    implicit none

    type(type_gl) :: gl


    INTEGER :: Max_Iterations, Iterations
    DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
    DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
    DOUBLE PRECISION :: Delta_allowed
    INTEGER :: errTBoyle,fluidnr,nrsubst
    type(type_additional_parameters) :: ResTBoyle_param
    DOUBLE PRECISION :: TBoyle

    nrsubst= 0
    TBoyle=0.d0
    !ResTBoyle_param = 0.d0
    Tstart_min = 0.d0
    Tstart_max = 0.d0
    Tmin_allowed = 0.d0
    Tmax_allowed = 0.d0
    Delta_allowed = 0.d0
    errTBoyle=0

    !Parameters for iteration:
    nrsubst = fluidnr
    Delta_allowed = 1.d-8
    Max_Iterations = 50
    ResTBoyle_param%a_p(1) = 0.d0
    ResTBoyle_param%a_p(2) = fluidnr
    !ResPsat_param(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration


    if (gl%ncomp > 1) then !property only available for pures

        TBoyle_calc = -9902

    else


        Tstart_min = gl%tc(nrsubst)
        Tmin_allowed = gl%tc(nrsubst)*0.5d0

        Tstart_max = 1.d4
        Tmax_allowed = 1.d6

        CALL Regula_Falsi(gl,Res_TBoyle,TBoyle,Tstart_min,Tstart_max,Delta_allowed,&
            Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errTBoyle,ResTBoyle_param)

        if (errTBoyle == 0) then
            TBoyle_calc=TBoyle
        else
            TBoyle_calc = -5662.d0
        end if


    end if

    End Function TBoyle_calc


    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Res_TBoyle(gl,T_akt,parameters)
    !--------------------------------------------------------------------------
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: T_akt, BTBoyle, TBoyle_akt
    type(type_additional_parameters) :: parameters
    Integer :: i


    BTBoyle=parameters%a_p(1)
    i=int(parameters%a_p(2))

    TBoyle_akt = B_CALC(gl,T_akt, i)

    Res_TBoyle=BTBoyle-TBoyle_akt

    RETURN
    END function




    !------------------------------------------
    !Function for the Joule Inversion temperature
    !criterion: dB/dT = 0
    !------------------------------------------
    DOUBLE PRECISION module FUNCTION TJT_calc(gl,fluidnr)

    implicit none

    type(type_gl) :: gl


    INTEGER :: Max_Iterations, Iterations
    DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
    DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
    DOUBLE PRECISION :: Delta_allowed
    INTEGER :: errTJT,fluidnr,nrsubst
    type(type_additional_parameters) :: ResTJT_param
    DOUBLE PRECISION :: TJT, TBoyle

    nrsubst= 0
    TJT=0.d0
    TBoyle=0.d0
    !ResTJT_param = 0.d0
    Tstart_min = 0.d0
    Tstart_max = 0.d0
    Tmin_allowed = 0.d0
    Tmax_allowed = 0.d0
    Delta_allowed = 0.d0
    errTJT=0

    !Parameters for iteration:
    nrsubst = fluidnr
    Delta_allowed = 1.d-8
    Max_Iterations = 500
    ResTJT_param%a_p(1) = 0.d0
    ResTJT_param%a_p(2) = fluidnr
    !ResPsat_param(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration



    if (gl%ncomp > 1) then !property only available for pures

        TJT_calc= -9902

    else

        TBoyle = TBoyle_calc(gl,fluidnr)

        Tstart_min = 1.1d0 * TBoyle
        Tmin_allowed = TBoyle

        !Tstart_min = gl%tc(nrsubst)*1.2d0
        !Tmin_allowed = gl%tc(nrsubst)*0.5d0

        !Tstart_max = 5.d4
        Tstart_max = 1.d4
        Tmax_allowed = 5.d5

        CALL Regula_Falsi(gl,Res_TJT,TJT,Tstart_min,Tstart_max,Delta_allowed,&
            Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errTJT,ResTJT_param)

        if (errTJT == 0) then
            TJT_calc=TJT
        else
            TJT_calc = -5662.d0
        end if

    end if

    End Function TJT_calc


    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Res_TJT(gl,T_akt,parameters)
    !--------------------------------------------------------------------------
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: T_akt, BTJT, TJT_akt
    type(type_additional_parameters) :: parameters
    Integer :: i

    BTJT=parameters%a_p(1)
    i=int(parameters%a_p(2))


    TJT_akt = DBDT_CALC(gl,T_akt, i)

    Res_TJT=BTJT-TJT_akt
    !sinnvoll für seewasser?
    RETURN
    END function



    !----------------------------------------------------
    !Function for the Joule-Thomson inversion temperature
    !criterion: dB/dT = B/T
    !----------------------------------------------------
    DOUBLE PRECISION module FUNCTION TJTINV_calc(gl,fluidnr)

    implicit none

    type(type_gl) :: gl


    INTEGER :: Max_Iterations, Iterations
    DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
    DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
    DOUBLE PRECISION :: Delta_allowed
    INTEGER :: errTJTINV,fluidnr,nrsubst
    type(type_additional_parameters):: ResTJTINV_param
    DOUBLE PRECISION :: TJTINV, TJT, TBoyle

    nrsubst= 0
    TJTINV=0.d0
    TJT = 0.d0
    !ResTJTINV_param = 0.d0
    Tstart_min = 0.d0
    Tstart_max = 0.d0
    Tmin_allowed = 0.d0
    Tmax_allowed = 0.d0
    Delta_allowed = 0.d0
    errTJTINV=0

    !Parameters for iteration:
    nrsubst = fluidnr
    Delta_allowed = 1.d-8
    Max_Iterations = 50
    ResTJTINV_param%a_p(1) = 0.d0
    ResTJTINV_param%a_p(2) = fluidnr
    !ResPsat_param(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration


    if (gl%ncomp > 1.d0) then !property only available for pures

        TJTINV_calc = -9902

    else
        TBoyle = TBoyle_calc(gl,fluidnr)

        Tstart_min = 1.1d0 * TBoyle
        Tmin_allowed = TBoyle

        TJT = TJT_calc(gl,fluidnr)
        Tstart_max = 0.9d0 * TJT
        if (TJT .gt. 1.d-10) then
            Tmax_allowed = TJT
        else
            Tmax_allowed = 1.d4
        end if

        CALL Regula_Falsi(gl,Res_TJTINV,TJTINV,Tstart_min,Tstart_max,Delta_allowed,&
            Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errTJTINV,ResTJTINV_param)

        if (errTJTINV == 0) then
            TJTINV_calc=TJTINV
        else
            TJTINV_calc = -5662.d0
        end if


    end if

    End Function TJTINV_calc


    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Res_TJTINV(gl,T_akt,parameters)
    !--------------------------------------------------------------------------
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: T_akt, BTJTINV, TJTINV_akt
    type(type_additional_parameters) :: parameters
    Integer :: i

    BTJTINV=parameters%a_p(1)
    i=int(parameters%a_p(2))

    TJTINV_akt =  DBDT_CALC(gl,T_akt, i) - B_CALC(gl,T_akt, i) / T_akt

    Res_TJTINV=BTJTINV-TJTINV_akt
    !sinnvoll für seewasser?
    RETURN
    END function



    !-----------------------------------------------
    !Function for the density maximum of H2O and D2O
    !-----------------------------------------------
    DOUBLE PRECISION module FUNCTION TDENSMAX_calc(gl,press,nrsubst)


    implicit none

    type(type_gl) :: gl


    INTEGER :: Max_Iterations, Iterations
    DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
    DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
    DOUBLE PRECISION :: Delta_allowed
    INTEGER :: errTDENMAX,nrsubst
    type(type_additional_parameters) :: ResTDENMAX_param
    DOUBLE PRECISION :: TDENMAX, press

    TDENMAX=0.d0
    !ResTDENMAX_param = 0.d0
    Tstart_min = 0.d0
    Tstart_max = 0.d0
    Tmin_allowed = 0.d0
    Tmax_allowed = 0.d0
    Delta_allowed = 0.d0
    errTDENMAX=0

    !Parameters for iteration:
    Delta_allowed = 1.d-8
    Max_Iterations = 500
    ResTDENMAX_param%a_p(1) = 0.d0
    ResTDENMAX_param%a_p(2) = nrsubst
    ResTDENMAX_param%a_p(3) = press
    !ResPsat_param(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration

    Tstart_min = 250.d0
    Tmin_allowed = 200.d0

    Tstart_max = 300.d0
    Tmax_allowed = 10000.d0

    CALL Regula_Falsi(gl,Res_TDENMAX,TDENMAX,Tstart_min,Tstart_max,Delta_allowed,&
        Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errTDENMAX,ResTDENMAX_param)

    if (errTDENMAX == 0) then
        TDENSMAX_calc=TDENMAX
    else
        TDENSMAX_calc = -5661.d0
    end if
    End Function TDENSMAX_calc


    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Res_TDENMAX(gl,T_akt,parameters)
    !--------------------------------------------------------------------------
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: T_akt, DDDT, DDDT_akt, press, D
    type(type_additional_parameters) :: parameters
    Integer :: nrsubst, IPHASE

    DDDT=parameters%a_p(1)
    nrsubst=int(parameters%a_p(2))
    press=parameters%a_p(3)

    IPHASE = 1
    D=rhomix_calc(gl,T_akt, press, 0.d0, 1, nrsubst)
    DDDT_akt = DDDT_CALC(gl,T_akt,D, nrsubst)

    Res_TDENMAX=DDDT-DDDT_akt
    !sinnvoll für seewasser?
    RETURN
    END function






    !--------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION Res_Dspin(gl,D_akt,parameters)
    !--------------------------------------------------------------------------
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters

    DOUBLE PRECISION :: D_akt, Dspin_slope, T_akt, DPDD_akt
    Integer :: i

    Dspin_slope=parameters%a_p(1)
    T_akt=parameters%a_p(2)
    i=int(parameters%a_p(3))

    DPDD_akt = DPDD_CALC(gl,T_akt,D_akt, i)

    Res_Dspin=Dspin_slope-DPDD_akt

    RETURN
    END function



    !----------------------------------------------------------------------------------------------------------------------------------
    !Function for calculating the fundamental derivative of gas dynamics according to Castier and Cabral
    !"Pure saturated gases with predicted negative fundamental derivative of gas dynamics"
    !Fluid Phase Equilibria, 334:128–136 (2012)   => unitless
    !----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION GAMMAGD_CALC(gl,T, D, nrsubst)




    implicit none

    type(type_gl) :: gl


    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    double precision:: T, D, gam1, gam2, gam3
    double precision:: DPDD, D2PDD2, DPDT, D2PDDT, D2PDT2, SND, DCVDT, CV, wmmix
    !double precision:: DPDD_CALC, D2PDD2_CALC, DPDT_CALC, D2PDDT_CALC, D2PDTT_CALC, WS_CALC, DCVDT_CALC, CV_CALC

    if (nrsubst == 0) then
        call wm_mix_calc(gl,wmmix)                 !subroutine calculates molmass of fluidmixture
    else
        wmmix = gl%wm(nrsubst)
    end if


    SND = WS_CALC(gl,T,D, nrsubst)          !m/s

    !for gam1
    DPDD = DPDD_CALC(gl,T,D, nrsubst)       ! MPa/(mol/m³)
    D2PDD2 = D2PDD2_CALC(gl,T,D, nrsubst)   ! [MPa/(mol/m³)]²

    !for gam2
    CV = CV_CALC(gl,T,D, nrsubst)           ! J/mol/K
    DPDT = DPDT_CALC(gl,T,D, nrsubst)       ! MPa/K
    D2PDDT = D2PDDT_CALC(gl,T,D, nrsubst)   ! J/mol/K

    !for gam3
    D2PDT2 = D2PDTT_CALC(gl,T,D, nrsubst)   !MPa/K²
    DCVDT = DCVDT_CALC(gl,T, D, nrsubst)    !J/mol/K²

    gam1 = 2.d0 * D**3 * DPDD + D**4 * D2PDD2       ! MPa*(mol/m³)²
    gam2 = -3.d0 * T / CV * DPDT * (-D**2 * D2PDDT)        ! MPa*(mol/m³)²
    gam3 = (T / CV * DPDT * 1.d6 / gl%factorpress)**2 * (3.d0 * D2PDT2 + 1.d0 / T * DPDT * (1.d0 - T / CV * DCVDT))  ! MPa*(mol/m³)²

    GAMMAGD_CALC = 0.D0
    GAMMAGD_CALC = D**(-3)/2.d0/snd**2/wmmix*(gam1 + gam2 + gam3)*1.d6 / gl%factorpress    !RB/TN factorpress nach Gefühl eingefügt


    END FUNCTION GAMMAGD_CALC


    !----------------------------------------------------------------------------------------------------------------------------------
    !Subroutine for calculating the thermodynamic factor according to Kooijman and Taylor
    ! "Composition derivatives of activity coefficient models",
    ! Chem. Eng. Comm. 1991, Vol. 102, pp. 87-106
    ! This property is unitless
    !
    ! M. Thol, October 2017
    !
    !----------------------------------------------------------------------------------------------------------------------------------
    module Subroutine GAMMATF_CALC(gl,T, D, GAMMATF)

    !!DEC$ ATTRIBUTES DLLEXPORT :: GAMMATF_CALC



    implicit none

    type(type_gl) :: gl


    double precision:: T, D
    double precision, dimension(30, 30), intent(out):: GAMMATF
    double precision, dimension(30, 30):: DCHPOTDX

    double precision::  Rmix, wmmix

    integer:: i, j, oir

    oir = 0                    !The overall chemical potential will be calculated (1 = ideal, 2 = residual)

    call R_mix_calc(gl,Rmix)
    call wm_mix_calc(gl,wmmix)

    call d2na_dnidxj_PT (gl,T, D, DCHPOTDX, OIR)

    GAMMATF = 0.d0
    Do j = 1, gl%NCOMP-1
        Do i = 1, gl%NCOMP
            GAMMATF(i,j) = DCHPOTDX(j,i) * gl%molfractions(i)
            !dchpotdx_res(j,i)
        end do
    end do

    END SUBROUTINE GAMMATF_CALC


    !----------------------------------------------------------------------------------------------------------------------------------
    !dCV/dT => numerical derivative of the isochoric heat capacity wrt temperature [J/mol/K²]
    !----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DCVDT_CALC(gl,T, D, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    double precision:: T, D, CVM, CVP, diff

    diff = 0.001d0

    CVM = CV_CALC(gl,T-diff,D, nrsubst)
    CVP = CV_CALC(gl,T+diff,D, nrsubst)

    DCVDT_CALC = (CVP-CVM)/(2.d0*diff)

    END FUNCTION DCVDT_CALC


    !----------------------------------------------------------------------------------------------------------------------------------
    !dCV/dT => numerical derivative of the isochoric heat capacity wrt temperature [J/mol/K²]
    !----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION DCVRDT_CALC(gl,T, D, nrsubst)

    implicit none

    type(type_gl) :: gl


    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    double precision:: T, D, CVRM, CVRP, diff

    diff = 0.001d0

    CVRM = CVR_CALC(gl,T-diff,D, nrsubst)
    CVRP = CVR_CALC(gl,T+diff,D, nrsubst)

    DCVRDT_CALC = (CVRP-CVRM)/(2.d0*diff)


    END FUNCTION DCVRDT_CALC

    !DPD2 DERIVATIVE WITH RESPECT TO D^2
    DOUBLE PRECISION module FUNCTION D4PD4_CALC(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    double precision:: T, D
    double precision:: num_delta     !numerical delta
    double precision, dimension(5):: num_dense,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_dense(1) = D
    num_dense(2) = D + num_delta
    num_dense(3) = D - num_delta
    num_dense(4) = D + 2d0*num_delta
    num_dense(5) = D - 2d0*num_delta

    DO i=1,5
        results(i) = D2PDD2_CALC(gl,T,num_dense(i),nrsubst)
    end do

    D4PD4_CALC = ( -1.d0 * results(4) + 16.d0 * results(2) - 30.d0 * results(1) + 16.d0 *results(3) - results(5) ) / (12.d0*num_delta**2)

    END FUNCTION D4PD4_CALC


    !####################################################################################################################################################################


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Second derivative of second virial coeff. B with respect to T^2
    !d2BdT2
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d2BdT2_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = B_CALC(gl,num_temp(i),nrsubst)
    end do

    d2BdT2_CALC = ( -1.d0 * results(4) + 16.d0 * results(2) - 30.d0 * results(1) + 16.d0 *results(3) - results(5) ) / (12.d0*num_delta**2)

    END FUNCTION d2BdT2_CALC


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Second derivative of third virial coeff. C with respect to T^2
    !d2CdT2
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d2CdT2_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = C_CALC(gl,num_temp(i),nrsubst)
    end do

    d2CdT2_CALC = ( -1.d0 * results(4) + 16.d0 * results(2) - 30.d0 * results(1) + 16.d0 *results(3) - results(5) ) / (12.d0*num_delta**2)

    END FUNCTION d2CdT2_CALC

    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Second derivative of fourth virial coeff. D with respect to T^2
    ! d2DdT2
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d2DdT2_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    double precision:: T
    double precision:: num_delta     !numerical delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = D_CALC(gl,num_temp(i),nrsubst)
    end do

    d2DdT2_CALC = ( -1.d0 * results(4) + 16.d0 * results(2) - 30.d0 * results(1) + 16.d0 *results(3) - results(5) ) / (12.d0*num_delta**2)

    END FUNCTION d2DdT2_CALC



    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Third derivative of second virial coeff. B with respect to T^3
    !d3BdT3
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d3BdT3_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(4):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T + num_delta
    num_temp(2) = T - num_delta
    num_temp(3) = T + 2d0*num_delta
    num_temp(4) = T - 2d0*num_delta

    DO i=1,4
        results(i) = B_CALC(gl,num_temp(i),nrsubst)
    end do

    d3BdT3_CALC = (results(3) - 2.d0 * results(1) + 2.d0 *results(2) - results(4) ) / (2.d0*num_delta**3)

    END FUNCTION d3BdT3_CALC


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Third derivative of third virial coeff. C with respect to T^3
    !d3CdT3
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d3CdT3_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(4):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T + num_delta
    num_temp(2) = T - num_delta
    num_temp(3) = T + 2d0*num_delta
    num_temp(4) = T - 2d0*num_delta

    DO i=1,4
        results(i) = C_CALC(gl,num_temp(i),nrsubst)
    end do

    d3CdT3_CALC = (results(3) - 2.d0 * results(1) + 2.d0 *results(2) - results(4) ) / (2.d0*num_delta**3)

    END FUNCTION d3CdT3_CALC


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Third derivative of fourth virial coeff. D with respect to T^3
    !d3DdT3
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d3DdT3_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(4):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_temp(1) = T + num_delta
    num_temp(2) = T - num_delta
    num_temp(3) = T + 2d0*num_delta
    num_temp(4) = T - 2d0*num_delta

    DO i=1,4
        results(i) = D_CALC(gl,num_temp(i),nrsubst)
    end do

    d3DdT3_CALC = (results(3) - 2.d0 * results(1) + 2.d0 *results(2) - results(4) ) / (2.d0*num_delta**3)

    END FUNCTION d3DdT3_CALC


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Fourth derivative of second virial coeff. B with respect to T^4
    !d4BdT4
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d4BdT4_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-6
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = B_CALC(gl,num_temp(i),nrsubst)
    end do

    d4BdT4_CALC = ( results(4) - 4.d0 * results(2) + 6.d0 * results(1) - 4.d0 *results(3) + results(5) ) / (num_delta**4)

    END FUNCTION d4BdT4_CALC



    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Fourth derivative of third virial coeff. C with respect to T^4
    !d4CdT4
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d4CdT4_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-2
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = C_CALC(gl,num_temp(i),nrsubst)
    end do

    d4CdT4_CALC = ( results(4) - 4.d0 * results(2) + 6.d0 * results(1) - 4.d0 *results(3) + results(5) ) / (num_delta**4)

    END FUNCTION d4CdT4_CALC


    ! ----------------------------------------------------------------------------------------------------------------------------------
    !Fourth derivative of fourth virial coeff. D with respect to T^4
    !d4DdT4
    ! ----------------------------------------------------------------------------------------------------------------------------------
    DOUBLE PRECISION module FUNCTION d4DdT4_CALC(gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    integer:: nrsubst
    double precision:: T
    double precision:: num_delta
    double precision, dimension(5):: num_temp,results
    integer:: i !loop variable

    num_delta  = 1.D-2
    num_temp(1) = T
    num_temp(2) = T + num_delta
    num_temp(3) = T - num_delta
    num_temp(4) = T + 2d0*num_delta
    num_temp(5) = T - 2d0*num_delta

    DO i=1,5
        results(i) = D_CALC(gl,num_temp(i),nrsubst)
    end do

    d4DdT4_CALC = ( results(4) - 4.d0 * results(2) + 6.d0 * results(1) - 4.d0 *results(3) + results(5) ) / (num_delta**4)

    END FUNCTION d4DdT4_CALC


    !Sven Pohl 28.11.2019
    !Numerical implementation of Virial coefficients
    !Second virial coefficient
    double precision  module function b_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T,D

    D=1D-8
    b_calc_num = A01_calc(gl,T,D,nrsubst)*(1d0/D)
    end function
    !Third virial coefficient
    double precision  module function c_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T,D
    D=1D-6
    c_calc_num = A02_calc(gl,T,D,nrsubst)*(1d0/D**2)
    end function

    !Fourth virial coefficient
    double precision  module function d_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T,D,diff
    double precision, dimension(500):: calcs
    double precision:: k
    logical:: converged
    converged = .false.

    D = 1d-5
    !calcs(1) = (A03_calc(gl,T,D,nrsubst)/abs(D)**3)/2.d0
    !D = 1d-4
    !calcs(2) = (A03_calc(gl,T,D,nrsubst)/D**3)/2.d0
    d_calc_num = (A03_calc(gl,T,D,nrsubst)/D**3)/2.d0
    !k=2
    !do while(.not.converged)
    !    D= D/k
    !    calcs(k) = A03_calc(gl,T,D,nrsubst)/D**3/2.d0
    !    diff = abs(calcs(k) - calcs(k-1))
    !    if(diff < 1.D-6) then
    !        d_calc_num=calcs(k)
    !        converged = .true.
    !    end if
    !    k = k + 1
    !end do
    end function

    ! Effective inverse-power-law exponent
    double precision module function n_eff_calc(gl,T,D,nrsubst)

    implicit none

    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D,Rmix
    DOUBLE PRECISION, dimension(nderivs)::FNRDER
    INTEGER,DIMENSION(nderivs),parameter :: GETDERR = (/0,1,0,0,1,1,0,0,0,0,0,0,0,0,0/)

    if (nrsubst == 0) then
        continue
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
        Rmix = gl%REQ(nrsubst)
    end if


    n_eff_calc = -3d0*(FNRDER(2)-FNRDER(6))/FNRDER(5)

    end function

    !function to prevent double loops
    ! ZETA = FNRD / log2(delta+1)
    double precision module function zeta_calc(gl,T,D,nrsubst)
    implicit none

    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D
    DOUBLE PRECISION, dimension(nderivs)::FNRDER
    INTEGER,DIMENSION(nderivs),parameter :: GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)

    if (nrsubst == 0) then
        continue
    else
        CALL FNRDERIVS(gl,T,D,GETDERR,FNRDER,nrsubst)   !subroutine calculates derivatives of residual part
    end if

    zeta_calc = FNRDER(2)/(log(D/gl%rhored(nrsubst)+1d0)/(log(2d0)))

    end function

    !function to prevent double loops
    ! d2ZETAdd2
    double precision module function d2zetadd2_calc(gl,T,D,nrsubst)
    implicit none

    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D,num_delta
    double precision, dimension(5):: num_d,results
    integer:: i !loop variable

    num_delta  = 1.D-4
    num_d(1) = D
    num_d(2) = D + num_delta
    num_d(3) = D - num_delta
    num_d(4) = D + 2d0*num_delta
    num_d(5) = D - 2d0*num_delta

    DO i=1,5
        results(i) =zeta_calc(gl,T,num_d(i),nrsubst)
    end do

    d2zetadd2_calc = ( -1.d0 * results(4) + 16.d0 * results(2) - 30.d0 * results(1) + 16.d0 *results(3) - results(5) ) / (12.d0*num_delta**2)

    end function

    end submodule impl