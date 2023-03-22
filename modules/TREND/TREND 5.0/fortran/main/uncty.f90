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

    ! module for file uncty.f90
    module uncty_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use rhomix_pt_module
    use flash_module
    use flash_pure_module


    contains




    !----------------------------------------------------------------------------------------------------------------------
    !                                                   CALL UNCERTAINTY
    !----------------------------------------------------------------------------------------------------------------------

    SUBROUTINE UNCERTAINTY (gl,fktProperty, dblT, dblD, blnTwophase, dblX_vap, calctype, errval, path)

    !**********************************************************************************************************************
    ! This subroutine returns estimated/known uncertainty.
    !
    !   INPUT:    fktProperty = calling property's function (prop_from_deriv) - > type function
    !		      dblT 	 	= temperature 						- > type double precision
    !		      dblD	 	= density							- > type double precision
    !		      blnTwophase = twophase yes/no 			    - > type logical
    !		      dblX_vap	= vapour amount						- > type double precision
    !
    !   OUTPUT:   estUncty_out = estimated or known uncertainty	- > type double precision
    !
    !**********************************************************************************************************************

    !USE module_eos_coefficients
    !USE module_uncertainty_parameters
    !USE module_fluid_parameters
    !USE module_ideal_gas_coefficients




    implicit none

    type(type_gl) :: gl

    ABSTRACT INTERFACE
    DOUBLE PRECISION FUNCTION FKTPROPERTY_PROTOTYPE(gl,T,D, nrsubst)
    use module_all_types
    type(type_gl) :: gl
    DOUBLE PRECISION :: T,D
    integer:: nrsubst
    END FUNCTION FKTPROPERTY_PROTOTYPE
    END INTERFACE
    PROCEDURE(FKTPROPERTY_PROTOTYPE), POINTER, intent(in) :: FKTPROPERTY

    double precision :: posvaried_Dens, negvaried_Dens, posvariedcoeff_Dens, negvariedcoeff_Dens, dblD_Gen
    DOUBLE PRECISION :: dblT, dblD, dblP
    DOUBLE PRECISION :: arEstUncty(5)
    DOUBLE PRECISION :: dblX_vap, dblRho_l, dblRho_v, dblUncty, dblUncty_vle
    DOUBLE PRECISION :: dblCpGen, dblRhoGen, dblWGen, dblPropGen, dblPropGen_vle, dblGen, dblGen_vle, dblPropGen_0
    DOUBLE PRECISION :: dblCpGen_vle, dblRhoGen_vle, dblWGen_vle, dblCVGen, dblCVGen_vle, dblPsat
    LOGICAL 		 :: blnTwophase, oncourse_pos, oncourse_neg, oncourse_out
    INTEGER  		 :: errPsat
    INTEGER			 :: i, intPhase, iter, iFlash, coeff_nr, variation, base_out, UNCTY_ERROR

    double precision :: PropMan_pos, PropMan_neg, estU_pos, estU_neg, estUncty, pos_varied, neg_varied
    double precision :: cpposvaried, cvposvaried, wposvaried, cpnegvaried, cvnegvaried, wnegvaried, unccp, unccv, uncd, uncw
    double precision :: coeff_variance_out
    double precision, dimension(100) :: nGen
    double precision :: rho_limit, cp_limit, cv_limit, ws_limit
    double precision :: rho_limit_varied, cp_limit_varied, cv_limit_varied, ws_limit_varied
    double precision :: rho_limit_plus, cp_limit_plus, cv_limit_plus, ws_limit_plus
    double precision :: rho_limit_minus, cp_limit_minus, cv_limit_minus, ws_limit_minus

    !Andreas April 2013
    character(255) :: path      !Path where the UNCTY files lie
    integer:: calctype          !Single Property uncertainty (H=1, U=2, CV=3, S=4, CP=5, WS=6, P=7, D=8)
    integer:: errval


    gl%uncty%blnCheck2Phase = .true.
    dblCPGen_vle = 0.D0
    dblRhoGen_vle = 0.D0
    dblWGen_vle = 0.D0
    dblCVGen_vle = 0.D0
    dblPropGen_vle = 0.D0
    arEstUncty = 0.D0
    dblPsat = 0.D0
    dblRho_v = 0.D0
    dblRho_l = 0.D0
    nGen = 0.D0
    estU_pos = 0.D0
    estU_neg = 0.D0
    estUncty = 0.D0
    gl%uncty%estUncty_out = 0.D0
    posvaried_Dens = 0.D0
    negvaried_Dens = 0.D0
    variation = 0
    UNCTY_ERROR = 0
    errval = 0

    !save genuine value
    dblD_Gen = dblD
    dblPropGen_0 = fktProperty(gl, dblT, dblD, 1)


    !T given, calculate psat
    iFlash = 1


    !get values for Psat, Rho_v and Rho_l
    !if(dblT < tc(1)) then
    !call VLEpure(dblT,dblPsat,dblRho_v,dblRho_l,errPsat,1)
    call Flash_Pure_PhaseBoundary (gl,dblPsat, dblT, dblRho_v, dblRho_l, iFlash, errPsat, iter, 1)
    if (errPsat /= 0) then
        dblPsat = vp_eq(gl,dblT, 1)
    end if
    !end if


    !from here on a check for twophase yes/no will lead to incorrect results near the phase boundary because of manipulated coefficients
    gl%uncty%blnCheck2Phase = .false.


    !if in twophase region p=vapour pressure | calculate properties on phase boundaries
    if (blnTwophase) then !in twophase region
        dblP = dblPsat
        dblCpGen = CP_CALC(gl,dblT, dblRho_l,1)
        dblCPGen_vle = CP_CALC(gl,dblT, dblRho_v, 1)
        dblRhoGen = dblRho_l
        dblRhoGen_vle = dblRho_v
        dblWGen = WS_CALC(gl,dblT, dblRho_l, 1)
        dblWGen_vle = WS_CALC(gl,dblT, dblRho_v, 1)
        dblCVGen = CV_CALC(gl,dblT, dblRho_l, 1)
        dblCVGen_vle = CV_CALC(gl,dblT, dblRho_v,1)
        dblPropGen = fktProperty(gl, dblT, dblRho_l, 1)
        dblPropGen_vle = fktProperty(gl, dblT, dblRho_v, 1)

    else !in singlephase region
        dblP = P_CALC(gl,dblT, dblD, 1)
        dblCPGen = CP_CALC(gl,dblT, dblD, 1)
        dblRhoGen = dblD
        dblWGen = WS_CALC(gl,dblT, dblD, 1)
        dblCVGen = CV_CALC(gl,dblT, dblD, 1)
        dblPropGen = fktProperty(gl, dblT, dblD, 1)

        !in singlephase only the liquid loop will be used, so dblRho_l is singlephase density
        dblRho_l = dblD
    end if


    !get uncertainties from UNCTY-File
    call uncertainty_area(gl,dblT, dblP, dblPsat, errval, path)
    !Andreas April 2013
    if (errval /= 0) then
        return
    end if

    !for consideration of equation's limits
    rho_limit = rhomix_calc(gl,gl%tmaxfluid(1), gl%pmaxfluid(1), 0.D0, 0, 0)
    rho_limit_plus = rho_limit*(1.D0 + gl%uncty%rho_limit_uncty/100.D0)
    rho_limit_minus = rho_limit*(1.D0 - gl%uncty%rho_limit_uncty/100.D0)

    if (gl%uncty%cp_limit_uncty /= 0.D0) then
        cp_limit = CP_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
        cp_limit_plus = cp_limit*(1.D0 + gl%uncty%cp_limit_uncty/100.D0)
        cp_limit_minus = cp_limit*(1.D0 - gl%uncty%cp_limit_uncty/100.D0)
    end if

    if (gl%uncty%cv_limit_uncty /= 0.D0) then
        cv_limit = CV_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
        cv_limit_plus = cv_limit*(1.D0 + gl%uncty%cv_limit_uncty/100.D0)
        cv_limit_minus = cv_limit*(1.D0 - gl%uncty%cv_limit_uncty/100.D0)
    end if

    if (gl%uncty%ws_limit_uncty /= 0.D0) then
        ws_limit = WS_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
        ws_limit_plus = ws_limit*(1.D0 + gl%uncty%ws_limit_uncty/100.D0)
        ws_limit_minus = ws_limit*(1.D0 - gl%uncty%ws_limit_uncty/100.D0)
    end if


    !if the user choses a property with known uncertainty - >  there is no need to calculate anything
    !this works only for cp and w, yet
    if ((fktProperty(gl, dblT,dblD,1) - CP_CALC(gl,dblT,dblD,1) == 0) .and. (gl%uncty%dblUncty_CP > 0)) then
        gl%uncty%estUncty_out = gl%uncty%dblUncty_CP
        gl%uncty%Uncty_cp = gl%uncty%dblUncty_CP / 100.D0 * dblCPGen
        !write(29,*) 'Known uncertainty for cp: ', estUncty_out, 'at T =', dblT
        gl%uncty%blnCheck2Phase = .true.
        RETURN

    else if ((fktProperty(gl, dblT,dblD,1) - CV_CALC(gl,dblT,dblD,1) == 0) .and. (gl%uncty%dblUncty_CV > 0)) then
        gl%uncty%estUncty_out = gl%uncty%dblUncty_CV
        gl%uncty%Uncty_cv = gl%uncty%dblUncty_CV  / 100.D0 * dblCVGen
        !write(29,*) 'Known uncertainty for cv: ', estUncty_out, 'at T =', dblT
        gl%uncty%blnCheck2Phase = .true.
        RETURN

    else if ((fktProperty(gl, dblT,dblD,1) - WS_CALC(gl,dblT,dblD,1) == 0) .and. (gl%uncty%dblUncty_WS > 0)) then
        gl%uncty%estUncty_out = gl%uncty%dblUncty_WS
        gl%uncty%Uncty_w = gl%uncty%dblUncty_ws  / 100.D0 * dblWGen
        !write(29,*) 'Known uncertainty for ws: ', estUncty_out, 'at T =', dblT
        gl%uncty%blnCheck2Phase = .true.
        RETURN
    end if


    !If the temperature AND Density has not changed, the coefficient manipulation does not have to be carried out but the last one is taken
    !Andreas April 2013
    if ((abs(gl%uncty%temp_last - dblT) > 1.D-8) .or. (abs(gl%uncty%dens_last - dblD) > 1.D-8)) then

        !get estimation of uncertainty for each base
        do i=1, 5
            arEstUncty(i) = 0.D0
            !        arEstUncty = 0.D0
            if (i == 1) then	        !- >  cp
                dblGen = dblCPGen
                dblGen_vle = dblCPGen_vle
                dblUncty = gl%uncty%dblUncty_CP/100
                dblUncty_vle = gl%uncty%dblUncty_CP_vle
            else if (i == 2) then	!- >  rho
                dblGen = dblRhoGen
                dblGen_vle = dblRhoGen_vle
                dblUncty = gl%uncty%dblUncty_D/100
                !Andreas April 2013
                !-------------------
                ! Absolute difference
                gl%uncty%Uncty_d = dblUncty * dblRhoGen
                !-------------------
                dblUncty_vle = gl%uncty%dblUncty_D_vle
            else if (i == 3) then	!- >  ws
                dblGen = dblWGen
                dblGen_vle = dblWGen_vle
                dblUncty = gl%uncty%dblUncty_WS/100
                dblUncty_vle = gl%uncty%dblUncty_WS_vle
            else if (i == 4) then	!- >  p
                dblGen = dblP
                dblGen_vle = dblP
                dblUncty = gl%uncty%dblUncty_P/100
                dblUncty_vle = gl%uncty%dblUncty_P_vle
            else if (i == 5) then	!- >  cv
                dblGen = dblCVGen
                dblGen_vle = dblCVGen_vle
                dblUncty = gl%uncty%dblUncty_CV/100
                dblUncty_vle = gl%uncty%dblUncty_CV_vle
            end if


            !only if there is a known uncertainty - >  continue | else, start with next base
            if (dblUncty /= 0.D0) then

                !these two function calls will return estimated upper/lower(+/- uncertainty) value of target property in single phase as well
                !as two phase. If in twhopase it will return saturated liquid estimation.
                intPhase = 1

                pos_varied = dblGen +(dblGen*dblUncty)            !highest expected Value for dblGen with the given uncertainty, ct
                neg_varied = dblGen -(dblGen*dblUncty)            !lowest expected Value for dblGen with the given uncertainty, ct

                call coeff_manipulation(gl,fktProperty, pos_varied, dblT, dblRho_l, dblP, i, dblPropGen, intPhase, dblPsat, dblGen, dblUncty)

                call coeff_manipulation(gl,fktProperty, neg_varied, dblT, dblRho_l, dblP, i, dblPropGen, intPhase, dblPsat, dblGen, dblUncty)


                !        !what results in a bigger deviation, upper Bound or lower Bound (regularly only little difference between them)
                !        if (abs(dblPropGen - dblUpEst) > abs(dblPropGen - dblLowEst)) then
                !          	dblEst = dblUpEst
                !       	   else
                !            dblEst = dblLowEst
                !        end if
                !
                !		!if in twophase do the same for saturated vapour
                !        if (blnTwophase) then
                !
                !            intPhase = 2
                !        	dblUpEst_vle = coeff_manipulation(fktProperty, (dblGen_vle+(dblGen_vle*dblUncty_vle)), dblT, dblRho_v, dblP, i,&
                !            dblPropGen_vle, intPhase, dblPsat)
                !
                !        	dblLowEst_vle = coeff_manipulation(fktProperty, (dblGen_vle-(dblGen_vle*dblUncty_vle)), dblT, dblRho_v, dblP, i,&
                !            dblPropGen_vle, intPhase, dblPsat)
                !
                !        	if (abs(dblPropGen_vle - dblUpEst_vle) > abs(dblPropGen_vle - dblLowEst_vle)) then
                !          		dblEst_vle = dblUpEst_vle
                !       	   		else
                !            	dblEst_vle = dblLowEst_vle
                !        	end if
                !
                !			!save the estimated uncertainties
                !			arEstUncty(i) = abs((dblPropGen_0 - (dblEst + dblX_vap*(dblEst_vle-dblEst))) / dblPropGen_0 * 100.D0)
                !
                !         else
                !           !save the estimated uncertainties
                !           arEstUncty(i) = abs((dblPropGen - dblEst) / dblPropGen * 100.D0)
                !       	end if

            end if
        end do

        !Andreas April 2013
        !Save the Temperature and density for which the pos and neg varied coefficients was already calculated
        gl%uncty%temp_last = dblT
        gl%uncty%dens_last = dblD

    end if

    !save genuine ni_u in nGen
    nGen(:) = gl%eos_coeff%ni(:,1)

    !das folgende wurde der Übersicht halber zunächst auskommentiert
    !            cp0coeff_Gen = cp0coeff
    !
    !    do i=1,20
    !
    !        !exchange genuine by POSITIVE varied ni_u
    !        if (ni_varied_pos_ideal(i) /= 0.d0) then
    !        	CP0COEFF(i,1) = ni_varied_pos_ideal(i)
    !
    !        dblD = rhomix_calc(gl,dblT, dblP, dblD_Gen, 0, 0)
    !
    !        !calculate target property with this manipulated ni_u
    !            PropMan_pos = fktProperty(dblT, dblD, 1)
    !
    !        !save delta between genuine and manipulated target property
    !            estU_pos = abs(dblPropGen - PropMan_pos)/dblPropGen * 100.D0
    !        end if
    !
    !
    !        !exchange genuine by NEGATIVE varied ni_u
    !        if (ni_varied_neg_ideal(i) /= 0.D0) then
    !        	CP0COEFF(i,1) = ni_varied_neg_ideal(i)
    !
    !        dblD = rhomix_calc(gl,dblT, dblP, dblD_Gen, 0, 0)
    !
    !        !calculate target property with this manipulated ni_u
    !            PropMan_neg = fktProperty(dblT, dblD, 1)
    !
    !        !save delta between genuine and manipulated target property
    !            estU_neg = abs(dblPropGen - PropMan_neg)/dblPropGen * 100.D0
    !        end if
    !
    !
    !        ! control if the estimated uncertainty of the positive varied ni_u (estU_pos) are bigger or the one of the negative varied ni_u (estU_neg)
    !        ! and save it as estUncty if it's bigger than the previous value
    !            if ((estU_pos > estU_neg) .and. (estU_pos > estUncty)) then
    !                estUncty = estU_pos
    !
    !            else if ((estU_neg > estU_pos) .and. (estU_neg > estUncty)) then
    !                estUncty = estU_neg
    !            else
    !                !do nothing (- >  leave estUncty as it was)
    !!                estUncty = 0.D0
    !            end if
    !
    !
    !        !reset ni_u
    !           cp0coeff = cp0coeff_Gen
    !            estU_pos = 0.D0
    !            estU_neg = 0.D0
    !
    !
    !    end do

    do i=1,100

        !exchange genuine by POSITIVE varied ni_u
        if (gl%uncty%ni_varied_pos(i) /= 0.d0) then
            gl%eos_coeff%ni(i,1) = gl%uncty%ni_varied_pos(i)

            dblD = rhomix_calc(gl,dblT, dblP, dblD_Gen, 0, 0)

            if (dblD == 0.D0) then
                !Da muss was schiefgelaufen sein!!!
                estU_pos = 0.D0

            else
                posvariedcoeff_Dens = dblD

                uncd =  abs(dblrhoGen - posvariedcoeff_Dens) / posvariedcoeff_Dens
                if (uncd >= gl%uncty%dblUNCTY_D) then
                    UNCTY_ERROR = 301
                end if

                if (gl%uncty%dblUncty_cp /= 0.D0) then
                    cpposvaried = CP_CALC(gl,dblT, dblD_Gen,1)
                    unccp = abs(dblCpGen - cpposvaried) / cpposvaried
                    if (unccp > gl%uncty%dblUncty_cp) then
                        UNCTY_ERROR = 302
                    end if
                end if

                if (gl%uncty%dblUncty_cv /= 0.D0) then
                    cvposvaried = CV_CALC(gl,dblT, dblD_Gen,1)
                    unccv = abs(dblCvGen - cvposvaried) / cvposvaried
                    if (unccv > gl%uncty%dblUncty_cv) then
                        UNCTY_ERROR = 303
                    end if
                end if


                if (gl%uncty%dblUncty_ws /= 0.D0) then
                    wposvaried = WS_CALC(gl,dblT, dblD_Gen,1)
                    uncw = abs(dblwGen - wposvaried) / wposvaried
                    if (uncw > gl%uncty%dblUncty_ws) then
                        UNCTY_ERROR = 304
                    end if
                end if

                !testing the limit of the equation, if there are too big discrepancies
                rho_limit_varied = rhomix_calc(gl,gl%tmaxfluid(1), gl%pmaxfluid(1), 0.D0, 0, 0)
                if ((rho_limit_varied > rho_limit_plus) .or. (rho_limit_varied < rho_limit_minus)) then
                    oncourse_pos = .false.
                end if

                oncourse_pos = .true.

                if ((gl%uncty%cp_limit_uncty /= 0.D0) .and. (oncourse_pos)) then
                    cp_limit_varied = CP_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((cp_limit_varied > cp_limit_plus) .or. (cp_limit_varied < cp_limit_minus)) then
                        oncourse_pos = .false.
                    end if
                end if

                if ((gl%uncty%cv_limit_uncty /= 0.D0) .and. (oncourse_pos)) then
                    cv_limit_varied = CV_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((cv_limit_varied > cv_limit_plus) .or. (cv_limit_varied < cv_limit_minus)) then
                        oncourse_pos = .false.
                    end if
                end if

                if ((gl%uncty%ws_limit_uncty /= 0.D0) .and. (oncourse_pos)) then
                    ws_limit_varied = WS_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((ws_limit_varied > ws_limit_plus) .or. (ws_limit_varied < ws_limit_minus)) then
                        oncourse_pos = .false.
                    end if
                end if


                !IMPORTANT: THE if statement was wrong!! if NO error occurs, calculate the uncertainty
                if (UNCTY_ERROR == 0) then
                    !calculate target property with this manipulated ni_u
                    PropMan_pos = fktProperty(gl, dblT, dblD, 1)
                    !save delta between genuine and manipulated target property
                    !Changed this to absolute differences, Andreas April 2013
                    estU_pos = abs(dblPropGen - PropMan_pos)!/dblPropGen * 100.D0

                else
                    estU_pos = 0.D0
                end if
            end if
        end if

        !exchange genuine by NEGATIVE varied ni_u
        if (gl%uncty%ni_varied_neg(i) /= 0.D0) then
            gl%eos_coeff%ni(i,1) = gl%uncty%ni_varied_neg(i)

            dblD = rhomix_calc(gl,dblT, dblP, dblD_Gen, 0, 0)

            if (dblD == 0.D0) then
                !Da muss was schiefgelaufen sein!!!
                estU_neg = 0.D0

            else
                negvariedcoeff_Dens = dblD

                uncd =  abs(dblrhoGen - negvariedcoeff_Dens) / negvariedcoeff_Dens
                if (uncd >= gl%uncty%dblUNCTY_D) then
                    UNCTY_ERROR = 301
                end if

                if ((gl%uncty%dblUncty_cp /= 0.D0) .and. (UNCTY_ERROR == 0)) then
                    cpnegvaried = CP_CALC(gl,dblT, dblD_Gen,1)
                    unccp = abs(dblCpGen - cpnegvaried) / cpnegvaried
                    if (unccp > gl%uncty%dblUncty_cp) then
                        UNCTY_ERROR = 302
                    end if
                end if

                if ((gl%uncty%dblUncty_cv /= 0.D0) .and. (UNCTY_ERROR == 0)) then
                    cvnegvaried = CV_CALC(gl,dblT, dblD_Gen,1)
                    unccv = abs(dblCvGen - cvnegvaried) / cvnegvaried
                    if (unccv > gl%uncty%dblUncty_cv) then
                        UNCTY_ERROR = 303
                    end if
                end if

                if ((gl%uncty%dblUncty_ws /= 0.D0) .and. (UNCTY_ERROR == 0)) then
                    wnegvaried = WS_CALC(gl,dblT, dblD_Gen,1)
                    uncw = abs(dblwGen - wnegvaried) / wnegvaried
                    if (uncw > gl%uncty%dblUncty_ws) then
                        UNCTY_ERROR = 304
                    end if
                end if

                !testing the limit of the equation, if there are too big discrepancies
                rho_limit_varied = rhomix_calc(gl,gl%tmaxfluid(1), gl%pmaxfluid(1), 0.D0, 0, 0)
                if ((rho_limit_varied > rho_limit_plus) .or. (rho_limit_varied < rho_limit_minus)) then
                    oncourse_neg = .false.
                end if

                oncourse_neg = .true.

                if ((gl%uncty%cp_limit_uncty /= 0.D0) .and. (oncourse_neg)) then
                    cp_limit_varied = CP_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((cp_limit_varied > cp_limit_plus) .or. (cp_limit_varied < cp_limit_minus)) then
                        oncourse_neg = .false.
                    end if
                end if

                if ((gl%uncty%cv_limit_uncty /= 0.D0) .and. (oncourse_neg)) then
                    cv_limit_varied = CV_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((cv_limit_varied > cv_limit_plus) .or. (cv_limit_varied < cv_limit_minus)) then
                        oncourse_neg = .false.
                    end if
                end if

                if ((gl%uncty%ws_limit_uncty /= 0.D0) .and. (oncourse_neg)) then
                    ws_limit_varied = WS_CALC(gl,gl%tmaxfluid(1), rho_limit, 1)
                    if ((ws_limit_varied > ws_limit_plus) .or. (ws_limit_varied < ws_limit_minus)) then
                        oncourse_neg = .false.
                    end if
                end if

                if (UNCTY_ERROR == 0) then
                    !calculate target property with this manipulated ni_u
                    PropMan_neg = fktProperty(gl, dblT, dblD, 1)
                    !save delta between genuine and manipulated target property
                    !Changed this to absolute differences, Andreas April 2013
                    estU_neg = abs(dblPropGen - PropMan_neg)!/dblPropGen * 100.D0

                else
                    estU_neg = 0.D0
                end if
            end if
        end if

        !reset ni_u
        gl%eos_coeff%ni(i,1) = nGen(i)

        ! control if the estimated uncertainty of the positive varied ni_u (estU_pos) are bigger or the one of the negative varied ni_u (estU_neg)
        ! and save it as estUncty if it's bigger than the previous value
        if ((estU_pos > estU_neg) .and. (estU_pos > estUncty)) then
            estUncty = estU_pos
            coeff_nr = i
            coeff_variance_out = gl%uncty%used_coeff_variance_pos(i)
            base_out = gl%uncty%used_base_pos(i)
            posvaried_Dens = posvariedcoeff_Dens
            variation = 1
            oncourse_out = oncourse_pos
        else if ((estU_neg > estU_pos) .and. (estU_neg > estUncty)) then
            estUncty = estU_neg
            coeff_nr = i
            coeff_variance_out = gl%uncty%used_coeff_variance_neg(i)
            base_out = gl%uncty%used_base_neg(i)
            negvaried_Dens = negvariedcoeff_Dens
            variation = 2
            oncourse_out = oncourse_neg
        else
            !leave estUncty as it was
            !estUncty = 0.D0
        end if

        estU_pos = 0.D0
        estU_neg = 0.D0
        dblD = dblD_Gen
    end do

    ! save the biggest found uncertainty
    ! Andreas April 2013
    ! (H=1, U=2, CV=3, S=4, CP=5, WS=6, P=7  or for "all props" =100?)
    ! estUncty_out = estUncty
    select case (calctype)
    case (1)
        gl%uncty%Uncty_h = estUncty !/ 100.D0 * dblPropGen
    case (2)
        gl%uncty%Uncty_u = estUncty !/ 100.D0 * dblPropGen
    case (3)
        gl%uncty%Uncty_cv = estUncty !/ 100.D0 * dblPropGen
    case (4)
        gl%uncty%Uncty_s = estUncty !/ 100.D0 * dblPropGen
    case (5)
        gl%uncty%Uncty_cp = estUncty !/ 100.D0 * dblPropGen
    case (6)
        gl%uncty%Uncty_w = estUncty !/ 100.D0 * dblPropGen
    case (7)
        gl%uncty%Uncty_p = estUncty !/ 100.D0 * dblPropGen
    end select


    !! output in .txt-File
    !! variation = 1  - >   the estimation of the uncertainty is caused by a positive varied coefficient
    !if (variation == 1) then
    !1001 FORMAT ('T[K]:', F6.1, 2x, 'p[MPa]:', F8.4, 2x, 'u_cp:', F8.5, 2x, 'u_rho:', F8.5, 2x, &
    !            &'u_w:', F8.5, 2x, 'u_cv:', F9.5, 2x, 'D+:', F17.9,2x, 'D_posvar:', F17.9, 2x, &
    !            &'D-:', F17.9,2x, 'i:', I3, 2x, 'coeff_variance:', f7.4, 2x, 'base:', I2, 2x, L2)
    !write (29,1001) dblT, dblP, dblUncty_CP, dblUncty_D, dblUncty_WS, estUncty_out, (dblD_Gen*(1.D0+(dblUncty_D/100))), &
    !                & posvaried_Dens, dblD_Gen*(1.D0-(dblUncty_D/100)), coeff_nr, coeff_variance_out, base_out, oncourse_out
    !
    !! variation = 2  - >   the estimation of the uncertainty is caused by a negative varied coefficient
    !else if (variation == 2) then
    !1002 FORMAT ('T[K]:', F6.1, 2x, 'p[MPa]:', F8.4, 2x, 'u_cp:', F8.5, 2x, 'u_rho:', F8.5, 2x, &
    !            &'u_w:', F8.5, 2x, 'u_cv:', F9.5, 2x, 'D+:', F17.9,2x, 'D_negvar:', F17.9, 2x, &
    !            &'D-:', F17.9,2x, 'i:', I3, 2x, 'coeff_variance:', f7.4, 2x, 'base:', I2, 2x, L2)
    !write (29,1002) dblT, dblP, dblUncty_CP, dblUncty_D, dblUncty_WS, estUncty_out, (dblD_Gen*(1.D0+(dblUncty_D/100))), &
    !                & negvaried_Dens, dblD_Gen*(1.D0-(dblUncty_D/100)), coeff_nr, coeff_variance_out, base_out, oncourse_out
    !
    !! there is no coefficient found
    !else
    !1003 FORMAT ('T[K]:', F6.1, 2x, 'p[MPa]:', F8.4, 2x, 'u_cp:', F8.5, 2x, 'u_rho:', F8.5, 2x, &
    !            &'u_w:', F8.5, 2x, 'u_cv:', F9.5, 2x, 'D+:', F17.9,2x, 'D_posvar:', F17.9, 2x, &
    !            &'D-:', F17.9,2x, 'i:', I3, 2x, 'coeff_variance:', f7.4, 2x, 'base:', I2, 2x, L2)
    !write (29,1003) dblT, dblP, dblUncty_CP, dblUncty_D, dblUncty_WS, estUncty_out, (dblD_Gen*(1.D0+(dblUncty_D/100))), &
    !                & posvaried_Dens, dblD_Gen*(1.D0-(dblUncty_D/100)), coeff_nr, coeff_variance_out, base_out, oncourse_out
    !end if

    gl%uncty%blnCheck2Phase = .true.

    END SUBROUTINE UNCERTAINTY


    !subroutine will be called by property function(prop_from_deriv) if third attibute, blnUncty = .true., then it will call subroutine UNCERTAINTY
    !this step ist necessary to pass on the target function
    subroutine callUncty(gl,intFunction, T, D, twophase, x_vap, errval, path)




    implicit none

    type(type_gl) :: gl


    INTEGER :: intFunction, errval
    DOUBLE PRECISION :: T, D, x_vap
    LOGICAL :: twophase

    !Andreas April 2013
    character (255) :: path

    if (intFunction == 1) then
        call UNCERTAINTY(gl,H_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 2) then
        call UNCERTAINTY(gl,U_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 3) then
        call UNCERTAINTY(gl,CV_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 4) then
        call UNCERTAINTY(gl,S_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 5) then
        call UNCERTAINTY(gl,CP_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 6) then
        call UNCERTAINTY(gl,WS_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 7) then
        call UNCERTAINTY(gl,P_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    else if (intFunction == 100) then
        call UNCERTAINTY(gl,H_CALC, T, D, twophase, x_vap, intFunction, errval, path)
    end if


    end subroutine callUncty


    !----------------------------------------------------------------------------------------------------------------------
    !                                             READING AREA OF KNOWN UNCERTAINTIES
    !----------------------------------------------------------------------------------------------------------------------

    SUBROUTINE uncertainty_area(gl,dblT, dblP, dblPsat, errval, path)

    !USE module_uncertainty_parameters
    !USE module_fluid_parameters


    implicit none

    type(type_gl) :: gl


    integer :: area_no, area_type, stat, errval1, nrsubst, errval, eqn
    double precision :: dblT, dblP, sum_vector_y
    double precision :: temp_begin, temp_end, pres_begin, pres_end
    double precision :: pointA_T, pointA_p, pointB_T, pointB_p, pointC_T, pointC_p, uncty_over, uncty_under
    double precision :: dblT_adjusted, T_transition, point_T, point_p, T_variation_downwards, T_variation_upwards
    double precision :: p_variation_downwards, p_variation_upwards
    double precision :: rhoL, rhoV, dblPsat, dblPsat_adjusted
    double precision, dimension(60, 60) :: matrix_A
    double precision, dimension(60) :: vector_y
    character(10) :: property, order, dummy, uncty_limit
    character(255) :: errval2

    character(255) :: path
    character(255) :: Uncty_File_path

    integer:: iFlash, iter, UNCTY_ERROR, uncty_unit

    nrsubst = 1
    gl%uncty%uncty = 0.D0
    matrix_A = 0.D0
    vector_y = 0.D0
    eqn = 2
    errval = 0
    uncty_unit = 0

    !changed for compatibility with gfortran
    !Uncty_File_path = trim(path) // "Uncty_Files\" // trim(components(1)) // ".uncty"
    Uncty_File_path = trim(path) // "Uncty_Files/" // trim(gl%components(1)) // ".uncty"
    !open uncertainty file
    open (newunit=uncty_unit, file=Uncty_File_path, iostat=stat, action='read')

    if (stat /= 0) then
        UNCTY_ERROR = -7001 !cannot find the UNCTY-FILE
        errval = UNCTY_ERROR
        RETURN
    end if


    !read uncertainty data
    read (uncty_unit,*) dummy   !Name of fluid and literature reference of the uncertainties
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy
    read (uncty_unit,*) dummy

    do
        read (uncty_unit,*) dummy
        read (uncty_unit,*) property

        do
            read (uncty_unit,*) area_no
            read (uncty_unit,*) area_type
            !case differentiation:
            !   1 = certain point with its surrounding
            !   3 = triangle
            !   4 = rectangle
            !   5 = rectangle with phase boundary line whithin
            !   6 = same as 5 but with a transition area beside the phase boundary line

            !   1 = certain point with its surrounding   ########################################################
            if (area_type == 1) then
                read (uncty_unit,*) point_T
                read (uncty_unit,*) point_p
                read (uncty_unit,*) T_variation_downwards
                read (uncty_unit,*) T_variation_upwards
                read (uncty_unit,*) p_variation_downwards
                read (uncty_unit,*) p_variation_upwards
                read (uncty_unit,*) gl%uncty

                temp_begin = point_T - T_variation_downwards
                temp_end = point_T + T_variation_upwards
                pres_begin = point_p - p_variation_downwards
                pres_end = point_p + p_variation_upwards

                !is the entered point in this area skip the next areas until the next property starts
                if ((dblT >= temp_begin) .and. (dblT <= temp_end) .and. (dblP >= pres_begin) .and. (dblP <= pres_end)) then
                    do
                        read(uncty_unit,*) order
                        if ((order == 'next_prop') .or. (order == 'end')) exit
                    end do
                else
                    read(uncty_unit,*) order
                end if

                !   3 = triangle   ##################################################################################
            else if (area_type == 3) then
                read (uncty_unit,*) pointA_T
                read (uncty_unit,*) pointA_p
                read (uncty_unit,*) pointB_T
                read (uncty_unit,*) pointB_p
                read (uncty_unit,*) pointC_T
                read (uncty_unit,*) pointC_p
                read (uncty_unit,*) gl%uncty

                !the following part describes a system of equations (and its solution), which can check if the entered point is in the triangle ABC or not:
                matrix_A(1, 1) = pointB_T - pointA_T
                matrix_A(2, 1) = pointC_T - pointA_T
                matrix_A(1, 2) = pointB_p - pointA_p
                matrix_A(2, 2) = pointC_p - pointA_p
                vector_y(1) = dblT - pointA_T
                vector_y(2) = dblP - pointA_p
                call LUdecomp(gl,eqn, matrix_A, vector_y, errval1, errval2)
                sum_vector_y = vector_y(1) + vector_y(2)

                !is the entered point in this area skip the next areas until the next property starts
                if ((vector_y(1) <= 1.D0) .and. (vector_y(1) >= 0.D0) .and. (vector_y(2) <= 1.D0) .and. (vector_y(2) >= 0.D0)&
                    & .and. (sum_vector_y <= 1.D0)) then
                    do
                        read(uncty_unit,*) order
                        if ((order == 'next_prop') .or. (order == 'end')) exit
                    end do
                else
                    read(uncty_unit,*) order
                end if

                !   4 = rectangle   #################################################################################
            else if (area_type == 4) then
                read (uncty_unit,*) temp_begin
                read (uncty_unit,*) temp_end
                read (uncty_unit,*) pres_begin
                read (uncty_unit,*) pres_end
                read (uncty_unit,*) gl%uncty

                !is the entered point in this area skip the next areas until the next property starts
                if ((dblT >= temp_begin) .and. (dblT <= temp_end) .and. (dblP >= pres_begin) .and. (dblP <= pres_end)) then
                    do
                        read(uncty_unit,*) order
                        if ((order == 'next_prop') .or. (order == 'end')) exit
                    end do
                else
                    read(uncty_unit,*) order
                end if

                !   5 = rectangle with phase boundary line whithin   ################################################
            else if (area_type == 5) then
                read (uncty_unit,*) temp_begin
                read (uncty_unit,*) temp_end
                read (uncty_unit,*) pres_begin
                read (uncty_unit,*) pres_end
                read (uncty_unit,*) uncty_over
                read (uncty_unit,*) uncty_under

                !is the entered point in this area skip the next areas until the next property starts
                if ((dblT >= temp_begin) .and. (dblT <= temp_end) .and. (dblP >= pres_begin) .and. (dblP <= pres_end)) then

                    !the differantiation, if it's over or under the phase boundary line:
                    if (dblP >= dblPsat) then
                        gl%uncty%uncty = uncty_over
                    else if (dblP < dblPsat) then
                        gl%uncty%uncty = uncty_under
                    end if

                    do
                        read(uncty_unit,*) order
                        if ((order == 'next_prop') .or. (order == 'end')) exit
                    end do
                else
                    read(uncty_unit,*) order
                end if

                !   6 = same as 5 but with a transition area beside the phase boundary line   #######################
            else if (area_type == 6) then
                read (uncty_unit,*) temp_begin
                read (uncty_unit,*) temp_end
                read (uncty_unit,*) pres_begin
                read (uncty_unit,*) pres_end
                read (uncty_unit,*) T_transition
                read (uncty_unit,*) uncty_over
                read (uncty_unit,*) uncty_under

                !is the entered point in this area skip the next areas until the next property starts
                if ((dblT >= temp_begin) .and. (dblT <= temp_end) .and. (dblP >= pres_begin) .and. (dblP <= pres_end)) then

                    dblT_adjusted = dblT - T_transition

                    !the differantiation, if it's over or under the phase boundary line:
                    !call VLEpure(dblT_adjusted,dblPsat_adjusted,RhoV,RhoL,errval,nrsubst)
                    iFlash = 1 !T given, calculate psat
                    RhoV = 0.D0
                    RhoL = 0.D0
                    call Flash_Pure_PhaseBoundary (gl,dblT_adjusted, dblPsat_adjusted, RhoV, RhoL, iFlash, errval, iter, nrsubst)

                    if (errval /= 0) then
                        dblPsat_adjusted = vp_eq(gl,dblT_adjusted, 1)
                    end if
                    if (dblP >= dblPsat_adjusted) then
                        gl%uncty%uncty = uncty_over
                    else if (dblP < dblPsat_adjusted) then
                        gl%uncty%uncty = uncty_under
                    end if

                    do
                        read(uncty_unit,*) order
                        if ((order == 'next_prop') .or. (order == 'end')) exit
                    end do
                else
                    read(uncty_unit,*) order
                end if

            else
                UNCTY_ERROR = 200 !no given Uncertainties for this point
            end if

            if ((order == 'next_prop') .or. (order == 'end')) exit
        end do

        if (property == 'd') then
            gl%uncty%dbluncty_d = gl%uncty%uncty
        else if (property == 'p') then
            gl%uncty%dbluncty_p = gl%uncty%uncty
        else if (property == 'cp') then
            gl%uncty%dbluncty_cp = gl%uncty%uncty
        else if (property == 'cv') then
            gl%uncty%dbluncty_cv = gl%uncty%uncty
        else if (property == 'ws') then
            gl%uncty%dbluncty_ws = gl%uncty%uncty
        end if

        if (order == 'end') then
            read (uncty_unit,*) dummy
            read (uncty_unit,*) dummy
            read (uncty_unit,*) dummy
            read (uncty_unit,*) dummy
            read (uncty_unit,*) dummy
            do
                read (uncty_unit,*) uncty_limit
                if (uncty_limit == 'd') then
                    read (uncty_unit,*) gl%uncty%rho_limit_uncty
                else if (uncty_limit == 'cp') then
                    read (uncty_unit,*) gl%uncty%cp_limit_uncty
                else if (uncty_limit == 'cv') then
                    read (uncty_unit,*) gl%uncty%cv_limit_uncty
                else if (uncty_limit == 'ws') then
                    read (uncty_unit,*) gl%uncty%ws_limit_uncty
                else if (uncty_limit == 'close') then
                    close(uncty_unit)
                    exit
                end if
            end do
            exit
        end if

    end do

    END SUBROUTINE uncertainty_area


    !----------------------------------------------------------------------------------------------------------------------
    !                                                COEFFICIENT MANIPULATION
    !----------------------------------------------------------------------------------------------------------------------


    subroutine coeff_manipulation(gl,fktProperty, dblBound, dblT, dblD, dblP, intBase, dblPropGen, intPhase, &
        & dblPsat, dblGen, dblUncty)

    !______________________________________________________________________________________________________________________
    !
    !    INPUT: fktProperty = target property function
    !           dblBound	= genuinie base +/- uncertainty
    ! 		    dblT		= temperature
    ! 		    dblD		= density
    ! 		    dblP		= pressure
    ! 		    intBase		= which base uncertainty? (1=cp,2=rho,3=w,4=p,5=cv)
    ! 	    	dblPropGen	= genuine target property value
    ! 		    intPhase	= which phase? relevant for Rho_TP
    !   	    dblPsat		= vapour pressure
    !
    !     OUTPUT: biggest, by coefficient manipulation calculated, target property
    !______________________________________________________________________________________________________________________




    implicit none

    type(type_gl) :: gl
    ABSTRACT INTERFACE
    DOUBLE PRECISION FUNCTION FKTPROPERTY_PROTOTYPE(gl,T,D, nrsubst)
    use module_all_types
    type(type_gl) :: gl
    DOUBLE PRECISION :: T,D
    integer:: nrsubst
    END FUNCTION FKTPROPERTY_PROTOTYPE
    END INTERFACE
    PROCEDURE(FKTPROPERTY_PROTOTYPE), POINTER, intent(in) :: FKTPROPERTY


    DOUBLE PRECISION :: dblXstart_Min, dblXstart_Max, dblDelta_allowed, dblXmax_allowed, dblXmin_allowed
    DOUBLE PRECISION :: dblNi, dblNGen, dblT, dblD, dblP, dblBound, dblPropGen, dblPsat, dblGen, dblUncty
    type(type_additional_parameters) :: parameters
    double precision :: vary_min, vary_max, variance
    integer :: coeff_nr
    INTEGER :: intIterations, intMaxIterations, errorflag, nr_ideal_coeff
    INTEGER	  :: i, intBase, intPhase

    double precision :: prop_plus, prop_minus, aux_slope, delta_prop, delta_ni, ni_factor, Dens

    !parameters = 0.D0

    !values for regula falsi
    parameters%a_p(2) = dblT
    parameters%a_p(3) = dblD
    parameters%a_p(4) = dblP
    parameters%a_p(5) = dblBound
    parameters%a_p(6) = intBase
    parameters%a_p(7) = intPhase
    parameters%a_p(8) = dblPsat

    !parameter(10) indicates, if ideal or real part coefficients are varied

    !the following part adds the coefficient of the ideal part of the heat capacity to the consideration, but it's not used by now
    !if ((intBase == 1) .or. (intBase  == 5)) then
    !    nr_ideal_coeff = ncp0poly(1)+ncp0pl(1)+ ncp0cosh(1)+ncp0sinh(1)   !the ideal coefficients are only relevant for the heat capacities
    !else
    nr_ideal_coeff = 0                                             !the ideal coefficients are not relevant for the other properties
    !end if

    coeff_nr = nr_ideal_coeff+gl%eos_coeff%nreg(1)+gl%eos_coeff%ncrt(1)

    !do for all regular and critical coefficients ni_u
    do i=1, (coeff_nr)

        if(nr_ideal_coeff >= i) then
            dblNGen = gl%CP0COEFF(i,1)		!preserves genuine value of ni(i)
        else
            dblNGen = gl%eos_coeff%ni(i-nr_ideal_coeff,1)		!preserves genuine value of ni(i)

            if (intBase == 1) then	        !- >  cp
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*1.0001D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_plus = CP_CALC(gl,dblT, Dens, 1)
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*0.9999D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_minus = CP_CALC(gl,dblT, Dens, 1)

            else if (intBase == 2) then	!- >  rho
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*1.0001D0
                prop_plus = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*0.9999D0
                prop_minus = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)

            else if (intBase == 3) then	!- >  ws
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*1.0001D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_plus = WS_CALC(gl,dblT, Dens, 1)
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*0.9999D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_minus = WS_CALC(gl,dblT, Dens, 1)

            else if (intBase == 4) then	!- >  p
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*1.0001D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_plus = P_CALC(gl,dblT, Dens, 1)
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*0.9999D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_minus = P_CALC(gl,dblT, Dens, 1)

            else if (intBase == 5) then	!- >  cv
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*1.0001D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_plus = CV_CALC(gl,dblT, Dens, 1)
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen*0.9999D0
                Dens = rhomix_calc(gl,dblT, dblP, dblD, 0, 0)
                prop_minus = CV_CALC(gl,dblT, Dens, 1)
            end if

            aux_slope = (prop_plus - prop_minus) / (0.0002D0*dblNGen)
            gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen
            delta_prop = dblBound - dblGen
            delta_ni = delta_prop / aux_slope

            ni_factor = 1.5D0

            variance = abs(delta_ni) / abs(dblNgen)

        end if

        if (variance < 0.3D0) then

            !variables for regula falsi
            if (delta_ni > 0.D0) then

                if(nr_ideal_coeff >= i) then
                    dblXstart_Min = vary_min * gl%CP0COEFF(i,1)
                    dblXstart_Max = vary_max * gl%CP0COEFF(i,1)
                    dblXmin_allowed = vary_min * gl%CP0COEFF(i,1)
                    dblXmax_allowed = vary_max * gl%CP0COEFF(i,1)
                    parameters%a_p(10) = 1.D0
                    parameters%a_p(1) = i
                else
                    dblXstart_Min = dblNGen
                    dblXstart_Max = dblNGen + (delta_ni * ni_factor)
                    dblXmin_allowed = dblNGen
                    dblXmax_allowed = dblNGen + (delta_ni * ni_factor)
                    parameters%a_p(10) = 2.D0
                    parameters%a_p(1) = i-nr_ideal_coeff
                end if

            else

                if(nr_ideal_coeff >= i) then
                    dblXstart_Min = vary_min * gl%CP0COEFF(i,1)
                    dblXstart_Max = vary_max * gl%CP0COEFF(i,1)
                    dblXmin_allowed = vary_min * gl%CP0COEFF(i,1)
                    dblXmax_allowed = vary_max * gl%CP0COEFF(i,1)
                    parameters%a_p(10) = 1.D0
                    parameters%a_p(1) = i
                else
                    dblXstart_Min = dblNGen + (delta_ni * ni_factor)
                    dblXstart_Max = dblNGen
                    dblXmin_allowed = dblNGen + (delta_ni * ni_factor)
                    dblXmax_allowed = dblNGen
                    parameters%a_p(10) = 2.D0
                    parameters%a_p(1) = i-nr_ideal_coeff
                end if

            end if


            !	!growing variances for the coefficients
            !
            !	do v = 1, 5
            !	!range of variance for regula falsi
            !	if (v == 1) then
            !	variance = 0.0001D0
            !	else if (v == 2) then
            !	variance = 0.001D0
            !	else if (v == 3) then
            !	variance = 0.01D0
            !	else if (v == 4) then
            !	variance = 0.1D0
            !	else if (v == 5) then
            !	variance = 0.5D0
            !	end if
            !
            !	vary_min = 1.D0 - variance
            !	vary_max = 1.D0 + variance
            !
            !
            !    if (dblNGen < 0.D0) then
            !        if(nr_ideal_coeff >= i) then
            !            dblXstart_Min = vary_max * CP0COEFF(i,1)
            ! 	        dblXstart_Max = vary_min * CP0COEFF(i,1)
            !	        dblXmin_allowed = vary_max * CP0COEFF(i,1)
            ! 	        dblXmax_allowed = vary_min * CP0COEFF(i,1)
            ! 	        parameters%a_p(10) = 1.D0
            ! 	        parameters%a_p(1) = i
            !        else
            !            dblXstart_Min = vary_max * ni(1,i-nr_ideal_coeff)
            ! 	        dblXstart_Max = vary_min * ni(1,i-nr_ideal_coeff)
            !	        dblXmin_allowed = vary_max * ni(1,i-nr_ideal_coeff)
            ! 	        dblXmax_allowed = vary_min * ni(1,i-nr_ideal_coeff)
            ! 	        parameters%a_p(10) = 2.D0
            ! 	        parameters%a_p(1) = i-nr_ideal_coeff
            ! 	    end if
            !    else
            !       !variables for regula falsi
            !        if(nr_ideal_coeff >= i) then
            !            dblXstart_Min = vary_min * CP0COEFF(i,1)
            ! 	        dblXstart_Max = vary_max * CP0COEFF(i,1)
            !	        dblXmin_allowed = vary_min * CP0COEFF(i,1)
            ! 	        dblXmax_allowed = vary_max * CP0COEFF(i,1)
            ! 	        parameters%a_p(10) = 1.D0
            ! 	        parameters%a_p(1) = i
            !        else
            !            dblXstart_Min = vary_min * ni(1,i-nr_ideal_coeff)
            ! 	        dblXstart_Max = vary_max * ni(1,i-nr_ideal_coeff)
            !	        dblXmin_allowed = vary_min * ni(1,i-nr_ideal_coeff)
            ! 	        dblXmax_allowed = vary_max * ni(1,i-nr_ideal_coeff)
            ! 	        parameters%a_p(10) = 2.D0
            ! 	        parameters%a_p(1) = i-nr_ideal_coeff
            ! 	    end if
            !	end if


            dblDelta_allowed = 1.D-7
            intMaxIterations = 30

            !do the manipulation for GBS- and NA- coefficients only in the immediate vacinity of the critical point
            if ((dblT < (gl%tc(1)*0.98D0) .or. (dblT > (gl%tc(1)*1.02D0))) .and. (i > (gl%eos_coeff%nreg(1)+nr_ideal_coeff))) EXIT
            if ((dblP < (gl%pc(1)*0.92D0) .or. (dblP > (gl%pc(1)*1.08D0))) .and. (i > (gl%eos_coeff%nreg(1)+nr_ideal_coeff))) EXIT

            call Regula_Falsi(gl,fktProp_diff, dblNi, dblXstart_Min, dblXstart_Max, dblDelta_allowed, dblXmin_allowed,&
                &	dblXmax_allowed, intMaxIterations, intIterations, errorflag, parameters)

            if(nr_ideal_coeff >= i) then
                gl%CP0COEFF(i,1) = dblNGen
            else
                gl%eos_coeff%ni(i-nr_ideal_coeff,1) = dblNGen
            end if

            if (errorflag /= 0) then
                errorflag = 0
                !do nothing, because no root found in reasonable vacinity of genuine ni_u or coefficient does not exist

            else
                ! save varied ni_u (but only if it's bigger than the previous one)
                if(nr_ideal_coeff >= i) then
                    if ((dblNi > dblNGen) .and. (gl%uncty%ni_varied_pos_ideal(i) == 0)) then
                        gl%uncty%ni_varied_pos_ideal(i) = dblNi
                        !used_coeff_variance_pos(i) = variance

                    else if ((dblNi > dblNGen) .and. (gl%uncty%ni_varied_pos_ideal(i) /= 0) .and. (dblNi < gl%uncty%ni_varied_pos_ideal(i))) then
                        gl%uncty%ni_varied_pos_ideal(i) = dblNi
                        !used_coeff_variance_pos(i) = variance

                    else if ((dblNi < dblNGen) .and. (gl%uncty%ni_varied_neg_ideal(i) == 0)) then
                        gl%uncty%ni_varied_neg_ideal(i) = dblNi
                        !used_coeff_variance_neg(i) = variance

                    else if ((dblNi < dblNGen) .and. (gl%uncty%ni_varied_neg_ideal(i) /= 0) .and. (dblNi > gl%uncty%ni_varied_neg_ideal(i))) then
                        gl%uncty%ni_varied_neg_ideal(i) = dblNi
                        !used_coeff_variance_neg(i) = variance

                    else
                        !do nothing
                    end if
                else
                    if ((dblNi > dblNGen) .and. (gl%uncty%ni_varied_pos(i-nr_ideal_coeff) == 0)) then
                        gl%uncty%ni_varied_pos(i-nr_ideal_coeff) = dblNi
                        !		            used_coeff_variance_pos(i-nr_ideal_coeff) = variance
                        !		            used_base_pos(i) = intBase

                    else if ((dblNi > dblNGen) .and. (gl%uncty%ni_varied_pos(i-nr_ideal_coeff) /= 0) .and. &
                        & (dblNi < gl%uncty%ni_varied_pos(i-nr_ideal_coeff))) then
                        gl%uncty%ni_varied_pos(i-nr_ideal_coeff) = dblNi
                        !		            used_coeff_variance_pos(i-nr_ideal_coeff) = variance
                        !		            used_base_pos(i) = intBase

                    else if ((dblNi < dblNGen) .and. (gl%uncty%ni_varied_neg(i-nr_ideal_coeff) == 0)) then
                        gl%uncty%ni_varied_neg(i-nr_ideal_coeff) = dblNi
                        !		            used_coeff_variance_neg(i-nr_ideal_coeff) = variance
                        !		            used_base_neg(i) = intBase

                    else if ((dblNi < dblNGen) .and. (gl%uncty%ni_varied_neg(i-nr_ideal_coeff) /= 0) .and. &
                        & (dblNi > gl%uncty%ni_varied_neg(i-nr_ideal_coeff))) then
                        gl%uncty%ni_varied_neg(i-nr_ideal_coeff) = dblNi
                        !		            used_coeff_variance_neg(i-nr_ideal_coeff) = variance
                        !		            used_base_neg(i) = intBase

                    else
                        !do nothing
                    end if
                end if
                !		    exit

            end if
        end if

        !end do

    end do

    !!find biggest delta                                                                !nun unnötig(?)
    !   j = maxloc(arDeltaPropMan, 1)
    !
    !!return manipulated target property which is responsible for the biggest delta
    !    if (arPropMan(j) == 0.D0)then
    !	    coeff_manipulation = dblPropGen
    !    else
    !	    coeff_manipulation = arPropMan(j)
    !    end if

    end subroutine coeff_manipulation


    !------------------------------------------------------------------------------------------------------------
    !                                                 ZEROFUNCTION
    !------------------------------------------------------------------------------------------------------------

    DOUBLE PRECISION FUNCTION fktProp_diff(gl,ni_u, parameters)
    !  function establishes target function for regula falsi, depending on whether cp, rho, w, p or cv has been
    !  chosen as base

    !		USE module_eos_coefficients
    !        USE module_fluid_parameters
    !        USE module_ideal_gas_coefficients




    implicit none

    type(type_gl) :: gl


    INTEGER :: j, intPhase
    DOUBLE PRECISION :: ni_u
    DOUBLE PRECISION :: dblBound, temp, dens, pres, dblD
    type(type_additional_parameters) :: parameters

    j = int(parameters%a_p (1))
    temp = parameters%a_p (2)
    dblD = parameters%a_p (3)
    pres = parameters%a_p(4)
    dblBound = parameters%a_p (5)
    intPhase = int(parameters%a_p(7))
    !parameters%a_p(8)- >  vapour pressure

    if(parameters%a_p(10) == 1.D0) then
        gl%CP0COEFF(j,1) = ni_u
    else
        gl%eos_coeff%ni(j,1) = ni_u
    End if

    !        Dens = rhomix_calc(gl,temp, pres, 0.D0, 0, 0)

    if (parameters%a_p(6) == 1) then					!1=cp
        Dens = rhomix_calc(gl,TEMP, PRES, dblD, 0, 1)
        fktProp_diff = dblBound - CP_CALC(gl,Temp,Dens, 1)

    else if (parameters%a_p(6) == 2) then				!2=rho
        fktProp_diff = dblBound - rhomix_calc(gl,TEMP, PRES, dblD, 0, 1)

    else if (parameters%a_p(6) == 3) then				!3=speed of sound
        Dens = rhomix_calc(gl,TEMP, PRES, dblD, 0, 1)
        fktProp_diff = dblBound - WS_CALC(gl,Temp,Dens, 1)

    else if (parameters%a_p(6) == 4) then				!4=p
        Dens = rhomix_calc(gl,TEMP, PRES, dblD, 0, 1)
        fktProp_diff = dblBound - P_CALC(gl,Temp,Dens, 1)

    else if (parameters%a_p(6) == 5) then				!5=cv
        Dens = rhomix_calc(gl,TEMP, PRES, dblD, 0, 1)
        fktProp_diff = dblBound - CV_CALC(gl,Temp,Dens, 1)
    end if


    END FUNCTION fktProp_diff



    end module uncty_module
