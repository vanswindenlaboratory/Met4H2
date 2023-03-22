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

    ! module for file fniderivs.f90
    submodule (fniderivs_module) impl
    !global use inclusion
    use module_all_types
    use calc_functions
    use rhomix_pt_module
    use flash_pure_module


    contains






    !****************************************************************************
    module SUBROUTINE FNIDERIVS(gl,T, D, GETDERIVS, SETDERIVS, nrsubst)
    !****************************************************************************
    ! FLORIAN DAUBER; DANMARK; 08.2009

    ! SUBROUTINE FOR THE CALCULATION OF ALL DERIVATIVES OF THE IDEAL PART
    ! OF THE HELMHOLTZ FREE ENERGY FOR PURE FLUIDS
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    ! GETDERIVS      - AN ARRAY WITH 6 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED IDEAL HELMHOLTZ ENERGY AS A FUNCTION OF D AND T
    !                2. 1ST DERIVATIVE WITH RESPECT TO D AT CONSTANT T
    !                3. 2ND DERIVATIVE WITH RESPECT TO D AT CONSTANT T
    !                4. 1ST MIXED DERIVATIVE WITH RESPECT TO D AND T
    !                5: 1ST DERIVATIVE WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                6: 2ND DERIVATIVE WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU**2
    !                7: 3RD DERIVATIVE WITH RESPECT TO TT AND D, MULTIPLIED BY TAU**2*D
    !                8: 3RD DERIVATIVE WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL**3
    !                9: 3RD DERIVATIVE WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU**3
    !               10: 3RD DERIVATIVE WITH RESPECT TO T AND DD, MULTIPLIED BY DEL**2 * TAU

    ! nrsubst - INTEGER THAT GIVES INFORMATION ON THE FLUID TO BE CONSIDERED
    !
    ! OUTPUT PARAMETERS:
    ! SETDERIVS      - AN ARRAY WITH 6 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !                  AS INDECATED IN "GETDERIVS"
    !-------------------------------------------------------------------------



    ! USE module_ideal_gas_coefficients !MODUL CALL
    ! USE module_fluid_parameters !MODUL CALL
    ! USE module_general_eos_parameters !MODUL CALL
    ! USE module_eos_coefficients !MODUL CALL



    implicit none

    type(type_gl) :: gl


    integer, dimension(nderivsi)::GETDERIVS
    DOUBLE PRECISION, dimension(nderivsi)::SETDERIVS
    integer::nrsubst

    double precision ::DELTA,T, D , FNI, FNID,FNIDD, FNIDT, FNITTAU, FNITTTAU2,TAU, FNIDTT, FNIDDD, FNITTT, FNITDD
    integer :: IPOL, KPE
    ! double precision :: RHOREF,TAUREF, DELREF
    integer :: K, I, J
    LOGICAL::C0KONST

    DOUBLE PRECISION, dimension(100) :: EXCH1, EXCH3, C, TPOTI, Mhyp
    DOUBLE PRECISION, dimension(100) :: EXCH2, EXCH4, THETA, M, THETAHYP



    !variables of SRK_FLUIDS
    !Thu February 2013
    double precision :: FNICV   !Ideal Helmholtz Energy by integration of CV0/R=cp0/R-1 two times.
    !Equation of CP0/R published by Michael Kleiber and Ralph Joh
    !and is called PPDS-Equation. In this equation is T=TR/TAU.
    double precision :: FNICVD, FNICVDD, FNICVDT, FNICVT, FNICVTT
    double precision::  EPSILON,C1_CV0,C2_CV0,TAU_CV0,DELTA_CV0,TR_SRK,DR_SRK

    !Term calculated before, last number equal to the exponent
    double precision::IOTA,IOTA2,IOTA3
    double precision:: ZETA
    double precision:: TR2, TR3
    double precision:: OMEGA_lc, OMEGA2, OMEGA3, OMEGA4, OMEGA5

    !variables of cp0-model by Joback
    !Thu 05/2014
    double precision:: BETA_1,BETA_2,BETA_3,BETA_4 !constant parameters e.g. BETA_1 =Delta_A-37.93
    double precision:: TAU_CV02, TAU_CV03,TAU_CV04 !Tau=Tc/T

    !Thu
    double precision:: EXCH2_2, T_1, T_2, TEXCH2, TEXCH2_2, TEXCH2_3, Theta2, Theta3, tau3

    double precision:: RMIX


    SETDERIVS = 0.d0 !ALL VALUES ARE SET TO 0 AND CAN BE OVERWRITEN

    DELTA = D / gl%RHORED(nrsubst) !REDUCED DENSITY

    !Monika, Februar 2015
    if (gl%tcp0red(nrsubst) == 1.d0) then
        TAU = gl%tc(nrsubst)/T
    else
        TAU = gl%tcp0red(nrsubst)/T
    end if

    !TAU = tredmix/T   !REDUCED TEMPERATURE FOR THE IDEAL GAS PART USING THE REDUCING TEMP OF THE RESIDUAL PART

    IPOL = gl%ncp0poly(nrsubst) !NUMBER OF CP0 COEFFICIENTS Ci, Mk ; FIRST VARIABLE IS C0
    KPE = gl%ncp0pl(nrsubst) !NUMBER OF EXPONENTS ti, THETAk

    EXCH1 = 0.D0
    EXCH2 = 0.D0
    EXCH3 = 0.D0
    EXCH4 = 0.D0
    C  = 0.D0
    TPOTI = 0.D0
    THETA = 0.D0
    M = 0.D0
    THETAHYP = 0.D0
    Mhyp = 0.D0
    TAU3 = TAU*TAU*TAU

    !Published fluid specific ideal gas equation or general treatment
    !Thu February 2013
    If ((gl%Eq_type(nrsubst) == 1) .or. (gl%Eq_type(nrsubst) == 7) .or. (gl%components(1) == "oil")) then

        !*********************************************************
        !      CONVERSION OF THE VARIABLES OF THE CP0-TERM FOR USING
        !      THEM IN THE DERIVATIVES OF THE HELMHOLTZ ENERGY OF THE IDEAL GAS
        !**********************************************************

        if (gl%cpmodel(nrsubst)) then
            J=0
            DO I = 1, IPOL
                J=J+1
                IF (gl%cp0exp(J,nrsubst) == 0.d0)THEN
                    C(I) = gl%CP0COEFF(J,nrsubst) - gl%REQ(nrsubst) / gl%cp0red(nrsubst)
                    if (gl%cp0red(nrsubst) == 1.d0) C(I) = gl%CP0COEFF(J,nrsubst) - 1.d0
                ELSE IF (gl%cp0exp(J,nrsubst) /= -1.D0) THEN ! for -1 the cp0-coeffs and exponents are used directly in FNI and its derivatives
                    C(I) = -gl%CP0COEFF(J,nrsubst)/(gl%cp0exp(J,nrsubst)*(gl%cp0exp(J,nrsubst)+1.d0)) &
                        & *gl%tred(nrsubst)**gl%cp0exp(J,nrsubst) / gl%tcp0red(nrsubst)**gl%cp0exp(J,nrsubst)
                    ! & *tcp0red(nrsubst)**cp0exp(J,nrsubst)
                END IF
            END DO

            DO K = 1, KPE
                J=J+1
                M(J) = gl%CP0COEFF(J,nrsubst) !FOR THE PLANCK-EINSTEIN TERMS
            END DO

            DO K = 1, gl%ncp0cosh(nrsubst) + gl%ncp0sinh(nrsubst)
                J=J+1
                Mhyp(J) = gl%CP0COEFF(J,nrsubst)/gl%cp0hyp1(J,nrsubst)**2 !FOR THE hyperbolic TERMS
            END DO

            J=0
            DO I = 1, IPOL
                J=J+1
                TPOTI(J) = -gl%cp0exp(J,nrsubst) !FIRST VARIABLE IS C0
            END DO

            DO K = 1, KPE
                J=J+1
                THETA(J) = gl%cp0exp(J,nrsubst)/gl%tred(nrsubst) !tredmix !Andreas August 2014  !FOR THE PLANCK-EINSTEIN TERMS
            END DO

            DO K = 1, gl%ncp0cosh(nrsubst) + gl%ncp0sinh(nrsubst)
                J=J+1
                THETAHYP(J) = gl%cp0hyp1(J,nrsubst)/gl%tred(nrsubst) !tredmix !Andreas August 2014/tredmix !FOR THE hyperbolic expressions
            END DO
        endif



        !*********************************************************
        !      Terms THAT ARE NEEDED IN SOME DERIVATIVES AND
        !      CAN BE CALCULATED IN ADVANCE TO SAVE TIME:
        !**********************************************************
        if (gl%cpmodel(nrsubst)) then
            J=0
            DO I = 1, IPOL
                J=J+1
                EXCH1(I) = TAU**TPOTI(I)
            END DO

            DO I = 1, KPE
                J=J+1
                EXCH2(J) = EXP(-THETA(J) * TAU) !FOR THE PLANCK-EINSTEIN TERM
            END DO
            !      else !if (cppcheck(nrsubst) == 'PHK') then
            !        J = ncp0log(nrsubst) + ncp0(nrsubst) + ncp0logexp(nrsubst) !count J to right positition for cosh and sinh
            !      endif

            IF(gl%ncp0cosh(nrsubst) /= 0) THEN
                DO I = 1,  gl%ncp0cosh(nrsubst)
                    J=J+1
                    EXCH3(J)=cosh(THETAHYP(J)*TAU)!FOR THE hyperbolic expressions
                END DO
            END IF

            IF(gl%ncp0SINh(nrsubst) /= 0) THEN
                DO I = 1, gl%ncp0sinh(nrsubst)
                    J=J+1
                    EXCH4(J)=sinh(THETAHYP(J)*TAU)!FOR THE hyperbolic expressions
                END DO
            END IF

        elseif (gl%phkmodel(nrsubst)) then
            J = 0
            J = gl%ncp0log(nrsubst) + gl%ncp0(nrsubst)
            IF(gl%ncp0cosh(nrsubst) /= 0) THEN
                DO I = 1,  gl%ncp0cosh(nrsubst)
                    J=J+1
                    EXCH3(J)=cosh(gl%cp0exp(J,nrsubst)*TAU)!FOR THE cos hyperbolic expressions
                END DO
            END IF

            IF(gl%ncp0SINh(nrsubst) /= 0) THEN
                DO I = 1, gl%ncp0sinh(nrsubst)
                    J=J+1
                    EXCH4(J)=sinh(gl%cp0exp(J,nrsubst)*TAU)!FOR THE sin hyperbolic expressions
                END DO
            END IF

            !ELSEIF (phmodel(nrsubst)) THEN
            !      J = IPOL
            !      DO I = 1, KPE
            !      J=J+1
            !       THETA(J) = -cp0exp(J,nrsubst) !FOR THE PLANCK-EINSTEIN EXPONENTS
            !      END DO

        End if




        !***************************************************
        !FNI       IDEAL GAS PART OF THE REDUCED
        !          HELMHOLTZ ENERGY [-]
        !***************************************************

        IF (GETDERIVS(1) == 1) THEN
            FNI = 0.D0
            J=0
            C0KONST = .FALSE.

            !**********************************************************************************************
            if(gl%cpmodel(nrsubst)) then
                !**********************************************************************************************

                !Correction for the different gas constants of components in a mixture model
                !The ideal part of the involved components might differ from each other. In this case
                ! the following correction of the gas constants applies R_ideal_part(component i)/R_mix
                !The correction of all derivatives of the ideal Helmholtz energy is
                ! done at the end of this subroutine in a general way.
                !Attention: the LN(delta) is not affected by the reducing gas constant of the ideal part
                ! but is a calculation of thermal properties and, therefore, affected by possibly different
                ! gas constants of the residual(!) part. Therfore, the correction with the reducing
                ! gas constant of the ideal part is reverted here:
                !  "/(gl%cp0red(nrsubst)/Rmix)"
                !and a correction with the gas constant of the residual part is done instead:
                ! "*gl%Req(nrsubst)/Rmix"

                if ((.not.gl%ref) .and. (gl%ncomp .gt. 1)) then
                    !The correction of the gas constants is only done when
                    ! not calculating the reference state

                    !FNI = DLOG(DELTA) / (gl%cp0red(nrsubst)/Rmix) * gl%Req(nrsubst)/Rmix
                    !in the original correction, the mixture gas constant cancels out:
                    FNI = DLOG(DELTA) * gl%Req(nrsubst) / gl%cp0red(nrsubst) 
                else
                    !This part applies when calculating the reference state
                    FNI = DLOG(DELTA)
                end if

                FNI = FNI + gl%C2(nrsubst) + gl%C1(nrsubst)*TAU
                DO  I = 1,IPOL  !FOR THE POLYNOMIAL TERM
                    J = J+1
                    IF (dabs(gl%cp0exp(J,nrsubst)) < 1.D-6) THEN
                        FNI =FNI+ C(I)*DLOG(TAU)
                        C0KONST = .TRUE.
                    ELSEIF (gl%cp0exp(J,nrsubst) == -1.D0) THEN
                        !Andreas August 2014
                        !FNI = FNI - cp0coeff(nrsubst, J) / tredmix * (TAU * dlog(TAU) - TAU)
                        FNI = FNI - gl%cp0coeff(J,nrsubst) / gl%tred(nrsubst) * (TAU * dlog(TAU) - TAU)
                    ELSE
                        FNI = FNI + C(J) * EXCH1(J)
                    END IF
                    if((.not.c0konst).AND.(i == ipol)) fni = fni - Dlog(tau)   !manual transformation from cv0 to cp0
                END DO

                IF(KPE  /= 0) THEN   !FOR THE PLANCK-EINSTEIN TERM
                    DO I = 1, KPE
                        J = J+1
                        FNI = FNI + M(J) * DLOG(1.d0 - EXCH2(J))
                    END DO
                END IF

                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNI = FNI - Mhyp(J) * DLOG(EXCH3(J))
                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINH(nrsubst)
                        J = J+1
                        FNI = FNI + Mhyp(J) * DLOG(EXCH4(J))
                    END DO
                END IF


                !**********************************************************************************************
            ELSEIF(gl%phkmodel(nrsubst)) then
                !**********************************************************************************************

                J = gl%ncp0log(nrsubst) + gl%ncp0(nrsubst)


                !Correction for the different gas constants of the GERG-2008
                !The ideal part of the GERG-2008 was developed with R = 3.31451
                ! but is used with R = 8.314472. The following correction applies
                ! for this difference.
                !Attention: the LN(delta) is not affected by this correction, because
                ! it is calculated ONLY with the residual part of the equation of state,
                ! which means with R = 8.314472. Therefore, no correction of this part is needed.
                !The correction of all derivatives of the ideal Helmholtz energy is
                ! done at the end of this subroutine in a general way. Therefore,
                ! this correction is reverted here for the LN(delta).

                if (.not.gl%ref) then
                    !The correction of the gas constants is only done when
                    ! not calculating the reference state
                    call R_mix_calc(gl,Rmix)
                    FNI = DLOG(DELTA) / (gl%cp0red(nrsubst)/Rmix)
                else
                    !This part applies when calculating the reference state
                    FNI = DLOG(DELTA)
                end if

                FNI = FNI + gl%C2(nrsubst) + gl%C1(nrsubst) * TAU + gl%CP0COEFF(1,nrsubst) * DLOG(TAU)

                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1

                        FNI = FNI + gl%CP0COEFF(J,nrsubst) * DLOG(EXCH3(J))

                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINH(nrsubst)
                        J = J+1

                        FNI = FNI + gl%CP0COEFF(J,nrsubst) * DLOG(EXCH4(J))

                    END DO
                END IF

                !**********************************************************************************************
            ELSEIF(gl%phmodel(nrsubst)) THEN
                !**********************************************************************************************

                FNI = DLOG (DELTA) + gl%cp0coeff(3,nrsubst) + gl%cp0coeff(2,nrsubst)*TAU + gl%cp0coeff(1,nrsubst)*DLOG(TAU)
                J=3

                IF(IPOL  > 3) THEN  !FOR THE POLYNOMIAL TERMS
                    DO  I = 4,IPOL
                        J=J+1
                        FNI = FNI + gl%CP0COEFF(I,nrsubst) * TAU**gl%cp0exp(I,nrsubst)
                    END DO
                END IF

                IF(KPE  /= 0) THEN   !FOR THE PLANCK-EINSTEIN TERMS
                    DO I = 1, KPE
                        J = J+1
                        FNI = FNI + gl%CP0COEFF(J,nrsubst) * DLOG(1.D0 - DEXP(gl%cp0exp(J,nrsubst) * TAU)) ! cp0exp(J,nrsubst)=-THETA
                    END DO
                END IF

            END IF
            SETDERIVS(1) = FNI
        END IF

        !***************************************************
        !FNID      FIRST DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY [-]
        !****************************************************
        IF (GETDERIVS(2) == 1) THEN
            SETDERIVS(2) = 1.d0
        END IF

        !***************************************************
        !FNIDD     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY [-]
        !****************************************************
        IF (GETDERIVS(3) == 1) THEN
            FNIDD = 0.D0
            FNIDD = -1.d0 / DELTA**2
            SETDERIVS(3) = FNIDD*DELTA**2
        END IF


        !***************************************************
        !FNIDT     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY AND THE REDUCED TEMPERATURE[-]
        !****************************************************
        IF (GETDERIVS(4) == 1) THEN
            SETDERIVS(4) = 0.d0
        END IF

        !***************************************************
        !FNITTAU     FIRST DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !         WITH RESPECT TO THE REDUCED TEMPERATURE [-]
        !         MULTIPLIED WITH TAU
        !****************************************************
        IF (GETDERIVS(5) == 1) THEN
            FNITTAU = 0.D0
            J=0
            C0KONST = .FALSE.

            !**********************************************************************************************
            IF(gl%cpmodel(nrsubst)) THEN
                !**********************************************************************************************

                FNITTAU = gl%C1(nrsubst)

                DO  I = 1, IPOL
                    J = J+1
                    IF (gl%cp0exp(J,nrsubst) == 0.D0) THEN
                        FNITTAU = FNITTAU+ C(I)/TAU
                        C0KONST = .TRUE.
                    ELSE IF (gl%cp0exp(J,nrsubst) == -1.D0) THEN
                        !Andreas August 2014
                        !FNITTAU = FNITTAU - cp0coeff(nrsubst, J) / tredmix * dlog(TAU)
                        FNITTAU = FNITTAU - gl%cp0coeff(J,nrsubst) / gl%tred(nrsubst) * dlog(TAU)
                    ELSE
                        FNITTAU = FNITTAU + C(J) * TPOTI(J) * EXCH1(J)/TAU
                    END IF
                    IF ((.not.c0konst).AND.(i == ipol)) FNITTAU = FNITTAU - 1.D0/TAU   !manual correction if cp0 constant is missing
                END DO

                IF(KPE  /= 0) THEN
                    DO I = 1, KPE
                        J = J+1
                        FNITTAU = FNITTAU + M(J) * THETA(J) * (1.D0/(1.D0-EXCH2(J))-1.D0)  !FOR THE PLANCK-EINSTEIN TERM
                    END DO
                END IF


                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTAU = FNITTAU - MHYP(J) * THETAHYP(J)* TANH(THETAHYP(J)*TAU)!TANH=SINH/COSH
                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTAU = FNITTAU + MHYP(J) * THETAHYP(J) / TANH(THETAHYP(J)*TAU)
                    END DO
                END IF

                !**********************************************************************************************
            ELSEIF(gl%phkmodel(nrsubst)) THEN
                !**********************************************************************************************
                J = gl%NCP0LOG(nrsubst) + gl%NCP0(nrsubst)
                !    FNITTAU = CP0COEFF(nrsubst,3) + CP0COEFF(nrsubst,1)/TAUIDEAL

                FNITTAU = gl%C1(nrsubst) + gl%CP0COEFF(1,nrsubst)/TAU

                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTAU = FNITTAU + gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)* TANH(gl%cp0exp(J,nrsubst)*TAU)
                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTAU = FNITTAU + gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst) / TANH(gl%cp0exp(J,nrsubst)*TAU)
                    END DO
                END IF

                !**********************************************************************************************
            ELSEIF(gl%phmodel(nrsubst)) THEN
                !**********************************************************************************************
                FNITTAU = gl%CP0COEFF(2,nrsubst) + gl%CP0COEFF(1,nrsubst)/TAU
                J=3

                IF(IPOL  > 3) THEN  !FOR THE POLYNOMIAL TERM
                    DO  I = 4,IPOL
                        J=J+1
                        FNITTAU = FNITTAU + gl%cp0coeff(I,nrsubst) * gl%cp0EXP(I,nrsubst) * tau**(gl%cp0exp(I,nrsubst) - 1.D0)
                    END DO
                END IF

                IF(KPE  /= 0) THEN   !FOR THE PLANCK-EINSTEIN TERM
                    DO I = 1, KPE
                        J = J+1
                        FNITTAU = FNITTAU + gl%CP0COEFF(J,nrsubst) * (-gl%cp0exp(J,nrsubst)) * (1.D0/(1.d0 - dexp(gl%cp0exp(J,nrsubst) * TAU))-1.D0)  !cp0exp(J,nrsubst)=-THETA
                    END DO
                END IF

            END IF

            SETDERIVS(5) = FNITTAU*TAU
        END IF

        !***************************************************
        !FNITTTAU2     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED TEMPERATURE [-]
        !          MULTIPLIED WITH TAU**2
        !****************************************************
        IF (GETDERIVS(6) == 1) THEN
            FNITTTAU2 = 0.D0
            J=0
            C0KONST = .FALSE.

            !**********************************************************************************************
            IF (gl%cpmodel(nrsubst)) THEN
                !**********************************************************************************************

                IF(IPOL /= 0)THEN
                    DO  I = 1, IPOL
                        J = J+1
                        IF (gl%cp0exp(J,nrsubst) == 0.d0) THEN
                            FNITTTAU2 = FNITTTAU2 - C(I)/TAU**2
                            C0KONST = .TRUE.
                        ELSEIF (gl%cp0exp(J,nrsubst) == -1.D0) THEN
                            !Andreas, August 2014
                            !FNITTTAU2 = FNITTTAU2 - cp0coeff(nrsubst, J) / tredmix / dlog(TAU)
                            FNITTTAU2 = FNITTTAU2 - gl%cp0coeff(J,nrsubst) / gl%tred(nrsubst) / TAU
                        ELSE
                            FNITTTAU2 = FNITTTAU2+ C(J) * TPOTI(J) * (TPOTI(J)-1.d0)* EXCH1(J)/TAU**2
                        END IF
                        IF ((.not.c0konst).AND.(i == ipol)) FNITTTAU2 = FNITTTAU2 + 1.D0/TAU**2   !manual correction if cp0 constant is missing
                    END DO
                END IF

                IF(int(KPE)  /= 0) THEN
                    DO I = 1, KPE
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 - M(J) * THETA(J)**2 * (EXCH2(J)/(1.D0-EXCH2(J))**2)  !FOR THE PLANCK-EINSTEIN TERM
                    END DO
                END IF
                !      else !if (cppcheck(nrsubst) == 'PHK') then
                !        J = ncp0log(nrsubst) + ncp0(nrsubst) + ncp0logexp(nrsubst) !count J to right positition for cosh and sinh
                !      endif
                !
                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 - Mhyp(J) * THETAHYP(J)**2 / EXCH3(J)**2
                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 - Mhyp(J) * THETAHYP(J)**2 / EXCH4(J)**2
                    END DO
                END IF

                !**********************************************************************************************
            ELSEIF(gl%phkmodel(nrsubst)) THEN
                !**********************************************************************************************
                J = gl%NCP0LOG(nrsubst) + gl%NCP0(nrsubst)
                FNITTTAU2 = - gl%CP0COEFF(1,nrsubst)/TAU**2

                IF(gl%ncp0cosh(nrsubst) /= 0) THEN       !FOR THE COShyperbolic expressions
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 + gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)**2 / EXCH3(J)**2
                    END DO
                END IF

                IF(gl%ncp0SINh(nrsubst) /= 0) THEN       !FOR THE SINhyperbolic expressions
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 - gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)**2 / EXCH4(J)**2
                    END DO
                END IF
                
                !**********************************************************************************************
            ELSEIF(gl%phmodel(nrsubst)) THEN
                !**********************************************************************************************
                FNITTTAU2 = - gl%CP0COEFF(1,nrsubst)/TAU**2
                J=3

                IF(IPOL  > 3) THEN  !FOR THE POLYNOMIAL TERM
                    DO  I = 4,IPOL
                        J=J+1
                        FNITTTAU2 = FNITTTAU2 + gl%cp0coeff(I,nrsubst) * gl%cp0EXP(I,nrsubst) * &
                            (gl%cp0EXP(I,nrsubst) - 1.D0) * tau**(gl%cp0exp(I,nrsubst) - 2.D0)
                    END DO
                END IF

                IF(KPE  /= 0) THEN   !FOR THE PLANCK-EINSTEIN TERM
                    DO I = 1, KPE
                        J = J+1
                        FNITTTAU2 = FNITTTAU2 - gl%CP0COEFF(J,nrsubst) * (-gl%cp0exp(J,nrsubst))**2 * &
                            dexp(gl%cp0exp(J,nrsubst) * TAU)/(1.d0 - dexp(gl%cp0exp(J,nrsubst) * TAU))**2  !cp0exp(J,nrsubst)=-THETA
                    END DO
                END IF

            END IF


            SETDERIVS(6) = FNITTTAU2*TAU**2

        END IF





        !***************************************************
        !FNIDTT    THIRD DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY AND TEMPERATURE [-]
        !          MULTIPLIED WITH DELTA*TAU**2
        !****************************************************
        IF (GETDERIVS(7)==1) THEN
            FNIDTT = 0.d0
            SETDERIVS(7)=FNIDTT*DELTA*TAU*TAU
        END IF

        !***************************************************
        !FNIDDD    THIRD DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY [-]
        !          MULTIPLIED WITH DELTA**3
        !****************************************************
        IF (GETDERIVS(8)==1) THEN
            !FNIDDD = 2.d0/DELTA3
            !SETDERIVS(8)=FNIDDD*DELTA3
            SETDERIVS(8)=2.d0
        END IF

        !***************************************************
        !FNITTT    THIRD DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED TEMPERATURE [-]
        !          MULTIPLIED WITH TAU**3
        !****************************************************
        IF (GETDERIVS(9)==1) THEN
            FNITTT=0.d0
            J=0
            C0KONST = .FALSE.

            !**********************************************************************************************
            IF(gl%cpmodel(nrsubst)) THEN
                !**********************************************************************************************
                !2014/10/02
                IF(IPOL/=0)THEN
                    DO  I = 1, IPOL
                        J = J+1
                        IF (gl%cp0exp(J,nrsubst) == 0.d0) THEN
                            FNITTT = FNITTT +2.d0*C(I)/TAU3
                            C0KONST = .TRUE.
                        ELSEIF (gl%cp0exp(J,nrsubst) == -1.D0) THEN !Ausnahme
                            T_1=1.d0/TAU
                            T_2=T_1*T_1
                            FNITTT = FNITTT + gl%cp0coeff(J,nrsubst) / gl%tred(nrsubst) * T_2
                        ELSE !Hauptteil
                            FNITTT = FNITTT+ C(J) * TPOTI(J) * (TPOTI(J)-1.d0)*(TPOTI(J)-2.d0)* EXCH1(J)/TAU3
                        END IF
                    END DO
                    IF ((.not.c0konst)) FNITTT = FNITTT - 2.D0/TAU3   !manual correction if cp0 constant is missing
                END IF

                !2014/10/02
                IF(int(KPE) /=0) THEN
                    DO I = 1, KPE
                        J = J+1
                        THETA2=THETA(J)*THETA(J)
                        THETA3=THETA2*THETA(J)
                        TEXCH2=1.d0-EXCH2(J)
                        TEXCH2_2=TEXCH2*TEXCH2
                        TEXCH2_3=TEXCH2_2*TEXCH2
                        EXCH2_2=EXCH2(J)*EXCH2(J)

                        FNITTT = FNITTT+ M(J) * THETA3 * ((EXCH2(J)/TEXCH2_2)+(2.d0*EXCH2_2/TEXCH2_3))  !FOR THE PLANCK-EINSTEIN TERM
                    END DO
                END IF

                !FOR THE hyperbolic expressions
                !*******************************
                IF(gl%ncp0cosh(nrsubst)/=0) THEN
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTT = FNITTT + 2.d0 * Mhyp(J) * THETAHYP(J)**3*TANH(THETAHYP(J)*TAU)/EXCH3(J)**2
                    END DO
                END IF
                IF(gl%ncp0SINh(nrsubst)/=0) THEN
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTT = FNITTT + 2.d0 * Mhyp(J) * THETAHYP(J)**3/(TANH(THETAHYP(J)*TAU)*EXCH4(J)**2)
                    END DO
                END IF
                !*******************************
                !2014/10/03
                !**********************************************************************************************
            ELSEIF(gl%phkmodel(nrsubst)) THEN
                !**********************************************************************************************
                J = gl%NCP0LOG(nrsubst) + gl%NCP0(nrsubst)
                FNITTT = 2.d0* gl%CP0COEFF(1,nrsubst)/TAU3

                !FOR THE hyperbolic expressions
                !*******************************
                IF(gl%ncp0cosh(nrsubst)/=0) THEN
                    DO I = 1, gl%ncp0cosh(nrsubst)
                        J = J+1
                        FNITTT = FNITTT -2.D0* gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)**3*TANH(gl%cp0exp(J,nrsubst)*tau) / EXCH3(J)**2
                    END DO
                END IF
                IF(gl%ncp0SINh(nrsubst)/=0) THEN
                    DO I = 1, gl%ncp0SINh(nrsubst)
                        J = J+1
                        FNITTT = FNITTT +2.D0* gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)**3 / (TANH(gl%cp0exp(J,nrsubst)*tau)*EXCH4(J)**2)
                    END DO
                END IF

                !**********************************************************************************************
            ELSEIF(gl%phmodel(nrsubst)) THEN
                !**********************************************************************************************
                FNITTT = 2.D0* gl%CP0COEFF(1,nrsubst)/TAU3
                J=3

                IF(IPOL >3) THEN  !FOR THE POLYNOMIAL TERM
                    DO  I = 4,IPOL
                        J=J+1
                        FNITTT = FNITTT + gl%cp0coeff(I,nrsubst) * gl%cp0EXP(I,nrsubst) * &
                            (gl%cp0EXP(I,nrsubst) - 1.D0) *(gl%cp0EXP(I,nrsubst) - 2.D0) * tau**(gl%cp0exp(I,nrsubst) - 3.D0)
                    END DO
                END IF

                IF(KPE /=0) THEN   !FOR THE PLANCK-EINSTEIN TERM
                    DO I = 1, KPE
                        J = J+1
                        FNITTT = FNITTT - gl%CP0COEFF(J,nrsubst) * gl%cp0exp(J,nrsubst)**3 * &
                            (dexp(gl%cp0exp(J,nrsubst) * tau)/(1.d0 - dexp(gl%cp0exp(J,nrsubst) * tau))**2+&  !cp0exp(J,nrsubst)=-THETA
                            (2.D0*dexp(gl%cp0exp(J,nrsubst) * tau)*dexp(gl%cp0exp(J,nrsubst) * tau))/(1.d0 - dexp(gl%cp0exp(J,nrsubst) * tau))**3)
                    END DO
                END IF

            END IF
            SETDERIVS(9)=FNITTT*TAU*TAU*TAU

        END IF

        !***************************************************
        !FNITDD    THIRD DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED TEMPERATURE AND DENSITY [-]
        !          MULTIPLIED WITH TAU*DELTA*DELTA
        !****************************************************
        IF (GETDERIVS(10)==1) THEN
            SETDERIVS(10)=0.d0
        END IF
        !**************************************************************

        !This if-case is important. The correction of the ratio of the gas constants
        ! must be considered for mixtures, but
        ! NOT for calculating the integration constants of the corresponding pure fluids
        if ((.not.gl%ref) .and. (gl%ncomp .gt. 1)) then
            call R_mix_calc(gl,Rmix)
            SETDERIVS = SETDERIVS * gl%cp0red(nrsubst)/Rmix
        end if

        continue




        !**********************************************************************************************
        !Equation of CP0/R published by Michael Kleiber and Ralph Joh (PPDS-Equation)
    else
        !**********************************************************************************************
        !constants of the helmoltz energy
        C1_CV0=gl%C1(nrsubst)
        C2_CV0=gl%C2(nrsubst)

        !main variables
        TR_SRK=gl%tc(nrsubst)
        TAU_CV0=TR_SRK/T
        DR_SRK=gl%rhoc(nrsubst)
        DELTA_CV0= D/DR_SRK
        EPSILON=gl%C_CV0(nrsubst)-gl%B_CV0(nrsubst)
        OMEGA_lc=gl%A_CV0(nrsubst)*TAU_CV0+TR_SRK


        !!terms that can be calculated in advance to save time: 25.03.2013
        OMEGA2=OMEGA_lc*OMEGA_lc
        OMEGA3=OMEGA2*OMEGA_lc
        OMEGA4=OMEGA3*OMEGA_lc
        OMEGA5=OMEGA4*OMEGA_lc

        TR2=TR_SRK*TR_SRK
        TR3=TR2*TR_SRK

        ZETA = gl%A_CV0(nrsubst)*TAU_CV0/(gl%A_CV0(nrsubst)*TAU_CV0+TR_SRK)

        IOTA =TR_SRK/(TAU_CV0*gl%A_CV0(nrsubst)+TR_SRK)
        IOTA2=IOTA*IOTA
        IOTA3=IOTA2*IOTA


        !***************************************************
        !FNI       IDEAL GAS PART OF THE REDUCED
        !          HELMHOLTZ ENERGY [-]
        !***************************************************

        IF (GETDERIVS(1) == 1) THEN
            !21.05.2013
            if (abs(gl%A_CV0(nrsubst)) .lt. 1.d-12) then    !noble gases
                FNICV = (1.d0 - gl%C_CV0(nrsubst)) * dlog(TAU_CV0) - C2_CV0 - C1_CV0 * TAU_CV0 - dlog(DELTA_CV0)
            else                                            !all other fluids
                !FNICV=-(gl%B_CV0(nrsubst)-1.D0)*dlog(TAU_CV0)+EPSILON*TR2*&
                !        (gl%A_CV0(nrsubst)*TAU_CV0*(gl%D_CV0(nrsubst)+gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+gl%E_CV0(nrsubst)+2.d0)*&
                !        (dlog(OMEGA_lc)-dlog(TAU_CV0))/TR3-&
                !        (gl%D_CV0(nrsubst)+gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+gl%E_CV0(nrsubst))*dlog(OMEGA_lc)/TR2+&
                !        (gl%D_CV0(nrsubst)+gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+gl%E_CV0(nrsubst)+2.d0)*dlog(OMEGA_lc)/TR2+&
                !        (gl%D_CV0(nrsubst)+gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+gl%E_CV0(nrsubst))/(2.d0*TR_SRK*OMEGA_lc)+&
                !        (gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst))*TR_SRK/(12.d0*OMEGA3)+&
                !        (gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+gl%E_CV0(nrsubst))/(6.d0*OMEGA2)+&
                !        gl%G_CV0(nrsubst)*TR2/(20.d0*OMEGA4)-&
                !        dlog(TAU_CV0)/TR2)&
                !        -C2_CV0-C1_CV0*TAU_CV0-dlog(DELTA_CV0)
                !        !C2_CV0*TAU_CV0+C1_CV0  changed by Theresa

                ! Numerical derivative of alpha_0 above did not yield dalpha0_dtau. alpha0 (and therefore the entropy) seem to be wrong.
                ! The integration of dalpha0_dtau was carried out again, see below
                ! Erik Mickoleit, Andreas Jäger, October 2019
                FNICV=-(gl%B_CV0(nrsubst)-1.D0) * dlog(TAU_CV0) + gl%A_CV0(nrsubst) * EPSILON * TR2 * ((- (TAU_CV0 * dlog(TAU_CV0) - TAU_CV0) + (TAU_CV0 + TR_SRK / gl%A_CV0(nrsubst)) * dlog(gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) - TAU_CV0) &
                    * (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst) + 2.d0) / TR3    &
                    - (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst) + 1.d0) * dlog(gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) / (gl%A_CV0(nrsubst) * TR2)   &
                    - (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst)) / (2.d0 * TR_SRK) * 1.D0 / (-1.D0) * (gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) ** (-1.D0) / gl%A_CV0(nrsubst)  &
                    - (gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst)) / 3.d0 * 1.D0 / (-2.D0) * (gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) ** (-2.D0) / gl%A_CV0(nrsubst)   &
                    - TR_SRK * (gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst)) / 4.d0 * 1.D0 / (-3.D0) * (gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) ** (-3.D0) / gl%A_CV0(nrsubst)    &
                    - gl%G_CV0(nrsubst) * TR2 / 5.d0 * 1.D0 / (-4.D0) * (gl%A_CV0(nrsubst) * TAU_CV0 + TR_SRK) ** (-4.D0) / gl%A_CV0(nrsubst) - 1.d0 / (gl%A_CV0(nrsubst) * TR2) * dlog(TAU_CV0)) - C2_CV0 - C1_CV0 * TAU_CV0 - dlog(DELTA_CV0)
            end if

            SETDERIVS(1) = -FNICV
        end if




        !***************************************************
        !FNID      FIRST DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY [-]
        !****************************************************
        if (GETDERIVS(2) == 1) then
            SETDERIVS(2) = 1.d0     !Monika, March 2017
        end if


        !***************************************************
        !FNIDD     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY [-]
        !****************************************************
        if (GETDERIVS(3) == 1) then
            !FNICVDD = -1.d0/(DELTA_CV0*DELTA_CV0)
            SETDERIVS(3) = -1.d0        !Monika, March 2017
        end if


        !***************************************************
        !FNIDT     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED DENSITY AND THE REDUCED TEMPERATURE[-]
        !****************************************************
        if (GETDERIVS(4) == 1) then
            FNICVDT=0.d0
            SETDERIVS(4) = FNICVDT
        end if


        !***************************************************
        !FNITTAU     FIRST DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !         WITH RESPECT TO THE REDUCED TEMPERATURE [-]
        !         MULTIPLIED WITH TAU
        !****************************************************
        if (GETDERIVS(5) == 1) then
            !16.05.2013
            if (abs(gl%A_CV0(nrsubst)) .lt. 1.d-12) then    !noble gases
                FNICVT = -((gl%C_CV0(nrsubst) - 1.d0) / TAU_CV0  + C1_CV0)
            else
                FNICVT = -(gl%B_Cv0(nrsubst) - 1.d0) / TAU_CV0 + gl%A_CV0(nrsubst) * EPSILON * TR2 * ((-dlog(TAU_CV0) + dlog(OMEGA_lc))  &
                    * (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst) + 2.d0)/TR3     &
                    - (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst)+gl%F_CV0(nrsubst)+gl%G_CV0(nrsubst)+1.d0)/(OMEGA_lc*TR2)   &
                    - (gl%D_CV0(nrsubst) + gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst)) / (2.d0*TR_SRK*OMEGA2)   &
                    - (gl%E_CV0(nrsubst) + gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst)) / (3.d0 * OMEGA3) - TR_SRK * (gl%F_CV0(nrsubst) + gl%G_CV0(nrsubst))  &
                    / (4.d0 * OMEGA4) - gl%G_CV0(nrsubst) * TR2 / (5.d0 * OMEGA5) - 1.d0 / (gl%A_CV0(nrsubst) * TR2*TAU_CV0)) - C1_CV0  !changed by Theresa
            end if

            SETDERIVS(5) = - FNICVT*TAU_CV0
        end if


        !***************************************************
        !FNITTTAU2     SECOND DERIVATIVE OF THE REDUCED HELMHOLTZ ENERGY
        !          WITH RESPECT TO THE REDUCED TEMPERATURE [-]
        !          MULTIPLIED WITH TAU**2
        !****************************************************
        if (GETDERIVS(6) == 1) then
            if (abs(gl%A_CV0(nrsubst)) .lt. 1.d-12) then    !noble gases
                FNICVTT = gl%C_CV0(nrsubst) - 1.d0
            else                                            !all other fluids
                FNICVTT=gl%B_CV0(nrsubst)-1.d0+EPSILON*IOTA2*&
                    (1-ZETA*(gl%D_CV0(nrsubst)+gl%E_CV0(nrsubst)*IOTA+&
                    gl%F_CV0(nrsubst)*IOTA2+gl%G_CV0(nrsubst)*IOTA3))
            end if
            SETDERIVS(6) = - FNICVTT
        end if

        !end if

    end if

    END SUBROUTINE FNIDERIVS



    module subroutine ref_calc(gl,nrsubst, error)
    !**********************************************************
    !   REFERENCE STATE PROPERTIES DETERMINATION
    !**********************************************************
    !----------------------------------------------------------
    ! J. Gernert, 09.2009, Denmark
    !----------------------------------------------------------
    ! Calculates the reference state for enthalpy, entropy, internal energy
    ! Calculates density for all reference states
    ! Calculates the pressure and temperature in the reference state, if needed






    implicit none

    type(type_gl) :: gl


    integer:: errPsat, nrsubst
    integer:: error, iter
    integer:: errTsat, Iphase         !Andreas Aug. 2010
    double precision :: Psat, RhoV, RhoL, rhoredmixorg, tredmixorg


    iter = 0
    rhoredmixorg = 0.d0
    tredmixorg = 0.d0
    error = 0
    rhov = 0.d0
    rhol = 0.d0
    psat = 0.d0

    gl%ref = .true.

    rhoredmixorg = gl%rhoredmix
    tredmixorg = gl%tredmix
    gl%rhoredmix = gl%rhored(nrsubst)
    gl%tredmix = gl%tred(nrsubst)
    errTsat = 0
    errPsat = 0
    Iphase = 0
    select case (gl%refstate(nrsubst))
    case ('OTH')            ! reference state: real gas (given in fluid file)
        gl%rhoref(nrsubst) = rhomix_calc(gl,gl%Tref(nrsubst), gl%pref(nrsubst), 0.d0, Iphase, nrsubst) ! rho_calc(Tref(i), pref(i), errorflag, i)    ! reference density calculated from EOS
    case ('OT0')            ! reference state: ideal gas
        !rhoref(i) = pref(i)*1.d6/(REQ(i)/WM(i)*Tref(i))            ! reference state density calculated with ideal gas equation
        gl%rhoref(nrsubst) = gl%pref(nrsubst)*1.d6/(gl%REQ(nrsubst)*gl%Tref(nrsubst))/gl%factorpress            ! reference state density calculated with ideal gas equation
    case ('NBP')            !reference state: normal boiling point
        !for calulation of VLEpure it is necessary to put into rhoredmix the value of rhored of the pure component
        !call VLEpure(Tref(i),Psat,RhoV,RhoL,errPsat, i)        ! reference state density calculated with Maxwell criteria
        psat = 0.101325d0                                                         !NBP pressure is atmospheric pressure Andreas Aug.2010
        call Flash_Pure_PhaseBoundary(gl,psat,gl%Tref(nrsubst),RhoV,RhoL, 2, errTsat, iter, nrsubst)
        !call VLEpurePres(psat,Tref(i),RhoV,RhoL,errTsat, i)   !Routine for temperature iteration Andreas Aug.2010
        if (errTsat == 0) then ! succesfull
            gl%rhoref(nrsubst) = RhoL
        else
            error = -6000
        end if
    case ('ASH')            ! reference state: ASHRAE definition
        !for calulation of VLEpure it is necessary to put into rhoredmix the value of rhored of the pure component
        call Flash_Pure_PhaseBoundary(gl,psat,gl%Tref(nrsubst),RhoV,RhoL, 1, errPsat, iter, nrsubst)
        !call VLEpure(Tref(i),Psat,RhoV,RhoL,errPsat,i)        ! reference state density and pressure calculated with Maxwell criteria
        if (errPsat == 0) then ! succesfull
            gl%rhoref(nrsubst) = RhoL
            gl%pref(nrsubst) = Psat
        else
            error = -6000
        end if
    case ('IIR')            ! reference state: ASHRAE definition
        !for calulation of VLEpure it is necessary to put into rhoredmix the value of rhored of the pure component
        call Flash_Pure_PhaseBoundary(gl,psat,gl%Tref(nrsubst),RhoV,RhoL, 1, errPsat, iter, nrsubst)
        !call VLEpure(Tref(i),Psat,RhoV,RhoL,errPsat,i)        ! reference state density and pressure calculated with Maxwell criteria
        if (errPsat == 0) then ! succesfull
            gl%rhoref(nrsubst) = RhoL
            gl%pref(nrsubst) = Psat
        else
            error = -6000
        end if
    end select
    !if density could not return a solution catch error here
    if (nrsubst /= 0) then
        if (gl%rhoref(nrsubst) < 0) error = -6000
    endif
    gl%rhoredmix= rhoredmixorg
    gl%tredmix = tredmixorg
    
    gl%ref = .false.
    
    end subroutine




    module SUBROUTINE IDEALCONSTS(gl,nrsubst)
    !!*************************************************
    !!CALCULATION OF THE CONSTANTS C1, C2
    !!***************************************************
    ! FLORIAN DAUBER; DANMARK; 08.2009

    ! SUBROUTINE FOR THE CALCULATION OF THE INTEGRATION CONSTANTS C1 AND C2
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! nrsubst - INTEGER THAT GIVES INFORMATION ON THE FLUID TO BE CONSIDERED
    !-------------------------------------------------------------------------

    !USE module_general_eos_parameters !MODUL CALL
    !USE module_fluid_parameters !MODUL CALL
    !USE module_ideal_gas_coefficients !MODUL CALL




    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION :: H0,S0,DH,DS,TAUREF,R,T_AIT,AI !H_REF,S_REF
    integer::nrsubst
    integer,dimension(nderivsi) :: GETDERI
    DOUBLE PRECISION,dimension(nderivsi) :: FNIDER                !Vector gives the values of the calculated derivatives


    GETDERI=(/1,0,0,0,1,0,0,0,0,0/)       !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    ! FIRST C1 AND C2 ARE ZERO IN THE MODUL TO CALCULATE H0 AND S0
    ! OVERWRITE C1 AND C2 LATER IN THIS ROUTINE WITH THE CALCULATED VALUES

    !USE REFERENCE ENTHALPY HREF AND ENTROPY SREF FROM THE MODULE
    gl%ref = .true.
    R=gl%Req(nrsubst)!/WM(nrsubst)
    TAUREF =  gl%tc(nrsubst)/gl%TREF(nrsubst)
    if (gl%refstate(nrsubst) == 'OT0') then   ! if the reference state is ideal gas then H0 and S0 must be ideal gas properties, too!
        CALL FNIDERIVS(gl,gl%TREF(nrsubst),gl%RHOREF(nrsubst),GETDERI,FNIDER,nrsubst)   !subroutine calculates derivatives of ideal part
        AI = FNIDER(1)
        T_AIT = FNIDER(5)
        H0 = (1.d0 + T_AIT ) * R * gl%TREF(nrsubst)
        S0 = (T_AIT - AI) * R
        !    H0 = H_REF(gl,TREF(nrsubst),1.D-12, nrsubst)  !   ideal gas enthalpy equals real gas enthalpy at very small density
        !    S0 = S_REF(gl,TREF(nrsubst),1.D-9, nrsubst) + R*DLOG(1.D-9 / RHOREF(nrsubst) )
    else
        H0 = H_REF(gl,gl%TREF(nrsubst),gl%RHOREF(nrsubst), nrsubst)
        DH =  gl%HREF(nrsubst)-H0
        gl%C1(nrsubst) = DH/(R*gl%TREF(nrsubst)*TAUREF)
        S0 = S_REF(gl,gl%TREF(nrsubst),gl%RHOREF(nrsubst), nrsubst)  !+ R*DLOG(rhored(nrsubst) / RHOREF(nrsubst))  !rhored(nrsubst)
    end if
    DH =  gl%HREF(nrsubst)-H0
    DS = gl%SREF(nrsubst)- S0

    gl%C1(nrsubst) = DH/(R*gl%TREF(nrsubst)*TAUREF)
    gl%C2(nrsubst) = -DS/R

    if (allocated(gl%eos_coeff)) then
        if (gl%eos_coeff%cppcheck(nrsubst) == 'PHK') then
            gl%C1(nrsubst) = gl%cp0coeff(3,nrsubst)
            gl%C2(nrsubst) = gl%cp0coeff(2,nrsubst)
        end if
    end if

    gl%ref = .false.

    END SUBROUTINE



    !****************************************************************************
    module SUBROUTINE MIXDERIVSFNI(gl,T, D, GETDERIVS, SETMIXDERIVS)
    !****************************************************************************
    ! FLORIAN DAUBER; DANMARK; 08.2009

    ! SUBROUTINE FOR THE CALCULATION OF ALL DERIVATIVES OF THE IDEAL PART
    ! OF THE HELMHOLTZ FREE ENERGY FOR MIXTURES
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    ! GETDERIVS      - AN ARRAY WITH 6 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED IDEAL HELMHOLTZ ENERGY AS A FUNCTION OF D AND T
    !                2. 1ST DERIVATIVE WITH RESPECT TO D AT CONSTANT T
    !                3. 2ND DERIVATIVE WITH RESPECT TO D AT CONSTANT T
    !                4. 1ST MIXED DERIVATIVE WITH RESPECT TO D AND T
    !                5: 1ST DERIVATIVE WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                6: 2ND DERIVATIVE WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU**2
    ! nrsubst - INTEGER THAT GIVES INFORMATION ON THE FLUID TO BE CONSIDERED
    !
    ! OUTPUT PARAMETERS:
    ! SETDERIVS      - AN ARRAY WITH 6 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !                  AS INDECATED IN "GETDERIVS"
    !-------------------------------------------------------------------------
    !USE MODULE_FLUID_PARAMETERS
    !USE module_general_eos_parameters       !Added Andreas Aug 2010


    implicit none

    type(type_gl) :: gl


    integer, dimension(nderivsi)::GETDERIVS
    DOUBLE PRECISION, dimension(nderivsi)::SETMIXDERIVS
    DOUBLE PRECISION, dimension(nderivsi)::SETDERIVS
    double precision ::T, D, Rmix!, rhored_org T.Wiens commented rhoredmix_org Jan 2012      !Andreas rhoredmix_org added Aug 2010
    integer :: I, J

    SETMIXDERIVS=(/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)

    DO I = 1, gl%NCOMP !LOOP TO CALCULATE THE DERIVATIVE OF ALL FLUIDS IN THE MIXTURE. EACH IS MULTIPLIED WITH ITS MOLAR FRACTION AND ALL ARE SUMMED UP
        SETDERIVS=(/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)

        !--------------------------------------------
        ! The ideal part must be calculated with rhoredmix instead of rhored(nrsubst)??? Andreas Aug 2010
        !    rhored_org = rhored(I)
        !    If (rhoredmix  /= 0) then
        !    rhored(I) = rhoredmix
        !    end if
        CALL FNIDERIVS(gl,T, D, GETDERIVS, SETDERIVS, I) ! CALL THE CALCULATION ROUTINE FOR THE DERIVATIVES OF TH nrsubst I
        !    rhored(I) = rhored_org
        !--------------------------------------------
        !call R_mix_calc(gl,Rmix)

        DO J=1, 10
            IF (GETDERIVS(j) == 1) THEN !CHEK WHICH DERIVATIVE IS NEEDED
                SETMIXDERIVS(j)=SETMIXDERIVS(j)+ gl%MOLFRACTIONS(I)*SETDERIVS(j) !* gl%REQ(i)/Rmix
                IF (j == 1)THEN !Corrected in Aug. 2010, Jäger, Gernert   old version: !if (GETDERIVS(1) == 1) then
                    SETMIXDERIVS(1)=SETMIXDERIVS(1)+ gl%MOLFRACTIONS(I)*DLOG(gl%MOLFRACTIONS(I))
                END IF
            END IF
        END DO
    END DO

    END SUBROUTINE MIXDERIVSFNI


    end submodule impl
