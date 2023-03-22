    module special_terms

    use module_all_types
    
    implicit none
    
    contains

    
    
    !**************************************************************************
    !**************************************************************************
    SUBROUTINE DELTADERIVS(gl,DELTACALC, i, del, tau, DeltaDeriv, TauDeriv, nrsubst)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DISTANCE FUNCTION DELTA^bi, WHICH IS
    ! A PART OF THE NONANALYTIC TERM. THE SUBROUTINE CALCULATES DELTA^bi AND ITS
    ! DERIVATIVES WITH RESPECT TO TEMPERATURE AND DENSITY.
    ! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
    ! AS PUBLISHED BY (PP. 33):
    !                   Span, R.
    !                   Multiparameter Equations of State
    !                   Springer, 2000
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! i         - THE INDEX NUMBER OF THE NONANALYTIC TERM
    ! del       - REDUCED DENSITY
    ! tau       - INVERSE REDUCED TEMPERATURE
    ! TauDeriv  - THE DERIVATIVE OF DELTA WITH RESPECT TO TEMPERATURE
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO TEMPERATURE
    !             2: 2ND DERIVATIVE WITH RESPECT TO TEMPERATURE
    ! DeltaDeriv- THE DERIVATIVE OF DELTA WITH RESPECT TO DENSITY
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO DENSITY
    !             2: 2ND DERIVATIVE WITH RESPECT TO DENSITY
    !             3: 3RD DERIVATIVE WITH RESPECT TO DENSITY
    ! nrsubst   - Number of the component the derivatives are calculated for
    !
    ! OUTPUT PARAMETERS:
    ! DELTACALC - 8 x 1 VECTOR CONTAINING DELTA^bi AND THE DERIVATIVES INDICATED
    !             BY TauDeriv AND DeltaDeriv
    ! DELTACALC = (/DELTA, DELTAD, DELTADD, DELTAT, DELTATT, DELTADT, DELTADTT, DELTADDD/)
    !-------------------------------------------------------------------------
    ! 3rd derivative w.r.t to del added - J.G., 9.2011






    implicit none

    type(type_gl) :: gl


    integer:: i, DeltaDeriv, TauDeriv, nrsubst
    double precision:: del, tau
    double precision, dimension(10):: DELTACALC

    double precision :: DELTA, DELTAD, DELTADD, DELTADDD, DELTAT, DELTATT, DELTADT, DELTADTT, DELTADDT, DELTATTT
    double precision :: TERM1, TERM2, deldiff
    double precision, dimension(4) :: THETACALC

    DELTACALC = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
    THETACALC = (/0.D0, 0.D0, 0.D0, 0.D0/)

    ! GET THE THETA FUNCTION AND ITS DERIVATIVES
    CALL THETADERIVS(gl,THETACALC, i, del, tau, DeltaDeriv, nrsubst)

    !*****************************************************
    ! DELTA         DISTANCE FUNCTION OF THE NA-TERM
    !*****************************************************
    ! CALCULATION OF DELTA IS PERFORMED REGARDLESS WHICH DERIVATIVE
    ! IS NEEDED BECAUSE IT IS PART OF ALL DERIVATIVES
    deldiff = del - 1.D0

    if (DABS(deldiff) < 1.D-10) then        ! if del - 1 is smaller than 1.D-10 then set it to 1.D-10, otherwise div. by zero!
        deldiff = DSIGN(1.D-10, deldiff)    ! keep the signum of del - 1
    end if

    TERM1 = deldiff**2.d0                      ! This term is used in the derivatives
    TERM2 = gl%eos_coeff%eidna(i,nrsubst)* (TERM1**gl%eos_coeff%eitna(i,nrsubst))     ! This term is used in the derivatives as well
    DELTA = THETACALC(1)**2 + TERM2
    DELTACALC(1) = DELTA**gl%eos_coeff%eta(i,nrsubst)

    !*****************************************************
    if ((DeltaDeriv == 1) .OR. (DeltaDeriv == 2) .OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! DELTAD        1ST DER. OF DELTA W.R.T. DENSITY
        !*****************************************************
        DELTAD = 2.D0*THETACALC(1)*THETACALC(2) + 2.d0 * gl%eos_coeff%eitna(i,nrsubst) * TERM2/deldiff
        !1st derivative of the distance function DELTA^bi w.r.t. density
        DELTACALC(2) = gl%eos_coeff%eta(i,nrsubst) * DELTACALC(1)/DELTA * DELTAD
    end if

    !*****************************************************
    if ((DeltaDeriv == 2) .OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! DELTADD       2ND DER. OF DELTA W.R.T. DENSITY
        !*****************************************************
        DELTADD = 2.D0*(THETACALC(3) * THETACALC(1) + THETACALC(2)**2) + &
            & 2.D0 * gl%eos_coeff%eitna(i,nrsubst)*(2.d0*gl%eos_coeff%eitna(i,nrsubst) - 1.D0)*TERM2/deldiff**2.d0

        DELTACALC(3) = gl%eos_coeff%eta(i,nrsubst) * ((gl%eos_coeff%eta(i,nrsubst) - 1.D0) * DELTACALC(1)/DELTA**2.d0 *DELTAD**2.d0 + DELTACALC(1)/DELTA * DELTADD)
    end if

    !*****************************************************
    if ((TauDeriv == 1) .OR. (TauDeriv == 2) .OR. (TauDeriv == 3)) then
        !*****************************************************
        ! DELTAT        1ST DER. OF DELTA W.R.T. TEMPERATURE
        !*****************************************************
        DELTAT = -2.D0 * THETACALC(1) * gl%eos_coeff%eta(i,nrsubst) * DELTACALC(1)/DELTA
        DELTACALC(4) = DELTAT
    end if

    !*****************************************************
    if ((TauDeriv == 2) .or. (TauDeriv == 3)) then
        !*****************************************************
        ! DELTATT       2ND DER. OF DELTA W.R.T. TEMPERATURE
        !*****************************************************
        DELTATT = 2.D0*gl%eos_coeff%eta(i,nrsubst)*DELTACALC(1)/DELTA*(1.d0+2.d0*THETACALC(1)**2*(gl%eos_coeff%eta(i,nrsubst) - 1.D0)/DELTA)

        DELTACALC(5) = DELTATT
    end if

    !*****************************************************
    if (((TauDeriv == 1) .AND. (DeltaDeriv == 1)) .or. &
        ((TauDeriv == 1) .AND. (DeltaDeriv == 2)) .or. &
        ((TauDeriv == 2) .AND. (DeltaDeriv == 1))) then
        !*****************************************************
        ! DELTADT       1ST MIXED DER. OF DELTA W.R.T.
        !               DENSITY AND TEMPERATURE
        !*****************************************************
        DELTADT = -2.D0 * gl%eos_coeff%eta(i,nrsubst) * DELTACALC(1) / DELTA * ((gl%eos_coeff%eta(i,nrsubst) - 1.D0) &
            & * THETACALC(1)/DELTA*DELTAD + THETACALC(2))

        DELTACALC(6) = DELTADT
    end if

    !*****************************************************
    if ((TauDeriv == 2) .AND. (DeltaDeriv == 1)) then
        !*****************************************************
        ! DELTADTT      2ND MIXED DER. OF DELTA W.R.T.
        !               DENSITY, TEMPERATURE AND TEMPERATURE
        !*****************************************************
        DELTADTT = 2.d0 * gl%eos_coeff%eta(i,nrsubst) * (gl%eos_coeff%eta(i,nrsubst) - 1.D0) * DELTACALC(1) / DELTA**2.d0 * (DELTAD * (1.d0 + 2.d0 * (gl%eos_coeff%eta(i,nrsubst) - 2.D0)/DELTA * THETACALC(1)**2) + 4.d0 * THETACALC(1) * THETACALC(2))

        DELTACALC(7) = DELTADTT
    end if

    !*****************************************************
    if (DeltaDeriv == 3) then
        !*****************************************************
        ! DELTADDD      3RD DER. OF DELTA W.R.T. DENSITY
        !*****************************************************
        DELTADDD = 6.D0*THETACALC(2)*THETACALC(3) + 2.D0*THETACALC(1)*THETACALC(4) &
            & + TERM2/(deldiff**3.d0) * 2.d0*gl%eos_coeff%eitna(i,nrsubst) * (2.d0*gl%eos_coeff%eitna(i,nrsubst) - 1.d0) * (2.d0*gl%eos_coeff%eitna(i,nrsubst) - 2.d0)
        DELTACALC(8) = gl%eos_coeff%eta(i,nrsubst)* (DELTA**(gl%eos_coeff%eta(i,nrsubst) - 3.D0) * (gl%eos_coeff%eta(i,nrsubst) - 1.D0) * (gl%eos_coeff%eta(i,nrsubst) - 2.D0) * DELTAD**3.d0 + &
            DELTA**(gl%eos_coeff%eta(i,nrsubst) - 2.D0) * 3.d0 * (gl%eos_coeff%eta(i,nrsubst) - 1.D0) * DELTAD * DELTADD + &
            DELTA**(gl%eos_coeff%eta(i,nrsubst) - 1.D0) * DELTADDD)
    end if

    !*****************************************************
    if (TauDeriv == 3) then
        !*****************************************************
        ! DELTATTT      3RD DER. OF DELTA W.R.T. TEMPERATURE
        !*****************************************************
        DELTATTT = - 4.d0 * THETACALC(1) * gl%eos_coeff%eta(i,nrsubst) * (gl%eos_coeff%eta(i,nrsubst) - 1.d0) * delta ** (gl%eos_coeff%eta(i,nrsubst) - 2.d0)  &
            & * (3.d0 + 2.d0 * THETACALC(1)**2 * (gl%eos_coeff%eta(i,nrsubst) - 2.d0) / delta )
        DELTACALC(9) = DELTATTT
    end if

    !*****************************************************
    if ((TauDeriv == 1) .AND. (DeltaDeriv == 2)) then
        !*****************************************************
        ! DELTADDT      2ND MIXED DER. OF DELTA W.R.T.
        !               DENSITY, DENSITY AND TEMPERATURE
        !*****************************************************
        DELTADDT = -2.d0 * gl%eos_coeff%eta(i,nrsubst) * (delta ** (gl%eos_coeff%eta(i,nrsubst) - 1.d0) * THETACALC(3) +  &
            delta ** (gl%eos_coeff%eta(i,nrsubst) - 2.d0) * (gl%eos_coeff%eta(i,nrsubst) - 1.d0) * (THETACALC(1) * DELTADD + 2.d0 * THETACALC(2) * DELTAD) +  &
            delta ** (gl%eos_coeff%eta(i,nrsubst) - 3.d0) * (gl%eos_coeff%eta(i,nrsubst) - 1.d0) * (gl%eos_coeff%eta(i,nrsubst) - 2.d0) * THETACALC(1) * DELTAD**2.d0)

        DELTACALC(10) = DELTADDT
    end if

    END SUBROUTINE DELTADERIVS

    !**************************************************************************
    !
    !**************************************************************************
    SUBROUTINE THETADERIVS(gl,THETACALC, i, del, tau, DeltaDeriv, nrsubst)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE THETA PART IN THE DISTANCE FUNCTION,
    ! WHICH IS A PART OF THE NONANALYTIC TERM. THE SUBROUTINE CALCULATES THETA
    ! AND ITS DERIVATIVES WITH RESPECT TO TEMPERATURE AND DENSITY.
    ! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
    ! AS PUBLISHED BY (PP. 33):
    !                   Span, R.
    !                   Multiparameter Equations of State
    !                   Springer, 2000
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! i         - THE INDEX NUMBER OF THE NONANALYTIC TERM
    ! del       - REDUCED DENSITY
    ! tau       - INVERSE REDUCED TEMPERATURE
    ! DeltaDeriv- THE DERIVATIVE OF THETA WITH RESPECT TO DENSITY
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO DENSITY
    !             2: 2ND DERIVATIVE WITH RESPECT TO DENSITY
    !             3: 3RD DERIVATIVE WITH RESPECT TO DENSITY
    ! nrsubst   - Number of the component the derivatives are calculated for
    !
    ! OUTPUT PARAMETERS:
    ! THETACALC - 3 x 1 VECTOR CONTAINING THETA AND THE DERIVATIVES INDICATED
    !             BY TauDeriv AND DeltaDeriv
    !-------------------------------------------------------------------------
    ! 3rd derivative w.r.t to del added - J.G., 9.2011






    implicit none

    type(type_gl) :: gl


    integer:: i, DeltaDeriv, nrsubst
    double precision:: del, tau
    double precision, dimension(4):: THETACALC

    double precision :: THETA, THETAD, THETADD, THETADDD, TERM1, TERM2
    double precision :: deldiff

    THETACALC = (/0.D0, 0.D0, 0.D0, 0.D0/)

    !*****************************************************
    ! THETA         THETA FUNCTION OF THE NA-TERM
    !*****************************************************
    ! CALCULATION OF THETA IS PERFORMED REGARDLESS WHICH DERIVATIVE
    ! IS NEEDED BECAUSE IT IS PART OF ALL DERIVATIVES
    deldiff = del - 1.D0

    if (DABS(deldiff) < 1.D-10) then        ! if del - 1 is smaller than 1.D-10 then set it to 1.D-10
        deldiff = DSIGN(1.D-10, deldiff)    ! keep the signum of del - 1
    end if

    TERM1 = deldiff**2.d0                              !this term can be used in the derivatives
    TERM2 = gl%eos_coeff%gam(i,nrsubst)*(TERM1**(1.D0/(2.D0*gl%eos_coeff%beta(i,nrsubst))))   !this term can be used in the derivatives
    THETA = 1.D0 - tau + TERM2
    THETACALC(1) = THETA

    !*****************************************************
    if ((DeltaDeriv == 1) .OR. (DeltaDeriv == 2) .OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! THETAD        1ST DER. OF THETA W.R.T. DENSITY
        !*****************************************************
        THETAD = TERM2/(deldiff * gl%eos_coeff%beta(i,nrsubst))
        THETACALC(2) = THETAD
    end if

    !*****************************************************
    if ((DeltaDeriv == 2).OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! THETADD       2ND DER. OF THETA W.R.T. DENSITY
        !*****************************************************
        THETADD = THETAD / deldiff * (1.D0 / gl%eos_coeff%beta(i,nrsubst) - 1.D0)
        THETACALC(3) = THETADD
    end if

    !*****************************************************
    if (DeltaDeriv == 3) then
        !*****************************************************
        ! THETADDD      3RD DER. OF THETA W.R.T. DENSITY
        !*****************************************************
        THETADDD = THETADD / deldiff * (1.D0 / gl%eos_coeff%beta(i,nrsubst) - 2.D0)
        THETACALC(4) = THETADDD
    end if

    ! The derivatives of THETA with respect to tau are not calculated here
    ! because they are just -1 and 0

    END SUBROUTINE THETADERIVS

    !**************************************************************************
    !
    !**************************************************************************
    SUBROUTINE PSIDERIVS(gl,PSICALC, i, del, tau, DeltaDeriv, TauDeriv, nrsubst)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE PSI PART IN THE NONANALYTIC TERMS.
    ! THE SUBROUTINE CALCULATES PSI AND ITS DERIVATIVES WITH RESPECT TO TEMPERATURE
    ! AND DENSITY. THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ
    ! FREE ENERGY AS PUBLISHED BY (PP. 33):
    !                   Span, R.
    !                   Multiparameter Equations of State
    !                   Springer, 2000
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! i         - THE INDEX NUMBER OF THE NONANALYTIC TERM
    ! del       - REDUCED DENSITY
    ! tau       - INVERSE REDUCED TEMPERATURE
    ! TauDeriv  - THE DERIVATIVE OF DELTA WITH RESPECT TO TEMPERATURE
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO TEMPERATURE
    !             2: 2ND DERIVATIVE WITH RESPECT TO TEMPERATURE
    ! DeltaDeriv- THE DERIVATIVE OF THETA WITH RESPECT TO DENSITY
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO DENSITY
    !             2: 2ND DERIVATIVE WITH RESPECT TO DENSITY
    !             3: 3RD DERIVATIVE WITH RESPECT TO DENSITY
    ! nrsubst   - Number of the component the derivatives are calculated for
    !
    ! OUTPUT PARAMETERS:
    ! PSICALC   - 8 x 1 VECTOR CONTAINING THETA AND THE DERIVATIVES INDICATED
    !             BY TauDeriv AND DeltaDeriv
    !-------------------------------------------------------------------------
    ! 3rd derivative w.r.t to del added - J.G., 9.2011






    implicit none

    type(type_gl) :: gl


    integer:: i, DeltaDeriv, TauDeriv, nrsubst
    double precision:: del, tau
    double precision, dimension(10):: PSICALC

    double precision :: PSI, PSID, PSIDD, PSIDDD, PSIT, PSITT, PSIDT, PSIDTT, TERM1, TERM2, PSIDDT, PSITTT
    double precision :: deldiff, taudiff

    PSICALC = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)

    !*****************************************************
    ! PSI           EXPONENTIAL FUNCTION OF THE NA-TERM
    !*****************************************************
    ! CALCULATION OF PSI IS PERFORMED REGARDLESS WHICH DERIVATIVE
    ! IS NEEDED BECAUSE IT IS PART OF ALL DERIVATIVES
    deldiff = del - 1.D0
    taudiff = tau - 1.D0

    if (DABS(deldiff) < 1.D-10) then        ! if del - 1 is smaller than 1.D-10 then set it to 1.D-10
        deldiff = DSIGN(1.D-10, deldiff)    ! keep the signum of del - 1
    end if
    if (DABS(taudiff) < 1.D-10) then        ! if tau - 1 is smaller than 1.D-10 then set it to 1.D-10
        taudiff = DSIGN(1.D-10, taudiff)    ! keep the signum of tau - 1
    end if

    TERM1 = gl%eos_coeff%eps(i,nrsubst)*deldiff**2.d0    ! exponential density term
    TERM2 = gl%eos_coeff%etana(i,nrsubst)*taudiff**2.d0    ! exponential temperature term

    PSI = dexp(- TERM1 - TERM2)

    !*****************************************************
    if (abs(PSI) < 1.D-10) then             ! checks, if psi is smaller than a very small value
        return                              ! if that's the case than the calculation of the nonanalytic
    end if                                  ! term and all of its derivatives is aborted
    !*****************************************************

    PSICALC(1) = PSI

    !*****************************************************
    If ((DeltaDeriv == 1) .OR. (DeltaDeriv == 2) .OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! PSID          1ST DER. OF PSI W.R.T. DENSITY
        !*****************************************************
        PSID = -2.D0*gl%eos_coeff%eps(i,nrsubst)*deldiff*PSI
        PSICALC(2) = PSID
    end if

    !*****************************************************
    if ((DeltaDeriv == 2) .OR. (DeltaDeriv == 3)) then
        !*****************************************************
        ! PSIDD         2ND DER. OF PSI W.R.T. DENSITY
        !*****************************************************
        PSIDD = PSI * 2.D0*gl%eos_coeff%eps(i,nrsubst) * (2.D0*TERM1 - 1.D0)
        PSICALC(3) = PSIDD
    end if

    !*****************************************************
    if ((TauDeriv == 1) .OR. (TauDeriv == 2) .OR. (TauDeriv == 3)) then
        !*****************************************************
        ! PSIT          1ST DER. OF PSI W.R.T. TEMPERATURE
        !*****************************************************
        PSIT = -2.D0*gl%eos_coeff%etana(i,nrsubst)*taudiff*PSI
        PSICALC(4) = PSIT
    end if

    !*****************************************************
    if ((TauDeriv == 2) .or. (TauDeriv == 3)) then
        !*****************************************************
        ! PSITT         2ND DER. OF PSI W.R.T. TEMPERATURE
        !*****************************************************
        PSITT = PSI * 2.D0*gl%eos_coeff%etana(i,nrsubst) * (2.D0*TERM2 - 1.D0)
        !    PSITT = PSIT**2.d0/PSI - 2.D0*etana(nrsubst, i)*PSI
        PSICALC(5) = PSITT
    end if

    !*****************************************************
    if (((DeltaDeriv == 1) .AND. (TauDeriv == 1)) .or. &
        ((DeltaDeriv == 2) .AND. (TauDeriv == 1)) .or. &
        ((DeltaDeriv == 1) .AND. (TauDeriv == 2))) then
        !*****************************************************
        ! PSIDT         1ST MIXED DER. OF PSI W.R.T.
        !               DENSITY AND TEMPERATURE
        !*****************************************************
        PSIDT = PSID*PSIT/PSI
        PSICALC(6) = PSIDT
    end if

    !*****************************************************
    if ((DeltaDeriv == 1) .AND. (TauDeriv == 2)) then
        !*****************************************************
        ! PSIDTT        2ND MIXED DER. OF PSI W.R.T.
        !               DENSITY, TEMPERATURE AND TEMPERATURE
        !*****************************************************
        PSIDTT = PSID*PSITT/PSI
        PSICALC(7) = PSIDTT
    end if

    !*****************************************************
    if (DeltaDeriv == 3) then
        !*****************************************************
        ! PSIDDD        3RD DER. OF PSI W.R.T. DENSITY
        !*****************************************************
        PSIDDD = PSI * 4.D0*gl%eos_coeff%eps(i,nrsubst)**2*(del - 1.d0)*(3.d0 - 2.d0*TERM1)
        PSICALC(8) = PSIDDD
    end if

    !*****************************************************
    if (TauDeriv == 3) then
        !*****************************************************
        ! PSITTT        3RD DER. OF PSI W.R.T. TEMPERATURE
        !*****************************************************
        !PSITTT = PSI * 4.D0*etana(nrsubst, i)**2*(tau - 1.d0)*(3.d0 - 2.d0*TERM2)
        PSITTT = PSI * 4.D0*gl%eos_coeff%etana(i,nrsubst)**2*(tau - 1.d0)*(3.d0 - 2.d0*TERM2)
        PSICALC(9) = PSITTT
    end if

    !*****************************************************
    if ((DeltaDeriv == 2) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! PSIDDT        2ND MIXED DER. OF PSI W.R.T.
        !               DENSITY, DENSITY AND TEMPERATURE
        !*****************************************************
        PSIDDT = PSIDD*PSIT/PSI
        PSICALC(10) = PSIDDT
    end if

    END SUBROUTINE PSIDERIVS

    
    !**************************************************************************
    subroutine SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DeltaDeriv, TauDeriv, nrsubst)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE SPECIAL GAUSSIAN BELL-SHAPED TERMS AND THEIR
    ! DERIVATIVES WITH RESPECT TO TAU AND DELTA
    ! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
    ! not published yet (20.03.2018)
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! i         - THE INDEX NUMBER OF THE SPECIAL GAUSSIAN BELL-SHAPED TERM
    ! del       - REDUCED DENSITY
    ! tau       - INVERSE REDUCED TEMPERATURE
    ! nrsubst   - Number of the component the derivatives are calculated for
    !
    ! OUTPUT PARAMETERS:
    ! SGBSCALC   - THE SPECIAL GAUSSIAN BELL-SHAPED TERM OR ITS DERIVATIVE
    !-------------------------------------------------------------------------
    !added by T.N. 03.2018



    implicit none

    type(type_gl) :: gl

    integer:: k,i, nrsubst, DeltaDeriv, TauDeriv
    double precision:: tau, del, expo
    double precision:: SGBSCALC, SGBS, SGBSD, SGBSDD, SGBSDDD, SGBST, SGBSTT, SGBSTTT, SGBSDT, SGBSDTT, SGBSDDT, part1, part2, part3, part4, part5, SGBSF, SGBSTT_1, SGBSTT_2,h, tau_org
    double precision, dimension (20):: TAUGAM, DELEPS



    !*****************************************************
    ! SGBS              Special GBS Terms
    !*****************************************************

    !!Summation over the special Gaussian bell-shaped terms

    k = i - gl%eos_coeff%nreg(nrsubst)
    DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
    TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)

    expo = dlog(del)*gl%eos_coeff%di(i,nrsubst) + dlog(tau)*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + 1.d0/(gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))

    SGBS = 0.d0
    if (expo < 200) SGBS = gl%eos_coeff%ni(i,nrsubst) * dexp (expo)   !prevents double precision overflow

    SGBSCALC = SGBS


    if ((DeltaDeriv == 1) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! SGBSD       1ST DER. OF THE Special GBS TERM W.R.T. DENSITY
        !*****************************************************

        SGBSD =  SGBS * (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))

        SGBSCALC = SGBSD
    end if


    if ((DeltaDeriv == 2) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! SGBSDD      2ND DER. OF THE Special GBS TERM W.R.T. DENSITY
        !*****************************************************

        SGBSDD =  SGBS * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))**2 &
            & - gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0)*(gl%eos_coeff%pli(k,nrsubst)-1.D0))

        SGBSCALC = SGBSDD
    end if


    if ((DeltaDeriv == 3) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! SGBSDDD     3RD DER. OF THE Special GBS TERM W.R.T. DENSITY

        !*****************************************************
        SGBSDDD =  SGBS * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 1.d0))**3 &
            & - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 + 2.d0*gl%eos_coeff%di(i,nrsubst) + 3.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*del*gl%eos_coeff%pli(k,nrsubst)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)* &
            & DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0) - 3.d0*gl%eos_coeff%di(i,nrsubst)*del*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0) &
            & + 3.D0*del*del*del*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)**2.d0*gl%eos_coeff%eta(k,nrsubst)**2*DELEPS(k)**(2.d0*gl%eos_coeff%pli(k,nrsubst) - 3.d0) &
            & + (gl%eos_coeff%pli(k,nrsubst) - 2.d0)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 3.d0))



        SGBSCALC = SGBSDDD
    end if


    if ((DeltaDeriv == 0) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! SGBST       1ST DER. OF THE Special GBS TERM W.R.T. TEMPERATURE
        !*****************************************************


        SGBST =  SGBS *(gl%eos_coeff%ti(i,nrsubst) - (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))

        SGBSCALC = SGBST
    end if



    if ((DeltaDeriv == 0) .AND. (TauDeriv == 2)) then
        !*****************************************************
        ! SGBSTT      2ND DER. OF THE Special GBS TERM W.R.T. TEMPERATURE
        !*****************************************************

        SGBSTT =  SGBS *(gl%eos_coeff%ti(i,nrsubst)**2.d0 - gl%eos_coeff%ti(i,nrsubst)*(((2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * &
            & ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))) + 1.d0) + &
            & ((tau*tau* (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0)) * (2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%etana(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
            & gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
            & gl%eos_coeff%tli(k,nrsubst)) - (gl%eos_coeff%tli(k,nrsubst) - 1.d0)*gl%eos_coeff%etana(k,nrsubst)**2.d0)) * ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-4.d0) )))

        SGBSCALC = SGBSTT
    end if



    if ((DeltaDeriv == 0) .AND. (TauDeriv == 3)) then
        !*****************************************************
        ! SGBSTTT     3RD DER. OF THE Special GBS TERM W.R.T. TEMPERATURE
        !*****************************************************

        !   SGBSF = gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst) !repeating Factor
        !
        !part1 = 3.d0*gl%eos_coeff%ti(i,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*tau*((gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(2.d0*gl%eos_coeff%tli(k,nrsubst) - 2.d0)*(2.d0*SGBSF+1.d0)) - (SGBSF**(2.d0)*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0)))
        !
        !part2 = - gl%eos_coeff%beta(k,nrsubst)**(3.d0)*gl%eos_coeff%tli(k,nrsubst)**(3.d0)*tau*tau*tau*TAUGAM(k)**(3.d0*gl%eos_coeff%tli(k,nrsubst) - 3.d0)*(6.d0*SGBSF**(2.d0)+6.d0*SGBSF+1.d0)
        !
        !part3 = 6.d0*gl%eos_coeff%beta(k,nrsubst)**(2.d0)*gl%eos_coeff%tli(k,nrsubst)**(2.d0)*tau*tau*tau*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(2.d0*gl%eos_coeff%tli(k,nrsubst) - 3.d0)*(SGBSF+5.d-1)
        !
        !part4 = - gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*tau*tau*(gl%eos_coeff%tli(k,nrsubst) - 2.d0)*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 3.d0)
        !
        !part5 = (gl%eos_coeff%ti(i,nrsubst) - 1.d0)*gl%eos_coeff%ti(i,nrsubst)*( - 3.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)+(gl%eos_coeff%ti(i,nrsubst) - 2.d0)*SGBSF**(2.d0))
        !
        !                SGBSTTT =  SGBS * (part1*SGBSF**(-4.d0) + part2*SGBSF**(-6.d0) + part3*SGBSF**(-4.d0) + part4*SGBSF**(-2.d0) + part5*SGBSF**(-2.d0))


        SGBSF = gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst) !repeating Factor


        SGBSTTT =  SGBS * ((3.d0*gl%eos_coeff%ti(i,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*tau*((gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(2.d0*gl%eos_coeff%tli(k,nrsubst) - 2.d0)*(2.d0*SGBSF+1.d0)) - &
            & (SGBSF**(2.d0)*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0))))*SGBSF**(-4.d0) + &
            & (- gl%eos_coeff%beta(k,nrsubst)**(3.d0)*gl%eos_coeff%tli(k,nrsubst)**(3.d0)*tau*tau*tau*TAUGAM(k)**(3.d0*gl%eos_coeff%tli(k,nrsubst) - 3.d0)*(6.d0*SGBSF**(2.d0)+6.d0*SGBSF+1.d0))*SGBSF**(-6.d0) + &
            & ( 6.d0*gl%eos_coeff%beta(k,nrsubst)**(2.d0)*gl%eos_coeff%tli(k,nrsubst)**(2.d0)*tau*tau*tau*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(2.d0*gl%eos_coeff%tli(k,nrsubst) - 3.d0)*(SGBSF+5.d-1))*SGBSF**(-4.d0) + &
            & ( - gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*tau*tau*(gl%eos_coeff%tli(k,nrsubst) - 2.d0)*(gl%eos_coeff%tli(k,nrsubst) - 1.d0)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 3.d0))*SGBSF**(-2.d0) + &
            & ((gl%eos_coeff%ti(i,nrsubst) - 1.d0)*gl%eos_coeff%ti(i,nrsubst)*( - 3.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)+(gl%eos_coeff%ti(i,nrsubst) - 2.d0)*SGBSF**(2.d0)))*SGBSF**(-2.d0))

        !
        !h=0.000000001
        !tau_org=tau
        !tau=tau+h
        !
        !k = i - gl%eos_coeff%nreg(nrsubst)
        !                DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
        !                TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
        !
        !                    expo = dlog(del)*gl%eos_coeff%di(i,nrsubst) + dlog(tau)*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + 1.d0/(gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))
        !
        !                SGBS = 0.d0
        !                if (expo < 200) SGBS = gl%eos_coeff%ni(i,nrsubst) * dexp (expo)   !prevents double precision overflow
        !
        !                    SGBSTT_1 =  SGBS*tau_org**(-2.d0) *(gl%eos_coeff%ti(i,nrsubst)**2.d0 - gl%eos_coeff%ti(i,nrsubst)*(((2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * &
        !                       & ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))) + 1.d0) + &
        !                       & ((tau*tau* (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0)) * (2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%etana(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
        !                       & gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
        !                       & gl%eos_coeff%tli(k,nrsubst)) - (gl%eos_coeff%tli(k,nrsubst) - 1.d0)*gl%eos_coeff%etana(k,nrsubst)**2.d0)) * ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-4.d0) )))
        !
        !                    tau=tau-(2.d0*h)
        !
        !                    k = i - gl%eos_coeff%nreg(nrsubst)
        !                DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
        !                TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
        !
        !                    expo = dlog(del)*gl%eos_coeff%di(i,nrsubst) + dlog(tau)*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + 1.d0/(gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))
        !
        !                SGBS = 0.d0
        !                if (expo < 200) SGBS = gl%eos_coeff%ni(i,nrsubst) * dexp (expo)   !prevents double precision overflow
        !
        !
        !                    SGBSTT_2 =  SGBS*tau_org**(-2.d0)  *(gl%eos_coeff%ti(i,nrsubst)**2.d0 - gl%eos_coeff%ti(i,nrsubst)*(((2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * &
        !                       & ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))) + 1.d0) + &
        !                       & ((tau*tau* (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0)) * (2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%etana(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
        !                       & gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
        !                       & gl%eos_coeff%tli(k,nrsubst)) - (gl%eos_coeff%tli(k,nrsubst) - 1.d0)*gl%eos_coeff%etana(k,nrsubst)**2.d0)) * ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-4.d0) )))
        !
        !
        !                    SGBSTTT=(SGBSTT_1-SGBSTT_2)*tau_org**(3.d0) /(2.d0*h)

        SGBSCALC = SGBSTTT
    end if


    if ((DeltaDeriv == 1) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! SGBSDT        1ST MIXED DER. OF THE Special GBS TERM W.R.T.
        !               DENSITY AND TEMPERATURE
        !*****************************************************

        SGBSDT =  SGBS * (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0)) *&
            & (gl%eos_coeff%ti(i,nrsubst) - (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))

        SGBSCALC = SGBSDT
    end if



    if ((DeltaDeriv == 2) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! SGBSDDT    2ND MIXED DER. OF THE Special GBS TERM W.R.T.
        !               DENSITY, DENSITY, AND TEMPERATURE
        !*****************************************************

        SGBSDDT =  SGBS * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))**2 &
            & - gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0)*(gl%eos_coeff%pli(k,nrsubst)-1.D0)) * &
            & (gl%eos_coeff%ti(i,nrsubst) - (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))

        SGBSCALC = SGBSDDT
    end if



    if ((DeltaDeriv == 1) .AND. (TauDeriv == 2)) then
        !*****************************************************
        ! SGBSDTT    2ND MIXED DER. OF THE Special GBS TERM W.R.T.
        !               DENSITY, TEMPERATURE AND TEMPERATURE
        !*****************************************************

        SGBSDTT =  SGBS * (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0)) * &
            & (gl%eos_coeff%ti(i,nrsubst)**2.d0 - gl%eos_coeff%ti(i,nrsubst)*(((2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 1.d0)) * &
            & ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-2.d0))) + 1.d0) + &
            & ((tau*tau* (gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst) - 2.d0)) * (2.d0*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%etana(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
            & gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) * (gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + &
            & gl%eos_coeff%tli(k,nrsubst)) - (gl%eos_coeff%tli(k,nrsubst) - 1.d0)*gl%eos_coeff%etana(k,nrsubst)**2.d0)) * ((gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst) + gl%eos_coeff%etana(k,nrsubst))**(-4.d0) )))


        SGBSCALC = SGBSDTT
    end if

    end subroutine SGBSTERMDERIVS
    !**************************************************************************


    !**************************************************************************
    subroutine NATERMDERIVS(gl,NACALC, i, del, tau, DeltaDeriv, TauDeriv, region_int, nrsubst)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE NONANALYTIC TERMS AND THEIR
    ! DERIVATIVES WITH RESPECT TO TAU AND DELTA
    ! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
    ! AS PUBLISHED BY (PP. 24):
    !                   Span, R.
    !                   Multiparameter Equations of State
    !                   Springer, 2000
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! i         - THE INDEX NUMBER OF THE NONANALYTIC TERM
    ! del       - REDUCED DENSITY
    ! tau       - INVERSE REDUCED TEMPERATURE
    ! TauDeriv  - THE DERIVATIVE OF THE NONANALYTIC TERM WITH RESPECT TO TEMPERATURE
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO TAU
    !             2: 2ND DERIVATIVE WITH RESPECT TO TAU
    ! DeltaDeriv- THE DERIVATIVE OF THE NONANALYTIC TERM WITH RESPECT TO DENSITY
    !             0: NO DERIVATIVE
    !             1: 1ST DERIVATIVE WITH RESPECT TO DELTA
    !             2: 2ND DERIVATIVE WITH RESPECT TO DELTA
    !             3: 3RD DERIVATIVE WITH RESPECT TO DELTA
    ! nrsubst   - Number of the component the derivatives are calculated for
    !
    ! OUTPUT PARAMETERS:
    ! NACALC   - THE NONANALYTIC TERM OR ITS DERIVATIVE
    ! region   - INTEGER, EITHER 0 OR 1;
    !            0: INDICATES, THAT THE NA-TERM IS NEGLIGIBLE (OUTSIDE CRIT. REGION)
    !            1: INDICATES, THAT THE NA-TERM IS NEEDED (CRIT. REGION)
    !-------------------------------------------------------------------------
    ! 3rd derivative w.r.t to del added - J.G., 9.2011






    implicit none

    type(type_gl) :: gl


    integer:: j,i, TauDeriv, DeltaDeriv, nrsubst
    double precision:: tau, del
    double precision:: NACALC
    integer :: region_int
    double precision :: NA, NAD, NADD, NAT, NATT, NADT, NADTT, NADDD, NADDT, NATTT
    double precision, dimension(10):: PSICALC, DELTACALC


    DELTACALC = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
    PSICALC = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
    region_int = 1
    j = i + gl%eos_coeff%nreg(nrsubst)

    CALL PSIDERIVS(gl,PSICALC, i, del, tau, DeltaDeriv, TauDeriv, nrsubst)

    !*****************************************************
    if (abs(PSICALC(1)) < 1.D-10) then             ! checks, if psi is smaller than a very small value
        NACALC = 0.D0                              ! if that's the case than the calculation of the nonanalytic
        region_int = 0                                 ! term and all of its derivatives is aborted
        return
    end if
    !*****************************************************

    CALL DELTADERIVS(gl,DELTACALC, i, del, tau, DeltaDeriv, TauDeriv, nrsubst)

    !*****************************************************
    if ((DeltaDeriv == 0) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! NA        NONANALYTIC TERM
        !*****************************************************
        NA = gl%eos_coeff%ni(j,nrsubst)*del*DELTACALC(1)*PSICALC(1)
        NACALC = NA

        !*****************************************************
    elseif ((DeltaDeriv == 1) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! NAD       1ST DER. OF THE NA TERM W.R.T. DENSITY
        !*****************************************************
        NAD = gl%eos_coeff%ni(j,nrsubst)*(DELTACALC(1)*(PSICALC(1) + del*PSICALC(2)) + DELTACALC(2)*del*PSICALC(1))
        NACALC = NAD

        !*****************************************************
    elseif ((DeltaDeriv == 2) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! NADD      2ND DER. OF THE NA TERM W.R.T. DENSITY
        !*****************************************************
        NADD = gl%eos_coeff%ni(j,nrsubst)*(DELTACALC(1)*(2.D0*PSICALC(2) + del*PSICALC(3)) &
            & + 2.D0*DELTACALC(2)*(PSICALC(1) + del*PSICALC(2)) + DELTACALC(3)*del*PSICALC(1))
        NACALC = NADD

        !*****************************************************
    elseif ((DeltaDeriv == 0) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! NAT       1ST DER. OF THE NA TERM W.R.T. TEMPERATURE
        !*****************************************************
        NAT = gl%eos_coeff%ni(j,nrsubst)*del*(DELTACALC(4)*PSICALC(1) + DELTACALC(1)*PSICALC(4))
        NACALC = NAT

        !*****************************************************
    elseif ((DeltaDeriv == 0) .AND. (TauDeriv == 2)) then
        !*****************************************************
        ! NATT      2ND DER. OF THE NA TERM W.R.T. TEMPERATURE
        !*****************************************************
        NATT = gl%eos_coeff%ni(j,nrsubst)*del*(DELTACALC(5)*PSICALC(1) + 2.D0*DELTACALC(4)*PSICALC(4) + DELTACALC(1)*PSICALC(5))
        NACALC = NATT

        !*****************************************************
    elseif ((DeltaDeriv == 1) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! NADT      1ST MIXED DER. OF THE NA TERM W.R.T.
        !           DENSITY AND TEMPERATURE
        !*****************************************************
        NADT = gl%eos_coeff%ni(j,nrsubst)*(DELTACALC(1)*(PSICALC(4) + del*PSICALC(6)) + del* DELTACALC(2)*PSICALC(4)&
            & + DELTACALC(4)*(PSICALC(1) + del*PSICALC(2)) + DELTACALC(6)*del*PSICALC(1))

        NACALC = NADT

        !*****************************************************
    elseif ((DeltaDeriv == 1) .AND. (TauDeriv == 2)) then
        !*****************************************************
        ! NADTT     2ND MIXED DER. OF THE NA TERM W.R.T.
        !           DENSITY, TEMPERATURE AND TEMPERATURE
        !*****************************************************
        NADTT = gl%eos_coeff%ni(j,nrsubst)*(DELTACALC(5)*PSICALC(1) + 2.D0*DELTACALC(4)*PSICALC(4) + DELTACALC(1)*PSICALC(5) &
            & + del*(DELTACALC(7)*PSICALC(1) + DELTACALC(5)*PSICALC(2) + 2.D0*DELTACALC(6)*PSICALC(4) &
            & + 2.d0*DELTACALC(4)*PSICALC(6) + DELTACALC(2)*PSICALC(5) + DELTACALC(1)*PSICALC(7)))

        NACALC = NADTT

        !*****************************************************
    elseif ((DeltaDeriv == 3) .AND. (TauDeriv == 0)) then
        !*****************************************************
        ! NADD3     3RD DER. OF THE NA TERM W.R.T. DENSITY
        !*****************************************************
        NADDD = gl%eos_coeff%ni(j,nrsubst)*(DELTACALC(1)*(3.D0*PSICALC(3) + del*PSICALC(8)) + DELTACALC(2)*(6.d0*PSICALC(2) + 3.D0*del*PSICALC(3)) &
            & + DELTACALC(3)*(3.d0*PSICALC(1) + 3.D0*del*PSICALC(2)) + DELTACALC(8)*del*PSICALC(1))
        NACALC = NADDD

        !*****************************************************
    elseif ((DeltaDeriv == 0) .AND. (TauDeriv == 3)) then
        !*****************************************************
        ! NADT3     3RD DER. OF THE NA TERM W.R.T. TEMPERATURE
        !*****************************************************
        NATTT = gl%eos_coeff%ni(j,nrsubst)* del * (DELTACALC(9) * PSICALC(1) + 3.d0 * (DELTACALC(5) * PSICALC(4) + DELTACALC(4) * PSICALC(5))  &
            + DELTACALC(1) * PSICALC(9))
        NACALC = NATTT

        !*****************************************************
    elseif ((DeltaDeriv == 2) .AND. (TauDeriv == 1)) then
        !*****************************************************
        ! NADD2T    3RD DER. OF THE NA TERM W.R.T. DENSITY,
        !           DENSITY AND TEMPERATURE
        !*****************************************************
        NADDT = gl%eos_coeff%ni(j,nrsubst)* (DELTACALC(1) * (2.d0 * PSICALC(6) + del * PSICALC(10)) + DELTACALC(4) * (2.d0 * PSICALC(2) + del * PSICALC(3)) &
            & + 2.d0 * DELTACALC(2) * (PSICALC(4) + del * PSICALC(6)) + DELTACALC(3) * del * PSICALC(4) + 2.d0 * DELTACALC(6) * (PSICALC(1) &
            & + del * PSICALC(2)) + DELTACALC(10) * del * PSICALC(1))

        NACALC = NADDT
    end if

    end subroutine NATERMDERIVS

    end module


    module fnrderivs_base

    use module_all_types
    use special_terms
    
    implicit none

    contains

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !  Helmholtz derivatives with respect to delta^n
    !
    ! ##########################################################################################################################################################################################################################
    ! First derivative of helmholtz energy with respect to delta multiplied by delta
    ! ##########################################################################################################################################################################################################################
    function FNRD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)

    double precision , dimension(nderivs) :: FNRD_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,DELEPS
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,tau,SGBSCALC,NACALC
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRD_FUNC = 0.d0

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the regular terms
    do i = 1, gl%eos_coeff%nreg(nrsubst)
        FNRD_FUNC(1) = FNRD_FUNC(1) + reg_term(i) * (gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))
    end do

    !!Summation over the Gaussian bell-shaped terms
    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)
        if (dabs(DELEPS(k)) < 1.d-10) then   ! for del = 1 and epsylon = 1 the exponent becomes unity
            FNRD_FUNC(2) = FNRD_FUNC(2) + gauss_term(k) * gl%eos_coeff%di(i,nrsubst)
        else if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
            !FNRD = FNRD + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + NINT(gl%eos_coeff%pli(k,nrsubst))*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**NINT((gl%eos_coeff%pli(k,nrsubst)-1.D0)))
            FNRD_FUNC(2) = FNRD_FUNC(2) + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + 2.d0*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k))
        else
            FNRD_FUNC(2) = FNRD_FUNC(2) + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))
        end if
    end do

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        DelDeriv = 1
        TauDeriv = 0
        SGBSCALC = 0d0
        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRD_FUNC(3) = FNRD_FUNC(3) + SGBSCALC

        end do

    end if

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the nonanalytic terms
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1.d0),gl%eos_coeff%ncrt(nrsubst)
        DelDeriv = 1
        TauDeriv = 0
        NACALC = 0d0
        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRD_FUNC(4) = FNRD_FUNC(4) + NACALC*del
        end if
    end do

    end function
    ! ##########################################################################################################################################################################################################################


    ! ##########################################################################################################################################################################################################################
    ! Second derivative of helmholtz energy with respect to delta multiplied by delta^2
    ! ##########################################################################################################################################################################################################################
    function FNRDD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)

    double precision , dimension(nderivs) :: FNRDD_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,DELEPS
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,del2,tau,SGBSCALC,NACALC
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRDD_FUNC = 0.D0
    del2 = del*del
    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the regualar terms
    do i = gl%eos_coeff%nreg(nrsubst),1, -1
        FNRDD_FUNC(1) = FNRDD_FUNC(1) + reg_term(i) * ((gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))* &
            (gl%eos_coeff%di(i,nrsubst) - 1.d0 - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst)) - &
            gl%eos_coeff%gama(i,nrsubst)*(gl%eos_coeff%p_i(i,nrsubst)**2)*(delpi(i,nrsubst)))
    end do

    !!Summation over the Gaussian bell-shaped terms
    do i = gl%eos_coeff%nreg(nrsubst) + 1, (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)
        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
            FNRDD_FUNC(2) = FNRDD_FUNC(2) + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + 2.d0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k))**2 &
                - gl%eos_coeff%di(i,nrsubst) + 2.d0*gl%eos_coeff%eta(k,nrsubst)*del2)
        else
            FNRDD_FUNC(2) = FNRDD_FUNC(2) + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))**2 &
                - gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0)*(gl%eos_coeff%pli(k,nrsubst)-1.D0))
        end if
    end do

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        DelDeriv = 2
        TauDeriv = 0
        SGBSCALC = 0d0
        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRDD_FUNC(3) = FNRDD_FUNC(3) + SGBSCALC

        end do

    end if

    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the nonanalytic terms
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
        DelDeriv = 2
        TauDeriv = 0
        NACALC = 0d0
        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRDD_FUNC(4) = FNRDD_FUNC(4) + NACALC*del*del
        end if
    end do


    end function
    ! ##########################################################################################################################################################################################################################

    ! ##########################################################################################################################################################################################################################
    ! Third derivative of helmholtz energy with respect to delta multiplied by delta^3
    ! ##########################################################################################################################################################################################################################
    function FNRDDD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)

    double precision , dimension(nderivs) :: FNRDDD_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,DELEPS
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,del2,del3,tau,SGBSCALC,NACALC
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRDDD_FUNC = 0.D0
    del2 = del*del
    del3 = del2*del
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !!Summation over the regular terms
    do i = 1, gl%eos_coeff%nreg(nrsubst)
        FNRDDD_FUNC(1) = FNRDDD_FUNC(1) + reg_term(i) * (gl%eos_coeff%di(i,nrsubst)*(gl%eos_coeff%di(i,nrsubst) - 1.d0)*(gl%eos_coeff%di(i,nrsubst) - 2.d0) + &
            gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst)*( - 2.d0 + 6.d0*gl%eos_coeff%di(i,nrsubst) &
            - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 - 3.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst) + 3.d0*gl%eos_coeff%p_i(i,nrsubst) - gl%eos_coeff%p_i(i,nrsubst)**2) &
            + 3.d0*gl%eos_coeff%gama(i,nrsubst)**2*gl%eos_coeff%p_i(i,nrsubst)**2*del**(2.d0*gl%eos_coeff%p_i(i,nrsubst))*(gl%eos_coeff%di(i,nrsubst) - 1.d0 + gl%eos_coeff%p_i(i,nrsubst))&
            - gl%eos_coeff%gama(i,nrsubst)**3*gl%eos_coeff%p_i(i,nrsubst)**3*del**(3*gl%eos_coeff%p_i(i,nrsubst)))
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over gaussian terms
    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)
        if (abs(DELEPS(k)) < 1.d-10) then   ! for del = 1 and epsylon = 1 the exponent becomes unity
            FNRDDD_FUNC(2) = FNRDDD_FUNC(2) + gauss_term(k)*gl%eos_coeff%di(i,nrsubst)*(gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%di(i,nrsubst) - 3.d0*gl%eos_coeff%di(i,nrsubst) + 2.d0)
        else if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then
            FNRDDD_FUNC = FNRDDD_FUNC + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + 2.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k))**3 - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 &
                & + 2.d0*gl%eos_coeff%di(i,nrsubst) + 6.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2 - 6.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)*(gl%eos_coeff%di(i,nrsubst) &
                & - 2.D0*gl%eos_coeff%eta(k,nrsubst)*del2))
        else
            FNRDDD_FUNC(2) = FNRDDD_FUNC(2) + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 1.d0))**3 &
                & - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 + 2.d0*gl%eos_coeff%di(i,nrsubst) + 3.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2*gl%eos_coeff%pli(k,nrsubst)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)* &
                & DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0) - 3.d0*gl%eos_coeff%di(i,nrsubst)*del*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0) &
                & + 3.D0*del3*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)**2*gl%eos_coeff%eta(k,nrsubst)**2*DELEPS(k)**(2.d0*gl%eos_coeff%pli(k,nrsubst) - 3.d0) &
                & + (gl%eos_coeff%pli(k,nrsubst) - 2.d0)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 3.d0))
        end if
    end do
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !!Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        SGBSCALC = 0.D0
        DelDeriv = 3
        TauDeriv = 0

        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRDDD_FUNC(3) = FNRDDD_FUNC(3) + SGBSCALC

        end do

    end if
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !!Summation over the nonanalytic terms
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)

        NACALC = 0d0
        DelDeriv = 3
        TauDeriv = 0

        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRDDD_FUNC(4) = FNRDDD_FUNC(4) + NACALC*del*del*del
        end if
    end do
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    end function
    ! ##########################################################################################################################################################################################################################


    ! -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !  Helmholtz derivatives with respect to tau^n
    !
    ! ##########################################################################################################################################################################################################################
    ! First derivative of helmholtz energy with respect to tau multiplied by tau
    ! ##########################################################################################################################################################################################################################
    function FNRT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,region)

    double precision , dimension(nderivs) :: FNRT_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,TAUGAM
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,tau,SGBSCALC,NACALC
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRT_FUNC = 0.D0
    !!Summation over the regular terms
    do i = 1, gl%eos_coeff%nreg(nrsubst)
        FNRT_FUNC(1) = FNRT_FUNC(1) + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !!Summation over the Gaussian bell-shaped terms
    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)
        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0.and.gl%eos_coeff%tli(k,nrsubst) == 2.D0) then       !if the exponents in the Gauss exponents are equal to 2
            FNRT_FUNC(2) = FNRT_FUNC(2) + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + 2.d0*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))
        else
            if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                FNRT_FUNC(2) = FNRT_FUNC(2) + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + (gl%eos_coeff%tli(k,nrsubst))*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT((gl%eos_coeff%tli(k,nrsubst)-1.D0)))
            else
                FNRT_FUNC(2) = FNRT_FUNC(2) + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + (gl%eos_coeff%tli(k,nrsubst))*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**((gl%eos_coeff%tli(k,nrsubst)-1.D0)))
            end if
        end if
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        DelDeriv = 0
        TauDeriv = 1
        SGBSCALC = 0d0
        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRT_FUNC(3) = FNRT_FUNC(3) + SGBSCALC

        end do

    end if

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Summation over the nonanalytic terms
    !do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
        DelDeriv = 0
        TauDeriv = 1
        NACALC = 0d0
        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRT_FUNC(4) = FNRT_FUNC(4) + NACALC*tau
        end if
    end do

    end function
    ! ##########################################################################################################################################################################################################################


    ! ##########################################################################################################################################################################################################################
    ! Second derivative of helmholtz energy with respect to tau multiplied by tau^2
    ! ##########################################################################################################################################################################################################################
    function FNRTT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,DELEPS,region)

    double precision , dimension(nderivs) :: FNRTT_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,TAUGAM,DELEPS
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,tau,SGBSCALC,NACALC,tau2part,delpart,tau2
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRTT_FUNC= 0.D0
    tau2 = tau*tau
    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the regular terms
    do i = 1, gl%eos_coeff%nreg(nrsubst)
        FNRTT_FUNC(1) = FNRTT_FUNC(1) + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%ti(i,nrsubst) - 1.d0)*(gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the Gaussian bell-shaped terms
    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)
        if ((gl%eos_coeff%pli(k,nrsubst) == 2.d0) .and. (gl%eos_coeff%tli(k,nrsubst) == 2.d0)) then       !if the exponents in the Gauss exponents are equal to 2
            FNRTT_FUNC(2) = FNRTT_FUNC(2) + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + 2.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)) &
                *((gl%eos_coeff%ti(i,nrsubst) + 2.D0*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k))**2 - gl%eos_coeff%ti(i,nrsubst) + 2.D0*gl%eos_coeff%beta(k,nrsubst)*tau2)
        else
            if ((gl%eos_coeff%tli(k,nrsubst) - 2.d0) == 0.d0) then
                FNRTT_FUNC(2) = FNRTT_FUNC(2) + gauss_term(k) * &
                    (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.d0))**2*&
                    ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.d0))**2&
                    - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*(gl%eos_coeff%tli(k,nrsubst)-1.d0)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-2.d0))
            else

                !für spätere Struktur der gesamten Routine
                tau2part = ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))**2 &
                    - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-2.D0)*(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                delpart = (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))

                FNRTT_FUNC(2) = FNRTT_FUNC(2) + gauss_term(k)*tau2part*delpart

            end if
        end if
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        DelDeriv = 1.d0
        TauDeriv = 2.d0

        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRTT_FUNC(3) = FNRTT_FUNC(3) + SGBSCALC

        end do

    end if

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the nonanalytic terms
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
        DelDeriv = 1
        TauDeriv = 2
        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRTT_FUNC(4) = FNRTT_FUNC(4) + NACALC*del*tau*tau
        end if
    end do


    end function


    ! ##########################################################################################################################################################################################################################
    ! Third derivative of helmholtz energy with respect to tau multiplied by tau^3
    ! ##########################################################################################################################################################################################################################
    function FNRTTT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,region)

    double precision , dimension(nderivs) :: FNRTTT_FUNC
    type(type_gl) ::  gl
    double precision, dimension(:) :: reg_term,gauss_term,TAUGAM
    double precision, dimension(:,:) :: delpi
    integer ,dimension(:) :: region
    double precision :: del,tau,SGBSCALC,NACALC,tau_inv,tau_inv2,tau_inv3,tau3
    integer :: DelDeriv,TauDeriv
    integer :: nrsubst
    integer :: i,k,nn

    FNRTTT_FUNC= 0.D0


    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the regular terms
    do i = 1, gl%eos_coeff%nreg(nrsubst)
        FNRTTT_FUNC(1) = FNRTTT_FUNC(1) + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%ti(i,nrsubst) - 1.d0)*(gl%eos_coeff%ti(i,nrsubst) - 2.d0)
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the Gaussian bell-shaped terms
    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
        k = i - gl%eos_coeff%nreg(nrsubst)

        !if the exponents in the Gauss exponents are equal to 2
        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then

            tau_inv = 1.d0/tau
            tau_inv2 = tau_inv*tau_inv
            tau_inv3 = tau_inv2*tau_inv
            tau3 = tau*tau*tau

            FNRTTT_FUNC(2) = FNRTTT_FUNC(2) + gauss_term(k)*tau3*&
                (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))**2&
                -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))+2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

        else if ((gl%eos_coeff%tli(k,nrsubst) - 2.d0) == 0.d0) then

            FNRTTT_FUNC(2) = FNRTTT_FUNC(2) + gauss_term(k)*tau3*&
                (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))**NINT(gl%eos_coeff%tli(k,nrsubst))&
                -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))&
                +2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

        else

            FNRTTT_FUNC(2) = FNRTTT_FUNC(2) + gauss_term(k)*tau3*&
                (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))**gl%eos_coeff%tli(k,nrsubst)&
                -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))+2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

        end if
    end do

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the special GBS terms
    if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

        DelDeriv = 0.d0
        TauDeriv = 3.d0

        do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

            call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

            FNRTTT_FUNC(3) = FNRTTT_FUNC(3) + SGBSCALC

        end do

    end if

    ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Summation over the nonanalytic terms
    do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
        DelDeriv = 0
        TauDeriv = 3
        if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
            call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
            FNRTTT_FUNC(4) = FNRTTT_FUNC(4) + NACALC*tau*tau*tau
        end if
    end do
    end function

    end  module