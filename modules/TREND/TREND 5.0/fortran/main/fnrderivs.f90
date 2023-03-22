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
#include "access_macros.fi"




    ! module for file fnrderivs.f90
    module fnrderivs_module
    !global use inclusion
    use setup_module
    use module_all_types
    use flash_module
    use association_module
    use cubic_eos_module
    use hard_sphere_module
    use pc_saft_module
    use reduced_parameters_calc_module
    use mixtures_AGA8_module

    use fnrderivs_base
    !    use crossover

    !DEC$ IF DEFINED(COMPARE)
    logical, parameter :: enable_cache = .false.
    !DEC$ ELSE
    logical, parameter :: enable_cache = .true.
    !DEC$ END IF


    integer, dimension(5), parameter ::upper_triangle_offset = (/0,0,1,3,6/)

    contains

    !function to set dimension of nreg and delpi
    pure integer function dim_set(gl,nrsubst)
    implicit none
    integer, intent(in):: nrsubst
    type(type_gl), intent(in) :: gl
    if(allocated(gl%eos_coeff))  then
        dim_set = gl%eos_coeff%nreg(nrsubst)
    else
        dim_set = 0
    end if
    end function




    !**************************************************************************
    !           --------------------------------------------------
    !           Routine for the calculation of the derivatives of
    !           the residual part of the helmholz energy
    !
    !           J. Gernert, Denmark, 09.2009
    !           --------------------------------------------------
    !**************************************************************************
    !
    !**************************************************************************
    subroutine FNRDERIVS (gl,TEMPERATURE, DENSITY, GETDER, FNRDER, nrsubst)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY
    ! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
    ! AS PUBLISHED BY (PP. 24):
    !                ****************************************
    !                *   Span, R.                           *
    !                *   Multiparameter Equations of State  *
    !                *   Springer, 2000                     *
    !                ****************************************
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED RESIDUAL HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
    !                2. 1ST DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del
    !                3. 2ND DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del^2
    !                4. 1ST DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU^2
    !                6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO DEL AND TAU, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO DEL, TAU, AND TAU, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    !               12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    !               13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    !               14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    !               15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    ! FLUID:      - INTEGER: THE NUMBER OF THE FLUID THE HELMHOLTZ ENERGY IS CALCULATED FOR IN THE ORDER GIVEN
    !                BY THE USER INPUT
    !
    ! OUTPUT PARAMETERS:
    ! FNRDER      - AN ARRAY WITH 10 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !                AS INDECATED IN "GETDER"
    !--------------------------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl


    double precision:: Temperature, Density
    integer, dimension(nderivs):: GETDER                         ! array specifier to indicate, which derivative is needed
    integer :: nrsubst                                      ! fluid numbers for the fluid F is calculated (FLUID1) and for the fluid FLUID1 is combined with (FLUID2)
    double precision, dimension(nderivs)::FNRDER                 ! array with the computed values for the derivatives
    double precision,dimension(nderivs) :: FNRASSODER            ! Vector gives the values of the calculated derivatives (association term)
    double precision,dimension(nderivs):: HSDER                  ! Vector gives the values of the calculated derivatives (hard sphere term)
    !double precision, dimension (0:6):: delp_i                ! delta to the power of i and 0 for index 0
    !double precision:: logdel, logtau                      ! reduced temperature and density
    double precision:: del, tau, expo                       ! reduced temperature and density
    double precision, dimension (dim_set(gl,nrsubst)) :: reg_term             ! array of the regular terms
    double precision, dimension (20):: gauss_term, TAUGAM, DELEPS      ! array of the gaussian bell-shaped terms
    double precision, dimension (dim_set(gl,nrsubst),gl%ncomp) :: DELPI             ! delta^pi in exponential terms
    double precision:: FNR, FNRD, FNRDD, FNRDDD, FNRDDT     ! residual helmholtz energy and its derivatives
    double precision:: FNRT, FNRTT, FNRDT, FNRDTT, FNRTTT   ! with respect to density and temperature
    double precision:: FNRDDDT,FNRDDTT,FNRDTTT
    double precision:: NACALC                               ! Return parameter for the calculation routine of the non-analytic term
    double precision:: SGBSCALC                             ! Return parameter for the calculation routine of the special GBS term
    !integer:: I_pol, I_exp(nrsubst), I_GBS(nrsubst), I_NA(nrsubst)                     ! number of polynomial, exponential, gaussian bell-shaped, and nonanalytic terms
    integer:: i, j, k, nn, o                                ! index variables
    integer:: DelDeriv, TauDeriv                            ! index variables for the calculation of the nonanalytic terms
    integer, dimension(3) :: region                         ! flag for NA-term-calc.: 1 - crit. region, NA-term is calculated, 0 - outside crit. region
    double precision:: diff_SM, NACALCNUM, delpart,taupart, del2part,tau2part, Rmix

    !SRK supportive variables. Lars Hüttermann, 14.02.2013
    double precision :: dT_dTau
    double precision :: d2T_dTau2
    double precision :: d3T_dTau3
    double precision :: da_dTau
    double precision :: d2a_dTau2
    double precision :: d3a_dTau3

    !PR supportive variables. Stefan, Feb 2014

    !LKP supportive variable. Stefan, May 2014

    double precision :: B_0
    double precision :: C_0
    double precision :: D_0
    double precision :: B_ref
    double precision :: C_ref
    double precision :: D_ref
    double precision :: part_id
    double precision :: part_ref
    double precision :: help_id1, help_id2
    double precision :: help_ref1, help_ref2
    double precision :: zcLKP2, zcLKP4, zcLKP5,zcLKP6, zcLKP8, tau2, tau3, del2, del3, del4, del5
    double precision :: DB_0_Dtau, DC_0_Dtau, DD_0_Dtau, DB_ref_Dtau, DC_ref_Dtau, DD_ref_Dtau, D2B_0_Dtau2, D2C_0_Dtau2, D2B_ref_Dtau2, D2C_ref_Dtau2
    double precision :: D3B_0_Dtau3, D3C_0_Dtau3, D3B_ref_Dtau3, D3C_ref_Dtau3
    double precision :: del_inv, del_inv2, tau_inv, tau_inv2, tau_inv3, tau4
    double precision :: srk_t, srk_t2, srk_t3, srk_t4, rhoc_del, bi_rhoc_del

    ! variables for calulation of 4th mixed derivatives
    double precision , dimension(nderivs) ::  FNRD_F,FNRDDD_F,FNRT_F,FNRTT_F,FNRDD_F,FNRDTTT_F,FNRTTT_F


    !SAFT
    double precision:: dens_orig, Temp_orig


    !saving derivatives
    logical :: complete
    integer, dimension(nderivs):: GETDERS             ! working array specifier to indicate, which derivative is still needed

    !numerical integration
    double precision:: integral

    integer:: nreg

    integral=0.d0

    SGBSCALC = 0.d0
    NACALC = 0.d0

    dT_dTau = 0.d0
    d2T_dTau2 = 0.d0
    d3T_dTau3 = 0.d0
    da_dTau = 0.d0
    d2a_dTau2 = 0.d0
    d3a_dTau3 = 0.d0

    tau = 0.d0
    tau2 = 0.d0
    tau3 = 0.d0
    tau4 = 0.d0
    tau_inv = 0.d0
    tau_inv2 = 0.d0
    tau_inv3 = 0.d0


    del = 0.d0
    del2 = 0.d0
    del3 = 0.d0
    del4 = 0.d0
    del5 = 0.d0
    del_inv = 0.d0
    del_inv2 = 0.d0

    srk_t = 0.d0
    srk_t2 = 0.d0
    srk_t3 = 0.d0
    srk_t4 = 0.d0
    rhoc_del = 0.d0
    bi_rhoc_del = 0.d0


    !set all array elements to default value 0
    FNR = 0.d0


    region = (/1, 1, 1/)

    FNRDER = 0.d0        ! initialize the return vector
    FNRASSODER = 0.d0    ! initialize the return vector
    HSDER = 0.d0         ! initialize the return vector


    GETDERS = GETDER
    !***************************************************************************************************************
    complete = .true.
    ! array specifier to indicate, which derivative is needed


    !if ( (dabs(TEMPERATURE-gl%TEMP_FNR_OLD) < 1.d-14) .and. (dabs(DENSITY-gl%DENS_FNR_OLD) < 1.d-14) &
    !        .and. (dabs(gl%tredmix-gl%tredmix_old) < 1.d-14) .and. (dabs(gl%rhoredmix-gl%rhoredmix_OLD) < 1.d-14) &
    !        .and. (nrsubst == gl%nrsubst_old) .and. (gl%Eq_type(nrsubst) == gl%Eq_type_old(nrsubst)) ) then
    !
    !        FNR = gl%FNR_OLD
    !        delpi = gl%delpi_old
    !        reg_term = gl%reg_term_old
    !        gauss_term = gl%gauss_term_old
    !        nacalc = gl%nacalc_old
    !        deleps = gl%deleps_old
    !        taugam = gl%taugam_old
    !
    !        do o = 1, 10
    !            if ( (GETDERS(o) == 1) .and. (dabs(gl%FNRDER_OLD(o)) > 1.d-14) ) then
    !                FNRDER(o) = gl%FNRDER_OLD (o)
    !                GETDERS(o) = 0
    !            else if (GETDERS(o) == 1) then
    !                complete = .false.
    !            end if
    !        end do
    !
    !    if (complete) then
    !        gl%zaehler = gl%zaehler +1
    !        return
    !    end if
    !
    !else
    !    same = .false.
    !end if

    !***************************************************************************************************************



    if (gl%NCOMP > 1) then
        if ((gl%bwr_read) .and. (gl%Req(nrsubst) == 1.d0) .and. (.not.gl%vir)) then      !when calculating Lennard-Jones-Fluid with MBWR
            del = Density
            tau = 1.d0/Temperature
        else
            del = Density / gl%rhoredmix      ! reduced density calculated with reducing function for mixtures
            tau = gl%tredmix / Temperature    ! inverse reduced temperature calculated with reducing function for mixtures
        endif
    else
        if ((gl%bwr_read) .and. (gl%Req(nrsubst) == 1.d0) .and. (.not.gl%vir)) then      !when calculating Lennard-Jones-Fluid with MBWR
            del = Density
            tau = 1.d0/Temperature
        else
            del = Density / gl%rhored(nrsubst)     ! reduced density for single fluid
            tau = gl%tred(nrsubst) / Temperature   ! inverse reduced temperature for single fluid
        endif
    end if


    If ((gl%Eq_type(nrsubst) == 1) .or. (gl%Eq_type(nrsubst) == 51) .or. (gl%Eq_type(nrsubst) == 52) .or. (gl%Eq_type(nrsubst) == 53)) then!EOS explicit in the Helmholtz energy is used

        ! code duplicated into both "if (del > 1.D-14) then"-branches for better cache-performance

        !!if (.not.allocated(reg_term))allocate(reg_term(gl%eos_coeff%nreg(nrsubst)))
        !reg_term = 0.D0
        !gauss_term = 0.D0
        !del2 = del*del
        !del3 = del2*del
        !tau2 = tau*tau
        !tau3 = tau2*tau
        !
        !
        !!calculate assoterms if exist for the fluid:
        !IF (gl%ASSOEXIST(nrsubst)) THEN
        !    CALL FNRASSODERIVS(gl,TEMPERATURE,DENSITY,GETDERS,FNRASSODER,nrsubst)
        !END IF
        !
        !!calculate the hard sphere term if exist for the fluid
        !IF (gl%eos_coeff%hard_sphere(nrsubst)) THEN
        !    CALL HSDERIVS (gl,TEMPERATURE, DENSITY, GETDERS, HSDER, nrsubst)
        !END IF
        !
        !! computes del and tau from reducing parameters
        !!! find the number of polynomial terms in the equation
        !!do l = 1, nreg(nrsubst)
        !!    if (p_i(nrsubst, l) > 0.D0) exit
        !!end do
        !!! find the number of gaussian bell-shaped terms in the equation
        !!do m = nreg(nrsubst) + 1, (nreg(nrsubst) + ncrt(nrsubst))
        !!    if (abs(etana(nrsubst, m)) > 0.D0) exit
        !!end do
        !!I_pol = l - 1               ! number of polynomial terms
        !!I_exp = nreg(nrsubst) - I_pol    ! number of exponential terms
        !!I_GBS(nrsubst) = m - 1 - nreg(nrsubst)    ! number of gaussian bell-shaped terms
        !!I_NA(nrsubst) = ncrt(nrsubst) - I_GBS(nrsubst)     ! number of nonanalytic terms
        !
        !!*****************************************************
        !!   FNR         RESIDUAL HELMHOLTZ FREE ENERGY F [-]
        !!*****************************************************
        !
        !! calculate delta^i (saves computing time)
        !delp_i(0) = 0.d0
        !delp_i(1)=del
        !do j =2, 6
        !    delp_i(j) = delp_i(j-1)*del
        !end do
        !
        !!if (.not. allocated(delpi)) then
        !!    allocate(delpi(gl%eos_coeff%nreg(nrsubst),gl%ncomp))
        !!end if
        !delpi = 0.d0
        !do i = 1, gl%eos_coeff%nreg(nrsubst)
        !    delpi(i,nrsubst) = delp_i(gl%eos_coeff%p_i(i,nrsubst))
        !    !            if (gl%eos_coeff%p_i(i,nrsubst) == 1) then
        !    !                delpi(i,nrsubst) = delp_i(1)
        !    !            elseif (gl%eos_coeff%p_i(i,nrsubst) == 2) then
        !    !                delpi(i,nrsubst) = delp_i(2)
        !    !            elseif (gl%eos_coeff%p_i(i,nrsubst) == 3) then
        !    !                delpi(i,nrsubst) = delp_i(3)
        !    !            elseif (gl%eos_coeff%p_i(i,nrsubst) == 4) then
        !    !                delpi(i,nrsubst) = delp_i(4)
        !    !            elseif (gl%eos_coeff%p_i(i,nrsubst) == 5) then
        !    !                delpi(i,nrsubst) = delp_i(5)
        !    !            elseif (gl%eos_coeff%p_i(i,nrsubst) == 6) then
        !    !                delpi(i,nrsubst) = delp_i(6)
        !    !            end if
        !end do

        !Andreas, September 2015
        !For the calculation of virial coefficients, the density needs to be set to 0. Elsewise, the virials could be inaccurate due to rounding errors
        !IMPORTANT!!!:  Furthermore, as the zero density limit of the delta-derivatives of alpha_r are needed for the calculation of the virial coefficients,
        !   the derivatives must not be multiplied by delta, delta², delta3!!!
        !ATTENTION!!    --> If the density is 0, then the derivatives are not longer multiplied by delta, delta² and so on!!!
        !Unfortunately, the system of calculating the derivatives in a clever manner cannot be applied then anymore. Thus, all derivatives are programmed in a different way for
        !the evaluation of the virial coefficients (see the following else-statement)
        !SH: Changed threshold from 1.d-12 to 1.d-14
#IF DEFINED(CON_FIT)
        if(del > 1.D-14 .or. gl%virial_num) then
#ELSE
            if (del > 1.D-14) then
#ENDIF
                ! check cache
                if (enable_cache .and. del == gl%fnrderivs_cache%del_old) then
                    continue ! nothing to do
                else
                    gl%fnrderivs_cache%del_old = del

                    gl%fnrderivs_cache%logdel = dlog(del)
                    ! calculate delta^i (saves computing time)
                    gl%fnrderivs_cache%delp_i(0) = 0.d0
                    gl%fnrderivs_cache%delp_i(1)=del
                    do j =2, 6
                        gl%fnrderivs_cache%delp_i(j) = gl%fnrderivs_cache%delp_i(j-1)*del
                    end do
                endif

                del2 = gl%fnrderivs_cache%delp_i(2)
                del3 = gl%fnrderivs_cache%delp_i(3)



                if (enable_cache .and. tau == gl%fnrderivs_cache%tau_old) then
                    continue ! nothing to do
                else
                    gl%fnrderivs_cache%tau_old = tau

                    gl%fnrderivs_cache%logtau = dlog(tau)
                    ! calculate tauta^i (saves computing time)
                    gl%fnrderivs_cache%tau2 = tau*tau
                    gl%fnrderivs_cache%tau3 = gl%fnrderivs_cache%tau2*tau
                endif

                tau2 = gl%fnrderivs_cache%tau2
                tau3 = gl%fnrderivs_cache%tau3




                !if (.not.allocated(reg_term))allocate(reg_term(gl%eos_coeff%nreg(nrsubst)))
                reg_term = 0.D0
                gauss_term = 0.D0


                !calculate assoterms if exist for the fluid:
                IF (gl%ASSOEXIST(nrsubst)) THEN
                    CALL FNRASSODERIVS(gl,TEMPERATURE,DENSITY,GETDERS,FNRASSODER,nrsubst)
                END IF

                !calculate the hard sphere term if exist for the fluid
                IF (gl%eos_coeff%hard_sphere(nrsubst)) THEN
                    CALL HSDERIVS (gl,TEMPERATURE, DENSITY, GETDERS, HSDER, nrsubst)
                END IF

                ! computes del and tau from reducing parameters
                !! find the number of polynomial terms in the equation
                !do l = 1, nreg(nrsubst)
                !    if (p_i(nrsubst, l) > 0.D0) exit
                !end do
                !! find the number of gaussian bell-shaped terms in the equation
                !do m = nreg(nrsubst) + 1, (nreg(nrsubst) + ncrt(nrsubst))
                !    if (abs(etana(nrsubst, m)) > 0.D0) exit
                !end do
                !I_pol = l - 1               ! number of polynomial terms
                !I_exp = nreg(nrsubst) - I_pol    ! number of exponential terms
                !I_GBS(nrsubst) = m - 1 - nreg(nrsubst)    ! number of gaussian bell-shaped terms
                !I_NA(nrsubst) = ncrt(nrsubst) - I_GBS(nrsubst)     ! number of nonanalytic terms

                !*****************************************************
                !   FNR         RESIDUAL HELMHOLTZ FREE ENERGY F [-]
                !*****************************************************

                !if (.not. allocated(delpi)) then
                !    allocate(delpi(gl%eos_coeff%nreg(nrsubst),gl%ncomp))
                !end if
                !delpi = 0.d0
                !do i = 1, gl%eos_coeff%nreg(nrsubst)
                !    delpi(i,nrsubst) = gl%fnrderivs_cache%delp_i(gl%eos_coeff%p_i(i,nrsubst))
                !end do





                !!Summation over the regular terms
                nreg = gl%eos_coeff%nreg(nrsubst)

                delpi(1:nreg,nrsubst) = gl%fnrderivs_cache%delp_i(gl%eos_coeff%p_i(1:nreg,nrsubst))
                reg_term(1:nreg) = gl%eos_coeff%ni(1:nreg,nrsubst) * dexp(gl%fnrderivs_cache%logdel*gl%eos_coeff%di(1:nreg,nrsubst) + gl%fnrderivs_cache%logtau*gl%eos_coeff%ti(1:nreg,nrsubst) - gl%eos_coeff%gama(1:nreg,nrsubst)*delpi(1:nreg,nrsubst))
                !do i = 1, gl%eos_coeff%nreg(nrsubst)
                !    reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * dexp(logdel*gl%eos_coeff%di(i,nrsubst) + logtau*gl%eos_coeff%ti(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*delpi(i,nrsubst))
                !    !  reg_term(:) = gl%eos_coeff%ni(:,nrsubst) * dexp(logdel*gl%eos_coeff%di(:,nrsubst) + logtau*gl%eos_coeff%ti(:,nrsubst) - gl%eos_coeff%gama(:,nrsubst)*delpi(:,nrsubst))
                !    FNR = FNR + reg_term(i)
                FNR = sum(reg_term)
                !end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                    TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                    if (gl%eos_coeff%pli(k,nrsubst) == 2.d0.and.gl%eos_coeff%tli(k,nrsubst) == 2.D0) then       !if the exponents in the Gauss exponents are equal to 2
                        expo = gl%fnrderivs_cache%logdel*gl%eos_coeff%di(i,nrsubst) + gl%fnrderivs_cache%logtau*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)
                    Else                                              !if the exponents in the Gauss exponents are not equal to 2
                        if(dabs(gl%eos_coeff%tli(k,nrsubst)-2.d0)<1.d-14) then
                            expo = gl%fnrderivs_cache%logdel*gl%eos_coeff%di(i,nrsubst) + gl%fnrderivs_cache%logtau*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**2
                        else
                            expo = gl%fnrderivs_cache%logdel*gl%eos_coeff%di(i,nrsubst) + gl%fnrderivs_cache%logtau*gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst)
                        end if
                    end if
                    gauss_term(k) = 0
                    if (expo < 200) gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * dexp (expo)   !prevents double precision overflow
                    FNR = FNR + gauss_term(k)
                end do

            else

                ! check cache
                if (enable_cache .and. del == gl%fnrderivs_cache%del_old) then
                    continue ! nothing to do
                else
                    gl%fnrderivs_cache%del_old = del

                    ! log not needed in this branch
                    !                gl%fnrderivs_cache%logdel = dlog(del)

                    ! calculate delta^i (saves computing time)
                    gl%fnrderivs_cache%delp_i(0) = 0.d0
                    gl%fnrderivs_cache%delp_i(1)=del
                    do j =2, 6
                        gl%fnrderivs_cache%delp_i(j) = gl%fnrderivs_cache%delp_i(j-1)*del
                    end do
                endif

                del2 = gl%fnrderivs_cache%delp_i(2)
                del3 = gl%fnrderivs_cache%delp_i(3)



                if (enable_cache .and. tau == gl%fnrderivs_cache%tau_old) then
                    continue ! nothing to do
                else
                    gl%fnrderivs_cache%tau_old = tau

                    ! log not needed in this branch
                    !                gl%fnrderivs_cache%logtau = dlog(tau)
                    ! calculate tauta^i (saves computing time)
                    gl%fnrderivs_cache%tau2 = tau*tau
                    gl%fnrderivs_cache%tau3 = gl%fnrderivs_cache%tau2*tau
                endif

                tau2 = gl%fnrderivs_cache%tau2
                tau3 = gl%fnrderivs_cache%tau3


                !if (.not.allocated(reg_term))allocate(reg_term(gl%eos_coeff%nreg(nrsubst)))
                reg_term = 0.D0
                gauss_term = 0.D0


                !calculate assoterms if exist for the fluid:
                IF (gl%ASSOEXIST(nrsubst)) THEN
                    CALL FNRASSODERIVS(gl,TEMPERATURE,DENSITY,GETDERS,FNRASSODER,nrsubst)
                END IF

                !calculate the hard sphere term if exist for the fluid
                IF (gl%eos_coeff%hard_sphere(nrsubst)) THEN
                    CALL HSDERIVS (gl,TEMPERATURE, DENSITY, GETDERS, HSDER, nrsubst)
                END IF

                ! computes del and tau from reducing parameters
                !! find the number of polynomial terms in the equation
                !do l = 1, nreg(nrsubst)
                !    if (p_i(nrsubst, l) > 0.D0) exit
                !end do
                !! find the number of gaussian bell-shaped terms in the equation
                !do m = nreg(nrsubst) + 1, (nreg(nrsubst) + ncrt(nrsubst))
                !    if (abs(etana(nrsubst, m)) > 0.D0) exit
                !end do
                !I_pol = l - 1               ! number of polynomial terms
                !I_exp = nreg(nrsubst) - I_pol    ! number of exponential terms
                !I_GBS(nrsubst) = m - 1 - nreg(nrsubst)    ! number of gaussian bell-shaped terms
                !I_NA(nrsubst) = ncrt(nrsubst) - I_GBS(nrsubst)     ! number of nonanalytic terms

                !*****************************************************
                !   FNR         RESIDUAL HELMHOLTZ FREE ENERGY F [-]
                !*****************************************************

                ! calculate delta^i (saves computing time)
                !delp_i(0) = 0.d0
                !delp_i(1)=del
                !do j =2, 6
                !    delp_i(j) = delp_i(j-1)*del
                !end do

                !if (.not. allocated(delpi)) then
                !    allocate(delpi(gl%eos_coeff%nreg(nrsubst),gl%ncomp))
                !end if
                delpi = 0.d0
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    delpi(i,nrsubst) = gl%fnrderivs_cache%delp_i(gl%eos_coeff%p_i(i,nrsubst))
                end do

                !IMPORTANT NOTE: Association terms, and not analytical terms not yet included!!! Andreas, September 2015

                if (GETDERS(1) == 1) then

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)
                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then
                            reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * del**gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(-gl%eos_coeff%gama(i,nrsubst)*delpi(i,nrsubst))
                        else
                            reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * del**gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst)
                        end if
                        FNR = FNR + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !!Summation over the Gaussian bell-shaped terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                            gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * del**gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                        Else                                              !if the exponents in the Gauss exponents are not equal to 2
                            if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * del**gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)))
                            else
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * del**gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst))
                            end if
                        end if

                        FNR = FNR + gauss_term(k)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    FNRDER(1) = FNR + HSDER(1)

                end if

                if (GETDERS(2) == 1) then
                    !*****************************************************
                    !    FNRD       1ST DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               del
                    !*****************************************************
                    FNRD = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)

                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 1 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3)  di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   1a) di + pi - 1 < 0 --> the exponential part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2a) di + pi - 1 = 0 --> delta vanishes, exponential term contributes to B
                            !                   3a) di + pi - 1 > 0 --> delta left --> exponential term does not contribute to B, since delta = 0! --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) -1.D0) == 0 .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * (-gl%eos_coeff%gama(i,nrsubst)* gl%eos_coeff%p_i(i,nrsubst)) * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 1 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3) di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if

                        end if
                        FNRD = FNRD + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !!Summation over the Gaussian bell-shaped terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        !Different cases:   1) di - 1 < 0 --> di = 0
                        !                   2) di - 1 = 0 --> delta vanishes, term contributes to B
                        !                   3) di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED

                        !case 1)
                        if(NINT(gl%eos_coeff%di(i,nrsubst)) == 0) then  !Monika 03/2018
                            if ((gl%eos_coeff%pli(k,nrsubst) - 2.d0) .lt. 1.d-14 .and. (gl%eos_coeff%tli(k,nrsubst) - 2.D0) .lt. 1.d-14) then       !if the exponents in the Gauss exponents are equal to 2
                                gauss_term(k) = 2.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                            else
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%pli(k,nrsubst) * gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst))
                            end if
                        end if

                        !case 2)
                        if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then  !modified by Monika 03/2018
                            if (gl%eos_coeff%pli(k,nrsubst) == 2.d0.and.gl%eos_coeff%tli(k,nrsubst) == 2.D0) then       !if the exponents in the Gauss exponents are equal to 2
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                            Else                                              !if the exponents in the Gauss exponents are not equal to 2
                                !these two cases are still wrong!!
                                if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                                    gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)))
                                else
                                    gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst))
                                end if
                            end if
                        end if

                        FNRD = FNRD + gauss_term(k)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    FNRDER(2) = FNRD + HSDER(2)
                end if

                if (GETDERS(3) == 1) then
                    !*****************************************************
                    !    FNRDD      2ND DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               del
                    !*****************************************************
                    FNRDD = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)
                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 2 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 2 = 0 --> delta vanishes, term contributes to B
                            !                   3)  di - 2 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   More different cases for the exponential terms required. See documentation

                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then   !that means di=1, p_i doesn't matter
                                !too complicated: reg_term(i) = ni(i,nrsubst) * di(i,nrsubst) * (di(i,nrsubst) - 1.D0) * tau**ti(i,nrsubst)
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * 2.d0 * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) - 2.D0) == 0 .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then   !that means di=1 and p_i=1
                                !too complicated, probably wrong: reg_term(i) = ni(i,nrsubst) * (-(2.D0 * di(i,nrsubst) - 1.D0 + gama(i,nrsubst)* p_i(i,nrsubst))) * gama(i,nrsubst)* p_i(i,nrsubst) * tau**ti(i,nrsubst)
                                reg_term(i) = -2.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + 2.D0 * gl%eos_coeff%p_i(i,nrsubst) - 2.D0) == 0 .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then     !that means di=0 and p_i=1
                                !too complicated: reg_term(i) = ni(i,nrsubst) * gama(i,nrsubst)**2 * p_i(i,nrsubst)**2 * tau**ti(i,nrsubst)
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)**2 * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 2 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 2 = 0 --> delta vanishes, term contributes to B
                            !                   3) di - 2 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                !too complicated: reg_term(i) = ni(i,nrsubst) * di(i,nrsubst) * (di(i,nrsubst) - 1.D0) * tau**ti(i,nrsubst)
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * 2.d0 * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if

                        end if
                        FNRDD = FNRDD + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !Gaussian bell-shaped terms
                    !modified by Monika and Theresa, Nov 2016

                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        !Different cases:   1) di - 2 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                        !                   2) di - 2 = 0 --> delta vanishes, term contributes to C
                        !                   3) di - 2 > 0 --> delta to the power di-1 left --> does not contribute to C, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED

                        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                gauss_term(k) = 2.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                                FNRDD = FNRDD + gauss_term(k)
                            elseif (NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) * (-4.d0) * gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst)
                                FNRDD = FNRDD + gauss_term(k)
                            end if
                        end if
                        !Attention: Fluids with exponents in the exponential part, which are not 2, are not considered here
                        !e.g., for R-125 (Lemmon and Jacobsen, 2005)
                    end do

                    FNRDER(3) = FNRDD + HSDER(3)
                end if

                if (GETDERS(4) == 1) then
                    !*****************************************************
                    !    FNRT       1ST DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau
                    !*****************************************************
                    !All derivatives of alpha_r with respect to ONLY tau
                    !are zero at zero density!!
                    FNRT = 0.D0
                    FNRDER(4) = FNRT + HSDER(4)
                end if

                if (GETDERS(5) == 1) then
                    !*****************************************************
                    !    FNRTT      2ND DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau
                    !*****************************************************
                    !All derivatives of alpha_r with respect to ONLY tau
                    !are zero at zero density!!
                    FNRTT = 0.D0
                    FNRDER(5) = FNRTT + HSDER(5)
                end if

                if (GETDERS(6) == 1) then
                    !*****************************************************
                    !    FNRDT      2ST DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau and del
                    !*****************************************************
                    FNRDT = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)

                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 1 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3)  di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   1a) di + pi - 1 < 0 --> the exponential part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2a) di + pi - 1 = 0 --> delta vanishes, exponential term contributes to B
                            !                   3a) di + pi - 1 > 0 --> delta left --> exponential term does not contribute to B, since delta = 0! --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst)* tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0)
                            end if
                            !if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) -1.D0) == 0) then
                            !    reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * (-gl%eos_coeff%gama(i,nrsubst)* gl%eos_coeff%p_i(i,nrsubst)) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0)
                            !end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 1 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3) di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0)
                            end if

                        end if
                        FNRDT = FNRDT + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !!Summation over the Gaussian bell-shaped terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        !Different cases:   1) di - 1 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                        !                   2) di - 1 = 0 --> delta vanishes, term contributes to B
                        !                   3) di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                        if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then

                            if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                                !gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                                !S.Pohl: Corrected 2019
                                gauss_term(k) =   dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) &
                                    & * gl%eos_coeff%ni(i,nrsubst)*tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0)*(-gl%eos_coeff%di(i,nrsubst)) &
                                    & * (2*gl%eos_coeff%beta(k,nrsubst)*tau*(gl%eos_coeff%gam(k,nrsubst)-tau) -gl%eos_coeff%ti(i,nrsubst))

                            Else                                              !if the exponents in the Gauss exponents are not equal to 2
                                if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                                    gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)))
                                    !gauss_term(k) = ( gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*tau)*tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0)* gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%gam(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*gl%eos_coeff%ni(i,nrsubst) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)))
                                else
                                    gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst)-1.D0) * dexp (gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**gl%eos_coeff%pli(k,nrsubst) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**gl%eos_coeff%tli(k,nrsubst))
                                end if
                            end if

                        end if

                        FNRDT = FNRDT + gauss_term(k)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    FNRDER(6) = FNRDT + HSDER(6)
                end if

                if (GETDERS(7) == 1) then
                    !*****************************************************
                    !    FNRDTT     3RD DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau, tau, and del
                    !*****************************************************
                    FNRDTT = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)

                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 1 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3)  di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   1a) di + pi - 1 < 0 --> the exponential part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2a) di + pi - 1 = 0 --> delta vanishes, exponential term contributes to B
                            !                   3a) di + pi - 1 > 0 --> delta left --> exponential term does not contribute to B, since delta = 0! --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * (gl%eos_coeff%ti(i,nrsubst)-1.D0) * tau**(gl%eos_coeff%ti(i,nrsubst)-2.D0)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) -1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * (-gl%eos_coeff%gama(i,nrsubst)* gl%eos_coeff%p_i(i,nrsubst)) * gl%eos_coeff%ti(i,nrsubst) * (gl%eos_coeff%ti(i,nrsubst)-1.D0) * tau**(gl%eos_coeff%ti(i,nrsubst)-2.D0)
                            end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 1 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 1 = 0 --> delta vanishes, term contributes to B
                            !                   3) di - 1 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * (gl%eos_coeff%ti(i,nrsubst)-1.D0) * tau**(gl%eos_coeff%ti(i,nrsubst)-2.D0)
                            end if

                        end if
                        FNRDTT = FNRDTT + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !Neglect Gaussian bell shaped terms! --> Check and discuss if valid.

                    FNRDER(7) = FNRDTT + HSDER(7)
                end if

                if (GETDERS(8) == 1) then
                    !*****************************************************
                    !    FNRDDD     3RD DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               del
                    !*****************************************************
                    FNRDDD = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)
                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 3 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 3 = 0 --> delta vanishes, term contributes to D
                            !                   3)  di - 3 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   More different cases for the exponential terms required. See documentation

                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 3.D0) == 0) then
                                !reg_term(i) = ni(i,nrsubst) * di(i,nrsubst) * (di(i,nrsubst) - 1.D0) * (di(i,nrsubst) - 2.D0) * tau**ti(i,nrsubst)
                                reg_term(i) = 6.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) - 3.D0) == 0  .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then   !di = 2, pi = 1  or di = 1, pi = 2
                                !reg_term(i) = ni(i,nrsubst) * (-((2.D0 * di(i,nrsubst) - 1.D0 + gama(i,nrsubst)* p_i(i,nrsubst)) * (di(i,nrsubst) + p_i(i,nrsubst) -2.D0) + di(i,nrsubst) * (di(i,nrsubst) - 1.D0)))* gama(i,nrsubst)* p_i(i,nrsubst) * tau**ti(i,nrsubst)
                                reg_term(i) = -6.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + 2.D0 * gl%eos_coeff%p_i(i,nrsubst) - 3.D0) == 0  .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then   !di = 1, pi = 1
                                !reg_term(i) = ni(i,nrsubst) * ((di(i,nrsubst) + 2.D0 * p_i(i,nrsubst)) + (2.D0 * di(i,nrsubst) - 1.D0 + gama(i,nrsubst) * p_i(i,nrsubst)) ) * gama(i,nrsubst)**2 * p_i(i,nrsubst)**2 * tau**ti(i,nrsubst)
                                reg_term(i) = 3.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)**2
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + 3.D0 * gl%eos_coeff%p_i(i,nrsubst) - 3.D0) == 0  .and. gl%eos_coeff%p_i(i,nrsubst) .ne. 0) then    !di = 0, pi = 1
                                !reg_term(i) = ni(i,nrsubst) * gama(i,nrsubst)**3 * p_i(i,nrsubst)**3 * tau**ti(i,nrsubst)
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * (-gl%eos_coeff%gama(i,nrsubst)**3)
                            end if
                            !if((NINT(gl%eos_coeff%di(i,nrsubst)) == 0) .and. (NINT(gl%eos_coeff%p_i(i,nrsubst)) == 2)) then    !di = 0, pi = 2 !TN The derivative has to vanich for this case.
                            !    reg_term(i) = -4.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)
                            !end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 3 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 3 = 0 --> delta vanishes, term contributes to D
                            !                   3) di - 3 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 3.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * (gl%eos_coeff%di(i,nrsubst) - 1.D0) * (gl%eos_coeff%di(i,nrsubst) - 2.D0) * tau**gl%eos_coeff%ti(i,nrsubst)
                            end if

                        end if
                        FNRDDD = FNRDDD + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------

                    !Gaussian bell-shaped terms
                    !modified by Monika and Theresa, Nov 2016

                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        !Different cases:   1) di - 3 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                        !                   2) di - 3 = 0 --> delta vanishes, term contributes to D
                        !                   3) di - 3 > 0 --> delta to the power di-1 left --> does not contribute to D, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED

                        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 3.D0) == 0) then
                                gauss_term(k) = 6.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k))
                                FNRDDD = FNRDDD + gauss_term(k)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                gauss_term(k) = -12.d0 * gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) * gl%eos_coeff%eta(k,nrsubst)* gl%eos_coeff%eps(k,nrsubst)
                                FNRDDD = FNRDDD + gauss_term(k)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 1.D0) == 0) then
                                !Fehler bei Eric !? Nein, Fehler bei uns. Korrigert von RB/TN
                                !alt, mit Fehler !gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) * (24.d0 * (gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst))**2 - 12 * gl%eos_coeff%eta(k,nrsubst))
                                gauss_term(k) = gl%eos_coeff%ni(i,nrsubst) * tau**gl%eos_coeff%ti(i,nrsubst) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) * (2.d0 * gl%eos_coeff%eta(k,nrsubst) * (6.d0 * gl%eos_coeff%eta(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) * gl%eos_coeff%eps(k,nrsubst) + 3.d0))
                                FNRDDD = FNRDDD + gauss_term(k)
                            end if
                            !if(NINT(di(i,nrsubst) - 0.D0) == 0) then    !Maybe this case does not exist, check this
                            !        gauss_term(k) = -72.d0 * ni(i,nrsubst) * tau**ti(i,nrsubst) * dexp(eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) * eta(k,nrsubst)**2 * eps(k,nrsubst)
                            !        FNRDDD = FNRDDD + gauss_term(k)
                            !end if
                        end if
                        !Attention: Fluids with exponents in the exponential part, which are not 2, are not considered here
                        !e.g., for R-125 (Lemmon and Jacobsen, 2005)
                    end do

                    FNRDER(8) = FNRDDD + HSDER(8)
                end if

                if (GETDERS(9) == 1) then
                    !*****************************************************
                    !    FNRTTT     3RD DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau
                    !*****************************************************
                    !All derivatives of alpha_r with respect to ONLY tau
                    !are zero at zero density!!
                    FNRTTT = 0.D0
                    FNRDER(9) = FNRTTT + HSDER(9)
                end if

                if (GETDERS(10) == 1) then
                    !*****************************************************
                    !    FNRDDT     3RD DERIVATIVE OF alpha_r WITH RESPECT TO
                    !               tau, del, and del
                    !*****************************************************
                    FNRDDT = 0.D0

                    !!Summation over the regular terms
                    !Andreas, September 2015
                    !------------------------------------------------------------------------------------------------------------------------
                    do i = 1, gl%eos_coeff%nreg(nrsubst)

                        if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                            reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * (gl%eos_coeff%di(i,nrsubst) - 1.D0) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                            !reg_term(i) = 2d0*gl%eos_coeff%ni(i,nrsubst)*gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                        end if
                        !
                        if (dabs(gl%eos_coeff%gama(i,nrsubst)) > 1.D-12) then !Exponential terms

                            !Different cases:   1)  di - 2 < 0 --> the polynomial part of the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2)  di - 2 = 0 --> delta vanishes, term contributes to B
                            !                   3)  di - 2 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            !                   More different cases for the exponential terms required. See documentation

                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                !reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * (gl%eos_coeff%di(i,nrsubst) - 1.D0) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                                reg_term(i) = 2d0*gl%eos_coeff%ni(i,nrsubst)*gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%p_i(i,nrsubst) - 2.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * (-(2.D0 * gl%eos_coeff%di(i,nrsubst) - 1.D0 + gl%eos_coeff%gama(i,nrsubst)* gl%eos_coeff%p_i(i,nrsubst))) * gl%eos_coeff%gama(i,nrsubst)* gl%eos_coeff%p_i(i,nrsubst) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                            end if
                            if(NINT(gl%eos_coeff%di(i,nrsubst) + 2.D0 * gl%eos_coeff%p_i(i,nrsubst) - 2.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%gama(i,nrsubst)**2 * gl%eos_coeff%p_i(i,nrsubst)**2 * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                            end if

                        else   !Polynomial terms:

                            !Different cases:   1) di - 2 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                            !                   2) di - 2 = 0 --> delta vanishes, term contributes to B
                            !                   3) di - 2 > 0 --> delta to the power di-1 left --> does not contribute to B, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                reg_term(i) = gl%eos_coeff%ni(i,nrsubst) * gl%eos_coeff%di(i,nrsubst) * (gl%eos_coeff%di(i,nrsubst) - 1.D0) * gl%eos_coeff%ti(i,nrsubst) * tau**(gl%eos_coeff%ti(i,nrsubst) - 1.D0)
                            end if

                        end if
                        FNRDDT = FNRDDT + reg_term(i)
                    end do
                    !------------------------------------------------------------------------------------------------------------------------
                    !Gaussian bell-shaped terms
                    !added by Sven, Nov 2019

                    do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                        k = i - gl%eos_coeff%nreg(nrsubst)
                        DELEPS(k) = DEL - gl%eos_coeff%eps(k,nrsubst)
                        TAUGAM(k) = TAU - gl%eos_coeff%gam(k,nrsubst)
                        gauss_term(k) = 0.D0
                        !Different cases:   1) di - 2 < 0 --> the derivative vanishes --> DOES NOT NEED TO BE EVALUATED
                        !                   2) di - 2 = 0 --> delta vanishes, term contributes to D
                        !                   3) di - 2 > 0 --> delta to the power di-1 left --> does not contribute to D, since delta = 0!!  --> DOES NOT NEED TO BE EVALUATED

                        if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                            if(NINT(gl%eos_coeff%di(i,nrsubst) - 2.D0) == 0) then
                                gauss_term(k) =  -2d0*gl%eos_coeff%ni(i,nrsubst) *  tau**(gl%eos_coeff%ti(i,nrsubst)-1d0) * dexp(gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)*DELEPS(k) + gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)*TAUGAM(k)) &
                                    & *(gl%eos_coeff%ti(i,nrsubst) + 2d0*gl%eos_coeff%beta(k,nrsubst)*TAU*(TAUGAM(k)))
                                FNRDDT = FNRDDT +gauss_term(k)
                            endif
                        end if
                        !Attention: Fluids with exponents in the exponential part, which are not 2, are not considered here
                        !e.g., for R-125 (Lemmon and Jacobsen, 2005)
                    end do
                    !Neglect Gaussian bell shaped terms! --> Check and discuss if valid. !NO

                    FNRDER(10) = FNRDDT + HSDER(10)
                end if

                return

            end if

            !!Summation over the special GBS terms
            if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                DelDeriv = 0.d0
                TauDeriv = 0.d0

                do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                    call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                    FNR = FNR + SGBSCALC

                end do

            end if

            !!Summation over the nonanalytic terms
            !do i = (gl%eos_coeff%nreg(nrsubst) + I_GBS(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + ncrt(nrsubst))
            do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)

                DelDeriv = 0.d0
                TauDeriv = 0.d0

                call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)

                FNR = FNR + NACALC
            end do



            !*****************************************************
            if (GETDERS(1) == 1) then
                FNRDER(1) = FNR + FNRASSODER(1) + HSDER(1)  !write the result into the output vector
            end if
            !*****************************************************
            if (GETDERS(2) == 1) then
                !*****************************************************
                !    FNRD       1ST DERIVATIVE OF H WITH RESPECT TO
                !               del MULTIPLIED WITH del [-]
                !*****************************************************
                FNRD = 0.d0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRD = FNRD + reg_term(i) * (gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if (dabs(DELEPS(k)) < 1.d-10) then   ! for del = 1 and epsylon = 1 the exponent becomes unity
                        FNRD = FNRD + gauss_term(k) * gl%eos_coeff%di(i,nrsubst)
                    else if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                        !FNRD = FNRD + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + NINT(gl%eos_coeff%pli(k,nrsubst))*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**NINT((gl%eos_coeff%pli(k,nrsubst)-1.D0)))
                        FNRD = FNRD + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + 2.d0*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k))
                    else
                        FNRD = FNRD + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))
                    end if
                end do


                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 1.d0
                    TauDeriv = 0.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRD = FNRD + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1.d0),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 1
                    TauDeriv = 0
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRD = FNRD + NACALC*del
                    end if
                end do

                FNRDER(2) = FNRD + FNRASSODER(2)+HSDER(2)     !write the result into the output vector
            end if
            !*****************************************************
            if (GETDERS(3) == 1) then
                !*****************************************************
                !   FNRDD       2ND DERIVATIVE OF H WITH RESPECT TO
                !               del MULTIPLIED WITH del^2 [-]
                !*****************************************************
                FNRDD = 0.D0

                !!Summation over the regualar terms
                do i = gl%eos_coeff%nreg(nrsubst),1, -1
                    FNRDD = FNRDD + reg_term(i) * ((gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))* &
                        (gl%eos_coeff%di(i,nrsubst) - 1.d0 - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst)) - &
                        gl%eos_coeff%gama(i,nrsubst)*(gl%eos_coeff%p_i(i,nrsubst)**2)*(delpi(i,nrsubst)))
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = gl%eos_coeff%nreg(nrsubst) + 1, (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then       !if the exponents in the Gauss exponents are equal to 2
                        FNRDD = FNRDD + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + 2.d0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k))**2 &
                            - gl%eos_coeff%di(i,nrsubst) + 2.d0*gl%eos_coeff%eta(k,nrsubst)*del2)
                    else
                        FNRDD = FNRDD + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))**2 &
                            - gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0)*(gl%eos_coeff%pli(k,nrsubst)-1.D0))
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 2.d0
                    TauDeriv = 0.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRDD = FNRDD + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 2
                    TauDeriv = 0
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRDD = FNRDD + NACALC*del*del
                    end if
                end do

                FNRDER(3) = FNRDD + FNRASSODER(3)+HSDER(3)     !write the result into the output vector
            end if


            !*****************************************************
            if (GETDERS(4) == 1) then
                !*****************************************************
                !  FNRT         1ST DERIVATIVE OF H WITH RESPECT TO
                !               TAU MULTIPLIED WITH TAU [-]
                !*****************************************************
                FNRT = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRT = FNRT + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if (gl%eos_coeff%pli(k,nrsubst) == 2.d0.and.gl%eos_coeff%tli(k,nrsubst) == 2.D0) then       !if the exponents in the Gauss exponents are equal to 2
                        FNRT = FNRT + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + 2.d0*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))
                    else
                        if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                            FNRT = FNRT + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + (gl%eos_coeff%tli(k,nrsubst))*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT((gl%eos_coeff%tli(k,nrsubst)-1.D0)))
                        else
                            FNRT = FNRT + gauss_term(k) * (gl%eos_coeff%ti(i,nrsubst) + (gl%eos_coeff%tli(k,nrsubst))*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**((gl%eos_coeff%tli(k,nrsubst)-1.D0)))
                        end if
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 0.d0
                    TauDeriv = 1.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRT = FNRT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 0
                    TauDeriv = 1
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRT = FNRT + NACALC*tau
                    end if
                end do

                FNRDER(4) = FNRT + FNRASSODER(4)+HSDER(4)     !write the result into the output vector
            end if
            !*****************************************************
            if (GETDERS(5) == 1) then
                !*****************************************************
                !  FNRTT        2ND DERIVATIVE OF H WITH RESPECT TO
                !               TAU MULTIPLIED WITH TAU^2 [-]
                !*****************************************************
                FNRTT = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRTT = FNRTT + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%ti(i,nrsubst) - 1.d0)
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if (gl%eos_coeff%pli(k,nrsubst) == 2.d0.and.gl%eos_coeff%tli(k,nrsubst) == 2.D0) then       !if the exponents in the Gauss exponents are equal to 2
                        FNRTT = FNRTT + gauss_term(k) * ((gl%eos_coeff%ti(i,nrsubst) + 2.d0*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k))**2 &
                            - gl%eos_coeff%ti(i,nrsubst) + 2.d0*gl%eos_coeff%beta(k,nrsubst)*tau2)
                    else
                        if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                            FNRTT = FNRTT + gauss_term(k) * ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))**2 &
                                - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-2.D0)*(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                        else
                            FNRTT = FNRTT + gauss_term(k) * ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))**2 &
                                - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-2.D0)*(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                        end if
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 0.d0
                    TauDeriv = 2.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRTT = FNRTT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 0
                    TauDeriv = 2
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRTT = FNRTT + NACALC*tau*tau
                    end if
                end do

                FNRDER(5) = FNRTT + FNRASSODER(5)+HSDER(5)     !write the result into the output vector
            end if
            !*****************************************************
            if (GETDERS(6) == 1) then
                !*****************************************************
                !   FNRDT       1ST MIXED DERIVATIVE OF F WITH RESPECT
                !               TO del AND tau MULTIPLIED WITH
                !               del*tau [-]
                !*****************************************************
                FNRDT = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRDT = FNRDT + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))
                end do


                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if ((gl%eos_coeff%pli(k,nrsubst) == 2.d0) .and. (gl%eos_coeff%tli(k,nrsubst) == 2.d0)) then       !if the exponents in the Gauss exponents are equal to 2
                        !FNRDT = FNRDT + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + 2.D0*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)) &
                        !    *(gl%eos_coeff%ti(i,nrsubst) + 2.D0*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))

                        !für spätere Struktur der gesamten Routine
                        taupart = (gl%eos_coeff%ti(i,nrsubst) + 2.d0*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))
                        delpart = (gl%eos_coeff%di(i,nrsubst) + 2.d0*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k))

                        FNRDT = FNRDT +gauss_term(k)*taupart*delpart
                    else
                        if ((gl%eos_coeff%tli(k,nrsubst)-2.d0) == 0.d0) then
                            FNRDT = FNRDT + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + del*gl%eos_coeff%eta(k,nrsubst)*gl%eos_coeff%pli(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0)) &
                                *(gl%eos_coeff%ti(i,nrsubst) + tau*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                        else
                            FNRDT = FNRDT + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + del*gl%eos_coeff%eta(k,nrsubst)*gl%eos_coeff%pli(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0)) &
                                *(gl%eos_coeff%ti(i,nrsubst) + tau*gl%eos_coeff%beta(k,nrsubst)*gl%eos_coeff%tli(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                        end if
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 1.d0
                    TauDeriv = 1.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRDT = FNRDT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 1
                    TauDeriv = 1
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRDT = FNRDT + NACALC*del*tau
                    end if
                end do

                FNRDER(6) = FNRDT + FNRASSODER(6) + HSDER(6)     !write the result into the output vector
            end if


            !*****************************************************
            if (GETDERS(7) == 1) then
                !*****************************************************
                !     FNRDTT    2ND MIXED DERIVATIVE OF F WITH RESPECT
                !               TO TAU AND del MULTIPLIED WITH
                !               del*TAU^2 [-]
                !*****************************************************
                FNRDTT = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRDTT = FNRDTT + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%ti(i,nrsubst) - 1.d0)*(gl%eos_coeff%di(i,nrsubst) - gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if ((gl%eos_coeff%pli(k,nrsubst) == 2.d0) .and. (gl%eos_coeff%tli(k,nrsubst) == 2.d0)) then       !if the exponents in the Gauss exponents are equal to 2
                        FNRDTT = FNRDTT + gauss_term(k) * (gl%eos_coeff%di(i,nrsubst) + 2.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)) &
                            *((gl%eos_coeff%ti(i,nrsubst) + 2.D0*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k))**2 - gl%eos_coeff%ti(i,nrsubst) + 2.D0*gl%eos_coeff%beta(k,nrsubst)*tau2)
                    else
                        if ((gl%eos_coeff%tli(k,nrsubst) - 2.d0) == 0.d0) then
                            FNRDTT = FNRDTT + gauss_term(k) * &
                                (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.d0))**2*&
                                ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.d0))**2&
                                - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*(gl%eos_coeff%tli(k,nrsubst)-1.d0)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-2.d0))
                        else

                            !für spätere Struktur der gesamten Routine
                            tau2part = ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))**2 &
                                - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-2.D0)*(gl%eos_coeff%tli(k,nrsubst)-1.D0))
                            delpart = (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*del*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))

                            FNRDTT = FNRDTT + gauss_term(k)*tau2part*delpart

                            !FNRDTT = FNRDTT + gauss_term(k) * &
                            !    (gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.d0))**2*&
                            !    ((gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.d0))**2&
                            !    - gl%eos_coeff%ti(i,nrsubst) + gl%eos_coeff%tli(k,nrsubst)*gl%eos_coeff%beta(k,nrsubst)*tau2*(gl%eos_coeff%tli(k,nrsubst)-1.d0)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-2.d0))
                        end if
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 1.d0
                    TauDeriv = 2.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRDTT = FNRDTT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 1
                    TauDeriv = 2
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRDTT = FNRDTT + NACALC*del*tau*tau
                    end if
                end do

                FNRDER(7) = FNRDTT + FNRASSODER(7) + HSDER(7)    !write the result into the output vector
            end if


            !*****************************************************
            if (GETDERS(8) == 1) then
                !*****************************************************
                !     FNRDDD    3RD DERIVATIVE OF F WITH RESPECT
                !               TO del MULTIPLIED WITH del^3 [-]
                !*****************************************************
                FNRDDD = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRDDD = FNRDDD + reg_term(i) * (gl%eos_coeff%di(i,nrsubst)*(gl%eos_coeff%di(i,nrsubst) - 1.d0)*(gl%eos_coeff%di(i,nrsubst) - 2.d0) + &    !Theresa added .d0
                        gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst)*( - 2.d0 + 6.d0*gl%eos_coeff%di(i,nrsubst) &
                        - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 - 3.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst) + 3.d0*gl%eos_coeff%p_i(i,nrsubst) - gl%eos_coeff%p_i(i,nrsubst)**2) &
                        + 3.d0*gl%eos_coeff%gama(i,nrsubst)**2*gl%eos_coeff%p_i(i,nrsubst)**2*del**(2.d0*gl%eos_coeff%p_i(i,nrsubst))*(gl%eos_coeff%di(i,nrsubst) - 1.d0 + gl%eos_coeff%p_i(i,nrsubst))&
                        - gl%eos_coeff%gama(i,nrsubst)**3*gl%eos_coeff%p_i(i,nrsubst)**3*del**(3*gl%eos_coeff%p_i(i,nrsubst)))
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)
                    if (abs(DELEPS(k)) < 1.d-10) then   ! for del = 1 and epsylon = 1 the exponent becomes unity
                        FNRDDD = FNRDDD + gauss_term(k)*gl%eos_coeff%di(i,nrsubst)*(gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%di(i,nrsubst) - 3.d0*gl%eos_coeff%di(i,nrsubst) + 2.d0) !Theresa added .d0
                    else if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then
                        FNRDDD = FNRDDD + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + 2.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k))**3 - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 &
                            & + 2.d0*gl%eos_coeff%di(i,nrsubst) + 6.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2 - 6.D0*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)*(gl%eos_coeff%di(i,nrsubst) &
                            & - 2.D0*gl%eos_coeff%eta(k,nrsubst)*del2))
                    else
                        FNRDDD = FNRDDD + gauss_term(k) * ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 1.d0))**3 &
                            & - 3.d0*gl%eos_coeff%di(i,nrsubst)**2 + 2.d0*gl%eos_coeff%di(i,nrsubst) + 3.d0*gl%eos_coeff%di(i,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2*gl%eos_coeff%pli(k,nrsubst)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)* &
                            & DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0) - 3.d0*gl%eos_coeff%di(i,nrsubst)*del*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0) &
                            & + 3.D0*del3*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)**2*gl%eos_coeff%eta(k,nrsubst)**2*DELEPS(k)**(2.d0*gl%eos_coeff%pli(k,nrsubst) - 3.d0) &
                            & + (gl%eos_coeff%pli(k,nrsubst) - 2.d0)*(gl%eos_coeff%pli(k,nrsubst) - 1.d0)*gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst) - 3.d0))
                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 3.d0
                    TauDeriv = 0.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRDDD = FNRDDD + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 3
                    TauDeriv = 0
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRDDD = FNRDDD + NACALC*del*del*del
                    end if
                end do

                FNRDER(8) = FNRDDD + FNRASSODER(8) + HSDER(8)    !write the result into the output vector
            end if

            !Thu
            !****************************************************************************************
            if (GETDERS(9) == 1) then
                !****************************************************************************************
                !  FNRTTT       3RD DERIVATIVE OF H WITH RESPECT TO
                !               TAU MULTIPLIED WITH TAU^3 [-]
                !****************************************************************************************
                FNRTTT = 0.D0

                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRTTT = FNRTTT + reg_term(i) * gl%eos_coeff%ti(i,nrsubst)*(gl%eos_coeff%ti(i,nrsubst) - 1.d0)*(gl%eos_coeff%ti(i,nrsubst) - 2.d0)
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)

                    !if the exponents in the Gauss exponents are equal to 2
                    if (gl%eos_coeff%pli(k,nrsubst) == 2.d0) then

                        tau_inv = 1.d0/tau
                        tau_inv2 = tau_inv*tau_inv
                        tau_inv3 = tau_inv2*tau_inv

                        FNRTTT = FNRTTT + gauss_term(k)*tau3*&
                            (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))**2&
                            -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                            (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))+2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

                    else if ((gl%eos_coeff%tli(k,nrsubst) - 2.d0) == 0.d0) then

                        FNRTTT = FNRTTT + gauss_term(k)*tau3*&
                            (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))**NINT(gl%eos_coeff%tli(k,nrsubst))&
                            -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                            (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**NINT(gl%eos_coeff%tli(k,nrsubst)-1.D0))&
                            +2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

                    else

                        FNRTTT = FNRTTT + gauss_term(k)*tau3*&
                            (((gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))**gl%eos_coeff%tli(k,nrsubst)&
                            -3.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv2+6.D0*gl%eos_coeff%beta(k,nrsubst))*&
                            (gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%tli(k,nrsubst)-1.D0))+2.D0*gl%eos_coeff%ti(i,nrsubst)*tau_inv3)

                    end if
                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 0.d0
                    TauDeriv = 3.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRTTT = FNRTTT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 0
                    TauDeriv = 3
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRTTT = FNRTTT + NACALC*tau*tau*tau
                    end if
                end do

                FNRDER(9) = FNRTTT + FNRASSODER(9) + HSDER(9)    !write the result into the output vector
            end if



            !*************************************************************************************
            if (GETDERS(10) == 1) then
                !*************************************************************************************
                !   FNRDDT      2ND  MIXED DERIVATIVE OF F WITH RESPECT
                !               TO del AND tau MULTIPLIED WITH
                !               del*del*tau [-]
                !*************************************************************************************

                FNRDDT = 0.D0
                !!Summation over the regular terms
                do i = 1, gl%eos_coeff%nreg(nrsubst)
                    FNRDDT = FNRDDT + reg_term(i)*((gl%eos_coeff%di(i,nrsubst)-gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))*&
                        (gl%eos_coeff%di(i,nrsubst)-1.d0-gl%eos_coeff%gama(i,nrsubst)*gl%eos_coeff%p_i(i,nrsubst)*delpi(i,nrsubst))-gl%eos_coeff%gama(i,nrsubst)**2*gl%eos_coeff%p_i(i,nrsubst)**2*delpi(i,nrsubst))*gl%eos_coeff%ti(i,nrsubst)
                end do

                !!Summation over the Gaussian bell-shaped terms
                do i = (gl%eos_coeff%nreg(nrsubst) + 1), (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst))
                    k = i - gl%eos_coeff%nreg(nrsubst)

                    del_inv = 1.d0/del
                    del_inv2 = del_inv*del_inv
                    tau_inv = 1.d0/tau

                    if (abs(DELEPS(k)) < 1.d-10) then       ! for del = 1 and epsylon = 1 the exponent becomes unity

                        FNRDDT = FNRDDT + gauss_term(k) * (del2*tau*(((gl%eos_coeff%di(i,nrsubst)*del_inv)**2 &
                            - gl%eos_coeff%di(i,nrsubst)*del_inv2 + 2.D0*gl%eos_coeff%eta(k,nrsubst))*(gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst))))


                    else if ((gl%eos_coeff%pli(k,nrsubst) == 2.d0) .and. (gl%eos_coeff%tli(k,nrsubst) == 2.d0)) then          !if the exponents in the Gauss exponents are equal to 2

                        FNRDDT = FNRDDT + gauss_term(k) * &
                            (del2*tau*(((gl%eos_coeff%di(i,nrsubst)*del_inv + 2.D0*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k))**2 &
                            - gl%eos_coeff%di(i,nrsubst)*del_inv2 + 2.D0*gl%eos_coeff%eta(k,nrsubst))*(gl%eos_coeff%ti(i,nrsubst)*tau_inv&
                            +2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k))))

                    else

                        !für spätere Struktur der gesamten Routine
                        taupart = (gl%eos_coeff%ti(i,nrsubst) + (gl%eos_coeff%tli(k,nrsubst))*tau*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**((gl%eos_coeff%tli(k,nrsubst)-1.D0)))
                        del2part = ((gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))**2 &
                            - gl%eos_coeff%di(i,nrsubst) + gl%eos_coeff%pli(k,nrsubst)*gl%eos_coeff%eta(k,nrsubst)*del2*DELEPS(k)**(gl%eos_coeff%pli(k,nrsubst)-2.D0)*(gl%eos_coeff%pli(k,nrsubst)-1.D0))

                        FNRDDT = FNRDDT + gauss_term(k)*taupart*del2part

                        !FNRDDT = FNRDDT + gauss_term(k) * &
                        !    (del2*tau*((gl%eos_coeff%di(i,nrsubst)*del_inv + 2.D0*gl%eos_coeff%eta(k,nrsubst)*DELEPS(k))**gl%eos_coeff%pli(k,nrsubst) &
                        !    - gl%eos_coeff%di(i,nrsubst)*del_inv2 + 2.D0*gl%eos_coeff%eta(k,nrsubst)*(gl%eos_coeff%ti(i,nrsubst)*tau_inv+2.D0*gl%eos_coeff%beta(k,nrsubst)*TAUGAM(k)**(gl%eos_coeff%pli(k,nrsubst)-1.D0))))

                    end if

                end do

                !!Summation over the special GBS terms
                if (gl%eos_coeff%nna(nrsubst) > 0.d0) then

                    DelDeriv = 2.d0
                    TauDeriv = 1.d0

                    do i = (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + 1.d0) , (gl%eos_coeff%nreg(nrsubst) + gl%eos_coeff%I_GBS(nrsubst) + gl%eos_coeff%nna(nrsubst))

                        call SGBSTERMDERIVS(gl,SGBSCALC, i, del, tau, DelDeriv, TauDeriv, nrsubst)

                        FNRDDT = FNRDDT + SGBSCALC

                    end do

                end if

                !!Summation over the nonanalytic terms
                !            do i = (nreg(nrsubst) + I_GBS(nrsubst) + 1), (nreg(nrsubst) + ncrt(nrsubst))
                do nn = (gl%eos_coeff%I_GBS(nrsubst) + 1),gl%eos_coeff%ncrt(nrsubst)
                    DelDeriv = 2
                    TauDeriv = 1
                    if (region(nn-gl%eos_coeff%I_GBS(nrsubst)) == 1) then
                        call NATERMDERIVS(gl,NACALC, nn, del, tau, DelDeriv, TauDeriv, region(nn-gl%eos_coeff%I_GBS(nrsubst)), nrsubst)
                        FNRDDT = FNRDDT + NACALC*del*del*tau
                    end if
                end do

                FNRDER(10) = FNRDDT + FNRASSODER(10) + HSDER(10)                           !write the result into the output vector
            end if


            ! fourth mixed derivatives of helmholtz energy
            !Sven Pohl, 14.09.2020

            !*****************************************************************************************
            if (GETDERS(12) == 1) then
                !*************************************************************************************
                !   FNRDDDT     4th  MIXED DERIVATIVE OF F WITH RESPECT
                !               TO del AND tau MULTIPLIED WITH
                !               tau*del*del*del [-]
                !*************************************************************************************
                FNRDDDT = 0.D0
                FNRDDD_F = FNRDDD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)
                FNRT_F   = FNRT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,region)
                FNRDDDT = sum(FNRDDD_F*FNRT_F)
                FNRDER(12) = FNRDDDT
            end if

            !*****************************************************************************************
            if (GETDERS(13) == 1) then !
                !*************************************************************************************
                !   FNRDDTT     4th  MIXED DERIVATIVE OF F WITH RESPECT
                !               TO del AND tau MULTIPLIED WITH
                !               tau*tau*del*del [-]
                !*************************************************************************************
                FNRDDTT = 0.D0
                FNRDD_F = FNRDD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)
                FNRTT_F = FNRTT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,DELEPS,region)
                FNRDDTT = sum(FNRTT_F*FNRDD_F)
                FNRDER(13) = FNRDDDT
            end if

            !*****************************************************************************************
            if (GETDERS(14) == 1) then !on work
                !*************************************************************************************
                !   FNRDDTT     4th  MIXED DERIVATIVE OF F WITH RESPECT
                !               TO del AND tau MULTIPLIED WITH
                !               tau*tau*tau*del [-]
                !*************************************************************************************
                FNRDTTT = 0.D0
                FNRD_F = FNRDD_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,DELEPS,region)
                FNRTTT_F = FNRTTT_FUNC(gl,tau,del,nrsubst,reg_term,delpi,gauss_term,TAUGAM,region)
                FNRDTTT = sum(FNRD_F*FNRTTt_F)
                FNRDER(14) = FNRDDDT
            end if

            !This if-case is important. The correction of the ratio of the gas constants
            ! must be considered for mixtures, but
            ! NOT for calculating the integration constants of the corresponding pure fluids
            if ((.not.gl%ref) .and. (gl%ncomp .gt. 1)) then
                call R_mix_calc(gl,Rmix)
                FNRDER = FNRDER * gl%REQ(nrsubst)/Rmix
            end if


            !deallocate(delpi)
            !deallocate(reg_term)


            !Crossover model: implemented by Sven Pohl
            !if(gl%is_crs) then
            !    call get_crs_derivs(gl,FNRDER,GETDER,tau,del)
            !end if
            !if(any(gl%eos_coeff%nr_nf .ne. 0)) then
            !    if(GETDERS(1) == 1)
            !    if(GETDERS(1) == 2)
            !    if(GETDERS(1) == 3)
            !    if(GETDERS(1) == 4)
            !    if(GETDERS(1) == 5)
            !    if(GETDERS(1) == 6)
            !    if(GETDERS(1) == 7)
            !    if(GETDERS(1) == 8)
            !    if(GETDERS(1) == 9)
            !    if(GETDERS(1) == 10)
            !    if(GETDERS(1) == 11)
            !end if




        elseif (gl%Eq_type(nrsubst) == 2) then !SRK pure fluid equation is used

            !Check whether the SRK was already initialized with the actual temperature. If it was already initialized, it does not have to be initialized again
            !This applies only for pure substances, in mixtures the SRK always has to be updated for every component
            !Andreas March 2013

            !OLD CODE REPLACED WITH NEW CODE FOR CUBIC EQUATIONS OF STATE (SRK AND PENG-ROBINSON)
            !Andreas Jäger, March 2017
            srk_t = Tau/gl%Tc(nrsubst) !1.d0 / (Tc(nrsubst) / Tau)
            srk_t2 = srk_t*srk_t
            srk_t3 = srk_t2*srk_t
            rhoc_del = gl%rhoc(nrsubst) * del
            bi_rhoc_del = 1.d0 +  gl%bi_SRK(nrsubst) * rhoc_del

            if ((dabs(Temperature - gl%Temp_init) > 1.D-14) .or. (gl%ncomp > 1)) then
                call update_parameters(gl,Temperature, nrsubst)
            end if

            !The MIXDERIVSFNR_HIGHER_CUBIC routine also works for pure substances, because
            !in the routine which is called before (update_parameters), the parameters a and b for the pure fluid
            !are calculated and saved in global variables which are used in this routine
            !Andreas Jäger, March 2017
            !call MIXDERIVSFNR_HIGHER_CUBIC (Temperature, Density, GETDERS, FNRDER)

            !OLD CODE REPLACED WITH NEW CODE FOR CUBIC EQUATIONS OF STATE (SRK AND PENG-ROBINSON)
            !Can be deleted if no problems occur
            !Andreas Jäger, March 2017
            !Old routines for now kept! There was indeed a problem when using cubics in the multi-fluid model.
            !This could maybe relatively easily be fixed, but for now, the old routines are used here instead of the new ones
            !Andreas Jäger, September 2017
            if (GETDERS(1) == 1) then
                FNRDER(1) = -dlog(1.D0-gl%bi_SRK(nrsubst)*rhoc_del)- &
                    & gl%ai_SRK(nrsubst)*srk_t/gl%R_SRK/gl%bi_SRK(nrsubst)*dlog(bi_rhoc_del)
            end if
            if (GETDERS(2) == 1) then
                FNRDER(2) = (gl%bi_SRK(nrsubst) * gl%rhoc(nrsubst) / (1.D0 - gl%bi_SRK(nrsubst) * rhoc_del) - &
                    & gl%ai_SRK(nrsubst) * gl%rhoc(nrsubst) * tau / (gl%R_SRK * gl%tc(nrsubst)*bi_rhoc_del))*del
            end if
            if (GETDERS(3) == 1) then
                FNRDER(3) = del*del*(((gl%bi_SRK(nrsubst)*gl%rhoc(nrsubst))/(1-gl%bi_SRK(nrsubst)*rhoc_del))**2 + &
                    & gl%ai_SRK(nrsubst) * srk_t / gl%R_SRK * gl%bi_SRK(nrsubst) * gl%rhoc(nrsubst)* gl%rhoc(nrsubst) / (bi_rhoc_del*bi_rhoc_del))
            end if
            if (GETDERS(4) == 1) then    !Lars H., 14.02.2013
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                dT_dTau = -gl%Tc(nrsubst) / (tau * tau)
                FNRDER(4) = -1.d0 / (gl%R_SRK * gl%bi_SRK(nrsubst)) * dlog(bi_rhoc_del) * &
                    & (da_dTau * (gl%Tc(nrsubst) / Tau) - gl%ai_SRK(nrsubst) * dT_dTau) * srk_t2 * Tau
            end if
            if (GETDERS(5) == 1) then    !Lars H., 14.02.2013
                dT_dTau = -gl%Tc(nrsubst) / (tau * tau)
                d2T_dTau2 = 2.d0 * gl%Tc(nrsubst) /  (tau * tau * tau)
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, nrsubst)
                FNRDER(5) = -1.d0 / (gl%R_SRK * gl%bi_SRK(nrsubst)) * dlog(bi_rhoc_del) * &
                    & (d2a_dTau2 * srk_t - gl%ai_SRK(nrsubst) * d2T_dTau2 * srk_t2 - &
                    & 2.d0 * dT_dTau * da_dTau * srk_t2 + &
                    & 2.d0 * gl%ai_SRK(nrsubst) * dT_dTau * dT_dTau * srk_t3) * Tau * Tau
            end if
            if (GETDERS(6) == 1) then
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                dT_dTau = -gl%Tc(nrsubst) / (tau * tau)
                !Andreas Nov. 2013
                FNRDER(6) = - (gl%rhoc(nrsubst) / gl%R_SRK / (1.D0 + gl%bi_SRK(nrsubst)*rhoc_del) * &
                    & (da_dTau * (gl%Tc(nrsubst) / Tau) - gl%ai_SRK(nrsubst) * dT_dTau) * srk_t2) * del * tau
            end if
            if (GETDERS(7) == 1) then    !Lars H., 14.02.2013
                tau2 = tau*tau
                tau3 = tau2*tau
                dT_dTau = -gl%Tc(nrsubst) / tau2
                d2T_dTau2 = 2.d0 * gl%Tc(nrsubst) / tau3
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, nrsubst)
                FNRDER(7) = -gl%rhoc(nrsubst) / (gl%R_SRK * (bi_rhoc_del)) * &
                    & (d2a_dTau2 * srk_t - gl%ai_SRK(nrsubst) * d2T_dTau2 * srk_t2 - &
                    & 2.d0 * dT_dTau * da_dTau * srk_t2 + &
                    & 2.d0 * gl%ai_SRK(nrsubst) * dT_dTau * dT_dTau * srk_t3) * tau2 * del
            end if
            if (GETDERS(8) == 1) then
                FNRDER(8) = del*del*del*(2.D0*(((gl%bi_SRK(nrsubst)*gl%rhoc(nrsubst))/(1-gl%bi_SRK(nrsubst)*rhoc_del))**3 - &
                    & gl%ai_SRK(nrsubst) *srk_t / gl%R_SRK * gl%bi_SRK(nrsubst)*gl%bi_SRK(nrsubst) * &
                    & gl%rhoc(nrsubst)*gl%rhoc(nrsubst)*gl%rhoc(nrsubst) / (bi_rhoc_del*bi_rhoc_del*bi_rhoc_del)))

            end if
            if (GETDERS(9) == 1) then    !Lars H., 14.02.2013
                srk_t4 = srk_t3 * srk_t
                tau2 = tau*tau
                tau3 = tau2*tau
                tau4 = tau3*tau
                dT_dTau = -gl%Tc(nrsubst) / tau2
                d2T_dTau2 = 2.d0 * gl%Tc(nrsubst) / tau3
                d3T_dTau3 = - 6.d0 * gl%Tc(nrsubst) / tau4
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, nrsubst)
                d3a_dTau3 = d3a_SRK_dtau3(gl,Temperature, nrsubst)
                FNRDER(9) = -1.d0 / (gl%R_SRK * gl%bi_SRK(nrsubst)) * dlog(bi_rhoc_del) * &
                    & (d3a_dTau3 * srk_t - gl%ai_SRK(nrsubst) * d3T_dTau3 * srk_t2 - &
                    & 3.d0 * d2T_dTau2 * da_dTau * srk_t2 - &
                    & 3.d0 * dT_dTau * d2a_dTau2 * srk_t2 + &
                    & 6.d0 * da_dTau * dT_dTau * dT_dTau * srk_t3 + &
                    & 6.d0 * gl%ai_SRK(nrsubst) * d2T_dTau2 * dT_dTau * srk_t4 - &
                    & 6.d0 * gl%ai_SRK(nrsubst) * dT_dTau* dT_dTau* dT_dTau * srk_t4) * tau3

            end if
            if (GETDERS(10) == 1) then    !Lars H., 14.02.2013
                dT_dTau = -gl%Tc(nrsubst) / (tau * tau)
                da_dTau = da_SRK_dtau(gl,Temperature, nrsubst)
                FNRDER(10) = gl%rhoc(nrsubst)*gl%rhoc(nrsubst) * gl%bi_SRK(nrsubst) / (gl%R_SRK * bi_rhoc_del*bi_rhoc_del) * &
                    & (da_dTau * srk_t - gl%ai_SRK(nrsubst) * dT_dTau * srk_t2) * tau * del * del
            end if

            !Stefan Feb 2014 (see Pit Podleschny)
        else If (gl%Eq_type(nrsubst) == 3) then  !PR pure fluid equation is used

            !Andreas, Nov 2015: Update PR parameters if temperature was changed (in mixtures always update)
            if ((dabs(Temperature - gl%Temp_init) > 1.D-14) .or. (gl%ncomp > 1)) then
                call update_parameters(gl,Temperature, nrsubst)
            end if

            !The MIXDERIVSFNR_HIGHER_CUBIC routine also works for pure substances, because
            !in the routine which is called before (update_parameters), the parameters a and b for the pure fluid
            !are calculated and saved in global variables which are used in this routine
            !Andreas Jäger, March 2017
            !call MIXDERIVSFNR_HIGHER_CUBIC (Temperature, Density, GETDERS, FNRDER)

            !OLD CODE REPLACED WITH NEW CODE FOR CUBIC EQUATIONS OF STATE (SRK AND PENG-ROBINSON)
            !Can be deleted if no problems occur
            !Andreas Jäger, March 2017
            !Old routines for now kept! There was indeed a problem when using cubics in the multi-fluid model.
            !This could maybe relatively easily be fixed, but for now, the old routines are used here instead of the new ones
            !Andreas Jäger, September 2017
            if (GETDER(1) == 1) then
                FNRDER(1) = -dlog(1.D0 - gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) * del) - gl%ai_PR(nrsubst) * tau * (1.d0 / &
                    & (2.d0 * 2.d0 ** 0.5 * gl%bi_PR(nrsubst) * gl%R_PR * gl%tc(nrsubst))) * dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * &
                    & del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst)) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst)))
            end if
            if (GETDER(2) == 1) then
                FNRDER(2) = (gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)) - gl%ai_PR(nrsubst) * tau / (2.d0 * &
                    & 2.d0 ** 0.5 * gl%R_PR * gl%tc(nrsubst) * gl%bi_PR(nrsubst)) * ((2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / &
                    & (1.d0 + (2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)) + (2.d0 ** 0.5 - 1.d0) * &
                    & gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))) * del
            end if
            if (GETDER(3) == 1) then
                FNRDER(3) = ((gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**2 - gl%ai_PR(nrsubst) * tau / (2.d0 * &
                    & 2.d0 ** 0.5 * gl%R_PR * gl%tc(nrsubst) * gl%bi_PR(nrsubst)) * (-((2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / &
                    & (1.d0 + (2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**2 + ((2.d0 ** 0.5 - 1.d0) * &
                    & gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**2)) * &
                    & del * del
            end if
            if (GETDER(4) == 1) then
                da_dtau = da_PR_dtau(gl,Temperature, nrsubst)
                FNRDER(4) = (-(tau * da_dtau + gl%ai_PR(nrsubst)) * dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) * del) / &
                    & (1.d0 - (2.d0 ** 0.5 - 1.d0) * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) * del)) / (2.d0 * 2.d0 ** 0.5 * gl%R_PR * &
                    & gl%tc(nrsubst) * gl%bi_PR(nrsubst))) * tau
            end if
            if (GETDER(5) == 1) then
                da_dtau = da_PR_dtau(gl,Temperature, nrsubst)
                d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, nrsubst)
                FNRDER(5) = (-((tau * d2a_dtau2 + 2.d0 * da_dtau) / (2.d0 * 2.d0 ** 0.5 * gl%bi_PR(nrsubst) * gl%R_PR * gl%tc(nrsubst)) * &
                    & dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst)) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * &
                    & del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst))))) * tau * tau
            end if
            if (GETDER(6) == 1) then
                da_dtau = da_PR_dtau(gl, Temperature, nrsubst)
                FNRDER(6) = ((gl%rhoc(nrsubst) * (tau * da_dtau + gl%ai_PR(nrsubst))) / (gl%R_PR * gl%tc(nrsubst) * ((del * gl%bi_PR(nrsubst) * &
                    & gl%rhoc(nrsubst))**2 - 2.d0 * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) - 1.d0))) * del * tau
            end if
            if (GETDER(7) == 1) then
                da_dtau = da_PR_dtau(gl, Temperature, nrsubst)
                d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, nrsubst)
                FNRDER(7) = ((gl%rhoc(nrsubst) * (tau * d2a_dtau2 + 2.d0 * da_dtau)) / (gl%R_PR * gl%tc(nrsubst) * ((del * gl%bi_PR(nrsubst) * &
                    & gl%rhoc(nrsubst))**2 - 2.d0 * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) - 1.d0))) * tau * tau * del
            end if
            if (GETDER(8) == 1) then
                FNRDER(8) = (2.d0 * (gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**3 - gl%ai_PR(nrsubst) * tau / &
                    & (2.d0 ** 0.5 * gl%R_PR * gl%tc(nrsubst) * gl%bi_PR(nrsubst)) * (((2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / &
                    & (1.d0 + (2.d0 ** 0.5 + 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**3 + ((2.d0 ** 0.5 - 1.d0) * &
                    & gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * gl%bi_PR(nrsubst) * del * gl%rhoc(nrsubst)))**3)) * &
                    & del * del * del
            end if
            if (GETDER(9) == 1) then
                d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, nrsubst)
                d3a_dtau3 = d3a_PR_dtau3(gl,Temperature, nrsubst)
                FNRDER(9) = (-((tau * d3a_dtau3 + 3.d0 * d2a_dtau2) / (2.d0 * 2.d0 ** 0.5 * gl%bi_PR(nrsubst) * gl%R_PR * gl%tc(nrsubst)) * &
                    & dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst)) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * &
                    & del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst))))) * tau * tau * tau
            end if
            if (GETDER(10) == 1) then
                da_dtau = da_PR_dtau(gl, Temperature, nrsubst)
                FNRDER(10) = (-(2.d0 * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst)**2 * (gl%ai_PR(nrsubst) + da_dtau * tau) * (del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) - 1.d0)) / &
                    & (gl%R_PR * gl%tc(nrsubst) * ((del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst))**2 - 2.d0 * del * gl%bi_PR(nrsubst) * gl%rhoc(nrsubst) - 1.d0)**2)) * &
                    &  tau * del * del
            end if

            !Stefan May 2014 (see Julian Hedtmann / Lars Hüttermann / Monika Thol)
        else If (gl%Eq_type(nrsubst) == 4) then  !LKP pure fluid equation is used

            !Hard coded "1" for every fluid replaced by "nrsubst". Andreas Jäger, March 2017
            !zcLKP=(pc(1)*1.d6)/(REQ(1)*rhored(1)*Tc(1))
            !wLKP=accen(1)
            gl%zcLKP=(gl%pc(nrsubst)*1.d6)/(gl%REQ(nrsubst)*gl%rhored(nrsubst)*gl%Tc(nrsubst))
            gl%wLKP=gl%accen(nrsubst)
            del2=del*del
            del3=del2*del
            del4=del*del3
            del5=del4*del
            zcLKP2=gl%zcLKP*gl%zcLKP
            zcLKP4=zcLKP2*zcLKP2
            zcLKP5=gl%zcLKP*zcLKP4
            zcLKP6=gl%zcLKP*zcLKP5
            zcLKP8=zcLKP6*zcLKP2
            tau2=tau*tau
            tau3=tau*tau2


            B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
            C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
            D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
            B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
            C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
            D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau



            if (GETDER(1) == 1) then

                !Ideal fluid contribution
                part_id = B_0 / gl%zcLKP * del + 0.5d0 * C_0 / zcLKP2 * del2 + 0.2d0 * D_0 / zcLKP5 * del5 -  &
                    & gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2)  &
                    & + gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_beta_0 + 1.d0)

                !Reference fluid contribution
                part_ref = B_ref / gl%zcLKP * del + 0.5d0 * C_ref / zcLKP2 * del2 + 0.2d0 * D_ref / zcLKP5 * del5 -  &
                    & gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2)  &
                    & + gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_beta_ref + 1.d0)

                FNRDER(1)=(1-gl%wLKP / gl%lkp_w_ref)*part_id + gl%wLKP / gl%lkp_w_ref*part_ref

            end if

            if (GETDER(2) == 1) then

                !Ideal fluid contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (1.d0 / del + B_0 / gl%zcLKP + C_0 * del / zcLKP2 + D_0 * del4 / zcLKP5 - gl%lkp_c4_0 * del * tau3 * help_id1 / zcLKP2 &
                    & + gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau3 * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                part_ref = gl%wLKP / gl%lkp_w_ref * (1.d0 / del + B_ref / gl%zcLKP + C_ref * del / zcLKP2 + D_ref * del4 / zcLKP5 - gl%lkp_c4_ref * del * tau3 * help_ref1 / zcLKP2 &
                    & + gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau3 * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
                FNRDER(2) = (part_id + part_ref - 1.d0 / del)*del

            end if

            if (GETDER(3) == 1) then

                !Ideal fluid contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
                help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (-1.d0 / del2 + C_0 / zcLKP2 + 4.d0 * D_0 * del3 / zcLKP5 &
                    & - gl%lkp_c4_0 * tau3 * help_id1 / zcLKP2 + 4.d0 * gl%lkp_c4_0 * tau3 * del2 * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                    & + gl%lkp_c4_0 * tau3 * gl%lkp_gamma_0 * help_id1 * help_id2 / zcLKP4 &
                    & - 2.d0 * gl%lkp_c4_0 * tau3 * del2 * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
                part_ref = gl%wLKP / gl%lkp_w_ref * (-1.d0 / del2 + C_ref / zcLKP2 + 4.d0 * D_ref * del3 / zcLKP5 &
                    & - gl%lkp_c4_ref * tau3 * help_ref1 / zcLKP2 + 4.d0 * gl%lkp_c4_ref * tau3 * del2 * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                    & + gl%lkp_c4_ref * tau3 * gl%lkp_gamma_ref * help_ref1 * help_ref2 / zcLKP4 &
                    & - 2.d0 * gl%lkp_c4_ref * tau3 * del2 * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
                FNRDER(3) = (part_id + part_ref + 1.d0 / del2)*del2

            end if


            if (GETDER(4) == 1) then

                DB_0_Dtau = -gl%lkp_b2_0 - 2.d0 * gl%lkp_b3_0 * tau - 3.d0 * gl%lkp_b4_0 * tau2
                DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
                DD_0_Dtau = gl%lkp_d2_0
                DB_ref_Dtau = -gl%lkp_b2_ref - 2.d0 * gl%lkp_b3_ref * tau - 3.d0 * gl%lkp_b4_ref * tau2
                DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
                DD_ref_Dtau = gl%lkp_d2_ref


                !Ideal fluid contribution
                part_id = DB_0_Dtau / gl%zcLKP * del + 0.5d0 * DC_0_Dtau / zcLKP2 * del2 + 0.2d0 * DD_0_Dtau / zcLKP5 * del5 &
                    & -1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                    & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
                !Reference fluid contribution
                part_ref = DB_ref_Dtau / gl%zcLKP * del + 0.5d0 * DC_ref_Dtau / zcLKP2 * del2 + 0.2d0 * DD_ref_Dtau / zcLKP5 * del5 &
                    & -1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                    & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_beta_ref + 1.d0)

                FNRDER(4)=((1-gl%wLKP / gl%lkp_w_ref)*part_id + gl%wLKP / gl%lkp_w_ref*part_ref) * tau

            end if

            if (GETDER(5) == 1) then

                D2B_0_Dtau2 = -2.d0 * gl%lkp_b3_0 - 6.d0 * gl%lkp_b4_0 * tau
                D2C_0_Dtau2 = 6.d0 * gl%lkp_c3_0 * tau
                D2B_ref_Dtau2 = -2.d0 * gl%lkp_b3_ref - 6.d0 * gl%lkp_b4_ref * tau
                D2C_ref_Dtau2 = 6.d0 * gl%lkp_c3_ref * tau

                !Ideal fluid contribution
                part_id = D2B_0_Dtau2 / gl%zcLKP * del + 0.5d0 * D2C_0_Dtau2 / zcLKP2 * del2 &
                    & -3.d0 * gl%lkp_c4_0 * tau / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                    & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_0 * tau / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
                !Reference fluid contribution
                part_ref = D2B_ref_Dtau2 / gl%zcLKP * del + 0.5d0 * D2C_ref_Dtau2 / zcLKP2 * del2 &
                    & -3.d0 * gl%lkp_c4_ref * tau / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                    & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_ref * tau / gl%lkp_gamma_ref * (gl%lkp_beta_ref + 1.d0)

                FNRDER(5)=((1-gl%wLKP / gl%lkp_w_ref)*part_id + gl%wLKP / gl%lkp_w_ref*part_ref) * tau2

            end if

            if (GETDER(6) == 1) then

                DB_0_Dtau = -gl%lkp_b2_0 - 2.d0 * gl%lkp_b3_0 * tau - 3.d0 * gl%lkp_b4_0 * tau2
                DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
                DD_0_Dtau = gl%lkp_d2_0
                DB_ref_Dtau = -gl%lkp_b2_ref - 2.d0 * gl%lkp_b3_ref * tau - 3.d0 * gl%lkp_b4_ref * tau2
                DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
                DD_ref_Dtau = gl%lkp_d2_ref

                !Ideal fluid contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (DB_0_Dtau / gl%zcLKP + DC_0_Dtau * del / zcLKP2 + DD_0_Dtau * del4 / zcLKP5 - 3.d0 * gl%lkp_c4_0 * del * tau2 * help_id1 / zcLKP2 &
                    & + 3.d0 * gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau2 * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                part_ref = gl%wLKP / gl%lkp_w_ref * (DB_ref_Dtau / gl%zcLKP + DC_ref_Dtau * del / zcLKP2 + DD_ref_Dtau * del4 / zcLKP5 - 3.d0 * gl%lkp_c4_ref * del * tau2 * help_ref1 / zcLKP2 &
                    & + 3.d0 * gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau2 * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
                FNRDER(6) = (part_id + part_ref)*tau*del

            end if

            if (GETDER(7) == 1) then

                D2B_0_Dtau2 = -2.d0 * gl%lkp_b3_0 - 6.d0 * gl%lkp_b4_0 * tau
                D2C_0_Dtau2 = 6.d0 * gl%lkp_c3_0 * tau
                D2B_ref_Dtau2 = -2.d0 * gl%lkp_b3_ref - 6.d0 * gl%lkp_b4_ref * tau
                D2C_ref_Dtau2 = 6.d0 * gl%lkp_c3_ref * tau

                !Ideal gas contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2/ zcLKP2)
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (D2B_0_Dtau2 / gl%zcLKP + D2C_0_Dtau2 * del / zcLKP2 - 6.d0 * gl%lkp_c4_0 * del * tau * help_id1 / zcLKP2 &
                    & + 6.d0 * gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                part_ref = gl%wLKP / gl%lkp_w_ref * (D2B_ref_Dtau2 / gl%zcLKP + D2C_ref_Dtau2 * del / zcLKP2 - 6.d0 * gl%lkp_c4_ref * del * tau * help_ref1 / zcLKP2 &
                    & + 6.d0 * gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
                FNRDER(7) = (part_id + part_ref)*del*tau2

            end if

            if (GETDER(8) == 1) then

                !Ideal fluid contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
                help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (2.d0 / del3 + 12.d0 * D_0 * del2 / zcLKP5 &
                    & - 12.d0 * gl%lkp_c4_0 * tau3 * del3 * gl%lkp_gamma_0 ** 2 * help_id1 / zcLKP6 &
                    & + 12.d0 * gl%lkp_c4_0 * tau3 * del * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                    & + 4.d0 * gl%lkp_c4_0 * tau3 * del3 * gl%lkp_gamma_0 ** 3 * help_id1 * help_id2 / zcLKP8 &
                    & - 6.d0 * gl%lkp_c4_0 * tau3 * del * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
                part_ref = gl%wLKP / gl%lkp_w_ref * (2.d0 / del3 + 12.d0 * D_ref * del2 / zcLKP5 &
                    & - 12.d0 * gl%lkp_c4_ref * tau3 * del3 * gl%lkp_gamma_ref ** 2 * help_ref1 / zcLKP6 &
                    & + 12.d0 * gl%lkp_c4_ref * tau3 * del * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                    & + 4.d0 * gl%lkp_c4_ref * tau3 * del3 * gl%lkp_gamma_ref ** 3 * help_ref1 * help_ref2 / zcLKP8 &
                    & - 6.d0 * gl%lkp_c4_ref * tau3 * del * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
                FNRDER(8) = (part_id + part_ref - 2.d0 / del3)*del3

            end if


            if (GETDER(9) == 1) then

                D3B_0_Dtau3 = -6.d0 * gl%lkp_b4_0
                D3C_0_Dtau3 = 6.d0 * gl%lkp_c3_0
                D3B_ref_Dtau3 = -6.d0 * gl%lkp_b4_ref
                D3C_ref_Dtau3 = 6.d0 * gl%lkp_c3_ref

                !Ideal fluid contribution
                part_id = D3B_0_Dtau3 / gl%zcLKP * del + 0.5d0 * D3C_0_Dtau3 / zcLKP2 * del2 &
                    & -3.d0 * gl%lkp_c4_0 / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                    & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_0 / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
                !Reference fluid contribution
                part_ref = D3B_ref_Dtau3 / gl%zcLKP * del + 0.5d0 * D3C_ref_Dtau3 / zcLKP2 * del2 &
                    & -3.d0 * gl%lkp_c4_ref / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                    & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_ref / gl%lkp_gamma_ref* (gl%lkp_beta_ref + 1.d0)

                FNRDER(9)=((1-gl%wLKP / gl%lkp_w_ref)*part_id + gl%wLKP / gl%lkp_w_ref*part_ref) * tau3

            end if


            if (GETDER(10) == 1) then

                DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
                DD_0_Dtau = gl%lkp_d2_0
                DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
                DD_ref_Dtau = gl%lkp_d2_ref

                !Ideal fluid contribution
                help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
                help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
                part_id = (1.d0 - gl%wLKP / gl%lkp_w_ref) * (DC_0_Dtau / zcLKP2 + 4.d0 * DD_0_Dtau * del3 / zcLKP5 &
                    & - 3.d0 * gl%lkp_c4_0 * tau2 * help_id1 / zcLKP2 + 12.d0 * gl%lkp_c4_0 * tau2 * del2 * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                    & + 3.d0 * gl%lkp_c4_0 * tau2 * gl%lkp_gamma_0 * help_id1 * help_id2 / zcLKP4 &
                    & - 6.d0 * gl%lkp_c4_0 * tau2 * del2 * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
                !Reference fluid contribution
                help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
                help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
                part_ref = gl%wLKP / gl%lkp_w_ref * (DC_ref_Dtau / zcLKP2 + 4.d0 * DD_ref_Dtau * del3 / zcLKP5 &
                    & - 3.d0 * gl%lkp_c4_ref * tau2 * help_ref1 / zcLKP2 + 12.d0 * gl%lkp_c4_ref * tau2 * del2 * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                    & + 3.d0 * gl%lkp_c4_ref * tau2 * gl%lkp_gamma_ref * help_ref1 * help_ref2 / zcLKP4 &
                    & - 6.d0 * gl%lkp_c4_ref * tau2 * del2 * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
                FNRDER(10) = (part_id + part_ref)*tau*del2

            end if

        elseif (gl%eq_type(nrsubst) == 6) then !PC-SAFT

            dens_orig = density
            Temp_orig = Temperature
            !Moni SAFT  Density = Density*N_A/1.d30
            !
            !if (GETDER(1) == 1) then
            !    continue
            !end if
            !
            !if (GETDER(2) == 1) then
            !
            !Andreas, Erik July 2019
            !Modification for using PC-SAFT in the multi-fluid mixture model, calculate the Temperature and density at which the pure fluid equation needs to be evaluated
            if (gl%ncomp > 1 .and. nrsubst == 0) then
                Temperature = gl%tc(nrsubst) / tau
                density = del * gl%rhoc(nrsubst)
            end if
            !eventuell nrsubst als Übergabe
            call ARDERIVS(gl,Temperature, Density, GETDER, FNRDER, nrsubst)

            !end if
            !
            !
            !if (GETDER(3) == 1) then
            !    !eventuell nrsubst als Übergabe
            !    call ARDERIVS(Temperature, Density, GETDER, FNRDER)
            !end if

            density=dens_orig
            Temperature = Temp_orig

            !AGA 8 -  thermodynamic properties of natural gas and related gas
        else if (gl%Eq_type(nrsubst) == 7) then  !aga8 pure fluid equation is used

            call mixture_parameters(gl,Temperature,Density, B_0)
            call residualpart_AGA8(Temperature, Density, GETDER, FNRDER)

        end if

        !! saving the derivatives and belonging inputs
        !gl%TEMP_FNR_OLD = TEMPERATURE
        !gl%DENS_FNR_OLD = DENSITY
        !gl%tredmix_old = gl%tredmix
        !gl%rhoredmix_old = gl%rhoredmix
        !gl%nrsubst_old = nrsubst
        !gl%Eq_type_old = gl%Eq_type
        !gl%FNR_OLD = FNR
        !gl%FNRDER_OLD = FNRDER
        !gl%delpi_old = delpi
        !gl%reg_term_old = reg_term
        !gl%gauss_term_old = gauss_term
        !gl%nacalc_old = nacalc
        !gl%deleps_old = deleps
        !gl%taugam_old = taugam

        !*****************************************************

    end subroutine FNRDERIVS
    !**************************************************************************







    !**************************************************************************
    !           --------------------------------------------------
    !           Routine for the calculation of the mixture
    !           derivatives of the residual part of the
    !           Helmholtz energy
    !
    !           J. Gernert, Denmark, 09.2009
    !           --------------------------------------------------
    !**************************************************************************

    !**************************************************************************
    subroutine MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY FOR MIXTURES
    ! THE ROUTINE CALLS THE ROUTINE 'FNRDERIVS' FOR THE SINGLE FLUIDS AND THE ROUTINE 'DEPFUNCFNR' FOR
    ! THE BINARY DEPARTURE FUNCTION AT THE REDUCING PARAMETERS FOR T_RED AND RHO_RED. IT RETURNS THE
    ! RESIDUAL MIXTURE HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    ! GETDER      - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X
    !                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
    !                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !
    ! OUTPUT PARAMETERS:
    ! MIXDERFNR   - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !                AS INDECATED IN "GETDER"
    !--------------------------------------------------------------------------------------------------










    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    integer, dimension(nderivs):: GETDER            ! array specifier to indicate, which derivative is needed
    double precision, dimension(nderivs)::MIXDERFNR  ! array with the computed values for the derivatives
    double precision, dimension(nderivs)::DEPFUNCDER              ! array with the computed values for the derivatives of the departure function
    double precision, dimension(nderivs)::DEPFUNCBIN              ! array with the computed values for the derivatives of the binary specific departure function
    double precision, dimension(nderivs):: FNRDER                 ! single fluid residual helmholtz energy return vector
    integer:: i, j, k, m,o                                   ! index variables
    integer:: FLD1, FLD2                                  ! Fluid numbers of the fluids that are part of the binary mixtures with departure functions
    double precision:: xi, xj

    !----------------------------------------------------
    !Variables needed for the SRK EOS
    !Andreas March 2012
    double precision:: tau, delta

    !SRK supportive variables. Lars Hüttermann, 14.02.2013
    double precision :: dT_dTau
    double precision :: d2T_dTau2
    double precision :: d3T_dTau3
    double precision :: da_dTau
    double precision :: d2a_dTau2
    double precision :: d3a_dTau3
    double precision :: da_SRK_dtau
    double precision :: d2a_SRK_dtau2
    double precision :: d3a_SRK_dtau3

    !PR supportive variables. Andreas Jäger, Nov. 2015
    double precision :: da_PR_dtau
    double precision :: d2a_PR_dtau2
    double precision :: d3a_PR_dtau3

    !LKP supportive variables. Stefan und Andreas (see Lars Hüttermann) June 14
    double precision :: B_0
    double precision :: C_0
    double precision :: D_0
    double precision :: B_ref
    double precision :: C_ref
    double precision :: D_ref
    double precision :: del
    double precision :: part_id
    double precision :: part_ref
    double precision :: help_id, help_id1, help_id2
    double precision :: help_ref, help_ref1, help_ref2
    double precision :: zcLKP2, zcLKP4, zcLKP5,zcLKP6, zcLKP8, tau2, tau3, del2, del3, del4, del5
    double precision :: DB_0_Dtau, DC_0_Dtau, DD_0_Dtau, DB_ref_Dtau, DC_ref_Dtau, DD_ref_Dtau, D2B_0_Dtau2, D2C_0_Dtau2, D2B_ref_Dtau2, D2C_ref_Dtau2
    double precision :: D3B_0_Dtau3, D3C_0_Dtau3, D3B_ref_Dtau3, D3C_ref_Dtau3

    double precision:: part_id0, part_ref0, delr, integral
    double precision:: test

    !saving derivatives
    logical :: complete
    integer, dimension(nderivs):: GETDERS             ! working array specifier to indicate, which derivative is still needed

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision, dimension(nderivs):: alpha_r_ij     !Combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij)) with all delta and tau derivatives
    double precision, dimension(30,nderivs):: alpha_r_i         !derivatives of all reduced residual Helmholtz energies of the pure fluids

    !New variables for excess based departure function
    double precision, dimension(nderivs)::DEPFUNCDER_GE              ! array with the computed values for the derivatives of the departure function

    !New variables for one-fluid model
    double precision, dimension(100):: reg_termi, reg_termj
    double precision:: deltai, deltaj, taui, tauj

    integer :: nrsubst

    dT_dTau = 0.d0
    d2T_dTau2 = 0.d0
    d3T_dTau3 = 0.d0
    da_dTau = 0.d0
    d2a_dTau2 = 0.d0
    d3a_dTau3 = 0.d0
    !----------------------------------------------------

    MIXDERFNR = 0.d0
    DEPFUNCDER = 0.d0
    DEPFUNCBIN = 0.d0



    GETDERS = GETDER

    !***************************************************************************************************************

    complete = .true.
    ! array specifier to indicate, which derivative is needed

    !Problems with the SRK here! tredmix and rhoredmix are always "1" for the SRK, so this cannot be used as
    !decision here. Preliminary solution: Exclude the SRK from saving here. This needs to be discussed!!
    !Andreas, September 2015
    !if ( (dabs(TEMPERATURE-TEMP_FNR_MIX_OLD) < 1.d-14) .and. (dabs(DENSITY-DENS_FNR_MIX_OLD) < 1.d-14) &
    !        .and. (dabs(tredmix-tredmix_MIX_old) < 1.d-14) .and. (dabs(rhoredmix-rhoredmix_MIX_OLD) < 1.d-14) ) then
    !
    !if ( (dabs(TEMPERATURE-gl%TEMP_FNR_MIX_OLD) < 1.d-14) .and. (dabs(DENSITY-gl%DENS_FNR_MIX_OLD) < 1.d-14) &
    !        .and. (dabs(gl%tredmix-gl%tredmix_MIX_old) < 1.d-14) .and. (dabs(gl%rhoredmix-gl%rhoredmix_MIX_OLD) < 1.d-14) &
    !        .and. (gl%mix_type == 1)    ) then
    !
    !        do o = 1, nderivs
    !            if ( (GETDERS(o) == 1) .and. (dabs(gl%FNRDER_MIX_OLD(o)) > 1.d-14) ) then
    !                MIXDERFNR(o) = gl%FNRDER_MIX_OLD (o)
    !                GETDERS(o) = 0
    !                gl%za_fnm = gl%za_fnm + 1
    !            else if (GETDERS(o) == 1) then
    !                complete = .false.
    !            end if
    !        end do
    !
    !    if (complete) then
    !        gl%zaehler_mix = gl%zaehler_mix +1
    !        return
    !    else
    !        gl%uncomp = gl%uncomp +1
    !    end if
    !
    !end if



    If (gl%Mix_type == 1) then     !Lorentz Berthelot or modified mixing rules used

        !sum of the residual Helmholtz energies and its derivatives of the pure fluids
        do i = 1, gl%NCOMP
            FNRDER = 0.D0
            call FNRDERIVS(gl,TEMPERATURE, DENSITY, GETDERS, FNRDER, i) ! call the calculation routine for the derivatives of the fluid i
            do j = 1, nderivs
                if (GETDERS(j) == 1) then     ! check which derivative is needed
                    MIXDERFNR(j) = MIXDERFNR(j) + gl%MOLFRACTIONS(i)*FNRDER(j) ! sum of derivatives over all pure fluids, multiplied with the respective mole fraction
                end if
            end do
        end do

        !sum of the binary departure functions
        do FLD1 = 1, gl%NCOMP - 1               !
            do FLD2 = FLD1 + 1, gl%NCOMP
                if (dabs(gl%Fij(FLD1, FLD2)) > 1.d-14) then    ! check if departure function exists for this binary combination
                    call DEPFUNCFNR (gl,TEMPERATURE, DENSITY, GETDERS, DEPFUNCBIN, FLD1, FLD2)   ! calls the routine for the calculation of the binary specific departure function and its derivatives
                    xi = gl%MOLFRACTIONS(FLD1)
                    xj = gl%MOLFRACTIONS(FLD2)
                    do k = 1, 10                      ! loop through all derivatives returned from DEPFUNCFNR
                        DEPFUNCDER(k) = DEPFUNCDER(k) + xi*xj*gl%Fij(FLD1, FLD2)*DEPFUNCBIN(k)     ! sum over all binary departure functions, multiplied with the mole fractions and the weighing factor Fij
                    end do
                    DEPFUNCBIN = 0.d0
                end if
            end do
        end do

        ! add the two parts of the residual Helmholtz enthalpy for mixtures
        do m = 1, nderivs
            MIXDERFNR(m) = MIXDERFNR(m) + DEPFUNCDER(m)
        end do


    elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then   !SRK with SRK mixing rules used

        !NEW PART FOR THE RESIDUAL HELMHOLTZ-ENERGY CALCULATED FROM THE SRK!!
        !Andreas March 2012

        !New routine for the calculation of the residual Helmholtz energy and its derivatives (up to 4th) wrt tau and delta for cubic EOSs
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_HIGHER_CUBIC (gl,Temperature, Density, GETDERS, MIXDERFNR)



        !OLD CODE REPLACED WITH NEW CODE FOR CUBIC EQUATIONS OF STATE (SRK AND PENG-ROBINSON)
        !Can be deleted if no problems occur
        !Andreas Jäger, March 2017

        !delta = Density / rhored_SRK
        !tau = Tred_SRK / Temperature
        !
        !
        !if (GETDERS(1) == 1) then
        !    MIXDERFNR(1) = -dlog(1.D0-b_SRK*rhored_SRK*delta)-a_SRK*tau/R_SRK/tred_SRK/b_SRK*dlog(1.D0+b_SRK*rhored_SRK*delta)
        !end if
        !
        !if (GETDERS(2) == 1) then
        !    MIXDERFNR(2) = (b_SRK * rhored_SRK / (1.D0 - b_SRK * rhored_SRK * delta) - &
        !                 & a_SRK * rhored_SRK * tau / (R_SRK * tred_SRK*(1.D0+b_SRK*rhored_SRK*delta)))*delta
        !end if
        !if (GETDERS(3) == 1) then
        !    MIXDERFNR(3) = delta*delta*(((b_SRK*rhored_SRK)/(1-b_SRK*rhored_SRK*delta))**2 + &
        !                 & a_SRK * tau / R_SRK / tred_SRK * b_SRK * rhored_SRK**2.d0 / (1.d0+b_SRK*rhored_SRK*delta)**2)
        !end if
        !!if (GETDERS(4) == 1) then
        !!    daSRK_dtau = da_SRK_dtau(gl,Temperature)
        !!    dT_dtau = - Temperature / tau
        !!    MIXDERFNR(4) = -(1.D0 / R_SRK / b_SRK * dlog(1.D0 + b_SRK*rhored_SRK*delta)* &
        !!                & (Temperature*daSRK_dtau-dT_dtau*a_SRK) / Temperature**2.d0) * tau
        !!end if
        !if (GETDERS(4) == 1) then    !Lars H., 14.02.2013
        !    da_dTau = da_SRK_dtau(gl,Temperature, 0)
        !    dT_dTau = -Tred_SRK / (tau * tau) !/ Tau ** 2
        !    MIXDERFNR(4) = -1.d0 / (R_SRK * b_SRK) * dlog(1.d0 +  b_SRK * rhored_SRK * delta) * &
        !              & (da_dTau * (Tred_SRK / Tau) - a_SRK * dT_dTau) / (Tred_SRK / Tau) ** 2 * Tau
        !end if
        !if (GETDERS(5) == 1) then    !Lars H., 14.02.2013
        !    dT_dTau = -Tred_SRK / (tau * tau) !/ Tau ** 2
        !    d2T_dTau2 = 2.d0 * Tred_SRK / (tau * tau * tau) !/ Tau ** 3.d0
        !    da_dTau = da_SRK_dtau(gl,Temperature, 0)
        !    d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, 0)
        !    MIXDERFNR(5) = -1.d0 / (R_SRK * b_SRK) * dlog(1.d0 +  b_SRK * rhored_SRK * delta) * &
        !              & (d2a_dTau2 / (Tred_SRK / Tau) - a_SRK * d2T_dTau2 / (Tred_SRK / Tau) ** 2 - &
        !              & 2.d0 * dT_dTau * da_dTau / (Tred_SRK / Tau) ** 2 + &
        !              & 2.d0 * a_SRK * dT_dTau ** 2 / (Tred_SRK / Tau) ** 3) * tau * tau
        !end if
        !if (GETDERS(6) == 1) then
        !    da_dtau = da_SRK_dtau(gl,Temperature, 0)
        !    dT_dtau = - Temperature / tau
        !    MIXDERFNR(6) = - (rhored_SRK / R_SRK / (1.D0 + b_SRK*rhored_SRK*delta) * &
        !                & (Temperature*da_dtau-dT_dtau*a_SRK) / Temperature**2.d0) * delta * tau
        !end if
        !if (GETDERS(7) == 1) then    !Lars H., 14.02.2013
        !    dT_dTau = -tred_SRK / (tau * tau) !/ Tau ** 2
        !    d2T_dTau2 = 2.d0 * tred_SRK / (tau * tau * tau) !/ Tau ** 3.d0
        !    da_dTau = da_SRK_dtau(gl,Temperature, 0)
        !    d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, 0)
        !    MIXDERFNR(7) = -rhored_SRK / (R_SRK * (1.d0 + b_SRK * rhored_SRK * delta)) * &
        !              & (d2a_dTau2 / (tred_SRK / Tau) - a_SRK * d2T_dTau2 / (tred_SRK / Tau) ** 2 - &
        !              & 2.d0 * dT_dTau * da_dTau / (tred_SRK / Tau) ** 2 + &
        !              & 2.d0 * a_SRK * dT_dTau ** 2 / (tred_SRK / Tau) ** 3) * &
        !              & tau * tau * delta
        !end if
        !if (GETDERS(8) == 1) then
        !    MIXDERFNR(8) = delta*delta*delta*(2.D0*(((b_SRK*rhored_SRK)/(1-b_SRK*rhored_SRK*delta))**3 - &
        !         & a_SRK * tau / R_SRK / tred_SRK * b_SRK**2.d0 * rhored_SRK**3.d0 / (1.d0+b_SRK*rhored_SRK*delta)**3))
        !end if
        !if (GETDERS(9) == 1) then    !Lars H., 14.02.2013
        !    dT_dTau = -TRed_SRK / (tau * tau)!/ Tau ** 2
        !    d2T_dTau2 = 2.d0 * TRed_SRK / (tau * tau * tau) !Tau ** 3.d0
        !    d3T_dTau3 = - 6.d0 * TRed_SRK / (tau * tau * tau * tau) !/ Tau ** 4.d0
        !    da_dTau = da_SRK_dtau(gl,Temperature, 0)
        !    d2a_dTau2 = d2a_SRK_dtau2(gl,Temperature, 0)
        !    d3a_dTau3 = d3a_SRK_dtau3(gl,Temperature, 0)
        !    MIXDERFNR(9) = -1.d0 / (R_SRK * b_SRK) * dlog(1.d0 +  b_SRK * rhored_SRK * delta) * &
        !              & (d3a_dTau3 / (TRed_SRK / Tau) - a_SRK * d3T_dTau3 / (TRed_SRK / Tau) ** 2 - &
        !              & 3.d0 * d2T_dTau2 * da_dTau / (TRed_SRK / Tau) ** 2 - &
        !              & 3.d0 * dT_dTau * d2a_dTau2 / (TRed_SRK / Tau) ** 2 + &
        !              & 6.d0 * da_dTau * dT_dTau ** 2 / (TRed_SRK / Tau) ** 3 + &
        !              & 6.d0 * a_SRK * d2T_dTau2 * dT_dTau / (TRed_SRK / Tau) ** 3 - &
        !              & 6.d0 * a_SRK * dT_dTau ** 3 / (TRed_SRK / Tau) ** 4) * tau * tau * tau
        !end if
        !if (GETDERS(10) == 1) then    !Lars H., 14.02.2013
        !    dT_dTau = -TRed_SRK / (tau * tau) !/ Tau ** 2
        !    da_dTau = da_SRK_dtau(gl,Temperature, 0)
        !    MIXDERFNR(10) = rhored_SRK ** 2 * b_SRK / (R_SRK * (1.d0 + b_SRK * rhored_SRK * delta) ** 2) * &
        !               & (da_dTau / (TRed_SRK / Tau) - a_SRK * dT_dTau / (TRed_SRK / Tau) ** 2) * tau * delta * delta
        !end if

    elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then   !PR with PR mixing rules used

        !THE RESIDUAL HELMHOLTZ-ENERGY CALCULATED FROM THE Peng-Robinson EOS!!
        !Andreas November 2015

        !New routine for the calculation of the residual Helmholtz energy and its derivatives (up to 4th) wrt tau and delta for cubic EOSs
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_HIGHER_CUBIC (gl,Temperature, Density, GETDERS, MIXDERFNR)



        !OLD CODE REPLACED WITH NEW CODE FOR CUBIC EQUATIONS OF STATE (SRK AND PENG-ROBINSON)
        !Can be deleted if no problems occur
        !Andreas Jäger, March 2017

        !delta = Density / rhored_PR
        !tau = Tred_PR / Temperature
        !
        !if (GETDER(1) == 1) then
        !    MIXDERFNR(1) = -dlog(1.D0 - b_PR * rhored_PR * delta) - a_PR * tau * (1.d0 / &
        !                & (2.d0 * 2.d0 ** 0.5 * b_PR * R_PR * Tred_PR)) * dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * &
        !                & delta * b_PR * rhored_PR) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * delta * b_PR * rhored_PR))
        !end if
        !if (GETDER(2) == 1) then
        !    MIXDERFNR(2) = (b_PR * rhored_PR / (1.d0 - b_PR * delta * rhored_PR) - a_PR * tau / (2.d0 * &
        !                & 2.d0 ** 0.5 * R_PR * Tred_PR * b_PR) * ((2.d0 ** 0.5 + 1.d0) * b_PR * rhored_PR / &
        !                & (1.d0 + (2.d0 ** 0.5 + 1.d0) * b_PR * delta * rhored_PR) + (2.d0 ** 0.5 - 1.d0) * &
        !                & b_PR * rhored_PR / (1.d0 - (2.d0 ** 0.5 - 1.d0) * b_PR * delta * rhored_PR))) * delta
        !end if
        !if (GETDER(3) == 1) then
        !    MIXDERFNR(3) = ((b_PR * rhored_PR / (1.d0 - b_PR * delta * rhored_PR))**2 - a_PR * tau / (2.d0 * &
        !                & 2.d0 ** 0.5 * R_PR * Tred_PR * b_PR) * (-((2.d0 ** 0.5 + 1.d0) * b_PR * rhored_PR / &
        !                & (1.d0 + (2.d0 ** 0.5 + 1.d0) * b_PR * delta * rhored_PR))**2 + ((2.d0 ** 0.5 - 1.d0) * &
        !                & b_PR * rhored_PR / (1.d0 - (2.d0 ** 0.5 - 1.d0) * b_PR * delta * rhored_PR))**2)) * &
        !                & delta * delta
        !end if
        !if (GETDER(4) == 1) then
        !    da_dtau = da_PR_dtau(gl, Temperature, 0)
        !    MIXDERFNR(4) = (-(tau * da_dtau + a_PR) * dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * b_PR * rhored_PR * delta) / &
        !                & (1.d0 - (2.d0 ** 0.5 - 1.d0) * b_PR * rhored_PR * delta)) / (2.d0 * 2.d0 ** 0.5 * R_PR * &
        !                & Tred_PR * b_PR)) * tau
        !end if
        !if (GETDER(5) == 1) then
        !    da_dtau = da_PR_dtau(gl, Temperature, 0)
        !    d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, 0)
        !    MIXDERFNR(5) = (-((tau * d2a_dtau2 + 2.d0 * da_dtau) / (2.d0 * 2.d0 ** 0.5 * b_PR * R_PR * Tred_PR) * &
        !                & dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * delta * b_PR * rhored_PR) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * &
        !                & delta * b_PR * rhored_PR)))) * tau * tau
        !end if
        !if (GETDER(6) == 1) then
        !    da_dtau = da_PR_dtau(gl, Temperature, 0)
        !    MIXDERFNR(6) = ((rhored_PR * (tau * da_dtau + a_PR)) / (R_PR * Tred_PR * ((delta * b_PR * &
        !                & rhored_PR)**2 - 2.d0 * delta * b_PR * rhored_PR - 1.d0))) * delta * tau
        !end if
        !if (GETDER(7) == 1) then
        !    da_dtau = da_PR_dtau(gl, Temperature, 0)
        !    d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, 0)
        !    MIXDERFNR(7) = ((rhored_PR * (tau * d2a_dtau2 + 2.d0 * da_dtau)) / (R_PR * Tred_PR * ((delta * b_PR * &
        !                & rhored_PR)**2 - 2.d0 * delta * b_PR * rhored_PR - 1.d0))) * tau * tau * delta
        !end if
        !if (GETDER(8) == 1) then
        !    MIXDERFNR(8) = (2.d0 * (b_PR * rhored_PR / (1.d0 - b_PR * delta * rhored_PR))**3 - a_PR * tau / &
        !                & (2.d0 ** 0.5 * R_PR * Tred_PR * b_PR) * (((2.d0 ** 0.5 + 1.d0) * b_PR * rhored_PR / &
        !                & (1.d0 + (2.d0 ** 0.5 + 1.d0) * b_PR * delta * rhored_PR))**3 + ((2.d0 ** 0.5 - 1.d0) * &
        !                & b_PR * rhored_PR / (1.d0 - (2.d0 ** 0.5 - 1.d0) * b_PR * delta * rhored_PR))**3)) * &
        !                & delta * delta * delta
        !end if
        !if (GETDER(9) == 1) then
        !    d2a_dtau2 = d2a_PR_dtau2(gl,Temperature, 0)
        !    d3a_dtau3 = d3a_PR_dtau3(gl,Temperature, 0)
        !    MIXDERFNR(9) = (-((tau * d3a_dtau3 + 3.d0 * d2a_dtau2) / (2.d0 * 2.d0 ** 0.5 * b_PR * R_PR * Tred_PR) * &
        !                & dlog((1.d0 + (2.d0 ** 0.5 + 1.d0) * delta * b_PR * rhored_PR) / (1.d0 - (2.d0 ** 0.5 - 1.d0) * &
        !                & delta * b_PR * rhored_PR)))) * tau * tau * tau
        !end if
        !if (GETDER(10) == 1) then
        !    da_dtau = da_PR_dtau(gl, Temperature, 0)
        !    MIXDERFNR(10) = (-(2.d0 * b_PR * rhored_PR**2.d0 * (a_PR + da_dtau * tau) * (delta * b_PR * rhored_PR - 1.d0)) / &
        !                & (R_PR * Tred_PR * ((delta * b_PR * rhored_PR)**2 - 2.d0 * delta * b_PR * rhored_PR - 1.d0)**2)) * &
        !                &  tau * delta * delta
        !end if

    elseif (gl%Mix_type  == 4) then    !LKP with LKP mixing rules used

        !edited by Stefan and Monika 10/2014 (alpha_r and tau-derivatives)

        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del3=del2*del
        del4=del*del3
        del5=del4*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP4=zcLKP2*zcLKP2
        zcLKP5=gl%zcLKP*zcLKP4
        zcLKP6=gl%zcLKP*zcLKP5
        zcLKP8=zcLKP6*zcLKP2
        tau2=tau*tau
        tau3=tau*tau2


        B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
        C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
        D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
        B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
        C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
        D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau


        if (GETDER(1) == 1) then

            !Ideal fluid contribution
            part_id = B_0 / gl%zcLKP * del + 0.5d0 * C_0 / zcLKP2 * del2 + 0.2d0 * D_0 / zcLKP5 * del5 -  &
                & gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2)  &
                & + gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_beta_0 + 1.d0)

            !Reference fluid contribution
            part_ref = B_ref / gl%zcLKP * del + 0.5d0 * C_ref / zcLKP2 * del2 + 0.2d0 * D_ref / zcLKP5 * del5 -  &
                & gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2)  &
                & + gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_beta_ref + 1.d0)

            MIXDERFNR(1)=(1-gl%accenLKPMix / gl%lkp_w_ref)*part_id + gl%accenLKPMix / gl%lkp_w_ref*part_ref
        end if

        if (GETDER(2) == 1) then

            !Ideal fluid contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (1.d0 / del + B_0 / gl%zcLKP + C_0 * del / zcLKP2 + D_0 * del4 / zcLKP5 - gl%lkp_c4_0 * del * tau3 * help_id1 / zcLKP2 &
                & + gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau3 * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (1.d0 / del + B_ref / gl%zcLKP + C_ref * del / zcLKP2 + D_ref * del4 / zcLKP5 - gl%lkp_c4_ref * del * tau3 * help_ref1 / zcLKP2 &
                & + gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau3 * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
            MIXDERFNR(2) = (part_id + part_ref - 1.d0 / del) *del

            !!Modified Stefan & Monika 2.12.14
            !!Ideal fluid contribution
            !help_id1 = dexp(-lkp_gamma_0 * del2 / zcLKP2)
            !part_id =  B_0 / zcLKP + C_0 * del / zcLKP2 + D_0 * del4 / zcLKP5 + lkp_c4_0 * del * tau3 * help_id1 / zcLKP2 *(lkp_beta_0 + lkp_gamma_0 * del2/zcLKP2)
            !!Reference fluid contribution
            !help_ref1 = dexp(-lkp_gamma_ref * del2 / zcLKP2)
            !part_ref = B_ref / zcLKP + C_ref * del / zcLKP2 + D_ref * del4 / zcLKP5 + lkp_c4_ref * del * tau3 * help_ref1 / zcLKP2 *(lkp_beta_ref + lkp_gamma_ref * del2/zcLKP2)
            !MIXDERFNR(2) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref) *del


        end if

        if (GETDER(3) == 1) then

            !Ideal fluid contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (-1.d0 / del2 + C_0 / zcLKP2 + 4.d0 * D_0 * del3 / zcLKP5 &
                & - gl%lkp_c4_0 * tau3 * help_id1 / zcLKP2 + 4.d0 * gl%lkp_c4_0 * tau3 * del2 * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                & + gl%lkp_c4_0 * tau3 * gl%lkp_gamma_0 * help_id1 * help_id2 / zcLKP4 &
                & - 2.d0 * gl%lkp_c4_0 * tau3 * del2 * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (-1.d0 / del2 + C_ref / zcLKP2 + 4.d0 * D_ref * del3 / zcLKP5 &
                & - gl%lkp_c4_ref * tau3 * help_ref1 / zcLKP2 + 4.d0 * gl%lkp_c4_ref * tau3 * del2 * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                & + gl%lkp_c4_ref * tau3 * gl%lkp_gamma_ref * help_ref1 * help_ref2 / zcLKP4 &
                & - 2.d0 * gl%lkp_c4_ref * tau3 * del2 * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
            MIXDERFNR(3) = (part_id + part_ref + 1.d0 / del2)*del2

            !!Modified Stefan & Monika 2.12.14
            !!Ideal fluid contribution
            !help_id = lkp_gamma_0 * del2 / zcLKP2
            !part_id = C_0 / zcLKP2 + 4.d0 * D_0 * del3 / zcLKP5 &
            !        & + lkp_c4_0 * tau3 / zcLKP2 * dexp(-help_id)*(lkp_beta_0 + 3.d0 * help_id - 2.d0 * help_id / del * (lkp_beta_0 + help_id))
            !!Reference fluid contribution
            !help_ref = lkp_gamma_ref * del2 / zcLKP2
            !part_ref = C_ref / zcLKP2 + 4.d0 * D_ref * del3 / zcLKP5 &
            !        & + lkp_c4_ref * tau3 / zcLKP2 * dexp(-help_ref)*(lkp_beta_ref + 3.d0 * help_ref - 2.d0 * help_ref / del * (lkp_beta_ref + help_ref))
            !MIXDERFNR(3) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref)*del2

        end if

        if (GETDER(4) == 1) then

            DB_0_Dtau = -gl%lkp_b2_0 - 2.d0 * gl%lkp_b3_0 * tau - 3.d0 * gl%lkp_b4_0 * tau2
            DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
            DD_0_Dtau = gl%lkp_d2_0
            DB_ref_Dtau = -gl%lkp_b2_ref - 2.d0 * gl%lkp_b3_ref * tau - 3.d0 * gl%lkp_b4_ref * tau2
            DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
            DD_ref_Dtau = gl%lkp_d2_ref


            !Ideal fluid contribution
            part_id = DB_0_Dtau / gl%zcLKP * del + 0.5d0 * DC_0_Dtau / zcLKP2 * del2 + 0.2d0 * DD_0_Dtau / zcLKP5 * del5 &
                & -1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
            !Reference fluid contribution
            part_ref = DB_ref_Dtau / gl%zcLKP * del + 0.5d0 * DC_ref_Dtau / zcLKP2 * del2 + 0.2d0 * DD_ref_Dtau / zcLKP5 * del5 &
                & -1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_beta_ref + 1.d0)

            MIXDERFNR(4)=((1-gl%accenLKPMix / gl%lkp_w_ref)*part_id + gl%accenLKPMix / gl%lkp_w_ref*part_ref) * tau

        end if

        if (GETDER(5) == 1) then

            D2B_0_Dtau2 = -2.d0 * gl%lkp_b3_0 - 6.d0 * gl%lkp_b4_0 * tau
            D2C_0_Dtau2 = 6.d0 * gl%lkp_c3_0 * tau
            D2B_ref_Dtau2 = -2.d0 * gl%lkp_b3_ref - 6.d0 * gl%lkp_b4_ref * tau
            D2C_ref_Dtau2 = 6.d0 * gl%lkp_c3_ref * tau

            !Ideal fluid contribution
            part_id = D2B_0_Dtau2 / gl%zcLKP * del + 0.5d0 * D2C_0_Dtau2 / zcLKP2 * del2 &
                & -3.d0 * gl%lkp_c4_0 * tau / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_0 * tau / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
            !Reference fluid contribution
            part_ref = D2B_ref_Dtau2 / gl%zcLKP * del + 0.5d0 * D2C_ref_Dtau2 / zcLKP2 * del2 &
                & -3.d0 * gl%lkp_c4_ref * tau / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_ref * tau / gl%lkp_gamma_ref * (gl%lkp_beta_ref + 1.d0)

            MIXDERFNR(5)=((1-gl%accenLKPMix / gl%lkp_w_ref)*part_id + gl%accenLKPMix / gl%lkp_w_ref*part_ref) * tau2

        end if

        if (GETDER(6) == 1) then

            DB_0_Dtau = -gl%lkp_b2_0 - 2.d0 * gl%lkp_b3_0 * tau - 3.d0 * gl%lkp_b4_0 * tau2
            DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
            DD_0_Dtau = gl%lkp_d2_0
            DB_ref_Dtau = -gl%lkp_b2_ref - 2.d0 * gl%lkp_b3_ref * tau - 3.d0 * gl%lkp_b4_ref * tau2
            DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
            DD_ref_Dtau = gl%lkp_d2_ref
            !
            !Ideal fluid contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (DB_0_Dtau / gl%zcLKP + DC_0_Dtau * del / zcLKP2 + DD_0_Dtau * del4 / zcLKP5 - 3.d0 * gl%lkp_c4_0 * del * tau2 * help_id1 / zcLKP2 &
                & + 3.d0 * gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau2 * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (DB_ref_Dtau / gl%zcLKP + DC_ref_Dtau * del / zcLKP2 + DD_ref_Dtau * del4 / zcLKP5 - 3.d0 * gl%lkp_c4_ref * del * tau2 * help_ref1 / zcLKP2 &
                & + 3.d0 * gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau2 * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
            MIXDERFNR(6) = (part_id + part_ref)*tau*del

            !!Modified Stefan & Monika 2.12.14
            !!Ideal fluid contribution
            !part_id = DB_0_Dtau / zcLKP + DC_0_Dtau * del / zcLKP2 + DD_0_Dtau * del4 / zcLKP5 + 3.d0 * lkp_c4_0 * del * tau2 / zcLKP2 &
            !     & * (lkp_beta_0 + lkp_gamma_0 * del2 / zcLKP2) * dexp(-lkp_gamma_0 * del2 / zcLKP2)
            !!Reference fluid contribution
            !part_ref = DB_ref_Dtau / zcLKP + DC_ref_Dtau * del / zcLKP2 + DD_ref_Dtau * del4 / zcLKP5 + 3.d0 * lkp_c4_ref * del * tau2 / zcLKP2 &
            !     & * (lkp_beta_ref + lkp_gamma_ref * del2 / zcLKP2) * dexp(-lkp_gamma_ref * del2 / zcLKP2)
            !MIXDERFNR(6) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref)*tau*del

        end if

        if (GETDER(7) == 1) then

            D2B_0_Dtau2 = -2.d0 * gl%lkp_b3_0 - 6.d0 * gl%lkp_b4_0 * tau
            D2C_0_Dtau2 = 6.d0 * gl%lkp_c3_0 * tau
            D2B_ref_Dtau2 = -2.d0 * gl%lkp_b3_ref - 6.d0 * gl%lkp_b4_ref * tau
            D2C_ref_Dtau2 = 6.d0 * gl%lkp_c3_ref * tau

            !Ideal gas contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2/ zcLKP2)
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (D2B_0_Dtau2 / gl%zcLKP + D2C_0_Dtau2 * del / zcLKP2 - 6.d0 * gl%lkp_c4_0 * del * tau * help_id1 / zcLKP2 &
                & + 6.d0 * gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau * help_id1 / zcLKP4 * (del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0))
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (D2B_ref_Dtau2 / gl%zcLKP + D2C_ref_Dtau2 * del / zcLKP2 - 6.d0 * gl%lkp_c4_ref * del * tau * help_ref1 / zcLKP2 &
                & + 6.d0 * gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau * help_ref1 / zcLKP4 * (del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref))
            MIXDERFNR(7) = (part_id + part_ref)*del*tau2

            !Modified Stefan & Monika 2.12.14
            !Ideal gas contribution

            !part_id = D2B_0_Dtau2 / zcLKP + D2C_0_Dtau2 * del / zcLKP2 + 6.d0 * lkp_c4_0 * del * tau  / zcLKP2 &
            !     & * (lkp_beta_0 + lkp_gamma_0 * del2 / zcLKP2) * dexp(-lkp_gamma_0 * del2 / zcLKP2)
            !!Reference fluid contribution
            !part_ref = D2B_ref_Dtau2 / zcLKP + D2C_ref_Dtau2 * del / zcLKP2 + 6.d0 * lkp_c4_ref * del * tau  / zcLKP2 &
            !     & * (lkp_beta_ref + lkp_gamma_ref * del2 / zcLKP2) * dexp(-lkp_gamma_ref * del2 / zcLKP2)
            !MIXDERFNR(7) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref)*del*tau2

        end if

        if (GETDER(8) == 1) then

            !Ideal fluid contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (2.d0 / del3 + 12.d0 * D_0 * del2 / zcLKP5 &
                & - 12.d0 * gl%lkp_c4_0 * tau3 * del3 * gl%lkp_gamma_0 ** 2 * help_id1 / zcLKP6 &
                & + 12.d0 * gl%lkp_c4_0 * tau3 * del * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                & + 4.d0 * gl%lkp_c4_0 * tau3 * del3 * gl%lkp_gamma_0 ** 3 * help_id1 * help_id2 / zcLKP8 &
                & - 6.d0 * gl%lkp_c4_0 * tau3 * del * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (2.d0 / del3 + 12.d0 * D_ref * del2 / zcLKP5 &
                & - 12.d0 * gl%lkp_c4_ref * tau3 * del3 * gl%lkp_gamma_ref ** 2 * help_ref1 / zcLKP6 &
                & + 12.d0 * gl%lkp_c4_ref * tau3 * del * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                & + 4.d0 * gl%lkp_c4_ref * tau3 * del3 * gl%lkp_gamma_ref ** 3 * help_ref1 * help_ref2 / zcLKP8 &
                & - 6.d0 * gl%lkp_c4_ref * tau3 * del * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
            MIXDERFNR(8) = (part_id + part_ref - 2.d0 / del3)*del3

            !!Modified Stefan & Monika 2.12.14
            !!Ideal fluid contribution
            !help_id = lkp_gamma_0 * del2 /zcLKP2
            !part_id = 12.d0 * D_0 * del2 / zcLKP5 + 2.d0 * lkp_c4_0 * tau3 * lkp_gamma_0 / zcLKP4 * exp(-help_id) * &
            !    & (-del * lkp_beta_0 + 3.d0 * del - lkp_beta_0 + help_id * (-3.d0 * del + 2.d0 * (lkp_beta_0 + help_id) - 3.d0))
            !!Reference fluid contribution
            !help_ref = lkp_gamma_ref * del2 /zcLKP2
            !part_ref = 12.d0 * D_ref * del2 / zcLKP5 + 2.d0 * lkp_c4_ref * tau3 * lkp_gamma_ref / zcLKP4 * exp(-help_ref) * &
            !    & (-del * lkp_beta_ref + 3.d0 * del - lkp_beta_ref + help_ref * (-3.d0 * del + 2.d0 * (lkp_beta_ref + help_ref) - 3.d0))
            !MIXDERFNR(8) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref)*del3

        end if


        if (GETDER(9) == 1) then

            D3B_0_Dtau3 = -6.d0 * gl%lkp_b4_0
            D3C_0_Dtau3 = 6.d0 * gl%lkp_c3_0
            D3B_ref_Dtau3 = -6.d0 * gl%lkp_b4_ref
            D3C_ref_Dtau3 = 6.d0 * gl%lkp_c3_ref

            !Ideal fluid contribution
            part_id = D3B_0_Dtau3 / gl%zcLKP * del + 0.5d0 * D3C_0_Dtau3 / zcLKP2 * del2 &
                & -3.d0 * gl%lkp_c4_0 / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_0 / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
            !Reference fluid contribution
            part_ref = D3B_ref_Dtau3 / gl%zcLKP * del + 0.5d0 * D3C_ref_Dtau3 / zcLKP2 * del2 &
                & -3.d0 * gl%lkp_c4_ref / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 3.d0 * gl%lkp_c4_ref / gl%lkp_gamma_ref* (gl%lkp_beta_ref + 1.d0)

            MIXDERFNR(9)=((1-gl%accenLKPMix / gl%lkp_w_ref)*part_id + gl%accenLKPMix / gl%lkp_w_ref*part_ref) * tau3

        end if


        if (GETDER(10) == 1) then

            DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
            DD_0_Dtau = gl%lkp_d2_0
            DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
            DD_ref_Dtau = gl%lkp_d2_ref

            !Ideal fluid contribution
            help_id1 = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            help_id2 = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
            part_id = (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * (DC_0_Dtau / zcLKP2 + 4.d0 * DD_0_Dtau * del3 / zcLKP5 &
                & - 3.d0 * gl%lkp_c4_0 * tau2 * help_id1 / zcLKP2 + 12.d0 * gl%lkp_c4_0 * tau2 * del2 * gl%lkp_gamma_0 * help_id1 / zcLKP4 &
                & + 3.d0 * gl%lkp_c4_0 * tau2 * gl%lkp_gamma_0 * help_id1 * help_id2 / zcLKP4 &
                & - 6.d0 * gl%lkp_c4_0 * tau2 * del2 * gl%lkp_gamma_0 ** 2 * help_id1 * help_id2 / zcLKP6)
            !Reference fluid contribution
            help_ref1 = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            help_ref2 = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
            part_ref = gl%accenLKPMix / gl%lkp_w_ref * (DC_ref_Dtau / zcLKP2 + 4.d0 * DD_ref_Dtau * del3 / zcLKP5 &
                & - 3.d0 * gl%lkp_c4_ref * tau2 * help_ref1 / zcLKP2 + 12.d0 * gl%lkp_c4_ref * tau2 * del2 * gl%lkp_gamma_ref * help_ref1 / zcLKP4 &
                & + 3.d0 * gl%lkp_c4_ref * tau2 * gl%lkp_gamma_ref * help_ref1 * help_ref2 / zcLKP4 &
                & - 6.d0 * gl%lkp_c4_ref * tau2 * del2 * gl%lkp_gamma_ref ** 2 * help_ref1 * help_ref2 / zcLKP6)
            MIXDERFNR(10) = (part_id + part_ref)*tau*del2

            !!Modified Monika & Stefan 2.12.14
            !!Ideal fluid contribution
            !help_id = lkp_gamma_0 * del2 / zcLKP2
            !part_id = DC_0_Dtau / zcLKP2 + 4.d0 * DD_0_Dtau * del3 / zcLKP5 &
            !     & + 3.d0 * lkp_c4_0 * tau2 / zcLKP2 * dexp(-help_id)*(lkp_beta_0 + 3.d0 * help_id - 2.d0 * help_id / del * (lkp_beta_0 + help_id))
            !!Reference fluid contribution
            !help_ref = lkp_gamma_ref * del2 / zcLKP2
            !part_ref = DC_ref_Dtau / zcLKP2 + 4.d0 * DD_ref_Dtau * del3 / zcLKP5 &
            !     & + 3.d0 * lkp_c4_ref * tau2 / zcLKP2 * dexp(-help_ref)*(lkp_beta_ref + 3.d0 * help_ref - 2.d0 * help_ref / del * (lkp_beta_ref + help_ref))
            !MIXDERFNR(10) = ((1-accenLKPMix / lkp_w_ref)*part_id + accenLKPMix / lkp_w_ref*part_ref)*tau*del2

        end if

    elseif (gl%Mix_type .eq. 6) then
        nrsubst = 0
        call ARDERIVS(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR, nrsubst)


        !AGA 8 -  thermodynamic properties of natural gas and related gas
    elseif (gl%mix_type == 7) then  !aga8

        call mixture_parameters(gl,Temperature,Density, B_0)
        call residualpart_AGA8(Temperature, Density, GETDER, MIXDERFNR)


    elseif (gl%Mix_type == 11) then

        !Andreas, February 2016. New quadratic mixing rules for the residual Helmholtz energy.
        !See ...

        !sum of the residual Helmholtz energies and its derivatives of the pure fluids
        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        FNRDER = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,TEMPERATURE, DENSITY, GETDERS, FNRDER, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_r_i(i,:) = FNRDER
        end do

        do i = 1, gl%ncomp
            do j = i, gl%ncomp
                if (i == j) then
                    alpha_r_ij(:) = alpha_r_i(i,:)
                end if
                if (i /= j) then
                    alpha_r_ij(:) = 0.5D0 * (alpha_r_i(i,:) + alpha_r_i(j,:)) * (1.D0 - ACCESS_KIJ_HELM(i,j)) * 2.D0    !Because of symmetry * 2.D0
                end if
                do k = 1, 10
                    if (GETDERS(k) == 1) then     ! check which derivative is needed
                        MIXDERFNR(k) = MIXDERFNR(k) + gl%MOLFRACTIONS(i)*gl%MOLFRACTIONS(j)*alpha_r_ij(k) ! sum of derivatives over all pure fluids, multiplied with the respective mole fraction
                    end if
                end do
            end do
        end do

        !sum of the binary departure functions
        do FLD1 = 1, gl%NCOMP - 1               !
            do FLD2 = FLD1 + 1, gl%NCOMP
                if (dabs(gl%Fij(FLD1, FLD2)) > 1.d-14) then    ! check if departure function exists for this binary combination
                    call DEPFUNCFNR (gl,TEMPERATURE, DENSITY, GETDERS, DEPFUNCBIN, FLD1, FLD2)   ! calls the routine for the calculation of the binary specific departure function and its derivatives
                    xi = gl%MOLFRACTIONS(FLD1)
                    xj = gl%MOLFRACTIONS(FLD2)
                    do k = 1, 10                      ! loop through all derivatives returned from DEPFUNCFNR
                        DEPFUNCDER(k) = DEPFUNCDER(k) + xi*xj*gl%Fij(FLD1, FLD2)*DEPFUNCBIN(k)     ! sum over all binary departure functions, multiplied with the mole fractions and the weighing factor Fij
                    end do
                    DEPFUNCBIN = 0.d0
                end if
            end do
        end do

        ! add the two parts of the residual Helmholtz enthalpy for mixtures
        do m = 1, nderivs
            MIXDERFNR(m) = MIXDERFNR(m) + DEPFUNCDER(m)
        end do


    elseif ((gl%Mix_type == 12) .or. (gl%Mix_type == 13)) then     !Excess based departure function

        !sum of the residual Helmholtz energies and its derivatives of the pure fluids
        do i = 1, gl%NCOMP
            FNRDER = 0.D0
            call FNRDERIVS(gl,TEMPERATURE, DENSITY, GETDERS, FNRDER, i) ! call the calculation routine for the derivatives of the fluid i
            do j = 1, nderivs
                if (GETDERS(j) == 1) then     ! check which derivative is needed
                    MIXDERFNR(j) = MIXDERFNR(j) + gl%MOLFRACTIONS(i)*FNRDER(j) ! sum of derivatives over all pure fluids, multiplied with the respective mole fraction
                end if
            end do
        end do

        !excess based departure function and derivatives
        call DEPFUNC_GE_BASED (gl,TEMPERATURE, DENSITY, GETDERS, DEPFUNCDER_GE)

        ! add the two parts of the residual Helmholtz energy for mixtures
        do m = 1, nderivs
            MIXDERFNR(m) = MIXDERFNR(m) + DEPFUNCDER_GE(m)
        end do


    elseif (gl%mix_type == 19) then     !non-corresponding states based mixture model (one-fluid model)

        !sum of the residual Helmholtz energies and its derivatives of the pure fluids
        MIXDERFNR = 0.D0
        do i = 1, gl%NCOMP
            FNRDER = 0.D0
            !Evaluate the pure fluids at the mixture temperature and density and not like in multifluid models at the same tau and delta
            gl%tredmix = gl%tc(i)
            gl%rhoredmix = gl%rhoc(i)
            call FNRDERIVS(gl,TEMPERATURE, DENSITY, GETDERS, FNRDER, i) ! call the calculation routine for the derivatives of the fluid i
            gl%tredmix = 1.D0
            gl%rhoredmix = 1.D0
            !The following summation still works and does not need to be adjusted! In the one fluid mixture model, the part that is calculated here is:
            ! alphar(T,rho,x) = x1^2 alphar1(T,rho) + x2^2 alphar2(T,rho) + ...
            ! Keeping the reduced variables delta and tau but now defining for a mixture: tau = 1 / T and delta = rho, for example the following derivatives are needed:
            ! tau * dalphar / dtau = x1^2 tau * dalphar1/dtau + tau * x2^2 dalphar2/dtau + ... = x1^2 tau * dalphar1/dtau1 * dtau1/dtau + tau * x2^2 dalphar2/dtau2 * dtau2/dtau + ...
            ! It is:  dtaui/dtau = d(tci/temp)/d(1/temp) = tci
            ! And consequently:
            ! tau * dalphar / dtau = x1^2 tau * dalphar1/dtau1 * tc1 + tau * x2^2 dalphar2/dtau2 * tc2 + ... = x1^2 tau1 * dalphar1/dtau1 + tau2 * x2^2 dalphar2/dtau2 + ...
            !
            ! It can be easily shown that this is also true for the delta-derivatives and higher order derivatives!
            !
            do j = 1, nderivs
                if (GETDERS(j) == 1) then     ! check which derivative is needed
                    MIXDERFNR(j) = MIXDERFNR(j) + gl%MOLFRACTIONS(i)**2 * FNRDER(j) ! sum of derivatives over all pure fluids, multiplied with the respective mole fraction to the power of 2
                end if
            end do
        end do


        !The following function can be regarded as the "departure function" of the new one-fluid model
        tau = 1.D0 / Temperature
        delta = Density
        DEPFUNCDER = 0.D0
        !sum over the binary functions
        do i = 1, (gl%ncomp-1)
            !deltai = density / rhoc(i)
            !taui = tc(i) / temperature

            do j = i+1, gl%ncomp
                !deltaj = density / rhoc(j)
                !tauj = tc(j) / temperature
                !
                !!!Summation over the regular terms
                !do k = 1, mix_nreg(i,j)
                !    reg_termi(k) = mix_nij(i,j,k,1) * taui**mix_tij(i,j,k) * deltai**mix_dij(i,j,k) * &
                !                &  dexp(-mix_gama(i,j,k) * deltai**mix_p_ij(i,j,k))  * (1.D0 - Helm_k_nij(i,j,k)) * 0.5D0
                !    reg_termj(k) = mix_nij(i,j,k,2) * tauj**mix_tij(i,j,k) * deltaj**mix_dij(i,j,k) * &
                !                &  dexp(-mix_gama(i,j,k) * deltaj**mix_p_ij(i,j,k))  * (1.D0 - Helm_k_nij(i,j,k)) * 0.5D0
                !    DEPFUNCDER(1) = DEPFUNCDER(1) + (reg_termi(k) + reg_termj(k)) * 2.D0 * molfractions(i) * molfractions(j)
                !end do
                !
                !!The departure function of the one-fluid model itself
                !if (GETDERS(1) .eq. 1) then
                !  DEPFUNCDER(1) = DEPFUNCDER(1)
                !end if
                !
                !!The derivative of the departure function with respect to delta, multiplied with delta
                !if (GETDERS(2) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(2)= DEPFUNCDER(2) &
                !                   & + (reg_termi(k) * (mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltai**mix_p_ij(i,j,k)) &
                !                   & + reg_termj(k) * (mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltaj**mix_p_ij(i,j,k)))  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The second derivative of the departure function with respect to delta^2, multiplied with delta^2
                !if (GETDERS(3) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(3)= DEPFUNCDER(3) &
                !                   & + (reg_termi(k) * ((mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltai**mix_p_ij(i,j,k)) * &
                !                   &   (mix_dij(i,j,k) - 1.D0 - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltai**mix_p_ij(i,j,k)) - &
                !                   &    mix_gama(i,j,k) * (mix_p_ij(i,j,k)**2) * deltai**mix_p_ij(i,j,k)  ) &
                !                   & + reg_termj(k) * ((mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltaj**mix_p_ij(i,j,k)) * &
                !                   &   (mix_dij(i,j,k) - 1.D0 - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltaj**mix_p_ij(i,j,k)) - &
                !                   &    mix_gama(i,j,k) * (mix_p_ij(i,j,k)**2) * deltaj**mix_p_ij(i,j,k)  ) )  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau, multiplied with tau
                !if (GETDERS(4) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(4)= DEPFUNCDER(4) + ((reg_termi(k) + reg_termj(k)) * mix_tij(i,j,k))  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau^2, multiplied with tau^2
                !if (GETDERS(5) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(5)= DEPFUNCDER(5) + ((reg_termi(k) + reg_termj(k)) * mix_tij(i,j,k) * (mix_tij(i,j,k) - 1.D0))  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau and del, multiplied with tau*del
                !if (GETDERS(6) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(6)= DEPFUNCDER(6) + &
                !            & ((reg_termi(k) * (mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltai**mix_p_ij(i,j,k)) + &
                !            & reg_termj(k) * (mix_dij(i,j,k) - mix_gama(i,j,k) * mix_p_ij(i,j,k) * deltaj**mix_p_ij(i,j,k))) * mix_tij(i,j,k))  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau^2 and del, multiplied with tau^2*del
                !if (GETDERS(7) .eq. 1) then
                !  DEPFUNCDER(7) = 0.D0
                !end if
                !
                !!The derivative of the departure function with respect to del^3, multiplied with del^3
                !if (GETDERS(8) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(8)= DEPFUNCDER(8) &
                !                & + (reg_termi(k) * (mix_dij(i,j,k)*(mix_dij(i,j,k) - 1.d0)*(mix_dij(i,j,k) - 2.d0) + &
                !                &   mix_gama(i,j,k)*mix_p_ij(i,j,k)*deltai**mix_p_ij(i,j,k)*(-2.d0 + 6.d0*mix_dij(i,j,k) &
                !                &   - 3.d0*mix_dij(i,j,k)**2 - 3.d0*mix_dij(i,j,k)*mix_p_ij(i,j,k) + 3.d0*mix_p_ij(i,j,k) - mix_p_ij(i,j,k)**2) &
                !                &   + 3.d0*mix_gama(i,j,k)**2*mix_p_ij(i,j,k)**2*deltai**(2.d0*mix_p_ij(i,j,k))*(mix_dij(i,j,k) - 1.d0 + mix_p_ij(i,j,k))&
                !                &   - mix_gama(i,j,k)**3*mix_p_ij(i,j,k)**3*deltai**(3*mix_p_ij(i,j,k)) ) &
                !                & + reg_termj(k) * (mix_dij(i,j,k)*(mix_dij(i,j,k) - 1.d0)*(mix_dij(i,j,k) - 2.d0) + &
                !                &   mix_gama(i,j,k)*mix_p_ij(i,j,k)*deltaj**mix_p_ij(i,j,k)*(-2.d0 + 6.d0*mix_dij(i,j,k) &
                !                &   - 3.d0*mix_dij(i,j,k)**2 - 3.d0*mix_dij(i,j,k)*mix_p_ij(i,j,k) + 3.d0*mix_p_ij(i,j,k) - mix_p_ij(i,j,k)**2) &
                !                &   + 3.d0*mix_gama(i,j,k)**2*mix_p_ij(i,j,k)**2*deltaj**(2.d0*mix_p_ij(i,j,k))*(mix_dij(i,j,k) - 1.d0 + mix_p_ij(i,j,k))&
                !                &   - mix_gama(i,j,k)**3*mix_p_ij(i,j,k)**3*deltaj**(3*mix_p_ij(i,j,k)) ) )  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau^3, multiplied with tau^3
                !if (GETDERS(9) .eq. 1) then
                !    !!Summation over the regular terms
                !    do k = 1, mix_nreg(i,j)
                !        DEPFUNCDER(9)= DEPFUNCDER(9) + ((reg_termi(k) + reg_termj(k)) * &
                !                     & mix_tij(i,j,k) * (mix_tij(i,j,k) - 1.D0) * (mix_tij(i,j,k) - 2.D0))  * 2.D0 * molfractions(i) * molfractions(j)
                !    end do
                !end if
                !
                !!The derivative of the departure function with respect to tau*del*del, multiplied with tau*del*del
                !if (GETDERS(10) .eq. 1) then
                !  DEPFUNCDER(10) = 0.D0
                !end if

                call DEPFUNC_NON_COR_STATE (gl,Temperature, Density, GETDERS, DEPFUNCBIN,i,j)
                do k = 1, 10                      ! loop through all derivatives returned from DEPFUNC_NON_COR_STATE
                    DEPFUNCDER(k) = DEPFUNCDER(k) + 2.D0*gl%molfractions(i)*gl%molfractions(j)*DEPFUNCBIN(k)     ! sum over all binary departure functions, multiplied with the mole fractions and the weighing factor Fij
                end do

            end do
        end do

        ! add the two parts of the residual Helmholtz energy for mixtures
        do m = 1, nderivs
            MIXDERFNR(m) = MIXDERFNR(m) + DEPFUNCDER(m)
        end do

    end if


    !! saving the derivatives and belonging inputs
    !gl%TEMP_FNR_MIX_OLD = TEMPERATURE
    !gl%DENS_FNR_MIX_OLD = DENSITY
    !gl%tredmix_MIX_old = gl%tredmix
    !gl%rhoredmix_MIX_old = gl%rhoredmix
    !comp_FNR_MIX_old = components
    !mol_FNR_MIX_old = molfractions
    !gl%FNRDER_MIX_OLD = MIXDERFNR


    end subroutine MIXDERIVSFNR




    !**************************************************************************
    !           --------------------------------------------------
    !           Routine for the calculation of the departure
    !           function (the residual mixture related part of
    !           the Helmholtz free energy function) and its
    !           derivatives
    !
    !           J. Gernert, Denmark, 09.2009
    !           --------------------------------------------------
    !**************************************************************************

    !**************************************************************************
    subroutine DEPFUNCFNR (gl,TEMPERATURE, DENSITY, GETDER, DEPFUNCDER, FLD1_arg, FLD2_arg)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE DEPARTURE FUNCTION DF, WHICH DESCRIBES THE
    ! RESIDUAL MIXTURE RELATED PART OF THE HELMHOLTZ FREE ENERGY
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    ! GETDER      - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED DEPARTURE FUNCTION DF AS A FUNCTION OF DEL AND TAU
    !                2. 1ST DERIVATIVE OF DF WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF DF WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
    !                4. 1ST DERIVATIVE OF DF WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF DF WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF DF WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !                7: 2ND MIXED DERIVATIVE OF DF WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF DF WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    ! FL1, FLD2   - BINARY SYSTEM, FOR WHICH THE BINARY DEPARTURE FUNCTION IS CALCULATED
    !
    ! OUTPUT PARAMETERS:
    ! DEPFUNCDER  - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !                AS INDECATED IN "GETDER"
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: Temperature, Density
    integer, dimension(nderivs):: GETDER            ! array specifier to indicate, which derivative is needed
    integer :: FLD1, FLD2, FLD1_arg, FLD2_arg                      ! line index of the array DFCOEFF which contains the coefficients of the binary departure functions
    double precision, dimension(nderivs)::DEPFUNCDER ! array with the computed values for the derivatives
    double precision:: del, tau                             ! reduced temperature and density
    double precision:: del2,del3,del4,del5,del6
    double precision, dimension (25):: pol_term             ! array of the polynomial terms
    double precision, dimension (25):: exp_term             ! array of the exponential terms
    double precision, dimension (25):: gauss_term           ! array of the gaussian terms
    double precision, dimension (25):: ni_depf, di_depf, ti_depf, li, pi      ! coefficients, density and temperature exponents, exponential term exponents
    double precision, dimension (25):: eta_depf, eps_depf, beta_depf, gama_depf                ! array of the parameters in the exponential terms
    double precision, dimension (25):: geta, geps, gbeta, ggama            ! array of the parameters in the gaussian terms
    double precision, dimension (25):: dli                  ! term for internal calc.
    double precision:: DEPF, DEPFD, DEPFDD, DEPFDDD                 ! departure function DF and its derivatives
    double precision:: DEPFT, DEPFTT, DEPFTD, DEPFDTT               ! with respect to density and temperature
    double precision:: DEPFTTT, DEPFDDT
    integer:: I_pol_mix, I_exp_mix, I_gau_mix             ! number of polynomial, exponential, gaussian bell-shaped, and nonanalytic terms
    integer:: i, j, m, o                                  ! index variables
    integer :: cache_index
    integer :: di_int, li_int, pi_int

    ! ensure FLD1 < FLD2
    if (FLD1_arg < FLD2_arg) then
        FLD1 = FLD1_arg
        FLD2 = FLD2_arg
    else
        FLD1 = FLD2_arg
        FLD2 = FLD1_arg
    endif

    !set all array elements to default value 0
    pol_term = 0.D0
    exp_term = 0.D0
    gauss_term = 0.d0
    ni_depf = 0.D0
    di_depf = 0.D0
    ti_depf = 0.D0
    li = 0.D0
    pi = 0.D0
    dli = 0.D0
    eta_depf = 0.D0
    eps_depf = 0.D0
    beta_depf = 0.D0
    gama_depf = 0.D0
    geta = 0.d0
    gbeta = 0.d0
    ggama = 0.d0
    geps = 0.d0
    I_pol_mix = 0
    I_exp_mix = 0
    I_gau_mix = 0
    li_int = 0
    pi_int = 0

    DEPFUNCDER = 0.D0

    ! computes del and tau from reducing parameters
    !Sebastian, September 2016
    !For the calculation of virial coefficients in mixtures density must no be 0
    !the evaluation of the virial coefficients (see the following else-statement)
    if (Density <= 1.d-12) then
        del = 1.d-12 / gl%rhoredmix
    else
        del = Density / gl%rhoredmix      ! reduced density calculated with reducing function for mixtures
    endif
    tau = gl%tredmix / Temperature    ! inverse reduced temperature calculated with reducing function for mixtures


    del2 = del*del
    del3 = del2*del
    del4 = del3*del
    del5 = del4*del
    del6 = del5*del

    I_pol_mix = int(gl%DFPOL(FLD1, FLD2))      ! number of polynomial terms
    I_exp_mix = int(gl%DFEXP(FLD1, FLD2))    ! number of exponential terms
    I_gau_mix = int(gl%DFGAU(FLD1, FLD2))    ! number of gaussian terms

    do m = 1, I_pol_mix + I_exp_mix + I_gau_mix
        ni_depf(m) = gl%DFN( m,FLD1, FLD2)          ! reading the coefficients of the polynomial terms from the module array into the local vector
        di_depf(m) = gl%DFD( m,FLD1, FLD2)          ! reading the density exponents of the polynomial terms from the module array into the local vector
        ti_depf(m) = gl%DFT( m,FLD1, FLD2)          ! reading the temperature exponents of the polynomial terms from the module array into the local vector
        li(m) = gl%DFL( m,FLD1, FLD2)          ! reading the exponential term density exponents from the module array into the local vector
        pi(m) = gl%DFP( m,FLD1, FLD2)          ! reading the exponential term density exponents from the module array into the local vector
    end do
    do o = 1, I_exp_mix
        eta_depf(o) = gl%DFETA( o,FLD1, FLD2)       ! reading the eta of the exponential terms from the module array into the local vector
        eps_depf(o) = gl%DFEPS( o,FLD1, FLD2)       ! reading the eps of the exponential terms from the module array into the local vector
        beta_depf(o) = gl%DFBETA( o,FLD1, FLD2)     ! reading the beta of the exponential terms from the module array into the local vector
        gama_depf(o) = gl%DFGAMMA( o,FLD1, FLD2)    ! reading the gamma of the exponential terms from the module array into the local vector
    end do

    do o = 1, I_gau_mix
        geta(o) = gl%DFGETA( o,FLD1, FLD2)       ! reading the beta of the gaussian terms from the module array into the local vector
        gbeta(o) = gl%DFGBETA( o,FLD1, FLD2)     ! reading the temperature exponents of the gaussian terms from the module array into the local vector
        ggama(o) = gl%DFGGAM( o,FLD1, FLD2)      ! reading the temperature exponents of the gaussian terms from the module array into the local vector
        geps(o) = gl%DFGEPS( o,FLD1, FLD2)       ! reading the density exponents of the gaussian terms from the module array into the local vector
    end do


    ! check cache
    if (enable_cache .and. FLD1 <= 5 .and. FLD2 <= 5) then
        cache_index = upper_triangle_offset(FLD2) + FLD1
    else
        cache_index = 0
    endif


    if (cache_index /= 0 .and. del == gl%depfuncfnr_cache%del_old(cache_index)) then
        continue ! nothing to do
    else
        gl%depfuncfnr_cache%del_old(cache_index) = del
        m = I_pol_mix + I_exp_mix + I_gau_mix

        select case(gl%dfd_coeff_structure(FLD1, FLD2))

        case (coeff_structure_all_zeros)
            gl%depfuncfnr_cache%del_pow_di_depf(1:m,cache_index) = 1d0

        case (coeff_structure_int)
            do j=1,m
                di_int = int(di_depf(j))
                select case(di_int)
                case (0)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = 1d0

                case (1)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del

                case (2)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**2

                case (3)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**3

                case (4)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**4

                case (5)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**5

                case (6)
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**6

                    case default
                    gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) = del**di_int
                end select
            end do
            !gl%depfuncfnr_cache%del_pow_di_depf(1:m,cache_index) = del**int(di_depf(1:m))
            case default
            gl%depfuncfnr_cache%del_pow_di_depf(1:m,cache_index) = del**di_depf(1:m)
        end select
    endif

    if (cache_index /= 0 .and. tau == gl%depfuncfnr_cache%tau_old(cache_index)) then
        continue ! nothing to do
    else
        gl%depfuncfnr_cache%tau_old(cache_index) = tau

        m = I_pol_mix + I_exp_mix + I_gau_mix
        gl%depfuncfnr_cache%tau_pow_ti_depf(1:m,cache_index) = tau**ti_depf(1:m)
    endif


    !*****************************************************
    !   DEPF         DEPARTURE FUNCTION DF [-]
    !*****************************************************
    DEPF = 0.D0

    !summation over the polynomial terms
    select case(gl%dfl_coeff_structure(FLD1, FLD2))

    case (coeff_structure_all_zeros)
        do i = 1, I_pol_mix
            dli(i) = 1d0
            !Andreas, December 2015
            !Changed comparison here, checking if a double value is equal to something can cause problems
            !if (li(i) == 0.D0) then     ! This checks if the term has an exponential term with the form exp(-delta^l)
            pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index)
            DEPF = DEPF + pol_term(i)
        end do

    case (coeff_structure_int)
        do i = 1, I_pol_mix
            li_int = int(li(i))
            pi_int = int(pi(i))
            if (li_int*pi_int == 0) li_int = 0
            !Andreas, December 2015
            !Changed comparison here, checking if a double value is equal to something can cause problems
            !if (li(i) == 0.D0) then     ! This checks if the term has an exponential term with the form exp(-delta^l)
            select case(li_int)
            case (0)
                dli(i) = 0.d0
                !2019/10/22 changed by Moni and Sebastian
                !dli(i) = 1.d0
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index)

            case (1)
                dli(i) = del
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))

            case (2)
                dli(i) = del*del
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))

            case (3)
                dli(i) = del*del*del
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))

                case default
                dli(i) = del**li_int
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))
            end select
            DEPF = DEPF + pol_term(i)
        end do

        case default
        do i = 1, I_pol_mix
            dli(i) = del**li(i)
            !Andreas, December 2015
            !Changed comparison here, checking if a double value is equal to something can cause problems
            !if (li(i) == 0.D0) then     ! This checks if the term has an exponential term with the form exp(-delta^l)
            if(dabs(li(i)) < 1.D-8) then
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index)
            else
                pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))
            end if
            DEPF = DEPF + pol_term(i)
        end do
    end select

    !do i = 1, I_pol_mix
    !    dli(i) = del**li(i)
    !    !Andreas, December 2015
    !    !Changed comparison here, checking if a double value is equal to something can cause problems
    !    !if (li(i) == 0.D0) then     ! This checks if the term has an exponential term with the form exp(-delta^l)
    !    if(dabs(li(i)) < 1.D-8) then
    !        pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index)
    !    else
    !        pol_term(i) = ni_depf(i) * gl%depfuncfnr_cache%del_pow_di_depf(i,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(i,cache_index) * dexp(-dli(i))
    !    end if
    !    DEPF = DEPF + pol_term(i)
    !end do

    !summation over the exponential terms
    do i = 1, I_exp_mix
        j = I_pol_mix + i
        exp_term(i) = ni_depf(j) * gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(j,cache_index) * dexp(-eta_depf(i)*(del - eps_depf(i))**2 - beta_depf(i)*(del - gama_depf(i)))
        DEPF = DEPF + exp_term(i)
    end do

    !summation over the gaussian terms
    do i = 1, I_gau_mix
        j = I_pol_mix + I_exp_mix + i
        gauss_term(i) = ni_depf(j) * gl%depfuncfnr_cache%del_pow_di_depf(j,cache_index) * gl%depfuncfnr_cache%tau_pow_ti_depf(j,cache_index) * dexp(geta(i)*(del - geps(i))**2 + gbeta(i)*(tau - ggama(i))**2)
        DEPF = DEPF + gauss_term(i)
    end do


    !*****************************************************
    if (GETDER(1) == 1) then
        DEPFUNCDER(1) = DEPF     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(2) == 1) then
        !*****************************************************
        !    DEPFD       1ST DERIVATIVE OF DF WITH RESPECT TO
        !               del MULTIPLIED WITH DEL [-]
        !*****************************************************
        DEPFD = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFD = DEPFD + pol_term(i) *( di_depf(i) - dli(i) * li(i))
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFD = DEPFD + exp_term(i) * (di_depf(j) + del*(-2.D0*eta_depf(i)*(del - eps_depf(i)) - beta_depf(i)))
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFD = DEPFD + gauss_term(i) * (di_depf(j) + 2.d0*del*geta(i)*(del - geps(i)))
        end do

        DEPFUNCDER(2) = DEPFD     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(3) == 1) then
        !*****************************************************
        !   DEPFDD       2ND DERIVATIVE OF DF WITH RESPECT TO
        !               del MULTIPLIED WITH DEL^2 [-]
        !*****************************************************
        DEPFDD = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFDD = DEPFDD + pol_term(i) * &
                & ((di_depf(i) - dli(i)*li(i))**2 - di_depf(i) - dli(i)*li(i)*(li(i) - 1.D0))
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFDD = DEPFDD + exp_term(i) * ((di_depf(j) + del*(-2.D0*eta_depf(i)*(del - eps_depf(i)) - beta_depf(i)))**2 &
                - di_depf(j) - 2.D0*eta_depf(i)*del2)
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFDD = DEPFDD + gauss_term(i) * ((di_depf(j) + 2.d0*del*geta(i)*(del - geps(i)))**2 - di_depf(j) + 2.d0*geta(i)*del**2.d0)
        end do

        DEPFUNCDER(3) = DEPFDD     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(4) == 1) then
        !*****************************************************
        !  DEPFT         1ST DERIVATIVE OF DF WITH RESPECT TO
        !               TAU MULTIPLIED WITH TAU [-]
        !*****************************************************
        DEPFT = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFT = DEPFT + pol_term(i) * ti_depf(i)
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFT = DEPFT + exp_term(i) * ti_depf(j)
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFT = DEPFT + gauss_term(i) * (ti_depf(j) + 2.d0*tau*gbeta(i)*(tau-ggama(i)))
        end do

        DEPFUNCDER(4) = DEPFT     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(5) == 1) then
        !*****************************************************
        !  DEPFTT        2ND DERIVATIVE OF DF WITH RESPECT TO
        !               TAU MULTIPLIED WITH TAU^2 [-]
        !*****************************************************
        DEPFTT = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFTT = DEPFTT + pol_term(i) * ti_depf(i)*(ti_depf(i) - 1.D0)
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFTT = DEPFTT + exp_term(i) * ti_depf(j)*(ti_depf(j) - 1.D0)
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFTT = DEPFTT + gauss_term(i) * ((ti_depf(j) + 2.d0*tau*gbeta(i)*(tau-ggama(i)))**2 - ti_depf(j) + 2.d0*gbeta(i)*tau**2.d0)
        end do

        DEPFUNCDER(5) = DEPFTT     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(6) == 1) then
        !*****************************************************
        !   DEPFTD       1ST MIXED DERIVATIVE OF DF WITH RESPECT
        !               TO TAU AND DEL MULTIPLIED WITH
        !               TAU*DEL [-]
        !*****************************************************
        DEPFTD = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFTD = DEPFTD + pol_term(i) * ti_depf(i) * (di_depf(i) - dli(i) * li(i))
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFTD = DEPFTD + exp_term(i) * ti_depf(j)*(di_depf(j) + del*(-2.D0*eta_depf(i)*(del - eps_depf(i)) - beta_depf(i)))
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFTD = DEPFTD + gauss_term(i) * (di_depf(j) + 2.d0*del*geta(i)*(del-geps(i))) * (ti_depf(j) + 2.d0*tau*gbeta(i)*(tau-ggama(i)))
        end do

        DEPFUNCDER(6) = DEPFTD     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(7) == 1) then
        !*****************************************************
        !     DEPFDTT    2ND MIXED DERIVATIVE OF DF WITH RESPECT
        !               TO TAU AND DEL MULTIPLIED WITH
        !               DEL*TAU^2 [-]
        !*****************************************************
        DEPFDTT = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFDTT = DEPFDTT + pol_term(i) * ti_depf(i)*(ti_depf(i) - 1.D0) * (di_depf(i) - dli(i) * li(i))
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFDTT = DEPFDTT + exp_term(i) * ti_depf(j)*(ti_depf(j) - 1.D0)*(di_depf(j) + del*(-2.D0*eta_depf(i) &
                *(del - eps_depf(i)) - beta_depf(i)))
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFDTT = DEPFDTT + gauss_term(i) * (di_depf(j) + 2.d0*del*geta(i)*(del-geps(i))) &
                * ((ti_depf(j) + 2.d0*tau*gbeta(i)*(tau-ggama(i)))**2 - ti_depf(j) + 2.d0*gbeta(i)*tau**2.d0)
        end do

        DEPFUNCDER(7) = DEPFDTT     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(8) == 1) then
        !*****************************************************
        !     DEPFDDD    3RD DERIVATIVE OF DF WITH RESPECT
        !               TO DEL MULTIPLIED WITH del^3 [-]
        !*****************************************************
        DEPFDDD = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            !        DEPFDDD = DEPFDDD + pol_term(i) * ((di(i) - dli(i)*li(i))**3 - di(i)*(3.D0*di(i) - 2.D0) &
            !                & - dli(i)*li(i)*(3.D0*di(i)*(li(i) - 2.D0) + li(i)*(li(i) - 3.D0) - 2.D0 &
            !                & + 3.D0*dli(i)*li(i)*(1.D0 - li(i))))
            DEPFDDD = DEPFDDD + pol_term(i) * (di_depf(i) ** 3 - 3.d0 * di_depf(i) ** 2 - 3.d0 * del ** li(i) * di_depf(i) ** 2 * li(i) + &
                & 2.d0 * di_depf(i) + 6.d0 * del ** li(i) * di_depf(i) * li(i) - 3.d0 * del ** li(i) * di_depf(i) * li(i) ** 2 + &
                & 3.d0 * del ** (2.d0 * li(i)) * di_depf(i) * li(i) ** 2 - del ** li(i) * li(i) ** 3 + &
                & 3.d0 * del ** li(i) * li(i) ** 2 + 3.d0 * del ** (2.d0 * li(i)) * li(i) ** 3 - &
                & 2.d0 * del ** li(i) * li(i) - 3.d0 * del ** (2.d0 * li(i)) * li(i) ** 2 - &
                & del ** (3.d0 * li(i)) * li(i) ** 3)
            !        DEPFDDD = DEPFDDD + pol_term(i) * (di(j) ** 3 - 3.d0 * di(j) ** 2 - 3.d0 * del ** li(i) * di(j) ** 2 * li(i) + &
            !                 & 2.d0 * di(j) + 6.d0 * del ** li(i) * di(j) * li(i) - 3.d0 * del ** li(i) * di(j) * li(i) ** 2 + &
            !                 & 3.d0 * del ** (2 * li(i)) * di(j) * li(i) ** 2 - del ** li(i) * li(i) ** 3 + &
            !                 & 3.d0 * del ** li(i) * li(i) ** 2 + 3.d0 * del ** (2 * li(i)) * li(i) ** 3 - &
            !                 & 2.d0 * del ** li(i) * li(i) - 3.d0 * del ** (2 * li(i)) * li(i) ** 2 - &
            !                 & del ** (3 * li(i)) * li(i) ** 3)
        end do

        !summation over the exponential terms   !TODO-EXP
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            !        DEPFDDD = DEPFDDD + exp_term(i) * ((di(j) + del*(2.D0*eta(i)*(del - eps(i)) - beta(i)))**3 &
            !                - 3.D0*(di(j)*del**2.d0 + del**3.d0*(-2.D0*eta(i)*(del - eps(i)) - beta(i)))*(di(j)/(del**2.d0) &
            !                + 2.D0*eta(i)) + 2.D0*di(j))
            DEPFDDD = DEPFDDD + exp_term(i) * ( 2.d0 * di_depf(j) - 6.d0 * di_depf(j) * eta_depf(i) * eps_depf(i) * del - &
                & 12.d0 * di_depf(j) * eta_depf(i) * eps_depf(i) * beta_depf(i) * del2 + 3.d0 * di_depf(j) * beta_depf(i) * del - &
                & 12.d0 * eta_depf(i) ** 2 * eps_depf(i) * del3 + 6.d0 * eta_depf(i) * beta_depf(i) * del3 - 3.d0 * di_depf(j) ** 2 + &
                & 12.d0 * eta_depf(i) ** 2 * del4 - beta_depf(i) ** 3 * del3 - 8.d0 * eta_depf(i) ** 3 * del6 + di_depf(j) ** 3 + &
                & 8.d0 * eta_depf(i) ** 3 * eps_depf(i) ** 3 * del3 + 12.d0 * del4 * di_depf(j) * eta_depf(i) ** 2 + &
                & 24.d0 * eta_depf(i) ** 3 * del5 * eps_depf(i) - 12.d0 * eta_depf(i) ** 2 * del5 * beta_depf(i) - &
                & 24.d0 * eta_depf(i) ** 3 * del4 * eps_depf(i) ** 2 - 6.d0 * eta_depf(i) * del4 * beta_depf(i) ** 2 - &
                & 6.d0 * di_depf(j) ** 2 * eta_depf(i) * del2 - 3.d0 * di_depf(j) **2 * beta_depf(i) * del + &
                & 3.d0 * di_depf(j) * beta_depf(i) ** 2 * del2 + &
                & 6.d0 * di_depf(j) ** 2 * eta_depf(i) * eps_depf(i) * del + 12.d0 * di_depf(j) * eta_depf(i) ** 2 * eps_depf(i) ** 2 * del2 - &
                & 24.d0 * di_depf(j) * eta_depf(i) ** 2 * eps_depf(i) * del3 + 12.d0 * di_depf(j) * eta_depf(i) * beta_depf(i) * del3 - &
                & 12.d0 * eta_depf(i) ** 2 * eps_depf(i) ** 2 * beta_depf(i) * del3 + 6.d0 * eta_depf(i) * eps_depf(i) * beta_depf(i) ** 2 * del3 + &
                & 24.d0 * eta_depf(i) ** 2 * del4 * eps_depf(i) * beta_depf(i))
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFDDD = DEPFDDD + gauss_term(i) * ((di_depf(j) + 2.d0*del*geta(i)*(del-geps(i)))**3 - 3.d0*di_depf(j)**2 &
                + 2.d0*di_depf(j) + 6.d0*di_depf(j)*geta(i)*del**2.d0 - 6.d0*geta(i)*del*(del-geps(i))*(di_depf(j) - 2.d0*geta(i)*del**2.d0))
        end do

        DEPFUNCDER(8) = DEPFDDD     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(9) == 1) then
        !*****************************************************
        !  DEPFTTT      3ND DERIVATIVE OF DF WITH RESPECT TO
        !               TAU MULTIPLIED WITH TAU^3 [-]
        ! Andreas Jäger, November 2016
        !*****************************************************
        DEPFTTT = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFTTT = DEPFTTT + pol_term(i) * ti_depf(i)*(ti_depf(i) - 1.D0)*(ti_depf(i) - 2.D0)
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFTTT = DEPFTTT + exp_term(i) * ti_depf(j)*(ti_depf(j) - 1.D0)*(ti_depf(j) - 2.D0)
        end do

        !summation over the gaussian terms
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            !wie Thu
            DEPFTTT = DEPFTTT + gauss_term(i)*tau**3.d0 * ( &
                ((ti_depf(j)/tau + 2.d0*gbeta(i)*(tau-ggama(i)))**2 - 3.d0*ti_depf(j)/tau**2.d0 + 6.d0*gbeta(i)) &
                * (ti_depf(j)/tau + 2.d0*gbeta(i)*(tau-ggama(i))) + 2.d0*ti_depf(j)/tau**3.d0)
        end do

        DEPFUNCDER(9) = DEPFTTT     !write the result into the output vector
    end if
    !*****************************************************
    if (GETDER(10) == 1) then
        !*****************************************************  !wirklich falsch, also doppelt DDD, oder nur falsch beschriftet?
        !   DEPFDDD     3ND DERIVATIVE OF DF WITH RESPECT TO
        !               del MULTIPLIED WITH DEL^3 [-]
        ! Andreas Jäger, November 2016
        !*****************************************************
        DEPFDDD = 0.D0
        !summation over the polynomial terms
        do i = 1, I_pol_mix
            DEPFDDD = DEPFDDD + pol_term(i) * ti_depf(i) *&
                & ((di_depf(i) - dli(i)*li(i))**2 - di_depf(i) - dli(i)*li(i)*(li(i) - 1.D0))
        end do

        !summation over the exponential terms
        do i = 1, I_exp_mix
            j = I_pol_mix + i
            DEPFDDD = DEPFDDD + exp_term(i) * ti_depf(j) * ((di_depf(j) + del*(-2.D0*eta_depf(i)*(del - eps_depf(i)) - beta_depf(i)))**2 &
                - di_depf(j) - 2.D0*eta_depf(i)*del2)
        end do

        !summation over the gaussian terms !DDT
        do i = 1, I_gau_mix
            j = I_pol_mix + I_exp_mix + i
            DEPFDDD = DEPFDDD + gauss_term(i) &
                * (del**2.d0*tau*(((di_depf(j)/del + 2.d0*geta(i)*(del-geps(i)))**2 - di_depf(j)/del**2.d0 + 2.d0*geta(i)) &
                * (ti_depf(j)/tau + 2.d0*gbeta(i)*(tau-ggama(i))) ))
        end do

        DEPFUNCDER(10) = DEPFDDD     !write the result into the output vector
    end if
    !*****************************************************


    if (Density <= 1.d-12) then
        depfuncder(2) = depfuncder(2) / del
        depfuncder(3) = depfuncder(3) / (del**2.d0)
        depfuncder(8) = depfuncder(8) / (del**3.d0)
    endif

    end subroutine DEPFUNCFNR


    !!!Subroutine for the calculation of the mixture virial coefficient. Andreas September 2015
    !!subroutine VIR_MIX_CALC(Temp, B_mix, C_mix, D_mix)
    !!

    !!
    !!implicit none
    !!
    !!!Variables for new virial based mixing rules. Ignore for now. Andreas September 2015
    !!double precision, dimension(30,30) :: B_ij
    !!double precision, dimension(30,30,30) :: C_ijk
    !!double precision, dimension(30,30,30,30) :: D_ijkl
    !!
    !!Double precision:: B_mix, C_mix, D_mix
    !!
    !!Double precision:: B_calc, C_calc, D_calc, temp
    !!integer:: i,j,k,l
    !!
    !!B_mix = 0.D0
    !!C_mix = 0.D0
    !!D_mix = 0.D0
    !!
    !!    !Calculate the second, third, (and fourth) virial coefficient for all pure components in the mixture
    !!    do i=1, ncomp
    !!        B_ij(i,i) = B_calc(gl,temp, i)
    !!        C_ijk(i,i,i) = C_calc(gl,temp, i)
    !!        D_ijkl(i,i,i,i) = D_calc(gl,temp, i)
    !!    end do
    !!
    !!    !Only for binary mixtures as first try:
    !!    !Use arithmetic combination rule
    !!    B_ij(1,2) = (B_ij(1,1) + B_ij(2,2)) * 0.5D0
    !!
    !!    C_ijk(1,1,2) = (2.D0 * C_ijk(1,1,1) + C_ijk(2,2,2)) / 3.D0
    !!    C_ijk(1,2,2) = (C_ijk(1,1,1) + 2.D0 * C_ijk(2,2,2)) / 3.D0
    !!
    !!    D_ijkl(1,1,1,2) = (3.D0 * D_ijkl(1,1,1,1) + D_ijkl(2,2,2,2)) * 0.25D0
    !!    D_ijkl(1,1,2,2) = (2.D0 * D_ijkl(1,1,1,1) + 2.D0 * D_ijkl(2,2,2,2)) * 0.25D0
    !!    D_ijkl(1,2,2,2) = (D_ijkl(1,1,1,1) + 3.D0 * D_ijkl(2,2,2,2)) * 0.25D0
    !!
    !!    B_ij(1,2) = B_ij(2,1)
    !!
    !!    C_ijk(1,2,1) =  C_ijk(1,1,2)
    !!    C_ijk(1,1,2) =  C_ijk(1,1,2)
    !!    C_ijk(2,1,2) =  C_ijk(1,2,2)
    !!    C_ijk(2,2,1) =  C_ijk(1,2,2)
    !!
    !!    D_ijkl(1,1,2,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(1,2,1,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(2,1,1,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(1,1,2,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(1,2,1,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,1,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(1,2,2,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,2,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,2,1,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,2,2) = D_ijkl(1,2,2,2)
    !!    D_ijkl(2,2,1,2) = D_ijkl(1,2,2,2)
    !!    D_ijkl(2,2,2,1) = D_ijkl(1,2,2,2)
    !!
    !!    !Calculate the mixture virial coefficients
    !!    do i = 1, ncomp
    !!        do j = 1, ncomp
    !!            B_mix = B_mix + molfractions(i) * molfractions(j) * B_ij(i,j)
    !!            do k = 1, ncomp
    !!                C_mix = C_mix + molfractions(i) * molfractions(j) * molfractions(k) * C_ijk(i,j,k)
    !!                do l = 1, ncomp
    !!                    D_mix = D_mix + molfractions(i) * molfractions(j) * molfractions(k) * molfractions(l) * D_ijkl(i,j,k,l)
    !!                end do
    !!            end do
    !!        end do
    !!    end do
    !!
    !!!Set D to 0 --> Not correct yet!
    !!D_mix = 0.D0
    !!
    !!end subroutine
    !!
    !!
    !!!Subroutine for the calculation of the first derivative of the mixture virial coefficient with respect to xi. Andreas September 2015
    !!subroutine d_VIR_MIX_d_xi(Temp, d_B_mix_dxi, d_C_mix_dxi, d_D_mix_dxi)
    !!

    !!
    !!implicit none
    !!
    !!!Variables for new virial based mixing rules. Ignore for now. Andreas September 2015
    !!double precision, dimension(30,30) :: B_ij
    !!double precision, dimension(30,30,30) :: C_ijk
    !!double precision, dimension(30,30,30,30) :: D_ijkl
    !!
    !!Double precision, dimension(30):: d_B_mix_dxi, d_C_mix_dxi, d_D_mix_dxi
    !!
    !!Double precision:: B_calc, C_calc, D_calc, temp
    !!integer:: i,j,k,l
    !!
    !!    d_B_mix_dxi = 0.D0
    !!    d_C_mix_dxi = 0.D0
    !!    d_D_mix_dxi = 0.D0
    !!
    !!    !Calculate the second, third, (and fourth) virial coefficient for all pure components in the mixture
    !!    do i=1, ncomp
    !!        B_ij(i,i) = B_calc(gl,temp, i)
    !!        C_ijk(i,i,i) = C_calc(gl,temp, i)
    !!        D_ijkl(i,i,i,i) = D_calc(gl,temp, i)
    !!    end do
    !!
    !!    !Only for binary mixtures as first try:
    !!    !Use arithmetic combination rule
    !!    B_ij(1,2) = (B_ij(1,1) + B_ij(2,2)) * 0.5D0
    !!
    !!    C_ijk(1,1,2) = (2.D0 * C_ijk(1,1,1) + C_ijk(2,2,2)) / 3.D0
    !!    C_ijk(1,2,2) = (C_ijk(1,1,1) + 2.D0 * C_ijk(2,2,2)) / 3.D0
    !!
    !!    D_ijkl(1,1,1,2) = (3.D0 * D_ijkl(1,1,1,1) + D_ijkl(2,2,2,2)) * 0.25D0
    !!    D_ijkl(1,1,2,2) = (2.D0 * D_ijkl(1,1,1,1) + 2.D0 * D_ijkl(2,2,2,2)) * 0.25D0
    !!    D_ijkl(1,2,2,2) = (D_ijkl(1,1,1,1) + 3.D0 * D_ijkl(2,2,2,2)) * 0.25D0
    !!
    !!    B_ij(1,2) = B_ij(2,1)
    !!
    !!    C_ijk(1,2,1) =  C_ijk(1,1,2)
    !!    C_ijk(1,1,2) =  C_ijk(1,1,2)
    !!    C_ijk(2,1,2) =  C_ijk(1,2,2)
    !!    C_ijk(2,2,1) =  C_ijk(1,2,2)
    !!
    !!    D_ijkl(1,1,2,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(1,2,1,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(2,1,1,1) = D_ijkl(1,1,1,2)
    !!    D_ijkl(1,1,2,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(1,2,1,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,1,2) = D_ijkl(1,1,2,2)
    !!    D_ijkl(1,2,2,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,2,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,2,1,1) = D_ijkl(1,1,2,2)
    !!    D_ijkl(2,1,2,2) = D_ijkl(1,2,2,2)
    !!    D_ijkl(2,2,1,2) = D_ijkl(1,2,2,2)
    !!    D_ijkl(2,2,2,1) = D_ijkl(1,2,2,2)
    !!
    !!    !Calculate the first derivative of the mixture virial coefficients with respect to xi
    !!    do l = 1, ncomp -1
    !!        do i = 1, ncomp
    !!            d_B_mix_dxi(l) = d_B_mix_dxi(l) + molfractions(i) * (B_ij(i,l) - B_ij(i,ncomp))
    !!            do j = 1, ncomp
    !!                d_C_mix_dxi(l) = d_C_mix_dxi(l) + molfractions(i) * molfractions(j) * (C_ijk(i,j,l) - C_ijk(i,j,ncomp))
    !!                do k = 1, ncomp
    !!                    d_D_mix_dxi(l) = d_D_mix_dxi(l) + molfractions(i) * molfractions(j) * molfractions(k) * (D_ijkl(i,j,k,l) - D_ijkl(i,j,k,ncomp))
    !!                end do
    !!            end do
    !!        end do
    !!    end do
    !!    d_B_mix_dxi = d_B_mix_dxi * 2.D0
    !!    d_C_mix_dxi = d_C_mix_dxi * 3.D0
    !!    d_D_mix_dxi = d_D_mix_dxi * 4.D0
    !!
    !!!Set D to 0 --> Not correct yet!
    !!d_D_mix_dxi = 0.D0
    !!
    !!end subroutine







    !Never used: commented
    ! also variables in gl are commented !

    !!--------------------------------------------------------------------------------------
    !subroutine gE_NRTL(gl,Temp, gE, ln_gamma)
    !!--------------------------------------------------------------------------------------
    !!Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !!components in the mixture according to the non-random-two-liquid model.
    !!
    !!The model is described in
    !!Renon, H.; Prausnitz, J.M.: "Local Compositions in Thermodynamic Excess Functions for
    !!                             Liquid Mixtures", AIChe J. 14(1), 135-144, 1968
    !!
    !!
    !!INPUT
    !!Temp       -    Temperature in K
    !!Note that the model makes use of the global variable "molfractions" defined in the module
    !!"module_fluid_parameters"
    !!
    !!OUTPUT
    !!gE         -    Molar excess Gibbs energy in J/molK
    !!ln_gamma   -    Vector containing the natural logarithm of the activity coefficients
    !!
    !!Andreas Jäger, June 2016
    !!--------------------------------------------------------------------------------------
    !
    !
    !
    !
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !
    !double precision, intent(in) :: Temp
    !double precision :: gE
    !double precision, dimension(30) :: ln_gamma
    !double precision, dimension(30) :: x
    !
    !double precision:: R_const
    !
    !!Help variables
    !double precision:: sum_xj_tauji_Gji
    !double precision:: sum_xk_Gki
    !
    !double precision:: sum_xk_Gkj
    !double precision:: sum_xm_taumj_Gmj
    !
    !integer:: i
    !integer:: j
    !integer:: k
    !integer:: m
    !
    !!double precision, dimension(2) :: test_lngamma
    !
    !!Value for the ideal gas constant
    !R_const = 8.3144598D0
    !
    !!Get mole fractions of the mixture from module variable
    !x = gl%molfractions
    !
    !do j = 1, gl%ncomp
    !    do i = 1, gl%ncomp
    !        gl%tau_ij_NRTL(i,j) = gl%tau_Aij_NRTL(i,j) + gl%tau_Bij_NRTL(i,j) / Temp + gl%tau_Cij_NRTL(i,j) / Temp**2.d0 &
    !            & + gl%tau_Dij_NRTL(i,j) * dlog(Temp) + gl%tau_Eij_NRTL(i,j) * Temp**gl%tau_Fij_NRTL(i,j)
    !        gl%alpha_ij_NRTL(i,j) = gl%alpha_ij_0_NRTL(i,j) + gl%alpha_ij_1_NRTL(i,j) * Temp
    !        gl%G_ij_NRTL(i,j) = dexp(-gl%tau_ij_NRTL(i,j) * gl%alpha_ij_NRTL(i,j))
    !    end do
    !end do
    !
    !
    !gE = 0.D0
    !do i = 1, gl%ncomp
    !
    !    !First part of the sum, ln(gammai) is composed of
    !    !----
    !    sum_xj_tauji_Gji = 0.D0
    !    do j = 1, gl%ncomp
    !        sum_xj_tauji_Gji = sum_xj_tauji_Gji + x(j) * gl%tau_ij_NRTL(j,i) * gl%G_ij_NRTL(j,i)
    !    end do
    !    sum_xk_Gki = 0.D0
    !    do k = 1, gl%ncomp
    !        sum_xk_Gki = sum_xk_Gki + x(k) * gl%G_ij_NRTL(k,i)
    !    end do
    !    ln_gamma(i) = sum_xj_tauji_Gji / sum_xk_Gki
    !    !----
    !
    !    !Second part of the sum, ln(gammai) is composed of
    !    do j = 1, gl%ncomp
    !        sum_xk_Gkj = 0.D0
    !        do k = 1, gl%ncomp
    !            sum_xk_Gkj = sum_xk_Gkj + x(k) * gl%G_ij_NRTL(k,j)
    !        end do
    !        sum_xm_taumj_Gmj = 0.D0
    !        do m = 1, gl%ncomp
    !            sum_xm_taumj_Gmj = sum_xm_taumj_Gmj + x(m) * gl%tau_ij_NRTL(m,j) * gl%G_ij_NRTL(m,j)
    !        end do
    !
    !        ln_gamma(i) = ln_gamma(i) + x(j) * gl%G_ij_NRTL(i,j) / sum_xk_Gkj * (gl%tau_ij_NRTL(i,j) - sum_xm_taumj_Gmj / sum_xk_Gkj)
    !
    !    end do
    !
    !    gE = gE + R_const * Temp * x(i) * ln_gamma(i)
    !end do
    !
    !!Test for a binary mixture
    !!test_lngamma(1) = x(2)**2 * (tau_ij_NRTL(2,1) * (G_ij_NRTL(2,1) / (x(1) + x(2) * G_ij_NRTL(2,1)))**2  + tau_ij_NRTL(1,2) * G_ij_NRTL(1,2) / (x(2) + x(1) * G_ij_NRTL(1,2))**2 )
    !!test_lngamma(2) = x(1)**2 * (tau_ij_NRTL(1,2) * (G_ij_NRTL(1,2) / (x(2) + x(1) * G_ij_NRTL(1,2)))**2  + tau_ij_NRTL(2,1) * G_ij_NRTL(2,1) / (x(1) + x(2) * G_ij_NRTL(2,1))**2 )
    !
    !end subroutine
    !!--------------------------------------------------------------------------------------



    !------------------------------------------------------------------------------------------------------------------------------
    !subroutine gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, ln_gamma_R, C_or_R, errval)
    subroutine gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the UNIFAC model as well as tau- and delta-
    !derivatives of the model.
    !
    !The model is described by
    !Fredenslund, A.; Jones, R. L.; Prausnitz, J. M.: "Group-Contribution Estimation of
    ! Activity Coefficients in Nonideal Liquid Mixtures", AIChE Journal 21(6), 1086-1099, 1975.
    !
    !
    !INPUT
    !Temp       -    Temperature in K
    !Note that the model makes use of the global variable "molfractions" defined in the module
    !"module_fluid_parameters"
    ! GETDER      - AN ARRAY WITH nderivs ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. gE from UNIFAC as a function of temperature and composition
    !                2. 1ST DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                3. 2ND DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                4. 1ST DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                5: 2ND DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                6: 2ND MIXED DERIVATIVE OF gE WITH RESPECT TO DEL AND TAU
    !                7: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO DEL, TAU, AND TAU
    !                8: 3RD DERIVATIVE OF gE WITH RESPECT TO DEL
    !                9: 3RD DERIVATIVE OF gE WITH RESPECT TO TAU
    !               10: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, AND DEL
    !               11: 4TH DERIVATIVE OF gE WITH RESPECT TO DEL
    !               12: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, DEL AND DEL
    !               13: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, DEL, AND DEL
    !               14: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, TAU AND DEL
    !               15: 4TH DERIVATIVE OF gE WITH RESPECT TO TAU
    !C_or_R       - Integer that specifies if the combinatorial, the residual or both parts shall be calculated
    !               C_or_R = 0: Calculate both parts
    !               C_or_R = 1: Calculate the combinatorial part only
    !               C_or_R = 2: Calculate the the residual part only
    !
    !OUTPUT
    !gE_C         -    Vector with molar excess Gibbs energy of the combinatorial part (C) in J/molK and derivatives with respect to tau and delta
    !gE_R         -    Vector with molar excess Gibbs energy of the residual part (R) in J/molK and derivatives with respect to tau and delta
    !ln_gamma_C   -    Matrix containing the natural logarithm of the activity coefficients of the combinatorial part for each component and derivatives with respect to tau and delta
    !ln_gamma_R   -    Matrix containing the natural logarithm of the activity coefficients of the residual part for each component and derivatives with respect to tau and delta
    !errval       -    Indicates if an error occurred during calculations
    !
    !Andreas Jäger, January 2017
    !------------------------------------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision :: Temp
    integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs) :: gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval

    !Test variables for overall excess Gibbs energy and activity coefficients
    double precision, dimension(nderivs) :: gE          !Excess Gibbs energy
    double precision, dimension(nderivs,30) :: ln_gamma !activity coefficients

    double precision, dimension(30) :: r_i              !Molecular vdW volume of components i
    double precision, dimension(30) :: q_i              !Molecular vdW surface area of components i

    double precision :: R_const                         !universal gas constant

    !Help variables
    double precision, dimension(30) :: x
    double precision, dimension(30) :: phi_i
    double precision, dimension(30) :: Theta_i
    double precision, dimension(100) :: lnGamma_k
    double precision, dimension(:,:), allocatable :: lnGamma_ki   !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: Theta_m
    double precision, dimension(:,:), allocatable :: Theta_mi     !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: X_m             !Mole fraction of groups in the mixture
    double precision, dimension(:,:), allocatable :: X_mi         !Mole fraction of groups in pure component i. !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(:,:), allocatable :: Psi_nm

    !More help variables
    double precision:: sum_rjxj
    double precision:: sum_qjxj
    double precision:: sum_QnXn
    double precision:: sumsum_vnxj
    !double precision:: sum_xjLj
    double precision:: sum_ThetamPsimk
    double precision:: sumsum_ThetamPsikm
    double precision:: sum_ThetanPsinm

    !Derivatives of the combinatorial and residual part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho

    !Derivatives of the residual part with respect to "natural" variables
    double precision, dimension(:), allocatable:: dlnGamma_k_dT
    double precision, dimension(:,:), allocatable :: dlnGamma_ki_dT
    double precision, dimension(:,:), allocatable :: dPsi_nm_dT
    double precision, dimension(:), allocatable :: d2lnGamma_k_dT2
    double precision, dimension(:,:), allocatable:: d2lnGamma_ki_dT2
    double precision, dimension(:,:), allocatable :: d2Psi_nm_dT2
    double precision:: sum_Thetam_dPsimk_dT
    double precision:: sum_Thetan_dPsinm_dT
    double precision:: sum_Thetam_d2Psimk_dT2
    double precision:: sum_Thetan_d2Psinm_dT2
    double precision:: sumsum_ThetamPsikm_dT        !Variable for double sum in first temperature derivative
    double precision:: sumsum_ThetamPsikm_dT2       !Variable for double sum in second temperature derivative

    !Corresponding states relevant variables
    double precision:: tau, del

    !Summation variables
    integer:: i, j, k, m, n, help_k

    !Variable that indicates if derivatives need to be recalculated
    logical:: recalc

    if (.not.(allocated(dlnGamma_k_dT)))      allocate(dlnGamma_k_dT(100))                                                                                    !    allocate(dlnGamma_k_dT((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(dlnGamma_ki_dT)))     allocate(dlnGamma_ki_dT(100,30))                                                                          !    allocate(dlnGamma_ki_dT((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix,gl%ncomp))
    if (.not.(allocated(dPsi_nm_dT)))         allocate(dPsi_nm_dT(100,100))          !    allocate(dPsi_nm_dT((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix,(gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(d2lnGamma_k_dT2)))    allocate(d2lnGamma_k_dT2(100))                                                                                  !    allocate(d2lnGamma_k_dT2((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(d2lnGamma_ki_dT2)))   allocate(d2lnGamma_ki_dT2(100,30))                                                                        !    allocate(d2lnGamma_ki_dT2((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix,gl%ncomp))
    if (.not.(allocated(d2Psi_nm_dT2)))       allocate(d2Psi_nm_dT2(100, 100))       !    allocate(d2Psi_nm_dT2((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix, (gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    !
    !
    if (.not.(allocated(lnGamma_ki)))       allocate(lnGamma_ki(100, 100))           !  allocate(lnGamma_ki((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix, (gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(Theta_mi)))         allocate(Theta_mi(100, 100))             !  allocate(Theta_mi((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix, (gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(X_mi)))             allocate(X_mi(100, 100))                 !  allocate(X_mi((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix, (gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))
    if (.not.(allocated(Psi_nm)))           allocate(Psi_nm(100, 100))               !  allocate(Psi_nm((gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix, (gl%ncomp * gl%nr_of_groups_mix - gl%nr_of_groups_mix) * gl%nr_of_groups_mix))

    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C = errval
        gl%ge%gE_R = errval
        gl%ge%ln_gamma_C = errval
        gl%ge%ln_gamma_R = errval
        return
    end if

    recalc = .true.

    !Check whether calculations have to be redone or not.
    !If the temperature, (density), composition, equation of state and, the required part of UNIFAC stayed the same, no need to recalculate
    if(dabs(gl%ge%Temp_prev - Temp) > 1.D-16) then
        recalc = .true.
    end if
    do i=1,gl%ncomp
        if(dabs(gl%ge%molfrac_prev(i) - gl%molfractions(i)) > 1.D-16) then
            recalc = .true.
        end if
        if(gl%ge%Eq_type_prev(i) .ne. gl%Eq_type(i)) then
            recalc = .true.
        end if
    end do
    if (gl%ge%mixtype_prev .ne. gl%mix_type) then
        recalc = .true.
    end if
    if (gl%ge%C_or_R_prev .ne. C_or_R) then
        recalc = .true.
    end if

    !Initialize variables
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.D0
    gl%ge%ln_gamma_C = 0.D0
    gl%ge%ln_gamma_R = 0.D0
    gl%ge%gE_C = 0.D0
    gl%ge%gE_R = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0

    !Get mole fractions of the mixture from module variable
    x = gl%molfractions

    !Calculate tau
    tau = gl%tredmix / Temp

    !Set dummy value for delta (Not needed because UNIFAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0

    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
        !Calculate the combinatorial part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_C is always needed for the following calculations, thus it is always calculated

        if ((gl%ge%GETDER_prev(1) == 1) .and. (recalc .eqv. .false.)) then

            gl%ge%ln_gamma_C(1,:) = gl%ge%ln_gamma_C_prev(1,:)

        else

            !Calculate help variables
            r_i = 0.D0
            q_i = 0.D0
            Do i = 1, gl%ncomp
                Do k = 1, gl%nr_of_groups_i(i)
                    r_i(i) = r_i(i) + gl%v_ik(i, k) * gl%R_ik(i, k)
                    q_i(i) = q_i(i) + gl%v_ik(i, k) * gl%Q_ik(i, k)
                End do
            End do

            sum_rjxj = 0.D0
            sum_qjxj = 0.D0
            Do j = 1, gl%ncomp
                sum_rjxj = sum_rjxj + r_i(j) * x(j)
                sum_qjxj = sum_qjxj + q_i(j) * x(j)
            End do
            Do i = 1, gl%ncomp
                phi_i(i) = r_i(i) / sum_rjxj
                Theta_i(i) = q_i(i) / sum_qjxj
            End do


            Do i = 1, gl%ncomp
                gl%ge%ln_gamma_C(1,i) = 1.D0 - phi_i(i) + dlog(phi_i(i)) - &
                    &  5.D0 * q_i(i) * (1 - phi_i(i) / Theta_i(i) + dlog(phi_i(i) / Theta_i(i)))
            End do

        end if

        !Theoretically not needed to copy the combinatorial gammas to the Trho derivatives variables (ln_gamma_C_Trho)
        !However, done anyways for safety and consistency
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)

        if (GETDER(1) .eq. 1) then
            if ((gl%ge%GETDER_prev(1) == 1) .and. (recalc .eqv. .false.)) then
                gl%ge%gE_C(1) = gl%ge%gE_C_prev(1)
            else
                Do i = 1, gl%ncomp
                    gl%ge%gE_C(1) = gl%ge%gE_C(1) + x(i) * gl%ge%ln_gamma_C(1,i) * R_const * Temp
                end do
            end if
            gE_C_Trho(1) = gl%ge%gE_C(1)  !Also theoretically not needed
        end if

        !All other derivatives of the combinatorial part are 0 (except for the tau-only derivatives of gE), because the combinatorial part is a function of composition only
        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_C_Trho(2,:) = 0.D0
            gE_C_Trho(2) = 0.d0
            gl%ge%ln_gamma_C(2,:) = 0.d0
            gl%ge%gE_C(2) = 0.D0
        end if
        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_C_Trho(3,:) = 0.D0
            gE_C_Trho(3) = 0.d0
            gl%ge%ln_gamma_C(3,:) = 0.d0
            gl%ge%gE_C(3) = 0.D0
        end if
        !Derivative wrt tau
        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then
            if ((gl%ge%GETDER_prev(4) == 1) .and. (recalc .eqv. .false.)) then   !Don't recalculate, take old solution
                gl%ge%gE_C(4) = gl%ge%gE_C_prev(4)
                gl%ge%ln_gamma_C(4,:) = gl%ge%ln_gamma_C_prev(4,:)
                ln_gamma_C_Trho(4,:) = 0.D0     !ln_gamma_C is not a function of temperature
                !Calculate the derivatives of gE_C with respect to the natural variables from the old derivatives with respect to tau and delta
                gE_C_Trho(4) = - tau**2.d0 / gl%tredmix * gl%ge%gE_C_prev(4)
            else
                ln_gamma_C_Trho(4,:) = 0.D0     !ln_gamma_C is not a function of temperature
                Do i = 1, gl%ncomp
                    gE_C_Trho(4) = gE_C_Trho(4) + x(i) * gl%ge%ln_gamma_C(1,i) * R_const
                end do
                gl%ge%ln_gamma_C(4,:) = - gl%tredmix / tau**2.d0 * ln_gamma_C_Trho(4,:)     !=0, because ln_gamma_C is not a function of temperature. For clearness sake nevertheless calculated here
                gl%ge%gE_C(4) = - gl%tredmix / tau**2.d0 * gE_C_Trho(4)
            end if
        end if
        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            if ((gl%ge%GETDER_prev(5) == 1) .and. (recalc .eqv. .false.)) then   !Don't recalculate, take old solution
                ln_gamma_C_Trho(5,:) = 0.D0
                gE_C_Trho(5) = 0.d0
                gl%ge%ln_gamma_C(5,:) = 0.d0
                gl%ge%gE_C(5) = gl%ge%gE_C_prev(5)
            else
                ln_gamma_C_Trho(5,:) = 0.D0
                gE_C_Trho(5) = 0.d0
                gl%ge%ln_gamma_C(5,:) = 0.d0
                gl%ge%gE_C(5) = 2.D0*gl%tredmix/tau**3.d0*gE_C_Trho(4)+gl%tredmix**2.d0/tau**4.d0*gE_C_Trho(5)
            end if
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_C_Trho(6,:) = 0.D0
            gE_C_Trho(6) = 0.d0
            gl%ge%ln_gamma_C(6,:) = 0.d0
            gl%ge%gE_C(6) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_C_Trho(7,:) = 0.D0
            gE_C_Trho(7) = 0.d0
            gl%ge%ln_gamma_C(7,:) = 0.d0
            gl%ge%gE_C(7) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_C_Trho(8,:) = 0.D0
            gE_C_Trho(8) = 0.d0
            gl%ge%ln_gamma_C(8,:) = 0.d0
            gl%ge%gE_C(8) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_C_Trho(9,:) = 0.D0
            gE_C_Trho(9) = 0.d0
            gl%ge%ln_gamma_C(9,:) = 0.d0
            gl%ge%gE_C(9) = 0.D0      !TO BE IMPLEMENTED
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_C_Trho(10,:) = 0.D0
            gE_C_Trho(10) = 0.d0
            gl%ge%ln_gamma_C(10,:) = 0.d0
            gl%ge%gE_C(10) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_C_Trho(11,:) = 0.D0
            gE_C_Trho(11) = 0.d0
            gl%ge%ln_gamma_C(11,:) = 0.d0
            gl%ge%gE_C(11) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_C_Trho(12,:) = 0.D0
            gE_C_Trho(12) = 0.d0
            gl%ge%ln_gamma_C(12,:) = 0.d0
            gl%ge%gE_C(12) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_C_Trho(13,:) = 0.D0
            gE_C_Trho(13) = 0.d0
            gl%ge%ln_gamma_C(13,:) = 0.d0
            gl%ge%gE_C(13) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_C_Trho(14,:) = 0.D0
            gE_C_Trho(14) = 0.d0
            gl%ge%ln_gamma_C(14,:) = 0.d0
            gl%ge%gE_C(14) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_C_Trho(15,:) = 0.D0
            gE_C_Trho(15) = 0.d0
            gl%ge%ln_gamma_C(15,:) = 0.d0
            gl%ge%gE_C(15) = 0.D0     !TO BE IMPLEMENTED
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end if



    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !Calculate the residual part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_R is always needed for the following calculations, thus it is always calculated

        !----------------------------------------------------------------------
        !Calculate the mole fractions of all groups in the mixture
        !----------------------------------------------------------------------
        sumsum_vnxj = 0.D0
        Do j = 1, gl%ncomp                     !sum over all components in the mixture
            Do n = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                sumsum_vnxj = sumsum_vnxj + gl%v_ik(j, n) * x(j)
            End do
        End do

        m = 1
        Do j = 1, gl%ncomp                     !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                X_m(m) = gl%v_ik(j, k) * x(j) / sumsum_vnxj
                m = m + 1
            End do
        End do
        !----------------------------------------------------------------------


        !----------------------------------------------------------------------
        !Calculate Theta_m for all groups in the mixture
        !----------------------------------------------------------------------
        m = 1
        sum_QnXn = 0.D0
        Do j = 1, gl%ncomp                      !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                sum_QnXn = sum_QnXn + gl%Q_ik(j, k) * X_m(m)
                m = m + 1
            End do
        End do

        m = 1
        Do j = 1, gl%ncomp                      !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                Theta_m(m) = gl%Q_ik(j, k) * X_m(m) / sum_QnXn
                m = m + 1
            End do
        End do
        !----------------------------------------------------------------------



        !----------------------------------------------------------------------
        !Calculate all group interaction parameters Psi_nm
        !----------------------------------------------------------------------
        Do n = 1, gl%nr_of_groups_mix
            Do m = 1, gl%nr_of_groups_mix
                Psi_nm(n, m) = dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2.d0) / Temp)
            End do
        End do
        !----------------------------------------------------------------------


        !Loop over all components in the mixture
        help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
        Do i = 1, gl%ncomp

            !----------------------------------------------------------------------
            !Calculate the mole fractions X_mi of the groups of PURE component i
            sumsum_vnxj = 0.D0
            Do n = 1, gl%nr_of_groups_i(i)
                sumsum_vnxj = sumsum_vnxj + gl%v_ik(i, n) * 1.D0     !For pure component i, the mole fraction is equal to 1!!
            End do

            m = 1
            Do k = 1, gl%nr_of_groups_i(i)
                X_mi(m,i) = gl%v_ik(i, k) * 1.D0 / sumsum_vnxj           !For pure component i, the mole fraction is equal to 1!!
                m = m + 1
            End do
            !----------------------------------------------------------------------

            !----------------------------------------------------------------------
            !Calculate Theta_mi for all groups of PURE component i
            m = 1
            sum_QnXn = 0.D0
            Do k = 1, gl%nr_of_groups_i(i)  !Sum over all groups of pure component i
                sum_QnXn = sum_QnXn + gl%Q_ik(i, k) * X_mi(m,i)
                m = m + 1
            End do

            m = 1
            Do k = 1, gl%nr_of_groups_i(i)  !Sum over all groups of pure component i
                Theta_mi(m,i) = gl%Q_ik(i, k) * X_mi(m,i) / sum_QnXn
                m = m + 1
            End do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            !Calculate lnGamma_k
            !----------------------------------------------------------------------
            Do k = 1, gl%nr_of_groups_i(i)
                sum_ThetamPsimk = 0.D0
                sumsum_ThetamPsikm = 0.D0
                help_k = help_k + 1
                Do m = 1, gl%nr_of_groups_mix
                    sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                    sum_ThetanPsinm = 0.D0
                    Do n = 1, gl%nr_of_groups_mix
                        sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                    End do
                    sumsum_ThetamPsikm = sumsum_ThetamPsikm + Theta_m(m) * Psi_nm(help_k, m) / sum_ThetanPsinm
                End do
                lnGamma_k(k) = gl%Q_ik(i, k) * (1.D0 - dlog(sum_ThetamPsimk) - sumsum_ThetamPsikm)
            End do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            !Calculate lnGamma_ki (PURE component i)
            !----------------------------------------------------------------------
            help_k = help_k - gl%nr_of_groups_i(i)     !Set back help_k for the calculation of lnGamma_ki
            Do k = 1, gl%nr_of_groups_i(i)
                sum_ThetamPsimk = 0.D0
                sumsum_ThetamPsikm = 0.D0
                help_k = help_k + 1
                Do m = 1, gl%nr_of_groups_i(i)
                    sum_ThetamPsimk = sum_ThetamPsimk + Theta_mi(m,i) * Psi_nm(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                    sum_ThetanPsinm = 0.D0
                    Do n = 1, gl%nr_of_groups_i(i)
                        sum_ThetanPsinm = sum_ThetanPsinm + Theta_mi(n,i) * Psi_nm(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                    End do
                    sumsum_ThetamPsikm = sumsum_ThetamPsikm + Theta_mi(m,i) * Psi_nm(help_k, m + gl%group_start_index(i) - 1) / sum_ThetanPsinm
                End do
                lnGamma_ki(k,i) = gl%Q_ik(i, k) * (1.D0 - dlog(sum_ThetamPsimk) - sumsum_ThetamPsikm)
            End do
            !----------------------------------------------------------------------

            Do k = 1, gl%nr_of_groups_i(i)
                gl%ge%ln_gamma_R(1,i) = gl%ge%ln_gamma_R(1,i) + gl%v_ik(i, k) * (lnGamma_k(k) - lnGamma_ki(k,i))
            End do
        End do

        !Theoretically not needed to copy the residual gammas to the Trho derivatives variables (ln_gamma_R_Trho)
        !However, done anyways for safety and consistency
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if (GETDER(1) .eq. 1) then
            if ((gl%ge%GETDER_prev(1) == 1) .and. (recalc .eqv. .false.)) then
                gl%ge%gE_R(1) = gl%ge%gE_R_prev(1)
            else
                Do i = 1, gl%ncomp
                    gl%ge%gE_R(1) = gl%ge%gE_R(1) + x(i) * gl%ge%ln_gamma_R(1,i) * R_const * Temp
                end do
            end if
            gE_R_Trho(1) = gl%ge%gE_R(1)  !Not needed, but for consistency written in this variable
        end if

        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho(2,:) = 0.D0
            gE_R_Trho(2) = 0.d0
            gl%ge%ln_gamma_R(2,:) = 0.d0
            gl%ge%gE_R(2) = 0.D0
        end if

        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho(3,:) = 0.D0
            gE_R_Trho(3) = 0.d0
            gl%ge%ln_gamma_R(3,:) = 0.d0
            gl%ge%gE_R(3) = 0.D0
        end if

        !Derivative wrt tau
        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then
            !---------------------------------------------------------------------------------------
            !Calculate the first temperature derivative of all group interaction parameters Psi_nm
            !Always calculate this, because it is needed for all following tau derivatives too
            !---------------------------------------------------------------------------------------
            Do n = 1, gl%nr_of_groups_mix
                Do m = 1, gl%nr_of_groups_mix
                    dPsi_nm_dT(n, m) = (gl%a_nm(n, m) / Temp**2.d0 - gl%c_nm(n, m)) * &
                        & dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2.d0) / Temp)
                End do
            End do
            !---------------------------------------------------------------------------------------

            if ((gl%ge%GETDER_prev(4) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%ln_gamma_R(4,:) = gl%ge%ln_gamma_R_prev(4,:)
                gl%ge%gE_R(4) = gl%ge%gE_R_prev(4)

                !Calculate the temperature derivatives from the known tau derivatives
                ln_gamma_R_Trho(4,:) = - tau**2.d0 / gl%tredmix * gl%ge%ln_gamma_R(4,:)
                gE_R_Trho(4) = - tau**2.d0 / gl%tredmix * gl%ge%gE_R(4)

            else

                !Loop over all components in the mixture
                help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
                Do i = 1, gl%ncomp

                    !----------------------------------------------------------------------
                    !Calculate d(lnGamma_k)/dT
                    !Here, it is theoretically not necessary to recalculate the following
                    !variables:
                    !sum_ThetamPsimk
                    !sum_ThetanPsinm
                    !However, for clearness sake and for better readability these sum are
                    !recalculated here. This could be changed by saving these sums in the
                    !calculations above and make these scalar sums vectors or matrices
                    !Andreas Jäger, April 2017
                    !----------------------------------------------------------------------
                    Do k = 1, gl%nr_of_groups_i(i)
                        sum_ThetamPsimk = 0.D0
                        sumsum_ThetamPsikm_dT = 0.D0
                        sum_Thetam_dPsimk_dT = 0.D0
                        help_k = help_k + 1
                        Do m = 1, gl%nr_of_groups_mix
                            sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                            sum_Thetam_dPsimk_dT = sum_Thetam_dPsimk_dT + Theta_m(m) * dPsi_nm_dT(m, help_k)
                            sum_ThetanPsinm = 0.D0
                            sum_Thetan_dPsinm_dT = 0.D0
                            Do n = 1, gl%nr_of_groups_mix
                                sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                                sum_Thetan_dPsinm_dT = sum_Thetan_dPsinm_dT + Theta_m(n) * dPsi_nm_dT(n, m)
                            End do
                            sumsum_ThetamPsikm_dT = sumsum_ThetamPsikm_dT + &
                                & ( Theta_m(m) * dPsi_nm_dT(help_k, m) * sum_ThetanPsinm &
                                & - Theta_m(m) * Psi_nm(help_k, m) * sum_Thetan_dPsinm_dT ) &
                                & / (sum_ThetanPsinm)**2
                        End do
                        dlnGamma_k_dT(k) = -gl%Q_ik(i, k) * sum_Thetam_dPsimk_dT / sum_ThetamPsimk &
                            & -gl%Q_ik(i, k) * sumsum_ThetamPsikm_dT
                    End do
                    !----------------------------------------------------------------------


                    !----------------------------------------------------------------------
                    !Calculate d(lnGamma_ki)/dT (PURE component i)
                    !----------------------------------------------------------------------
                    help_k = help_k - gl%nr_of_groups_i(i)     !Set back help_k for the calculation of lnGamma_ki
                    Do k = 1, gl%nr_of_groups_i(i)
                        sum_ThetamPsimk = 0.D0
                        sumsum_ThetamPsikm_dT = 0.D0
                        sum_Thetam_dPsimk_dT = 0.D0
                        help_k = help_k + 1
                        Do m = 1, gl%nr_of_groups_i(i)
                            sum_ThetamPsimk = sum_ThetamPsimk + Theta_mi(m,i) * Psi_nm(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                            sum_Thetam_dPsimk_dT = sum_Thetam_dPsimk_dT + Theta_mi(m,i) * dPsi_nm_dT(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                            sum_ThetanPsinm = 0.D0
                            sum_Thetan_dPsinm_dT = 0.D0
                            Do n = 1, gl%nr_of_groups_i(i)
                                sum_ThetanPsinm = sum_ThetanPsinm + Theta_mi(n,i) * Psi_nm(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                                sum_Thetan_dPsinm_dT = sum_Thetan_dPsinm_dT + Theta_mi(n,i) * dPsi_nm_dT(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                            End do
                            sumsum_ThetamPsikm_dT = sumsum_ThetamPsikm_dT + &
                                & ( Theta_mi(m,i) * dPsi_nm_dT(help_k, m + gl%group_start_index(i) - 1) * sum_ThetanPsinm &
                                & - Theta_mi(m,i) * Psi_nm(help_k, m + gl%group_start_index(i) - 1) * sum_Thetan_dPsinm_dT ) &
                                & / (sum_ThetanPsinm)**2
                        End do
                        dlnGamma_ki_dT(k,i) = -gl%Q_ik(i, k) * sum_Thetam_dPsimk_dT / sum_ThetamPsimk &
                            & -gl%Q_ik(i, k) * sumsum_ThetamPsikm_dT
                    End do
                    !----------------------------------------------------------------------


                    Do k = 1, gl%nr_of_groups_i(i)
                        ln_gamma_R_Trho(4,i) = ln_gamma_R_Trho(4,i) + gl%v_ik(i, k) * (dlnGamma_k_dT(k) - dlnGamma_ki_dT(k,i))
                    End do
                End do

                Do i = 1, gl%ncomp
                    gE_R_Trho(4) = gE_R_Trho(4) + R_Const * x(i) * gl%ge%ln_gamma_R(1,i) + R_Const * x(i) * Temp * ln_gamma_R_Trho(4,i)
                end do

                !Deriatives with respect to tau at constant delta and x
                gl%ge%ln_gamma_R(4,:) = - gl%tredmix / tau**2.d0 * ln_gamma_R_Trho(4,:)
                gl%ge%gE_R(4) = - gl%tredmix / tau**2.d0 * gE_R_Trho(4)

            end if
        end if


        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            !---------------------------------------------------------------------------------------
            !Calculate the second temperature derivative of all group interaction parameters Psi_nm
            !Always calculate this, because it is needed for all following tau derivatives too
            !---------------------------------------------------------------------------------------
            Do n = 1, gl%nr_of_groups_mix
                Do m = 1, gl%nr_of_groups_mix
                    d2Psi_nm_dT2(n, m) = (-2.D0 * gl%a_nm(n, m) / Temp**3.d0 + (gl%a_nm(n, m) / Temp**2.d0 - gl%c_nm(n, m))**2) * &
                        & dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2.d0) / Temp)
                End do
            End do
            !---------------------------------------------------------------------------------------

            if ((gl%ge%GETDER_prev(5) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%ln_gamma_R(5,:) = gl%ge%ln_gamma_R_prev(5,:)
                gl%ge%gE_R(5) = gl%ge%gE_R_prev(5)

                !Calculate the temperature derivatives from the known tau derivatives
                ln_gamma_R_Trho(5,:) = tau**4.d0 / gl%tredmix**2.d0 * gl%ge%ln_gamma_R(5,:) - 2.D0 * tau / gl%tredmix * ln_gamma_R_Trho(4,:)
                gE_R_Trho(5) = tau**4.d0 / gl%tredmix**2.d0 * gl%ge%gE_R(5) - 2.D0 * tau / gl%tredmix * gE_R_Trho(4)

            else


                !Loop over all components in the mixture
                help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
                Do i = 1, gl%ncomp

                    !----------------------------------------------------------------------
                    !Calculate d2(lnGamma_k)/dT2
                    !Here, it is theoretically not necessary to recalculate the following
                    !variables:
                    !sum_ThetamPsimk
                    !sum_ThetanPsinm
                    !sum_Thetam_dPsimk_dT
                    !sum_Thetan_dPsinm_dT
                    !However, for clearness sake and for better readability these sum are
                    !recalculated here. This could be changed by saving these sums in the
                    !calculations above and make these scalar sums vectors or matrices
                    !Andreas Jäger, April 2017
                    !----------------------------------------------------------------------
                    Do k = 1, gl%nr_of_groups_i(i)
                        sum_ThetamPsimk = 0.D0
                        sumsum_ThetamPsikm_dT = 0.D0
                        sum_Thetam_dPsimk_dT = 0.D0
                        sum_Thetam_d2Psimk_dT2 = 0.D0
                        sumsum_ThetamPsikm_dT2 = 0.D0
                        help_k = help_k + 1
                        Do m = 1, gl%nr_of_groups_mix
                            sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                            sum_Thetam_dPsimk_dT = sum_Thetam_dPsimk_dT + Theta_m(m) * dPsi_nm_dT(m, help_k)
                            sum_Thetam_d2Psimk_dT2 = sum_Thetam_d2Psimk_dT2 + Theta_m(m) * d2Psi_nm_dT2(m, help_k)
                            sum_ThetanPsinm = 0.D0
                            sum_Thetan_dPsinm_dT = 0.D0
                            sum_Thetan_d2Psinm_dT2 = 0.D0
                            Do n = 1, gl%nr_of_groups_mix
                                sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                                sum_Thetan_dPsinm_dT = sum_Thetan_dPsinm_dT + Theta_m(n) * dPsi_nm_dT(n, m)
                                sum_Thetan_d2Psinm_dT2 = sum_Thetan_d2Psinm_dT2 + Theta_m(n) * d2Psi_nm_dT2(n, m)
                            End do
                            sumsum_ThetamPsikm_dT2 = sumsum_ThetamPsikm_dT2 + &
                                & (( Theta_m(m) * d2Psi_nm_dT2(help_k, m) * sum_ThetanPsinm &
                                & - Theta_m(m) * Psi_nm(help_k, m) * sum_Thetan_d2Psinm_dT2 ) &
                                & / (sum_ThetanPsinm)**2) &
                                & -(2.D0 * sum_Thetan_dPsinm_dT * (Theta_m(m) * dPsi_nm_dT(help_k, m) * sum_ThetanPsinm &
                                & - Theta_m(m) * Psi_nm(help_k, m) * sum_Thetan_dPsinm_dT ) &
                                & / (sum_ThetanPsinm)**3)
                        End do
                        d2lnGamma_k_dT2(k) = gl%Q_ik(i, k) * sum_Thetam_dPsimk_dT**2.d0 / sum_ThetamPsimk**2.d0 &
                            & - gl%Q_ik(i, k) * sum_Thetam_d2Psimk_dT2 / sum_ThetamPsimk &
                            & - gl%Q_ik(i, k) * sumsum_ThetamPsikm_dT2
                    End do
                    !----------------------------------------------------------------------


                    !----------------------------------------------------------------------
                    !Calculate d2(lnGamma_ki)/dT2 (PURE component i)
                    !----------------------------------------------------------------------
                    help_k = help_k - gl%nr_of_groups_i(i)     !Set back help_k for the calculation of lnGamma_ki
                    Do k = 1, gl%nr_of_groups_i(i)
                        sum_ThetamPsimk = 0.D0
                        sumsum_ThetamPsikm_dT = 0.D0
                        sum_Thetam_dPsimk_dT = 0.D0
                        sum_Thetam_d2Psimk_dT2 = 0.D0
                        sumsum_ThetamPsikm_dT2 = 0.D0
                        help_k = help_k + 1
                        Do m = 1, gl%nr_of_groups_i(i)
                            sum_ThetamPsimk = sum_ThetamPsimk + Theta_mi(m,i) * Psi_nm(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                            sum_Thetam_dPsimk_dT = sum_Thetam_dPsimk_dT + Theta_mi(m,i) * dPsi_nm_dT(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                            sum_Thetam_d2Psimk_dT2 = sum_Thetam_d2Psimk_dT2 + Theta_mi(m,i) * d2Psi_nm_dT2(m + gl%group_start_index(i) - 1, help_k)  !The index "m + group_start_index(i) - 1" of Psi makes sure that ONLY THE PURE FLUID interactions are used. This is necessary because Psi contains the interactions of ALL GROUPS IN THE MIXTURE
                            sum_ThetanPsinm = 0.D0
                            sum_Thetan_dPsinm_dT = 0.D0
                            sum_Thetan_d2Psinm_dT2 = 0.D0
                            Do n = 1, gl%nr_of_groups_i(i)
                                sum_ThetanPsinm = sum_ThetanPsinm + Theta_mi(n,i) * Psi_nm(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                                sum_Thetan_dPsinm_dT = sum_Thetan_dPsinm_dT + Theta_mi(n,i) * dPsi_nm_dT(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                                sum_Thetan_d2Psinm_dT2 = sum_Thetan_d2Psinm_dT2 + Theta_mi(n,i) * d2Psi_nm_dT2(n + gl%group_start_index(i) - 1, m + gl%group_start_index(i) - 1)
                            End do
                            sumsum_ThetamPsikm_dT = sumsum_ThetamPsikm_dT + &
                                & (( Theta_mi(m,i) * d2Psi_nm_dT2(help_k, m) * sum_ThetanPsinm &
                                & - Theta_mi(m,i) * Psi_nm(help_k, m) * sum_Thetan_d2Psinm_dT2 ) &
                                & / (sum_ThetanPsinm)**2) &
                                & -(2.D0 * sum_Thetan_dPsinm_dT * (Theta_mi(m,i) * dPsi_nm_dT(help_k, m) * sum_ThetanPsinm &
                                & - Theta_mi(m,i) * Psi_nm(help_k, m) * sum_Thetan_dPsinm_dT ) &
                                & / (sum_ThetanPsinm)**3)
                        End do
                        d2lnGamma_ki_dT2(k,i) = gl%Q_ik(i, k) * sum_Thetam_dPsimk_dT**2.d0 / sum_ThetamPsimk**2.d0 &
                            & - gl%Q_ik(i, k) * sum_Thetam_d2Psimk_dT2 / sum_ThetamPsimk &
                            & - gl%Q_ik(i, k) * sumsum_ThetamPsikm_dT2
                    End do
                    !----------------------------------------------------------------------


                    Do k = 1, gl%nr_of_groups_i(i)
                        ln_gamma_R_Trho(5,i) = ln_gamma_R_Trho(5,i) + gl%v_ik(i, k) * (d2lnGamma_k_dT2(k) - d2lnGamma_ki_dT2(k,i))
                    End do
                End do


                Do i = 1, gl%ncomp
                    gE_R_Trho(5) = gE_R_Trho(5) + 2.D0 * R_Const * x(i) * ln_gamma_R_Trho(4,i) + R_Const * Temp * x(i) * ln_gamma_R_Trho(5,i)
                end do
                gl%ge%ln_gamma_R(5,:) = 2.D0*gl%tredmix/tau**3.d0*ln_gamma_R_Trho(4,:)+gl%tredmix**2.d0/tau**4.d0*ln_gamma_R_Trho(5,:)
                gl%ge%gE_R(5) = 2.D0*gl%tredmix/tau**3.d0*gE_R_Trho(4)+gl%tredmix**2.d0/tau**4.d0*gE_R_Trho(5)

            end if

        end if



        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho(6,:) = 0.D0
            gE_R_Trho(6) = 0.d0
            gl%ge%ln_gamma_R(6,:) = 0.d0
            gl%ge%gE_R(6) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho(7,:) = 0.D0
            gE_R_Trho(7) = 0.d0
            gl%ge%ln_gamma_R(7,:) = 0.d0
            gl%ge%gE_R(7) = 0.D0
        end if

        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho(8,:) = 0.D0
            gE_R_Trho(8) = 0.d0
            gl%ge%ln_gamma_R(8,:) = 0.d0
            gl%ge%gE_R(8) = 0.D0
        end if

        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho(9,:) = 0.D0
            gE_R_Trho(9) = 0.d0
            gl%ge%ln_gamma_R(9,:) = 0.d0
            gl%ge%gE_R(9) = 0.D0      !TO BE IMPLEMENTED
        end if

        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho(10,:) = 0.D0
            gE_R_Trho(10) = 0.d0
            gl%ge%ln_gamma_R(10,:) = 0.d0
            gl%ge%gE_R(10) = 0.D0
        end if

        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho(11,:) = 0.D0
            gE_R_Trho(11) = 0.d0
            gl%ge%ln_gamma_R(11,:) = 0.d0
            gl%ge%gE_R(11) = 0.D0
        end if

        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho(12,:) = 0.D0
            gE_R_Trho(12) = 0.d0
            gl%ge%ln_gamma_R(12,:) = 0.d0
            gl%ge%gE_R(12) = 0.D0
        end if

        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho(13,:) = 0.D0
            gE_R_Trho(13) = 0.d0
            gl%ge%ln_gamma_R(13,:) = 0.d0
            gl%ge%gE_R(13) = 0.D0
        end if

        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho(14,:) = 0.D0
            gE_R_Trho(14) = 0.d0
            gl%ge%ln_gamma_R(14,:) = 0.d0
            gl%ge%gE_R(14) = 0.D0
        end if

        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho(15,:) = 0.D0
            gE_R_Trho(15) = 0.d0
            gl%ge%ln_gamma_R(15,:) = 0.d0
            gl%ge%gE_R(15) = 0.D0     !TO BE IMPLEMENTED
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    End if

    !Calculate the overall activity coefficient and the excess Gibbs energy (not needed, only for checking results)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !gE = 0.d0
    !ln_gamma = 0.D0
    !Do j = 1,nderivs
    !    Do i = 1, ncomp
    !        ln_gamma(j,i) = ln_gamma_C(j,i) + ln_gamma_R(j,i)
    !        gE(j) = gE(j) + x(i) * ln_gamma(j,i) * R_const * Temp
    !    end do
    !End do
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !Write new values in save variables
    !gl%ge%Temp_prev = Temp
    !gl%ge%GETDER_prev = GETDER
    !gl%ge%molfrac_prev = gl%molfractions
    !gl%ge%Eq_type_prev = gl%Eq_type
    !gl%ge%mixtype_prev = gl%mix_type
    !gl%ge%C_or_R_prev = C_or_R
    !gl%ge%gE_C_prev = gl%ge%gE_C
    !gl%ge%gE_R_prev = gl%ge%gE_R
    !gl%ge%ln_gamma_C_prev = gl%ge%ln_gamma_C
    !gl%ge%ln_gamma_R_prev = gl%ge%ln_gamma_R
    deallocate(dlnGamma_k_dT)
    deallocate(dlnGamma_ki_dT)
    deallocate(dPsi_nm_dT)
    deallocate(d2lnGamma_k_dT2)
    deallocate(d2lnGamma_ki_dT2)
    deallocate(d2Psi_nm_dT2)


    deallocate(lnGamma_ki)
    deallocate(Theta_mi)
    deallocate(X_mi)
    deallocate(Psi_nm)


    end subroutine
    !--------------------------------------------------------------------------------------




    !The following functions and subroutines should be moved to other FORTRAN files, if the concept works

    !--------------------------------------------------------------------------------------
    subroutine DEPFUNC_GE_BASED (gl,Temp, Dens, getder, depfuncderivs)
    !--------------------------------------------------------------------------------------
    !Excess model based departure function, see
    ! Jäger and ... (XXXX)
    !
    !INPUT:
    ! Temp      -       Temperature in K
    ! Dens      -       Density in mol / m³
    ! getder    -       specify the derivatives that need to be calculated
    !                       INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                       1. NORMALIZED RESIDUAL HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
    !                       2. 1ST DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del
    !                       3. 2ND DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del^2
    !                       4. 1ST DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU
    !                       5: 2ND DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU^2
    !                       6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO DEL AND TAU, MULTIPLIED BY TAU*DEL
    !                       7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO DEL, TAU, AND TAU, MULTIPLIED BY TAU*TAU*DEL
    !                       8: 3RD DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^3
    !                       9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !                       10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !                       11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    !                       12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    !                       13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    !                       14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    !                       15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    !
    !   NOTE THAT THE EXCESS BASED DEPARTURE FUNCTION IS NOT A FUNCTION OF DENSITY. THUS, ALL DERIVATIVES WITH RESPECT TO DELTA ARE 0
    !
    ! OUTPUT
    ! depfuncderivs     -   array containing the derivatives of the excess based departure function with respect to delta and tau
    !
    !
    ! Andreas Jäger, June 2016
    ! Erik Mickolait, December 2020
    !--------------------------------------------------------------------------------------


    !use module_nderivs
    !use module_general_eos_parameters
    !use module_fluid_parameters
    !use module_cubic
    !use module_HelmgE

    use module_all_types
    implicit none

    type(type_gl) :: gl


    double precision:: Temp, Dens
    integer, dimension(nderivs):: GETDER                         ! array specifier to indicate, which derivative is needed
    double precision, dimension(nderivs)::depfuncderivs          ! array with the computed values for the derivatives

    double precision, dimension(30):: x

    double precision:: aE                           !Excess Helmholtz energy from the gE-model
    double precision, dimension(30) :: ln_gammai    !Activity coefficients from the gE-model

    double precision:: tredmix_orig, rhoredmix_orig

    double precision, dimension(30):: alpha_oi_r_mix, alpha_oi_r_i
    integer, dimension(nderivs):: GETDER_i                         ! array specifier to indicate, which derivative is needed for the pure fluid residual Helmholtz energies
    double precision, dimension(nderivs)::FNRDER_i                 ! array with the computed values for the derivatives of the residual Helmholtz energy of the pure fluids

    double precision:: help_f                               !This term is inspired by the Helmholtz translation of cubic EOS, specifically of the SRK. For the combination of the SRK
    double precision:: dhelpf_ddel, d2helpf_ddel2, d3helpf_ddel3, d4helpf_ddel4
    double precision:: help_h
    double precision:: dhelph_dtau, d2helph_dtau2, d3helph_dtau3, d4helph_dtau4
    double precision:: daE_dtau, d2aE_dtau2
    double precision:: const_A1
    double precision, dimension(30):: dalpha_oi_r_mix_dtau, d2alpha_oi_r_mix_dtau2
    double precision, dimension(30):: dalpha_oi_r_i_dtaui, d2alpha_oi_r_i_dtaui2

    double precision:: Dep_func_T_part                              !Part of departure function that only depends on temperature and not on density

    double precision:: tau, tau2, tau3, tau4
    double precision:: delta, delta2, delta3, delta4

    integer:: i

    !Variables for UNIFAC
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval

    !Variables for numerical derivatives
    double precision:: tau_p, tau_m, delta_tau
    double precision:: Dep_func_T_part_taup, Dep_func_T_part_taum
    double precision:: Temp_p, Temp_m

    double precision:: rhomix_ref_copy
    double precision, dimension(30):: rho_i_ref_copy

    depfuncderivs = 0.D0

    tredmix_orig = gl%tredmix
    rhoredmix_orig = gl%rhoredmix

    !Calculate delta and tau, as well as powers of delta and tau
    delta = Dens / gl%rhoredmix
    delta2 = delta * delta
    delta3 = delta2 * delta
    delta4 = delta2 * delta2
    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2

    !write mole fractions to other variable
    x = gl%molfractions

    !!Calculate the density of the mixture with the SRK-mixture covolume and the constant packing fraction
    !rho_mix_ref = 1.D0 / u_pack / b_HelmgE
    !!Calculate the densities of the pure fluids with the SRK covolumes and the constant packing fraction
    !do i = 1, ncomp
    !    rho_i_ref(i) = 1.D0 / u_pack / bi_HelmgE(i)
    !end do

    GETDER_i = 0
    !GETDER_i(1) = 1 !Only get the reduced Helmholtz energy itself at the moment. All derivatives with respect to tau and delta are done numerically
    GETDER_i(1) = 1 !The Helmholtz energy itself is needed
    GETDER_i(4) = 1 !The first derivative with respect to tau is needed
    GETDER_i(5) = 1 !The second derivative with respect to tau is needed
    do i = 1, gl%ncomp
        !Get the reduced Helmholtz energies of the pure components at mixture conditions
        rhomix_ref_copy = gl%rho_mix_ref
        call FNRDERIVS(gl,Temp, rhomix_ref_copy, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
        alpha_oi_r_mix(i) = FNRDER_i(1)
        dalpha_oi_r_mix_dtau(i) = FNRDER_i(4) / tau
        d2alpha_oi_r_mix_dtau2(i) = FNRDER_i(5) / tau**2.d0
    end do
    do i = 1, gl%ncomp
        !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
        gl%tredmix = gl%tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
        gl%rhoredmix = gl%rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
        rho_i_ref_copy(i) = gl%rho_i_ref(i)
        call FNRDERIVS(gl,Temp, rho_i_ref_copy(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
        alpha_oi_r_i(i) = FNRDER_i(1)
        dalpha_oi_r_i_dtaui(i)  = FNRDER_i(4) / (gl%tc(i)/temp)
        d2alpha_oi_r_i_dtaui2(i)  = FNRDER_i(5) / (gl%tc(i)/temp)**2
    end do
    gl%tredmix = tredmix_orig
    gl%rhoredmix = rhoredmix_orig

    !Get gE (or aE) from the gE-model
    !IMPLEMENTED UNIFAC AS STANDARD GE MODEL, Andreas Jäger, June 2017
    !call gE_NRTL(Temp, aE, ln_gammai)
    C_or_R = 2  !Get residual part of UNIFAC
    GETDER_i = 0
    
    !Andreas, August 2020
    !-------------
    if ((GETDER(1) == 1) .or. (GETDER(2) == 1) .or. (GETDER(3) == 1) .or. (GETDER(4) == 1) .or. (GETDER(6) == 1) .or. (GETDER(8) == 1) .or. (GETDER(11) == 1)) then !Only calculate the derivative when it is needed
        GETDER_i(1) = 1
    end if
    if ((GETDER(4) == 1) .or. (GETDER(6) == 1) .or. (GETDER(10) == 1) .or. (GETDER(12) == 1) .or. (GETDER(5) == 1) .or. (GETDER(7) == 1) .or. (GETDER(13) == 1)) then !Only calculate the derivative when it is needed
        GETDER_i(4) = 1
    end if
    if ((GETDER(5) == 1) .or. (GETDER(7) == 1) .or. (GETDER(13) == 1)) then !Only calculate the derivative when it is needed
        GETDER_i(5) = 1
    end if
    !-------------
    !for speed test
    !GETDER_i(1) = 1
    !GETDER_i(4) = 1
    !GETDER_i(5) = 1
    
    if (gl%mix_type == 12) then
        !call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_i, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, ln_gamma_R, C_or_R, errval)
        call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_i, C_or_R, errval)
    elseif (gl%mix_type == 13) then
        !ERIK, April 2018
        call gE_COSMO_SAC_MIXDERIVS(gl,Temp, GETDER_i, C_or_R, errval)
    end if
    aE = gl%ge%gE_C(1) + gl%ge%gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)
    daE_dtau = gl%ge%gE_C(4) + gl%ge%gE_R(4)
    d2aE_dtau2 = gl%ge%gE_C(5) + gl%ge%gE_R(5)

    !The density dependent term of the departure function is inspired by the combination of the SRK with gE models.
    !For the SRK, the reduced Helmholtz energy alpha^r gets:
    !   alpha^r = -ln(1-b*rho) - a/RT ln(b*rho + 1) / b
    !The mixing rule for a for the PSRK is: a = b*gE/A1 + b * sum(xi * ai / bi) + RTb / A1 * sum(xi ln(b/bi)), where A1 = -0.64663 = -ln(1 + 1/u) with packing fraction u = 1.1
    !Thus, alpha^r can be split in these parts:
    !   alpha^r = -ln(1-b*rho) - ln(b*rho + 1)/A1 * (gE/RT + A1/RT * sum (xi ai / bi) + sum(xi ln(b/bi) )
    !Thus, the PSRK inspired term is f = - ln(b*rho + 1)/A1
    const_A1 = -dlog(1.D0 / gl%u_pack + 1.D0)
    if ((GETDER(1) == 1) .or. (GETDER(4) == 1) .or. (GETDER(5) == 1) .or. (GETDER(9) == 1) .or. (GETDER(15) == 1)) then !Only calculate the derivative when it is needed
        help_f = -dlog(gl%b_HelmgE * Dens + 1.D0) / const_A1
    end if
    !Calculate derivatives of the help function f(del,x) with respect to delta
    if ((GETDER(2) == 1) .or. (GETDER(6) == 1) .or. (GETDER(7) == 1) .or. (GETDER(14) == 1)) then !Only calculate the derivative when it is needed
        dhelpf_ddel = - 1.D0 / const_A1 * gl%b_HelmgE * gl%rhoredmix / (gl%b_HelmgE * gl%rhoredmix * delta + 1.D0)
    end if
    if ((GETDER(3) == 1) .or. (GETDER(10) == 1) .or. (GETDER(13) == 1)) then !Only calculate the derivative when it is needed
        d2helpf_ddel2 = 1.D0 / const_A1 * gl%b_HelmgE**2.d0 * gl%rhoredmix**2.d0 / (gl%b_HelmgE * gl%rhoredmix * delta + 1.D0)**2
    end if
    if ((GETDER(8) == 1) .or. (GETDER(12) == 1)) then !Only calculate the derivative when it is needed
        d3helpf_ddel3 = - 2.D0 / const_A1 * gl%b_HelmgE**3.d0 * gl%rhoredmix**3.d0 / (gl%b_HelmgE * gl%rhoredmix * delta + 1.D0)**3
    end if
    if (GETDER(11) == 1) then !Only calculate the derivative when it is needed
        d4helpf_ddel4  = 6.D0 / const_A1 * gl%b_HelmgE**4.d0 * gl%rhoredmix**4.d0 / (gl%b_HelmgE * gl%rhoredmix * delta + 1.D0)**4
    end if


    !Calculate the auxiliary term h(tau,x) and its tau derivatives
    if ((GETDER(1) == 1) .or. (GETDER(2) == 1) .or. (GETDER(3) == 1) .or. (GETDER(8) == 1) .or. (GETDER(11) == 1)) then !Only calculate the derivative when it is needed
        help_h = aE / R_HelmgE / Temp
        do i = 1, gl%ncomp
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !help_h = help_h - x(i) * dlog(bi_HelmgE(i) / b_HelmgE) - x(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
            help_h = help_h - x(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
        end do
    end if
    if ((GETDER(4) == 1) .or. (GETDER(6) == 1) .or. (GETDER(10) == 1) .or. (GETDER(12) == 1)) then !Only calculate the derivative when it is needed
        dhelph_dtau = (tau * daE_dtau + aE) / R_HelmgE / gl%tredmix
        do i = 1, gl%ncomp
            dhelph_dtau = dhelph_dtau - x(i) * (dalpha_oi_r_mix_dtau(i) - dalpha_oi_r_i_dtaui(i) * gl%tc(i) / gl%tredmix)
        end do
    end if
    if ((GETDER(5) == 1) .or. (GETDER(7) == 1) .or. (GETDER(13) == 1)) then !Only calculate the derivative when it is needed
        d2helph_dtau2 = (2.D0 * daE_dtau + tau * d2aE_dtau2) / R_HelmgE / gl%tredmix
        do i = 1, gl%ncomp
            d2helph_dtau2 = d2helph_dtau2 - x(i) * (d2alpha_oi_r_mix_dtau2(i) - d2alpha_oi_r_i_dtaui2(i) * (gl%tc(i) / gl%tredmix)**2)
        end do
    end if
    if ((GETDER(9) == 1) .or. (GETDER(14) == 1)) then !Only calculate the derivative when it is needed
        d3helph_dtau3 = 0.D0    !NOT YET IMPLEMENTED
    end if
    if (GETDER(15) == 1) then !Only calculate the derivative when it is needed
        d4helph_dtau4 = 0.D0    !NOT YET IMPLEMENTED
    end if


    !!FOR NUMERICAL DERIVATIVES
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!Calculate the temperature dependent part separately, because it is needed in different derivatives
    !Dep_func_T_part = aE / R_HelmgE / Temp
    !do i = 1, ncomp
    !    Dep_func_T_part = Dep_func_T_part - x(i) * dlog(bi_HelmgE(i) / b_HelmgE) - x(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
    !end do
    !
    !!For all numerical derivatives with respect to tau, the temperature dependent part of the departure function is evaluated at higher and lower temperature
    !!-----------------------------------------------------------------------------------------------------------------------------------------------------------
    !    delta_tau = 1.D-4
    !
    !    !Increased tau
    !    !----------------------------------------------------------
    !    tau_p = tau + delta_tau
    !    Temp_p = tredmix/tau_p
    !
    !    GETDER_i = 0
    !    GETDER_i(1) = 1
    !    do i = 1, ncomp
    !        !Get the reduced Helmholtz energies of the pure components at mixture conditions
    !        call FNRDERIVS(Temp_p, rho_mix_ref, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
    !        alpha_oi_r_mix(i) = FNRDER_i(1)
    !    end do
    !    do i = 1, ncomp
    !        !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
    !        tredmix = tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
    !        rhoredmix = rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
    !        call FNRDERIVS(Temp_p, rho_i_ref(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
    !        alpha_oi_r_i(i) = FNRDER_i(1)
    !    end do
    !    tredmix = tredmix_orig
    !    rhoredmix = rhoredmix_orig
    !
    !    !Get gE (or aE) from the gE-model
    !    !IMPLEMENTED UNIFAC AS STANDARD GE MODEL, Andreas Jäger, June 2017
    !    !call gE_NRTL(Temp_p, aE, ln_gammai)
    !    call gE_UNIFAC_MIXDERIVS(Temp_p, GETDER, gE_C, gE_R, ln_gamma_C, ln_gamma_R, C_or_R, errval)
    !    aE = gE_C(1) + gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)
    !
    !    !Calculate the temperature dependent part separately, because it is needed in different derivatives
    !    Dep_func_T_part_taup = aE / R_HelmgE / Temp_p
    !    do i = 1, ncomp
    !        Dep_func_T_part_taup = Dep_func_T_part_taup - x(i) * dlog(bi_HelmgE(i) / b_HelmgE) - x(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
    !    end do
    !    !----------------------------------------------------------
    !
    !    !Decreased tau
    !    !----------------------------------------------------------
    !    tau_m = tau - delta_tau
    !    Temp_m = tredmix/tau_m
    !
    !    GETDER_i = 0
    !    GETDER_i(1) = 1
    !    do i = 1, ncomp
    !        !Get the reduced Helmholtz energies of the pure components at mixture conditions
    !        call FNRDERIVS(Temp_m, rho_mix_ref, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
    !        alpha_oi_r_mix(i) = FNRDER_i(1)
    !    end do
    !    do i = 1, ncomp
    !        !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
    !        tredmix = tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
    !        rhoredmix = rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
    !        call FNRDERIVS(Temp_m, rho_i_ref(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
    !        alpha_oi_r_i(i) = FNRDER_i(1)
    !    end do
    !    tredmix = tredmix_orig
    !    rhoredmix = rhoredmix_orig
    !
    !    !Get gE (or aE) from the gE-model
    !    !IMPLEMENTED UNIFAC AS STANDARD GE MODEL, Andreas Jäger, June 2017
    !    !call gE_NRTL(Temp_m, aE, ln_gammai)
    !    call gE_UNIFAC_MIXDERIVS(Temp_m, GETDER, gE_C, gE_R, ln_gamma_C, ln_gamma_R, C_or_R, errval)
    !    aE = gE_C(1) + gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)
    !
    !    !Calculate the temperature dependent part separately, because it is needed in different derivatives
    !    Dep_func_T_part_taum = aE / R_HelmgE / Temp_m
    !    do i = 1, ncomp
    !        Dep_func_T_part_taum = Dep_func_T_part_taum - x(i) * dlog(bi_HelmgE(i) / b_HelmgE) - x(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
    !    end do
    !    !----------------------------------------------------------
    !
    !!-----------------------------------------------------------------------------------------------------------------------------------------------------------
    !!END OF NUMERICAL DERIVATIVE WITH RESPECT TO TAU
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !Get the excess based departure function (reduced Helmholtz energy)
    if (GETDER(1) == 1) then
        !depfuncderivs(1) = Dep_func_T_part * help_f     !NUMERICAL
        depfuncderivs(1) = help_f * help_h              !ANALYTICAL
    end if

    !Get the first derivative of the excess based departure function with respect to delta, multiplied with delta
    if (GETDER(2) == 1) then
        !depfuncderivs(2) = Dep_func_T_part * dhelpf_ddel * delta    !NUMERICAL
        depfuncderivs(2) = dhelpf_ddel * help_h * delta             !ANALYTICAL
    end if

    !Get the second derivative of the excess based departure function with respect to delta, multiplied with delta^2
    if (GETDER(3) == 1) then
        !depfuncderivs(3) = Dep_func_T_part * d2helpf_ddel2 * delta2     !NUMERICAL
        depfuncderivs(3) = d2helpf_ddel2 * help_h * delta2                !ANALYTICAL
    end if

    !Get the first derivative of the excess based departure function with respect to tau, multiplied with tau
    if (GETDER(4) == 1) then
        !Numerical derivative with respect to tau
        !depfuncderivs(4) = (Dep_func_T_part_taup - Dep_func_T_part_taum) / 2.D0 / delta_tau * help_f * tau  !NUMERICAL
        depfuncderivs(4) = help_f * dhelph_dtau * tau                !ANALYTICAL
    end if

    !Get the second derivative of the excess based departure function with respect to tau, multiplied with tau^2
    if (GETDER(5) == 1) then
        !Second numerical derivative with respect to tau^2
        !depfuncderivs(5) = (Dep_func_T_part_taup - 2.D0 * Dep_func_T_part + Dep_func_T_part_taum) / delta_tau**2.d0 * help_f * tau2  !NUMERICAL
        depfuncderivs(5) = help_f * d2helph_dtau2 * tau2                !ANALYTICAL
    end if

    !Get the second derivative of the excess based departure function with respect to tau and delta, multiplied with tau*del
    if (GETDER(6) == 1) then
        !Numerical tau and analytical delta derivative
        !depfuncderivs(6) = (Dep_func_T_part_taup - Dep_func_T_part_taum) / 2.D0 / delta_tau * dhelpf_ddel * tau *delta  !NUMERICAL
        depfuncderivs(6) = dhelpf_ddel * dhelph_dtau * tau * delta               !ANALYTICAL
    end if

    !Get the third derivative of the excess based departure function with respect to tau, tau and delta, multiplied with tau^2*del
    if (GETDER(7) == 1) then
        !Numerical tau and analytical delta derivative
        !depfuncderivs(7) = (Dep_func_T_part_taup - 2.D0 * Dep_func_T_part + Dep_func_T_part_taum) / delta_tau**2.d0 * dhelpf_ddel * tau2 * delta     !NUMERICAL
        depfuncderivs(7) = dhelpf_ddel * d2helph_dtau2 * tau2 * delta               !ANALYTICAL
    end if

    !Get the third derivative of the excess based departure function with respect to delta, multiplied with delta^3
    if (GETDER(8) == 1) then
        !depfuncderivs(8) = Dep_func_T_part * d3helpf_ddel3 * delta3     !NUMERICAL
        depfuncderivs(8) = d3helpf_ddel3 * help_h * delta3               !ANALYTICAL
    end if

    !Get the third derivative of the excess based departure function with respect to tau, multiplied with tau^3
    if (GETDER(9) == 1) then
        !NOT YET IMPLEMENTED
        !depfuncderivs(9) = 0.D0
        depfuncderivs(9) = help_f * d3helph_dtau3 * tau3               !ANALYTICAL
    end if

    !Get the third derivative of the excess based departure function with respect to tau, delta and delta, multiplied with tau*delta^2
    if (GETDER(10) == 1) then
        !depfuncderivs(10) = d2helpf_ddel2 * (Dep_func_T_part_taup - Dep_func_T_part_taum) / 2.D0 * delta2 * tau   !NUMERICAL
        depfuncderivs(10) = d2helpf_ddel2 * dhelph_dtau * delta2 * tau               !ANALYTICAL
    end if

    !Get the fourth derivative of the excess based departure function with respect to delta, multiplied with delta^4
    if (GETDER(11) == 1) then
        !depfuncderivs(11) = Dep_func_T_part * d4helpf_ddel4 * delta4    !NUMERICAL
        depfuncderivs(11) = d4helpf_ddel4 * help_h * delta4    !ANALYTICAL
    end if

    !Get the fourth derivative of the excess based departure function with respect to delta, delta, delta and tau, multiplied with del^3*tau
    if (GETDER(12) == 1) then
        !depfuncderivs(12) = d3helpf_ddel3 * (Dep_func_T_part_taup - Dep_func_T_part_taum) / 2.D0 * delta3 * tau   !NUMERICAL
        depfuncderivs(12) = d3helpf_ddel3 * dhelph_dtau * delta3 * tau               !ANALYTICAL
    end if

    !Get the fourth derivative of the excess based departure function with respect to delta, delta, tau and tau, multiplied with del^2*tau^2
    if (GETDER(13) == 1) then
        !NOT YET IMPLEMENTED
        !depfuncderivs(13) = d2helpf_ddel2 * (Dep_func_T_part_taup - 2.D0 * Dep_func_T_part + Dep_func_T_part_taum) / delta_tau**2.d0 * delta2 * tau2   !NUMERICAL
        depfuncderivs(13) = d2helpf_ddel2 * d2helph_dtau2 * delta2 * tau2               !ANALYTICAL
    end if

    !Get the fourth derivative of the excess based departure function with respect to tau, tau, tau and delta, multiplied with tau^3*del
    if (GETDER(14) == 1) then
        !NOT YET IMPLEMENTED
        !depfuncderivs(14) = 0.D0
        depfuncderivs(14) = dhelpf_ddel * d3helph_dtau3 * delta * tau3               !ANALYTICAL
    end if

    !Get the fourth derivative of the excess based departure function with respect to tau, multiplied with tau^4
    if (GETDER(15) == 1) then
        !NOT YET IMPLEMENTED
        !depfuncderivs(15) = 0.D0
        depfuncderivs(15) = help_f * d4helph_dtau4 * tau4              !ANALYTICAL
    end if


    end subroutine
    !--------------------------------------------------------------------------------------



    !--------------------------------------------------------------------------------------
    subroutine DEPFUNC_NON_COR_STATE (gl,Temp, Dens, getder, depfuncder, FLD1, FLD2)
    !--------------------------------------------------------------------------------------
    !Departure function for mixture model that is not any more corresponding states based
    !(mix_type = 19)
    ! Jäger and ... (XXXX)
    !
    !INPUT:
    ! Temp      -       Temperature in K
    ! Dens      -       Density in mol / m³
    ! getder    -       specify the derivatives that need to be calculated
    !                       INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                       1. NORMALIZED RESIDUAL HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
    !                       2. 1ST DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del
    !                       3. 2ND DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del^2
    !                       4. 1ST DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU
    !                       5: 2ND DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU^2
    !                       6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO DEL AND TAU, MULTIPLIED BY TAU*DEL
    !                       7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO DEL, TAU, AND TAU, MULTIPLIED BY TAU*TAU*DEL
    !                       8: 3RD DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^3
    !                       9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !                       10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !                       11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    !                       12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    !                       13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    !                       14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    !                       15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    !
    !
    ! OUTPUT
    ! depfuncder     -   array containing the derivatives of the excess based departure function with respect to delta and tau
    !
    !
    ! Andreas Jäger, October 2016
    !--------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: Temp, Dens
    integer, dimension(nderivs):: GETDER                         ! array specifier to indicate, which derivative is needed
    double precision, dimension(nderivs)::depfuncder          ! array with the computed values for the derivatives
    integer, intent(in):: FLD1, FLD2

    double precision, dimension(100):: reg_termi, reg_termj
    double precision:: deltai, deltaj, taui, tauj

    integer:: k, j, i

    !The following function can be regarded as the "departure function" of the new one-fluid model
    DEPFUNCDER = 0.D0
    !sum over the binary functions

    i = FLD1
    j = FLD2
    !Binary departure function for this model is symmetric. Only implemented case is i<j, so in case i>j this is resorted here
    if (i .gt. j) then
        k = i
        i = j
        j = k
    end if

    deltai = dens / gl%rhoc(i)
    taui = gl%tc(i) / Temp

    deltaj = dens / gl%rhoc(j)
    tauj = gl%tc(j) / Temp



    !!Summation over the regular terms
    do k = 1, gl%mix_nreg(i,j)
        reg_termi(k) = gl%mix_nij(i,j,k,1) * taui**gl%mix_tij(i,j,k) * deltai**gl%mix_dij(i,j,k) * &
            &  dexp(-gl%mix_gama(i,j,k) * deltai**gl%mix_p_ij(i,j,k))   * 0.5D0 !* (1.D0 - gl%Helm_k_nij(i,j,k))
        reg_termj(k) = gl%mix_nij(i,j,k,2) * tauj**gl%mix_tij(i,j,k) * deltaj**gl%mix_dij(i,j,k) * &
            &  dexp(-gl%mix_gama(i,j,k) * deltaj**gl%mix_p_ij(i,j,k))   * 0.5D0 !* (1.D0 - gl%Helm_k_nij(i,j,k))
        DEPFUNCDER(1) = DEPFUNCDER(1) + (reg_termi(k) + reg_termj(k))
    end do

    !The departure function of the one-fluid model itself
    if (GETDER(1) .eq. 1) then
        DEPFUNCDER(1) = DEPFUNCDER(1)
    end if

    !The derivative of the departure function with respect to delta, multiplied with delta
    if (GETDER(2) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(2)= DEPFUNCDER(2) &
                & + (reg_termi(k) * (gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltai**gl%mix_p_ij(i,j,k)) &
                & + reg_termj(k) * (gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltaj**gl%mix_p_ij(i,j,k)))
        end do
    end if

    !The second derivative of the departure function with respect to delta^2, multiplied with delta^2
    if (GETDER(3) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(3)= DEPFUNCDER(3) &
                & + (reg_termi(k) * ((gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltai**gl%mix_p_ij(i,j,k)) * &
                &   (gl%mix_dij(i,j,k) - 1.D0 - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltai**gl%mix_p_ij(i,j,k)) - &
                &    gl%mix_gama(i,j,k) * (gl%mix_p_ij(i,j,k)**2) * deltai**gl%mix_p_ij(i,j,k)  ) &
                & + reg_termj(k) * ((gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltaj**gl%mix_p_ij(i,j,k)) * &
                &   (gl%mix_dij(i,j,k) - 1.D0 - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltaj**gl%mix_p_ij(i,j,k)) - &
                &    gl%mix_gama(i,j,k) * (gl%mix_p_ij(i,j,k)**2) * deltaj**gl%mix_p_ij(i,j,k)  ) )
        end do
    end if

    !The derivative of the departure function with respect to tau, multiplied with tau
    if (GETDER(4) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(4)= DEPFUNCDER(4) + ((reg_termi(k) + reg_termj(k)) * gl%mix_tij(i,j,k))
        end do
    end if

    !The derivative of the departure function with respect to tau^2, multiplied with tau^2
    if (GETDER(5) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(5)= DEPFUNCDER(5) + ((reg_termi(k) + reg_termj(k)) * gl%mix_tij(i,j,k) * (gl%mix_tij(i,j,k) - 1.D0))
        end do
    end if

    !The derivative of the departure function with respect to tau and del, multiplied with tau*del
    if (GETDER(6) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(6)= DEPFUNCDER(6) + &
                & ((reg_termi(k) * (gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltai**gl%mix_p_ij(i,j,k)) + &
                & reg_termj(k) * (gl%mix_dij(i,j,k) - gl%mix_gama(i,j,k) * gl%mix_p_ij(i,j,k) * deltaj**gl%mix_p_ij(i,j,k))) * gl%mix_tij(i,j,k))
        end do
    end if

    !The derivative of the departure function with respect to tau^2 and del, multiplied with tau^2*del
    if (GETDER(7) .eq. 1) then
        DEPFUNCDER(7) = 0.D0
    end if

    !The derivative of the departure function with respect to del^3, multiplied with del^3
    if (GETDER(8) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(8)= DEPFUNCDER(8) &
                & + (reg_termi(k) * (gl%mix_dij(i,j,k)*(gl%mix_dij(i,j,k) - 1.d0)*(gl%mix_dij(i,j,k) - 2.d0) + &
                &   gl%mix_gama(i,j,k)*gl%mix_p_ij(i,j,k)*deltai**gl%mix_p_ij(i,j,k)*(-2.d0 + 6.d0*gl%mix_dij(i,j,k) &
                &   - 3.d0*gl%mix_dij(i,j,k)**2 - 3.d0*gl%mix_dij(i,j,k)*gl%mix_p_ij(i,j,k) + 3.d0*gl%mix_p_ij(i,j,k) - gl%mix_p_ij(i,j,k)**2) &
                &   + 3.d0*gl%mix_gama(i,j,k)**2*gl%mix_p_ij(i,j,k)**2*deltai**(2.d0*gl%mix_p_ij(i,j,k))*(gl%mix_dij(i,j,k) - 1.d0 + gl%mix_p_ij(i,j,k))&
                &   - gl%mix_gama(i,j,k)**3*gl%mix_p_ij(i,j,k)**3*deltai**(3*gl%mix_p_ij(i,j,k)) ) &
                & + reg_termj(k) * (gl%mix_dij(i,j,k)*(gl%mix_dij(i,j,k) - 1.d0)*(gl%mix_dij(i,j,k) - 2.d0) + &
                &   gl%mix_gama(i,j,k)*gl%mix_p_ij(i,j,k)*deltaj**gl%mix_p_ij(i,j,k)*(-2.d0 + 6.d0*gl%mix_dij(i,j,k) &
                &   - 3.d0*gl%mix_dij(i,j,k)**2 - 3.d0*gl%mix_dij(i,j,k)*gl%mix_p_ij(i,j,k) + 3.d0*gl%mix_p_ij(i,j,k) - gl%mix_p_ij(i,j,k)**2) &
                &   + 3.d0*gl%mix_gama(i,j,k)**2*gl%mix_p_ij(i,j,k)**2*deltaj**(2.d0*gl%mix_p_ij(i,j,k))*(gl%mix_dij(i,j,k) - 1.d0 + gl%mix_p_ij(i,j,k))&
                &   - gl%mix_gama(i,j,k)**3*gl%mix_p_ij(i,j,k)**3*deltaj**(3*gl%mix_p_ij(i,j,k)) ) )
        end do
    end if

    !The derivative of the departure function with respect to tau^3, multiplied with tau^3
    if (GETDER(9) .eq. 1) then
        !!Summation over the regular terms
        do k = 1, gl%mix_nreg(i,j)
            DEPFUNCDER(9)= DEPFUNCDER(9) + ((reg_termi(k) + reg_termj(k)) * &
                & gl%mix_tij(i,j,k) * (gl%mix_tij(i,j,k) - 1.D0) * (gl%mix_tij(i,j,k) - 2.D0))
        end do
    end if

    !The derivative of the departure function with respect to tau*del*del, multiplied with tau*del*del
    if (GETDER(10) .eq. 1) then
        DEPFUNCDER(10) = 0.D0
    end if


    end subroutine
    !--------------------------------------------------------------------------------------

    !*****************************************************************************************************
    !END OF PRELIMINARY CODE
    !*****************************************************************************************************

    !*****************************************************************************************************
    !START OF COSMO-SAC CODE
    !*****************************************************************************************************


    !ERIK, April 2018
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine gE_COSMO_SAC_MIXDERIVS(gl, Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the COSMO-SAC model as well as tau- and delta-
    !derivatives of the model.
    !
    !The model is described by
    !....
    !....
    !
    !
    !INPUT
    !Temp       -    Temperature in K
    !Note that the model makes use of the global variable "molfractions" defined in the module
    !"module_fluid_parameters"
    ! GETDER      - AN ARRAY WITH nderivs ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. gE from UNIFAC as a function of temperature and composition
    !                2. 1ST DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                3. 2ND DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                4. 1ST DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                5: 2ND DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                6: 2ND MIXED DERIVATIVE OF gE WITH RESPECT TO DEL AND TAU
    !                7: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO DEL, TAU, AND TAU
    !                8: 3RD DERIVATIVE OF gE WITH RESPECT TO DEL
    !                9: 3RD DERIVATIVE OF gE WITH RESPECT TO TAU
    !               10: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, AND DEL
    !               11: 4TH DERIVATIVE OF gE WITH RESPECT TO DEL
    !               12: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, DEL AND DEL
    !               13: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, DEL, AND DEL
    !               14: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, TAU AND DEL
    !               15: 4TH DERIVATIVE OF gE WITH RESPECT TO TAU
    !C_or_R       - Integer that specifies if the combinatorial, the residual or both parts shall be calculated
    !               C_or_R = 0: Calculate both parts
    !               C_or_R = 1: Calculate the combinatorial part only
    !               C_or_R = 2: Calculate the the residual part only
    !
    !OUTPUT
    !gE_C         -    Vector with molar excess Gibbs energy of the combinatorial part (C) in J/molK and derivatives with respect to tau and delta
    !gE_R         -    Vector with molar excess Gibbs energy of the residual part (R) in J/molK and derivatives with respect to tau and delta
    !ln_gamma_C   -    Matrix containing the natural logarithm of the activity coefficients of the combinatorial part for each component and derivatives with respect to tau and delta
    !ln_gamma_R   -    Matrix containing the natural logarithm of the activity coefficients of the residual part for each component and derivatives with respect to tau and delta
    !errval       -    Indicates if an error occurred during calculations
    !
    !Andreas Jäger, January 2017
    !------------------------------------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl
    !type(type_COSMO_SAC) :: COSMO_SAC

    double precision :: Temp
    integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs) :: gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval

    !Test variables for overall excess Gibbs energy and activity coefficients
    double precision, dimension(nderivs) :: gE          !Excess Gibbs energy
    double precision, dimension(nderivs,30) :: ln_gamma !activity coefficients

    double precision, dimension(30) :: r_i              !Molecular vdW volume of components i
    double precision, dimension(30) :: q_i              !Molecular vdW surface area of components i

    double precision :: R_const, R_const_cal            !universal gas constant
    double precision, dimension(30) :: x
    double precision:: tau, del

    !Derivatives of the combinatorial and residual part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho

    !New variables for COSMO-SAC
    double precision :: delta_tau, Temp_p, Temp_m
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_C_act, gE_R_p, gE_R_m, gE_R_act
    double precision, dimension(nderivs, 30) :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_C_act, ln_gamma_R_p, ln_gamma_R_m, ln_gamma_R_act
    integer :: i, j, k, m, n, o

    !Analytical derivation of COSMO-SAC based gE
    double precision, dimension(51,30) :: dseggammadt_pure, dseggammadt_pure_gauss, dseggammadt_pure_test
    double precision, dimension(51) :: dseggammadt_mix, dseggammadt_mix_gauss
    double precision, dimension (51,30) :: Fi_T, Fi_T_pure
    double precision :: Det_seggamma, Det_T
    !double precision, dimension(60,60):: MatrixA, MatrixB, dMatrixAdT, dMatrixBdT, MatrixC, MatrixD, MatrixC_mix!, MatrixA_p, MatrixA_m, MatrixB_p, MatrixB_m,
    double precision, dimension(:,:), allocatable:: MatrixA
    double precision, dimension(:,:), allocatable:: MatrixB
    double precision, dimension(:,:), allocatable:: MatrixC
    double precision, dimension(:,:), allocatable:: MatrixC_mix
    double precision, dimension(:,:), allocatable:: MatrixD
    double precision, dimension(:,:), allocatable:: dMatrixBdT
    double precision, dimension(:,:), allocatable:: dMatrixAdT
    !November 2018, erik, allocatable
    !double precision, dimension(51,60,60):: MatrixA_mix, MatrixA_mix_p, MatrixA_mix_m
    double precision, dimension(:,:,:), allocatable:: MatrixA_mix
    !double precision, dimension(60,60) ::MatrixB_mix!, MatrixB_mix_p, MatrixB_mix_m
    double precision, dimension(:,:), allocatable:: MatrixB_mix
    !November 2018, erik, allocatable
    !double precision, dimension(51,30,60,60) :: MatrixA_pure, MatrixA_pure_p, MatrixA_pure_m
    double precision, dimension(:,:,:,:), allocatable :: MatrixA_pure
    !double precision, dimension(:,:,:,:), allocatable :: MatrixA_pure_p
    !double precision, dimension(:,:,:,:), allocatable :: MatrixA_pure_m
    !double precision, dimension(60,60,30) :: MatrixB_pure!, MatrixB_pure_p, MatrixB_pure_m
    double precision, dimension(:,:,:), allocatable :: MatrixB_pure
    integer:: rankA
    double precision, dimension(60):: vectorb, vectorb_gauss
    double precision, dimension(60):: vectorx
    integer:: errorflag
    double precision, dimension(15,30) :: ln_gamma_R_test
    double precision, dimension(51) :: dlnseggammadt_mix, dlnseggammadtau_mix
    double precision, dimension(51,30) :: dlnseggammadt_pure, dlnseggammadtau_pure
    double precision, dimension(51) :: sum_d2Fi_dT2, sum2_d2Fi_dT2, sum_d2Fi_dseggammadT
    double precision, dimension(51,30) :: d2Fi_dT2_pure, d2Fi_dT2_pure_test, sum3_d2Fi_dT2
    double precision, dimension(51) :: d2Fi_dT2_mix
    double precision, dimension(51,30) :: DetMatrixA_pure! , DetMatrixA_pure_p, DetMatrixA_pure_m, DetMatrixA_pure_test
    double precision, dimension(30) :: DetMatrixB_pure!, DetMatrixB_pure_p, DetMatrixB_pure_m
    double precision, dimension(51) :: DetMatrixA_mix!, DetMatrixA_mix_p, DetMatrixA_mix_m
    double precision :: DetMatrixB_mix!, DetMatrixB_mix_p, DetMatrixB_mix_m

    !double precision, dimension(60,60):: adj_A, adj_B
    double precision, dimension(:,:), allocatable:: adj_A
    double precision, dimension(:,:), allocatable:: adj_B
    double precision, dimension(51,30) :: dDetAdT_pure, dDetAdT_pure_num, dDetBdT_pure_num
    double precision, dimension(30) :: dDetBdT_pure
    double precision, dimension(51) :: dDetAdT_mix, dDetAdT_mix_num, dDetBdT_mix_num
    double precision :: dDetBdT_mix
    double precision, dimension(51,30) :: d2seggammadT2_pure, d2lnseggammadT2_pure
    double precision, dimension(51) :: d2seggammadT2_mix, d2lnseggammadT2_mix
    double precision, dimension(51,30) :: d2lnseggammadtau2_pure
    double precision, dimension(51) :: d2lnseggammadtau2_mix

    double precision, dimension(51,30) :: seggamma_pure_p, seggamma_pure_m, d2FdT2_pure!, Fi_T_pure_p, Fi_T_pure_m
    !double precision, dimension(51,51,30) :: d2Fi_dseggammadT_pure!,Fi_seggamma_pure_p, Fi_seggamma_pure_m, d2FdTdseggamma_pure
    double precision, dimension(:,:,:), allocatable :: d2Fi_dseggammadT_pure!,Fi_seggamma_pure_p, Fi_seggamma_pure_m, d2FdTdseggamma_pure
    double precision, dimension(51) :: seggamma_p, seggamma_m, d2FdT2_mix!, Fi_T_p, Fi_T_m
    !double precision, dimension(51,51) :: d2Fi_dseggammadT_mix!, d2FdTdseggamma_mix!, Fi_seggamma_mix_p, Fi_seggamma_mix_m
    double precision, dimension(:,:), allocatable :: d2Fi_dseggammadT_mix
    !integer, dimension(3) :: start_sigma, end_sigma
    integer, dimension(51) :: pos
    double precision, dimension(51) :: sum2_v1


    !for test
    !double precision, dimension(60,60,51) :: MatrixC_mix
    double precision, dimension(15) :: gE_R_Trho_p, gE_R_Trho_m, gE_R_trho_num, gE_R_num

    !November 2018, erik, allocatable
    if(.not. allocated(MatrixA_pure)) allocate(MatrixA_pure(51,30,60,60))
    if(.not. allocated(MatrixB_pure)) allocate(MatrixB_pure(60,60,30))
    !if(.not. allocated(MatrixA_pure_p)) allocate(MatrixA_pure_p(51,30,60,60))
    !if(.not. allocated(MatrixA_pure_m)) allocate(MatrixA_pure_m(51,30,60,60))
    if(.not. allocated(MatrixA_mix)) allocate(MatrixA_mix(51,60,60))
    if(.not. allocated(dMatrixBdT)) allocate(dMatrixBdT(60,60))
    if(.not. allocated(dMatrixAdT)) allocate(dMatrixAdT(60,60))
    if(.not. allocated(MatrixA)) allocate(MatrixA(60,60))
    if(.not. allocated(MatrixB)) allocate(MatrixB(60,60))
    if(.not. allocated(MatrixB_mix)) allocate(MatrixB_mix(60,60))
    if(.not. allocated(MatrixC)) allocate(MatrixC(60,60))
    if(.not. allocated(MatrixC_mix)) allocate(MatrixC_mix(60,60))
    if(.not. allocated(MatrixD)) allocate(MatrixD(60,60))
    if(.not. allocated(adj_A)) allocate(adj_A(60,60))
    if(.not. allocated(adj_B)) allocate(adj_B(60,60))
    if(.not. allocated(d2Fi_dseggammadT_pure)) allocate(d2Fi_dseggammadT_pure(51,51,30))
    if(.not. allocated(d2Fi_dseggammadT_mix)) allocate(d2Fi_dseggammadT_mix(51,51))



    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C = errval
        gl%ge%gE_R = errval
        gl%ge%ln_gamma_C = errval
        gl%ge%ln_gamma_R = errval
        return
    end if

    !Initialize
    gl%ge%gE_C = 0.D0
    gl%ge%gE_R = 0.D0
    gl%ge%ln_gamma_C = 0.D0
    gl%ge%ln_gamma_R = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.D0
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0
    R_const_cal = 8.3144598D0 / 4184.D0

    !Get mole fractions of the mixture from module variable
    x = gl%molfractions

    !Calculate tau
    tau = gl%tredmix / Temp

    !Set dummy value for delta (Not needed because COSMO_SAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0
    if ((GETDER(1) .eq. 1) .or. (GETDER(5) .eq. 1)) then
        call gE_COSMO_SAC_CALC(gl, Temp, C_or_R, errval)
        gE_C_act = gl%ge%gE_C
        gE_R_act = gl%ge%gE_R
        ln_gamma_C_act = gl%ge%ln_gamma_C
        ln_gamma_R_act = gl%ge%ln_gamma_R
    end if

    if ((gl%cosmo%analytical == .false.) .or. (gl%cosmo%COSMO_ver /= 1)) then

        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then
            delta_tau = 1.d-4
            Temp_p = gl%tredmix / (tau + delta_tau)
            Temp_m = gl%tredmix / (tau - delta_tau)
            call gE_COSMO_SAC_CALC(gl, Temp_p, C_or_R, errval)
            gE_C_p = gl%ge%gE_C
            gE_R_p = gl%ge%gE_R
            ln_gamma_C_p = gl%ge%ln_gamma_C
            ln_gamma_R_p = gl%ge%ln_gamma_R
            call gE_COSMO_SAC_CALC(gl, Temp_m, C_or_R, errval)
            gE_C_m = gl%ge%gE_C
            gE_R_m = gl%ge%gE_R
            ln_gamma_C_m = gl%ge%ln_gamma_C
            ln_gamma_R_m = gl%ge%ln_gamma_R
        end if
    end if

    !Erik, September 2018
    !for now combinatorial part not needed
    !If it will be needed in the future, the analytical derivations are the same as in UNIFAC, has to be implemented...
    !if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
    !    !Calculate the combinatorial part of the activity coefficients
    !    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !    !ln_gamma_C is always needed for the following calculations, thus it is always calculated
    !
    !    if (GETDER(1) .eq. 1) then
    !        !call gE_COSMO_SAC_CALC(Temp, gE_C, gE_R, ln_gamma_C, ln_gamma_R, C_or_R, errval)
    !        gl.gE_C(1) = gE_C_act(1)
    !        gl.ln_gamma_C(1,:) = ln_gamma_C_act(1,:)
    !        if (C_or_R == 1) then
    !            gl.gE_R(1) = 0.D0
    !            gl.ln_gamma_R(1,:) =0.D0
    !        end if
    !    end if
    !
    !    !All other derivatives of the combinatorial part are 0 (except for the tau-only derivatives of gE), because the combinatorial part is a function of composition only
    !    !Derivative wrt delta
    !    if (GETDER(2) .eq. 1) then
    !        ln_gamma_C_Trho(2,:) = 0.D0
    !        gE_C_Trho(2) = 0.d0
    !        gl.ln_gamma_C(2,:) = 0.d0
    !        gl.gE_C(2) = 0.D0
    !    end if
    !    !Second derivative wrt delta
    !    if (GETDER(3) .eq. 1) then
    !        ln_gamma_C_Trho(3,:) = 0.D0
    !        gE_C_Trho(3) = 0.d0
    !        gl.ln_gamma_C(3,:) = 0.d0
    !        gl.gE_C(3) = 0.D0
    !    end if
    !    !Derivative wrt tau
    !    if (GETDER(4) .eq. 1) then
    !        gl.gE_C(4) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_tau)
    !        !gl.gE_R(4) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_tau)
    !        do i = 1, gl.ncomp
    !            gl.ln_gamma_C(4,i) = (ln_gamma_C_p(1,i) - ln_gamma_C_m(1,i)) / (2.D0 * delta_tau)
    !            !gl.ln_gamma_R(4,i) = (ln_gamma_R_p(1,i) - ln_gamma_R_m(1,i)) / (2.D0 * delta_tau)
    !        end do
    !    end if
    !    !Second derivative wrt tau
    !    if (GETDER(5) .eq. 1) then
    !        gl.gE_C(5) = (gE_C_p(1) - 2.D0 * gE_C_act(1) + gE_C_m(1)) / (delta_tau ** 2)
    !        !gl.gE_R(5) = (gE_R_p(1) - 2.D0 * gl.gE_R(1) + gE_R_m(1)) / (delta_tau ** 2)
    !        do i = 1, gl.ncomp
    !            gl.ln_gamma_C(5,i) = (ln_gamma_C_p(1,i) - 2.D0 * ln_gamma_C_act(1,i) + ln_gamma_C_m(1,i)) / (delta_tau ** 2)
    !            !gl.ln_gamma_R(5,i) = (ln_gamma_R_p(1,i) - 2.D0 * gl.ln_gamma_R(1,i) + ln_gamma_R_m(1,i)) / (delta_tau ** 2)
    !        end do
    !    end if
    !    !Second derivative wrt delta and tau
    !    if (GETDER(6) .eq. 1) then
    !        ln_gamma_C_Trho(6,:) = 0.D0
    !        gE_C_Trho(6) = 0.d0
    !        gl.ln_gamma_C(6,:) = 0.d0
    !        gl.gE_C(6) = 0.D0
    !    end if
    !    !Third derivative wrt delta, tau, and tau
    !    if (GETDER(7) .eq. 1) then
    !        ln_gamma_C_Trho(7,:) = 0.D0
    !        gE_C_Trho(7) = 0.d0
    !        gl.ln_gamma_C(7,:) = 0.d0
    !        gl.gE_C(7) = 0.D0
    !    end if
    !    !Third derivative wrt delta
    !    if (GETDER(8) .eq. 1) then
    !        ln_gamma_C_Trho(8,:) = 0.D0
    !        gE_C_Trho(8) = 0.d0
    !        gl.ln_gamma_C(8,:) = 0.d0
    !        gl.gE_C(8) = 0.D0
    !    end if
    !    !Third derivative wrt tau
    !    if (GETDER(9) .eq. 1) then
    !        ln_gamma_C_Trho(9,:) = 0.D0
    !        gE_C_Trho(9) = 0.d0
    !        gl.ln_gamma_C(9,:) = 0.d0
    !        gl.gE_C(9) = 0.D0      !TO BE IMPLEMENTED
    !    end if
    !    !Third derivative wrt tau, delta, and delta
    !    if (GETDER(10) .eq. 1) then
    !        ln_gamma_C_Trho(10,:) = 0.D0
    !        gE_C_Trho(10) = 0.d0
    !        gl.ln_gamma_C(10,:) = 0.d0
    !        gl.gE_C(10) = 0.D0
    !    end if
    !    !Fourth derivative wrt delta
    !    if (GETDER(11) .eq. 1) then
    !        ln_gamma_C_Trho(11,:) = 0.D0
    !        gE_C_Trho(11) = 0.d0
    !        gl.ln_gamma_C(11,:) = 0.d0
    !        gl.gE_C(11) = 0.D0
    !    end if
    !    !Fourth derivative wrt delta, delta, delta and tau
    !    if (GETDER(12) .eq. 1) then
    !        ln_gamma_C_Trho(12,:) = 0.D0
    !        gE_C_Trho(12) = 0.d0
    !        gl.ln_gamma_C(12,:) = 0.d0
    !        gl.gE_C(12) = 0.D0
    !    end if
    !    !Fourth derivative wrt delta, delta, tau, tau
    !    if (GETDER(13) .eq. 1) then
    !        ln_gamma_C_Trho(13,:) = 0.D0
    !        gE_C_Trho(13) = 0.d0
    !        gl.ln_gamma_C(13,:) = 0.d0
    !        gl.gE_C(13) = 0.D0
    !    end if
    !    !Fourth derivative wrt delta, tau, tau, tau
    !    if (GETDER(14) .eq. 1) then
    !        ln_gamma_C_Trho(14,:) = 0.D0
    !        gE_C_Trho(14) = 0.d0
    !        gl.ln_gamma_C(14,:) = 0.d0
    !        gl.gE_C(14) = 0.D0
    !    end if
    !    !Fourth derivative wrt tau
    !    if (GETDER(15) .eq. 1) then
    !        ln_gamma_C_Trho(15,:) = 0.D0
    !        gE_C_Trho(15) = 0.d0
    !        gl.ln_gamma_C(15,:) = 0.d0
    !        gl.gE_C(15) = 0.D0     !TO BE IMPLEMENTED
    !    end if
    !    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !end if

    !TODO !!!
    !    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
    !Calculate the combinatorial part of the activity coefficients
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !ln_gamma_C is always needed for the following calculations, thus it is always calculated

    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !TODO !!!

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if (GETDER(1) .eq. 1) then
            !call gE_COSMO_SAC_CALC(Temp, gE_C, gE_R, ln_gamma_C, ln_gamma_R, C_or_R, errval)
            gl%ge%gE_R(1) = gE_R_act(1)
            gl%ge%ln_gamma_R(1,:) = ln_gamma_R_act(1,:)
            if (C_or_R == 2) then
                gl%ge%gE_C(1) = 0.D0
                gl%ge%ln_gamma_C(1,:) =0.D0
            end if
        end if

        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho(2,:) = 0.D0
            gE_R_Trho(2) = 0.d0
            gl%ge%ln_gamma_R(2,:) = 0.d0
            gl%ge%gE_R(2) = 0.D0
        end if

        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho(3,:) = 0.D0
            gE_R_Trho(3) = 0.d0
            gl%ge%ln_gamma_R(3,:) = 0.d0
            gl%ge%gE_R(3) = 0.D0
        end if

        !alternative analytical method to calculate Derivative wrt tau
        if (GETDER(4) .eq. 1) then
            if ((gl%cosmo%COSMO_ver == 1) .and. (gl%cosmo%analytical)) then

                gl%ge%ln_gamma_R(4,:) = 0.D0
                gl%ge%gE_R(4) = 0.D0
                gE_R_Trho(4) = 0.D0
                ln_gamma_R_Trho(4,:) = 0.D0

                !set all COSMO-SAC values to orig
                !call gE_COSMO_SAC_CALC(gl, Temp, C_or_R, errval)

                !calculate all derivatives of seggamma and seggamma_pure

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !pure fluid
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                !In case segment activity coefficients were not calculated using Newton-Raphson, dFi_pure_v1_gl needs to be calculated here
                if (maxval(gl.cosmo.dFi_pure_v1_gl) < 1.D-14) then
                    !calculate gl.cosmo.dFi_pure_v1_gl
                    do i = 1, gl.ncomp
                        sum2_v1 = 0.0
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                                sum2_v1(k) = sum2_v1(k) + (gl.cosmo.sigma(j,i) / gl.cosmo.ACOSMO(i)) * gl.cosmo.seggamma_pure_gl(j,i) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            end do
                            !derivatives of Fi after each seggamma results in two different equations
                            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                                if (k == j) then
                                    gl.cosmo.dFi_pure_v1_gl(k,j,i) = sum2_v1(k) + gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                                else
                                    gl.cosmo.dFi_pure_v1_gl(k,j,i) = gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.sigma(j,i) / gl.cosmo.ACOSMO(i)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                                end if
                            end do
                        end do
                    end do
                end if

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                Fi_T_pure = 0.D0
                do i = 1, gl.ncomp

                    !Fi_T_pure_p = 1.D-14
                    !Fi_T_pure_m = 1.D-14
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)    !all Fi_T_pure's
                        !if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)     !sum for Fi_T_pure
                            !partial derivatives in respect to Temperature
                            Fi_T_pure(j,i) = Fi_T_pure(j,i) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !!nummerical derivatives test
                            !Fi_T_pure_p(j,i) = Fi_T_pure_p(j,i) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * seggamma_pure_p(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_p ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_p))
                            !Fi_T_pure_m(j,i) = Fi_T_pure_m(j,i) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * seggamma_pure_m(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_m ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_m))
                            !d2FdTdseggamma_pure(j,k,i) = (Fi_seggamma_pure_p(j,k,i) - Fi_seggamma_pure_m(j,k,i)) / (Temp_p - Temp_m)
                        end do
                        Fi_T_pure(j,i) = gl.cosmo.seggamma_pure_gl(j,i) * Fi_T_pure(j,i)
                        !!nummerical test
                        !Fi_T_pure_p(j,i) = seggamma_pure_p(j,i) * Fi_T_pure_p(j,i)
                        !Fi_T_pure_m(j,i) = seggamma_pure_m(j,i) * Fi_T_pure_m(j,i)
                        !d2FdT2_pure(j,i) = (Fi_T_pure_p(j,i) - Fi_T_pure_m(j,i)) / (Temp_p - Temp_m)
                    end do

                    !Calculation Determinante of Matrix B (only needs to be calculated ones for each component)
                    !Matrix partial derivatives in respect to seggamma only
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            MatrixB(m,n) = gl.cosmo.dFi_pure_v1_gl(j,k,i)
                            MatrixB_pure(j,k,i) = MatrixB(m,n)
                            !!nummerical derivatives test
                            !MatrixB_pure_p(i,j,k) = Fi_seggamma_pure_p(j,k,i)
                            !MatrixB_p(m,n) = Fi_seggamma_pure_p(j,k,i)
                            !MatrixB_m(m,n) = Fi_seggamma_pure_m(j,k,i)
                            !MatrixB_pure_m(i,j,k) = Fi_seggamma_pure_m(j,k,i)
                        end do
                    end do

                    rankA = m

                    call Gauss_algorithm(gl,MatrixB, rankA, vectorb, vectorx, Det_seggamma, errval)
                    DetMatrixB_pure(i) = Det_seggamma
                    !!nummerical derivatives test
                    !call Gauss_algorithm(gl,MatrixB_p, rankA, vectorb, vectorx, Det_seggamma, errval)
                    !DetMatrixB_pure_p(i) = Det_seggamma
                    !!nummerical derivatives test
                    !call Gauss_algorithm(gl,MatrixB_m, rankA, vectorb, vectorx, Det_seggamma, errval)
                    !DetMatrixB_pure_m(i) = Det_seggamma
                    !
                    !dDetBdT_pure_num(:,i) = (DetMatrixB_pure_p(i) - DetMatrixB_pure_m(i)) / (Temp_p - Temp_m)

                    pos = 0

                    do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all derivatives of seggamma_pure in respect to T
                        if (abs(gl.cosmo.sigma(o,i)) < 1.D-14) then
                            cycle
                        end if
                        !save position of occupied intervalls of sigma profile to match dseggammadxa_mix to the respective intervall
                        !as this information is lost using gauss algorithm to solve dseggammadxa_mix
                        pos(o) = 1
                        m = 0
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile
                            if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                                cycle
                            end if
                            m = m + 1
                            n = 0
                            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung)
                                if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                    cycle
                                end if
                                n = n + 1
                                if (k == o) then
                                    MatrixA(m,n) = Fi_T_pure(j,i)
                                    MatrixA_pure(o,i,j,k) = Fi_T_pure(j,i)
                                    vectorb_gauss(m) = - Fi_T_pure(j,i)
                                    !!nummerical derivatives test
                                    !MatrixA_pure_p(o,i,j,k) = Fi_T_pure_p(j,i)
                                    !MatrixA_pure_m(o,i,j,k) = Fi_T_pure_m(j,i)
                                    !MatrixA_p(m,n) = Fi_T_pure_p(j,i)
                                    !MatrixA_m(m,n) = Fi_T_pure_m(j,i)
                                else
                                    !partial derivatives in respect to seggamma
                                    MatrixA(m,n) = gl.cosmo.dFi_pure_v1_gl(j,k,i)
                                    MatrixA_pure(o,i,j,k) = gl.cosmo.dFi_pure_v1_gl(j,k,i)
                                    !!nummerical derivatives test
                                    !MatrixA_pure_p(o,i,j,k) = Fi_seggamma_pure_p(j,k,i)
                                    !MatrixA_pure_m(o,i,j,k) = Fi_seggamma_pure_m(j,k,i)
                                    !MatrixA_p(m,n) = Fi_seggamma_pure_p(j,k,i)
                                    !MatrixA_m(m,n) = Fi_seggamma_pure_m(j,k,i)
                                end if
                            end do
                        end do

                        rankA = n
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !replaced with gauss algorithm
                        !vectorb = 0.D0
                        !call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_T, errval)
                        !DetMatrixA_pure(o,i) = Det_T
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!nummerical derivatives test
                        !call Gauss_algorithm(gl,MatrixA_p, rankA, vectorb, vectorx, Det_T, errval)
                        !DetMatrixA_pure_p(o,i) = Det_T
                        !!nummerical derivatives test
                        !call Gauss_algorithm(gl,MatrixA_m, rankA, vectorb, vectorx, Det_T, errval)
                        !DetMatrixA_pure_m(o,i) = Det_T
                        !!nummerical derivatives test
                        !dDetAdT_pure_num(o,i) = (DetMatrixA_pure_p(o,i) - DetMatrixA_pure_m(o,i)) / (Temp_p - Temp_m)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !replaced with gauss algorithm
                        !calculate dseggammadt_pure
                        !dseggammadt_pure(o,i) = - DetMatrixA_pure(o,i) / DetMatrixB_pure(i)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end do

                    !Gauss Algorithm is faster than Cramer's rule to calculate derivatives of system of equations
                    call Gauss_algorithm(gl,MatrixB, rankA, vectorb_gauss, vectorx, Det_T, errval)
                    dseggammadt_pure_gauss(:,i) = vectorx(1:51)

                    m = 0
                    do o = 1, 51
                        If (pos(o) == 0) then
                        else
                            m = m + 1
                            dseggammadt_pure(o,i) = dseggammadt_pure_gauss(m,i)
                        end if
                    end do
                    !for second Derivative wrt tau
                    DetMatrixA_pure(:,i) = - dseggammadt_pure(:,i) * DetMatrixB_pure(i)
                end do




                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !mixture
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                !In case segment activity coefficients were not calculated using Newton-Raphson, dFi_pure_v1_gl needs to be calculated here
                if (maxval(gl.cosmo.dFi_pure_v1_gl) < 1.D-14) then
                    !calculate gl.cosmo.dFi_pure_v1_gl
                    sum2_v1 = 0.0
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            sum2_v1(k) = sum2_v1(k) + gl.cosmo.sigma_profile_mix_gl(j) * gl.cosmo.seggamma_gl(j) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        end do
                        !derivatives of Fi after each seggamma results in two different equations
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (k == j) then
                                gl.cosmo.dFi_mix_v1_gl(k,j) = sum2_v1(k) + gl.cosmo.seggamma_gl(k) * gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            else
                                gl.cosmo.dFi_mix_v1_gl(k,j) = gl.cosmo.seggamma_gl(k) * gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            end if
                        end do
                    end do
                end if

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                ! Calculation Determinante of Matrix B (only needs to be calculated ones)
                !Matrix partial derivatives in respect to seggamma only
                m = 0
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    n = 0
                    if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        cycle
                    end if
                    m = m + 1
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            cycle
                        end if
                        n = n + 1
                        MatrixB(m,n) = gl.cosmo.dFi_mix_v1_gl(j,k)
                        MatrixB_mix(j,k) = MatrixB(m,n)
                        !MatrixB_mix_p(j,k) = Fi_seggamma_mix_p(j,k)
                        !MatrixB_p(m,n) = Fi_seggamma_mix_p(j,k)
                        !MatrixB_m(m,n) = Fi_seggamma_mix_m(j,k)
                        !MatrixB_mix_m(j,k) = Fi_seggamma_mix_m(j,k)
                    end do
                end do
                call Gauss_algorithm(gl,MatrixB, rankA, vectorb, vectorx, Det_seggamma, errval)
                DetMatrixB_mix = Det_seggamma
                ! !nummerical derivatives test
                !call Gauss_algorithm(gl,MatrixB_p, rankA, vectorb, vectorx, Det_seggamma, errval)
                !DetMatrixB_mix_p = Det_seggamma
                !!nummerical derivatives test
                !call Gauss_algorithm(gl,MatrixB_m, rankA, vectorb, vectorx, Det_seggamma, errval)
                !DetMatrixB_mix_m = Det_seggamma
                !
                !dDetBdT_mix_num(:) = (DetMatrixB_mix_p - DetMatrixB_mix_m) / (Temp_p - Temp_m)

                pos = 0

                do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all derivatives of seggamma in respect to T
                    if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                        cycle
                    end if
                    Fi_T = 0.D0
                    !Fi_T_p = 1.D-14
                    !Fi_T_m = 1.D-14
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all Fi_T's
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !sum for Fi_T
                            !partial derivatives in respect to Temperature
                            Fi_T(j,1) = Fi_T(j,1) + gl.cosmo.sigma_profile_mix_gl(k) * gl.cosmo.seggamma_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !!test for second derivative of Fi_T in respect to Temperature, can be removed after test
                            !Fi_T_p(j) = Fi_T_p(j) + gl.cosmo.sigma_profile_mix_gl(k) * seggamma_p(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_p ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_p))
                            !Fi_T_m(j) = Fi_T_m(j) + gl.cosmo.sigma_profile_mix_gl(k) * seggamma_m(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_m ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_m))
                            !d2FdTdseggamma_mix(j,k) = (Fi_seggamma_mix_p(j,k) - Fi_seggamma_mix_m(j,k)) / (Temp_p - Temp_m)
                        end do
                        Fi_T(j,1) = gl.cosmo.seggamma_gl(j) * Fi_T(j,1)
                        !!test for second derivative of Fi_T in respect to Temperature, can be removed after test
                        !Fi_T_p(j) = seggamma_p(j) * Fi_T_p(j)
                        !Fi_T_m(j) = seggamma_m(j) * Fi_T_m(j)
                        !d2FdT2_mix(j) = (Fi_T_p(j) - Fi_T_m(j)) / (Temp_p - Temp_m)
                    end do
                    !save position of occupied intervalls of sigma profile to match dseggammadxa_mix to the respective intervall
                    !as this information is lost using gauss algorithm to solve dseggammadxa_mix
                    pos(o) = 1
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)     !Zeile MatrixA
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        n = 0
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)     !Spalte (Ableitung) MatrixA
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            if (k == o) then
                                MatrixA(m,n) = Fi_T(j,1)
                                MatrixA_mix(o,j,k) = MatrixA(m,n)
                                vectorb_gauss(m) = - Fi_T(j,1)
                                !!nummerical derivatives test
                                !MatrixA_mix_p(o,j,k) = Fi_T_p(j)
                                !MatrixA_mix_m(o,j,k) = Fi_T_m(j)
                                !MatrixA_p(m,n) = Fi_T_p(j)
                                !MatrixA_m(m,n) = Fi_T_m(j)
                            else
                                !partial derivatives in respect to seggamma
                                MatrixA(m,n) = gl.cosmo.dFi_mix_v1_gl(j,k)
                                MatrixA_mix(o,j,k) = MatrixA(m,n)
                                !!nummerical derivatives test
                                !MatrixA_mix_p(o,j,k) = Fi_seggamma_mix_p(j,k)
                                !MatrixA_mix_m(o,j,k) = Fi_seggamma_mix_m(j,k)
                                !MatrixA_p(m,n) = Fi_seggamma_mix_p(j,k)
                                !MatrixA_m(m,n) = Fi_seggamma_mix_m(j,k)
                            end if
                        end do
                    end do

                    rankA = n
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !replaced with gauss algorithm
                    !vectorb = 0.D0
                    !
                    !call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_T, errval)
                    !DetMatrixA_mix(o) = Det_T
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    !!nummerical derivatives test
                    !call Gauss_algorithm(gl,MatrixA_p, rankA, vectorb, vectorx, Det_T, errval)
                    !DetMatrixA_mix_p(o) = Det_T
                    !!nummerical derivatives test
                    !call Gauss_algorithm(gl,MatrixA_m, rankA, vectorb, vectorx, Det_T, errval)
                    !DetMatrixA_mix_m(o) = Det_T
                    !!nummerical derivatives test
                    !dDetAdT_mix_num(o) = (DetMatrixA_mix_p(o) - DetMatrixA_mix_m(o)) / (Temp_p - Temp_m)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !replaced with gauss algorithm
                    !calculate dseggammadt_pure
                    !dseggammadt_mix(o) = - DetMatrixA_mix(o) / DetMatrixB_mix
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                end do
                !Gauss Algorithm is faster than Cramer's rule to calculate derivatives of system of equations
                call Gauss_algorithm(gl,MatrixB, rankA, vectorb_gauss, vectorx, Det_T, errval)
                dseggammadt_mix_gauss = vectorx(1:51)

                m = 0
                do o = 1, 51
                    If (pos(o) == 0) then
                    else
                        m = m + 1
                        dseggammadt_mix(o) = dseggammadt_mix_gauss(m)
                    end if
                end do
                !for second Derivative wrt tau
                DetMatrixA_mix = - dseggammadt_mix * DetMatrixB_mix

                !save for other derivations
                gl.cosmo.dseggamma_dT_mix = dseggammadt_mix

                !Derivations are dy/dx, need to be transformed to dlny/dx; dlny/dx = 1/y dy/dx
                do i = 1, gl.ncomp
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        dlnseggammadt_pure(j,i) = dseggammadt_pure(j,i) / gl.cosmo.seggamma_pure_gl(j,i)
                    end do
                end do
                do i = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    dlnseggammadt_mix(i) = dseggammadt_mix(i) / gl.cosmo.seggamma_gl(i)
                end do

                !!also T has to be tranformed to tau
                !dlnseggammadtau_mix = dlnseggammadt_mix * (-Temp ** 2.D0 / gl.tredmix)
                !do i = 1, gl.ncomp
                !    dlnseggammadtau_pure(:,i) = dlnseggammadt_pure(:,i) * (-Temp ** 2.D0 / gl.tredmix)
                !end do

                do i = 1, gl.ncomp
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        ln_gamma_R_Trho(4,i) = ln_gamma_R_Trho(4,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadt_mix(j) - dlnseggammadt_pure(j,i))
                    end do
                end do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !test
                !gl.cosmo.ln_gamma_R_Trho_test = ln_gamma_R_Trho
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                Do i = 1, gl.ncomp
                    gE_R_Trho(4) = gE_R_Trho(4) + R_const * x(i) * gl%ge%ln_gamma_R(1,i) + R_const * x(i) * Temp * ln_gamma_R_Trho(4,i)
                end do

                !Deriatives with respect to tau at constant delta and x
                gl%ge%ln_gamma_R(4,:) = - gl.tredmix / tau**2 * ln_gamma_R_Trho(4,:)
                gl%ge%gE_R(4) = - gl.tredmix / tau**2 * gE_R_Trho(4)


                !do i = 1, gl.ncomp
                !    do j = 1, 51
                !        !ln_gamma_R_test(4,i) = ln_gamma_R_test(4,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadtau_mix(j) - dlnseggammadtau_pure(j,i))
                !        gl.ln_gamma_R(4,i) = gl.ln_gamma_R(4,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadtau_mix(j) - dlnseggammadtau_pure(j,i))
                !    end do
                !end do

                !delete after testing
                !do i = 1, gl.ncomp
                !    do j = 1, 51
                !        ln_gamma_R_test(3,i) = ln_gamma_R_test(3,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadt_mix(j) - dlnseggammadt_pure(j,i))
                !    end do
                !end do
            else
                !Derivative wrt tau
                !gl%ge%gE_C(4) = (gE_C_p(1) - gE_C_m(1)) / (2.d0 * delta_tau)
                gl%ge%gE_R(4) = (gE_R_p(1) - gE_R_m(1)) / (2.d0 * delta_tau)
                do i = 1, gl%ncomp
                    !gl%ge%ln_gamma_C(4,i) = (ln_gamma_C_p(1,i) - ln_gamma_C_m(1,i)) / (2.d0 * delta_tau)
                    gl%ge%ln_gamma_R(4,i) = (ln_gamma_R_p(1,i) - ln_gamma_R_m(1,i)) / (2.d0 * delta_tau)
                end do
            end if
        end if

        !alternative analytic derivative wrt tau
        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then


                !initialize variables
                ln_gamma_R_Trho(5,:) = 0.D0
                gE_R_Trho(5) = 0.D0
                gl%ge%ln_gamma_R(5,:) = 0.D0
                gl%ge%gE_R(5) = 0.D0

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !pure substance
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                do i = 1, gl.ncomp
                    sum_d2Fi_dT2 = 0.D0
                    sum2_d2Fi_dT2 = 0.D0
                    sum_d2Fi_dseggammadT = 0.D0
                    d2Fi_dT2_pure = 0.D0
                    !D2Fi_dT2 for all derivatives of Fi
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            !if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                            !    cycle
                            !end if
                            !Sums for d2Fi_dT2
                            !sum_d2Fi_dT2(j) = sum_d2Fi_dT2(j) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !sum2_d2Fi_dT2(j) = sum2_d2Fi_dT2(j) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_pure(k,i) - gl.cosmo.seggamma_pure_gl(k,i) * 2.D0 / Temp + (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * gl.cosmo.seggamma_pure_gl(k,i))
                            d2Fi_dT2_pure(j,i) = d2Fi_dT2_pure(j,i) +  (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_pure(j,i) * gl.cosmo.seggamma_pure_gl(k,i) + gl.cosmo.seggamma_pure_gl(j,i) * (dseggammadt_pure(k,i) + gl.cosmo.seggamma_pure_gl(k,i) * (-2 / Temp) + gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0))))
                            !Sum for d2Fi_dseggammadT
                            sum_d2Fi_dseggammadT(j) = sum_d2Fi_dseggammadT(j) + (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_pure(k,i) + gl.cosmo.seggamma_pure_gl(k,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))
                        end do
                        !calculate d2Fi_dseggammadT_pure
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            !if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                            !    cycle
                            !end if
                            if (j == k) then
                                d2Fi_dseggammadT_pure(j,k,i) = sum_d2Fi_dseggammadT(j) + (gl.cosmo.sigma(j,i) / gl.cosmo.ACOSMO(i)) * exp(-gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp)) * (dseggammadt_pure(j,i) + gl.cosmo.seggamma_pure_gl(j,i) * (gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp ** 2.D0)))
                            else
                                d2Fi_dseggammadT_pure(j,k,i) = (gl.cosmo.sigma(k,i) / gl.cosmo.ACOSMO(i)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_pure(j,i) + gl.cosmo.seggamma_pure_gl(j,i) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))
                            end if
                        end do
                    end do

                    !Matrix dB/dT, second partial derivatives of B in respect to T
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            dMatrixBdT(m,n) = d2Fi_dseggammadT_pure(j,k,i)
                        end do
                    end do

                    !Adjugate of MatrixB
                    MatrixB = 0.D0
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            MatrixB(m,n) = MatrixB_pure(j,k,i)
                        end do
                    end do

                    rankA = m

                    call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                    !Multiply adj(MatrixA) * dMatrixAdT
                    call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)

                    !calculate dDetA/dT = Trace(Adjugate(A) * dMatrixA/dT
                    call Trace(gl,MatrixD, rankA, dDetBdT_pure(i), errval)


                    do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)   !all derivatuves of seggamma
                        if (abs(gl.cosmo.sigma(o,i)) < 1.D-14) then
                            cycle
                        end if
                        !matrix dA/dT, second partial derivatives of A in respect to T
                        m = 0
                        dMatrixAdT = dMatrixBdT
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)     !Zeile seggamma
                            if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                                cycle
                            end if
                            m = m + 1
                            n = 0
                            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung)
                                if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                    cycle
                                end if
                                n = n + 1
                                if (k == o) then
                                    dMatrixAdT(m,n) = d2Fi_dT2_pure(j,i)
                                else
                                end if
                            end do
                        end do



                        !Adjugates of MatrixA and MatrixB
                        MatrixA = 0.D0
                        !MatrixB = 0.D0
                        m = 0
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            n = 0
                            if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                                cycle
                            end if
                            m = m + 1
                            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                                if (abs(gl.cosmo.sigma(k,i)) < 1.D-14) then
                                    cycle
                                end if
                                n = n + 1
                                MatrixA(m,n) = MatrixA_pure(o,i,j,k)
                                !MatrixB(m,n) = MatrixB_pure(j,k,i)
                            end do
                        end do

                        rankA = m

                        call Adjugate(gl,MatrixA, rankA, adj_A, errval)
                        !call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                        !Multiply adj(MatrixA) * dMatrixAdT
                        call Mat_mult(gl,adj_A, rankA, rankA, dMatrixAdT, rankA, rankA, MatrixC, errval)
                        !call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)

                        !calculate dDetA/dT = Trace(Adjugate(A) * dMatrixA/dT
                        call Trace(gl,MatrixC, rankA, dDetAdT_pure(o,i), errval)
                        !call Trace(gl,MatrixD, rankA, dDetBdT_pure(o,i), errval)

                        !dDetAdT_pure(o,i) = dDetAdT_pure(o,i)

                    end do
                end do

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !mixture
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                !sum_d2Fi_dT2 = 0.D0
                !sum2_d2Fi_dT2 = 0.D0
                sum_d2Fi_dseggammadT = 0.D0
                d2Fi_dT2_pure_test = 0.D0
                sum3_d2Fi_dT2 = 0.D0
                d2Fi_dT2_mix = 0.D0
                !D2Fi_dT2 for all derivatives of Fi
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    !    cycle
                    !end if
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        !    cycle
                        !end if
                        !Sums for d2Fi_dT2
                        !sum_d2Fi_dT2(j) = sum_d2Fi_dT2(j) + gl.cosmo.sigma_profile_mix_gl(k) * gl.cosmo.seggamma_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        !sum2_d2Fi_dT2(j) = sum2_d2Fi_dT2(j) + gl.cosmo.sigma_profile_mix_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(k) - gl.cosmo.seggamma_gl(k) * 2.D0 / Temp + (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * gl.cosmo.seggamma_gl(k))
                        d2Fi_dT2_mix(j) = d2Fi_dT2_mix(j) +  gl.cosmo.sigma_profile_mix_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(j) * gl.cosmo.seggamma_gl(k) + gl.cosmo.seggamma_gl(j) * (dseggammadt_mix(k) + gl.cosmo.seggamma_gl(k) * (-2.D0 / Temp) + gl.cosmo.seggamma_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0))))
                        !sum3_d2Fi_dT2(j,1) = sum3_d2Fi_dT2(j,1) + gl.cosmo.sigma_profile_mix_gl(k) * gl.cosmo.seggamma_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        !sum3_d2Fi_dT2(j,2) = sum3_d2Fi_dT2(j,2) + gl.cosmo.sigma_profile_mix_gl(k) * dseggammadt_mix(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        !sum3_d2Fi_dT2(j,3) = sum3_d2Fi_dT2(j,3) + gl.cosmo.sigma_profile_mix_gl(k) * gl.cosmo.seggamma_gl(k) * (- 2.D0 * gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 3.D0)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        !sum3_d2Fi_dT2(j,4) = sum3_d2Fi_dT2(j,4) + gl.cosmo.sigma_profile_mix_gl(k) * gl.cosmo.seggamma_gl(k) * ((gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0))) ** 2.D0 * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        !Sum for d2Fi_dseggammadT
                        sum_d2Fi_dseggammadT(j) = sum_d2Fi_dseggammadT(j) + gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(k) + gl.cosmo.seggamma_gl(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))

                    end do
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        !    cycle
                        !end if
                        if (j == k) then
                            d2Fi_dseggammadT_mix(j,k) = sum_d2Fi_dseggammadT(j) + gl.cosmo.sigma_profile_mix_gl(j) * exp(-gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp)) * (dseggammadt_mix(j) + gl.cosmo.seggamma_gl(j) * (gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp ** 2.D0)))
                        else
                            d2Fi_dseggammadT_mix(j,k) = gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(j) + gl.cosmo.seggamma_gl(j) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))
                        end if
                    end do
                end do
                !Matrix dB/dT, second partial derivatives of B in respect to T
                m = 0
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    n = 0
                    if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        cycle
                    end if
                    m = m + 1
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            cycle
                        end if
                        n = n + 1
                        dMatrixBdT(m,n) = d2Fi_dseggammadT_mix(j,k)
                    end do
                end do

                !Adjugates MatrixB
                MatrixB = 0.D0
                m = 0
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    n = 0
                    if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        cycle
                    end if
                    m = m + 1
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            cycle
                        end if
                        n = n + 1
                        MatrixB(m,n) = MatrixB_mix(j,k)
                    end do
                end do

                rankA = m

                call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                !Multiply adj(MatrixA) * dMatrixAdT
                call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)

                !calculate dDetA/dT = Trace(Adjugate(A) dMatrixA/dT)
                call Trace(gl,MatrixD, rankA, dDetBdT_mix, errval)

                !dDetAdT_mix = 0.D0
                !matrix dA/dT, partial derivatives of A in respect to T
                do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)   !all derivatives of seggamma
                    if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                        cycle
                    end if
                    m = 0
                    dMatrixAdT = dMatrixBdT
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        n = 0
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            if (k == o) then
                                dMatrixAdT(m,n) = d2Fi_dT2_mix(j)
                            else
                            end if
                        end do
                    end do

                    !Adjugates of MatrixA and MatrixB
                    MatrixA = 0.D0
                    !MatrixB = 0.D0
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            MatrixA(m,n) = MatrixA_mix(o,j,k)
                            !MatrixB(m,n) = MatrixB_mix(j,k)
                        end do
                    end do

                    rankA = m

                    call Adjugate(gl,MatrixA, rankA, adj_A, errval)
                    !call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                    !Multiply adj(MatrixA) * dMatrixAdT
                    call Mat_mult(gl,adj_A, rankA, rankA, dMatrixAdT, rankA, rankA, MatrixC_mix, errval)
                    !call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)

                    !calculate dDetA/dT = Trace(Adjugate(A) dMatrixA/dT)
                    call Trace(gl,MatrixC_mix, rankA, dDetAdT_mix(o), errval)
                    !call Trace(gl,MatrixD, rankA, dDetBdT_mix(o), errval)

                    !dDetAdT_mix(o) = dDetAdT_mix(o)

                end do

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !Calculate second derivative of segment gamma wrt T for pure substance and mixture
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                !pure substance
                do i = 1, gl.ncomp
                    do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        d2seggammadT2_pure(o,i) = - ( dDetAdT_pure(o,i) * detMatrixB_pure(i) - dDetBdT_pure(i) * detMatrixA_pure(o,i)) / (detMatrixB_pure(i) ** 2.D0)
                    end do
                end do

                !mixture
                do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    d2seggammadT2_mix(o) = - ( dDetAdT_mix(o) * detMatrixB_mix - dDetBdT_mix * detMatrixA_mix(o)) / (detMatrixB_mix ** 2.D0)
                end do

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !Derivatives need to be for tau not Temp
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                !Derivations are d2y/dx2, need to be transformed to d2lny/dx2; d2lny/dx2 = -1/y^2 * (dy/dx)^2 + 1/y d2y/dx^2
                do i = 1, gl.ncomp
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma(j,i)) < 1.D-14) then
                        !    cycle
                        !end if
                        d2lnseggammadT2_pure(j,i) = - 1.D0 / gl.cosmo.seggamma_pure_gl(j,i) ** 2.D0 * dseggammadt_pure(j,i) ** 2.D0 + d2seggammadT2_pure(j,i) / gl.cosmo.seggamma_pure_gl(j,i)
                    end do
                end do
                do i = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !if (abs(gl.cosmo.sigma_profile_mix_gl(i)) < 1.D-14) then
                    !    cycle
                    !end if
                    d2lnseggammadT2_mix(i) =  - 1.D0 / gl.cosmo.seggamma_gl(i) ** 2.D0 * dseggammadt_mix(i) ** 2.D0 + d2seggammadT2_mix(i) / gl.cosmo.seggamma_gl(i)
                end do

                !T has to be tranformed to tau
                !d2lnseggammadtau2_mix = 2 * Temp ** 3.D0 / gl.tredmix ** 2.D0 * dlnseggammadt_mix + Temp ** 4.D0 / gl.tredmix ** 2.D0 * d2lnseggammadT2_mix
                !
                !do i = 1, gl.ncomp
                !    d2lnseggammadtau2_pure(:,i) = 2 * Temp ** 3.D0 / gl.tredmix ** 2.D0 * dlnseggammadt_pure(:,i) + Temp ** 4.D0 / gl.tredmix ** 2.D0 * d2lnseggammadT2_pure(:,i)
                !end do

                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !Calculate Second derivative wrt tau
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                do i = 1, gl.ncomp
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !constant Tau, rho
                        !gl.ln_gamma_R(5,i) = gl.ln_gamma_R(5,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (d2lnseggammadtau2_mix(j) - d2lnseggammadtau2_pure(j,i))
                        !constant T, rho
                        ln_gamma_R_Trho(5,i) = ln_gamma_R_Trho(5,i) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (d2lnseggammadT2_mix(j) - d2lnseggammadT2_pure(j,i))
                    end do
                end do


                Do i = 1, gl.ncomp
                    gE_R_Trho(5) = gE_R_Trho(5) + 2.D0 * R_Const * x(i) * ln_gamma_R_Trho(4,i) + R_Const * Temp * x(i) * ln_gamma_R_Trho(5,i)
                end do

                !Deriatives with respect to tau at constant delta and x
                gl%ge%ln_gamma_R(5,:) = 2.D0*gl.tredmix/tau**3*ln_gamma_R_Trho(4,:)+gl.tredmix**2/tau**4*ln_gamma_R_Trho(5,:)
                gl%ge%gE_R(5) = 2.D0*gl.tredmix/tau**3*gE_R_Trho(4)+gl.tredmix**2/tau**4*gE_R_Trho(5)

                !save for other derivations
                gl.cosmo.d2Fi_dseggamma_dT_mix = d2Fi_dseggammadT_mix

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !test, remove after test
                !gl.gE_R(5) = gE_R_num(5)
                !gl.cosmo.gE_R_Trho_num = gE_R_Trho_num
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else
                !Second derivative wrt tau
                !gl.gE_C(5) = (gE_C_p(1) - 2.d0 * gl.gE_C(1) + gE_C_m(1)) / (delta_tau ** 2)
                gl%ge%gE_R(5) = (gE_R_p(1) - 2.d0 * gE_R_act(1) + gE_R_m(1)) / (delta_tau ** 2)
                do i = 1, gl.ncomp
                    !   gl.ln_gamma_C(5,i) = (ln_gamma_C_p(1,i) - 2.d0 * gl.ln_gamma_C(1,i) + ln_gamma_C_m(1,i)) / (delta_tau ** 2)
                    gl%ge%ln_gamma_R(5,i) = (ln_gamma_R_p(1,i) - 2.d0 * ln_gamma_R_act(1,i) + ln_gamma_R_m(1,i)) / (delta_tau ** 2)
                end do
            end if
        end if


        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho(6,:) = 0.D0
            gE_R_Trho(6) = 0.d0
            gl%ge%ln_gamma_R(6,:) = 0.d0
            gl%ge%gE_R(6) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho(7,:) = 0.D0
            gE_R_Trho(7) = 0.d0
            gl%ge%ln_gamma_R(7,:) = 0.d0
            gl%ge%gE_R(7) = 0.D0
        end if

        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho(8,:) = 0.D0
            gE_R_Trho(8) = 0.d0
            gl%ge%ln_gamma_R(8,:) = 0.d0
            gl%ge%gE_R(8) = 0.D0
        end if

        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho(9,:) = 0.D0
            gE_R_Trho(9) = 0.d0
            gl%ge%ln_gamma_R(9,:) = 0.d0
            gl%ge%gE_R(9) = 0.D0      !TO BE IMPLEMENTED
        end if

        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho(10,:) = 0.D0
            gE_R_Trho(10) = 0.d0
            gl%ge%ln_gamma_R(10,:) = 0.d0
            gl%ge%gE_R(10) = 0.D0
        end if

        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho(11,:) = 0.D0
            gE_R_Trho(11) = 0.d0
            gl%ge%ln_gamma_R(11,:) = 0.d0
            gl%ge%gE_R(11) = 0.D0
        end if

        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho(12,:) = 0.D0
            gE_R_Trho(12) = 0.d0
            gl%ge%ln_gamma_R(12,:) = 0.d0
            gl%ge%gE_R(12) = 0.D0
        end if

        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho(13,:) = 0.D0
            gE_R_Trho(13) = 0.d0
            gl%ge%ln_gamma_R(13,:) = 0.d0
            gl%ge%gE_R(13) = 0.D0
        end if

        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho(14,:) = 0.D0
            gE_R_Trho(14) = 0.d0
            gl%ge%ln_gamma_R(14,:) = 0.d0
            gl%ge%gE_R(14) = 0.D0
        end if

        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho(15,:) = 0.D0
            gE_R_Trho(15) = 0.d0
            gl%ge%ln_gamma_R(15,:) = 0.d0
            gl%ge%gE_R(15) = 0.D0     !TO BE IMPLEMENTED
        end if

        if (C_or_R == 2) then
            gl%ge%gE_C = 0.D0
            gl%ge%ln_gamma_C =0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    End if

    !Calculate the overall activity coefficient and the excess Gibbs energy (not needed, only for checking results)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !gE = 0.d0
    !ln_gamma = 0.D0
    !Do j = 1,nderivs
    !    Do i = 1, ncomp
    !        ln_gamma(j,i) = ln_gamma_C(j,i) + ln_gamma_R(j,i)
    !        gE(j) = gE(j) + x(i) * ln_gamma(j,i) * R_const_cal * Temp
    !    end do
    !End do
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    end subroutine
    !--------------------------------------------------------------------------------------


    !Erik, April 2018
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine gE_COSMO_SAC_CALC(gl, Temp, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Calculate the activity coefficients and the excess Gibbs energy from COSMO SAC

    implicit none

    type(type_gl) :: gl

    double precision :: Temp
    integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed
    integer:: C_or_R
    integer:: errval
    double precision, dimension(51) :: Sigma_area_mix           !total surface area from all segments with same surface charge of mixture
    double precision, dimension(51) :: sigma_profile_mix        !Sigma-profile of the mixture
    double precision, dimension(51,gl%ncomp) :: sigma_profile_pure    !sigma-profile of the pure components
    double precision :: ACOSMO_mix                      !COSMO Area of mixture
    double precision :: sigmaacc, sigmadon              !Acceptor and Donator (hydrogen bonding interactions) in calculation of electrostatic interaction
    double precision, dimension(:, :), allocatable :: delta_w
    double precision :: aeff
    double precision :: coord
    double precision, dimension(51) :: seggamma
    double precision :: r_gas, vnorm, anorm
    double precision, dimension(51, gl%ncomp) :: seggamma_pure
    double precision, dimension(30) :: rnorm, qnorm, phi_cosmo, theta_cosmo, l_cosmo
    double precision :: rnorm_mix, qnorm_mix, sum_l_cosmo
    double precision, dimension(30) :: ln_gamma_SG, sum_gamma_res, gamma_cosmo, ln_gamma_cosmo
    integer :: i, j, k, l, m, n, o, p
    double precision :: R_cosmo
    double precision :: A_ES, B_ES, c_es
    double precision, dimension(:,:,:,:), allocatable :: delta_w_v23
    double precision, dimension(51,3) :: Sigma_area_mix_v23, seggamma_old_v23, converg_v23, seggamma_v23, sigma_profile_mix_v23
    double precision, dimension(:,:,:), allocatable :: seggamma_pure_v23
    double precision, dimension(:,:,:), allocatable :: sigma_profile_pure_v23
    logical :: recalc
    double precision, dimension(30,30) :: A_dsp
    double precision :: w_dsp, sum_xi_Aki, sum_Ajk_xj, sum_Aji_xj, sum_xi_sum_xj_Aji
    double precision :: k_b_cosmo   !Boltzmann constant
    logical :: faster

    integer, dimension(3) :: start_sigma, end_sigma
    !dispersive interactions for version 4 (Xiong); this model does not work
    double precision, dimension(17) :: eps_element = (/ 42976.D0,	36380.D0,	19708.D0,	6753.D0,	16031.D0,	8252.D0,	7799.D0,	4979.D0,	4497.D0,	13501.D0,	7766.D0,    67738.D0,   59665.D0,   46518.D0,   52577.D0,   115568.D0,  0.D0/)
    double precision, dimension(30,17) :: eps_dsp, m_tau_dsp
    integer, dimension(30) :: n_a_dsp
    double precision :: V_mix, sum_xi_xj_mk_ml_epskl
    double precision, dimension(30) :: sum_xi_mk_ml_epskl
    double precision, dimension(30) :: V_molar_cosmo
    double precision :: T_red_cosmo
    double precision, dimension(51,3) :: test_cosmo
    integer, dimension(30) :: check_seggamma
    integer :: interval
    integer :: cycle_count

    !November 2018, Erik
    !speed improvements
    double precision, dimension(51,51) :: AA
    character(len=255) :: fileout
    double precision :: conv_criteria
    logical :: NR

    !speed test
    integer :: start, ende, rate, cmax
    real :: time

    if(.not. allocated(delta_w_v23)) allocate(delta_w_v23(51,51,3,3))
    if(.not. allocated(delta_w)) allocate(delta_w(51,51))
    if(.not. allocated(seggamma_pure_v23)) allocate(seggamma_pure_v23(51,gl%ncomp,3))
    if(.not. allocated(sigma_profile_pure_v23)) allocate(sigma_profile_pure_v23(51,gl%ncomp,3))

    !convergence criteria
    conv_criteria = 1.d-8
    interval = 51

    !new universal gas constant
    R_cosmo = 1.38064903D-23 * 6.022140758D23

    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3)) then    !Hsieh
        aeff = 7.25D0    !surface area of a standard segment in Angstrom^2
    elseif (gl%cosmo%COSMO_ver == 1) then    !Mullins
        aeff = 7.5D0    !surface area of a standard segment in Angstrom^2
    elseif (gl%cosmo%COSMO_ver == 4) then    !Xiong
        aeff = 7.9D0    !surface area of a standard segment in Angstrom^2
    end if

    r_gas= R_cosmo / 4184.D0    !universal gas constant in kcal/(mol*K)

    vnorm= 66.69D0
    anorm = 79.53D0
    coord=10    !coordination number, Klamt used 7.2
    A_ES = 6525.69D0 * 4184.D0    !4184 to convert from (kcal/mol)*(A^4/e^2) to (J/mol)*(A^4/e^2); for electrostatic interaction coefficient
    B_ES = 1.4859D0 * 10 ** 8 * 4184.D0
    k_b_cosmo = 1.380649*1.D-23  !Boltzmann constant

    !Intialize
    Sigma_area_mix = 0.D0
    Sigma_area_mix_v23 = 0.D0
    ACOSMO_mix = 0.D0
    check_seggamma = 0
    seggamma_v23 = 0.d0
    seggamma_pure_v23 = 0.d0

    recalc =.true.
    faster = gl.cosmo.solv_method
    !molescheck = .false.
    NR = .false.

    if(dabs(gl.cosmo.temp_cosmo_prev - Temp) > 1.D-16) then
        recalc = .true.
    end if
    do i=1,gl.ncomp
        if(dabs(gl.cosmo.molefraction_cosmo_prev(i) - gl.molfractions(i)) > 1.D-16) then
            recalc = .true.
        end if
        if(gl%ge%Eq_type_prev(i) .ne. gl.Eq_type(i)) then
            recalc = .true.
        end if
    end do

    if (recalc == .false.) then
        gl%ge%ln_gamma_C = gl.cosmo.ln_gamma_C_cosmo_prev
        gl%ge%ln_gamma_R = gl.cosmo.ln_gamma_R_cosmo_prev
        gl%ge%gE_C = gl.cosmo.gE_C_cosmo_prev
        gl%ge%gE_R = gl.cosmo.gE_R_cosmo_prev
        return
    end if

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !CALCULATE THE MIXTURE SIGMA PROFILE
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (gl%cosmo%COSMO_ver == 1) then
        !Calculate total surface area from all segments with same surface charge of mixture
        do j = 1, interval
            !Sigma_area_mix(j)=0
            do i = 1, gl%ncomp
                Sigma_area_mix(j) = Sigma_area_mix(j) + gl%molfractions(i) * gl%cosmo%sigma(j,i)
            end do
        end do
        !Calculate COSMO area of the mixture
        do i=1,gl%ncomp
            ACOSMO_mix=ACOSMO_mix+gl%molfractions(i)*gl%cosmo%ACOSMO(i)
        end do
        !Calculate the sigma-profile of the mixture
        do j=1,interval
            sigma_profile_mix(j)=Sigma_area_mix(j)/ACOSMO_mix
        end do

        !Calculate the sigma-profiles of the pure components
        do i = 1, gl%ncomp
            do j = 1, 51
                sigma_profile_pure(j,i) = gl.cosmo.sigma(j,i) / gl.cosmo.ACOSMO(i)
            end do
        end do

        start_sigma(1) = gl.cosmo.interval_min(1)
        end_sigma(1) = gl.cosmo.interval_max(1)
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    elseif ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
        !Calculate total surface area from all segments with same surface charge of mixture
        do k = 1, 3
            do j = 1, interval
                do i = 1, gl%ncomp
                    Sigma_area_mix_v23(j,k) = Sigma_area_mix_v23(j,k) + gl%molfractions(i) * gl%cosmo%sigma_v23(j,i,k)
                end do
            end do
        end do
        !Calculate COSMO area of the mixture; for versions 2 and 3 of COSMO-SAC
        do i=1,gl%ncomp
            ACOSMO_mix = ACOSMO_mix + gl%molfractions(i) * gl%cosmo%ACOSMO(i)
        end do
        !Calculate the sigma-profile of the mixture
        do k = 1, 3
            do j=1,interval
                sigma_profile_mix_v23(j,k) = Sigma_area_mix_v23(j,k) / ACOSMO_mix
            end do
        end do

        !Calculate the sigma-profiles of the pure components
        do i = 1, gl%ncomp
            do k = 1, 3
                do j = 1, 51
                    sigma_profile_pure_v23(j,i,k) = gl.cosmo.sigma_v23(j,i,k) / gl.cosmo.ACOSMO(i)
                end do
            end do
        end do

        start_sigma = gl.cosmo.interval_min
        end_sigma = gl.cosmo.interval_max
    end if


    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Electrostatic interaction for resdiual part of activity coefficient
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (gl.cosmo.COSMO_ver == 1) then
        !Calculate electrostatic interaction COSMO-SAC Version 1
        delta_w(start_sigma(1):end_sigma(1),start_sigma(1):end_sigma(1)) = gl.cosmo.delta_w_gl(start_sigma(1):end_sigma(1),start_sigma(1):end_sigma(1))
    elseif ((gl.cosmo.COSMO_ver == 2) .or. (gl.cosmo.COSMO_ver == 3) .or. (gl.cosmo.COSMO_ver == 4)) then
        !Calculate electrostatic interaction COSMO-SAC Version 2 and 3
        !Electrostatic interaction coefficient c_es
        if ((gl.cosmo.COSMO_ver == 2) .or. (gl.cosmo.COSMO_ver == 3)) then
            c_es = A_ES + B_ES / Temp ** 2
        elseif (gl.cosmo.COSMO_ver == 4) then
            c_es = 9254.D0 * 4184.D0
        end if
        delta_w_v23 = c_es * gl.cosmo.sum_square_counter - gl.cosmo.hb_inter
    end if

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Segment activity coefficient COSMO version 1
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (gl.cosmo.COSMO_ver == 1) then
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !Newton-Raphson for mixtures
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if (faster == .true.) then
            call Newton_Raphson_cosmo1(gl, Temp, delta_w, r_gas, sigma_profile_mix, seggamma)
        end if

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !successive substitution method for mixtures
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if ((faster == .false.) .or. (errval == -12000)) then
            if (errval == -12000) then
                errval = 0
                seggamma = 1.0D0
            end if

            call Succ_Sub_cosmo1(gl, Temp, delta_w, r_gas, sigma_profile_mix, NR, seggamma)

        end if

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES)
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !start values for segment activity coefficients for Pure Component
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        !Use standard initial values of 1, May 2018, changed to 1.1 due to occurance of solving the equations with negative seggamma
        seggamma_pure = 1.0D0

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !Newton-Raphson for pure fluids
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if (faster == .true.) then
            do p = 1, gl.ncomp
                call Newton_Raphson_cosmo1(gl, Temp, delta_w, r_gas, sigma_profile_pure(:,p), seggamma_pure(:,p))
            end do
        end if

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !successive substitution method for pure fluids
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if ((faster == .false.) .or. (errval == -12000)) then
            if (errval == -12000) then
                errval = 0
                seggamma_pure = 1.0D0
            end if
            do i = 1, gl.ncomp
                call Succ_Sub_cosmo1(gl, Temp, delta_w, r_gas, sigma_profile_pure(:,i), NR, seggamma_pure(:,i))
            end do
        end if
    end if

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Segment activity coefficient COSMO Version 2,3 and 4
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !ITERATION FOR SEGMENT ACTIVITY COEF. (MIXTURE)
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !Newton-Raphson for mixtures
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if (faster == .true.) then
            call Newton_Raphson_cosmo3(gl, Temp, delta_w_v23, R_cosmo, sigma_profile_mix_v23, gl%cosmo%seggamma_prev, seggamma_v23)
        end if

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !successive substitution method for mixtures
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if ((faster == .false.) .or. (errval == -12000)) then
            if (errval == -12000) then
                errval = 0
                seggamma_v23 = 1.D0
            end if
            call Succ_Sub_cosmo3(gl, Temp, delta_w_v23, R_cosmo, sigma_profile_mix_v23, NR, gl%cosmo%seggamma_prev, seggamma_v23, errval)
        end if

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES)
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !Newton-Raphson method
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if (faster == .true.) then
            do p = 1, gl.ncomp
                call Newton_Raphson_cosmo3(gl, Temp, delta_w_v23, R_cosmo, sigma_profile_pure_v23(:,p,:), gl%cosmo%seggamma_pure_prev(:,p,:), seggamma_pure_v23(:,p,:))
            end do
        end if
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !Successive substitution method for pure fluids in case of not converging, negative seggamma or as "slower" option
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if ((faster == .false.) .or. (maxval(check_seggamma) > 0)) then
            do i = 1, gl%ncomp
                call Succ_Sub_cosmo3(gl, Temp, delta_w_v23, R_cosmo, sigma_profile_pure_v23(:,i,:), NR, gl%cosmo%seggamma_pure_prev(:,i,:), seggamma_pure_v23(:,i,:), errval)
            end do
        end if
    end if

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !calculation of residual contribution term
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (gl%cosmo%COSMO_ver == 1) then
        sum_gamma_res = 0.0
        do i = 1, gl%ncomp
            do j = start_sigma(1), end_sigma(1)
                sum_gamma_res(i) = sum_gamma_res(i) + ((gl%cosmo%sigma(j,i) / aeff) * (log(seggamma(j) / (seggamma_pure(j,i)))))
            end do
        end do
    end if

    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
        sum_gamma_res = 0.D0
        do i = 1, gl%ncomp
            do k = 1, 3
                if (start_sigma(k) == 0) then
                    cycle
                end if
                do j = start_sigma(k), end_sigma(k)
                    sum_gamma_res(i) = sum_gamma_res(i) + ((gl.cosmo.sigma_v23(j,i,k) / aeff) * (log(seggamma_v23(j,k) / (seggamma_pure_v23(j,i,k)))))
                end do
            end do
        end do
    end if

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Combinatorial term
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    !THE STAVERMAN-GUGGENHEIM EQUATION
    do i = 1, gl%ncomp
        rnorm(i) = gl%cosmo%VCOSMO(i)/vnorm
        qnorm(i) = gl%cosmo%ACOSMO(i)/anorm
    end do

    qnorm_mix = 0.d0
    rnorm_mix = 0.d0

    do i = 1, gl%ncomp
        qnorm_mix = qnorm_mix + gl%molfractions(i) * qnorm(i)
        rnorm_mix = rnorm_mix + gl%molfractions(i) * rnorm(i)
    end do

    do i = 1, gl%ncomp
        theta_cosmo(i) = gl%molfractions(i) * qnorm(i) / qnorm_mix
        phi_cosmo(i) = gl%molfractions(i) * rnorm(i) / rnorm_mix
        l_cosmo(i) = (coord/2.0) * (rnorm(i) - qnorm(i)) - (rnorm(i) - 1.0)
    end do

    sum_l_cosmo = 0.D0

    do i = 1, gl%ncomp
        sum_l_cosmo = sum_l_cosmo + gl%molfractions(i) * l_cosmo(i)
    end do

    do i = 1, gl%ncomp
        if (gl%molfractions(i) < 1.D-14) then
            ln_gamma_SG(i) = log(rnorm(i) / rnorm_mix) + (coord / 2)*qnorm(i) * log(qnorm(i) / qnorm_mix / (rnorm(i) / rnorm_mix)) + l_cosmo(i) - (rnorm(i) / rnorm_mix) * sum_l_cosmo
        else
            ln_gamma_SG(i) = log(phi_cosmo(i) / gl%molfractions(i)) + (coord / 2)*qnorm(i) * log(theta_cosmo(i) / phi_cosmo(i)) + l_cosmo(i) - (phi_cosmo(i) / gl%molfractions(i)) * sum_l_cosmo
        end if
    end do

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !dispersive interactions
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    !Calculate dispersive interactions
    gl%cosmo%ln_gamma_dsp = 0.D0
    if (gl%cosmo%COSMO_ver == 3) then
        !only do this if the eps_molecule values are not the specified error value eps_molecule = -7000D0 or below 0
        if (minval(gl%cosmo%eps_molecule(1:gl%ncomp)) > -1.D-10) then
            A_dsp = 0.D0
            do i = 1, gl%ncomp - 1
                do j = i + 1, gl%ncomp
                    !do j = 1, gl%ncomp
                    !calculate w_dsp and A_dsp
                    if ((gl%cosmo%molecule_type(i) == 2 .and. gl%cosmo%molecule_type(j) == 3) .or. (gl%cosmo%molecule_type(i) == 3 .and. gl%cosmo%molecule_type(j) == 2) .or. &    !water + hb-only-acceptor
                        (gl%cosmo%molecule_type(i) == 1 .and. (gl%cosmo%molecule_type(j) == 5 .or. gl%cosmo%molecule_type(j) == 4)) .or. ((gl%cosmo%molecule_type(i) == 5 .or. gl%cosmo%molecule_type(i) == 4) .and. gl%cosmo%molecule_type(j) == 1) .or. &    !COOH + nhb or hb-donor-acceptor
                        (gl%cosmo%molecule_type(i) == 2 .and. gl%cosmo%molecule_type(j) == 1) .or. (gl%cosmo%molecule_type(i) == 1 .and. gl%cosmo%molecule_type(j) == 2)) then  !water + COOH
                        w_dsp = -0.27027D0
                    else
                        w_dsp = 0.27027D0
                    end if
                    A_dsp(i,j) = w_dsp * (0.5D0 * (gl%cosmo%eps_molecule(i) + gl%cosmo%eps_molecule(j)) - sqrt(gl%cosmo%eps_molecule(i) * gl%cosmo%eps_molecule(j)))
                    A_dsp(j,i) = A_dsp(i,j) !mirror matrix A_dsp as A_dsp12 equals A_dsp21
                end do
            end do
            !Margules equation of multi component mixtures
            do k = 1, gl%ncomp
                sum_xi_Aki = 0.D0
                sum_Ajk_xj = 0.D0
                sum_xi_sum_xj_Aji = 0.D0
                do i = 1, gl%ncomp
                    sum_xi_Aki = sum_xi_Aki + gl%molfractions(i) * A_dsp(k,i)
                    sum_Ajk_xj = sum_xi_Aki
                    sum_Aji_xj = 0.D0
                    do j = 1, gl%ncomp
                        sum_Aji_xj = sum_Aji_xj + gl%molfractions(j) * A_dsp(j,i)
                    end do
                    sum_xi_sum_xj_Aji = sum_xi_sum_xj_Aji + gl%molfractions(i) * sum_Aji_xj
                end do
                gl%cosmo%ln_gamma_dsp(k) = 1.D0 / 2.D0 * (sum_Ajk_xj + sum_xi_Aki - sum_xi_sum_xj_Aji)
            end do

            !exeption for CO and SO3 as their eps_molecule is <0, set ln_gamma_dsp to 0, eps_molecule is now set to 0 for these components
            do i = 1, gl%ncomp
                if (gl%cosmo%eps_molecule(i) == 0.D0) then
                    gl%cosmo%ln_gamma_dsp = 0.D0
                end if
            end do
        else
            gl%cosmo%ln_gamma_dsp = 0.D0
        end if


    elseif (gl%cosmo%COSMO_ver == 4) then
        !Funktion zur epsilon bestimmung im Anhang---> funktion wird nicht verwendet. einfach array mit allen epsilons und über m_tau raussuchen welche element_dsp besetzt werden, um rechenzeit zu sparen
        !select atom types, effective area of atom types and corresponding epsilons which actually exist in given components
        do i = 1, gl%ncomp
            k = 0
            do j = 1, 17
                if (gl%cosmo%m_tau(i,j) > 1.D-16) then
                    k = k + 1
                    eps_dsp(i,k) = eps_element(j)
                    m_tau_dsp(i,k) = gl%cosmo%m_tau(i,j)
                end if
                !save number of element types for each component
                if (j == 17) then
                    n_a_dsp(i) = k
                end if
            end do
        end do
        !Molar volume of mixture asuming ideal mixing
        do i = 1, gl%ncomp
            T_red_cosmo = Temp / gl%cosmo%NBT_cosmo(i)
            V_molar_cosmo(i) = 0.000971 * T_red_cosmo ** 2 * gl%cosmo%VCOSMO(i) ** 2 + (0.233 * T_red_cosmo ** 2 + 1.14) * gl%cosmo%VCOSMO(i) - 4.42
        end do

        V_mix = 0.D0
        do i = 1, gl%ncomp
            V_mix = V_mix + V_molar_cosmo(i) * gl%molfractions(i)
        end do

        gl%cosmo%ln_gamma_dsp = 0.D0
        do i = 1, gl%ncomp
            sum_xi_xj_mk_ml_epskl = 0.D0
            do j = 1, gl%ncomp
                do k = 1, n_a_dsp(i)
                    do l = 1, n_a_dsp(j)
                        sum_xi_xj_mk_ml_epskl = sum_xi_xj_mk_ml_epskl + gl%molfractions(i) * gl%molfractions(j) * m_tau_dsp(i,k) * m_tau_dsp(j,l) * (eps_dsp(i,k) * eps_dsp(j,l)) ** 0.5D0
                    end do
                end do
            end do
        end do

        sum_xi_mk_ml_epskl = 0.D0
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                do k = 1, n_a_dsp(i)
                    do l = 1, n_a_dsp(j)
                        sum_xi_mk_ml_epskl(i) = sum_xi_mk_ml_epskl(i) + gl%molfractions(j) * m_tau_dsp(i,k) * m_tau_dsp(j,l) * (eps_dsp(i,k) * eps_dsp(j,l)) ** 0.5D0
                    end do
                end do
            end do
        end do

        do i = 1, gl%ncomp
            gl%cosmo%ln_gamma_dsp(i) = 1 / (R_cosmo * Temp * V_mix) * (gl%cosmo%VCOSMO(i) / V_mix * sum_xi_xj_mk_ml_epskl - 2 * sum_xi_mk_ml_epskl(i))
        end do
    end if


    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !CALCULATION OF GAMMAS for COSMO-SAC
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ! maybe for stand alone program, needs to be refreshed as the dispersion interaction is not included


    do i = 1, gl%ncomp
        gamma_cosmo(i) = exp(sum_gamma_res(i) + ln_gamma_SG(i) + gl%cosmo%ln_gamma_dsp(i))
        ln_gamma_cosmo(i) = log(gamma_cosmo(i))
    end do

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Calculation of ln_gamma's and gE's
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    gl%ge%gE_C(1)=0.D0
    gl%ge%gE_R(1)=0.D0
    do i = 1, gl%ncomp
        gl%ge%ln_gamma_C(1,i) = ln_gamma_SG(i)
        gl%ge%ln_gamma_R(1,i) = sum_gamma_res(i) + gl%cosmo%ln_gamma_dsp(i)
        !gl%ge%ln_gamma_R(1,i) = sum_gamma_res(i) !+ gl%cosmo%ln_gamma_dsp(i)    !for Hsieh test
        gl%ge%gE_C(1) = gl%ge%gE_C(1) + gl%molfractions(i) * gl%ge%ln_gamma_C(1,i) * R_cosmo * Temp
        gl%ge%gE_R(1) = gl%ge%gE_R(1) + gl%molfractions(i) * gl%ge%ln_gamma_R(1,i) * R_cosmo * Temp
    end do

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Save previous values
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    gl%cosmo%temp_cosmo_prev = Temp
    gl%cosmo%molefraction_cosmo_prev = gl%molfractions
    gl%cosmo%seggamma_pure_prev(:,1:gl%ncomp,:) = seggamma_pure_v23(:,1:gl%ncomp,:)
    gl%cosmo%seggamma_pure_prev_v1(:,1:gl%ncomp) = seggamma_pure(:,1:gl%ncomp)
    gl%ge%Eq_type_prev = gl%Eq_type
    gl%cosmo%seggamma_prev = seggamma_v23
    gl%cosmo%seggamma_prev_v1 =seggamma
    gl%cosmo%ln_gamma_R_cosmo_prev = gl%ge%ln_gamma_R
    gl%cosmo%gE_C_cosmo_prev = gl%ge%gE_C
    gl%cosmo%gE_R_cosmo_prev = gl%ge%gE_R

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !Save values for analytical derivations
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    !only save those values if analytical derivatives are on
    if (gl.cosmo.analytical == .true.) then

        gl.cosmo.seggamma_gl = seggamma
        gl.cosmo.seggamma_pure_gl(:,1:gl%ncomp) = seggamma_pure(:,1:gl%ncomp)
        gl.cosmo.seggamma_v23_gl = seggamma_v23
        gl.cosmo.seggamma_pure_v23_gl(:,1:gl%ncomp,:) = seggamma_pure_v23(:,1:gl%ncomp,:)
        gl.cosmo.aeff_gl = aeff
        gl.cosmo.sigma_profile_mix_gl = sigma_profile_mix
        gl.cosmo.sigma_profile_mix_v23_gl = sigma_profile_mix_v23
        gl.cosmo.delta_w_gl = delta_w
        gl.cosmo.delta_w_v23_gl = delta_w_v23
        gl.cosmo.start_sigma = start_sigma
        gl.cosmo.end_sigma = end_sigma

    end if
    end subroutine
    !------------------------------------------------------------------------------------------------------------------------------
    !**************************************************************************
    subroutine Gauss_algorithm_cosmo(MatrixA, rankA, vectorb, vectorx, Det_A, errorflag)
    !Gauß-algorithm to bring a matrix into upper diagonal form
    !This algorithm can be used to solve a linear system of equation of the
    !form A*x = b, where A is the matrix, x is the solution vector and b is
    !the right side vector of the system of equations.
    !Furthermore, this algorithm can be used to calculate the determinant of
    !a matrix
    !INPUT:
    !MatrixA: A quadratic matrix containing a maximum of 60 rows and columns
    !vectorb: The right hand side vector, containing up to 60 entries
    !RankA: How many columns and rows does the matrix have? (means: how many
    !       equations need to be solved
    !OUTPUT:
    !vectorx: Vector containing the solutions

    !Erik, Januar 2018
    !**************************************************************************

    implicit none

    double precision, dimension(160,160):: MatrixA
    integer:: rankA
    double precision, dimension(160):: vectorb
    double precision, dimension(160):: vectorx
    double precision:: det_A
    integer:: errorflag

    integer:: i,j,k,l
    double precision:: eps_zero
    double precision:: factor
    double precision, dimension(160,160):: MatA_Gauss

    eps_zero = 1.D-16   !Tolerance for values close to 0
    det_A = 0.D0


    !A rank smaller than 0 does not exist. Rank > 60 is not allowed because the matrix is restricted to 60 entries -> quit with error in both cases
    if ((rankA <= 0) .or. (rankA > 160)) then
        errorflag = -3322
        return
    end if

    !If the rank is 1, the solution is simple and treated separately
    if (rankA == 1) then
        Det_A = MatrixA(1,1)
        if (dabs(MatrixA(1,1)) > eps_zero) then
            vectorx(1) = vectorb(1) / MatrixA(1,1)
            return
        end if
    end if

    !Copy Matrix A to another matrix for Gauss transformations
    do i= 1, rankA
        do j = 1, rankA
            MatA_Gauss(i,j) = MatrixA(i,j)
        end do
    end do


    !If the rank is larger than 1:
    !Bring the matrix in form of an upper triangular matrix
    do j = 1, rankA-1     !Loop over the colums of matrix A
        !When starting to look at a new column, the element on the diagonal must have a value which is not equal to 0
        !This is ensured here:
        !1) Check whether the entry is 0, if yes, go through other rows k > j and subtract a row which has not a value of 0 in column j . If the entry is not equal to 0, nothing needs to be done
        If (dabs(MatA_Gauss(j,j)) < eps_zero) then
            Do k = j+1, rankA
                if (dabs(MatA_Gauss(k,j)) > eps_zero) then
                    do l = j, rankA
                        MatA_Gauss(j,l) = MatA_Gauss(j,l) - MatA_Gauss(k,l)
                    end do
                    vectorb(j) = vectorb(j) - vectorb(k)
                    exit
                end if
                if (k == rankA) then
                    !No row was found that has an entry in the needed column --> quit with error
                    det_A = 0.D0
                    errorflag = -3322
                    return
                end if
            end do
        end if

        Do i = j+1, rankA     !Loop over the rows of matrix A
            !For all rows i>j, make entry i,j 0 by subtracting row j
            if (dabs(MatA_Gauss(i,j)) > eps_zero) then
                factor = MatA_Gauss(i,j)/MatA_Gauss(j,j)
                do l = j, rankA
                    MatA_Gauss(i,l) = MatA_Gauss(i,l) - factor * MatA_Gauss(j,l)
                end do
                vectorb(i) = vectorb(i) - factor * vectorb(j)
            end if
        end do
    end do


    !Calculate the determinant
    Det_A = 1.D0
    Do i = rankA, 1, -1
        !Calculate determinant
        Det_A = Det_A * MatA_Gauss(i,i)
    end do

    !Calculate the solution of the system of linear equations, if the determinant is not equal (close to) 0
    if (abs(Det_A) > eps_zero) then
        Do i = rankA, 1, -1
            !Calculate the solution x(i)
            vectorx(i) = vectorb(i)/MatA_Gauss(i,i)
            do k= i+1, rankA
                vectorx(i) = vectorx(i) - MatA_Gauss(i,k) * vectorx(k) / MatA_Gauss(i,i)
            end do
        end do
    end if

    end subroutine Gauss_algorithm_cosmo
    !**************************************************************************

    !**************************************************************************
    SUBROUTINE LUdecomp_cosmo (n,aMatrix,cMatrix,ierr,herr)
    !**************************************************************************
    ! subroutine to solve a linear system of quations
    !                               X*A = C
    ! for the unknown vector X
    ! written by E.Lemmon
    ! adapted by J. Gernert, Jan. 2011
    !--------------------------------------------------------------------------
    ! Variables:
    !   n       - Number of components in the mixture
    !   aMatrix - the matrix A
    !   cMatrix - the vector C
    !--------------------------------------------------------------------------
    implicit none
    !implicit double precision (a-h,o-z)
    !implicit integer (i-n)
    !max number of components in mixture
    double precision, dimension(160, 160):: amatrix
    double precision, dimension(160):: cmatrix,  ctemp, sdecomp
    integer, dimension(160)::iord
    character(255):: herr
    integer:: ierr, i, j, k, n
    double precision:: sum
    !dimension i,j,k,sum
    ierr=0
    do i=1,n
        iord(i)=i
        sdecomp(i)=Abs(aMatrix(1,i))
        do j=2,n
            If (Abs(aMatrix(j,i)) > sdecomp(i)) sdecomp(i)=Abs(aMatrix(j,i))
        end do
        !Singular matrix
        If (sdecomp(i) == 0) Then
            ierr=1
            herr='Singular matrix'
            RETURN
        end if
    end do

    j=1
    Call Pivot_cosmo(n,j,iord,aMatrix,sdecomp)
    If (aMatrix(1,iord(1)) == 0) Then
        ierr=1
        RETURN
    end if

    do j=2,n
        aMatrix(j,iord(1))=aMatrix(j,iord(1))/aMatrix(1,iord(1))
    end do
    do j=2,n-1
        do i=j,n
            sum=0
            do k=1,j-1
                sum=sum+aMatrix(k,iord(i))*aMatrix(j,iord(k))
            end do
            aMatrix(j,iord(i))=aMatrix(j,iord(i))-sum
        end do
        Call Pivot_cosmo(n,j,iord,aMatrix,sdecomp)
        do k=j+1,n
            sum=0
            do i=1,j-1
                sum=sum+aMatrix(i,iord(j))*aMatrix(k,iord(i))
            end do
            If (aMatrix(j,iord(j)) == 0) aMatrix(j,iord(j))=1E+20
            aMatrix(k,iord(j))=(aMatrix(k,iord(j))-sum)/aMatrix(j,iord(j))
        end do
    end do
    sum=0
    do k=1,n-1
        sum=sum+aMatrix(k,iord(n))*aMatrix(n,iord(k))
    end do
    aMatrix(n,iord(n))=aMatrix(n,iord(n))-sum
    If (aMatrix(n,iord(n)) == 0) aMatrix(n,iord(n))=1E+20

    cMatrix(iord(1))=cMatrix(iord(1))/aMatrix(1,iord(1))
    do i=2,n
        sum=0
        do j=1,i-1
            sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        end do
        cMatrix(iord(i))=(cMatrix(iord(i))-sum)/aMatrix(i,iord(i))
    end do

    do i=n-1,1,-1
        sum=0
        do j=i+1,n
            sum=sum+aMatrix(j,iord(i))*cMatrix(iord(j))
        end do
        cMatrix(iord(i))=cMatrix(iord(i))-sum
    end do
    do i=1,n
        ctemp(i)=cMatrix(iord(i))
    end do
    do i=1,n
        cMatrix(i)=ctemp(i)
    end do

    END subroutine LUdecomp_cosmo
    !**************************************************************************

    !**************************************************************************
    SUBROUTINE Pivot_cosmo (n,j,iord,aMatrix,sdecomp)
    !**************************************************************************
    implicit double precision (a-h,o-z)
    implicit integer (i-n)

    double precision, dimension(160,160)::amatrix
    integer, dimension(160):: iord
    double precision, dimension(160):: sdecomp


    ipivt=j
    big=Abs(aMatrix(j,iord(j))/sdecomp(iord(j)))
    do i=j+1,n
        dummy=Abs(aMatrix(j,iord(i))/sdecomp(iord(i)))
        If (dummy > big) Then
            big=dummy
            ipivt=i
        end if
    end do
    idummy=iord(ipivt)
    iord(ipivt)=iord(j)
    iord(j)=idummy

    END SUBROUTINE Pivot_cosmo

    !*****************************************************************************************************
    subroutine Succ_Sub_cosmo1(gl, Temp, delta_w, R_cosmo, sigma_profile, NR, seggamma)
    !*****************************************************************************************************


    implicit none

    type(type_gl) :: gl

    double precision :: Temp, R_cosmo
    double precision, dimension(51) :: sigma_profile
    double precision, dimension(51,51) :: delta_w
    logical :: NR
    integer :: m,o,i,j,k,l
    integer :: iter_counter
    double precision, dimension(51) :: seggamma, seggamma_old
    double precision :: sum
    double precision, dimension(51) :: converg
    double precision :: conv_crit
    double precision, dimension(51,51) :: AA
    logical :: molescheck
    integer :: start_sigma, end_sigma
    integer :: errval

    !convergence criteria
    conv_crit = 1.d-8

    !start and end values for occupied sigma-profile intervals
    start_sigma = gl%cosmo%interval_min(1)
    end_sigma = gl%cosmo%interval_max(1)

    !Iteration for segment activity coefficient
    !sigma_profile(k) * exp(-delta_w(j,k)/(r_gas*Temp)) can be calculated before iterating; only seggamma_old() is changing during iteration
    do j = start_sigma, end_sigma
        do k = start_sigma, end_sigma
            AA(j,k) =  sigma_profile(k) * exp(-delta_w(j,k)/(R_cosmo*Temp))
        end do
    end do

    !If temperature is very close to the previous temperature, use the previous segment activity coefficients as initial estimate
    molescheck = .false.
    do i = 1, gl%ncomp
        if (dabs(gl%cosmo%molefraction_cosmo_prev(i) - gl.molfractions(i)) < 1.D-2) then
        else
            molescheck = .false.
            exit
        end if
    end do

    if ((dabs(gl%cosmo%temp_cosmo_prev - Temp) < 1.D-2) .and. (molescheck == .true.)) then
        seggamma = gl%cosmo%seggamma_prev_v1
    else
        !Use standard initial values of 1, May 2018, changed to 1.1 due to occurance of solving the equations with negative seggamma
        seggamma = 1.0D0

    end if

    converg = 0.D0
    iter_counter = 0
    do
        iter_counter = iter_counter + 1
        seggamma_old = seggamma
        do j = start_sigma, end_sigma
            sum= 0.D0
            do k = start_sigma, end_sigma
                sum = sum + seggamma_old(k) * AA(j,k)
            end do
            seggamma(j) = 1 / sum
            seggamma(j) = (seggamma(j)+seggamma_old(j))/2.D0
        end do
        do j = start_sigma, end_sigma
            converg(j)=abs((seggamma(j)-seggamma_old(j))/seggamma_old(j))
        end do
        if (maxval(converg(start_sigma:end_sigma)) <=conv_crit) exit
        !if there is no convergence
        if (iter_counter == 1000) then
            errval = -7899
            exit
        end if
        if ((NR == .true.) .and. (iter_counter == 2)) exit
    end do

    end subroutine
    !*****************************************************************************************************

    !*****************************************************************************************************
    subroutine Succ_Sub_cosmo3(gl, Temp, delta_w3, R_cosmo, sigma_profile3, NR, seggamma_prev3, seggamma3, errval)
    !*****************************************************************************************************


    implicit none

    type(type_gl) :: gl

    double precision :: Temp, R_cosmo
    double precision, dimension(51,3) :: sigma_profile3
    double precision, dimension(51,51,3,3) :: delta_w3
    logical :: NR
    integer :: m,o,i,j,k,l
    integer :: iter_counter
    double precision, dimension(51,3) :: seggamma3, seggamma_old3, seggamma_prev3
    double precision :: sum
    double precision, dimension(51,3) :: converg3
    double precision :: conv_crit
    double precision, dimension(:,:,:,:), allocatable :: AA3
    integer, dimension(3) :: start_sigma, end_sigma
    logical :: molescheck
    integer :: errval

    if(.not. allocated(AA3)) allocate(AA3(51,51,3,3))

    !convergence criteria
    conv_crit = 1.d-8

    !start and end values for occupied sigma-profile intervals
    start_sigma = gl%cosmo%interval_min
    end_sigma = gl%cosmo%interval_max

    !Iteration for segment activity coefficient
    !sigma_profile_mix_v23(j,l) * exp(-delta_w_v23(i,j,k,l) / (R_cosmo*Temp)) can be calculated before iterating; only seggamma_old() is changing during iteration
    do k = 1, 3
        if (start_sigma(k) == 0) cycle
        do i = start_sigma(k), end_sigma(k)
            do l = 1, 3
                if (start_sigma(l) == 0) then
                    cycle
                end if
                do j = start_sigma(l), end_sigma(l)
                    AA3(i,j,k,l) = sigma_profile3(j,l) * exp(-delta_w3(i,j,k,l) / (R_cosmo*Temp))
                end do
            end do
        end do
    end do

    !If temperature is very close to the previous temperature, use the previous segment activity coefficients as initial estimate
    molescheck = .false.
    do i = 1, gl%ncomp
        if (dabs(gl%cosmo%molefraction_cosmo_prev(i) - gl%molfractions(i)) < 1.D-2) then
        else
            molescheck = .false.
            exit
        end if
    end do

    if ((dabs(gl%cosmo%temp_cosmo_prev - Temp) < 1.D-1) .and. (molescheck == .true.)) then
        seggamma3 = seggamma_prev3
    else
        !Use standard initial values of 1, May 2018, changed to 1.1 due to occurance of solving the equations with negative seggamma
        seggamma3 = 1.0D0

    end if

    converg3 = 0.D0
    iter_counter = 0
    do
        iter_counter = iter_counter + 1
        seggamma_old3 = seggamma3
        do k = 1, 3
            if (start_sigma(k) == 0) then
                cycle
            end if
            do i = start_sigma(k), end_sigma(k)
                sum= 0.D0
                do l = 1, 3
                    if (start_sigma(l) == 0) then
                        cycle
                    end if
                    do j = start_sigma(l), end_sigma(l)
                        sum = sum + seggamma_old3(j,l) * AA3(i,j,k,l)
                    end do
                end do
                seggamma3(i,k) = 1.D0 / sum
                seggamma3(i,k) = (seggamma3(i,k) + seggamma_old3(i,k))/2.D0
            end do
            do i=start_sigma(k), end_sigma(k)
                converg3(i,k)=abs((seggamma3(i,k) - seggamma_old3(i,k)) / seggamma_old3(i,k))
            end do
        end do
        if (maxval(converg3) <=conv_crit) exit
        !if iteration does not converge
        if (iter_counter == 1000) then
            !müssen noch aus der routine übergeben werden
            errval = -7899
            exit
        end if
        if ((NR == .true.) .and. (iter_counter == 2)) exit
    end do





    end subroutine
    !*****************************************************************************************************

    !*****************************************************************************************************
    subroutine Newton_Raphson_cosmo1(gl, Temp, delta_w, R_cosmo, sigma_profile, seggamma)
    !*****************************************************************************************************


    implicit none

    type(type_gl) :: gl

    double precision :: Temp, R_cosmo
    double precision, dimension(51) :: sigma_profile
    double precision, dimension(51,51) :: delta_w
    integer :: m,n,o,i,j,k,l
    integer :: iter_counter
    double precision, dimension(51) :: seggamma, seggamma_old
    double precision :: sum
    double precision, dimension(51) :: converg
    double precision :: conv_crit
    double precision, dimension(51,51) :: AA
    logical :: molescheck
    integer :: start_sigma, end_sigma
    integer :: errval
    logical :: speedup, NR

    double precision, dimension(51) :: Fi0_mix
    double precision, dimension(51,51) :: dFi_mix
    integer :: rankJ
    double precision, dimension(160) :: F_vector
    double precision, dimension(:,:), allocatable :: Jacobi_cosmo
    character(255) :: herr

    if(.not. allocated(Jacobi_cosmo)) allocate(Jacobi_cosmo(160,160))


    !speed up
    speedup = .true.

    !convergence criteria
    conv_crit = 1.d-8

    !start and end values for occupied sigma-profile intervals
    start_sigma = gl%cosmo%interval_min(1)
    end_sigma = gl%cosmo%interval_max(1)

    !Iteration for segment activity coefficient
    !sigma_profile(k) * exp(-delta_w(j,k)/(r_gas*Temp)) can be calculated before iterating; only seggamma_old() is changing during iteration
    do j = start_sigma, end_sigma
        do k = start_sigma, end_sigma
            AA(j,k) =  sigma_profile(k) * exp(-delta_w(j,k)/(R_cosmo*Temp))
        end do
    end do

    !If temperature is very close to the previous temperature, use the previous segment activity coefficients as initial estimate
    molescheck = .false.
    do i = 1, gl%ncomp
        if (dabs(gl%cosmo%molefraction_cosmo_prev(i) - gl.molfractions(i)) < 1.D-2) then
        else
            molescheck = .false.
            exit
        end if
    end do

    if (speedup == .false.) then
        if ((dabs(gl%cosmo%temp_cosmo_prev - Temp) < 1.D-2) .and. (molescheck == .true.)) then
            seggamma = gl%cosmo%seggamma_prev_v1
        else
            !Use standard initial values of 1, May 2018, changed to 1.1 due to occurance of solving the equations with negative seggamma
            seggamma = 1.0D0
        end if
    else
        NR = .true.
        call Succ_Sub_cosmo1(gl, Temp, delta_w, R_cosmo, sigma_profile, NR, seggamma)
    end if

    converg = 0.D0
    iter_counter = 0
    Fi0_mix = 0.D0
    !seggamma = 1.0D0
    dFi_mix = 0.D0
    do o = 1, 30    !Do 30 iterations maximum, there needs to be an errorval when it happens
        seggamma_old = seggamma
        !seggamma_old = 1.0D0
        m = 0
        n = 0
        F_vector = 0.D0
        Jacobi_cosmo = 0.D0
        do i = start_sigma, end_sigma
            sum = 0.D0
            m = m + 1
            n = 0
            do j = start_sigma, end_sigma
                !sum = sum + sigma_profile_mix(j) * seggamma_old(j) * exp(-delta_w(i,j) / (r_gas*Temp))
                sum = sum + AA(i,j) * seggamma_old(j)
            end do
            !Vector Fi0 for Taylor approximation
            Fi0_mix(i) = seggamma_old(i) * sum - 1.D0
            F_vector(m) = - Fi0_mix(i)
            !derivatives of Fi after each seggamma results in two different equations. One
            do j = start_sigma, end_sigma
                n = n + 1
                if (i == j) then
                    !dFi_mix(i,j) = sum + seggamma_old(i) * sigma_profile_mix(i) * exp(-delta_w(i,j) / (r_gas*Temp))
                    dFi_mix(i,j) = sum + seggamma_old(i) * AA(i,j)
                    Jacobi_cosmo(n,m) = dFi_mix(i,j)
                else
                    !dFi_mix(i,j) = seggamma_old(i) * sigma_profile_mix(j) * exp(-delta_w(i,j) / (r_gas*Temp))
                    dFi_mix(i,j) = seggamma_old(i) * AA(i,j)
                    Jacobi_cosmo(n,m) = dFi_mix(i,j)
                end if
            end do
        end do
        j = 0
        !To use Ludecomp_cosmo function F_vector and Jacobi_cosmo matrix have to be defined from Fi0_mix and dFi_mix vectors with suitable dimensions
        !Jacobi_cosmo = 0.D0
        !F_vector = 0.D0
        !do i = start_sigma(1), end_sigma(1)
        !    j = j + 1
        !    F_vector(j) = - Fi0_mix(i)
        !end do
        !m = 0
        !do i = start_sigma(1), end_sigma(1)
        !    m = m + 1
        !    n = 0
        !    do j = start_sigma(1), end_sigma(1)
        !        n = n + 1
        !        Jacobi_cosmo(n,m) = dFi_mix(i,j)
        !    end do
        !end do
        !RankJ has to be recalculated for all the actual occupied intervals of the sigma profile
        RankJ = end_sigma - start_sigma + 1
        !For RankJ start_sigma(k) was not yet accounted for
        call LUdecomp_cosmo (rankJ,Jacobi_cosmo,F_vector,errval,herr)
        j = 0
        do i = start_sigma, end_sigma
            j = j + 1
            seggamma(i) = F_vector(j) + seggamma(i)
        end do
        !Stop Iteration when Residual is very low
        if (maxval(dabs(Fi0_mix)) < conv_crit) then
            exit
        else
            if (o == 30) then
                errval = -12000 !check whether this value is already being used)
                !dFi_mix will be needed for analytical derivatives, needs to be set to 0 to avoid using wrong values
                dFi_mix = 0.D0
                !for next iterations of this system of compounds set solver to successive substitution
                gl%cosmo%solv_method = .false.
            end if
        end if

    end do
    !There is also the possibility of getting convergence, but the seggammas being negative. In that case the successive substitution method has to be used
    if (minval(seggamma) < -1.D-14) then
        !errval = -12000
        !use successive substitution methode in case of negative seggammas
        NR = .false.
        call Succ_Sub_cosmo1(gl, Temp, delta_w, R_cosmo, sigma_profile, NR, seggamma)
        !dFi_mix will be needed for analytical derivatives, needs to be set to 0 to avoid using wrong values
        dFi_mix = 0.D0
        !for next iterations of this system of compounds set solver to successive substitution
        gl%cosmo%solv_method = .false.
    end if





    end subroutine
    !*****************************************************************************************************

    !*****************************************************************************************************
    subroutine Newton_Raphson_cosmo3(gl, Temp, delta_w3, R_cosmo, sigma_profile3, seggamma_prev3, seggamma3)
    !*****************************************************************************************************


    implicit none

    type(type_gl) :: gl

    double precision :: Temp, R_cosmo
    double precision, dimension(51,3) :: sigma_profile3
    double precision, dimension(51,51,3,3) :: delta_w3
    integer :: m,n,o,i,j,k,l
    integer :: iter_counter
    double precision, dimension(51,3) :: seggamma3, seggamma_old3, seggamma_prev3
    double precision :: sum
    double precision, dimension(51,3) :: converg3
    double precision :: conv_crit
    double precision, dimension(:,:,:,:), allocatable :: AA3
    integer, dimension(3) :: start_sigma, end_sigma
    logical :: molescheck
    integer :: errval

    double precision, dimension(51,3) :: Fi0_mix
    double precision, dimension(51,51,3,3) :: dFi_mix
    integer :: rankJ
    double precision, dimension(160) :: F_vector
    double precision, dimension(:,:), allocatable :: Jacobi_cosmo
    character(255) :: herr
    logical :: speedup, NR

    if(.not. allocated(Jacobi_cosmo)) allocate(Jacobi_cosmo(160,160))
    if(.not. allocated(AA3)) allocate(AA3(51,51,3,3))

    !convergence criteria
    !nach test wieder auf 1.d-8 stellen
    conv_crit = 1.d-8

    !start and end values for occupied sigma-profile intervals
    start_sigma = gl%cosmo%interval_min
    end_sigma = gl%cosmo%interval_max

    !speed up
    speedup = .true.

    !Iteration for segment activity coefficient
    !sigma_profile_mix_v23(j,l) * exp(-delta_w_v23(i,j,k,l) / (R_cosmo*Temp)) can be calculated before iterating; only seggamma_old() is changing during iteration
    do k = 1, 3
        if (start_sigma(k) == 0) cycle
        do i = start_sigma(k), end_sigma(k)
            do l = 1, 3
                if (start_sigma(l) == 0) then
                    cycle
                end if
                do j = start_sigma(l), end_sigma(l)
                    AA3(i,j,k,l) = sigma_profile3(j,l) * exp(-delta_w3(i,j,k,l) / (R_cosmo*Temp))
                end do
            end do
        end do
    end do

    !If temperature is very close to the previous temperature, use the previous segment activity coefficients as initial estimate
    molescheck = .true.
    do i = 1, gl%ncomp
        if (dabs(gl%cosmo%molefraction_cosmo_prev(i) - gl.molfractions(i)) < 1.D-2) then
        else
            molescheck = .false.
            exit
        end if
    end do

    if (speedup == .false.) then
        if ((dabs(gl%cosmo%temp_cosmo_prev - Temp) < 1.D-1) .and. (molescheck == .true.)) then
            seggamma3 = seggamma_prev3
        else
            !Use standard initial values of 1, May 2018, changed to 1.1 due to occurance of solving the equations with negative seggamma
            seggamma3 = 1.0D0
        end if
    else
        NR = .true.
        call Succ_Sub_cosmo3(gl, Temp, delta_w3, R_cosmo, sigma_profile3, NR, seggamma_prev3, seggamma3, errval)
    end if

    converg3 = 0.D0
    iter_counter = 0
    Fi0_mix = 0.D0
    !seggamma3 = 1.0D0
    dFi_mix = 0.D0
    do o = 1, 30    !Do 30 iterations maximum, there needs to be an errorval when it happens
        seggamma_old3 = seggamma3
        !seggamma_old = 1.0D0
        m = 0
        n = 0
        F_vector = 0.D0
        Jacobi_cosmo = 0.D0
        do k = 1, 3
            if (start_sigma(k) == 0) cycle
            do i = start_sigma(k), end_sigma(k)
                sum = 0.D0
                m = m + 1
                n = 0
                do l = 1, 3
                    if (start_sigma(l) == 0) cycle
                    do j = start_sigma(l), end_sigma(l)
                        !sum = sum + sigma_profile_mix(j) * seggamma_old(j) * exp(-delta_w(i,j) / (r_gas*Temp))
                        sum = sum + AA3(i,j,k,l) * seggamma_old3(j,l)
                    end do
                end do
                !Vector Fi0 for Taylor approximation
                Fi0_mix(i,k) = seggamma_old3(i,k) * sum - 1.D0
                F_vector(m) = - Fi0_mix(i,k)
                !derivatives of Fi after each seggamma results in two different equations. One
                do l = 1, 3
                    if (start_sigma(l) == 0) cycle
                    do j = start_sigma(l), end_sigma(l)
                        n = n + 1
                        if ((i == j) .and. (k == l)) then
                            !dFi_mix(i,j) = sum + seggamma_old(i) * sigma_profile_mix(i) * exp(-delta_w(i,j) / (r_gas*Temp))
                            dFi_mix(i,j,k,l) = sum + seggamma_old3(i,l) * AA3(i,j,k,l)
                            Jacobi_cosmo(n,m) = dFi_mix(i,j,k,l)
                        else
                            !dFi_mix(i,j) = seggamma_old(i) * sigma_profile_mix(j) * exp(-delta_w(i,j) / (r_gas*Temp))
                            dFi_mix(i,j,k,l) = seggamma_old3(i,k) * AA3(i,j,k,l)
                            Jacobi_cosmo(n,m) = dFi_mix(i,j,k,l)
                        end if
                    end do
                end do
            end do
        end do
        !j = 0
        !To use Ludecomp_cosmo function F_vector and Jacobi_cosmo matrix have to be defined from Fi0_mix and dFi_mix vectors with suitable dimensions
        !Jacobi_cosmo = 0.D0
        !F_vector = 0.D0
        !do i = start_sigma(1), end_sigma(1)
        !    j = j + 1
        !    F_vector(j) = - Fi0_mix(i)
        !end do
        !m = 0
        !do i = start_sigma(1), end_sigma(1)
        !    m = m + 1
        !    n = 0
        !    do j = start_sigma(1), end_sigma(1)
        !        n = n + 1
        !        Jacobi_cosmo(n,m) = dFi_mix(i,j)
        !    end do
        !end do
        !RankJ has to be recalculated for all the actual occupied intervals of the sigma profile
        RankJ = end_sigma(1) - start_sigma(1) + end_sigma(2) - start_sigma(2) + end_sigma(3) - start_sigma(3)
        !For RankJ begin_sigma_mix(k) was not yet accounted for
        do k = 1,3
            if (start_sigma(k) > 0.D0) then
                RankJ = RankJ + 1
            end if
        end do

        call LUdecomp_cosmo (rankJ,Jacobi_cosmo,F_vector,errval,herr)
        j = 0
        do k = 1, 3
            if (start_sigma(k) == 0) cycle
            do i = start_sigma(k), end_sigma(k)
                j = j + 1
                seggamma3(i,k) = F_vector(j) + seggamma3(i,k)
            end do
        end do
        !Stop Iteration when Residual is very low
        if (maxval(dabs(Fi0_mix)) < conv_crit) then
            exit
        else
            if (o == 30) then
                errval = -12000 !check whether this value is already being used)
                !dFi_mix will be needed for analytical derivatives, needs to be set to 0 to avoid using wrong values
                dFi_mix = 0.D0
                !for next iterations of this system of compounds set solver to successive substitution
                !nach TEST wieder einkommentieren
                gl%cosmo%solv_method = .false.
            end if
        end if

    end do
    !There is also the possibility of getting convergence, but the seggammas being negative. In that case the successive substitution method has to be used
    if (minval(seggamma3) < -1.D-14) then
        !errval = -12000
        !use successive substitution methode in case of negative seggammas
        !nach test wieder einkommentieren
        NR = .false.
        call Succ_Sub_cosmo3(gl, Temp, delta_w3, R_cosmo, sigma_profile3, NR, seggamma_prev3, seggamma3, errval)
        !dFi_mix will be needed for analytical derivatives, needs to be set to 0 to avoid using wrong values
        dFi_mix = 0.D0
        !for next iterations of this system of compounds set solver to successive substitution
        gl%cosmo%solv_method = .false.
    end if





    end subroutine
    !*****************************************************************************************************

    !*****************************************************************************************************
    !END OF COSMO-SAC CODE
    !*****************************************************************************************************




    end module fnrderivs_module
