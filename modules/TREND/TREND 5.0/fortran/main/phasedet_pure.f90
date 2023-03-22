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
    !Cite as: Span, R.; Beckmuller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.;
    !          Jager, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020):
    !          TREND. Thermodynamic Reference and Engineering Data 5.0.
    !          Lehrstuhl fur Thermodynamik, Ruhr-Universitat Bochum.

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

    ! module for file phasedet_pure.f90
    module phasedet_pure_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use rhomix_pt_module
    use flash_pure_module
    use ancillary_equations_mix_module


    contains



    !************************************************************************************
    subroutine PhaseDet_pure(gl,press, Temp, rho, d_vap, d_liq, phasetype, vapfrac, nrofphases, errval)
    !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "PhaseDet_pure" :: PhaseDet_pure
    
    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   rho         - density [mol/m設
    !   Temp        - Temperature [K]
    !
    ! OUTPUT:
    !   press       - Pressure [MPa]
    !   d_vap       - Vapor phase density [mol/m設
    !   d_liq       - liquid phase density [mol/m設
    !   phasetype   - phase identifier
    !   vapfrac        - Molar vapor fraction
    !   nrofphases  - number of phases
    !   errval      - Error value
    !               - -2211: VLE-Algorithm failed
    !************************************************************************************
    !Theresa Wiens Dez. 2011








    implicit none

    type(type_gl) :: gl


    double precision :: press, temp, rho
    double precision :: d_liq, d_vap
    double precision :: vapfrac
    integer :: phasetype
    integer :: nrofphases
    integer :: errval
    !double precision :: p_calc
    double precision :: psat,pmelt
    integer :: errPsat
    integer :: nrsubst, iter

    d_liq = 0.d0
    d_vap = 0.d0
    Psat = 0.d0
    nrofphases = 0
    errval = 0
    errPsat = 0
    nrsubst = 1
    phasetype = 0

    press = p_calc(gl,Temp, rho, 1)
    if ((Temp .ge. (gl%ttp(nrsubst) - 1.D-6)) .and. (Temp < gl%tc(nrsubst))) then                   !can be twophase, vapor or liquid
        call Flash_Pure_PhaseBoundary(gl,Psat, Temp, d_vap, d_liq, 1, errPsat, iter, 1)
        !call VLEpure (Temp, Psat, d_vap, d_liq, errPsat, nrsubst)
        if (errPsat == 0) then                                                                                  !successful
            if ((d_liq - rho) > -1e-14 .and. (rho - d_vap) > -1e-14) then                                                         !twophase point
                nrofphases = 2
                phasetype = 0
                vapfrac = (1.d0/rho - 1.d0/d_liq) / (1.d0/d_vap - 1.d0/d_liq)
                press = psat                                                                                    !no need to calculate the pressur again
            elseif (rho > d_liq) then                                                                          !liquid
                phasetype = 1
                nrofphases = 1
                vapfrac = -2.d0
            elseif (rho < d_vap) then                                                                          !vapor
                phasetype = 2
                nrofphases = 1
                vapfrac = -2.d0
            end if
        else                                                                                                    !not successful
            errval = -2211
        end if
    elseif ((Temp < gl%tc(nrsubst)) .and. (press > gl%pc(nrsubst))) then                                              !liquid-like supercrit
        phasetype = 3
        nrofphases = 1
        vapfrac = -3.d0
    elseif ((Temp > gl%tc(nrsubst)) .and. (press < gl%pc(nrsubst))) then                                              !vapor-like supercrit
        phasetype = 4
        nrofphases = 1
        vapfrac = -3.d0
    elseif ((Temp > gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                            !vapor phase                                                                                                                         !supercrit
        phasetype = 2
        nrofphases = 1
        vapfrac = -1.d0
    elseif ((Temp < gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                            !maybe solid phase -- >  safety first                                                                                                                       !supercrit
        phasetype = 6
        nrofphases = 1
        vapfrac = -4.d0
    else                                                                                                        !supercrit
        phasetype = 5
        nrofphases = 1
        vapfrac = -3.d0
    end if

    if ((phasetype == 0) .or. (phasetype == 1) .or. (phasetype == 3) .or. (phasetype == 5)) then
        pmelt = pmelt_eq(gl,Temp,nrsubst)
        if (trim(gl%components(nrsubst)) == "water" .and. Temp < 273.31d0) then                                              !solid phase?
            if (press < pmelt .or. press > gl%pmelt_high) then
                phasetype = 6
                nrofphases = 1
                vapfrac = -4.d0
                !DANGER
            end if
        else if (press > pmelt) then
            if ((temp > gl%ttp(nrsubst)) .and. (pmelt > 1.d-12)) then    !Monika, 23/04/2018, this is only reasonable if the melting pressure is greater than 0 and if it is not exactly at the triple point
                phasetype = 6
                vapfrac = -4.d0
            end if
            nrofphases = 1
            !DANGER
        end if
    end if

    gl%phase_id(1) = phasetype

    end subroutine



    !************************************************************************************
    subroutine PhaseDet_pure_tp (gl,press, Temp, rho, d_vap, d_liq, phasetype, vapfrac, nrofphases, errval)
    !DEC$ ATTRIBUTES DLLEXPORT :: PhaseDet_pure_tp

    !************************************************************************************
    !--------------------------------------------------------------------------
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   Temp        - Temperature [K]
    !
    ! OUTPUT:
    !   rho         - density [mol/m設
    !   d_vap       - Vapor phase density [mol/m設
    !   d_liq       - liquid phase density [mol/m設
    !   phasetype   - phase identifier
    !   vapfrac        - Molar vapor fraction
    !   nrofphases  - number of phases
    !   errval      - Error value
    !               - -2211: VLE-Algorithm failed
    !               - -8877: density iteration failed
    !************************************************************************************
    !Theresa Wiens Dez. 2011





    implicit none

    type(type_gl) :: gl


    double precision :: press, temp
    double precision :: rho, d_vap, d_liq
    double precision :: vapfrac
    integer :: phasetype !1: liquid / 2: vapor
    integer :: nrofphases
    integer :: errval
    double precision :: psat,pmelt
    integer :: errPsat
    integer :: nrsubst
    integer :: errorflag, iter

    d_vap = 0.d0
    d_liq = 0.d0
    iter = 0
    Psat = 0.d0
    nrofphases = 0
    errval = 0
    errPsat = 0
    nrsubst = 1
    errorflag = 0
    phasetype = 0
    rho = 0.D0


    if ((Temp > gl%ttp(nrsubst)) .and. (Temp < gl%tc(nrsubst)) .and. (press > gl%ptp(nrsubst)) .and. (press < gl%pc(nrsubst))) then          !can be twophase, vapor or liquid
        call Flash_Pure_PhaseBoundary(gl,Psat,Temp,d_vap,d_liq, 1, errPsat, iter, 1)
        !call VLEpure(Temp,Psat,d_vap,d_liq,errPsat,nrsubst)
        if (errPsat == 0) then                                                                                                      !succesfull
            if (dabs(psat-press)/Psat*100.d0 < 1.d-7) then                                                                          !comparision of absolute, relative, procentual deviation
                nrofphases = 2                                                                                                      !T,p combination twophase
                phasetype = 0
                vapfrac = -1.d0                                                                                                     !twophase
                rho = -1.d0                                                                                                         !impossible to determine this density
                errval = -898968
            elseif (press > psat) then                                                                                              !liquid
                phasetype = 1
                nrofphases = 1
                vapfrac = -2.d0
                rho = rhomix_calc(gl,Temp, press, 0.d0, 1, nrsubst)
                !Andreas June 2014
                if (dabs(rho) < 1.D-14) errval = -8888
            elseif (press < psat) then                                                                                              !vapor
                phasetype = 2
                nrofphases = 1
                vapfrac = -2.d0
                rho = rhomix_calc(gl,Temp, press, 0.d0, 2, nrsubst)
                !Andreas June 2014
                if (dabs(rho) < 1.D-14) errval = -8888
            end if
        else                                                                                                                        !not successful
            errval = -2211
        end if
    elseif ((Temp < gl%tc(nrsubst)) .and. (press > gl%pc(nrsubst))) then                                                                  !liquid-like supercrit
        phasetype = 3
        nrofphases = 1
        vapfrac = -3.d0
        rho = rhomix_calc(gl,Temp, press, 0.d0, 0, nrsubst)
        !Andreas June 2014
        if (dabs(rho) < 1.D-14) errval = -8888
        !if ((errorflag /= 0) .and. (errorflag /= 3)) then
        !    rho = -4444.d0
        !    errval = -8877
        !end if
    elseif ((Temp > gl%tc(nrsubst)) .and. (press < gl%pc(nrsubst))) then                                                                  !vapor-like supercrit
        phasetype = 4
        nrofphases = 1
        vapfrac = -3.d0
        rho = rhomix_calc(gl,Temp, press, 0.d0, 0, nrsubst)
        !Andreas June 2014
        if (dabs(rho) < 1.D-14) errval = -8888
        !if ((errorflag /= 0) .and. (errorflag /= 3)) then
        !    rho = -4444.d0
        !    errval = -8877
        !end if
    elseif ((Temp > gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                               !vapor phase                                                                                                                         !supercrit
        phasetype = 2
        nrofphases = 1
        vapfrac = -1.d0
        rho = rhomix_calc(gl,Temp, press, 0.d0, 0, nrsubst)
        !Andreas June 2014
        if (dabs(rho) < 1.D-14) errval = -8888
        !if ((errorflag /= 0) .and. (errorflag /= 3)) then
        !    rho = -4444.d0
        !    errval = -8877
        !end if
    elseif ((Temp < gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                                !maybe solid phase -- >  safety first                                                                                                                       !supercrit
        phasetype = 62
        nrofphases = 1
        vapfrac = -4.d0
        rho = rhomix_calc(gl,Temp, press, 0.d0, 0, nrsubst)
        !Andreas June 2014
        if (dabs(rho) < 1.D-14) errval = -8888
        !if ((errorflag /= 0) .and. (errorflag /= 3)) then
        !    rho = -4444.d0
        !    errval = -8877
        !end if
    else                                                                                                                            !supercrit
        phasetype = 5
        nrofphases = 1
        vapfrac = -3.d0
        rho = rhomix_calc(gl,Temp, press, 0.d0, 0, nrsubst)
        !Andreas June 2014
        if (dabs(rho) < 1.D-14) errval = -8888
        !if ((errorflag /= 0) .and. (errorflag /= 3)) then
        !    rho = -4444.d0
        !    errval = -8877
        !end if
    end if

    !Andreas Jun 2014
    if (Temp < gl%pmeltmintemp(nrsubst)) then
        phasetype = 6
        nrofphases = 1
        vapfrac = -4.d0
    else    !Special case for water
        if (trim(gl%components(nrsubst)) == "water" .and. Temp < 273.31d0) then                                              !solid state
            pmelt = pmelt_eq(gl,Temp,nrsubst)
            if (pmelt > 0.d0) then ! pmelt successfully calculated
                if (press < pmelt .or. press > gl%pmelt_high) then
                    phasetype = 6
                    nrofphases = 1
                    vapfrac = -4.d0
                    !DANGER
                end if
            end if
        end if
    end if

    gl%phase_id(1) = phasetype

    end subroutine



    !************************************************************************************
    subroutine PhaseDet_ps_pure(gl,press, Temp, rho, d_vap, d_liq, phasetype, vapfrac, s_spec, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified fluid
    ! For this calculation the pressure and the entropy are given
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   s_spec      - Entropy [J / mol K]
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - homogeneous or mixture density [mol/m設
    !   d_vap       - Vapor phase density [mol/m設
    !   d_liq       - Liquid phase density [mol/m設
    !   phasetype   - Phase identifier
    !   vapfrac     - Molar vapor fraction
    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************
    !Theresa Wiens Dez. 2011








    implicit none

    type(type_gl) :: gl


    double precision :: press, s_spec, vapfrac, Temp
    double precision :: d_vap, d_liq, rho
    integer :: phasetype
    integer :: nrofphases
    integer :: errval

    integer :: errTsat
    integer :: nrsubst

    !  double precision :: S_CALC, P_CALC
    double precision :: Tsat
    double precision :: s_vap
    double precision :: s_liq


    ! needed to transmit T, p to Regula Falsi:
    double precision :: t_min
    double precision :: t_max
    double precision :: t_min_allowed
    double precision :: t_max_allowed
    double Precision :: DeltaS_allowed,pmelt
    type(type_additional_parameters) :: parameters
    integer :: max_iterations
    integer :: iterations
    integer :: errorflag , iter

    nrofphases = 0
    errval = 0
    errTsat = 0
    nrsubst = 1  ! this routine is exclusively used for pure fluids
    Tsat = 0.d0
    s_vap = 0.d0
    s_liq = 0.d0
    DeltaS_allowed = 1.d-8
    !parameters = 0.d0
    max_iterations = 50
    iterations = 0
    errorflag = 0
    phasetype = 0
    rho = 0.D0


    !default boundaries
    t_min = gl%tminfluid(nrsubst)
    t_max = gl%tmaxfluid(nrsubst)
    t_min_allowed = gl%tminfluid(nrsubst)
    t_max_allowed = gl%tmaxfluid(nrsubst)

    parameters%a_p(1) = press
    parameters%a_p(2) = s_spec


    if ((press < gl%pc(nrsubst)) .and. (press > gl%ptp(nrsubst))) then                                                         !can be twophase, vapor, or liquid
        call Flash_Pure_PhaseBoundary(gl,press,Tsat,d_vap,d_liq, 2, errTsat, iter, 1)
        !call VLEpurePres(press,Tsat,d_vap,d_liq, errTsat,nrsubst)
        if (errTsat == 0) then                                                                                              !successful
            s_vap = S_CALC(gl,Tsat, d_vap, nrsubst)
            s_liq = S_CALC(gl,Tsat, d_liq, nrsubst)
            if ((s_liq <= s_spec).and.(s_spec <= s_vap)) then                                                               !twophase point
                nrofphases = 2
                phasetype = 0
                vapfrac = (s_spec-s_liq)/(s_vap-s_liq)
                rho = 1.d0 / (1.d0/d_liq + vapfrac *(1.d0/d_vap - 1.d0/d_liq))
                Temp = Tsat
            else
                if (s_spec < s_liq) then                                                                                    !liquid
                    vapfrac = 0.d0
                    nrofphases = 1
                    phasetype = 1
                    t_max = tsat
                    t_max_allowed = tsat
                elseif (s_spec > s_vap) then                                                                                !vapor
                    vapfrac = 1.d0
                    nrofphases = 1
                    phasetype = 2
                    t_min = tsat
                    t_min_allowed = tsat
                end if

                d_vap = 0.d0
                d_liq = 0.d0

                !determine temperature from p and s
                call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, Errorflag, Parameters)

                if ((Errorflag == 0) .or. (Errorflag == 3)) then  !successful
                    gl%p_it = p_calc(gl,temp, gl%rho_it, nrsubst) ! control value checking whether the two iterations were successful
                    !density already determined from t and p
                    if (dabs(press-gl%p_it)/press*100.d0 < 1.d-5) then
                        rho = gl%rho_it
                    else ! rho_it doesn't belong to this ps-configuration, for whatever reason
                        errval = -1234
                    end if
                else    ! iteration failed
                    rho = -4444.d0
                    errval = -8877
                end if
            end if
        else                                                                                                                !not successful
            errval = -2211
        end if
    else                                                                                                                    !can be supercrit, liquid-like supercrit or vapor-like supercrit
        !determine temperature from p and s
        call Regula_Falsi(gl,entropy_diff, temp, t_min, t_max, DeltaS_allowed, t_min_allowed, t_max_allowed, &
            &                       Max_iterations, Iterations, Errorflag, Parameters)

        if ((Errorflag == 0) .or. (Errorflag == 3)) then                                                                !successful
            if ((Temp < gl%tc(nrsubst)) .and. (press > gl%pc(nrsubst))) then                                                      !liquid-like supercrit
                phasetype = 3
                nrofphases = 1
                vapfrac = -3.d0
            elseif ((Temp > gl%tc(nrsubst)) .and. (press < gl%pc(nrsubst))) then                                                  !vapor-like supercrit
                phasetype = 4
                nrofphases = 1
                vapfrac = -3.d0
            elseif ((Temp > gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                !vapor phase                                                                                                                         !supercrit
                phasetype = 2
                nrofphases = 1
                vapfrac = -1.d0
            elseif ((Temp < gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                !maybe solid phase -- >  safety first                                                                                                                       !supercrit
                phasetype = 62
                nrofphases = 1
                vapfrac = -4.d0
            else                                                                                                            !supercrit
                phasetype = 5
                nrofphases = 1
                vapfrac = -3.d0
            end if
            gl%p_it = p_calc(gl,temp, gl%rho_it, nrsubst) !control value checking whether the two iterations were successful
            !density already determined from t and p
            if (dabs(press-gl%p_it)/press*100.d0 < 1.d-6) then
                rho = gl%rho_it
            else !rho_it doesn't belong to this ps-configuration, for whatever reason
                errval = -1234
            end if
        else    !iteration failed
            rho = -4444.d0
            errval = -8877
        end if
    end if

    if ((phasetype == 1) .or. (phasetype == 3) .or. (phasetype == 5).or. (phasetype == 62)) then
        pmelt = pmelt_eq(gl,Temp,nrsubst)
        if (trim(gl%components(nrsubst)) == "water" .and. Temp < 273.31d0) then                                              !solid state
            if (press < pmelt .or. press > gl%pmelt_high) then
                phasetype = 6
                nrofphases = 1
                vapfrac = -4.d0
                !DANGER
            end if
        else if (press > pmelt) then
            phasetype = 6
            nrofphases = 1
            vapfrac = -4.d0
            !DANGER
        end if
    end if

    gl%phase_id(1) = phasetype

    end subroutine



    Double Precision Function entropy_diff(gl,temp, Parameters)



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: temp                                       ! Inputvariablen
    Double Precision :: p, S, S_TP
    ! Double Precision :: S_CALC_TP                                ! function allocating calculated pressure
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    S = parameters%a_p(2)
    S_TP = S_CALC_TP(gl,temp, p)
    if (dabs(S_TP + 44444.d0) <= 1.d-12) then
        Parameters%a_p(65)= S_TP
    end if
    entropy_diff = S-S_TP ! using internal function

    End Function



    !************************************************************************************
    subroutine PhaseDet_ph_pure(gl,press, Temp, rho, d_vap, d_liq, phasetype, vapfrac, h_spec, nrofphases, errval)
    !************************************************************************************
    ! Subroutine for determining how many phases are present for the specified fluid
    ! For this calculation the pressure and the enthalpy are given
    ! Variables:
    ! INPUT:
    !   press       - Pressure [MPa]
    !   h_spec      - Enthalpy [J / mol]
    !
    ! OUTPUT:
    !   Temp        - Temperature [K]
    !   rho         - homogeneous or mixture density [mol/m設
    !   d_vap       - Vapor phase density [mol/m設
    !   d_liq       - Liquid phase density [mol/m設
    !   phasetype   - Phase identifier
    !   vapfrac     - Molar vapor fraction
    !   nrofphases  - number of phases
    !   errval      - Error value
    !************************************************************************************
    !Theresa Wiens Dez. 2011








    implicit none

    type(type_gl) :: gl


    double precision :: press, h_spec, vapfrac, Temp
    double precision :: d_vap, d_liq, rho
    integer :: phasetype
    integer :: nrofphases
    integer :: errval

    integer :: errTsat
    integer :: nrsubst

    !double precision :: H_CALC, P_CALC
    double precision :: Tsat
    double precision :: h_vap
    double precision :: h_liq


    ! needed to transmit T, p to Regula Falsi:
    double precision :: t_min
    double precision :: t_max
    double precision :: t_min_allowed
    double precision :: t_max_allowed
    double Precision :: DeltaH_allowed,pmelt
    type(type_additional_parameters) :: parameters
    integer :: max_iterations
    integer :: iterations
    integer :: errorflag, iter

    nrofphases = 0
    errval = 0
    errTsat = 0
    nrsubst = 1  ! this routine is exclusively used for pure fluids
    Tsat = 0.d0
    h_vap = 0.d0
    h_liq = 0.d0
    DeltaH_allowed = 1.d-8
    max_iterations = 200
    iterations = 0
    errorflag = 0
    !parameters = type_additional_parameters()
    phasetype = 0
    rho = 0.D0


    !default boundaries
    t_min = gl%tminfluid(nrsubst)
    t_max = gl%tmaxfluid(nrsubst)
    t_min_allowed = gl%tminfluid(nrsubst)
    t_max_allowed = gl%tmaxfluid(nrsubst)

    parameters%a_p(1) = press
    parameters%a_p(2) = h_spec

    if ((press < gl%pc(nrsubst)) .and. (press > gl%ptpmod(nrsubst))) then                                                         !can be twophase, vapor, or liquid
        call Flash_Pure_PhaseBoundary(gl,press,Tsat,d_vap,d_liq, 2, errTsat, iter, 1)
        !call VLEpurePres(press,Tsat,d_vap,d_liq, errTsat,nrsubst)
        if (errTsat == 0) then                                                                                              !successful
            h_vap = H_CALC(gl,Tsat, d_vap, nrsubst)
            h_liq = H_CALC(gl,Tsat, d_liq, nrsubst)
            if ((h_liq <= h_spec).and.(h_spec <= h_vap)) then                                                               !twophase point
                nrofphases = 2
                phasetype = 0
                vapfrac = (h_spec-h_liq)/(h_vap-h_liq)
                rho = 1.d0 / (1.d0/d_liq + vapfrac *(1.d0/d_vap - 1.d0/d_liq))
                Temp = Tsat
            else
                if (h_spec < h_liq) then                                                                                    !liquid
                    vapfrac = 0.d0
                    nrofphases = 1
                    phasetype = 1
                    t_max = tsat
                    t_max_allowed = tsat
                elseif (h_spec > h_vap) then                                                                                !vapor
                    vapfrac = 1.d0
                    nrofphases = 1
                    phasetype = 2
                    t_min = tsat
                    t_min_allowed = tsat
                end if

                d_vap = 0.d0
                d_liq = 0.d0

                !determine temperature from p and h
                call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
                    &                       Max_iterations, Iterations, Errorflag, Parameters)

                if ((Errorflag == 0) .or. (Errorflag == 3)) then  !successful
                    gl%p_it = p_calc(gl,temp, gl%rho_it, nrsubst) ! control value checking whether the two iterations were successful
                    !density already determined from t and p
                    if (dabs(press-gl%p_it)/press*100.d0 < 1.d-4) then
                        rho = gl%rho_it
                    else ! rho_it doesn't belong to this ph-configuration, for whatever reason
                        errval = -1234
                    end if
                else    ! iteration failed
                    rho = -4444.d0
                    errval = -1235
                end if
            end if
        else                                                                                                                !not successful
            errval = -2211
        end if
    else                                                                                                                    !can be supercrit, liquid-like supercrit or vapor-like supercrit
        !determine temperature from p and h
        call Regula_Falsi(gl,enthalpy_diff, temp, t_min, t_max, DeltaH_allowed, t_min_allowed, t_max_allowed, &
            &                       Max_iterations, Iterations, Errorflag, Parameters)

        if ((Errorflag == 0) .or. (Errorflag == 3)) then                                                                !successful
            if ((Temp < gl%tc(nrsubst)) .and. (press > gl%pc(nrsubst))) then                                                      !liquid-like supercrit
                phasetype = 3
                nrofphases = 1
                vapfrac = -3.d0
            elseif ((Temp > gl%tc(nrsubst)) .and. (press < gl%pc(nrsubst))) then                                                  !vapor-like supercrit
                phasetype = 4
                nrofphases = 1
                vapfrac = -3.d0
            elseif ((Temp > gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                !vapor phase                                                                                                                         !supercrit
                phasetype = 2
                nrofphases = 1
                vapfrac = -1.d0
            elseif ((Temp < gl%ttp(nrsubst)) .and. (press < gl%ptp(nrsubst))) then                                                !maybe solid phase -- >  safety first                                                                                                                       !supercrit
                phasetype = 62
                nrofphases = 1
                vapfrac = -4.d0
            else                                                                                                            !supercrit
                phasetype = 5
                nrofphases = 1
                vapfrac = -3.d0
            end if
            gl%p_it = p_calc(gl,temp, gl%rho_it, nrsubst) !control value checking whether the two iterations were successful
            !density already determined from t and p
            if (dabs(press-gl%p_it)/press*100.d0 < 1.d-4) then
                rho = gl%rho_it
            else !rho_it doesn't belong to this ph-configuration, for whatever reason
                errval = -1234
            end if
        else    !iteration failed
            rho = -4444.d0
            errval = -8877
        end if
    end if

    if ((phasetype == 1) .or. (phasetype == 3) .or. (phasetype == 5).or. (phasetype == 62)) then
        pmelt = pmelt_eq(gl,Temp,nrsubst)
        if (trim(gl%components(nrsubst)) == "water" .and. Temp < 273.31d0) then                                              !solid state
            if (press < pmelt .or. press > gl%pmelt_high) then
                phasetype = 6
                nrofphases = 1
                vapfrac = -4.d0
                !DANGER
            end if
        else if (press > pmelt) then
            phasetype = 6
            nrofphases = 1
            vapfrac = -4.d0
            !DANGER
        end if
    end if

    gl%phase_id(1) = phasetype

    end subroutine



    Double Precision Function enthalpy_diff(gl,temp, Parameters)



    implicit none

    type(type_gl) :: gl


    ! Variable declaration:
    !  --------------------------------------------------
    Double Precision :: temp                                    ! Inputvariablen
    Double Precision :: p, H, H_TP
    !Double Precision :: H_CALC_TP                                ! function allocating calculated pressure
    type(type_additional_parameters) :: parameters
    !  --------------------------------------------------

    p = parameters%a_p(1)
    H = parameters%a_p(2)
    H_TP = H_CALC_TP(gl,temp, p)
    if (dabs(H_TP + 44444.d0) <= 1.d-12) then
        Parameters%a_p(65)= H_TP
    end if
    enthalpy_diff = H-H_CALC_TP(gl,temp, p) ! using internal function

    End Function






    end module phasedet_pure_module
