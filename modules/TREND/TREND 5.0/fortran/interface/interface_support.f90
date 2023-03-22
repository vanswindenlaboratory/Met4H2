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


    ! module for file interface_support.f90
    module interface_support_module
    !global use inclusion
    use module_all_types
    use module_trend_info
    use calc_functions
    use dryice_module
    use waterice_module
    use seawater_module
    use transport_module
    use hdrt_chem_pot_module
    use phasedet_mix_module
    use phasedet_pure_module
    use phasedet_sol_module
    use phasedet_sol_pure_module
    use costald_module
    use electrolytes
    use gibbsderivs_module

    contains



    subroutine  TREND_RELEASE_CHECK_STDCALL (version)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "TREND_RELEASE_CHECK_STDCALL" :: TREND_RELEASE_CHECK_STDCALL
    !DEC$ ATTRIBUTES stdcall, reference :: TREND_RELEASE_CHECK_STDCALL
    implicit none



    !define input variables
    character (255), intent(out) :: version

    call TREND_RELEASE_CHECK (version)

    end subroutine

    subroutine TREND_RELEASE_CHECK(version)


    character (255), intent(out) :: version
    !DEC$ ATTRIBUTES DLLEXPORT :: TREND_RELEASE_CHECK
    version = Version_parameter

    end subroutine TREND_RELEASE_CHECK



    !subroutine moles_incheck(gl,moles,  errorflag)
    !!Subroutine to check validity of input compositions (for example of FLASH_EXPERT or FLASH_NC_3P_EXPERT)
    !!convert if neccessary
    !
    !implicit none
    !
    !type(type_gl) :: gl
    !
    !
    !double precision, dimension(30) :: moles, molesave
    !integer:: l, e, j, t, i, c, b, k
    !integer:: errorflag, converttype
    !double precision:: sum
    !character(12) :: dd
    !double precision :: sum_weight, wm_mix
    !double precision, dimension (30) :: x_molar
    !
    !errorflag = 0
    !i=0
    !
    !! check if all entries are within 0 < x_i < 1 and fill the vector to the size of 30 with zeros
    !sum = 0.d0
    !
    !do k = 1, gl%ncomp
    !    if (((0.d0 >= moles(k)) .OR. (moles(k) > 1.d0)) .AND. (k == 1)) then    !check if 0 < xi <= 1
    !        errorflag = -9951   ! there is a valid molefraction given in the first entry
    !        return
    !    end if
    !    if (((0.d0 > moles(k)) .OR. (moles(k) >= 1.d0)) .AND. (k > 1)) then    !check if 0 < xi < 1
    !        errorflag = -9951
    !        return
    !    end if
    !
    !    if (moles(k) == 0.d0) then  ! check whether there are any further molefractions given -- >  error
    !        do c = k, 30
    !            if (moles(c) /= 0.d0) then
    !                errorflag = -9951
    !                return
    !            end if
    !        end do
    !        exit
    !    end if
    !    sum = sum + moles(k)
    !end do
    !
    !if (abs(sum - 1.d0) > 1.d-10) then  !check if sum /= 1
    !    errorflag = -9952
    !    return
    !end if
    !
    !!if(gl%unitin == 2) then
    !!
    !!    molesave = gl%molfractions
    !!    gl%molfractions = moles
    !!    converttype = 1
    !!    call wm_mix_calc(gl, wm_mix)
    !!    call convert_fractions(gl, converttype, wm_mix, moles)
    !!
    !!    gl%molfractions = molesave
    !!
    !!end if
    !
    !
    !
    !End subroutine moles_incheck




    !**************************************************************************************************
    !Core-Routine that calculates the requested values
    !**************************************************************************************************
    double precision function PROP_EOS (gl,prop, input, prop1, prop2, fluids, moles, EOS_indicator, mixtype, path, errorflag)






    implicit none

    type(type_gl) :: gl



    !define input variables
    integer, intent(in) :: prop
    !   prop    property
    !    1       T
    !    2       D
    !    3       P
    !    4       U
    !    5       H
    !    6       S
    !    7       G
    !    8       A
    !    9       CP
    !   10       CV
    !   11       WS
    !   12       B
    !   13       C
    !   14       CP0
    !   15       Q
    !   16       Z
    !   17       PNUMER
    !   18       WSNUMER
    !   19       A00
    !   20       A01
    !   21       A02
    !   22       A03
    !   23       A10
    !   24       A11
    !   25       A12
    !   26       A20
    !   27       A21
    !   28       A30
    !   29       UR
    !   30       HR
    !   31       SR
    !   32       CVR
    !   33       CPOTR
    !   34       GRUEN
    !   35       PIP
    !   36       DE
    !   37       ST
    !   38       VISDYN
    !   39       TCX
    !   40       VOLEXP
    !   41       DVIR
    !   42       RIEM
    !   43       COMPT
    !   44       COMPS
    !   45       THROT
    !   46       EXPANS
    !   47       EXPANT
    !   48       DUDV
    !   49       JTCO
    !   50       DPDT
    !   51       VE
    !   52       HE
    !   53       GE
    !   54       B12
    !   55       GAMMAGD
    !   56       Residual pressure
    !   57       dBdT (derivative of the second virial coeffcient)
    !   58       dCdT (derivative of the third virial coeffcient)
    !   59       dPdD (1st derivative of the pressure with respect to density at constant temperature)
    !   60       d2PdD2 (2nd derivative of the pressure with respect to density at constant temperature)
    !   61       d2PdDdT (1st mixed derivatives of the pressure with respect to density and temperature)
    !   62       d2PdT2 (2nd derivative of the pressure with respect to temperature at constant density)
    !   63       dDdT (1st derivative of the density with respect to temperature at constant pressure)
    !   64       DSPIN (Calculation of Spinodals)
    !   65       PHASEEOS (Phase Determination)
    !   66       AF Acentric factor
    !   67       Density maximum of water and heavy water
    !   68       Boyle-Temperature
    !   69       Joule-Thomson-Temperature
    !   70       Joule-Thomson-Inversion-Temperature
    !   71       First PIP derivative with respect to density (at constant temperature) (numerical derivative)
    !   72       Second PIP derivative with respect to density^2 (at constant temperature) (numerical derivative)
    !   73       Second derivative of second virial coeff. B with respect to T^2
    !   74       Second derivative of third virial coeff. C with respect to T^2
    !   75       Second derivative of fourth virial coeff. D with respect to T^2
    !   76       Third derivative of second virial coeff. B with respect to T^3
    !   77       Third derivative of third virial coeff. C with respect to T^3
    !   78       Third derivative of fourth virial coeff. D with respect to T^3
    !   79       Fourth derivative of second virial coeff. B with respect to T^4
    !   80       Fourth derivative of third virial coeff. C with respect to T^4
    !   81       Fourth derivative of fourth virial coeff. D with respect to T^4

    character (12), intent(inout) :: input
    double precision, intent(inout) :: prop1, prop2
    character (30), dimension(30), intent(inout) :: fluids
    double precision, dimension(30), intent(inout) :: moles
    integer, dimension(30), intent(inout) :: EOS_indicator
    integer, intent(inout) :: mixtype
    character (255), intent(inout) :: path

    !double precision :: U_CALC, H_CALC, S_CALC, G_CALC, A_CALC, CP_CALC, CV_CALC, WS_CALC, B_CALC, C_CALC, CP0_CALC, P_CALC, Z_CALC
    !double precision :: A00_CALC, A01_CALC, A02_CALC, A03_CALC, A10_CALC, A11_CALC, A12_CALC, A20_CALC, A21_CALC, A30_CALC, D_calc
    !double precision :: UR_CALC,HR_CALC,SR_CALC,CVR_CALC,CPOTR_CALC,GRUEN_CALC,PIP_CALC,DE_CALC,ST_CALC,VISDYN_CALC,TCX_CALC,VOLEXP_CALC
    !double precision :: RIEM_CALC, COMPT_CALC, COMPS_CALC, THROT_CALC, EXPANS_CALC, EXPANT_CALC, DUDV_CALC, JTCO_CALC, DPDT_calc, GAMMAGD_calc
    !double precision :: PR_calc, DBDT_calc, DCDT_calc, LNG_density, DPDD_CALC, D2PDD2_CALC, D2PDDT_CALC, D2PDTT_CALC, DDDT_CALC, DSPIN_CALC
    !define variables for calculating
    double precision :: t, d, p, h, s
    double precision :: dvap
    double precision :: dliq
    double precision :: propvap
    double precision :: propliq
    double precision :: x_liq(30)
    double precision :: x_vap(30)

    !double precision, dimension(30) :: molev

    integer :: errorflag
    integer :: nrsubst

    double precision:: rho(5)           ! vector containing the densities of all phases
    double precision:: x_Phase(30,5)    ! vector conaining the compositions of all phases
    double precision:: vapfrac             ! molar vapor fraction
    integer:: nrofphases                ! at this point: only 1 or 2 implemented
    !Indicate which phases are present
    integer:: phasetype(5)              !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    integer :: ILIMITS
    integer :: i,j

    double precision, dimension(5) :: phasefrac
    double precision, dimension(5) :: prop_phase, dens_reac

    !variable for spinodals
    integer :: spinprop

    !Variables needed for hydrate enthalpy and entropy calculation
    double precision :: d_phase, d_fluid, ChemPot_hyd, fug_CO2, prop_hyd
    double precision, dimension(30) :: x_fluid, x_hyd, chem_pot, fug_g, Chempot
    integer:: pos_DryIce

    double precision :: wm_mix, p1,p2,d1,d2,D_EOS,DP_DRHO_S,dens_phase,dens, salb

    logical:: costald_eq_type
    logical:: RKM_eq_type, ERKM_eq_type
    double precision :: rhoredmixorg, tredmixorg, PROP_MIX_IDEAL, PROP_PURE, BVIR1, BVIR2
    double precision :: ttrip, ptrip, pressb, pressf, addpart, rmix
    double precision, dimension(30) :: chempot_num_vap, chempot_num_liq
    double precision, dimension(200) :: tb, pb
    integer :: phase

    PROP_MIX_IDEAL = 0.d0
    dens = 0.d0

    !phase variables
    dvap = 0.d0
    dliq = 0.d0
    propvap = 0.d0
    propliq = 0.d0
    vapfrac = 0.d0
    x_liq = 0.d0
    x_vap = 0.d0
    errorflag = 0
    ILIMITS = 0
    phasefrac = 0.D0
    !molev = 0.d0
    prop_hyd = 0.d0
    t = 0.d0; p = 0.d0; d = 0.d0
    x_phase = 0.d0
    rho = 0.d0
    prop_phase = 0.d0


    call setup(gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, mixtype, errorflag)
    PROP_EOS = 0.d0
    !irregular input, or problems with fluid file or reference state
    if (errorflag /= 0) then
        PROP_EOS = errorflag
        input=gl%inptorig
        return
    end if

    call reduced_parameters_calc(gl,300.d0)


    !special case for density calculation with RKM, COSTALD, or enhanced RKM equation of state (eq_type 8, 9, or 81)
    if ((gl%Mix_type == 8) .or. (gl%Mix_type == 9) .or. (gl%Mix_type == 81)) then
        PROP_EOS = LNG_density(gl,input, prop1, prop2, prop,errorflag)
        input=gl%inptorig
        return
    end if

    !special case for spinodals: dspin not implemented for mixtures, yet
    if ((prop == 64) .and. (gl%ncomp > 1)) then
        PROP_EOS = -5668.d0
        errorflag = PROP_EOS
        input=gl%inptorig
        return
    end if

    !special case for density maximum: input is always pressure
    if (prop == 67) then
        if ((gl%ncomp > 1) .or. (trim(input) .ne. 'tp')) then
            PROP_EOS = -5660.d0
            errorflag = PROP_EOS
            input=gl%inptorig
            return
        else
            PROP_EOS = TDENSMAX_calc(gl,prop1,1)
            return
        end if
    end if

    !test if model for dielectric constant exists
    if (prop == 36) then
        if(.not.allocated(gl%de)) then
            PROP_EOS = -5234
            errorflag = -5234
            input=gl%inptorig
            return
        else
            if (.not.gl%de%de_read) then
                PROP_EOS = -5234
                errorflag = -5234
                input=gl%inptorig
                return
            end if
        end if
    end if

    !test if model for surface tension exists, different to other transport properties because no individual subtype exists for the surface tension, yet
    if ((prop == 37) .and. (.not.gl%stn_read)) then
        if ((gl%ncomp == 1) .and. (gl%substcasnr(1)=="7722-84-1"))  then   !hydrogen peroxide, special case implemented for Vaclav Vins, Monika 2019/02/13
            PROP_EOS = ST_CALC(gl,prop1,1)
            input=gl%inptorig
            return
        else
            PROP_EOS = -5234
            errorflag = PROP_EOS
            input=gl%inptorig
            return
        end if
    end if

    !test if model for viscosity exists
    if (prop == 38) then
        if (.not.gl%visco%eta_read) then
            if ((gl%ncomp == 1) .and. (gl%substcasnr(1)=="7722-84-1"))  then   !hydrogen peroxide, special case implemented for Vaclav Vins, Monika 2019/02/13
                PROP_EOS = VISDYN_CALC(gl,prop1,prop2,1)
                input=gl%inptorig
                return
            else
                PROP_EOS = -5234
                errorflag = PROP_EOS
                input=gl%inptorig
                return
            end if
        end if
    end if

    !test if model for thermal conductivity exists
    if (prop == 39) then
        if (.not.gl%tcx%tcx_read) then
            if ((gl%ncomp == 1) .and. (gl%substcasnr(1)=="7722-84-1"))  then   !hydrogen peroxide, special case implemented for Vaclav Vins, Monika 2019/02/13
                PROP_EOS = TCX_CALC(gl,prop1,prop2,1)
                input=gl%inptorig
                return
            else
                PROP_EOS = -5234
                errorflag = PROP_EOS
                input=gl%inptorig
                return
            end if
        end if
    end if


    !If the SRK, PR, LKP, or PC-SAFT is used check whether a model for the ideal gas exists
    !Andreas November 2013
    !Only neccessary for caloric properties
    if(((prop > 3).and. (prop < 12)).or. (prop == 14).or. (prop == 18) .or. (prop == 34).or. (prop == 35).or. (prop == 38).or. (prop == 39) .or. (prop == 42) .or. (prop == 49) .or. (prop == 55)) then
        Do i = 1, gl%ncomp
            if ((gl%Eq_type(i) > 1) .and. (gl%Eq_type(i) .ne. 7))  then  !Eq_type 7 => AGA8
                !---------------------------------------------------------
                !If constants B and C for the SRK cp0 model are 0, no model is implemented
                if (((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8))) then
                    errorflag = -2908
                    PROP_EOS = errorflag
                    input=gl%inptorig
                    return
                end if
            end if
        end do
    end if
    PROP_EOS = 0.d0


    if (.not.(((prop ==12) .or. (prop ==13) .or. (prop ==14) .or. (prop == 17) .or. (prop ==41) .or. (prop ==57) .or. (prop ==58)) .and. (input(1:1) == 't'))) then !input of these properties is only T so that no phase determination is necessary
        !Handle different input parameters and do all calculations
        call inpt_handle (gl,input, prop1, prop2, p, t, d, dvap, dliq, moles, rho, x_Phase, phasetype, &
            & phasefrac, nrofphases, nrsubst, errorflag)


    elseif ((prop == 17)) then !for PNUMER_EOS there are no phase determination calculations needed
        nrofphases = 0
        if (input(1:1) == 't') then
            t = prop1
            d = 1.d0
            dens = 1.d0
        end if
    else
        nrofphases=1
        rho(1)=1.d0
        phasetype(1)=1
        x_phase(:,phasetype(1)) = gl%molfractions
        phasefrac(1) = 1.d0
        t=prop1
        if (gl%ncomp >1) then
            nrsubst = 0
        else
            nrsubst = 1
        end if
    end if


    if (prop == 68) then
        PROP_EOS = TBoyle_calc(gl,nrsubst)
        if(PROP_EOS.lt.1d-12) errorflag = PROP_EOS
        input=gl%inptorig
        return
    end if

    if (prop == 69) then
        PROP_EOS = TJT_calc(gl,nrsubst)
        if(PROP_EOS.lt.1d-12) errorflag = PROP_EOS
        input=gl%inptorig
        return
    end if

    if (prop == 70) then
        PROP_EOS = TJTINV_calc(gl,nrsubst)
        if(PROP_EOS.lt.1d-12) errorflag = PROP_EOS
        input=gl%inptorig
        return
    end if

    if (prop == 71) then
        PROP_EOS = D_PIP_DD_CALC(gl,T,D,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 72) then
        PROP_EOS = D2_PIP_D2D_CALC(gl,T,D,nrsubst)
        input=gl%inptorig
        return
    end if

    !vanessa
    if (prop == 73) then
        PROP_EOS = d2BdT2_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 74) then
        PROP_EOS = d2CdT2_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 75) then
        PROP_EOS = d2DdT2_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 76) then
        PROP_EOS = d3BdT3_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 77) then
        PROP_EOS = d3CdT3_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 78) then
        PROP_EOS = d3DdT3_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 79) then
        PROP_EOS = d4BdT4_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 80) then
        PROP_EOS = d4CdT4_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 81) then
        PROP_EOS = d4DdT4_CALC(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    !Numerical calculation of virial coeff (temporary Sven)
    if (prop == 82) then
        PROP_EOS = b_calc_num(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    !Numerical calculation of virial coeff (temporary Sven)
    if (prop == 83) then
        PROP_EOS = c_calc_num(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    !Numerical calculation of virial coeff (temporary Sven)
    if (prop == 84) then
        PROP_EOS = d_calc_num(gl,T,nrsubst)
        input=gl%inptorig
        return
    end if

    if (prop == 85) then
        PROP_EOS = n_eff_calc(gl,T,D,nrsubst)
        input=gl%inptorig
        return
    end if

    if (errorflag /= 0) then
        PROP_EOS = errorflag
        input=gl%inptorig
        return
    end if

    !check range of validity
    if (gl%hold_limits) then   !hold_limits = true => range of validity has to be checked, hold_limits = false => extrapolation to arbitrary state points
        if ((prop /= 36) .and. (prop /= 38) .and. (prop /= 39)) then  !check limits for dielectric constant, viscosity, and thermal conductivity
            if ((prop == 12) .or. (prop == 13) .or. (prop == 57) .or. (prop ==58)) d = 0.d0 !density for virials and derivatives is not calculated -> initialize d
            call check_limits (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)
        else
            call check_limits_trans (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, prop, ILIMITS)
        endif

        if (ILIMITS /= 0) then
            PROP_EOS = ILIMITS
            input=gl%inptorig
            return
        end if
    end if

    ! CP,CV,WS,GRUEN,DE,VISDYN,TCX only available in homogenious region
    if ((nrofphases > 1) .and. (.not. gl%phaseboundary) .and. ((prop == 9) .or. (prop == 10) .or. (prop == 11) &
        .or. (prop == 34) .or. (prop == 35) .or. (prop == 36) .or. (prop == 38) .or. (prop == 39) .or. (prop == 40) .or. (prop == 49) )) then
        PROP_EOS = -6666
        errorflag = PROP_EOS
        input=gl%inptorig
        return
    end if

    if (prop == 65) then
        if (phasetype(2) == 0) then  !single phase
            if (phasetype(1) == 1) then
                PROP_EOS = 2
            else if (phasetype(1) == 3) then
                PROP_EOS = 1
            else if (phasetype(1) == 4) then
                PROP_EOS = 6
            else
                PROP_EOS = -111
                errorflag = PROP_EOS
            end if
        else    !twophase
            PROP_EOS = 12
        end if
    else
        !Start of property calculation
        do i = 1, nrofphases

            !-----------------------------------------------------------------------------------------------------------------------
            ! Fluid phase
            gl%gecarrier = .false.
            if (phasetype(i) < 4) then
                if(((gl%seawater) .or. (gl%el_present)) .and. ((phasetype(i) .eq. 3) .or. (phasetype(i) .eq. 2)) ) gl%gecarrier = .true.
                if (gl%ncomp > 1) then !mixture
                    if ((prop == 37).or.(prop == 38).or.(prop == 39)) then    !ST,VISDYN,TCX: no mixture model available, yet
                        PROP_EOS = -5223
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                    !If mixture, recalculate tredmix and rhoredmix
                    if (gl%ncomp > 1) then
                        gl%molfractions = x_phase(:,phasetype(i))
                        call reduced_parameters_calc(gl,t) !Dummy temperature 300 K for the SRK
                    end if
                end if

                d_phase = rho(phasetype(i)) ! properties will be calculated with this phase specific density

                dens_phase = 1.d0/d_phase       ! needed for the summation of the overall density

                if (gl%gecarrier) then
                    dens_reac(phasetype(i)) = brine_dens(gl, t, p, d_phase)
                end if

                select case (prop)

                case (1) !Temperature
                    PROP_EOS = t

                    !input=inptorig
                    !return

                case (3) !Pressure

                    ! AGA8
                    if (gl%mix_type .eq. 7) then
                        PROP_EOS = P_calc(gl,t, d_phase, nrsubst)
                        ! AGA8: DETAIL EOS can be applied for pressures up to 280MPa only
                        if (PROP_EOS > 280) then
                            errorflag = -9999
                            PROP_EOS = errorflag
                            input=gl%inptorig
                            return
                        end if
                    else
                        PROP_EOS = p
                    end if


                    !input=inptorig
                    !return

                case (4) !Internal Energy
                    prop_phase(i) = U_calc(gl,t, d_phase, nrsubst)
                case (5) !Enthalpy
                    prop_phase(i) = H_calc(gl,t, d_phase, nrsubst)
                case (6) !Entropy
                    prop_phase(i) = S_calc(gl,t, d_phase, nrsubst)
                case (7) !Gibbs Free Energy
                    prop_phase(i) = G_calc(gl,t, d_phase, nrsubst)
                case (8) !Helmholtz Free Energy
                    prop_phase(i) = A_calc(gl,t, d_phase, nrsubst)
                case (9) !Isobaric Heat Capacity
                    prop_phase(i) = CP_calc(gl,t, d_phase, nrsubst)
                case (10) !Isochoric Heat Capacity
                    prop_phase(i) = CV_calc(gl,t, d_phase, nrsubst)
                case (11) !Speed of Sound
                    prop_phase(i) = WS_calc(gl,t, d_phase, nrsubst)
                case (12) !Second Virial Coefficient
                    !Saving causes a problem here: if the SRK is used to calculate B for a mixture, the parameters a and b are not updated. reduced_parameters_calc(gl) was not called.
                    !Same happens when using Helmholtz EOS! So, reduced_paramters_calc(gl) is called here. This is a priliminary solution which needs to be discussed!! Andreas, September 2015
                    call reduced_parameters_calc(gl,t)
                    prop_phase(i) = B_calc(gl,t, nrsubst)
                case (13) !Third Virial Coefficient
                    !Saving causes a problem here: if the SRK is used to calculate C for a mixture, the parameters a and b are not updated. reduced_parameters_calc(gl) was not called.
                    !Same happens when using Helmholtz EOS! So, reduced_paramters_calc(gl) is called here. This is a priliminary solution which needs to be discussed!! Andreas, September 2015
                    call reduced_parameters_calc(gl,t)
                    prop_phase(i) = C_calc(gl,t, nrsubst)
                case (14) !Ideal Gas Heat Capacity
                    prop_phase(i) = CP0_calc(gl,t, nrsubst)
                case (15) !Quality (Molar Vapor Fraction)
                    !only depending on temperature, so only one calculation necessary
                    !this way, equation boundaries will not be checked
                    if (nrofphases == 2) then
                        PROP_EOS = phasefrac(phasetype(i))
                    else
                        PROP_EOS = -6667.d0
                        errorflag = PROP_EOS
                    end if
                    input=gl%inptorig
                    return
                case (16) !Compressibility factor
                    prop_phase(i) = Z_calc(gl,t, d_phase, nrsubst)
                case (29) !Residual Internal Energy
                    prop_phase(i) = UR_calc(gl,t, d_phase, nrsubst)
                case (30) !Residual Enthalpy
                    prop_phase(i) = HR_calc(gl,t, d_phase, nrsubst)
                case (31) !Residual Entropy
                    prop_phase(i) = SR_calc(gl,t, d_phase, nrsubst)
                case (32) !Residual Isochoric Heat Capacity
                    prop_phase(i) = CVR_calc(gl,t, d_phase, nrsubst)
                case (33) !Residual Chemical Potential
                    prop_phase(i) = CPOTR_calc(gl,t, d_phase, nrsubst)
                case (34) !Grueneisen Coefficient
                    prop_phase(i) = GRUEN_calc(gl,t, d_phase, nrsubst)
                case (35) !PIP
                    prop_phase(i) = PIP_calc(gl,t, d_phase, nrsubst)
                case (36) !Dielectric Constant
                    prop_phase(i) = DE_calc(gl,t, d_phase, nrsubst)
                case (37) !Surface Tension
                    !only depending on temperature, so only one calculation necessary
                    !this way, equation boundaries will not be checked
                    if (nrofphases == 2) then
                        PROP_EOS = ST_calc(gl,t, nrsubst)
                    else
                        PROP_EOS = -5222.d0
                        errorflag = PROP_EOS
                    end if
                    input=gl%inptorig
                    return
                case (38) !Dynamic Viscosity
                    prop_phase(i) = VISDYN_calc(gl,t, d_phase, nrsubst)
                case (39) !Thermal Conductivity
                    prop_phase(i) = TCX_calc(gl,t, d_phase, nrsubst)
                case (40) !volume expansivity at constant pressure
                    prop_phase(i) = VOLEXP_calc(gl,t, d_phase, nrsubst)
                case (41) !Fourth Virial Coefficient
                    !Saving causes a problem here: if the SRK is used to calculate C for a mixture, the parameters a and b are not updated. reduced_parameters_calc(gl) was not called.
                    !Same happens when using Helmholtz EOS! So, reduced_paramters_calc(gl) is called here. This is a priliminary solution which needs to be discussed!! Andreas, September 2015
                    call reduced_parameters_calc(gl,t)
                    prop_phase(i) = D_calc(gl,t, nrsubst)
                case (42) !RIEM
                    if (gl%ncomp /= 1) then
                        PROP_EOS = -5223d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                    prop_phase(i) = RIEM_calc(gl,t, d_phase, nrsubst)
                case (43) !isothermal compressibility
                    prop_phase(i) = COMPT_calc(gl,t, d_phase, nrsubst)
                case (44) !isentropic compressibility
                    prop_phase(i) = COMPS_calc(gl,t, d_phase, nrsubst)
                case (45) !isothermal throttling coefficient
                    prop_phase(i) = THROT_calc(gl,t, d_phase, nrsubst)
                case (46) !isentropic expansion coefficient
                    prop_phase(i) = EXPANS_calc(gl,t, d_phase, nrsubst)
                case (47) !isothermal expansion coefficient
                    prop_phase(i) = EXPANT_calc(gl,t, d_phase, nrsubst)
                case (48) !derivative of the internal energy with respect to the volume at constant temperature
                    prop_phase(i) = DUDV_calc(gl,t, d_phase, nrsubst)
                case (49) !joule - thomson coefficient
                    prop_phase(i) = JTCO_calc(gl,t, d_phase, nrsubst)
                case (50) !dp/dT @ v=const
                    prop_phase(i) = DPDT_calc(gl,t, d_phase, nrsubst)
                case (51) !vE
                    if (gl%ncomp == 1) then
                        PROP_EOS = -5248d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                case (52) !hE
                    if (gl%ncomp == 1) then
                        PROP_EOS = -5248d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                    prop_phase(i) = H_calc(gl,t, d_phase, nrsubst)
                case (53) !gE
                    if (gl%ncomp == 1) then
                        PROP_EOS = -5248d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                    prop_phase(i) = G_calc(gl,t, d_phase, nrsubst)
                case(54) !Cross Virial Coefficient
                    !only depending on temperature, so only one calculation necessary
                    !this way, equation boundaries will not be checked
                    if (gl%ncomp /= 2) then
                        PROP_EOS = -5249d0
                        errorflag = PROP_EOS
                    else
                        !Saving causes a problem here: if the SRK is used to calculate C for a mixture, the parameters a and b are not updated. reduced_parameters_calc(gl) was not called.
                        !Same happens when using Helmholtz EOS! So, reduced_paramters_calc(gl) is called here. This is a priliminary solution which needs to be discussed!! Andreas, September 2015
                        call reduced_parameters_calc(gl,t)
                        PROP_EOS = B_calc(gl,t, nrsubst)
                        BVIR1 = B_calc(gl,t, 1)
                        BVIR2 = B_calc(gl,t, 2)
                        PROP_EOS = (PROP_EOS-gl%molfractions(1)*gl%molfractions(1)*BVIR1-gl%molfractions(2)*gl%molfractions(2)*BVIR2)/(2.d0*gl%molfractions(1)*gl%molfractions(2))
                    end if

                    input=gl%inptorig
                    return
                case (55) !fundamental derivative of gas dynamics
                    prop_phase(i) = GAMMAGD_calc(gl,t, d_phase, nrsubst)
                case (56) !residual pressure
                    prop_phase(i) = PR_calc(gl,t, d_phase, nrsubst)
                case (57) !derivative of the second virial coefficient with respect to the temperature
                    prop_phase(i) = DBDT_calc(gl,t, nrsubst)
                case (58) !derivative of the third virial coefficient with respect to the temperature
                    prop_phase(i) = DCDT_calc(gl,t, nrsubst)
                case (59) !dp/dT @ v=const
                    prop_phase(i) = DPDD_calc(gl,t, d_phase, nrsubst)
                case (60) !d2p/dD2 @ T=const
                    prop_phase(i) = D2PDD2_calc(gl,t, d_phase, nrsubst)
                case (61) !d2p/dDdT
                    prop_phase(i) = D2PDDT_CALC(gl,t, d_phase, nrsubst)
                case (62) !d2p/dT2
                    prop_phase(i) = D2PDTT_CALC(gl,t, d_phase, nrsubst)
                case (63) !dD/dT @p=const
                    prop_phase(i) = DDDT_CALC(gl,t, d_phase, nrsubst)
                case (64) !spinodalds
                    !Only 4 input codes valid
                    if ((input == trim('tliq'))) then
                        spinprop = 1
                    elseif ((input == trim('tvap'))) then
                        spinprop = 2
                    else
                        PROP_EOS = -5666.d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if

                    if (((spinprop == 1) .or. (spinprop == 2)) .and. (prop1 == gl%tc(1))) then
                        PROP_EOS=gl%rhoc(1)
                        input=gl%inptorig
                        return
                    elseif (((spinprop == 3) .or. (spinprop == 4)) .and. (prop1 == gl%pc(1))) then
                        PROP_EOS=gl%rhoc(1)
                        input=gl%inptorig
                        return
                    end if

                    PROP_EOS = DSPIN_CALC (gl,prop1,spinprop,1)
                    input=gl%inptorig
                    if ((PROP_EOS + 5667.d0) < 1.d-12) return
                    return

                    !case (65) not applicable for loop over phases
                    !case 65 is placed before the loop over phases

                case (66) !acentric factor
                    prop_phase(i) = OMEGA_CALC(gl,nrsubst)

                    !case 67,68,69,70 is placed before the loop over phases


                end select
                !
                !gl%seacalc = .false.
                !gl%gecarrier = .false.

                !-----------------------------------------------------------------------------------------------------------------------
                !Solid Phase
            elseif (phasetype(i) == 4) then

                if (gl%solidtype_akt_phase == 1) then         !Water
                    !Andreas Feb 2015, Get the molar volume of solid water
                    dens_phase = 1.D0 / rho(phasetype(i))
                    select case (prop)
                    case (4) !Internal Energy
                        prop_phase(i) = u_WaterIce(gl,t, p)
                    case (5) !Enthalpy
                        prop_phase(i) = h_WaterIce(gl,t, p)
                    case (6) !Entropy
                        prop_phase(i) = s_WaterIce(gl,t, p)
                    case (7) !Gibbs Energy
                        prop_phase(i) = g_WaterIce(gl,t, p)
                    case (8) !Helmholtz Energy
                        prop_phase(i) = f_WaterIce(gl,t, p)
                    case (9) !Isobaric Heat Capacity
                        prop_phase(i) = cp_WaterIce(gl,t, p)
                    case (10:) !e.g. Isochoric Heat Capacity and Speed of Sound
                        if (dabs(phasefrac(phasetype(i))) > 1.D-14) then
                            PROP_EOS = -9904d0
                            errorflag = PROP_EOS
                            input=gl%inptorig
                            return
                        else
                            prop_phase(i) = 0.D0  !cv, ws of solids cannot be calculated yet, but if e.g. prop of vapor in sublimation equilibrium is asked for, set prop_sol to 0.D0
                        end if
                    end select

                elseif (gl%solidtype_akt_phase == 2) then     !CO2
                    !Andreas Feb 2015, Get the molar volume of solid CO2
                    dens_phase = 1.D0 / rho(phasetype(i))
                    select case (prop)
                    case (4) !Internal Energy
                        prop_phase(i) = u_DryIce(gl,t, p)
                    case (5) !Enthalpy
                        prop_phase(i) = h_DryIce(gl,t, p)
                    case (6) !Entropy
                        prop_phase(i) = s_DryIce(gl,t, p)
                    case (7) !Gibbs Energy
                        prop_phase(i) = g_DryIce(gl,t, p)
                    case (8) !Helmholtz Energy
                        prop_phase(i) = f_DryIce(gl,t, p)
                    case (9) !Isobaric Heat Capacity
                        prop_phase(i) = cp_DryIce(gl,t, p)
                    case (10:) !e.g. Isochoric Heat Capacity and Speed od Sound
                        if (dabs(phasefrac(phasetype(i))) > 1.D-14) then
                            PROP_EOS = -9904d0
                            errorflag = PROP_EOS
                            input=gl%inptorig
                            return
                        else
                            prop_phase(i) = 0.D0  !cv, ws of solids cannot be calculated yet, but if e.g. prop of vapor in sublimation equilibrium is asked for, set prop_sol to 0.D0
                        end if
                    end select

                else                                !No solid equation available
                    !melting/sublimation pressure/temperature for pure fluids
                    !=> ancillary equations will only be evaluated in the fluid phase (phasetype(i) == 3);
                    !a solid equation is not available for non-hydrate components
                    if ((gl%ncomp ==1) .and. ((prop == 1) .or. (prop == 3))) then
                        !if a pure fluid is calculated, the sublimation pressure extrapolated from the vapor pressure curve will be returned;
                        !the boolean variable p_fake enables to check if it is an extrapolated or a real property
                        if (gl%p_fake) then
                            PROP_EOS = -5234d0
                            errorflag = PROP_EOS
                        end if
                        gl%p_fake = .false.
                        input=gl%inptorig
                        return
                    else
                        PROP_EOS = -9904d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if
                end if



                !-----------------------------------------------------------------------------------------------------------------------
                !Hydrate
            elseif (phasetype(i) == 5) then

                !Andreas Feb 2015, Get the molar volume of gas hydrate
                dens_phase = 1.D0 / rho(phasetype(i))

                select case (prop)  !no hydrates for e.g. cp,cv,ws
                case (9:)
                    PROP_EOS =-12900d0
                    errorflag = PROP_EOS
                    input=gl%inptorig
                    return
                end select

                !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
                !calculate the entropy and enthalpy of hydrates
                if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives

                    x_fluid = x_phase(:,phasetype(1))
                    d_fluid = rho(phasetype(1))
                    x_hyd = x_phase(:,phasetype(i))
                    call hdrt_ancillary_hs(gl,t, d_fluid, p*1.D6, x_fluid, chem_pot, fug_g, errorflag)
                    if (errorflag == 0) then
                        select case (prop)
                        case (4,5) !Internal Energy and Enthalpy
                            call hdrt_enthalpy(gl,t, p*1.D6, x_hyd, chem_pot, fug_g,  prop_hyd, errorflag) !prop_hyd = h_hyd
                        case (6) !Entropy
                            call hdrt_entropy(gl,t, p*1.D6, x_hyd, fug_g,  prop_hyd, errorflag) !prop_hyd = s_hyd
                        case (7,8) !Gibbs Energy and Helmholtz Energy
                            call hdrt_Gibbs_energy(gl,x_hyd, chem_pot, prop_hyd, errorflag) !prop_hyd = g_hyd
                        end select
                    else
                        input=gl%inptorig
                        return
                    end if

                else

                    if ((gl%solidtype_akt_phase == 2) .and. (gl%ncomp == 2)) then !dry ice in equilibrium with hydrate
                        pos_DryIce = gl%solidpos_akt_phase
                        fug_g(1) = fug_DryIce(gl,t,p, pos_dryIce) *1.D6
                        if (fug_g(1) < 1.D-12) then
                            errorflag = -7779
                        else
                            !pure hydrates
                            !fug_CO2 = fug_g(1)
                            chem_pot(2) = g_DryIce(gl,t,p)
                            !pure hydrates
                            !call hdrt_chem_potent_w(t,p*1.d6,fug_CO2,ChemPot_hyd)
                            call hdrt_chem_potent_w(gl,t,p*1.d6,fug_g,ChemPot_hyd)
                            chem_pot(1) = ChemPot_hyd

                            x_hyd = x_phase(:,phasetype(i))

                            select case (prop)
                            case (4,5) !Internal Energy and Enthalpy
                                call hdrt_enthalpy(gl,t, p*1.D6, x_hyd, chem_pot, fug_g, prop_hyd, errorflag)
                            case (6) !Entropy
                                call hdrt_entropy(gl,t, p*1.D6, x_hyd, fug_g, prop_hyd, errorflag)
                            case (7,8) !Gibbs Energy and Helmholtz energy
                                call hdrt_Gibbs_energy(gl,x_hyd, chem_pot, prop_hyd, errorflag) !Gibbs Energy
                            end select
                        end if

                    else
                        PROP_EOS = -12900d0
                        errorflag = PROP_EOS
                        input=gl%inptorig
                        return
                    end if

                end if

                if (prop == 4) then !internal energy
                    prop_phase(i) = prop_hyd - p * 1.D6 / rho(phasetype(i))
                else if (prop == 8) then !Helmholtz energy
                    prop_phase(i) = prop_hyd - p * 1.D6 / rho(phasetype(i))
                else
                    prop_phase(i) = prop_hyd
                end if
            end if
            !-----------------------------------------------------------------------------------------------------------------------

            !Density might be needed afterwards for further calculations
            if(gl%gecarrier) then
                dens_phase = 1.d0/dens_reac(phasetype(i))
            end if

            dens = dens + phasefrac(phasetype(i)) * dens_phase

            if ((prop /= 1) .and. (prop /= 3)) then ! exception for temperature and pressure
                PROP_EOS = PROP_EOS + phasefrac(phasetype(i)) * prop_phase(i)
            end if

        end do

        dens = 1.d0/dens



        !-----------------------------------------------------------------------------------------------------------------------
        ! exception for non extensive properties that don't follow the addtion of saturated states values
        select case (prop)

        case(2) !Density
            PROP_EOS = dens
            !input=inptorig
            !return

        case (17) !Numerical Pressure Solution (in the VLE Region)
            if (input == 'td') then
                t = prop1
                d = prop2
                if (gl%ncomp >1) then
                    nrsubst = 0
                else
                    nrsubst = 1
                end if
                call reduced_parameters_calc(gl,t)
                PROP_EOS = P_calc(gl,t, d, nrsubst)
                return
            else
                errorflag = -5700
            endif
        case (18) !Numerical Speed of Sound Solution (in the VLE Region)  sqrt[(dp/drho)_s]
            !Moni Sven und Sebastian waren hier und haben auskommentiert
            !BEGIN Kommentar
            !if (input == 'ps') then
            !    call wm_mix_calc(gl,wm_mix)
            !    p1 = prop1
            !    p2 = p1 + 1.d-3
            !    s = prop2
            !    d1 = D_EOS(input, p1, s, fluids, moles, EOS_indicator,path)
            !    d2 = D_EOS(input, p2, s, fluids, moles, EOS_indicator,path)
            !    if ((d2 > d1) .AND. (d1 > 0.d0) .AND. (d2 > 0.d0)) then
            !        dp_drho_s = (p2-p1)*1.d+6/(d2-d1)/wm_mix
            !        PROP_EOS = dp_drho_s**0.5d0
            !    else
            !        PROP_EOS = -2222
            !    end if
            !else
            PROP_EOS = -9955
            errorflag = PROP_EOS
            !end if
            !ENDE Kommentar

        case (19) !residual reduced Helmholtz energy ar
            PROP_EOS = A00_calc(gl,t, dens, nrsubst)
        case (20) !(1st derivative of the residual reduced Helmholtz energy with respect to delta) * delta
            PROP_EOS = A01_calc(gl,t, dens, nrsubst)
        case (21) !( 2nd Derivative of the Residual Reduced Helmholtz Energy with Respect to Delta)*Delta^2
            PROP_EOS = A02_calc(gl,t, dens, nrsubst)
        case (22) !( 3rd Derivative of the Residual Reduced Helmholtz Energy with Respect to Delta)*Delta^3
            PROP_EOS = A03_calc(gl,t, dens, nrsubst)
        case (23) !( 1st Derivative of the Residual Reduced Helmholtz Energy with Respect to Tau)*Tau
            PROP_EOS = A10_calc(gl,t, dens, nrsubst)
        case (24) !( 2nd Derivative of the Residual Reduced Helmholtz Energy with Respect to Delta and Tau)*Delta*Tau
            PROP_EOS = A11_calc(gl,t, dens, nrsubst)
        case (25) !( 3rd Derivative of the Residual Reduced Helmholtz Energy with Respect to Delta^2 and Tau)*Delta^2*Tau
            PROP_EOS = A12_calc(gl,t, dens, nrsubst)
        case (26) !( 2nd Derivative of the Residual Reduced Helmholtz Energy with Respect to Tau)*Tau^2
            PROP_EOS = A20_calc(gl,t, dens, nrsubst)
        case (27) !( 3rd Derivative of the Residual Reduced Helmholtz Energy with Respect to Delta^2 and Tau)*Delta^2*Tau
            PROP_EOS = A21_calc(gl,t, dens, nrsubst)
        case (28) !( 3rd Derivative of the Residual Reduced Helmholtz Energy with Respect to Tau)*Tau^3
            PROP_EOS = A30_calc(gl,t, dens, nrsubst)
        case(51) !VE
            PROP_EOS = 1.d0/dens
        end select


        !Excess Properties VE and HE and GE
        if ((prop == 51) .or. (prop == 52) .or. (prop == 53)) then
            rhoredmixorg = gl%rhoredmix
            tredmixorg = gl%tredmix
            do j = 1, gl%ncomp
                gl%rhoredmix = gl%rhored(j)
                gl%tredmix = gl%tred(j)
                !call PhaseDet_pure_tp(p, t, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
                d = rhomix_calc(gl,t, p, 0.d0, 0, j)
                if (prop == 51) PROP_PURE = 1.d0/d
                if (prop == 52) PROP_PURE = H_CALC(gl,T,d,j)
                if (prop == 53) PROP_PURE = G_CALC(gl,T,d,j)
                PROP_MIX_IDEAL = PROP_MIX_IDEAL + gl%molfractions(j) * PROP_PURE
            end do

            PROP_EOS = PROP_EOS - PROP_MIX_IDEAL
            if (prop == 53) then
                do j = 1, gl%ncomp
                    PROP_EOS = PROP_EOS -gl%Req(j)*t *gl%molfractions(j)*dlog(gl%molfractions(j))
                end do
            end if
            gl%rhoredmix = rhoredmixorg
            gl%tredmix = tredmixorg
        end if
    end if

    if(gl%seawater) then
        if(gl%ncomp .eq. 1) then
            rho(3) = d_sea_calc(gl, t, p,rho(3)) / gl%sea%wm_sea
            if(prop == 2) PROP_EOS = rho(3) !/ gl%sea%wm_sea
        end if
    end if

    if(gl%seawater) then
        call seawater_remap(gl, fluids, moles, gl%sea%dummyarray, gl%sea%dummyarray1, phasetype)
    end if

    if(errorflag .ne. 0) PROP_EOS = errorflag

    input=gl%inptorig
    gl%molfractions = moles


    end function PROP_EOS






    !****************************************************************************************************************************************************
    subroutine inpt_handle(gl,input, prop1, prop2, p, t, d, dvap, dliq, molev, rho, x_Phase, phasetype, phasefrac, nrofphases, nrsubst, errorflag)
    !****************************************************************************************************************************************************
    !New routine to handle all input types. This routine should be called from all interface functions because
    !the handle of the different input codes is for most interface functions the same
    !Andreas June 2014
    !****************************************************************************************************************************************************



    implicit none

    type(type_gl) :: gl


    character (12), intent(in) :: input
    double precision, intent(in) :: prop1, prop2
    double precision, intent(out) :: p, t, d                    !p - pressure, t - temperature, d - density at requested state (Density of phasetype(1) in case of multiple phases)
    double precision, intent(out) :: dvap, dliq
    double precision, intent(inout) , dimension(30) :: molev
    double precision, intent(out) :: rho(5)                     ! vector containing the densities of all phases
    double precision, intent(out) :: x_Phase(30,5)              ! vector containing the compositions of all phases
    integer, dimension(5), intent(out) :: phasetype             !phasetype contains the phase indicator number
    !E.g.: 2 phases are present: liquid and liquid equilibrium
    !-- >  nrofphases = 2
    !-- >  phasetype(1) = 2 (light liquid)
    !-- >  phasetype(2) = 3 (heavy liquid)
    double precision, dimension(5), intent(out) :: phasefrac    ! Phase amounts for all phases
    integer, intent(out) :: nrofphases                          ! at this point: only 1 or 2 implemented
    integer, intent(out) :: nrsubst
    integer, intent(out) :: errorflag


    double precision :: h, s, vapfrac
    integer :: iFlash, iphase, iter, phase

    double precision, dimension(30) :: x_vap, x_liq
    integer :: h_or_s

    p = 0.D0
    t = 0.D0
    d = 0.D0
    h = 0.D0
    s = 0.D0
    dvap = 0.D0
    dliq = 0.D0
    vapfrac = 0.D0
    phasetype = 0
    phasefrac = 0.d0

    x_vap = 0.D0
    x_liq = 0.D0
    errorflag = 0
    nrsubst = 1
    iter=0

    rho = 0.d0

    !call parameters
    ! ---------------------------------------------------------------------
    !  TD-FLASH
    ! ---------------------------------------------------------------------
    if (input == 'td') then
        t = prop1
        d = prop2
        !Hier muss wenn auf die Wasserdichte zurückgerechnet werden!!!!
        !if(gl%seawater) then
        !d = d_w_new(gl, T, P, d)!dw oder pw? d_w new gibt wasserdichte zueück, pw_new den druck!
        !gl%sea%seap = p
        !d = rhomix_calc(gl, t, p, d, 1,1)
        !end if
        if (gl%ncomp == 1) then
            nrsubst = 1     ! pure fluid
            if (gl%check_solid) then
                call PhaseDet_pure_td_sol(gl,p, t, d, rho, phasetype, phasefrac, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
                if (errorflag /= 0) return
                if ((phase == 1) .or. (phase == 3)) then        !   liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif ((phase == 2) .or. (phase == 4)) then    !   vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif (phase == 5) then                        !   Supercrit
                    !Treat supercritical phase as vapor
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif (phase == 0) then                        !   VLE
                    phasetype(1) = 1
                    phasetype(2) = 3
                    x_phase(1,phasetype(2)) = 1.D0
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                    rho(phasetype(1)) = dvap
                    rho(phasetype(2)) = dliq
                else                                            !   Assume solid for all other phaseIDs (except if hold limits was deactivated)
                    if (gl%hold_limits .eqv. .true.) then
                        phasetype(1) = 4
                        phasefrac(phasetype(1)) = 1.D0
                    else
                        phasetype(1) = 1
                        phasefrac(phasetype(1)) = 1.D0
                        rho(phasetype(1)) = d
                    End if
                end if
            end if
            x_phase(1,phasetype(1)) = 1.D0
        else
            nrsubst = 0     ! mixture
            if (gl%check_solid) then
                !Not implemented yet
                errorflag = -12901
            else
                call PhaseDet_td(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, d, nrofphases, errorflag)
                if (errorflag /= 0) return
                if (nrofphases == 1) then
                    phasefrac(phasetype(1)) = 1.D0
                else
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                end if
            end if
        end if


        ! ---------------------------------------------------------------------
        !  TP-FLASH
        ! ---------------------------------------------------------------------
    elseif (input == 'tp') then
        t = prop1
        p = prop2
        if (gl%ncomp == 1) then
            nrsubst = 1     ! pure fluid
            if (gl%check_solid) then
                call PhaseDet_pure_tp_sol(gl,p, t, d, rho, phasetype, phasefrac, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet_pure_tp(gl,p, t, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
                if (errorflag /= 0) return
                if ((phase == 1) .or. (phase == 3)) then        !   liquid phase
                    phasetype(1) = 3
                elseif ((phase == 2) .or. (phase == 4)) then    !   vapor phase
                    phasetype(1) = 1
                elseif (phase == 5) then                        !   Supercrit
                    !Treat supercritical phase as vapor
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                else                                            !   Any other phase, assume solid
                    if (gl%hold_limits .eqv. .true.) then
                        phasetype(1) = 4
                        phasefrac(phasetype(1)) = 1.D0
                    else
                        phasetype(1) = 1
                        phasefrac(phasetype(1)) = 1.D0
                    End if
                end if
                rho(phasetype(1)) = d
            end if
            phasefrac(phasetype(1)) = 1.D0
            x_phase(1,phasetype(1)) = 1.D0
        else
            nrsubst = 0     ! mixture
            if (gl%check_solid) then
                call PhaseDet_sol(gl,p, t, molev, rho, x_Phase, phasetype, Phasefrac, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, nrofphases, errorflag)
                if (errorflag /= 0) return
                if (nrofphases == 1) then
                    phasefrac(phasetype(1)) = 1.D0
                else
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                end if
            end if
        end if


        ! ---------------------------------------------------------------------
        !  PH-FLASH
        ! ---------------------------------------------------------------------
    elseif (input == 'ph') then
        p = prop1
        h = prop2

        if (gl%ncomp == 1) then
            nrsubst = 1     ! pure fluid
            if (gl%check_solid) then
                call PhaseDet_pure_ph_sol(gl,p, t, d, rho, phasetype, phasefrac, h, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet_ph_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, h, nrofphases, errorflag)
                if (errorflag /= 0) return
                if ((phase == 1) .or. (phase == 3)) then        !   liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif ((phase == 2) .or. (phase == 4)) then    !   vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif (phase == 0) then                        !   VLE
                    phasetype(1) = 1
                    phasetype(2) = 3
                    x_phase(1,phasetype(2)) = 1.D0
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                    rho(phasetype(1)) = dvap
                    rho(phasetype(2)) = dliq
                elseif (phase == 5) then                        !   Supercrit
                    !Treat supercritical phase as vapor
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                else                                            ! Any other phase, assume solid
                    if (gl%hold_limits .eqv. .true.) then
                        phasetype(1) = 4
                        phasefrac(phasetype(1)) = 1.D0
                    else
                        phasetype(1) = 1
                        phasefrac(phasetype(1)) = 1.D0
                        rho(phasetype(1)) = d
                    End if
                end if
            end if
            x_phase(1,phasetype(1)) = 1.D0
            if (nrofphases == 2) then
                x_phase(1,phasetype(2)) = 1.D0
            end if
        else
            nrsubst = 0     ! mixture
            if (gl%check_solid) then
                h_or_s = 1  !Enthalpy
                call PhaseDet_ph_ps_sol(gl,p, t, molev, rho, x_Phase, phasetype, Phasefrac, h, h_or_s, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet_ph(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, h, nrofphases, errorflag)
                if (errorflag /= 0) return
                if (nrofphases == 1) then
                    phasefrac(phasetype(1)) = 1.D0
                else
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                end if
            end if
        end if


        ! ---------------------------------------------------------------------
        !  PS-FLASH
        ! ---------------------------------------------------------------------
    elseif (input == 'ps') then
        p = prop1
        s = prop2

        if (gl%ncomp == 1) then
            nrsubst = 1     ! pure fluid
            if (gl%check_solid) then
                call PhaseDet_pure_ps_sol(gl,p, t, d, rho, phasetype, phasefrac, s, nrofphases, errorflag)
                if (errorflag /= 0) return
            else
                call PhaseDet_ps_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, s, nrofphases, errorflag)
                if (errorflag /= 0) return
                if ((phase == 1) .or. (phase == 3)) then        !   liquid phase
                    phasetype(1) = 3
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif ((phase == 2) .or. (phase == 4)) then    !   vapor phase
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                elseif (phase == 0) then                        !   VLE
                    phasetype(1) = 1
                    phasetype(2) = 3
                    x_phase(1,phasetype(2)) = 1.D0
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                    rho(phasetype(1)) = dvap
                    rho(phasetype(2)) = dliq
                elseif (phase == 5) then                        !   Supercrit
                    !Treat supercritical phase as vapor
                    phasetype(1) = 1
                    rho(phasetype(1)) = d
                    phasefrac(phasetype(1)) = 1.D0
                else                                            ! Any other phase, assume solid
                    if (gl%hold_limits .eqv. .true.) then
                        phasetype(1) = 4
                        phasefrac(phasetype(1)) = 1.D0
                    else
                        phasetype(1) = 1
                        phasefrac(phasetype(1)) = 1.D0
                        rho(phasetype(1)) = d
                    End if
                end if
            end if
            x_phase(1,phasetype(1)) = 1.D0
            if (nrofphases == 2) then
                x_phase(1,phasetype(2)) = 1.D0
            end if
        else
            nrsubst = 0     ! mixture
            if (gl%check_solid) then
                h_or_s = 2  !Entropy
                call PhaseDet_ph_ps_sol(gl,p, t, molev, rho, x_Phase, phasetype, Phasefrac, s, h_or_s, nrofphases, errorflag)
            else
                call PhaseDet_ps(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, s, nrofphases, errorflag)
                if (errorflag /= 0) return
                if (nrofphases == 1) then
                    phasefrac(phasetype(1)) = 1.D0
                else
                    phasefrac(phasetype(1)) = vapfrac
                    phasefrac(phasetype(2)) = 1.D0 - vapfrac
                end if
            end if
        end if


        ! ---------------------------------------------------------------------
        !  OTHER INPUT CODES
        ! ---------------------------------------------------------------------
    elseif ((input == 'tliq') .or. (input == 'tvap') .or. (input == 'pliq') .or. (input == 'pvap') .or. (input == 'tq') .or. (input == 'pq')) then
        !VLE
        gl%phaseboundary=.true.
        phasetype(1) = 1
        phasetype(2) = 3
        nrofphases = 2
        if (gl%ncomp == 1) then
            nrsubst = 1     ! pure fluid
            x_liq = molev
            x_vap = molev
            x_phase(1, phasetype(1)) = 1.D0
            x_phase(1, phasetype(2)) = 1.D0
            if (input(1:1) == 'p') then   ! vapor pressure given
                p = prop1
                t = 0.d0; dvap = 0.d0; dliq = 0.d0
                call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 2, errorflag, iter, nrsubst)
                if (errorflag /= 0) return
            else                        ! Sat. Temperature given
                t = prop1
                if (abs(prop1 - gl%tc(nrsubst)) .lt. 1.d-14) then  !critical point
                    dvap = gl%rhoc(nrsubst)
                    dliq = gl%rhoc(nrsubst)
                    p = p_calc(gl,prop1,dvap,nrsubst)
                else
                    p = 0.d0; dvap = 0.d0; dliq = 0.d0
                    call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 1, errorflag, iter, nrsubst)
                    if (errorflag /= 0) return
                end if
            end if
            rho(phasetype(1)) = dvap
            rho(phasetype(2)) = dliq
            if (input(2:4) == 'liq') then
                phasefrac(phasetype(2)) = 1.d0
            else if (input(2:4) == 'vap') then
                phasefrac(phasetype(1)) = 1.d0
            else
                phasefrac(phasetype(1)) = prop2 !q
                phasefrac(phasetype(2)) = 1.d0 - prop2
            end if
        else
            if ( (index(input,'q') /= 0) .and. (len(trim(input)) < 3) ) then   ! this input is not implemented for mixtures
                errorflag = -5250
                return
            end if
            nrsubst = 0     ! mixture
            !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
            !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
            !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
            !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
            iFlash = 0
            if (input == 'pliq') then
                iFlash = 1
                p = prop1
                t = 0.d0
                phasefrac(phasetype(2)) = 1.d0
            end if
            if (input == 'pvap') then
                iFlash = 2
                p = prop1
                t = 0.d0
                phasefrac(phasetype(1)) = 1.d0
            end if
            if (input == 'tliq') then
                iFlash = 3
                t = prop1
                p = 0.d0
                phasefrac(phasetype(2)) = 1.d0
            end if
            if (input == 'tvap') then
                iFlash = 4
                t = prop1
                p = 0.d0
                phasefrac(phasetype(1)) = 1.d0
            end if
            call Flash_PhaseBoundary(gl,p, t, molev, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash, 0, errorflag, iter)
            if (errorflag /= 0) return
            dliq = gl%rho_liq
            dvap = gl%rho_vap
            rho(phasetype(1)) = dvap
            rho(phasetype(2)) = dliq
            x_phase(:, phasetype(1)) = x_vap
            x_phase(:, phasetype(2)) = x_liq
        end if
    end if


    !New inputs: Sublimation and Melting equilibria (for CO2 and water)
    !Andreas Feb 2014
    !--------------------------------------------------------------------------------------
    if (gl%check_solid) then
        !input = PSUBV
        if (input == 'psubv') then
            p = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 1 !p given
                iPhase = 2 !Fluidphase is vapor
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(1))   !Vapor density
                phasefrac(phasetype(1)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                iFlash = 1 !p given
                call Flash_PhaseBoundary_sol (gl,p, t, molev, rho, x_Phase, phasetype, Phasefrac, nrofphases, iFlash, errorflag)
                return
            end if


            !input = PSUBS
        elseif (input == 'psubs') then
            p = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 1 !p given
                iPhase = 2 !Fluidphase is vapor
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(2))   !Solid density
                phasefrac(phasetype(2)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if


            !input = TSUBV
        elseif (input == 'tsubv') then
            t = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 2 !T given
                iPhase = 2 !Fluidphase is vapor
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(1))   !Vapor density
                phasefrac(phasetype(1)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if


            !input = TSUBS
        elseif (input == 'tsubs') then
            t = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 2 !T given
                iPhase = 2 !Fluidphase is vapor
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(2))   !Solid density
                phasefrac(phasetype(2)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if


            !input = PMLTL
        elseif (input == 'pmltl') then
            p = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 1 !p given
                iPhase = 1 !Fluidphase is liquid
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(1))   !liquid density
                phasefrac(phasetype(1)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                iFlash = 1 !p given
                call Flash_PhaseBoundary_sol (gl,p, t, molev, rho, x_Phase, phasetype, Phasefrac, nrofphases, iFlash, errorflag)
                return
            end if


            !input = PMLTS
        elseif (input == 'pmlts') then
            p = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 1 !p given
                iPhase = 1 !Fluidphase is liquid
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(2))   !Solid density
                phasefrac(phasetype(2)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if


            !input = TMLTL
        elseif (input == 'tmltl') then
            t = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 2 !T given
                iPhase = 1 !Fluidphase is liquid
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(1))   !liquid density
                phasefrac(phasetype(1)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if


            !input = TMLTS
        elseif (input == 'tmlts') then
            t = prop1
            gl%phaseboundary=.true.
            if (gl%ncomp == 1) then
                nrsubst = 1
                iFlash = 2 !T given
                iPhase = 1 !Fluidphase is liquid
                call Flash_Pure_PhaseBoundary_sol (gl,p, t, d, rho, phasetype, phasefrac, nrofphases, iFlash, iPhase, errorflag)
                if (errorflag /= 0) return
                d = rho(phasetype(2))   !Solid density
                phasefrac(phasetype(2)) = 1.D0
                x_phase(1, phasetype(1)) = 1.D0
                x_phase(1, phasetype(2)) = 1.D0
            else
                nrsubst = 0
                !Not implemented yet
                errorflag = -9902
                return
            end if
        end if
    end if
    !--------------------------------------------------------------------------------------


    end subroutine




    !****************************************************************************************
    ! Subroutine looks for a Trivial Solution, whether the unknwon is given with the input
    !****************************************************************************************
    subroutine check_triv_input(gl,input, propID, found_prop)

    implicit none

    type(type_gl) :: gl


    character (12), intent(inout) :: input  !TD,TP,PH,PS,PVAP... : Given Input
    character (1), intent(in) :: propID     !T,D,P,H,S: Property that is to be calculated
    integer, intent(out) :: found_prop      !states whether trivial solution was found, 0: no, 1: first property, 2: second property

    character (5) :: lower, upper
    integer :: i, j, noc


    found_prop = 1

    noc=len(trim(input))

    !lower = 'dhpst'
    !upper = 'DHPST'
    !
    !!write lower cases in upper cases for the input variables
    !do i = 1,noc
    !  do j = 1,5
    !    if (input(i:i) == lower(j:j)) then
    !    input(i:i) = upper(j:j)
    !    exit
    !    end if
    !  end do
    !end do


    i = index(input, propID)

    ! extra condition for VAP,SUB...  -- >    < 4
    if ((i == 2).and.(len(trim(input)) < 4)) then
        found_prop=2
    else if ((i > 1).OR.(i == 0)) then
        found_prop=0
    end if

    end subroutine check_triv_input


    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    Double precision function LNG_density(gl,input, prop1, prop2, prop, errorflag)
    !M. Thol

    implicit none

    type(type_gl) :: gl

    character (12), intent(inout) :: input
    double precision, intent(inout) :: prop1, prop2
    integer :: prop, errorflag, i

    if (gl%mix_type == 9) then       ! special case for density calculation with the Costald equation of state
        !M. Thol

        if (((trim(input) /= 'tp') .and. (trim(input) /= 'tliq')) .or. (prop /= 2)) then
            LNG_density = -5778.d0
            errorflag = -5778
            input=gl%inptorig
            return
        end if

        do i = 1, gl%ncomp
            if (gl%eq_type(i) /= 9) then
                LNG_density = -5777.d0
                errorflag = -5777
                input=gl%inptorig
                return
            end if
        end do

        LNG_density = 1.d0/vtaiteq(gl,input, prop1, prop2)
        input=gl%inptorig
        return

    elseif (gl%mix_type == 8) then       ! special case for density calculation with the revised Klosek-McKinley equation of state (eq_type 8)

        ! Tietz Feb 2016

        if ((trim(input) /= 'tliq') .or. (prop /= 2)) then
            LNG_density = -5997.d0
            errorflag = -5997
            input=gl%inptorig
            return
        end if

        do i = 1, gl%ncomp
            if (gl%eq_type(i) /= 8) then
                LNG_density = -5997.d0
                errorflag = -5997
                input=gl%inptorig
                return
            end if
        end do

        LNG_density = D_RKM(gl,prop1, prop2)
        input=gl%inptorig
        return

    elseif (gl%mix_type == 81) then          ! special case for density calculation with the enhanced revised Klosek-McKinley equation of state (eq_type 81)

        ! Tietz Apr 2016

        if (((trim(input) /= 'tp') .and. (trim(input) /= 'tliq')) .or. (prop /= 2)) then
            LNG_density = -5998.d0
            input=gl%inptorig
            return
        end if

        do i = 1, gl%ncomp
            if (gl%eq_type(i) /= 81) then
                LNG_density = -5998.d0
                input=gl%inptorig
                return
            end if
        end do

        !prop2=0.D0

        LNG_density = D_RKM(gl,prop1, prop2)
        return
        input=gl%inptorig

    end if

    End function LNG_density


    end module interface_support_module

    module interface_helper

    contains
    ! HELPER FUNCTIONS ---------------------------------------------------------------------------------------------------------------------------------
    ! Unit definition transformation from character array to integer -> unitin
    integer function get_unitin(unit_char)
    character(*), intent(in) :: unit_char

    if(trim(unit_char) == "specific") then
        get_unitin=2
    elseif(trim(unit_char) == "molar") then
        get_unitin=1
    elseif(trim(unit_char) == "reduced") then
        get_unitin=3
    else
        get_unitin=-123456
    end if

    end function get_unitin
    end module